#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
coalescence_check.py — 좌굴 모드 병합(coalescence) 진단 도구
================================================================================
배경 (arXiv:2609.19603, Prabha & Kumar 2026)
  좌굴 고유치 λ₁·λ₂가 접근하는 설계점(병합점)에서는 "어느 모드로 좌굴되는가"가
  구조가 아니라 불완전성·메시·솔버 경로가 결정한다. 그 결과 포스트버클링(HF) 응답은
  그 지점에서 임의적이 되고, LF-HF 상관 가설도 붕괴할 수 있다.
  이 스크립트는 고유치 상대 간극 g = (λ₂-λ₁)/λ₁ 을 그 위험의 **진단 지표**로 쓴다.
  (논문의 modality 판정은 최적화 기반이지만, 여기서는 진단용 근사만 쓴다. §한계 참조)

두 가지 모드
  1) record  — Abaqus 실행 직후 좌굴 고유치를 읽어 이력(eig_history.jsonl)에 1줄 추가.
               run_abaqus.py 의 Buckle_Analysis 직후 훅에서 자동 호출된다([N-6]).
               표준 라이브러리만 사용 → Abaqus Python 에서도 안전.
  2) analyze — 이력과 HF 결과를 조인해 간극 vs LF-HF 잔차(또는 반복 산포)를 분석.
               Abaqus 불필요, pn40/Windows 어디서나 실행 가능.

사용 예
  # (1) Abaqus 측 기록 (수동)
  abaqus python coalescence_check.py record --dat Buckle_Analysis.dat --x 0.3520 --d 0.01230

  # 콘솔에 찍힌 고유치를 직접 붙여넣기
  python coalescence_check.py record --values 4.3575 4.4677 --x 0.3520 --d 0.01230

  # (2) 회고 분석 (CSV 경로)
  python coalescence_check.py analyze --history code/aba/eig_history.jsonl \
         --hf-csv hf_points.csv --out coalescence_report.csv --plot coalescence.png

  # (2') 회고 분석 (체크포인트 경로 — 정규화 좌표를 mfbo.py 규약으로 역변환)
  python coalescence_check.py analyze --history code/aba/eig_history.jsonl \
         --checkpoint mfbo_checkpoint.pt --plot coalescence.png

  # (3) 동작 확인 (Abaqus/데이터 불필요)
  python coalescence_check.py selftest

HF CSV 형식 (헤더 필수, 열 이름은 자동 인식)
  x_c,d_c,lf,hf          ← lf/hf 는 추력손실(양수) 또는 LF 지표
  0.3520,0.01230,0.0041,0.0097
  반복 실행 산포를 쓰려면 같은 (x_c,d_c) 행을 여러 개 넣으면 된다.
"""

from __future__ import print_function

import argparse
import json
import math
import os
import re
import sys
import time

# ──────────────────────────────────────────────────────────────────────────────
# 공통 유틸
# ──────────────────────────────────────────────────────────────────────────────

HISTORY_DEFAULT = "eig_history.jsonl"

# run_abaqus.py / mfbo.py 규약과 동일해야 한다 (mfbo.py:unnormalize_params)
X1_MIN, X1_MAX = 0.05, 0.95
X2_MIN, X2_MAX = 1e-4, 1.0
KEY_TOL = 6  # (x_c,d_c) 조인 키 반올림 자릿수


def _fnum(tok):
    """Fortran D 지수/유니코드 마이너스까지 허용하는 float 파서."""
    t = tok.strip().replace("D", "E").replace("d", "e")
    t = t.replace("\u2212", "-").replace("\u2013", "-")
    t = t.rstrip(",;")
    return float(t)


def parse_eigenvalues(text):
    """Abaqus .dat / .msg 텍스트에서 좌굴 고유치(하중계수) 목록을 뽑는다.

    지원 형식 (Abaqus 버전별로 헤더가 다름):
        E I G E N V A L U E   O U T P U T
         MODE NO.    EIGENVALUE
             1        4.3575
    또는
        B U C K L I N G   F A C T O R   O U T P U T
         MODE NO.   BUCKLING FACTOR
             1        4.3575
    부분적으로 수렴한 경우 음수/0 모드가 섞여 있을 수 있다(부정정 기저상태) →
    호출측에서 양수만 쓰거나 별도 판단해야 한다. 여기서는 **발견 순서대로 전부** 반환한다.
    """
    if not text:
        return []

    lines = text.splitlines()
    vals = []
    in_table = False
    seen_header = False

    hdr_re = re.compile(r"(e\s*i\s*g\s*e\s*n\s*v\s*a\s*l\s*u\s*e|b\s*u\s*c\s*k\s*l\s*i\s*n\s*g\s+f\s*a\s*c\s*t\s*o\s*r)", re.I)
    row_re = re.compile(r"^\s*(\d{1,6})\s+([-+0-9.][0-9eEdD+\-. ]*)\s*$")

    def _is_banner(line):
        """'E I G E N V A L U E   O U T P U T' 처럼 자간 공백이 섞여도 인식한다."""
        compact = re.sub(r"\s+", "", line.upper())
        return ("EIGENVALUE" in compact) or ("BUCKLINGFACTOR" in compact)

    for raw in lines:
        line = raw.rstrip("\n")
        if (hdr_re.search(line) or _is_banner(line)) and _is_banner(line):
            in_table, seen_header = True, True
            continue
        if not in_table:
            continue
        if not line.strip():
            # 표 사이의 빈 줄은 허용, 두 줄 이상 연속 공백이면 표 종료로 본다
            continue
        if hdr_re.search(line) or re.match(r"^\s*(MODE|NO\.)", line, re.I):
            continue  # 열 제목 줄
        m = row_re.match(line)
        if m:
            try:
                vals.append(_fnum(m.group(2).split()[0]))
            except ValueError:
                pass
            continue
        # 숫자 행이 아니면 표 종료
        if vals:
            break
        if len(line.strip()) > 0 and not re.match(r"^[\s\-_=*]+$", line):
            in_table = False

    if not vals:
        # 폴백: .msg 에 "EIGENVALUE n = x" 혹은 "BUCKLING FACTOR n = x" 형태로 찍히는 경우
        pat = re.compile(r"(?:EIGENVALUE|BUCKLING\s+FACTOR)\s*(?:NO\.?)?\s*\d*\s*[=:]\s*([-+0-9eEdD.]+)", re.I)
        for m in pat.finditer(text):
            try:
                vals.append(_fnum(m.group(1)))
            except ValueError:
                pass

    return vals


def gaps(eigvals):
    """고유치 목록 → (gap_rel, gap_abs). 양수 모드만 사용한다.

    gap_rel = (λ₂-λ₁)/λ₁ : 다음 좌굴 모드가 얼마나 가까운가 (상대).
    부정정 기저상태로 음수 모드가 섞이면 양수만 취해 순서를 재정렬한다.
    """
    pos = sorted([v for v in eigvals if v > 0.0])
    if len(pos) < 2:
        return float("nan"), float("nan")
    l1, l2 = pos[0], pos[1]
    return (l2 - l1) / l1, (l2 - l1)


def unnormalize(nx1, nx2):
    """mfbo.py:unnormalize_params 와 동일한 변환 (x1 선형, x2 지수)."""
    rx1 = X1_MIN + nx1 * (X1_MAX - X1_MIN)
    rx2 = X2_MIN * (X2_MAX / X2_MIN) ** nx2
    return rx1, rx2


def normalize(rx1, rx2):
    """역변환 (체크포인트의 정규화 좌표로 조인할 때 필요)."""
    nx1 = (rx1 - X1_MIN) / (X1_MAX - X1_MIN)
    nx2 = math.log(rx2 / X2_MIN) / math.log(X2_MAX / X2_MIN)
    return nx1, nx2


def key(x_c, d_c):
    """(x_c, d_c) 조인 키. 부동소수 잔차를 흡수하도록 반올림한다."""
    return (round(float(x_c), KEY_TOL), round(float(d_c), KEY_TOL))


# ──────────────────────────────────────────────────────────────────────────────
# record — Abaqus 측 기록기 (표준 라이브러리만)
# ──────────────────────────────────────────────────────────────────────────────

def cmd_record(args):
    vals = []
    src = "manual"

    if args.values:
        vals = [float(v) for v in args.values]
        src = "cli-values"
    else:
        path = args.dat or ""
        if not path:
            print("[record] --dat 또는 --values 중 하나가 필요합니다.", file=sys.stderr)
            return 2
        if not os.path.exists(path):
            if args.msg and os.path.exists(args.msg):
                path = args.msg          # .dat 미생성 환경 대비: .msg 로 폴백
                print("[record] .dat 없음 → .msg 로 폴백: %s" % path)
            else:
                print("[record] 파일 없음: %s (좌굴 해석이 실패했을 수 있음)" % path, file=sys.stderr)
                return 1
        with open(path, "r", errors="replace") as f:
            text = f.read()
        vals = parse_eigenvalues(text)
        src = os.path.basename(path)
        if not vals and args.msg and os.path.exists(args.msg):
            with open(args.msg, "r", errors="replace") as f:
                vals = parse_eigenvalues(f.read())
            src = os.path.basename(args.msg)

    if not vals:
        print("[record] 고유치를 찾지 못했습니다 — Abaqus 버전의 .dat 헤더 형식을 확인하세요.", file=sys.stderr)
        # 실패 원인을 그 자리에서 알 수 있도록 파일 앞부분과 관련 줄을 함께 찍는다.
        #   (a) 좌굴이 실패하면 .dat 에 고유치 표가 아예 없다 -> 여기서 바로 구분된다.
        #   (b) 좌굴은 성공했는데 표 형식이 다르면 아래 'E|' 줄에 그 형식이 보인다.
        try:
            with open(path, 'r', errors='replace') as f:
                lines = f.read().splitlines()
            print("[record] %s 앞 12줄:" % os.path.basename(path))
            for ln in [l for l in lines[:12] if l.strip()]:
                print("    | %s" % ln.strip()[:120])
            hits = [l.strip() for l in lines
                    if ('eigen' in l.lower() or 'buckling factor' in l.lower())]
            for ln in hits[:8]:
                print("    E| %s" % ln[:120])
            if not hits:
                print("    (EIGEN/BUCKLING FACTOR 포함 줄 없음 -> 좌굴이 모드를 못 낸 상태)")
        except Exception as _e:
            print("[record] (진단 출력 실패: %s)" % _e)
        return 1

    g_rel, g_abs = gaps(vals)
    rec = {
        "ts": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "x_c": None if args.x is None else float(args.x),
        "d_c": None if args.d is None else float(args.d),
        "fidelity": (args.fidelity or "").upper() or None,
        "n_modes": len(vals),
        "eigenvalues": [float(v) for v in vals],
        "gap_rel": None if g_rel != g_rel else float(g_rel),
        "gap_abs": None if g_abs != g_abs else float(g_abs),
        "source": src,
        "dat_mtime": None,
    }
    if args.dat and os.path.exists(args.dat):
        rec["dat_mtime"] = os.path.getmtime(args.dat)

    hist = args.history or HISTORY_DEFAULT
    # 멱등성: 같은 (x_c,d_c,dat_mtime) 가 이미 있으면 건너뛴다 (재실행/재호출 대비)
    if os.path.exists(hist):
        try:
            with open(hist, "r", errors="replace") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        old = json.loads(line)
                    except ValueError:
                        continue
                    same_pt = (old.get("x_c") == rec["x_c"] and old.get("d_c") == rec["d_c"])
                    same_src = (old.get("dat_mtime") == rec["dat_mtime"] and rec["dat_mtime"] is not None)
                    if same_pt and (same_src or rec["dat_mtime"] is None):
                        print("[record] 이미 기록됨 — 건너뜀 (x_c=%s, d_c=%s)" % (rec["x_c"], rec["d_c"]))
                        return 0
        except IOError:
            pass

    try:
        with open(hist, "a") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    except IOError as e:
        print("[record] 기록 실패: %s" % e, file=sys.stderr)
        return 1

    _pos = sorted([v for v in vals if v > 0])
    print("[record] %s | n_modes=%d | lambda1=%s lambda2=%s | gap_rel=%.4g"
          % (os.path.basename(hist), rec["n_modes"],
             ("%.6g" % _pos[0]) if _pos else "n/a",
             ("%.6g" % _pos[1]) if len(_pos) > 1 else "n/a", g_rel))
    return 0


# ──────────────────────────────────────────────────────────────────────────────
# analyze — 회고 분석
# ──────────────────────────────────────────────────────────────────────────────

def load_history(paths):
    recs = []
    for p in paths:
        if os.path.isdir(p):
            for root, _dirs, files in os.walk(p):
                for fn in files:
                    if fn.endswith(".jsonl") and ("eig" in fn.lower() or "coales" in fn.lower()):
                        recs.extend(load_history([os.path.join(root, fn)]))
            continue
        if not os.path.exists(p):
            print("[analyze] 이력 파일 없음: %s" % p, file=sys.stderr)
            continue
        with open(p, "r", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    recs.append(json.loads(line))
                except ValueError:
                    pass
    return recs


def load_hf_csv(path, colmap=None):
    """CSV → {(x_c,d_c): [hf, ...]}, {(x_c,d_c): [lf, ...]} (같은 점 반복 실행 허용)."""
    import csv
    hf, lf, rows = {}, {}, []
    with open(path, "r", errors="replace") as f:
        rd = csv.DictReader(f)
        fnames = rd.fieldnames
        if not fnames:
            raise ValueError("CSV 헤더가 없습니다.")
        names = [(n or "").strip().lower() for n in fnames]

        def pick(cands):
            for c in cands:
                if c in names:
                    return fnames[names.index(c)]
            return None
        cx = (colmap or {}).get("x") or pick(["x_c", "x1", "xc", "clamp_x", "x"])
        cd = (colmap or {}).get("d") or pick(["d_c", "x2", "dc", "clamp_d", "d", "disp"])
        chf = (colmap or {}).get("hf") or pick(["hf", "hf_loss", "thrust_loss", "hf_val", "f_hf"])
        clf = (colmap or {}).get("lf") or pick(["lf", "lf1", "lf_val", "compression_area_ratio", "f_lf"])
        if cx is None or cd is None or chf is None:
            raise ValueError("CSV에서 x_c/d_c/hf 열을 찾지 못했습니다. 헤더: %s" % fnames)
        for r in rd:
            try:
                k = key(float(r[cx]), float(r[cd]))
                v = float(r[chf])
            except (TypeError, ValueError):
                continue
            hf.setdefault(k, []).append(v)
            rows.append(k)
            if clf is not None and r.get(clf) not in (None, ""):
                try:
                    lf.setdefault(k, []).append(float(r[clf]))
                except ValueError:
                    pass
    return hf, lf, rows


def load_checkpoint(path, hf_only=True):
    """mfbo_checkpoint.pt → {(x_c,d_c): [y, ...]}, {(x_c,d_c): [lf, ...]}.

    train_x 는 [정규화 x1, 정규화 x2, fidelity(0=LF,1=HF)], train_y 는 부호 반전된 값
    (mfbo.py 에서 최소화를 위해 -값 저장) → 여기서 다시 양수로 되돌린다.
    """
    try:
        import torch
    except ImportError:
        print("[analyze] torch 가 없어 체크포인트를 읽을 수 없습니다. --hf-csv 를 쓰세요.", file=sys.stderr)
        return None, None
    ck = torch.load(path, map_location="cpu", weights_only=False)
    tx = ck["train_x"].detach().cpu().numpy().tolist()
    ty = ck["train_y"].detach().cpu().numpy().ravel().tolist()
    return _ckpt_rows_to_points(zip(tx, ty))


def _ckpt_rows_to_points(rows):
    """[정규화 x1, 정규화 x2, fidelity], 저장값 → ({pt: [hf...]}, {pt: [lf...]}).

    mfbo.py 규약을 그대로 되돌린다:
      - train_x 좌표는 단위 초화면(unit hypercube) → unnormalize_params 로 실좌표 환원
      - train_y 는 '최소화'를 위해 부호가 반전되어 저장됨 → 다시 -를 붙여 원래 값 복원
      - 3번째 열이 1.0 이면 HF, 0.0 이면 LF
    torch 없이도 검증 가능하도록 순수 파이썬만 쓴다 (selftest 에서 호출).
    """
    hf, lf = {}, {}
    for row in rows:
        tx_i, y_i = row
        nx1, nx2, fid = float(tx_i[0]), float(tx_i[1]), float(tx_i[2])
        rx1, rx2 = unnormalize(nx1, nx2)
        k = key(rx1, rx2)
        val = -float(y_i)            # 저장 시 부호 반전 → 복원
        if fid > 0.5:
            hf.setdefault(k, []).append(val)
        else:
            lf.setdefault(k, []).append(val)
    return hf, lf


def spearman(xs, ys):
    """스피어만 순위상관 (동순위 평균 순위). scipy 불필요."""
    n = len(xs)
    if n < 3:
        return float("nan")
    def rank(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for t in range(i, j + 1):
                r[order[t]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)))
    dy = math.sqrt(sum((ry[i] - my) ** 2 for i in range(n)))
    return float("nan") if dx == 0 or dy == 0 else num / (dx * dy)


def _median(v):
    s = sorted(v)
    n = len(s)
    if n == 0:
        return float("nan")
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def cmd_analyze(args):
    hist = load_history(args.history)
    if not hist:
        print("[analyze] 이력이 비어 있습니다. run_abaqus.py 훅 또는 record 모드로 먼저 채우세요.")
        return 1

    # 이력 → {(x_c,d_c): 간극}. 같은 점이 여러 번이면 마지막 기록 사용.
    gap_by_pt, eig_by_pt = {}, {}
    skipped = 0
    for r in hist:
        if r.get("x_c") is None or r.get("d_c") is None:
            skipped += 1
            continue
        k = key(r["x_c"], r["d_c"])
        gap_by_pt[k] = r.get("gap_rel")
        eig_by_pt[k] = r.get("eigenvalues") or []
    if skipped:
        print("[analyze] x_c/d_c 가 없는 기록 %d건은 제외했습니다." % skipped)

    # HF 결과 로드
    hf_by_pt = lf_by_pt = None
    if args.hf_csv:
        hf_by_pt, lf_by_pt, _rows = load_hf_csv(args.hf_csv, args.colmap)
    elif args.checkpoint:
        hf_by_pt, lf_by_pt = load_checkpoint(args.checkpoint)
        if hf_by_pt is None:
            return 1
    else:
        print("[analyze] --hf-csv 또는 --checkpoint 중 하나가 필요합니다.")
        return 2

    # 조인
    rows, unmatched = [], 0
    for k, g in sorted(gap_by_pt.items()):
        if k not in hf_by_pt:
            unmatched += 1
            continue
        hf_vals = hf_by_pt[k]
        lf_vals = (lf_by_pt or {}).get(k, [])
        hf_mean = sum(hf_vals) / len(hf_vals)
        spread = (max(hf_vals) - min(hf_vals)) if len(hf_vals) > 1 else None
        resid = None
        if lf_vals:
            lf_mean = sum(lf_vals) / len(lf_vals)
            resid = abs(lf_mean - hf_mean)
        rows.append({
            "x_c": k[0], "d_c": k[1], "gap_rel": g,
            "n_pos_modes": sum(1 for v in (eig_by_pt.get(k) or []) if v > 0),
            "lam1": (sorted([v for v in (eig_by_pt.get(k) or []) if v > 0]) + [None, None])[0],
            "lam2": (sorted([v for v in (eig_by_pt.get(k) or []) if v > 0]) + [None, None, None])[1],
            "hf_mean": hf_mean, "hf_repeat_spread": spread,
            "lf_mean": (sum(lf_vals) / len(lf_vals)) if lf_vals else None,
            "lf_hf_resid": resid,
            "in_band": (g is not None and g <= args.gap_thresh),
        })

    if not rows:
        print("[analyze] 조인된 점이 없습니다 (이력 %d점, HF %d점). "
              "x_c/d_c 좌표계가 같은지 확인하세요." % (len(gap_by_pt), len(hf_by_pt)))
        return 1

    # 지표: 같은 설계점 반복 실행 산포가 있으면 그것을, 없으면 |LF-HF| 잔차를 대용으로 쓴다
    if any(r["hf_repeat_spread"] is not None for r in rows):
        metric_name = "hf_repeat_spread"
    elif any(r["lf_hf_resid"] is not None for r in rows):
        metric_name = "lf_hf_resid"
    else:
        metric_name = None
    metric = [(r.get(metric_name) if metric_name else None) for r in rows]
    metric = [(m if m is not None else float("nan")) for m in metric]

    use = [(r["gap_rel"], m) for r, m in zip(rows, metric) if r["gap_rel"] is not None and m == m]
    inb = [m for r, m in zip(rows, metric) if r["in_band"] and m == m]
    outb = [m for r, m in zip(rows, metric) if (not r["in_band"]) and m == m]

    # CSV 출력
    if args.out:
        import csv
        cols = ["x_c", "d_c", "lam1", "lam2", "gap_rel", "in_band", "hf_mean",
                "hf_repeat_spread", "lf_mean", "lf_hf_resid"]
        with open(args.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            for r in rows:
                w.writerow({c: r.get(c) for c in cols})
        print("[analyze] CSV 저장: %s (%d행)" % (args.out, len(rows)))

    # 요약 출력
    n_in = sum(1 for r in rows if r["in_band"])
    print("")
    print("=" * 74)
    print(" 좌굴 모드 병합(coalescence) 진단 요약")
    print("=" * 74)
    print(" 기록 점수            : %d (조인 %d, 미조인 %d)" % (len(hist), len(rows), unmatched))
    print(" 임계값 gap_thresh    : %.3g  → 밴드 내 %d점 / 밖 %d점" % (args.gap_thresh, n_in, len(rows) - n_in))
    if metric_name:
        print(" 신뢰도 대용 지표      : %s%s" % (metric_name,
              "  (같은 점 반복 실행 산포)" if metric_name == "hf_repeat_spread"
              else "  (|LF-HF| 잔차 — 반복 실행 데이터가 없을 때의 대용)"))
        if inb:
            print("   밴드 내 중앙값      : %.6g  (n=%d)" % (_median(inb), len(inb)))
        if outb:
            print("   밴드 밖 중앙값      : %.6g  (n=%d)" % (_median(outb), len(outb)))
        if inb and outb:
            ratio = _median(inb) / _median(outb) if _median(outb) else float("inf")
            print("   비율 (내/외)        : %.3g  → %s" % (ratio,
                  "가설 지지(밴드 내가 더 불안정)" if ratio >= args.ratio_thresh
                  else "가설 미지지 — 밴드 내가 특별히 나쁘지 않음"))
        rho = spearman([g for g, _m in use], [_m for _g, _m in use])
        print("   스피어만 상관 ρ(gap, 지표): %.4f  (n=%d, gap이 작을수록 지표가 커야 하므로 음수 기대)"
              % (rho, len(use)))
    else:
        print(" 신뢰도 대용 지표      : 없음 (CSV에 lf 열이 없고 반복 실행도 없음)")
        print("   → LF 열을 추가하거나, 같은 설계점을 imperfection seed/메시를 바꿔 2회 이상 돌려")
        print("     같은 (x_c,d_c) 행을 여러 개 넣어 주세요.")

    if len(rows) < args.min_points:
        print("")
        print(" ⚠ 표본 %d점 < 권장 %d점 — 이 결과로 가설을 채택/기각하지 마세요 (탐색적 지표)."
              % (len(rows), args.min_points))
    print("=" * 74)

    if args.plot:
        try:
            make_plot(rows, metric, metric_name, args.gap_thresh, args.plot)
            print("[analyze] 그림 저장: %s" % args.plot)
        except Exception as e:      # matplotlib 없음 등 — 분석 결과에는 영향 없음
            print("[analyze] 그림 생략: %s" % e)

    return 0


def make_plot(rows, metric, metric_name, thresh, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    g = [r["gap_rel"] for r in rows if r["gap_rel"] is not None]
    ax[0].hist(g, bins=min(20, max(5, len(g) // 3)), color="0.7", edgecolor="0.3")
    ax[0].axvline(thresh, color="crimson", ls="--", lw=1.2, label="gap_thresh=%.3g" % thresh)
    ax[0].set_xlabel(r"relative eigenvalue gap  $(\lambda_2-\lambda_1)/\lambda_1$")
    ax[0].set_ylabel("count")
    ax[0].set_title("gap distribution")
    ax[0].legend(fontsize=8)
    if g and max(g) > 4 * thresh:
        ax[0].set_xlim(left=0, right=min(max(g), max(10 * thresh, 0.05)))

    xs = [r["gap_rel"] for r in rows if metric_name and r.get("gap_rel") is not None
          and r.get(metric_name) is not None]
    ys = [r.get(metric_name) for r in rows if metric_name and r.get("gap_rel") is not None
          and r.get(metric_name) is not None]
    if xs:
        ax[1].scatter(xs, ys, s=18, color="0.25")
        ax[1].axvline(thresh, color="crimson", ls="--", lw=1.2)
        ax[1].set_xlabel(r"gap  $(\lambda_2-\lambda_1)/\lambda_1$")
        ax[1].set_ylabel(metric_name)
        ax[1].set_title("gap vs HF reliability proxy")
    else:
        ax[1].text(0.5, 0.5, "no reliability proxy available", ha="center", va="center")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────────
# selftest — 실제 Abaqus 산출물 없이 동작 확인
# ──────────────────────────────────────────────────────────────────────────────

SELFTEST_DAT = """\
                            A B A Q U S
   STEP     1  BUCKLE ANALYSIS

                  B U C K L I N G   F A C T O R   O U T P U T

 MODE NO.   BUCKLING FACTOR
        1       4.3575
        2       4.4677
        3       5.1230
        4       5.9001

   T H E   A N A L Y S I S   H A S   B E E N   C O M P L E T E D
"""

SELFTEST_DAT_NEG = """\
                  E I G E N V A L U E    O U T P U T

 MODE NO.    EIGENVALUE
     1      -2.10000
     2      -0.45000
     3       3.75000
     4       3.90000
     5       4.50000
"""


def cmd_selftest(args):
    ok = True
    print("=== selftest: 파서 ===")
    v = parse_eigenvalues(SELFTEST_DAT)
    print("  정상 .dat → %s" % v)
    ok &= (v == [4.3575, 4.4677, 5.1230, 5.9001])
    v2 = parse_eigenvalues(SELFTEST_DAT_NEG)
    print("  음수모드 .dat → %s" % v2)
    ok &= (v2 == [-2.1, -0.45, 3.75, 3.9, 4.5])
    g, ga = gaps(v2)
    print("  양수만 사용 gap_rel=%.6g gap_abs=%.6g (3.75→3.90)" % (g, ga))
    ok &= abs(g - (3.90 - 3.75) / 3.75) < 1e-12
    v3 = parse_eigenvalues("EIGENVALUE  1 = 1.234E+00\nBUCKLING FACTOR 2 = 1.500D+00")
    print("  .msg 폴백 → %s" % v3)
    ok &= (len(v3) == 2)
    nx1, nx2 = normalize(0.5, 0.01)
    rx1, rx2 = unnormalize(nx1, nx2)
    print("  좌표 왕복: x1=%.10f x2=%.10f" % (rx1, rx2))
    ok &= (abs(rx1 - 0.5) < 1e-12 and abs(rx2 - 0.01) < 1e-12)

    print("=== selftest: 체크포인트 규약(부호·정규화) ===")
    hf_pts, lf_pts = _ckpt_rows_to_points([
        ([0.0, 0.0, 1.0], -0.5),     # HF, 정규화 최소 → 실좌표 (0.05, 1e-4), 값 +0.5
        ([1.0, 1.0, 0.0], -0.25),    # LF, 정규화 최대 → 실좌표 (0.95, 1.0), 값 +0.25
        ([0.5, 0.5, 1.0], -0.5),     # HF 중앙
    ])
    hk_lo = key(0.05, 1e-4)
    hk_hi = key(0.95, 1.0)
    print("  HF 키 %s → %s / LF 키 %s → %s" % (hk_lo, hf_pts.get(hk_lo), hk_hi, lf_pts.get(hk_hi)))
    ok &= (abs(hf_pts.get(hk_lo, [None])[0] - 0.5) < 1e-12)
    ok &= (abs(lf_pts.get(hk_hi, [None])[0] - 0.25) < 1e-12)
    ok &= (hk_hi not in hf_pts)

    print("=== selftest: record → analyze 왕복 (기능 검증) ===")
    import csv as _csv
    import random
    import tempfile

    tmp = tempfile.mkdtemp(prefix="coalescence_selftest_")
    hist = os.path.join(tmp, "eig_history.jsonl")
    csvp = os.path.join(tmp, "hf.csv")
    outcsv = os.path.join(tmp, "report.csv")

    def _ns(**kw):
        class _O(object):
            pass
        o = _O()
        for k, v in kw.items():
            setattr(o, k, v)
        return o

    # 두 설계점: A = 간극 큼(위험 밖), B = 간극 작음(위험 밴드)
    pts = {
        "A": dict(x=0.3520, d=0.01230, eig=[4.3575, 4.4677, 5.1230]),   # gap_rel ~ 2.5e-2
        "B": dict(x=0.8000, d=0.50000, eig=[4.3575, 4.3620, 4.9000]),   # gap_rel ~ 1.0e-3
    }
    for tag, p in pts.items():
        rc = cmd_record(_ns(values=[str(v) for v in p["eig"]], dat=None, msg=None,
                            history=hist, x=p["x"], d=p["d"], fidelity="HF"))
        ok &= (rc == 0)
    n_lines = sum(1 for _ in open(hist))
    print("  기록 줄 수 = %d (2 기대)" % n_lines)
    ok &= (n_lines == 2)

    # 멱등성: 같은 점 재기록 → 줄 수 그대로
    rc = cmd_record(_ns(values=[str(v) for v in pts["A"]["eig"]], dat=None, msg=None,
                        history=hist, x=pts["A"]["x"], d=pts["A"]["d"], fidelity="HF"))
    ok &= (rc == 0 and sum(1 for _ in open(hist)) == 2)
    print("  멱등성 확인 (재기록 후에도 2줄)")

    # 합성 HF: 밴드 내(B)의 반복 산포가 밴드 밖(A)보다 크도록 (가설 지지 데이터)
    random.seed(7)
    with open(csvp, "w") as f:
        f.write("x_c,d_c,lf,hf\n")
        for tag, p in pts.items():
            noise = 0.5 if tag == "A" else 3.0
            for _ in range(4):
                base = 0.010
                f.write("%.6f,%.6f,%.6f,%.6f\n"
                        % (p["x"], p["d"], base, base + noise * random.uniform(0.5, 1.5)))

    rc = cmd_analyze(_ns(history=[hist], hf_csv=csvp, checkpoint=None, colmap=None,
                         out=outcsv, plot=None, gap_thresh=0.01, ratio_thresh=1.5, min_points=2))
    ok &= (rc == 0)

    # 결과 CSV 검증: in_band 플래그와 산포 순서
    got = {}
    with open(outcsv) as f:
        for r in _csv.DictReader(f):
            got[round(float(r["x_c"]), 4)] = (r["in_band"], float(r["hf_repeat_spread"]))
    a = got.get(round(pts["A"]["x"], 4))
    b = got.get(round(pts["B"]["x"], 4))
    print("  A(간극 큼): in_band=%s spread=%.4f" % a if a else "  A: 조인 실패")
    print("  B(간극 작음): in_band=%s spread=%.4f" % b if b else "  B: 조인 실패")
    ok &= bool(a and b)
    if a and b:
        ok &= (a[0] == "False" and b[0] == "True")      # 밴드 판정
        ok &= (b[1] > a[1])                              # 밴드 내 산포가 더 큼
        print("  판정·순서 검증: %s" % ("OK" if (a[0] == "False" and b[0] == "True" and b[1] > a[1]) else "FAIL"))

    print("")
    print("SELFTEST: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="좌굴 모드 병합(coalescence) 진단 — 고유치 간극 기록/분석",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="예: python coalescence_check.py analyze --history code/aba/eig_history.jsonl --hf-csv hf.csv --out report.csv --plot gap.png")
    sub = ap.add_subparsers(dest="cmd")

    r = sub.add_parser("record", help="좌굴 고유치를 이력(JSONL)에 추가 (Abaqus 측)")
    r.add_argument("--dat", help="Buckle_Analysis.dat 경로")
    r.add_argument("--msg", help="폴백으로 읽을 .msg 경로 (선택)")
    r.add_argument("--values", nargs="*", help="고유치를 직접 지정 (예: --values 4.3575 4.4677)")
    r.add_argument("--x", type=float, help="설계변수 x_c (실수, run_abaqus.py 인자와 동일)")
    r.add_argument("--d", type=float, help="설계변수 d_c (실수)")
    r.add_argument("--fidelity", default="HF", help="HF/LF 표기 (기본 HF)")
    r.add_argument("--history", default=None, help="이력 파일 (기본 ./%s)" % HISTORY_DEFAULT)
    r.set_defaults(func=cmd_record)

    a = sub.add_parser("analyze", help="이력 + HF 결과 조인 분석 (회고)")
    a.add_argument("--history", nargs="+", required=True, help="eig_history.jsonl 파일 또는 디렉터리")
    a.add_argument("--hf-csv", default=None, help="HF 결과 CSV (x_c,d_c,hf[,lf])")
    a.add_argument("--checkpoint", default=None, help="mfbo_checkpoint.pt (대안 입력)")
    a.add_argument("--colmap", default=None, help="열 이름 강제 지정 JSON (예: '{\"hf\":\"thrust\"}')")
    a.add_argument("--gap-thresh", type=float, default=1e-2, help="병합 위험 구간 임계값 (기본 1e-2)")
    a.add_argument("--ratio-thresh", type=float, default=1.5, help="밴드 내/외 중앙값 비율 채택 기준 (기본 1.5)")
    a.add_argument("--min-points", type=int, default=8, help="권장 최소 표본 수 (기본 8)")
    a.add_argument("--out", default=None, help="결과 CSV 경로")
    a.add_argument("--plot", default=None, help="그림 PNG 경로 (matplotlib 필요)")
    a.set_defaults(func=cmd_analyze)

    s = sub.add_parser("selftest", help="Abaqus 없이 파서·기록·분석 동작 확인")
    s.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    if not getattr(args, "cmd", None):
        ap.print_help()
        return 0
    if args.cmd == "analyze" and args.colmap:
        args.colmap = json.loads(args.colmap)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
