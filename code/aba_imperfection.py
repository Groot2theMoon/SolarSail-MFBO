#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""임퍼펙션(좌굴모드) 소스 스테이징 — 2-모델 레시피 (2026-09-22)

무엇을 하는가
  좌굴모드를 '클램프 없는' 모델에서 추출해, '클램프가 있는' 포스트버클링 모델의
  기하학적 임퍼펙션(*IMPERFECTION)으로 주입한다(Galhofo 2022 와 같은 2-모델 구조).

    모드 추출 : run_abaqus_cable.py   (클램프 없음, 실측상 항상 성공, 모드 4개)
    모드 소비 : run_abaqus.py -- HF   (클램프 포함 포스트버클링)

  이 모듈은 소스 .fil 을 HF 잡 디렉터리로 '스테이징'하고, 그 전에 모드 개수를 검증한다.

왜 이름을 바꿔 복사하는가 (원장 C-1)
  run_abaqus.py 는 '자기' 좌굴 잡도 Buckle_Analysis 라는 이름으로 제출한다(클램프 모델,
  실측상 음수 고유값 598~2897 개 / CONVERGED=0 -> 모드 0개). 같은 이름을 소비하면
  0-모드 .fil 이 4-모드 .fil 을 덮어써 임퍼펙션이 '조용히' 사라진다.
  -> 별도 이름(IMPERFECTION_NAME)으로 복사해 소비하면 이름 충돌이 원천 차단된다.
  -> 스테이징 단계에서 모드 개수를 검증하므로, 실패를 조용히 넘기지 않는다.

왜 노드 매핑이 성립하는가 (라이선스 0 정적 검증, 2026-09-22)
  create_rigid_patch 는 face partition 없이 getByBoundingSphere 로 '기존 노드'를 골라
  Coupling(KINEMATIC) 만 건다 -> 클램프가 메쉬를 바꾸지 않는다.
  막 메쉬 결정 인자(형상 BASE/HEIGHT/THICKNESS, seedPart(BASE/200.0, deviationFactor=0.1),
  QUAD_DOMINATED/FREE/MEDIAL_AXIS, 막 요소코드 S4)가 케이블 런과 HF 에서 동일하다
  -> 노드 좌표/라벨이 동일 -> *IMPERFECTION 의 노드 라벨 매핑이 성립한다.

abaqus 모듈을 쓰지 않는다(stdlib only) -> 하네스에서 단위시험 가능.
"""

from __future__ import print_function

import os
import re
import shutil
import time

try:
    # 기존 검증된 고유치(.dat MODE NO 표) 파서 재사용 — 중복 구현 금지
    from coalescence_check import parse_eigenvalues
except Exception:                                        # pragma: no cover
    parse_eigenvalues = None

_ROW_RE = re.compile(r"^\s*(\d{1,6})\s+([-+0-9.][0-9eEdD+\-. ]*)\s*$")


class ImperfectionSourceError(RuntimeError):
    """임퍼펙션 소스(.fil)를 스테이징할 수 없다. HF 잡 제출 '전에' 발생해야 한다."""


def _count_rows_after_header(text):
    """parse_eigenvalues 를 쓸 수 없을 때의 최소 폴백: MODE NO 표의 행 수."""
    n, seen = 0, False
    for line in text.splitlines():
        up = line.upper()
        if "MODE NO" in up or "BUCKLING FACTOR" in up:
            seen = True
            continue
        if not seen:
            continue
        if _ROW_RE.match(line):
            n += 1
        elif line.strip() and not re.match(r"^[\s\-_=*]+$", line):
            break
    return n


def _count_from_msg(text):
    """'.msg' 의 CONVERGED / REQUESTED 줄에서 모드 개수 추정 (최대값)."""
    best = 0
    for line in text.splitlines():
        up = line.upper()
        if ("CONVERGED UPTO THIS ITERATION" in up) or ("REQUESTED BY THE USER" in up):
            tail = up.split("=")[-1] if "=" in up else up.split(":")[-1]
            tok = tail.strip().split()[0] if tail.strip() else ""
            if tok.isdigit():
                best = max(best, int(tok))
    return best


def count_modes(dat_path=None, msg_path=None):
    """소스 좌굴 런의 모드 개수. (.dat MODE NO 표 우선, 없으면 .msg 폴백) -> int."""
    if dat_path and os.path.exists(dat_path):
        with open(dat_path, "r", errors="replace") as f:
            text = f.read()
        if parse_eigenvalues is not None:
            vals = parse_eigenvalues(text)
            if vals:
                return len(vals)
        n = _count_rows_after_header(text)
        if n:
            return n
    if msg_path and os.path.exists(msg_path):
        with open(msg_path, "r", errors="replace") as f:
            return _count_from_msg(f.read())
    return 0


# =====================================================================
# ODB 모드표 -> 기하 섭동 (2026-09-22)
#   `abaqus python aba_mode_from_odb.py <src.odb> <고유스텝> <out.txt> [n]` 이 만든 표를 읽는다.
#   표의 u3 은 모드별 max|u3|=1 정규화 + 부호 고정이므로, 진폭만 곱해 더하면
#   *IMPERFECTION 의 'mode, scale' 합 규약과 같은 기하 섭동이 된다(셸에서는 노드 좌표 섭동 = 임퍼펙션).
#   .fil 경로와 달리 *BUCKLE 스텝의 파일출력 금지에 걸리지 않는다(원본 docstring 의 계획 경로).
# =====================================================================
MODE_TABLE_HEADER = 'MODE'


def load_mode_table(path, modes=(1, 2, 3, 4)):
    """모드표를 {mode: {node_label: u3(정규화)}} 로 읽는다. 문제가 있으면 ImperfectionSourceError."""
    if not os.path.exists(path):
        raise ImperfectionSourceError(
            "%s 를 찾지 못했습니다.\n  먼저 모드표를 만드세요:\n"
            "    abaqus python aba_mode_from_odb.py <소스.odb> <고유스텝> \"%s\" %d"
            % (path, path, len(modes)))
    want = set(int(m) for m in modes)
    table = {}
    cur = None
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            if s.upper().startswith(MODE_TABLE_HEADER):
                parts = s.split()
                if len(parts) < 2 or not parts[1].isdigit():
                    raise ImperfectionSourceError("모드표 %s: 잘못된 MODE 행 %r" % (path, s))
                cur = int(parts[1])
                if cur in want:
                    table[cur] = {}
                continue
            if cur is None or cur not in want:
                continue
            parts = s.split()
            if len(parts) < 2:
                raise ImperfectionSourceError("모드표 %s: 잘못된 행 %r" % (path, s))
            table[cur][int(parts[0])] = float(parts[1])
    missing = sorted(want - set(table.keys()))
    if missing:
        raise ImperfectionSourceError("모드표 %s 에 모드 %s 가 없습니다 (있는 것: %s)"
                                      % (path, missing, sorted(table.keys())))
    for m in sorted(table):
        if not table[m]:
            raise ImperfectionSourceError("모드표 %s: 모드 %d 가 비어 있습니다" % (path, m))
    return table


def mode_table_info(path):
    """모드표 헤더(# key=value)를 dict 로 돌려준다(추적용: source/step/instance).

    출처를 상수로 중복 기재하지 않고 표에서 읽어 로그에 남긴다(교란 차단: 엉뚱한 ODB 를 쓰면 로그가 말해준다).
    """
    info = {}
    if not os.path.exists(path):
        return info
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s.startswith('#'):
                break
            body = s.lstrip('#').strip()
            if '=' in body:
                k, v = body.split('=', 1)
                info[k.strip()] = v.strip()
    return info


def build_perturbation(table, amplitude, modes=(1, 2, 3, 4)):
    """{label: Δz} = Σ_modes amplitude * φ_mode(label). (*IMPERFECTION 의 합 규약과 동일)"""
    out = {}
    for m in modes:
        for lab, u3 in table[int(m)].items():
            out[lab] = out.get(lab, 0.0) + amplitude * u3
    return out


def perturbation_report(pert):
    """로그용 한 줄 요약(섭동의 실현 크기를 정직하게 남긴다)."""
    if not pert:
        return "빈 섭동"
    vals = [abs(v) for v in pert.values()]
    lab_max = max(pert.items(), key=lambda kv: abs(kv[1]))[0]
    return ("노드 %d개, max|dz|=%.3e m (노드 %s), 평균|dz|=%.3e m"
            % (len(pert), max(vals), lab_max, sum(vals) / len(vals)))


def resolve_source(hint, run_dir=None, here=None):
    """소스 .fil 의 실제 경로를 찾는다. 못 찾으면 후보 목록과 함께 예외."""
    cands = []
    if hint:
        cands.append(hint if os.path.isabs(hint) else os.path.join(run_dir or os.getcwd(), hint))
        cands.append(hint)
        if here:
            cands.append(os.path.join(here, hint))
            cands.append(os.path.join(here, os.path.basename(hint)))
    seen, uniq = set(), []
    for c in cands:
        c = os.path.abspath(c)
        if c not in seen:
            seen.add(c)
            uniq.append(c)
    for c in uniq:
        if os.path.exists(c) and os.path.getsize(c) > 0:
            return c
    raise ImperfectionSourceError(
        "임퍼펙션 소스 .fil 을 찾을 수 없습니다.\n"
        "  찾아본 경로: %s\n"
        "  -> 먼저 클램프 없는 모드 추출 런을 실행하세요: abaqus cae noGUI=run_abaqus_cable.py"
        % "\n               ".join(uniq))


def stage(source_hint, target_name, run_dir, required_modes,
          dat_hint=None, msg_hint=None, here=None):
    """소스 .fil 을 <run_dir>/<target_name>.fil 로 스테이징하고 검증 결과를 dict 로 반환.

    required_modes : 주입할 모드 번호 목록(예: (1,2,3,4)). max() 개수 이상이어야 통과.
    """
    src = resolve_source(source_hint, run_dir, here)
    sib = os.path.dirname(src)
    dat = dat_hint if (dat_hint and os.path.exists(dat_hint)) else (
        os.path.join(sib, os.path.basename(dat_hint)) if dat_hint else None)
    msg = msg_hint if (msg_hint and os.path.exists(msg_hint)) else (
        os.path.join(sib, os.path.basename(msg_hint)) if msg_hint else None)

    n = count_modes(dat, msg)
    need = max(required_modes) if required_modes else 0
    if n < need:
        raise ImperfectionSourceError(
            "소스 .fil 에 모드가 부족합니다: %d개 (필요 %d개)\n"
            "  source = %s\n  dat    = %s\n  msg    = %s\n"
            "  -> 좌굴 런이 CONVERGED=0 으로 끝났을 가능성이 큽니다"
            "(클램프가 있는 base state 에서 뽑으려 했는지 확인)."
            % (n, need, src, dat, msg))

    dst = os.path.join(run_dir, target_name + ".fil")
    copied = os.path.abspath(src) != os.path.abspath(dst)
    if copied:
        shutil.copyfile(src, dst)
    return {
        "source": src, "target": dst, "name": target_name, "copied": copied,
        "modes": n, "src_size": os.path.getsize(src),
        "src_mtime": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(src))),
        "dat": dat, "msg": msg,
    }


def imperfection_text(name, step_no, modes, amplitude_m):
    """*IMPERFECTION 키워드 텍스트. 각 줄은 '모드번호, 진폭[길이단위]' (Abaqus 규약)."""
    lines = ["*IMPERFECTION, FILE=%s, STEP=%d" % (name, int(step_no))]
    for m in modes:
        lines.append("%d, %.6e" % (int(m), float(amplitude_m)))
    return "\n".join(lines)


def report(staged):
    """로그 한 줄용 요약 문자열."""
    return ("source=%s (%.1f KB, mtime=%s, 모드 %d개) -> target=%s [%s]"
            % (staged["source"], staged["src_size"] / 1024.0, staged["src_mtime"],
               staged["modes"], staged["target"],
               "복사됨" if staged["copied"] else "제자리"))
