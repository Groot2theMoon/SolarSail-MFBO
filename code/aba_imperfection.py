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
  막 메쉬 결정 인자(형상 BASE/HEIGHT/THICKNESS, seedPart(BASE/SEED_DIV, deviationFactor=0.05),
  QUAD_DOMINATED/FREE/ADVANCING_FRONT, 막 요소코드 S4)가 케이블 런과 HF 에서 동일하다
  [2026-10-02] MEDIAL_AXIS -> ADVANCING_FRONT, deviationFactor 0.1 -> 0.05 로 변경(메쉬 품질 시험).
  노드 라벨이 바뀌므로 임퍼펙션 모드표와 모드 소스 ODB 를 반드시 함께 재생성해야 한다.
  -> 노드 좌표/라벨이 동일 -> *IMPERFECTION 의 노드 라벨 매핑이 성립한다.

abaqus 모듈을 쓰지 않는다(stdlib only) -> 하네스에서 단위시험 가능.
"""

from __future__ import print_function

import hashlib
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

# ---- 모드표 지문 (2026-09-24) ----------------------------------------------------
# 좌굴모드는 '클램프 없는' 모델에서 한 번만 계산해 두고 HF 가 재사용한다. 그런데 모드 소스의
# 상수(프리텐션/패치/안정화 등)를 바꾼 뒤 다시 뽑지 않으면 HF 는 '다른 형상에서 나온 모드'를
# 조용히 주입하게 된다. 그래서 표 헤더에 모드 소스의 지문을 박아 두고 HF 가 대조한다.
FINGERPRINT_PREFIX = '# FINGERPRINT'
FINGERPRINT_KEYS = ('BASE', 'HEIGHT', 'THK', 'SIGMA0', 'PRETENSION_MODE', 'DISP_GLOBAL',
                    'PATCH_RADIUS', 'COUPLING_TYPE', 'MODE_STABILIZATION', 'PERTURBATION', 'PATTERN_SIGN',
                    'N_EIG_BUCKLE', 'BUCKLE_VECTORS', 'N_MODES', 'JOB_NAME', 'MODE_STEP_NAME',
                    'ELEM_CODE_QUAD', 'ELEM_CODE_TRI', 'SEED_DIV')


def model_fingerprint(source_path, keys=FINGERPRINT_KEYS):
    """모드 소스 스크립트에서 모델/모드 정의 상수만 뽑아 sha1 지문(16자)을 만든다.

    반환: (지문 or None, {상수: 값}). 소스를 못 읽거나 상수를 하나도 못 찾으면 (None, {}).
    """
    vals = []
    try:
        # [2026-10-05] encoding 을 반드시 명시한다. 이 파일에는 한글 주석이 있고, Windows 의
        #   기본 인코딩(cp949)으로 UTF-8 파일을 열면 UnicodeDecodeError 가 나는데, 아래
        #   `except Exception: return None, {}` 가 그것을 **조용히 삼켜** 지문이 항상 사라졌다
        #   (실측: 모드 소스 로그에 '[MODE] 재사용 지문 없음(생략)').
        #   지문은 '모드표가 다른 모델에서 만들어졌다'를 잡는 유일한 장치이므로 조용히 꺼지면 안 된다.
        try:
            fh = open(source_path, 'r', encoding='utf-8', errors='replace')
        except TypeError:                      # encoding 인자를 모르는 파이썬 폴백
            fh = open(source_path, 'r', encoding='utf-8', errors='replace')
        with fh as f:
            for line in f.read().splitlines():
                code = line.split('#')[0].strip()
                for k in keys:
                    m = re.match(r'^%s\s*=\s*(.+)$' % re.escape(k), code)
                    if m:
                        vals.append('%s=%s' % (k, m.group(1).strip()))
    except Exception:
        return None, {}
    if not vals:
        return None, {}
    digest = hashlib.sha1('|'.join(vals).encode('utf-8')).hexdigest()[:16]
    return digest, dict(v.split('=', 1) for v in vals)


def model_fingerprint_line(source_path):
    """모드표 헤더에 넣을 '# FINGERPRINT sha1=... src=...' 줄. 못 만들면 None."""
    digest, _vals = model_fingerprint(source_path)
    if not digest:
        return None
    return '%s sha1=%s src=%s' % (FINGERPRINT_PREFIX, digest, os.path.basename(source_path))


def table_fingerprint(table_path):
    """모드표 헤더에서 지문을 읽는다. 없으면 None(모드 데이터 행을 만나면 중단)."""
    try:
        with open(table_path, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                s = line.strip()
                if s.startswith(FINGERPRINT_PREFIX):
                    for tok in s.split():
                        if tok.startswith('sha1='):
                            return tok.split('=', 1)[1]
                if s.startswith('MODE ') or s.upper().startswith('# MODE'):
                    break
    except Exception:
        return None
    return None


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
        with open(dat_path, "r", encoding="utf-8", errors="replace") as f:
            text = f.read()
        if parse_eigenvalues is not None:
            vals = parse_eigenvalues(text)
            if vals:
                return len(vals)
        n = _count_rows_after_header(text)
        if n:
            return n
    if msg_path and os.path.exists(msg_path):
        with open(msg_path, "r", encoding="utf-8", errors="replace") as f:
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
    with open(path, "r", encoding="utf-8", errors="replace") as f:
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
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s.startswith('#'):
                break
            body = s.lstrip('#').strip()
            if body.upper().startswith('MODE '):
                break          # '# MODE n frame=... lambda=...' = 모드 데이터 시작 -> 헤더 끝
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


def _norm_elem_code(x):
    """요소코드를 비교 가능한 문자열로 정규화한다.

    왜 필요한가(실측 2026-10-02): ELEM_CODE_QUAD 는 abaqusConstants 의 **SymbolicConstant** 이고
    SymbolicConstant 는 plain str 과 `==` 비교가 성립하지 않는다. 그래서 .inp 에서 파싱한
    문자열 'M3D4' 와 `SymbolicConstant(M3D4) in {...}` 가 False 가 되어 **가드가 오탐**했다
    (요소코드는 정상이었는데 '무시됐다'고 보고). .name 을 먼저 쓰고, 없으면 str() 로 떨어뜨린 뒤
    따옴표/공백을 벗긴다.
    """
    n = getattr(x, 'name', x)
    n = str(n).strip()
    return n.strip("'\"").strip()


#   [§19.112] RELEASE_EDGE_U3=True 는 자유변 면외 구속을 '하드 u3=0' 에서 to-ground 스프링으로
#     바꾼다 -> 덱에 **노드별 스프링 요소**(n = All_Edges_NoClamp 노드수)가 정당하게 생긴다.
#     CAE 는 그 요소를 `*Element, type=Spring1` 로 쓴다(고전 키워드 표기는 SPRINGA — 같은 스프링).
#     종전 가드는 'S' 로 시작하는 타입을 전부 셸 폴백 후보로 봤기 때문에 `Spring1` 이
#     extra=['Spring1'] 로 잡혀 **정상 덱을 제출 전에 막았다**(실측 2026-10-19).
EDGE_SPRING_ELEM_TYPES = ('Spring1', 'SPRINGA')


def verify_inp_element_types(inp_path, want_codes, allow_extra=None,
                             expect_extra_total=None):
    """방금 쓴 .inp 에서 `*Element, type=` 블록을 읽어 **요소 수까지** 세고, 요청 요소코드가
    실제로 들어갔는지 확인한다.

    왜 필요한가(실측 2026-10-02): ELEM_CODE_QUAD='M3D4' 로 바꿨는데도 CAE 가 setElementType 을
    조용히 무시하고 기본 셸(S4R/S3)로 남겼다 -> `***ERROR: 12086 elements have missing property
    definitions` (막 섹션 vs 셸 요소 불일치). 콘솔에는 단서가 없었다. '요소코드를 바꿨다'는
    선언은 검증이 아니다.

    allow_extra  이번 덱에서 **정당하게** 있을 수 있는 비(非)셸 요소 타입(예: §19.111/§19.112 의
                 `Spring1`). 기본 None = 종전과 **완전히 같은** 판정(하위호환). 여기 넣은 타입만
                 extra 목록에서 빠지고, 다른 셸/막 계열 타입은 그대로 실패시킨다 — 허용은 기능
                 스위치가 켜졌을 때만 호출자가 넘긴다.
    expect_extra_total  allow_extra 요소 수의 기대값(보통 스프링 대상 경계 노드수).
                 None 이면 개수 검증을 하지 않는다. 값이 있으면 합계와 비교해 **불일치면 실패**한다.
                 왜: 스프링이 절반만 생성되거나(셋이 비었거나) 다른 영역에 붙으면 BC 대체가
                 조용히 어긋난다 — '타입이 있으니 통과' 로는 잡히지 않는다.

    반환: (ok, counts, err_str)
      ok      = 요청 코드가 전부 .inp 에 존재하는가 (+ allow_extra 개수 기대값 일치)
      counts  = {요소코드: **요소 수**}   (블록 수가 아니라 요소 수 — 메쉬가 바뀌었는지 판정용)
    """
    counts, cur = {}, None
    try:
        with open(inp_path, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('*'):
                    if line.startswith('*Element, type='):
                        cur = _norm_elem_code(line.split('type=', 1)[1].split(',')[0])
                        counts.setdefault(cur, 0)
                    else:
                        cur = None
                    continue
                if cur is not None and line.strip():
                    counts[cur] += 1
    except Exception as e:
        return None, {}, '%s: %s' % (type(e).__name__, e)
    want = [_norm_elem_code(w) for w in want_codes]
    #   [§19.112] 허용목록(allow_extra)은 **비어 있으면 종전 판정과 문자 단위로 같다**.
    #     `t not in allow` 는 allow 가 빈 리스트일 때 항상 True 이므로 extra 계산이 동일하다.
    allow = [_norm_elem_code(t) for t in (allow_extra or ())]
    # 판정 규칙 (실측 2026-10-05 3차 교정):
    #   처음엔 `all(w in got)` 로 썼는데 **메쉬가 순수 사각형이라 S3 블록이 없는 정상 덱**을
    #   'setElementType 이 무시됐다'고 오탐했다(요청 (S4,S3) 인데 .inp 에 S4 만 존재).
    #   요소 타입은 **모양별로** 배정되므로, 어떤 모양의 요소가 없으면 그 코드는 안 나오는 게 정상이다.
    #   => 실패 조건은 두 가지다:
    #      (1) 요청 코드가 **하나도** 없다                      -> 전면 폴백
    #      (2) 요청 집합에 없는 셸/막 계열 타입이 .inp 에 있다   -> 부분 폴백(예: 사각형만 S4R)
    got = set(counts)
    _extra = [t for t in got
              if (t.startswith('S') or t.startswith('M3D') or t.startswith('STRI'))
              and t not in want and t not in allow]
    ok = any(w in got for w in want) and not _extra
    err = '' if ok else 'extra=%s want=%s' % (_extra, want)
    #   [§19.112] 개수 검증: 허용한 스프링이 **기대한 노드수만큼** 실제로 나갔는가.
    #     ('타입이 있으니 통과' 는 셋이 비거나 일부만 걸린 경우를 놓친다.)
    if expect_extra_total is not None:
        _n_extra = sum(counts.get(t, 0) for t in allow)
        if _n_extra != int(expect_extra_total):
            ok = False
            _detail = ('extra_count=%d != expected %d (허용=%s, 실제=%s)'
                       % (_n_extra, int(expect_extra_total), allow,
                          dict((t, counts[t]) for t in allow if t in counts)))
            err = (err + '; ' if err else '') + _detail
    return ok, counts, err


def job_completed_from_logs(job_name):
    """잡이 실제로 완주했는지 .msg/.sta 의 완주 문자열로 판정한다.

    왜 필요한가 (실측 2026-10-04)
      Abaqus/CAE `noGUI` 에서 `job.status` 는 **None** 으로 돌아온다(실측: 'Status: None').
      그래서 `job.status != COMPLETED` 로 쓰면 **완주한 잡도 실패로 오판**한다
      (실측: .msg 에 THE ANALYSIS HAS BEEN COMPLETED 가 있는데 ok=False 보고).
      반대로 `== ABORTED` 만 보면 Abaqus 가 실패해도 **부분 .odb 를 남기므로** 실패를 놓친다.
      => 판정은 산출물이 아니라 **완주 문자열**로 한다.
    """
    for f in (job_name + '.msg', job_name + '.sta'):
        try:
            with open(f, 'r', encoding='utf-8', errors='replace') as fh:
                t = fh.read()
        except Exception:
            continue
        if ('THE ANALYSIS HAS BEEN COMPLETED' in t) or ('HAS COMPLETED SUCCESSFULLY' in t):
            return True
    return False


# =====================================================================
# 설계점별 시드표 (2026-10-12) — x_c/d_c 로 표 경로를 결정하고, 지문으로 드리프트를 막는다.
#   HF(run_abaqus.py) 와 make_seed.py 가 **같은 규칙**을 써야 하므로 stdlib 전용인 여기에 둔다.
#   결정 근거: code/notes/plan-hf-seed-pipeline-2026-10-12.md (Q1=a, Q2=a, Q3=b)
# =====================================================================
BUCKLE_FINGERPRINT_KEYS = FINGERPRINT_KEYS + (
    'THICKNESS', 'USE_GRID_MESH', 'MESH_ORDER', 'CLAMP_DC', 'CLAMP_CABLE_LEN', 'CLAMP_EXCL_R',
    'CLAMP_DIR', 'BUCKLE_SOLVER', 'CORNER_ANGLE_DEG')
BUCKLE_SCRIPT_NAME = 'run_abaqus_buckle.py'
SEED_DIR_NAME = 'seeds'


def value_tag(v):
    """[§19.76] 수치 -> 파일명 태그. **반올림하지 않는다**.

    종전 `int(round(v*100))` 는 0 과 0.001 을, x_c=0.1234 와 0.1249 를 같은 이름으로 만들었다
    -> 서로 다른 설계점이 같은 시드표/ODB 를 덮어쓰거나 재사용해 조용한 오염이 가능했다.
    '%g' 는 값을 복원 가능한 문자열로 만든다(0 -> '0', 0.5 -> '0p5', 0.0556 -> '0p0556',
    1.8e-05 -> '1p8em05'). 파일명 안전: 영숫자와 p/m 만 쓴다.
    """
    return ('%g' % float(v)).replace('+', '').replace('-', 'm').replace('.', 'p')


def seed_table_name(x_c, d_c, prefix='modes'):
    """설계점 -> 시드표 파일명. 태그는 value_tag(HF 잡 이름과 같은 규칙)."""
    return '%s_xc%s_dc%s.txt' % (prefix, value_tag(x_c), value_tag(d_c))


def composite_table_name(x_c, prefix='modes'):
    """[§19.77] 합성 시드표 파일명 — **d_c 를 넣지 않는다**.

    모드 내용이 d_c 와 무관하기 때문이다(시드 클램프 당김은 0.01 고정, §19.70). 종전에는 HF 의
    d_c 가 이름에 들어가 새 d_c 마다 합성을 다시 해야 했는데, 그 d_c 는 모드에 아무 영향이 없다.
    """
    return '%s_xc%s.txt' % (prefix, value_tag(x_c))


def composite_table_path(seeds_dir, x_c, env=None):
    """합성 시드표 절대경로. env(SOLARSAIL_MODE_TABLE) > 신형 이름 > 구형(d_c 포함) 최신 폴백.

    구형 폴백: §19.77 이전에 만든 'modes_xc<x>_dc<d>.txt' 가 있으면 그중 **최신**을 쓴다.
    내용이 같으므로 재생성 없이 그대로 쓸 수 있다(어느 쪽을 썼는지는 basename 으로 보인다).
    """
    if env:
        return os.path.abspath(env)
    import glob
    new = os.path.join(seeds_dir, composite_table_name(x_c))
    if os.path.exists(new):
        return new
    hits = sorted(glob.glob(os.path.join(seeds_dir, 'modes_xc%s_dc*.txt' % value_tag(x_c))),
                  key=os.path.getmtime, reverse=True)
    return hits[0] if hits else new


def seed_table_path(seeds_dir, x_c, d_c, env=None):
    """합성 시드표 경로. [§19.77] d_c 는 이름·지문에서 빠졌다 — 인자는 호출부 호환용."""
    return composite_table_path(seeds_dir, x_c, env=env)


def buckle_fingerprint(script_path, keys=BUCKLE_FINGERPRINT_KEYS):
    """buckle 스크립트의 모델 상수 지문(16자). 못 읽으면 None."""
    return model_fingerprint(script_path, keys)[0]


def point_tag(x_c, d_c):
    """설계점 태그 — 산출물 이름에 넣어 짝 어긋남(mispairing)을 구조적으로 막는다."""
    return '_xc%s_dc%s' % (value_tag(x_c), value_tag(d_c))


def lf_odb_name(x_c, d_c):
    """설계점별 LF ODB 파일명. 고정 'LF_Analysis.odb' 를 대체한다(§19.66)."""
    return 'LF_Analysis%s.odb' % point_tag(x_c, d_c)


def lf_job_name(x_c, d_c):
    """설계점별 LF 잡 이름(run_job_safely 에 넘긴다)."""
    return 'LF_Analysis%s' % point_tag(x_c, d_c)


# ---------------------------------------------------------------------------
# [paper-match] LF ODB model sidecar (2026-10-21)
#   LF ODB 이름에는 (x_c, d_c) 만 있고 **모델 정의 버전**은 없다. LF 물리가 바뀌면
#   기존 ODB 가 조용히 재사용된다(예: 변위구동 LF -> 초기응력 프리스트레스 LF).
#   성공한 LF 잡 옆에 .model.txt 를 남기고, 재사용 전에 tag 를 대조한다.
# ---------------------------------------------------------------------------
LF_SIDECAR_SUFFIX = '.model.txt'
#   paper physical default: clamps=1, full-membrane u3=1. The tag includes the sub-mode so
#   switching PAPER_INCLUDE_CLAMPS / PAPER_HOLD_U3_GLOBAL invalidates old LF ODBs.
#   run_abaqus.py computes the same string from those two flags; keep them in sync.
LF_MODEL_TAG = 'globaltension-v1-c1-u1'


def lf_sidecar_path(odb_path):
    return odb_path + LF_SIDECAR_SUFFIX


def lf_sidecar_write(odb_path, tag=None, extra=None):
    """LF ODB 옆에 model sidecar 를 쓴다. 실패는 호출자가 처리한다."""
    kv = {'tag': tag or LF_MODEL_TAG,
          'written': time.strftime('%Y-%m-%d %H:%M:%S')}
    kv.update(extra or {})
    path = lf_sidecar_path(odb_path)
    with open(path, 'w') as f:
        for k in sorted(kv):
            f.write('%s=%s\n' % (k, kv[k]))
    return path


def lf_sidecar_read(odb_path):
    out = {}
    try:
        with open(lf_sidecar_path(odb_path), 'r') as f:
            for line in f:
                if '=' in line:
                    k, v = line.strip().split('=', 1)
                    out[k] = v
    except Exception:
        return {}
    return out


def lf_odb_is_current(odb_path, tag=None):
    """LF ODB 가 현재 모델 정의(tag)로 만들어진 것인가. sidecar 없으면 False."""
    if not odb_path or not os.path.exists(odb_path):
        return False
    return lf_sidecar_read(odb_path).get('tag') == (tag or LF_MODEL_TAG)


def file_token(path):
    """파일 내용 sha1 앞 16자. 없으면 None."""
    try:
        h = hashlib.sha1()
        with open(path, 'rb') as f:
            for chunk in iter(lambda: f.read(65536), b''):
                h.update(chunk)
        return h.hexdigest()[:16]
    except Exception:
        return None


def canon_seed_inputs(x_c, d_c, buckle_fp, seed_sha, clamp_sha):
    """지문 입력의 정본 문자열 — HF 와 make_seed 가 **같은 문자열**을 만들어야 한다.

    [§19.77] **d_c 를 뺀다**: 모드가 d_c 무관이므로(§19.70) d_c 가 다르면 같은 표를 쓰는 게 맞다.
    인자 d_c 는 호출부 호환을 위해 받기만 하고 지문에는 쓰지 않는다.
    """
    return 'x_c=%.4f|buckle=%s|seed=%s|clamp=%s' % (
        float(x_c), buckle_fp or 'none', seed_sha or 'none', clamp_sha or 'none')


def seed_fingerprint(x_c, d_c, buckle_fp, seed_sha, clamp_sha):
    return hashlib.sha1(canon_seed_inputs(x_c, d_c, buckle_fp, seed_sha, clamp_sha)
                        .encode('utf-8')).hexdigest()[:16]


def seed_header_lines(x_c, d_c, buckle_fp, seed_tbl, seed_modes, clamp_tbl, clamp_modes,
                      step='Step-Buckle', instance='MEMBRANE-1'):
    """합성 시드표 맨 앞에 붙일 헤더 줄들. (지문, [줄...]) 반환."""
    seed_sha = file_token(seed_tbl)
    clamp_sha = file_token(clamp_tbl)
    fp = seed_fingerprint(x_c, d_c, buckle_fp, seed_sha, clamp_sha)
    return fp, [
        '# FINGERPRINT sha1=%s src=%s' % (fp, BUCKLE_SCRIPT_NAME),
        '# SEED x_c=%g d_c_ref=%g buckle_fp=%s fingerprint=%s'
        % (float(x_c), float(d_c), buckle_fp or 'none', fp),
        '# SOURCE seed=%s seed_sha=%s modes=%s'
        % (os.path.basename(seed_tbl), seed_sha or 'none', seed_modes),
        '# SOURCE clamp=%s clamp_sha=%s modes=%s'
        % (os.path.basename(clamp_tbl), clamp_sha or 'none', clamp_modes),
        '# STEP %s' % step,
        '# INSTANCE %s' % instance,
    ]


def read_seed_header(table_path):
    """시드표 헤더(# 의 k=v 토큰)를 dict 로. 'MODE ' 줄을 만나면 멈춘다. 없으면 {}."""
    out = {}
    try:
        with open(table_path, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                if not s.startswith('#'):
                    break
                body = s.lstrip('#').strip()
                if body.upper().startswith('MODE '):
                    break
                for tok in body.split():
                    if '=' in tok:
                        k, v = tok.split('=', 1)
                        out.setdefault(k.strip(), v.strip())
    except Exception:
        return {}
    return out


def seed_gate(x_c, d_c, table_path, buckle_script, keys=BUCKLE_FINGERPRINT_KEYS):
    """이 시드표가 **이 설계점·이 buckle 소스**와 맞는가. (ok, 사유, 헤더dict) 반환.

    차단 사유는 사람이 읽고 바로 조치할 수 있게 쓴다(재생성 명령은 호출자가 붙인다).
    """
    info = read_seed_header(table_path)
    if not info:
        return False, '헤더를 읽지 못함(구버전 표이거나 파일 없음)', info
    if 'x_c' not in info or 'fingerprint' not in info:
        return False, '헤더에 x_c/fingerprint 가 없음(구버전 표)', info
    try:
        _xc = float(info['x_c'])
        _dc = float(info['d_c']) if 'd_c' in info else float('nan')
    except ValueError:
        return False, '헤더의 x_c/d_c 를 숫자로 읽지 못함', info
    if abs(_xc - float(x_c)) > 1e-9:
        return False, 'x_c 불일치(표 %.4f vs 런 %.4f)' % (_xc, float(x_c)), info
    #   [§19.77] d_c 는 비교하지 않는다 — 모드 내용이 d_c 무관이므로(§19.70) 다른 d_c 의 표를
    #     그대로 쓰는 것이 맞다. 표가 어느 d_c 로 만들어졌는지는 d_c_ref 로 참고만 한다.
    _now = buckle_fingerprint(buckle_script, keys)
    if _now != (info.get('buckle_fp') or None):
        return False, ('buckle 소스가 표 생성 이후 바뀜(현재 %s vs 표 %s) -> 재생성 필요'
                       % (_now, info.get('buckle_fp'))), info
    _exp = seed_fingerprint(_xc, _dc, info.get('buckle_fp'),
                            info.get('seed_sha'), info.get('clamp_sha'))
    if _exp != info.get('fingerprint'):
        return False, ('지문 재계산 불일치(%s vs 표 %s) -> 표가 손상/수정됨'
                       % (_exp, info.get('fingerprint'))), info
    #   [2026-10-12] 소스 모드표가 표 생성 뒤 바뀌면(재추출했는데 재합성 안 함) 조용히 낡은 표를 쓴다.
    #   같은 폴더에 소스가 남아 있으면 지금 해시와 대조한다(없으면 건너뜀 — 이동/삭제는 정상 흐름).
    _tdir = os.path.dirname(os.path.abspath(table_path))
    for _k in ('seed', 'clamp'):
        _nm, _sha = info.get(_k), info.get(_k + '_sha')
        if not (_nm and _sha):
            continue
        _p = os.path.join(_tdir, _nm)
        if os.path.exists(_p):
            _cur = file_token(_p)
            if _cur != _sha:
                return False, ('소스 모드표 %s 가 표 생성 이후 바뀜(%s vs %s) -> 재합성 필요'
                               % (_nm, _cur, _sha)), info
    return True, 'OK', info
