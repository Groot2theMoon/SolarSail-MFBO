# -*- coding: utf-8 -*-
"""ctx.py - run_abaqus.py 의 실행 문맥: 로그 경로와 emit.

왜 별도 모듈인가
  `abaqus cae noGUI=script.py` 는 스크립트를 execfile 로 로드하므로 `__file__` 이 없고,
  stdout 도 콘솔로 오지 않는다(history §1). 그래서 로그는 **파일이 주 채널**이고,
  경로는 스크립트가 있는 폴더(code/)를 기준으로 잡아야 한다.
  모듈은 정상 import 되므로 `__file__` 이 있지만, 호출자가 `init(here)` 로 확정해 주면
  MFBO 가 주입한 경로(`MFBO_CODE_DIR`)까지 그대로 따른다.

계약
  - `emit()` 은 **어떤 경우에도 예외를 올리지 않는다** (로그 실패가 해석을 죽이면 안 된다).
  - `emit_start()` 는 로그를 truncate 하고 canary 한 줄을 쓴다.
"""

import os
import sys
import time

#   init() 전에도 emit 이 동작하도록 모듈 위치를 기본값으로 둔다(같은 code/ 폴더).
_HERE = os.path.dirname(os.path.abspath(__file__))
_EMIT_PATHS = None


def init(here):
    """호출자가 확정한 스크립트 폴더를 등록한다. 경로 캐시를 비운다."""
    global _HERE, _EMIT_PATHS
    if here:
        _HERE = here
    _EMIT_PATHS = None


def here():
    return _HERE


def emit_paths():
    """로그 파일 후보 경로 전부. __file__ 이 없는 execfile 환경을 전제로 한다."""
    cand = [_HERE, os.getcwd()]
    try:
        if sys.argv and sys.argv[0]:
            cand.append(os.path.dirname(os.path.abspath(sys.argv[0])))
    except Exception:
        pass
    out, seen = [], set()
    for d in cand:
        if not d:
            continue
        p = os.path.join(d, 'hf_run_log.txt')
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def emit(*a):
    """콘솔 + hf_run_log.txt(후보 경로 전부) 기록. **어떤 경우에도 예외를 올리지 않는다**."""
    global _EMIT_PATHS
    try:
        msg = ' '.join(str(x) for x in a)
    except Exception:
        msg = '<emit: 인자 변환 실패>'
    try:
        print(msg)
        sys.stdout.flush()
    except Exception:
        pass
    try:
        if _EMIT_PATHS is None:
            _EMIT_PATHS = emit_paths()
        for p in _EMIT_PATHS:
            try:
                try:
                    f = open(p, 'a', encoding='utf-8', errors='replace')
                except TypeError:          # encoding 인자를 모르는 파이썬이면 폴백
                    f = open(p, 'a')
                f.write(msg + '\n')
                f.close()
            except Exception:
                pass
    except Exception:
        pass


def first_existing(*cands):
    """후보 경로 중 처음 존재하는 것을 돌려준다(없으면 None).

    (2026-09-28) HF 는 os.chdir('aba') 로 돌기 때문에 'Buckle_Analysis.odb' 같은 상대 경로가
    code\\aba 에서는 안 찾히고, R-13/N-6 진단이 조용히 건너뛰어졌다(로그: odb=False).
    소스 잡의 산출물은 code\\ 에 있으므로 '..' 후보를 함께 본다.
    """
    for c in cands:
        if c and os.path.exists(c):
            return c
    return None


def probe_keep_line(_ls):
    """[§19.97] R-13 프로브 stdout 필터 - 로그(hf_run_log.txt)에 남길 줄만 True.

    probe 는 40여 줄을 쏟는다(요소 좌표·단면별 통계·압축지도·RF 상위 8개·인스턴스별 max|u3|).
    로그에 필요한 것은 (i) '>>>' 요약, (ii) ODB 신원/경고 줄, (iii) 중앙 영역 개수,
    (iv) 스텝 못박기 줄('스텝=...')이다.

    [§19.94] '>>>' 만 통과시키면 신원/경고 줄이 통째로 사라져 stale 판정 자체가 불가능해진다.
    [§19.97] 스텝 접미사([R-13:CT] / [R-13:GT])를 **떼고** 판정한다 -> 스텝이 늘어도 필터를
      고칠 필요가 없고, 접미사를 못 떼면 CT 판독이 통째로 로그에서 사라진다(같은 유형의 사고).
    """
    tag = _ls
    if tag.startswith('[R-13:'):
        _i = tag.find(']')
        if _i > 0:
            tag = '[R-13]' + tag[_i + 1:]
    return ('>>>' in tag) or tag.startswith(('[R-13] ODB 신원', '[R-13] !!!',
                                             '[R-13] 중앙 영역', '[R-13] 스텝='))


def emit_start():
    """런 시작 표식(파일을 truncate). 호출되면 이 파일이 실행된 것이 증명된다."""
    global _EMIT_PATHS
    _EMIT_PATHS = emit_paths()
    for p in _EMIT_PATHS:
        try:
            with open(p, 'w') as f:
                f.write('[canary] run_abaqus.py start %s  cwd=%s  argv=%s\n'
                        % (time.strftime('%Y-%m-%d %H:%M:%S'), os.getcwd(), list(sys.argv)))
        except Exception:
            pass
