# -*- coding: utf-8 -*-
"""solver.py — Abaqus 잡 실행 보조: 덱 패처, 완주 판정, 실패 진단.

run_abaqus.py 에서 **함수만** 옮겼다(2026-11-04, P2). `run_job_safely`/`run_job_with_ladder`
는 `mdb`/`abq` 와 30여 개 모델 상수에 직접 묶여 있어 이 모듈로 옮기지 않았다 — 그 둘은
스크립트에 남고, 여기 있는 함수는 경로(`here`)만 받는다.
"""

import io
import os
import re
import shutil
import subprocess
import sys
import traceback

from ctx import emit


def run_patcher(here, deck, step='Step-Postbuckle', ctrl=('--line-search', '--ran=5e-3',
                                                          '--can=1e-2', '--i0=8', '--ir=10'),
                stabilize=None):
    """패치 전 덱에 *Controls 를 주입해 같은 파일명으로 덮어쓴다.

    riks_patch_input.py 는 SRC 를 전부 메모리로 읽은 뒤 DST 를 열어 쓰므로
    SRC == DST (제자리 덮어쓰기) 가 안전하다. 덕분에 잡 이름/입력 파일명이
    그대로 유지되어 아래 job.submit() 흐름을 손대지 않아도 된다.

    반환: True(성공) / False(실패 — 호출자가 RESULTS:FAIL 로 끊는다)
    """
    #   [§19.85] stabilize 가 주어지면 패치된 덱의 `*Static, stabilize=...` 값만 교체한다.
    #     (패처 자체의 안정화 기본값에 의존하지 않고 사다리를 돌리기 위한 최소 개입)
    def _set_stabilize(_p, _v):
        with open(_p, 'r', encoding='utf-8', errors='replace') as _f:
            _t = _f.read()
        _new, _n = re.subn(r'(\*Static,\s*stabilize=)[0-9.eE+-]+', r'\g<1>%g' % _v, _t)
        if _n != 1:
            emit("[LADDER] !!! stabilize 치환 실패(패턴 %d건) — 덱을 확인하세요: %s" % (_n, _p))
            return False
        with open(_p, 'w', encoding='utf-8') as _f:
            _f.write(_new)
        emit("[LADDER] stabilize=%g 로 덱 수정 (%s)" % (_v, os.path.basename(_p)))
        return True

    _rp = os.path.join(here, 'riks_patch_input.py')
    if not os.path.exists(_rp):
        emit("[PATCHER] !!! riks_patch_input.py 없음: %s" % _rp)
        return False
    if not os.path.exists(deck):
        emit("[PATCHER] !!! 덱 없음: %s" % deck)
        return False
    # 멱등성: 이미 *Controls 가 들어 있으면 먼저 걷어낸다.
    #   덱은 보통 run_abaqus.py 가 writeInput 으로 새로 쓰므로 초기화되지만,
    #   'write' 로 만든 덱을 'auto' 로 다시 돌리면 주입이 두 번 되어
    #   Abaqus 가 모르는 파라미터를 읽는다. 걷어내고 항상 같은 결과를 만든다.
    _lines = io.open(deck, 'r', errors='replace').read().splitlines()
    if any(ln.lstrip().lower().startswith('*controls') for ln in _lines):
        _keep, _skip = [], False
        for ln in _lines:
            _ls = ln.lstrip()
            if _ls.startswith('*'):
                # 새 키워드 줄 — *Controls 면 블록 시작, 아니면 블록 종료
                _skip = _ls.lower().startswith('*controls')
            if not _skip:
                _keep.append(ln)
        io.open(deck, 'w', encoding='ascii', errors='replace',
                newline='\n').write('\n'.join(_keep) + '\n')
        emit("[PATCHER] 기존 *Controls %d 줄 제거 (멱등성)"
             % (len(_lines) - len(_keep)))
    before = os.path.getsize(deck)
    #   [2026-10-08] `abaqus python <script>` 가 Windows 에서 **조용히 아무 일도 안 하는**
    #   사례가 나왔다: exit 0, stdout 없음, 덱 무변경 -> 크기 검사에 걸려 HF 가 멈췄다.
    #   riks_patch_input.py 는 io/os/sys 만 쓰는 순수 텍스트 변환기라 Abaqus 인터프리터가
    #   필요 없다. 그래서 런처를 거치지 않고 **이 프로세스 안에서** 실행한다.
    #   (로컬 검증: 가짜 덱에 *Controls/line search 마커가 실제로 붙고 8줄 -> 15줄)
    #   실패하면 셸(cmd.exe) 경로로 한 번 더, 그래도 안 되면 손으로 돌릴 명령을 찍어 준다.
    _MARKERS = ('*Controls', 'line search', 'discontinuous', 'field=displacement')

    def _deck_markers():
        _t = io.open(deck, 'r', errors='replace').read().lower()
        return [m for m in _MARKERS if m.lower() in _t]

    def _run_inproc():
        _argv = [os.path.basename(_rp), deck, deck, step] + list(ctrl)
        _src = io.open(_rp, 'r', errors='replace').read()
        _save = list(sys.argv)
        _o, _e = io.StringIO(), io.StringIO()
        _so, _se = sys.stdout, sys.stderr
        sys.argv = _argv
        sys.stdout, sys.stderr = _o, _e
        try:
            exec(compile(_src, _rp, 'exec'), {'__name__': '__main__', '__file__': _rp})
            return True, _o.getvalue(), _e.getvalue()
        except SystemExit as _ex:
            return (not getattr(_ex, 'code', 0)), _o.getvalue(), \
                _e.getvalue() + ('\nSystemExit(%s)' % getattr(_ex, 'code', None))
        except Exception:
            return False, _o.getvalue(), _e.getvalue() + '\n' + traceback.format_exc()
        finally:
            sys.argv, sys.stdout, sys.stderr = _save, _so, _se

    _ok, _out, _err = _run_inproc()
    for _ln in (_out or '').splitlines()[-14:]:
        emit("[PATCHER] %s" % _ln.strip())
    for _ln in (_err or '').splitlines()[-14:]:
        emit("[PATCHER] (stderr) %s" % _ln.strip())
    emit("[PATCHER] in-process ok=%s / 마커=%s / 크기 %d -> %d"
         % (_ok, _deck_markers() or '없음', before, os.path.getsize(deck)))

    if not _deck_markers():
        _abq2 = shutil.which('abaqus') or shutil.which('abaqus.bat') or 'abaqus'
        _cmd = '"%s" python "%s" "%s" "%s" %s %s' % (
            _abq2, _rp, deck, deck, step, ' '.join(ctrl))
        emit("[PATCHER] 재시도(셸): %s" % _cmd)
        try:
            _p = subprocess.run(_cmd, shell=True, cwd=os.getcwd(),
                                capture_output=True, text=True)
            for _ln in (_p.stdout or '').splitlines()[-10:]:
                emit("[PATCHER] %s" % _ln.strip())
            for _ln in (_p.stderr or '').splitlines()[-10:]:
                emit("[PATCHER] (stderr) %s" % _ln.strip())
            emit("[PATCHER] 셸 rc=%s / 마커=%s" % (_p.returncode, _deck_markers() or '없음'))
        except Exception as _e2:
            emit("[PATCHER] 셸 실행 실패: %s" % _e2)

    if not _deck_markers():
        emit("[PATCHER] !!! 패치 마커가 덱에 없다 — 손으로 실행해 보세요:")
        emit('[PATCHER]     abaqus python riks_patch_input.py "%s" "%s" %s %s'
             % (deck, deck, step, ' '.join(ctrl)))
        return False
    emit("[PATCHER] OK — %s 패치 확인 (%.1f KB -> %.1f KB, 마커 %s)"
         % (deck, before / 1024.0, os.path.getsize(deck) / 1024.0, _deck_markers()))
    return True


def job_completed_ok(job_name):
    """R-12: .sta/.msg 에 'HAS COMPLETED SUCCESSFULLY' 가 있는지로 완주를 판정한다.

    odb 파일 존재만 보면 '중도 중단된 odb'를 성공으로 오판한다
    (2026-09-21 사례: HF_Postbuckle.odb 는 존재하나 Step-Postbuckle 프레임 0개).
    """
    for ext in ('sta', 'msg'):
        fn = '%s.%s' % (job_name, ext)
        if os.path.exists(fn):
            with open(fn, 'r', errors='replace') as f:
                if 'HAS COMPLETED SUCCESSFULLY' in f.read().upper():
                    return True
    return False


def print_job_diag(job_name):
    """잡 산출물(.msg/.dat)의 원인 판별용 핵심 줄만 콘솔에 찍는다.

    로그를 통째로 붙여넣지 않아도 (a) 좌굴모드가 나왔는지 (b) 왜 죽었는지를
    한 화면에서 볼 수 있게 하기 위한 것. 실패해도 해석에는 영향 없음.
    """
    keys = ('NEGATIVE EIGENVALUES', 'CONVERGED', 'EIGENVALUES CANNOT BE FOUND',
            'HAS COMPLETED SUCCESSFULLY', 'THE ANALYSIS HAS BEEN COMPLETED',
            'misplaced', 'STIFFNESS MATRIX IS SINGULAR', 'TOO MANY ATTEMPTS',
            #   [2026-10-08] 실제로 '적용된' 솔루션 컨트롤을 로그에 싣는다. 문서:
            #   "The controls in effect for an analysis are listed in the .dat and .msg files.
            #    Nondefault controls are marked by ***".  *CONTROLS 블록을 덱에 넣었는데
            #   효과가 없을 때(라인서치 0회 등) 이 목록이 유일한 판정 근거다.
            'CONVERGENCE TOLERANCE PARAMETERS', 'CRIT. FOR', 'LINE SEARCH',
            '***ERROR')
    for ext in ('msg', 'dat'):
        fn = '%s.%s' % (job_name, ext)
        if not os.path.exists(fn):
            emit("[DIAG:%s] %s 없음" % (job_name, fn))
            continue
        size = os.path.getsize(fn)
        with open(fn, 'r') as f:
            if size > 2000000:
                f.seek(size - 2000000)
            text = f.read()
        #   [2026-10-08] 매치 줄만 찍으면 원인을 가린다: Abaqus 는 오류문을 줄바꿈하므로
        #   "AN INITIAL CONDITION HAS BEEN SPECIFIED ON ELEMENT SET" 다음 줄에 **집합 이름**이
        #   온다. 그 줄은 키워드에 안 걸려 버려졌고, 매번 원인을 추측해야 했다(HF 가 여기서 막혔다).
        #   그래서 매치 지점마다 **뒤 2줄까지** 줄번호와 함께 남긴다(덱 대조용).
        #   [2026-10-08] 먼저 **매치 줄 번호만** 모으고(정수 리스트 = 메모리 무시 가능),
        #   출력은 앞 4개 + 끝 11개를 쓴다. 매치마다 문자열을 쌓으며 상한으로 끊으면
        #   8MB .msg 에서 **마지막(=중단 원인)** 을 잃는다 — 그 실수를 두 번 했다.
        _lines = text.splitlines()
        _idx = [i for i, _ln in enumerate(_lines)
                if _ln.strip() and any(k.lower() in _ln.lower() for k in keys)]
        emit("[DIAG:%s] %s (%.0f KB) 핵심줄 %d개 (문맥 포함)"
              % (job_name, fn, size / 1024.0, len(_idx)))
        _keep = _idx if len(_idx) <= 16 else (_idx[:4] + [None] + _idx[-11:])
        for _i in _keep:
            if _i is None:
                emit("      | ... (중략 %d개) ..." % (len(_idx) - 15))
                continue
            for _j in range(_i, min(_i + 3, len(_lines))):
                if _lines[_j].strip():
                    emit("      | %d| %s" % (_j + 1, _lines[_j].strip()[:170]))
