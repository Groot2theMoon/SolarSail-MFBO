# -*- coding: utf-8 -*-
"""run_abaqus.py 를 감싸서 **출력이 삼켜져도 흔적을 남긴다** (2026-10-08)

왜 필요한가
    CAE noGUI 는 스크립트 stdout 을 콘솔로 보내지 않고 버리는 경우가 있고, 프로세스
    종료코드는 0 으로 돌아온다. 그래서 화면에는 "Exit code: 0" 만 남는다 — HF 가 실제로
    이렇게 여러 번 죽었고 원인을 볼 수 없었다(§19.33, §19.34).

무엇을 하는가
    run_abaqus.py 를 exec 로 실행하되 **모든 stdout/stderr 를 로그 파일로 복제**한다.
    예외/SystemExit 도 잡아 traceback 을 같은 파일에 남긴다. 끝나면 한 줄 요약을 출력한다.
    -> 어떤 경우에도 '어디까지 갔는지'가 파일에 남는다.
    또한 noGUI 실행은 __file__ 을 정의하지 않는데, 여기서 넣어 주므로
    run_abaqus.py 의 스크립트-폴더 기준 경로 해석(모드 소스)도 정상 동작한다.

사용 (인자는 그대로 넘어간다)
    abaqus cae noGUI=run_hf_capture.py -- HF 0.5 1.0
    -> code 폴더에 hf_capture.log 가 생긴다. 마지막 40줄을 보면 된다.

인자(-- HF 0.5 1.0)는 sys.argv 로 그대로 전달된다.
"""
from __future__ import print_function

import io
import os
import sys
import time
import traceback


class _Tee(object):
    """여러 스트림에 동시에 쓰는 최소 래퍼(flush 포함)."""

    def __init__(self, *streams):
        self._s = streams

    def write(self, text):
        for s in self._s:
            try:
                s.write(text)
            except Exception:
                pass

    def flush(self):
        for s in self._s:
            try:
                s.flush()
            except Exception:
                pass


def main():
    try:
        here = os.path.dirname(os.path.abspath(__file__))
    except Exception:
        here = os.getcwd()
    target = os.path.join(here, 'run_abaqus.py')
    log = os.path.join(here, 'hf_capture.log')
    t0 = time.time()

    f = io.open(log, 'w', encoding='utf-8', errors='replace')
    old_out, old_err = sys.stdout, sys.stderr
    sys.stdout = _Tee(old_out, f)
    sys.stderr = _Tee(old_err, f)
    try:
        f.write('[capture] start %s\n' % time.strftime('%Y-%m-%d %H:%M:%S'))
        f.write('[capture] cwd=%s\n[capture] argv=%s\n[capture] target=%s\n'
                % (os.getcwd(), list(sys.argv), target))
        f.flush()
        if not os.path.exists(target):
            raise IOError('run_abaqus.py 를 찾을 수 없습니다: %s' % target)
        src = io.open(target, encoding='utf-8', errors='replace').read()
        g = {'__name__': '__main__', '__file__': target}
        exec(compile(src, target, 'exec'), g)
        f.write('[capture] 정상 종료 (%.1f s)\n' % (time.time() - t0))
    except SystemExit as e:
        f.write('[capture] SystemExit(code=%s) (%.1f s)\n'
                % (getattr(e, 'code', None), time.time() - t0))
    except Exception:
        f.write('[capture] === 예외로 중단 ===\n')
        traceback.print_exc(file=f)
        f.write('[capture] (%.1f s)\n' % (time.time() - t0))
    finally:
        f.write('[capture] 로그 파일: %s\n' % log)
        f.flush()
        try:
            f.close()
        except Exception:
            pass
        sys.stdout, sys.stderr = old_out, old_err
        print('[capture] 완료 — 로그: %s' % log)


main()
