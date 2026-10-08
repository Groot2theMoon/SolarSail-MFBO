"""
Solar Sail의 삼각형 멤브레인 모델을 생성하고, 입력 변수(Clamp 위치, 변위량)에 따라
구조 해석(Linear/Post-buckling)을 수행하는 Abaqus Python Script.

Usage:
    abaqus cae noGUI=run_abaqus.py -- [Fidelity] [x_c] [d_c]

Arguments:
    - Fidelity (str): 'LF' (Low-Fidelity, Linear/Static) or 'HF' (High-Fidelity, Post-buckling)
    - x_c (float): Clamp Cable의 부착 위치 파라미터 (0.0 ~ 1.0)
                    0에 가까울수록 위쪽 꼭짓점, 1에 가까울수록 변의 아래 끝점에 위치.
    - d_c (float): Clamp Cable을 당기는 변위 비율 (Displacement Ratio).
                    각 step에서 꼭짓점 변위량에 대한 배율로 작용.

Flow (Simulation Logic):
    1. 형상(Geometry) & 물성(Material) & 메쉬(Mesh) 생성 (공통)
    2. 초기 장력(Global Tension) 단계 정의
    3. Fidelity에 따른 분기:
       [LF 모드]  (클램프 있음 + 완전 형상)
         - 'Step-GlobalTension' -> 'Step-ClampTension' -> 'Step-HighTension' 순차 진행.
         - 주름(Wrinkle) 거동을 무시하고 선형적(혹은 단순 기하 비선형) 강성을 빠르게 계산.
         - 임퍼펙션(*IMPERFECTION) 없음 -> 좌굴(Buckle) 잡도 돌리지 않는다(모드 불필요).
       [HF 모드]  (클램프 있음 + 임퍼펙션 주입 = 2-모델 레시피, 2026-09-22)
         - 1차: 모드 소스는 '클램프 없는' 별도 모델에서 **선형 좌굴해석(*BUCKLE)** 으로 모드를
                뽑는다(run_abaqus_mode.py). 클램프가 있는 base state 는 선형 좌굴 스펙트럼을
                주지 못한다(실측: 음수고유값 598~2897 / CONVERGED=0 / EIGENVALUES CANNOT BE FOUND).
                *BUCKLE 은 .fil 출력이 금지되므로(실측 7025) 모드는 ODB 모드 프레임에만 있다
                -> aba_mode_from_odb.py 가 그 ODB 를 모드표(modes_ClampFree_Buckle.txt)로 만든다.
         - 2차: 그 모드를 노드 좌표 섭동으로 주입한다(IMPERFECTION_MODE='odb_direct' = 기본,
                ODB 에서 HF 가 직접 읽는다. 'odb_table' = 모드표 txt 경유).
                aba_imperfection.py 가 모드 개수를 검증하고, 부족하면 HF 제출 전에 중단한다.
                자기 좌굴 잡(Buckle_Analysis)은 클램프 base state 증거용으로만 유지한다.
                (옵션 IMPERFECTION_MODE='file' = .fil 스테이징 + *IMPERFECTION, FILE=, STEP=n:
                 사용자 원본 08d4cbc 방식. 단 좌굴모드에는 쓸 수 없다 — *BUCKLE 은 .fil 금지.)
         - 3차: 'Step-Postbuckle' 수행 (Riks/Stabilization). 주름 거동을 포함한 비선형 해석.
    4. Post-Processing:
       - 'eval_abaqus.py'를 subprocess로 호출하여 ODB에서 필요한 값(추력, 면적 등)만 추출.
       - 추출된 값을 stdout으로 출력하여 상위 프로세스(MFBO)에 전달.

solar-sail 해석 방법론은 대부분 Galhofo2022 논문을 레퍼런스로 하였음.
(물성치, 크기, 각도, step의 종류와 순서 등등)
 - displacementBC 값들은 경험적으로 찾은 값들. 더 좋은 숫자가 있을 수 있음.
"""

from abaqus import *
from abaqusConstants import *
from step import *
from mesh import ElemType
import regionToolset
import interaction
import sys
import os
import io
import shutil
import traceback
import subprocess
import time
import numpy as np

# ---- 스크립트 위치(_HERE) / Abaqus 작업 디렉터리(_RUN) 결정 ----
# 주의: `abaqus cae noGUI=script.py` 로 실행하면 스크립트가 execfile 로 로드되어
#       __file__ 이 정의되지 않는다(NameError). 아래처럼 우선순위를 둔다.
#   (1) mfbo.py 가 주입한 MFBO_CODE_DIR  -> 가장 확실
#   (2) __file__ (abaqus python 등으로 직접 실행될 때)
#   (3) sys.argv[0] / cwd / cwd\code 중 eval_abaqus.py 가 실제로 있는 곳
def _resolve_here():
    def _ok(d):
        return bool(d) and os.path.exists(os.path.join(d, "eval_abaqus.py"))
    try:
        d = os.path.dirname(os.path.abspath(__file__))
        if _ok(d):
            return d
    except NameError:
        pass
    cand = []
    try:
        if sys.argv and sys.argv[0]:
            cand.append(os.path.dirname(os.path.abspath(sys.argv[0])))
    except Exception:
        pass
    cand.append(os.getcwd())
    # cwd = code/aba (mfbo.py 구동 시) 이면 그 부모가 code 다.
    cand.append(os.path.dirname(os.getcwd()))
    cand.append(os.path.join(os.getcwd(), "code"))
    for c in cand:
        if _ok(c):
            return os.path.abspath(c)
    emit("!!! WARNING: eval_abaqus.py 위치를 찾지 못했습니다. cwd=%s" % os.getcwd())
    return os.getcwd()

_HERE = _resolve_here()
# 산출물 디렉터리 = code/aba. mfbo.py 의 RUN_DIR_NAME 과 같은 값이어야 한다.
RUN_DIR_NAME = "aba"
# Abaqus 병렬 스레드/도메인 수 (1 = 단일 CPU, 4 = 대략 1.5~2.5배 빠름).
#   라이선스 토큰이 없으면 잡이 라이선스 오류로 즉시 죽는다 -> 그때는 1 로.
NUMCPUS = 4
BASE_PROBE = True                # 각 잡 직후 base_state_probe.py (R-13) 자동 실행
EIG_RECORD = True                # 좌굴 고유치 이력 기록 (coalescence_check.py record)
# ---- 코드 상수 끝 ----
_RUN = os.path.join(_HERE, RUN_DIR_NAME)
os.makedirs(_RUN, exist_ok=True)
os.chdir(_RUN)

# ---- 임퍼펙션 소스 스테이징 모듈 (2-모델 레시피, stdlib only) ----
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from aba_imperfection import (ImperfectionSourceError, stage,           # noqa: E402
                              imperfection_text, report, load_mode_table,
                              build_perturbation, perturbation_report, mode_table_info,
                              verify_inp_element_types, job_completed_from_logs)
# aba_grid_mesh: boundary_node_labels / clamp_exclude_labels / make_set (면외 z 구속 노드셋). §9
#   격자는 이제 사용 중이다: USE_GRID_MESH=True -> aba_grid_mesh.import_grid_part
#   (.inp -> PartFromInputFile). fill_part 는 Part.addNodes 가 없어 폐기했다(§19.15).
import aba_grid_mesh                                                            # noqa: E402
# 로깅 — cae noGUI 는 스크립트 stdout 을 콘솔로 보내지 않고, execfile 이라 __file__ 도 없다(§1).
try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass

_EMIT_PATHS = None


def _emit_paths():
    """로그 파일 후보 경로 전부. __file__ 이 없는 execfile 환경을 전제로 한다."""
    cand = [globals().get('_HERE'), os.getcwd()]
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
            _EMIT_PATHS = _emit_paths()
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


def _first_existing(*cands):
    """후보 경로 중 처음 존재하는 것을 돌려준다(없으면 None).

    (2026-09-28) HF 는 os.chdir('aba') 로 돌기 때문에 'Buckle_Analysis.odb' 같은 상대 경로가
    code\aba 에서는 안 찾히고, R-13/N-6 진단이 조용히 건너뛰어졌다(로그: odb=False).
    소스 잡의 산출물은 code\ 에 있으므로 '..' 후보를 함께 본다.
    """
    for c in cands:
        if c and os.path.exists(c):
            return c
    return None


def emit_start():
    """런 시작 표식(파일을 truncate). 호출되면 이 파일이 실행된 것이 증명된다."""
    global _EMIT_PATHS
    _EMIT_PATHS = _emit_paths()
    for p in _EMIT_PATHS:
        try:
            with open(p, 'w') as f:
                f.write('[canary] run_abaqus.py start %s  cwd=%s  argv=%s\n'
                        % (time.strftime('%Y-%m-%d %H:%M:%S'), os.getcwd(), list(sys.argv)))
        except Exception:
            pass


emit_start()


try:
    # abaqus cae noGUI=run_abaqus.py -- [HF/LF] x_c d_c
    fidelity = sys.argv[-3].upper()
    x_c = float(sys.argv[-2])
    d_c = float(sys.argv[-1])
    # 범위 가드: 두 클램프 패치(반경 0.2 m)는 사선변 위에서 거리 20*x_c 만큼 떨어진다.
    #   x_c < 0.05 이면 두 패치가 겹치고, x_c = 1.0 이면 꼭짓점 RP 패치와 겹친다.
    #   상한은 mfbo.py 의 설계공간(x1_max=0.95)과 반드시 일치해야 한다 — 좁게 잡으면
    #   MFBO 가 그 구간을 탐색할 때 조용히 sys.exit 으로 죽는다.
    if not (0.05 <= x_c <= 0.95):
        emit('!!! ERROR: x_c=%.4g 는 유효범위 [0.05, 0.95] 밖이다 (패치 겹침).' % x_c)
        sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
        sys.exit(1)
    if not (0.0 < d_c <= 2.0):
        emit('!!! ERROR: d_c=%.4g 는 유효범위 (0, 2.0] 밖이다.' % d_c)
        sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
        sys.exit(1)
    
except:
    emit("Error: Invalid arguments. Usage: abaqus cae noGUI=run_abaqus.py -- HF x_c d_c")
    sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
    sys.exit(1)

# P1-3: 잡 제출 전에 지워야 하는 이전 실행 산출물
_JOB_ARTIFACTS = ('odb', 'fil', 'sta', 'msg', 'lck', 'com', 'prt', 'sim', 'log',
                  'dat', 'res', 'abq', 'ipm', 'mdl', 'stt', 'cid')

def run_job_safely(job_name, model_name=None):
    """
    job 실행 중 .odb, .lck 파일 충돌을 방지하고, 완료까지 대기하는 함수

    P1-3: 제출 전에 이전 실행 산출물을 지운다.
      - 성공 판정(not os.path.exists(odb))이 '지난 실행의 odb'를 보고 오판하는 것을 막고,
      - *IMPERFECTION, FILE=Buckle_Analysis 가 stale .fil 을 조용히 읽는 것을 막는다.
    """
    model_name = model_name or MODEL_NAME
    for _ext in _JOB_ARTIFACTS:
        _f = '%s.%s' % (job_name, _ext)
        if os.path.exists(_f):
            try:
                os.remove(_f)
            except OSError as _e:
                emit("  (warning) could not remove %s: %s" % (_f, _e))

    lck_file = job_name + '.lck'
    odb_file = job_name + '.odb'
    # 기존 Lock 파일이 있다면 삭제 시도 < 이전 실행 강제 종료 시 남게 됨
    if os.path.exists(lck_file):
        emit("Detected old lock file: %s. Removing it..." % lck_file)
        try:
            os.remove(lck_file)
        except OSError:
            # 만약 삭제가 안 된다면 다른 프로세스가 실제로 사용 중
            emit("!!! FATAL ERROR: Cannot remove lock file. Is another Abaqus process running?")
            sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
            sys.exit(1)

    if job_name in mdb.jobs:
        del mdb.jobs[job_name]
    
    # 병렬 실행 (속도). 기본 1 = 기존과 완전히 동일한 거동.
    #   단일 CPU 로 373 증분/30분 수준이면 4 스레드로 대략 1.5~2.5배 단축 여지가 있다.
    #   주의: 라이선스 토큰이 부족하면 잡이 라이선스 오류로 죽는다 -> 그때는 1 로 되돌린다.
    #   라이선스 토큰 오류가 나면 위 NUMCPUS 상수를 1 로 바꾼다.
    _ncp = NUMCPUS
    emit("[run_abaqus] 케이스: x_c=%.4g  d_c=%.4g  클램프=%s  RUN_MODE=%s  NUMCPUS=%d"
         % (x_c, d_c, "OFF" if NO_CLAMP else "ON", RUN_MODE, _ncp))
    emit("   DISP_GLOBAL=%.6g m  SEED_DIV=%.6g  THICKNESS=%.3g m  stabilize=0.003"
         % (DISP_GLOBAL, SEED_DIV, THICKNESS))
    emit("   GLOBAL_FINAL=%.6g m  PRETENSION_SCALE=%.6g  비율=%.1f"
         % (GLOBAL_FINAL, PRETENSION_SCALE, GLOBAL_FINAL / (0.0001 * PRETENSION_SCALE)))
    job = mdb.Job(name=job_name, model=model_name, numCpus=_ncp, numDomains=_ncp)
    emit("Submitting Job: %s" % job_name)
    job.writeInput(consistencyChecking=OFF)
    # 제출 전 요소 타입 검증 — 조용한 폴백을 Standard 토큰을 태우기 전에 잡는다(§3).
    _ok_et, _cnt_et, _err_et = verify_inp_element_types(
        job_name + '.inp', (ELEM_CODE_QUAD, ELEM_CODE_TRI))
    emit("[GUARD] .inp 요소 타입 블록: %s" % (_cnt_et or _err_et))
    if _err_et or not _ok_et:
        emit("!!! ERROR: .inp 요소 타입 검증 실패 (요청=%s): %s"
             % ((ELEM_CODE_QUAD, ELEM_CODE_TRI), _err_et or _cnt_et))
        emit("    잡을 제출하지 않는다. ELEM_CODE_* 지정 방식을 고쳐라.")
        return False
    # 'write' 는 덱만 생성하고 종료한다(수동 패처 워크플로우).
    if RUN_MODE == 'write':
        emit("[run_abaqus] RUN_MODE='write' : 덱만 생성(%s.inp)하고 제출하지 않는다." % job_name)
        emit("             다음: python riks_patch_input.py aba\\%s.inp <out>.inp Step-Postbuckle "
             "--line-search --ran=2e-2 --can=1e-2 --i0=8 --ir=10" % job_name)
        sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
        sys.exit(0)
    # 'auto' 는 제출 전에 *Controls 를 주입한다. 제자리 덮어쓰기라 잡 이름/입력명이 유지된다.
    if RUN_MODE == 'auto':
        if not run_patcher(job_name + '.inp'):
            emit("!!! ERROR: 패처 실패 — 패치되지 않은 덱을 제출하지 않는다(§14).")
            emit("RESULTS:FAIL")
            sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
            sys.exit(1)
    else:
        emit("[run_abaqus] RUN_MODE='plain' : 패처 없이 원본 덱을 제출한다(진단용).")
    job.submit(consistencyChecking=OFF)
    
    job.waitForCompletion()

    time.sleep(1.0)

    
    # 잡 성공 판정 = 산출물 + 완주 문자열. 'ODB 존재'도 'Exit code 0'도 증거가 아니다(§2).
    if (job.status == ABORTED or not os.path.exists(odb_file)
            or not job_completed_from_logs(job_name)):
        print_job_diag(job_name)   # 원인 판별용 — 실패했을 때만(완주 잡에는 노이즈)
        emit("!!! ERROR: Job %s failed. Actual Status: %s" % (job_name, str(job.status)))
        sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
        sys.exit(1)
        
    if not job_completed_ok(job_name):
        emit("!!! WARNING: %s — .sta/.msg 에 'HAS COMPLETED SUCCESSFULLY' 없음 "
              "(중도 중단 의심; odb 존재만으로는 판정 불가 — R-12)" % job_name)
    emit("Job %s completed successfully (Status: %s)." % (job_name, str(job.status)))
    return True

def run_patcher(deck, step='Step-Postbuckle', ctrl=('--line-search', '--ran=2e-2',
                                                     '--can=1e-2', '--i0=8', '--ir=10')):
    """패치 전 덱에 *Controls 를 주입해 같은 파일명으로 덮어쓴다.

    riks_patch_input.py 는 SRC 를 전부 메모리로 읽은 뒤 DST 를 열어 쓰므로
    SRC == DST (제자리 덮어쓰기) 가 안전하다. 덕분에 잡 이름/입력 파일명이
    그대로 유지되어 아래 job.submit() 흐름을 손대지 않아도 된다.

    반환: True(성공) / False(실패 — 호출자가 RESULTS:FAIL 로 끊는다)
    """
    _rp = os.path.join(_HERE, 'riks_patch_input.py')
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

sqrt2 = 1.414
N_EIG = 4          # 임퍼펙션에 쓸 좌굴모드 수 (*IMPERFECTION / *NODE FILE)

# 2-모델 레시피 — 모드는 클램프-프리 모델에서 뽑아 클램프 모델에 임퍼펙션으로 주입한다(§5).
#   클램프 패치는 메쉬를 바꾸지 않으므로 두 모델의 막 노드 라벨이 같다 -> *IMPERFECTION 매핑이 성립한다.
MODE_SOURCE = 'external'     # 'external' = 클램프 없는 외부 .fil(run_abaqus_mode.py) | 'self' = 자기 좌굴 잡
#   모드 소스 = run_abaqus_mode.py (클램프 없음, HF 파라미터 정렬). 산출물은 code\ 에 쌓인다.
#   (이전 소스였던 run_abaqus_cable.py 의 Buckle_Analysis.fil 을 쓰려면 아래 이름만 교체한다.)
MODE_SOURCE_FIL = os.path.join('..', 'ClampFree_Buckle.fil')   # code\aba -> code\
MODE_SOURCE_DAT = os.path.join('..', 'ClampFree_Buckle.dat')
MODE_SOURCE_MSG = os.path.join('..', 'ClampFree_Buckle.msg')
MODE_SOURCE_STEP = 2         # 소스 .fil 안의 스텝 번호 (케이블 런: 1=GlobalTension 2=Buckle)
IMPERFECTION_NAME = 'ClampFree_Buckle'   # 모드 소스 잡 이름과 '같은' 이름 -> 원장 C-1(조용한 0-모드 소비) 차단
                                         # (MODE_SOURCE_FIL 의 basename 과 함께 움직여야 한다)
IMPERFECTION_MODES = (1, 2, 3, 4)   # 기본 = 논문(Galhofo) 주입 모드 1~4.
IMPERFECTION_MAX_MODES = 4          # 상한. 모드표에 있는 모드가 더 적으면 그만큼만 쓴다(아래에서 자동 조정).
                                    #  모드 3개를 주입한다(lambda1 탈락). 논문 lambda 스프레드 0.036% 라 모드 3.4 의 기여가 작다(§5a).
IMPERFECTION_AMPL_T = 1.0    # 막 두께 배수. Galhofo 채택값은 0.10 t 이지만 우리는 1.0 을 쓴다.
#   [2026-10-06 실측] 0.10 이면 증분이 눌려 2215 증분 / 약 90 분이 걸린다. 1.0 이면 271 증분 / 11 분이다
#   (HF_x025k). 진폭이 작으면 초기 결함이 약해 분기가 늦게 갈라지고 수렴이 어려워진다.
#   진폭 민감도를 보려면 이 상수를 0.50 등으로 바꿔 재실행한다.
RUN_SELF_BUCKLE_JOB = False  # 자기(클램프) 좌굴 잡을 끈다 (2026-09-24). 이유는 측정된 실패 경로다:
                             #   이 잡은 _EIGENSOLVER='LANCZOS' 로 제출되는데(아래), 이 모델의 좌굴 base state 는
                             #   이미 분기하중을 넘은 부정정 상태다(시스템 음수 고유값 실측 56~88개). 매뉴얼이 열거한
                             #   LANCZOS 금지 조건 중 'preloaded above the bifurcation load' 에 정확히 해당하므로
                             #   ***ERROR: THE LANCZOS SOLVER CANNOT BE USED FOR BUCKLING ANALYSIS 로 거부된다.
                             #   run_job_safely 는 ABORTED 를 보면 sys.exit(1) 하므로, 이 잡이 켜져 있으면
                             #   HF_Postbuckle 이 제출되기도 전에 스크립트가 끝난다(제출 순서: LF 840 -> Buckle 862 -> HF 1010).
                             # SUBSPACE 로 바꿔도 클램프 모델은 좌굴모드 0개로 죽어 같은 경로다. 모드 추출은 이제
                             #   ODB 모드표 경로(run_abaqus_mode.py -> modes_ClampFree_Buckle.txt)가 담당하므로
                             #   이 진단 잡은 불필요하다. 클램프 base state 증거가 다시 필요하면 True 로 되돌리되
                             #   그때는 _EIGENSOLVER 와 제출 순서를 함께 손봐야 한다.
if MODE_SOURCE not in ('external', 'self'):
    raise RuntimeError("MODE_SOURCE 는 'external' 또는 'self' 여야 합니다 (현재 %r)" % (MODE_SOURCE,))

# 임퍼펙션 주입 방식 — 모드 소스의 스텝 타입이 방식을 결정한다(§6).
#   *BUCKLE   -> 'odb_direct'(기본). *BUCKLE 은 .fil 출력이 금지되므로 모드 프레임은 ODB 에만 있다.
#   *FREQUENCY -> 'file' (.fil + *IMPERFECTION, FILE=).
IMPERFECTION_MODE = 'odb_direct'   # 'odb_direct'(기본) | 'odb_table' | 'file'
MODE_TABLE = os.path.join('..', 'modes_ClampFree_Buckle.txt')   # code\aba -> code\
MODE_SOURCE_ODB = os.path.join('..', 'ClampFree_Buckle.odb')    # 'odb_direct' 가 읽는 ODB
MODE_SOURCE_STEP_NAME = 'Step-Buckle'                           # 그 ODB 안의 좌굴 스텝 이름
MODE_INSTANCE = 'MEMBRANE-1'                                    # 막 인스턴스 이름

# ---- 모드 소스 오버라이드 (2026-10-07, C-route) ----
#   기본 소스는 **클램프 없는** 모델(ClampFree_Buckle)이다 — 클램프 모델의 좌굴이 오래 실패했기
#   때문(음수 고유값 598~2897 / CONVERGED=0 -> 모드 0개). 이제 C-route(run_abaqus_buckle.py)가
#   클램프를 **포함한** 좌굴모드를 성공적으로 낸다. 주름이 클램프 기인이라는 54배 규명과 맞물려
#   그쪽이 더 나은 시드다. 경로는 CLI 숫자 파서와 성격이 달라 **환경변수**로 받는다
#   (숫자 인자 파싱을 건드리지 않는다).
#     PowerShell:  $env:SOLARSAIL_MODE_ODB='buckle\Buckle_xc050_d018um_dc006_f050_s4_grid_clamp_lf_normal_ex210.odb'
#   .dat(λ 표 출처)는 **ODB 옆 동명 파일**이 기본이고 SOLARSAIL_MODE_DAT 로 덮을 수 있다.
#   경로는 abspath 로 고정한다(스크립트가 상대경로를 자기 기준으로 해석하는 혼동 방지).
_SRC_ENV = os.environ.get('SOLARSAIL_MODE_ODB')
if _SRC_ENV:
    #   [2026-10-08 실측] 이 스크립트의 cwd 는 **code\aba** 다(canary: cwd=...\code\aba).
    #   그래서 사용자가 code 에서 준 상대경로(buckle\...)를 그대로 abspath 하면
    #   code\aba\buckle\... 이 되어 'ODB 없음' 으로 죽는다.
    #   => 후보를 (1) 준 그대로 (2) **스크립트가 있는 폴더(=code)** (3) cwd 순으로 만들어
    #      실제로 존재하는 것을 쓴다. 절대경로면 그대로 통과한다.
    _srcdirs = []
    try:
        _srcdirs.append(os.path.dirname(os.path.abspath(__file__)))      # code\
    except Exception:
        pass
    _srcdirs.append(os.getcwd())                                         # 안전망
    _pick = None
    for _d in [None] + _srcdirs:
        _c = os.path.abspath(_SRC_ENV if _d is None else os.path.join(_d, _SRC_ENV))
        if os.path.exists(_c):
            _pick = _c
            break
    MODE_SOURCE_ODB = _pick if _pick else os.path.abspath(_SRC_ENV)
    if _pick is None:
        print("[IMPERFECTION] 경고: 모드 소스 ODB 를 찾지 못했습니다 -> %s" % MODE_SOURCE_ODB)
    MODE_SOURCE_DAT = os.path.abspath(os.environ.get(
        'SOLARSAIL_MODE_DAT', os.path.splitext(MODE_SOURCE_ODB)[0] + '.dat'))
    IMPERFECTION_NAME = os.path.splitext(os.path.basename(MODE_SOURCE_ODB))[0]
    print("[IMPERFECTION] 모드 소스 override (env SOLARSAIL_MODE_ODB)")
    print("               ODB = %s" % MODE_SOURCE_ODB)
    print("               DAT = %s" % MODE_SOURCE_DAT)
    print("               IMPERFECTION_NAME = %s" % IMPERFECTION_NAME)
    sys.stdout.flush()
#   (출처 ODB/스텝은 상수로 중복 기재하지 않고 모드표 헤더에서 읽어 로그에 남긴다)
if IMPERFECTION_MODE not in ('odb_direct', 'odb_table', 'file'):
    raise RuntimeError("IMPERFECTION_MODE 는 'odb_direct' / 'odb_table' / 'file' 중 하나여야 합니다"
                       " (현재 %r)" % (IMPERFECTION_MODE,))

# 모드표는 모델을 만들기 전에 읽는다(없으면 HF 잡을 제출하지 않고 즉시 중단).
_pert = None
_mode_table = None
_labels = ()
def _available_mode_numbers(path):
    """모드표에 실제로 들어 있는 모드 번호를 오름차순으로 돌려준다(없으면 빈 리스트)."""
    nums = []
    if not os.path.exists(path):
        return nums
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            if s.upper().startswith('MODE '):
                try:
                    nums.append(int(s.split()[1]))
                except (IndexError, ValueError):
                    pass
    return sorted(nums)


if IMPERFECTION_MODE == 'odb_table':
    # 모드표가 정본이다: 표에 있는 모드만 주입할 수 있으므로 상수를 표에 맞춰 조정한다.
    _avail = _available_mode_numbers(MODE_TABLE)
    if _avail:
        _want = tuple(_avail[:IMPERFECTION_MAX_MODES])
        if _want != tuple(IMPERFECTION_MODES):
            emit("[IMPERFECTION] 모드표의 모드 %s -> 주입 모드 %s 로 조정 (상수 %s, 상한 %d)"
                  % (_avail, list(_want), list(IMPERFECTION_MODES), IMPERFECTION_MAX_MODES))
        IMPERFECTION_MODES = _want
    try:
        _mode_table = load_mode_table(MODE_TABLE, IMPERFECTION_MODES)
    except ImperfectionSourceError as _imp_err_tab:
        emit("!!! ERROR: 임퍼펙션 모드표를 읽지 못했습니다 -> HF 잡을 제출하지 않고 중단합니다.")
        emit("[ERROR-EN] imperfection mode table missing/unreadable -> aborting BEFORE HF job submit.")
        emit(str(_imp_err_tab))
        sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
        sys.exit(1)
    _tinfo = mode_table_info(MODE_TABLE)
    emit("[IMPERFECTION] 모드표 %s 로드 완료: 모드 %s / 모드당 노드 %d개"
          % (MODE_TABLE, sorted(_mode_table.keys()), len(_mode_table[IMPERFECTION_MODES[0]])))
    emit("[IMPERFECTION] 모드표 출처: source=%s / step=%s / instance=%s"
          % (_tinfo.get('source', '?'), _tinfo.get('step', '?'), _tinfo.get('instance', '?')))
    # 모드표 재사용 게이트 — 소스 상수를 바꾼 뒤 재추출하지 않으면 '다른 형상의 모드'를 조용히 주입한다(§5b).
    #   차단하지 않고 기록만 한다(mode_table_check.txt).
    try:
        from aba_imperfection import model_fingerprint, table_fingerprint
        _src_mode = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 'run_abaqus_mode.py')
        _fp_now = model_fingerprint(_src_mode)[0]
        _fp_tbl = table_fingerprint(MODE_TABLE)
        if _fp_tbl is None:
            _fp_msg = ("[MODE-TABLE] 지문 없음(표가 구버전) -> 지금 표가 이 모드 소스에서 나온 "
                       "것인지 확인할 수 없다. 의심되면 다시 뽑을 것: "
                       "abaqus cae noGUI=run_abaqus_mode.py")
        elif _fp_now == _fp_tbl:
            _fp_msg = ("[MODE-TABLE] 지문 일치 OK sha1=%s -> 모드 소스가 표 생성 이후 바뀌지 "
                       "않았다(재사용 안전)" % _fp_tbl)
        else:
            _fp_msg = ("[MODE-TABLE] !! STALE: 지문 불일치(현재 소스 %s / 표 %s) -> 모드 소스가 "
                       "바뀌었다. 좌굴모드를 다시 뽑아야 한다: "
                       "abaqus cae noGUI=run_abaqus_mode.py" % (_fp_now, _fp_tbl))
        emit(_fp_msg)
        try:
            with open('mode_table_check.txt', 'w') as _fpc:
                _fpc.write('%s\ntable=%s\nsource_sha1=%s\ntable_sha1=%s\ntime=%s\n'
                           % (_fp_msg, os.path.abspath(MODE_TABLE), _fp_now, _fp_tbl,
                              time.strftime('%Y-%m-%d %H:%M:%S')))
        except Exception:
            pass
    except Exception:
        pass
elif IMPERFECTION_MODE == 'odb_direct':
    # (2026-09-28, 사용자 제안) 모드표 txt 를 거치지 않고 ODB 에서 직접 모드를 읽는다.
    #   물리는 odb_table 과 완전히 같다(노드 좌표 섭동). 파일이 하나 줄 뿐이다.
    try:
        from aba_mode_from_odb import extract_modes
        _tbl, _meta, _picks, _lam, _msgs = extract_modes(
            MODE_SOURCE_ODB, MODE_SOURCE_STEP_NAME, MODE_INSTANCE,
            len(IMPERFECTION_MODES), verbose=True,
            dat_hint=MODE_SOURCE_DAT)   # λ 는 .dat MODE NO 표에서 읽는다(frameValue 는 모드번호)
        #   MODE_SOURCE_DAT 는 코드\aba 기준 '..\ClampFree_Buckle.dat' (모드 소스와 같은 폴더)
    except Exception as _e_odb:
        emit("!!! ERROR: ODB 에서 좌굴모드를 읽지 못했습니다 -> HF 잡을 제출하지 않고 중단합니다.")
        emit("[ERROR-EN] cannot read modes directly from ODB (%s) -> aborting BEFORE HF job submit."
             % MODE_SOURCE_ODB)
        emit("%s: %s" % (type(_e_odb).__name__, _e_odb))
        sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
        sys.exit(1)
    _want = [m for m in IMPERFECTION_MODES if m in _tbl]
    if not _want:
        emit("!!! ERROR: ODB 에서 요청 모드 %s 를 얻지 못했습니다 (ODB 모드: %s) -> 중단."
             % (list(IMPERFECTION_MODES), sorted(_tbl)))
        emit("[ERROR-EN] requested modes not available in ODB -> aborting BEFORE HF job submit.")
        sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
        sys.exit(1)
    if tuple(_want) != tuple(IMPERFECTION_MODES):
        emit("[IMPERFECTION] ODB 모드 %s -> 주입 모드 %s 로 조정"
             % (sorted(_tbl), _want))
    IMPERFECTION_MODES = tuple(_want)
    _mode_table = dict((m, _tbl[m]) for m in IMPERFECTION_MODES)
    emit("[IMPERFECTION] ODB 직접 읽기 완료: %s / %s / %s -> 모드 %s"
         % (MODE_SOURCE_ODB, MODE_SOURCE_STEP_NAME, MODE_INSTANCE, list(IMPERFECTION_MODES)))
    emit("[IMPERFECTION] 모드당 노드 %d개 / λ %s"
         % (len(_mode_table[IMPERFECTION_MODES[0]]),
            ['%.6e' % _lam[_picks[m - 1]] for m in IMPERFECTION_MODES]))
# A: 추출 요청 고유값 수. base state 가 부정정이면 줄이는 것이 subspace 수렴에 유리하다(§8a).
N_EIG_BUCKLE = 100
# subspace 반복의 기저 벡터 수 (구 MFBO_VECTORS 환경변수는 2026-09-21 제거됨)
BUCKLE_VECTORS = 250

MODEL_NAME = 'SailModel_Triangle'
INSTANCE_NAME = 'MEMBRANE-1'

BASE = 20.0   # m
HEIGHT = 10.0 # m
# 막 요소는 1차(S4/S3) 유지 — 2차는 base state 를 더 나쁘게 만들어 기각됐다(§4).
#   요소코드는 모드 소스와 같아야 한다(*IMPERFECTION 이 노드 라벨로 주입된다).
ELEM_CODE_QUAD = S4         # [2026-10-04] 막(M3D4/M3D3) 실험 철회 -> 1차 셸 복귀.
ELEM_CODE_TRI = S3          # [2026-10-04] 위와 같은 이유로 3절점 1차 셸 복귀
SEED_DIV = 200.0           # seed = BASE/SEED_DIV -> 약 1.82만 요소 (실측 2026-09-28)
THICKNESS = 5.0e-6 # F2: 2.5e-6 -> 5.0e-6 (cable 변형, 2.5um는 수렴 매우 어려움)
TARGET_STRESS = 7000.0 # Pa   # R-13: 목표 운용점 - 실제 도달 응력 미검증(측정 필요)

CLAMP_EXCL_R = 0.0   # m — 0.0 이면 제외 노드가 없어 All_Edges_NoClamp == All_Edges 다.
#   [2026-10-06 철회] 0.2 로 두면 클램프 부착 구간의 경계 노드가 면외 자유로워져 base state 가
#   불안정해진다(실측: lambda 전부 음수, CONVERGED 0, 스프레드 0.406 %). 이전에 CONVERGED=100
#   이었던 설정과 같게 0.0 으로 되돌린다. 근거와 실측은 §9.
                     #     이 반경 내 경계 노드는 면외 z 구속에서 제외한다(§9).
                     #     근거: 논문 (c) 는 클램프 없는 모델의 조건이고, 클램프가 붙은 선을 면외 고정하면 클램프의 물리가 왜곡된다.

# 케이블 파라미터 (Galhofo reference)
CABLE_RADIUS = 5.0e-4 # m
CABLE_AREA = np.pi * (CABLE_RADIUS**2)
LEN_TOP = 0.280 # m
LEN_BOT = 0.689 # m

CLAMP_CABLE_LEN = 0.5 # m

# 좌표 정의
V1 = (BASE/2.0, HEIGHT, 0.0) # Top
V2 = (BASE, 0.0, 0.0)        # Right
V3 = (0.0, 0.0, 0.0)         # Left

def clamp_coord_L(x): return (10-10*x, 10-10*x, 0)
def clamp_coord_R(x): return (10+10*x, 10-10*x, 0)

V_CL = clamp_coord_L(x_c)
V_CR = clamp_coord_R(x_c)

# ---- 
# ---- 사전 장력(prestrain) 캘리브레이션 ----
# 사전 장력 크기. base state 장력을 키워 시스템행렬 부정정(음수 고유값) 완화를 시도한다.
# 1.0 = 기존값. 값을 바꾸려면 이 상수를 직접 수정한다(구 MFBO_PRETENSION_SCALE 환경변수는 제거됨).
PRETENSION_SCALE = 10.0          # 프리텐션 변위 = 5e-6 m * 이 값 = 5e-5 m
# PRETENSION_SCALE=10 -> 면내응력 2122 Pa = 목표 7000 Pa 의 0.303배. 운용점 조정은 DISP_GLOBAL 로 한다(§8).
DISP_GLOBAL = 165e-6
CLAMP_PULL = DISP_GLOBAL * d_c
# 좌굴 스텝 perturbation — K_delta 를 만든다. 너무 작으면 고유값 분리가 나빠져 추출이 실패한다(§8b).
PERTURBATION = 0.01
CLAMP_PERT = PERTURBATION * d_c
GLOBAL_FINAL = 1e-3  # 최종 하중: 코너 당김 1e-3 m
CLAMP_FINAL = GLOBAL_FINAL * d_c


# 하중 각도 (28.6도)
angle_deg = 28.6
angle_rad = np.deg2rad(angle_deg)
cos_val = float(np.cos(angle_rad))
sin_val = float(np.sin(angle_rad))

# 모델 초기화 및 재질
def _apply_imperfection_to_part(part):
    """모드표 -> 파트 노드 좌표 섭동(*IMPERFECTION 과 같은 합 규약). (pert, labels) 반환.

    [2026-10-08] 원래 모델 구축부의 `if not USE_GRID_MESH:` **안에 인라인**으로 있었다.
    격자(aba_grid_mesh.PartFromInputFile) 경로로 바꾼 뒤 그 블록이 죽은 코드가 되어
    섭동이 조용히 사라졌다: _pert=None -> 로그 '빈 섭동' -> HF 가 임퍼펙션 없이 제출될 뻔했다.
    두 경로가 같은 함수를 부르도록 뽑아냈다(중복 금지).
    """
    _amp = THICKNESS * IMPERFECTION_AMPL_T
    _pert = build_perturbation(_mode_table, _amp, IMPERFECTION_MODES)
    _labels = tuple(sorted(_pert.keys()))
    _apply, _new = [], []
    for _nd in part.nodes.sequenceFromLabels(labels=_labels):
        _dz = _pert.get(_nd.label)
        if _dz is None:
            continue
        _cx, _cy, _cz = _nd.coordinates
        _apply.append(_nd)
        _new.append((_cx, _cy, _cz + _dz))
    #   MeshNode 에는 setValues(coordinates=...) 가 없다 — Part.editNode 로 한 번에 넘긴다(§6b).
    part.editNode(nodes=tuple(_apply), coordinates=tuple(_new))
    emit("[IMPERFECTION] 기하 섭동 적용: %d/%d 노드, %s (진폭 %.2f t = %.3e m)"
         % (len(_apply), len(_labels), perturbation_report(_pert),
            IMPERFECTION_AMPL_T, _amp))
    if len(_apply) != len(_labels):
        emit("!!! WARNING: 모드표 노드 %d개 중 %d개만 적용 -> 메쉬/라벨 불일치 의심"
             % (len(_labels), len(_apply)))
    #   말이 아니라 실측: 좌표에 실제로 들어갔는지 샘플 3개를 로그에 남긴다.
    for _nd in part.nodes.sequenceFromLabels(
            labels=(_labels[0], _labels[len(_labels) // 2], _labels[-1])):
        emit("               파트 노드 %d z=%.9e (dz=%.3e m)"
             % (_nd.label, _nd.coordinates[2], _pert.get(_nd.label, 0.0)))
    return _pert, _labels


if MODEL_NAME in mdb.models: del mdb.models[MODEL_NAME]
my_model = mdb.Model(name=MODEL_NAME)

# 멤브레인 (Kapton)
mat = my_model.Material(name='Kapton')
mat.Density(table=((1420.0, ),))
mat.Elastic(table=((2.5e9, 0.34),))
my_model.HomogeneousShellSection(name='Section-Membrane', material='Kapton', thickness=THICKNESS)   # [2026-10-04] 막 실험 철회로 셸 섹션 복귀

# 케이블 (Kevlar)
mat_cable = my_model.Material(name='Kevlar')
mat_cable.Density(table=((1440.0, ),)) 
mat_cable.Elastic(table=((62.0e9, 0.36),))
my_model.TrussSection(name='Section-Cable', material='Kevlar', area=CABLE_AREA)

# 파트 생성: 멤브레인
# ---- 메쉬 종류 (2026-10-07) ----
#   True  = aba_grid_mesh 균일 격자(.inp -> PartFromInputFile, orphan mesh). 좌우 완전 대칭.
#   False = 기존 자유 메쉬(seedPart + QUAD_DOMINATED/FREE/MEDIAL_AXIS, 약 1.82만 요소).
#   !! 네 스크립트(buckle/mode/HF/cable)가 **같은 값**을 써야 한다.
USE_GRID_MESH = True
#   orphan(격자) 경로는 한 이름이 노드용/요소용을 겸할 수 없다(2026-10-07 실측:
#     ***ERROR: Unknown part instance node set MEMBRANE-1.ALL).
#     All      = **절점**집합 -> 변위 BC.  All_Elem = **요소**집합 -> 섹션 + 초기응력.
_ALL_ELEM = 'All_Elem' if USE_GRID_MESH else 'All'
if USE_GRID_MESH:
    # 격자: `.inp` 로 **orphan mesh part** 를 import 한다.
    #   Part.addNodes / addElements / deleteMesh 는 odb.Part 전용이라 mdb Part 에는 없다
    #   (실측 AttributeError 2026-10-06/07). CAE 스크립트로 격자를 넣는 경로는
    #   mdb.models[].PartFromInputFile 뿐이다.
    p, _nn, _n4, _n3 = aba_grid_mesh.import_grid_part(
        my_model, part_name='Membrane', base=BASE, height=HEIGHT, seed_div=SEED_DIV)
    print("[mesh] grid(.inp import, orphan): 노드 %d / S4 %d + S3 %d = %d"
          % (_nn, _n4, _n3, _n4 + _n3))
    #   orphan mesh: 섹션은 **요소 기반**으로, 요소가 있는 **뒤**에 준다(§10a).
    #   셋 이름 'All' 유지 — inst_memb.sets['All'] 참조가 남아 있다.
    #   orphan mesh 는 요소/절점 겸용 집합이 없다 -> 두 개를 만든다.
    p.Set(elements=p.elements, name='All_Elem')   # 섹션+초기응력(요소)
    p.Set(nodes=p.nodes, name='All')              # 변위 BC(절점) — 이름 유지
    p.SectionAssignment(region=p.sets['All_Elem'], sectionName='Section-Membrane')
else:
    s = my_model.ConstrainedSketch(name='triangle_profile', sheetSize=BASE*2)
    s.Line(point1=V3[:2], point2=V2[:2])
    s.Line(point1=V2[:2], point2=V1[:2])
    s.Line(point1=V1[:2], point2=V3[:2])
    p = my_model.Part(name='Membrane', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseShell(sketch=s)
    p.SectionAssignment(region=p.Set(faces=p.faces, name='All'), sectionName='Section-Membrane')

# 케이블 생성
def create_cable_part(name, length):
    p_c = my_model.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p_c.WirePolyLine(points=((0.0, 0.0, 0.0), (length, 0.0, 0.0)), mergeType=IMPRINT, meshable=ON)
    p_c.SectionAssignment(region=p_c.Set(edges=p_c.edges, name='Wire'), sectionName='Section-Cable')
    # 메쉬 (T3D2)
    p_c.seedPart(size=length) # 요소 1개
    elemTypeTruss = ElemType(elemCode=T3D2, elemLibrary=STANDARD)
    p_c.setElementType(regions=(p_c.edges,), elemTypes=(elemTypeTruss,))
    p_c.generateMesh()
    return p_c

p_cable_top = create_cable_part('Cable_Top', LEN_TOP)
p_cable_bot_l = create_cable_part('Cable_BotL', LEN_BOT)
p_cable_bot_r = create_cable_part('Cable_BotR', LEN_BOT)
p_cable_cl = create_cable_part("Cable_CL", CLAMP_CABLE_LEN)
p_cable_cr = create_cable_part("Cable_CR", CLAMP_CABLE_LEN)

a = my_model.rootAssembly
a.DatumCsysByDefault(CARTESIAN)
inst_memb = a.Instance(name=INSTANCE_NAME, part=p, dependent=ON)

# --- Rigid Patch 생성  ---
def create_rigid_patch(name, coord, radius):
    """
    지정된 좌표 기준 radius 내의 노드들을 묶어 '분산(distributing)' 결합한다
    실제 solar sail 에서 케이블이나 클램프를 설치하기 위해 sail 면에 테이프 등을 설치하는 과정을 모사
    """
    # RP 생성
    rp = a.ReferencePoint(point=coord)

    rp_key = a.referencePoints[rp.id]
    rp_region = a.Set(name=name+'_RP_Set', referencePoints=(rp_key,))
    
    # 노드 찾기
    nodes = inst_memb.nodes.getByBoundingSphere(center=coord, radius=radius)
    if not nodes:
        nodes = (inst_memb.nodes.getClosest(coordinates=coord),)
    
    patch_set = a.Set(name=name+'_Nodes', nodes=nodes)
    
    # Coupling — DISTRIBUTING 은 표면 절점을 강체로 묶지 않고 하중/변위만 가중 분배한다(§7).
    #   KINEMATIC+WHOLE_SURFACE 는 패치 경계에 응력 특이점(52배)을 만들었다. influenceRadius 는 API 필수 인자.
    my_model.Coupling(
        name=name+'_Coupling', controlPoint=rp_region, surface=patch_set, 
        influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, weightingMethod=UNIFORM,
        u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON
    )
    return rp, rp_region

# --- 케이블 배치 및 연결 ---
def connect_cable(name, part, coord, vector_dir):
    """케이블 part를 불러와 회전/이동 후 원하는 좌표에 배치하는 함수"""
    inst_name = 'Inst_' + name
    inst = a.Instance(name=inst_name, part=part, dependent=ON)
    
    # 회전축 계산 (Z축 회전)
    target_vec = np.array(vector_dir)
    target_vec = target_vec / np.linalg.norm(target_vec) # Normalize
    
    # 2D 회전각 계산 (atan2)
    rot_angle_deg = np.degrees(np.arctan2(target_vec[1], target_vec[0]))
    
    # 회전 및 이동
    a.rotate(instanceList=(inst_name,), axisPoint=(0,0,0), axisDirection=(0,0,1), angle=rot_angle_deg)
    a.translate(instanceList=(inst_name,), vector=coord)
    
    a.regenerate()

    region_start = a.Set(name=name.upper()+'_STARTSET', nodes=inst.nodes[0:1])
    region_end   = a.Set(name=name.upper()+'_ENDSET',   nodes=inst.nodes[1:2])

    return region_start, region_end

if not USE_GRID_MESH:
    p.seedPart(size=BASE/SEED_DIV, deviationFactor=0.1)  # 약 1.82만개 (1차 요소)
    p.setMeshControls(regions=p.faces, elemShape=QUAD_DOMINATED, technique=FREE, algorithm=MEDIAL_AXIS)
    elemTypeQuad = ElemType(elemCode=ELEM_CODE_QUAD, elemLibrary=STANDARD)
    elemTypeTri = ElemType(elemCode=ELEM_CODE_TRI, elemLibrary=STANDARD)
    p.setElementType(regions=(p.faces,), elemTypes=(elemTypeQuad, elemTypeTri))
    p.generateMesh()

    a.regenerate()

# 의존 인스턴스는 파트 메쉬를 공유한다 -> 어셈블리 쪽에서도 좌표가 같아야 한다(전달 확인).
# 격자/자유 **두 경로 모두** 메쉬가 준비된 뒤에 한 번만 적용한다(2026-10-08 수정).
#   이 호출이 없으면 _pert=None 이 되어 로그에 '빈 섭동'만 남고 HF 가 평탄한 막으로 돌아간다.
if _mode_table:
    _pert, _labels = _apply_imperfection_to_part(p)
else:
    _pert, _labels = None, ()
    emit("[IMPERFECTION] 기하 섭동 없음 (모드 없음: IMPERFECTION_MODE=%s, table=%s)"
         % (IMPERFECTION_MODE, type(_mode_table).__name__))

a.regenerate()   # 파트 메쉬 좌표 변경을 (의존) 인스턴스로 전파(§6a)

if _pert is not None and _labels:
    _i0 = _labels[0]
    emit("[IMPERFECTION] 어셈블리 인스턴스 확인: 노드 %d z=%.9e (파트와 같아야 함)"
          % (_i0, inst_memb.nodes.sequenceFromLabels(labels=(_i0,))[0].coordinates[2]))

# 꼭짓점 RP
rp1_obj, rp1_reg = create_rigid_patch('Top', V1, radius=0.2)
rp2_obj, rp2_reg = create_rigid_patch('Right', V2, radius=0.2)
rp3_obj, rp3_reg = create_rigid_patch('Left', V3, radius=0.2)

# 클램프 RP: 우측 빗변 중점 (15,5), 좌측 빗변 중점 (5,5). NO_CLAMP=True 면 클램프를 아예 만들지 않는다.
NO_CLAMP = False          # True = 클램프 생략 진단 모델 (run_abaqus_cable 대조용)
# 실행 모드.
#   'auto'  : 덱 생성 -> riks_patch_input.py 로 *Controls 주입 -> 제출 -> 결과 추출까지 한 번에.
#             mfbo.py 가 이 경로로 HF 를 부른다(덱 생성만 하면 RESULTS: 가 없어 실패로 읽힌다).
#   'write' : 덱만 생성하고 종료. 수동으로 패처/job 을 돌리는 워크플로우(latex/디버깅용).
#   'plain' : 패처 없이 원본 덱을 그대로 제출. 패처 자체를 의심할 때만 쓴다(진단용, 40분 낭비 주의).
RUN_MODE = 'auto'
# 초기 가짜 응력(수렴 보조). 케이블 변형=700 Pa, 우리=500 Pa -> 정렬 노브
SIGMA0 = 500.0                   # 수렴 보조용 초기응력 [Pa]
if NO_CLAMP:
    emit("[run_abaqus] NO_CLAMP=True : 클램프(cable_CL/CR + 강체패치 + BC) 없이 모델링 (대조 실험)")
    rp_cl_obj = rp_cr_obj = rp_cl_reg = rp_cr_reg = None
else:
    rp_cl_obj, rp_cl_reg = create_rigid_patch('CL', V_CL, radius=0.2)
    rp_cr_obj, rp_cr_reg = create_rigid_patch('CR', V_CR, radius=0.2)

# 케이블 연결
# V1(Top): 위로 (0, 1)
start_c1, end_c1 = connect_cable('Cable_Top', p_cable_top, V1, (0.0, 1.0, 0.0))
# V2(Right): 우하향 (+cos, -sin)
start_c2, end_c2 = connect_cable('Cable_Right', p_cable_bot_r, V2, (cos_val, -sin_val, 0.0))
# V3(Left): 좌하향 (-cos, -sin)
start_c3, end_c3 = connect_cable('Cable_Left', p_cable_bot_l, V3, (-cos_val, -sin_val, 0.0))
# V_CL : 좌상향 (-1, 1)
if NO_CLAMP:
    start_cl = end_cl = start_cr = end_cr = None
else:
    start_cl, end_cl = connect_cable('cable_CL', p_cable_cl, V_CL, (-1.0, 1.0, 0.0))
    # V_CR : 우상향 (1, 1)
    start_cr, end_cr = connect_cable('cable_CR', p_cable_cr, V_CR, (1.0, 1.0, 0.0))

a.Set(name='RP_Top_Set', referencePoints=(a.referencePoints[rp1_obj.id],))
a.Set(name='RP_Right_Set', referencePoints=(a.referencePoints[rp2_obj.id],))
a.Set(name='RP_Left_Set', referencePoints=(a.referencePoints[rp3_obj.id],))

if not NO_CLAMP:
    a.Set(name='RP_CL_Set', referencePoints=(a.referencePoints[rp_cl_obj.id],))
    a.Set(name='RP_CR_Set', referencePoints=(a.referencePoints[rp_cr_obj.id],))

my_model.Tie(name='Tie_Top', main=a.sets['RP_Top_Set'], secondary=start_c1, positionToleranceMethod=COMPUTED)
my_model.Tie(name='Tie_Right', main=a.sets['RP_Right_Set'], secondary=start_c2, positionToleranceMethod=COMPUTED)
my_model.Tie(name='Tie_Left', main=a.sets['RP_Left_Set'], secondary=start_c3, positionToleranceMethod=COMPUTED)

if not NO_CLAMP:
    my_model.Tie(name='Tie_CL', main=a.sets['RP_CL_Set'], secondary=start_cl, positionToleranceMethod=COMPUTED)
    my_model.Tie(name='Tie_CR', main=a.sets['RP_CR_Set'], secondary=start_cr, positionToleranceMethod=COMPUTED)

a.regenerate()

# Step Global Tension : 3개 꼭짓점에 변위 가하기
my_model.StaticStep(
    name='Step-GlobalTension',
    previous='Initial', 
    nlgeom=ON,
    stabilizationMagnitude=0.0002,      # Galhofo Reference
    stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
    initialInc=0.0001, minInc=1e-8, maxNumInc=1000    # P1-2: 1e-15 는 발산 시 증분 폭주
)

# Step Clamp Tension : 클램프에 변위 가하기
my_model.StaticStep(
    name='Step-ClampTension',
    previous='Step-GlobalTension',
    nlgeom=ON,
    stabilizationMagnitude=0.0002,      # Galhofo Reference
    stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
    initialInc=0.0001, minInc=1e-8, maxNumInc=1000    # P1-2: 1e-15 는 발산 시 증분 폭주
)

# Step Buckle : 버클 모드 찾기. LF에서는 직접적으로 쓰이지 않음. 
# HF 포스트버클링의 imperfection으로 사용됨.
#Galhofo 참조 BUCKLE+subspace, 최초 4모드.
# BuckleStep도 아래 주석처리로 남겨놓았음.

"""if 'Step-Buckle' in my_model.steps: del my_model.steps['Step-Buckle']
my_model.BuckleStep(
    name='Step-Buckle',
    previous='Step-ClampTension',
    numEigen=20,             
    eigensolver=LANCZOS,     
    maxBlocks=DEFAULT        
)"""

if 'Step-Buckle' in my_model.steps: del my_model.steps['Step-Buckle']
# A: 좌굴 고유값 추출 솔버. **기본 LANCZOS** — Abaqus 의 *BUCKLE 기본 솔버이자 이 코드의
#    원래 설정(LANCZOS + maxBlocks=DEFAULT)이다. 요소 ~1만 개 / DOF ~4만 개 + 두께 5um 라
#    강성 조건수가 극단적이어서 SUBSPACE(작은 모델·다수 모드용)로는 100개 모드 추출이
#    'EIGENVALUES CANNOT BE FOUND' (0 CONVERGED) 로 실패했다.
#    SUBSPACE 로 강제하려면 _EIGENSOLVER 상수를 'SUBSPACE' 로 수정 (구 MFBO_EIGENSOLVER 환경변수는 제거됨)
_EIGENSOLVER = "LANCZOS"         # "SUBSPACE" 로 바꾸면 좌굴 추출 재시도
emit("[run_abaqus] buckle eigensolver=%s numEigen=%d vectors=%s"
      % (_EIGENSOLVER, N_EIG_BUCKLE, (BUCKLE_VECTORS if _EIGENSOLVER != 'LANCZOS' else 'n/a')))
if _EIGENSOLVER == 'LANCZOS':
    my_model.BuckleStep(
        name='Step-Buckle',
        # (B) base state = '클램프 없는(GlobalTension 말단)' 상태.
        #   ClampTension 말단을 base 로 쓰면 시스템행렬이 부정정(음수 고유값 수천 개)이 되어
        #   좌굴모드 0개 -> 임퍼펙션 seed 없음 -> 포스트버클 불가.
        previous='Step-GlobalTension',
        numEigen=N_EIG_BUCKLE,
        eigensolver=LANCZOS,
        maxBlocks=DEFAULT,
    )
else:
    my_model.BuckleStep(
        name='Step-Buckle',
        previous='Step-GlobalTension',
        numEigen=N_EIG_BUCKLE,
        eigensolver=SUBSPACE,
        vectors=BUCKLE_VECTORS,      # A: 기본 8 -> 250
        maxIterations=5000,
    )
# *IMPERFECTION, STEP=n 의 n 은 'Buckle_Analysis.fil 안의 스텝 번호'
# (현재 2 = GlobalTension/Buckle/ClampTension). 스텝 순서가 바뀌어도 자동 추종.
#   n 은 'Initial 을 제외한' 1-based 스텝 번호 (len() 은 Initial 포함해 1 크다 — NEW-1).
_steps_in_order = [s for s in my_model.steps.keys() if s != 'Initial']
_BUCKLE_STEP_NO = _steps_in_order.index('Step-Buckle') + 1
emit("[run_abaqus] steps=%s  _BUCKLE_STEP_NO=%d" % (_steps_in_order, _BUCKLE_STEP_NO))

# *IMPERFECTION, FILE= 은 results file(.fil) 을 읽는다 -> 좌굴 모드를 .fil 에 기록해야 임퍼펙션이 실제로 주입된다 (미요청 시 조용히 무시됨)

my_model.Stress(
    name='Initial_Stiffness',
    region=inst_memb.sets[_ALL_ELEM],
    distributionType=UNIFORM,
    sigma11=SIGMA0, sigma22=SIGMA0, sigma33=0.0,
    sigma12=0.0, sigma13=0.0, sigma23=0.0
)

# 시작하면서 전체 면 z 변위 고정
my_model.DisplacementBC(
    name='BC_Stabilize_Z', 
    createStepName='Initial', 
    region=inst_memb.sets['All'],
    u3=SET
)

# 구속조건 생성
my_model.DisplacementBC(name='BC_Anchor', createStepName='Initial', region=end_c1, 
                        u1=0, u2=0, u3=0, ur1=0, ur2=0, ur3=0
                        )
my_model.DisplacementBC(name='BC_Right', createStepName='Initial', region=end_c2, u3=0, ur1=0, ur2=0, ur3=0)
my_model.DisplacementBC(name='BC_Left', createStepName='Initial', region=end_c3, u3=0, ur1=0, ur2=0, ur3=0)
if not NO_CLAMP:
    my_model.DisplacementBC(name='BC_CL', createStepName='Initial', region=end_cl, u3=0)
    my_model.DisplacementBC(name='BC_CR', createStepName='Initial', region=end_cr, u3=0)

my_model.DisplacementBC(name='Disp_Control_Right', createStepName='Initial', region=end_c2, 
    u1=0, u2=0
)
my_model.DisplacementBC(name='Disp_Control_Left',createStepName='Initial', region=end_c3, 
    u1=0, u2=0
)
if not NO_CLAMP:
    my_model.DisplacementBC(name='Disp_Control_CL', createStepName='Initial', region=end_cl, 
        u1=0, u2=0
    )
    my_model.DisplacementBC(name='Disp_Control_CR', createStepName='Initial', region=end_cr, 
        u1=0, u2=0
    )

# Right 케이블: 우하향 당기기
my_model.boundaryConditions['Disp_Control_Right'].setValuesInStep(stepName='Step-GlobalTension', 
    u1=DISP_GLOBAL * cos_val,
    u2=-DISP_GLOBAL * sin_val
)
# Left Cable: 좌하향 당기기
my_model.boundaryConditions['Disp_Control_Left'].setValuesInStep(stepName='Step-GlobalTension', 
    u1=-DISP_GLOBAL * cos_val,
    u2=-DISP_GLOBAL * sin_val
)

if not NO_CLAMP:
    # 좌측 클램프: (-1, +1) 방향
    my_model.boundaryConditions['Disp_Control_CL'].setValuesInStep(stepName='Step-ClampTension', 
        u1=-(CLAMP_PULL/sqrt2),
        u2=(CLAMP_PULL/sqrt2)
    )
    # 우측 클램프: (+1, +1) 방향
    my_model.boundaryConditions['Disp_Control_CR'].setValuesInStep(stepName='Step-ClampTension', 
        u1=(CLAMP_PULL/sqrt2),
        u2=(CLAMP_PULL/sqrt2)
    )

# Buckle Perturbation
my_model.boundaryConditions['Disp_Control_Right'].setValuesInStep(stepName='Step-Buckle', 
    u1=PERTURBATION * cos_val,
    u2=-PERTURBATION * sin_val
)
my_model.boundaryConditions['Disp_Control_Left'].setValuesInStep(stepName='Step-Buckle', 
    u1=-PERTURBATION * cos_val,
    u2=-PERTURBATION * sin_val
)
if not NO_CLAMP:
    my_model.boundaryConditions['Disp_Control_CL'].setValuesInStep(stepName='Step-Buckle', 
        u1=-(CLAMP_PERT / sqrt2),
        u2=(CLAMP_PERT / sqrt2)
    )
    my_model.boundaryConditions['Disp_Control_CR'].setValuesInStep(stepName='Step-Buckle', 
        u1=(CLAMP_PERT / sqrt2),
        u2=(CLAMP_PERT/ sqrt2)
    )

# Step Buckle 에서는 전체 z 구속에서 모서리 z 구속으로 교체 Galhofo reference
my_model.boundaryConditions['BC_Stabilize_Z'].deactivate('Step-Buckle')

# all_edges = inst_memb.edges 제거 — orphan mesh 에서는 인스턴스의 기하 edge 가 비어 셋이 무효다(§10a).

# 면외 z 구속용 노드셋(All_Edges_NoClamp) — 클램프 부착 구간을 제외한다(§9).
#   세 변의 방정식은 y=0 / y=x / y=BASE-x. 좌표 판정은 aba_grid_mesh 로 일원화(세 스크립트가 같은 함수).
_excl_edges = (aba_grid_mesh.clamp_exclude_labels(inst_memb, (V_CL, V_CR), CLAMP_EXCL_R)
               if not NO_CLAMP else set())
_keep_labels = aba_grid_mesh.boundary_node_labels(
    inst_memb, base=BASE, height=HEIGHT, exclude_labels=_excl_edges)
emit("[BC] All_Edges_NoClamp 후보: 경계 노드 %d개 (클램프 반경 %.3g m 내 %d개 제외, 전체 막 노드 %d개)"
     % (len(_keep_labels), CLAMP_EXCL_R, len(_excl_edges), len(inst_memb.nodes)))
aba_grid_mesh.make_set(a, 'All_Edges_NoClamp', inst_memb, _keep_labels)
emit("[BC] All_Edges_NoClamp 셋 생성 완료 (노드 %d개)" % len(_keep_labels))

my_model.DisplacementBC(
    name='BC_Edges_Only_Z',
    createStepName='Step-Buckle',
    region=a.sets['All_Edges_NoClamp'],
    u3=0
)

BUCKLE_ODB = 'Buckle_Analysis.odb'
LF_ODB = 'LF_Analysis.odb'
HF_ODB = 'HF_Postbuckle.odb'

cmd = ""
# 버클링 준비 (모델 완성 후)
# R-7: *IMPERFECTION, FILE= 은 results file(.fil) 을 읽으므로 좌굴모드를 .fil 에 기록해야 한다.
# 결함 A/NEW-2: 그 *NODE FILE 은 Step-Buckle 안에서만 유효한데 포스트버클 잡은 그 스텝을
#   삭제하므로 키워드가 스텝 밖으로 밀려나 'the keyword is misplaced' 로 입력처리 직사한다.
#   -> *NODE FILE 은 좌굴 잡 전용 '모델 복사본'에만 넣고, 원본 모델은 건드리지 않는다.
BUCKLE_MODEL = MODEL_NAME
if fidelity == 'HF':
    _noderef = '*NODE FILE, GLOBAL=YES, LAST MODE=%d\nU' % N_EIG
    try:
        BUCKLE_MODEL = 'Model-Buckle'
        if BUCKLE_MODEL in mdb.models:
            del mdb.models[BUCKLE_MODEL]
        mdb.Model(name=BUCKLE_MODEL, objectToCopy=my_model)
        _bkm = mdb.models[BUCKLE_MODEL]
        _bkm.keywordBlock.synchVersions(storeNodesAndElements=False)
        _found = False
        for _i, _b in enumerate(_bkm.keywordBlock.sieBlocks):
            if _b.strip().upper().startswith(('*BUCKLE', '*FREQUENCY')):
                _bkm.keywordBlock.insert(_i + 1, _noderef)
                _found = True
                break
        emit("[run_abaqus] 좌굴 잡 모델=%s, *NODE FILE 삽입=%s" % (BUCKLE_MODEL, _found))
    except Exception as _e:
        emit("!!! WARNING: 모델 복사 실패(%s) -> 원본에 삽입 (결함 A 재발 가능)" % _e)
        BUCKLE_MODEL = MODEL_NAME
        my_model.keywordBlock.synchVersions(storeNodesAndElements=False)
        for _i, _b in enumerate(my_model.keywordBlock.sieBlocks):
            if _b.strip().upper().startswith(('*BUCKLE', '*FREQUENCY')):
                my_model.keywordBlock.insert(_i + 1, _noderef)
                break

if fidelity == 'LF':

    # LF 는 좌굴모드를 쓰지 않는다(HF 분기가 같은 설계점에서 자체 실행) -> 비용 절감
    # run_job_safely('Buckle_Analysis')
    
    if 'Initial_Stiffness' in my_model.predefinedFields:
        del my_model.predefinedFields['Initial_Stiffness']
    
    # 낮은 값으로 다시 생성 (region 은 반드시 요소집합 _ALL_ELEM — 절점집합 'All' 을
    #   주면 'not an element set' 으로 입력처리기가 죽는다. 2026-10-08)
    my_model.Stress(
        name='Initial_Stiffness',
        region=inst_memb.sets[_ALL_ELEM],
        distributionType=UNIFORM,
        sigma11=SIGMA0, sigma22=SIGMA0, sigma33=0.0, 
        sigma12=0.0, sigma13=0.0, sigma23=0.0
    )
    if 'Step-Buckle' in my_model.steps:
        del my_model.steps['Step-Buckle']

    my_model.StaticStep(
        name='Step-HighTension',
        previous='Step-ClampTension',
        nlgeom=ON,
        initialInc=0.01, 
        maxNumInc=100
    )

    my_model.boundaryConditions['BC_Stabilize_Z'].setValuesInStep(
        stepName='Step-HighTension', 
        u3=0.0        
    )

    for name, sign in [('Disp_Control_Right', 1), ('Disp_Control_Left', -1)]:
        my_model.boundaryConditions[name].setValuesInStep(
            stepName='Step-HighTension',
            u1=sign * GLOBAL_FINAL * cos_val,
            u2=-GLOBAL_FINAL * sin_val
        )
    for name, sign in ([] if NO_CLAMP else [('Disp_Control_CL', -1), ('Disp_Control_CR', 1)]):
        my_model.boundaryConditions[name].setValuesInStep(
            stepName='Step-HighTension',
            u1=sign * CLAMP_FINAL / sqrt2,
            u2=CLAMP_FINAL / sqrt2
        )

    my_model.fieldOutputRequests['F-Output-1'].setValues(
        variables=('S', 'E', 'U', 'COORD', 'EVOL'), 
        frequency=10   # R-10: odb 크기 절감
    )

    run_job_safely('LF_Analysis')

    cmd = 'abaqus python "%s" %s LF' % (os.path.join(_HERE, "eval_abaqus.py"), LF_ODB)

elif fidelity == 'HF':

    # 좌굴모드와 포스트버클 해석이 '같은 초기응력 상태'에서 계산되도록
    # 초기응력 재생성(500 -> 100 Pa)을 좌굴 잡 제출 '앞'으로 이동.
    if 'Initial_Stiffness' in my_model.predefinedFields:
        del my_model.predefinedFields['Initial_Stiffness']

    my_model.Stress(
        name='Initial_Stiffness',
        region=inst_memb.sets[_ALL_ELEM],
        distributionType=UNIFORM,
        sigma11=SIGMA0, sigma22=SIGMA0, sigma33=0.0, 
        sigma12=0.0, sigma13=0.0, sigma23=0.0
    )

    if RUN_SELF_BUCKLE_JOB:
        # 자기(클램프) 좌굴 잡: 모드 소스가 아니라 '클램프 base state 증거 + 실패 진단'용이다.
        #   (모드 추출은 MODE_SOURCE='cable' 경로가 담당 -> 클램프 모델은 CONVERGED=0.)
        run_job_safely('Buckle_Analysis', model_name=BUCKLE_MODEL)   # P0-2: 500 Pa 상태에서 좌굴모드 산출
    else:
        emit("[IMPERFECTION] 자기 좌굴 잡 건너뜀 (RUN_SELF_BUCKLE_JOB=False)")

    # 좌굴 모드 병합(coalescence) 진단용 고유치 기록
    #   - 이유: Buckle_Analysis.dat 는 '다음 설계점'의 좌굴 잡이 시작될 때 삭제되므로
    #           여기서 고유치를 뽑아 이력(JSONL)에 남기지 않으면 회고 분석이 불가능하다.
    #   - 기록 전용: 실패해도 해석에는 전혀 영향을 주지 않는다 (예외 전부 삼킴).
    #   - 끄려면 EIG_RECORD 상수를 False 로 (구 MFBO_EIG_RECORD 환경변수는 제거됨)
    try:
        if EIG_RECORD:
            _rec = os.path.join(_HERE, 'coalescence_check.py')
            _dat = _first_existing(os.path.join(os.getcwd(), 'Buckle_Analysis.dat'),
                                   os.path.join(_HERE, 'Buckle_Analysis.dat'),
                                   os.path.join('..', 'Buckle_Analysis.dat'))
            _msg = _first_existing(os.path.join(os.getcwd(), 'Buckle_Analysis.msg'),
                                   os.path.join(_HERE, 'Buckle_Analysis.msg'),
                                   os.path.join('..', 'Buckle_Analysis.msg'))
            _hist = os.path.join(os.getcwd(), 'eig_history.jsonl')
            if os.path.exists(_rec) and (_dat or _msg):
                emit("[N-6] 고유치 소스: dat=%s msg=%s" % (_dat, _msg))
                _cmd = ('abaqus python "%s" record --dat "%s" --msg "%s" --history "%s" '
                        '--x %s --d %s --fidelity HF' % (_rec, _dat, _msg, _hist, x_c, d_c))
                _n6 = subprocess.run(_cmd, shell=True, capture_output=True, text=True)
                # 하위 프로세스 출력을 emit 으로 회수한다. 그냥 두면 콘솔에만 새고
                # hf_run_log.txt 에는 남지 않아 로그와 콘솔이 어긋난다.
                for _ln in ((_n6.stdout or "") + (_n6.stderr or "")).splitlines():
                    _ls = _ln.strip()
                    if _ls.startswith("[record]") and (
                            "n_modes=" in _ls or "건너뜀" in _ls or "실패" in _ls
                            or "없음" in _ls or "찾지 못" in _ls):
                        emit("   %s" % _ls)
                emit("[N-6] coalescence record rc=%d (x_c=%s, d_c=%s)" % (_n6.returncode, x_c, d_c))
            else:
                emit("[N-6] 고유치 기록 건너뜀 (script=%s, dat=%s, msg=%s)"
                      % (os.path.exists(_rec), bool(_dat), bool(_msg)))
    except Exception as _eig_err:
        emit("[N-6] coalescence record 실패(무시): %s" % _eig_err)

    # [R-13 / 스파이크 10-0] base state(GlobalTension 말단) 응력·반력 측정
    #   좌굴 성공/실패와 무관하게 Buckle_Analysis.odb 의 GlobalTension 프레임에서 읽는다.
    #   -> "프리텐션 운용점이 물리적인가"를 실측으로 답하기 위한 것 (미해결 최우선 1건).
    #   끄려면 BASE_PROBE 상수를 False 로 (구 MFBO_BASE_PROBE 환경변수는 제거됨)
    try:
        if BASE_PROBE:
            _probe = os.path.join(_HERE, 'base_state_probe.py')
            _bodb = _first_existing(BUCKLE_ODB,
                                    os.path.join(_HERE, os.path.basename(BUCKLE_ODB)),
                                    os.path.join('..', os.path.basename(BUCKLE_ODB)))
            if os.path.exists(_probe) and _bodb:
                emit("[R-13] base state ODB: %s" % _bodb)
                _pcmd = ('abaqus python "%s" "%s" Step-GlobalTension'
                         % (_probe, _bodb))
                _pr = subprocess.run(_pcmd, shell=True, capture_output=True, text=True)
                # probe 는 40여 줄을 쏟는다(요소 좌표·단면별 통계·RF 상위 8개·인스턴스별 max|u3|).
                # 로그에 필요한 것은 ">>>" 가 붙은 요약뿐이므로 그것만 남긴다.
                _pr_all = ((_pr.stdout or "") + (_pr.stderr or "")).splitlines()
                _n_kept = 0
                for _ln in _pr_all:
                    if ">>>" in _ln:
                        emit("   %s" % _ln.strip()); _n_kept += 1
                emit("[R-13] 요약 %d줄 / 원본 %d줄 (상세 생략)" % (_n_kept, len(_pr_all)))
            else:
                emit("[R-13] base state 측정 건너뜀 (probe=%s, odb=%s : 후보 %s / %s)"
                      % (os.path.exists(_probe), bool(_bodb), BUCKLE_ODB, _HERE))
    except Exception as _probe_err:
        emit("[R-13] base state 측정 실패(무시): %s" % _probe_err)

    # 기존 Step 정리: Post-buckling은 GlobalTension 직후에서 시작하며,
    # 중간 단계(ClampTension)를 건너뛰고 바로 최종 하중으로 Ramping함
    if 'Step-Buckle' in my_model.steps: del my_model.steps['Step-Buckle']
    if 'Step-ClampTension' in my_model.steps: del my_model.steps['Step-ClampTension']

    # Static, General 포스트버클링 수행    
    my_model.StaticStep(
        name='Step-Postbuckle', 
        previous='Step-GlobalTension',  
        nlgeom=ON, 
        stabilizationMagnitude=0.003,      # Galhofo 참조 2e-4 / [2026-09-29] 0.001 -> 0.003 (살짝만)
        #   근거: 실측 ALLSD/ALLIE = 0.50%% 인데 allsdtol = 5%% -> 인공감쇠가 약 10배의 여유를 남기고
        #   거의 작동하지 않았다. 즉 "감쇠가 컷백을 만든다"는 과거 반증은 "감쇠를 늘리면 도움이 될 것"을
        #   반증하지 않는다(다른 주장). allsdtol 이 상한을 걸어 자체 제한되므로 위험은 유계다.
        #   판정: 초기 100~200 증분의 증분 상한이 3e-4 이상으로 회복되는가 (10~20분) + energy_check.py
        stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
        continueDampingFactors=False,
        adaptiveDampingRatio=0.05,
        initialInc=1e-4,
        minInc=1e-8,          # R-9: 1e-15 는 발산 시 증분 소진까지 수시간
        maxInc=0.1,
        maxNumInc=10000       # R-9 / [2026-09-28] 1000 -> 10000
        # 상한 상향은 해결책이 아니라 연장이다 — 실제 원인은 국소 불안정으로 증분이 눌리는 것(§13).
        #   ALLSD/ALLIE 가 0.50%(허용 5%)라 allsdtol 상향은 binding 이 아니다(폐기된 레버).
    )
    my_model.keywordBlock.synchVersions(storeNodesAndElements=False)

    my_model.boundaryConditions['BC_Stabilize_Z'].deactivate('Step-Postbuckle')

    # Step-Buckle 삭제 시 BC_Edges_Only_Z(prescribed condition)가 함께 삭제되므로
    # Postbuckle 스텝에 모서리 z 구속을 재생성한다 (Galhofo 참조: 3개 모서리 u3=0).
    # 셋은 위에서 만든 All_Edges_NoClamp(클램프 구간 제외)를 그대로 쓴다.
    my_model.DisplacementBC(
        name='BC_Edges_Only_Z',
        createStepName='Step-Postbuckle',
        region=a.sets['All_Edges_NoClamp'],
        u3=0
    )

    # 초기 상태에서 최종 장력(목표 7000 Pa, 미검증)까지 증가
    for name, sign in [('Disp_Control_Right', 1), ('Disp_Control_Left', -1)]:
        my_model.boundaryConditions[name].setValuesInStep(
            stepName='Step-Postbuckle',
            u1=sign * GLOBAL_FINAL * cos_val,
            u2=-GLOBAL_FINAL * sin_val
        )

    for name, sign in ([] if NO_CLAMP else [('Disp_Control_CL', -1), ('Disp_Control_CR', 1)]):
        my_model.boundaryConditions[name].setValuesInStep(
            stepName='Step-Postbuckle',
            u1=sign * CLAMP_FINAL / sqrt2,
            u2=CLAMP_FINAL / sqrt2
        )
    

    my_model.fieldOutputRequests['F-Output-1'].setValues(
        variables=('S', 'E', 'U', 'RF', 'COORD', 'EVOL'),   # R-8: LF 지표(lf1~lf3) 산출에 필요
        frequency=10   # R-10: odb 크기 절감
    )
    
    # Imperfection Injection (Keyword Editing)
    # Buckle_Analysis.odb 파일의 결과(고유모드)를 초기 결함으로 주입
    my_model.keywordBlock.synchVersions(storeNodesAndElements=False)
    
    # 2-모델 레시피: 모드 형상은 클램프 없는 모델에서 온다. 소스 .fil 을 별도 이름으로 스테이징해
    #   자기 좌굴 잡의 0-모드 .fil 이 조용히 소비되는 사고를 차단한다.
    imp_scale = THICKNESS * IMPERFECTION_AMPL_T   # Galhofo 채택값 0.10 t
    if IMPERFECTION_MODE in ('odb_table', 'odb_direct'):
        # 이미 모델 빌드 단계에서 노드 좌표를 섭동했다(모드표 txt 또는 ODB 직접). 키워드 경로는 쓰지 않는다.
        imp_name, imp_step = None, None
        _src_desc = MODE_TABLE if IMPERFECTION_MODE == 'odb_table' else (
            "%s / %s" % (MODE_SOURCE_ODB, MODE_SOURCE_STEP_NAME))
        emit("[IMPERFECTION] 주입=기하 섭동 (출처 %s, 모드 %s, 진폭 %.2f t = %.3e m)"
              % (_src_desc, list(IMPERFECTION_MODES), IMPERFECTION_AMPL_T, imp_scale))
        emit("               %s" % perturbation_report(_pert))
    elif MODE_SOURCE == 'self':
        imp_name, imp_step = 'Buckle_Analysis', _BUCKLE_STEP_NO
        emit("[IMPERFECTION] 소스=self 좌굴 잡(model=%s) STEP=%d amplitude=%.3e m"
              % (BUCKLE_MODEL, imp_step, imp_scale))
    else:
        try:
            _st = stage(MODE_SOURCE_FIL, IMPERFECTION_NAME, os.getcwd(),
                        IMPERFECTION_MODES, dat_hint=MODE_SOURCE_DAT,
                        msg_hint=MODE_SOURCE_MSG, here=_HERE)
        except ImperfectionSourceError as _imp_err:
            emit("!!! ERROR: 임퍼펙션 소스 스테이징 실패 -> HF 잡을 제출하지 않고 중단합니다.")
            emit(str(_imp_err))
            sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
            sys.exit(1)
        imp_name, imp_step = _st['name'], MODE_SOURCE_STEP
        emit("[IMPERFECTION] %s STEP=%d amplitude=%.3e m (모드 %s)"
              % (report(_st), imp_step, imp_scale, list(IMPERFECTION_MODES)))
        if _st['modes'] < N_EIG:
            emit("!!! WARNING: 소스 모드 %d개 < N_EIG=%d -> 요청 모드 일부만 주입됩니다."
                  % (_st['modes'], N_EIG))

    if IMPERFECTION_MODE in ('odb_table', 'odb_direct'):
        emit("[IMPERFECTION] 키워드 삽입 생략 (기하 섭동으로 이미 주입됨)")
    else:
        imp_text = imperfection_text(imp_name, imp_step, IMPERFECTION_MODES, imp_scale)

        # 키워드 삽입 위치 찾기 : *STEP 블록 직전에 삽입하는 것이 안전함
        inserted = False
        for i, block in enumerate(my_model.keywordBlock.sieBlocks):
            # Step 정의 시작 부분 찾기
            if block.lower().strip().startswith('*step'):
                my_model.keywordBlock.insert(max(i - 1, 0), imp_text)
                inserted = True
                break

        if not inserted:
            my_model.keywordBlock.insert(len(my_model.keywordBlock.sieBlocks)-1, imp_text)

    run_job_safely('HF_Postbuckle')


    cmd = 'abaqus python "%s" %s %s HF' % (os.path.join(_HERE, "eval_abaqus.py"), LF_ODB, HF_ODB)   # P0-C: 짝지은 LF odb

try:
    # shell=True로 eval_abaqus.py 실행
    # stdout 을 회수해 로그에 남기되 "RESULTS:" 줄은 버린다 — 상위에는 아래
    # extraction.txt 경로로만 전달해야 mfbo.py 가 같은 신호를 두 번 읽지 않는다.
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    rc = p.returncode
    for _ln in ((p.stdout or "") + (p.stderr or "")).splitlines():
        _ls = _ln.strip()
        if _ls.startswith("RESULTS:") or not _ls:
            continue
        if _ls.startswith("!!!") or "ERROR" in _ls or "FAIL" in _ls or "Traceback" in _ls:
            emit("   [eval] %s" % _ls)
    if rc != 0:
        emit("!!! ERROR: eval_abaqus.py exited with code %d" % rc)
    
    # 추출 스크립트가 출력한 "RESULTS:..." 라인을 찾아 전달
    if os.path.exists('extraction.txt'):
        with open('extraction.txt', 'r') as f:
            emit("RESULTS:" + f.read().strip())
    else:
        emit("!!! ERROR: Extraction failed. 'extraction.txt' not found.") 
        emit("RESULTS:FAIL")   # 상위(mfbo)가 원인을 식별하도록 명시적 실패 신호

except Exception as err:
    emit("Error during data extraction: %s" % str(err))
    sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
    sys.exit(1)
