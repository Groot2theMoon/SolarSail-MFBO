"""
run_abaqus_buckle.py — 좌굴(선형 고유값) 해석 전용. 1회 실행, 스윕 없음.

목적
    클램프 패치가 있는 삼각 막 모델의 선형 좌굴 고유모드를 한 번에 추출한다.
    추출한 모드는 HF 포스트버클링의 초기결함(*IMPERFECTION, FILE=..., STEP=2)으로 쓴다.

모델 (run_abaqus_new.py 의 HF 모델과 같은 빌드 블록을 공유한다)
    삼각 막(BASE=20 m, HEIGHT=10 m, 두께 5e-6 m)
    + 꼭짓점 강체패치 3개 + 클램프 강체패치 2개 (위치 x_c).
    케이블은 없다 — 논문 §3.2 의 좌굴 모델과 같은 방식으로, 케이블이 당기던
    지점을 직접 prescribed 변위로 구속한다.
    스텝: Step-GlobalTension(프리텐션) -> Step-Buckle(SUBSPACE, numEigen=100,
    vectors=250). 좌굴 스텝의 perturbation 은 0.01 m (대조 스크립트와 동일).

사용법
    abaqus cae noGUI=run_abaqus_buckle.py -- <x_c> [disp_m]
        x_c     클램프 위치 파라미터. 0.5 -> 좌(5,5) / 우(15,5)
        disp_m  GlobalTension 코너 당김 [m]. 생략하면 기본 5e-5 m.
                (논문 좌굴 모델의 값: 모서리 5e-4 m, 정점 1e-3 m)
    예)  abaqus cae noGUI=run_abaqus_buckle.py -- 0.5
         abaqus cae noGUI=run_abaqus_buckle.py -- 0.5 1e-3

산출물 (code/buckle/)
    <job>.odb / .dat / .msg / .sta / .fil / .diag.txt
    job 이름은 인자+요소+클램프모드에서 자동 생성한다: Buckle_xc<NNN>_d<NNN>um_<elem>_<clamp>
      (예: Buckle_xc050_d0050um_s4r_driven) — 케이스가 바뀌어도 산출물이 서로 덮이지 않는다.
    CLAMP_MODE 는 상단 상수: none | passive | driven | fixed (클램프 존재/작동 분리).
    고유값 표는 .dat 의 MODE NO / EIGENVALUE 블록에 있고, 스크립트가 콘솔에도 덤프한다.

HF 에서 모드를 쓸 때 — 경로 주의 (HF 잡의 작업 디렉터리는 code/aba 다)
    *IMPERFECTION, FILE=..\\buckle\\<job>, STEP=2
    <job> 은 이 스크립트가 마지막에 그대로 찍어 준다.
"""

from abaqus import *
from abaqusConstants import *
from step import *
from mesh import ElemType
import regionToolset
import interaction
import sys
import os
import re
import subprocess
import math
import time
import traceback
import numpy as np

TAG = '[run_abaqus_buckle]'

# ---- 스크립트 위치(_HERE) / Abaqus 작업 디렉터리(_RUN) 결정 ----
# 주의: `abaqus cae noGUI=script.py` 로 실행하면 스크립트가 execfile 로 로드되어
#       __file__ 이 정의되지 않는다(NameError). 아래처럼 우선순위를 둔다.
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
    cand.append(os.path.dirname(os.getcwd()))
    cand.append(os.path.join(os.getcwd(), "code"))
    for c in cand:
        if _ok(c):
            return os.path.abspath(c)
    print("!!! WARNING: eval_abaqus.py 위치를 찾지 못했습니다. cwd=%s" % os.getcwd())
    return os.getcwd()

_HERE = _resolve_here()
# 작업 디렉터리: 좌굴 산출물 전용.
#   HF/LF 산출물(run_abaqus_new.py -> code/aba)과 섞으면 두 가지가 조용히 망가진다:
#     (a) *IMPERFECTION, FILE=... 이 stale .fil 을 읽을 수 있다
#     (b) .dat/.msg/.sta 진단 로그가 어느 해석 것인지 구분되지 않는다
RUN_DIR_NAME = "buckle"
NUMCPUS = 4               # run_abaqus_new.py 와 동일
_RUN = os.path.join(_HERE, RUN_DIR_NAME)
os.makedirs(_RUN, exist_ok=True)
os.chdir(_RUN)
print("%s _HERE = %s" % (TAG, _HERE))
print("%s _RUN  = %s" % (TAG, _RUN))

# aba_grid_mesh: BC 노드셋 헬퍼만 쓴다(clamp_exclude_labels/boundary_node_labels/make_set).
#   격자(fill_part)는 현재 미사용 — Part.addNodes 가 없어 되돌렸다(§10).
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import aba_grid_mesh                                                            # noqa: E402



# ============================================================================
# 모델 상수 — run_abaqus_new.py 의 HF 모델과 같은 값이어야 한다.
#   검사: python code/check_model_consistency.py
# ============================================================================
MODEL_PREFIX = 'SailModel_Buckle'
INSTANCE_NAME = 'MEMBRANE-1'

BASE = 20.0   # m
HEIGHT = 10.0 # m
# 막 요소는 1차(S4/S3) 유지. 2차는 base state 를 더 나쁘게 만들어 기각됐다(§4).
#   요소코드는 HF/mode 와 같아야 한다 — *IMPERFECTION 이 노드 라벨로 주입된다.
ELEM_CODE_QUAD = S4         # [2026-10-04] 막(M3D4/M3D3) 실험 철회 -> 1차 셸 복귀.
ELEM_CODE_TRI = S3          # [2026-10-04] 위와 같은 이유로 3절점 1차 셸 복귀
SEED_DIV = 200.0           # seed = BASE/SEED_DIV -> 약 1.82만 요소 (실측 2026-09-28)
THICKNESS = 5.0e-6

# ---- 케이블 (mode.py / HF 와 동일) ----
CABLE_RADIUS = 5.0e-4 # m
CABLE_AREA = np.pi * (CABLE_RADIUS**2)
LEN_TOP = 0.280 # m
LEN_BOT = 0.689 # m

# ---- 프리텐션 모드 (mode.py 와 동일하게 맞춘다) ----
#   'corner2' = Top 정점 완전고정 + 아래 두 꼭짓점만 당김. mode.py 가 4개 전역 모드를 낸 설정.
#   'paper3'  = 꼭짓점 3개를 모두 당김(논문 좌굴모델). mode.py 실측: 양수 2개/0개.
#   이것이 mode.py 와의 마지막 비(非)클램프 차이였다 — buckle 은 3꼭짓점 구동이었다.
PRETENSION_MODE = 'corner2'
if PRETENSION_MODE not in ('corner2', 'paper3'):
    raise RuntimeError("PRETENSION_MODE must be 'corner2' or 'paper3' (got %r)" % (PRETENSION_MODE,))
DISP_TOP_OVER_CORNER = 1.4142135623   # paper3 전용

CLAMP_EXCL_R = 0.0   # m — 0.0 이면 제외 노드가 없어 All_Edges_NoClamp == All_Edges 다.
#   [2026-10-06 철회] 0.2 로 두면 클램프 부착 구간의 경계 노드가 면외 자유로워져 base state 가
#   불안정해진다(실측: lambda 전부 음수, CONVERGED 0, 스프레드 0.406 %). 이전에 CONVERGED=100
#   이었던 설정과 같게 0.0 으로 되돌린다. 근거와 실측은 §9.
#     이 반경 내 경계 노드는 면외 z 구속에서 제외한다(§9). HF 와 같은 값.

V1 = (BASE/2.0, HEIGHT, 0.0) # Top
V2 = (BASE, 0.0, 0.0)        # Right
V3 = (0.0, 0.0, 0.0)         # Left

def clamp_coord_L(x): return (10-10*x, 10-10*x, 0)
def clamp_coord_R(x): return (10+10*x, 10-10*x, 0)

# ---- 프리텐션 (CLI 두 번째 인자로 덮어쓴다) ----
#   [mode 정합] mode.py 는 1.8e-5 m 를 쓴다. buckle 은 5e-5(=5e-6*10) 였다.
#   이 차이가 클램프와 교란되어 '클램프가 원인'인지 판별할 수 없었다.
DISP_GLOBAL = 1.8e-5    # m — mode.py(ClampFree_Buckle) 와 동일한 운용점

# ---- 클램프 처리 모드 (CLAMP_MODE) ----
#   A-route(run_abaqus.py)의 클램프 '존재'(패치+케이블 Tie, u3=0)와 '작동'(법선 방향
#   CLAMP_PULL 당김)을, C-route 는 케이블 없이 RP 직접 구속/구동으로 재현한다. 4종은 §11.
CLAMP_MODE = 'passive'
CLAMP_DC = 0.5   # d_c — 클램프 당김 비율 (CLAMP_PULL = 코너 당김 * d_c, A-route 와 동일)
# 오타로 조용히 다른 케이스가 되는 것을 막는다 (값 검증은 메쉬 생성 전에).
if CLAMP_MODE not in ('none', 'passive', 'driven', 'fixed'):
    raise RuntimeError("CLAMP_MODE must be 'none'|'passive'|'driven'|'fixed' (got %r)"
                       % (CLAMP_MODE,))

# ---- 좌굴 스텝 (run_abaqus_cable.py 에서 완주가 확인된 설정과 동일) ----
PERTURBATION = 0.01  # m — 좌굴 스텝의 prescribed 변위(증분 응력 -> K_delta).
#   [mode 정합] mode.py 와 HF 는 0.01 을 쓴다. buckle 만 1e-4 였다.
#   [2026-10-06 이력] 1e-4 에서 CONVERGED=100, 0.01 에서 CONVERGED=0 이었던 관측은
#   클램프가 있는 상태에서 얻은 것이다. 이제 나머지 조건을 mode.py 와 맞췄으므로
#   0.01 로 두고 다시 판정한다 — 클램프만 다른 상태에서 비교하려면 이 값도 같아야 한다.
# ---- 추출 창 (래더 L4) ----
#   규칙(mode.py L26): '요청 고유값 수 > base state 의 음수 고유값 수' 여야 양수
#   좌굴모드가 subspace 창에 들어온다. 그런데 클램프(passive)를 켜면 음수 고유값이
#   48 개로 늘고, N_EIG=100 / vectors=250 에서는 ***ERROR: THE EIGENVALUES CANNOT BE
#   FOUND (INSTABILITIES IN THE BASE STATE) 로 죽는다. 클램프를 끄면(none) 같은
#   조건에서 CONVERGED=4 다 — 즉 창이 아니라 음수 개수가 문제다.
#   그래서 창을 넓혀(D) 음수 48 개를 넘겨 본다. 실패 비용이 40 초라 판정이 빠르다.
N_EIG_BUCKLE = 200      # 추출 요청 고유값 수  [2026-10-07] 100 -> 200 (래더 L4)
BUCKLE_VECTORS = 500    # subspace 기저 벡터 수  [2026-10-07] 250 -> 500 (래더 L4)
BUCKLE_MAXITER = 5000
BUCKLE_SOLVER = 'SUBSPACE'   # 'SUBSPACE' | 'LANCZOS' — 제어 흐름용 문자열
BUCKLE_BLOCK_SIZE = 8           # LANCZOS 전용
BUCKLE_MIN_EIGEN = 0.0          # LANCZOS 전용 (음의 고유값도 보고 싶으면 -1e30 등)
BUCKLE_MAX_EIGEN = None         # LANCZOS 전용 (None 이면 인자를 아예 넘기지 않는다)

# Abaqus API 는 심볼릭 상수를 요구한다 — 문자열을 넘기면 즉시 죽는다(§7).
if BUCKLE_SOLVER not in ('SUBSPACE', 'LANCZOS'):
    # 원장 P3: 값 검증이 없으면 dict KeyError 로 죽어 원인이 안 보인다.
    raise RuntimeError("BUCKLE_SOLVER must be 'SUBSPACE' or 'LANCZOS' (got %r)"
                       % (BUCKLE_SOLVER,))
EIGENSOLVER_CONST = {'SUBSPACE': SUBSPACE, 'LANCZOS': LANCZOS}[BUCKLE_SOLVER]

MODE_STABILIZATION = 0.0005   # GlobalTension 안정화 계수. [mode 정합] mode.py 와 동일.
PATTERN_SIGN = 1.0            # 좌굴 '하중 패턴' 부호: 1.0 = 바깥으로 더 당김 (mode.py 동일)

SIGMA0 = 700.0          # 초기응력 [Pa] — 수렴 보조. [mode 정합] mode.py=700 (buckle 은 500 이었다)


# 인자 — 한 번에 한 케이스만(x_c [disp_m]). 스윕 없음.
#   러너가 자기 플래그와 그 값을 앞에 붙이고 셸의 '--' 는 전달되지 않는다.
#   => 뒤에서부터 훑어 러너 토큰을 건너뛰고 숫자가 끊길 때까지 모은다(§1).
_LAUNCHER_WORDS = ('abaqus', 'cae', 'cae.exe', 'abq', 'standard', 'explicit')


def _is_launcher_token(tok):
    """러너/구분자/스크립트 토큰이면 True (사용자 인자가 아니다)."""
    t = tok.strip()
    if t in ('', '--'):
        return True
    if t.startswith('-'):                     # -cae, -noGUI, -tmpdir ...
        return True
    if t.lower() in _LAUNCHER_WORDS:
        return True
    if t.lower().endswith('.py') or t.lower().startswith('nogui='):
        return True
    if len(t) > 1 and t[1] == ':':            # C:\... 같은 Windows 경로
        return True
    return False


def parse_args(argv):
    """맨 뒤의 숫자 1~2개를 <x_c> [disp_m] 로 읽는다.

    조용한 치환 금지: 개수가 맞지 않거나 숫자가 아닌 토큰을 만나면 추측하지 않고
    사용법과 argv 전문을 찍고 예외로 끝낸다.
    """
    toks = list(argv)
    print("%s argv=%s" % (TAG, toks))
    if not any(t.lower().endswith('.py') for t in toks):
        print("%s WARNING: no script-name token in argv (unexpected launcher layout)." % TAG)

    # 레거시 `-- HF x_c d_c` 형식은 두 번째 숫자의 의미가 다르다(d_c 비율 vs disp_m [m]).
    # 조용히 재해석하면 1.0 m 같은 엉뚱한 변위가 되므로 거부한다.
    legacy = [t for t in toks if t.strip().upper() in ('LF', 'HF')]
    if legacy:
        raise RuntimeError(
            'Legacy fidelity form %r is not accepted: the second number now means '
            'disp_m [m], not the d_c ratio.\n'
            'Usage: abaqus cae noGUI=run_abaqus_buckle.py -- <x_c> [disp_m] '
            '(e.g. -- 0.5 1e-3)' % (legacy[0],))

    nums, stop = [], None
    for t in reversed(toks):
        if _is_launcher_token(t):
            continue
        try:
            nums.append(float(t))
        except ValueError:
            stop = t
            break
    nums.reverse()
    if not (1 <= len(nums) <= 3):
        raise RuntimeError(
            'Expected 1 to 3 trailing numeric arguments (<x_c> [disp_m] [clamp_pull_m]), got %r '
            '(first non-numeric token from the end: %r).\n'
            'Usage: abaqus cae noGUI=run_abaqus_buckle.py -- <x_c> [disp_m] [clamp_pull_m]'
            % (nums, stop))
    return (nums[0],
            (nums[1] if len(nums) > 1 else None),
            (nums[2] if len(nums) > 2 else None))


x_c, _disp_arg, _clamp_arg = parse_args(sys.argv)
DISP = DISP_GLOBAL if _disp_arg is None else _disp_arg
# 클램프 법선 당김 — 기본은 코너 당김과 같은 비율(A-route 의 d_c), CLI 3번째 인자로 직접 지정 가능.
#   예: -- 0.5 0 5e-5  -> 코너 미구동 + 클램프만 5e-5 m (클램프 단독 구동 진단)
CLAMP_PULL = (DISP * CLAMP_DC if _clamp_arg is None else _clamp_arg)
CLAMP_PERT = PERTURBATION * CLAMP_DC    # 좌굴 스텝 클램프 섭동 [m] (A-route: PERTURBATION*d_c)

if not (0.03 <= x_c <= 0.95):
    # 원장 M-8: x_c < 0.028 이면 클램프 패치가 정점 패치와 겹치고, x_c ~ 1 이면
    # 모서리 패치와 겹친다 -> 조용히 다른 모델이 된다. 경고만 찍고 진행한다.
    print("%s WARNING: x_c=%g is outside the safe band (0.03, 0.95) — "
          "the clamp patch may overlap a vertex patch (ledger M-8)." % (TAG, x_c))
if DISP < 0.0:
    raise RuntimeError('disp_m must be >= 0 (got %r).' % (DISP,))
if DISP == 0.0:
    # 코너를 구동하지 않는다 = 클램프 단독 구동 케이스 (진단용).
    # Initial 단계의 Disp_Control_Right/Left(u1=u2=0) 값이 그대로 유지된다.
    print("%s DISP=0 : 코너 미구동 — 클램프 단독 구동 케이스로 진행한다." % (TAG,))
if CLAMP_MODE == 'driven' and CLAMP_PULL <= 0.0:
    raise RuntimeError('CLAMP_MODE=driven 인데 CLAMP_PULL=%r <= 0 이다.' % (CLAMP_PULL,))
print("%s CLAMP_MODE=%s d_c=%.3g -> CLAMP_PULL=%.4e m (buckle pert %.3e m)"
      % (TAG, CLAMP_MODE, CLAMP_DC, CLAMP_PULL, CLAMP_PERT))

V_CL = clamp_coord_L(x_c)
V_CR = clamp_coord_R(x_c)

print("%s x_c=%g -> V_CL=%s V_CR=%s" % (TAG, x_c, V_CL, V_CR))
print("%s corner pull DISP=%.4e m (= alpha %.4g x 5e-5 m)  perturbation=%.3e m"
      % (TAG, DISP, DISP / 5.0e-5, PERTURBATION))
print("%s buckle step: solver=%s numEigen=%d vectors=%d"
      % (TAG, BUCKLE_SOLVER, N_EIG_BUCKLE, BUCKLE_VECTORS))

# 하중 각도 (28.6도) — 케이블 방향과 동일하게 유지 (하중 경로 동일화)
angle_deg = 28.6
angle_rad = np.deg2rad(angle_deg)
cos_val = float(np.cos(angle_rad))
sin_val = float(np.sin(angle_rad))

_JOB_ARTIFACTS = ('odb', 'fil', 'sta', 'msg', 'lck', 'com', 'prt', 'sim', 'log',
                  'dat', 'res', 'abq', 'ipm', 'mdl', 'stt', 'cid')


def run_job_safely(job_name, model_name=None):
    if model_name is None:
        raise RuntimeError('run_job_safely: model_name is required.')
    for _ext in _JOB_ARTIFACTS:
        _f = '%s.%s' % (job_name, _ext)
        if os.path.exists(_f):
            print("Removing stale artifact: %s" % _f)
            try:
                os.remove(_f)
            except OSError as _e:
                print("  (warning) could not remove %s: %s" % (_f, _e))

    lck_file = job_name + '.lck'
    odb_file = job_name + '.odb'
    if os.path.exists(lck_file):
        print("Detected old lock file: %s. Removing it..." % lck_file)
        try:
            os.remove(lck_file)
        except OSError:
            raise RuntimeError('락 파일을 지울 수 없습니다: %s '
                               '(다른 Abaqus 프로세스가 돌고 있습니까?)' % lck_file)

    if job_name in mdb.jobs:
        del mdb.jobs[job_name]

    _ncp = NUMCPUS
    print("%s numCpus=%d numDomains=%d (코드 상수 NUMCPUS)" % (TAG, _ncp, _ncp))
    job = mdb.Job(name=job_name, model=model_name, numCpus=_ncp, numDomains=_ncp)
    print("Submitting Job: %s" % job_name)
    job.writeInput(consistencyChecking=OFF)
    job.submit(consistencyChecking=OFF)

    job.waitForCompletion()

    time.sleep(1.0)

    # 진단 출력은 호출측에서 report_job() 한 번으로 끝낸다.

# 잡 성공 판정 = 산출물 + 완주 문자열. 'ODB 존재'도 'Exit code 0'도 증거가 아니다(§2).
    if (job.status == ABORTED or not os.path.exists(odb_file)
            or not job_completed_ok(job_name)):
        raise RuntimeError('Job %s 실패 (Status=%s). sys.exit 대신 예외로 올린다: '
                           'CAE noGUI 러너에서 sys.exit 은 종료코드 0으로 보인다.'
                           % (job_name, str(job.status)))

    if not job_completed_ok(job_name):
        print("!!! WARNING: %s — .sta/.msg 에 완주 문자열 없음 (%s) "
              "(중도 중단 의심; odb 존재만으로는 판정 불가 — R-12)"
              % (job_name, ' / '.join(_COMPLETION_STRINGS)))
    print("Job %s completed successfully (Status: %s)." % (job_name, str(job.status)))
    return True


_COMPLETION_STRINGS = ('HAS COMPLETED SUCCESSFULLY', 'THE ANALYSIS HAS BEEN COMPLETED')


def job_completed_ok(job_name):
    """R-12: .sta/.msg 의 완주 문자열로 완주를 판정한다.

    odb 파일 존재만 보면 '중도 중단된 odb'를 성공으로 오판한다
    (2026-09-21 사례: HF_Postbuckle.odb 는 존재하나 Step-Postbuckle 프레임 0개).

    완주 문자열은 두 가지다. 좌굴 스텝은 증분 블록이 없어 .sta 에
    'THE ANALYSIS HAS COMPLETED SUCCESSFULLY' 를 남기지 않고 .msg 에
    'THE ANALYSIS HAS BEEN COMPLETED' 만 남긴다 -> 한 문자열만 보면 성공한
    좌굴 잡을 '중도 중단 의심'으로 잘못 보고한다 (실측 2026-09-22: 좌굴형
    .msg + .sta 조합에서 False, 이 스크립트의 α 요약표 '완주판정' 열이 전부 False).
    중단된 잡은 'HAS NOT BEEN COMPLETED' 이므로 두 문자열 어느 것과도 일치하지 않는다.
    """
    for ext in ('sta', 'msg'):
        fn = '%s.%s' % (job_name, ext)
        if os.path.exists(fn):
            with open(fn, 'r', errors='replace') as f:
                _up = f.read().upper()
            if any(_s in _up for _s in _COMPLETION_STRINGS):
                return True
    return False


def report_job(job_name):
    """결과 판정에 필요한 줄만 콘솔에 찍고 <job>.diag.txt 로도 남긴다.

    좌굴 잡의 결론은 두 곳에만 있다.
      .dat -> MODE NO / EIGENVALUE 표 (고유값 원문)
      .msg -> CONVERGED / REQUESTED BY THE USER / CANNOT BE FOUND / REDUCED TO
    전체 로그를 다 찍으면 정작 필요한 줄이 묻히므로 그 줄만 뽑는다.
    """
    out = ['===== job_completed_ok = %s =====' % job_completed_ok(job_name)]

    msg_keys = ('CONVERGED', 'REQUESTED BY THE USER', 'CANNOT BE FOUND',
                'REDUCED TO', 'NEGATIVE EIGENVALUES', 'HAS BEEN COMPLETED',
                'HAS NOT BEEN COMPLETED', '***ERROR')
    # .log 는 '잡이 시작조차 못 한' 실패(라이선스 거부/입력 거부)에만 단서가 있다.
    _logfn = '%s.log' % job_name
    if os.path.exists(_logfn):
        out.append('===== %s (tail 20) =====' % _logfn)
        try:
            with open(_logfn, 'r', errors='replace') as _lf:
                for _l2 in _lf.read().splitlines()[-20:]:
                    out.append('      | ' + _l2[:150])
        except Exception as _le:
            out.append('      | (읽기 실패: %s)' % _le)

    fn = '%s.msg' % job_name
    if os.path.exists(fn):
        with open(fn, 'r', errors='replace') as f:
            lines = f.read().splitlines()
        hits = [ln.rstrip() for ln in lines
                if any(k in ln.upper() for k in msg_keys)]
        out.append('===== %s : %d lines, %d key hits (last 12) ====='
                   % (fn, len(lines), len(hits)))
        out.extend(hits[-12:])

    fn = '%s.dat' % job_name
    if os.path.exists(fn):
        with open(fn, 'r', errors='replace') as f:
            lines = f.read().splitlines()
        idx = [i for i, ln in enumerate(lines)
               if 'MODE NO' in ln.upper() or 'BUCKLING FACTOR' in ln.upper()]
        if idx:
            i0, i1 = max(0, idx[0] - 3), min(len(lines), idx[-1] + 14)
            out.append('===== %s : eigen table (lines %d-%d) ====='
                       % (fn, i0 + 1, i1))
            out.extend(lines[i0:i1])
        else:
            out.append('===== %s : NO eigen table -> the step produced no modes =====' % fn)

    if len(out) == 1:
        out.append('(no .msg/.dat found — the job died before writing them)')

    # 파일을 **먼저** 쓴다. 콘솔 출력은 예외로 프로세스가 죽으면 버퍼가 flush 되지 않고
    #   사라진다 — 그러면 사용자에게 남는 것이 없다(Status=None 사례에서 실제로 겪었다).
    dst = os.path.join(_RUN, '%s.diag.txt' % job_name)
    import io as _io
    with _io.open(dst, 'w', encoding='utf-8', errors='replace') as f:
        f.write('\n'.join(out) + '\n')
    for ln in out:
        print('      | %s' % ln[:170], flush=True)
    print('%s diag saved: %s (%d lines)' % (TAG, dst, len(out)), flush=True)
    return dst
    return dst
def create_rigid_patch(a, inst_memb, name, coord, radius):
    """
    지정된 좌표 기준 radius 내의 노드들을 묶어 강체운동을 하도록 Tie 설정
    실제 solar sail 에서 케이블이나 클램프를 설치하기 위해 sail 면에 테이프 등을 설치하는 과정을 모사
    (run_abaqus_new.py 와 동일한 정의 — 좌표/반지름/결합 방식 모두 그대로)
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

    # Coupling (RP <-> Membrane Nodes)
    my_model.Coupling(
        name=name+'_Coupling', controlPoint=rp_region, surface=patch_set,
        influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING, weightingMethod=UNIFORM,
        u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON
    )
    return rp, rp_region



def connect_cable(a, name, part, sail_corner, vector_dir, radius=1e-4):
    """케이블 인스턴스를 꼭짓점에 붙이고 (시작단, 끝단) 노드 영역을 돌려준다.

    mode.py 와 같은 정의. 케이블은 'RP <-> 구동 끝단' 사이의 하중 전달 로드이므로,
    이것을 빼고 RP 를 직접 구동하면 하중 경로가 같다고 볼 수 없다 — 그래서 복원한다.
    """
    inst_name = 'Inst_' + name
    inst = a.Instance(name=inst_name, part=part, dependent=ON)
    target_vec = np.array(vector_dir)
    target_vec = target_vec / np.linalg.norm(target_vec)
    rot_angle_deg = np.degrees(np.arctan2(target_vec[1], target_vec[0]))
    a.rotate(instanceList=(inst_name,), axisPoint=(0, 0, 0), axisDirection=(0, 0, 1),
             angle=rot_angle_deg)
    a.translate(instanceList=(inst_name,), vector=sail_corner)
    node_start = inst.nodes.getByBoundingSphere(center=sail_corner, radius=radius)
    if len(node_start) == 0:
        print("%s Warning: Node not found by sphere, trying closest for %s" % (TAG, name))
        node_start = inst.nodes.getClosest(coordinates=sail_corner)
    region_start = regionToolset.Region(nodes=node_start)
    cable_len = LEN_TOP if 'Top' in name else LEN_BOT
    end_coord = (sail_corner[0] + target_vec[0] * cable_len,
                 sail_corner[1] + target_vec[1] * cable_len, 0.0)
    node_end = inst.nodes.getByBoundingSphere(center=end_coord, radius=radius)
    region_end = regionToolset.Region(nodes=node_end)
    return region_start, region_end



def build_model(disp):
    """run_abaqus_new.py 의 모델 생성을 '케이블만 제외' 하고 재현한다."""
    global my_model
    model_name = MODEL_PREFIX
    if model_name in mdb.models:
        del mdb.models[model_name]
    my_model = mdb.Model(name=model_name)

    # 멤브레인 (Kapton) — 동일
    mat = my_model.Material(name='Kapton')
    mat.Density(table=((1420.0, ),))
    mat.Elastic(table=((2.5e9, 0.34),))
    my_model.HomogeneousShellSection(name='Section-Membrane', material='Kapton', thickness=THICKNESS)   # [2026-10-04] 막 실험 철회로 셸 섹션 복귀

    # 케이블 재질 (Kevlar) — mode.py 와 동일. 클램프와 무관하게 케이블 경로를 복원한다.
    mat_cable = my_model.Material(name='Kevlar')
    mat_cable.Elastic(table=((62.0e9, 0.36),))
    my_model.TrussSection(name='Section-Cable', material='Kevlar', area=CABLE_AREA)

    # 파트 생성: 멤브레인 — 동일
    s = my_model.ConstrainedSketch(name='triangle_profile', sheetSize=BASE*2)
    s.Line(point1=V3[:2], point2=V2[:2])
    s.Line(point1=V2[:2], point2=V1[:2])
    s.Line(point1=V1[:2], point2=V3[:2])
    p = my_model.Part(name='Membrane', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseShell(sketch=s)
    p.SectionAssignment(region=p.Set(faces=p.faces, name='All'), sectionName='Section-Membrane')

    def create_cable_part(name, length):
        """mode.py 와 동일한 케이블 1요소 파트 (T3D2)."""
        p_c = my_model.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        p_c.WirePolyLine(points=((0.0, 0.0, 0.0), (length, 0.0, 0.0)), mergeType=IMPRINT, meshable=ON)
        p_c.SectionAssignment(region=p_c.Set(edges=p_c.edges, name='Wire'), sectionName='Section-Cable')
        p_c.seedPart(size=length)  # 요소 1개
        elemTypeTruss = ElemType(elemCode=T3D2, elemLibrary=STANDARD)
        p_c.setElementType(regions=(p_c.edges,), elemTypes=(elemTypeTruss,))
        p_c.generateMesh()
        return p_c

    p_cable_top = create_cable_part('Cable_Top', LEN_TOP)
    p_cable_bot = create_cable_part('Cable_Bot', LEN_BOT)

    a = my_model.rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    inst_memb = a.Instance(name=INSTANCE_NAME, part=p, dependent=ON)

    p.seedPart(size=BASE/SEED_DIV, deviationFactor=0.1)  # 약 1.82만개 (1차 요소)
    p.setMeshControls(regions=p.faces, elemShape=QUAD_DOMINATED, technique=FREE, algorithm=MEDIAL_AXIS)
    elemTypeQuad = ElemType(elemCode=ELEM_CODE_QUAD, elemLibrary=STANDARD)
    elemTypeTri = ElemType(elemCode=ELEM_CODE_TRI, elemLibrary=STANDARD)
    p.setElementType(regions=(p.faces,), elemTypes=(elemTypeQuad, elemTypeTri))
    p.generateMesh()
    a.regenerate()

    # 꼭짓점 RP — 동일 (radius=0.2)
    rp1_obj, rp1_reg = create_rigid_patch(a, inst_memb, 'Top', V1, radius=0.2)
    rp2_obj, rp2_reg = create_rigid_patch(a, inst_memb, 'Right', V2, radius=0.2)
    rp3_obj, rp3_reg = create_rigid_patch(a, inst_memb, 'Left', V3, radius=0.2)

    # 클램프 RP : 우측 빗변 중점 (15, 5), 좌측 빗변 중점 (5, 5) — 동일
    #   CLAMP_MODE='none' 이면 패치를 아예 만들지 않는다 (논문 재현 기준선 P0).
    if CLAMP_MODE == 'none':
        rp_cl_obj = rp_cr_obj = None
    else:
        rp_cl_obj, rp_cl_reg = create_rigid_patch(a, inst_memb, 'CL', V_CL, radius=0.2)
        rp_cr_obj, rp_cr_reg = create_rigid_patch(a, inst_memb, 'CR', V_CR, radius=0.2)

    # 케이블 연결: 정점은 위로, 두 아래 모서리는 각 케이블 축(28.6도) 방향 — mode.py 와 동일
    start_c1, end_c1 = connect_cable(a, 'Cable_Top', p_cable_top, V1, (0.0, 1.0, 0.0))
    start_c2, end_c2 = connect_cable(a, 'Cable_Right', p_cable_bot, V2, (cos_val, -sin_val, 0.0))
    start_c3, end_c3 = connect_cable(a, 'Cable_Left', p_cable_bot, V3, (-cos_val, -sin_val, 0.0))

    a.Set(name='RP_Top_Set', referencePoints=(a.referencePoints[rp1_obj.id],))
    a.Set(name='RP_Right_Set', referencePoints=(a.referencePoints[rp2_obj.id],))
    a.Set(name='RP_Left_Set', referencePoints=(a.referencePoints[rp3_obj.id],))
    if CLAMP_MODE != 'none':
        a.Set(name='RP_CL_Set', referencePoints=(a.referencePoints[rp_cl_obj.id],))
        a.Set(name='RP_CR_Set', referencePoints=(a.referencePoints[rp_cr_obj.id],))
    # Tie (RP <-> 케이블 시작단) — mode.py 와 동일
    my_model.Tie(name='Tie_Top', main=a.sets['RP_Top_Set'], secondary=start_c1,
                 positionToleranceMethod=COMPUTED)
    my_model.Tie(name='Tie_Right', main=a.sets['RP_Right_Set'], secondary=start_c2,
                 positionToleranceMethod=COMPUTED)
    my_model.Tie(name='Tie_Left', main=a.sets['RP_Left_Set'], secondary=start_c3,
                 positionToleranceMethod=COMPUTED)
    a.regenerate()

    # ---- Step 1: GlobalTension (프리텐션) — run_abaqus_new.py 와 동일 ----
    my_model.StaticStep(
        name='Step-GlobalTension',
        previous='Initial',
        nlgeom=ON,
        stabilizationMagnitude=MODE_STABILIZATION,   # mode.py 와 동일 (0.0005)
        stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
        continueDampingFactors=False,
        adaptiveDampingRatio=0.05,
        initialInc=0.0001, minInc=1e-8, maxNumInc=1000
    )

    # ---- Step 2: Buckle — Trigger 스텝을 대체한다 ----
    # 대조 스크립트(run_abaqus_cable.py)에서 완주한 설정을 그대로 사용:
    #   numEigen=100 (음수 모드 건너뛰기), SUBSPACE, vectors=250, maxIterations=5000
    if 'Step-Buckle' in my_model.steps:
        del my_model.steps['Step-Buckle']
    # EVOL 필수 — 없으면 base_state_probe 가 면적가중 대신 균등가중으로 떨어진다.
    my_model.FieldOutputRequest(name='F-Output-1',
                                createStepName='Step-GlobalTension',
                                variables=('S', 'E', 'U', 'COORD', 'EVOL', 'RF'))
    # RF = 앵커 반력. 각 앵커가 당김을 얼마나 흡수하는가(하중 경로 분담)를 읽는다.
    #   출력 요청은 해석 결과를 바꾸지 않는다.

    # ---- 좌굴 스텝: 솔버는 코드 상수 하나로 교체 (SUBSPACE <-> LANCZOS) ----
    _eig = dict(name='Step-Buckle', previous='Step-GlobalTension',
                numEigen=N_EIG_BUCKLE, eigensolver=EIGENSOLVER_CONST)
    if BUCKLE_SOLVER == 'SUBSPACE':
        _eig.update(vectors=BUCKLE_VECTORS, maxIterations=BUCKLE_MAXITER)
    else:
        _eig.update(blockSize=BUCKLE_BLOCK_SIZE, maxIterations=BUCKLE_MAXITER,
                    minEigen=BUCKLE_MIN_EIGEN)
        if BUCKLE_MAX_EIGEN is not None:
            _eig.update(maxEigen=BUCKLE_MAX_EIGEN)
    print("%s eigensolver=%s numEigen=%d : %s"
          % (TAG, BUCKLE_SOLVER, N_EIG_BUCKLE,
             ', '.join('%s=%s' % (_k, _eig[_k]) for _k in sorted(_eig)
                       if _k not in ('name', 'previous'))))
    my_model.BuckleStep(**_eig)

    # ---- 초기응력 (수렴 보조) — 동일 ----
    my_model.Stress(
        name='Initial_Stiffness',
        region=inst_memb.sets['All'],
        distributionType=UNIFORM,
        sigma11=SIGMA0, sigma22=SIGMA0, sigma33=0.0,
        sigma12=0.0, sigma13=0.0, sigma23=0.0
    )

    # ---- 구속조건 ----
    # 시작 시 전체 면 z 고정 (평탄) — 동일
    my_model.DisplacementBC(
        name='BC_Stabilize_Z',
        createStepName='Step-GlobalTension',   # [mode 정합] mode.py 와 동일 (buckle 은 Initial 이었다)
        region=inst_memb.sets['All'],
        u3=SET
    )
    # Top 정점: 완전 고정 — mode.py 와 동일하게 케이블 끝단(end_c1)에 건다.
    #   RP_Top_Set 은 Tie 로만 연결된다(mode.py 도 그렇다).
    my_model.DisplacementBC(
        name='BC_Anchor_Top',
        createStepName='Initial',
        region=end_c1,
        u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET
    )

    # Right/Left 정점: mode.py 와 동일하게 u3+회전만 고정하고, 구동은 케이블 끝단에 건다.
    my_model.DisplacementBC(
        name='BC_Right_Z', createStepName='Initial', region=end_c2,
        u3=SET, ur1=SET, ur2=SET, ur3=SET
    )
    my_model.DisplacementBC(
        name='BC_Left_Z', createStepName='Initial', region=end_c3,
        u3=SET, ur1=SET, ur2=SET, ur3=SET
    )
    disp_a = disp
    my_model.DisplacementBC(name='Disp_Control_Right', createStepName='Initial',
                            region=end_c2, u1=SET, u2=SET)
    my_model.DisplacementBC(name='Disp_Control_Left', createStepName='Initial',
                            region=end_c3, u1=SET, u2=SET)
    my_model.boundaryConditions['Disp_Control_Right'].setValuesInStep(
        stepName='Step-GlobalTension',
        u1=disp_a * cos_val,
        u2=-disp_a * sin_val
    )
    my_model.boundaryConditions['Disp_Control_Left'].setValuesInStep(
        stepName='Step-GlobalTension',
        u1=-disp_a * cos_val,
        u2=-disp_a * sin_val
    )
    # [paper3] 위 꼭짓점도 케이블 축(+y)으로 당긴다 — mode.py 와 동일.
    #   corner2(기본)에서는 BC_Anchor_Top 의 완전고정이 그대로 유지된다.
    if PRETENSION_MODE == 'paper3':
        my_model.boundaryConditions['BC_Anchor_Top'].setValuesInStep(
            stepName='Step-GlobalTension',
            u1=0.0, u2=disp_a * DISP_TOP_OVER_CORNER
        )
    print("%s 프리텐션 = %s : 아래 %.4e m%s"
          % (TAG, PRETENSION_MODE, disp_a,
             (' / 위 %.4e m' % (disp_a * DISP_TOP_OVER_CORNER)) if PRETENSION_MODE == 'paper3'

             else ' (Top 완전고정)'))
    # ---- CLAMP_MODE 별 RP 처리 (4종 의미는 §11) ----
    _sq2 = 2.0 ** 0.5
    if CLAMP_MODE == 'none':
        pass
    elif CLAMP_MODE == 'fixed':
        my_model.DisplacementBC(name='BC_Clamp_CL', createStepName='Initial',
                                region=a.sets['RP_CL_Set'], u1=0, u2=0, u3=0)
        my_model.DisplacementBC(name='BC_Clamp_CR', createStepName='Initial',
                                region=a.sets['RP_CR_Set'], u1=0, u2=0, u3=0)
    else:
        my_model.DisplacementBC(name='BC_Clamp_CL', createStepName='Initial',
                                region=a.sets['RP_CL_Set'], u3=0)
        my_model.DisplacementBC(name='BC_Clamp_CR', createStepName='Initial',
                                region=a.sets['RP_CR_Set'], u3=0)
        if CLAMP_MODE == 'driven':
            # A-route 와 동일 패턴: u3 전용 BC 와 in-plane 구동 BC 를 분리해 만든다.
            my_model.DisplacementBC(name='Disp_Clamp_CL', createStepName='Initial',
                                    region=a.sets['RP_CL_Set'], u1=0, u2=0)
            my_model.DisplacementBC(name='Disp_Clamp_CR', createStepName='Initial',
                                    region=a.sets['RP_CR_Set'], u1=0, u2=0)
            my_model.boundaryConditions['Disp_Clamp_CL'].setValuesInStep(
                stepName='Step-GlobalTension',
                u1=-CLAMP_PULL/_sq2, u2=+CLAMP_PULL/_sq2)
            my_model.boundaryConditions['Disp_Clamp_CR'].setValuesInStep(
                stepName='Step-GlobalTension',
                u1=+CLAMP_PULL/_sq2, u2=+CLAMP_PULL/_sq2)

    # ---- Buckle 스텝: 논문 (c) 3변 u3=0 유지 + 클램프 구간은 면외 구속에서 제외(§9) ----
    _excl_edges = (aba_grid_mesh.clamp_exclude_labels(inst_memb, (V_CL, V_CR), CLAMP_EXCL_R)
                   if CLAMP_MODE != 'none' else set())
    _keep_labels = aba_grid_mesh.boundary_node_labels(
        inst_memb, base=BASE, height=HEIGHT, exclude_labels=_excl_edges)
    print("%s All_Edges_NoClamp 후보: 경계 노드 %d개 (클램프 반경 %.3g m 내 %d개 제외, 전체 막 노드 %d개)"
          % (TAG, len(_keep_labels), CLAMP_EXCL_R, len(_excl_edges), len(inst_memb.nodes)))
    aba_grid_mesh.make_set(a, 'All_Edges_NoClamp', inst_memb, _keep_labels)
    print("%s All_Edges_NoClamp 셋 생성 완료 (노드 %d개)" % (TAG, len(_keep_labels)))

    my_model.boundaryConditions['BC_Stabilize_Z'].deactivate('Step-Buckle')
    my_model.DisplacementBC(
        name='BC_Edges_Only_Z',
        createStepName='Step-Buckle',
        region=a.sets['All_Edges_NoClamp'],
        u3=0
    )

    # ---- Buckle 스텝의 perturbation 하중 (검증된 대조 스크립트와 동일한 크기) ----
    my_model.boundaryConditions['Disp_Control_Right'].setValuesInStep(
        stepName='Step-Buckle',
        u1=PERTURBATION * cos_val,
        u2=-PERTURBATION * sin_val
    )
    my_model.boundaryConditions['Disp_Control_Left'].setValuesInStep(
        stepName='Step-Buckle',
        u1=-PERTURBATION * cos_val,
        u2=-PERTURBATION * sin_val
    )

    # 클램프 섭동 (driven 모드만): A-route 의 CLAMP_PERT = PERTURBATION*d_c 와 동일.
    #   클램프도 구동점이면 좌굴 스텝에서 같은 방식으로 섭동을 줘야 K_delta 가 일관된다.
    if CLAMP_MODE == 'driven':
        my_model.boundaryConditions['Disp_Clamp_CL'].setValuesInStep(
            stepName='Step-Buckle',
            u1=-CLAMP_PERT/_sq2, u2=+CLAMP_PERT/_sq2)
        my_model.boundaryConditions['Disp_Clamp_CR'].setValuesInStep(
            stepName='Step-Buckle',
            u1=+CLAMP_PERT/_sq2, u2=+CLAMP_PERT/_sq2)

    return model_name

# job 이름에 요소 태그 + 클램프 모드를 넣는다. 같은 이름이면 run_job_safely 가 stale
#   산출물을 지워 직전 증거가 사라진다(§11).
ELEM_TAG = 's4'
JOB_NAME = 'Buckle_xc%03d_d%03dum_%s_%s' % (int(round(x_c * 100.0)),
                                            int(round(DISP * 1.0e6)),
                                            ELEM_TAG, CLAMP_MODE)

print("")
print("=" * 78)
print("%s single case: x_c=%g  DISP=%.4e m  model=%s  job=%s"
      % (TAG, x_c, DISP, MODEL_PREFIX, JOB_NAME))
print("=" * 78)

try:
    _model_name = build_model(DISP)
    run_job_safely(JOB_NAME, _model_name)
except Exception:
    print("%s BUILD/RUN FAILED — traceback:" % TAG)
    print(traceback.format_exc())
    # 실패해도 진단을 찍는다. 원인은 .msg/.dat/.log 에만 있고, 여기서 안 찍으면
    #   사용자에게는 traceback 한 줄만 남아 아무것도 판정할 수 없다.
    for _fn in (JOB_NAME + ".log", JOB_NAME + ".msg", JOB_NAME + ".dat"):
        _p = os.path.join(os.getcwd(), _fn)
        if os.path.exists(_p):
            print("%s [실패진단] %s (%d bytes)" % (TAG, _fn, os.path.getsize(_p)), flush=True)
        else:
            print("%s [실패진단] %s 없음 — 잡이 시작조차 못 했다" % (TAG, _fn), flush=True)
        try:
            sys.stdout.flush()
        except Exception:
            pass
    _diag_dst = os.path.join(_RUN, "%s.diag.txt" % JOB_NAME)
    print("%s [실패진단] report_job 호출 -> %s" % (TAG, _diag_dst), flush=True)
    try:
        report_job(JOB_NAME)
    except Exception as _re:
        print("%s report_job 실패(무시): %s" % (TAG, _re), flush=True)
    raise

report_job(JOB_NAME)

# base state 의 압축 정도를 lambda 판정과 함께 보기 위한 선택 단계
try:
    _probe = os.path.join(_HERE, 'base_state_probe.py')
    _odb = '%s.odb' % JOB_NAME
    if os.path.exists(_probe) and os.path.exists(_odb):
        print("%s base-state probe: %s / Step-GlobalTension" % (TAG, _odb))
        subprocess.call('abaqus python "%s" "%s" Step-GlobalTension'
                        % (_probe, _odb), shell=True)
    else:
        print("%s base-state probe skipped (probe=%s, odb=%s)"
              % (TAG, os.path.exists(_probe), os.path.exists(_odb)))
except Exception as _e:
    print("%s base-state probe failed (ignored): %s" % (TAG, _e))

_rel_fil = '..' + os.sep + RUN_DIR_NAME + os.sep + JOB_NAME
print("")
print("=" * 78)
print("%s PASS CRITERION: CONVERGED > 0 and the first eigenvalues positive." % TAG)
print("%s eigen table: %s.dat" % (TAG, JOB_NAME))
print("%s HF imperfection keyword for this run:" % TAG)
print("     *IMPERFECTION, FILE=%s, STEP=2" % _rel_fil)
print("=" * 78)
