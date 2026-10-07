"""
run_abaqus_buckle.py — 좌굴(선형 고유값) 해석 전용. 1회 실행, 스윕 없음.

목적
    클램프 패치가 있는 삼각 막 모델의 선형 좌굴 고유모드를 한 번에 추출한다.
    추출한 모드는 HF 포스트버클링의 초기결함(*IMPERFECTION, FILE=..., STEP=2)으로 쓴다.

모델 (run_abaqus_new.py 의 HF 모델과 같은 빌드 블록을 공유한다)
    삼각 막(BASE=20 m, HEIGHT=10 m, 두께 5e-6 m)
    + 꼭짓점 강체패치 3개 + 클램프 강체패치 2개 (위치 x_c).
    스텝 (2개 — 클램프를 별도 텐션 스텝으로 분리하지 않는다):
        Step-GlobalTension      꼭짓점 + 클램프 하중 (DEAD = base state)
        Step-Buckle             고유값 추출 (LIVE = perturbation 0.01 m)
    창은 N_EIG_BUCKLE / BUCKLE_VECTORS 상수로 정한다(요청 수 <= 실제 subspace 차원).

케이스 (CASE 상수 — 2026-10-07). 예전의 LOAD_MODE x CLAMP_MODE 조합 표면을 대체한다.
    CASE          DEAD (base state)                 LIVE (lambda 의 대상)      답하는 질문
    seed          클램프 없음, 코너 당김            코너 당김                  코너 하중의 좌굴 배수
    control_none  클램프 패치 u3=0, 당김 0          코너 당김                  (클램프 구속이 무해한가)
    clamp_lf      DEAD_FRAC x P + 클램프 u3=0       (1-DEAD_FRAC) x P          운용 램프의 좌굴 배수
    paper_s1      논문 baseline(케이블/클램프 없음) 논문 0.10 m                우리 추출 셋업이 맞는가
      P = 꼭짓점 + 클램프 하중 세트(CORNER_F0 / CLAMP_F0). clamp_lf 는 케이블 0개다.
      BUCKLING LOAD = DEAD + lambda*LIVE 이므로 **clamp_lf 의 lambda = 1 이 운용점**이다.
      무효 조합(예: 코너는 케이블, 클램프는 CLOAD)은 상수 검증이 런타임 예외로 막는다.

사용법
    abaqus cae noGUI=run_abaqus_buckle.py -- <x_c> [disp_m] [clamp_pull_m] [clamp_excl_r_m]
        x_c            클램프 위치 파라미터. 0.5 -> 좌(5,5) / 우(15,5)
        disp_m         GlobalTension 코너 당김 [m]. 생략하면 DISP_GLOBAL (1.8e-5).
        clamp_pull_m   클램프 당김 [m]. 생략하면 disp_m * CLAMP_DC.
        clamp_excl_r_m 클램프 반경 내 노드를 면외 구속에서 제외 [m]. 기본 0.0.
    예)  abaqus cae noGUI=run_abaqus_buckle.py -- 0.5
         abaqus cae noGUI=run_abaqus_buckle.py -- 0.25 1.8e-5 9.0e-6
    CASE 는 CLI 가 아니라 **상수부에서 직접 수정한다**. 값이 잡 이름에 들어가므로 산출물이
    서로 덮이지 않는다. 케이스가 요구하는 상수(PRETENSION_MODE / CORNER_ANGLE_DEG)가 어긋나면
    실행 전에 예외로 막힌다.

    추출 창(N_EIG_BUCKLE / BUCKLE_VECTORS)은 **상수부에서 직접 수정한다** (CLI 로 받지 않는다).
    ⚠️ N_EIG_BUCKLE 은 '실제 subspace 차원' 이하여야 한다. 실측(2026-10-07): vectors=500 을
       요청해도 'VECTORS IN SUBSPACE IS REDUCED TO 14' 가 찍히면 200 개는 원리적으로 못 찾는다.
       창을 키우는 대신 요청 수를 줄이는 쪽이 맞다.
    disp_m 에는 1 mm 상한 가드가 있다(인자 순서가 어긋나 조용히 다른 모델이 만들어지는 것을 막는다).

산출물 (code/buckle/)
    <job>.odb / .dat / .msg / .sta / .fil / .diag.txt
    job 이름은 인자+요소+클램프모드에서 자동 생성한다: Buckle_xc<NNN>_d<NNN>um_<elem>_<clamp>
      (예: Buckle_xc050_d0050um_s4r_driven) — 케이스가 바뀌어도 산출물이 서로 덮이지 않는다.
    CASE / CORNER_F0 / CLAMP_F0 / DEAD_FRAC 는 상단 상수.
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
#   [2026-10-07] CLI 4번째 인자로 덮어쓸 수 있다(d_c 스윕과 함께 클램프 구간 면외 자유를 본다).
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

# 코너 하중 각도 [deg] — CASE 가 요구한다('hf' = 28.6 = 우리 HF/케이블 방향, 'paper45' = 논문).
CORNER_ANGLE_DEG = 28.6

# ============================================================================
# 케이스 (CASE) — '무엇을 DEAD, 무엇을 LIVE 로 두는가' 만 고른다.   (§18)
#   왜: 예전에는 LOAD_MODE(2) x CLAMP_MODE(5) = 10 조합이 독립이라 **무효 조합이 조용히
#       만들어졌다**. 실제로 실패 조합(cable + cload)이 커밋된 기본값이었다.
#   원칙: 수치는 이름 있는 리터럴 상수로 유지한다 — check_model_consistency.py 의 SHARED 가
#         리터럴 문자열을 비교하므로 값을 dict 로 감추면 게이트가 깨진다. CASE 는 배선만
#         고르고, 요구되는 상수 조합은 아래에서 **런타임 검증**한다(조용한 오답 금지).
# ============================================================================
CASE = 'clamp_lf'   # 'seed' | 'control_none' | 'clamp_lf' | 'paper_s1'

_CASE_TABLE = {
    #  케이스          케이블  클램프패치  클램프하중      코너각
    'seed':         dict(cables=True,  clamps=False, clamp_load='none', corner_dir='hf'),
    'control_none': dict(cables=True,  clamps=True,  clamp_load='none', corner_dir='hf'),
    'clamp_lf':     dict(cables=False, clamps=True,  clamp_load='live', corner_dir='hf'),
    'paper_s1':     dict(cables=False, clamps=False, clamp_load='none', corner_dir='paper45'),
}


def resolve_case(case):
    """케이스 -> 배선 dict. **순수 함수** — Abaqus 없이 로컬에서 전 케이스를 실행 검증한다."""
    if case not in _CASE_TABLE:
        raise RuntimeError("CASE must be one of %s (got %r)"
                           % ('|'.join(sorted(_CASE_TABLE)), case))
    return dict(_CASE_TABLE[case])


_ROUTE = resolve_case(CASE)
HAS_CABLES = _ROUTE['cables']       # 꼭짓점을 케이블(T3D2)+Tie 로 잇고 끝단을 구동하는가
HAS_CLAMPS = _ROUTE['clamps']       # 클램프 강체패치를 만드는가
CLAMP_LOAD = _ROUTE['clamp_load']   # 'none' | 'live' (live = 클램프를 집중하중으로 LIVE 에)
CORNER_DIR = _ROUTE['corner_dir']   # 'hf' | 'paper45'


def _require_constant(name, have, want):
    """케이스가 요구하는 상수 조합. 무효 조합은 조용한 오답 대신 즉시 예외로 막는다."""
    if have != want:
        raise RuntimeError("CASE=%r 은 %s = %r 을 요구합니다 (현재 %r).\n"
                           "  -> 상수부를 고치거나 CASE 를 바꾸세요." % (CASE, name, want, have))


if CASE == 'paper_s1':
    _require_constant('PRETENSION_MODE', PRETENSION_MODE, 'paper3')
    _require_constant('CORNER_ANGLE_DEG', CORNER_ANGLE_DEG, 45.0)
else:
    _require_constant('PRETENSION_MODE', PRETENSION_MODE, 'corner2')
    _require_constant('CORNER_ANGLE_DEG', CORNER_ANGLE_DEG, 28.6)

CLAMP_DC = 0.5   # d_c — 클램프 당김 비율 (CLAMP_PULL = 코너 당김 * d_c, A-route 와 동일)

# ---- 클램프 당김 **방향** (driven 과 cload 가 함께 쓴다) ----
#   교수님 지적(2026-10-07): 클램프 없이 케이블만 쓸 때는 각 케이블에 걸리는 변위가
#   힘평형을 이루도록 각도를 맞췄다. 클램프가 있는 모델도 마찬가지로 맞춰야 base state 가
#   안정하고 선형 모드를 찾을 수 있다.
#   실측한 현재 불평형 (DISP=1.8e-5, d_c=0.5):
#       꼭짓점 2개 y합 = -1.723e-05 m (아래로)      <- u2=-DISP*sin(28.6)
#       클램프  2개 y합 = +1.273e-05 m (위로)        <- u2=+CLAMP_PULL/√2  (normal)
#     => 두 하중계가 **반대 방향**이다. Top 완전고정이 순합(-4.5e-06 m)을 반력으로
#        흡수하지만, 막 내부는 '위로 잡아당겨진' 면내 응력 상태가 된다.
#   세 후보 (실측 계산):
#       'normal'   : CL=(-1,+1)/√2, CR=(+1,+1)/√2 , d_c=0.5 그대로
#                    클램프 y합 +1.273e-05 m / 꼭짓점 y합 -1.723e-05 m
#                    -> 남는 -4.505e-06 m 를 Top 완전고정이 반력으로 흡수.
#                    **HF 와 동일한 실제 부품 형상**이다(기본값).
#       'balanced' : 'normal' 방향 그대로, **크기만** 힘평형으로 자동 조정한다.
#                    목표: 클램프 y합 + 꼭짓점 y합 = 0  =>  Top 반력 = 0
#                    필요 CLAMP_PULL = DISP*sin(28.6)/DIR_y  =>  d_c_eff = √2*sin(28.6) = 0.677
#                    교수님이 말한 '케이블이 거는 힘의 합 = 0' 에 가장 부합한다.
#       'inward'   : CL=(+1,0), CR=(-1,0)  -> y 성분 0
#                    y 불평형은 없애지만 사선변 중점을 **안쪽으로 눌러** 면내 압축을 키운다
#                    => 좌굴 하중을 낮춰 오히려 더 불안정해질 수 있다. 케이블 장력(인장)과 어긋난다.
CLAMP_DIR = 'normal'
if CLAMP_DIR not in ('normal', 'balanced', 'inward'):
    raise RuntimeError("CLAMP_DIR must be 'normal'|'balanced'|'inward' (got %r)" % (CLAMP_DIR,))

# ---- [2026-10-07 삭제] 클램프 케이블 경유 구동('driven') ----
#   기록상 막다른 길이라 CASE 에서 제거했다: T3D2 Truss 는 축방향 강성만 주고 Tie 는 위치
#   공차만 정하므로 RP 의 저강성 모드가 subspace 를 삼킨다(500 -> 45~151, CONVERGED=0,
#   실측 888~1256 SYSTEM 음수). RP 에 u3=0 을 더해도 카운터가 변하지 않았다.
#   클램프는 이제 케이블 없이 집중하중으로 준다(CASE='clamp_lf'). 근거는 노트 §18.

# ---- CLOAD 라우트 (CLAMP_MODE='cload') ----
#   [2026-10-07] 사용자 제안: 클램프를 케이블로 구동하는 대신 **집중하중**으로 준다.
#   근거 1 — 진단: 실측에서 passive(클램프 하중 없음 + u3=0)는 48/성공인데
#     driven(클램프 케이블 추가)은 d_c=0.5/0.25 모두 SYSTEM 888, DIFFERENTIAL 27461~27761,
#     lambda1 ~ 4e-08 로 **d_c 와 무관하게** 같은 값이 나온다. RP 에 u3=0 을 더해도 변하지 않았다.
#     그리고 mode.py 는 꼭짓점 케이블 3개로 성공한다 -> 꼭짓점 케이블은 검증됐고
#     변수는 **클램프 케이블(Cable_CL/Cable_CR + Tie_CL/Tie_CR)** 하나만 남는다.
#   근거 2 — Galhofo 2022 재현: 이 논문은 '케이블 없는 only membrane' 에서 좌굴모드를 뽑아
#     케이블 포함 모델에 임퍼펙션으로 주입했다(2-모델 레시피). 모드 추출에는 케이블이 필요 없다.
#   근거 3 — *BUCKLE 의 본질: BUCKLING LOAD ESTIMATE = ("DEAD" LOADS) + lambda * ("LIVE" LOADS).
#     즉 **Buckle 스텝의 LIVE 하중 크기는 lambda 에 흡수된다.** 따라서 클램프 하중은
#     임의의 기준값 CLAMP_F0 로 주면 되고, 나온 lambda 가 곧 '클램프 하중 계수' 다.
#     이 lambda 가 x_c/d_c 에 따라 변하는 것이 MFBO 목적함수의 재료가 된다.
#   구성: GlobalTension(꼭짓점)은 **그대로** 둔다(검증된 경로, base state 보존).
#         클램프는 passive 처럼 RP 면외만 구속하고, Step-Buckle 에서 CLOAD 로 LIVE 를 준다.
CLAMP_F0 = 1.0   # N — (구 CLOAD 라우트의 값) **현행 clamp_lf 는 변위로 구동한다(§19)**:
#   (lambda 가 전체 배수를 흡수한다) 절대값 자체는 결과를 바꾸지 않는다. 비율 CLAMP_F0/CORNER_F0
#   을 로그에 찍는다.

# ---- 꼭짓점 CLOAD 크기 [N] (CASE='clamp_lf') ----
#   실측 RF(꼭짓점 케이블 3개 반력, DISP_GLOBAL=1.8e-5): 0.03032 / 0.03032 / 0.02903 N
#   => 약 0.03 N. 'cable' 라우트의 프리텐션과 같은 크기에서 출발한다.
CORNER_F0 = 0.03

# ---- DEAD/LIVE 배분 (CASE='clamp_lf') ----
DEAD_FRAC = 0.018   # [-] 프리텐션 단계가 운용 하중 세트 P 의 어느 비율까지 올리는가.
#   BUCKLING LOAD = DEAD + lambda*LIVE 이고 DEAD = DEAD_FRAC x P, LIVE = (1-DEAD_FRAC) x P 이므로
#   **lambda = 1 이면 적용 하중 = 정확히 P = 운용점**이다 (lambda > 1 = 좌굴 전 = 설계 여유).
#   기준점은 HF 의 실제 램프 비율: DISP_GLOBAL / GLOBAL_FINAL = 1.8e-5 / 1.0e-3.
#   0.0 으로 두면 base state 가 무응력이 되어 면외 강성이 SIGMA0 만 남고, 굽힘강성 없는 막의
#   영에너지 모드가 subspace 창을 먹는다(기록된 실패 유형). 최종값은 base_state_probe 의
#   압축 면적비·max|u3|·전막 평균 응력을 보고 정한다(= 퇴화하지 않는 최소값). 1열 스윕 대상.

# ---- 스텝 체인: GlobalTension 하나에 꼭짓점 + 클램프를 모두 넣는다 ----
#   [2026-10-07 사용자 결정] 버클모드 추출에서는 클램프를 **별도 텐션 스텝으로 분리하지
#   않는다**. 글로벌 텐션 스텝에 클램프도 함께 포함시킨다(HF 처럼 3스텝으로 늘리지 않는다).
#       Step-GlobalTension   꼭짓점 + 클램프 (DEAD = base state)
#       Step-Buckle          고유값 추출   (LIVE = perturbation)
#   주의: 클램프 하중을 Break 스텝(=perturbation)에 걸면 LIVE 로만 작용해 base state 에
#   압축을 만들지 못한다. 그래서 GlobalTension 쪽에 있어야 한다(이번 수정의 핵심).
#
# ---- base state 진단에 쓸 스텝 (base_state_probe 의 인자) ----
#   클램프 하중이 GlobalTension 에 있으므로 이 스텝을 보면 클램프 기여가 진단에 잡힌다.
BASE_STATE_STEP = 'Step-GlobalTension'

# ---- 판정 임계값 (계획 §4: 3번 런 결과 뒤에 확정한다 — 지금은 잠정값) ----
SPREAD_MIN_PCT = 2.0   # lambda 스프레드 하한(중복 근 배제).
#   실측 대조: 클램프-프리 성공 런 26.9~31.3% (구별되는 모드) /
#              클램프 포함 국소 중복 근 0.28% (사실상 같은 값 100개).
#   2.0 은 그 사이의 잠정값이고, 3번 런의 실측 분포를 보고 확정한다.
U3_FRAC_MIN = 0.10     # 모드의 면외 비율 하한(면내 모드 배제). O(1) vs O(1e-16).
#   buckle_mode_report.U3_FRAC_MIN 과 같은 값이어야 한다(형상 도구와 게이트가 같은 기준).
CLAMP_ZONE_R = 1.0     # |u3| 무게중심이 클램프 부착점에서 이 거리 안이면 '압축영역 분포' 로 본다.

# ---- 좌굴 스텝 (run_abaqus_cable.py 에서 완주가 확인된 설정과 동일) ----
PERTURBATION = 0.01  # m — 좌굴 스텝의 prescribed 변위(증분 응력 -> K_delta).
#   [mode 정합] mode.py 와 HF 는 0.01 을 쓴다. buckle 만 1e-4 였다.
#   ⚠️ CASE='paper_s1' 에는 이 값을 쓰지 않는다 — 아래 PAPER_S1_LIVE_M 참조
#      (논문 자신의 좌굴 스텝 크기에서만 lambda 비교가 성립한다).
#   [2026-10-06 이력] 1e-4 에서 CONVERGED=100, 0.01 에서 CONVERGED=0 이었던 관측은
#   클램프가 있는 상태에서 얻은 것이다. 이제 나머지 조건을 mode.py 와 맞췄으므로
#   0.01 로 두고 다시 판정한다 — 클램프만 다른 상태에서 비교하려면 이 값도 같아야 한다.

# ---- 논문 baseline 좌굴 스텝 크기 (CASE='paper_s1' 전용) ----
PAPER_S1_LIVE_M = 0.10   # m — Galhofo 2022 의 좌굴 스텝은 꼭짓점 변위를 0.10 m 까지 올린다.
#   lambda 는 LIVE 크기에 반비례하므로, **논문 자신의 크기**에서만 논문의 lambda
#   (대역 3.18e-4)와 직접 비교가 성립한다. PERTURBATION(0.01, mode.py 와의 SHARED 계약)과
#   분리해 둔 이유: 계약 항목을 케이스별로 바꾸면 HF<->buckle<->mode 게이트가 깨진다.
#   [주의] paper_s1 의 LIVE 는 **prescribed 변위**(꼭짓점 3개)이지 CLOAD 가 아니다.
# ---- 추출 창 (래더 L4) ----
#   규칙(mode.py L26): '요청 고유값 수 > base state 의 음수 고유값 수' 여야 양수
#   좌굴모드가 subspace 창에 들어온다. 그런데 클램프(passive)를 켜면 음수 고유값이
#   48 개로 늘고, N_EIG=100 / vectors=250 에서는 ***ERROR: THE EIGENVALUES CANNOT BE
#   FOUND (INSTABILITIES IN THE BASE STATE) 로 죽는다. 클램프를 끄면(none) 같은
#   조건에서 CONVERGED=4 다 — 즉 창이 아니라 음수 개수가 문제다.
#   그래서 창을 넓혀(D) 음수 48 개를 넘겨 본다. 실패 비용이 40 초라 판정이 빠르다.
#   창을 정하는 규칙은 이제 두 조건의 **교집합**이다(parse_buckle_log._window_hint 가 같은 판정을 찍는다):
#     (a) N_EIG_BUCKLE > base state 의 SYSTEM 음수 고유값 개수
#     (b) N_EIG_BUCKLE <= 실제 subspace 차원 ("REDUCED TO n" 의 n)
#   (a)는 음수 모드 뒤에 양수 모드가 들어오게 하는 조건, (b)는 원리적 상한이다.
N_EIG_BUCKLE = 10      # 추출 요청 고유값 수  [2026-10-07 이력] 100 -> 200 -> (REDUCED TO 14) -> 10
BUCKLE_VECTORS = 20    # subspace 기저 벡터 수  [2026-10-07 이력] 250 -> 500 -> (REDUCED TO 14) -> 20
#   ⚠️ 실측 반증(2026-10-07, cload 라우트): 요청 500 -> "VECTORS IN SUBSPACE IS REDUCED TO 14".
#      즉 실제 기저는 14 차원인데 200 개를 요청했다 -> ITERATION 마다 수렴 수가 출렁이고
#      (4,5,2,4,4,2) 결국 ***ERROR: THE EIGENVALUES CANNOT BE FOUND.
#      "N_EIG <= 실제 subspace 차원" 이 조건이다. 창을 키우는 것은 역효과다.
#      => 창을 키우는 대신 **이 두 상수를 직접 줄여서** 시험한다 (CLI 로 받지 않는다).
#         예: N_EIG_BUCKLE = 10 / BUCKLE_VECTORS = 20   (우리에게 필요한 것은 양수 모드 몇 개뿐)
#      [이력] 잠시 `neig=`/`vec=` CLI 로 열었다가 되돌렸다(2026-10-07). Windows 런처가
#      `key=value` 의 `=` 를 소비해 값만 도착했고, 파서가 그것을 disp_m/clamp_pull_m 으로 받아
#      조용히 다른 모델(10 m 당김)을 만들었다. 인자 표면을 늘리지 않는 편이 안전하다.
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
    """맨 뒤의 숫자 1~4개를 <x_c> [disp_m] [clamp_pull_m] [clamp_excl_r_m] 로 읽는다.

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
    if not (1 <= len(nums) <= 4):
        raise RuntimeError(
            'Expected 1 to 4 trailing numeric arguments '
            '(<x_c> [disp_m] [clamp_pull_m] [clamp_excl_r_m]), got %r '
            '(first non-numeric token from the end: %r).\n'
            'Usage: abaqus cae noGUI=run_abaqus_buckle.py -- <x_c> [disp_m] [clamp_pull_m] '
            '[clamp_excl_r_m]'
            % (nums, stop))
    return (nums[0],
            (nums[1] if len(nums) > 1 else None),
            (nums[2] if len(nums) > 2 else None),
            (nums[3] if len(nums) > 3 else None))


x_c, _disp_arg, _clamp_arg, _excl_arg = parse_args(sys.argv)
DISP = DISP_GLOBAL if _disp_arg is None else _disp_arg
# 클램프 법선 당김 — 기본은 코너 당김과 같은 비율(A-route 의 d_c), CLI 3번째 인자로 직접 지정 가능.
#   예: -- 0.5 0 5e-5  -> 코너 미구동 + 클램프만 5e-5 m (클램프 단독 구동 진단)
CLAMP_PULL = (DISP * CLAMP_DC if _clamp_arg is None else _clamp_arg)
# [2026-10-07 수정] CLAMP_PERT 를 CLAMP_DC 가 아니라 **실제 CLAMP_PULL** 에 비례시킨다.
#   CLI 3번째 인자로 CLAMP_PULL 을 직접 주면 예전 코드는 섭동만 CLAMP_DC(0.5) 기준으로
#   남아 하중과 섭동이 어긋났다(d_c 스윕이 조용히 오염된다).
if _excl_arg is not None:
    CLAMP_EXCL_R = _excl_arg
    if CLAMP_EXCL_R < 0.0:
        raise RuntimeError('CLAMP_EXCL_R=%r < 0' % (CLAMP_EXCL_R,))
CLAMP_DC_EFF = CLAMP_PULL / DISP if DISP > 0.0 else 0.0
CLAMP_PERT = PERTURBATION * CLAMP_DC_EFF   # 좌굴 스텝 클램프 섭동 [m]

if not (0.03 <= x_c <= 0.95):
    # 원장 M-8: x_c < 0.028 이면 클램프 패치가 정점 패치와 겹치고, x_c ~ 1 이면
    # 모서리 패치와 겹친다 -> 조용히 다른 모델이 된다. 경고만 찍고 진행한다.
    print("%s WARNING: x_c=%g is outside the safe band (0.03, 0.95) — "
          "the clamp patch may overlap a vertex patch (ledger M-8)." % (TAG, x_c))
if DISP < 0.0:
    raise RuntimeError('disp_m must be >= 0 (got %r).' % (DISP,))
# [2026-10-07] 상한 가드. DISP 는 이 모델의 크기(20 m)에 비해 작아야 한다.
#   정상 운용점은 DISP_GLOBAL = 1.8e-05 m (18 um), 논문 좌굴모델도 모서리 5e-4 m / 정점 1e-3 m 다.
#   실사고: 런처가 `neig=10` 의 `=` 를 소비해 값 `10` 이 disp_m 으로 들어왔고, 10 m 당김으로
#   잡이 돌았다(잡 이름 d10000000um). 오류가 아니라 **엉뚱한 해석**이라 눈치채기 어렵다.
DISP_MAX = 1.0e-3   # m
if DISP > DISP_MAX:
    raise RuntimeError(
        'disp_m=%r m is above the sanity cap DISP_MAX=%r m. '
        'This model is 20 m x 10 m with a 5e-6 m membrane; the operating point is 1.8e-05 m. '
        'Check the argument order: <x_c> [disp_m] [clamp_pull_m] [clamp_excl_r_m]'
        % (DISP, DISP_MAX))
if DISP == 0.0:
    # 코너를 구동하지 않는다 = 클램프 단독 구동 케이스 (진단용).
    # Initial 단계의 Disp_Control_Right/Left(u1=u2=0) 값이 그대로 유지된다.
    print("%s DISP=0 : 코너 미구동 — 클램프 단독 구동 케이스로 진행한다." % (TAG,))
print("%s CASE=%s (cables=%s clamps=%s clamp_load=%s corner_dir=%s)"
      % (TAG, CASE, HAS_CABLES, HAS_CLAMPS, CLAMP_LOAD, CORNER_DIR))
print("%s d_c=%.3g -> CLAMP_PULL=%.4e m (buckle pert %.3e m)"
      % (TAG, CLAMP_DC, CLAMP_PULL, CLAMP_PERT))

V_CL = clamp_coord_L(x_c)
V_CR = clamp_coord_R(x_c)

# 하중 각도 — HF/케이블 방향(28.6도). 아래 첫 줄은 HF↔buckle SHARED 계약 항목이라
#   **리터럴로 유지한다**(check_model_consistency 가 이 문자열을 직접 비교한다).
angle_deg = 28.6
if CORNER_DIR == 'paper45':
    angle_deg = 45.0    # 논문 좌굴모델(45도) 재현 케이스 전용 — HF 와 같은 각도가 아니다.
angle_rad = np.deg2rad(angle_deg)
cos_val = float(np.cos(angle_rad))
sin_val = float(np.sin(angle_rad))

# 클램프 당김 단위벡터 (DIR_CL = 좌측 클램프, DIR_CR = 우측 클램프)
_s2 = 2.0 ** 0.5
if CLAMP_DIR == 'inward':
    DIR_CL, DIR_CR = (+1.0, 0.0), (-1.0, 0.0)   # 사선변 안쪽 수평 -> y 성분 0
else:
    DIR_CL, DIR_CR = (-1.0 / _s2, +1.0 / _s2), (+1.0 / _s2, +1.0 / _s2)

# 'balanced': 클램프 y합 + 꼭짓점 y합 = 0 이 되도록 크기를 맞춘다 (Top 반력 0).
#   꼭짓점 y합 = -2*DISP*sin(theta)  이므로  2*CLAMP_PULL*DIR_y = 2*DISP*sin(theta).
if CLAMP_DIR == 'balanced':
    _theta = np.deg2rad(angle_deg)
    CLAMP_PULL = DISP * np.sin(_theta) / DIR_CL[1]
    CLAMP_PERT = PERTURBATION * np.sin(_theta) / DIR_CL[1]
    print("%s [balanced] CLAMP_PULL %.4e -> %.4e m (d_c_eff=%.4f)"
          % (TAG, DISP * CLAMP_DC, CLAMP_PULL, CLAMP_PULL / DISP))

print("%s x_c=%g -> V_CL=%s V_CR=%s" % (TAG, x_c, V_CL, V_CR))
print("%s d_c_eff=%.4f (CLAMP_PULL=%.4e, CLAMP_PERT=%.4e)  CLAMP_EXCL_R=%.4f m"
      % (TAG, CLAMP_DC_EFF, CLAMP_PULL, CLAMP_PERT, CLAMP_EXCL_R))
print("%s CLAMP_DIR=%s -> DIR_CL=%s DIR_CR=%s  (y성분 합=%+.3e m)"
      % (TAG, CLAMP_DIR, DIR_CL, DIR_CR, CLAMP_PULL * (DIR_CL[1] + DIR_CR[1])))
print("%s corner pull DISP=%.4e m (= alpha %.4g x 5e-5 m)  perturbation=%.3e m"
      % (TAG, DISP, DISP / 5.0e-5, PERTURBATION))
print("%s buckle step: solver=%s numEigen=%d vectors=%d"
      % (TAG, BUCKLE_SOLVER, N_EIG_BUCKLE, BUCKLE_VECTORS))


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
    # [2026-10-07] 예전에는 같은 조건을 두 번 검사해 두 번째 WARNING 이 **도달 불가**였다
    #   (첫 검사가 먼저 raise 한다). 한 번만 검사하고 실패 이유를 메시지에 구분해 넣는다.
    if job.status == ABORTED or not os.path.exists(odb_file):
        print_failure_cause(job_name)
        raise RuntimeError('Job %s 실패 (Status=%s): ABORTED 이거나 .odb 가 없다. '
                           'sys.exit 대신 예외로 올린다 — CAE noGUI 러너에서 sys.exit 은 '
                           '종료코드 0으로 보인다.' % (job_name, str(job.status)))
    if not job_completed_ok(job_name):
        print_failure_cause(job_name)
        raise RuntimeError('Job %s 실패 (Status=%s): .sta/.msg 에 완주 문자열(%s)이 없다 '
                           '(중도 중단 의심; odb 존재만으로는 판정 불가 — R-12).'
                           % (job_name, str(job.status), ' / '.join(_COMPLETION_STRINGS)))
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


def msg_cause_from_text(msg_text, sta_text='', limit=8):
    """실패한 잡의 원인 줄을 뽑는다. **순수 함수** — 합성 .msg/.sta 로 Abaqus 없이 검증한다.

    왜: 원인은 .msg 에만 있는데 예전에는 '실패 원인 미상' 만 찍혀 왕복이 늘었다.
    ***ERROR 를 WARNING 보다 우선한다(원인에 가깝다). ERROR 가 없으면 .sta 마지막 줄로 대체한다.
    """
    _keys = ('***ERROR', 'TOO MANY ATTEMPTS', 'HAS NOT BEEN COMPLETED',
             'THE EIGENVALUES CANNOT BE FOUND', 'HAS BEEN TERMINATED',
             'NUMERICAL SINGULARITY', 'EXCESSIVE DISTORTION', '***WARNING')
    out, seen = [], set()
    for _raw in (msg_text or '').splitlines():
        _s = _raw.strip()
        if not _s or _s in seen:
            continue
        for _k in _keys:
            if _k.upper() in _s.upper():
                seen.add(_s)
                out.append(_s[:200])
                break
    _errs = [s for s in out if '***ERROR' in s.upper() or 'ERROR:' in s.upper()]
    if _errs:
        return _errs[:limit], '***ERROR %d줄 (전체 후보 %d줄)' % (len(_errs), len(out))
    if out:
        return out[:limit], '***ERROR 없음 — 경고/기타 후보 %d줄' % len(out)
    _sta = [s.strip() for s in (sta_text or '').splitlines() if s.strip()]
    if _sta:
        return _sta[-3:], 'ERROR 줄 없음 — .sta 마지막 %d줄로 대체' % min(3, len(_sta))
    return [], '.msg/.sta 가 비었거나 없다'


def print_failure_cause(job_name, limit=8):
    """실패 원인을 콘솔에 찍는다. .msg 원문은 ASCII 라 **콘솔 한글 깨짐과 무관하게** 읽힌다."""
    _read = {}
    for _ext in ('msg', 'sta'):
        _fn = '%s.%s' % (job_name, _ext)
        _read[_ext] = ''
        if os.path.exists(_fn):
            try:
                with open(_fn, 'r', errors='replace') as _f:
                    _read[_ext] = _f.read()
            except (IOError, OSError) as _e:
                print('%s [FAIL] %s 를 읽지 못함: %s' % (TAG, _fn, _e))
    _pick, _summary = msg_cause_from_text(_read['msg'], _read['sta'], limit=limit)
    print('%s [FAIL] 원인 판정: %s' % (TAG, _summary))
    for _s in _pick:
        print('%s   | %s' % (TAG, _s))
    print('%s [FAIL] 전문 파일: %s.msg / %s.sta / %s.dat'
          % (TAG, job_name, job_name, job_name))
    return _pick


def report_job(job_name):
    """결과 판정에 필요한 줄만 콘솔에 찍고 <job>.diag.txt 로도 남긴다.

    좌굴 잡의 결론은 두 곳에만 있다.
      .dat -> MODE NO / EIGENVALUE 표 (고유값 원문)
      .msg -> CONVERGED / REQUESTED BY THE USER / CANNOT BE FOUND / REDUCED TO
    전체 로그를 다 찍으면 정작 필요한 줄이 묻히므로 그 줄만 뽑는다.
    """
    out = ['===== job_completed_ok = %s =====' % job_completed_ok(job_name)]

    msg_keys = ('CONVERGED', 'REQUESTED BY THE USER', 'CANNOT BE FOUND', 'ERROR MESSAGES',
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



def connect_cable(a, name, part, sail_corner, vector_dir, radius=1e-4, cable_len=None):
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
    if cable_len is None:
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

    if HAS_CABLES:
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
    if not HAS_CLAMPS:
        # 'cload' 라우트는 CLAMP_MODE 와 무관하게 클램프 패치가 필요하다(집중하중 대상).
        rp_cl_obj = rp_cr_obj = None
    else:
        rp_cl_obj, rp_cl_reg = create_rigid_patch(a, inst_memb, 'CL', V_CL, radius=0.2)
        rp_cr_obj, rp_cr_reg = create_rigid_patch(a, inst_memb, 'CR', V_CR, radius=0.2)

    # 케이블 연결: 정점은 위로, 두 아래 모서리는 각 케이블 축(28.6도) 방향 — mode.py 와 동일
    if HAS_CABLES:
        start_c1, end_c1 = connect_cable(a, 'Cable_Top', p_cable_top, V1, (0.0, 1.0, 0.0))
        start_c2, end_c2 = connect_cable(a, 'Cable_Right', p_cable_bot, V2, (cos_val, -sin_val, 0.0))
        start_c3, end_c3 = connect_cable(a, 'Cable_Left', p_cable_bot, V3, (-cos_val, -sin_val, 0.0))

    a.Set(name='RP_Top_Set', referencePoints=(a.referencePoints[rp1_obj.id],))
    a.Set(name='RP_Right_Set', referencePoints=(a.referencePoints[rp2_obj.id],))
    a.Set(name='RP_Left_Set', referencePoints=(a.referencePoints[rp3_obj.id],))
    if HAS_CLAMPS:
        a.Set(name='RP_CL_Set', referencePoints=(a.referencePoints[rp_cl_obj.id],))
        a.Set(name='RP_CR_Set', referencePoints=(a.referencePoints[rp_cr_obj.id],))
    # Tie (RP <-> 케이블 시작단) — mode.py 와 동일. cload 라우트에는 케이블이 없다.
    if HAS_CABLES:
        my_model.Tie(name='Tie_Top', main=a.sets['RP_Top_Set'], secondary=start_c1,
                     positionToleranceMethod=COMPUTED)
        my_model.Tie(name='Tie_Right', main=a.sets['RP_Right_Set'], secondary=start_c2,
                     positionToleranceMethod=COMPUTED)
        my_model.Tie(name='Tie_Left', main=a.sets['RP_Left_Set'], secondary=start_c3,
                     positionToleranceMethod=COMPUTED)
    a.regenerate()

    # ---- 구동점 추상화 ----
    #   'cable' : 케이블 끝단(end_c*)에 변위를 준다.
    #   'cload' : 케이블이 없으므로 **강체패치 RP** 에 하중/구속을 직접 건다.
    if HAS_CABLES:
        A1, A2, A3 = end_c1, end_c2, end_c3
    else:
        A1, A2, A3 = (a.sets['RP_Top_Set'], a.sets['RP_Right_Set'], a.sets['RP_Left_Set'])
    print("%s HAS_CABLES=%s -> 구동점 %s" % (TAG, HAS_CABLES,
          'CABLE ENDS' if HAS_CABLES else 'RIGID PATCH RPs'))

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
    #   클램프 하중은 GlobalTension 스텝에 이미 들어 있다(별도 ClampTension 스텝 없음).
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
    # Top 정점: 완전 고정. 'cable' 이면 케이블 끝단, 'cload' 면 강체패치 RP 에 건다.
    my_model.DisplacementBC(
        name='BC_Anchor_Top',
        createStepName='Initial',
        region=A1,
        u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET
    )

    # Right/Left 정점: u3+회전만 고정하고, 구동은 A2/A3 에 건다.
    my_model.DisplacementBC(
        name='BC_Right_Z', createStepName='Initial', region=A2,
        u3=SET, ur1=SET, ur2=SET, ur3=SET
    )
    my_model.DisplacementBC(
        name='BC_Left_Z', createStepName='Initial', region=A3,
        u3=SET, ur1=SET, ur2=SET, ur3=SET
    )
    disp_a = disp
    if HAS_CABLES or CORNER_DIR == 'paper45':
        # 변위 구동: 'cable' 은 케이블 끝단, 'paper45' 는 꼭짓점 강체패치 RP (논문 레시피).
        my_model.DisplacementBC(name='Disp_Control_Right', createStepName='Initial',
                                region=A2, u1=SET, u2=SET)
        my_model.DisplacementBC(name='Disp_Control_Left', createStepName='Initial',
                                region=A3, u1=SET, u2=SET)
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
    # 'cload' 라우트의 꼭짓점 하중은 아래 CLAMP_MODE 분기(LOAD_MODE=='cload')에서
    #   일괄 생성한다. 여기서 또 만들면 같은 이름이 두 createStepName 으로 등록되어
    #   나중 것이 0 으로 덮어쓴다(스텁 하네스로 적발: CF_Corner_* 가 Initial/GlobalTension
    #   양쪽에 나타나 최종값이 0 이 됐다).
    # [paper3] 위 꼭짓점도 케이블 축(+y)으로 당긴다 — mode.py 와 동일.
    #   corner2(기본)에서는 BC_Anchor_Top 의 완전고정이 그대로 유지된다.
    if PRETENSION_MODE == 'paper3':
        my_model.boundaryConditions['BC_Anchor_Top'].setValuesInStep(
            stepName='Step-GlobalTension',
            u1=0.0, u2=disp_a * DISP_TOP_OVER_CORNER
        )
    _top_note = (' / 위 %.4e m' % (disp_a * DISP_TOP_OVER_CORNER)
                 if PRETENSION_MODE == 'paper3' else ' (Top 완전고정)')
    print("%s 프리텐션 = %s : 아래 %.4e m%s"
          % (TAG, PRETENSION_MODE, disp_a, _top_note))
    # ---- 하중/구속 배선 (CASE 가 결정한다. 삭제된 구 의미는 §18) ----
    #   HAS_CLAMPS=False      : 클램프 패치가 없다 -> 할 일이 없다.
    #   CLAMP_LOAD='none'     : 패치 RP 의 면외(u3)만 막는다. in-plane 자유, 하중 없음.
    #                           실측: none 과 lambda 상대차 6e-6 (구속은 물리를 안 바꾼다).
    #   CLAMP_LOAD='live'     : 꼭짓점과 클램프를 **prescribed 변위**로 구동하고 DEAD/LIVE 로 나눈다.
    #                           (검증된 관용구 = BC. CLOAD 라우트는 완주 기록이 없다 — §18.1,
    #                            CAE 도 스텝 간 하중 수정을 거부했다 — §19)
    if CLAMP_LOAD == 'live':
        # =====================================================================
        # clamp_lf — 꼭짓점과 클램프를 **prescribed 변위**로 구동한다.
        #   [왜 변위인가 — 2026-10-07 실측] 이 저장소에서 **완주한** 구동 경로는 prescribed 변위
        #     BC 뿐이다: HF(run_abaqus.py) · 대조(run_abaqus_cable.py) · 모드소스
        #     (run_abaqus_mode.py) 가 전부 `Disp_Control_*` 를 Initial 에 만들고
        #     `boundaryConditions[...].setValuesInStep(stepName=...)` 로 스텝마다 값을 바꾼다.
        #     CLOAD 는 이 저장소에서 **완주한 적이 없고**(기록된 유일한 CLOAD 조합이 실패 조합,
        #     §18.1), CAE 는 스텝 간 하중 수정 자체를 거부했다:
        #       ValueError: The load does not exist in the specified step or is suppressed...
        #     => 검증된 관용구로 되돌린다. **회계는 그대로다**:
        #        DEAD 를 base state(Step-GlobalTension)에, LIVE 를 좌굴 스텝에 둔다.
        #        setValuesInStep 은 그 스텝의 값을 대체하므로 lambda=1 에서
        #        총 구동 = DEAD + LIVE = 운용 변위 P  =>  lambda = 1 이 운용점.
        #   RP 면외(u3)는 BC_Clamp_* 가 따로 막는다(in-plane 은 아래 변위가 잡는다).
        # =====================================================================
        my_model.DisplacementBC(name='BC_Clamp_CL', createStepName='Initial',
                                region=a.sets['RP_CL_Set'], u3=0)
        my_model.DisplacementBC(name='BC_Clamp_CR', createStepName='Initial',
                                region=a.sets['RP_CR_Set'], u3=0)
        #   크기 출처는 HF 와 동일: 꼭짓점 = DISP, 클램프 = CLAMP_PULL(= CLAMP_DC x DISP).
        _DRV = (('Disp_Control_Right', A2, cos_val, -sin_val, DISP),
                ('Disp_Control_Left', A3, -cos_val, -sin_val, DISP),
                ('Disp_Control_CL', a.sets['RP_CL_Set'], DIR_CL[0], DIR_CL[1], CLAMP_PULL),
                ('Disp_Control_CR', a.sets['RP_CR_Set'], DIR_CR[0], DIR_CR[1], CLAMP_PULL))
        for _nm, _reg, _d1, _d2, _mag in _DRV:
            my_model.DisplacementBC(name=_nm, createStepName='Initial', region=_reg,
                                    u1=0.0, u2=0.0)
        for _nm, _reg, _d1, _d2, _mag in _DRV:
            # DEAD = DEAD_FRAC x P — base state. 0 을 생략하지 않고 **명시**한다.
            my_model.boundaryConditions[_nm].setValuesInStep(
                stepName='Step-GlobalTension',
                u1=DEAD_FRAC * _mag * _d1, u2=DEAD_FRAC * _mag * _d2)
            # LIVE = (1-DEAD_FRAC) x P x PATTERN_SIGN — 좌굴 스텝 = 섭동 패턴(lambda 의 대상).
            my_model.boundaryConditions[_nm].setValuesInStep(
                stepName='Step-Buckle',
                u1=PATTERN_SIGN * (1.0 - DEAD_FRAC) * _mag * _d1,
                u2=PATTERN_SIGN * (1.0 - DEAD_FRAC) * _mag * _d2)
        print("%s [clamp_lf] 변위 구동: 꼭짓점 %.4g m + 클램프 %.4g m (CLAMP_DC=%.3g, 비율 %.3g)"
              % (TAG, DISP, CLAMP_PULL, CLAMP_DC, (CLAMP_PULL / DISP) if DISP else 0.0))
        print("%s [clamp_lf] DEAD_FRAC=%.4g -> DEAD=%.4g x P (GlobalTension) / "
              "LIVE=%.4g x P x PATTERN_SIGN(%+.2g) (Buckle)  => lambda=1 이 운용점"
              % (TAG, DEAD_FRAC, DEAD_FRAC, 1.0 - DEAD_FRAC, PATTERN_SIGN))
    elif HAS_CLAMPS:
        # 구속만 (CLAMP_LOAD='none'): 패치 RP 의 면외(u3)만 막고 in-plane 은 자유. 하중 없음.
        my_model.DisplacementBC(name='BC_Clamp_CL', createStepName='Initial',
                                region=a.sets['RP_CL_Set'], u3=0)
        my_model.DisplacementBC(name='BC_Clamp_CR', createStepName='Initial',
                                region=a.sets['RP_CR_Set'], u3=0)

    # ---- Buckle 스텝: 논문 (c) 3변 u3=0 유지 + 클램프 구간은 면외 구속에서 제외(§9) ----
    _excl_edges = (aba_grid_mesh.clamp_exclude_labels(inst_memb, (V_CL, V_CR), CLAMP_EXCL_R)
                   if HAS_CLAMPS else set())
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

    # ---- Buckle 스텝의 perturbation (케이블 라우트 전용) ----
    #   [2026-10-07 결함수정] PATTERN_SIGN 을 여기에도 곱한다. mode.py 는 곱하는데 buckle 은
    #   빠져 있어 PATTERN_SIGN 이 1.0 을 벗어나면 두 스크립트가 조용히 갈라졌다(검사기는
    #   상수값만 비교하고 '사용 여부'는 보지 않는다). CLAMP_LOAD='live' 의 LIVE CLOAD 는
    #   위에서 이미 정의했으므로 여기서 건드리지 않는다.
    #   [함정] 라우트에 없는 BC 를 참조하면 KeyError 로 모델 생성이 죽는다(스텁 하네스가 적발).
    if HAS_CABLES:
        _ns = PATTERN_SIGN * PERTURBATION
        my_model.boundaryConditions['Disp_Control_Right'].setValuesInStep(
            stepName='Step-Buckle', u1=_ns * cos_val, u2=-_ns * sin_val)
        my_model.boundaryConditions['Disp_Control_Left'].setValuesInStep(
            stepName='Step-Buckle', u1=-_ns * cos_val, u2=-_ns * sin_val)
    if CORNER_DIR == 'paper45':
        # 논문 baseline 은 **논문 자신의 좌굴 스텝 크기**로 올린다 — lambda 가 LIVE 크기에
        #   반비례하므로 그 크기에서만 논문의 3.18e-4 와 직접 비교가 성립한다. 꼭짓점 3개 전부.
        #   (CORNER_DIR=='paper45' 와 HAS_CABLES 는 CASE 로 상호배타 — _require_constant.)
        _ns = PATTERN_SIGN * PAPER_S1_LIVE_M
        my_model.boundaryConditions['Disp_Control_Right'].setValuesInStep(
            stepName='Step-Buckle', u1=_ns * cos_val, u2=-_ns * sin_val)
        my_model.boundaryConditions['Disp_Control_Left'].setValuesInStep(
            stepName='Step-Buckle', u1=-_ns * cos_val, u2=-_ns * sin_val)
        my_model.boundaryConditions['BC_Anchor_Top'].setValuesInStep(
            stepName='Step-Buckle', u1=0.0, u2=_ns * DISP_TOP_OVER_CORNER)
        print("%s [paper_s1] 꼭짓점 3개 prescribed 구동: DEAD=%.4e m (Top x%.4g) / "
              "LIVE=%.4e m — 논문 좌굴 스텝 크기"
              % (TAG, disp_a, DISP_TOP_OVER_CORNER, PAPER_S1_LIVE_M))

    return model_name

# job 이름에 요소 태그 + 클램프 모드를 넣는다. 같은 이름이면 run_job_safely 가 stale
#   산출물을 지워 직전 증거가 사라진다(§11).
# ELEM_TAG 는 요소 상수에서 **유도**한다 — 하드코딩이면 요소를 바꿔도 잡 이름이 그대로여서
#   '산출물 덮임 방지' 장치가 무력해진다(실측: S4 모델에 's4r' 이 박혀 있었다).
ELEM_TAG = {S4: 's4', S4R: 's4r', S3: 's3', M3D4: 'm3d4', M3D3: 'm3d3'}.get(ELEM_CODE_QUAD)
if ELEM_TAG is None:
    raise RuntimeError('ELEM_TAG 를 유도할 수 없는 요소코드입니다: %r. 잡 이름이 케이스를 '
                       '구분하지 못하면 run_job_safely 가 직전 증거를 지운다.' % (ELEM_CODE_QUAD,))
#   d_c_eff 도 이름에 넣는다 — CLAMP_PULL 을 CLI 로 바꿔 d_c 스윕을 하면 DISP 만으로는
#   이름이 겹쳐 run_job_safely 가 직전 증거를 지운다(§11).
_excl_tag = '' if CLAMP_EXCL_R == 0.0 else '_ex%03d' % int(round(CLAMP_EXCL_R * 1000.0))
JOB_NAME = 'Buckle_xc%03d_d%03dum_dc%03d_%s_%s_%s%s' % (int(round(x_c * 100.0)),
                                                       int(round(DISP * 1.0e6)),
                                                       int(round(CLAMP_DC_EFF * 100.0)),
                                                       ELEM_TAG, CASE, CLAMP_DIR, _excl_tag)

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

# ============================================================================
# 판정 재료 수집 — lambda 표 + base state 진단 + ledger 한 행.   (§18)
#   규칙: "완주했는가" 와 "그 base state 가 우리가 주장하는 상태인가" 는 별개 질문이다.
#   lambda 표가 나와도 base state 진단(압축 면적비 / max|u3| / RF 상위)을 함께 읽는다.
# ============================================================================
_M = {}
_HAS_PBL = False
try:
    import parse_buckle_log as _pbl                  # 같은 폴더, 표준 라이브러리만 (DRY)
    _HAS_PBL = True
    _msgp = os.path.join(_RUN, '%s.msg' % JOB_NAME)
    if os.path.exists(_msgp):
        import io as _io1
        with _io1.open(_msgp, encoding='utf-8', errors='replace') as _f:
            _M = _pbl.parse(_f.read())
        print("%s [DIAG] CONVERGED series %s / SYSTEM 음수 %s / subspace %s"
              % (TAG, _M.get('conv'), _M.get('sysNeg'), _M.get('reduced_to')))
except Exception as _e:
    print("%s [DIAG] msg 파싱 실패(무시): %s" % (TAG, _e))

_BASE = {'probe_ok': False}
try:
    _probe = os.path.join(_HERE, 'base_state_probe.py')
    _odb = os.path.join(_RUN, '%s.odb' % JOB_NAME)
    if os.path.exists(_probe) and os.path.exists(_odb):
        print("%s [BASE] probe: %s / %s" % (TAG, _odb, BASE_STATE_STEP))
        _out = subprocess.check_output('abaqus python "%s" "%s" %s'
                                       % (_probe, _odb, BASE_STATE_STEP), shell=True)
        if not isinstance(_out, str):
            _out = _out.decode('utf-8', 'replace')
        _keep = [l.strip() for l in _out.splitlines()
                 if ('압축' in l or 'max|u3|' in l or 'CENTRE' in l or 'mean|u3|' in l)]
        for _l in _keep[:12]:
            print("      | %s" % _l[:170])
        _BASE.update({'probe_ok': True, 'lines': _keep[:40]})
    else:
        print("%s [BASE] skip (probe=%s, odb=%s)"
              % (TAG, os.path.exists(_probe), os.path.exists(_odb)))
except Exception as _e:
    print("%s [BASE] probe 실패(무시): %s" % (TAG, _e))

# ---------------------------------------------------------------------------
# 모드 형상 판정 (C3/C4) + 4조건 게이트.   (§18)
#   왜 형상인가: 임퍼펙션은 z 섭동이므로 면내 모드(u3~0)는 아무것도 주입하지 못하고,
#   사선변에 붙은 경계조건 아티팩트는 클램프 물리가 아니라 틀린 초기결함이 된다.
#   판정 산술은 buckle_mode_report.py 의 순수 함수를 쓴다(같은 저장소, DRY — 로컬 단위검증됨).
# ---------------------------------------------------------------------------
def judge_run(conv_max, lams, modes_ok, c3, c4, spread_min=SPREAD_MIN_PCT,
              u3_min=U3_FRAC_MIN, zone_r=CLAMP_ZONE_R):
    """4조건의 (이름, 통과, 실측근거) 목록을 돌려준다. 각 조건의 근거를 함께 찍기 위해서다."""
    res = []
    _pos = [v for v in (lams or []) if v > 0.0]
    res.append(('C1 수렴/양수lambda',
                bool(conv_max) and len(_pos) >= 1,
                'CONVERGED=%s, 양수 lambda %d개' % (conv_max, len(_pos))))
    if len(_pos) >= 2:
        _mn, _mx = min(_pos), max(_pos)
        _mean = sum(_pos) / float(len(_pos))
        _sp = 100.0 * (_mx - _mn) / abs(_mean) if _mean else 0.0
        res.append(('C2 스프레드', _sp >= spread_min,
                    '%.3f%% (하한 %.1f%%)' % (_sp, spread_min)))
    else:
        res.append(('C2 스프레드', False, '양수 lambda 가 2개 미만'))
    res.append(('C3 면외 성분', bool(c3 is not None and c3 >= u3_min),
                'u3_frac=%s (하한 %.2g)'
                % (('%.3g' % c3) if c3 is not None else '-', u3_min)))
    res.append(('C4 압축영역 분포', bool(c4 is not None and c4 <= zone_r),
                '클램프까지 %s m (상한 %.2g m)'
                % (('%.3f' % c4) if c4 is not None else '-', zone_r)))
    return res


_MET, _VD, _RES = [], {'c3': None, 'c4': None, 'modes_ok': []}, []
try:
    import buckle_mode_report as _bmr                 # _HERE 는 sys.path 에 있다
    _datp = os.path.join(_RUN, '%s.dat' % JOB_NAME)
    _nodes, _disps, _lams, _msgs = _bmr.read_odb_modes(
        os.path.join(_RUN, '%s.odb' % JOB_NAME), 'Step-Buckle',
        n_modes=max(1, min(8, N_EIG_BUCKLE)), instance=INSTANCE_NAME,
        dat_hint=_datp if os.path.exists(_datp) else None)
    for _m in _msgs:
        print("%s [MODE] %s" % (TAG, _m))
    for _d in _disps:
        _m3 = max([abs(v[2]) for v in _d.values()] or [0.0])
        if _m3 > 0:
            _d = dict((k, (v[0] / _m3, v[1] / _m3, v[2] / _m3)) for k, v in _d.items())
        _MET.append(_bmr.mode_metrics(_nodes, _d))
    _ok, _drop = _bmr.judge_modes(_MET)
    _VD['modes_ok'] = _ok
    for _i, _why in _drop:
        print("%s [MODE] 버림 M%-3d %s" % (TAG, _i, _why))
    if _ok:
        _mi = min(_ok) - 1
        _mm = _MET[_mi]
        _VD['c3'] = _mm['u3_frac']
        _VD['c4'] = min(((_mm['cx'] - V_CL[0]) ** 2 + (_mm['cy'] - V_CL[1]) ** 2) ** 0.5,
                        ((_mm['cx'] - V_CR[0]) ** 2 + (_mm['cy'] - V_CR[1]) ** 2) ** 0.5)
        print("%s [MODE] 후보 M%d : u3_frac=%.3g centroid=(%.3f,%.3f) rw=%.3f 클램프까지 %.3f m"
              % (TAG, _mi + 1, _mm['u3_frac'], _mm['cx'], _mm['cy'], _mm['rw'], _VD['c4']))
except Exception as _e:
    print("%s [MODE] 형상 판정 불가(무시): %s" % (TAG, _e))

# ledger: 런 1개 = 1행. 원인 판정은 '한 열만 다른 두 행' 에서만 한다(규칙 2).
try:
    import io as _io
    import json as _js
    import datetime as _dt
    _row = {
        'ts': _dt.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ'),
        'case': CASE, 'x_c': x_c, 'DISP': DISP, 'CLAMP_DC': CLAMP_DC,
        'DEAD_FRAC': DEAD_FRAC, 'CLAMP_DIR': CLAMP_DIR,
        'CORNER_F0': CORNER_F0, 'CLAMP_F0': CLAMP_F0,
        'clamp_load': CLAMP_LOAD, 'corner_dir': CORNER_DIR,
        'window': [N_EIG_BUCKLE, BUCKLE_VECTORS], 'PERTURBATION': PERTURBATION,
        'ELEM_TAG': ELEM_TAG, 'job': JOB_NAME,
        'sysNeg': _M.get('sysNeg'), 'reduced_to': _M.get('reduced_to'),
        'convMax': (max(_M['conv']) if _M.get('conv') else None),
        'wallclock': _M.get('wall'), 'lambdas': [], 'base_probe_ok': _BASE['probe_ok'],
    }
    # lambda 는 **.dat 의 MODE NO / EIGENVALUE 표에서만** 읽는다(요청 수만큼 CONVERGED 한 런의
    #   표가 스펙트럼이고, 실패 런의 ITERATION 목록은 레일리 몫 스냅샷이다).
    #   파서는 parse_buckle_log.dat_lambdas() 에 있다(순수 함수 -> 로컬 단위검증 가능, DRY).
    _datf = os.path.join(_RUN, '%s.dat' % JOB_NAME)
    if os.path.exists(_datf) and _HAS_PBL:
        _txt = _io.open(_datf, encoding='utf-8', errors='replace').read()
        _row['lambdas'] = _pbl.dat_lambdas(_txt, limit=12)
    _RES = judge_run(_row['convMax'], _row['lambdas'], _VD.get('modes_ok'),
                     _VD.get('c3'), _VD.get('c4'))
    _row['verdict'] = {
        'pass': all(_ok for _, _ok, _ in _RES),
        'c3_u3_frac': _VD.get('c3'), 'c4_dist_m': _VD.get('c4'),
        'modes_ok': _VD.get('modes_ok'),
    }
    _led = os.path.join(_RUN, 'buckle_ledger.jsonl')
    with _io.open(_led, 'a', encoding='utf-8') as _f:
        _f.write(_js.dumps(_row, ensure_ascii=False) + '\n')
    print("%s [LEDGER] += %s   (lambda %d개: %s)"
          % (TAG, _led, len(_row['lambdas']), _row['lambdas'][:5]))
except Exception as _e:
    print("%s [LEDGER] 기록 실패(무시): %s" % (TAG, _e))

print("")
print("=" * 78)
for _n, _ok, _why in _RES:
    print("%s [JUDGE] %-4s %-18s %s" % (TAG, 'PASS' if _ok else 'FAIL', _n, _why))
print("%s [JUDGE] 종합: %s   (조건이 비어 있으면 앞의 [DIAG]/[MODE] 줄의 실패 이유를 보라)"
      % (TAG, 'PASS' if (_RES and all(_ok for _, _ok, _ in _RES)) else 'FAIL'))
print("")
print("%s PASS CRITERION (4조건, 위 [JUDGE] 가 실측값과 함께 판정한다):" % TAG)
print("   C1 CONVERGED >= 1 그리고 양수 lambda >= 1")
print("   C2 lambda 스프레드 >= SPREAD_MIN_PCT   (중복 근 배제)")
print("   C3 모드가 면외 성분을 가짐 (u3 비율)   (면내 모드 배제)")
print("   C4 |u3| 무게중심이 클램프 부착점 근처  (경계조건 아티팩트 배제)")
print("%s eigen table : %s.dat" % (TAG, JOB_NAME))
print("%s ledger      : %s%sbuckle_ledger.jsonl" % (TAG, RUN_DIR_NAME, os.sep))
print("%s 모드 형상 판정: abaqus python buckle_mode_report.py %s.odb Step-Buckle"
      % (TAG, JOB_NAME))
print("%s HF 소비(주의): HF(run_abaqus.py)는 기본 IMPERFECTION_MODE='odb_direct' 로" % TAG)
print("                 ..\\ClampFree_Buckle.odb 를 읽는다 — 이 잡을 쓰려면 MODE_SOURCE_ODB 교체.")
print("=" * 78)
