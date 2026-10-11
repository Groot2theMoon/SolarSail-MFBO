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

케이스 (CASE 상수). 예전의 LOAD_MODE x CLAMP_MODE 조합 표면을 대체한다.
    CASE          DEAD (base state)                 LIVE (lambda 의 대상)      답하는 질문
    seed          코너 전량 DISP (케이블 끝단 구동)  코너 PERTURBATION(끝단)    케이블 라우트 corner 모드 — 전역 4모드 재현
    clamp_lf      DEAD_FRAC x P + 클램프 u3=0       (1-DEAD_FRAC) x P          운용 램프의 좌굴 배수(lambda=1=운용점)
    corner_lf     코너 전량 DISP, 클램프 없음        코너 PERTURBATION          그리드에서의 코너 하중 모드(합성 씨드 1/2)
    clamp_live    코너 전량 DISP + 클램프 u3=0,      클램프 DIR x PERTURBATION  클램프 구동 모드(합성 씨드 2/2)
                  클램프 당김 0(명시)
    ⚠️ corner_lf/clamp_live 의 lambda 는 운용점 회계(위 clamp_lf 의 lambda=1)가 **아니다** —
       프로브 크기의 배수다. 모드 형상은 프로브 크기에 불변.
      seed 는 케이블 3개(Cable_Top/Right/Left)를 Tie 로 꼭짓점 패치에 연결하고 끝단을 구동한다
      — 이 케이스만 HAS_CABLES=True 이며 buckle<->mode 정합 계약의 적용 대상이다.
      P = 꼭짓점 + 클램프 하중 세트(CORNER_F0 / CLAMP_F0). clamp_lf 는 케이블 0개다.
      BUCKLING LOAD = DEAD + lambda*LIVE 이므로 **clamp_lf 의 lambda = 1 이 운용점**이다.
      무효 조합(예: 코너는 케이블, 클램프는 CLOAD)은 상수 검증이 런타임 예외로 막는다.

사용법
    abaqus cae noGUI=run_abaqus_buckle.py -- <x_c> [disp_m] [clamp_pull_m] [clamp_excl_r_m] [dead_frac]
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
#   격자는 이제 사용 중이다: USE_GRID_MESH=True -> aba_grid_mesh.import_grid_part
#   (.inp -> PartFromInputFile). fill_part 는 Part.addNodes 가 없어 폐기했다(§19.15).
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
ELEM_CODE_QUAD = S8R5         # 1차 셸(S4/S3). 막 계열(M3D4/M3D3)은 기각 — 재시도 금지(notes §4·§24).
ELEM_CODE_TRI = STRI65          # 위와 같은 이유로 3절점 1차 셸

# ---- 메쉬 종류 ----
#   True = aba_grid_mesh 균일 격자 (S4 9900 + S3 200 = 10100 요소, 노드 10201, 미러쌍 누락 0).
#           논문 Table A.1 토폴로지 일치. 자유 메쉬와 달리 좌우대칭이라 한쪽에만 모드가 생기지 않는다.
#   False = 자유 메쉬(약 1.82만 요소, 구 경로).
USE_GRID_MESH = True
#   orphan(격자) 경로는 **한 이름이 노드용/요소용을 겸할 수 없다**(2026-10-07 실측):
#     입력 처리기 오류: ***ERROR: Unknown part instance node set MEMBRANE-1.ALL
#                       NODE SET ASSEMBLY_MEMBRANE-1_ALL HAS NOT BEEN DEFINED
#     All      = **절점**집합 -> 변위 BC(BC_Stabilize_Z 등).  free 경로에서는 면집합이라 겸용 가능.
#     All_Elem = **요소**집합 -> 섹션 할당 + 초기응력(Stress).  free 경로에서는 All 그대로 쓴다.
_ALL_ELEM = 'All_Elem' if USE_GRID_MESH else 'All'

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

CLAMP_EXCL_R = 0.21   # m — 클램프 패치 노드를 경계 u3 구속에서 **뺀다**.
#   CLAMP_EXCL_R > 0 인 이유: 클램프 패치 노드는 Coupling 으로 RP(u3=0)에 묶여 있어, 같은
#   노드에 BC_Edges_Only_Z(u3=0)이 겹치면 Abaqus 가 이중 구속으로 죽는다(.dat FTL).
#   격자 메쉬는 경계 노드가 규칙적으로 놓여 패치가 그 노드들을 확실히 잡는다(자유 메쉬는 잠복).
#   물리 변화는 없다 — 제외 노드는 Coupling(u3=0)으로 그대로 면외 고정이다.
#   0.0 으로 두면 check_model_consistency 가 실패한다(이중 구속 게이트). CLI 4번째 인자로 덮는다.

V1 = (BASE/2.0, HEIGHT, 0.0) # Top
V2 = (BASE, 0.0, 0.0)        # Right
V3 = (0.0, 0.0, 0.0)         # Left

def clamp_coord_L(x): return (10-10*x, 10-10*x, 0)
def clamp_coord_R(x): return (10+10*x, 10-10*x, 0)

# ---- 프리텐션 (CLI 두 번째 인자로 덮어쓴다) ----
#   [mode 정합] mode.py 는 1.8e-5 m 를 쓴다. buckle 은 5e-5(=5e-6*10) 였다.
#   이 차이가 클램프와 교란되어 '클램프가 원인'인지 판별할 수 없었다.
DISP_GLOBAL = 5e-6    # m — mode.py(ClampFree_Buckle) 와 동일한 운용점

# 코너 하중 각도 [deg] — CASE 가 요구한다('hf' = 28.6 = 우리 HF/케이블 방향, 'paper45' = 논문).
CORNER_ANGLE_DEG = 28.6

# 케이스 (CASE) — '무엇을 DEAD, 무엇을 LIVE 로 두는가'만 고른다.
#   원칙: 수치는 모듈 최상위 리터럴 상수로 유지한다(check_model_consistency 의 SHARED 가 리터럴
#   문자열을 비교한다). CASE 는 배선만 고르고, 요구 상수 조합은 아래에서 런타임 검증한다.
CASE = 'seed'   # 현재 런 대상 — 케이블 구동 corner-family 재현 (은퇴/복원 이력: §19.35, §19.38)

_CASE_TABLE = {
    #  seed: [2026-10-11 복원] 꼭짓점을 케이블(T3D2)+Tie 로 잇고 **케이블 끝단**을 구동한다 —
    #    mode.py(ClampFree_Buckle) 와 같은 하중 경로. 원래 전역 4모드의 재현용(은퇴 경위 §19.35).
    'seed':         dict(cables=True,  clamps=False, clamp_load='none', corner_dir='hf'),
    'clamp_lf':     dict(cables=False, clamps=True,  clamp_load='live', corner_dir='hf'),
    #  corner_lf: 그리드·무클램프. base=코너 전량(DISP), LIVE=코너 프로브(PERTURBATION).
    #    합성 씨드(코너 패밀리) 재료 — clamp_lf 의 lambda=1 운용점 회계가 아니다(lambda=프로브 배수).
    'corner_lf':    dict(cables=False, clamps=False, clamp_load='none', corner_dir='hf'),
    #  clamp_live: base=초기 인장만(코너 전량, 클램프 당김 0 명시), LIVE=클램프 단독 프로브.
    #    합성 씨드(클램프 패밀리) 재료. lambda=프로브 배수. live_driver='clamp' 가 배선을 가른다.
    'clamp_live':   dict(cables=False, clamps=True,  clamp_load='none', corner_dir='hf',
                         live_driver='clamp'),
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
LIVE_DRIVER = _ROUTE.get('live_driver', 'corner')   # 'corner'(기존 케이스) | 'clamp'(clamp_live)


def _require_constant(name, have, want):
    """케이스가 요구하는 상수 조합. 무효 조합은 조용한 오답 대신 즉시 예외로 막는다."""
    if have != want:
        raise RuntimeError("CASE=%r 은 %s = %r 을 요구합니다 (현재 %r).\n"
                           "  -> 상수부를 고치거나 CASE 를 바꾸세요." % (CASE, name, want, have))


#   paper_s1 은퇴(2026-10-08) -> 단일 조건. 되살릴 때는 paper3 / 45.0 분기도 함께 복원한다.
_require_constant('PRETENSION_MODE', PRETENSION_MODE, 'corner2')
_require_constant('CORNER_ANGLE_DEG', CORNER_ANGLE_DEG, 28.6)

CLAMP_DC = 0.0556  # d_c — 클램프 당김 비율 (CLAMP_PULL = 코너 당김 * d_c).
#   [실측] 이 비율이 '모드가 존재하는 창' 안의 값이다 — 5.6%: 4모드 정상 /
#   50%·100%: THE EIGENVALUES CANNOT BE FOUND(음수 고유값 234개, 대각 24528개).
#   f 를 0.5/0.10/0.05 로 낮춰도 실패한다 -> 클램프가 큰 기저에는 선형 좌굴모드가 없다.
#   CLI 3번째 인자로 덮을 수 있고, 잡 이름 dc%03d 가 이 값을 인코딩한다(0.0556 -> dc006).

CLAMP_RP_FIX_ROT = False   # True 면 클램프 RP 의 **회전 3성분도** 고정한다(패치 기울기 차단).
#   현행은 RP 의 u3 만 고정 -> 강체패치가 **기울 수 있다**(회전 자유). 그 틈이 클램프 옆
#   국소 면외모드의 통로가 된다(§19.12). 실제 클램프는 패치 방향을 잡아 주므로 회전 고정이
#   더 물리적일 수 있다. 기본 False = HF 배선과 동일(출력 이름이 _rotfix 로 갈린다).
_CLAMP_ROT_KW = dict(ur1=0.0, ur2=0.0, ur3=0.0) if CLAMP_RP_FIX_ROT else {}

# ---- CLEAN_PERTURBATION ----
#   True = 좌굴 스텝에서 구속집합을 바꾸지 않는다(자유경계 u3=0 은 Initial 에서, BC_Stabilize_Z 는
#   생성 안 함). 선형 섭동 스텝에서 **새 구속의 반력·구속 해제의 스프링백은 하중**이라 LIVE 패턴을
#   오염시킨다 — 실측: lambda 가 f 에 따라 부호까지 뒤집혔다(+8.2e-3 @f=0.5 / -0.16 @f=0.25).
#   False = 종전 배선(비교/복구용).
CLEAN_PERTURBATION = True

# ---- 클램프 당김 방향 (CLAMP_DIR) ----
#   normal(기본, HF 와 동일 부품 형상): CL=(-1,+1)/√2, CR=(+1,+1)/√2, 크기 CLAMP_PULL=d_c x DISP.
#     힘평형은 아니다 — 남는 y 불평형(-4.5e-6 m)을 Top 고정이 반력으로 흡수한다.
#   balanced: 같은 방향, 크기만 힘평형으로(d_c_eff=√2*sin 28.6도=0.677). along/inward 는 미시험·기각.
#   상세 계산과 기각 근거: notes/history.md §12.
CLAMP_DIR = 'normal'
if CLAMP_DIR not in ('normal', 'balanced', 'inward'):
    raise RuntimeError("CLAMP_DIR must be 'normal'|'balanced'|'inward' (got %r)" % (CLAMP_DIR,))


CLAMP_F0 = 1.0   # N — 구 CLOAD 라우트용. 현행 clamp_lf 는 변위 구동(§19):
#   (lambda 가 전체 배수를 흡수한다) 절대값 자체는 결과를 바꾸지 않는다. 비율 CLAMP_F0/CORNER_F0
#   을 로그에 찍는다.

# ---- 꼭짓점 CLOAD 크기 [N] (CASE='clamp_lf') ----
#   실측 RF(꼭짓점 케이블 3개 반력, DISP_GLOBAL=1.8e-5): 0.03032 / 0.03032 / 0.02903 N
#   => 약 0.03 N. 'cable' 라우트의 프리텐션과 같은 크기에서 출발한다.
CORNER_F0 = 0.03

# ---- DEAD/LIVE 배분 (CASE='clamp_lf') ----
DEAD_FRAC = 0.5   # [-] 프리텐션 단계가 운용 하중 세트 P 의 어느 비율까지 올리는가.
#   [실측] f 창: 0.018 -> 무응력 퇴화(음수 126) / 0.900 -> 개시 초과(음수 671) / 0.5 -> 정상.
#   즉 '응력이 살아 있고 개시 앞'인 f 가 필요하다.
#   자기일관성 검증: 개시 = f + (1-f)*lambda1*운용 은 f 와 무관해야 한다. lambda=1 <=> 운용점 도.

# ---- 스텝 체인: GlobalTension 하나에 꼭짓점 + 클램프를 모두 넣는다 ----
#   Step-GlobalTension = 꼭짓점+클램프 (DEAD = base state) / Step-Buckle = 고유값 (LIVE = perturbation).
#   클램프를 perturbation 스텝에 걸면 LIVE 로만 작용해 base state 에 기압축을 못 만든다 — 그래서
#   GlobalTension 쪽에 있어야 한다. base_state_probe 는 이 스텝을 본다(클램프 기여 진단).
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
#   [mode 정합] PERTURBATION=0.01: mode.py/HF 와 동일해야 한다(클램프만 다른 상태로 비교하려면).
#   paper_s1 은 예외로 PAPER_S1_LIVE_M 을 쓴다(논문 자신의 좌굴 스텝 크기에서만 lambda 비교 성립).

# ---- 논문 baseline 좌굴 스텝 크기 (CASE='paper_s1' 전용) ----
PAPER_S1_LIVE_M = 0.10   # m — Galhofo 2022 의 좌굴 스텝은 꼭짓점 변위를 0.10 m 까지 올린다.
#   논문 lambda 비교는 **논문 자신의 좌굴 스텝 크기**(PAPER_S1_LIVE_M)에서만 성립한다 — lambda 는
#   LIVE 크기에 반비례. PERTURBATION(0.01)은 mode.py 와의 SHARED 계약이라 케이스별로 바꾸지 않는다.
# ---- 추출 창: N_EIG_BUCKLE 은 두 조건의 교집합 (parse_buckle_log._window_hint 와 같은 판정) ----
#   (a) N_EIG_BUCKLE > base state 의 SYSTEM 음수 고유값 개수   (b) N_EIG_BUCKLE <= 실제 subspace 차원
N_EIG_BUCKLE = 100      # 추출 요청 고유값 수 (mode.py 와 같은 값이어야 한다 — 계약)
                        #   [2026-10-11 복원] 4 -> 100: 과거 성공 조합(100/250). base state 의 음수
                        #   클러스터를 넘겨야 양수 모드가 subspace 창에 들어온다 — 케이블 런 실측
                        #   "음수 52 < 100 -> CONVERGED=4".
#   [이력] N_EIG=4 / 기저 250: 대조군이 '요청 4 / 기저 250' 으로 수렴했다. 기저 250 = 천장 272 아래.
BUCKLE_VECTORS = 250    # subspace 기저 벡터 수. 실측 천장: 300 을 요청하면
#   ⚠️ "N_EIG <= 실제 subspace 차원" 이 조건이다: 기저를 키워 요청하면 REDUCED 되고 결국
#   EIGENVALUES CANNOT BE FOUND. 창을 키우지 말고 이 두 상수를 직접 줄여 시험한다.
#   (neig=/vec= CLI 는 런처가 '=' 를 삼켜 조용한 오답 — 쓰지 말 것.)
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
    """맨 뒤의 숫자 1~5개를 <x_c> [disp_m] [clamp_pull_m] [clamp_excl_r_m] [dead_frac] 로 읽는다.

    dead_frac 을 주면 DEAD_FRAC(기저 상태가 운용 세트의 몇 %까지 올라가는가)을 덮는다.
    x_c 스윕에서 필요하다: 개시(좌굴) 하중비는 x_c 마다 다르므로 x_c=0.5 에서 통한 0.5 가
    다른 x_c 에서는 **개시 위**가 되어 고유값 추출이 실패한다(§19.24).

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
    if not (1 <= len(nums) <= 5):
        raise RuntimeError(
            'Expected 1 to 5 trailing numeric arguments '
            '(<x_c> [disp_m] [clamp_pull_m] [clamp_excl_r_m] [dead_frac]), got %r '
            '(first non-numeric token from the end: %r).\n'
            'Usage: abaqus cae noGUI=run_abaqus_buckle.py -- <x_c> [disp_m] [clamp_pull_m] '
            '[clamp_excl_r_m] [dead_frac]'
            % (nums, stop))
    return (nums[0],
            (nums[1] if len(nums) > 1 else None),
            (nums[2] if len(nums) > 2 else None),
            (nums[3] if len(nums) > 3 else None),
            (nums[4] if len(nums) > 4 else None))


x_c, _disp_arg, _clamp_arg, _excl_arg, _dead_arg = parse_args(sys.argv)
DISP = DISP_GLOBAL if _disp_arg is None else _disp_arg
# 클램프 법선 당김 — 기본은 코너 당김과 같은 비율(A-route 의 d_c), CLI 3번째 인자로 직접 지정 가능.
#   예: -- 0.5 0 5e-5  -> 코너 미구동 + 클램프만 5e-5 m (클램프 단독 구동 진단)
CLAMP_PULL = (DISP * CLAMP_DC if _clamp_arg is None else _clamp_arg)
# CLAMP_PERT 를 CLAMP_DC 가 아니라 **실제 CLAMP_PULL** 에 비례시킨다.
#   CLI 3번째 인자로 CLAMP_PULL 을 직접 주면 예전 코드는 섭동만 CLAMP_DC(0.5) 기준으로
#   남아 하중과 섭동이 어긋났다(d_c 스윕이 조용히 오염된다).
if _excl_arg is not None:
    CLAMP_EXCL_R = _excl_arg
    if CLAMP_EXCL_R < 0.0:
        raise RuntimeError('CLAMP_EXCL_R=%r < 0' % (CLAMP_EXCL_R,))
CLAMP_DC_EFF = CLAMP_PULL / DISP if DISP > 0.0 else 0.0
#   DEAD_FRAC (5번째 인자) — 기저 상태의 하중비. 창이 좁다: 0.018 -> 무응력 퇴화 /
#   0.5 는 x_c=0.5(개시 61%)에서만 성공, 0.35·0.65 에서는 개시 위라 추출 실패. 스윕은 x_c 별로 낮춘다.
#   => f 는 [조건화될 만큼 크고 개시 아래일 만큼 작은] 창에서 고른다.
if _dead_arg is not None:
    DEAD_FRAC = float(_dead_arg)
CLAMP_PERT = PERTURBATION * CLAMP_DC_EFF   # 좌굴 스텝 클램프 섭동 [m]

if not (0.03 <= x_c <= 0.95):
    # 원장 M-8: x_c < 0.028 이면 클램프 패치가 정점 패치와 겹치고, x_c ~ 1 이면
    # 모서리 패치와 겹친다 -> 조용히 다른 모델이 된다. 경고만 찍고 진행한다.
    print("%s WARNING: x_c=%g is outside the safe band (0.03, 0.95) — "
          "the clamp patch may overlap a vertex patch (ledger M-8)." % (TAG, x_c))
if DISP < 0.0:
    raise RuntimeError('disp_m must be >= 0 (got %r).' % (DISP,))
# 상한 가드. DISP 는 이 모델의 크기(20 m)에 비해 작아야 한다.
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
elif CLAMP_DIR == 'along':
    # 접선(변을 따라 코너 쪽) — 클램프가 사선변을 **인장**한다. 두 클램프 모두 아래-바깥으로
    #   당기므로 꼭짓점 케이블 당김과 **같은 부호**다(현행 'normal' 은 위-바깥 = 반대 부호).
    DIR_CL, DIR_CR = (-1.0 / _s2, -1.0 / _s2), (+1.0 / _s2, -1.0 / _s2)
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
    # 예전에는 같은 조건을 두 번 검사해 두 번째 WARNING 이 **도달 불가**였다
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
    my_model.HomogeneousShellSection(name='Section-Membrane', material='Kapton', thickness=THICKNESS)   # 1차 셸 섹션(막 계열 기각)

    # 케이블 재질 (Kevlar) — mode.py 와 동일. 클램프와 무관하게 케이블 경로를 복원한다.
    mat_cable = my_model.Material(name='Kevlar')
    mat_cable.Elastic(table=((62.0e9, 0.36),))
    my_model.TrussSection(name='Section-Cable', material='Kevlar', area=CABLE_AREA)

    # 파트 생성: 멤브레인 — 동일
    if USE_GRID_MESH:
        # 격자: `.inp` 로 **orphan mesh part** 를 import 한다.
        #   Part.addNodes / addElements / deleteMesh 는 odb.Part 전용이라 mdb Part 에는 없다
        #   (실측 AttributeError 2026-10-06/07). CAE 스크립트로 격자를 넣는 경로는
        #   mdb.models[].PartFromInputFile 뿐이다.
        p, _nn, _n4, _n3 = aba_grid_mesh.import_grid_part(
            my_model, part_name='Membrane', base=BASE, height=HEIGHT, seed_div=SEED_DIV)
        print("%s [mesh] grid(.inp import, orphan): 노드 %d / S4 %d + S3 %d = %d"
              % (TAG, _nn, _n4, _n3, _n4 + _n3))
        #   orphan mesh: 섹션은 **요소 기반**으로, 요소가 있는 **뒤**에 준다(§10a).
        #   셋 이름 'All' 유지 — inst_memb.sets['All'] 참조가 남아 있다.
        #   orphan mesh 는 요소/절점 겸용 집합이 없다 -> 두 개를 만든다.
        p.Set(elements=p.elements, name='All_Elem')      # 섹션 할당 + 초기응력(요소)
        p.Set(nodes=p.nodes, name='All')                 # 변위 BC 용(절점) — 이름 'All' 유지
        p.SectionAssignment(region=p.sets['All_Elem'], sectionName='Section-Membrane')
    else:
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

    print("%s [mesh] %s : 파트 요소 %d개 / 인스턴스 노드 %d개"
          % (TAG, ('grid(격자, orphan)' if USE_GRID_MESH else 'free(자유, seedPart)'),
             len(p.elements), len(inst_memb.nodes)))
    if not USE_GRID_MESH:
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
        region=inst_memb.sets[_ALL_ELEM],
        distributionType=UNIFORM,
        sigma11=SIGMA0, sigma22=SIGMA0, sigma33=0.0,
        sigma12=0.0, sigma13=0.0, sigma23=0.0
    )

    # ---- 구속조건 ----
    # 시작 시 전체 면 z 고정 (평탄) — 동일
    #   [CLEAN_PERTURBATION] 이 면 이 BC 를 만들지 않는다: 좌굴 스텝에서 해제되는 순간
    #   스프링백 하중이 되어 섭동을 오염시킨다(위 상수 설명 참조).
    if not CLEAN_PERTURBATION:
        my_model.DisplacementBC(
            name='BC_Stabilize_Z',
            createStepName='Step-GlobalTension',   # [mode 정합] mode.py 와 동일
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
    #   corner_lf/clamp_live 도 여기서 구동 BC 를 만든다 — 새 케이스의 base(전량 프리텐션).
    if HAS_CABLES or CORNER_DIR == 'paper45' or CASE in ('corner_lf', 'clamp_live'):
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
    # 'cload' 라우트용 꼭짓점 하중은 아래 CLAMP_MODE 분기에서 일괄 생성한다 — 여기서 또 만들면
    #   같은 이름이 두 스텝에 등록되어 최종값이 0 이 된다(스텁 하네스로 적발).
    # [paper3] 꼭짓점도 케이블 축(+y)으로 당긴다(corner2 는 BC_Anchor_Top 고정 유지).
    if PRETENSION_MODE == 'paper3':
        my_model.boundaryConditions['BC_Anchor_Top'].setValuesInStep(
            stepName='Step-GlobalTension',
            u1=0.0, u2=disp_a * DISP_TOP_OVER_CORNER
        )
    _top_note = (' / 위 %.4e m' % (disp_a * DISP_TOP_OVER_CORNER)
                 if PRETENSION_MODE == 'paper3' else ' (Top 완전고정)')
    print("%s 프리텐션 = %s : 아래 %.4e m%s"
          % (TAG, PRETENSION_MODE, disp_a, _top_note))
    # ---- 하중/구속 배선: clamp_lf 는 꼭짓점·클램프를 prescribed 변위로 구동하고 DEAD/LIVE 로 나눈다 ----
    if CLAMP_LOAD == 'live':
        # =====================================================================
        # clamp_lf — 꼭짓점과 클램프를 prescribed 변위로 구동한다.
        #   이 저장소에서 **완주한** 구동 경로는 prescribed 변위 BC 뿐이다(HF/케이블/모드소스 모두
        #   Disp_Control_* + setValuesInStep). CLOAD 는 완주 기록이 없고 CAE 도 스텝 간 하중 수정을
        #   거부했다. 회계: DEAD=base state(GlobalTension) / LIVE=좌굴 스텝 -> lambda=1 이 운용점.
        # =====================================================================
        my_model.DisplacementBC(name='BC_Clamp_CL', createStepName='Initial',
                                region=a.sets['RP_CL_Set'], u3=0, **_CLAMP_ROT_KW)
        my_model.DisplacementBC(name='BC_Clamp_CR', createStepName='Initial',
                                region=a.sets['RP_CR_Set'], u3=0, **_CLAMP_ROT_KW)
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
            # LIVE = (1-DEAD_FRAC) x P — 남은 하중(운용점까지의 나머지)을 패턴으로 준다.
            #   [실측 확정] *BUCKLE 의 LIVE 는 '그 스텝에 지정한 값'이며 이전 스텝과의 차분이 아니다.
            #   증거: 스텝 값을 2배로 하니 lambda 가 정확히 절반(4모드 전부), 모드 형상은 불변.
            #   => lambda x (좌굴 스텝 값) = 좌굴까지의 추가 변위(불변량), 총 구동 = DEAD + lambda*LIVE.
            #   setValuesInStep 은 '그 스텝의 총값'을 대체하므로 lambda=1 이 운용점(설계 의도).
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
                                region=a.sets['RP_CL_Set'], u3=0, **_CLAMP_ROT_KW)
        my_model.DisplacementBC(name='BC_Clamp_CR', createStepName='Initial',
                                region=a.sets['RP_CR_Set'], u3=0, **_CLAMP_ROT_KW)

    if LIVE_DRIVER == 'clamp':
        # =====================================================================
        # clamp_live — base 는 **초기 인장만**(코너 전량, 클램프 당김 0 명시).
        #   좌굴 스텝의 패턴은 클램프 단독 — 합성 씨드(클램프 패밀리)용. lambda=프로브 배수.
        # =====================================================================
        my_model.DisplacementBC(name='Disp_Control_CL', createStepName='Initial',
                                region=a.sets['RP_CL_Set'], u1=0.0, u2=0.0)
        my_model.DisplacementBC(name='Disp_Control_CR', createStepName='Initial',
                                region=a.sets['RP_CR_Set'], u1=0.0, u2=0.0)
        for _nm, _d1, _d2 in (('Disp_Control_CL', DIR_CL[0], DIR_CL[1]),
                              ('Disp_Control_CR', DIR_CR[0], DIR_CR[1])):
            #   base: 클램프 당김 0 을 생략하지 않고 **명시**(초기 인장만). 좌굴 스텝: 클램프 단독.
            my_model.boundaryConditions[_nm].setValuesInStep(
                stepName='Step-GlobalTension', u1=0.0, u2=0.0)
            my_model.boundaryConditions[_nm].setValuesInStep(
                stepName='Step-Buckle',
                u1=PATTERN_SIGN * PERTURBATION * _d1,
                u2=PATTERN_SIGN * PERTURBATION * _d2)
        #   코너는 좌굴 스텝에서 0 명시 — LIVE = 클램프 단독 (E3 프로브 블록과 역할이 다르다).
        for _nm in ('Disp_Control_Right', 'Disp_Control_Left'):
            my_model.boundaryConditions[_nm].setValuesInStep(
                stepName='Step-Buckle', u1=0.0, u2=0.0)
        print("%s [clamp_live] base=코너 %.4g m + 클램프 당김 0 / LIVE=클램프 %.4g m x DIR (프로브)"
              % (TAG, DISP, PERTURBATION))

    # ---- Buckle 스텝: 논문 (c) 3변 u3=0 유지 + 클램프 구간은 면외 구속에서 제외(§9) ----
    _excl_edges = (aba_grid_mesh.clamp_exclude_labels(inst_memb, (V_CL, V_CR), CLAMP_EXCL_R)
                   if HAS_CLAMPS else set())
    _keep_labels = aba_grid_mesh.boundary_node_labels(
        inst_memb, base=BASE, height=HEIGHT, exclude_labels=_excl_edges)
    print("%s All_Edges_NoClamp 후보: 경계 노드 %d개 (클램프 반경 %.3g m 내 %d개 제외, 전체 막 노드 %d개)"
          % (TAG, len(_keep_labels), CLAMP_EXCL_R, len(_excl_edges), len(inst_memb.nodes)))
    aba_grid_mesh.make_set(a, 'All_Edges_NoClamp', inst_memb, _keep_labels)
    print("%s All_Edges_NoClamp 셋 생성 완료 (노드 %d개)" % (TAG, len(_keep_labels)))

    if CLEAN_PERTURBATION:
        #   자유경계 u3=0 을 **Initial** 에서 건다(패치 반경 내 노드는 이미 제외 — 이중구속 방지).
        #   => 스텝 간 구속집합이 동일해지고, LIVE 는 prescribed 변위 패턴 하나만 남는다.
        my_model.DisplacementBC(
            name='BC_Edges_Only_Z',
            createStepName='Initial',
            region=a.sets['All_Edges_NoClamp'],
            u3=0
        )
    else:
        my_model.boundaryConditions['BC_Stabilize_Z'].deactivate('Step-Buckle')
        my_model.DisplacementBC(
            name='BC_Edges_Only_Z',
            createStepName='Step-Buckle',
            region=a.sets['All_Edges_NoClamp'],
            u3=0
        )

    # ---- Buckle 스텝의 perturbation (케이블 라우트 전용) ----
    #   PATTERN_SIGN 을 여기에도 곱한다. mode.py 는 곱하는데 buckle 은
    #   빠져 있어 PATTERN_SIGN 이 1.0 을 벗어나면 두 스크립트가 조용히 갈라졌다(검사기는
    #   상수값만 비교하고 '사용 여부'는 보지 않는다). CLAMP_LOAD='live' 의 LIVE CLOAD 는
    #   위에서 이미 정의했으므로 여기서 건드리지 않는다.
    #   [함정] 라우트에 없는 BC 를 참조하면 KeyError 로 모델 생성이 죽는다(스텁 하네스가 적발).
    if HAS_CABLES or CASE == 'corner_lf':
        #   corner_lf: 좌굴 스텝 값 = 프로브(PERTURBATION). 모드 형상은 크기에 불변, lambda 만 스케일.
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
ELEM_TAG = {S4: 's4', S4R: 's4r', S3: 's3', S8R5: 's8r5', STRI65: 'stri65'}.get(ELEM_CODE_QUAD)
if ELEM_TAG is None:
    raise RuntimeError('ELEM_TAG 를 유도할 수 없는 요소코드입니다: %r. 잡 이름이 케이스를 '
                       '구분하지 못하면 run_job_safely 가 직전 증거를 지운다.' % (ELEM_CODE_QUAD,))
#   d_c_eff 도 이름에 넣는다 — CLAMP_PULL 을 CLI 로 바꿔 d_c 스윕을 하면 DISP 만으로는
#   이름이 겹쳐 run_job_safely 가 직전 증거를 지운다(§11).
_excl_tag = '' if CLAMP_EXCL_R == 0.0 else '_ex%03d' % int(round(CLAMP_EXCL_R * 1000.0))
#   새 케이스의 표기: base=전량 프리텐션(f100) / base 클램프 당김 없음(dc000). CASE 세그먼트가 달라 유일성은 보장.
_f_tag = int(round(DEAD_FRAC * 100.0)) if CASE == 'clamp_lf' else 100
_dc_tag = int(round(CLAMP_DC_EFF * 100.0)) if CASE == 'clamp_lf' else 0
#   clamp_live 는 프로브 크기가 이름에 들어간다 — 프로브를 바꾸면 산출물 이름도 갈린다.
_live_tag = ('_plv%05d' % int(round(PERTURBATION * 1.0e5))) if CASE == 'clamp_live' else ''
_rot_tag = '_rotfix' if CLAMP_RP_FIX_ROT else ''   # 스위치를 켜면 출력 이름이 갈린다
MESH_TAG = 'grid' if USE_GRID_MESH else 'free'   # 메쉬가 바뀌면 산출물 이름도 갈린다(§11)
JOB_NAME = 'Buckle_xc%03d_d%03dum_dc%03d_f%03d_%s_%s_%s_%s%s%s%s' % (int(round(x_c * 100.0)),
                                                       int(round(DISP * 1.0e6)),
                                                       _dc_tag,
                                                       _f_tag,
                                                       ELEM_TAG, MESH_TAG, CASE, CLAMP_DIR,
                                                       _excl_tag, _rot_tag, _live_tag)

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
