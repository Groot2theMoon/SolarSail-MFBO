"""
Solar Sail의 삼각형 멤브레인 모델을 생성하고, 입력 변수(Clamp 위치, 변위량)에 따라
구조 해석(Linear/Post-buckling)을 수행하는 Abaqus Python Script.

Usage:
    abaqus cae noGUI=run_abaqus.py -- [Fidelity] [x_c] [d_c]

Arguments:
    - Fidelity (str): 'LF' (Low-Fidelity: paper 경로는 초기응력 프리스트레스 정적해석,
                     구 기준선은 변위구동 정적해석) 또는 'HF' (High-Fidelity, Post-buckling)
    - x_c (float): Clamp Cable의 부착 위치 파라미터 (0.0 ~ 1.0)
                    0에 가까울수록 위쪽 꼭짓점, 1에 가까울수록 변의 아래 끝점에 위치.
    - d_c (float): paper 경로 = 클램프/꼭짓점 케이블 **초기응력 비율**
                    (CLAMP_PRESTRESS = d_c * CABLE_PRESTRESS).
                    구 disp 경로 = 클램프 케이블의 변위 비율.

Flow (Simulation Logic):
    1. 형상(Geometry) & 물성(Material) & 메쉬(Mesh) 생성 (공통)
    2. 초기 장력(Global Tension) 단계 정의
    3. Fidelity에 따른 분기:
       [LF 모드]  (클램프 있음 + 완전 형상)
         - paper 경로(PAPER_MATCH=True): 'Step-GlobalTension' 하나로 끝난다.
           케이블/클램프 초기응력 + 전막/3변 u3=0 인 **프리스트레스 상태**가 LF 최종 상태다.
           (ClampTension/HighTension 빈 스텝 없음. HF 의 Postbuckle base 와 같은 상태.)
         - 구 MFBO 기준선(PAPER_MATCH=False): 'Step-GlobalTension' -> 'Step-ClampTension'
           -> 'Step-HighTension' (1 mm 변위구동). 이 경로는 보존된다.
         - 임퍼펙션(*IMPERFECTION) 없음 -> 좌굴(Buckle) 잡도 돌리지 않는다(모드 불필요).
         - eval 은 LF ODB 에서 lf1/lf2/lf3 + ZK(Zhang & Kiendl) 지표를 낸다.
       [HF 모드]  (클램프 있음 + 임퍼펙션 주입 = 2-모델 레시피, 2026-09-22)
         - 1차: 모드 소스는 '클램프 없는' 별도 모델에서 **선형 좌굴해석(*BUCKLE)** 으로 모드를
                뽑는다(run_abaqus_mode.py). 클램프가 있는 base state 는 선형 좌굴 스펙트럼을
                주지 못한다(실측: 음수고유값 598~2897 / CONVERGED=0 / EIGENVALUES CANNOT BE FOUND).
                *BUCKLE 은 .fil 출력이 금지되므로(실측 7025) 모드는 ODB 모드 프레임에만 있다
                -> aba_mode_from_odb.py 가 그 ODB 를 모드표(modes_ClampFree_Buckle.txt)로 만든다.
         - 2차: 그 모드를 노드 좌표 섭동으로 주입한다(IMPERFECTION_MODE='odb_table' = 기본,
                seeds/modes_xc<tag>.txt 모드표를 읽는다. 'odb_direct' = ODB 직접 읽기).
                aba_imperfection.py 가 모드 개수를 검증하고, 부족하면 HF 제출 전에 중단한다.
                자기 좌굴 잡(Buckle_Analysis)은 클램프 base state 증거용으로만 유지한다.
                (옵션 IMPERFECTION_MODE='file' = .fil 스테이징 + *IMPERFECTION, FILE=, STEP=n:
                 사용자 원본 08d4cbc 방식. 단 좌굴모드에는 쓸 수 없다 — *BUCKLE 은 .fil 금지.)
         - 3차: 'Step-Postbuckle' 수행 (Riks/Stabilization). 주름 거동을 포함한 비선형 해석.
    4. Post-Processing:
       - 'eval_abaqus.py'를 subprocess로 호출하여 ODB에서 필요한 값(추력, 면적 등)만 추출.
       - 추출된 값을 stdout으로 출력하여 상위 프로세스(MFBO)에 전달.
       - LF/HF(짝 LF ODB 가 있을 때) 공통: `ZK:` 한 줄(별도 계약) + zk_metrics.txt.
         RESULTS 3/4필드는 불변.

solar-sail 해석 방법론은 대부분 Galhofo2022 논문을 레퍼런스로 하였음.
(물성치, 크기, 각도, step의 종류와 순서 등등)
 - displacementBC 값들은 경험적으로 찾은 값들. 더 좋은 숫자가 있을 수 있음.
"""

from abaqus import *
from abaqusConstants import *
from step import *
from mesh import ElemType
import re
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

# ---- explicit CLI config (no environment variables) ----
for _cp in (os.getcwd(), os.path.dirname(os.getcwd()),
            os.path.dirname(os.path.dirname(os.getcwd()))):
    if os.path.exists(os.path.join(_cp, 'analysis_config.py')):
        if _cp not in sys.path:
            sys.path.insert(0, _cp)
        break
import analysis_config as _AC
_CLI_REST, _CLI_CONFIG, _CLI_CASE = _AC.split_cli(sys.argv[1:])
_CFG = _AC.load_config(_CLI_CONFIG)
_CFG_HASH = _AC.config_hash(_CFG)

# ===== SETTINGS =====
#   [§19.110] 실행/모델 설정값은 전부 이 블록에 있다. 값을 고칠 때는 여기만 본다.
#     종전에는 설정이 2259줄 본문에 흩어져 있었다 — NUMCPUS 97 / N_EIG 836 / MODE_SOURCE 840 /
#     STABILIZE_LADDER 862 / NO_CLAMP(정의 1028 + 인라인 재대입 1421) / PULL_TOP 1032 /
#     TOP_PULL_SCALE 1039 / SIGMA0 1429 / DISP_GLOBAL 1205 / GLOBAL_FINAL 1210 /
#     _EIGENSOLVER 1523 / BUCKLE_ODB 1703. 값/주석/**그룹 내 정의 순서는 그대로** 두고
#     위치만 모았다(값 변경 0건 -> 생성 덱과 콘솔 로그가 종전과 바이트 동일해야 한다).
#   이 아래의 파생부는 이 값을 읽어 계산한다 (손대지 않는다):
#     _HERE/_RUN(경로) -> CLI 파싱(fidelity/x_c/d_c) -> V_CL/V_CR
#     -> CLAMP_PRESTRESS = d_c * CABLE_PRESTRESS -> MODE_TABLE/IMPERFECTION_MODES(모드표 게이트)
#     -> LF_ODB/LF_JOB -> _ALL_ELEM -> _steps_in_order/_BUCKLE_STEP_NO -> BUCKLE_MODEL
#   그룹: ① 실행 모드/스위치 ② 형상/격자 ③ 재료 ④ 하중/프리텐션 ⑤ 좌굴/임퍼펙션
#         ⑥ HF 솔버 ⑦ 경로/식별자
#   ⚠️ 런타임 인자(Fidelity, x_c, d_c)는 상수가 아니다 — 아래 CLI 파싱부가 sys.argv[-3:] 로 읽는다.
# ============================================================================

# ---- ⓪ 논문 정합 프리셋 (Galhofo 2022, 5PCS) ----
#   True  = 논문 조건: 3케이블 등응력 프리텐션(초기응력), 3변 u3=0(Initial부터),
#           포스트버클 안정화 2e-4, 2차 요소(S8R5+STRI65), 지정변위 0(변위구동 끔).
#   False = 종전 MFBO 기준선: 지정변위 구동, 1차 S4/S3, 안정화 사다리(0.005~0.05).
#   !! PAPER_MATCH=True 는 **시드표 재생성**이 필요하다:
#        make_seed.py 로 modes_seed.txt 를 2차 격자에서 다시 만든 뒤 HF 를 돌린다.
#      표가 1차 라벨이면 seed_family_gate 가 HF 제출 전에 중단한다.
PAPER_MATCH = True
PAPER_MESH_ORDER = 1             # 기본 1차 S4/S3. 2차(S8R5/STRI65)는 실험용 — 모드 추출 base state 불안정
MESH_ORDER = PAPER_MESH_ORDER if PAPER_MATCH else 1
LOAD_MODE = 'cable_stress' if PAPER_MATCH else 'disp'   # 'disp' | 'cable_stress'
CABLE_PRESTRESS = 2.69e7         # Pa — 논문 5PCS Fig. 11, 3케이블 동일
PAPER_STABILIZE = 2.0e-4         # 논문 stability factor
PAPER_EDGE_BC_FROM_INITIAL = bool(PAPER_MATCH)   # 3변 u3=0 을 Initial 에서 (HF 전용)
#   [2026-10-20] 논문 재현과 실제 설계를 분리한다.
#     False = 논문 5PCS 그대로 (클램프 없음, NO_CLAMP=True)
#     True  = 논문 프리텐션·BC·1차 메쉬는 유지하되 **클램프 + 꼭짓점 케이블을 모두 포함**
#             (실제 포스트버클링 모델. 클램프 프리텐션 = 초기응력 d_c*CABLE_PRESTRESS.)
PAPER_INCLUDE_CLAMPS = True
#   [2026-10-21] u3 정책은 **fidelity 종속**이다:
#     LF  = Initial 부터 전막 u3=0 (평평한 프리스트레스)
#     HF  = paper 경로는 3변만 u3=0 (처음부터 자유; 주름은 GlobalTension 에서 발생)
#           old disp 경로는 기존처럼 GlobalTension hold -> Postbuckle release
#   PAPER_HOLD_U3_GLOBAL 상수는 제거했다. _HOLD_U3_ALL 이 런타임 fidelity 로 분기한다.
#   [paper-match] 모드표 = **seed + clamp_only 합성표** (사용자 결정 2026-10-20).
#     True  = seeds/modes_xc<xc>.txt (2+2 합성) 를 쓰고 seed_gate 지문 검사를 적용한다.
#     False = 종전 NO_CLAMP 순수 seed 가족 표(modes_seed.txt) — clamp 모드 배제.
#     합성표를 쓰면 clamp 가족 모드가 임퍼펙션에 섞인다(설계 의도). 모드 추출은 make_seed.py.
USE_COMPOSITE_MODE_TABLE = True

# ---- ① 실행 모드/스위치 ----
BASE_PROBE = True                # R-13/CT 라이브. MFBO에서는 False로 끌 수 있다.
EIG_RECORD = True                # 좌굴 고유치 이력 기록 (coalescence_check.py record)
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
#   [§19.104] NO_CLAMP 상수를 이 블록으로 올렸다(값/의미 불변 - 종전에는 클램프 RP 조립부에 있었고,
#     모드표 결정보다 **뒤**라 여기서 볼 수 없었다). True = 클램프 생략 baseline(run_abaqus_cable 대조).
NO_CLAMP = not PAPER_INCLUDE_CLAMPS
#   [paper-match] 논문 5PCS 재현은 클램프가 없다. 실제 설계는 PAPER_INCLUDE_CLAMPS=True 로 켠다.
if PAPER_MATCH and not NO_CLAMP and not PAPER_INCLUDE_CLAMPS:
    raise RuntimeError('PAPER_MATCH=True 의 논문 재현은 NO_CLAMP=True 를 요구합니다. '
                       '실제 클램프 모델은 PAPER_INCLUDE_CLAMPS=True 로 켜십시오.')
#   [§19.105] PULL_TOP: True = 윗 꼭짓점(Top)도 케이블 방향 지정변위로 당긴다(Galhofo 2022 3-tendon 정합).
#     False(기본) = 종전 그대로 Top 완전 고정 - 기존 런 재현성이 우선이다(BC/스텝 호출 바이트 동일).
#     NO_CLAMP 와 **독립**이다: 네 조합(PULL_TOP x NO_CLAMP) 모두 유효하고, 클램프 유무와 서로 간섭하지 않는다.
PULL_TOP = False
#   [§19.109] TOP_PULL_SCALE: 세 당김의 힘 평형 보정 계수 - **Top 지정변위에만** 곱한다.
#     케이블 힘 F = E*A*u/L 이므로 대칭 3-tendon 의 평형은 F_top = 2*sin(theta)*F_bot,
#     즉 u_top = 2*sin(theta)*(LEN_TOP/LEN_BOT)*u_bot (theta = angle_deg = 28.6 deg).
#     1.0 (기본) = 종전과 동일(보정 없음, 덱/로그 바이트 동일). 0.389 = 힘 평형 조건.
#     유도: 2*sin(28.6 deg)*0.280/0.689 = 0.389 (길이를 무시한 0.9574 는 틀린 값이다).
#     Right/Left/클램프 변위에는 곱하지 않는다(세 당김 중 Top 만 평형 보정 대상).
TOP_PULL_SCALE = 1.0
# 실행 모드.
#   'auto'  : 덱 생성 -> riks_patch_input.py 로 *Controls 주입 -> 제출 -> 결과 추출까지 한 번에.
#             mfbo.py 가 이 경로로 HF 를 부른다(덱 생성만 하면 RESULTS: 가 없어 실패로 읽힌다).
#   'write' : 덱만 생성하고 종료. 수동으로 패처/job 을 돌리는 워크플로우(latex/디버깅용).
#   'plain' : 패처 없이 원본 덱을 그대로 제출. 패처 자체를 의심할 때만 쓴다(진단용, 40분 낭비 주의).
RUN_MODE = 'auto'

# ---- ② 형상/격자 ----
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

#   [§19.111] RELEASE_EDGE_U3: 자유변의 면외 구속을 '하드 u3=0' 에서 '무른 to-ground 스프링' 으로 바꾼다.
#   True = 자유변 BC_Edges_Only_Z 를 제거하고 '무른 면외 스프링'으로 대체 (경계 아티팩트 검증용).
#   False = 종전 그대로(u3=0). 덱/로그 바이트 동일 유지.
RELEASE_EDGE_U3 = False
#   스프링 강성 = ALPHA * (σt)_ref,  σt = 7 kPa * 5 µm = 0.035 N/m (막 장력/단위폭 = 스프링 [N/m])
EDGE_SPRING_ALPHA = 1.0e-3     # 1e-3(거의 자유) ~ 1(자연 강성) 스윕용
EDGE_SPRING_K = EDGE_SPRING_ALPHA * 0.035

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
# ---- 메쉬 종류 (2026-10-07) ----
#   True  = aba_grid_mesh 균일 격자(.inp -> PartFromInputFile, orphan mesh). 좌우 완전 대칭.
#   False = 기존 자유 메쉬(seedPart + QUAD_DOMINATED/FREE/MEDIAL_AXIS, 약 1.82만 요소).
#   !! 네 스크립트(buckle/mode/HF/cable)가 **같은 값**을 써야 한다.
USE_GRID_MESH = True

# ---- ③ 재료 ----
#   승격한 것: 케이블(Kevlar) 물성 3개 + 클램프 패치 반경. 값은 종전 리터럴과 **동일**하다.
#   승격하지 않은 것(검사기 계약): 막(Kapton) `1420.0` / `2.5e9` / `0.34` 와 꼭짓점 패치
#     `radius=0.2` 는 check_model_consistency.py 의 SHARED 항목이
#     `mat.Density(table=((1420.0,),))` / `mat.Elastic(table=((2.5e9, 0.34),))` / `radius=0.2`
#     같은 **리터럴 문자열**로 buckle 스크립트와 대조한다 -> 이름을 붙이면 그 게이트가 FAIL 한다.
#     (값 이중화를 만들지 않으려고 리터럴을 그대로 뒀다. notes §19.110 참조.)
CABLE_RHO = 1440.0                 # kg/m^3 — Kevlar 밀도 (종전 mat_cable.Density 리터럴)
CABLE_E = 62.0e9                   # Pa — Kevlar 영률 (종전 mat_cable.Elastic 리터럴)
CABLE_NU = 0.36                    # - — Kevlar 푸아송비 (종전 mat_cable.Elastic 리터럴)
# 초기 가짜 응력(수렴 보조). 케이블 변형=700 Pa, 우리=500 Pa -> 정렬 노브
SIGMA0 = 500.0                   # 수렴 보조용 초기응력 [Pa]
PATCH_R_CLAMP = 0.21               # m — 클램프 강체패치 반경 (종전 radius=0.21 리터럴 2곳)

# ---- ④ 하중/프리텐션 ----
sqrt2 = 1.414
# ---- 사전 장력(prestrain) 캘리브레이션 ----
# 사전 장력 크기. base state 장력을 키워 시스템행렬 부정정(음수 고유값) 완화를 시도한다.
# 1.0 = 기존값. 값을 바꾸려면 이 상수를 직접 수정한다(구 MFBO_PRETENSION_SCALE 환경변수는 제거됨).
PRETENSION_SCALE = 10.0          # 프리텐션 변위 = 5e-6 m * 이 값 = 5e-5 m
# PRETENSION_SCALE=10 -> 면내응력 2122 Pa = 목표 7000 Pa 의 0.303배. 운용점 조정은 DISP_GLOBAL 로 한다(§8).
DISP_GLOBAL = 165e-6
# 좌굴 스텝 perturbation — K_delta 를 만든다. 너무 작으면 고유값 분리가 나빠져 추출이 실패한다(§8b).
PERTURBATION = 0.01
GLOBAL_FINAL = 1e-3  # 최종 하중: 코너 당김 1e-3 m
# 하중 각도 (28.6도)
angle_deg = 28.6
angle_rad = np.deg2rad(angle_deg)
cos_val = float(np.cos(angle_rad))
sin_val = float(np.sin(angle_rad))
#   [§19.105] Top 케이블 방향 코사인. 케이블은 connect_cable('Cable_Top', ..., V1, (0.0, 1.0, 0.0)) 로
#     +y 축에 놓인다 -> 방향 (u1, u2) = (0, +1). Right/Left 와 **같은 대입식**을 쓰기 위해
#     (u1 = v*cos, u2 = -v*sin) 규약에 맞춘 값이 top_cos=0.0, top_sin=-1.0 이다(즉 u2 = +v).
top_cos = 0.0
top_sin = -1.0

# ---- ⑤ 좌굴/임퍼펙션 ----
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
IMPERFECTION_AMPL_T = 0.1    # 막 두께 배수. 현재 0.10 t (= 5.0e-7 m). 과거 1.0 t 실험은 수렴/증분 문제로 폐기.
#   [§19.81] LF 는 **무결함** 모델로 계산한다(사용자 결정).
#     LF 는 MFBO 의 "싼 대리 모델"이다 — 그 정의가 시드표/진폭에 의존하면 같은 설계점의 LF 값이
#     시드에 따라 흔들려 다중충실도 구조가 약해진다. LF 분기에서는 섭동 크기를 0 으로 만들어
#     메쉬를 완전하게 둔다(주름은 HF 에서만 발달한다는 해석과도 일치).
#   [2026-10-21] u3-free HF: GlobalTension 에서 주름이 생기므로 사다리는 **내리는 방향**이다.
#     physical LF 는 기존 GlobalTension 안정화(1.0)를 유지하고, HF 만 2e-4 에서 시작한다.
#     Postbuckle 은 정착 스텝(paper) 또는 최종 하중(disp)이며 5e-4 에서 시작한다.
#     (구 0.005-0.05 단일 사다리는 Postbuckle 전용 실험이었다.)
_PHYSICAL_PAPER = bool(PAPER_MATCH and PAPER_INCLUDE_CLAMPS)
STABILIZE_GLOBAL_LADDER = ((2.0e-4, 1.0e-3, 1.0e-2) if _PHYSICAL_PAPER
                           else ((PAPER_STABILIZE,) if PAPER_MATCH else (0.0002,)))
STABILIZE_POST_LADDER = ((5.0e-4, 2.0e-3, 1.0e-2) if _PHYSICAL_PAPER
                         else ((PAPER_STABILIZE,) if PAPER_MATCH else (0.01,)))
#   LF 는 u3=0 고정이라 기존 GlobalTension 감쇠를 그대로 쓴다(회귀 기준선 보존).
STABILIZE_GLOBAL_LF = (1.0 if _PHYSICAL_PAPER else 0.0002)
STABILIZE_GLOBAL_HF = STABILIZE_GLOBAL_LADDER[0]
STABILIZE_GLOBAL = STABILIZE_GLOBAL_LF
STABILIZE_POST = STABILIZE_POST_LADDER[0]
#   [2026-10-21] u3-free HF 에서는 임퍼펙션 진폭이 **주름 패턴 자체를 선택**한다.
#     기본 0.10 t. 0.1 vs 0.5 민감도를 반드시 비교한다(PAPER_MATCH_BRANCH.md).
#     과거 1.0 t 는 수렴/증분 문제로 폐기했고, 0.10 t 는 평평한 base 의 tie-breaker 였다.
IMPERFECTION_MODE = 'odb_table'   # 'odb_direct'(기본) | 'odb_table' | 'file'
#   [2026-10-12] 설계점별 시드표 - 경로 결정과 시드 게이트는 아래 `_available_mode_numbers` def
#     뒤 "[§19.104] 임퍼펙션 표 선택" 블록으로 **옮겼다**. 이유: 클램프 없는 baseline(NO_CLAMP=True)
#     은 순수 seed 가족 표를 써야 하는데 표 선택 시점에 NO_CLAMP 상수가 필요하고, 그 상수는 종전
#     모델 조립부(클램프 RP 직전)에 있어 여기서 볼 수 없었다(값/의미 불변, 위치만 이동).
MODE_SOURCE_ODB = os.path.join('..', 'ClampFree_Buckle.odb')    # 'odb_direct' 가 읽는 ODB
MODE_SOURCE_STEP_NAME = 'Step-Buckle'                           # 그 ODB 안의 좌굴 스텝 이름
MODE_INSTANCE = 'MEMBRANE-1'                                    # 막 인스턴스 이름
SEED_FAMILY_NAME = 'modes_seed.txt'   # seed 가족 표(꼭짓점 케이블만, x_c 무관, 모드 4개)
# A: 추출 요청 고유값 수. base state 가 부정정이면 줄이는 것이 subspace 수렴에 유리하다(§8a).
N_EIG_BUCKLE = 100
# subspace 반복의 기저 벡터 수 (구 MFBO_VECTORS 환경변수는 2026-09-21 제거됨)
BUCKLE_VECTORS = 250

# ---- ⑥ HF 솔버 ----
# Abaqus 병렬 스레드/도메인 수 (1 = 단일 CPU, 4 = 대략 1.5~2.5배 빠름).
#   라이선스 토큰이 없으면 잡이 라이선스 오류로 즉시 죽는다 -> 그때는 1 로.
NUMCPUS = 4
# A: 좌굴 고유값 추출 솔버. **기본 LANCZOS** — Abaqus 의 *BUCKLE 기본 솔버이자 이 코드의
#    원래 설정(LANCZOS + maxBlocks=DEFAULT)이다. 요소 ~1만 개 / DOF ~4만 개 + 두께 5um 라
#    강성 조건수가 극단적이어서 SUBSPACE(작은 모델·다수 모드용)로는 100개 모드 추출이
#    'EIGENVALUES CANNOT BE FOUND' (0 CONVERGED) 로 실패했다.
#    SUBSPACE 로 강제하려면 _EIGENSOLVER 상수를 'SUBSPACE' 로 수정 (구 MFBO_EIGENSOLVER 환경변수는 제거됨)
_EIGENSOLVER = "LANCZOS"         # "SUBSPACE" 로 바꾸면 좌굴 추출 재시도

# ---- ⑦ 경로/식별자 ----
#   SEEDS_DIR / SEED_FAMILY_TABLE 은 여기 없다 — _HERE 와 aba_imperfection.SEED_DIR_NAME 의
#   파생이라 import 뒤(모드표 선택 블록)에서 계산한다. LF_ODB/LF_JOB 도 (x_c,d_c) 파생이다.
# 산출물 디렉터리 = code/aba. mfbo.py 의 RUN_DIR_NAME 과 같은 값이어야 한다.
RUN_DIR_NAME = "aba"
BUCKLE_ODB = 'Buckle_Analysis.odb'   # legacy self-buckle name; R-13 does NOT use it (point tag 없음)
HF_JOB = 'HF_Postbuckle_u3free'
HF_ODB = HF_JOB + '.odb'
# ============================================================================
#   Explicit config overrides the default settings above.  No env vars.
#   [CONFIG] 다른 프로세스가 값에 의존하는 키는 config 로 바꿀 수 없다.
#     RUN_DIR_NAME: mfbo.py 가 "aba" 를 하드코딩한다 -> 바꾸면 mfbo 가 ODB 를 못 찾는다.
_CFG_FORBIDDEN = ('RUN_DIR_NAME',)
_CFG_APPLY = _AC.apply_globals(globals(), _CFG, forbidden=_CFG_FORBIDDEN)
#   config 신원(sha1). LF sidecar 와 점별 results txt 에 남겨 "이 점이 어느 config 였나"를 추적한다.
_CFG_SHA1 = _AC.config_hash(_CFG)

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
#   RUN_DIR_NAME/NUMCPUS/BASE_PROBE/EIG_RECORD 는 SETTINGS ⑦/⑥/① 에 있다.
_RUN = os.path.join(_HERE, RUN_DIR_NAME)
os.makedirs(_RUN, exist_ok=True)
os.chdir(_RUN)

# ---- 임퍼펙션 소스 스테이징 모듈 (2-모델 레시피, stdlib only) ----
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from aba_imperfection import (ImperfectionSourceError, stage,           # noqa: E402
                              imperfection_text, report, load_mode_table,
                              build_perturbation, perturbation_report, mode_table_info,
                              verify_inp_element_types,
                              EDGE_SPRING_ELEM_TYPES,
                              point_tag, lf_odb_name, lf_job_name,
                              lf_sidecar_write, lf_odb_is_current,
                              lf_sidecar_read, model_fingerprint,
                              seed_table_path, seed_gate, BUCKLE_SCRIPT_NAME,
                              SEED_DIR_NAME, job_completed_from_logs, value_tag)
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


def _probe_keep_line(_ls):
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
    _EMIT_PATHS = _emit_paths()
    for p in _EMIT_PATHS:
        try:
            with open(p, 'w') as f:
                f.write('[canary] run_abaqus.py start %s  cwd=%s  argv=%s\n'
                        % (time.strftime('%Y-%m-%d %H:%M:%S'), os.getcwd(), list(sys.argv)))
        except Exception:
            pass


# [§19.94] 런 시작 시각. R-13 이 읽은 ODB 의 mtime 이 이보다 오래되었으면 '이 런이 만든 파일이
#   아니다' -> stale 가능 경고를 찍는다(그 숫자를 이 설계점의 base state 로 읽으면 안 된다).
_RUN_T0 = time.time()
# [§19.97] R-13 Step-ClampTension 판독 결과(점별 results_*.txt 아카이빙용). BASE_PROBE 블록이 채운다.
#   비어 있으면 점별 txt 에는 'nan' + 이유가 들어간다(_r13_ct_keys 참조).
_R13_CT = {}
emit_start()
_AC.log_config(emit, _CFG, _CLI_CASE)
#   [CONFIG] 적용/무시된 키를 남긴다. 오타 키가 조용히 사라지면 "설정했는데 안 먹었다"가 된다.
#     unknown 은 이 스크립트 globals 에 없는 키다(다른 스크립트용 키일 수 있으므로 경고만).
if _CFG:
    emit('[CONFIG] applied=%d keys=%s'
         % (len(_CFG_APPLY.get('applied') or []),
            ','.join(_CFG_APPLY.get('applied') or []) or '-'))
if _CFG_APPLY.get('unknown'):
    emit('[CONFIG] !!! unknown(이 스크립트에 없는 키 -> 무시): %s'
         % ', '.join(sorted(_CFG_APPLY['unknown'])))


try:
    # abaqus cae noGUI=run_abaqus.py -- [HF/LF] x_c d_c [--config ...]
    fidelity = _CLI_REST[-3].upper()
    x_c = float(_CLI_REST[-2])
    d_c = float(_CLI_REST[-1])
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

#   [paper-match] 논문 하중/BC 경로는 **LF/HF 공통**이다 (2026-10-21).
#     `_PAPER_PATH` = 케이블/클램프 초기응력 + 3변/패치 u3=0 + 변위구동 끔.
#     `fidelity` 는 해석 종류(LF=정적 프리스트레스, HF=포스트버클)만 제어한다.
#     `_PAPER_HF` 는 HF 전용 분기(예: Postbuckle 에서 BC 재생성 회피)에만 남긴다.
_PAPER_PATH = bool(PAPER_EDGE_BC_FROM_INITIAL)
_PAPER_HF = _PAPER_PATH and (fidelity == 'HF')
#   [2026-10-23] paper HF의 no-load Step-Postbuckle은 제거한다.
#     실측: GlobalTension은 완주하는데 Postbuckle은 TOO MANY ATTEMPTS로 실패(no load,
#     allsdtol=0.002, line search 0회). 주름은 이미 GlobalTension에서 형성되므로 hf/DN은
#     GlobalTension 마지막 프레임에서 읽는다. old disp HF는 BC 해제/최종 변위 때문에 유지.
_HF_USE_POSTBUCKLE = (LOAD_MODE == 'disp') or (not _PAPER_HF)
#   [2026-10-21] u3 정책을 fidelity 로 분기한다.
#     LF  = 항상 전막 u3=0 (깨끗한 프리스트레스 base)
#     HF  = paper 경로는 처음부터 자유(3변만 u3=0). old disp 는 기존 hold 유지.
_LF_HOLD_U3_ALL = (LOAD_MODE == 'disp') or bool(_PAPER_PATH)
_HOLD_U3_ALL = (LOAD_MODE == 'disp') or (_PAPER_PATH and (fidelity == 'LF'))
_DISP_PRE = 0.0 if (LOAD_MODE == 'cable_stress' and _PAPER_PATH) else DISP_GLOBAL
_DISP_FINAL = 0.0 if (LOAD_MODE == 'cable_stress' and _PAPER_PATH) else GLOBAL_FINAL
#   [paper-match] LF ODB stale guard: 모델 정의가 바뀌면 기존 LF ODB 를 재사용하지 않는다.
#     sidecar(.model.txt)의 tag 를 현재 경로와 대조한다(aba_imperfection 참조).
#     tag 에 clamp/u3 sub-mode 와 field-output 버전(e2: ZK 전용 F-ZK E 요청)을 포함한다.
#     기본값은 aba_imperfection.LF_MODEL_TAG 와 같아야 한다(mfbo.py 가 그 값을 쓴다).
#     v3-e2 로 올린 이유: v2-e1 ODB 에도 E 가 없어 zk_Iw/ZK 가 FAIL 할 수 있다.
_LF_MODEL_TAG = (('globaltension-v3-e2-c%d-u%d'
                  % (int(PAPER_INCLUDE_CLAMPS), int(_LF_HOLD_U3_ALL)))
                 if (LOAD_MODE == 'cable_stress' and _PAPER_PATH)
                 else 'disp-v3-e2')

#   [CONFIG] 모델 정의값 지문(model_fp). 하드코딩 tag 만으로는 config 로 바뀐
#     THICKNESS/SEED_DIV/MESH_ORDER/프리텐션/재료가 재사용 가드에 안 잡혔다
#     -> 같은 (x_c,d_c) 의 옛 LF ODB 가 조용히 재사용되어 ZK/lf 가 오염된다(실측 위험:
#     configs/paper_sd632.json 은 SEED_DIV=632 로 요소 ~10배인데 태그/파일명이 같다).
#     이 지문을 sidecar 에 남기고 재사용 전에 대조한다(모르는 호출자는 config sha 로 대조).
_LF_MODEL_FP = model_fingerprint({
    'SEED_DIV': globals().get('SEED_DIV'),
    'MESH_ORDER': globals().get('MESH_ORDER'),
    'THICKNESS': globals().get('THICKNESS'),
    'SIGMA0': globals().get('SIGMA0'),
    'CABLE_PRESTRESS': globals().get('CABLE_PRESTRESS'),
    'CABLE_AREA': globals().get('CABLE_AREA'),
    'BASE': globals().get('BASE'),
    'HEIGHT': globals().get('HEIGHT'),
    'ELEM_CODE_QUAD': globals().get('ELEM_CODE_QUAD'),
    'ELEM_CODE_TRI': globals().get('ELEM_CODE_TRI'),
    'USE_GRID_MESH': globals().get('USE_GRID_MESH'),
    'PAPER_INCLUDE_CLAMPS': globals().get('PAPER_INCLUDE_CLAMPS'),
    'PAPER_MATCH': globals().get('PAPER_MATCH'),
    'LOAD_MODE': globals().get('LOAD_MODE'),
    'LF_HOLD_U3_ALL': _LF_HOLD_U3_ALL,
})
emit('[CONFIG] lf_model_fp=%s config_sha1=%s lf_tag=%s'
     % (_LF_MODEL_FP, _CFG_SHA1, _LF_MODEL_TAG))

# P1-3: 잡 제출 전에 지워야 하는 이전 실행 산출물
_JOB_ARTIFACTS = ('odb', 'fil', 'sta', 'msg', 'lck', 'com', 'prt', 'sim', 'log',
                  'dat', 'res', 'abq', 'ipm', 'mdl', 'stt', 'cid')

def _archive_result(x_c, d_c, fidelity, results_csv, extra=None):
    """[§19.88] 설계점별 결과 아카이빙.

    왜: HF 잡 이름이 `HF_Postbuckle` 고정이라 다음 런이 **직전 산출물을 지운다**
    (`run_job_safely` 가 odb/dat/msg 를 제거). 그래서 (0.5,0.3)의 ODB 를 이미 잃었다 —
    결과를 설계점 이름의 텍스트 + 누적 CSV 로 남겨 스윕 분석이 한 파일로 가능하게 한다.

    남기는 것: RESULTS 4필드(lf1,lf2,lf3,hf) + 사다리가 쓴 stabilize_used/에너지 비/판정 + 시각.
    """
    try:
        import time
        _dir = os.path.join(_HERE, 'results')
        if not os.path.isdir(_dir):
            os.makedirs(_dir)
        _p = [str(v).strip() for v in (results_csv or '').split(',')]
        while len(_p) < 4:
            _p.append('nan')
        _kv = {'x_c': x_c, 'd_c': d_c, 'fidelity': fidelity,
               'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
               'lf1': _p[0], 'lf2': _p[1], 'lf3': _p[2], 'hf': _p[3]}
        _kv.update(extra or {})
        #   [CONFIG] 이 점이 **어느 config/모델**로 계산됐는지 남긴다. HF 는 짝 LF ODB 를
        #     재사용할 수 있으므로 "이 런의 config"(_CFG_SHA1)와 "LF ODB sidecar 의 config"를
        #     따로 적는다(둘이 다르면 재사용이 섞였다는 신호다). CSV 열은 늘리지 않는다.
        _kv.setdefault('config_sha1', _CFG_SHA1)
        try:
            _lfsd = lf_sidecar_read(globals().get('LF_ODB') or '')
        except Exception:
            _lfsd = {}
        _kv.setdefault('lf_config_sha1', _lfsd.get('config_sha1', 'nan'))
        _kv.setdefault('lf_model_fp', _lfsd.get('model_fp', 'nan'))
        _tag = 'results_xc%s_dc%s.txt' % (value_tag(x_c), value_tag(d_c))
        with io.open(os.path.join(_dir, _tag), 'w', encoding='utf-8') as _f:
            for _k in sorted(_kv):
                _f.write('%s=%s\n' % (_k, _kv[_k]))
        #   [§19.91] 열 스키마를 **고정**한다 — 진단키(stabilize/에너지/flip_bias)도 CSV 에 실제로
        #     실려야 스윕 분석이 한 파일로 된다(헤더가 1회만 쓰이므로 열 목록은 상수여야 한다).
        _cols = ['x_c', 'd_c', 'fidelity', 'timestamp', 'stabilize_used',
                 'allds_allie_last_pct', 'allds_allie_max_pct', 'energy_verdict',
                 'backward_area_pct', 'backward_elems', 'total_elems',
                 'hf_flip', 'hf_true_normal', 'delta_pct',
                 'lf1', 'lf2', 'lf3', 'hf']
        _csv = os.path.join(_dir, 'all_results.csv')
        _new = not os.path.exists(_csv)
        with io.open(_csv, 'a', encoding='utf-8') as _f:
            if _new:
                _f.write(','.join(_cols) + '\n')
            _f.write(','.join(str(_kv.get(_c, '')) for _c in _cols) + '\n')
        emit("[ARCHIVE] %s + all_results.csv (1행 추가)" % _tag)
    except Exception as _e:
        emit("[ARCHIVE] 실패(무시): %s" % _e)


def _archive_slope(x_c, d_c):
    """[§19.102] **성공한 HF 런에서만** 호출: (1) HF ODB 설계점 태깅 복사,
    (2) slope_stats.compute_stats 를 인프로세스로 돌려 점별 results txt 에 slope_ 접두로 아카이빙.

    모두 예외 안전 - 실패해도 런을 죽이지 않고 경고 1줄만 남긴다. RESULTS 4필드 계약과
    all_results.csv 열 스키마는 건드리지 않는다(점별 txt append + 별도 JSON 만).
    스위치: SOLARSAIL_NO_SLOPE=1 -> 전체 생략, SOLARSAIL_NO_ODB_COPY=1 -> ODB 복사만 생략.
    로그 접두 [SLOPE]/[ODBCOPY] 는 기존 stdout 필터가 버리지 않는 형식이다.
    """
    _tag = 'xc%s_dc%s_u3free' % (value_tag(x_c), value_tag(d_c))
    _src = os.path.join(_RUN, HF_ODB)
    _src = _src if os.path.exists(_src) else (HF_ODB if os.path.exists(HF_ODB) else None)

    if str(os.environ.get('SOLARSAIL_NO_SLOPE', '')).strip() == '1':
        emit("[SLOPE] SKIP (SOLARSAIL_NO_SLOPE=1): 설계점 %s slope_stats 아카이빙 생략" % _tag)
        return

    # ---- (1) ODB 설계점 태깅 복사: results/odb/HF_<tag>.odb (태그 규칙 = value_tag, 기존 관례)
    _odb_used = _src
    if str(os.environ.get('SOLARSAIL_NO_ODB_COPY', '')).strip() == '1':
        emit("[ODBCOPY] SKIP (SOLARSAIL_NO_ODB_COPY=1): 설계점 %s ODB 복사 생략" % _tag)
    else:
        try:
            if _src is None:
                emit("[ODBCOPY] !!! HF ODB 를 찾지 못했다(%s) - 복사 생략" % HF_ODB)
            else:
                _odir = os.path.join(_HERE, 'results', 'odb')
                if not os.path.isdir(_odir):
                    os.makedirs(_odir)
                _dst = os.path.join(_odir, 'HF_%s.odb' % _tag)
                shutil.copy2(_src, _dst)
                _odb_used = _dst
                _sz = os.path.getsize(_dst)
                emit("[ODBCOPY] %s -> %s (%.2f MB, %d B)"
                     % (os.path.abspath(_src), os.path.abspath(_dst), _sz / 1048576.0, _sz))
        except Exception as _e:
            emit("[ODBCOPY] 실패(무시): %s" % _e)

    # ---- (2) slope_stats 인프로세스 실행 (Abaqus 안이므로 서브프로세스 금지)
    try:
        if _odb_used is None:
            emit("[SLOPE] !!! 읽을 ODB 가 없다(%s) - 통계 생략" % HF_ODB)
            return
        if _HERE not in sys.path:
            sys.path.insert(0, _HERE)
        import slope_stats
        _S = slope_stats.compute_stats(_odb_used)
        if _S is None:
            emit("[SLOPE] !!! compute_stats 실패(None) - 아카이빙 생략")
            return
        # 점별 results txt 에 slope_ 접두 key=value append (RESULTS 4필드 계약은 그대로)
        _lines = slope_stats.scalar_items(_S, prefix='slope_')
        _ptxt = os.path.join(_HERE, 'results',
                             'results_xc%s_dc%s.txt' % (value_tag(x_c), value_tag(d_c)))
        with io.open(_ptxt, 'a', encoding='utf-8') as _f:
            for _l in _lines:
                _f.write(_l + '\n')
        emit("[SLOPE] %d개 key(slope_ 접두)를 %s 에 append | f1=%.10e f2=%.10e f3=%.10e f4=%.10e"
             % (len(_lines), os.path.basename(_ptxt), _S.get('f1_hf_pred'),
                _S.get('f2_p95'), _S.get('f3_worst_patch'), _S.get('f4_cvar90')))
        # JSON 은 slope_stats 기존 관례: results/slope_<stem>_<step>.json
        try:
            _stem = os.path.splitext(os.path.basename(_odb_used))[0]
            _jpath = os.path.join(_HERE, 'results', 'slope_%s_%s.json' % (_stem, _S['step']))
            _rec = slope_stats.build_json_record(_S, _odb_used)
            slope_stats.write_json_record(_rec, _jpath)
            emit("[SLOPE] JSON 저장: %s" % os.path.abspath(_jpath))
        except Exception as _e:
            emit("[SLOPE] !!! JSON 저장 실패(무시): %s" % _e)
    except Exception as _e:
        emit("[SLOPE] 실패(무시, 런은 계속): %s" % _e)


def _archive_dn(x_c, d_c, dn_lines):
    """[19.106] eval_abaqus.py 가 낸 `DN:` 줄을 점별 results txt 에 `dn_` 접두로 append.

    계약 보존:
      - `RESULTS:` 4필드(lf1,lf2,lf3,hf)는 **읽지도 쓰지도 않는다**.
      - `all_results.csv` 는 이 함수에서 **열지 않는다**(헤더 1회 규칙 - 열 추가 금지).
      - 실패해도 런을 죽이지 않는다(경고 1줄). DN 줄이 아예 없으면 조용히 넘기지 않고 경고한다.
    호출 위치: _archive_result('w') **뒤** (그 함수가 txt 를 다시 쓰므로).
    """
    try:
        _ptxt = os.path.join(_HERE, 'results',
                             'results_xc%s_dc%s.txt' % (value_tag(x_c), value_tag(d_c)))
        if not os.path.exists(_ptxt):
            emit('[DN] !!! 점별 txt 없음(%s) - dn_ 아카이빙 생략' % os.path.basename(_ptxt))
            return
        _kv = {}
        _fail = None
        for _ln in (dn_lines or []):
            _s = str(_ln).strip()
            if not _s.startswith('DN:'):
                continue
            _body = _s[3:].strip()
            if _body.startswith('FAIL'):
                _fail = _body
                continue
            for _tok in _body.split():
                if '=' in _tok:
                    _k, _v = _tok.split('=', 1)
                    _kv[_k if _k.startswith('dn_') else ('dn_' + _k)] = _v

        if _fail is not None and not _kv:
            #   계산 불가: 0 으로 뭉개지 않는다(nan + 사유). CSV 는 그대로.
            with io.open(_ptxt, 'a', encoding='utf-8') as _f:
                _f.write('dn_deg=nan\n')
                _f.write('dn_note=%s\n' % _fail)
            emit('[DN] !!! DN 계산 실패 -> dn_deg=nan (%s)' % _fail)
            return

        if not _kv:
            emit('[DN] !!! 줄 없음 (eval_abaqus.py 출력에 DN: 이 없다) - dn_ 아카이빙 생략')
            return

        _order = ['dn_deg', 'dn_psi_deg', 'dn_n_eq_x', 'dn_n_eq_y', 'dn_n_eq_z',
                  'dn_n_elem', 'dn_area_ref']
        _keys = [k for k in _order if k in _kv] + sorted(k for k in _kv if k not in _order)
        with io.open(_ptxt, 'a', encoding='utf-8') as _f:
            for _k in _keys:
                _f.write('%s=%s\n' % (_k, _kv[_k]))
        emit('[DN] %d개 key(dn_ 접두)를 %s 에 append | dn_deg=%s deg dn_psi_deg=%s '
             'n_eq=(%s,%s,%s) n_elem=%s'
             % (len(_keys), os.path.basename(_ptxt), _kv.get('dn_deg', 'nan'),
                _kv.get('dn_psi_deg', 'nan'), _kv.get('dn_n_eq_x', 'nan'),
                _kv.get('dn_n_eq_y', 'nan'), _kv.get('dn_n_eq_z', 'nan'),
                _kv.get('dn_n_elem', 'nan')))
    except Exception as _e:
        emit('[DN] 실패(무시, 런은 계속): %s' % _e)


def _archive_zk(x_c, d_c, zk_lines, zkiw_lines=None):
    """[paper-match] eval_abaqus.py 가 낸 `ZK:` 줄을 점별 results txt 에 `zk_` 접두로 append.

    `zk_metrics.txt` sidecar 의 추가 키(`zk_Iw`/기준면적/가중모드)도 함께 append 한다.
    [A2] `ZKIW:` 줄이 stdout 의 **정본**이다 -> sidecar 병합이 실패해도 값이 남는다.

    계약 보존:
      - `RESULTS:` 3/4필드는 읽지도 쓰지도 않는다.
      - `all_results.csv` 는 열지 않는다(헤더 1회 규칙 - 열 추가 금지).
    호출 위치: _archive_result('w') 뒤 (그 함수가 txt 를 다시 쓰므로).
    """
    try:
        _ptxt = os.path.join(_HERE, 'results',
                             'results_xc%s_dc%s.txt' % (value_tag(x_c), value_tag(d_c)))
        if not os.path.exists(_ptxt):
            emit('[ZK] !!! 점별 txt 없음(%s) - zk_ 아카이빙 생략' % os.path.basename(_ptxt))
            return
        _kv = {}
        _fail = None
        for _ln in (zk_lines or []):
            _s = str(_ln).strip()
            if not _s.startswith('ZK:'):
                continue
            _body = _s[3:].strip()
            if _body.startswith('FAIL'):
                _fail = _body
                continue
            for _tok in _body.split():
                if '=' in _tok:
                    _k, _v = _tok.split('=', 1)
                    if _k == 'state_wrinkled_frac':
                        _k = 'zk_state_wrinkled_frac'
                    elif not _k.startswith('zk_'):
                        _k = 'zk_' + _k
                    _kv[_k] = _v

        #   [A2] ZKIW: 는 stdout 이 정본이다(sidecar 병합이 실패해도 값이 남는다).
        #     ZK: 4키와 겹치지 않으므로 그대로 덮어쓴다.
        for _ln in (zkiw_lines or []):
            _s = str(_ln).strip()
            if not _s.startswith('ZKIW:'):
                continue
            _body2 = _s[5:].strip()
            if _body2.startswith('FAIL'):
                continue
            for _tok in _body2.split():
                if '=' not in _tok:
                    continue
                _k, _v = _tok.split('=', 1)
                if not _k.startswith('zk_'):
                    _k = 'zk_' + _k      #   sidecar 병합과 같은 키 이름(zk_strain_field 등)
                _kv[_k] = _v

        #   `zk_metrics.txt` sidecar 의 추가 키(`zk_Iw` 등)도 점별 txt 에 남긴다.
        #     ZK: 줄은 public 4키만 갖고, I_w/기준면적/가중모드는 sidecar 에만 있다.
        #   [A1] ZK:FAIL 이어도 sidecar 를 병합한다 — 원인 단서(fields/frame_value)가
        #     거기에만 있고, 다음 런이 zk_metrics.txt 를 덮으면 영구 소실된다.
        #     stale 차단은 "FAIL 에서는 진단 키만 통과"로 대신한다(값 키는 병합 금지).
        _FAILDX = ('fields', 'frame_value', 'reason', 'odb')
        if _fail is not None or _kv:
            try:
                _zp = os.path.join(os.getcwd(), 'zk_metrics.txt')
                if os.path.exists(_zp):
                    with io.open(_zp, 'r', encoding='utf-8', errors='replace') as _zf:
                        for _ln in _zf:
                            _s = _ln.strip()
                            if not _s or '=' not in _s:
                                continue
                            _k, _v = _s.split('=', 1)
                            if _k == 'status':
                                continue
                            if _fail is not None and _k not in _FAILDX:
                                continue
                            if _k == 'state_wrinkled_frac':
                                _k = 'zk_state_wrinkled_frac'
                            elif _k in ('step', 'frame', 'odb', 'weights', 'note'):
                                _k = 'zk_' + _k
                            elif not _k.startswith('zk_'):
                                _k = 'zk_' + _k
                            if _k not in _kv:      #   stdout 우선
                                _kv[_k] = _v
            except Exception as _e:
                emit('[ZK] zk_metrics.txt sidecar merge 실패(무시): %s' % _e)

        with io.open(_ptxt, 'a', encoding='utf-8') as _f:
            if _fail is not None:
                #   계산 불가: 0 으로 뭉개지 않는다(nan + 사유). CSV 는 그대로.
                _f.write('zk_state_wrinkled_frac=nan\n')
                _f.write('zk_note=%s\n' % _fail.replace(' ', '_'))
                for _k in sorted(_kv):
                    _f.write('%s=%s\n' % (_k, _kv[_k]))
                emit('[ZK] !!! ZK 계산 실패 -> zk_state_wrinkled_frac=nan (%s) / 진단키 %d개'
                     % (_fail, len(_kv)))
                return
            if not _kv:
                emit('[ZK] !!! 줄 없음 (eval_abaqus.py 출력에 ZK: 이 없다) - zk_ 아카이빙 생략')
                return
            for _k in sorted(_kv):
                _f.write('%s=%s\n' % (_k, _kv[_k]))
        emit('[ZK] %d개 key를 %s 에 append | state_wrinkled_frac=%s dir=%s'
             % (len(_kv), os.path.basename(_ptxt),
                _kv.get('zk_state_wrinkled_frac', 'nan'), _kv.get('zk_dir_deg', 'nan')))
    except Exception as _e:
        emit('[ZK] 실패(무시, 런은 계속): %s' % _e)


def _r13_ct_keys():
    """[§19.97] 점별 results_*.txt 에 넣을 R-13 Step-ClampTension 판독 키.

    **CSV(all_results.csv) 열은 건드리지 않는다** - 헤더가 1회만 쓰이므로 열을 추가하면 기존 행과
    어긋난다(§19.91). 그래서 진단값은 점별 txt(key=value)로만 남긴다.
    읽지 못한 항목은 'nan' + r13_ct_note 에 이유를 남긴다(조용한 0 금지).
    """
    _d = {'r13_ct_centre_mean_pa': 'nan', 'r13_ct_centre_minp_pa': 'nan',
          'r13_ct_area_frac_neg': 'nan', 'r13_ct_minp_frac_neg': 'nan',
          'r13_ct_uz_max_m': 'nan', 'r13_ct_uz_ave_m': 'nan',
          'r13_ct_odb': 'nan', 'r13_ct_odb_mtime': 'nan', 'r13_ct_stale': 'nan',
          'r13_ct_note': 'probe 미실행(BASE_PROBE=False, ODB 없음, 또는 비-HF 경로)'}
    _d.update(_R13_CT or {})
    return _d


def run_job_safely(job_name, model_name=None, stabilize=None, stabilize_global=None,
                   hard_fail=True, patch_step=None):
    """
    job 실행 중 .odb, .lck 파일 충돌을 방지하고, 완료까지 대기하는 함수

    P1-3: 제출 전에 이전 실행 산출물을 지운다.
      - 성공 판정(not os.path.exists(odb))이 '지난 실행의 odb'를 보고 오판하는 것을 막고,
      - *IMPERFECTION, FILE=Buckle_Analysis 가 stale .fil 을 조용히 읽는 것을 막는다.
    """
    model_name = model_name or MODEL_NAME
    patch_step = patch_step or 'Step-Postbuckle'
    _locked = []
    for _ext in _JOB_ARTIFACTS:
        _f = '%s.%s' % (job_name, _ext)
        if os.path.exists(_f):
            try:
                os.remove(_f)
            except OSError as _e:
                emit("  (warning) could not remove %s: %s" % (_f, _e))
                _locked.append(_f)
    if _locked:
        emit("!!! ERROR: 이전 Abaqus 프로세스가 %s 파일을 잠그고 있습니다." % job_name)
        emit("    잠긴 파일: %s" % ', '.join(_locked))
        emit("    Abaqus/CAE·Standard·pre 프로세스를 종료한 뒤 다시 실행하세요.")
        emit("    (이 상태로 제출하면 새 job이 이전 파일/프로세스와 섞여 결과를 신뢰할 수 없습니다.)")
        sys.stdout.flush()
        sys.exit(1)

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
    emit("   DISP_GLOBAL=%.6g m  SEED_DIV=%.6g  THICKNESS=%.3g m  stabilize(post)=%s  stabilize(global)=%s"
         % (DISP_GLOBAL, SEED_DIV, THICKNESS,
            ('%.4g' % stabilize) if stabilize is not None else 'ladder/기본값',
            ('%.4g' % stabilize_global) if stabilize_global is not None else 'ladder/기본값'))
    emit("   DISP_PRE=%.6g m  DISP_FINAL=%.6g m  u3_policy=%s"
         % (_DISP_PRE, _DISP_FINAL, 'HOLD_ALL' if _HOLD_U3_ALL else 'FREE'))
    emit("   GLOBAL_FINAL=%.6g m  PRETENSION_SCALE=%.6g  (참고: 실제 DISP_FINAL=%.6g m)"
         % (GLOBAL_FINAL, PRETENSION_SCALE, _DISP_FINAL))
    job = mdb.Job(name=job_name, model=model_name, numCpus=_ncp, numDomains=_ncp)
    emit("Submitting Job: %s" % job_name)
    job.writeInput(consistencyChecking=OFF)
    # 제출 전 요소 타입 검증 — 조용한 폴백을 Standard 토큰을 태우기 전에 잡는다(§3).
    #   [§19.112] RELEASE_EDGE_U3=True 는 자유변 면외 구속을 to-ground 스프링으로 바꾼다 ->
    #     덱에 노드별 `Spring1` 요소(n = All_Edges_NoClamp 노드수)가 **정당하게** 들어간다.
    #     종전 가드는 'S' 로 시작하는 타입을 전부 셸 폴백 후보로 봤기 때문에 이 이름을
    #     extra=['Spring1'] 으로 오탐해 잡을 제출 전에 막았다(실측 2026-10-19).
    #     -> True 일 때만 허용목록에 추가하고 개수(경계 노드수)까지 검증한다.
    #     False 는 추가 인자 없이 호출 = 종전 판정/문구 그대로(스프링이 있는데 False 면 종전처럼 실패).
    _edge_kw = {}
    _edge_n = 0
    _allow = [ELEM_CODE_QUAD, ELEM_CODE_TRI]
    if RELEASE_EDGE_U3:
        #   스프링 1개 = 경계 노드 1개. 이 수를 모르면 '타입만 맞다' 로 통과시킬 수 없으므로
        #   검증 불가를 명시하고 제출을 거부한다(조용한 통과 금지).
        _edge_n = len(globals().get('_keep_labels') or ())
        if _edge_n <= 0:
            emit("!!! ERROR: RELEASE_EDGE_U3=True 인데 All_Edges_NoClamp 경계 노드수를 모른다"
                 " (_keep_labels 미정의/빈 셋) -> 스프링 개수를 검증할 수 없어 제출하지 않는다.")
            return False
        _edge_kw = dict(allow_extra=EDGE_SPRING_ELEM_TYPES, expect_extra_total=_edge_n)
        _allow += list(EDGE_SPRING_ELEM_TYPES)
    emit("[GUARD] 허용목록=%s  RELEASE_EDGE_U3=%s%s"
         % ([getattr(_c, 'name', _c) for _c in _allow], RELEASE_EDGE_U3,
            ('  (스프링 기대 개수 %d = All_Edges_NoClamp 노드수)' % _edge_n)
            if RELEASE_EDGE_U3 else '  (스프링 타입 불허)'))
    _ok_et, _cnt_et, _err_et = verify_inp_element_types(
        job_name + '.inp', (ELEM_CODE_QUAD, ELEM_CODE_TRI), **_edge_kw)
    emit("[GUARD] .inp 요소 타입 블록: %s" % (_cnt_et or _err_et))
    if _err_et or not _ok_et:
        emit("!!! ERROR: .inp 요소 타입 검증 실패 (요청=%s): %s"
             % ((ELEM_CODE_QUAD, ELEM_CODE_TRI), _err_et or _cnt_et))
        emit("    허용목록=%s  RELEASE_EDGE_U3=%s  -> Spring1 은 RELEASE_EDGE_U3=True 일 때만 허용된다."
             % ([getattr(_c, 'name', _c) for _c in _allow], RELEASE_EDGE_U3))
        emit("    잡을 제출하지 않는다. ELEM_CODE_* 지정 방식을 고쳐라"
             " (스프링이 정당하면 RELEASE_EDGE_U3 와 스프링 개수를 확인).")
        return False
    # 'write' 는 덱만 생성하고 종료한다(수동 패처 워크플로우).
    if RUN_MODE == 'write':
        emit("[run_abaqus] RUN_MODE='write' : 덱만 생성(%s.inp)하고 제출하지 않는다." % job_name)
        emit("             다음: python riks_patch_input.py aba\\%s.inp <out>.inp Step-Postbuckle "
             "--line-search --ran=5e-3 --can=1e-2 --i0=8 --ir=10" % job_name)
        sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
        sys.exit(0)
    # 'auto' 는 제출 전에 *Controls 를 주입한다. 제자리 덮어쓰기라 잡 이름/입력명이 유지된다.
    if RUN_MODE == 'auto':
        #   [§19.80] 패처는 `Step-Postbuckle`(=HF 전용 스텝)에 *Controls 를 주입한다.
        #     LF 덱에는 그 스텝이 없다 — 실측: "ERROR: *Step, name=Step-Postbuckle 을(를) 찾을 수
        #     없습니다" 로 패처가 실패하고, §14 규칙(패치 안 된 덱은 제출 금지) 때문에 LF 잡이
        #     아예 제출되지 않았다. LF 경로는 그동안 밟지 않아 드러나지 않던 버그다.
        _needs_patch = True
        try:
            with open(job_name + '.inp', 'r', encoding='utf-8', errors='replace') as _fh:
                _needs_patch = (patch_step in _fh.read())
        except Exception:
            _needs_patch = True          # 읽지 못하면 종전 동작(패처 시도) 유지
        if not _needs_patch:
            emit("[PATCHER] 이 덱에 %s 이 없다 -> 패처 생략(§19.80, LF 경로)"
                 % patch_step)
        if _needs_patch and not run_patcher(job_name + '.inp', step=patch_step,
                                             stabilize=stabilize):
            emit("!!! ERROR: 패처 실패 — 패치되지 않은 덱을 제출하지 않는다(§14).")
            emit("RESULTS:FAIL")
            sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
            sys.exit(1)
        #   [§19.92] 패처(riks_patch_input.py)는 **자기 기본값(stabilize=0.01)** 으로 `*Static` 줄을
        #     써 넣는다. 사다리가 넘긴 레벨을 실제 덱에 반영하는 단계가 빠져 있어서, 지금까지
        #     "0.005 시도"는 사실상 전부 0.01 로 돌았다(그래서 (0.8,1.0)이 수동 시험(0.01)과
        #     같은 결과를 냈다). 여기서 덱의 값을 확정한다 — 기록(stabilize_used.txt)과
        #     실제 실행이 일치해야 점 간 비교가 성립한다.
        if _needs_patch and (stabilize is not None or stabilize_global is not None):
            #   [2026-10-21] 두 스텝을 **스텝 범위 한정**으로 바꾼다. 전역 정규식은
            #     Step-ClampTension 등의 `*Static` 까지 잡으므로 금지.
            #     GlobalTension 은 이제 주름이 생기는 스텝이라 사다리 대상이다.
            with open(job_name + '.inp', 'r', encoding='utf-8', errors='replace') as _fh2:
                _ls2 = _fh2.read().splitlines(True)
            _in_gt = _in_ps = False
            _done_gt = _done_ps = False
            for _i2, _l2 in enumerate(_ls2):
                _ll = _l2.lstrip().lower()
                if _ll.startswith('*step'):
                    _in_gt = ('step-globaltension' in _ll)
                    _in_ps = ('step-postbuckle' in _ll)
                    continue
                if (_in_gt and stabilize_global is not None
                        and _ll.startswith('*static') and 'stabilize=' in _ll):
                    _nw, _n2 = re.subn(r'(stabilize=)[0-9.eE+-]+',
                                       r'\g<1>%g' % stabilize_global, _l2)
                    if _n2 == 1:
                        _ls2[_i2] = _nw
                        _done_gt = True
                    _in_gt = False
                    continue
                if (_in_ps and stabilize is not None
                        and _ll.startswith('*static') and 'stabilize=' in _ll):
                    _nw, _n2 = re.subn(r'(stabilize=)[0-9.eE+-]+',
                                       r'\g<1>%g' % stabilize, _l2)
                    if _n2 == 1:
                        _ls2[_i2] = _nw
                        _done_ps = True
                    _in_ps = False
                    continue
            if _done_gt or _done_ps:
                with io.open(job_name + '.inp', 'w', encoding='utf-8', newline='') as _fh2:
                    _fh2.write(''.join(_ls2))
                with open(job_name + '.inp', 'r', encoding='utf-8', errors='replace') as _vf:
                    _vt = _vf.read()
                if stabilize_global is not None:
                    _chk = re.search(
                        r'\*Step,\s*name=Step-GlobalTension.*?\*Static,\s*stabilize=([0-9.eE+-]+)',
                        _vt, re.S)
                    emit("[LADDER] 덱 확인: Step-GlobalTension stabilize=%s (요청 %g)"
                         % (_chk.group(1) if _chk else '?', stabilize_global))
                if stabilize is not None:
                    _chk = re.search(
                        r'\*Step,\s*name=Step-Postbuckle.*?\*Static,\s*stabilize=([0-9.eE+-]+)',
                        _vt, re.S)
                    emit("[LADDER] 덱 확인: Step-Postbuckle stabilize=%s (요청 %g)"
                         % (_chk.group(1) if _chk else '?', stabilize))
            else:
                emit("[LADDER] !!! stabilize 줄을 찾지 못함 "
                     "(GT=%s, PS=%s)" % (_done_gt, _done_ps))
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
        if not hard_fail:
            #   [§19.92] 사다리(run_job_with_ladder)가 이 실패를 **보고 다음 감쇠 레벨로 올라갈 수**
            #     있게 bool 을 돌려준다. 종전에는 여기서 sys.exit(1) 로 프로세스가 끝나
            #     사다리가 실패를 볼 기회 자체가 없었다 → §19.85 자동 상향이 한 번도 작동하지 않음.
            emit("    -> 사다리 모드: 실패를 반환하고 다음 레벨로 넘긴다(§19.92)")
            return False
        sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
        sys.exit(1)
        
    if not job_completed_ok(job_name):
        emit("!!! WARNING: %s — .sta/.msg 에 'HAS COMPLETED SUCCESSFULLY' 없음 "
              "(중도 중단 의심; odb 존재만으로는 판정 불가 — R-12)" % job_name)
    emit("Job %s completed successfully (Status: %s)." % (job_name, str(job.status)))
    return True

def run_job_with_ladder(job_name, model_name=None):
    """[2026-10-21] GlobalTension + Postbuckle **paired** 안정화 사다리.

    u3-free HF 는 GlobalTension 에서 주름이 생기므로, 그 스텝의 감쇠를 내리는 방향으로
    시도한다(2e-4 -> 1e-3 -> 1e-2). Postbuckle 은 정착/최종 하중 스텝으로 5e-4 -> 2e-3 -> 1e-2.
    각 시도는 두 스텝을 같은 인덱스의 레벨로 패치하고, 완주 후 **두 스텝 모두**
    ALLSD/ALLIE 를 읽어 5% 초과 여부를 기록한다.
    """
    def _probe_step(step_name):
        _pe = os.path.join(_HERE, 'probe_energy.py')
        if not os.path.exists(_pe):
            return None, None
        try:
            _pp = subprocess.run(
                'abaqus python "%s" %s %s' % (_pe, job_name + '.odb', step_name),
                shell=True, capture_output=True, text=True)
            _mm = re.search(
                r'ALLSD/ALLIE\s+last=([0-9.eE+-]+)%\s+max=([0-9.eE+-]+)%',
                (_pp.stdout or '') + (_pp.stderr or ''))
            if _mm:
                return float(_mm.group(1)), float(_mm.group(2))
            emit("[LADDER] (%s 에너지 이력 없음 -> probe_energy 출력 참조)" % step_name)
        except Exception as _pe_err:
            emit("[LADDER] %s 에너지 프로브 실패(무시): %s" % (step_name, _pe_err))
        return None, None

    #   paper HF는 Postbuckle이 없으므로 GlobalTension 자체를 line-search/안정화 대상으로
    #   삼는다. old disp HF는 기존처럼 Postbuckle을 대상으로 한다.
    _patch_step = ('Step-Postbuckle' if _HF_USE_POSTBUCKLE else 'Step-GlobalTension')
    for _g_lv, _p_lv in zip(STABILIZE_GLOBAL_LADDER, STABILIZE_POST_LADDER):
        #   paper HF에는 Postbuckle이 없으므로 post stabilize 인자를 넘기지 않는다.
        _stab_post = _p_lv if _HF_USE_POSTBUCKLE else None
        emit("[LADDER] GlobalTension stabilize=%g / Postbuckle stabilize=%s 시도 (%s)"
             % (_g_lv, ('%g' % _p_lv) if _stab_post is not None else 'n/a', job_name))
        if run_job_safely(job_name, model_name=model_name, stabilize=_stab_post,
                          stabilize_global=_g_lv, hard_fail=False,
                          patch_step=_patch_step):
            emit("[LADDER] 완주 (GlobalTension=%g, Postbuckle=%g)" % (_g_lv, _p_lv))
            #   [§19.86] 두 스텝의 ALLSD/ALLIE 를 모두 읽는다. GlobalTension 은 이제
            #     주름이 생기는 스텝이므로 이 값이 선택의 핵심이다.
            _r_last = _r_max = None
            if _HF_USE_POSTBUCKLE:
                _r_last, _r_max = _probe_step('Step-Postbuckle')
            _g_last, _g_max = _probe_step('Step-GlobalTension')
            _vals = [v for v in (_r_max, _g_max) if v is not None]
            _all_max = max(_vals) if _vals else None
            _verdict = 'unknown'
            if _all_max is not None:
                _verdict = ('OK' if _all_max < 1.0 else
                            ('ACCEPT' if _all_max < 5.0 else 'SUSPECT'))
                emit("[LADDER] ALLSD/ALLIE PS(last=%.4g%%,max=%.4g%%) "
                     "GT(last=%.4g%%,max=%.4g%%) -> %s"
                     % (_r_last if _r_last is not None else float('nan'),
                        _r_max if _r_max is not None else float('nan'),
                        _g_last if _g_last is not None else float('nan'),
                        _g_max if _g_max is not None else float('nan'), _verdict))
                if _verdict == 'SUSPECT':
                    emit("!!! WARNING: 안정화 에너지가 5%% 이상 — 이 점의 손실은 "
                         "감쇠가 만든 과소평가일 수 있다")
            try:
                with open(os.path.join(os.getcwd(), 'stabilize_used.txt'), 'w',
                          encoding='utf-8') as _f:
                    _f.write('stabilize=%g\n' % _p_lv)
                    _f.write('stabilize_used=%g\n' % _p_lv)
                    _f.write('stabilize_global=%g\n' % _g_lv)
                    _f.write('allds_allie_last_pct=%s\n'
                             % ('%.6g' % _r_last if _r_last is not None else 'unknown'))
                    _f.write('allds_allie_max_pct=%s\n'
                             % ('%.6g' % _r_max if _r_max is not None else 'unknown'))
                    _f.write('allds_allie_global_last_pct=%s\n'
                             % ('%.6g' % _g_last if _g_last is not None else 'unknown'))
                    _f.write('allds_allie_global_max_pct=%s\n'
                             % ('%.6g' % _g_max if _g_max is not None else 'unknown'))
                    _f.write('allds_allie_all_max_pct=%s\n'
                             % ('%.6g' % _all_max if _all_max is not None else 'unknown'))
                    _f.write('energy_verdict=%s\n' % _verdict)
                    _f.write('job=%s\nx_c=%s\nd_c=%s\n'
                             % (job_name, globals().get('x_c', '?'), globals().get('d_c', '?')))
            except Exception:
                pass
            return True
        emit("[LADDER] 실패 (GlobalTension=%g, Postbuckle=%g) -> 다음 레벨로"
             % (_g_lv, _p_lv))
    emit("!!! ERROR: 사다리 소진 (global=%s, post=%s) — 이 설계점은 완주 불가"
         % (STABILIZE_GLOBAL_LADDER, STABILIZE_POST_LADDER))
    return False


def run_patcher(deck, step='Step-Postbuckle', ctrl=('--line-search', '--ran=5e-3',
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


#   [§19.110] LF 진폭 override — 런타임 인자(fidelity) 파생이라 SETTINGS 밖에 남긴다.
#     근거([§19.81] LF 는 무결함 모델)는 SETTINGS ⑤ 의 IMPERFECTION_AMPL_T 주석 참조.
try:
    _IS_HF_RUN = (fidelity == 'HF')
except NameError:
    _IS_HF_RUN = True                 # 알 수 없으면 종전 동작(HF 가정) 유지
if not _IS_HF_RUN:
    IMPERFECTION_AMPL_T = 0.0

if MODE_SOURCE not in ('external', 'self'):
    raise RuntimeError("MODE_SOURCE 는 'external' 또는 'self' 여야 합니다 (현재 %r)" % (MODE_SOURCE,))

# 임퍼펙션 주입 방식 — 모드 소스의 스텝 타입이 방식을 결정한다(§6).
#   *BUCKLE   -> 'odb_direct'(기본). *BUCKLE 은 .fil 출력이 금지되므로 모드 프레임은 ODB 에만 있다.
#   *FREQUENCY -> 'file' (.fil + *IMPERFECTION, FILE=).

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
        _srcdirs.append(_HERE)      # code\
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
    with open(path, "r", encoding="utf-8", errors="replace") as f:
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


# =====================================================================
# [§19.104] 임퍼펙션 표 선택 - 클램프 없는 baseline 은 **순수 seed 가족**만 쓴다
#   발견된 오염(2026-10-16 실런): `NO_CLAMP=True`(클램프 없는 모델)인데 임퍼펙션은 합성표
#   `modes_xc<x>.txt` 였고 그 내용이 `modes_seed.txt[1,2] + modes_clamp_xc0p5_dc0.txt[1,2]` 였다.
#   클램프 케이스의 좌굴모드는 클램프 없는 모델의 모드가 아니므로 **오염**이다(과거 로그의
#   "성공한 모드가 클램프 근처에만 몰린다"도 같은 원인으로 읽힌다).
#   근거: 2026-10-01 문헌(Galhofo) 재현 성공값(u_z,ave 97.7 t vs 논문 91.1 t = 1.072배)은
#   **합성표(2+2) 도입 전**의 값이다 -> 그 성공은 seed 가족만으로 얻어졌다. Galhofo 도 "케이블
#   없는 only membrane 으로 모드를 뽑아 실제 모델에 주입"했다 -> seed 가족 사용이 문헌 방식과 일치.
#   그래서 NO_CLAMP=True 는 (a) 합성표 대신 seed 가족 표를 쓰고 (b) 그 표의 상위 4모드를 쓴다.
#   NO_CLAMP=False 경로(합성표 + 설계점 지문 게이트)는 종전과 **동일**하다.
# =====================================================================


def seed_family_gate(table_path, max_modes):
    """[§19.104] 클램프 없는 baseline 이 쓸 **순수 seed 가족 표**인지 검증한다.

    반환 (ok, 사유, 요약dict). 요약은 실패해도 채워서 '무엇을 쓰려 했는지'를 로그에 남긴다.
      ① 파일이 있는가  ② 파일명이 seed 가족 표인가(합성표/클램프표 차단)
      ③ 헤더에 합성/클램프 흔적이 없는가  ④ 모드가 max_modes 개 이상인가
    stdlib 전용이라 Abaqus 없이 단위시험할 수 있다(code/seed_family_stub.py).
    """
    summary = {'table': table_path, 'modes': [], 'source': '?', 'header': ''}
    if not table_path or not os.path.exists(table_path):
        return False, 'seed 가족 표 파일이 없음: %s' % table_path, summary
    if os.path.basename(table_path) != SEED_FAMILY_NAME:
        return False, ('seed 가족 표가 아니다(합성표/클램프표): %s (기대 %s)'
                       % (os.path.basename(table_path), SEED_FAMILY_NAME)), summary
    nums = _available_mode_numbers(table_path)
    summary['modes'] = nums[:max_modes]
    summary['source'] = mode_table_info(table_path).get('source', '?')
    try:
        hdr = []
        with open(table_path, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                if not line.startswith('#'):
                    break
                hdr.append(line)
        summary['header'] = ''.join(hdr).strip()
    except Exception as e:
        return False, '표 헤더를 읽지 못함(%s: %s)' % (type(e).__name__, e), summary
    bad = [t for t in ('modes_clamp', 'combined', 'clamp_sha') if t in summary['header'].lower()]
    if bad:
        return False, ('표 헤더에 클램프/합성 흔적 %s -> 오염(seed 가족이 아니다)' % bad), summary
    #   [paper-match] 2차 격자(MESH_ORDER=2)는 노드 라벨이 1차와 다르다. 표의 source ODB 이름에
    #     'grid2' 가 없으면 1차 격자에서 뽑은 표이므로 주입하면 안 된다(조용한 오염 차단).
    if MESH_ORDER == 2:
        _src = str(summary.get('source', ''))
        if 'grid2' not in os.path.basename(_src):
            return False, ('MESH_ORDER=2 인데 시드표 source 가 grid2 가 아니다: %s '
                           '-> make_seed.py 로 2차 격자 모드를 다시 뽑아라' % _src), summary
    if len(nums) < max_modes:
        return False, ('모드 부족: %d개(모드 %s) < %d개 필요(상위 %d개를 주입한다)'
                       % (len(nums), nums, max_modes, max_modes)), summary
    return True, 'OK', summary


#   [2026-10-12] 설계점별 시드표 - 기본은 결정적 경로(seeds/modes_xc<NNN>.txt),
#     실험은 env SOLARSAIL_MODE_TABLE 로 덮는다. 종전 상수 '..\modes_2p2.txt' 는 x_c=0.5 산출물
#     하나만 가리켜, MFBO(x_c 스윕)에서 조용히 같은 시드를 재사용했다(계획 Q1=a / Q2=a).
SEEDS_DIR = os.path.join(_HERE, SEED_DIR_NAME)   # [§19.73] cae noGUI 는 execfile 이라 __file__ 이 없다
SEED_FAMILY_TABLE = os.path.join(SEEDS_DIR, SEED_FAMILY_NAME)   # [§19.104] 순수 seed 가족 표
_MODE_TABLE_ENV = os.environ.get('SOLARSAIL_MODE_TABLE')
if NO_CLAMP and not USE_COMPOSITE_MODE_TABLE:
    # 순수 seed 가족 표(클램프 모드 배제) — 종전 §19.104 경로.
    MODE_TABLE = SEED_FAMILY_TABLE
else:
    # [paper-match] 기본: seed + clamp_only **합성표**. 클램프가 없는 HF 라도 모드 형상 기저에는
    #   두 가족을 함께 쓴다(사용자 결정 2026-10-20). 지문 게이트(seed_gate)가 설계점/소스 드리프트를 막는다.
    MODE_TABLE = seed_table_path(SEEDS_DIR, x_c, d_c, env=_MODE_TABLE_ENV)


if IMPERFECTION_MODE == 'odb_table':
    #   [§19.104] 클램프 없는 baseline -> 순수 seed 가족 표. 여기서 검증/로그를 끝낸다:
    #     표가 없거나 모드가 부족하면 **크게 찍고 중단**한다(조용히 합성표로 되돌아가면 오염이 재발).
    if NO_CLAMP and not USE_COMPOSITE_MODE_TABLE:
        _sf_ok, _sf_why, _sf = seed_family_gate(MODE_TABLE, IMPERFECTION_MAX_MODES)
        emit("[SEED] NO_CLAMP=True -> 순수 seed 가족 사용(clamp 모드 제외) / "
             "table=%s / modes=%s / source=%s / x_c=%g d_c=%g"
             % (_sf.get('table', MODE_TABLE), _sf.get('modes', []), _sf.get('source', '?'),
                x_c, d_c))
        emit("[SEED] seed 가족 판정: %s (%s)" % ('OK' if _sf_ok else '거부', _sf_why))
        if not _sf_ok:
            emit("!!! ERROR: NO_CLAMP=True baseline 은 순수 seed 가족 표가 필요합니다 -> 중단: %s"
                 % _sf_why)
            emit("[ERROR-EN] NO_CLAMP baseline requires the pure seed-family mode table"
                 " -> aborting BEFORE HF job submit.")
            emit("    기대 파일: %s" % SEED_FAMILY_TABLE)
            emit("    재생성: python3 make_seed.py --x_c %g --d_c %g   ([1/3] 단계가 seed 가족 표를 만든다)"
                 % (x_c, d_c))
            sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
            sys.exit(1)
        # 표가 가진 모드 중 상위 IMPERFECTION_MAX_MODES 개([1,2,3,4] 규약). 클램프 모드는
        #   seed 가족 표에 없으므로 구조적으로 섞일 수 없다.
        IMPERFECTION_MODES = tuple(_sf['modes'])
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
    #   [2026-10-12] 설계점 시드 게이트 (계획 Q2=a: 불일치면 **제출 전 중단**).
    #     판정: ① 표가 이 x_c/d_c 에서 나왔는가 ② 만든 뒤 buckle 소스가 바뀌지 않았는가
    #           ③ 소스 모드표가 그 뒤 바뀌지 않았는가.
    #     로직은 aba_imperfection.seed_gate(stdlib) — 단위시험 code/test_seed_path.py.
    #   [§19.104] 이 게이트는 **합성표 전용**이다(합성표 헤더의 x_c/fingerprint 를 본다). 클램프
    #     없는 baseline 의 순수 seed 가족 표에는 그 헤더가 없으므로 여기서는 검증하지 않고,
    #     위 seed_family_gate 가 대신한다(검증 없는 통과가 아니라 대상이 다르다).
    if NO_CLAMP and not USE_COMPOSITE_MODE_TABLE:
        emit("[SEED] 합성 지문 게이트 생략 (NO_CLAMP=True 순수 seed 가족 표, §19.104)")
    else:
        _buckle_script = os.path.join(_HERE, BUCKLE_SCRIPT_NAME)
        _gate_ok, _gate_why, _sinfo = seed_gate(x_c, d_c, MODE_TABLE, _buckle_script)
        emit("[SEED] x_c=%.4g d_c=%.4g table=%s" % (x_c, d_c, MODE_TABLE))
        emit("[SEED] 모드 %s / source_seed=%s / source_clamp=%s / 지문=%s"
             % (sorted(_mode_table.keys()), _sinfo.get('seed', '?'), _sinfo.get('clamp', '?'),
                _sinfo.get('fingerprint', '?')))
        if not _gate_ok:
            emit("!!! ERROR: 시드표 게이트 실패 - %s" % _gate_why)
            emit("[ERROR-EN] seed table gate failed -> aborting BEFORE HF job submit.")
            emit("    재생성: python3 make_seed.py --x_c %g --d_c %g --force" % (x_c, d_c))
            try:
                with open('mode_table_check.txt', 'w') as _fpc:
                    _fpc.write('GATE FAIL\n%s\ntable=%s\ntime=%s\n'
                               % (_gate_why, os.path.abspath(MODE_TABLE),
                                  time.strftime('%Y-%m-%d %H:%M:%S')))
            except Exception:
                pass
            sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
            sys.exit(1)
        emit("[SEED] 게이트 OK - 이 표는 이 설계점·이 buckle 소스와 일치한다")
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


def clamp_coord_L(x): return (10-10*x, 10-10*x, 0)
def clamp_coord_R(x): return (10+10*x, 10-10*x, 0)

V_CL = clamp_coord_L(x_c)
V_CR = clamp_coord_R(x_c)

# ---- 런타임 파생 하중: SETTINGS ④ 의 값 x CLI d_c (정의는 여기 — d_c 가 필요하다) ----
#   [2026-10-20] 클램프 프리텐션 = 초기응력. d_c = 클램프/꼭짓점 케이블 **응력비**.
#     단, 종전 MFBO 기준선(LOAD_MODE='disp')은 변위 당김을 그대로 유지한다.
if LOAD_MODE == 'cable_stress':
    CLAMP_PRESTRESS = CABLE_PRESTRESS * d_c
    CLAMP_PULL = 0.0             # 변위 당김 없음
    CLAMP_FINAL = 0.0            # 변위 당김 없음
    CLAMP_PERT = PERTURBATION    # 자기좌굴 프로브용(설계하중 아님, d_c 무관)
else:
    CLAMP_PRESTRESS = None       # 이 경로는 초기응력을 쓰지 않는다
    CLAMP_PULL = DISP_GLOBAL * d_c
    CLAMP_FINAL = GLOBAL_FINAL * d_c
    CLAMP_PERT = PERTURBATION * d_c


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
        #   [2026-10-12 계획 T3] 경고로 두면 드리프트가 조용한 오주입으로 남는다 -> 중단.
        emit("!!! ERROR: 모드표 노드 %d개 중 %d개만 적용 -> 메쉬/라벨 불일치(오주입 방지)"
             % (len(_labels), len(_apply)))
        emit("[ERROR-EN] imperfection node mapping incomplete -> aborting.")
        sys.stdout.flush()
        sys.exit(1)
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
mat_cable.Density(table=((CABLE_RHO, ),)) 
mat_cable.Elastic(table=((CABLE_E, CABLE_NU),))
my_model.TrussSection(name='Section-Cable', material='Kevlar', area=CABLE_AREA)

# 파트 생성: 멤브레인
#   [§19.110] 메쉬 종류 스위치(USE_GRID_MESH)의 정의와 근거는 SETTINGS ② 에 있다.
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
        my_model, part_name='Membrane', base=BASE, height=HEIGHT, seed_div=SEED_DIV,
        elem_order=MESH_ORDER)
    print("[mesh] grid(.inp import, orphan): order=%d, 노드 %d / quad %d + tri %d = %d"
          % (MESH_ORDER, _nn, _n4, _n3, _n4 + _n3))
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
    #   [paper-match] 케이블 초기응력(*INITIAL CONDITIONS, TYPE=STRESS)용 요소셋.
    #     절점셋을 주면 'not an element set' 으로 죽는다(막의 All_Elem 과 같은 이유).
    p_c.Set(elements=p_c.elements, name='All_Elem')
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
#   [§19.110] NO_CLAMP 정의는 SETTINGS ① 한 곳이다. 여기 있던 인라인 재대입
#     (`NO_CLAMP = True`)은 같은 값의 중복이라 지웠다 — 정의(1028줄)와 재대입(1421줄)
#     사이에 NO_CLAMP 를 쓰는 문장이 없어 유효값은 종전과 동일하다(상수 대조로 확인).
if NO_CLAMP:
    emit("[run_abaqus] NO_CLAMP=True : 클램프(cable_CL/CR + 강체패치 + BC) 없이 모델링 (대조 실험)")
    rp_cl_obj = rp_cr_obj = rp_cl_reg = rp_cr_reg = None
else:
    #   [2026-10-12 §19.58] 0.2 -> 0.21: 거리 0.2 노드에서 부동소수 경계로 좌우 패치 노드 수가
    #   달라져(L6/R7) 클램프 부근 최소 주응력이 비대칭이었다. 0.21 은 8/8 대칭(실측)이고
    #   CLAMP_EXCL_R 과 같아 제외가 패치를 덮는다. 빌드 시점에 대칭을 강제 검사한다.
    rp_cl_obj, rp_cl_reg = create_rigid_patch('CL', V_CL, radius=PATCH_R_CLAMP)
    rp_cr_obj, rp_cr_reg = create_rigid_patch('CR', V_CR, radius=PATCH_R_CLAMP)
    _clm = aba_grid_mesh.patch_coords_from_set(a.sets['CL_Nodes'])
    _crm = aba_grid_mesh.patch_coords_from_set(a.sets['CR_Nodes'])
    if _clm is None or _crm is None:
        print("[대칭] 경고: Set 노드좌표 API 없음 -> 클램프 패치 좌우 대칭 검사 건너뜀(스텁 전용 경고).")
    else:
        _okm, _msgm = aba_grid_mesh.patch_nodes_mirror_ok(_clm, _crm)
        if not _okm:
            raise RuntimeError("[클램프 패치 대칭] %s — 좌우 클램프 부근 응력이 달라진다(§19.58)."
                               % _msgm)
        print("[대칭] 클램프 패치 미러검사 %s" % _msgm)

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
#   [2026-10-21] HF 는 여기서 주름이 생긴다. initialInc 를 0.001 로 낮추고
#     maxNumInc 를 10000 으로 늘린다. LF 는 기존 값(0.01/1000, GlobalTension 1.0)을 유지한다.
my_model.StaticStep(
    name='Step-GlobalTension',
    previous='Initial', 
    nlgeom=ON,
    stabilizationMagnitude=(STABILIZE_GLOBAL_HF if _IS_HF_RUN else STABILIZE_GLOBAL),
    stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
    initialInc=((0.001 if _IS_HF_RUN else 0.01) if PAPER_MATCH else 0.0001),
    minInc=1e-8,
    maxInc=0.01,
    maxNumInc=(10000 if _IS_HF_RUN else 1000)
    #   [paper-match] 초기응력은 스텝 시작에 전부 들어간다(램프 없음) -> 첫 증분이 너무 작으면
    #     큰 잔차를 작은 증분으로 풀다 컷백이 반복된다. 논문 모드에서는 0.001(HF)/0.01(LF) 로 시작한다.
)

# Step Clamp Tension : 클램프에 변위 가하기
my_model.StaticStep(
    name='Step-ClampTension',
    previous='Step-GlobalTension',
    nlgeom=ON,
    stabilizationMagnitude=0.0002,      # Galhofo Reference
    stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
    initialInc=0.0001, 
    minInc=1e-8, 
    maxInc=0.01,
    maxNumInc=1000    # P1-2: 1e-15 는 발산 시 증분 폭주
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
#   [§19.110] 솔버 선택(_EIGENSOLVER)의 정의/근거는 SETTINGS ⑥ 에 있다.
emit("[run_abaqus] buckle eigensolver=%s numEigen=%d vectors=%s"
      % (_EIGENSOLVER, N_EIG_BUCKLE, (BUCKLE_VECTORS if _EIGENSOLVER != 'LANCZOS' else 'n/a')))
if _EIGENSOLVER == 'LANCZOS':
    my_model.BuckleStep(
        name='Step-Buckle',
        # (B) base state = '클램프 당김 없는(GlobalTension 말단)' 상태.
        #   [§19.96 정정] 클램프 자체는 Initial 에서 생성되어 GlobalTension 내내 구속이 살아 있다
        #     (Tie_CL/Tie_CR + BC_CL/BC_CR u3=0 + Disp_Control_CL/CR u1=u2=0 — 위 구속조건 블록).
        #     이 상태에서 빠지는 것은 '클램프 변위 당김'뿐이다. 클램프 프리텐션은 초기응력으로 바뀌었다.
        #     따라서 이 값은 클램프 구속(하중 분담)을 포함한다 — Galhofo 7000 Pa 같은
        #     '클램프 없는 프리텐션'과 대조하려면 NO_CLAMP=True 모델의 GlobalTension 말단을 써야 한다.
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

#   [paper-match] LF/HF 공통: 케이블 초기응력 + 클램프 초기응력(d_c*CABLE_PRESTRESS).
#     HF 는 이 상태를 Postbuckle 의 base 로 쓰고, LF 는 이 상태를 **최종 프리스트레스**로 쓴다.
if LOAD_MODE == 'cable_stress' and _PAPER_PATH:
    for _nm in ('Cable_Top', 'Cable_Right', 'Cable_Left'):
        my_model.Stress(
            name='PreTension_' + _nm,
            region=a.instances['Inst_' + _nm].sets['All_Elem'],
            distributionType=UNIFORM,
            sigma11=CABLE_PRESTRESS,
            sigma22=0.0, sigma33=0.0,
            sigma12=0.0, sigma13=0.0, sigma23=0.0)
    emit("[LOAD] cable_stress: sigma11=%.4g Pa x 3 cables (A=%.4g m2 -> F=%.4g N/cable)"
         % (CABLE_PRESTRESS, CABLE_AREA, CABLE_PRESTRESS * CABLE_AREA))
    #   [2026-10-20] 클램프도 초기응력. d_c = 클램프/꼭짓점 케이블 응력비.
    #     인스턴스 이름은 connect_cable('cable_CL', ...) -> 'Inst_cable_CL' (소문자 c).
    if not NO_CLAMP:
        for _nm in ('cable_CL', 'cable_CR'):
            my_model.Stress(
                name='PreTension_' + _nm,
                region=a.instances['Inst_' + _nm].sets['All_Elem'],
                distributionType=UNIFORM,
                sigma11=CLAMP_PRESTRESS,
                sigma22=0.0, sigma33=0.0,
                sigma12=0.0, sigma13=0.0, sigma23=0.0)
        emit("[LOAD] clamp_stress: sigma11=%.4g Pa (= d_c %.3f x cable %.4g Pa)"
             % (CLAMP_PRESTRESS, d_c, CABLE_PRESTRESS))

# 시작하면서 전체 면 z 변위 고정
#   [paper-match] 논문은 정적 스텝에서 **세 변만** u3=0 이다.
#     - 종전 disp / physical HF(PAPER_HOLD_U3_GLOBAL=True): 전막 u3=0 (BC_Stabilize_Z)
#     - LF(NO_CLAMP) / paper-only: 전막 고정 없이 BC_Edges_Only_Z 만
if _HOLD_U3_ALL:
    my_model.DisplacementBC(
        name='BC_Stabilize_Z', 
        createStepName='Initial', 
        region=inst_memb.sets['All'],
        u3=SET
    )

# 구속조건 생성
#   [§19.105] BC_Anchor 는 케이블 '바깥 끝'(end_c1)의 구속이다. 종전에는 u1~ur3 완전 고정이었다.
#     PULL_TOP=True 면 이 점도 케이블 방향 지정변위(Disp_Control_Top)를 받아야 하므로 면내 u1/u2 를 푼다.
#     남기는 것은 u3(면외) + ur1~ur3(회전): 평면성과 강체 면외/회전을 막고, 면내 강체운동은
#     Top/Right/Left 세 점(비공선)의 u1=u2=0(Initial, Disp_Control_*)이 이미 막는다 -> 특이행렬 아님.
if PULL_TOP:
    my_model.DisplacementBC(name='BC_Anchor', createStepName='Initial', region=end_c1,
                            u3=0, ur1=0, ur2=0, ur3=0
                            )
else:
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
if PULL_TOP:
    #   [§19.105] Right/Left 와 완전히 같은 패턴: Initial 에서 u1=u2=0, 스텝에서 케이블 방향 지정변위.
    #     하중은 지정 변위만 쓴다(CLOAD 없음). 방향은 위 top_cos/top_sin (+y 축).
    my_model.DisplacementBC(name='Disp_Control_Top', createStepName='Initial', region=end_c1, 
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
#   [paper-match] cable_stress 모드에서는 _DISP_PRE=0 -> BC 가 Initial 의 0 에 머문다.
my_model.boundaryConditions['Disp_Control_Right'].setValuesInStep(stepName='Step-GlobalTension', 
    u1=_DISP_PRE * cos_val,
    u2=-_DISP_PRE * sin_val
)
# Left Cable: 좌하향 당기기
my_model.boundaryConditions['Disp_Control_Left'].setValuesInStep(stepName='Step-GlobalTension', 
    u1=-_DISP_PRE * cos_val,
    u2=-_DISP_PRE * sin_val
)
if PULL_TOP:
    #   [§19.105] Top Cable: 상향(+y) 당기기 - Right/Left 와 같은 대입식.
    my_model.boundaryConditions['Disp_Control_Top'].setValuesInStep(stepName='Step-GlobalTension', 
        u1=_DISP_PRE * top_cos * TOP_PULL_SCALE,
        u2=-_DISP_PRE * top_sin * TOP_PULL_SCALE
    )
    emit("[BC] PULL_TOP=True -> Top 꼭짓점도 케이블 방향 변위(프리텐션 %.6g m / 최종 %.6g m), BC_Anchor 는 u3 만"
         % (DISP_GLOBAL, GLOBAL_FINAL))
    if TOP_PULL_SCALE != 1.0:
        #   [§19.109] 1.0 이 아니면 실제 적용값을 1줄 더 남긴다(1.0 은 종전 로그와 바이트 동일).
        emit("[BC] PULL_TOP=True TOP_PULL_SCALE=%.6g -> Top(u1=%.6e,u2=%.6e)"
             % (TOP_PULL_SCALE,
                DISP_GLOBAL * top_cos * TOP_PULL_SCALE,
                -DISP_GLOBAL * top_sin * TOP_PULL_SCALE))
    emit("   [BC] 케이블 방향 (u1=v*cos, u2=-v*sin): Top cos=%.6f sin=%.6f / Right cos=%.6f sin=%.6f / Left cos=%.6f sin=%.6f"
         % (top_cos, top_sin, cos_val, sin_val, -cos_val, sin_val))

#   [2026-10-20] cable_stress 경로는 클램프 프리텐션을 초기응력으로 준다 -> 변위 당김 없음.
#     종전 MFBO 기준선(LOAD_MODE='disp')만 아래 변위 당김을 유지한다.
if not NO_CLAMP and LOAD_MODE == 'disp':
    my_model.boundaryConditions['Disp_Control_CL'].setValuesInStep(
        stepName='Step-ClampTension', u1=-(CLAMP_PULL/sqrt2), u2=(CLAMP_PULL/sqrt2))
    my_model.boundaryConditions['Disp_Control_CR'].setValuesInStep(
        stepName='Step-ClampTension', u1=(CLAMP_PULL/sqrt2), u2=(CLAMP_PULL/sqrt2))

# Buckle Perturbation
my_model.boundaryConditions['Disp_Control_Right'].setValuesInStep(stepName='Step-Buckle', 
    u1=PERTURBATION * cos_val,
    u2=-PERTURBATION * sin_val
)
my_model.boundaryConditions['Disp_Control_Left'].setValuesInStep(stepName='Step-Buckle', 
    u1=-PERTURBATION * cos_val,
    u2=-PERTURBATION * sin_val
)
if PULL_TOP:
    #   [§19.105] 좌굴 perturbation 도 세 꼭짓점에 같은 크기로 준다(한 점만 다른 크기면
    #     좌굴 base state 가 비대칭이 되어 모드가 인위적 비대칭을 담는다).
    my_model.boundaryConditions['Disp_Control_Top'].setValuesInStep(stepName='Step-Buckle', 
        u1=PERTURBATION * top_cos * TOP_PULL_SCALE,
        u2=-PERTURBATION * top_sin * TOP_PULL_SCALE
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

#   [2026-10-21] Step-Buckle 은 어차피 삭제되므로 BC_Stabilize_Z 해제는 죽은 코드였다.
#     paper HF 는 BC 가 없고, old disp HF 는 Step-Postbuckle 에서 해제한다.

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

#   [§19.111] 자유변 면외 구속: 기본은 하드 u3=0(BC_Edges_Only_Z). RELEASE_EDGE_U3=True 면 이 BC 를
#     만들지 않고 **노드별 to-ground 스프링(dof=3)** 으로 대체한다.
#     왜: 가장 긴 변을 면외 고정하면 국부 압축이 주름으로 이완되지 못해 압축 무게중심이 바닥변
#       중앙(10.03, 0.24)에 몰리는 경계 아티팩트가 생긴다(§19.111).
#     스프링은 model-level engineeringFeature 다 -> 스텝 삭제(Step-Buckle -> Step-Postbuckle)
#       와 무관하게 모델에 남는다(BC_Edges_Only_Z 만 스텝과 함께 사라진다).
#     Abaqus 덱: *ELEMENT, TYPE=SPRINGA (노드-지면 1개씩) + *SPRING, DOF=3, ELSET=<셋>
if RELEASE_EDGE_U3:
    my_model.rootAssembly.engineeringFeatures.SpringDashpotToGround(
        name='EdgeZ_Spring',
        region=a.sets['All_Edges_NoClamp'],
        orientation=None,
        dof=3,
        springBehavior=ON,
        springStiffness=EDGE_SPRING_K,
        dashpotBehavior=OFF,
    )
    #   면외 강체모드 방지(3점 u3): u3 하드 구속이 사라지면 z 강체병진이 특이모드가 된다.
    #     꼭짓점 강체패치 **제어점(RP)** 에 u3=0 을 건다 - 좌굴 스크립트의 BC_Right_Z/BC_Left_Z 와
    #     같은 역할/같은 대상(RP)이다. HF 의 BC_Right/BC_Left 는 케이블 **바깥 끝**(end_c2/end_c3)
    #     이라 코너 RP 의 u3 를 잡지 못한다(트러스 T3D2 는 횡강성이 0).
    #     Top 도 포함한다: BC_Anchor 는 end_c1(케이블 바깥 끝)이라 같은 이유로 RP_Top_Set 의
    #     u3 를 잡지 못한다 -> 세 꼭짓점 모두 RP u3=0 이어야 3점 지지가 성립한다.
    #     앵커/클램프 자체의 u3(BC_Anchor/BC_CL/BC_CR)는 건드리지 않는다.
    _edge_u3_rps = ['RP_Top_Set', 'RP_Right_Set', 'RP_Left_Set']
    if not NO_CLAMP:
        #   클램프 패치도 경계 u3 구속에서 빠지므로(CLAMP_EXCL_R=0 이면 원래 겹쳤다) 같은 이유로
        #   제어점 u3 를 잡는다. 클램프는 물리적으로 면외 고정 부품이다.
        _edge_u3_rps += ['RP_CL_Set', 'RP_CR_Set']
    for _rp_name in _edge_u3_rps:
        my_model.DisplacementBC(
            name='BC_%s_Corner_Z' % _rp_name.replace('RP_', '').replace('_Set', ''),
            createStepName='Initial',
            region=a.sets[_rp_name],
            u3=0
        )
    emit("[BC] RELEASE_EDGE_U3=True -> BC_Edges_Only_Z 제거, EdgeZ_Spring k=%.3e N/m (alpha=%.1e, n=%d), Right/Left u3 지지 확인"
         % (EDGE_SPRING_K, EDGE_SPRING_ALPHA, len(_keep_labels)))
else:
    my_model.DisplacementBC(
        name='BC_Edges_Only_Z',
        #   [paper-match] 논문 정적 스텝은 세 변 u3=0 을 처음부터 건다.
        createStepName=('Initial' if _PAPER_PATH else 'Step-Buckle'),
        region=a.sets['All_Edges_NoClamp'],
        u3=0
    )

#   [2026-10-12 §19.66] LF ODB 는 **설계점별 이름**이다(LF_Analysis_xc<NNN>_dc<NNN>.odb).
#     고정 이름은 이전 설계점의 LF 가 이번 점의 지표로 둔갑할 수 있었다(짝 오염).
#     HF ODB 는 고정 이름을 유지한다 — 즉시 eval 로 소비되고 남는 것은 수치이며, 점마다
#     새 이름을 주면 수 GB 가 쌓인다.
LF_ODB = lf_odb_name(x_c, d_c)
LF_JOB = lf_job_name(x_c, d_c)

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

    if LOAD_MODE == 'cable_stress':
        #   [paper-match] LF/HF 물리 정합: 초기응력 프리스트레스 상태가 **최종 상태**다.
        #     GlobalTension 뒤에 빈 ClampTension/HighTension 스텝을 두지 않는다.
        #     HF 도 Postbuckle 직전에 ClampTension 을 삭제하고 GlobalTension 에서 시작한다.
        #     -> LF ODB 의 마지막(그리고 유일한) 해석 스텝 = Step-GlobalTension.
        if 'Step-ClampTension' in my_model.steps:
            del my_model.steps['Step-ClampTension']
        emit("[LF] paper 경로: Step-GlobalTension 이 최종 프리스트레스 스텝 "
             "(ClampTension/HighTension 없음)")
    else:
        #   ---- 구 MFBO 기준선(LOAD_MODE='disp') 경로: 바이트/거동 보존 ----
        my_model.StaticStep(
            name='Step-HighTension',
            previous='Step-ClampTension',
            nlgeom=ON,
            initialInc=0.01,
            maxNumInc=100
        )

        if 'BC_Stabilize_Z' in my_model.boundaryConditions:
            my_model.boundaryConditions['BC_Stabilize_Z'].setValuesInStep(
                stepName='Step-HighTension',
                u3=0.0
            )

        #   [2026-10-20] 구 disp 경로만 GLOBAL_FINAL 1 mm 를 유지한다.
        for name, sign in [('Disp_Control_Right', 1), ('Disp_Control_Left', -1)]:
            my_model.boundaryConditions[name].setValuesInStep(
                stepName='Step-HighTension',
                u1=sign * _DISP_FINAL * cos_val,
                u2=-_DISP_FINAL * sin_val
            )
        if PULL_TOP:
            #   [§19.105] 최종 하중도 세 꼭짓점 동일 - Top 은 +y 방향 _DISP_FINAL.
            my_model.boundaryConditions['Disp_Control_Top'].setValuesInStep(
                stepName='Step-HighTension',
                u1=_DISP_FINAL * top_cos * TOP_PULL_SCALE,
                u2=-_DISP_FINAL * top_sin * TOP_PULL_SCALE
            )
        for name, sign in ([] if NO_CLAMP else [('Disp_Control_CL', -1), ('Disp_Control_CR', 1)]):
            my_model.boundaryConditions[name].setValuesInStep(
                stepName='Step-HighTension',
                u1=sign * CLAMP_FINAL / sqrt2,
                u2=CLAMP_FINAL / sqrt2)

    #   [paper-match] 변형률 = Green-Lagrange E. Abaqus NLGEOM=ON 은 ODB 에 E 대신
    #     LE(logarithmic)를 쓰는 경우가 있다(실측: Lean field 만 존재). eval_abaqus 가
    #     LE 주성분을 E=0.5*(exp(2*LE)-1) 로 변환하므로 LE 를 명시적으로 요청한다.
    #   [2026-10-23] F-Output-1 에 E 를 setValues 로 덧붙이는 방식이 ODB 에 반영되지 않는
    #     사례가 확인됐다(tag v2-e1 인데도 E_field_missing). ZK 전용 **별도 요청 F-ZK** 를
    #     선택 LF 스텝에 만들고 frequency=1 로 고정한다. F-Output-1 은 S/U/COORD/EVOL/RF 만 남긴다.
    my_model.fieldOutputRequests['F-Output-1'].setValues(
        variables=('S', 'U', 'COORD', 'EVOL', 'RF'),
        frequency=10   # R-10: odb 크기 절감
    )
    _zk_output_step = ('Step-HighTension'
                       if LOAD_MODE != 'cable_stress' else 'Step-GlobalTension')
    try:
        if 'F-ZK' in my_model.fieldOutputRequests:
            del my_model.fieldOutputRequests['F-ZK']
        my_model.FieldOutputRequest(
            name='F-ZK',
            createStepName=_zk_output_step,
            variables=('LE',),
            frequency=1,   # 마지막 프레임 포함 보장
        )
        emit("[LF] field output F-ZK: step=%s variables=('LE',) frequency=1"
             % _zk_output_step)
    except Exception as _fze:
        #   API/버전 차이 폴백: 기존 요청에 LE 를 넣고 주기를 1 로.
        emit("[LF] !!! F-ZK 생성 실패(%s) -> F-Output-1 폴백" % _fze)
        try:
            my_model.fieldOutputRequests['F-Output-1'].setValues(
                variables=('S', 'LE', 'U', 'COORD', 'EVOL', 'RF'),
                frequency=1,
            )
        except Exception as _fbe:
            emit("[LF] !!! F-Output-1 폴백도 실패: %s" % _fbe)
    #   [2026-10-21] 덱에 실제로 E 가 들어갔는지 로그로 남긴다. E 가 없으면 ZK/zk_Iw 가
    #     `E_field_missing` 으로 실패하므로, 재사용/stale ODB 와 출력요청 누락을 구분한다.
    for _fon in ('F-Output-1', 'F-ZK'):
        try:
            _fo = my_model.fieldOutputRequests[_fon]
            emit("[LF] field output %s: variables=%s frequency=%s"
                 % (_fon, tuple(_fo.variables), _fo.frequency))
        except Exception as _foe:
            emit("[LF] field output %s 확인 실패: %s" % (_fon, _foe))

    #   [§19.66] 같은 설계점의 LF 는 재사용한다(잡 재제출 ~3분 절약).
    #   [paper-match] 단, 모델 정의가 바뀐 뒤의 **stale LF ODB** 는 재사용하지 않는다.
    #     sidecar(.model.txt)의 tag 가 현재 경로와 일치할 때만 재사용한다.
    if os.path.exists(LF_ODB) and lf_odb_is_current(LF_ODB, tag=_LF_MODEL_TAG,
                                                    model_fp=_LF_MODEL_FP):
        emit("[LF] 재사용: %s (설계점 %s, sidecar tag=%s) - 잡 재제출 생략"
             % (LF_ODB, point_tag(x_c, d_c), _LF_MODEL_TAG))
    else:
        if os.path.exists(LF_ODB):
            emit("[LF] !!! stale LF ODB (sidecar 없음/불일치) -> 재제출: %s" % LF_ODB)
        run_job_safely(LF_JOB)
        try:
            lf_sidecar_write(LF_ODB, tag=_LF_MODEL_TAG, extra={
                'x_c': x_c, 'd_c': d_c, 'load_mode': LOAD_MODE,
                'paper_path': int(_PAPER_PATH), 'paper_hold_u3': int(_LF_HOLD_U3_ALL),
                'paper_include_clamps': int(PAPER_INCLUDE_CLAMPS),
                'mesh_order': MESH_ORDER, 'thickness': THICKNESS,
                'cable_prestress': CABLE_PRESTRESS,
                'clamp_prestress': CLAMP_PRESTRESS,
                #   [CONFIG] 기계 판독용: 재사용 전에 대조하는 값들(sidecar extra 는 사람용).
                'model_fp': _LF_MODEL_FP, 'config_sha1': _CFG_SHA1,
                'seed_div': SEED_DIV, 'sigma0': SIGMA0,
                'lf_steps': ','.join([s for s in my_model.steps.keys() if s != 'Initial']),
            })
            emit("[LF] sidecar 기록: %s.model.txt (tag=%s)" % (LF_ODB, _LF_MODEL_TAG))
        except Exception as _sc_err:
            emit("[LF] !!! sidecar 기록 실패(무시): %s" % _sc_err)

    cmd = 'abaqus python "%s" %s LF' % (os.path.join(_HERE, "eval_abaqus.py"), LF_ODB)

elif fidelity == 'HF':

    #   [§19.87] 짝 LF 유무와 **무관하게** HF 를 수행한다(사용자 지시 — §19.79 의 중단 로직 삭제).
    #     짝 LF 가 있으면 eval 이 lf1~lf3 + hf 를 모두 내고, 없으면 HF 단독 평가(HFONLY)로 hf 만 낸다.
    #     주의: LF 가 없으면 그 점의 lf 지표는 nan 이므로 **MFBO 의 LF 자료로는 못 쓴다**
    #     (LF↔HF 상관 계산에는 LF 산출물이 필요하다). HF 단독은 완주·hf 확인용 경로다.
    _have_lf = os.path.exists(LF_ODB)
    if not _have_lf:
        emit("[eval] 짝 LF 없음(%s) -> HF 단독 평가(HFONLY). lf 지표는 nan 이 된다."
             % os.path.basename(LF_ODB))
    elif not lf_odb_is_current(LF_ODB, tag=_LF_MODEL_TAG):
        emit("[eval] !!! 짝 LF ODB sidecar 불일치(stale 가능) - lf/ZK 가 현재 모델 것이 아닐 수 있다. "
             "LF 먼저 재실행을 권장한다: %s" % os.path.basename(LF_ODB))

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
    #   기본 BASE_PROBE=True (live). False 로 바꾸면 이 런에서 자동 실행하지 않는다.
    #   소스는 **짝 LF ODB(u3=0)** 우선; 없으면 skip, 디버깅은 SOLARSAIL_BASE_ODB.
    #   수동 진단:
    #     abaqus python base_state_probe.py aba\LF_Analysis_xc0p5_dc0p5.odb Step-GlobalTension 1.0 1.0
    #     abaqus python base_state_probe.py results\odb\HF_xc0p5_dc0p5_u3free.odb Step-GlobalTension 1.0 1.0
    #     abaqus python base_state_probe.py aba\HF_Postbuckle_u3free.odb Step-Postbuckle 1.0 1.0
    #   고정 이름 Buckle_Analysis.odb 는 포인트 태그가 없어 per-point 진단에 쓰지 않는다.
    #   좌굴 성공/실패와 무관하게 ODB 의 GlobalTension 프레임에서 읽는다.
    try:
        if BASE_PROBE:
            _probe = os.path.join(_HERE, 'base_state_probe.py')
            #   [2026-10-21] HF 가 처음부터 u3-free 이므로 HF 의 GlobalTension 말단은 이미 주름진
            #     상태다. R-13 은 **짝 LF ODB(u3=0)를 첫 후보**로 읽는다. 고정 이름
            #     Buckle_Analysis.odb 는 포인트 태그가 없어 폴백으로 쓰지 않는다.
            #     디버깅은 SOLARSAIL_BASE_ODB 로만 명시한다. 없으면 skip.
            _bodb = None
            _bodb_src = ''
            _base_override = os.environ.get('SOLARSAIL_BASE_ODB')
            if _base_override:
                _bodb = os.path.abspath(_base_override)
                _bodb_src = 'explicit SOLARSAIL_BASE_ODB'
            elif os.path.exists(LF_ODB) and lf_odb_is_current(LF_ODB, tag=_LF_MODEL_TAG):
                _bodb = os.path.abspath(LF_ODB)
                _bodb_src = '짝 LF ODB (same x_c,d_c; u3=0)'
            if _bodb is None or not os.path.exists(_bodb):
                emit("[R-13] skip: no point-consistent base ODB "
                     "(override missing or LF sidecar missing/mismatch)")
            elif os.path.exists(_probe):
                emit("[R-13] base state ODB: %s  [%s]" % (_bodb, _bodb_src))
                if _bodb_src.startswith('짝 LF'):
                    #   LF 는 HF 선행 런이므로 mtime < run start 가 정상. sidecar tag 가 신선도 증거다.
                    emit("[R-13] source=LF prerequisite; mtime<run start is normal; sidecar tag=%s"
                         % _LF_MODEL_TAG)
                else:
                    try:
                        _mt = os.path.getmtime(_bodb)
                        _age = _RUN_T0 - _mt
                        _mt_txt = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(_mt))
                        emit("[R-13] explicit override mtime=%s (age %.1f s; stale check not enforced)"
                             % (_mt_txt, _age))
                    except Exception as _mt_err:
                        emit("[R-13] ODB mtime 확인 실패(무시): %s" % _mt_err)
                _pcmd = ('abaqus python "%s" "%s" Step-GlobalTension'
                         % (_probe, _bodb))
                _pr = subprocess.run(_pcmd, shell=True, capture_output=True, text=True)
                # probe 는 40여 줄을 쏟는다(요소 좌표·단면별 통계·RF 상위 8개·인스턴스별 max|u3|).
                # 로그에 필요한 것은 ">>>" 가 붙은 요약 + [§19.94] 신원/경고/영역 개수 줄이다.
                _pr_all = ((_pr.stdout or "") + (_pr.stderr or "")).splitlines()
                _n_kept = 0
                #   [§19.94] '>>>' 요약만 통과시키면 **새로 추가한 ODB 신원/경고 줄이 통째로 사라진다**
                #     (신원 로그가 hf_run_log.txt 에 없으면 stale 판정 자체가 불가능하다). 그래서
                #     신원·경고·영역 개수 줄도 함께 남긴다.
                #   [§19.97] 스텝 접미사([R-13:CT] / [R-13:GT])를 **떼고** 판정하는 공용 필터로 교체:
                #     `_pass_pf` 를 스텝마다 늘리는 방식은 스텝이 추가될 때마다 고쳐야 하고, 한 곳을
                #     빠뜨리면 CT 판독이 로그에서 통째로 사라진다(§19.94 와 같은 유형의 사고).
                for _ln in _pr_all:
                    _ls = _ln.strip()
                    if _probe_keep_line(_ls):
                        emit("   %s" % _ls); _n_kept += 1
                emit("[R-13] 요약 %d줄 / 원본 %d줄 (상세 생략)" % (_n_kept, len(_pr_all)))
            else:
                emit("[R-13] probe missing: %s" % _probe)
    except Exception as _probe_err:
        emit("[R-13] base state 측정 실패(무시): %s" % _probe_err)

    # ============================================================================================
    # [§19.97] R-13 2차 판독: `Step-ClampTension` 말단 = 꼭짓점 당김 + 클램프 당김이 **모두** 가해진
    #   후좌굴 직전 = **이 설계점의 실제 프리텐션**(§19.96 정정: GlobalTension 말단은 클램프 당김이
    #   빠진 상태라 '실제 프리텐션'이 아니다). 위 1차 판독(GlobalTension)은 그대로 둔다.
    #   배경 질문: 이 상태의 중앙 응력 Pa 스케일이 7000 Pa 운용점과 얼마나 다른가.
    #
    #   [2026-10-21] 이 시점에 이 런의 HF ODB 는 아직 없다 - 포스트버클 잡은 아래에서
    #     제출된다. HF 모델은 포스트버클 제출 **직전에 Step-ClampTension 을 삭제**하므로
    #     HF ODB 에는 그 스텝이 **영원히 없다**. -> 짝 LF ODB(또는 명시 override)만 읽는다.
    #   우선순위: (1) 짝 LF ODB(sidecar 일치, u3=0) -> (2) SOLARSAIL_BASE_ODB override -> (3) skip.
    #     고정 Buckle_Analysis.odb 폴백은 제거했다.
    try:
        if BASE_PROBE:
            _probe2 = os.path.join(_HERE, 'base_state_probe.py')
            _ct_odb = None
            _ct_src = ''
            _base_override = os.environ.get('SOLARSAIL_BASE_ODB')
            if _base_override:
                if os.path.exists(_base_override):
                    _ct_odb = os.path.abspath(_base_override)
                    _ct_src = 'explicit SOLARSAIL_BASE_ODB'
                else:
                    _ct_odb = None
                    _ct_src = 'override_missing'
            elif os.path.exists(LF_ODB) and lf_odb_is_current(LF_ODB, tag=_LF_MODEL_TAG):
                _ct_odb = os.path.abspath(LF_ODB)
                _ct_src = '짝 LF ODB (same x_c,d_c)'
            if _ct_odb is None:
                emit("[R-13:CT] skip: no point-consistent base ODB "
                     "(override missing or LF sidecar missing/mismatch)")
                _R13_CT['r13_ct_note'] = 'no point-consistent base ODB'
            elif os.path.exists(_probe2):
                #   [paper-match] paper LF ODB 에는 Step-ClampTension 이 없다. 그 경로의
                #     프리스트레스 상태는 Step-GlobalTension 이므로 그 스텝을 읽는다.
                _r13_step = ('Step-GlobalTension'
                             if (LOAD_MODE == 'cable_stress' and _PAPER_PATH)
                             else 'Step-ClampTension')
                emit("[R-13:CT] --- 프리스트레스 말단 판독 (step=%s) ---" % _r13_step)
                emit("[R-13:CT] ODB 선택: %s  [%s]" % (_ct_odb, _ct_src))
                emit("[R-13:CT] !!! 이 런의 HF ODB(%s)는 아직 없다(포스트버클 제출 전) -> 위 파일로 "
                     "**대체 판독**했다. 이 숫자를 '이 런의 HF 값'으로 읽지 마라." % HF_ODB)
                if not _ct_src.startswith('짝 LF'):
                    emit("[R-13:CT]     (주의) explicit override 사용 -> 다른 설계점/다른 런의 "
                         "산출물일 수 있다(아래 mtime 확인).")
                #   [§19.94] 신선도 판정 - mtime 이 런 시작보다 오래되면 이 런이 만든 파일이 아니다.
                _mt2_txt = 'nan'
                _stale2 = 'nan'
                try:
                    _mt2 = os.path.getmtime(_ct_odb)
                    _age2 = _RUN_T0 - _mt2
                    _mt2_txt = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(_mt2))
                    if _ct_src.startswith('짝 LF'):
                        #   LF 는 HF 선행 런이므로 mtime<run start 가 정상. sidecar tag 가 증거다.
                        _stale2 = '0'
                        emit("[R-13:CT] source=LF prerequisite; sidecar tag=%s; "
                             "mtime<run start is normal (mtime %s)"
                             % (_LF_MODEL_TAG, _mt2_txt))
                    elif _age2 > 0.0:
                        _stale2 = '1'
                        emit("[R-13:CT] !!! explicit override is older than run start: "
                             "mtime(%s) vs run(%s), age %.1f s -> 이 런의 산출물이 아니다."
                             % (_mt2_txt, time.strftime('%Y-%m-%d %H:%M:%S',
                                                        time.localtime(_RUN_T0)), _age2))
                    else:
                        _stale2 = '0'
                        emit("[R-13:CT] 신선도 OK: 런 시작 후 %.1f 초에 생성됨 (mtime %s)"
                             % (-_age2, _mt2_txt))
                except Exception as _mt2_err:
                    emit("[R-13:CT] mtime 확인 실패(무시): %s" % _mt2_err)
                #   사후(authoritative) 재판독 안내. **스텝을 실제로 보유한 파일**로 안내한다 -
                #     HF_Postbuckle.odb 를 안내하면 '스텝 없음' 경고만 돌아온다(위 **주의** 참조).
                emit("[R-13:CT] 사후 재판독: abaqus python \"%s\" \"%s\" %s"
                     % (_probe2, _ct_odb, _r13_step))
                emit("[R-13:CT] (한계) 이 런의 HF ODB(%s)에는 이 프리스트레스 스텝이 없다 - "
                     "HF 런이 포스트버클 직전에 ClampTension 을 삭제하고 GlobalTension 에서 시작한다. "
                     "프리스트레스 확정은 위 ODB(짝 LF/override)로 한다(이 진단은 전용)." % HF_ODB)
                _pcmd2 = 'abaqus python "%s" "%s" %s' % (_probe2, _ct_odb, _r13_step)
                _pr2 = subprocess.run(_pcmd2, shell=True, capture_output=True, text=True)
                _pr2_all = ((_pr2.stdout or "") + (_pr2.stderr or "")).splitlines()
                _n_kept2 = 0
                _kv2 = {}
                for _ln in _pr2_all:
                    _ls = _ln.strip()
                    if not _ls:
                        continue
                    if 'KV:' in _ls:            # 프로브가 찍는 기계 판독용 1줄
                        for _tok in _ls.split('KV:', 1)[1].split():
                            if '=' in _tok:
                                _k2, _v2 = _tok.split('=', 1)
                                _kv2[_k2] = _v2
                    if _probe_keep_line(_ls):
                        emit("   %s" % _ls); _n_kept2 += 1
                emit("[R-13:CT] 요약 %d줄 / 원본 %d줄 (rc=%d, 상세 생략)"
                     % (_n_kept2, len(_pr2_all), _pr2.returncode))
                _R13_CT.update({
                    'r13_ct_centre_mean_pa': _kv2.get('centre_mean_pa', 'nan'),
                    'r13_ct_centre_minp_pa': _kv2.get('centre_minp_pa', 'nan'),
                    'r13_ct_area_frac_neg': _kv2.get('area_frac_neg', 'nan'),
                    'r13_ct_minp_frac_neg': _kv2.get('minp_frac_neg', 'nan'),
                    'r13_ct_uz_max_m': _kv2.get('uz_max_m', 'nan'),
                    'r13_ct_uz_ave_m': _kv2.get('uz_ave_m', 'nan'),
                    'r13_ct_odb': _ct_odb,
                    'r13_ct_odb_mtime': _mt2_txt,
                    'r13_ct_stale': _stale2,
                    'r13_ct_note': '%s / %s%s' % (
                        _ct_src,
                        ('sidecar-validated LF (mtime<run start is normal)'
                         if _ct_src.startswith('짝 LF') else
                         ('stale(런 시작 이전 파일)' if _stale2 == '1' else
                          ('fresh(런 시작 이후 생성)' if _stale2 == '0' else 'mtime 미확인'))),
                        '' if _kv2 else ' / KV 줄 없음: 스텝 없음 또는 S/프레임 없음(프로브 rc=%d)'
                        % _pr2.returncode)})
                if not _kv2:
                    emit("[R-13:CT] !!! KV 줄이 없다 -> 판독 실패(스텝 없음/프레임 0/S 없음). "
                         "점별 txt 는 nan 으로 남기고 이유는 r13_ct_note 에 남긴다.")
            else:
                emit("[R-13:CT] probe missing: %s" % _probe2)
                _R13_CT['r13_ct_note'] = 'base_state_probe.py missing'
    except Exception as _probe2_err:
        emit("[R-13:CT] 판독 실패(무시): %s" % _probe2_err)
        #   예외로 끝나도 이유를 남긴다(기본 note 'probe 미실행' 은 이 경우 사실이 아니다).
        _R13_CT.setdefault('r13_ct_note', 'CT 판독 예외: %s' % _probe2_err)

    # 기존 Step 정리: Post-buckling은 GlobalTension 직후에서 시작하며,
    # 중간 단계(ClampTension)를 건너뛰고 바로 최종 하중으로 Ramping함
    if 'Step-Buckle' in my_model.steps: del my_model.steps['Step-Buckle']
    if 'Step-ClampTension' in my_model.steps: del my_model.steps['Step-ClampTension']

    # Static, General 포스트버클링 수행    
    my_model.StaticStep(
        name='Step-Postbuckle', 
        previous='Step-GlobalTension',  
        nlgeom=ON, 
        stabilizationMagnitude=STABILIZE_POST,      # paper 2e-4 / physical 5e-4 / old disp 0.01
        #   근거: 실측 ALLSD/ALLIE = 0.50%% 인데 allsdtol = 5%% -> 인공감쇠가 약 10배의 여유를 남기고
        #   거의 작동하지 않았다. 즉 "감쇠가 컷백을 만든다"는 과거 반증은 "감쇠를 늘리면 도움이 될 것"을
        #   반증하지 않는다(다른 주장). allsdtol 이 상한을 걸어 자체 제한되므로 위험은 유계다.
        #   판정: 초기 100~200 증분의 증분 상한이 3e-4 이상으로 회복되는가 (10~20분) + energy_check.py
        stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
        continueDampingFactors=False,
        adaptiveDampingRatio=0.002,
        initialInc=1e-4,
        minInc=1e-8,
        maxInc=0.1,
        maxNumInc=1000
        # 상한 상향은 해결책이 아니라 연장이다 — 실제 원인은 국소 불안정으로 증분이 눌리는 것(§13).
        #   ALLSD/ALLIE 가 0.50%(허용 5%)라 allsdtol 상향은 binding 이 아니다(폐기된 레버).
    )
    my_model.keywordBlock.synchVersions(storeNodesAndElements=False)

    #   paper HF 는 BC_Stabilize_Z 가 없다(u3-free). old disp HF 만 여기서 전막 u3 를 해제한다.
    if 'BC_Stabilize_Z' in my_model.boundaryConditions:
        my_model.boundaryConditions['BC_Stabilize_Z'].deactivate('Step-Postbuckle')

    # Step-Buckle 삭제 시 BC_Edges_Only_Z(prescribed condition)가 함께 삭제되므로
    # Postbuckle 스텝에 모서리 z 구속을 재생성한다 (Galhofo 참조: 3개 모서리 u3=0).
    # 셋은 위에서 만든 All_Edges_NoClamp(클램프 구간 제외)를 그대로 쓴다.
    #   [§19.111] RELEASE_EDGE_U3=True 면 재생성하지 않는다: EdgeZ_Spring 은 model-level
    #     engineeringFeature 라 스텝 삭제와 무관하게 모델에 남아 있고, 같은 이름으로 다시 만들면
    #     중복 이름 오류다. 3점 RP u3(BC_*_Corner_Z, Initial)도 그대로 살아 있다.
    #   [paper-match] _PAPER_HF 면 Initial 에서 만든 BC 가 모든 스텝에 살아 있으므로
    #     여기서 다시 만들면 중복 이름 오류다.
    if not RELEASE_EDGE_U3 and not _PAPER_HF:
        my_model.DisplacementBC(
            name='BC_Edges_Only_Z',
            createStepName='Step-Postbuckle',
            region=a.sets['All_Edges_NoClamp'],
            u3=0
        )

    #   old disp: 최종 장력(목표 7000 Pa, 미검증)까지 증가.
    #   paper: _DISP_FINAL=0 -> 하중 변화가 없는 **정착 스텝**이다. 주름은 GlobalTension 에서
    #     이미 발달했고, get_hf/DN 은 이 Postbuckle 마지막 프레임을 읽는다.
    for name, sign in [('Disp_Control_Right', 1), ('Disp_Control_Left', -1)]:
        my_model.boundaryConditions[name].setValuesInStep(
            stepName='Step-Postbuckle',
            u1=sign * _DISP_FINAL * cos_val,
            u2=-_DISP_FINAL * sin_val
        )

    if PULL_TOP:
        #   [§19.105] 포스트버클 최종 하중도 세 꼭짓점 동일 - Top 은 +y 방향 GLOBAL_FINAL.
        my_model.boundaryConditions['Disp_Control_Top'].setValuesInStep(
            stepName='Step-Postbuckle',
            u1=_DISP_FINAL * top_cos * TOP_PULL_SCALE,
            u2=-_DISP_FINAL * top_sin * TOP_PULL_SCALE
        )

    #   [2026-10-20] cable_stress 경로는 클램프 초기응력 -> Postbuckle 변위 없음.
    #     종전 MFBO 기준선(disp)만 CLAMP_FINAL 변위를 유지한다.
    if LOAD_MODE == 'disp':
        for name, sign in ([] if NO_CLAMP else [('Disp_Control_CL', -1), ('Disp_Control_CR', 1)]):
            my_model.boundaryConditions[name].setValuesInStep(
                stepName='Step-Postbuckle',
                u1=sign * CLAMP_FINAL / sqrt2,
                u2=CLAMP_FINAL / sqrt2)
    

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

    if not _HF_USE_POSTBUCKLE:
        #   paper HF: no-load 정착 스텝을 모델에서 제거한다. F-Output-1/임퍼펙션은 이미
        #   GlobalTension에 적용되어 있으므로 hf/DN을 그 마지막 프레임에서 읽을 수 있다.
        if 'Step-Postbuckle' in my_model.steps:
            del my_model.steps['Step-Postbuckle']
            emit("[HF] paper no-load Step-Postbuckle 제거 -> GlobalTension에서 hf/DN")

    if not run_job_with_ladder(HF_JOB):
        emit("!!! ERROR: HF 완주 실패(안정화 사다리 소진) -> RESULTS:FAIL")
        _archive_result(x_c, d_c, 'HF_FAIL', 'nan,nan,nan,nan',
                        dict({'status': 'LADDER_EXHAUSTED',
                              'stabilize_global_ladder': str(STABILIZE_GLOBAL_LADDER),
                              'stabilize_post_ladder': str(STABILIZE_POST_LADDER)},
                             **_r13_ct_keys()))
        emit("RESULTS:FAIL")
        sys.stdout.flush()
        sys.exit(1)


    if _have_lf:
        cmd = 'abaqus python "%s" %s %s HF' % (os.path.join(_HERE, "eval_abaqus.py"),
                                               LF_ODB, HF_ODB)
    else:
        cmd = 'abaqus python "%s" %s HFONLY' % (os.path.join(_HERE, "eval_abaqus.py"),
                                                HF_ODB)

try:
    #   [paper-match] ZK 도 stale 파일을 남기지 않는다(extraction.txt 와 같은 처리).
    for _stale in ('extraction.txt', 'zk_metrics.txt'):
        try:
            if os.path.exists(_stale):
                os.remove(_stale)
        except Exception:
            pass
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    rc = p.returncode
    _eval_bad = (rc != 0)
    #   [19.106] DN: 줄 수집. 아래 FAIL 검사가 'DN:FAIL reason=...' 의 부분문자열 'FAIL' 을
    #     잡아 성공한 HF 런을 실패로 뒤집기 때문에, DN 줄은 **먼저** 걷어내고 continue 한다.
    #     (RESULTS: 4필드 계약과 무관한 진단 줄이다.)
    _dn_lines = []
    _zk_lines = []
    _zkiw_lines = []
    for _ln in ((p.stdout or "") + (p.stderr or "")).splitlines():
        _ls = _ln.strip()
        if not _ls:
            continue
        if _ls.startswith("ZKIW:"):
            #   `zk_Iw` sidecar 진단 줄. FAIL 문자열이 있어도 정상 런을 뒤집지 않는다.
            _zkiw_lines.append(_ls)
            emit("   [ZKIW] %s" % _ls)
            continue
        if _ls.startswith("DN:"):
            _dn_lines.append(_ls)
            continue
        if _ls.startswith("ZK:"):
            #   [paper-match] ZK 도 DN 과 같은 이유로 FAIL 검사보다 먼저 걷어낸다.
            #     'ZK:FAIL reason=...' 의 FAIL 이 정상 LF/HF 런을 실패로 뒤집으면 안 된다.
            _zk_lines.append(_ls)
            emit("   [ZK] %s" % _ls)
            continue
        if _ls.startswith("RESULTS:FAIL") or "Target step not found" in _ls or "FAIL" in _ls:
            _eval_bad = True
            emit("   [eval] %s" % _ls)
            continue
        if _ls.startswith("RESULTS:"):
            continue
        if _ls.startswith("!!!") or "ERROR" in _ls or "Traceback" in _ls:
            emit("   [eval] %s" % _ls)
    if rc != 0:
        emit("!!! ERROR: eval_abaqus.py exited with code %d" % rc)

    # 추출 스크립트가 출력한 "RESULTS:..." 라인을 찾아 전달
    if _eval_bad:
        emit("!!! ERROR: eval_abaqus.py 가 실패했습니다(LF 스텝/입력 확인) -> RESULTS:FAIL")
        emit("RESULTS:FAIL")
    elif os.path.exists('extraction.txt'):
        with open('extraction.txt', 'r', encoding='utf-8', errors='replace') as f:
            _rc = f.read().strip()
        emit("RESULTS:" + _rc)
        _extra = {}
        try:
            for _fn in (('stabilize_used.txt', 'flip_bias.txt')
                        if fidelity == 'HF' else ()):
                if os.path.exists(_fn):
                    with open(_fn, 'r', encoding='utf-8', errors='replace') as _sf:
                        for _ln in _sf:
                            if '=' in _ln:
                                _k2, _v2 = _ln.strip().split('=', 1)
                                _extra[_k2] = _v2
        except Exception:
            pass
        _extra.update(_r13_ct_keys())
        _archive_result(x_c, d_c, fidelity, _rc, _extra)
        # [19.106] DN: 줄 -> 점별 txt 에 dn_ 접두 append (_archive_result 의 'w' 뒤).
        #     RESULTS 4필드와 all_results.csv 열 스키마는 그대로다.
        _archive_dn(x_c, d_c, _dn_lines)
        # [paper-match] ZK: 줄 -> 점별 txt 에 zk_ 접두 append (CSV 열은 그대로).
        #   HFONLY 처럼 짝 LF 가 없으면 ZK 줄 자체가 없다(경고 생략).
        if _zk_lines or fidelity == 'LF':
            _archive_zk(x_c, d_c, _zk_lines, _zkiw_lines)
        # 성공한 HF 런에서만: HF ODB 설계점 태깅 + slope_stats 아카이빙(예외 안전).
        #     _archive_result 가 점별 txt 를 'w' 로 다시 쓰므로 **그 뒤에** append 해야 한다.
        _archive_slope(x_c, d_c)
    else:
        emit("!!! ERROR: Extraction failed. 'extraction.txt' not found.") 
        emit("RESULTS:FAIL")   # 상위(mfbo)가 원인을 식별하도록 명시적 실패 신호

except Exception as err:
    emit("Error during data extraction: %s" % str(err))
    sys.stdout.flush()   # CAE noGUI: sys.exit 앞에서 버퍼 비움
    sys.exit(1)
