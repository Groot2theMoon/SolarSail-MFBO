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
         - 1차: 모드 소스는 '클램프 없는' 별도 모델의 좌굴 런이다(run_abaqus_cable.py 또는
                전용 모드 스크립트). 클램프가 있는 base state 는 선형 좌굴 스펙트럼을 주지
                못한다(실측: 음수고유값 598~2897 / CONVERGED=0 / EIGENVALUES CANNOT BE FOUND).
         - 2차: 그 .fil 을 IMPERFECTION_NAME 으로 스테이징해 초기 결함으로 주입한다.
                aba_imperfection.py 가 모드 개수를 검증하고, 부족하면 HF 제출 전에 중단한다.
                자기 좌굴 잡(Buckle_Analysis)은 클램프 base state 증거용으로만 유지한다.
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
    print("!!! WARNING: eval_abaqus.py 위치를 찾지 못했습니다. cwd=%s" % os.getcwd())
    return os.getcwd()

_HERE = _resolve_here()
# Abaqus 작업 디렉터리(=산출물 위치) = code/aba.
# ---- 코드 상수 (2026-09-21: 환경변수 전부 제거, 값은 이 파일에 고정) ----
# 산출물 디렉터리 이름. mfbo.py 의 RUN_DIR_NAME 과 반드시 같은 값이어야 한다.
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
                              build_perturbation, perturbation_report, mode_table_info)
print("[run_abaqus] _HERE = %s" % _HERE)
print("[run_abaqus] _RUN  = %s" % _RUN)
 
print("DEBUG: All sys.argv: " + str(sys.argv))

try:
    # abaqus cae noGUI=run_abaqus.py -- [HF/LF] x_c d_c
    fidelity = sys.argv[-3].upper()
    x_c = float(sys.argv[-2])
    d_c = float(sys.argv[-1])
    
except:
    print("Error: Invalid arguments. Usage: abaqus cae noGUI=run_abaqus.py -- HF x_c d_c")
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
            print("Removing stale artifact: %s" % _f)
            try:
                os.remove(_f)
            except OSError as _e:
                print("  (warning) could not remove %s: %s" % (_f, _e))

    lck_file = job_name + '.lck'
    odb_file = job_name + '.odb'
    # 기존 Lock 파일이 있다면 삭제 시도 < 이전 실행 강제 종료 시 남게 됨
    if os.path.exists(lck_file):
        print("Detected old lock file: %s. Removing it..." % lck_file)
        try:
            os.remove(lck_file)
        except OSError:
            # 만약 삭제가 안 된다면 다른 프로세스가 실제로 사용 중
            print("!!! FATAL ERROR: Cannot remove lock file. Is another Abaqus process running?")
            sys.exit(1)

    if job_name in mdb.jobs:
        del mdb.jobs[job_name]
    
    # 병렬 실행 (속도). 기본 1 = 기존과 완전히 동일한 거동.
    #   단일 CPU 로 373 증분/30분 수준이면 4 스레드로 대략 1.5~2.5배 단축 여지가 있다.
    #   주의: 라이선스 토큰이 부족하면 잡이 라이선스 오류로 죽는다 -> 그때는 1 로 되돌린다.
    #   라이선스 토큰 오류가 나면 위 NUMCPUS 상수를 1 로 바꾼다.
    _ncp = NUMCPUS
    print("[run_abaqus] numCpus=%d numDomains=%d (코드 상수 NUMCPUS)" % (_ncp, _ncp))
    job = mdb.Job(name=job_name, model=model_name, numCpus=_ncp, numDomains=_ncp)
    print("Submitting Job: %s" % job_name)
    job.writeInput(consistencyChecking=OFF)
    job.submit(consistencyChecking=OFF)
    
    job.waitForCompletion()

    time.sleep(1.0)

    # 잡 결과 핵심 줄을 콘솔에 직접 찍는다 (로그 일부만 붙여넣어도 원인 판별 가능)
    print_job_diag(job_name)
    
    # ABORTED가 아니면서, ODB 파일이 실제로 존재하면 성공으로 간주
    if job.status == ABORTED or not os.path.exists(odb_file):
        print("!!! ERROR: Job %s failed. Actual Status: %s" % (job_name, str(job.status)))
        sys.exit(1)
        
    if not job_completed_ok(job_name):
        print("!!! WARNING: %s — .sta/.msg 에 'HAS COMPLETED SUCCESSFULLY' 없음 "
              "(중도 중단 의심; odb 존재만으로는 판정 불가 — R-12)" % job_name)
    print("Job %s completed successfully (Status: %s)." % (job_name, str(job.status)))
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
            '***ERROR')
    for ext in ('msg', 'dat'):
        fn = '%s.%s' % (job_name, ext)
        if not os.path.exists(fn):
            print("[DIAG:%s] %s 없음" % (job_name, fn))
            continue
        size = os.path.getsize(fn)
        with open(fn, 'r') as f:
            if size > 2000000:
                f.seek(size - 2000000)
            text = f.read()
        hits = [ln.strip() for ln in text.splitlines()
                if ln.strip() and any(k.lower() in ln.lower() for k in keys)]
        print("[DIAG:%s] %s (%.0f KB) 핵심줄 %d개"
              % (job_name, fn, size / 1024.0, len(hits)))
        for ln in hits[-5:]:
            print("      | %s" % ln[:150])

sqrt2 = 1.414
N_EIG = 4          # 임퍼펙션에 쓸 좌굴모드 수 (*IMPERFECTION / *NODE FILE)

# ---- 2-모델 레시피: 임퍼펙션(좌굴모드) 소스 분리 (2026-09-22) ----------------
#   모드 추출 = run_abaqus_cable.py (클램프 없음 -> 실측상 항상 성공, 모드 4개)
#     abaqus cae noGUI=run_abaqus_cable.py          (code\ 에 Buckle_Analysis.fil 생성)
#     abaqus cae noGUI=run_abaqus.py -- HF <x_c> <d_c>   (이 파일을 스테이징해 소비)
#   근거: 클램프 패치는 기존 노드에 Coupling 만 걸어 메쉬를 바꾸지 않으므로 두 모델의
#         막 노드 좌표/라벨이 동일하다 -> *IMPERFECTION 의 노드 라벨 매핑이 성립한다.
#         (클램프가 있는 base state 는 음수 고유값 598~2897/CONVERGED=0 -> 모드 추출 불가)
#         상세 근거·실패 이력: aba_imperfection.py, references/buckle-failure-triage.md
MODE_SOURCE = 'external'     # 'external' = 클램프 없는 외부 .fil(run_abaqus_mode.py) | 'self' = 자기 좌굴 잡
#   모드 소스 = run_abaqus_mode.py (클램프 없음, HF 파라미터 정렬). 산출물은 code\ 에 쌓인다.
#   (이전 소스였던 run_abaqus_cable.py 의 Buckle_Analysis.fil 을 쓰려면 아래 이름만 교체한다.)
MODE_SOURCE_FIL = os.path.join('..', 'ClampFree_Buckle.fil')   # code\aba -> code\
MODE_SOURCE_DAT = os.path.join('..', 'ClampFree_Buckle.dat')
MODE_SOURCE_MSG = os.path.join('..', 'ClampFree_Buckle.msg')
MODE_SOURCE_STEP = 2         # 소스 .fil 안의 스텝 번호 (케이블 런: 1=GlobalTension 2=Buckle)
IMPERFECTION_NAME = 'ClampFree_Buckle'   # 모드 소스 잡 이름과 '같은' 이름 -> 원장 C-1(조용한 0-모드 소비) 차단
                                         # (MODE_SOURCE_FIL 의 basename 과 함께 움직여야 한다)
IMPERFECTION_MODES = (1, 2, 3, 4)
IMPERFECTION_AMPL_T = 0.10   # 막 두께 배수(Galhofo 채택값 0.10 t). 진폭 민감도 = 0.50 으로 바꿔 재실행
RUN_SELF_BUCKLE_JOB = True   # 자기(클램프) 좌굴 잡도 계속 돌린다 -> 클램프 base state probe
                             # (코너 당김 하중분담) + 음수고유값 진단 증거를 함께 얻는다.
                             # 이 잡이 sys.exit(1) 로 스크립트를 끊으면 False 로 두고,
                             # 그때는 base state probe 가 '건너뜀' 으로 출력된다.
if MODE_SOURCE not in ('external', 'self'):
    raise RuntimeError("MODE_SOURCE 는 'external' 또는 'self' 여야 합니다 (현재 %r)" % (MODE_SOURCE,))

# ---- 임퍼펙션 주입 방식 (2026-09-22) --------------------------------------
#   모드 소스의 스텝 타입이 선택을 결정한다:
#     좌굴모드를 '선형 좌굴해석(*BUCKLE)' 으로 찾는 경우 -> 'odb_table'  [기본, 정본]
#        *BUCKLE 은 .fil 출력이 금지된다(실측: ClampFree_Buckle.dat:7025
#          "FILE OUTPUT IS NOT AVAILABLE FOR BUCKLING ANALYSIS"). 모드 프레임은 ODB 에만 있다.
#          -> abaqus python aba_mode_from_odb.py <job>.odb <step> modes_ClampFree_Buckle.txt 4
#          -> 그 표를 노드 좌표 섭동으로 주입(셸에서 *IMPERFECTION 과 등가).
#     '진동 고유모드(*FREQUENCY)' 를 모드 소스로 쓰는 경우 -> 'file'
#        .fil 에 모드가 기록되므로 스테이징 + *IMPERFECTION, FILE=, STEP=n (사용자 원본 08d4cbc 방식).
IMPERFECTION_MODE = 'odb_table'
MODE_TABLE = os.path.join('..', 'modes_ClampFree_Buckle.txt')   # code\aba -> code\
#   (출처 ODB/스텝은 상수로 중복 기재하지 않고 모드표 헤더에서 읽어 로그에 남긴다)
if IMPERFECTION_MODE not in ('odb_table', 'file'):
    raise RuntimeError("IMPERFECTION_MODE 는 'odb_table' 또는 'file' 여야 합니다 (현재 %r)"
                       % (IMPERFECTION_MODE,))

# 모드표는 모델을 만들기 전에 읽는다(없으면 HF 잡을 제출하지 않고 즉시 중단).
_pert = None
_mode_table = None
_labels = ()
if IMPERFECTION_MODE == 'odb_table':
    try:
        _mode_table = load_mode_table(MODE_TABLE, IMPERFECTION_MODES)
    except ImperfectionSourceError as _imp_err_tab:
        print("!!! ERROR: 임퍼펙션 모드표를 읽지 못했습니다 -> HF 잡을 제출하지 않고 중단합니다.")
        print(str(_imp_err_tab))
        sys.exit(1)
    _tinfo = mode_table_info(MODE_TABLE)
    print("[IMPERFECTION] 모드표 %s 로드 완료: 모드 %s / 모드당 노드 %d개"
          % (MODE_TABLE, sorted(_mode_table.keys()), len(_mode_table[IMPERFECTION_MODES[0]])))
    print("[IMPERFECTION] 모드표 출처: source=%s / step=%s / instance=%s"
          % (_tinfo.get('source', '?'), _tinfo.get('step', '?'), _tinfo.get('instance', '?')))
# ---- 2-모델 레시피 끝 ----------------------------------------------------
# A: 추출 요청 고유값 수 (음수모드 우회; run_abaqus_cable 과 동일)
#   base state 가 부정정이면 요청 개수를 줄이는 것이 subspace 수렴에 유리하다.
#   (구 MFBO_N_EIG_BUCKLE 환경변수 스윕은 2026-09-21 제거됨 -> 값을 바꾸려면 이 상수를 직접 수정)
N_EIG_BUCKLE = 100
# subspace 반복의 기저 벡터 수 (구 MFBO_VECTORS 환경변수는 2026-09-21 제거됨)
BUCKLE_VECTORS = 250

MODEL_NAME = 'SailModel_Triangle'
INSTANCE_NAME = 'MEMBRANE-1'

BASE = 20.0   # m
HEIGHT = 10.0 # m
THICKNESS = 5.0e-6 # F2: 2.5e-6 -> 5.0e-6 (cable 변형, 2.5um는 수렴 매우 어려움)
TARGET_STRESS = 7000.0 # Pa   # R-13: 목표 운용점 - 실제 도달 응력 미검증(측정 필요)

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
# R-13 실측(2026-09-21): PRETENSION_SCALE=10 -> 평균 면내응력 2122 Pa = 목표 7000 Pa 의 0.303배.
#   운용점을 목표에 맞추려면 약 33배(= DISP_GLOBAL 165um)가 필요하다.
#   단 PRETENSION_SCALE 는 최종 하중까지 함께 키우므로(포스트버클 변위 1mm -> 3.3mm),
#   운용점만 따로 맞추려면 DISP_GLOBAL / GLOBAL_FINAL 상수를 직접 수정한다(절대값).
DISP_GLOBAL = 0.000005 * PRETENSION_SCALE    # 운용점: 코너 당김 5e-5 m
CLAMP_PULL = DISP_GLOBAL * d_c
# 좌굴 스텝의 perturbation 변위 (K_delta 를 만드는 항).
#   Abaqus 문서 §6.2.3: 좌굴 스텝의 nonzero prescribed BC 는 '증분 응력'에 기여하고,
#   그 증분이 미분 초기응력 강성 K_delta 를 만든다. 크기 자체는 lambda 로 스케일되어
#   사라지지만(CONVERGED 수에는 영향 없음), K_delta 가 K0 대비 너무 작으면 고유값 분리가
#   나빠져 subspace 반복이 'EIGENVALUES CANNOT BE FOUND' 로 실패한다.
#   성공한 run_abaqus_cable.py 는 같은 솔버 설정(numEigen=100/SUBSPACE/vectors=250)에서
#   0.01 m 를 쓴다 -> 2026-09-22 부로 우리도 0.01 m 로 맞췄다(구 5e-4 = 1/20 이었다).
#   값을 바꾸려면 이 상수를 직접 수정한다(구 MFBO_PERT_MAG 환경변수는 제거됨).
PERTURBATION = 0.01
CLAMP_PERT = PERTURBATION * d_c
GLOBAL_FINAL = 0.0001 * PRETENSION_SCALE     # 최종 하중: 코너 당김 1e-3 m
CLAMP_FINAL = GLOBAL_FINAL * d_c


# 하중 각도 (28.6도)
angle_deg = 28.6
angle_rad = np.deg2rad(angle_deg)
cos_val = float(np.cos(angle_rad))
sin_val = float(np.sin(angle_rad))

# 모델 초기화 및 재질
if MODEL_NAME in mdb.models: del mdb.models[MODEL_NAME]
my_model = mdb.Model(name=MODEL_NAME)

# 멤브레인 (Kapton)
mat = my_model.Material(name='Kapton')
mat.Density(table=((1420.0, ),))
mat.Elastic(table=((2.5e9, 0.34),))
my_model.HomogeneousShellSection(name='Section-Membrane', material='Kapton', thickness=THICKNESS)

# 케이블 (Kevlar)
mat_cable = my_model.Material(name='Kevlar')
mat_cable.Density(table=((1440.0, ),)) 
mat_cable.Elastic(table=((62.0e9, 0.36),))
my_model.TrussSection(name='Section-Cable', material='Kevlar', area=CABLE_AREA)

# 파트 생성: 멤브레인
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
    지정된 좌표 기준 radius 내의 노드들을 묶어 강체운동을 하도록 Tie 설정
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
    
    # Coupling (RP <-> Membrane Nodes)
    my_model.Coupling(
        name=name+'_Coupling', controlPoint=rp_region, surface=patch_set, 
        influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
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

p.seedPart(size=BASE/200.0, deviationFactor=0.1) # 약 1만개
p.setMeshControls(regions=p.faces, elemShape=QUAD_DOMINATED, technique=FREE, algorithm=MEDIAL_AXIS)
# S4R-> S4 (cable 변형과 동일)
elemTypeQuad = ElemType(elemCode=S4, elemLibrary=STANDARD)
elemTypeTri = ElemType(elemCode=S3, elemLibrary=STANDARD)  
p.setElementType(regions=(p.faces,), elemTypes=(elemTypeQuad, elemTypeTri))
p.generateMesh()

# ---- 임퍼펙션(기하 섭동) 주입 (2026-09-22) --------------------------------
#   원본(업스트림) docstring 의 계획: "1차 해석 결과(.odb)에서 고유모드를 추출하여 초기 결함으로 주입".
#   셸 요소에서 *IMPERFECTION 은 결국 노드 좌표를 모드 형상만큼 옮기는 것이므로, 모드표
#   (aba_mode_from_odb.py)를 진폭 0.10 t 로 합산해 좌표를 직접 섭동한다(.fil/스텝타입 제약 우회).
#   위치: generateMesh 직후 + 어셈블리 regenerate 전 -> 의존 인스턴스가 이 좌표를 물려받는다.
if _pert is None:
    print("[IMPERFECTION] 기하 섭동 없음 (IMPERFECTION_MODE=%s)" % IMPERFECTION_MODE)
else:
    _amp = THICKNESS * IMPERFECTION_AMPL_T
    _pert = build_perturbation(_mode_table, _amp, IMPERFECTION_MODES)
    _labels = tuple(sorted(_pert.keys()))
    _seq = p.nodes.sequenceFromLabels(labels=_labels)
    _n = 0
    for _nd in _seq:
        _dz = _pert.get(_nd.label)
        if _dz is None:
            continue
        _cx, _cy, _cz = _nd.coordinates
        _nd.setValues(coordinates=(_cx, _cy, _cz + _dz))
        _n += 1
    print("[IMPERFECTION] 기하 섭동 적용: %d/%d 노드, %s (진폭 %.2f t = %.3e m)"
          % (_n, len(_labels), perturbation_report(_pert), IMPERFECTION_AMPL_T, _amp))
    if _n != len(_labels):
        print("!!! WARNING: 모드표 노드 %d개 중 %d개만 적용 -> 메쉬/라벨 불일치 의심"
              % (len(_labels), _n))
    # 말이 아니라 실측: 섭동이 실제 좌표에 들어갔는지 3개 노드를 찍어 로그에 남긴다.
    for _nd in p.nodes.sequenceFromLabels(
            labels=(_labels[0], _labels[len(_labels) // 2], _labels[-1])):
        print("               파트 노드 %d z=%.9e (Δz=%.3e m)"
              % (_nd.label, _nd.coordinates[2], _pert.get(_nd.label, 0.0)))
a.regenerate()

# 의존 인스턴스는 파트 메쉬를 공유한다 -> 어셈블리 쪽에서도 좌표가 같아야 한다(전달 확인).
if _pert is not None and _labels:
    _i0 = _labels[0]
    print("[IMPERFECTION] 어셈블리 인스턴스 확인: 노드 %d z=%.9e (파트와 같아야 함)"
          % (_i0, inst_memb.nodes.sequenceFromLabels(labels=(_i0,))[0].coordinates[2]))

# 꼭짓점 RP
rp1_obj, rp1_reg = create_rigid_patch('Top', V1, radius=0.2)
rp2_obj, rp2_reg = create_rigid_patch('Right', V2, radius=0.2)
rp3_obj, rp3_reg = create_rigid_patch('Left', V3, radius=0.2)

# 클램프 RP : 우측 빗변 중점 (15, 5), 좌측 빗변 중점 (5, 5)
# [대조 실험] 아래 NO_CLAMP 상수를 True 로 바꾸면 클램프(강체패치 + cable_CL/CR + BC)를
#   아예 만들지 않는다 (환경변수 아님 — 2026-09-21 부로 환경변수는 전부 제거됨).
#   run_abaqus_cable.py(성공)는 클램프가 없다. 우리만 클램프가 base state 하중의 33%를
#   받아 sigma2<0 영역(21.5%)을 만들고, 그 때문에 좌굴 고유값 추출이 실패한다는 가설을
#   클램프만 제거해 직접 검증한다.
NO_CLAMP = False                 # True = 클램프 생략 진단 모델 (run_abaqus_cable 대조용)
# 초기 가짜 응력(수렴 보조). 케이블 변형=700 Pa, 우리=500 Pa -> 정렬 노브
SIGMA0 = 500.0                   # 수렴 보조용 초기응력 [Pa]
if NO_CLAMP:
    print("[run_abaqus] NO_CLAMP=True : 클램프(cable_CL/CR + 강체패치 + BC) 없이 모델링 (대조 실험)")
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
print("[run_abaqus] buckle eigensolver=%s numEigen=%d vectors=%s"
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
print("[run_abaqus] steps=%s  _BUCKLE_STEP_NO=%d" % (_steps_in_order, _BUCKLE_STEP_NO))

# *IMPERFECTION, FILE= 은 results file(.fil) 을 읽는다 -> 좌굴 모드를 .fil 에 기록해야 임퍼펙션이 실제로 주입된다 (미요청 시 조용히 무시됨)

my_model.Stress(
    name='Initial_Stiffness',
    region=inst_memb.sets['All'],
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

all_edges = inst_memb.edges
a.Set(name='All_Edges', edges=all_edges)

my_model.DisplacementBC(
    name='BC_Edges_Only_Z', 
    createStepName='Step-Buckle',
    region=a.sets['All_Edges'], 
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
        print("[run_abaqus] 좌굴 잡 모델=%s, *NODE FILE 삽입=%s" % (BUCKLE_MODEL, _found))
    except Exception as _e:
        print("!!! WARNING: 모델 복사 실패(%s) -> 원본에 삽입 (결함 A 재발 가능)" % _e)
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
    
    # 낮은 값으로 다시 생성 
    my_model.Stress(
        name='Initial_Stiffness',
        region=inst_memb.sets['All'],
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
        region=inst_memb.sets['All'],
        distributionType=UNIFORM,
        sigma11=SIGMA0, sigma22=SIGMA0, sigma33=0.0, 
        sigma12=0.0, sigma13=0.0, sigma23=0.0
    )

    if RUN_SELF_BUCKLE_JOB:
        # 자기(클램프) 좌굴 잡: 모드 소스가 아니라 '클램프 base state 증거 + 실패 진단'용이다.
        #   (모드 추출은 MODE_SOURCE='cable' 경로가 담당 -> 클램프 모델은 CONVERGED=0.)
        run_job_safely('Buckle_Analysis', model_name=BUCKLE_MODEL)   # P0-2: 500 Pa 상태에서 좌굴모드 산출
    else:
        print("[IMPERFECTION] 자기 좌굴 잡 건너뜀 (RUN_SELF_BUCKLE_JOB=False)")

    # 좌굴 모드 병합(coalescence) 진단용 고유치 기록
    #   - 이유: Buckle_Analysis.dat 는 '다음 설계점'의 좌굴 잡이 시작될 때 삭제되므로
    #           여기서 고유치를 뽑아 이력(JSONL)에 남기지 않으면 회고 분석이 불가능하다.
    #   - 기록 전용: 실패해도 해석에는 전혀 영향을 주지 않는다 (예외 전부 삼킴).
    #   - 끄려면 EIG_RECORD 상수를 False 로 (구 MFBO_EIG_RECORD 환경변수는 제거됨)
    try:
        if EIG_RECORD:
            _rec = os.path.join(_HERE, 'coalescence_check.py')
            _dat = os.path.join(os.getcwd(), 'Buckle_Analysis.dat')
            _msg = os.path.join(os.getcwd(), 'Buckle_Analysis.msg')
            _hist = os.path.join(os.getcwd(), 'eig_history.jsonl')
            if os.path.exists(_rec) and (os.path.exists(_dat) or os.path.exists(_msg)):
                _cmd = ('abaqus python "%s" record --dat "%s" --msg "%s" --history "%s" '
                        '--x %s --d %s --fidelity HF' % (_rec, _dat, _msg, _hist, x_c, d_c))
                _rc = subprocess.call(_cmd, shell=True)
                print("[N-6] coalescence record rc=%d (x_c=%s, d_c=%s)" % (_rc, x_c, d_c))
            else:
                print("[N-6] 고유치 기록 건너뜀 (script=%s, dat=%s, msg=%s)"
                      % (os.path.exists(_rec), os.path.exists(_dat), os.path.exists(_msg)))
    except Exception as _eig_err:
        print("[N-6] coalescence record 실패(무시): %s" % _eig_err)

    # [R-13 / 스파이크 10-0] base state(GlobalTension 말단) 응력·반력 측정
    #   좌굴 성공/실패와 무관하게 Buckle_Analysis.odb 의 GlobalTension 프레임에서 읽는다.
    #   -> "프리텐션 운용점이 물리적인가"를 실측으로 답하기 위한 것 (미해결 최우선 1건).
    #   끄려면 BASE_PROBE 상수를 False 로 (구 MFBO_BASE_PROBE 환경변수는 제거됨)
    try:
        if BASE_PROBE:
            _probe = os.path.join(_HERE, 'base_state_probe.py')
            if os.path.exists(_probe) and os.path.exists(BUCKLE_ODB):
                _pcmd = ('abaqus python "%s" "%s" Step-GlobalTension'
                         % (_probe, BUCKLE_ODB))
                subprocess.call(_pcmd, shell=True)
            else:
                print("[R-13] base state 측정 건너뜀 (probe=%s, odb=%s)"
                      % (os.path.exists(_probe), os.path.exists(BUCKLE_ODB)))
    except Exception as _probe_err:
        print("[R-13] base state 측정 실패(무시): %s" % _probe_err)

    # 기존 Step 정리: Post-buckling은 GlobalTension 직후에서 시작하며,
    # 중간 단계(ClampTension)를 건너뛰고 바로 최종 하중으로 Ramping함
    if 'Step-Buckle' in my_model.steps: del my_model.steps['Step-Buckle']
    if 'Step-ClampTension' in my_model.steps: del my_model.steps['Step-ClampTension']

    # Static, General 포스트버클링 수행    
    my_model.StaticStep(
        name='Step-Postbuckle', 
        previous='Step-GlobalTension',  
        nlgeom=ON, 
        stabilizationMagnitude=0.0002,      # Galhofo 참조 2e-4 
        stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
        continueDampingFactors=False,
        adaptiveDampingRatio=0.05,
        initialInc=1e-4,
        minInc=1e-8,          # R-9: 1e-15 는 발산 시 증분 소진까지 수시간
        maxInc=0.1,
        maxNumInc=1000        # R-9
    )
    my_model.keywordBlock.synchVersions(storeNodesAndElements=False)

    my_model.boundaryConditions['BC_Stabilize_Z'].deactivate('Step-Postbuckle')

    # Step-Buckle 삭제 시 BC_Edges_Only_Z(prescribed condition)가 함께 삭제되므로
    # Postbuckle 스텝에 모서리 z 구속을 재생성한다 (Galhofo 참조: 3개 모서리 u3=0)
    my_model.DisplacementBC(
        name='BC_Edges_Only_Z',
        createStepName='Step-Postbuckle',
        region=a.sets['All_Edges'],
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
    
    # --- 2-모델 레시피: 임퍼펙션 소스 (2026-09-22) ---------------------------
    #   모드 형상은 '클램프 없는' 모델에서 나온 것을 쓴다(Galhofo 와 동일한 구조).
    #   소스 .fil 을 IMPERFECTION_NAME 으로 스테이징 -> 자기 좌굴 잡의 0-모드 .fil 이
    #   조용히 소비되는 사고(원장 C-1)를 이름 분리로 차단하고, 모드 개수를 검증한다.
    imp_scale = THICKNESS * IMPERFECTION_AMPL_T   # Galhofo 채택값 0.10 t
    if IMPERFECTION_MODE == 'odb_table':
        # 이미 모델 빌드 단계에서 노드 좌표를 섭동했다(ODB 모드표). 키워드 경로는 쓰지 않는다.
        imp_name, imp_step = None, None
        print("[IMPERFECTION] 주입=기하 섭동 (모드표 %s, 모드 %s, 진폭 %.2f t = %.3e m)"
              % (MODE_TABLE, list(IMPERFECTION_MODES), IMPERFECTION_AMPL_T, imp_scale))
        print("               %s" % perturbation_report(_pert))
    elif MODE_SOURCE == 'self':
        imp_name, imp_step = 'Buckle_Analysis', _BUCKLE_STEP_NO
        print("[IMPERFECTION] 소스=self 좌굴 잡(model=%s) STEP=%d amplitude=%.3e m"
              % (BUCKLE_MODEL, imp_step, imp_scale))
    else:
        try:
            _st = stage(MODE_SOURCE_FIL, IMPERFECTION_NAME, os.getcwd(),
                        IMPERFECTION_MODES, dat_hint=MODE_SOURCE_DAT,
                        msg_hint=MODE_SOURCE_MSG, here=_HERE)
        except ImperfectionSourceError as _imp_err:
            print("!!! ERROR: 임퍼펙션 소스 스테이징 실패 -> HF 잡을 제출하지 않고 중단합니다.")
            print(str(_imp_err))
            sys.exit(1)
        imp_name, imp_step = _st['name'], MODE_SOURCE_STEP
        print("[IMPERFECTION] %s STEP=%d amplitude=%.3e m (모드 %s)"
              % (report(_st), imp_step, imp_scale, list(IMPERFECTION_MODES)))
        if _st['modes'] < N_EIG:
            print("!!! WARNING: 소스 모드 %d개 < N_EIG=%d -> 요청 모드 일부만 주입됩니다."
                  % (_st['modes'], N_EIG))

    if IMPERFECTION_MODE == 'odb_table':
        print("[IMPERFECTION] 키워드 삽입 생략 (기하 섭동으로 이미 주입됨)")
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
    print("Calling extraction script: %s" % cmd)
    # shell=True로 eval_abaqus.py 실행
    p = subprocess.Popen(cmd, shell=True)
    rc = p.wait()
    if rc != 0:
        print("!!! ERROR: eval_abaqus.py exited with code %d" % rc)
    
    # 추출 스크립트가 출력한 "RESULTS:..." 라인을 찾아 전달
    if os.path.exists('extraction.txt'):
        with open('extraction.txt', 'r') as f:
            print("RESULTS:" + f.read().strip())
    else:
        print("!!! ERROR: Extraction failed. 'extraction.txt' not found.") 
        print("RESULTS:FAIL")   # 상위(mfbo)가 원인을 식별하도록 명시적 실패 신호

except Exception as err:
    print("Error during data extraction: %s" % str(err))
    sys.exit(1)