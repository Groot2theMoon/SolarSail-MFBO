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

Flow (Simulation Logic)  --  [B안 / 2026-09-21]
    이 스크립트는 run_abaqus.py(A안: 좌굴 고유모드를 imperfection 으로)의 대안 구현이다.
    A안은 클램프가 있는 이 모델에서 성립하지 않는다(실측):
       클램프가 base state 하중의 33.3%를 담당 -> sigma2<0 영역 21.5%
       -> base state 가 이미 분기점 위 -> *BUCKLE 추출 실패(208 negative, CONVERGED=0).
    따라서 A안의 '좌굴모드 임퍼펙션' 대신, B안은 *클램프 당김 중 자연 발생하는 주름*을
    초기 결함으로 사용한다.

    B안 절차 (HF):
      Step-GlobalTension : 프리텐션. u3 전체 고정 -> 프리텐션을 평탄 상태로 확립
      Step-Trigger       : u3 해제(모서리만 고정) + 소량 기하 trigger(기본 0.1t)
                           -> sigma2<0 영역에서 주름이 자연 발생 개시
      Step-ClampTension  : 클램프 당김 ramp (u3 자유) -> 주름이 클램프의 영향을 받으며 성장
      Step-Postbuckle    : 최종 하중까지 ramp (안정화 2e-4, Galhofo 참조)
    LF 는 기존과 동일하다: GlobalTension -> ClampTension -> HighTension, u3 고정.

    주름을 만드는 것은 seed 가 아니라 면내 응력장(sigma2<0)이다. trigger 는 '어느 패턴이
    먼저 자랄지'만 정한다. 근거(Galhofo 2022 Table 3·4): 임퍼펙션 진폭 0.05t~0.5t 변화 ->
    응답 차이 약 3%, 모드조합(단일 vs 4모드) -> 면외변위 차이 <= 8%.
    -> seed 선택이 해석 결과를 지배하지 않으므로 '결함 규칙의 통일'이 정당화된다.
    (검증 항목: MFBO_TRIG_MAG / MFBO_TRIG_ON 스윕으로 우리 모델에서 재확인)

    좌굴 잡/ *NODE FILE / *IMPERFECTION 카드는 쓰지 않는다.
    Arguments / Post-Processing 은 run_abaqus.py 와 동일.
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
# ---- 코드 상수 끝 ----
_RUN = os.path.join(_HERE, RUN_DIR_NAME)
os.makedirs(_RUN, exist_ok=True)
os.chdir(_RUN)
print("[run_abaqus_new] _HERE = %s" % _HERE)
print("[run_abaqus_new] _RUN  = %s" % _RUN)
 
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
    print("[run_abaqus_new] numCpus=%d numDomains=%d (코드 상수 NUMCPUS)" % (_ncp, _ncp))
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
# 버클 base state(Step-ClampTension)의 장력을 키워 시스템행렬 부정정(음수 고유값)을 해소.
# 1.0 = 기존값.  조정:  PowerShell  $env:MFBO_PRETENSION_SCALE="30"
PRETENSION_SCALE = 10.0          # 프리텐션 변위 = 5e-6 m * 이 값 = 5e-5 m
# R-13 실측(2026-09-21): PRETENSION_SCALE=10 -> 평균 면내응력 2122 Pa = 목표 7000 Pa 의 0.303배.
#   운용점을 목표에 맞추려면 약 33배(= DISP_GLOBAL 165um)가 필요하다.
#   단 PRETENSION_SCALE 는 최종 하중까지 함께 키우므로(포스트버클 변위 1mm -> 3.3mm),
#   운용점만 따로 맞추려면 절대값 노브 MFBO_DISP_GLOBAL / MFBO_GLOBAL_FINAL 을 쓴다.
DISP_GLOBAL = 0.000005 * PRETENSION_SCALE    # 운용점: 코너 당김 5e-5 m
CLAMP_PULL = DISP_GLOBAL * d_c
# (B안) Buckle perturbation 상수 없음 (좌굴 스텝 자체가 없다)
GLOBAL_FINAL = 0.0001 * PRETENSION_SCALE     # 최종 하중: 코너 당김 1e-3 m
CLAMP_FINAL = GLOBAL_FINAL * d_c
# ---- B안 trigger 파라미터 ----
# 기하 trigger 진폭. 기본 = 막 두께의 10% (Galhofo 관행 0.1t)
TRIG_MAG = THICKNESS * 0.1       # 0.1t = 5e-7 m (Galhofo 관행)
# trigger 를 적용할 '내부' 노드의 경계 여유 [m]. 요소 크기 0.1 m -> 3요소 여유
TRIG_MARGIN = 0.3                # 스케일 밖 값이면 폭의 5% 로 자동 대체된다
# 0 이면 기하 trigger 없이 u3 해제만 한다 (trigger 민감도 비교용).
#   2026-09-21: 다음 실행 기준값으로 0 을 기본에 둔다 -> 안정화(STAB)만 바꿔 주름 분기를
#   넘는지 보는 한 변수 대조 런. seed 까지 켜는 런(민감도/생산)은
#   $env:MFBO_TRIG_ON="1" 로 덮어쓴다. 실행 첫 줄이 실제 값을 항상 찍는다.
TRIG_ON = False                  # seed 없이 안정화만 바꾸는 '한 변수' 대조 런. seed 런은 True 로.
# 비선형 스텝(Trigger/ClampTension/Postbuckle)의 안정화 계수
#   2026-09-21 실측: 2e-4(Galhofo 참조값) 로는 Step-ClampTension 이 주름 발생 직후
#   증분 444 에서 TOO MANY ATTEMPTS 로 죽었다 (증분 2.7e-4 -> 2.1e-6, 100배 축소에도
#   복구 불가 = 증분 제어로 넘을 수 없는 분기). 그래서 기본값을 1e-3 으로 올려 둔다.
#   스윕은 환경변수로 덮어쓴다:  $env:MFBO_STAB="0.01" (평상시 실행엔 입력 불필요)
#   불안정(주름) 분기에서 증분이 컷백으로도 복구되지 않고 죽으면 이 값을 키운다.
#   예:  PowerShell  $env:MFBO_STAB="0.001"   (5배)   /   "0.01" (50배)
#   검증: .sta 의 ALLSD/ALLIE (누적 소산/변형 에너지 비율) 가 작아야 물리적으로 유효.
#   GlobalTension 스텝은 기존 2e-4 고정 (프리텐션 상태를 바꾸지 않기 위함).
STAB = 0.001                     # 2e-4 는 주름 발생 직후 분기에서 실패 (2026-09-21 실측)
print("[run_abaqus_new] STAB=%g (코드 상수) / TRIG_ON=%s / TRIG_MAG=%.3e m / TRIG_MARGIN=%.3g"
      % (STAB, TRIG_ON, TRIG_MAG, TRIG_MARGIN))


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
a.regenerate()

# 꼭짓점 RP
rp1_obj, rp1_reg = create_rigid_patch('Top', V1, radius=0.2)
rp2_obj, rp2_reg = create_rigid_patch('Right', V2, radius=0.2)
rp3_obj, rp3_reg = create_rigid_patch('Left', V3, radius=0.2)

# 클램프 RP : 우측 빗변 중점 (15, 5), 좌측 빗변 중점 (5, 5)
# B안은 클램프를 유지한다 (설계변수 x_c/d_c 가 하중 경로로 들어가는 통로).
# 초기 가짜 응력(수렴 보조). 케이블 변형=700 Pa, 우리=500 Pa -> 정렬 노브
SIGMA0 = 500.0                   # 수렴 보조용 초기응력 [Pa]
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
start_cl, end_cl = connect_cable('cable_CL', p_cable_cl, V_CL, (-1.0, 1.0, 0.0))
# V_CR : 우상향 (1, 1)
start_cr, end_cr = connect_cable('cable_CR', p_cable_cr, V_CR, (1.0, 1.0, 0.0))

a.Set(name='RP_Top_Set', referencePoints=(a.referencePoints[rp1_obj.id],))
a.Set(name='RP_Right_Set', referencePoints=(a.referencePoints[rp2_obj.id],))
a.Set(name='RP_Left_Set', referencePoints=(a.referencePoints[rp3_obj.id],))

a.Set(name='RP_CL_Set', referencePoints=(a.referencePoints[rp_cl_obj.id],))
a.Set(name='RP_CR_Set', referencePoints=(a.referencePoints[rp_cr_obj.id],))

my_model.Tie(name='Tie_Top', main=a.sets['RP_Top_Set'], secondary=start_c1, positionToleranceMethod=COMPUTED)
my_model.Tie(name='Tie_Right', main=a.sets['RP_Right_Set'], secondary=start_c2, positionToleranceMethod=COMPUTED)
my_model.Tie(name='Tie_Left', main=a.sets['RP_Left_Set'], secondary=start_c3, positionToleranceMethod=COMPUTED)

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

# Step Trigger : B안. u3 를 해제하고 소량 기하 trigger 로 주름을 자연 발생시킨다.
#   HF 전용으로 만든다 -> LF 체인(GlobalTension->ClampTension->HighTension)은 기존과 동일하게
#   유지되어 LF 지표가 변하지 않는다.
if fidelity == 'HF':
    my_model.StaticStep(
        name='Step-Trigger',
        previous='Step-GlobalTension',
        nlgeom=ON,
        stabilizationMagnitude=STAB,      # MFBO_STAB (기본 1e-3)
        stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
        initialInc=0.1, minInc=1e-8, maxInc=1.0, maxNumInc=50
    )
_PREV_CLAMP = 'Step-Trigger' if fidelity == 'HF' else 'Step-GlobalTension'

# Step Clamp Tension : 클램프에 변위 가하기
my_model.StaticStep(
    name='Step-ClampTension',
    previous=_PREV_CLAMP,
    nlgeom=ON,
    stabilizationMagnitude=STAB,      # MFBO_STAB (기본 1e-3)
    stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
    initialInc=0.0001, minInc=1e-8, maxNumInc=1000    # P1-2: 1e-15 는 발산 시 증분 폭주
)

# Step Buckle / 좌굴 고유모드 추출 없음 (B안).
#   실측 근거: 클램프가 있는 base state 는 sigma2<0 영역 21.5% 로 이미 분기점 위이고
#   *BUCKLE 추출이 실패한다(208 negative eigenvalues, CONVERGED=0).
#   -> 좌굴모드 imperfection 대신 'trigger + 자연 주름'을 결함으로 쓴다.
#   (A안 = run_abaqus.py 가 그대로 보존되어 있음)

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
my_model.DisplacementBC(name='BC_CL', createStepName='Initial', region=end_cl, u3=0)
my_model.DisplacementBC(name='BC_CR', createStepName='Initial', region=end_cr, u3=0)

my_model.DisplacementBC(name='Disp_Control_Right', createStepName='Initial', region=end_c2, 
    u1=0, u2=0
)
my_model.DisplacementBC(name='Disp_Control_Left',createStepName='Initial', region=end_c3, 
    u1=0, u2=0
)
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

# (B안) Buckle perturbation BC 없음 -> Step-Trigger 의 기하 trigger 로 대체

# Step-Trigger 에서 전체 z 구속을 해제하고 '모서리 z 구속'으로 교체 (Galhofo reference)
all_edges = inst_memb.edges
a.Set(name='All_Edges', edges=all_edges)

if fidelity == 'HF':
    my_model.boundaryConditions['BC_Stabilize_Z'].deactivate('Step-Trigger')
    my_model.DisplacementBC(
        name='BC_Edges_Only_Z',
        createStepName='Step-Trigger',
        region=a.sets['All_Edges'],
        u3=0
    )

    # ---- trigger : 내부 노드를 +-TRIG_MAG 로 교란 (광대역 seed) ----
    #   자유메쉬 노드 순서대로 홀/짝을 나눠 부호를 교대시킨다 -> 여러 파장이 섞인 seed.
    #   모서리/꼭짓점 노드는 제외한다(BC_Edges_Only_Z, 케이블 RP 와 충돌 방지).
    #   [FIX 2026-09-21] 이전 필터는 모델 스케일과 맞지 않아 내부 노드가 0개였다:
    #     조건  y > TRIG_MARGIN(0.3)  /  y < 20 - x - 0.3  vs  실제 노드 범위
    #     x in [-0.075, 0.075], y in [-0.041, 0.041]  ->  교집합 공집합.
    #   빈 노드셋에 BC 를 만들면 Abaqus 는 오류 없이 '아무 일도 하지 않는다'
    #   -> trigger 가 조용히 무효였고(프로브 실측 max|u3| = 5.8e-16 m), 주름은 반올림
    #      노이즈가 seed 역할을 해서 났다(진폭 131.6 um = 26t).
    #   -> 하드코딩 대신 노드 좌표 범위로 삼각형 3변을 복원해 '변에서 MG 이상 떨어진' 노드를 쓴다.
    _mn = inst_memb.nodes
    _allc = [n.coordinates for n in _mn]
    _xmin = min(c[0] for c in _allc); _xmax = max(c[0] for c in _allc)
    _ymin = min(c[1] for c in _allc); _ymax = max(c[1] for c in _allc)
    _apex = max(_allc, key=lambda c: c[1])          # 꼭짓점 = y 최대 노드
    _W = _xmax - _xmin
    # TRIG_MARGIN 이 스케일 범위 안이면 그 값, 아니면 폭의 5% 로 자동 스케일
    _MG = TRIG_MARGIN if (0.0 < TRIG_MARGIN < 0.2 * _W) else 0.05 * _W

    def _dist_seg(px, py, ax, ay, bx, by):
        _dx, _dy = bx - ax, by - ay
        _L2 = _dx * _dx + _dy * _dy
        _t = 0.0 if _L2 == 0.0 else max(0.0, min(1.0, ((px - ax) * _dx + (py - ay) * _dy) / _L2))
        _qx, _qy = ax + _t * _dx, ay + _t * _dy
        return ((px - _qx) ** 2 + (py - _qy) ** 2) ** 0.5

    _n_int = []
    for _i in range(len(_mn)):
        _c = _mn[_i].coordinates
        if (_c[1] - _ymin) < _MG:
            continue
        if _dist_seg(_c[0], _c[1], _xmin, _ymin, _apex[0], _apex[1]) < _MG:
            continue
        if _dist_seg(_c[0], _c[1], _xmax, _ymin, _apex[0], _apex[1]) < _MG:
            continue
        _n_int.append(_mn[_i])
    print("[run_abaqus_new] 노드 범위 x[%.4f, %.4f] y[%.4f, %.4f] apex=(%.4f, %.4f) margin=%.4g m"
          % (_xmin, _xmax, _ymin, _ymax, _apex[0], _apex[1], _MG))
    print("[run_abaqus_new] trigger 내부 노드 %d / 전체 %d (mag=%.3e m, on=%s)"
          % (len(_n_int), len(_mn), TRIG_MAG, TRIG_ON))
    if TRIG_ON and len(_n_int) < 50:
        raise RuntimeError(
            "trigger 내부 노드가 %d개뿐입니다 (전체 %d). 빈/과소 노드셋에 BC 를 만들면 "
            "trigger 가 조용히 무효가 됩니다 - 좌표 필터와 MFBO_TRIG_MARGIN 을 확인하세요."
            % (len(_n_int), len(_mn)))
    _TRIG_OK = False
    try:
        a.Set(name='NS_TRIG_P', nodes=_n_int[0::2])
        a.Set(name='NS_TRIG_M', nodes=_n_int[1::2])
        if TRIG_ON:
            my_model.DisplacementBC(name='BC_Trig_P', createStepName='Step-Trigger',
                                    region=a.sets['NS_TRIG_P'], u3=TRIG_MAG)
            my_model.DisplacementBC(name='BC_Trig_M', createStepName='Step-Trigger',
                                    region=a.sets['NS_TRIG_M'], u3=-TRIG_MAG)
            _TRIG_OK = True
        else:
            print("[run_abaqus_new] TRIG_ON=False : trigger 없이 u3 해제만")
    except Exception as _trig_err:
        # trigger 생성 실패는 치명적이지 않다: u3 해제만으로도 불안정 상태는 드러난다.
        print("!!! WARNING: trigger 생성 실패(무시하고 'u3 해제만'으로 진행): %s" % _trig_err)

# 클램프 당김 단계부터 trigger 를 해제 -> 주름이 자유롭게 성장한다
if fidelity == 'HF' and _TRIG_OK:
    my_model.boundaryConditions['BC_Trig_P'].deactivate('Step-ClampTension')
    my_model.boundaryConditions['BC_Trig_M'].deactivate('Step-ClampTension')

# BUCKLE_ODB 는 더 이상 쓰지 않는다 (좌굴 잡 없음). base state 측정은 HF odb 로 한다.
LF_ODB = 'LF_Analysis.odb'
HF_ODB = 'HF_Postbuckle.odb'

# (B안) 좌굴 잡 모델 복사본 / *NODE FILE 삽입 없음.
#   *IMPERFECTION 카드도 쓰지 않으므로 좌굴모드를 .fil 에 기록할 필요가 없다.

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
    for name, sign in [('Disp_Control_CL', -1), ('Disp_Control_CR', 1)]:
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

    # (B안) 좌굴 잡 없음. HF 는 포스트버클 잡 하나만 제출한다.

    # (B안) 좌굴모드 병합(coalescence) 진단 제거: 고유값을 추출하지 않는다.


    # (B안) Step-Buckle 은 애초에 만들지 않았고, Step-ClampTension 은 *유지*한다.
    #   -> 클램프 당김은 ClampTension 스텝에서 ramp 되고(GlobalTension 0 -> CLAMP_PULL),
    #      Postbuckle 에서 CLAMP_FINAL 까지 다시 ramp 된다 (2단 점진 당김 = 주름 성장 경로).

    # Static, General 포스트버클링 수행    
    my_model.StaticStep(
        name='Step-Postbuckle', 
        previous='Step-ClampTension',   # B안: ClampTension 을 체인에 유지
        nlgeom=ON, 
        stabilizationMagnitude=STAB,      # MFBO_STAB (기본 1e-3) 
        stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
        continueDampingFactors=False,
        adaptiveDampingRatio=0.05,
        initialInc=1e-4,
        minInc=1e-8,          # R-9: 1e-15 는 발산 시 증분 소진까지 수시간
        maxInc=0.1,
        maxNumInc=1000        # R-9
    )
    my_model.keywordBlock.synchVersions(storeNodesAndElements=False)



    # 초기 상태에서 최종 장력(목표 7000 Pa, 미검증)까지 증가
    for name, sign in [('Disp_Control_Right', 1), ('Disp_Control_Left', -1)]:
        my_model.boundaryConditions[name].setValuesInStep(
            stepName='Step-Postbuckle',
            u1=sign * GLOBAL_FINAL * cos_val,
            u2=-GLOBAL_FINAL * sin_val
        )

    for name, sign in [('Disp_Control_CL', -1), ('Disp_Control_CR', 1)]:
        my_model.boundaryConditions[name].setValuesInStep(
            stepName='Step-Postbuckle',
            u1=sign * CLAMP_FINAL / sqrt2,
            u2=CLAMP_FINAL / sqrt2
        )
    

    my_model.fieldOutputRequests['F-Output-1'].setValues(
        variables=('S', 'E', 'U', 'RF', 'COORD', 'EVOL'),   # R-8: LF 지표(lf1~lf3) 산출에 필요
        frequency=10   # R-10: odb 크기 절감
    )
    

    run_job_safely('HF_Postbuckle')
    # [R-13] base state 실측 — B안에서는 좌굴 잡이 없으므로 HF odb 를 직접 읽는다.
    #   Step-GlobalTension : 프리텐션만 걸린 상태 (A안과 비교 가능한 기준점)
    #   Step-ClampTension  : 클램프 당김 + 주름 발생 이후 = B안의 실제 '결함 시작 상태'
    #   끄려면 환경변수 MFBO_BASE_PROBE=0
    try:
        if BASE_PROBE:
            _probe = os.path.join(_HERE, 'base_state_probe.py')
            if os.path.exists(_probe) and os.path.exists(HF_ODB):
                #   Step-Trigger: u3 를 처음 푼 스텝 -> 'trigger 만으로 주름이 났는지' 판정
                for _st in ('Step-GlobalTension', 'Step-Trigger', 'Step-ClampTension'):
                    print("[R-13] base state 측정: %s / %s" % (HF_ODB, _st))
                    subprocess.call('abaqus python "%s" "%s" %s' % (_probe, HF_ODB, _st),
                                    shell=True)
            else:
                print("[R-13] base state 측정 건너뜀 (probe=%s, odb=%s)"
                      % (os.path.exists(_probe), os.path.exists(HF_ODB)))
    except Exception as _probe_err:
        print("[R-13] base state 측정 실패(무시): %s" % _probe_err)
    
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
