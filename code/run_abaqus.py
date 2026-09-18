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
       [LF 모드]
         - 'Step-GlobalTension' -> 'Step-ClampTension' -> 'Step-HighTension' 순차 진행.
         - 주름(Wrinkle) 거동을 무시하고 선형적(혹은 단순 기하 비선형) 강성을 빠르게 계산.
         - 'Buckle_Analysis'는 모드 확인용으로만 돌리고 실제 결과엔 반영 안 함.
       [HF 모드]
         - 1차: 'Step-GlobalTension' -> 'Step-ClampTension' 상태에서 고유치 해석(Buckle) 수행.
         - 2차: 1차 해석 결과(.odb)에서 고유모드를 추출하여 초기 결함(Imperfection)으로 주입.
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
    env_dir = os.environ.get("MFBO_CODE_DIR")
    if _ok(env_dir):
        return os.path.abspath(env_dir)
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
    cand.append(os.path.join(os.getcwd(), "code"))
    for c in cand:
        if _ok(c):
            return os.path.abspath(c)
    print("!!! WARNING: eval_abaqus.py 위치를 찾지 못했습니다. cwd=%s" % os.getcwd())
    return os.getcwd()

_HERE = _resolve_here()
# Abaqus 작업 디렉터리(=산출물 위치) = code/aba.
# mfbo.py 가 MFBO_RUN_DIR 를 넘겨주면 그 값을 그대로 사용 (직접 실행해도 동일하게 동작)
_RUN = os.path.abspath(os.environ.get("MFBO_RUN_DIR") or os.path.join(_HERE, "aba"))
os.makedirs(_RUN, exist_ok=True)
os.chdir(_RUN)
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

def run_job_safely(job_name):
    """
    job 실행 중 .odb, .lck 파일 충돌을 방지하고, 완료까지 대기하는 함수

    P1-3: 제출 전에 이전 실행 산출물을 지운다.
      - 성공 판정(not os.path.exists(odb))이 '지난 실행의 odb'를 보고 오판하는 것을 막고,
      - *IMPERFECTION, FILE=Buckle_Analysis 가 stale .fil 을 조용히 읽는 것을 막는다.
    """
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
    
    job = mdb.Job(name=job_name, model=MODEL_NAME, numCpus=1, numDomains=1)
    print("Submitting Job: %s" % job_name)
    job.writeInput(consistencyChecking=OFF)
    job.submit(consistencyChecking=OFF)
    
    job.waitForCompletion()
    
    time.sleep(1.0) 
    
    # ABORTED가 아니면서, ODB 파일이 실제로 존재하면 성공으로 간주
    if job.status == ABORTED or not os.path.exists(odb_file):
        print("!!! ERROR: Job %s failed. Actual Status: %s" % (job_name, str(job.status)))
        sys.exit(1)
        
    print("Job %s completed successfully (Status: %s)." % (job_name, str(job.status)))
    return True

sqrt2 = 1.414
N_EIG = 4          # 임퍼펙션에 쓸 좌굴모드 수 (*IMPERFECTION / *NODE FILE)
N_EIG_BUCKLE = 40  # A: 추출 요청 고유값 수 (음수모드 우회 여유; run_abaqus_cable 은 100)

MODEL_NAME = 'SailModel_Triangle'
INSTANCE_NAME = 'MEMBRANE-1'

BASE = 20.0   # m
HEIGHT = 10.0 # m
THICKNESS = 2.5e-6 # Galhofo Reference 에서는 2.5e-6
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
# ---- C: 사전 장력(prestrain) 캘리브레이션 ----
# 버클 base state(Step-ClampTension)의 장력을 키워 시스템행렬 부정정(음수 고유값)을 해소.
# 1.0 = 기존값.  조정:  PowerShell  $env:MFBO_PRETENSION_SCALE="30"
PRETENSION_SCALE = float(os.environ.get("MFBO_PRETENSION_SCALE", "10"))
DISP_GLOBAL = 0.000005 * PRETENSION_SCALE
CLAMP_PULL = DISP_GLOBAL * d_c
PERTURBATION = 0.0005
CLAMP_PERT = PERTURBATION * d_c
GLOBAL_FINAL = 0.0001 * PRETENSION_SCALE
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
elemTypeQuad = ElemType(
    elemCode=S4R,
    elemLibrary=STANDARD,
    kinematicSplit=AVERAGE_STRAIN,
    hourglassControl=ENHANCED
)
elemTypeTri = ElemType(elemCode=S3, elemLibrary=STANDARD)  
p.setElementType(regions=(p.faces,), elemTypes=(elemTypeQuad, elemTypeTri))
p.generateMesh()
a.regenerate()

# 꼭짓점 RP
rp1_obj, rp1_reg = create_rigid_patch('Top', V1, radius=0.2)
rp2_obj, rp2_reg = create_rigid_patch('Right', V2, radius=0.2)
rp3_obj, rp3_reg = create_rigid_patch('Left', V3, radius=0.2)

# 클램프 RP : 우측 빗변 중점 (15, 5), 좌측 빗변 중점 (5, 5)
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
# Abacus BuckleStep을 사용하는 것이 정석이고 옳으나, 매우 얇은 solar-sail 자체의 불안정성에 의해 
# 부적합했음이 확인되어 BuckleStep + subspace iteration 으로 복원 (R-6). LANCZOS 는 강성행렬이
# indefinite 한 좌굴 해석에서 금지이며(예비하중 초과 시 해석 종료), FREQUENCY 모드는 질량
# 정규화라 '0.1t' 임퍼펙션 관행이 성립하지 않는다. Galhofo 참조도 BUCKLE+subspace, 최초 4모드.
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
my_model.BuckleStep(
    name='Step-Buckle',          
    previous='Step-ClampTension',  
    numEigen=N_EIG_BUCKLE,        # A: 음수모드 건너뛰기 위해 여유있게
    eigensolver=LANCZOS,          # A: SUBSPACE -> LANCZOS (부정정 대응)
    maxBlocks=DEFAULT,
)
# *IMPERFECTION, STEP=n 의 n 은 'Buckle_Analysis.fil 안의 스텝 번호'
# (현재 3 = GlobalTension/ClampTension/Buckle). 스텝 구성이 바뀌면 자동 추종.
_BUCKLE_STEP_NO = len(my_model.steps)

# R-7: *IMPERFECTION, FILE= 은 results file(.fil) 을 읽는다 -> 좌굴 모드를 .fil 에 기록해야
#      임퍼펙션이 실제로 주입된다 (미요청 시 조용히 무시됨)
my_model.keywordBlock.synchVersions(storeNodesAndElements=False)
for _i, _b in enumerate(my_model.keywordBlock.sieBlocks):
    if _b.strip().upper().startswith(('*BUCKLE', '*FREQUENCY')):
        my_model.keywordBlock.insert(_i + 1, '*NODE FILE, GLOBAL=YES, LAST MODE=%d\nU' % N_EIG)
        break
# ----

my_model.Stress(
    name='Initial_Stiffness',
    region=inst_memb.sets['All'],
    distributionType=UNIFORM,
    sigma11=500.0, sigma22=500.0, sigma33=0.0,
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

# Buckle Perturbation
my_model.boundaryConditions['Disp_Control_Right'].setValuesInStep(stepName='Step-Buckle', 
    u1=PERTURBATION * cos_val,
    u2=-PERTURBATION * sin_val
)
my_model.boundaryConditions['Disp_Control_Left'].setValuesInStep(stepName='Step-Buckle', 
    u1=-PERTURBATION * cos_val,
    u2=-PERTURBATION * sin_val
)
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
if fidelity == 'LF':

    # R-17: LF 는 좌굴모드를 쓰지 않는다(HF 분기가 같은 설계점에서 자체 실행) -> 비용 절감
    # run_job_safely('Buckle_Analysis')
    
    if 'Initial_Stiffness' in my_model.predefinedFields:
        del my_model.predefinedFields['Initial_Stiffness']
    
    # 낮은 값으로 다시 생성 
    my_model.Stress(
        name='Initial_Stiffness',
        region=inst_memb.sets['All'],
        distributionType=UNIFORM,
        sigma11=500.0, sigma22=500.0, sigma33=0.0, 
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

    # P0-2: 좌굴모드와 포스트버클 해석이 '같은 초기응력 상태'에서 계산되도록
    #       초기응력 재생성(500 -> 100 Pa)을 좌굴 잡 제출 '앞'으로 이동.
    if 'Initial_Stiffness' in my_model.predefinedFields:
        del my_model.predefinedFields['Initial_Stiffness']

    my_model.Stress(
        name='Initial_Stiffness',
        region=inst_memb.sets['All'],
        distributionType=UNIFORM,
        sigma11=500.0, sigma22=500.0, sigma33=0.0, 
        sigma12=0.0, sigma13=0.0, sigma23=0.0
    )

    run_job_safely('Buckle_Analysis')   # P0-2: 500 Pa 상태에서 좌굴모드 산출

    # 기존 Step 정리: Post-buckling은 GlobalTension 직후에서 시작하며,
    # 중간 단계(ClampTension)를 건너뛰고 바로 최종 하중으로 Ramping함 (수렴성 향상 전략)
    if 'Step-Buckle' in my_model.steps: del my_model.steps['Step-Buckle']
    if 'Step-ClampTension' in my_model.steps: del my_model.steps['Step-ClampTension']

    # Static, General 포스트버클링 수행    
    my_model.StaticStep(
        name='Step-Postbuckle', 
        previous='Step-GlobalTension',  
        nlgeom=ON, 
        stabilizationMagnitude=0.0002,      # R-4: Galhofo 참조 2e-4 (기존 0.05 = 250배)
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

    # R-2: Step-Buckle 삭제 시 BC_Edges_Only_Z(prescribed condition)가 함께 삭제되므로
    #      Postbuckle 스텝에 모서리 z 구속을 재생성한다 (Galhofo 참조: 3개 모서리 u3=0)
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
    
    # Imperfection Injection (Keyword Editing)
    # Buckle_Analysis.odb 파일의 결과(고유모드)를 초기 결함으로 주입
    my_model.keywordBlock.synchVersions(storeNodesAndElements=False)
    
    imp_scale = THICKNESS * 0.1 # Galhofo Reference

    imp_text = ("*IMPERFECTION, FILE=Buckle_Analysis, STEP=%d\n" % _BUCKLE_STEP_NO +
                "\n".join("%d, %e" % (m, imp_scale) for m in range(1, N_EIG + 1)))

    # 키워드 삽입 위치 찾기 : *STEP 블록 직전에 삽입하는 것이 안전함
    inserted = False
    for i, block in enumerate(my_model.keywordBlock.sieBlocks):
        # Step 정의 시작 부분 찾기
        if block.lower().strip().startswith('*step'):
            my_model.keywordBlock.insert(i-1, imp_text)
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