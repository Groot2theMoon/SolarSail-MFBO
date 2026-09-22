"""
run_abaqus.py 의 클램프 없이 solar-sail 과 꼭짓점 케이블만 있는 버전의 스크립트.
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
 
print("DEBUG: All sys.argv: " + str(sys.argv))

try:
    # abaqus cae noGUI=run_abaqus_cable.py -- [HF/LF]
    fidelity = sys.argv[-1].upper()
    
except:
    print("Error: Invalid arguments. Usage: abaqus cae noGUI=run_abaqus_cable.py -- [fidelity]")
    sys.exit(1)

def run_job_safely(job_name):
    lck_file = job_name + '.lck'
    odb_file = job_name + '.odb'
    # 기존 Lock 파일이 있다면 삭제 시도
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
    
    job = mdb.Job(name=job_name, model=MODEL_NAME)
    print("Submitting Job: %s" % job_name)
    job.submit(consistencyChecking=OFF)
    
    job.waitForCompletion()
    
    time.sleep(1.0) 
    
    # ABORTED가 아니면서, ODB 파일이 실제로 존재하면 성공으로 간주
    if job.status == ABORTED or not os.path.exists(odb_file):
        print("!!! ERROR: Job %s failed. Actual Status: %s" % (job_name, str(job.status)))
        sys.exit(1)
        
    print("Job %s completed successfully (Status: %s)." % (job_name, str(job.status)))
    return True


MODEL_NAME = 'SailModel_Triangle'
INSTANCE_NAME = 'MEMBRANE-1'

BASE = 20.0   # m
HEIGHT = 10.0 # m
# [수정] 두께를 5um로 설정하여 초기 강성 확보 (논문 2.5um는 수렴 매우 어려움)
THICKNESS = 5.0e-6 
TARGET_STRESS = 7000.0 # Pa

# 케이블 파라미터 (논문 참조)
CABLE_RADIUS = 5.0e-4 # m
CABLE_AREA = np.pi * (CABLE_RADIUS**2)
LEN_TOP = 0.280 # m
LEN_BOT = 0.689 # m

# 좌표 정의
V1 = (BASE/2.0, HEIGHT, 0.0) # Top
V2 = (BASE, 0.0, 0.0)        # Right
V3 = (0.0, 0.0, 0.0)         # Left

# 하중 각도 (28.6도)
angle_deg = 28.6
angle_rad = np.deg2rad(angle_deg)
cos_val = float(np.cos(angle_rad))
sin_val = float(np.sin(angle_rad))

# -------------------------------------------------------------
# 2. 모델 초기화 및 재질
# -------------------------------------------------------------
if MODEL_NAME in mdb.models: del mdb.models[MODEL_NAME]
my_model = mdb.Model(name=MODEL_NAME)

# (1) 멤브레인 재질 (Kapton)
mat = my_model.Material(name='Kapton')
mat.Elastic(table=((2.5e9, 0.34),))
my_model.HomogeneousShellSection(name='Section-Membrane', material='Kapton', thickness=THICKNESS)

# (2) 케이블 재질 (Kevlar)
mat_cable = my_model.Material(name='Kevlar')
mat_cable.Elastic(table=((62.0e9, 0.36),))
my_model.TrussSection(name='Section-Cable', material='Kevlar', area=CABLE_AREA)

# -------------------------------------------------------------
# 3. 파트 생성: 멤브레인
# -------------------------------------------------------------
s = my_model.ConstrainedSketch(name='triangle_profile', sheetSize=BASE*2)
s.Line(point1=V3[:2], point2=V2[:2])
s.Line(point1=V2[:2], point2=V1[:2])
s.Line(point1=V1[:2], point2=V3[:2])
p = my_model.Part(name='Membrane', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p.BaseShell(sketch=s)
p.SectionAssignment(region=p.Set(faces=p.faces, name='All'), sectionName='Section-Membrane')

# -------------------------------------------------------------
# 4. 파트 생성: 케이블 (함수화)
# -------------------------------------------------------------
def create_cable_part(name, length):
    # 길이만큼의 Truss Wire 생성
    p_c = my_model.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    # (0,0,0)에서 (length, 0, 0)으로 생성 (나중에 회전/이동)
    p_c.WirePolyLine(points=((0.0, 0.0, 0.0), (length, 0.0, 0.0)), mergeType=IMPRINT, meshable=ON)
    
    # 섹션 할당
    p_c.SectionAssignment(region=p_c.Set(edges=p_c.edges, name='Wire'), sectionName='Section-Cable')
    
    # 메쉬 (T3D2)
    p_c.seedPart(size=length) # 요소 1개
    elemTypeTruss = ElemType(elemCode=T3D2, elemLibrary=STANDARD)
    p_c.setElementType(regions=(p_c.edges,), elemTypes=(elemTypeTruss,))
    p_c.generateMesh()
    return p_c

# 케이블 파트 생성
p_cable_top = create_cable_part('Cable_Top', LEN_TOP)
p_cable_bot = create_cable_part('Cable_Bot', LEN_BOT)

# -------------------------------------------------------------
# 5. 어셈블리 및 연결
# -------------------------------------------------------------
a = my_model.rootAssembly
a.DatumCsysByDefault(CARTESIAN)
inst_memb = a.Instance(name=INSTANCE_NAME, part=p, dependent=ON)

# --- Rigid Patch 생성 함수 (기존 로직 유지) ---
def create_rigid_patch(name, coord, radius=0.4):
    # RP 생성
    rp = a.ReferencePoint(point=coord)
    rp_key = a.referencePoints[rp.id]
    rp_region = regionToolset.Region(referencePoints=(rp_key,))
    
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

# --- 케이블 배치 및 연결 함수 ---
def connect_cable(name, part, sail_corner, vector_dir, radius=1e-4):
    # 1. Instance 생성
    inst_name = 'Inst_' + name
    inst = a.Instance(name=inst_name, part=part, dependent=ON)
    
    # 2. 회전 및 이동
    # 케이블은 기본적으로 X축(1,0,0)으로 생성됨. 이를 vector_dir 방향으로 회전.
    # 회전축 계산 (Z축 회전)
    base_vec = np.array([1.0, 0.0, 0.0])
    target_vec = np.array(vector_dir)
    target_vec = target_vec / np.linalg.norm(target_vec) # Normalize
    
    # 2D 회전각 계산 (atan2)
    rot_angle_deg = np.degrees(np.arctan2(target_vec[1], target_vec[0]))
    
    # 회전 (0,0,0 기준 Z축 회전)
    a.rotate(instanceList=(inst_name,), axisPoint=(0,0,0), axisDirection=(0,0,1), angle=rot_angle_deg)
    
    # 이동 (시작점을 돛의 꼭짓점으로)
    a.translate(instanceList=(inst_name,), vector=sail_corner)
    
    # 3. 끝점(Loading Point) 찾기
    # 시작점은 sail_corner, 끝점은 sail_corner + vector_dir * length
    # 하지만 회전/이동했으므로 노드 좌표로 찾는 게 안전
    # 시작점 노드 (Constraint용)
    node_start = inst.nodes.getByBoundingSphere(center=sail_corner, radius=radius)
    if len(node_start) == 0:
        # 혹시 못 찾을 경우를 대비한 안전장치 (getClosest 후 배열로 변환 시도)
        # 하지만 getByBoundingSphere가 정석입니다.
        print("Warning: Node not found by sphere, trying closest for " + name)
        # Abaqus 버전에 따라 슬라이싱으로 MeshSequence를 만들 수 있음
        node_start = inst.nodes.getClosest(coordinates=sail_corner)

    region_start = regionToolset.Region(nodes=node_start)
    
    # 끝점 노드 (Load/BC용)
    # 끝점 좌표 계산
    cable_len = LEN_TOP if 'Top' in name else LEN_BOT
    end_coord = (sail_corner[0] + target_vec[0]*cable_len, sail_corner[1] + target_vec[1]*cable_len, 0.0)
    node_end = inst.nodes.getByBoundingSphere(center=end_coord, radius=radius)
    region_end = regionToolset.Region(nodes=node_end)

    return region_start, region_end

# 1. 메쉬 생성 (먼저 해야 노드 찾기 가능)
# [핵심] 아워글래싱 방지를 위해 S4(Full Integration) 사용
p.seedPart(size=BASE/200.0, deviationFactor=0.1) # 약 5~8천개 요소 목표
p.setMeshControls(regions=p.faces, elemShape=QUAD_DOMINATED, technique=FREE, algorithm=MEDIAL_AXIS)
elemTypeQuad = ElemType(elemCode=S4, elemLibrary=STANDARD) # Full Integration
elemTypeTri = ElemType(elemCode=S3, elemLibrary=STANDARD)  # Full Integration
p.setElementType(regions=(p.faces,), elemTypes=(elemTypeQuad, elemTypeTri))
p.generateMesh()
a.regenerate()

# 2. RP 생성
rp1_obj, rp1_reg = create_rigid_patch('Top', V1)
rp2_obj, rp2_reg = create_rigid_patch('Right', V2)
rp3_obj, rp3_reg = create_rigid_patch('Left', V3)

# 3. 케이블 연결
# V1(Top): 위로 (0, 1)
start_c1, end_c1 = connect_cable('Cable_Top', p_cable_top, V1, (0.0, 1.0, 0.0))
# V2(Right): 우하향 (+cos, -sin)
start_c2, end_c2 = connect_cable('Cable_Right', p_cable_bot, V2, (cos_val, -sin_val, 0.0))
# V3(Left): 좌하향 (-cos, -sin)
start_c3, end_c3 = connect_cable('Cable_Left', p_cable_bot, V3, (-cos_val, -sin_val, 0.0))

# 4. Tie Constraint (RP <-> Cable Start)
# Master: Sail RP, Slave: Cable Node (Node-to-Surface or TIE)
# RP는 Geometry가 아니므로 Tie에 직접 쓰기 까다로울 수 있음.
# RP를 Set으로 만들어서 Master로 지정.
a.Set(name='RP_Top_Set', referencePoints=(a.referencePoints[rp1_obj.id],))
a.Set(name='RP_Right_Set', referencePoints=(a.referencePoints[rp2_obj.id],))
a.Set(name='RP_Left_Set', referencePoints=(a.referencePoints[rp3_obj.id],))

# Cable Start Node는 Region이므로 Set으로 변환 필요 (create_cable_connection 리턴값 활용)
# 여기서는 편의상 Tie 정의
my_model.Tie(name='Tie_Top', main=a.sets['RP_Top_Set'], secondary=start_c1, positionToleranceMethod=COMPUTED)
my_model.Tie(name='Tie_Right', main=a.sets['RP_Right_Set'], secondary=start_c2, positionToleranceMethod=COMPUTED)
my_model.Tie(name='Tie_Left', main=a.sets['RP_Left_Set'], secondary=start_c3, positionToleranceMethod=COMPUTED)

# -------------------------------------------------------------
# 6. 경계 조건
# -------------------------------------------------------------
DISP_MAG = 0.000018

# Step 1: Global Tension
my_model.StaticStep(
    name='Step-GlobalTension', previous='Initial', nlgeom=ON,
    stabilizationMagnitude=0.0005, # 케이블과 두께 덕분에 작아도 됨
    stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
    initialInc=0.0001, minInc=1e-15, maxNumInc=1000
)

my_model.DisplacementBC(
    name='BC_Anchor_Top', 
    createStepName='Initial', 
    region=end_c1,  # Top Cable의 끝점
    u1=SET, u2=SET, u3=SET, 
    ur1=SET, ur2=SET, ur3=SET 
)

# Right: Z 고정
my_model.DisplacementBC(
    name='BC_Right_Z', 
    createStepName='Initial', 
    region=end_c2, 
    u3=SET, ur1=SET, ur2=SET, ur3=SET
)
# Left: Z 고정
my_model.DisplacementBC(
    name='BC_Left_Z', 
    createStepName='Initial',
    region=end_c3, 
    u3=SET, ur1=SET, ur2=SET, ur3=SET
)

my_model.DisplacementBC(
    name='Disp_Control_Right', 
    createStepName='Initial', 
    region=end_c2, 
    u1=SET, u2=SET
)

my_model.DisplacementBC(
    name='Disp_Control_Left',
    createStepName='Initial', 
    region=end_c3, 
    u1=SET, u2=SET
)

my_model.DisplacementBC(
    name='BC_Stabilize_Z', 
    createStepName='Step-GlobalTension', 
    region=inst_memb.sets['All'],  # 전체 면
    u3=SET
)


# Load: 케이블 방향으로 힘 적용 (Vector Math 필요 없음, 케이블 축방향이 곧 힘방향)
# 하지만 ConcentratedForce는 Global 축 기준이므로 분해 필요.
# 1. Top Cable: +Y 방향 당기기
# Top Cable: 위로(+Y) 당기기
my_model.boundaryConditions['Disp_Control_Right'].setValuesInStep(
    stepName='Step-GlobalTension', 
    u1=DISP_MAG * cos_val,
    u2=-DISP_MAG * sin_val
)

# 2. Left Cable: 좌하향 당기기 (-X, -Y)
# 좌우 대칭이므로 X부호만 반대
my_model.boundaryConditions['Disp_Control_Left'].setValuesInStep(
    stepName='Step-GlobalTension', 
    u1=-DISP_MAG * cos_val,
    u2=-DISP_MAG * sin_val
)

# 초기 가짜 응력 (수렴 보조용)
my_model.Stress(
    name='Initial_Stiffness',
    region=inst_memb.sets['All'],
    distributionType=UNIFORM,
    sigma11=700.0, sigma22=700.0, sigma33=0.0,
    sigma12=0.0, sigma13=0.0, sigma23=0.0
)

# Step 2: Buckle
if 'Step-Buckle' in my_model.steps: del my_model.steps['Step-Buckle']

my_model.BuckleStep(
    name='Step-Buckle',
    previous='Step-GlobalTension',
    numEigen=100,           # [핵심] 76개의 음수 모드를 건너뛰기 위해 100개 요청
    eigensolver=SUBSPACE,   # [핵심] Subspace 사용
    vectors=250,            # 모드 수의 2배 이상 (100 * 2.5)
    maxIterations=5000      # 끈기 있게 찾도록 횟수 증가
)

PERTURBATION_MAG = 0.01 

my_model.boundaryConditions['Disp_Control_Right'].setValuesInStep(
    stepName='Step-Buckle', 
    u1=PERTURBATION_MAG * cos_val,
    u2=-PERTURBATION_MAG * sin_val
)

my_model.boundaryConditions['Disp_Control_Left'].setValuesInStep(
    stepName='Step-Buckle', 
    u1=-PERTURBATION_MAG * cos_val,
    u2=-PERTURBATION_MAG * sin_val
)

my_model.boundaryConditions['BC_Stabilize_Z'].deactivate('Step-Buckle')

all_edges = inst_memb.edges
a.Set(name='All_Edges', edges=all_edges)

my_model.DisplacementBC(
    name='BC_Edges_Only_Z', 
    createStepName='Step-Buckle',  # Buckle Step에서 새로 생성
    region=a.sets['All_Edges'], 
    u3=SET
)

BUCKLE_ODB = 'Buckle_Analysis.odb'
LF_ODB = 'LF_Analysis.odb'
HF_ODB = 'HF_Postbuckle.odb'
FINAL_DISP = 0.00015


# -------------------------------------------------------------
# 7. 실행 로직 (HF/LF)
# -------------------------------------------------------------

cmd = ""
if fidelity == 'LF':
    run_job_safely('Buckle_Analysis')

    my_model.StaticStep(
        name='Step-HighTension',
        previous='Step-GlobalTension',
        nlgeom=ON
    )

    my_model.boundaryConditions['BC_Stabilize_Z'].setValuesInStep(
        stepName='Step-HighTension', 
        u3=0.0
    )

    for name, sign in [('Disp_Control_Right', 1), ('Disp_Control_Left', -1)]:
        my_model.boundaryConditions[name].setValuesInStep(
            stepName='Step-HighTension',
            u1=sign * FINAL_DISP * cos_val,
            u2=-FINAL_DISP * sin_val
        )

    my_model.fieldOutputRequests['F-Output-1'].setValues(
        variables=('S', 'U', 'COORD', 'EVOL'), 
        frequency=1
    )
    if 'Step-Buckle' in my_model.steps:
        del my_model.steps['Step-Buckle']

    run_job_safely('LF_Analysis')

    cmd = "abaqus python eval_abaqus.py %s LF" % LF_ODB

elif fidelity == 'HF':
    run_job_safely('Buckle_Analysis')

    my_model.StaticStep(name='Step-HighTension', previous='Step-GlobalTension', nlgeom=ON)
    my_model.boundaryConditions['BC_Stabilize_Z'].setValuesInStep(stepName='Step-HighTension', u3=0.0)
    
    my_model.StaticStep(
        name='Step-Postbuckle', 
        previous='Step-HighTension',  
        nlgeom=ON, 
        stabilizationMagnitude=0.0005, # 댐핑 계수 Galhofo Reference
        stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
        continueDampingFactors=False,
        adaptiveDampingRatio=0.05,
        initialInc=1e-5,
        minInc=1e-12,
        maxInc=0.05,
        maxNumInc=5000
    )

    my_model.boundaryConditions['BC_Stabilize_Z'].deactivate('Step-Postbuckle')

    # 꼭짓점 하중 증가
    # 초기  상태에서 최종 장력(7000 Pa)까지 증가 시나리오
    my_model.boundaryConditions['Disp_Control_Right'].setValuesInStep(
        stepName='Step-Postbuckle',
        u1=FINAL_DISP * cos_val, 
        u2=-FINAL_DISP * sin_val
    )
    
    # Left Cable
    my_model.boundaryConditions['Disp_Control_Left'].setValuesInStep(
        stepName='Step-Postbuckle',
        u1=-FINAL_DISP * cos_val,
        u2=-FINAL_DISP * sin_val
    )

    my_model.keywordBlock.synchVersions(storeNodesAndElements=False)
    
    imp_scale = THICKNESS / 10.0 # Galhofo Reference

    imp_text = "*IMPERFECTION, FILE=Buckle_Analysis, STEP=2\n1, %e\n2, %e" % (imp_scale, imp_scale)

    # 키워드 삽입 위치 찾기
    inserted = False
    for i, block in enumerate(my_model.keywordBlock.sieBlocks):
        # Step 정의 시작 부분 찾기
        if block.lower().strip().startswith('*step'):
            my_model.keywordBlock.insert(i-1, imp_text)
            inserted = True
            break
            
    if not inserted:
        my_model.keywordBlock.insert(len(my_model.keywordBlock.sieBlocks)-1, imp_text)

    # 3. HF (Post-buckling) 실행
    run_job_safely('HF_Postbuckle')
    
    # 결과 추출
    cmd = "abaqus python eval_abaqus.py %s %s HF" % (LF_ODB, HF_ODB)

try:
    print("Calling extraction script: %s" % cmd)
    # shell=True로 eval_abaqus.py 실행
    p = subprocess.Popen(cmd, shell=True)
    p.wait()
    
    # 추출 스크립트가 출력한 "RESULTS:..." 라인을 찾아 전달
    if os.path.exists('extraction.txt'):
        with open('extraction.txt', 'r') as f:
            print("RESULTS:" + f.read().strip())
    else:
        print("!!! ERROR: Extraction failed. 'extraction.txt' not found.") 

except Exception as err:
    print("Error during data extraction: %s" % str(err))
    sys.exit(1)