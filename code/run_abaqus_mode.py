#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""클램프 없는 모델에서 '고유모드'를 추출해 .fil 로 넘기는 전용 스크립트 (2026-09-22)

목적 — 2-모델 레시피의 '모드 소스'
    run_abaqus.py -- HF 는 클램프가 있는 모델이라 선형 좌굴 스펙트럼을 주지 못한다
    (실측: 시스템행렬 음수 고유값 598~2897개 / CONVERGED=0 / EIGENVALUES CANNOT BE FOUND).
    그래서 모드는 '클램프 없는' 이 모델에서 뽑아 .fil 로 넘긴다(Galhofo2022 의 2-모델 구조).

    [중요, 2026-09-22 실측] '*BUCKLE' 스텝은 .fil 출력이 금지된다:
      ClampFree_Buckle.dat:7025 -> "***WARNING: FILE OUTPUT IS NOT AVAILABLE FOR BUCKLING ANALYSIS"
    즉 좌굴 스텝에 *NODE FILE 을 넣어 모드를 .fil 로 보내는 경로는 구조적으로 막혀 있다.
    그래서 기본 모드 스텝은 '*FREQUENCY'(같은 base state 위의 주파수 추출, .fil 허용)다.
    λ(좌굴 고유치) 비교가 필요하면 MODE_STEP_TYPE='buckle' 로 바꾼다(그때는 .dat 의 MODE NO 표를 읽는다).

run_abaqus_cable.py 와 다른 점
    1) 모드 추출만 한다 — 포스트버클링/LF/HF 잡을 돌리지 않으므로 라이선스 소모가 최소다.
    2) 베이스 스테이트 파라미터를 HF 와 정렬한다(아래 '정렬' 주석 참조):
         코너 당김 DISP_GLOBAL = 5e-5 m / SIGMA0 = 500 Pa / 안정화 = 2e-4
       -> 모드 소스와 HF 가 '클램프 유무' 만 다른 비교가 된다(교란 제거).
       실패하면 케이블 런 값(1.8e-5 / 700 / 5e-4)으로 되돌린다 — 상수 3개만 수정.
    3) 잡 이름이 ClampFree_Buckle 이고, HF 는 이 파일을 IMPERFECTION_NAME 으로 스테이징한다
       (이름 분리 -> 자기 좌굴 잡의 0-모드 .fil 이 조용히 소비되는 사고를 차단, 원장 C-1).
    4) 메쉬/형상/BC/프리텐션 정의 줄은 run_abaqus_new.py(HF)·run_abaqus_buckle.py 와
       '문자 그대로' 같아야 한다 -> check_model_consistency.py 가 제출 전에 정적으로 검사한다.
       (모드는 HF 와 같은 노드 라벨에 정의되어야 *IMPERFECTION 으로 이식된다.)

실행
    abaqus cae noGUI=run_abaqus_mode.py
    (이어서) abaqus cae noGUI=run_abaqus.py -- HF <x_c> <d_c>

종료 코드
    0 = 모드 1개 이상 추출 / 1 = 0개(실패) 또는 잡 이상 -> 상위 스크립트가 판정 가능
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

# ---- 스크립트 위치(_HERE): base_state_probe.py / aba_imperfection.py 를 찾기 위해 ----
def _resolve_here():
    cand = []
    try:
        cand.append(os.path.dirname(os.path.abspath(__file__)))
    except NameError:
        pass
    cand.append(os.getcwd())
    if sys.argv and sys.argv[0]:
        cand.append(os.path.dirname(os.path.abspath(sys.argv[0])))
    for c in cand:
        if c and os.path.exists(os.path.join(c, 'aba_imperfection.py')):
            return c
    return os.getcwd()

_HERE = _resolve_here()
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from aba_imperfection import count_modes                                    # noqa: E402
print("[run_abaqus_mode] _HERE = %s" % _HERE)


def run_job_safely(job_name, model_name=None):
    """잡 제출 + 완료 대기. 제출 전 이전 산출물을 지워 stale .fil/.odb 오판을 막는다."""
    model_name = model_name or MODEL_NAME
    for _ext in ('odb', 'fil', 'sta', 'msg', 'lck', 'com', 'prt', 'sim', 'log',
                 'dat', 'res', 'abq', 'ipm', 'mdl', 'stt', 'cid'):
        _f = '%s.%s' % (job_name, _ext)
        if os.path.exists(_f):
            print("Removing stale artifact: %s" % _f)
            try:
                os.remove(_f)
            except OSError as _e:
                print("  (warning) could not remove %s: %s" % (_f, _e))
    if job_name in mdb.jobs:
        del mdb.jobs[job_name]
    job = mdb.Job(name=job_name, model=model_name)
    print("Submitting Job: %s" % job_name)
    job.writeInput(consistencyChecking=OFF)
    job.submit(consistencyChecking=OFF)
    job.waitForCompletion()
    time.sleep(1.0)
    if job.status == ABORTED or not os.path.exists(job_name + '.odb'):
        print("!!! ERROR: Job %s failed. Actual Status: %s" % (job_name, str(job.status)))
        return False
    print("Job %s finished (Status: %s)." % (job_name, str(job.status)))
    return True


def print_job_diag(job_name):
    """원인 판별용 핵심 줄만 콘솔에 찍는다(로그 전체를 붙여넣지 않아도 되도록)."""
    keys = ('NEGATIVE EIGENVALUES', 'CONVERGED', 'EIGENVALUES CANNOT BE FOUND',
            'HAS COMPLETED SUCCESSFULLY', 'THE ANALYSIS HAS BEEN COMPLETED',
            'STIFFNESS MATRIX IS SINGULAR', 'TOO MANY ATTEMPTS', '***ERROR',
            'USER INPUT PROCESSING', 'misplaced', 'NOT A VALID', 'UNKNOWN PARAMETER')
    for ext in ('msg', 'dat'):
        fn = '%s.%s' % (job_name, ext)
        if not os.path.exists(fn):
            print("[DIAG:%s] %s 없음" % (job_name, fn))
            continue
        size = os.path.getsize(fn)
        with open(fn, 'r', errors='replace') as f:
            # 키워드/입력처리 경고는 파일 '앞'에, 수렴 실패는 '뒤'에 있다 -> 양쪽을 본다
            head = f.read(600000)
            tail = ''
            if size > 2000000:
                f.seek(size - 2000000)
                tail = f.read()
        hits = ['[앞] ' + ln.strip() for ln in head.splitlines()
                if ln.strip() and any(k.lower() in ln.lower() for k in keys)]
        hits += ['[뒤] ' + ln.strip() for ln in tail.splitlines()
                 if ln.strip() and any(k.lower() in ln.lower() for k in keys)]
        print("[DIAG:%s] %s (%.0f KB) 핵심줄 %d개" % (job_name, fn, size / 1024.0, len(hits)))
        for ln in hits[-10:]:
            print("      | %s" % ln[:150])

    fn = '%s.fil' % job_name
    if os.path.exists(fn):
        with open(fn, 'rb') as f:
            data = f.read()
        n_eig = data.count(b'EIGENVALUE') + data.count(b'EIGEN')
        print("[DIAG:%s] .fil %.0f KB, 'EIGEN' 마커 %d개 (휴리스틱 — 0 이면 "
              "임퍼펙션 주입이 조용히 무시될 수 있다)"
              % (job_name, os.path.getsize(fn) / 1024.0, n_eig))
    else:
        print("[DIAG:%s] .fil 없음 -> 임퍼펙션 주입이 조용히 무시된다" % job_name)

    # 설계 규칙 (run_abaqus_cable.py 의 좌굴 스텝 주석에서 확인):
    #   '요청 고유값 수 > base state 의 음수 고유값 수' 여야 양수 좌굴모드가 subspace 창에 들어온다.
    #   (케이블 런: 52 < 100 -> 성공 / 클램프 C-route: 598~2897 > 100 -> 0모드)
    fn = '%s.msg' % job_name
    if os.path.exists(fn):
        with open(fn, 'r', errors='replace') as f:
            txt = f.read()
        key = 'SYSTEM MATRIX HAS'
        idx = txt.find(key)
        neg = None
        if idx >= 0:
            tok = txt[idx + len(key):].strip().split()
            if tok and tok[0].isdigit():
                neg = int(tok[0])
        if neg is None:
            print("[DIAG:%s] 'SYSTEM MATRIX HAS ... NEGATIVE EIGENVALUES' 를 못 찾음" % job_name)
        else:
            print("[DIAG:%s] 음수 고유값 %d개 vs 요청 %d개 -> %s"
                  % (job_name, neg, N_EIG_BUCKLE,
                     'OK (창에 양수 모드가 들어온다)' if N_EIG_BUCKLE > neg
                     else '!!! 요청 수를 늘려야 한다 (N_EIG_BUCKLE <= 음수 개수)'))


# =====================================================================
# 모델 정의 — HF 와 '문자 그대로' 같은 줄 (check_model_consistency.py 가 검사)
# =====================================================================
MODEL_NAME = 'SailModel_Triangle'
INSTANCE_NAME = 'MEMBRANE-1'

BASE = 20.0   # m
HEIGHT = 10.0 # m
THICKNESS = 5.0e-6
TARGET_STRESS = 7000.0 # Pa (참고용: 운용점 목표. 모드 추출에는 쓰이지 않는다)

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

# ---- 정렬: 프리텐션 정의 (HF 와 동일한 줄) ----
#   HF(run_abaqus.py) 의 코너 당김 5e-5 m·SIGMA0 500 Pa·안정화 2e-4 와 맞춘다.
#   실패하면 케이블 런 값으로 되돌린다: PRETENSION_SCALE=3.6(=1.8e-5/5e-6) / SIGMA0=700.0 /
#   안정화 0.0005. (그 경우 '클램프 유무' 외에 파라미터도 달라진다는 사실을 논문에 명시할 것.)
PRETENSION_SCALE = 10.0
DISP_GLOBAL = 0.000005 * PRETENSION_SCALE    # = 5e-5 m
SIGMA0 = 500.0                               # 초기 가짜 응력(수렴 보조) [Pa]
MODE_STABILIZATION = 0.0002                  # GlobalTension 안정화 계수 (HF 정렬)
PERTURBATION = 0.01                          # 좌굴 스텝 섭동 (케이블 런·HF 와 동일)
N_EIG_BUCKLE = 100                           # subspace 요청 고유값 수
BUCKLE_VECTORS = 250                         # subspace 기저 벡터 수
MODE_FREQ_NUM_EIGEN = 10                     # 주파수 스텝 요청 모드 수 (업스트림 원본과 동일)
N_MODE_FILE = 4                              # .fil 에 기록할 모드 수 (*NODE FILE, LAST MODE)
INSERT_NODE_FILE = True                      # 좌굴모드 .fil 기록을 키워드로 '명시 요청'할지.
                                             #   케이블 런(성공)에는 이 요청이 없다 -> 2026-09-22 실패
                                             #   (모드 0)의 비물리적 후보 1순위. False 래더(L1)로 검증.
JOB_NAME = 'ClampFree_Buckle'
BUCKLE_STEP_NO = 2                           # 이 모델의 스텝 순서: 1=GlobalTension 2=Buckle

# -------------------------------------------------------------
# 2. 모델 초기화 및 재질
# -------------------------------------------------------------
if MODEL_NAME in mdb.models: del mdb.models[MODEL_NAME]
my_model = mdb.Model(name=MODEL_NAME)

# (1) 멤브레인 재질 (Kapton)
mat = my_model.Material(name='Kapton')
mat.Density(table=((1420.0,),))
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
# 4. 파트 생성: 케이블
# -------------------------------------------------------------
def create_cable_part(name, length):
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

# -------------------------------------------------------------
# 5. 어셈블리 및 연결
# -------------------------------------------------------------
a = my_model.rootAssembly
a.DatumCsysByDefault(CARTESIAN)
inst_memb = a.Instance(name=INSTANCE_NAME, part=p, dependent=ON)


def create_rigid_patch(name, coord, radius=0.2):
    """강체 패치 = 기존 노드에 Coupling(KINEMATIC). face partition 을 하지 않으므로
    메쉬(노드 좌표/라벨)를 바꾸지 않는다 -> HF 와 노드 라벨이 일치한다."""
    rp = a.ReferencePoint(point=coord)
    rp_key = a.referencePoints[rp.id]
    rp_region = regionToolset.Region(referencePoints=(rp_key,))
    nodes = inst_memb.nodes.getByBoundingSphere(center=coord, radius=radius)
    if not nodes:
        nodes = (inst_memb.nodes.getClosest(coordinates=coord),)
    patch_set = a.Set(name=name+'_Nodes', nodes=nodes)
    my_model.Coupling(
        name=name+'_Coupling', controlPoint=rp_region, surface=patch_set,
        influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
        u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON
    )
    return rp, rp_region


def connect_cable(name, part, sail_corner, vector_dir, radius=1e-4):
    inst_name = 'Inst_' + name
    inst = a.Instance(name=inst_name, part=part, dependent=ON)
    base_vec = np.array([1.0, 0.0, 0.0])
    target_vec = np.array(vector_dir)
    target_vec = target_vec / np.linalg.norm(target_vec)
    rot_angle_deg = np.degrees(np.arctan2(target_vec[1], target_vec[0]))
    a.rotate(instanceList=(inst_name,), axisPoint=(0,0,0), axisDirection=(0,0,1), angle=rot_angle_deg)
    a.translate(instanceList=(inst_name,), vector=sail_corner)
    node_start = inst.nodes.getByBoundingSphere(center=sail_corner, radius=radius)
    if len(node_start) == 0:
        print("Warning: Node not found by sphere, trying closest for " + name)
        node_start = inst.nodes.getClosest(coordinates=sail_corner)
    region_start = regionToolset.Region(nodes=node_start)
    cable_len = LEN_TOP if 'Top' in name else LEN_BOT
    end_coord = (sail_corner[0] + target_vec[0]*cable_len, sail_corner[1] + target_vec[1]*cable_len, 0.0)
    node_end = inst.nodes.getByBoundingSphere(center=end_coord, radius=radius)
    region_end = regionToolset.Region(nodes=node_end)
    return region_start, region_end


# 1. 메쉬 생성 (노드를 찾기 전에 필요)
p.seedPart(size=BASE/200.0, deviationFactor=0.1)
p.setMeshControls(regions=p.faces, elemShape=QUAD_DOMINATED, technique=FREE, algorithm=MEDIAL_AXIS)
elemTypeQuad = ElemType(elemCode=S4, elemLibrary=STANDARD)
elemTypeTri = ElemType(elemCode=S3, elemLibrary=STANDARD)
p.setElementType(regions=(p.faces,), elemTypes=(elemTypeQuad, elemTypeTri))
p.generateMesh()
a.regenerate()

# 2. 꼭짓점 강체 패치 3개 (클램프 패치는 만들지 않는다 — 이 모델의 요점)
rp1_obj, rp1_reg = create_rigid_patch('Top', V1, radius=0.2)
rp2_obj, rp2_reg = create_rigid_patch('Right', V2, radius=0.2)
rp3_obj, rp3_reg = create_rigid_patch('Left', V3, radius=0.2)

# 3. 케이블 연결: 정점은 위로, 두 아래 모서리는 각 케이블 축(28.6도) 방향
start_c1, end_c1 = connect_cable('Cable_Top', p_cable_top, V1, (0.0, 1.0, 0.0))
start_c2, end_c2 = connect_cable('Cable_Right', p_cable_bot, V2, (cos_val, -sin_val, 0.0))
start_c3, end_c3 = connect_cable('Cable_Left', p_cable_bot, V3, (-cos_val, -sin_val, 0.0))

# 4. Tie (RP <-> Cable Start)
a.Set(name='RP_Top_Set', referencePoints=(a.referencePoints[rp1_obj.id],))
a.Set(name='RP_Right_Set', referencePoints=(a.referencePoints[rp2_obj.id],))
a.Set(name='RP_Left_Set', referencePoints=(a.referencePoints[rp3_obj.id],))
my_model.Tie(name='Tie_Top', main=a.sets['RP_Top_Set'], secondary=start_c1, positionToleranceMethod=COMPUTED)
my_model.Tie(name='Tie_Right', main=a.sets['RP_Right_Set'], secondary=start_c2, positionToleranceMethod=COMPUTED)
my_model.Tie(name='Tie_Left', main=a.sets['RP_Left_Set'], secondary=start_c3, positionToleranceMethod=COMPUTED)

# -------------------------------------------------------------
# 6. 경계 조건
# -------------------------------------------------------------
# Step 1: Global Tension (= 좌굴의 base state)
my_model.StaticStep(
    name='Step-GlobalTension', previous='Initial', nlgeom=ON,
    stabilizationMagnitude=0.0002,
    stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
    continueDampingFactors=False,
    adaptiveDampingRatio=0.05,
    initialInc=0.0001, minInc=1e-8, maxNumInc=1000
)

my_model.DisplacementBC(
    name='BC_Anchor_Top',
    createStepName='Initial',
    region=end_c1,
    u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET
)
my_model.DisplacementBC(
    name='BC_Right_Z', createStepName='Initial', region=end_c2,
    u3=SET, ur1=SET, ur2=SET, ur3=SET
)
my_model.DisplacementBC(
    name='BC_Left_Z', createStepName='Initial', region=end_c3,
    u3=SET, ur1=SET, ur2=SET, ur3=SET
)
my_model.DisplacementBC(
    name='Disp_Control_Right', createStepName='Initial', region=end_c2, u1=SET, u2=SET
)
my_model.DisplacementBC(
    name='Disp_Control_Left', createStepName='Initial', region=end_c3, u1=SET, u2=SET
)
# 평탄 강제(강제 평탄) — 좌굴 스텝에서 해제한다(케이블 런·HF 와 동일한 순서)
my_model.DisplacementBC(
    name='BC_Stabilize_Z', createStepName='Step-GlobalTension',
    region=inst_memb.sets['All'], u3=SET
)

my_model.boundaryConditions['Disp_Control_Right'].setValuesInStep(
    stepName='Step-GlobalTension', u1=DISP_GLOBAL * cos_val, u2=-DISP_GLOBAL * sin_val
)
my_model.boundaryConditions['Disp_Control_Left'].setValuesInStep(
    stepName='Step-GlobalTension', u1=-DISP_GLOBAL * cos_val, u2=-DISP_GLOBAL * sin_val
)

# 초기 가짜 응력 (수렴 보조)
my_model.Stress(
    name='Initial_Stiffness',
    region=inst_memb.sets['All'],
    distributionType=UNIFORM,
    sigma11=SIGMA0, sigma22=SIGMA0, sigma33=0.0,
    sigma12=0.0, sigma13=0.0, sigma23=0.0
)

# Step 2: 모드 추출 스텝
#   실측(ClampFree_Buckle.dat:7025): '*BUCKLE' 스텝은 .fil 출력이 금지된다
#     "***WARNING: FILE OUTPUT IS NOT AVAILABLE FOR BUCKLING ANALYSIS"
#   -> *NODE FILE 로 모드를 .fil 에 넣는 경로가 구조적으로 막혀 있다(요청이 무의미).
#   같은 base state 위의 '주파수 추출'은 .fil 출력이 허용된다 -> 모드 소스는 이쪽을 쓴다.
#   λ(좌굴 고유치)만 필요하면 MODE_STEP_TYPE='buckle' -> .dat 의 MODE NO 표를 읽는다.
MODE_STEP_TYPE = 'frequency'   # 'frequency' = 모드 소스(.fil 기록 가능) | 'buckle' = λ 비교용
if MODE_STEP_TYPE not in ('frequency', 'buckle'):
    raise RuntimeError("MODE_STEP_TYPE 은 'frequency' 또는 'buckle' 여야 합니다 (현재 %r)"
                       % (MODE_STEP_TYPE,))
MODE_STEP_NAME = 'Step-Mode' if MODE_STEP_TYPE == 'frequency' else 'Step-Buckle'
if MODE_STEP_NAME in my_model.steps:
    del my_model.steps[MODE_STEP_NAME]
if MODE_STEP_TYPE == 'frequency':
    # 업스트림 원본과 같은 설정(numEigen=10, LANCZOS). 원본 주석 그대로:
    #   "BuckleStep을 사용하는 것이 정석이고 옳으나, 매우 얇은 solar-sail 자체의 불안정성에 의해
    #    negative eigenvalue만 찾는 경우가 대부분이라, mfbo 적용을 위해 Abaqus FrequencyStep으로 대체."
    # CAE FrequencyStep 은 vectors/maxIterations 를 받지 않는 형태가 원본에서 검증됐다.
    my_model.FrequencyStep(
        name=MODE_STEP_NAME,
        previous='Step-GlobalTension',
        numEigen=MODE_FREQ_NUM_EIGEN,
        eigensolver=LANCZOS
    )
else:
    my_model.BuckleStep(
        name=MODE_STEP_NAME,
        previous='Step-GlobalTension',
        numEigen=N_EIG_BUCKLE,
        eigensolver=SUBSPACE,
        vectors=BUCKLE_VECTORS,
        maxIterations=5000
    )
    # 좌굴 스텝에만 있는 '하중 패턴'(섭동) -> λ 를 그 패턴 기준으로 얻는다.
    # (주파수 추출에는 하중 패턴 개념이 없으므로 이 섭동을 주지 않는다)
    my_model.boundaryConditions['Disp_Control_Right'].setValuesInStep(
        stepName=MODE_STEP_NAME, u1=PERTURBATION * cos_val, u2=-PERTURBATION * sin_val
    )
    my_model.boundaryConditions['Disp_Control_Left'].setValuesInStep(
        stepName=MODE_STEP_NAME, u1=-PERTURBATION * cos_val, u2=-PERTURBATION * sin_val
    )
my_model.boundaryConditions['BC_Stabilize_Z'].deactivate(MODE_STEP_NAME)
a.Set(name='All_Edges', edges=inst_memb.edges)
my_model.DisplacementBC(
    name='BC_Edges_Only_Z',
    createStepName=MODE_STEP_NAME,
    region=a.sets['All_Edges'],
    u3=SET
)

# 필드출력: base_state_probe 가 면적가중(EVOL)·반력(RF)을 읽는다
my_model.fieldOutputRequests['F-Output-1'].setValues(
    variables=('S', 'E', 'U', 'RF', 'COORD', 'EVOL'),
    frequency=1
)

# -------------------------------------------------------------
# 7. *NODE FILE 삽입 -> 좌굴모드를 .fil 에 기록 (없으면 *IMPERFECTION 이 조용히 무시된다)
#    위치: *BUCKLE 키워드 '직후' = 좌굴 스텝 '안' (스텝 밖으로 밀리면 misplaced 로 죽는다)
# -------------------------------------------------------------
if INSERT_NODE_FILE:
    my_model.keywordBlock.synchVersions(storeNodesAndElements=False)
    _noderef = '*NODE FILE, GLOBAL=YES, LAST MODE=%d\nU' % N_MODE_FILE
    _found = False
    for _i, _b in enumerate(my_model.keywordBlock.sieBlocks):
        if _b.strip().upper().startswith(('*BUCKLE', '*FREQUENCY')):
            my_model.keywordBlock.insert(_i + 1, _noderef)
            _found = True
            break
    print("[run_abaqus_mode] *NODE FILE 삽입=%s (LAST MODE=%d)" % (_found, N_MODE_FILE))
else:
    print("[run_abaqus_mode] *NODE FILE 삽입 생략 (INSERT_NODE_FILE=False)"
          " -> *BUCKLE 자체 기록에 의존. HF 쪽 stage() 가 모드 0개를 잡으므로 조용한 주입은 없다.")

# -------------------------------------------------------------
# 8. 실행
# -------------------------------------------------------------
print("=" * 74)
print("모드 소스 런: %s  (스텝=%s/%s, 코너 당김 %.3e m / SIGMA0 %.1f Pa / 안정화 %.4f)"
      % (JOB_NAME, MODE_STEP_TYPE, MODE_STEP_NAME, DISP_GLOBAL, SIGMA0, MODE_STABILIZATION))
print("  모델: 클램프 없음 / 케이블 3개 / 꼭짓점 앵커 / 두 아래 모서리 구동 (논문 좌굴 모델 구성)")
print("=" * 74)

ok = run_job_safely(JOB_NAME)
print_job_diag(JOB_NAME)

# base state 실측(클램프 없는 상태의 압축 영역 비율 등) — 라이선스 추가 소모 없음
try:
    _probe = os.path.join(_HERE, 'base_state_probe.py')
    if os.path.exists(_probe) and os.path.exists(JOB_NAME + '.odb'):
        print("[MODE] base state 측정: %s.odb" % JOB_NAME)
        subprocess.call('abaqus python "%s" "%s" Step-GlobalTension'
                        % (_probe, JOB_NAME + '.odb'), shell=True)
    else:
        print("[MODE] base state 측정 건너뜀 (probe=%s)" % os.path.exists(_probe))
except Exception as _e:
    print("[MODE] base state 측정 실패(무시): %s" % _e)

n_modes = count_modes(JOB_NAME + '.dat', JOB_NAME + '.msg')
print("=" * 74)
if not ok or n_modes <= 0:
    print("RESULT:MODE_FAIL — 모드 %d개 (job_ok=%s, 스텝=%s)." % (n_modes, ok, MODE_STEP_TYPE))
    print("  래더(1줄씩, 1회 ~61초):")
    print("    L1  MODE_STEP_TYPE='frequency' (기본값) — .fil 출력이 허용되는 유일한 스텝")
    print("    L2  SIGMA0 = 700.0 + MODE_STABILIZATION = 0.0005   (수치 보조, 더 보수적)")
    print("    L3  PRETENSION_SCALE = 3.6   (코너 당김 1.8e-5 = 케이블 런과 같은 인장)")
    print("    L4  N_EIG_BUCKLE = 200 + BUCKLE_VECTORS = 500   (요청 수 > 음수 고유값 수)")
    print("  라이선스 0: 위 [DIAG] 의 '음수 고유값 N개 vs 요청 M개' 판정을 먼저 본다.")
    sys.exit(1)

print("RESULT:MODE_OK — 모드 %d개 계산. 파일: %s" % (n_modes, os.path.abspath(JOB_NAME + '.fil')))
print("  스텝=%s(%s) / *NODE FILE 요청=%s (주파수 스텝이므로 .fil 기록이 허용된다)"
      % (MODE_STEP_TYPE, MODE_STEP_NAME, INSERT_NODE_FILE))
print("  [기본 경로] run_abaqus.py 가 이 .fil 을 IMPERFECTION_NAME=%s 로 스테이징해"
      " *IMPERFECTION, FILE= 로 주입한다" % 'ClampFree_Buckle')
print("  다음 단계: abaqus cae noGUI=run_abaqus.py -- HF <x_c> <d_c>")
print("  [대체 경로] .fil 에 모드가 0개면(위 DIAG 가 MODE_FAIL) ODB 모드표로 우회 - code\\ 에서:")
print("    abaqus python aba_mode_from_odb.py %s %s modes_ClampFree_Buckle.txt %d"
      % (JOB_NAME + '.odb', MODE_STEP_NAME, N_MODE_FILE))
print("    -> run_abaqus.py 에서 IMPERFECTION_MODE='odb_table' 로 바꾼다")
print("=" * 74)
