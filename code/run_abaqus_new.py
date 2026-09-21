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
import math
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
RUN_DIR_NAME = "aba"
NUMCPUS = 4
BASE_PROBE = True
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
    dump_job_diag(job_name)
    
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


def dump_job_diag(job_name):
    """실패 진단에 필요한 것만 파일 하나로 모아 둔다 (사용자가 그 파일만 보내면 되도록).

    내용: ① .sta 마지막 25줄 ② .msg/.dat 의 원인 판별 키 줄 ③ 완주 판정 결과.
    파일: <RUN_DIR>/<job>.diag.txt  (예: aba/HF_Postbuckle.diag.txt)
    """
    out = ['===== job_completed_ok = %s =====' % job_completed_ok(job_name)]
    sta = '%s.sta' % job_name
    if os.path.exists(sta):
        with open(sta, 'r', errors='replace') as f:
            out.append('===== %s (last 25 lines) =====' % sta)
            out.extend(f.read().splitlines()[-25:])
    KEYS = ('***ERROR', '***WARNING', 'TOO MANY', 'DISTORT', 'NEGATIVE EIGENVALUE',
            'HAS COMPLETED', 'NOT BEEN COMPLETED', 'EXCESSIVE', 'CUT BACK',
            'CANNOT BE', 'ATTEMPT NUMBER  2',
            # 2026-09-21: 아래 5개가 빠져 있어 원인 판별이 한 라운드 지연되었다
            'CONSTANT DAMPING', 'OVERCONSTRAINT', 'INACTIVE DOF',
            'ALLSDTOL', 'SEVERE ELEMENT')
    for ext in ('msg', 'dat'):
        fn = '%s.%s' % (job_name, ext)
        if not os.path.exists(fn):
            continue
        out.append('===== %s (key lines, max 500) =====' % fn)
        n_hit = 0
        with open(fn, 'r', errors='replace') as f:
            for _ln, _line in enumerate(f, 1):
                _u = _line.upper()
                if any(_k in _u for _k in KEYS):
                    out.append('%d: %s' % (_ln, _line.rstrip()))
                    n_hit += 1
                    if n_hit >= 500:
                        out.append('... (truncated at 500)')
                        break
    if len(out) <= 2:
        out.append('(no .sta/.msg/.dat found)')
    dst = os.path.join(_RUN, '%s.diag.txt' % job_name)
    # ASCII 내용 + 명시적 UTF-8 저장: 한글 라벨을 쓰면 Windows 기본 인코딩(CP949)으로
    # 저장되어 다른 도구에서 읽히지 않는다 (2026-09-21 실제 발생).
    import io as _io
    with _io.open(dst, 'w', encoding='utf-8', errors='replace') as f:
        f.write('\n'.join(out) + '\n')
    print('[run_abaqus_new] 진단 요약 저장: %s (%d줄)' % (dst, len(out)))
    return dst


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
THICKNESS = 5.0e-6
TARGET_STRESS = 7000.0 # Pa   # 목표 운용점 - 실제 도달 응력 미검증(측정 필요)

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

# ---- 사전 장력(prestrain) 캘리브레이션 ----

PRETENSION_SCALE = 10.0          # 프리텐션 변위 = 5e-6 m * 이 값 = 5e-5 m
# R-13 실측(2026-09-21): PRETENSION_SCALE=10 -> 평균 면내응력 2122 Pa = 목표 7000 Pa 의 0.303배.
DISP_GLOBAL = 0.000005 * PRETENSION_SCALE    # 운용점: 코너 당김 5e-5 m
CLAMP_PULL = DISP_GLOBAL * d_c
GLOBAL_FINAL = 0.0001 * PRETENSION_SCALE     # 최종 하중: 코너 당김 1e-3 m
CLAMP_FINAL = GLOBAL_FINAL * d_c
# ---- trigger 파라미터 ----
TRIG_MAG = THICKNESS * 0.1       # 0.1t = 5e-7 m (Galhofo 관행)
TRIG_MARGIN = 0.05               # 삼각형 빗변에서 띄울 여유 (모델 폭 W 대비 비율)
TRIG_ON = True
STAB = 0.0002
ADAPT_DAMP_MAX = 0.15
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
#   HF 전용으로 만든다 -> LF 체인(GlobalTension->ClampTension->HighTension)은 기존과 동일하게 유지되어 LF 지표가 변하지 않는다.
if fidelity == 'HF':
    my_model.StaticStep(
        name='Step-Trigger',
        previous='Step-GlobalTension',
        nlgeom=ON,
        stabilizationMagnitude=STAB,      # 코드 상수 STAB (기본 2e-4)
        stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
        initialInc=0.1, minInc=1e-8, maxInc=1.0, maxNumInc=50
    )
_PREV_CLAMP = 'Step-Trigger' if fidelity == 'HF' else 'Step-GlobalTension'

# Step Clamp Tension : 클램프에 변위 가하기
#   u3 는 Step-Trigger 에서 해제되므로 이 스텝부터가 사실상 포스트버클링 구간이다 (프리텐션 상태가 이미 분기 위 -> 해제 즉시 주름 발생).
#   물리적 스테이징 : GlobalTension 이 '운용 프리텐션' 수준(DISP_GLOBAL)을 만들고,
#         이 스텝에서 설계변수 클램프(d_c)를 결합하고, Postbuckle 이 최종 하중까지 올린다.
#   수치적 스테이징 : 분기 핵생성과 20배 램프를 한 스텝에 몰면 수렴이 나빠진다
#   LF 도 동일 스텝 구조를 쓴다 (하중경로 동일화 -> LF->HF 보정이 '주름 효과'만 학습).

my_model.StaticStep(
    name='Step-ClampTension',
    previous=_PREV_CLAMP,
    nlgeom=ON,
    stabilizationMagnitude=STAB,      # 코드 상수: 2e-4 (주름을 죽이지 않는 기저값)
    stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
    continueDampingFactors=False,     # 스텝마다 감쇠 초기화
    adaptiveDampingRatio=ADAPT_DAMP_MAX,        # 적응 감쇠: 수렴이 어려울 때만 Abaqus 가 자동으로 키운다
    initialInc=0.0001, minInc=1e-8, maxNumInc=20000   # 2026-09-21: 1000 소진 -> 말단 속도 1.25e-4/inc 기준 ~2300 필요
)


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

    # ---- trigger : 삼각형 내부를 해석적 사인 패턴으로 교란 (명시적 seed) ----
    _mn = inst_memb.nodes
    _x0 = min(n.coordinates[0] for n in _mn); _x1 = max(n.coordinates[0] for n in _mn)
    _y0 = min(n.coordinates[1] for n in _mn); _y1 = max(n.coordinates[1] for n in _mn)
    _W, _H = _x1 - _x0, _y1 - _y0
    _xc = 0.5 * (_x0 + _x1)

    def _trig_sign(c):
        _t = (c[1] - _y0) / _H                    # 0 = 밑변, 1 = 꼭짓점
        if not (0.30 < _t < 0.65):
            return None                           # 밑변/꼭짓점에서 충분히 안쪽
        _half = 0.5 * _W * (1.0 - _t)             # 삼각형 반폭
        if abs(c[0] - _xc) > _half - TRIG_MARGIN * _W:
            return None                           # 빗변에서 TRIG_MARGIN*W 이상 안쪽
        return math.sin(math.pi * (c[0] - _x0) / _W) * math.sin(2.0 * math.pi * (c[1] - _y0) / _H) >= 0.0

    _n_p = [n for n in _mn if _trig_sign(n.coordinates) is True]
    _n_m = [n for n in _mn if _trig_sign(n.coordinates) is False]
    _n_int = _n_p + _n_m
    print("[run_abaqus_new] 노드 범위 x[%.4f, %.4f] y[%.4f, %.4f] W=%.4f H=%.4f"
          % (_x0, _x1, _y0, _y1, _W, _H))
    print("[run_abaqus_new] trigger 내부 노드 %d / 전체 %d  (+, -)=(%d, %d)  mag=%.3e m  on=%s"
          % (len(_n_int), len(_mn), len(_n_p), len(_n_m), TRIG_MAG, TRIG_ON))
    if TRIG_ON and len(_n_int) < 50:
        raise RuntimeError(
            "trigger 내부 노드가 %d개뿐입니다 (전체 %d). 빈/과소 노드셋에 BC 를 만들면 "
            "trigger 가 조용히 무효가 됩니다 - 영역 조건(_trig_sign)을 확인하세요."
            % (len(_n_int), len(_mn)))
    _TRIG_OK = False
    try:
        a.Set(name='NS_TRIG_P', nodes=_n_p)
        a.Set(name='NS_TRIG_M', nodes=_n_m)
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

LF_ODB = 'LF_Analysis.odb'
HF_ODB = 'HF_Postbuckle.odb'

try:
    my_model.HistoryOutputRequest(
        name='H-Energy', createStepName='Step-GlobalTension',
        variables=('ALLIE', 'ALLSD', 'ALLKE'))
    print("[run_abaqus_new] H-Energy 추가: ALLIE/ALLSD/ALLKE (감쇠 가격 측정)")
except Exception as _e:
    print("[run_abaqus_new] (warning) H-Energy 추가 실패: %s" % _e)

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

    # Static, General 포스트버클링 수행    
    my_model.StaticStep(
        name='Step-Postbuckle', 
        previous='Step-ClampTension',
        nlgeom=ON, 
        stabilizationMagnitude=STAB,      # 코드 상수 STAB (기본 2e-4) 
        stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
        continueDampingFactors=False,
        adaptiveDampingRatio=ADAPT_DAMP_MAX,
        initialInc=1e-4,
        minInc=1e-8,
        maxInc=0.1,
        maxNumInc=20000
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
