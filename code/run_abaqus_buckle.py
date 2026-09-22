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
    job 이름은 인자에서 자동 생성한다: Buckle_xc<NNN>_d<NNN>um
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



# ============================================================================
# 모델 상수 — run_abaqus_new.py 의 HF 모델과 같은 값이어야 한다.
#   검사: python code/check_model_consistency.py
# ============================================================================
MODEL_PREFIX = 'SailModel_Buckle'
INSTANCE_NAME = 'MEMBRANE-1'

BASE = 20.0   # m
HEIGHT = 10.0 # m
THICKNESS = 5.0e-6

V1 = (BASE/2.0, HEIGHT, 0.0) # Top
V2 = (BASE, 0.0, 0.0)        # Right
V3 = (0.0, 0.0, 0.0)         # Left

def clamp_coord_L(x): return (10-10*x, 10-10*x, 0)
def clamp_coord_R(x): return (10+10*x, 10-10*x, 0)

# ---- 프리텐션 기본값 (CLI 두 번째 인자로 덮어쓴다) ----
PRETENSION_SCALE = 10.0
DISP_GLOBAL = 0.000005 * PRETENSION_SCALE    # 기본 코너 당김 5e-5 m

# ---- 좌굴 스텝 (run_abaqus_cable.py 에서 완주가 확인된 설정과 동일) ----
PERTURBATION = 0.01     # m — 좌굴 스텝의 prescribed 변위(증분 응력 -> K_delta)
N_EIG_BUCKLE = 100      # 추출 요청 고유값 수
BUCKLE_VECTORS = 250    # subspace 기저 벡터 수 (요청 수의 2.5배)
BUCKLE_MAXITER = 5000
BUCKLE_SOLVER = 'SUBSPACE'   # 'SUBSPACE' | 'LANCZOS' — 제어 흐름용 문자열
BUCKLE_BLOCK_SIZE = 8           # LANCZOS 전용
BUCKLE_MIN_EIGEN = 0.0          # LANCZOS 전용 (음의 고유값도 보고 싶으면 -1e30 등)
BUCKLE_MAX_EIGEN = None         # LANCZOS 전용 (None 이면 인자를 아예 넘기지 않는다)

# Abaqus API 는 **심볼릭 상수**를 요구한다. 문자열을 그대로 넘기면
#   "eigensolver; found string, expecting SUBSPACE, LANCZOS or AMS"
# 로 즉시 죽는다(2026-09-22 실측: 전 케이스 BUILD_FAIL). 여기서 변환해서 넘긴다.
EIGENSOLVER_CONST = {'SUBSPACE': SUBSPACE, 'LANCZOS': LANCZOS}[BUCKLE_SOLVER]

SIGMA0 = 500.0          # 초기응력 [Pa] — 수렴 보조 (run_abaqus_new.py 와 동일)


# ============================================================================
# 인자 — 한 번에 한 케이스만 받는다. 스윕 기능은 없다.
#   abaqus cae noGUI=run_abaqus_buckle.py -- <x_c> [disp_m]
# ============================================================================
# Abaqus CAE 러너가 scripts 에 넘기는 sys.argv 전문(2026-09-22 Windows, Abaqus 2026 실측):
#   ['C:\...\win_b64\code\bin\ABQcaeK.exe', '-cae', '-noGUI', 'run_abaqus_buckle.py',
#    '-academic', 'RESEARCH', '-tmpdir', 'C:\Users\...\Temp', '-lmlog', 'ON', '0.5', '1e-3']
# 즉 (a) 러너가 자기 플래그와 **그 값**(RESEARCH, Temp 경로, ON)을 앞에 붙이고,
#    (b) 셸의 '--' 구분자는 스크립트까지 전달되지 않으며,
#    (c) 사용자 인자는 **맨 뒤에 숫자로** 온다.
# => 뒤에서부터 훑어 러너 토큰은 건너뛰고 숫자가 끊길 때까지 모은다.
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
    if not (1 <= len(nums) <= 2):
        raise RuntimeError(
            'Expected 1 or 2 trailing numeric arguments (<x_c> [disp_m]), got %r '
            '(first non-numeric token from the end: %r).\n'
            'Usage: abaqus cae noGUI=run_abaqus_buckle.py -- <x_c> [disp_m]'
            % (nums, stop))
    return nums[0], (nums[1] if len(nums) > 1 else None)


x_c, _disp_arg = parse_args(sys.argv)
DISP = DISP_GLOBAL if _disp_arg is None else _disp_arg

if not (0.03 <= x_c <= 0.95):
    # 원장 M-8: x_c < 0.028 이면 클램프 패치가 정점 패치와 겹치고, x_c ~ 1 이면
    # 모서리 패치와 겹친다 -> 조용히 다른 모델이 된다. 경고만 찍고 진행한다.
    print("%s WARNING: x_c=%g is outside the safe band (0.03, 0.95) — "
          "the clamp patch may overlap a vertex patch (ledger M-8)." % (TAG, x_c))
if DISP <= 0.0:
    raise RuntimeError('disp_m must be > 0 (got %r).' % (DISP,))

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

    # ABORTED가 아니면서, ODB 파일이 실제로 존재하면 성공으로 간주
    if job.status == ABORTED or not os.path.exists(odb_file):
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

    for ln in out:
        print('      | %s' % ln[:170])

    dst = os.path.join(_RUN, '%s.diag.txt' % job_name)
    import io as _io
    with _io.open(dst, 'w', encoding='utf-8', errors='replace') as f:
        f.write('\n'.join(out) + '\n')
    print('%s diag saved: %s (%d lines)' % (TAG, dst, len(out)))
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
        influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
        u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON
    )
    return rp, rp_region



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
    my_model.HomogeneousShellSection(name='Section-Membrane', material='Kapton', thickness=THICKNESS)

    # [제외] 케이블 재질(Kevlar) / TrussSection — 케이블이 없으므로 만들지 않는다.

    # 파트 생성: 멤브레인 — 동일
    s = my_model.ConstrainedSketch(name='triangle_profile', sheetSize=BASE*2)
    s.Line(point1=V3[:2], point2=V2[:2])
    s.Line(point1=V2[:2], point2=V1[:2])
    s.Line(point1=V1[:2], point2=V3[:2])
    p = my_model.Part(name='Membrane', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseShell(sketch=s)
    p.SectionAssignment(region=p.Set(faces=p.faces, name='All'), sectionName='Section-Membrane')

    # [제외] create_cable_part(...) 5개 — 케이블 part 자체를 만들지 않는다.

    a = my_model.rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    inst_memb = a.Instance(name=INSTANCE_NAME, part=p, dependent=ON)

    # ---- 메쉬: run_abaqus_new.py 와 완전히 동일 ----
    p.seedPart(size=BASE/200.0, deviationFactor=0.1) # 약 1만개
    p.setMeshControls(regions=p.faces, elemShape=QUAD_DOMINATED, technique=FREE, algorithm=MEDIAL_AXIS)
    # ★ 요소 타입은 run_abaqus_new.py(HF) 와 반드시 같아야 한다.
    #   이유: 좌굴 모드는 HF 와 **같은 노드**에 정의되어야 *IMPERFECTION 으로 이식된다.
    #   HF 쪽 요소 타입을 바꾸면 이 두 줄도 함께 바꿔야 한다.
    #   현재값은 HF 기준 S4/S3 (커밋 6e0ef6e "cable deformation compatibility" 이후).
    #   참고: Galhofo 검증모델은 S4R(s4R) 을 썼다 — 남은 차이는 요소 종류 하나다.
    elemTypeQuad = ElemType(elemCode=S4, elemLibrary=STANDARD)
    elemTypeTri = ElemType(elemCode=S3, elemLibrary=STANDARD)
    p.setElementType(regions=(p.faces,), elemTypes=(elemTypeQuad, elemTypeTri))
    p.generateMesh()
    a.regenerate()

    # 꼭짓점 RP — 동일 (radius=0.2)
    rp1_obj, rp1_reg = create_rigid_patch(a, inst_memb, 'Top', V1, radius=0.2)
    rp2_obj, rp2_reg = create_rigid_patch(a, inst_memb, 'Right', V2, radius=0.2)
    rp3_obj, rp3_reg = create_rigid_patch(a, inst_memb, 'Left', V3, radius=0.2)

    # 클램프 RP : 우측 빗변 중점 (15, 5), 좌측 빗변 중점 (5, 5) — 동일
    rp_cl_obj, rp_cl_reg = create_rigid_patch(a, inst_memb, 'CL', V_CL, radius=0.2)
    rp_cr_obj, rp_cr_reg = create_rigid_patch(a, inst_memb, 'CR', V_CR, radius=0.2)

    # [제외] connect_cable(...) 5개 — 케이블 배치/부착 없음.
    # 케이블 버전에서 케이블은 'RP <-> 구동 끝단' 사이의 하중 전달 로드였다.
    # 따라서 케이블을 뺀 뒤에는 RP 를 직접 구동/구속하면 하중 경로가 동일하다.

    a.Set(name='RP_Top_Set', referencePoints=(a.referencePoints[rp1_obj.id],))
    a.Set(name='RP_Right_Set', referencePoints=(a.referencePoints[rp2_obj.id],))
    a.Set(name='RP_Left_Set', referencePoints=(a.referencePoints[rp3_obj.id],))
    a.Set(name='RP_CL_Set', referencePoints=(a.referencePoints[rp_cl_obj.id],))
    a.Set(name='RP_CR_Set', referencePoints=(a.referencePoints[rp_cr_obj.id],))
    a.regenerate()

    # ---- Step 1: GlobalTension (프리텐션) — run_abaqus_new.py 와 동일 ----
    my_model.StaticStep(
        name='Step-GlobalTension',
        previous='Initial',
        nlgeom=ON,
        stabilizationMagnitude=0.0002,      # Galhofo Reference
        stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
        initialInc=0.0001, minInc=1e-8, maxNumInc=1000
    )

    # ---- Step 2: Buckle — Trigger 스텝을 대체한다 ----
    # 대조 스크립트(run_abaqus_cable.py)에서 완주한 설정을 그대로 사용:
    #   numEigen=100 (음수 모드 건너뛰기), SUBSPACE, vectors=250, maxIterations=5000
    if 'Step-Buckle' in my_model.steps:
        del my_model.steps['Step-Buckle']
    # ---- 필드출력: HF(run_abaqus_new.py:724) 와 동일. EVOL 이 없으면 base_state_probe 가
    #      면적가중을 못 하고 균등가중으로 떨어진다(2026-09-22 실측으로 발견).
    my_model.FieldOutputRequest(name='F-Output-1',
                                createStepName='Step-GlobalTension',
                                variables=('S', 'E', 'U', 'COORD', 'EVOL', 'RF'))
    # 'RF' = 앵커 반력. base_state_probe 가 "각 앵커(정점/클램프)가 당김을 얼마나
    # 흡수하는가"(하중 경로 분담)를 읽는 데 쓴다 — 2026-09-22 첫 좌굴 런에서
    # RF 미출력으로 그 계측이 비어 있었다. 출력 요청은 해석 결과를 바꾸지 않는다.

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
        createStepName='Initial',
        region=inst_memb.sets['All'],
        u3=SET
    )

    # Top 정점: 고정 (케이블 버전의 BC_Anchor(end_c1: u1,u2,u3=0) 대응)
    my_model.DisplacementBC(
        name='BC_Anchor_Top',
        createStepName='Initial',
        region=a.sets['RP_Top_Set'],
        u1=0, u2=0, u3=0
    )

    # Right/Left 정점: GlobalTension 에서 α*DISP_GLOBAL 만큼 당긴다.
    #   케이블 버전의 Disp_Control_Right/Left 와 동일한 방향/크기.
    disp_a = disp
    my_model.DisplacementBC(name='Disp_Control_Right', createStepName='Initial',
                            region=a.sets['RP_Right_Set'], u1=0, u2=0)
    my_model.DisplacementBC(name='Disp_Control_Left', createStepName='Initial',
                            region=a.sets['RP_Left_Set'], u1=0, u2=0)
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

    # 클램프 RP: 고정 (케이블 버전은 클램프 케이블 끝단을 u3=0 으로 잡았고 u1,u2 는
    #   ClampTension 에서 구동했다. 좌굴 모델은 논문 (b) 레시피대로 클램프를 구동하지
    #   않으므로 클램프 위치를 그대로 고정한다.)
    my_model.DisplacementBC(name='BC_Clamp_CL', createStepName='Initial',
                            region=a.sets['RP_CL_Set'], u1=0, u2=0, u3=0)
    my_model.DisplacementBC(name='BC_Clamp_CR', createStepName='Initial',
                            region=a.sets['RP_CR_Set'], u1=0, u2=0, u3=0)

    # ---- Buckle 스텝: 논문 (c) "z-displacement is fixed in the three edges" ----
    all_edges = inst_memb.edges
    a.Set(name='All_Edges', edges=all_edges)
    my_model.boundaryConditions['BC_Stabilize_Z'].deactivate('Step-Buckle')
    my_model.DisplacementBC(
        name='BC_Edges_Only_Z',
        createStepName='Step-Buckle',
        region=a.sets['All_Edges'],
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

    return model_name

# ============================================================================
# 실행 — 한 케이스: 모델 빌드 -> 잡 제출 -> 결과 덤프 -> base state 측정(선택)
# ============================================================================
JOB_NAME = 'Buckle_xc%03d_d%03dum' % (int(round(x_c * 100.0)), int(round(DISP * 1.0e6)))

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
