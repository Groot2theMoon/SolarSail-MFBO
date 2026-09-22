"""
선형 좌굴(고유치) 해석 전용 스크립트 — 케이블만 제외하면 run_abaqus_new.py 와 동일한 모델.

목적
    포스트버클링이 ClampTension step time 0.466 에서 멈추는 원인 = **준축퇴 다중 분기**.
    (실측: A안에서 208 negative eigenvalues, 논문 첫 4개 고유값 차이가 0.04%)
    분기 경로를 '지정'해 주려면 먼저 그 경로(고유모드)를 얻어야 한다.
    이 스크립트가 그 모드를 뽑는다.

논문(Galhofo 2022) §3.2 좌굴 레시피 중 우리에게 필요한 부분
    (a) 케이블 없는 별도 모델          <- 이 스크립트 (케이블 제외)
    (b) 초기 작은 프리텐션(3정점 변위)  <- Step-GlobalTension
    (c) 3개 변의 z 고정                <- BC_Edges_Only_Z
    (d) buckle 스텝 + SUBSPACE 솔버     <- Step-Buckle
    클램프는 **유지**한다 (우리 모델의 주 하중원이며 MFBO 설계변수 x_c/d_c 의 통로다).

run_abaqus_new.py 와의 차이 (그 외는 동일해야 한다)
    제외 : 케이블 part 생성(create_cable_part), Cable 재질/TrussSection,
           connect_cable, Tie/Join 부착, 케이블 끝단 BC(end_c1..end_cr)
    대체 : 케이블이 'RP <-> 구동 끝단'을 잇는 하중 전달 로드였으므로,
           케이블을 빼면 **RP 를 직접 구동/구속**하면 하중 경로가 동일해진다.
             R  : RP_Right_Set 에 (DISP*cos, -DISP*sin)
             L  : RP_Left_Set  에 (-DISP*cos, -DISP*sin)
             Top: RP_Top_Set   을 고정 (케이블 버전의 BC_Anchor 대응)
             CL/CR: RP_CL/CR_Set 을 고정 (클램프는 구동하지 않는다 — 논문 (b) 레시피)
    교체 : Step-Trigger -> Step-Buckle (BuckleStep, SUBSPACE, numEigen=100, vectors=250)
    유지 : 메쉬(seedPart/setMeshControls/setElementType: S4+S3 — HF 와 같은 요소),
           create_rigid_patch(radius=0.2), Coupling KINEMATIC, SIGMA0,
           Step-GlobalTension 의 증분/안정화 설정, NUMCPUS.
    분리 : RUN_DIR_NAME = "buckle" (HF/LF 는 "aba"). 산출물을 섞지 않기 위함이다.
           -> 이 때문에 HF 의 *IMPERFECTION 경로에 ..\buckle\ 를 붙여야 한다 (아래 참조).

α 스윕 (프리텐션 수준)
    λ>0 인 base state 를 찾기 위해 프리텐션을 α 배로 바꿔가며 한 번에 전부 돌린다.
    α = ALPHA_LIST (0.05, 0.10, 0.25, 0.50, 1.00) x DISP_GLOBAL(5e-5 m)
    판정: λ1..λ4 > 0 이며 CONVERGED 인 **최대 α** 를 고르고, 그 모드를 HF 초기결함으로 쓴다.

LF/HF 구분은 없다 (의도적)
    이 스크립트는 순수 선형 좌굴(고유치) 해석만 한다. 스텝은 GlobalTension -> Buckle 둘뿐이다.
    LF 지표/Postbuckle/추출(eval_abaqus.py)/감쇠/trigger 는 이 스크립트에 없다 —
    그것들은 run_abaqus_new.py 의 책무이고, 여기서 중복되면 역할이 흐려진다.

Usage
    abaqus cae noGUI=run_abaqus_buckle.py -- x_c
    abaqus cae noGUI=run_abaqus_buckle.py -- HF 0.5 1.0
        (두 번째 형식은 기존 호출 습관 호환용. fidelity/d_c 는 쓰이지 않으며 그 사실을 출력한다.)

산출물
    code/buckle/Buckle_a<alpha>.odb / .dat / .msg / .sta / .fil / .diag.txt
    .fil 은 *BUCKLE 스텝이 자동 기록한다 (run_abaqus_cable.py 로 검증된 사실).
    -> HF 쪽에서 읽을 때는 **경로를 붙여야 한다**. HF 잡의 작업 디렉터리는 code/aba 이고
       좌굴 산출물은 code/buckle 이므로, 그냥 FILE=Buckle_a025 라고 쓰면 못 찾거나
       (더 나쁘게) code/aba 에 남은 stale .fil 을 조용히 읽는다.
         *IMPERFECTION, FILE=..\buckle\Buckle_a025, STEP=2   (STEP=2 = Buckle 스텝)
       ..\buckle\ 로 시작하는 파일명을 피하려면 절대경로를 써도 된다.
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
print("DEBUG: All sys.argv: " + str(sys.argv))

# ================= CLI =================
#  ★ 2026-09-22 실측: "아무것도 안 돌고 exit code 0" 이 나왔다. 원인은 두 취약점의 겹침이다.
#   (1) 위치 기반 인덱싱(sys.argv[1:])은 Abaqus 가 넘기는 argv 형태에 따라 빈 리스트가 된다.
#       검증된 기존 스크립트들(run_abaqus_new.py 등)이 sys.argv[-3:] (접미 기반)을 쓴 이유가 이것이다.
#   (2) 그 실패를 sys.exit(1) 로 알리려 했더니 CAE 의 noGUI 러너에서 종료코드가 0으로 보였다.
#       => "실패했는데 성공처럼 보인다" 는 최악의 형태. 그래서 이 스크립트는 실패를
#          반드시 **예외**로 올린다 (Abaqus 가 'cae exited with an error' 로 표면화한다).
print("DEBUG: All sys.argv: " + str(sys.argv))


def parse_x_c(argv):
    """argv 에서 x_c 하나만 뽑는다. Abaqus 가 어떤 형태로 넘겨도 동작해야 한다.

    허용 형태:
        -- 0.5
        -- HF 0.5 1.0
        HF 0.5 1.0
        run_abaqus_buckle.py -- 0.5
        C:\...\run_abaqus_buckle.py HF 0.5 1.0
    반환: (x_c 문자열 또는 None, 오류 메시지 또는 None)
    """
    toks = [a for a in argv if a != '--']
    # 스크립트 파일명/경로 토큰 제거
    toks = [a for a in toks
            if not a.lower().endswith('.py') and not re.match(r'^[A-Za-z]:', a)]
    if not toks:
        return None, '인자가 하나도 없습니다'
    if toks[0].upper() in ('LF', 'HF'):
        if len(toks) < 2:
            return None, 'fidelity 뒤에 x_c 가 없습니다: %s' % toks
        if len(toks) >= 3:
            print("%s 참고: fidelity='%s' 및 d_c='%s' 는 쓰이지 않습니다 "
                  "(순수 선형 좌굴 전용)." % (TAG, toks[0], toks[2]))
        return toks[1], None
    if len(toks) >= 3:
        print("%s 참고: 첫 토큰 '%s' 를 fidelity 로 보지 않았습니다. 마지막 값 %s 을 x_c 로 씁니다."
              % (TAG, toks[0], toks[-1]))
    return toks[-1], None


_x_c_raw, _cli_err = parse_x_c(sys.argv)
if _cli_err:
    print("Error: %s" % _cli_err)
    print("Usage: abaqus cae noGUI=run_abaqus_buckle.py -- x_c        (예: -- 0.5)")
    print("       abaqus cae noGUI=run_abaqus_buckle.py -- HF x_c d_c (호환 형식)")
    raise RuntimeError('CLI 인자를 해석할 수 없습니다: %s (위 Usage 참조). '
                       'sys.exit 를 쓰지 않는 이유: CAE noGUI 러너에서 종료코드가 0으로 보인다.'
                       % _cli_err)
try:
    x_c = float(_x_c_raw)
except Exception:
    raise RuntimeError('x_c 를 float 로 변환할 수 없습니다: %r' % (_x_c_raw,))
print("%s x_c = %g" % (TAG, x_c))

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
    # MODEL_NAME 전역 상수를 쓰지 않는다: 이 스크립트는 alpha 마다 모델 이름이 달라진다.
    if model_name is None:
        raise RuntimeError('run_job_safely: model_name 을 반드시 인자로 넘겨야 합니다 '
                           '(이 스크립트의 모델 이름은 alpha 마다 동적으로 생성된다).')
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

    # 잡 결과 핵심 줄을 콘솔에 직접 찍는다 (로그 일부만 붙여넣어도 원인 판별 가능)
    print_job_diag(job_name)
    dump_job_diag(job_name)

    # ABORTED가 아니면서, ODB 파일이 실제로 존재하면 성공으로 간주
    if job.status == ABORTED or not os.path.exists(odb_file):
        raise RuntimeError('Job %s 실패 (Status=%s). sys.exit 대신 예외로 올린다: '
                           'CAE noGUI 러너에서 sys.exit 은 종료코드 0으로 보인다.'
                           % (job_name, str(job.status)))

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
    파일: <RUN_DIR>/<job>.diag.txt
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
            'CONSTANT DAMPING', 'OVERCONSTRAINT', 'INACTIVE DOF',
            'ALLSDTOL', 'SEVERE ELEMENT',
            'DIVERG', 'MINIMUM SPECIFIED', 'TIME INCREMENT REQUIRED',
            'LINE SEARCH', 'ANALYSIS SUMMARY', 'CUTBACKS IN AUTOMATIC',
            # 좌굴 전용 키 (이 스크립트의 핵심 산출물)
            'EIGEN', 'CONVERGED', 'BUCKLING FACTOR', 'MODE NO')
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
    for ext in ('msg', 'dat'):
        fn = '%s.%s' % (job_name, ext)
        if not os.path.exists(fn):
            continue
        with open(fn, 'r', errors='replace') as f:
            _tail = f.read().splitlines()[-45:]
        out.append('===== %s (TAIL 45 - 진짜 사망 원인/요약은 여기에만 있다) =====' % fn)
        out.extend(_tail)
    if len(out) <= 2:
        out.append('(no .sta/.msg/.dat found)')
    dst = os.path.join(_RUN, '%s.diag.txt' % job_name)
    # ASCII 내용 + 명시적 UTF-8 저장: 한글 라벨을 쓰면 Windows 기본 인코딩(CP949)으로
    # 저장되어 다른 도구에서 읽히지 않는다 (2026-09-21 실제 발생).
    import io as _io
    with _io.open(dst, 'w', encoding='utf-8', errors='replace') as f:
        f.write('\n'.join(out) + '\n')
    print('%s 진단 요약 저장: %s (%d줄)' % (TAG, dst, len(out)))
    return dst


def print_job_diag(job_name):
    """잡 산출물(.msg/.dat)의 원인 판별용 핵심 줄만 콘솔에 찍는다."""
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


def print_job_eigen(job_name):
    """좌굴 고유값을 .dat/.msg 에서 '원문 그대로' 출력한다 — 파서 없음.

    과거 A안에서 .dat 헤더 형식을 추정한 파서가 '고유치를 찾지 못했습니다' 로 실패했다.
    여기서는 EIGEN 류 키워드가 나온 줄 주변을 통째로 덤프해서 사람이 직접 판독한다.
    """
    KEYS = ('EIGEN', 'CONVERGED', 'BUCKLING FACTOR', 'MODE NO')
    for ext in ('dat', 'msg'):
        fn = '%s.%s' % (job_name, ext)
        if not os.path.exists(fn):
            print('%s [EIGEN:%s] %s 없음' % (TAG, job_name, fn))
            continue
        with open(fn, 'r', errors='replace') as f:
            lines = f.read().splitlines()
        hit = [i for i, ln in enumerate(lines)
               if any(k in ln.upper() for k in KEYS)]
        if not hit:
            print('%s [EIGEN:%s] %s 에 고유값 키워드가 없습니다 (총 %d줄). '
                  '이 경우 .msg 의 ***ERROR 를 보세요.' % (TAG, job_name, fn, len(lines)))
            continue
        keep = set()
        for i in hit:
            for j in range(max(0, i - 2), min(len(lines), i + 25)):
                keep.add(j)
        print('%s [EIGEN:%s] %s 원문 덤프 — 히트 %d개 / 출력 %d줄'
              % (TAG, job_name, fn, len(hit), len(keep)))
        prev = None
        for j in sorted(keep):
            if prev is not None and j != prev + 1:
                print('      ...')
            print('      | %s' % lines[j][:160])
            prev = j


MODEL_PREFIX = 'SailModel_Buckle'
INSTANCE_NAME = 'MEMBRANE-1'

BASE = 20.0   # m
HEIGHT = 10.0 # m
THICKNESS = 5.0e-6

# 좌표 정의 (run_abaqus_new.py 와 동일)
V1 = (BASE/2.0, HEIGHT, 0.0) # Top
V2 = (BASE, 0.0, 0.0)        # Right
V3 = (0.0, 0.0, 0.0)         # Left

def clamp_coord_L(x): return (10-10*x, 10-10*x, 0)
def clamp_coord_R(x): return (10+10*x, 10-10*x, 0)

V_CL = clamp_coord_L(x_c)
V_CR = clamp_coord_R(x_c)

# ---- 사전 장력(prestrain) 캘리브레이션 (run_abaqus_new.py 와 동일) ----
PRETENSION_SCALE = 10.0          # 프리텐션 변위 = 5e-6 m * 이 값 = 5e-5 m
DISP_GLOBAL = 0.000005 * PRETENSION_SCALE    # 운용점: 코너 당김 5e-5 m

# ---- 좌굴 전용 상수 ----
# α 스윕: 프리텐션 수준을 바꿔가며 λ>0 인 base state 를 찾는다.
ALPHA_LIST = (0.05, 0.10, 0.25, 0.50, 1.00)
# 좌굴 스텝의 perturbation 크기.
#   근거(A안 주석 + run_abaqus_cable.py 실측): 좌굴 스텝의 nonzero prescribed BC 는
#   '증분 응력' 을 만들고 그 증분이 미분 초기응력 강성 K_delta 를 만든다.
#   K_delta 가 K0 대비 너무 작으면 고유값 분리가 나빠져 subspace 가
#   'EIGENVALUES CANNOT BE FOUND' 로 실패한다. 성공한 대조 스크립트는 0.01 m 를 썼다.
# ---- 좌굴 솔버 (논문 :247 "The solver subspace interaction is selected" => 기본 SUBSPACE) ----
#   A/B 시험은 이 상수 한 줄만 바꾼다 (환경변수 금지 원칙).
#   두 솔버는 인자 이름이 다르다:
#     SUBSPACE -> vectors        (검증된 대조 run_abaqus_cable.py 와 동일 조합)
#     LANCZOS  -> blockSize / minEigen / maxEigen
BUCKLE_SOLVER = 'SUBSPACE'      # 'SUBSPACE' | 'LANCZOS'
BUCKLE_BLOCK_SIZE = 8           # LANCZOS 전용
BUCKLE_MIN_EIGEN = 0.0          # LANCZOS 전용 (음의 고유값도 보고 싶으면 -1e30 등)
BUCKLE_MAX_EIGEN = None         # LANCZOS 전용 (None 이면 인자를 아예 넘기지 않는다)

PERTURBATION = 0.01     # m
N_EIG_BUCKLE = 100      # 추출 요청 고유값 수 (음수 모드 건너뛰기 위해 100 — 대조 스크립트와 동일)
BUCKLE_VECTORS = 250    # subspace 기저 벡터 수 (numEigen 의 2.5배 — 대조 스크립트와 동일)
SIGMA0 = 500.0          # 수렴 보조용 초기응력 [Pa] (run_abaqus_new.py 와 동일)

print("%s x_c=%.3g -> V_CL=%s V_CR=%s" % (TAG, x_c, V_CL, V_CR))
print("%s DISP_GLOBAL=%.3e m  PERTURBATION=%.3e m  ALPHA_LIST=%s"
      % (TAG, DISP_GLOBAL, PERTURBATION, str(ALPHA_LIST)))
print("%s N_EIG_BUCKLE=%d BUCKLE_VECTORS=%d SUBSPACE" % (TAG, N_EIG_BUCKLE, BUCKLE_VECTORS))

# 하중 각도 (28.6도) — 케이블 방향과 동일하게 유지 (하중 경로 동일화)
angle_deg = 28.6
angle_rad = np.deg2rad(angle_deg)
cos_val = float(np.cos(angle_rad))
sin_val = float(np.sin(angle_rad))


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


def alpha_tag(alpha):
    """0.05 -> 'a005'  (모델/잡 이름용, 파일명에 점이 들어가지 않게)"""
    return 'a%03d' % int(round(alpha * 100.0))


def build_model(alpha):
    """run_abaqus_new.py 의 모델 생성을 '케이블만 제외' 하고 재현한다."""
    global my_model
    model_name = '%s_%s' % (MODEL_PREFIX, alpha_tag(alpha))
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
                                variables=('S', 'E', 'U', 'COORD', 'EVOL'))

    # ---- 좌굴 스텝: 솔버는 코드 상수 하나로 교체 (SUBSPACE <-> LANCZOS) ----
    #   논문 :247 은 SUBSPACE 를 명시 선택했다. 바꾸면 논문 사양에서 이탈한다.
    #   두 솔버는 인자 이름이 다르므로 분기해서 넘긴다 (LANCZOS 는 vectors 를 받지 않는다).
    _eig = dict(name='Step-Buckle', previous='Step-GlobalTension',
                numEigen=N_EIG_BUCKLE, eigensolver=BUCKLE_SOLVER)
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
    disp_a = alpha * DISP_GLOBAL
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


# ================= α 스윕 실행 =================
SUMMARY = []
for _alpha in ALPHA_LIST:
    _tag = alpha_tag(_alpha)
    _model = None
    _job = 'Buckle_' + _tag
    print("")
    print("=" * 78)
    print("%s ===== alpha=%.4g  (프리텐션 %.4e m)  model=%s_%s  job=%s ====="
          % (TAG, _alpha, _alpha * DISP_GLOBAL, MODEL_PREFIX, _tag, _job))
    print("=" * 78)
    try:
        _model = build_model(_alpha)
    except Exception as _e:
        print("!!! ERROR: alpha=%.4g 모델 생성 실패: %s" % (_alpha, _e))
        print("----- traceback (실패한 줄을 확인하세요) -----")
        print(traceback.format_exc())
        print("---------------------------------------------")
        SUMMARY.append((_alpha, _job, 'BUILD_FAIL', str(_e)))
        continue

    # 스윕 중 한 케이스가 죽어도 나머지를 계속 돌린다 (run_job_safely 는 실패 시
    # RuntimeError 를 올린다 — sys.exit 은 CAE 러너에서 종료코드 0으로 보이므로 쓰지 않는다).
    _ok = True
    try:
        run_job_safely(_job, _model)
    except Exception as _e:
        _ok = False
        print("!!! alpha=%.4g : %s 실패 -> 다음 alpha 로 계속: %s" % (_alpha, _job, _e))

    try:
        print_job_eigen(_job)
    except Exception as _e:
        print("!!! [EIGEN] %s 고유값 덤프 실패: %s" % (_job, _e))

    # base state 의 압축 정도 측정 (λ 판정을 물리와 함께 보기 위해)
    try:
        _probe = os.path.join(_HERE, 'base_state_probe.py')
        _odb = '%s.odb' % _job
        if os.path.exists(_probe) and os.path.exists(_odb):
            print("%s base state 측정: %s / Step-GlobalTension" % (TAG, _odb))
            subprocess.call('abaqus python "%s" "%s" Step-GlobalTension' % (_probe, _odb),
                            shell=True)
        else:
            print("%s base state 측정 건너뜀 (probe=%s, odb=%s)"
                  % (TAG, os.path.exists(_probe), os.path.exists(_odb)))
    except Exception as _e:
        print("%s base state 측정 실패(무시): %s" % (TAG, _e))

    SUMMARY.append((_alpha, _job, 'OK' if _ok else 'JOB_FAIL',
                    job_completed_ok(_job)))

    # 스윕은 모델을 5개 만든다. 다 쓰면 지운다 (스윕 도중 메모리 부족으로 남은 alpha 가
    # 죽는 것을 막는다. odb/.dat/.fil 은 이미 디스크에 있으므로 잃는 것이 없다).
    _freed = False
    if _model is not None and _model in mdb.models:
        del mdb.models[_model]
        _freed = True
    print("%s alpha=%.4g 종료 (model 해제=%s, 남은 모델 %d개)"
          % (TAG, _alpha, _freed, len(mdb.models.keys())))

print("")
print("=" * 78)
print("%s ===== α 스윕 요약 =====" % TAG)
print("     alpha        프리텐션[m]     job            완주판정")
for _row in SUMMARY:
    _a = _row[0]
    print("     %-12.4g %-15.4e %-14s %s"
          % (_a, _a * DISP_GLOBAL, _row[1], _row[-1]))
if not any(r[2] == 'OK' for r in SUMMARY):
    raise RuntimeError('전 alpha 에서 잡이 하나도 완주하지 못했습니다 (SUMMARY=%s). '
                       '위 로그의 첫 ERROR/traceback 을 보세요. '
                       '이 상태를 exit 0 으로 끝내면 성공처럼 보이므로 예외로 올린다.' % str(SUMMARY))

print("%s 선택 규칙: 위 .dat 원문 덤프에서 λ1..λ4 > 0 이며 CONVERGED 인 최대 alpha." % TAG)
print("%s 그 alpha 의 모드를 HF 초기결함으로 쓴다: *IMPERFECTION, FILE=%s, STEP=2"
      % (TAG, r'..\buckle\Buckle_a<alpha>'))
print("=" * 78)
