#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""클램프 없는 모델에서 **선형 좌굴해석(*BUCKLE)** 으로 좌굴모드를 구하는 전용 스크립트 (2026-09-22)

무엇을 하는가 (한 줄)
    클램프가 '없는' 형상에서 좌굴모드를 **선형 좌굴해석**으로 계산해 ODB 에 남기고,
    HF(run_abaqus.py)가 먹을 모드표(modes_ClampFree_Buckle.txt)까지 이 스크립트 안에서 만든다.

왜 클램프 없는 모델인가 (2-모델 레시피, Galhofo2022)
    run_abaqus.py -- HF 는 클램프(당김)가 있는 모델이라 base state 가 이미 분기점을 넘어
    선형 좌굴 스펙트럼을 주지 못한다(실측: 음수 고유값 598~2897개 / CONVERGED=0 /
    "THE EIGENVALUES CANNOT BE FOUND"). 그래서 논문의 좌굴 모델처럼 '클램프 없이 당긴'
    상태에서 모드를 뽑아, 클램프가 있는 HF 의 초기 결함으로 주입한다.

왜 ODB 인가 (실측 2026-09-22)
    ClampFree_Buckle.dat:7025 -> "***WARNING: FILE OUTPUT IS NOT AVAILABLE FOR BUCKLING ANALYSIS"
    *BUCKLE 스텝은 .fil 출력이 금지된다 -> 좌굴모드의 유일한 통로는 **ODB 의 모드 프레임**이다.
    (진동 *FREQUENCY 스텝은 .fil 이 허용되지만 좌굴모드가 아니라 다른 물리량이므로 쓰지 않는다.)

스텝 구성
    Step 1  Step-GlobalTension : 코너 당김(=base state) + 평탄 강제(u3=SET)
    Step 2  Step-Buckle       : 선형 좌굴해석(*BUCKLE, SUBSPACE). 평탄 강제 해제 + 모서리만 u3 고정,
                                코너 당김을 PERTURBATION(0.01 m)로 키워 '하중 패턴'(λ 기준)으로 준다.

요청 수 규칙 (run_abaqus_cable.py 좌굴 스텝 주석에서 확인)
    '요청 고유값 수(N_EIG_BUCKLE) > base state 의 음수 고유값 수' 여야 양수 좌굴모드가 subspace 창에 들어온다.
    성공 실적: 케이블 런 = 음수 52 < 100 -> CONVERGED=4.

메쉬/형상/BC/프리텐션 정의 줄은 run_abaqus.py(HF)·run_abaqus_buckle.py 와 '문자 그대로' 같아야 한다
(check_model_consistency.py 가 제출 전에 정적으로 검사). 모드는 HF 와 같은 노드 라벨에 정의되어
기하 섭동(IMPERFECTION_MODE='odb_table')으로 이식된다.

실행
    abaqus cae noGUI=run_abaqus_mode.py          # 잡 1회 -> ODB + modes_ClampFree_Buckle.txt
    (이어서) abaqus cae noGUI=run_abaqus.py -- HF <x_c> <d_c>

종료 코드
    0 = 좌굴모드 1개 이상 + 모드표 저장 / 1 = 0개(실패) 또는 잡 이상
"""

# (2026-09-24 제거) from matplotlib.image import LANCZOS — IDE 자동삽입. matplotlib 의 LANCZOS 는
#   이미지 필터이고 Abaqus Python 에 matplotlib 이 없으면 ImportError 로 죽는다(§16).
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

# 로깅 — cae noGUI 는 스크립트 stdout 을 콘솔로 보내지 않고 execfile 이라 __file__ 도 없다(§1).
#   파일이 주 채널이다. canary 로 이 파일이 실제 실행되는지부터 확인한다.
try:
    with open('mode_run_log.txt', 'w') as _cf:
        _cf.write('[canary] run_abaqus_mode.py start %s  cwd=%s  argv=%s\n'
                  % (time.strftime('%Y-%m-%d %H:%M:%S'), os.getcwd(), list(sys.argv)))
except Exception:
    pass

try:
    sys.stdout.reconfigure(line_buffering=True)      # 되면 좋고, 안 되면 파일 채널이 담당한다
except Exception:
    pass

_EMIT_PATHS = None


def _emit_paths():
    """로그 파일 후보 경로(첫 번째가 주 경로). __file__ 이 없는 execfile 환경을 전제로 한다."""
    cand = []
    here = globals().get('_HERE')                    # _resolve_here() 가 고른 스크립트 폴더
    if here:
        cand.append(here)
    try:
        if sys.argv and sys.argv[0]:
            cand.append(os.path.dirname(os.path.abspath(sys.argv[0])))
    except Exception:
        pass
    cand.append(os.getcwd())                         # Abaqus 실행 CWD ( = code\ )
    out, seen = [], set()
    for d in cand:
        if not d:
            continue
        p = os.path.join(d, 'mode_run_log.txt')
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def emit(*a):
    """콘솔 + mode_run_log.txt(후보 경로 전부) 기록. **어떤 경우에도 예외를 올리지 않는다**."""
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

emit("DEBUG: All sys.argv: " + str(sys.argv))

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
from aba_imperfection import (count_modes, parse_eigenvalues, verify_inp_element_types,  # noqa: E402
                              job_completed_from_logs)
# aba_grid_mesh: BC 노드셋 헬퍼만 쓴다(boundary_node_labels/make_set). 격자는 미사용(§10).
import aba_grid_mesh                                                            # noqa: E402
emit("[run_abaqus_mode] _HERE = %s" % _HERE)
# 코드 지문 — 스크립트·공용 모듈 md5 를 로그 첫머리에 찍어 'pull 누락'과 '옛 코드'를 구분한다(§1).
try:
    import hashlib as _hl
    for _fn in ('run_abaqus_mode.py', 'aba_imperfection.py'):
        _p = os.path.join(_HERE, _fn)
        with open(_p, 'rb') as _fh:
            emit("[CODE] %s md5=%s" % (_fn, _hl.md5(_fh.read()).hexdigest()[:12]))
except Exception as _ce:
    emit("[CODE] 지문 계산 실패(무시): %s" % _ce)


def run_job_safely(job_name, model_name=None):
    """잡 제출 + 완료 대기. 제출 전 이전 산출물을 지워 stale .fil/.odb 오판을 막는다."""
    model_name = model_name or MODEL_NAME
    for _ext in ('odb', 'fil', 'sta', 'msg', 'lck', 'com', 'prt', 'sim', 'log',
                 'dat', 'res', 'abq', 'ipm', 'mdl', 'stt', 'cid'):
        _f = '%s.%s' % (job_name, _ext)
        if os.path.exists(_f):
            emit("Removing stale artifact: %s" % _f)
            try:
                os.remove(_f)
            except OSError as _e:
                emit("  (warning) could not remove %s: %s" % (_f, _e))
    if job_name in mdb.jobs:
        del mdb.jobs[job_name]
    job = mdb.Job(name=job_name, model=model_name)
    emit("Submitting Job: %s" % job_name)
    job.writeInput(consistencyChecking=OFF)
    # 제출 전 요소 타입 검증 — 조용한 폴백을 Standard 토큰 전에 잡는다(§3).
    _ok_et, _cnt_et, _err_et = verify_inp_element_types(
        job_name + '.inp', (ELEM_CODE_QUAD, ELEM_CODE_TRI))
    emit("[GUARD] .inp 요소 타입 블록: %s" % (_cnt_et or _err_et))
    if _err_et:
        emit("!!! ERROR: .inp 요소 타입 검증 실패: %s" % _err_et)
        return False
    if not _ok_et:
        emit("!!! ERROR: 요소코드 폴백 의심 (요청=%s / .inp=%s / %s)"
             % ((ELEM_CODE_QUAD, ELEM_CODE_TRI), _cnt_et, _err_et))
        emit("    잡을 제출하지 않는다(라이선스 절약). ELEM_CODE_* 지정 방식을 고쳐라.")
        return False
    job.submit(consistencyChecking=OFF)
    job.waitForCompletion()
    time.sleep(1.0)
    # 잡 성공 판정 = 산출물 + 완주 문자열(§2).
    if (job.status == ABORTED or not os.path.exists(job_name + '.odb')
            or not job_completed_from_logs(job_name)):
        emit("!!! ERROR: Job %s failed (status=%s, odb=%s, .msg 완주=%s)"
             % (job_name, str(job.status), os.path.exists(job_name + '.odb'),
                job_completed_from_logs(job_name)))
        return False
    emit("Job %s finished (Status: %s)." % (job_name, str(job.status)))
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
            emit("[DIAG:%s] %s 없음" % (job_name, fn))
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
        emit("[DIAG:%s] %s (%.0f KB) 핵심줄 %d개" % (job_name, fn, size / 1024.0, len(hits)))
        for ln in hits[-10:]:
            emit("      | %s" % ln[:150])

    # 좌굴모드의 실제 산출물은 ODB 다(*BUCKLE 은 .fil 출력이 금지된다 — 실측 7025).
    #   ODB 프레임/모드표 생성은 잡 후 8단계에서 한 번만 한다(여기서는 파일 존재만 본다).
    fn = '%s.odb' % job_name
    emit("[DIAG:%s] ODB %s (%.0f KB)"
          % (job_name, 'OK' if os.path.exists(fn) else '없음',
             (os.path.getsize(fn) / 1024.0) if os.path.exists(fn) else 0.0))

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
            emit("[DIAG:%s] 'SYSTEM MATRIX HAS ... NEGATIVE EIGENVALUES' 를 못 찾음" % job_name)
        else:
            emit("[DIAG:%s] base state 음수 고유값 %d개 vs 요청 %d개" % (job_name, neg, N_EIG_BUCKLE))
            emit("[DIAG:%s]   -> %s" % (job_name,
                  '요청 수가 충분하다(창 안에 양수 모드가 들어온다)'
                  if N_EIG_BUCKLE > neg else
                  '요청 수 부족 (N_EIG_BUCKLE <= 음수 개수) — 래더 L4'))
            if neg > 0:
                emit("[DIAG:%s]   !! 음수 고유값 %d개 = base state 자체가 이미 좌굴/압축 상태다."
                      % (job_name, neg))
                emit("[DIAG:%s]      (Abaqus 오류문의 'INSTABILITIES IN THE BASE STATE' 와 같은 진단)"
                      % job_name)
                emit("[DIAG:%s]      창을 넓히는 것만으로는 안 풀린다 — 4회차 실측: 88개 / 양수 모드 0개."
                      % job_name)
                emit("[DIAG:%s]      위 [R-13] 의 '면내 압축(<0) 면적비' 를 먼저 확인해라(프리텐션 부족 여부)."
                      % job_name)


# =====================================================================
# 모델 정의 — HF 와 '문자 그대로' 같은 줄 (check_model_consistency.py 가 검사)
# =====================================================================
MODEL_NAME = 'SailModel_Triangle'
INSTANCE_NAME = 'MEMBRANE-1'

BASE = 20.0   # m
HEIGHT = 10.0 # m
# 막 요소는 1차(S4/S3) 유지 — 2차는 기각됐다(§4). 요소코드는 HF 와 같아야 한다
#   (*IMPERFECTION 이 노드 라벨로 주입된다).
ELEM_CODE_QUAD = S4         # [2026-10-04] 막(M3D4/M3D3) 실험 철회 -> 1차 셸 복귀.
ELEM_CODE_TRI = S3          # [2026-10-04] 위와 같은 이유로 3절점 1차 셸 복귀
SEED_DIV = 200.0           # seed = BASE/SEED_DIV -> 약 1.82만 요소 (실측 2026-09-28)
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

# 프리텐션 정의 — HF 와 다르다(모드 소스는 3꼭짓점 prescribed displacement). 근거는 §8.
#   메쉬/형상/요소/재료는 동일하므로 모드 노드 라벨 매핑은 그대로 성립한다.
# CHECKER-DIVERGENCE: DISP_GLOBAL, SIGMA0, PRETENSION_SCALE, stabilizationMagnitude, radius
#   이 선언이 남아 있는 동안에는 PATCH_RADIUS 를 바꿔도 검사기가 통과시킨다 —
#   PATCH_RADIUS 를 바꿀 때는 사람이 HF 값(0.2)과 직접 대조할 것.
PRETENSION_MODE = 'corner2'  # [2026-09-24 오라클 재현] 케이블 런(run_abaqus_cable.py)과 동일한 구동:
                             #   정점(BC_Anchor_Top)은 완전 고정한 채, 아래 두 꼭짓점만 케이블 축으로 당긴다.
                             #   성공 실적: 이 설정에서 CONVERGED=4 (양수 λ 4개).
                             # 'paper3' = 논문 좌굴모델: 꼭짓점 3개를 모두 당김(3중 대칭).
                             #   5·7회차 실측: 양수 2개 / 0개 -> 재현 실험 뒤 bisect 단계에서 하나씩 되돌린다.
DISP_GLOBAL = 1.8e-5         # 꼭짓점 당김 [m] — 모드 소스의 운용점. HF 와 **다르다(의도된 선언 분기)**.
# 모드 소스는 저프리텐션에 둔다 — HF 운용점(165e-6)은 base state 가 불안정해 모드 추출이 실패했다(§8).
#   양수 모드가 나온 회차는 전부 저프리텐션(1e-7 / 5e-5 m)이고 1e-3 m 에서는 0개였다.
#   1e-7 m 대는 모드가 영에너지로 퇴화하므로 그보다 위에 둔다. 안정 상한은 40e-6/80e-6 각 1회로 좁힌다.
DISP_TOP_OVER_CORNER = 1.4142135623
SIGMA0 = 700.0                               # 초기 가짜 응력(수렴 보조) [Pa] — 프리텐션의 대체물이 아니다
                                             #   [2026-09-24 오라클 재현] 케이블 런 값 = 700.0 (HF 는 500.0)
MODE_STABILIZATION = 0.0005                  # GlobalTension 안정화 계수
PERTURBATION = 0.01                          # 좌굴 스텝 섭동 크기 [m] (케이블 런·HF 와 동일)
PATTERN_SIGN = 1.0                           # 좌굴 '하중 패턴'의 부호: 1.0 = 바깥으로 더 당김 (오라클과 동일)
                                             # PATTERN_SIGN 은 바꾸지 않는다 — λ 부호를 정하는 것은 패턴이 아니라 base state 크기다(§8).
N_EIG_BUCKLE = 150                           # [2026-10-07] 100 -> 200 (buckle 정합)
#   근거: buckle 은 클램프가 있으면 100/250 에서 ***ERROR: THE EIGENVALUES CANNOT BE
#   FOUND (INSTABILITIES IN THE BASE STATE) 로 죽고, 200/500 에서 CONVERGED=4 로 성공했다.
#   실측: passive 와 none 의 lambda 는 상대차 6e-6(수치 잡음)으로 동일 -> 창만 문제였다.
#   두 스크립트의 창을 통일한다. 요청 수는
                                             #   **base state 의 음수 고유값 개수보다 커야 한다**
                                             #   (규칙: N_EIG_BUCKLE > 음수 개수).
                                             #   실측: 구 메쉬(18,200)에서 음수 48개 -> 요청 8 로는
                                             #   창이 좁아 CONVERGED=0 / THE EIGENVALUES CANNOT BE
                                             #   FOUND 였다. 8 로 줄였던 근거(음수 2개)는 **다른 메쉬**
                                             #   (12,086)에서 얻은 값이라 잘못된 진단이었다.
                                             #   => 메쉬를 먼저 고정하고, [DIAG] 의 음수 개수를 보고
                                             #      그보다 큰 값을 고른다(과거 성공값 100).
BUCKLE_VECTORS = 300                         # [2026-10-05 복귀] 40 -> 250. 기저 벡터는 요청 수와
                                             #   함께 움직인다(과거 성공 조합 100/250).
N_MODES = 4                                  # HF 에 주입할 모드 수(= ODB 모드 프레임에서 뽑는 개수)
PATCH_RADIUS = 0.2                           # 꼭짓점 강체패치 반경 [m] — HF 와 동일값(0.2)으로 복원.
COUPLING_TYPE = DISTRIBUTING     # 패치 절점 결합 방식: 'DISTRIBUTING' | 'KINEMATIC'.
# Coupling — DISTRIBUTING + weightingMethod=UNIFORM. KINEMATIC+WHOLE_SURFACE 는 응력 특이점(52배)을
#   만들었다(§7). couplingType 에 문자열을 넘기면 죽는다 — abaqusConstants 심볼을 쓴다.
#   influenceRadius 는 API 필수 인자다. 모델 지문에 포함되므로 값을 바꾸면 모드 재추출이 필요하다.
MODE_STEP_NAME = 'Step-Buckle'               # 이 모델의 2번째 스텝(HF 의 MODE_SOURCE_STEP=2 와 짝)
JOB_NAME = 'ClampFree_Buckle'
MODE_TABLE = 'modes_%s.txt' % JOB_NAME       # HF 는 ..\modes_ClampFree_Buckle.txt 를 읽는다
WRITE_MODE_TABLE = True                      # 잡 성공 후 ODB 에서 모드표를 직접 만든다
RUN_BASE_STATE_PROBE = True                  # base state 압축영역 측정(토큰 소량, 결과 해석용)
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
my_model.HomogeneousShellSection(name='Section-Membrane', material='Kapton', thickness=THICKNESS)   # [2026-10-04] 막 실험 철회로 셸 섹션 복귀

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


def create_rigid_patch(name, coord, radius=PATCH_RADIUS):
    """강체 패치 = 기존 노드에 Coupling(COUPLING_TYPE). face partition 을 하지 않으므로
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
        influenceRadius=WHOLE_SURFACE, couplingType=COUPLING_TYPE, weightingMethod=UNIFORM,
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
        emit("Warning: Node not found by sphere, trying closest for " + name)
        node_start = inst.nodes.getClosest(coordinates=sail_corner)
    region_start = regionToolset.Region(nodes=node_start)
    cable_len = LEN_TOP if 'Top' in name else LEN_BOT
    end_coord = (sail_corner[0] + target_vec[0]*cable_len, sail_corner[1] + target_vec[1]*cable_len, 0.0)
    node_end = inst.nodes.getByBoundingSphere(center=end_coord, radius=radius)
    region_end = regionToolset.Region(nodes=node_end)
    return region_start, region_end


p.seedPart(size=BASE/SEED_DIV, deviationFactor=0.1)  # 약 1.82만개 (1차 요소)
p.setMeshControls(regions=p.faces, elemShape=QUAD_DOMINATED, technique=FREE, algorithm=MEDIAL_AXIS)
elemTypeQuad = ElemType(elemCode=ELEM_CODE_QUAD, elemLibrary=STANDARD)
elemTypeTri = ElemType(elemCode=ELEM_CODE_TRI, elemLibrary=STANDARD)
p.setElementType(regions=(p.faces,), elemTypes=(elemTypeQuad, elemTypeTri))
p.generateMesh()
# 요소 타입 읽기 검증 — 선언이 아니라 실제로 무엇이 붙었는가를 본다(§3). 실패해도 계속 간다.
for _shp, _nm in ((QUAD, 'QUAD'), (TRI, 'TRI')):
    try:
        _et = p.getElementType(region=regionToolset.Region(faces=p.faces), elemShape=_shp)
        emit("[ET-CHECK] %s -> elemCode=%s" % (_nm, getattr(_et, 'elemCode', _et)))
    except Exception as _e2:
        emit("[ET-CHECK] %s 읽기 생략 (%s) — 실제 검증은 제출 전 .inp 가드가 한다" % (_nm, _e2))
emit("[ET-CHECK] 메쉬 요소 %d개 / 요청 elemCode=(%s, %s)"
     % (len(p.elements), getattr(ELEM_CODE_QUAD, 'name', ELEM_CODE_QUAD),
        getattr(ELEM_CODE_TRI, 'name', ELEM_CODE_TRI)))
a.regenerate()

# 2. 꼭짓점 강체 패치 3개 (클램프 패치는 만들지 않는다 — 이 모델의 요점)
rp1_obj, rp1_reg = create_rigid_patch('Top', V1)
rp2_obj, rp2_reg = create_rigid_patch('Right', V2)
rp3_obj, rp3_reg = create_rigid_patch('Left', V3)

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
    stabilizationMagnitude=MODE_STABILIZATION,
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
# [paper3] 위 꼭짓점도 바깥(+y = 위 케이블 축)으로 당긴다 -> 3중 대칭 프리텐션 (논문 좌굴모델).
#   BC_Anchor_Top 은 Initial 에서 u1=u2=u3=SET(완전고정) 이므로 여기서 u1/u2 를 값으로 덮어쓴다.
#   u3/회전은 고정 유지 -> 면외 구속 조건은 그대로다.
if PRETENSION_MODE == 'paper3':
    my_model.boundaryConditions['BC_Anchor_Top'].setValuesInStep(
        stepName='Step-GlobalTension', u1=0.0, u2=DISP_GLOBAL * DISP_TOP_OVER_CORNER
    )
    emit("[run_abaqus_mode] 프리텐션 = paper3 (3꼭짓점 당김): 아래 %.4e m / 위 %.4e m"
          % (DISP_GLOBAL, DISP_GLOBAL * DISP_TOP_OVER_CORNER))
else:
    emit("[run_abaqus_mode] 프리텐션 = corner2 (아래 두 꼭짓점만 %.4e m, 비대칭)" % DISP_GLOBAL)

# 초기 가짜 응력 (수렴 보조)
my_model.Stress(
    name='Initial_Stiffness',
    region=inst_memb.sets['All'],
    distributionType=UNIFORM,
    sigma11=SIGMA0, sigma22=SIGMA0, sigma33=0.0,
    sigma12=0.0, sigma13=0.0, sigma23=0.0
)

# Step 2: **선형 좌굴해석(*BUCKLE)** — 이 스크립트의 본체
#   실측(ClampFree_Buckle.dat:7025): *BUCKLE 스텝은 .fil 출력이 금지된다
#     "***WARNING: FILE OUTPUT IS NOT AVAILABLE FOR BUCKLING ANALYSIS"
#   -> 좌굴모드(고유벡터)는 ODB 모드 프레임에만 있다. 잡 후 8단계에서 ODB 를 읽어 모드표를 만든다.
#   요청 수 규칙: N_EIG_BUCKLE > base state 의 음수 고유값 수 (케이블 런 주석에서 확인)
if MODE_STEP_NAME in my_model.steps:
    del my_model.steps[MODE_STEP_NAME]
my_model.BuckleStep(
    name=MODE_STEP_NAME,
    previous='Step-GlobalTension',
    numEigen=N_EIG_BUCKLE,
    eigensolver=SUBSPACE,
    vectors=BUCKLE_VECTORS,
    maxIterations=5000
)
# 좌굴 스텝의 '하중 패턴': 꼭짓점 당김을 PERTURBATION * PATTERN_SIGN 만큼 준다(λ 를 이 패턴 기준으로 얻는다).
#   4회차 스펙트럼이 전부 λ<0 이었지만, 그 원인은 패턴 부호가 아니라 base state 프리텐션 과대였다
#   (성공한 케이블 런도 같은 '바깥 당김' 패턴으로 λ 4개가 양수다 — 위 PATTERN_SIGN 주석 참조).
#   하중 패턴은 오라클·논문과 같은 방향(바깥 당김)으로 유지한다.
#   논문 좌굴 BC 의 '방향비'만 가져온다: 정점 u_y / 아래 두 꼭짓점 |(5e-4,-5e-4)| = 1e-3/7.071e-4 = √2
#   (절대 크기는 논문값이 아니라 DISP_GLOBAL = 1.8e-5 m 를 쓴다 — 위 상수 근거 참조)
#   (패턴이 비대칭이면 좌굴 스펙트럼이 뭉개진다 — 3회차 실측의 교훈)
my_model.boundaryConditions['Disp_Control_Right'].setValuesInStep(
    stepName=MODE_STEP_NAME,
    u1=PATTERN_SIGN * PERTURBATION * cos_val,
    u2=-PATTERN_SIGN * PERTURBATION * sin_val
)
my_model.boundaryConditions['Disp_Control_Left'].setValuesInStep(
    stepName=MODE_STEP_NAME,
    u1=-PATTERN_SIGN * PERTURBATION * cos_val,
    u2=-PATTERN_SIGN * PERTURBATION * sin_val
)
if PRETENSION_MODE == 'paper3':
    my_model.boundaryConditions['BC_Anchor_Top'].setValuesInStep(
        stepName=MODE_STEP_NAME, u1=0.0, u2=PATTERN_SIGN * PERTURBATION * DISP_TOP_OVER_CORNER
    )
emit("[MODE] 좌굴 하중 패턴 부호 PATTERN_SIGN=%+0.1f (%s) — 크기 %g m"
      % (PATTERN_SIGN, '안쪽(당김을 푸는 방향)' if PATTERN_SIGN < 0 else '바깥 당김', PERTURBATION))
my_model.boundaryConditions['BC_Stabilize_Z'].deactivate(MODE_STEP_NAME)
# all_edges = inst_memb.edges -> 좌표 판정으로 교체. orphan mesh 에서는 기하 edge 가 비고,
#   HF/buckle 과 같은 함수를 써야 노드 집합이 일치한다(§10a). mode 는 클램프가 없어 제외가 없다.
_keep_labels = aba_grid_mesh.boundary_node_labels(inst_memb, base=BASE, height=HEIGHT)
aba_grid_mesh.make_set(a, 'All_Edges_NoClamp', inst_memb, _keep_labels)
emit("[BC] All_Edges_NoClamp (mode): 경계 노드 %d개 (클램프 없음, 전체 막 노드 %d개)"
     % (len(_keep_labels), len(inst_memb.nodes)))
my_model.DisplacementBC(
    name='BC_Edges_Only_Z',
    createStepName=MODE_STEP_NAME,
    region=a.sets['All_Edges_NoClamp'],
    u3=SET
)

# 필드출력: base_state_probe 가 면적가중(EVOL)·반력(RF)을 읽는다
my_model.fieldOutputRequests['F-Output-1'].setValues(
    variables=('S', 'E', 'U', 'RF', 'COORD', 'EVOL'),
    frequency=1
)

# 스텝 순서 가드 — HF 는 MODE_SOURCE_STEP 번째 스텝의 프레임을 모드 소스로 읽는다.
#   순서가 밀리면 엉뚱한 스텝을 조용히 읽으므로 여기서 중단한다.
_step_seq = [s for s in my_model.steps.keys() if s != 'Initial']
if len(_step_seq) != BUCKLE_STEP_NO or _step_seq[BUCKLE_STEP_NO - 1] != MODE_STEP_NAME:
    raise RuntimeError("스텝 순서 불일치: %s (기대: %s) — HF 의 MODE_SOURCE_STEP=%d 와 어긋난다."
                       % (_step_seq, ['Step-GlobalTension', MODE_STEP_NAME], BUCKLE_STEP_NO))
_hf_step_no = None
try:
    with open(os.path.join(_HERE, 'run_abaqus.py'), 'r', errors='replace') as _f:
        for _ln in _f:
            if _ln.lstrip().startswith('MODE_SOURCE_STEP'):
                _hf_step_no = int(_ln.split('=', 1)[1].split('#')[0].strip())
                break
except Exception as _e:
    emit("[run_abaqus_mode] HF 상수 교차확인 건너뜀: %s" % _e)
if _hf_step_no is not None and _hf_step_no != BUCKLE_STEP_NO:
    raise RuntimeError("HF 의 MODE_SOURCE_STEP=%d <> 모드 소스의 BUCKLE_STEP_NO=%d — 모드표 스텝이 어긋난다."
                       % (_hf_step_no, BUCKLE_STEP_NO))
emit("[run_abaqus_mode] 스텝 순서 확인: %s (BUCKLE_STEP_NO=%d, HF MODE_SOURCE_STEP=%s)"
      % (_step_seq, BUCKLE_STEP_NO, _hf_step_no))

# -------------------------------------------------------------
# 7. (삭제) *NODE FILE 삽입 — *BUCKLE 은 .fil 출력이 금지되므로(실측 7025) 요청 자체가 무효다.
#    좌굴모드는 ODB 모드 프레임에만 있으며, 아래 8단계에서 모드표로 뽑는다(키워드 조작 불필요).
# -------------------------------------------------------------
emit("[run_abaqus_mode] *NODE FILE 삽입 없음 — *BUCKLE 좌굴모드는 ODB 에서 뽑는다")

# -------------------------------------------------------------
# 8. 실행 -> 좌굴모드 확인 -> 모드표 생성 (이 스크립트가 끝나면 HF 가 바로 돌 수 있다)
# -------------------------------------------------------------
emit("=" * 74)
emit("모드 소스 런: %s  (선형 좌굴해석 %s, 프리텐션 %s, 꼭짓점 당김 %.3e m / SIGMA0 %.1f Pa / 안정화 %.4f)"
      % (JOB_NAME, MODE_STEP_NAME, PRETENSION_MODE, DISP_GLOBAL, SIGMA0, MODE_STABILIZATION))
emit("  모델: 클램프 없음 / 케이블 3개 / 꼭짓점 강체패치 — 구동 = %s"
      % ('3꼭짓점 = 논문 좌굴모델' if PRETENSION_MODE == 'paper3' else '아래 2꼭짓점'))
emit("  요청 고유값 %d개 / 기저벡터 %d / 최대반복 5000 (SUBSPACE)" % (N_EIG_BUCKLE, BUCKLE_VECTORS))
emit("=" * 74)

ok = run_job_safely(JOB_NAME)
print_job_diag(JOB_NAME)

# base state 실측(클램프 없는 상태의 압축 영역 비율) — 결과 해석용, 토큰 소량
if RUN_BASE_STATE_PROBE:
    try:
        _probe = os.path.join(_HERE, 'base_state_probe.py')
        if os.path.exists(_probe) and os.path.exists(JOB_NAME + '.odb'):
            emit("[MODE] base state 측정: %s.odb" % JOB_NAME)
            subprocess.call('abaqus python "%s" "%s" Step-GlobalTension'
                            % (_probe, JOB_NAME + '.odb'), shell=True)
        else:
            emit("[MODE] base state 측정 건너뜀 (probe=%s)" % os.path.exists(_probe))
    except Exception as _e:
        emit("[MODE] base state 측정 실패(무시): %s" % _e)

# --- 모드 확인: 두 경로를 교차 확인한다 -------------------------------------
#   (1) .dat MODE NO 표  = Abaqus 가 수렴시킨 고유값 개수(λ 포함)
#   (2) ODB 모드 프레임  = HF 가 실제로 먹는 데이터 -> 이걸 모드표로 만든다
n_dat = count_modes(JOB_NAME + '.dat', JOB_NAME + '.msg')
lams = []
try:
    _dat = os.path.join(os.getcwd(), JOB_NAME + '.dat')
    if parse_eigenvalues is not None and os.path.exists(_dat):
        with open(_dat, 'r', errors='replace') as _f:
            lams = parse_eigenvalues(_f.read())
except Exception as _e:
    emit("[MODE] λ 읽기 실패(무시): %s" % _e)
emit("[MODE] .dat 고유값 %d개 %s" % (len(lams), ['%.6e' % v for v in lams[:6]]))

n_modes = 0
# 왜 이 줄이 필요한가: 표가 안 만들어졌을 때 '잡이 실패해서(ok=False)'인지 '표 생성이 예외로 죽어서'인지
#   콘솔/로그로 즉시 구분해야 한다(실측: 표가 없는 채로 HF 가 중단되는 일이 있었다).
emit("[MODE] 잡 완료 ok=%s / WRITE_MODE_TABLE=%s / ODB=%s.odb 존재=%s / .dat 고유값=%d개"
     % (ok, WRITE_MODE_TABLE, JOB_NAME, os.path.exists(JOB_NAME + '.odb'), n_dat))
if ok and WRITE_MODE_TABLE:
    try:
        from aba_mode_from_odb import write_mode_table
        from aba_imperfection import model_fingerprint_line
        _fp_line = model_fingerprint_line(os.path.join(_HERE, 'run_abaqus_mode.py'))
        n_modes, _msgs = write_mode_table(JOB_NAME + '.odb', MODE_STEP_NAME, MODE_TABLE, N_MODES,
                                          INSTANCE_NAME,
                                          extra_header=[_fp_line] if _fp_line else None)
        emit("[MODE] 재사용 지문 %s" % (_fp_line or '없음(생략)'))
    except Exception as _e:
        emit("!!! 모드표 생성 실패: %s" % _e)
        emit("[ERROR-EN] mode table write FAILED: %s: %s" % (type(_e).__name__, _e))
        emit("    (수동 확인: abaqus python aba_mode_from_odb.py %s.odb %s %s %d %s)"
              % (JOB_NAME, MODE_STEP_NAME, MODE_TABLE, N_MODES, INSTANCE_NAME))
if n_dat != n_modes:
    emit("!!! WARNING: .dat 고유값 %d개 <> ODB 모드 프레임 %d개 -> 어느 쪽이 맞는지 확인 필요"
          % (n_dat, n_modes))

emit("=" * 74)
if not ok or n_modes <= 0:
    emit("RESULT:MODE_FAIL — ODB 모드 %d개 / .dat 고유값 %d개 (job_ok=%s)." % (n_modes, n_dat, ok))
    emit("  원인 판정 순서(라이선스 0):")
    emit("    1) ITERATION 2 이후에 양수 고유값이 하나도 없는가? -> base state 프리텐션 크기 문제다.")
    emit("       실측 추이: 양수 4개(5e-5 m) -> 3~4개(1e-7 m) -> 0개(1e-3 m, 면내 68.7 kPa).")
    emit("       [R-13] 의 '목표 7000 Pa 대비 배율'로 DISP_GLOBAL 을 1.8e-5 m 근처로 맞춘다(래더 L1).")
    emit("    2) [DIAG] 의 base state 음수 고유값 개수는 참고용이다 — 오라클은 76개여도 성공했다.")
    emit("       다만 개수가 100 에 근접하면 subspace 창이 부족하므로 N_EIG_BUCKLE 을 올린다(래더 L4).")
    emit("  래더(1줄씩, 1회 ~2분):")
    emit("    L1  DISP_GLOBAL 스케일 — 1.8e-5(오라클 검증) <-> base_state_probe 배율 보정값")
    emit("    L2  SIGMA0 = 700.0 (케이블 런 값) + MODE_STABILIZATION = 0.0005")
    emit("    L3  PRETENSION_MODE='corner2' 로 되돌려 3꼭짓점 대칭 효과와 교란 분리")
    emit("    L4  N_EIG_BUCKLE = 200 + BUCKLE_VECTORS = 500  (요청 수 > 음수 고유값 수)")
    emit("    L5  0.4 m 강체 패치 + KINEMATIC 커플링으로 로드 분산 — 응력집중 완화")
    emit("        (오라클 run_abaqus_cable.py:144 create_rigid_patch 가 쓰는 장치. 실측 응력비 maxP/mean = 52배)")
    sys.exit(1)

emit("RESULT:MODE_OK — 좌굴모드 %d개. ODB=%s" % (n_modes, os.path.abspath(JOB_NAME + '.odb')))
emit("  모드표: %s" % os.path.abspath(MODE_TABLE))
emit("  스텝=%s / λ(하중계수) %s" % (MODE_STEP_NAME, ['%.6e' % v for v in lams[:N_MODES]]))
emit("  참고: Galhofo2022 λ1..4 = 3.18260/3.18295/3.18341/3.18374e-4 (스프레드 0.036%)")
emit("  다음 단계: abaqus cae noGUI=run_abaqus.py -- HF <x_c> <d_c>"
      "   (IMPERFECTION_MODE='odb_table' 가 이 모드표를 노드 좌표 섭동으로 주입)")
emit("=" * 74)
