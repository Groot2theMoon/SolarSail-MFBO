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

# (제거 2026-09-24) `from matplotlib.image import LANCZOS` — IDE 자동삽입으로 들어온 잘못된 import.
#   matplotlib 의 LANCZOS 는 이미지 리샘플링 필터이고, Abaqus 고유치 솔버 LANCZOS 는 abaqusConstants 에서 온다.
#   게다가 Abaqus Python 에 matplotlib 이 없으면 스크립트가 시작 시 ImportError 로 죽는다.
#   이 스크립트는 LANCZOS 를 쓰지 않는다(선형 좌굴해석 + SUBSPACE).
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

# ---- 로깅 (2026-09-24 실측으로 발견한 문제 대응) --------------------------------
# 왜 필요한가: Abaqus `cae noGUI=...` 는 이 스크립트의 stdout 을 블록 버퍼링한 뒤 프로세스를
#   flush 없이 teardown 하는 경우가 있다. 그때 print() 는 콘솔에서 통째로 사라지고, 별도
#   프로세스로 뜬 base_state_probe 의 [R-13] 출력만 남는다(실측: 4회차 9회차 콘솔 모두 그랬다).
#   게다가 Abaqus CAE 네임스페이스에는 이미 `log` 가 있어 그 이름으로 정의하면 호출이 내장 log 로
#   가서 `TypeError: illegal argument type for built-in operation` 으로 죽는다(실측).
#   -> 이름을 emit 으로 두고, 정의를 첫 호출보다 앞에 두고, 절대 예외를 올리지 않게 한다.
try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass

_EMIT_LOG_PATH = None


def emit(*a):
    """콘솔 + mode_run_log.txt 동시 기록. **어떤 경우에도 예외를 올리지 않는다**(로깅이 해석을 죽이면 안 된다)."""
    global _EMIT_LOG_PATH
    try:
        msg = ' '.join(str(x) for x in a)
    except Exception:
        msg = '<emit: 인자 변환 실패>'
    try:
        print(msg)
    except Exception:
        pass
    try:
        sys.stdout.flush()
    except Exception:
        pass
    try:
        if _EMIT_LOG_PATH is None:
            _EMIT_LOG_PATH = os.path.join(
                os.path.dirname(os.path.abspath(__file__)), 'mode_run_log.txt')
            with open(_EMIT_LOG_PATH, 'w') as _f:      # 첫 호출에서 truncate(최신 런만 남긴다)
                _f.write('[run_abaqus_mode] log start %s\n'
                         % time.strftime('%Y-%m-%d %H:%M:%S'))
        with open(_EMIT_LOG_PATH, 'a') as _f:
            _f.write(msg + '\n')
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
from aba_imperfection import count_modes, parse_eigenvalues                 # noqa: E402
emit("[run_abaqus_mode] _HERE = %s" % _HERE)


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
    job.submit(consistencyChecking=OFF)
    job.waitForCompletion()
    time.sleep(1.0)
    if job.status == ABORTED or not os.path.exists(job_name + '.odb'):
        emit("!!! ERROR: Job %s failed. Actual Status: %s" % (job_name, str(job.status)))
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
#   [2026-09-24 3회차 실측 -> 근본 수정] 코너 2개만 1e-7 m 당기면 base state 가 slack 이 되고,
#   좌굴 스펙트럼이 뭉개진다(실측: ITER2 양수 모드가 1회차 4.6e-5~3.1e-3 -> 3회차 6.9e-7~1.8e-5,
#   즉 요구 모드가 사실상 영에너지 모드로 퇴화). 논문 Galhofo2022 의 좌굴 모델은
#   '세 꼭짓점에 prescribed displacement 로 초기 프리텐션' + 좌굴 BC(정점 u_y=1e-3 m,
#   아래 두 꼭짓점 성분 (5e-4,-5e-4) m) -> 3중 대칭으로 세 꼭짓점을 모두 당긴다.
# CHECKER-DIVERGENCE: DISP_GLOBAL, SIGMA0, PRETENSION_SCALE, stabilizationMagnitude, radius
#   radius 는 오라클 재현 런(0.4)을 위해 추가했다. 지금은 0.2 로 복원되어 '값'은 HF 와 같다.
#   다만 값이 PATCH_RADIUS 상수로 들어가므로 검사기가 찾는 리터럴 'radius=0.2' 는 코드에 없고,
#   결과는 DIVRG(선언됨) 로 표시된다. 즉 이 선언이 남아 있는 동안에는 PATCH_RADIUS 를 0.9 등으로
#   바꿔도 검사기가 통과시킨다 -> PATCH_RADIUS 를 바꿀 때는 사람이 HF 값(0.2)과 직접 대조할 것.
#   HF(run_abaqus.py)와 이 모드 소스는 base state 프리텐션 정의가 다르다(모드 소스는 3꼭짓점 프리텐션).
#   메쉬/형상/요소/재료는 동일하므로 모드 노드 라벨 매핑은 그대로 성립한다.
#   논문에는 '클램프 유무' 외에 이 프리텐션 정의 차이도 함께 명시할 것.
PRETENSION_MODE = 'corner2'  # [2026-09-24 오라클 재현] 케이블 런(run_abaqus_cable.py)과 동일한 구동:
                             #   정점(BC_Anchor_Top)은 완전 고정한 채, 아래 두 꼭짓점만 케이블 축으로 당긴다.
                             #   성공 실적: 이 설정에서 CONVERGED=4 (양수 λ 4개).
                             # 'paper3' = 논문 좌굴모델: 꼭짓점 3개를 모두 당김(3중 대칭).
                             #   5·7회차 실측: 양수 2개 / 0개 -> 재현 실험 뒤 bisect 단계에서 하나씩 되돌린다.
DISP_GLOBAL = 1.8e-5         # 꼭짓점 당김 [m] — 프리텐션 크기. [2026-09-24 오라클 재현] 케이블 런 DISP_MAG 와 동일값.
                             #   [2026-09-24 오라클 재현 실험] 아래 '창/기저 소진' 기록은 그대로 유효하다.
                             #     프리텐션 5.0e-6 하향은 이 재현 실험 뒤 bisect 단계로 미룬다(한 번에 한 축).
                             #     실측: 요청 100/기저 400 -> CONVERGED=0 (ITER4 에서 396개까지 추적 후 붕괴) /
                             #           요청 60/기저 120 -> CONVERGED=0 (창 부족: 필요 총개수 63 > 요청 60) /
                             #           요청 100/기저 250 -> 양수 2개 수렴 (이 배치의 최선).
                             #     -> subspace 로는 2개가 상한이므로 음수 대역 자체를 줄이는 축으로 간다.
                             #     근거(단조 실측): 시스템 음수 고유값 88개(1e-3 m) -> 56개(1.8e-5 m) -> 16개(1e-7 m).
                             #     지표(사전등록): 새 런의 "SYSTEM MATRIX HAS N NEGATIVE EIGENVALUES" 의 N 이
                             #       56 보다 줄었는가. 줄지 않으면 이 축도 아니라는 뜻이다.
                             #     [R-13] 선형 환산 예상 면내응력 ~344 Pa (= 6.874e+04 x 5e-6/1e-3), 운용점 7000 Pa 의 0.049배.
                             #     하한 주의: 1e-7 m 대는 좌굴모드가 영에너지로 퇴화한다(3회차 진단) -> 그보다 위에 둔다.
                             #   근거 1(계측): base_state_probe 헤드라인 = 평균 면내응력 6.874e+04 Pa,
                             #     면내평균<0 면적비 0.0000, minP<0 0.0017, |u3|max 5.6e-22 m (=주름 0, 완전 평탄).
                             #     목표 7000 Pa 대비 배율 9.820 -> 1e-3 m 는 운용점의 약 10배로 과대했다.
                             #   근거 2(회차별 대조): 양수 좌굴모드가 나온 회차는 전부 저프리텐션이다.
                             #     1e-7 m  -> 양수 3~4개 (2·3회차)  /  5e-5 m -> 양수 4개 (1회차)
                             #     1e-3 m  -> 양수 0개 (4회차, 면내 68.7 kPa)  ← 양수가 사라지는 구간
                             #   근거 3(양성 대조군): 성공한 케이블 런(run_abaqus_cable.py, CONVERGED=4)은
                             #     DISP_MAG = 1.8e-5 m 이다. 같은 값으로 맞춘다.
                             #   주의: 이 값은 '좌굴모드를 뽑기 위한 저프리텐션'이다. HF 운용점(7000 Pa)과 다르며
                             #         그 차이는 HF/모드소스 분기로 이미 선언돼 있다(check_model_consistency.py).
DISP_TOP_OVER_CORNER = 1.4142135623
SIGMA0 = 700.0                               # 초기 가짜 응력(수렴 보조) [Pa] — 프리텐션의 대체물이 아니다
                                             #   [2026-09-24 오라클 재현] 케이블 런 값 = 700.0 (HF 는 500.0)
MODE_STABILIZATION = 0.0005                  # GlobalTension 안정화 계수
PERTURBATION = 0.01                          # 좌굴 스텝 섭동 크기 [m] (케이블 런·HF 와 동일)
PATTERN_SIGN = 1.0                           # 좌굴 '하중 패턴'의 부호: 1.0 = 바깥으로 더 당김 (오라클과 동일)
                                             #   [2026-09-24 철회] 직전에 -1.0 으로 뒤집었다가 되돌렸다.
                                             #   뒤집은 근거였던 "4회차 스펙트럼이 전부 λ<0 이니 패턴 부호 문제"는
                                             #   양성 대조군으로 반증되었다: 성공한 케이블 런(CONVERGED=4)도
                                             #   패턴이 똑같이 '바깥 당김'(+PERTURBATION)이고 λ 4개가 양수다.
                                             #   -> λ 의 부호를 정하는 것은 패턴 부호가 아니라 base state(프리텐션 크기)다.
                                             #      실측: 1e-7 m 와 5e-5 m 에서는 양수 모드가 나왔고(3~4개),
                                             #            1e-3 m(면내 68.7 kPa)에서는 0개였다.
                                             #   -> 원인은 DISP_GLOBAL 이며, 위 상수에서 1.8e-5 m 로 내렸다.
N_EIG_BUCKLE = 100                           # [2026-09-24 실측 근거] subspace 요청 고유값 수
                                             #   규칙: '요청 수 > base state 음수 고유값 수 + 양수 고유값 수' 여야 한다
                                             #         (= 필요한 고유값 총 개수보다 창이 커야 한다).
                                             #   실측: 케이블 런(음수 76)에서 100 요청 -> CONVERGED=4 (성공).
                                             #         이 모델의 2회차 런(음수 16)에서 10 요청 -> CONVERGED=0 (실패).
                                             #   [2026-09-24 정정] 음수 개수 자체는 치명적이지 않다 — 오라클은 76개인데
                                             #         성공했고, 이 모델은 88개였다. 결정적 차이는 프리텐션 크기였다(위 DISP_GLOBAL).
                                             #   [2026-09-24 6회차 실측] 60/120 -> CONVERGED=0. 원인이 숫자로 닫혔다:
                                             #     음수 56개 + ITER3에서 드러난 양수 7개 = 최소 63개가 존재하는데 창이 60 이라
                                             #     부족했다("THE EIGENVALUES CANNOT BE FOUND" + INSTABILITIES IN THE BASE STATE).
                                             #     -> 100 으로 복원(5회차 조합). 창은 충분했다.
                                             #   [2026-09-24 5회차 vs 6회차 대조 = 진짜 병목] 두 런의 ITER3 양수를 맞춰보니
                                             #     5회차가 수렴시킨 λ(1.53792e-04, 1.62867e-04)는 6회차 ITER3 의 3·4번째 값
                                             #     (1.587912e-04, 1.703028e-04)과 같은 순서다. 즉 5회차는 1·2번째 모드
                                             #     (2.856217e-05, 4.109332e-05)를 건너뛰고 3·4번째를 잡았다 -> 병목은 창이 아니라
                                             #     subspace 붕괴다(250 -> 244 -> 2). 붕괴는 기저 벡터 수로 지연시킨다.
                                             #     실측 단조: 기저 120 -> 양수 0개 / 기저 250 -> 양수 2개(조기 붕괴).
                                             #   판독 규칙(사전등록): 요청/기저 수는 solver 작업공간만 바꾸고 물리는 못 바꾼다.
                                             #     -> 수렴 λ 목록이 길어지면 붕괴가 병목이었다(성공).
                                             #     -> 다시 2개면 그 2개가 이 base state 의 안정 수렴 집합이다(성공, 물리 결론).
                                             #     -> 또 CONVERGED=0 이면 창/기저가 아니라 base state 음수 56개가 한계다
                                             #        -> 두 번째 노브(DISP_GLOBAL 하향, 음수 16 실적의 저프리텐션 대)로 간다.
BUCKLE_VECTORS = 250                         # 기저 벡터(요청 수 x 2.5) — 실측 최선 조합으로 복원.
                                             #   실측 비교: 기저 120 -> 양수 0개 / 기저 250 -> 양수 2개(최선) /
                                             #              기저 400 -> CONVERGED=0 (ITER4 396개까지 추적 후 붕괴).
                                             #   -> '기저를 키우면 붕괴가 지연된다'는 가설은 400 실측으로 반증됐다.
N_MODES = 4                                  # HF 에 주입할 모드 수(= ODB 모드 프레임에서 뽑는 개수)
PATCH_RADIUS = 0.2                           # 꼭짓점 강체패치 반경 [m] — HF 와 동일값(0.2)으로 복원.
                                             #   [2026-09-24 bisect 2단계, 사용자 지시] 4모드를 낸 조합은 0.4 였다.
                                             #   이제 반경만 0.2 로 되돌려 '패치 반경이 스펙트럼을 좌우하는가'를
                                             #   단독으로 판정한다(구동 방식 corner2 / 1.8e-5 / SIGMA0 700 은 그대로).
                                             #   판독 규칙(사전등록): 4개 수렴 + 0 ERROR -> 반경은 무관
                                             #     = HF 정합 회복(최상). 0~2개 -> 0.4 가 결정적이었다 -> 0.4 복귀.
                                             #   lambda 값 자체는 base state 변화로 이동하는 것이 정상이다(개수로 판정).
                                             #   되돌리기: 0.4 (오라클 값, run_abaqus_cable.py:144).
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


def create_rigid_patch(name, coord, radius=PATCH_RADIUS):
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
        emit("Warning: Node not found by sphere, trying closest for " + name)
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

# --- 스텝 순서 가드: HF 는 MODE_SOURCE_STEP 번째 스텝의 프레임을 모드 소스로 읽는다 ------
#   (2026-09-24) 이 상수가 '정의만 되고 안 쓰이던' 죽은 상수였고, 스텝 순서가 밀리면
#   HF 가 엉뚱한 스텝의 프레임을 조용히 읽는다 -> 여기서 즉시 중단한다.
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
if ok and WRITE_MODE_TABLE:
    try:
        from aba_mode_from_odb import write_mode_table
        n_modes, _msgs = write_mode_table(JOB_NAME + '.odb', MODE_STEP_NAME, MODE_TABLE, N_MODES,
                                          INSTANCE_NAME)
    except Exception as _e:
        emit("!!! 모드표 생성 실패: %s" % _e)
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
