# -*- coding: utf-8 -*-
"""core/constants.py — HF와 좌굴 스크립트가 **문자 그대로** 공유해야 하는 모델 상수.

왜 한 곳에 두나
  `check_model_consistency.py` 가 두 스크립트의 같은 문자열을 비교해 왔다. 값이 두 곳에
  있으면 한쪽만 고쳐도 조용히 어긋나고, 그 사실은 라이선스를 쓴 뒤에야 드러난다.
  여기 값이 **유일한 출처**이고 두 스크립트는 import 한다.

지문(fingerprint) 주의
  `aba_imperfection.source_model_fingerprint` 는 이 파일의 `KEY = 값` **텍스트**를 해싱한다
  (buckle 시드표 헤더의 `buckle_fp`). 값뿐 아니라 **우변 표기**가 바뀌면 지문이 바뀐다.
  예: `5.0e-6` -> `5e-6` 은 같은 값이지만 지문은 달라진다. 표기를 바꾸려면 시드표를 재생성해야 한다.
"""

#   요소 코드는 Abaqus 심볼(S4/S3)이다. 이 모듈은 Abaqus 프로세스 안에서만 import 된다.
from abaqusConstants import *          # noqa: F401,F403

# ---- 형상 ----
BASE = 20.0                            # m
HEIGHT = 10.0                          # m

# ---- 요소/격자 ----
ELEM_CODE_QUAD = S4                    # 1차 셸(S4). 막 계열(M3D4/M3D3)은 기각(notes §4·§24).
ELEM_CODE_TRI = S3                     # 3절점 1차 셸
SEED_DIV = 200.0                       # seed = BASE/SEED_DIV -> 약 1.82만 요소 (실측 2026-09-28)

# ---- 두께 ----
THICKNESS = 5.0e-6                     # m. 2.5e-6 은 cable 변형으로 수렴이 매우 어렵다.

# ---- 클램프 케이블 ----
CLAMP_CABLE_LEN = 0.5                  # m - HF/좌굴/모드가 같은 값이어야 한다(검사기 계약).
