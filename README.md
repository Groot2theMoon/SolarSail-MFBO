# Solar Sail MFBO Framework

이 저장소는 **Solar Sail의 형상 최적설계**를 위한 **Multi-Fidelity Bayesian Optimization (MFBO)** 프레임워크를 포함하고 있습니다.

계산 비용이 매우 높은 **비선형 포스트버클링 해석(High-Fidelity)** 횟수를 줄이기 위해, 계산이 저렴한 **선형 Static 해석**을 활용하여 최적의 설계 변수(클램프 부착 위치 및 변위량)를 탐색합니다.

## Key Concepts

### 1. 최적화 문제 정의
*   **목표:** 태양광 압력(SRP) 모델 기반 **추력 손실(Thrust Loss) 최소화**.
*   **설계 변수:**
    *   `x_c` ($x_1$): 클램프 케이블 부착 위치 (0.05 $\le$ $x_c$ $\le$ 0.95).
    *   `d_c` ($x_2$): 클램프 케이블 인장 변위 비율. (0 $<$ $d_c$ $\le$ 1.0)
*   **Minimization Strategy:** BoTorch는 기본적으로 최대화를 수행하므로, 코드 내부에서 결과값(추력 손실량)에 `-1`을 곱하여 최적화를 수행합니다.

### 2. Multi-Fidelity 구성
| Fidelity | Abaqus Method | Output Metric | Cost (Relative) |
| :--- | :--- | :--- | :--- |
| **LF (Low)** | Linear / Nonlinear Static | **Compression Area Ratio** (압축 응력 발생 면적비) | 1.0 |
| **HF (High)** | Post-Buckling (Static, General) | **Thrust Loss** (실제 형상 기반 추력 감소량) | 10.0 |

> **가설:** "선형 해석(LF)에서 압축 응력이 넓게 분포하는 디자인은, 실제 포스트버클링 해석(HF)에서도 주름(Wrinkle)이 크게 발생하여 추력 성능이 떨어질 것이다."

* eval_abaqus.py docstring에서도 밝히듯이, 실제로는 3가지의 lf 지표를 제안하므로, 3개 중 원하는 지표를 선택하여 mfbo를 진행할 수 있습니다.
* 수정을 위해서는 mfbo.py 의 get_abaqus 함수내에서 parsed_data[0] 를 1 or 2 로 수정하면 됩니다.
---

## System Architecture

이 프로젝트는 서로 다른 두 개의 Python 환경이 `subprocess`를 통해 통신하도록 작성되었습니다.
1.  **Host Process (MFBO)**: `mfbo.py`
    *   **Environment:** Python 3.x (Anaconda recommended).
    *   **Libraries:** `PyTorch`, `BoTorch`, `GPyTorch`, `WandB`.
    *   **Role:** 베이지안 최적화 알고리즘 구동, GPU 연산, 다음 해석 지점 추천.
2.  **Solver Process (Abaqus)**: `run_abaqus.py` & `eval_abaqus.py`
    *   **Environment:** Abaqus Internal Python (Python 2.7).
    *   **Role:** FEM 모델 생성, 해석 수행, ODB 결과 파싱.
    *   **Communication:** `subprocess` 호출 -> `stdout` 및 텍스트 파일(`extraction.txt`)을 통해 결과값 반환.

---

## Quick Start

### 1. 필수 요구 사항
*   **Abaqus:** 설치되어 있어야 하며, 터미널에서 `abaqus` 명령어로 실행 가능해야 함. (Simulia License 필요)
*   **Python 3.x 환경:**
    ```bash
    conda create -n mfbo_env python=3.9
    conda activate mfbo_env
    pip install torch botorch gpytorch scipy wandb
    ```

### 2. 실행 방법
* 터미널에서 `mfbo.py`만 실행하면 됩니다. 나머지는 자동으로 호출됩니다.
* 각각의 파일을 별개로 실행하고 싶을 경우 (run_abaqus.py 등) 각 파일의 docstring을 참고 바랍니다.

```bash
# 1. 가상환경 활성화
conda activate mfbo_env

# 2. 실행
python mfbo.py
```

### 3. WandB 설정
*   `mfbo.py` 내의 `WANDB_PROJECT` 변수를 본인의 프로젝트명으로 수정하세요.
*   실행 시 WandB 로그인을 요구할 수 있습니다.

---

## 파일 상세 설명

### 1. `mfbo.py` 
*   **역할:** 전체 최적화 루프 제어.
*   **주요 로직:**
    *   `load_checkpoint()`: 중단된 실험 자동 재개.
    *   `get_abaqus()`: `run_abaqus.py`를 CLI 명령어로 실행하고 결과를 파싱. **값의 부호 반전(-)** 로직이 여기에 포함됨.
    *   **Acquisition Function:** GIBBON (`qMultiFidelityLowerBoundMaxValueEntropy`) 사용.
*   **수정 포인트:** `N_ITERATIONS`, `LF_INIT`, `HF_INIT` 변수로 실험 규모 조절.

### 2. `run_abaqus.py` (ABAQUS Simulation)
*   **역할:** Abaqus CAE 모델링 및 해석 실행.
*   **실행 방식:** `abaqus cae noGUI=run_abaqus.py -- [Mode] [x1] [x2]`
*   **주요 로직:**
    *   **LF 모드:** `Step-GlobalTension` -> `Step-ClampTension` -> `Step-HighTension`.
    *   **HF 모드:** `Step-GlobalTension` -> `Step-ClampTension` -> `Step-Buckle` (모드 추출)
    *                  & `Step-GlobalTension` -> Imperfection 주입 -> `Step-Postbuckle` (비선형 해석).
*   **주의:** `*IMPERFECTION` 키워드 삽입 로직(삽입 위치 등)은 Abaqus 버전에 따라 민감할 수 있음.

### 3. `eval_abaqus.py` (Post-Processing)
*   **역할:** run_abaqus.py 로 생성된 ODB 파일을 열어 수치적 성능 지표 계산.
*   **주요 로직:**
    *   **LF Metrics:** 압축 응력 부피비(`lf1`), wrinkle 에너지 밀도(`lf2`), 응력 이방성(`lf3`).
    *   **HF Metric:** SRP 모델을 적용하여 변형된 메쉬의 법선 벡터를 적분, wrinkle에 의한 **Thrust Loss** 계산.
    *   **성능:** `bulkDataBlocks`와 `numpy`를 사용하여 Python 루프 속도 문제 개선.

---

## Optical Properties Calculation (Optional)

`eval_abaqus.py`에서 사용되는 추력 모델(SRP Model)의 상수인 반사율($R_0$)과 흡수율($A_0$)은 `RA_calc.py`를 통해 계산되었습니다.

*   **Source Data:**
    *   `nk_data.csv`: 해당 멤브레인(Kapton 등)의 파장별 굴절률($n$) 및 소멸 계수($k$).
    *   `sun_data.csv`: 태양 복사 스펙트럼 (Solar Irradiance Spectrum).
*   **Methodology:**
    *   Fresnel 방정식을 이용하여 파장별 반사율/흡수율 계산.
    *   태양 스펙트럼 강도(Irradiance)를 가중치로 하여 전체 파장 대역에 대해 적분(Weighted Integration).
*   **Usage:**
    ```bash
    python RA_calc.py
    ```
    *   코드 파일과 동일 폴더에 nk_data.csv, sun_data.csv가 있음을 확인한 후, 코드 실행으로 출력된 $R_0, A_0$ 값을 `eval_abaqus.py`의 상수로 입력하여 사용합니다.

---

## 주요 Troubleshooting

1.  **"abaqus command not found" 에러**
    *   실행 전 시스템 PATH에 Abaqus 실행 파일 경로가 등록되어 있는지 확인이 필요합니다.
    *   `mfbo.py`의 `command` 변수에서 `abaqus` 대신 전체 경로(예: `C:\SIMULIA\Commands\abaqus.bat`)를 입력해야 할 수도 있습니다.

2.  **Lock File 에러 (.lck)**
    *   Abaqus 해석이 비정상 종료되면 `.lck` 파일이 남아 다음 해석을 방해할 수 있습니다.
    *   `run_abaqus.py`에 자동 삭제 로직이 있지만, 권한 문제로 실패할 경우 수동으로 삭제해야 합니다.

3.  **라이센스 문제**
    *   `run_abaqus.py`는 `abaqus cae` 라이센스를 점유합니다. 사용 중인 라이센스 서버 상태를 확인하세요.

4.  **HF 해석 수렴 실패**
    *   포스트버클링 해석은 수렴이 어려울 수 있습니다. `eval_abaqus.py`는 시간이 `0.99` 미만일 경우 경고 메세지를 띄우지만, 마지막 프레임을 사용해서 값은 계산합니다. 필요시 `mfbo.py`에서 `1e6` 패널티 처리를 강화.

---

## References
*   **Simulation Model:** "참고논문" 폴더 참조.

*   **Libraries:** [BoTorch Documentation](https://botorch.org/), [Abaqus Scripting Reference](http://130.149.89.49:2080/v2016/books/ker/default.htm)


