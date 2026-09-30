# Riks 전환 계획 (SolarSail-MFBO 포스트버클링)

2026-09-30 작성. 현행 정적 전략을 Riks 기반으로 바꿀 때 **무엇을 왜 어떻게 고치는가** 를 실제 코드 위치와 함께 정리한다.
현행 정적 조합은 이미 전 램프를 완주했고(990 증분 / 5972 s / 오류 0, 추력손실 4.59 %) 검증 4종을 통과했다.
따라서 Riks 전환은 **필요가 아니라 선택**이며, 정당화되는 축은 하나뿐이다: **속도**.

## 0. 판정 요약

| 축 | 판정 | 근거 |
|---|---|---|
| **필요성**(경로 추적) | **반증** | 변위 제어 정적 런이 전 램프를 완주 = 그 경로에 스냅백/폴드 없음. 게다가 장력 이력이 두 런 모두 단조 |
| **정확도** | **철회** | Riks 는 경로 추적 도구지 정확도 기구가 아니다. 정확도는 수렴기준·이산화가 결정 |
| **속도** | **미확정 - 유일하게 남은 축** | 아래 측정 참조 |

### 속도에 대한 기존 측정 (주목할 만하다)

| | 정적 (현행, 완주) | Riks 시도 (완화 없음, 2026-09-29) |
|---|---|---|
| 5 % 도달까지 증분 | 약 90 | **38** |
| 5 % 도달까지 시간 | 약 9 분 | **4 분 15 초** |
| 증분당 비용 | 6.03 s | 6.71 s |
| 전 구간 | **완주** (990 증분 / 5972 s) | **미완주** (λ 0.05 에서 아크길이 붕괴) |

즉 **5 % 까지는 Riks 가 약 2.2 배 빠르다**. 그런데 그 Riks 런은 **완화를 전혀 주지 않았고**(기본 C_n^a, 감쇠 없음),
정지 원인은 정적 런을 막았던 것과 **같은 C_n^a 함정**이다. 따라서 "Riks 가 빠르다"는 **아직 주장이 아니라 가설**이다.

## 1. 상실과 위험 (먼저 봐야 할 것)

1. **감쇠를 쓸 수 없다** - `*STABILIZE` 는 RIKS 와 양립 불가(실측). 그런데 감쇠 0.003 은 **현재 완주 조합의 일부**다.
   전환은 검증된 조합을 **덜 검증된 상태로 되돌린다**.
2. **논문 방법과 어긋난다** - Galhofo 는 원문에 "the geometrically nonlinear incremental analysis was carried out,
   using the **Newton-Raphson method** ... The **stability factor** was used in order to overcome numerical
   instabilities" 라고 명시한다. Riks 로 가면 **참조와 다른 솔버**가 되어 재현성 주장에서 불리해진다.
3. **line search / C_n^a 완화의 합법성이 미검증** - Riks 스텝 안에서 `*Controls, parameters=field` 와
   `parameters=line search` 가 허용되는지 문서에 명시가 없다(입력단계에서 수 초만에 판정 가능).
4. **평가 코드 재검증 필요** - 미달 상태 처리와 종료 조건이 λ 기준으로 바뀐다.

## 2. 코드 수정 목록

| 파일 | 위치 | 무엇을 |
|---|---|---|
| `run_abaqus.py` | 1126-1156 (`my_model.StaticStep(name='Step-Postbuckle', ...)`) | Riks 스텝으로 교체. **`stabilizationMagnitude` / `stabilizationMethod` / `continueDampingFactors` / `adaptiveDampingRatio` 네 인자를 반드시 제거**(RIKS 와 양립 불가) |
| `run_abaqus.py` | 1138-1141 (`initialInc=1e-4, minInc=1e-8, maxInc=0.1, maxNumInc=10000`) | 증분 인자 → Riks 8 항목: `Δl_in, l_period, Δl_min, Δl_max, λ_end, node, dof, value`. **λ_end 는 5번째**(3번이 아니다 - 입력단계 실패로 확인된 함정). `maxNumInc` 는 유지 |
| `run_abaqus.py` | 1254-1269 (`[STEP-CHECK]` 되읽기) | 지금은 `*Static` 줄만 찾는다. `*Static, riks` 도 잡도록 필터 확장(리크스 전환이 실제 반영됐는지 확인하는 유일한 장치) |
| `eval_abaqus.py` | 246-252 (`current_time = last_frame.frameValue; if current_time < 0.99`) | Riks 에서 `frameValue` 는 **LPF** 이므로 의미는 살아있지만, `λ_end` 를 직접 읽어 **그 값과 비교**하도록 바꾼다(부분 종료를 설계로 쓸 때 필수) |
| `check_model_consistency.py` | 파일 문자열 검사 목록 | **변경 불필요** - 검사 대상 `stabilizationMagnitude=0.0002,` 는 GlobalTension/ClampTension 스텝의 값이라 Postbuckle 만 바꾸면 그대로 남는다(코드로 확인) |
| `aba_imperfection.py`, `run_abaqus_mode.py` | - | **변경 불필요** - 모드 소스는 여전히 `*BUCKLE` 이고 임퍼펙션은 노드 라벨 매핑이라 무관 |

## 3. 구현 경로 두 가지

**A. `.inp` 패처 경로 (검증됨, 즉시 사용 가능)**
`riks_patch_input.py` 가 이미 3 줄 교체를 수행하고 프로브와 실모델에서 모두 통과했다.
```
python riks_patch_input.py aba\HF_Postbuckle.inp aba\HF_riks.inp Step-Postbuckle --riks --relax-corr --line-search
```
MFBO 루프에 넣으려면 `run_abaqus.py` 가 잡을 제출하기 **전에** 같은 로직을 파이썬으로 호출하면 된다(패처의 함수를 import).

**B. CAE API 경로 (`StaticRiksStep`)** - 자동화의 정공법이지만 **인자 목록이 미검증**이다.
추측 금지: CAE 커맨드 레퍼런스로 확인하거나 소형 CAE 프로브로 확정할 것. 확인 없이 인자를 늘리거나 빼면
정확도 검증 없이 죽는다(이 프로젝트에서 세 번 반복된 실패 유형).

## 4. 검증 사다리

| 단계 | 비용 | 명령 / 방법 | 판정 |
|---|---|---|---|
| **V0 합법성** | 수 초 | 패처로 `--riks --relax-corr --line-search` 생성 → `abaqus job=<n> input=<path> interactive` | **PASS (2026-09-30)** - 입력 프로세서가 오류 없이 끝나고 `.dat` 에 `***ERROR` 0. 즉 **Riks 스텝에서 `*Controls, parameters=field` + `parameters=line search` 둘 다 합법**이다(문서에 명시가 없던 사실). 분석단계 실패는 `.dat` 가 아니라 **`.msg`** 를 봐야 한다. |
| **V1 속도** | 수 분 | `--lpf=0.05` 로 5 % 도달 (`0.01, 1.0, 1e-07, 0.025, 0.05`) | 정적 대비 **증분 수·시간** 비교. 유의미하게 적으면 V2, 아니면 **여기서 중단**. 주의: V0 는 lambda_end=1.0 / dl_min=1e-5 로 만들어져 이전 실패 모드(아크길이 바닥)에 그대로 노출된다 - V0 의 분석 실패는 속도 판정이 아니므로, 속도는 반드시 V1 으로 본다(dl_min=1e-7 로 100배 낮다). |
| **V2 완주** | 30~90 분 | `--lpf=1.0` | 완주 여부 + 도달 λ |
| **V3 물리** | 수 분 | ① 같은 λ 에서 `base_state_probe.py <odb> Step-Postbuckle <λ>` 로 정적 완주 해와 대조 ② `eval_abaqus.py` HF 지표 ③ `energy_check.py`(감쇠가 없으니 ALLSD≈0 예상) ④ `node_probe.py` / `wrinkle_spectrum.py` | 정적 해와 **적분량 ≤2 %** 면 대안으로 성립 |

**사전 고정 판정 규칙**: V1 에서 증분·시간이 줄지 않으면 **Riks 는 채택하지 않는다**. 채택 조건은 "더 빠르고, V3 에서 정적 해와 일치" 둘 다다.

## 5. 참고: Riks 데이터 라인과 종료 조건

- 데이터 라인 8 항목 순서: `Δl_in, l_period, Δl_min, Δl_max, λ_end, node, dof, value`.
- 종료 3 방식: 모니터 절점·자유도·도달 값(6,7,8) 또는 `λ_end`(5). 후자는 **부분 종료를 설계로 지정**할 수 있어
  정적 스텝의 `timePeriod` 보다 의미가 분명하다.
- `self-explanatory` 원칙: Riks 는 `STABILIZE` 와 양립 불가하고 **마지막 스텝**이어야 한다(가드 구현됨).
