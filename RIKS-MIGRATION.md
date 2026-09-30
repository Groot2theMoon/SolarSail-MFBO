# Riks 전환 계획 (SolarSail-MFBO 포스트버클링)

2026-09-30 작성. 현행 정적 전략을 Riks 기반으로 바꿀 때 **무엇을 왜 어떻게 고치는가** 를 실제 코드 위치와 함께 정리한다.
현행 정적 조합은 이미 전 램프를 완주했고(990 증분 / 5972 s / 오류 0, 추력손실 4.59 %) 검증 4종을 통과했다.
따라서 Riks 전환은 **필요가 아니라 선택**이며, 정당화되는 축은 하나뿐이다: **속도**.

## 0. 판정 요약

| 축 | 판정 | 근거 |
|---|---|---|
| **필요성**(경로 추적) | **확인 - 경로에 폴드가 있다** | Riks 런의 `.sta` LPF 열에 **음수 증분**이 수락된 증분에서 나타난다(`-0.01994`, `-0.002500`, `-0.0009962`). λ 가 이 모델에서 지정 변위 BC 값을 스케일하므로 구동 변위가 되돌아간다는 뜻이고, 정적 런이 완주한 것은 **감쇠가 폴드를 통과시킨 결과**다. "완주 = 스냅백 없음" 추론은 **철회** |
| **정확도** | **철회** | Riks 는 경로 추적 도구지 정확도 기구가 아니다. 정확도는 수렴기준·이산화가 결정 |
| **속도** | **측정됨 - 그러나 완주 불가** | 5 % 까지 약 2 배 빠르지만(아래) 두 런 모두 아크길이 바닥에서 정지했다. 2 배 이득은 완주 불가를 정당화하지 않는다 |
| **종합** | **Riks 종료 (문서상 부적합)** | 공식문서 Restrictions: 접촉 상실이 있는 포스트버클링에서 "usually will not work", 대신 점성 감쇠를 정적/동적 해석에 도입하라고 명시 = **현행 경로**. Riks 는 감쇠를 쓸 수 없으므로(§1-1) 문서가 권한 도구를 스스로 버린다 |

### 속도에 대한 기존 측정 (주목할 만하다)

| | 정적 (현행, 완주) | Riks 시도 (완화 없음, 2026-09-29) |
|---|---|---|
| 5 % 도달까지 증분 | 약 90 | **38** |
| 5 % 도달까지 시간 | 약 9 분 | **4 분 15 초** |
| 증분당 비용 | 6.03 s | 6.71 s |
| 전 구간 | **완주** (990 증분 / 5972 s) | **미완주** (λ 0.05 에서 아크길이 붕괴) |

그 뒤 **완화를 모두 주고**(`--riks --relax-corr --line-search`) 두 번 더 돌렸다:

| | 정적(완주) | Riks `Δl_max=0.1` | Riks `Δl_max=0.025` |
|---|---|---|---|
| 도달 λ | 1.00 (완주) | 0.0503 | 0.0129 |
| 증분 / 시간 | 990 / 5972 s | **73 / 278 s** | 67 / 264 s |
| line search 실효 | 833 + 647 평가 | **0 / 0** | **0 / 0** |
| 정지 사유 | - | `TIME INCREMENT REQUIRED IS LESS THAN THE MINIMUM SPECIFIED` | 같음 |

**C_n^a 함정 가설은 반증됐다**: 완화를 줘도 같은 λ 구간에서 멈춘다. 움직인 것은 **아크길이 창**뿐이고,
**더 좁게 조인 `Δl_max` 가 더 못 갔다**. 문서 기본값이 "상한 없음"이므로 상한을 걸면 자동 증분이 아크길이를
*늘릴* 여지가 사라진다(§6). 속도 이득(5 % 까지 약 2 배)은 실재하지만 **완주하지 못하므로 채택 근거가 되지 않는다**.

## 1. 상실과 위험 (먼저 봐야 할 것)

1. **감쇠를 쓸 수 없다** - `*STABILIZE` 는 RIKS 와 양립 불가(실측). 그런데 감쇠 0.003 은 **현재 완주 조합의 일부**다.
   전환은 검증된 조합을 **덜 검증된 상태로 되돌린다**.
2. **논문 방법과 어긋난다** - Galhofo 는 원문에 "the geometrically nonlinear incremental analysis was carried out,
   using the **Newton-Raphson method** ... The **stability factor** was used in order to overcome numerical
   instabilities" 라고 명시한다. Riks 로 가면 **참조와 다른 솔버**가 되어 재현성 주장에서 불리해진다.
3. **line search 는 Riks 에서 실효가 없다** - 문법은 합법으로 확정됐지만(V0 통과), `.msg` 의
   `ADDITIONAL RESIDUAL/OPERATOR EVALUATIONS FOR LINE SEARCHES` 가 **0/0** 이다(두 런 모두). `*CONTROLS` 와
   convergence-control 문서에 "Riks" 가 **0회** 등장하고 `Nls` 기본값이 Newton 스텝에서 0인 것과 일치한다.
   즉 완주 조합의 카드 하나를 Riks 는 실제로 쓰지 못한다.
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

## 공식 문서 원문 (2017 Abaqus Analysis User's Guide, MIT mirror에서 확인)

### 1) Restrictions - **이 문제 부류에 대한 결정적 진술**

> "For postbuckling problems involving loss of contact, the Riks method will usually not work;
> inertia or viscous damping forces (such as those provided by dashpots) must be introduced in a
> dynamic or static analysis to stabilize the solution."

주름(wrinkling)은 막의 국소 좌굴 = 접촉 상실이므로 **이 부류에 정확히 해당**한다. 문서가 권하는
방법은 점성 감쇠를 넣은 정적 해석 = `*STABILIZE`. 즉 **현행 정적 + stabilize 경로가 문서 권장**이다.

### 2) Bifurcation - 어떤 문제에서 잘 작동하는가

> "The Riks method works well in snap-through problems - those in which the equilibrium path in
> load-displacement space is smooth and does not branch."

우리는 국소 분기(주름)가 계속 생기는 문제다. 임퍼펙션 주입 자체는 문서가 요구하는 올바른 조치
("the exact postbuckling problem cannot be analyzed directly due to the discontinuous response at
the point of buckling ... must be turned into a problem with continuous response").

### 3) 시간/속도 의존 효과 - 감쇠 금지

> "any effects involving time or strain rate (such as viscous damping or rate-dependent plasticity)
> are no longer treated correctly and should not be used." / "Dashpots should not be used."

Riks에서는 감쇠를 못 쓴다. 그리고 그 감쇠가 1)에서 문서가 권한 도구다 - **구조적 모순**.

### 4) 증분 - 1% 외삽 한계

> "The Riks procedure uses only a 1% extrapolation of the strain increment."

일반 정적 해석보다 외삽이 약하다 = 심한 비선형에서 코렉터가 약하다.
`STATIC, RIKS, DIRECT`(고정 아크길이)는 "not recommended ... prevents Abaqus/Standard from reducing
the arc length when a severe nonlinearity is encountered".

### 5) 종료 조건

> "You can specify a maximum value of the load proportionality factor, lambda_end, or a maximum
> displacement value at a specified degree of freedom." / "The Riks algorithm cannot obtain a
> solution at a given load or displacement value since these are treated as unknowns."

둘 다 없으면 `inc=` 개수까지 계속 간다. **lambda_end 는 '종료 트리거'일 뿐 정확한 해를 보장하지 않는다**
-> 다른 ODB와 대조할 때는 프레임 보간 + 시각 불일치 가드가 필요하다.

### 6) `*STATIC` 데이터줄 기본값 - **여기서 우리가 틀렸다**

| 항목 | 기본값(문서) |
|---|---|
| 1 dl_in | 0/미지정이면 **스텝의 총 아크길이** (= 첫 증분에서 스텝 전체) |
| 2 lperiod | 1.0. 첫 증분의 lambda 증가는 `dlambda_in = dl_in / lperiod` |
| 3 dl_min | 0이면 `min(제안 dl_in, 1e-5 x 총 아크길이)` |
| 4 **dl_max** | **미지정이면 상한 없음(no upper limit is imposed)** |
| 5 lambda_end | lambda 상한 |

**교훈**: dl_max 를 지정하면 자동 증분이 아크길이를 *늘릴* 여지를 잘라먹는다. 기존 실패 런들은
0.1/V0, 0.025/V1 로 조여놨고 **더 좁게 조인 V1 이 더 못 갔다**(lambda 0.0503 -> 0.0129).
문서 기본으로 되돌리려면 4번째 항목을 **비운다**: `0.005, 1.0, 1e-06, , 1`

### 7) `*CONTROLS` 와 line search - 문서에 Riks 언급이 **0회**

`simakey-r-controls.htm` 과 `simaanl-c-convergecontrol.htm` 어디에도 "Riks" 가 없다.
그리고 line search 데이터줄: "Nls, maximum number of line search iterations. **Default Nls=0 for
steps that use the Newton method** and Nls=5 for steps that use the quasi-Newton method."

실측과 일치: Riks 스텝에서 `*Controls, parameters=line search` + `5,` 는 **문법은 통과하나
`ADDITIONAL RESIDUAL/OPERATOR EVALUATIONS FOR LINE SEARCHES = 0`** (V0/V1 두 런 모두) -
**실효 없음**. Riks 는 자기 아크길이 코렉터만 쓴다.

### 8) 도구

`riks_patch_input.py --arc=<dl_in>` (신규): dl_in 직접 지정 + **dl_max 항목을 비워 상한 없음** +
dl_min=1e-6. `--arc-max=<v>` 로만 상한을 건다. `--lpf=` 는 종전 동작 유지(하위 호환).
검증: `--arc=0.005` -> `0.005, 1.0, 1e-06, , 1` / `--arc=0.005 --arc-max=0.05` -> `..., 0.05, 1` /
`--lpf=0.05` -> `0.01, 1.0, 1e-07, 0.025, 0.05` (변화 없음).

## V2 실패의 진짜 원인 - 실험이 아니었다 (2026-10-01)

`--arc=0.005` 로 만든 줄 알았던 `HF_riks_v2` 는 **V0 과 같은 덱**이었다. 증거: 두 `.sta` 의
마지막 6행이 **소수점까지 완전 일치**(0.002126 / -0.0009402 / 1.301e-06 / -0.001231 / -0.001183).
즉 lambda 0.0503 정지가 재현된 것이지, 상한 제거 효과가 없는 것이 아니다.

**원인 2개**
1. 사용자 로컬 사본이 옛 리비전이었다(git pull 누락).
2. **옛 패처가 모르는 `--arc=` 를 조용히 무시**하고 기본값으로 진행했다 - 이게 진짜 결함.

**도구 수정 (재사용 교훈)**
- **모르는 `--` 옵션은 즉시 중단(rc=2)** 한다. 오타 하나로 다른 실험이 조용히 돌아가는 것을 막는다.
- **출력 첫 줄에 빌드 마커**를 찍는다: `[패처] build <날짜><태그>`. 사용자가 새 리비전을 돌렸는지
  로그만 보고 판정할 수 있다(= "도구를 배포했으면 사용자가 새 리비전을 돌렸는지 확인할 장치를
  도구 안에 넣어라").
- 검증은 로컬에서 실제 실행: `--arc=0.005` -> `0.005, 1.0, 1e-06, , 1` /
  오타 `--arcx=0.005` -> rc=2 중단 / `--lpf=0.05` -> 종전과 동일(하위 호환) /
  플래그 없음 -> `0.05, 1.0, 1e-07, 0.1, 1`.
  이 검증 중 `_argv` 이름이 깨진 버그를 잡았다 - 로컬 실행 검증이 아니었으면 배포될 뻔했다.
