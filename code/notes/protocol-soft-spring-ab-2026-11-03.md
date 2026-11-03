# soft-spring A/B 프로토콜 — 3자유변 u3=0 이 δn 채널을 닫았는가

작성 2026-11-03. 대상 커밋: `a008425` 이후.

## 목적

우리 HF는 Galhofo 2022의 **좌굴 모델**(§3.3, Fig. 12)에서 `u3=0`(세 자유변)을 상속했다. 경계 전체가 면외 고정이면 발산정리에 의해

    ⟨∇u⟩ = (1/A)∮ u n ds = 0   ⇒   δn ≡ 0

이 되어 **자세 오차 채널이 구조적으로 닫힌다**. 이 A/B는 그것이 "상속한 가정의 결과"인지 "이 형상의 물리"인지 가른다.

## 사전 조건

| 항목 | 값 | 이유 |
|:---|:---|:---|
| config | `configs/seed.json` (`DISP_GLOBAL` 없음) | HF 운용점 165 µm 유지 |
| 운용점 | 165 µm | 코드 주석의 선형 외삽: 2122 Pa @ 50 µm → **7.00 kPa @ 165 µm**. 스프링 기준 σt = 7 kPa × 5 µm 와 일치. **외삽이며 측정 아님**(R-13) |
| seed 표 | 기존 `seeds/modes_xc<xc>.txt` 재사용 (`make_seed: false`) | 모드 기저를 고정해 변수를 하나로 |
| x_c, d_c | A와 동일 | — |

## Arm A — 기준 (완료)

`RELEASE_EDGE_U3 = False` (양 스크립트). 2026-11-03 런이 이것이다.

기록할 값: `hf`, `dn_deg`, `u_z,ave`, LF 지표.

## Arm B — 스프링 해제

**`run_abaqus.py:190`만 `True`로 바꾼다.** `run_abaqus_buckle.py`는 건드리지 않는다 — 모드 소스까지 바꾸면 임퍼펙션 기저와 BC가 동시에 변해 교란된다.

```python
RELEASE_EDGE_U3 = True
EDGE_SPRING_ALPHA = 1.0e-3   # 기본값. k = 3.5e-05 N/m = 1e-3 × (7 kPa × 5 µm)
```

```powershell
abaqus cae noGUI=run_abaqus.py -- HF 0.2 0.0001 --config configs/seed.json
```

## 기록 (양 arm 공통)

- 로그의 `[BC]` 줄 — 이제 **양 분기 모두** 출력한다
  - A: `[BC] RELEASE_EDGE_U3=False -> BC_Edges_Only_Z 유지 (세 자유변 u3=0, n=..., step=...)`
  - B: `[BC] RELEASE_EDGE_U3=True -> BC_Edges_Only_Z 제거, EdgeZ_Spring k=...`
- `hf`, `dn_deg`, `u_z,ave`

## 판정

| `dn_deg` (B) | 해석 | 조치 |
|:---|:---|:---|
| O(1°) | BC가 자세 축을 닫고 있었다 → **상속 가정의 artifact** | §3.5.6의 "δn ≡ 0" 서술을 조건부로 수정 |
| ≤ O(1e−4°) | 이 형상의 **물리** | δn=0을 결과로 확정, 문헌 비교 시 명시 |
| 그 사이 | 애매 | 모드별 기여 분해 또는 α 스윕(1e-2, 1e-1)로 경향 확인 |

`hf`는 A/B에서 크게 변하지 않아야 정상이다(BC는 경계 근처 국부 주름에 주로 작용). 크게 변하면 경계 artifact가 지배적이었다는 뜻이므로 그것 자체가 결과다.

## 주의

- 스프링은 **model-level** engineeringFeature다 → LF 스텝에도 작용한다. A와 LF 지표를 비교하려면 B에서 LF도 재측정해야 한다.
- 스프링 강성은 **7 kPa 기준 절대값**이다. 운용점을 165 µm에서 옮기면 α를 그 운용점 장력으로 재계산해야 한다.
- B에서 수렴 실패나 음수 고유값이 다수 나오면 α를 1e-2로 올려 1회 재시도한다. 그래도 실패면 "스프링으로는 해제 불가"로 기록하고 soft-spring 경로를 접는다.

## 복원

```python
RELEASE_EDGE_U3 = False
```
로 되돌리고 `git diff`가 clean인지 확인한다. 다음 런 로그에 `[BC] ... False`가 다시 찍히면 복원 완료다.
