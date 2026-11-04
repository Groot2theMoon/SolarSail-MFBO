# attic/ — 은퇴한 구현 (동결)

P1(2026-11-03)에서 `code/` 밖으로 내린 것 중 **더 이상 맞지 않는** 파일만 여기 있다.
`tools/` 와 다르다: tools 는 지금도 돌아가고 의미가 있는 장비, attic 은 **대체되었거나 깨진** 구현이다.

## 원칙

- **고치지 않는다.** 여기 파일은 현재 `code/` 계약(상수 이름·함수 시그니처)을 따라가지 않는다.
  `sys.path` 보정도 하지 않았다. 그래서 import 하면 실패하는 것이 정상이다.
- 되살리려면 **증거를 먼저 만들 것**: 현재 `code/` 에 맞게 고치고, 무엇을 검증했는지 남긴 뒤 `git mv` 로 올린다.
- 여기서 배우려는 사람은 `git log --follow -- <파일>` 로 이력부터 볼 것.

## 왜 여기 있나 (2026-11-03 실측)

| 파일 | 이유 |
|:---|:---|
| `torque_stats.py` + `torque_stats_stub.py` | 토크 손실 지표 폐기. `code/slope_stats.py` 가 대체(면적가중 기울기 통계) |
| `eval_abaqus_avg.py` + `eval_abaqus_avg_stub.py` | 평균방식 평가 폐기. `code/eval_abaqus.py` 가 대체 |
| `check_odb_fields.py` | `tools/compare_odb_fields.py` 가 대체(필드 단위 노드 비교) |
| `eval_lf_stub.py` | tag `pre-refactor-2026-11-03` 에서 이미 **FAIL** (LF 계약 변경 미반영) |
| `slope_stats_stub.py` | 같은 tag 에서 **106/110 FAIL** (추출 배선 금지문자 검사 등) |
| `pull_top_stub.py` | `PULL_TOP` 경로 폐기 → 배선 검사 **FAIL** |
| `seed_family_stub.py` | `NameError: SEED_FAMILY_NAME` — 그 상수는 `SEED_FAMILY_TABLE` 로 바뀌었다 |
| `run_abaqus_mode.py` | 모드 소스 비교(run_abaqus.py ↔ 이 스크립트)가 §19.64 로 은퇴. HF 임퍼펙션은 `run_abaqus_buckle.py` → `aba_mode_from_odb.py` 경로를 쓴다 |
| `RA_calc.py` | 1회성 광학 계산. 입력 `nk_data.csv` / `sun_data.csv` 가 저장소에 없다 |

"지금도 통과한다"는 이유로 `tools/` 로 올린 것: `buckle_stub_harness.py`(`seed 0.5` exit 0),
`buckle_compare_stub.py`(24 checks PASS), `buckle_compare.py`, `compare_states.py`.
