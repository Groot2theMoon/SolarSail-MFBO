# tools/ — 수동 실행 장비 (런타임 경로 밖)

여기 있는 것은 `run_analysis.py` 가 부르지 않는다. 사람이 직접 돌려서 **측정·검증**하는 스크립트다.
P1(2026-11-03)에서 `code/` 밖으로 옮겼다. 런타임 경로는 `code/*.py` 20개뿐이다.

## 경로 규칙

옮기면서 각 스크립트의 `HERE` 를 **부모(`code/`)** 로 바꿨다. 이유: 이들이 읽는 것은
`run_abaqus.py`·`aba_imperfection.py`·`aba/`·`results/` 처럼 전부 `code/` 아래에 있다.
그래서 출력 경로(`results/`, `sweep_out/`, `buckle/`)도 **이동 전과 같은 곳**에 쓰인다.

부모 모듈을 import 하는 스크립트는 `_CODE` 부트스트랩을 갖는다:

```python
_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _CODE)
```

## 목록

| 파일 | 무엇을 하나 | Abaqus 필요 |
|:---|:---|:---:|
| `check_model_consistency.py` | HF ↔ 좌굴 ↔ 모드소스 세 스크립트의 모델 정의 일치 검사. **제출 전 게이트** | 아니오 |
| `static_leak_check.py` | 실행 없이 `NameError` 부류(빠진 import / 복사가 끌고 온 이름 / 절점집합을 요소자리에) 정적 검사 | 아니오 |
| `buckle_stub_harness.py` | Abaqus 없이 `run_abaqus_buckle.py` 의 `build_model()` 실행 → BC/CLOAD 배선 추적 | 아니오 |
| `buckle_compare_stub.py` | `buckle_compare.py` 를 가짜 입력으로 실행하는 하네스 | 아니오 |
| `test_seed_path.py` / `test_fingerprint_api.py` | 설계점 시드표 경로·지문 순수함수 단위시험 | 아니오 |
| `ab_softspring.py` | `RELEASE_EDGE_U3` A/B 실행기 (상수 1줄 편집 → 런 → 원본 바이트 복원 → sha1 검증) | 예 |
| `probe_grid_import.py` | `.inp` → `PartFromInputFile` 경로가 이 기계에서 되는지 30초 검증 | 예(CAE) |
| `hf_by_state.py` | hf 를 LF 상태(taut/wrinkled/slack) 별로 분해 | 예 |
| `slope_stats.py`(code) 짝 | — | — |
| `zk_map_extract.py` / `zk_map_check.py` / `zk_map_plot.py` | 요소별 ZK/I_w/손실 지도 추출·검증·도시 | 추출만 예 |
| `compare_odb_fields.py` / `compare_states.py` | 두 ODB 의 변위장·상태 비교 | 예 |
| `buckle_compare.py` / `buckle_mirror_check.py` / `mode_symmetry_probe.py` | 좌굴 런 비교 / 미러 대칭 / 모드 좌우대칭 진단 | 예 |
| `energy_check.py` / `sta_increment_report.py` / `probe_partial_hf.py` | 에너지 이력 / `.sta` 완주 판정 / 미완주 잡 상태 판독 | 예 |
| `rf_history.py` / `node_probe.py` / `mesh_quality_probe.py` | 반력 이력 / 국소 노드 진단 / 왜곡 요소 위치 | 예 |
| `probe_prestress_sweep.py` | 케이블 프리스트레스 스윕 (중앙 von Mises 운용점 탐색) | 예 |
| `sweep_buckle_xc.py` / `sweep_clamp_preset.py` | 좌굴모드가 성공하는 `x_c` 범위 / 클램프 프리셋 OFAT | 예 |
| `screen_report.py` / `inp_mesh_compare.py` / `wrinkle_spectrum.py` | `.msg`/`.sta` 판정 / 두 `.inp` 메쉬 정합 / 주름 파장 측정 | 아니오(파일) |
| `run_abaqus_cable.py` | 클램프 없는 대조 모델(좌굴 오라클). **모드 소스는 `run_abaqus_mode.py`(code/)** | 예 |
| `run_hf_capture.py` | `run_abaqus.py` 를 감싸 출력이 삼켜져도 흔적을 남긴다 | 예 |

## 검증 (2026-11-03, P1 직후)

```
python tools/static_leak_check.py selftest   # OK
python tools/buckle_compare_stub.py          # PASS 24 checks
python tools/buckle_stub_harness.py seed 0.5 # exit 0
python tools/mode_symmetry_probe.py selftest # 12/12
python tools/test_fingerprint_api.py         # 10/10
python tools/test_seed_path.py               # 47/49 (2건은 이동 전 tag 에서도 동일 실패)
```
