# HF 임퍼펙션 시드 연결 파이프라인 — 구현 계획 (2026-10-12)

**목표**: 설계점 (x_c, d_c) 마다 좌굴 시드 모드표를 자동 준비해 HF 포스트버클링에 주입하고,
그 전 과정을 MFBO 루프에서 호출 가능하게 만든다. 사람이 손으로 경로를 넘기지 않는다.

**전제(실측 확인)**
- `seed`(꼭짓점 케이블만) 모델은 **x_c 무관** — `HAS_CLAMPS=False` 라 클램프 패치/RP/케이블이
  만들어지지 않고, 경계 제외집합도 빈 집합이다(`run_abaqus_buckle.py:906-907`). → 한 번 만들어 재사용.
- `clamp_only`(클램프 케이블만) 모델은 **x_c 의존** → 설계점마다 재추출(런 약 40초).
- 주입 기구 자체는 이미 검증됨(§19.55: 기하 섭동 10201/10201 노드, max|dz|=9.16e-7 m).
- 사용자 실측: x_c 를 바꿔도 좌굴 해석은 문제없이 수렴 → 창(10/25) 문제는 없음.
  남은 미검증 가설은 "시드 패턴이 클램프 위치와 어긋나면 후좌굴 결과가 편향되는가" (→ §8 A/B).

---

## 1. 현황 — 되는 것 / 빠진 것

| 구성 | 상태 | 근거 |
|---|---|---|
| 모드 추출 (2케이스) | ✅ | `run_abaqus_buckle.py` `_CASE_TABLE = {seed, clamp_only}` |
| 모드표 생성 | ✅ | `aba_mode_from_odb.py <odb> Step-Buckle <out> [n]` |
| 합성 (2+2) | ✅ | `combine_mode_tables.py <out> <w>:<tbl>:<m1,m2> …` |
| HF 기하 섭동 주입 | ✅ | `run_abaqus.py:527` `IMPERFECTION_MODE='odb_table'` + `aba_imperfection.build_perturbation` |
| **설계점별 표 전달 통로** | ❌ | `run_abaqus.py:528` `MODE_TABLE = os.path.join('..','modes_2p2.txt')` — **하드코딩**. env 는 `SOLARSAIL_MODE_ODB`(odb_direct 전용)만 있음 |
| **지문 검증이 현행 소스와 불일치** | ❌ | `run_abaqus.py:620-645` 이 `run_abaqus_mode.py`(§19.64 로 비교 은퇴된 스크립트)와 지문을 비교하고, **차단하지 않고 기록만** 함 |
| **시드 출처 로그** | ⚠️ | 로그가 `source=? / step=? / instance=?` 로 찍힘 — 합성표 헤더가 `mode_table_info` 형식과 어긋남 |
| MFBO → HF 호출 | ⚠️ | `mfbo.py:79 get_abaqus()` 가 `-- {HF\|LF} {x1} {x2}` 로 호출, **"env 주입 없음"** 주석대로 시드 통로가 없음 |

---

## 2. 설계 결정 (권고안)

**D1. 시드 구성** — `seed`(x_c 무관, 고정) + `clamp_only`(x_c별) 2+2 합성 유지. (현행 확정 설계)
- 대안: `seed` 4모드만(가장 단순, 클램프 국소 기여 없음) / `clamp_only` 4모드만(전역 모드 없음).
- **권고: 현행 2+2 유지.** `IMPERFECTION_MAX_MODES=4` 와도 맞는다.

**D2. 표 전달 통로** — **결정적 경로 규약 + env 오버라이드** 이중.
- 기본: `<repo>/code/seeds/modes_xc<NNN>_dc<NNN>.txt` (스크립트가 x_c·d_c 로 조립, 없으면 **즉시 중단**).
- 실험용: `SOLARSAIL_MODE_TABLE` env 로 덮기(기존 `SOLARSAIL_MODE_ODB` 와 같은 방식).
- 대안: CLI 인자 추가(✗ 위치기반 파서 `argv[-3:]` 와 충돌) / HF 안에서 자동 합성(✗ HF 안에서 Abaqus 재귀 호출).

**D3. 지문(fingerprint) 정책** — 시드 소스 = buckle 스크립트 + (x_c, d_c) + 입력표 지문.
- 합성표 헤더에 `fingerprint=` 를 쓰고, HF 가 **불일치면 중단**(현행은 경고만 ✗ → 조용한 오주입 차단).
- `seed`/`clamp` 반쪽 지문을 각각 남기려면 헤더에 `source_seed=` / `source_clamp=` 두 줄을 둔다.

**D4. 노드 수 게이트** — 현행 `!!! WARNING: 모드표 노드 N개 중 M개만 적용`(WARN)을 **중단**으로 승격.
  근거: 메쉬/라벨 드리프트는 조용한 오주입이고, 게이트가 없으면 사용자가 못 본다.

**D5. d_c 반영 여부** (⚠️ 결정 필요) —
- (a) 시드의 클램프 당김을 HF 의 d_c 와 일치시킨다: 정확, 설계점마다 buckle 1런 추가.
- (b) 시드를 `d_c=CLAMP_DC=0.0556` 로 고정: 저렴, "씨드는 guidance" 논리에 부합.
- **권고 (b)** 로 시작하고 §8 A/B 로 (a) 필요성을 판정한다.

---

## 3. `run_abaqus.py` 변경 명세 (5항목)

### 3-1. 시드 표 경로 결정 (신규, `:528` 대체)
```python
#   [2026-10-12] 설계점별 시드표. 기본은 결정적 경로, 실험은 env 로 덮는다.
_SEED_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'seeds')
_TBL_ENV = os.environ.get('SOLARSAIL_MODE_TABLE')
if _TBL_ENV:
    MODE_TABLE = os.path.abspath(_TBL_ENV)
else:
    MODE_TABLE = os.path.join(_SEED_DIR,
                              'modes_xc%03d_dc%03d.txt' % (int(round(x_c * 100.0)),
                                                           int(round(d_c * 100.0))))
```
- `x_c`/`d_c` 는 `:203-205` 에서 이미 파싱된 값 → **이 블록은 그 뒤**(모드 로딩 전)에 둔다.
- 없으면 기존 `load_mode_table` 이 `ImperfectionSourceError` → `sys.exit(1)`(현행 동작 유지).
- **검증**: `python3 -c` 로 모듈 상수만 뽑는 단위테스트는 불가(Abaqus import) → 경로 조립만
  순수 함수 `seed_table_path(x_c, d_c, env=None)` 로 빼고 라이선스 0 테스트(§7 Task 1).

### 3-2. 시드 provenance 로그 (수정)
- 현행 `:616-619` 의 `source=?` 를 고치려면 합성표 헤더 형식을 `mode_table_info` 규약에 맞춘다(§4).
- HF 쪽에서 한 줄 추가:
  `[SEED] x_c=%.4g d_c=%.4g table=%s modes=%s source_seed=%s source_clamp=%s`

### 3-3. 지문 게이트 교체 (수정 `:620-645`)
- 비교 대상: `run_abaqus_mode.py` → **시드 산출물의 지문**(합성표 헤더 + buckle 스크립트 상수 해시).
- 정책: 불일치 → `!!!` + `sys.exit(1)`(제출 전 중단). `mode_table_check.txt` 기록은 유지.

### 3-4. 노드 수 게이트 승격 (수정 `_apply_imperfection_to_part`)
- `len(_apply) != len(_labels)` → `raise RuntimeError`(현행 WARN).

### 3-5. (선택) `IMPERFECTION_MODES` 기본값
- 현행 `(1,2,3,4)` + 표 자동 조정 유지. 합성표에 4모드가 있으면 그대로 4.

---

## 4. 신규 스크립트 `code/make_seed.py`

```
python3 make_seed.py --x_c 0.5 --d_c 0.0556 [--force] [--dry-run] [--n-clamp 2]
```
| 단계 | 명령 | 산출물 |
|---|---|---|
| 1 | (캐시 확인) `seeds/modes_seed.txt` 없거나 `--force` | `seed` 계열 고정표(x_c 무관) |
| 1a | `abaqus cae noGUI=run_abaqus_buckle.py -- 0.5 1.8e-5 0` (CASE=`seed`) | `buckle/Buckle_…_seed_….odb` |
| 1b | `abaqus python aba_mode_from_odb.py <odb> Step-Buckle seeds/modes_seed.txt 4` | 모드 4개 |
| 2a | `abaqus cae noGUI=run_abaqus_buckle.py -- <x_c> 1.8e-5 <clamp_pull>` (CASE=`clamp_only`) | `buckle/Buckle_…_clamp_only_….odb` |
| 2b | `abaqus python aba_mode_from_odb.py <odb> Step-Buckle seeds/modes_clamp_xc<NNN>.txt 2` | 모드 2개 |
| 3 | `abaqus python combine_mode_tables.py seeds/modes_xc<NNN>_dc<NNN>.txt 1.0:seeds/modes_seed.txt:1,2 1.0:seeds/modes_clamp_xc<NNN>.txt:1,2` | **HF 가 읽는 표** |
| 4 | 검증: 노드수 == 10201 · 모드 == 4 · 모드별 peak\|u3\| == 1 · 헤더에 지문/x_c/d_c | `seeds/<table>.meta.txt` |

- `clamp_pull = DISP_GLOBAL(1.8e-5) * d_c` (§D5).
- 멱등: 표 + meta 가 있고 지문이 맞으면 **건너뜀**. `--force` 로 재생성.
- `--dry-run`: 실행할 `abaqus …` 명령만 출력(라이선스 0 검증용).
- `CASE` 는 상수라 스크립트가 소스를 수정하지 않는다 → **`--case` 인자를 신설**하거나
  두 케이스를 연속 실행하기 위해 `run_abaqus_buckle.py` 에 `SOLARSAIL_BUCKLE_CASE` env 오버라이드를
  추가한다(권고: env, CLI 파서를 건드리지 않음).

---

## 5. MFBO 연결 (`code/mfbo.py` 수정)

- `get_abaqus()` HF 분기에서 HF 실행 **직전**:
```python
subprocess.run(f'python3 "{os.path.join(_HERE,"make_seed.py")}" '
               f'--x_c {x1} --d_c {x2}', shell=True, check=False,
               capture_output=True, text=True, cwd=_HERE)
```
  → 실패 시 `return FAIL, FAIL`(기존 LF 선행 실패 처리와 동일 패턴).
- HF 명령은 **그대로**(`-- HF {x1} {x2}`) — 표 경로는 `run_abaqus.py` 가 결정적으로 조립한다.
- 캐시 덕분에 같은 x_c 재방문 시 재추출 비용 0.
- (선택) 초기 스윕 격자를 **사전 배치 생성**(`make_seed.py --batch`)해 루프 중 라이선스 대기를 제거.

---

## 6. 작업 순서 (각 task = 커밋 1개, 검증 명령 포함)

**Task 1 — `seed_table_path()` 순수 함수 + 테스트**
- Modify: `code/run_abaqus.py` (경로 조립만 함수로 추출, `:528` 대체)
- Test: `code/test_seed_path.py` (신규, Abaqus 불필요)
- Run: `python3 test_seed_path.py` → 기대: 4케이스 PASS (env 우선 / 기본 경로 / 자리수 / 대소문자)

**Task 2 — `SOLARSAIL_MODE_TABLE` env + 결정적 경로 적용**
- Run: `python3 -m py_compile run_abaqus.py` + `python3 static_leak_check.py` + `python3 check_model_consistency.py`(PASS)

**Task 3 — 노드 수 게이트 승격(WARN→중단)**
- Verify: 하네스 `python3 buckle_stub_harness.py seed 0.5` EXIT=0 (HF 하네스는 없으므로 정적 확인)

**Task 4 — `make_seed.py` 골격 + `--dry-run`**
- Create: `code/make_seed.py`
- Run: `python3 make_seed.py --x_c 0.5 --d_c 0.0556 --dry-run` → 기대: 명령 4줄, 실행 없음

**Task 5 — 실제 생성 1회 + 검증 블록**
- Run(사용자, Abaqus): `python3 make_seed.py --x_c 0.5 --d_c 0.0556 --force`
- 기대: `seeds/modes_xc050_dc006.txt` 생성, meta 노드 10201 / 모드 4

**Task 6 — 지문 헤더 규약 정합 (`combine_mode_tables.py` + `aba_imperfection.mode_table_info`)**
- Verify: HF `RUN_MODE='write'` 로그에서 `source_seed=`/`source_clamp=` 가 `?` 가 아닌 값으로 찍힘

**Task 7 — 지문 게이트 교체(HF 측)**
- Verify: 표를 손으로 옛 지문으로 바꾸면 HF 가 **제출 전 중단**(로그로 확인)

**Task 8 — `mfbo.py` 연결**
- Verify: `MFBO_MOCK=1` 로 배관 테스트(기존 목 오라클) → `make_seed` 호출 경로만 확인

---

## 7. 검증 계획

| 단계 | 방법 | 라이선스 |
|---|---|---|
| V1 | `py_compile` / `static_leak_check` / `check_model_consistency` / 하네스 | 0 |
| V2 | `make_seed.py --dry-run` 명령 생성 확인 | 0 |
| V3 | 경로·지문 순수 함수 단위테스트 | 0 |
| V4 | HF `RUN_MODE='write'` → `.inp` 노드 좌표 확인(max\|dz\| ≈ `THICKNESS*IMPERFECTION_AMPL_T`) | cae (제출 없음) |
| V5 | x_c=0.5 한 점 E2E: LF → seed → HF → `eval_abaqus.py` | standard (~40분) |
| V6 | **시드 감도 A/B**(같은 x_c, 다른 x_c 의 clamp 시드) → 목적함수 비교 | standard ×2 |

V6 가 §D5/미검증 가설을 판정한다: 두 결과가 오차범위 내면 고정 시드로 충분(파이프라인 단순화),
유의하게 다르면 x_c별 재추출이 필수(현 계획 유지).

---

## 8. 리스크 / 열린 질문

1. **d_c 의존성(§D5)** — 미검증. V6 로 판정.
2. **지문 정의** — 합성표는 소스가 2개다. 헤더에 2줄로 나눠 남기는 안(§D3)이 최소 침습.
3. **기존 실행 호환** — 현재 HF 는 `code/modes_2p2.txt`(x_c=0.5 산출물)를 읽는다.
   Task 5 에서 이 파일을 `seeds/modes_xc050_dc006.txt` 로 옮기면 기존 런과 동일해진다.
4. **라이선스** — buckle 준비는 cae+standard 를 쓰고 설계점당 ~1분. 사전 배치 생성을 권고.
5. **`CASE` 상수** — 스크립트가 자기 소스를 고치지 않게 env 오버라이드(`SOLARSAIL_BUCKLE_CASE`)를
   쓸지, `make_seed.py` 가 소스를 패치할지 결정 필요(권고: env).
6. **`seed` 표의 재현성** — x_c 무관이지만 메쉬/상수 변경 시 지문으로 무효화된다(§D3).

---

## 9. 결정이 필요한 것 (사용자)

- **Q1** §D5: d_c 를 시드에 반영할지 — (a) 반영(정확·비쌈) / **(b) 고정(권고)**
- **Q2** §D3: 지문 불일치 시 (a) 중단 / (b) 경고만 — **(a) 중단 권고**
- **Q3** §5: `CASE` 전환 방식 — **(a) env `SOLARSAIL_BUCKLE_CASE`(권고)** / (b) 소스 패치
- **Q4** V6 A/B 를 지금 돌릴지(standard 2런 ≈ 80분) 아니면 V5 이후로 미룰지
