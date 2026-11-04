# baseline/ — 리팩토링 회귀 기준선

`code/refactor_gate.py` 가 읽고 쓰는 디렉터리. 리팩토링 **전** 상태를 여기에 얼려 두고,
각 단계 후에 같은 명령으로 대조한다.

| 게이트 | 대상 | 판정 |
|:---|:---|:---|
| G1 deck | `aba/*.inp` (job name·날짜·절대경로 제거 후) | sha256 동일 |
| G2 results | `results/results_*.txt` (`key=value`) | 상대 1e-5 이내, `timestamp`/`r13_ct_odb_mtime` 제외 |
| G3 probe | `hf_run_log.txt` 의 `key=value` 쌍 | 다중집합 동일 |

## 기준선 목록

| label | 운용점 | 용도 |
|:---|:---|:---|
| `p0` | `CABLE_PRESTRESS=2.69e7` | 2026-11-04 이전 리팩토링 회귀 기준. **역사적 기록** |
| `p5` | `CABLE_PRESTRESS=1.71e7` (중심 von Mises 1 psi) | 2026-11-04 이후 **모든 비교의 기준** |

운용점은 §19.122 의 5점 LF 스윕(14–23 MPa, 단조)으로 확정했다. 이후 `check` 는 `--label p5` 를 쓴다.
`p0` 대조는 운용점이 달라 실패하는 것이 정상이다.

## P0 캡처 절차 (Abaqus 머신에서 1회)

운용점은 **x_c=0.2, d_c=0.0001** 이다(2026-11-03 실런과 같은 점). `--config` 를 주지 않는다 =
코드 기본값(`SEED_DIV=200`, `MESH_ORDER=1`, `N_EIG_BUCKLE=4`, `DISP_GLOBAL=165e-6`)으로 돌아
그 실런과 같은 조건이 된다.

```bat
cd code
ren aba aba_pre_p0          & rem G1 은 aba/*.inp 전부를 기준선에 넣는다 - 과거 덱을 먼저 치운다
abaqus cae noGUI=run_abaqus.py -- LF 0.2 0.0001
abaqus cae noGUI=run_abaqus.py -- HF 0.2 0.0001

python refactor_gate.py capture --label p0 ^
    --deck "aba/*.inp" --results "results/results_*.txt" --probe "hf_run_log.txt"
```

### 왜 0.5 0.5 가 아닌가

- HF 는 제출 **전에** 시드표를 검사한다(`run_abaqus.py:1564`). 표가 없으면
  `ImperfectionSourceError` → `... aborting BEFORE HF job submit`. x_c=0.2 의 합성표
  (`seeds/modes_xc0p2.txt`)는 2026-11-03 런으로 존재가 확인됐지만, `0.5/0.5` 표는 저장소에 없다.
- `0.5/0.5` 는 어느 설정의 점도 아니다(`paper_sd632.json`=0.5/1.5, `seed.json`=0.2/1e-4).

### `--config configs/seed.json` 을 쓰면 안 된다

`seed.json` 은 `"make_seed": true, "force_seed": true` 다 — 런마다 좌굴 모드를 **재생성**한다.
재생성된 모드의 부호 기준 노드가 달라지면 임퍼펙션 장이 뒤집혀 `hf` 가 바뀐다. 코드가 한 줄도
안 변해도 게이트가 FAIL 한다. `--config` 없이 도는 것이 이 게이트의 전제다.

### 측정 사이클은 반드시 초기화한다 (2026-11-04 실측)

- LF 는 sidecar(`LF_*.odb.model.txt`)의 tag 가 같으면 **ODB 를 재사용하고 잡을 재제출하지 않는다**
  (실측: 두 번째 LF 런이 7.8초에 끝남). 그러면 덱도 다시 쓰지 않아 G1 이 낡은 파일을 비교한다.
- LF 런도 포인트 아카이브 꼬리에서 `results_*.txt` 에 key 를 **append** 한다. 파일을 지우지 않으면
  같은 키가 두 번 들어가 G2 가 중복으로 FAIL 한다.
- 그래서 각 측정 사이클 전에 LF 산출물과 결과를 지운다:

```bat
del /q aba\LF_Analysis_*.odb aba\LF_Analysis_*.odb.model.txt results\results_*.txt
```

(HF 잡은 매 런 제출되고, 제출 전에 자기 stale 산출물을 스스로 지운다.)

### 노이즈 바닥 (캡처 직후, 같은 두 런을 한 번 더)

```bat
abaqus cae noGUI=run_abaqus.py -- LF 0.2 0.0001
abaqus cae noGUI=run_abaqus.py -- HF 0.2 0.0001
python refactor_gate.py check --label p0 ^
    --deck "aba/*.inp" --results "results/results_*.txt" --probe "hf_run_log.txt"
```

두 번째 런이 PASS 면 기계 노이즈가 0 이라는 뜻이고, FAIL 이면 그 diff 가 **허용오차의 근거**다.
이 단계를 건너뛰면 P2 이후 모든 단계에서 원인 불명 FAIL 을 만난다.

**측정 결과 (2026-11-04, labpc, 위 절차 그대로 재해석)**

| 게이트 | 결과 |
|:---|:---|
| G1 deck | **3/3 바이트 동일** (LF/HF/grid) |
| G2 results | 51키 불일치 — 최대 상대차 **1.44e-6**(`slope_f2_p95`), `hf` 는 1.5e-9 |
| G3 probe | 12kv 불일치 — 같은 잡음(λ/f/zk) |

절대오차 0 은 이 기계에서 재현 불가능하다. 그래서 `REL_TOL = 1e-5`(실측 최대의 약 7배)를
G2/G3 에 적용하고 `r13_ct_odb_mtime`(파일 수정시각)을 제외했다. 같은 diff 를 다시 돌리면 **PASS**.
G1 은 여전히 오차 0 이며, 이 숫자 허용오차는 덱 비교에 적용되지 않는다.

## 이후 단계마다

```bat
python refactor_gate.py check --label p0 ^
    --deck "aba/*.inp" --results "results/results_*.txt" --probe "hf_run_log.txt"
```

종료코드 0 = PASS, 1 = FAIL. FAIL 이면 게이트가 어느 파일 몇 번째 줄이 어떻게 달라졌는지
출력한다(G1 은 정규화 덱 unified diff, G2/G3 는 `키 kind baseline=.. now=..`).

## 주의

- G1 은 **숫자를 반올림하지 않는다.** 부동소수 표기 차이는 FAIL 로 뜨고 사람이 diff 로 본다.
  반올림 허용은 "조용한 모델 변경"을 통과시키므로 게이트가 아니다.
- `hf_run_log.txt` 는 런마다 truncate 된다(`emit_start`). 캡처는 **HF 런 직후**에 한다.
- `--label` 은 단계별로 나눈다: `p0`, `p2`, `p4` … 비교는 항상 `p0` 기준.
- G2/G3 의 **절대경로 값**(`r13_ct_odb`, `table`, `cwd`)은 **basename 으로 비교**한다. baseline 은
  임시 사본에서, check 는 정본 체크아웃에서 돌 수 있기 때문이다(G1 이 덱 경로를 정규화하는 것과 같은 이유).
- `code/` 루트의 `Buckle_Analysis.*`(자기 좌굴 잔재)가 있으면 N-6 λ 로그 줄이 추가되어 G3 다중집합이
  달라진다. 측정 전에 `stale_pre_p2/` 같은 곳으로 옮긴다(삭제하지 않는다).
- `source_manifest.txt` 는 `snapshot --label p0` 산출물이다(파일별 sha256 + 줄 수).
  해시는 **CRLF -> LF 정규화 후** 계산한다 - Windows/Linux 체크아웃이 같은 값을 낸다(§19.117).
  코드 변경 자체는 게이트 대상이 아니다 — 덱/결과만 본다.
