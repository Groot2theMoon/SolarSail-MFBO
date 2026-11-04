# baseline/ — 리팩토링 회귀 기준선

`code/refactor_gate.py` 가 읽고 쓰는 디렉터리. 리팩토링 **전** 상태를 여기에 얼려 두고,
각 단계 후에 같은 명령으로 대조한다.

| 게이트 | 대상 | 판정 |
|:---|:---|:---|
| G1 deck | `aba/*.inp` (job name·날짜·절대경로 제거 후) | sha256 동일 |
| G2 results | `results/results_*.txt` (`key=value`) | 전 키 값 동일 (오차 0, `timestamp` 제외) |
| G3 probe | `hf_run_log.txt` 의 `key=value` 쌍 | 다중집합 동일 |

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

### 노이즈 바닥 (캡처 직후, 같은 두 런을 한 번 더)

```bat
abaqus cae noGUI=run_abaqus.py -- LF 0.2 0.0001
abaqus cae noGUI=run_abaqus.py -- HF 0.2 0.0001
python refactor_gate.py check --label p0 ^
    --deck "aba/*.inp" --results "results/results_*.txt" --probe "hf_run_log.txt"
```

두 번째 런이 PASS 면 기계 노이즈가 0 이라는 뜻이고, FAIL 이면 그 diff 가 **허용오차의 근거**다
(G2 는 오차 0, G3 는 다중집합 동일인데 로그에 `λ=-1.586670e-03` 같은 고유값이 들어간다).
이 단계를 건너뛰면 P2 이후 모든 단계에서 원인 불명 FAIL 을 만난다. 허용오차를 넣을 때는
제외 키와 근거를 이 문서에 적는다.

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
- `source_manifest.txt` 는 `snapshot --label p0` 산출물이다(파일별 sha256 + 줄 수).
  코드 변경 자체는 게이트 대상이 아니다 — 덱/결과만 본다.
