# baseline/ — 리팩토링 회귀 기준선

`code/refactor_gate.py` 가 읽고 쓰는 디렉터리. 리팩토링 **전** 상태를 여기에 얼려 두고,
각 단계 후에 같은 명령으로 대조한다.

| 게이트 | 대상 | 판정 |
|:---|:---|:---|
| G1 deck | `aba/*.inp` (job name·날짜·절대경로 제거 후) | sha256 동일 |
| G2 results | `results/results_*.txt` (`key=value`) | 전 키 값 동일 (오차 0, `timestamp` 제외) |
| G3 probe | `hf_run_log.txt` 의 `key=value` 쌍 | 다중집합 동일 |

## P0 캡처 절차 (Abaqus 머신에서 1회)

`code/` 에서 기본 설정(x_c=0.5, d_c=0.5) 2런. `--config` 를 주지 않는다 =
`SEED_DIV=200` 기본 격자(약 1.8만 요소). config 를 주면 설계점이 바뀌어 게이트가 다른 것을 잰다.

```bat
cd code
abaqus cae noGUI=run_abaqus.py -- LF 0.5 0.5
abaqus cae noGUI=run_abaqus.py -- HF 0.5 0.5

python refactor_gate.py capture --label p0 ^
    --deck "aba/*.inp" --results "results/results_*.txt" --probe "hf_run_log.txt"
```

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
