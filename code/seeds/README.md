# HF 임퍼펙션 시드표 (`code/seeds/`)

HF 포스트버클링(`run_abaqus.py`)이 읽는 좌굴모드 임퍼펙션 표가 여기 쌓인다.
경로 규약과 게이트는 `aba_imperfection.py`(stdlib) + `make_seed.py` 가 정본이다.

## 파일

| 파일 | 의미 | x_c 의존 |
|---|---|---|
| `modes_seed.txt` | `seed` 가족(꼭짓점 케이블만) 좌굴 4모드 | **아니오** — 한 번 만들면 재사용 |
| `modes_clamp_xc<NNN>_dc<NNN>.txt` | `clamp_only` 가족(클램프 케이블만) 모드 | 예 |
| `modes_xc<NNN>.txt` | **HF 가 읽는 합성표**(2+2) + 설계점/지문 헤더 | 예 |
| `modes_xc<NNN>.meta.txt` | 생성 시각·노드수·모드수·peak u3·게이트 결과 | 예 |

`<NNN>` 은 `value_tag` 표기다(0.5 -> `0p5`, 0.0556 -> `0p0556`). HF 는 자기 `x_c`·`d_c` 로
`modes_xc<NNN>.txt` 를 조립해 읽고, **없거나 지문이 어긋나면 제출 전에 중단**한다.

## 커밋된 표 (2026-11-04, 운용점 1 psi 확정 시점)

5×5 결정변수 격자(x_c = 0.1/0.2/0.5/0.7/0.8)가 쓰는 표를 **저장소에 커밋했다**.
HF 결과를 만드는 입력이므로, 표가 없으면 같은 그리드를 재현할 수 없다.

```
modes_seed.txt
modes_clamp_xc{0p1,0p2,0p5,0p7,0p8}_dc0p0556.txt
modes_xc{0p1,0p2,0p5,0p7,0p8}.txt (+ .meta.txt)
```

## 생성 (커밋된 표를 다시 만들 때)

`modes_seed.txt` 는 x_c 무관이지만, **합성표 헤더가 이 파일의 sha 를 참조**한다.
따라서 seed 가족을 재생성하면 그 뒤에 합성표를 **전부 재합성**해야 게이트가 통과한다.

```bat
:: 1) seed 가족 + 첫 설계점 (force)
python make_seed.py --x_c 0.1 --d_c 1e-4 --force

:: 2) 나머지 x_c (seed 가족은 건너뛰고 clamp + 합성만)
python make_seed.py --x_c 0.2 --d_c 1e-4
python make_seed.py --x_c 0.5 --d_c 1e-4
python make_seed.py --x_c 0.7 --d_c 1e-4
python make_seed.py --x_c 0.8 --d_c 1e-4

:: 3) 전부 게이트 확인 (Abaqus 불필요)
python -c "import aba_imperfection as a; src=('core/constants.py','run_abaqus_buckle.py'); [print(x, a.seed_gate(x,1e-4,'seeds/modes_xc%s.txt'%a.value_tag(x),src)[:2]) for x in (0.1,0.2,0.5,0.7,0.8)]"
```

`apply` 는 `--force` 로 불필요하게 seed 가족을 다시 쓰지 않는다. `--force` 는 clamp 표가
2차 격자(40401 노드)처럼 오염됐거나 buckle ODB 가 실패했을 때만 쓴다.

⚠️ **모드 부호/순서는 재생성마다 뒤집힐 수 있다.** 좌굴 고유값이 쌍으로 거의 겹치기
때문이다(실측 λ −1.58667e-3 / −1.58678e-3). 그래서 표를 재생성하면 HF 임퍼펙션 장이
바뀔 수 있고, 같은 그리드를 "다시 만든다"는 보장이 없다. **논문/비교 재현에는 커밋된 표를
그대로 쓴다.** 새로 만들었다면 `seed_gate` 통과와 `hf` 재현을 모두 확인한다.

## 게이트가 막는 것 (계획 Q2=a)

1. 표가 **다른 x_c/d_c** 에서 나온 것
2. 표를 만든 뒤 **buckle 소스 상수**가 바뀐 것(`aba_imperfection.model_fingerprint`)
3. **소스 모드표가 그 뒤 바뀐 것**(재추출했는데 재합성 안 함)
4. 모드표 노드 수 != 실제 적용 노드 수(메쉬/라벨 드리프트)

차단되면 HF 로그에 `!!! ERROR: 시드표 게이트 실패 - <사유>` 와 재생성 명령이 찍히고 잡은 제출되지 않는다.

## 커밋 정책

표는 텍스트(모드 4개 x 10,201 노드, 수백 KB)이고 HF 가 직접 의존하는 산출물이라
**저장소에 커밋한다**. 재생성 비용(Abaqus)이 크고 모드 부호가 바뀔 수 있으므로 이력에
남기는 편이 낫다. `.meta.txt` 도 함께 커밋해 출처를 남긴다.
