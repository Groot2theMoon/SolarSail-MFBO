# HF 임퍼펙션 시드표 (`code/seeds/`)

HF 포스트버클링(`run_abaqus.py`)이 읽는 좌굴모드 임퍼펙션 표가 여기 쌓인다.
경로 규약과 게이트는 `aba_imperfection.py`(stdlib) + `make_seed.py` 가 정본이다.

## 파일

| 파일 | 의미 | x_c 의존 |
|---|---|---|
| `modes_seed.txt` | `seed` 가족(꼭짓점 케이블만) 좌굴 4모드 | **아니오** — 한 번 만들면 재사용 |
| `modes_clamp_xc<NNN>_dc<NNN>.txt` | `clamp_only` 가족(클램프 케이블만) 모드 | 예 |
| `modes_xc<NNN>_dc<NNN>.txt` | **HF 가 읽는 합성표**(2+2) + 설계점/지문 헤더 | 예 |
| `modes_xc<NNN>_dc<NNN>.meta.txt` | 생성 시각·노드수·모드수·peak u3·게이트 결과 | 예 |

`<NNN>` 은 백분율이다(0.5 -> 050, 0.0556 -> 006). HF 는 자기 `x_c`·`d_c` 로
`modes_xc<NNN>_dc<NNN>.txt` 를 조립해 읽고, **없거나 지문이 어긋나면 제출 전에 중단**한다.

## 생성

```
python3 make_seed.py --x_c 0.5 --d_c 0.0556            # 캐시 있으면 건너뜀
python3 make_seed.py --x_c 0.25 --d_c 0.0556 --force   # 강제 재생성
python3 make_seed.py --x_c 0.5 --d_c 0.0556 --dry-run  # 명령만 출력(라이선스 0)
```

내부적으로 `run_abaqus_buckle.py`(2회: seed / clamp_only) -> `aba_mode_from_odb.py`(2회)
-> `combine_mode_tables.py`(1회) 순으로 돌린다. `seed` 가족은 x_c 무관이라 표가 있으면 건너뛴다.

## 게이트가 막는 것 (계획 Q2=a)

1. 표가 **다른 x_c/d_c** 에서 나온 것
2. 표를 만든 뒤 **buckle 소스 상수**가 바뀐 것(`aba_imperfection.model_fingerprint`)
3. **소스 모드표가 그 뒤 바뀐 것**(재추출했는데 재합성 안 함)
4. 모드표 노드 수 != 실제 적용 노드 수(메쉬/라벨 드리프트)

차단되면 HF 로그에 `!!! ERROR: 시드표 게이트 실패 - <사유>` 와 재생성 명령이 찍히고 잡은 제출되지 않는다.

## 커밋 정책

표는 텍스트(모드 4개 x 10,201 노드, 수백 KB)이고 HF 가 직접 의존하는 산출물이라
**저장소에 커밋한다**. 재생성 비용(Abaqus)이 크므로 이력에 남기는 편이 낫다.
`.meta.txt` 도 함께 커밋해 출처를 남긴다.
