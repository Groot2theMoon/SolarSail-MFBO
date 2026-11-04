#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""지문(fingerprint) API 단위시험 - Abaqus 불필요(라이선스 0).

    python3 test_fingerprint_api.py        # PASS/FAIL, 실패 시 EXIT=1

배경(2026-11-02): `aba_imperfection.py` 에 `model_fingerprint` 가 두 번 정의돼
뒤쪽 dict 버전이 앞쪽 스크립트 버전을 가렸다. 그 결과

    buckle_fingerprint()      -> TypeError  (2 인자 -> 1 인자 함수)
    model_fingerprint_line()  -> AttributeError ('str' object has no attribute 'items')

가 되어 `make_seed.py` 3단계(합성)가 항상 죽었다. 이름을 `source_model_fingerprint`
로 분리해 고쳤고, 같은 사고가 다시 나지 않도록 **중복 top-level 정의 자체**를 검사한다.
"""
import ast
import os
import sys

_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # [P1] tools/ -> code/
sys.path.insert(0, _CODE)
import aba_imperfection as A                                                   # noqa: E402

FAIL = []
N = [0]


def check(name, cond, detail=None):
    N[0] += 1
    ok = bool(cond)
    print("  %-4s %s%s" % ('OK' if ok else 'FAIL', name,
                           ('   <- ' + str(detail)) if (detail and not ok) else ''))
    if not ok:
        FAIL.append(name)


def main():
    here = _CODE          # [P1] tools/ -> code/ (검사 대상은 부모에 있다)
    src = os.path.join(here, 'aba_imperfection.py')
    buckle = os.path.join(here, 'run_abaqus_buckle.py')

    print("[1] 중복 top-level 정의 (섀도잉)")
    tree = ast.parse(open(src, encoding='utf-8').read())
    names = [n.name for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))]
    dup = sorted({n for n in names if names.count(n) > 1})
    check('aba_imperfection.py 중복 정의 없음', not dup, dup)

    print("\n[2] 소스 스크립트 지문 (16자)")
    fp = A.buckle_fingerprint(buckle)
    check('buckle_fingerprint() 가 문자열 반환', isinstance(fp, str), type(fp).__name__)
    check('16자 hex', isinstance(fp, str) and len(fp) == 16
          and all(c in '0123456789abcdef' for c in fp), fp)

    line = A.model_fingerprint_line(buckle)
    check('model_fingerprint_line() 문자열', isinstance(line, str), type(line).__name__)
    check('헤더 접두사', isinstance(line, str) and line.startswith(A.FINGERPRINT_PREFIX), line)

    print("\n[3] dict 지문 (12자)")
    d1 = A.model_fingerprint({'SEED_DIV': 200, 'THICKNESS': 5e-6})
    d2 = A.model_fingerprint({'THICKNESS': 5e-6, 'SEED_DIV': 200.0})
    d3 = A.model_fingerprint({'SEED_DIV': 632, 'THICKNESS': 5e-6})
    check('12자', isinstance(d1, str) and len(d1) == 12, d1)
    check('키 순서/표기 무관', d1 == d2, (d1, d2))
    check('값 변화에 민감', d1 != d3, (d1, d3))

    print("\n[4] 두 API 가 서로를 가리지 않음")
    check('source_model_fingerprint 존재', hasattr(A, 'source_model_fingerprint'))
    check('dict 버전과 소스 버전이 다른 함수',
          A.model_fingerprint is not A.source_model_fingerprint)

    print()
    if FAIL:
        print("RESULT: FAIL %d / %d - %s" % (len(FAIL), N[0], ', '.join(FAIL)))
        return 1
    print("RESULT: PASS %d/%d" % (N[0], N[0]))
    return 0


if __name__ == '__main__':
    sys.exit(main())
