#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""설계점 시드표 경로/지문 순수함수 단위시험 — Abaqus 불필요(라이선스 0).

    python3 test_seed_path.py        # PASS/FAIL, 실패 시 EXIT=1

검증 대상은 `aba_imperfection.py` 의 설계점 시드 API(2026-10-12)다.
계획: code/notes/plan-hf-seed-pipeline-2026-10-12.md
"""
import io
import os
import re
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aba_imperfection as A                                                   # noqa: E402

FAIL = []
N = [0]


def check(name, cond, detail=''):
    N[0] += 1
    ok = bool(cond)
    print("  %-4s %s%s" % ('OK' if ok else 'FAIL', name,
                           ('   <- ' + str(detail)) if (detail and not ok) else ''))
    if not ok:
        FAIL.append(name)


def _write(path, text):
    with io.open(path, 'w', encoding='utf-8') as f:
        f.write(text)


def main():
    tmp = tempfile.mkdtemp(prefix='seedpath_')

    print("--- 1) 파일명 규칙 ---")
    check('name 0.5/0.0556 -> xc050_dc006',
          A.seed_table_name(0.5, 0.0556) == 'modes_xc0p5_dc0p0556.txt',
          A.seed_table_name(0.5, 0.0556))
    check('name 0.25/0.1 -> xc025_dc010',
          A.seed_table_name(0.25, 0.1) == 'modes_xc0p25_dc0p1.txt',
          A.seed_table_name(0.25, 0.1))
    check('name 0.3/0.5 -> xc030_dc050',
          A.seed_table_name(0.3, 0.5) == 'modes_xc0p3_dc0p5.txt',
          A.seed_table_name(0.3, 0.5))

    print("--- 2) 경로 · env 우선 ---")
    check('path 기본 = seeds/modes_xc050_dc006.txt',
          A.seed_table_path(tmp, 0.5, 0.0556) == os.path.join(tmp, 'modes_xc0p5_dc0p0556.txt'),
          A.seed_table_path(tmp, 0.5, 0.0556))
    check('path env 오버라이드 우선',
          A.seed_table_path(tmp, 0.5, 0.0556, env='other/t.txt') == os.path.abspath('other/t.txt'))

    print("--- 3) 지문 민감도 ---")
    f0 = A.seed_fingerprint(0.5, 0.0556, 'abc', 's1', 'c1')
    check('지문 결정적', f0 == A.seed_fingerprint(0.5, 0.0556, 'abc', 's1', 'c1'))
    check('지문 16자', len(f0) == 16, f0)
    check('지문 x_c 민감', f0 != A.seed_fingerprint(0.6, 0.0556, 'abc', 's1', 'c1'))
    check('지문 d_c 민감', f0 != A.seed_fingerprint(0.5, 0.10, 'abc', 's1', 'c1'))
    check('지문 buckle 민감', f0 != A.seed_fingerprint(0.5, 0.0556, 'XYZ', 's1', 'c1'))
    check('지문 소스표 민감', f0 != A.seed_fingerprint(0.5, 0.0556, 'abc', 's2', 'c1'))

    print("--- 4) buckle 소스 지문(상수 드리프트 감지) ---")
    sc = os.path.join(tmp, 'fake_buckle.py')
    base = ("SIGMA0 = 700.0\nTHICKNESS = 5.0e-6\nUSE_GRID_MESH = True\nCLAMP_DC = 0.0556\n"
            "N_EIG_BUCKLE = 10\nBUCKLE_VECTORS = 25\n")
    _write(sc, base)
    b1 = A.buckle_fingerprint(sc)
    check('buckle 지문 생성(16자)', b1 and len(b1) == 16, b1)
    _write(sc, base.replace('BUCKLE_VECTORS = 25', 'BUCKLE_VECTORS = 26'))
    check('buckle 지문 값변경 감지', A.buckle_fingerprint(sc) != b1)
    _write(sc, base)
    check('buckle 지문 복원 시 동일', A.buckle_fingerprint(sc) == b1)

    print("--- 5) 헤더 왕복 + 게이트 ---")
    seed_t = os.path.join(tmp, 'modes_seed.txt')
    clamp_t = os.path.join(tmp, 'modes_clamp_xc050.txt')
    _write(seed_t, 'MODE 1\n1 1.0\nMODE 2\n1 -1.0\n')
    _write(clamp_t, 'MODE 1\n1 1.0\n')
    fp, hdr = A.seed_header_lines(0.5, 0.0556, b1, seed_t, '1,2', clamp_t, '1')
    tbl = os.path.join(tmp, 'modes_xc0p5_dc0p0556.txt')
    _write(tbl, '\n'.join(hdr) + '\nMODE 1\n1 1.0\nMODE 2\n1 -1.0\n')
    check('헤더에 지문 기록', fp in '\n'.join(hdr))
    check('table_fingerprint 호환(기존 리더)', A.table_fingerprint(tbl) == fp,
          A.table_fingerprint(tbl))
    ok, why, info = A.seed_gate(0.5, 0.0556, tbl, sc)
    check('게이트 통과', ok, why)
    check('헤더 x_c 읽기', info.get('x_c') == '0.5', info.get('x_c'))
    check('헤더 source 표 이름', info.get('seed') == 'modes_seed.txt', info.get('seed'))
    ok2, why2, _ = A.seed_gate(0.6, 0.0556, tbl, sc)
    check('게이트 x_c 불일치 차단', (not ok2) and ('x_c' in why2), why2)
    ok3, why3, _ = A.seed_gate(0.5, 0.10, tbl, sc)
    check('게이트 d_c 불일치 차단', (not ok3) and ('d_c' in why3), why3)
    _write(sc, base.replace('N_EIG_BUCKLE = 10', 'N_EIG_BUCKLE = 11'))
    ok4, why4, _ = A.seed_gate(0.5, 0.0556, tbl, sc)
    check('게이트 buckle 드리프트 차단', (not ok4) and ('바뀜' in why4), why4)
    _write(sc, base)
    # 소스 표가 바뀌면(내용 변경) 지문 재계산이 어긋나야 한다
    _write(seed_t, 'MODE 1\n1 2.0\nMODE 2\n1 -1.0\n')
    ok5, why5, _ = A.seed_gate(0.5, 0.0556, tbl, sc)
    check('게이트 소스표 변조 차단', (not ok5) and ('modes_seed.txt' in why5), why5)
    # 구버전 표(헤더 없음)
    old = os.path.join(tmp, 'old.txt')
    _write(old, 'MODE 1\n1 1.0\n')
    ok6, why6, _ = A.seed_gate(0.5, 0.0556, old, sc)
    check('게이트 구버전표 차단', not ok6, why6)

    print("--- 10) 로케일 C(cp949 상당)에서도 UTF-8 모드표를 읽는가 (§19.71) ---")
    _d = tempfile.mkdtemp(prefix='utf8tab_')
    _p = os.path.join(_d, 't.txt')
    _write(_p, '# 헤더 — 한글 주석\nMODE 1\n1 1.0\n2 -0.5\n')
    _code = ("import sys,aba_imperfection as A;"
             "t=A.load_mode_table(r'%s',(1,));print(sorted(t.keys()),len(t[1]))" % _p)
    _r = subprocess.run([sys.executable, '-c', _code], capture_output=True, text=True,
                        env=dict(os.environ, LC_ALL='C', LANG='C', PYTHONUTF8='0'),
                        cwd=os.path.dirname(os.path.abspath(__file__)))
    check('locale=C 에서 UTF-8 모드표 읽기', _r.returncode == 0 and '[1] 2' in (_r.stdout or ''),
          (_r.stderr or '')[-140:])

    print("--- 9) 자식 프로세스 인코딩 고정 (§19.68) ---")
    try:
        import make_seed as MS
        _rc, out = MS.sh("python3 -c \"import os;print(os.environ.get('PYTHONIOENCODING'))\"")
        check('자식 PYTHONIOENCODING=utf-8 전달', 'utf-8' in out, out.strip()[-60:])
        _rc, out = MS.sh("python3 -c \"import sys;print(str(chr(8212)))\"")   # U+2014
        check('자식이 em-dash 출력해도 정상 종료', _rc == 0 and '\u2014' in out, out.strip()[-60:])
    except ImportError as e:
        check('make_seed import', False, e)

    print("--- 8) combine_mode_tables.parse_spec (드라이브 콜론 회귀) ---")
    try:
        import combine_mode_tables as CMT
        check('상대경로 스펙', CMT.parse_spec('1.0:seeds/a.txt:1,2') == (1.0, 'seeds/a.txt', [1, 2]),
              CMT.parse_spec('1.0:seeds/a.txt:1,2'))
        check('Windows 절대경로 스펙', CMT.parse_spec(r'1.0:C:\x\a.txt:1,2') == (1.0, r'C:\x\a.txt', [1, 2]),
              CMT.parse_spec(r'1.0:C:\x\a.txt:1,2'))
        check('유닉스 절대경로 스펙', CMT.parse_spec('0.5:/tmp/a.txt:3') == (0.5, '/tmp/a.txt', [3]),
              CMT.parse_spec('0.5:/tmp/a.txt:3'))
    except ImportError as e:
        check('combine_mode_tables import', False, e)

    print("--- 7) 설계점 태그 · LF 산출물 이름 (§19.66) ---")
    check("point_tag(0.5, 0.5)", A.point_tag(0.5, 0.5) == '_xc0p5_dc0p5', A.point_tag(0.5, 0.5))
    check("point_tag(0.25, 0.0556)", A.point_tag(0.25, 0.0556) == '_xc0p25_dc0p0556', A.point_tag(0.25, 0.0556))
    # ---- [§19.76] 태그 충돌 가드: 서로 다른 설계점이 같은 이름을 쓰면 안 된다 ----
    check("[§19.76] 0 vs 0.001 충돌 없음",
          A.seed_table_name(0.5, 0.0) != A.seed_table_name(0.5, 0.001),
          A.seed_table_name(0.5, 0.0) + ' / ' + A.seed_table_name(0.5, 0.001))
    check("[§19.76] x_c 0.1234 vs 0.1249 충돌 없음",
          A.seed_table_name(0.1234, 0.5) != A.seed_table_name(0.1249, 0.5),
          A.seed_table_name(0.1234, 0.5) + ' / ' + A.seed_table_name(0.1249, 0.5))
    check("[§19.76] d_c 1e-04 vs 0 충돌 없음",
          A.seed_table_name(0.5, 1.0e-4) != A.seed_table_name(0.5, 0.0),
          A.seed_table_name(0.5, 1.0e-4))
    _dir = os.path.dirname(os.path.abspath(__file__))
    _bs = io.open(os.path.join(_dir, 'run_abaqus_buckle.py'), encoding='utf-8').read()
    _ms = io.open(os.path.join(_dir, 'make_seed.py'), encoding='utf-8').read()
    _fmt = 'Buckle_xc%s_d%03dum_dc%s_'
    check("[§19.76] 잡이름 형식 buckle↔make_seed 동일",
          (_fmt in _bs) and (_fmt in _ms), 'buckle=%s make_seed=%s' % (_fmt in _bs, _fmt in _ms))
    check("lf_odb_name(0.5, 0.5)", A.lf_odb_name(0.5, 0.5) == 'LF_Analysis_xc050_dc050.odb',
          A.lf_odb_name(0.5, 0.5))
    check("lf_job_name(0.5, 0.5)", A.lf_job_name(0.5, 0.5) == 'LF_Analysis_xc050_dc050',
          A.lf_job_name(0.5, 0.5))
    check("LF 이름 x_c 로 구분", A.lf_odb_name(0.5, 0.5) != A.lf_odb_name(0.6, 0.5))
    check("LF 이름 d_c 로 구분", A.lf_odb_name(0.5, 0.5) != A.lf_odb_name(0.5, 0.6))
    check("구 고정명과 충돌 없음", A.lf_odb_name(0.5, 0.5) != 'LF_Analysis.odb')

    print("--- 6b) run_abaqus.py 는 execfile 안전해야 (__file__ 금지) ---")
    _here = os.path.dirname(os.path.abspath(__file__))
    _lines = io.open(os.path.join(_here, 'run_abaqus.py'), encoding='utf-8').read().split('\n')
    _i_def = next((k for k, l in enumerate(_lines) if l.startswith('def _resolve_here')), -1)
    _i_here = next((k for k, l in enumerate(_lines) if l.startswith('_HERE = ')), -1)
    _bad = [k + 1 for k, l in enumerate(_lines)
            if 'abspath(__file__)' in l and not (_i_def < k < _i_here)]
    check('abspath(__file__) 는 _resolve_here 안에만', not _bad, _bad)

    print("--- 6) run_abaqus.py 의 aba_imperfection import 이름 실재성 ---")
    here = os.path.dirname(os.path.abspath(__file__))
    raw = io.open(os.path.join(here, 'run_abaqus.py'), encoding='utf-8').read()
    m = re.search(r'from aba_imperfection import \((.*?)\)', raw, re.S)
    check('run_abaqus.py import 블록 발견', bool(m))
    if m:
        names = [n.strip().split('#')[0].strip() for n in m.group(1).split(',')]
        names = [n for n in names if n]
        missing = [n for n in names if not hasattr(A, n)]
        check('import 이름 전부 실재(%d개)' % len(names), not missing, missing)

    print()
    if FAIL:
        print("RESULT: FAIL %d / %d — %s" % (len(FAIL), N[0], ', '.join(FAIL)))
        return 1
    print("RESULT: PASS %d/%d" % (N[0], N[0]))
    return 0


if __name__ == '__main__':
    sys.exit(main())
