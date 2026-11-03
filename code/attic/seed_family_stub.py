#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""[§19.104] NO_CLAMP baseline 임퍼펙션 표 선택 하네스 - Abaqus 없이 **실제 실행**한다.

무엇을 검증하나 (실제 Abaqus 런이 아니라 '배선/선택 로직' 검증이다)
    run_abaqus.py 에서 임퍼펙션 표 선택 배선(`_available_mode_numbers` def ~
    `if IMPERFECTION_MODE == 'odb_table':` 블록 끝)을 **소스에서 떼어** 임시 디렉터리에서
    그대로 실행한다. 검증 항목:

      (i)   NO_CLAMP=False -> 합성표 이름(modes_xc<...>.txt)과 기존 모드 목록([1,2,3,4])이
            **바이트 동일**하다. §19.104 를 도입한 커밋의 **부모 소스**(없으면 HEAD)에서 같은
            블록을 떼어 실행한 결과와 로그/표/모드를 문자열 비교한다(회귀).
      (ii)  NO_CLAMP=True -> **seed 가족 표**(seeds/modes_seed.txt) + 그 표의 상위 4모드만 쓴다.
            합성표도 env 오버라이드(SOLARSAIL_MODE_TABLE)도 무시되고, 로그에 modes_clamp 가 없다.
      (iii) NO_CLAMP=True 인데 seed 가족 표가 없으면 **명확한 에러 + 종료(1)** - 조용히 합성표로
            되돌아가지 않는다(오염 재발 차단). 모드가 4개 미만이어도 같다.
      (iv)  로그 줄에 해결된 표/모드/출처가 들어간다.
      (v)   seed_family_gate 단위 검증(합성표 이름 / 헤더 오염 / 모드 부족 / 파일 없음).
      (vi)  새로 쓴 줄에 금지문자(em-dash/en-dash/화살표/smart quote/NBSP/bullet/ellipsis) 0,
            출력 로그 cp949 인코딩 가능(§19.95 교훈), 기존 스크립트 git diff 무수정.

사용
    python3 seed_family_stub.py        # PASS/FAIL, 실패 시 EXIT=1
Abaqus/라이선스 불필요(numpy 도 불필요).
종료코드: 0 OK / 1 FAIL
"""
from __future__ import print_function

import io
import os
import re
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import aba_imperfection as A                                                     # noqa: E402

RA = os.path.join(HERE, 'run_abaqus.py')
# 금지 문자(§19.95: cp949 콘솔에서 UnicodeEncodeError 로 도구 자체가 죽는다)
FORBIDDEN_RE = re.compile(u'[\u2010\u2011\u2012\u2013\u2014\u2015'
                          u'\u2190\u2191\u2192\u2193\u21d2'
                          u'\u2018\u2019\u201c\u201d'
                          u'\u00a0\u2022\u2026]')

FAIL = []
N = [0]


def check(name, cond, detail=''):
    N[0] += 1
    ok = bool(cond)
    print("  %-4s %s%s" % ('OK' if ok else 'FAIL', name,
                           ('   <- ' + str(detail)) if (detail and not ok) else ''))
    if not ok:
        FAIL.append(name)
    return ok


# ------------------------------------------------------------------ 소스에서 배선 떼기
def read_source(path):
    with io.open(path, encoding='utf-8') as f:
        return f.read()


def extract_region(src):
    """`_available_mode_numbers` def 부터 odb_table 블록 끝까지. 해석부가 블록 **밖**(앞)이면 붙인다.

    (구버전/신버전 모두 같은 추출기로 처리해 '구버전 vs 신버전'을 같은 조건에서 실행한다.)
    """
    m = re.search(r'\ndef _available_mode_numbers\(.*?'
                  r'(?=\nelif IMPERFECTION_MODE == .odb_direct.:)', src, re.S)
    if not m:
        return None
    blk = m.group(0)
    r = re.search(r'^SEEDS_DIR = os\.path\.join\(_HERE, SEED_DIR_NAME\)[^\n]*\n'
                  r'^MODE_TABLE = [^\n]*\n', src, re.M)
    if r and r.start() < m.start():
        blk = '\n' + r.group(0) + blk          # 구버전: 해석부가 def 앞에 있다 -> 앞에 붙인다
    return blk


def run_region(region_src, code_dir, seeds_dir, x_c=0.5, d_c=0.0556, no_clamp=False,
               env_table=None, imperf_modes=(1, 2, 3, 4), cwd=None):
    """뗀 배선 블록을 실제로 실행하고 (ns, logs, exc) 반환. SystemExit 는 잡아서 exc 로 돌려준다."""
    src, nflip = region_src, 0
    if no_clamp:
        src, nflip = re.subn(r'^NO_CLAMP = False$', 'NO_CLAMP = True',
                             region_src, count=1, flags=re.M)
    logs = []

    def _emit(*a):
        logs.append(' '.join(str(x) for x in a))

    ns = {'__name__': 'run_abaqus_wiring', 'os': os, 'sys': sys, 'time': time,
          '_HERE': code_dir, 'SEED_DIR_NAME': A.SEED_DIR_NAME,
          'BUCKLE_SCRIPT_NAME': A.BUCKLE_SCRIPT_NAME,
          'x_c': x_c, 'd_c': d_c,
          'IMPERFECTION_MODE': 'odb_table', 'IMPERFECTION_MODES': tuple(imperf_modes),
          'IMPERFECTION_MAX_MODES': 4, 'NO_CLAMP': no_clamp,
          'emit': _emit, 'load_mode_table': A.load_mode_table, 'seed_gate': A.seed_gate,
          'seed_table_path': A.seed_table_path, 'mode_table_info': A.mode_table_info,
          'ImperfectionSourceError': A.ImperfectionSourceError}
    old = os.getcwd()
    old_env = os.environ.get('SOLARSAIL_MODE_TABLE')
    try:
        if env_table is None:
            os.environ.pop('SOLARSAIL_MODE_TABLE', None)
        else:
            os.environ['SOLARSAIL_MODE_TABLE'] = env_table
        os.chdir(cwd or code_dir)
        exc = None
        try:
            exec(compile(src, '<run_abaqus_wiring>', 'exec'), ns)
        except SystemExit as e:
            exc = e
    finally:
        os.chdir(old)
        if old_env is None:
            os.environ.pop('SOLARSAIL_MODE_TABLE', None)
        else:
            os.environ['SOLARSAIL_MODE_TABLE'] = old_env
    ns['_nflip'] = nflip
    return ns, logs, exc


# ------------------------------------------------------------------ 픽스처
def write(path, text):
    d = os.path.dirname(path)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with io.open(path, 'w', encoding='utf-8') as f:
        f.write(text)


def make_buckle_script(code_dir):
    """buckle 지문 키를 가진 가짜 buckle 소스(합성표 게이트 판정에 필요)."""
    write(os.path.join(code_dir, A.BUCKLE_SCRIPT_NAME),
          "SIGMA0 = 700.0\nTHICKNESS = 5.0e-6\nUSE_GRID_MESH = True\nCLAMP_DC = 0.0556\n"
          "N_EIG_BUCKLE = 10\nBUCKLE_VECTORS = 25\n")


def seed_family_text(n_modes=4, odb_name='Buckle_xc0p5_d018um_dc0_f100_s4_grid_seed_normal_ex210_plv01000.odb',
                     extra_header=()):
    """aba_mode_from_odb.write_mode_table 가 쓰는 형식의 seed 가족 표."""
    lines = ['# source=/work/buckle/%s' % odb_name,
             '# step=Step-Buckle',
             '# instance=MEMBRANE-1',
             '# u3 = per-mode max|u3|=1 normalized, sign fixed at max node',
             '# table_modes=%d / source_frames=%d / lambdas=1.0e-04,2.0e-04' % (n_modes, n_modes + 1)]
    lines += list(extra_header)
    out = ['\n'.join(lines), '']
    for m in range(1, n_modes + 1):
        out.append('MODE %d' % m)
        for lab in (1, 2, 3):
            out.append('%d %.9e' % (lab, (1.0 if lab == 1 else -0.25) / m))
    return '\n'.join(out) + '\n'


def write_seed_family(seeds_dir, n_modes=4, extra_header=()):
    p = os.path.join(seeds_dir, 'modes_seed.txt')
    write(p, seed_family_text(n_modes, extra_header=extra_header))
    return p


def write_source_tables(seeds_dir):
    """합성표의 **소스**가 되는 seed/clamp 가족 표(게이트가 해시를 대조한다)."""
    seed_t = os.path.join(seeds_dir, 'modes_seed.txt')
    clamp_t = os.path.join(seeds_dir, 'modes_clamp_xc0p5_dc0.txt')
    write(seed_t, seed_family_text())
    write(clamp_t, 'MODE 1\n1 1.0\n2 -0.5\nMODE 2\n1 0.5\n2 1.0\n')
    return seed_t, clamp_t


def write_composite(seeds_dir, code_dir, x_c, d_c, seed_tbl, clamp_tbl):
    """합성표(2+2) - HF 가 NO_CLAMP=False 에서 읽는 표. 진짜 지문 게이트를 통과해야 한다."""
    fp = A.buckle_fingerprint(os.path.join(code_dir, A.BUCKLE_SCRIPT_NAME))
    _fp2, hdr = A.seed_header_lines(x_c, d_c, fp, seed_tbl, '1,2', clamp_tbl, '1,2')
    body = []
    for m in range(1, 5):
        body.append('MODE %d' % m)
        body.append('1 %.9e' % (1.0 / m))
        body.append('2 %.9e' % (-0.5 / m))
    p = os.path.join(seeds_dir, A.composite_table_name(x_c))
    write(p, '\n'.join(hdr) + '\n' + '\n'.join(body) + '\n')
    return p


def _norm(line, code_dir):
    """임시 디렉터리 경로만 지운다(로그 내용 비교용)."""
    return line.replace(code_dir, '<CODE>')


def fixture(tag):
    """(code_dir, seeds_dir, work_dir) - code_dir 에 seeds/ 가 있고 가짜 buckle 소스가 있다."""
    base = tempfile.mkdtemp(prefix=tag + '_')
    code_dir = os.path.join(base, 'code')
    seeds_dir = os.path.join(code_dir, 'seeds')
    os.makedirs(seeds_dir)
    make_buckle_script(code_dir)
    return code_dir, seeds_dir, base


# ------------------------------------------------------------------ 회귀 기준(구버전 소스)
def baseline_source():
    """§19.104 를 도입한 커밋이 있으면 그 **부모** 소스, 없으면 HEAD 소스. (출처 문자열, 텍스트)."""
    try:
        h = subprocess.run(['git', '-C', ROOT, 'log', '--format=%H', '-S19.104', '--',
                            'code/run_abaqus.py'],
                           capture_output=True, text=True).stdout.split()
    except Exception:
        h = []
    if h:
        c = h[-1]                              # 가장 오래된 도입 커밋
        r = subprocess.run(['git', '-C', ROOT, 'show', '%s^:code/run_abaqus.py' % c],
                           capture_output=True, text=True)
        if r.returncode == 0 and r.stdout:
            return '%s^ (도입 커밋 %s)' % (c[:8], c[:8]), r.stdout
    r = subprocess.run(['git', '-C', ROOT, 'show', 'HEAD:code/run_abaqus.py'],
                       capture_output=True, text=True)
    if r.returncode == 0 and r.stdout:
        return 'HEAD', r.stdout
    return None, None


def added_lines(path='code/run_abaqus.py'):
    """이번 변경(§19.104)으로 **새로 추가된** 줄. 커밋 전이면 worktree vs HEAD, 커밋 후면 HEAD vs HEAD~1."""
    st = subprocess.run(['git', '-C', ROOT, 'status', '--porcelain'],
                        capture_output=True, text=True).stdout
    if st.strip():
        args = ['git', '-C', ROOT, 'diff', '-U0', 'HEAD', '--', path]
    else:
        msg = subprocess.run(['git', '-C', ROOT, 'log', '-1', '--format=%s'],
                             capture_output=True, text=True).stdout
        if '19.104' not in msg:
            return []
        args = ['git', '-C', ROOT, 'diff', '-U0', 'HEAD~1', 'HEAD', '--', path]
    d = subprocess.run(args, capture_output=True, text=True).stdout
    return [l[1:] for l in d.splitlines() if l.startswith('+') and not l.startswith('+++')]


def changed_py_files():
    """이번 변경(§19.104)으로 바뀐 code/*.py - 커밋 전이면 worktree, 커밋 후면 HEAD 커밋 대조."""
    files = set()
    p = subprocess.run(['git', '-C', ROOT, 'status', '--porcelain'],
                       capture_output=True, text=True).stdout
    for ln in p.splitlines():
        path = ln[3:].strip()
        if path.endswith('.py'):
            files.add(path)
    if not files:
        msg = subprocess.run(['git', '-C', ROOT, 'log', '-1', '--format=%s'],
                             capture_output=True, text=True).stdout
        if '19.104' in msg:
            d = subprocess.run(['git', '-C', ROOT, 'diff', '--name-only', 'HEAD~1', 'HEAD'],
                               capture_output=True, text=True).stdout
            files = set(l for l in d.splitlines() if l.endswith('.py'))
    return sorted(files)


# ------------------------------------------------------------------ 본체
def main():
    ra = read_source(RA)
    region = extract_region(ra)
    print("--- 0) 소스/추출 ---")
    if not check('run_abaqus.py 에서 배선 블록 추출', bool(region)):
        return 1
    check('배선 블록에 seed_family_gate 포함', 'def seed_family_gate' in region)
    check('배선 블록에 NO_CLAMP 표 분기 포함', 'NO_CLAMP = False' in region)
    check("SEED_FAMILY_NAME = 'modes_seed.txt'", "SEED_FAMILY_NAME = 'modes_seed.txt'" in ra)
    check('IMPERFECTION_MODES 기본 = (1, 2, 3, 4)', 'IMPERFECTION_MODES = (1, 2, 3, 4)' in ra)
    check('IMPERFECTION_MAX_MODES 기본 = 4', 'IMPERFECTION_MAX_MODES = 4' in ra)
    check('NO_CLAMP 정의가 파일에 정확히 1회', len(re.findall(r'^NO_CLAMP = ', ra, re.M)) == 1,
          len(re.findall(r'^NO_CLAMP = ', ra, re.M)))
    check('NO_CLAMP 오버라이드 대상이 SEED_FAMILY_TABLE',
          re.search(r'if NO_CLAMP:\n(?:[ \t]*#[^\n]*\n)*[ \t]*MODE_TABLE = SEED_FAMILY_TABLE\n',
                    region) is not None)
    check('NO_CLAMP=False 는 seed_table_path(합성표) 그대로',
          re.search(r'else:\n    MODE_TABLE = seed_table_path\(SEEDS_DIR, x_c, d_c, env=_MODE_TABLE_ENV\)\n',
                    region) is not None)

    # --------------------------------------------------------------- (i) 회귀
    print("--- (i) NO_CLAMP=False 회귀: 합성표 이름 + 기존 모드 목록 ---")
    code_dir, seeds_dir, base = fixture('seedreg')
    seed_t, clamp_t = write_source_tables(seeds_dir)
    comp = write_composite(seeds_dir, code_dir, 0.5, 0.0556, seed_t, clamp_t)
    ns, logs, exc = run_region(region, code_dir, seeds_dir, no_clamp=False)
    check('NO_CLAMP=False: 예외 없이 통과', exc is None, exc)
    check('NO_CLAMP=False: 표 = seeds/modes_xc0p5.txt (합성표)',
          ns.get('MODE_TABLE') == comp, ns.get('MODE_TABLE'))
    check('NO_CLAMP=False: 모드 = (1, 2, 3, 4)', tuple(ns.get('IMPERFECTION_MODES', ())) == (1, 2, 3, 4),
          ns.get('IMPERFECTION_MODES'))
    check('NO_CLAMP=False: 합성표 4모드 로드', sorted((ns.get('_mode_table') or {}).keys()) == [1, 2, 3, 4],
          sorted((ns.get('_mode_table') or {}).keys()))
    check('NO_CLAMP=False: 지문 게이트 통과 로그', any('게이트 OK' in l for l in logs), logs)
    check('NO_CLAMP=False: NO_CLAMP 전용 로그 없음',
          not any('NO_CLAMP=True' in l for l in logs), logs)

    bsrc_label, bsrc_txt = baseline_source()
    print("--- (i-b) 회귀 기준 소스와 **바이트 비교** (기준=%s) ---" % bsrc_label)
    if not bsrc_txt:
        check('회귀 기준 소스 확보(git)', False, 'git show 실패')
    else:
        breg = extract_region(bsrc_txt)
        if not check('기준 소스에서 배선 블록 추출', bool(breg)):
            pass
        elif breg.strip() == region.strip():
            print("      (참고) 기준 소스가 이미 수정본이다 -> 회귀 판정은 위 기대값으로 한다")
        else:
            c2, s2, b2 = fixture('seedreg_base')
            st2, cl2 = write_source_tables(s2)
            comp2 = write_composite(s2, c2, 0.5, 0.0556, st2, cl2)
            ns2, logs2, exc2 = run_region(breg, c2, s2, no_clamp=False)
            check('기준 소스: 예외 없이 통과', exc2 is None, exc2)
            _n1 = [_norm(l, code_dir) for l in logs]
            _n2 = [_norm(l, c2) for l in logs2]
            check('(i) 합성표 이름 바이트 동일(%s)' % os.path.basename(comp2),
                  os.path.basename(ns.get('MODE_TABLE', '')) == os.path.basename(ns2.get('MODE_TABLE', '')),
                  (ns.get('MODE_TABLE'), ns2.get('MODE_TABLE')))
            check('(i) 주입 모드 목록 바이트 동일',
                  tuple(ns.get('IMPERFECTION_MODES', ())) == tuple(ns2.get('IMPERFECTION_MODES', ())),
                  (ns.get('IMPERFECTION_MODES'), ns2.get('IMPERFECTION_MODES')))
            check('(i) 로드된 모드 키 동일',
                  sorted((ns.get('_mode_table') or {}).keys()) == sorted((ns2.get('_mode_table') or {}).keys()),
                  (sorted((ns.get('_mode_table') or {}).keys()), sorted((ns2.get('_mode_table') or {}).keys())))
            check('(i) 로그 줄 수 동일', len(_n1) == len(_n2), (len(_n1), len(_n2)))
            check('(i) 로그 내용 바이트 동일(임시경로 정규화 후)', _n1 == _n2,
                  [l for l in _n1 if l not in _n2] + [l for l in _n2 if l not in _n1])

    # --------------------------------------------------------------- (ii) NO_CLAMP=True
    print("--- (ii) NO_CLAMP=True: 순수 seed 가족 표 ---")
    code_dir, seeds_dir, base = fixture('seedfam')
    seed_t, clamp_t = write_source_tables(seeds_dir)
    comp = write_composite(seeds_dir, code_dir, 0.5, 0.0556, seed_t, clamp_t)   # 합성표도 함께 존재
    ns, logs, exc = run_region(region, code_dir, seeds_dir, no_clamp=True)
    check('NO_CLAMP=True: 치환(NO_CLAMP = False -> True)이 정확히 1회', ns.get('_nflip') == 1, ns.get('_nflip'))
    check('NO_CLAMP=True: 예외 없이 통과', exc is None, exc)
    check('NO_CLAMP=True: 표 = seeds/modes_seed.txt (seed 가족)',
          ns.get('MODE_TABLE') == seed_t, ns.get('MODE_TABLE'))
    check('NO_CLAMP=True: 합성표(modes_xc0p5.txt)를 집지 않는다',
          os.path.basename(ns.get('MODE_TABLE', '')) != os.path.basename(comp), ns.get('MODE_TABLE'))
    check('NO_CLAMP=True: 모드 = (1, 2, 3, 4) (seed 표 상위 4개)',
          tuple(ns.get('IMPERFECTION_MODES', ())) == (1, 2, 3, 4), ns.get('IMPERFECTION_MODES'))
    check('NO_CLAMP=True: seed 표 4모드만 로드', sorted((ns.get('_mode_table') or {}).keys()) == [1, 2, 3, 4],
          sorted((ns.get('_mode_table') or {}).keys()))
    check('NO_CLAMP=True: 로그/표 어디에도 modes_clamp 없음(문자열 검사)',
          not any('modes_clamp' in l for l in logs) and 'modes_clamp' not in (ns.get('MODE_TABLE') or ''),
          [l for l in logs if 'modes_clamp' in l])
    check('NO_CLAMP=True: 합성 지문 게이트를 돌리지 않는다',
          not any('게이트 OK' in l for l in logs) and any('지문 게이트 생략' in l for l in logs), logs)
    #   (iv) 로그 줄에 해결된 표/모드/출처
    _sflog = [l for l in logs if l.startswith('[SEED] NO_CLAMP=True')]
    check('(iv) [SEED] NO_CLAMP=True 로그 1줄 존재', len(_sflog) == 1, _sflog)
    if _sflog:
        _l = _sflog[0]
        check('(iv) 로그에 해결된 표 경로', ('table=%s' % seed_t) in _l, _l)
        check('(iv) 로그에 주입 모드 [1, 2, 3, 4]', 'modes=[1, 2, 3, 4]' in _l, _l)
        check('(iv) 로그에 출처(source ODB)',
              ('source=/work/buckle/Buckle_xc0p5_d018um_dc0_f100_s4_grid_seed_normal_ex210_plv01000.odb' in _l),
              _l)
        check('(iv) 로그에 "순수 seed 가족 사용(clamp 모드 제외)"',
              '순수 seed 가족 사용(clamp 모드 제외)' in _l, _l)
    check('(iv) 모드표 출처 로그에 seed ODB 가 찍힌다',
          any(l.startswith('[IMPERFECTION] 모드표 출처') and 'source=/work/buckle/' in l for l in logs), logs)

    # env 오버라이드가 있어도 무시된다
    ns2, logs2, exc2 = run_region(region, code_dir, seeds_dir, no_clamp=True, env_table=comp)
    check('NO_CLAMP=True: env SOLARSAIL_MODE_TABLE(합성표) 무시',
          ns2.get('MODE_TABLE') == seed_t, ns2.get('MODE_TABLE'))
    check('NO_CLAMP=True: env 무시 경로도 예외 없음', exc2 is None, exc2)

    # --------------------------------------------------------------- (iii) 실패 경로
    print("--- (iii) seed 가족 표가 없거나 모드가 부족하면 중단 ---")
    code_dir, seeds_dir, base = fixture('seedmiss')
    seed_t, clamp_t = write_source_tables(seeds_dir)
    comp = write_composite(seeds_dir, code_dir, 0.5, 0.0556, seed_t, clamp_t)
    os.remove(seed_t)                                   # seed 가족 표만 없앤다
    ns, logs, exc = run_region(region, code_dir, seeds_dir, no_clamp=True)
    check('seed 가족 표 없음: SystemExit', isinstance(exc, SystemExit), exc)
    check('seed 가족 표 없음: 종료코드 1', getattr(exc, 'code', None) == 1, getattr(exc, 'code', None))
    check('seed 가족 표 없음: 명확한 에러 로그',
          any(l.startswith('!!! ERROR:') and 'seed 가족' in l for l in logs), logs)
    check('seed 가족 표 없음: ASCII 에러 라인([ERROR-EN])', any('[ERROR-EN]' in l for l in logs), logs)
    check('seed 가족 표 없음: 재생성 명령 안내(make_seed.py)',
          any('make_seed.py --x_c 0.5 --d_c 0.0556' in l for l in logs), logs)
    check('seed 가족 표 없음: 합성표로 조용히 되돌아가지 않음',
          ns.get('MODE_TABLE') == os.path.join(seeds_dir, 'modes_seed.txt') and comp not in logs,
          ns.get('MODE_TABLE'))
    check('seed 가족 표 없음: 모드 로드에 도달하지 않는다(_mode_table 없음)',
          ns.get('_mode_table') is None, ns.get('_mode_table'))

    # 모드 부족
    code_dir, seeds_dir, base = fixture('seedshort')
    write_seed_family(seeds_dir, n_modes=2)
    ns, logs, exc = run_region(region, code_dir, seeds_dir, no_clamp=True)
    check('seed 표 모드 2개: SystemExit(1)', isinstance(exc, SystemExit) and exc.code == 1, exc)
    check('seed 표 모드 2개: "모드 부족" 에러', any('모드 부족' in l for l in logs), logs)

    # 헤더 오염(합성/클램프 흔적이 있는 modes_seed.txt)
    code_dir, seeds_dir, base = fixture('seeddirty')
    write_seed_family(seeds_dir, n_modes=4,
                      extra_header=('# FINGERPRINT sha1=deadbeef src=run_abaqus_buckle.py',
                                    '# source=combined: modes_seed.txt[1,2] + '
                                    'modes_clamp_xc0p5_dc0.txt[1,2]'))
    ns, logs, exc = run_region(region, code_dir, seeds_dir, no_clamp=True)
    check('seed 표 헤더 오염(combined/modes_clamp): SystemExit(1)',
          isinstance(exc, SystemExit) and exc.code == 1, exc)
    check('seed 표 헤더 오염: 사유에 클램프/합성 흔적',
          any('클램프/합성 흔적' in l for l in logs), logs)

    # --------------------------------------------------------------- (v) gate 단위
    print("--- (v) seed_family_gate 단위 ---")
    code_dir, seeds_dir, base = fixture('seedgate')
    good = write_seed_family(seeds_dir)
    ns = {'os': os, 'SEED_FAMILY_NAME': 'modes_seed.txt', 'mode_table_info': A.mode_table_info}
    exec(compile(re.search(r'\ndef _available_mode_numbers\(.*?(?=\n# ====|\nSEED_FAMILY_NAME = )',
                           ra, re.S).group(0), '<gate>', 'exec'), ns)
    exec(compile(re.search(r'\ndef seed_family_gate\(.*?(?=\n\n#   \[2026-10-12\] 설계점별 시드표)',
                           ra, re.S).group(0), '<gate>', 'exec'), ns)
    gate = ns['seed_family_gate']
    ok, why, s = gate(good, 4)
    check('gate: 정상 seed 가족 표 통과', ok, why)
    check('gate: 요약에 표/모드/출처', s.get('table') == good and s.get('modes') == [1, 2, 3, 4]
          and s.get('source', '').endswith('.odb'), s)
    #   합성표 이름의 **실제 파일**을 만들어 '이름 게이트'만 걸리게 한다(파일 없음 사유와 구분).
    write(os.path.join(seeds_dir, 'modes_xc0p5.txt'), seed_family_text())
    ok, why, _ = gate(os.path.join(seeds_dir, 'modes_xc0p5.txt'), 4)
    check('gate: 합성표 이름은 거부', (not ok) and ('seed 가족 표가 아니다' in why), why)
    ok, why, _ = gate(os.path.join(seeds_dir, 'nope.txt'), 4)
    check('gate: 파일 없음 거부', (not ok) and ('없음' in why), why)
    ok, why, _ = gate(good, 5)
    check('gate: 모드 부족 거부(5개 요구)', (not ok) and ('모드 부족' in why), why)

    # --------------------------------------------------------------- (vi) 문자/스키마/무수정
    print("--- (vi) 금지문자 / 인코딩 / 기존 파일 무수정 ---")
    _new = added_lines()
    check('금지문자: 이번에 **추가된** 줄에 0 (기존 줄은 불변)',
          len(_new) > 0 and len(FORBIDDEN_RE.findall('\n'.join(_new))) == 0,
          sorted(set(FORBIDDEN_RE.findall('\n'.join(_new))))[:4])
    _ptr = "\n".join(l for l in ra.splitlines() if '19.104' in l)
    check('금지문자: §19.104 로 표시된 새 줄 0', len(FORBIDDEN_RE.findall(_ptr)) == 0,
          sorted(set(FORBIDDEN_RE.findall(_ptr)))[:4])
    _enc_bad = []
    for _l in logs:
        try:
            _l.encode('cp949')
        except Exception:
            _enc_bad.append(_l)
    check('로그 cp949 인코딩 가능(Windows 콘솔, §19.95)', not _enc_bad, _enc_bad[:2])
    check('RESULTS 4필드 계약 문자열 불변(코드 내 RESULTS: 3종)',
          ra.count('RESULTS:FAIL') >= 3 and 'emit("RESULTS:" + _rc)' in ra)
    check('all_results.csv 열 스키마 불변(hf 마지막, tau/slope 미포함)',
          "'lf1', 'lf2', 'lf3', 'hf']" in ra and 'tau' not in _csv_cols(ra) and 'slope' not in _csv_cols(ra))
    _chg = changed_py_files()
    check('변경된 code/*.py 는 run_abaqus.py (+ harness) 뿐', set(_chg) <= set(
        ['code/run_abaqus.py', 'code/seed_family_stub.py']), _chg)
    check('run_abaqus.py 가 실제로 변경 목록에 있음', 'code/run_abaqus.py' in _chg, _chg)

    print()
    if FAIL:
        print("RESULT: FAIL %d / %d - %s" % (len(FAIL), N[0], ', '.join(FAIL)))
        return 1
    print("RESULT: PASS %d/%d" % (N[0], N[0]))
    return 0


def _csv_cols(src):
    m = re.search(r"_cols = \[(.*?)\]", src, re.S)
    return m.group(0) if m else ''


if __name__ == '__main__':
    sys.exit(main())
