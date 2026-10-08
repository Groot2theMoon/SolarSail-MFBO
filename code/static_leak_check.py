#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""정적 누출 검사 — 실행 없이 잡히는 NameError 두 종류 (2026-10-08)

왜 필요한가
    좌굴 스크립트(run_abaqus_buckle.py)는 스텁 하네스로 **실행** 검증이 된다. 그런데
    나머지 세 스크립트(run_abaqus.py / run_abaqus_mode.py / run_abaqus_cable.py)에는
    그런 검사가 없다. 그래서 격자 블록을 복사할 때 딸려온 이름이나 빠진 임포트가 조용히
    통과해 실제 Abaqus 런에서 죽는다 — 실제로 HF 가 그렇게 두 번 죽었다(§19.33, §19.34):
        NameError: name 'TAG' is not defined          <- 복사가 끌고 온 좌굴 전용 상수
        NameError: name 'io'  is not defined          <- 임포트 누락(덱 검사 코드 3곳)

무엇을 잡나 (Abaqus 없이, 표준 라이브러리만)
    R1  모듈처럼 쓰였는데(`io.open(...)`) 임포트가 없는 **표준 라이브러리 이름**.
        `from abaqus import *` 가 채워 주지 않는 이름들이라 반드시 잡아야 한다.
    R2  다른 스크립트에서 **모듈 수준으로 정의됐는데, 이 파일에서는 정의 없이 쓰이는** 이름.
        예: 좌굴의 `TAG = '[run_abaqus_buckle]'` 를 복사한 세 스크립트.
        양쪽 모두에서 '정의되지 않은' 이름(예: abaqusConstants 가 주는 SET)은 잡지 않는다.
        그래야 오탐이 0 이다.

사용
    python3 static_leak_check.py              # code/ 안의 대상 스크립트 전부
    python3 static_leak_check.py selftest     # 순수 로직 자체 점검 (파일 없이)

종료 코드: 0 이상 없음 / 1 누출 발견 / 2 사용법 / 3 selftest 실패

한계 (정직하게)
    `from abaqus import *` 가 실제로 제공하는 이름은 정적으로 알 수 없다. R2 는
    '어떤 스크립트가 **스스로** 정의한 이름'만 보므로 그 애매함을 피한다. 따라서
    별칭/동적 접근(getattr)으로 숨긴 누출은 잡지 못한다.
"""
import io
import os
import re
import sys

TARGETS = ('run_abaqus.py', 'run_abaqus_mode.py', 'run_abaqus_cable.py', 'run_abaqus_buckle.py')

# R1 대상: 흔한 표준 라이브러리 모듈 이름. `from abaqus import *` 와 겹치지 않는다.
MODS = ('io', 'os', 'sys', 're', 'math', 'time', 'json', 'csv', 'shutil', 'glob', 'subprocess',
        'itertools', 'collections', 'tempfile', 'datetime', 'platform', 'struct', 'hashlib',
        'random', 'copy', 'warnings')

IMPORT_RE = re.compile(r'^\s*import\s+(.+)$', re.M)
FROM_RE = re.compile(r'^\s*from\s+[\w.]+\s+import\s+(.+)$', re.M)
MODULE_ASSIGN_RE = re.compile(r'^([A-Za-z_][A-Za-z0-9_]*)\s*=', re.M)
ANY_ASSIGN_RE = re.compile(r'^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=', re.M)
FOR_RE = re.compile(r'\bfor\s+([A-Za-z_][A-Za-z0-9_]*)')


def strip_comments(text):
    """줄 끝 주석을 지운다. `import aba_grid_mesh  # noqa: E402` 같은 줄이 임포트
    인식에서 빠지던 실수를 막는다(1차 스캔이 오탐 8개를 낸 원인이었다)."""
    return '\n'.join(l.split('#')[0] for l in text.split('\n'))


TRIPLE_RE = re.compile(r'""".*?"""|\'\'\'.*?\'\'\'', re.S)
STRING_RE = re.compile(r'"[^"\n]*"|\'[^\'\n]*\'')


def strip_strings(code):
    """문자열 리터럴을 지운다. docstring/도움말에 이름이 적혀 있으면 이름 '사용'이 아니다 —
    1차 스캔의 오탐 7건 중 다수가 이 때문이었다(예: mode.py 도움말 속 IMPERFECTION_MODE)."""
    return STRING_RE.sub('', TRIPLE_RE.sub('', code))


def used_modules(code):
    """`X.` 형태로 쓰인 이름 중 표준 라이브러리 목록에 있는 것.
    앞에 점이 붙은(`obj.datetime.`) 것은 제외한다 — 그건 속성 접근이지 모듈이 아니다."""
    return set(re.findall(r'(?<![\w.])([a-z_][a-z0-9_]*)\.', code)) & set(MODS)


def imported_names(code):
    out = set()
    for m in IMPORT_RE.finditer(code):
        for part in m.group(1).split(','):
            out.add(part.strip().split(' as ')[-1].strip())
    for m in FROM_RE.finditer(code):
        for part in m.group(1).split(','):
            part = part.strip()
            if part and part != '*':
                out.add(part.split(' as ')[-1].strip())
    return out


def module_level_names(code):
    return set(MODULE_ASSIGN_RE.findall(code))


DEF_RE = re.compile(r'^\s*def\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(([^)]*)\)', re.M | re.S)
AS_RE = re.compile(r'\bas\s+([A-Za-z_][A-Za-z0-9_]*)')
LHS_ID_RE = re.compile(r'(?<![\w.])([A-Za-z_][A-Za-z0-9_]*)\s*(?:,|=(?!=))')


def locally_bound_names(code):
    """이 파일 안에서 '묶이는' 모든 이름 — 오탐 억제용으로 넓게 잡는다.
    모듈수준/함수내 대입, for 대상, def 인자, with/except ... as, 튜플 대입 좌변까지."""
    out = set(ANY_ASSIGN_RE.findall(code)) | set(FOR_RE.findall(code))
    for m in DEF_RE.finditer(code):
        for a in m.group(2).split(','):
            a = a.strip().split('=')[0].strip().lstrip('*')
            if re.match(r'^[A-Za-z_][A-Za-z0-9_]*$', a):
                out.add(a)
    out |= set(AS_RE.findall(code))
    # 튜플 대입 좌변: 'a, b = ...' / 'x, y, z = ...'
    for m in re.finditer(r'^\s*([A-Za-z_][\w]*(?:\s*,\s*[A-Za-z_][\w]*)+)\s*=', code, re.M):
        for a in m.group(1).split(','):
            out.add(a.strip())
    return out


ELEM_NEEDING_CALL = re.compile(r'(\.Stress|\.SectionAssignment|InitialCondition)\s*\(')


def grid_node_set_misuse(codes):
    """격자 모드에서 'All'(절점집합)을 **요소가 필요한 호출**에 넘기면 잡는다.

    격자(USE_GRID_MESH=True)에서 파트의 `All` 은 **절점**집합이다. 초기응력(Stress)이나
    SectionAssignment 처럼 요소/면 집합이 필요한 곳에 `All` 을 주면 Abaqus 입력처리기가
    'not an element set' 으로 죽는다. 실제로 HF 가 이걸로 죽었다 — 초기응력을 만드는
    Stress 호출이 3곳인데 격자 전환에서 1곳만 _ALL_ELEM 으로 바뀌고 나머지 2곳이
    'All' 로 남아, 마지막 재생성이 절점집합 위에 초기응력을 걸었다(2026-10-08).
    자유메쉬에서는 All 이 면집합이라 정상이므로 격자일 때만 검사한다.
    (DisplacementBC 는 절점집합이 맞으므로 대상이 아니다.)
    """
    out = []
    for p in sorted(codes):
        code = codes[p]
        if not re.search(r'^USE_GRID_MESH\s*=\s*True', code, re.M):
            continue
        for m in ELEM_NEEDING_CALL.finditer(code):
            seg = code[m.start():m.start() + 700]
            if re.search(r"sets\['All'\]", seg):
                ln = code[:m.start()].count('\n') + 1
                out.append((p, 'R3', '%s @%d 에 sets[All]' % (m.group(1), ln)))
    return out


def scan_files(paths):
    """반환: (problems, notes). problems = [(path, rule, name), ...]"""
    codes = {}
    for p in paths:
        with io.open(p, encoding='utf-8', errors='replace') as f:
            codes[p] = strip_comments(f.read())

    defined_somewhere = set()
    for p in paths:
        defined_somewhere |= module_level_names(codes[p])

    problems = []
    for p in paths:
        code = codes[p]
        imp = imported_names(code)
        loc = locally_bound_names(code)
        for name in sorted(used_modules(code)):
            if name not in imp and name not in loc:
                problems.append((p, 'R1', name))
        for name in sorted(defined_somewhere):
            if name in imp or name in loc:
                continue
            if re.search(r'(?<![\w.])%s(?![\w.])' % re.escape(name), strip_strings(code)):
                problems.append((p, 'R2', name))
    problems.extend(grid_node_set_misuse(codes))
    return problems, codes


def selftest():
    """순수 로직만 검증한다(파일 불필요)."""
    ok = True
    code = 'import os\nio.open("x")\nos.path\n_dt.datetime.now()\n'
    ok &= ('io' in used_modules(code))
    ok &= ('os' not in used_modules(code) or 'os' in imported_names(code))
    ok &= ('datetime' not in used_modules(code))          # `_dt.datetime.` 는 제외돼야 한다
    ok &= (imported_names('import os, sys as s\nfrom x import a, b\n') == {'os', 's', 'a', 'b'})
    ok &= (module_level_names('TAG = 1\n_x = 2\n') == {'TAG', '_x'})
    ok &= ('io' in used_modules(strip_comments('io.open(1)  # io.open(2)\n')))
    ok &= (imported_names(strip_comments('import aba_grid_mesh  # noqa: E402\n')) == {'aba_grid_mesh'})
    ok &= (used_modules(strip_strings('"IMPERFECTION_MODE"\n')) == set())      # 문자열 언급은 제외
    ok &= ('p' in locally_bound_names('def f(p, q=1):\n    return p\n'))
    ok &= ('a' in locally_bound_names('a, b = 1, 2\n'))
    ok &= ('e' in locally_bound_names('try:\n    pass\nexcept X as e:\n    pass\n'))
    _grid_bad = {'x.py': 'USE_GRID_MESH = True\nmy_model.Stress(\n    region=inst.sets[\'All\'],\n)\n'}
    _grid_ok = {'x.py': 'USE_GRID_MESH = True\nmy_model.Stress(\n    region=inst.sets[_ALL_ELEM],\n)\n'}
    _free = {'x.py': 'USE_GRID_MESH = False\nmy_model.Stress(\n    region=inst.sets[\'All\'],\n)\n'}
    ok &= (len(grid_node_set_misuse(_grid_bad)) == 1)      # 잡아야 한다
    ok &= (len(grid_node_set_misuse(_grid_ok)) == 0)
    ok &= (len(grid_node_set_misuse(_free)) == 0)          # 자유메쉬는 정상
    print('[selftest] static_leak_check %s' % ('OK' if ok else 'FAIL'))
    return 0 if ok else 3


def main(argv):
    here = os.path.dirname(os.path.abspath(__file__))
    if len(argv) > 1 and argv[1] == 'selftest':
        return selftest()
    if len(argv) > 1:
        paths = argv[1:]
    else:
        paths = [os.path.join(here, t) for t in TARGETS if os.path.exists(os.path.join(here, t))]
    if not paths:
        print('usage: python3 static_leak_check.py [<script.py> ...] | selftest')
        return 2
    problems, _ = scan_files(paths)
    if not problems:
        print('[leak] 이상 없음 (%d 파일 검사)' % len(paths))
        return 0
    print('[leak] 누출 %d건 — 실행 전에 고치세요:' % len(problems))
    for p, rule, name in problems:
        why = ('모듈처럼 쓰였는데 임포트 없음' if rule == 'R1'
               else '다른 스크립트의 모듈수준 이름을 정의 없이 사용')
        print('  %-24s %s  %-20s  %s' % (os.path.basename(p), rule, name, why))
    print('\n  R1 = 빠진 임포트 / R2 = 복사가 끌고 온 이름 / R3 = 격자에서 절점집합을 요소자리에 사용')
    return 1


if __name__ == '__main__':
    sys.exit(main(sys.argv))
