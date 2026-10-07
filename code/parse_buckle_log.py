#!/usr/bin/env python3
"""Digest Abaqus buckling logs (`.msg` or pasted console text) into one comparison table.

Usage:
    python parse_buckle_log.py <file-or-dir> [<file-or-dir> ...]
    python parse_buckle_log.py --attachments [N]      # newest N files in ~/.hermes/attachments

Why: a buckling round is judged from `numEigen`/`vectors`, the two negative-eigenvalue counts, the per-iteration
positive/negative split, the `CONVERGED` series and the wallclock. Eyeballing a 100 KB log gets that wrong (an
all-negative iteration-1 ladder reads like progress), and the cross-round table is what shows whether a knob
actually moved. Stdlib only; ASCII labels on purpose, because these logs come from a CP949 Windows shell.
"""
from __future__ import print_function

import glob
import os
import re
import sys

NUMVAL = re.compile(r'[-+]?\d\.\d{6}e[-+]\d{2}')
ERR = re.compile(r'\*\*\*ERROR[^\n]*')
ITER_SPLIT = re.compile(r'\bITERATION\s+(\d+)\s*\n')


def read_text(path):
    """Decode a log that may be UTF-8, CP949 or latin-1 — never report it unreadable."""
    with open(path, 'rb') as fh:
        data = fh.read()
    for enc in ('utf-8', 'cp949', 'latin-1'):
        try:
            return data.decode(enc)
        except UnicodeDecodeError:
            continue
    return data.decode('utf-8', 'replace')


def _int(pat, text):
    m = re.search(pat, text, re.M)
    return int(m.group(1)) if m else None


def _float(pat, text):
    m = re.search(pat, text, re.M)
    try:
        return float(m.group(1)) if m else None
    except ValueError:
        return None


def parse(text):
    d = {
        'numEigen': _int(r'NUMBER OF EIGENVALUES\s+(\d+)', text),
        'vectors': _int(r'NUMBER OF VECTORS IN ITERATION\s+(\d+)', text),
        'maxIter': _int(r'MAXIMUM NUMBER OF ITERATIONS\s+(\d+)', text),
        'eqs': _int(r'NUMBER OF EQUATIONS:\s+(\d+)', text),
        'sysNeg': _int(r'SYSTEM MATRIX HAS (\d+) NEGATIVE', text),
        # [2026-10-07] subspace 실제 차원 — 창(window) 추천의 상한이다.
        #   실측 문구: "THE NUMBER OF VECTORS IN SUBSPACE IS REDUCED TO 14 DUE TO LINEAR DEPENDENCY"
        'reduced_to': _int(r'VECTORS IN SUBSPACE IS REDUCED TO\s+(\d+)', text),
        'diffNeg': _int(r'DIFFERENTIAL MATRIX HAS (\d+) NEGATIVE DIAGONAL', text),
        'wall': _float(r'WALLCLOCK TIME \(SEC\)\s*=\s*([\d.eE+]+)', text),
        'cpu': _float(r'TOTAL CPU TIME \(SEC\)\s*=\s*([\d.eE+]+)', text),
        'mem': _float(r'MEMORY PEAK \(GB\)\s*=\s*([\d.eE+]+)', text),
        'warnAnalysis': _int(r'(\d+)\s+WARNING MESSAGES DURING ANALYSIS', text),
        'nErr': _int(r'(\d+)\s+ERROR MESSAGES', text),
        'iters': [],
        'conv': [int(v) for v in re.findall(r'CONVERGED UPTO THIS ITERATION =\s*(\d+)', text)],
        'errors': ERR.findall(text),
    }
    parts = ITER_SPLIT.split(text)
    for i in range(1, len(parts) - 1, 2):
        body = parts[i + 1]
        cut = body.find('NUMBER OF EIGENVALUES CONVERGED')
        if cut >= 0:
            body = body[:cut]
        vals = [float(v) for v in NUMVAL.findall(body)]
        if not vals:
            continue
        pos = sorted(v for v in vals if v > 0.0)
        d['iters'].append({
            'n': int(parts[i]), 'vals': len(vals), 'pos': len(pos),
            'min': min(vals), 'max': max(vals), 'posFirst': pos[:4],
        })
    return d


MODE_TABLE_HDR = re.compile(r'MODE\s+NO\s+.*EIGEN', re.I)


def dat_lambdas(text, limit=12):
    """`.dat` 의 MODE NO / EIGENVALUE 표에서 lambda 를 읽는다 (순수 함수 — 단위검증 가능).

    왜 표에서만 읽는가: 실패 런의 `ITERATION n` 목록은 레일리 몫 스냅샷이지 스펙트럼이 아니다.
    요청 수만큼 `CONVERGED` 한 런의 표만 lambda 로 인용할 수 있다. 표가 없으면 빈 리스트
    (= 모드를 못 냈다)이고, 그것이 정상적인 '실패' 표현이다.
    """
    lines = text.splitlines()
    hdr = [i for i, l in enumerate(lines) if MODE_TABLE_HDR.search(l)]
    if not hdr:
        return []
    out = []
    for l in lines[hdr[-1] + 1:]:
        p = l.split()
        if not p:
            if out:
                break
            continue
        if not (p[0].isdigit() and len(p) >= 2):
            break
        try:
            out.append(float(p[-1]))
        except ValueError:
            break
        if len(out) >= limit:
            break
    return out


def _window_hint(d):
    """창(window) 추천 — 두 조건의 교집합이다.

      (a) N_EIG_BUCKLE > base state 의 SYSTEM 음수 고유값 개수
      (b) N_EIG_BUCKLE <= 실제 subspace 차원 ("VECTORS IN SUBSPACE IS REDUCED TO n" 의 n)

    (a)는 '음수 모드 뒤에 양수 모드가 들어오게' 하는 조건이고 (b)는 '원리적으로 찾을 수 있는 수'의
    상한이다. 둘 다 만족해야 하며, 교집합이 비면 **창으로 못 푸는 실패**다 — 그때는 창을 키우지 말고
    base state 를 의심한다(ITERATION 1 이 전부 음수면 그 점이 좌굴 후다).
    """
    neg, red = d.get('sysNeg'), d.get('reduced_to')
    if neg is None:
        return
    if red is None:
        print('  window hint  : N_EIG_BUCKLE > %d  (SYSTEM 음수 %d개)' % (neg, neg))
        print('                 로그에 "REDUCED TO n" 이 있으면 n 이 상한이다(그때는 n 이하로).')
    elif neg + 1 <= red:
        print('  window hint  : N_EIG_BUCKLE in (%d, %d]   (음수 %d < 요청 <= subspace %d)'
              % (neg, red, neg, red))
    else:
        print('  window hint  : 공집합 — 음수 %d >= subspace %d  ==> 창으로 못 푼다.' % (neg, red))
        print('                 base state 가 좌굴 후이거나(ITERATION 1 전부 음수) 다른 원인이다.')


def report(path, d):
    print('=' * 78)
    try:
        size = os.path.getsize(path)
    except OSError:
        size = 0
    print('FILE  %s   (%d bytes)' % (path, size))
    if d['numEigen'] is None and not d['iters']:
        print('  no buckling echo found in this file')
        return
    print('  request      : numEigen=%s  vectors=%s  maxIterations=%s  equations=%s'
          % (d['numEigen'], d['vectors'], d['maxIter'], d['eqs']))
    print('  base state   : SYSTEM MATRIX HAS %s NEGATIVE EIGENVALUES   |   DIFFERENTIAL MATRIX HAS %s NEGATIVE DIAGONALS'
          % (d['sysNeg'], d['diffNeg']))
    print('  iterations   : %d listed   CONVERGED series %s' % (len(d['iters']), d['conv']))
    _window_hint(d)
    for it in d['iters']:
        note = ('positives: ' + ', '.join('%+.4e' % v for v in it['posFirst'])) if it['posFirst'] else ''
        print('    ITER %-4d n=%-4d pos=%-3d neg=%-4d  min=%+.4e  max=%+.4e  %s'
              % (it['n'], it['vals'], it['pos'], it['vals'] - it['pos'], it['min'], it['max'], note))
    print('  wallclock    : %s s  cpu=%s s  mem=%s GB  analysis warnings=%s  errors=%s'
          % (d['wall'], d['cpu'], d['mem'], d['warnAnalysis'], d['nErr']))
    for e in d['errors'][:4]:
        print('  ERROR        : %s' % e.rstrip())


def table(rows):
    if len(rows) < 2:
        return
    print('=' * 78)
    print('COMPARISON')
    hdr = ('%-32s %-13s %-7s %-8s %-6s %-7s %-8s %s'
           % ('file', 'num/vec', 'sysNeg', 'diffNeg', 'iters', 'posSum', 'convMax', 'wall'))
    print(hdr)
    print('-' * len(hdr))
    for path, d in rows:
        pos_sum = sum(it['pos'] for it in d['iters']) if d['iters'] else '-'
        conv_max = max(d['conv']) if d['conv'] else '-'
        print('%-32s %-13s %-7s %-8s %-6s %-7s %-8s %s'
              % (os.path.basename(path)[-32:],
                 '%s/%s' % (d['numEigen'], d['vectors']), d['sysNeg'], d['diffNeg'],
                 len(d['iters']) or '-', pos_sum, conv_max, d['wall']))
    print('Read it as: a non-zero convMax is the only pass line; posSum==0 with a long CONVERGED-zero series means the'
          ' perturbation pattern pushes the wrong way or the base state is unstable (check [R-13] in the console).')


def expand(args):
    if args and args[0] == '--attachments':
        n = int(args[1]) if len(args) > 1 else 5
        home = os.path.expanduser('~/.hermes/attachments')
        files = sorted(glob.glob(os.path.join(home, '*')), key=os.path.getmtime)
        return files[-n:]
    out = []
    for a in args:
        if os.path.isdir(a):
            out += sorted(glob.glob(os.path.join(a, '*')))
        else:
            hit = sorted(glob.glob(a))
            out += hit or [a]
    return out


def main(argv):
    args = argv[1:]
    if not args:
        print(__doc__)
        return 2
    rows = []
    for path in expand(args):
        if not os.path.isfile(path):
            print('SKIP %s (not a file)' % path)
            continue
        try:
            if os.path.getsize(path) > 8 * 1024 * 1024:
                print('SKIP %s (too large)' % path)
                continue
            d = parse(read_text(path))
        except Exception as exc:  # one odd file must not kill the digest
            print('SKIP %s (%s)' % (path, exc))
            continue
        report(path, d)
        rows.append((path, d))
    table(rows)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
