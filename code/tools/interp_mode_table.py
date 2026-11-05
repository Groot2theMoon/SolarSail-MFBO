# -*- coding: utf-8 -*-
"""Interpolate a 1st-order mode table (10,201 nodes) to the 2nd-order grid.

Mid-side values are the average of the two corner values; corner labels are
unchanged.  Run from code/ with normal Python 3.
"""
import argparse
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import aba_grid_mesh
from core.constants import BASE, HEIGHT, SEED_DIV


def parse(path):
    header = []
    tables = {}
    cur = None
    with io.open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            s = line.strip()
            if s.startswith('#'):
                if cur is None:
                    header.append(s)
                continue
            if s.upper().startswith('MODE '):
                cur = int(s.split()[1])
                tables[cur] = {}
                continue
            if cur is None or not s:
                continue
            p = s.split()
            if len(p) >= 2:
                tables[cur][int(p[0])] = float(p[1])
    return header, tables


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', default='seeds/modes_xc0p2.txt')
    ap.add_argument('--out', default='seeds/modes_xc0p2_o2.txt')
    a = ap.parse_args()

    nodes, s4, s3 = aba_grid_mesh.generate_grid(BASE, HEIGHT, SEED_DIV)
    nodes2, s8, s6 = aba_grid_mesh._second_order(nodes, s4, s3)
    n1 = len(nodes)
    n2 = len(nodes2)
    mid = {}
    for e in s8:
        a0, b0, c0, d0, mab, mbc, mcd, mda = e
        mid[mab + 1] = (a0 + 1, b0 + 1)
        mid[mbc + 1] = (b0 + 1, c0 + 1)
        mid[mcd + 1] = (c0 + 1, d0 + 1)
        mid[mda + 1] = (d0 + 1, a0 + 1)
    for e in s6:
        a0, b0, c0, mab, mbc, mca = e
        mid[mab + 1] = (a0 + 1, b0 + 1)
        mid[mbc + 1] = (b0 + 1, c0 + 1)
        mid[mca + 1] = (c0 + 1, a0 + 1)

    header, tables = parse(a.src)
    missing = [m for m in (1, 2, 3, 4) if m not in tables]
    if missing:
        raise SystemExit('missing modes %s' % missing)
    if len(tables[1]) != n1:
        raise SystemExit('source count %d != first-order nodes %d' % (len(tables[1]), n1))
    if len(mid) != n2 - n1:
        raise SystemExit('mid count %d != added nodes %d' % (len(mid), n2 - n1))

    with io.open(a.out, 'w', encoding='utf-8') as f:
        for h in header:
            f.write(h + '\n')
        f.write('# interp=2nd-order-midpoint from %s\n' % a.src)
        for m in (1, 2, 3, 4):
            f.write('MODE %d\n' % m)
            src = tables[m]
            for lab in range(1, n2 + 1):
                if lab <= n1:
                    v = src[lab]
                else:
                    i, j = mid[lab]
                    v = 0.5 * (src[i] + src[j])
                f.write('%d %.10e\n' % (lab, v))
    print('wrote %s: modes=%d nodes/mode=%d' % (a.out, 4, n2))


if __name__ == '__main__':
    main()
