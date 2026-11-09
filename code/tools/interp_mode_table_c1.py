# -*- coding: utf-8 -*-
"""C1 control: map a 1st-order mode table to the 2nd-order grid with Clough-Tocher.

The production path (`interp_mode_table.py`) uses the exact Q1 embedding: corners
are copied and mid-side nodes get 0.5*(a+b).  That field is C0 at element edges.
This tool instead evaluates a globally C1 piecewise-cubic surface
(scipy.interpolate.CloughTocher2DInterpolator) through the same corner values.
Only the seed field changes; the 2nd-order HF model is identical.

Run with normal Python 3 + numpy/scipy (not Abaqus):
    python tools/interp_mode_table_c1.py --src seeds/modes_xc0p2.txt \
        --out seeds/modes_xc0p2_o2_c1.txt
"""
from __future__ import print_function

import argparse
import io
import os
import sys

import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator

_HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import aba_grid_mesh  # noqa: E402


def parse(path):
    header, tables, cur = [], {}, None
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
    ap.add_argument('--out', default='seeds/modes_xc0p2_o2_c1.txt')
    ap.add_argument('--base', type=float, default=20.0)
    ap.add_argument('--height', type=float, default=10.0)
    ap.add_argument('--seed-div', type=float, default=200.0)
    a = ap.parse_args()

    header, tables = parse(a.src)
    nodes, s4, s3 = aba_grid_mesh.generate_grid(a.base, a.height, a.seed_div)
    nodes2, s8, s6 = aba_grid_mesh._second_order(nodes, s4, s3)
    if len(nodes) != 10201:
        raise RuntimeError('expected 10201 1st nodes, got %d' % len(nodes))
    if len(nodes2) != 30501:
        raise RuntimeError('expected 30501 2nd nodes, got %d' % len(nodes2))

    pts = np.asarray(nodes, dtype=np.float64)
    q = np.asarray(nodes2, dtype=np.float64)
    out = []
    out.append('# C1 control: CloughTocher2DInterpolator through the 1st-order nodal values')
    out.extend(header)
    for m in sorted(tables):
        vals = np.asarray([tables[m][i + 1] for i in range(len(nodes))], dtype=np.float64)
        interp = CloughTocher2DInterpolator(pts, vals)
        z = interp(q)
        bad = ~np.isfinite(z)
        if bad.any():
            raise RuntimeError('C1 interpolation produced %d non-finite 2nd-grid values' % bad.sum())
        out.append('MODE %d' % m)
        for i, v in enumerate(z):
            out.append('%d %.16e' % (i + 1, v))
    with io.open(a.out, 'w', encoding='utf-8') as f:
        f.write('\n'.join(out) + '\n')
    print('wrote %s: modes=%d nodes/mode=%d' % (a.out, len(tables), len(nodes2)))


if __name__ == '__main__':
    main()
