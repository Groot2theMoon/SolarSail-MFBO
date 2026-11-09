# -*- coding: utf-8 -*-
"""Interpolate a 0.1 m mode table to a finer first-order triangle grid.

Usage (from code/):
    python tools/interp_fine_mode.py --src seeds/modes_xc0p2.txt --div 400 --out seeds/modes_xc0p2_fine400.txt

Requires numpy/scipy (normal Python 3, not Abaqus).
"""
import argparse
import io
import os
import sys

import numpy as np
from scipy.interpolate import griddata

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import aba_grid_mesh

BASE, HEIGHT = 20.0, 10.0


def read_table(path):
    header, modes, cur = [], {}, None
    with io.open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith('#'):
                if cur is None:
                    header.append(s)
                continue
            if s.upper().startswith('MODE '):
                cur = int(s.split()[1]); modes[cur] = {}
                continue
            if cur is None:
                continue
            p = s.split()
            if len(p) >= 2:
                modes[cur][int(p[0])] = float(p[1])
    return header, modes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', default='seeds/modes_xc0p2.txt')
    ap.add_argument('--div', type=int, default=400)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    src_nodes, _, _ = aba_grid_mesh.generate_grid(BASE, HEIGHT, 200.0)
    dst_nodes, _, _ = aba_grid_mesh.generate_grid(BASE, HEIGHT, float(a.div))
    src_xy = np.asarray(src_nodes, dtype=np.float64)
    dst_xy = np.asarray(dst_nodes, dtype=np.float64)
    header, modes = read_table(a.src)
    if sorted(modes) != [1, 2, 3, 4]:
        raise SystemExit('source must have modes 1..4')
    if len(modes[1]) != len(src_nodes):
        raise SystemExit('source nodes %d != 0.1 grid %d' % (len(modes[1]), len(src_nodes)))

    with io.open(a.out, 'w', encoding='utf-8') as f:
        for h in header:
            f.write(h + '\n')
        f.write('# interp_fine: 0.1m -> 0.05m (div=%d), linear+nearest boundary\n' % a.div)
        for m in (1, 2, 3, 4):
            vals = np.array([modes[m][i + 1] for i in range(len(src_nodes))], dtype=np.float64)
            out = griddata(src_xy, vals, dst_xy, method='linear')
            bad = ~np.isfinite(out)
            if np.any(bad):
                near = griddata(src_xy, vals, dst_xy[bad], method='nearest')
                out[bad] = near
            f.write('MODE %d\n' % m)
            for i, v in enumerate(out):
                f.write('%d %.10e\n' % (i + 1, v))
    print('wrote %s: %d nodes, div=%d' % (a.out, len(dst_nodes), a.div))


if __name__ == '__main__':
    main()
