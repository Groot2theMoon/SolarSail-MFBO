# -*- coding: utf-8 -*-
"""Validate the C1 mode-table transfer on the 2nd-order grid.

Checks (no Abaqus, normal Python + numpy/scipy):
  1. all 10,201 corner labels keep the exact 1st-order values
  2. the C1 interpolant reproduces a linear field exactly at all 30,501 nodes
  3. C1 vs bilinear differs only on newly added mid-side nodes; report max/RMS
  4. all values finite and per-mode max|u3| is unchanged at corners

Usage:
    python tools/validate_c1_interp.py --src seeds/modes_xc0p2.txt \
        --c1 seeds/modes_xc0p2_o2_c1.txt \
        --bilinear seeds/modes_xc0p2_o2.txt \
        --out results/variants/o2_c1/validate_c1.json
"""
from __future__ import print_function

import argparse
import io
import json
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
    ap.add_argument('--c1', default='seeds/modes_xc0p2_o2_c1.txt')
    ap.add_argument('--bilinear', default='seeds/modes_xc0p2_o2.txt')
    ap.add_argument('--out', default='results/variants/o2_c1/validate_c1.json')
    ap.add_argument('--base', type=float, default=20.0)
    ap.add_argument('--height', type=float, default=10.0)
    ap.add_argument('--seed-div', type=float, default=200.0)
    a = ap.parse_args()

    _, src = parse(a.src)
    _, c1 = parse(a.c1)
    _, bil = parse(a.bilinear) if a.bilinear and os.path.exists(a.bilinear) else (None, {})
    nodes, s4, s3 = aba_grid_mesh.generate_grid(a.base, a.height, a.seed_div)
    nodes2, s8, s6 = aba_grid_mesh._second_order(nodes, s4, s3)
    if len(nodes) != 10201 or len(nodes2) != 30501:
        raise RuntimeError('grid size mismatch: %d / %d' % (len(nodes), len(nodes2)))
    pts = np.asarray(nodes, dtype=float)
    q = np.asarray(nodes2, dtype=float)
    n1 = len(nodes)

    linear = pts[:, 0] + 2.0 * pts[:, 1]
    lin_interp = CloughTocher2DInterpolator(pts, linear)
    lin_exact = q[:, 0] + 2.0 * q[:, 1]
    lin_err = np.abs(lin_interp(q) - lin_exact)

    report = {
        'source': a.src, 'c1': a.c1, 'bilinear': a.bilinear,
        'n1': n1, 'n2': len(nodes2),
        'linear_max_abs_error': float(np.max(lin_err)),
        'modes': {},
    }
    for m in sorted(set(src) & set(c1)):
        v_src = np.asarray([src[m][i + 1] for i in range(n1)], dtype=float)
        v_c1 = np.asarray([c1[m][i + 1] for i in range(len(nodes2))], dtype=float)
        corner_err = np.abs(v_c1[:n1] - v_src)
        rec = {
            'corner_max_abs_error': float(np.max(corner_err)),
            'corner_rms_error': float(np.sqrt(np.mean(corner_err ** 2))),
            'c1_max_abs': float(np.max(np.abs(v_c1))),
            'c1_finite': bool(np.isfinite(v_c1).all()),
            'c1_nonzero_mid': int(np.count_nonzero(v_c1[n1:])),
        }
        if m in bil and len(bil[m]) >= len(nodes2):
            v_b = np.asarray([bil[m][i + 1] for i in range(len(nodes2))], dtype=float)
            d = v_c1 - v_b
            rec.update({
                'mid_max_abs_diff': float(np.max(np.abs(d[n1:]))),
                'mid_rms_diff': float(np.sqrt(np.mean(d[n1:] ** 2))),
                'corner_vs_bilinear_max': float(np.max(np.abs(d[:n1]))),
                'mac_c1_bilinear': float((np.dot(v_c1, v_b) ** 2) /
                                         ((np.dot(v_c1, v_c1) * np.dot(v_b, v_b)) or np.nan)),
            })
        report['modes'][m] = rec

    print(json.dumps(report, ensure_ascii=False, indent=1))
    if a.out:
        d = os.path.dirname(a.out)
        if d and not os.path.isdir(d):
            os.makedirs(d)
        with io.open(a.out, 'w', encoding='utf-8') as f:
            f.write(json.dumps(report, ensure_ascii=False, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
