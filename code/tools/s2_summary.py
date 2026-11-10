# -*- coding: utf-8 -*-
"""Summarize S2 seed runs (hf/cp/torque; seed ICC; sign consistency)."""
from __future__ import print_function

import argparse
import csv
import io
import json
import math
import os
import sys

import numpy as np


def val(path, key):
    if not os.path.exists(path):
        return None
    with io.open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith(key + '='):
                try:
                    return float(line.strip().split('=', 1)[1])
                except Exception:
                    return None
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', default='results/variants/uq_s2/s2_runs.csv')
    ap.add_argument('--out', default='results/variants/uq_s2/s2_summary.json')
    a = ap.parse_args()
    rows = list(csv.DictReader(io.open(a.csv, encoding='utf-8', errors='replace')))
    groups = {}
    for r in rows:
        try:
            x = float(r['x_c']); d = float(r['d_c']); seed = int(r['seed'])
        except Exception:
            continue
        key = (r['kind'], x, d)
        hf = val(r['results'], 'hf') if r.get('results') else None
        cp = val(r['results'], 'cp_offset_m') if r.get('results') else None
        tau = val(r['results'], 'torque_nm') if r.get('results') else None
        tx = val(r['results'], 'tau_x_nm') if r.get('results') else None
        ty = val(r['results'], 'tau_y_nm') if r.get('results') else None
        # run_abaqus currently writes cp components, not tau components.
        if tx is None:
            tx = val(r['results'], 'cp_x_m') if r.get('results') else None
        if ty is None:
            ty = val(r['results'], 'cp_y_m') if r.get('results') else None
        groups.setdefault(key, []).append({'seed': seed, 'hf': hf, 'cp': cp,
                                           'tau': tau, 'tx': tx, 'ty': ty})
    out = {'groups': {}}
    for key, vals in sorted(groups.items()):
        vals.sort(key=lambda z: z['seed'])
        rec = {'kind': key[0], 'x_c': key[1], 'd_c': key[2], 'n': len(vals)}
        for field in ('hf', 'cp', 'tau', 'tx', 'ty'):
            arr = np.asarray([v[field] for v in vals if v[field] is not None], dtype=float)
            if len(arr):
                rec[field] = {'n': int(len(arr)), 'mean': float(arr.mean()),
                              'std': float(arr.std(ddof=1)) if len(arr) > 1 else 0.0,
                              'min': float(arr.min()), 'max': float(arr.max()),
                              'rel_std': (float(arr.std(ddof=1) / arr.mean())
                                          if len(arr) > 1 and arr.mean() != 0 else None)}
        s = []
        for f in ('tx', 'ty'):
            signs = [math.copysign(1.0, v[f]) for v in vals if v[f] is not None and abs(v[f]) > 1e-12]
            if len(signs) >= 2:
                s.append((f, max(signs.count(1), signs.count(-1)) / float(len(signs))))
        rec['sign_consistency'] = dict(s)
        out['groups']['%s_x%g_d%g' % (key[0], key[1], key[2])] = rec
    with io.open(a.out, 'w', encoding='utf-8') as f:
        f.write(json.dumps(out, ensure_ascii=False, indent=1))
    print(json.dumps(out, ensure_ascii=False, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
