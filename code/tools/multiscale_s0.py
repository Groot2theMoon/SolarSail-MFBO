# -*- coding: utf-8 -*-
"""S0: HF ODB -> theta^2, eps_m (LE), hf consistency.

Usage (from code/):
    abaqus python tools/multiscale_s0.py --all --out results/multiscale/s0_metrics.csv

theta^2 = area-weighted mean of angle(n,z)^2  (same definition as slope_stats)
eps_m   = area-weighted mean(LE11+LE22)       (LE because E is not in these ODBs)
hf_code = points.csv hf
check   hf_code ~= (1+eps_m) * k * theta^2, k=1.461716
"""
import argparse
import csv
import glob
import math
import os

import numpy as np
from odbAccess import openOdb

K = 1.461716


def _tag(path):
    name = os.path.basename(path)
    if name.startswith('HF_'):
        name = name[3:]
    return name.replace('_u3free.odb', '').replace('.odb', '')


def _step(odb):
    for nm in ('Step-GlobalTension', 'Step-Postbuckle'):
        if nm in odb.steps:
            return odb.steps[nm]
    return odb.steps[list(odb.steps.keys())[-1]]


def _strains(frame, inst):
    key = 'LE' if 'LE' in frame.fieldOutputs else ('E' if 'E' in frame.fieldOutputs else None)
    if key is None:
        return {}, None
    sub = frame.fieldOutputs[key].getSubset(region=inst)
    acc = {}
    for v in sub.values:
        lab = v.elementLabel
        if lab is None:
            continue
        d = v.data
        val = float(d[0] + d[1])  # LE11 + LE22 (small-strain: areal)
        acc.setdefault(lab, []).append(val)
    return {k: sum(v) / len(v) for k, v in acc.items()}, key


def one(path, points):
    odb = openOdb(path=path, readOnly=True)
    try:
        step = _step(odb)
        fr = step.frames[-1]
        inst = odb.rootAssembly.instances['MEMBRANE-1']
        coord = {}
        for v in fr.fieldOutputs['COORD'].getSubset(region=inst).values:
            if v.nodeLabel is not None:
                coord[v.nodeLabel] = v.data
        eps, ekey = _strains(fr, inst)
        th2 = []
        area = []
        epsv = []
        for e in inst.elements:
            conn = list(e.connectivity)
            if len(conn) not in (3, 4):
                continue
            if any(n not in coord for n in conn):
                continue
            cur = np.array([coord[n] for n in conn])
            if len(conn) == 4:
                cr = np.cross(cur[2] - cur[0], cur[3] - cur[1])
            else:
                cr = np.cross(cur[1] - cur[0], cur[2] - cur[0])
            nrm = float(np.linalg.norm(cr))
            if nrm <= 1e-14:
                continue
            n = cr / nrm
            if n[2] < 0.0:
                n = -n
            nz = abs(float(n[2]))
            if nz < 1e-9:
                continue
            a = 0.5 * nrm
            ev = eps.get(e.label)
            if ev is None:
                continue
            th2.append(math.atan2(math.hypot(n[0], n[1]), nz) ** 2)
            area.append(a)
            epsv.append(ev)
        th2 = np.array(th2); area = np.array(area); epsv = np.array(epsv)
        w = area / area.sum() if area.sum() > 0 else None
        theta2 = float(np.sum(w * th2)) if w is not None else float('nan')
        eps_m = float(np.sum(w * epsv)) if w is not None else float('nan')
        p = points.get(_tag(path), {})
        hf = float(p['hf']) if p.get('hf') not in (None, '') else float('nan')
        pred = (1.0 + eps_m) * K * theta2
        return {
            'tag': _tag(path), 'x_c': p.get('x_c', ''), 'd_c': p.get('d_c', ''),
            'theta2': '%.10e' % theta2, 'theta_rms_rad': '%.10e' % math.sqrt(theta2),
            'eps_m': '%.10e' % eps_m, 'strain_field': ekey or '',
            'hf_code': '%.10e' % hf, 'hf_pred': '%.10e' % pred,
            'hf_ratio': '%.6f' % (pred / hf if hf else float('nan')),
            'f_proxy_theta2_eps': '%.6f' % (theta2 / eps_m if eps_m else float('nan')),
            'zk_state_wrinkled_frac': p.get('zk_state_wrinkled_frac', ''),
        }
    finally:
        odb.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--odb')
    ap.add_argument('--out', default='results/multiscale/s0_metrics.csv')
    ap.add_argument('--points', default='results/decision_grid/points.csv')
    args = ap.parse_args()
    points = {}
    with open(args.points, newline='') as f:
        for r in csv.DictReader(f):
            points[r['tag']] = r
    paths = sorted(glob.glob('results/odb/HF_*.odb')) if args.all else [args.odb]
    rows = []
    for p in paths:
        try:
            r = one(p, points)
            rows.append(r)
            print('[ok] %s theta2=%.4e eps=%.4e hf=%.4e ratio=%.4f'
                  % (r['tag'], float(r['theta2']), float(r['eps_m']),
                     float(r['hf_code']), float(r['hf_ratio'])))
        except Exception as ex:
            print('[fail] %s %s' % (p, ex))
    fields = list(rows[0].keys()) if rows else ['tag']
    d = os.path.dirname(args.out)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with open(args.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    print('wrote %s (%d rows)' % (args.out, len(rows)))


if __name__ == '__main__':
    main()
