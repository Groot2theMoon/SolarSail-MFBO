# -*- coding: utf-8 -*-
"""Postprocess one subgrid RVE ODB -> theta^2, s0_meas, f_out.

Usage (from code/):
    abaqus python tools/subgrid_post.py --odb subgrid/RVE....odb --s0 0.01 --mode biaxial --seed 0
"""
import argparse
import csv
import math
import os

import numpy as np
from odbAccess import openOdb

L = 0.5


def _area_theta(path):
    odb = openOdb(path=path, readOnly=True)
    try:
        step = odb.steps['Step-Tension']
        fr = step.frames[-1]
        fr0 = step.frames[0]
        inst = odb.rootAssembly.instances['RVE-1']
        c = {}
        for v in fr.fieldOutputs['COORD'].getSubset(region=inst).values:
            if v.nodeLabel is not None:
                c[v.nodeLabel] = v.data
        c0 = {}
        for v in fr0.fieldOutputs['COORD'].getSubset(region=inst).values:
            if v.nodeLabel is not None:
                c0[v.nodeLabel] = v.data
        th2 = ar = ainit = 0.0
        for e in inst.elements:
            conn = list(e.connectivity)
            cur = np.array([c[n] for n in conn])
            cur0 = np.array([c0[n] for n in conn])
            cr = np.cross(cur[2] - cur[0], cur[3] - cur[1])
            cr0 = np.cross(cur0[2] - cur0[0], cur0[3] - cur0[1])
            nrm = float(np.linalg.norm(cr)); nrm0 = float(np.linalg.norm(cr0))
            if nrm <= 1e-14 or nrm0 <= 1e-14:
                continue
            a = 0.5 * nrm
            ainit += 0.5 * nrm0
            n = cr / nrm
            if n[2] < 0:
                n = -n
            nz = abs(float(n[2]))
            if nz < 1e-9:
                continue
            th2 += a * math.atan2(math.hypot(n[0], n[1]), nz) ** 2
            ar += a
        return th2 / ar if ar else float('nan'), ainit / (L * L) - 1.0
    finally:
        odb.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--odb', required=True)
    ap.add_argument('--s0', type=float, required=True)
    ap.add_argument('--mode', required=True)
    ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--eps', type=float, default=0.005)
    ap.add_argument('--out', default='results/subgrid/subgrid_runs.csv')
    a = ap.parse_args()
    t2, s0m = _area_theta(a.odb)
    row = {'job': os.path.basename(a.odb).replace('.odb', ''),
           's0_target': '%.6e' % a.s0, 'mode': a.mode, 'seed': a.seed,
           'eps': '%.6e' % a.eps, 'status': 'COMPLETED',
           'theta2': '%.10e' % t2, 's0_meas': '%.10e' % s0m,
           'f_out': '%.10e' % (t2 / (2.0 * s0m) if s0m else float('nan'))}
    d = os.path.dirname(a.out)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    exists = os.path.exists(a.out)
    with open(a.out, 'a', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()))
        if not exists:
            w.writeheader()
        w.writerow(row)
    print('theta2=%.6e s0_meas=%.6e f_out=%.6f' % (t2, s0m, t2 / (2.0 * s0m)))


if __name__ == '__main__':
    main()
