# -*- coding: utf-8 -*-
"""P0/P1: LF torque proxy vs HF torque on the 25-point clamp grid.

Run with abaqus python from code/:
  abaqus python tools/lf_torque_proxy.py --out results/lf_proxy/lf_proxy.csv
"""
from __future__ import print_function

import argparse
import csv
import io
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from odbAccess import openOdb
from abaqusConstants import CENTROID
import eval_abaqus as ev

C_SLOPE = 4.0
P_SRP = 4.56e-6
R0 = np.array([10.0, 10.0 / 3.0, 0.0])


def node_coords(frame, max_label):
    nc = np.zeros((max_label + 1, 3), dtype=np.float64)
    cf = frame.fieldOutputs['COORD']
    for b in cf.bulkDataBlocks:
        if b.nodeLabels is not None:
            nc[b.nodeLabels] = b.data
    return nc


def lf_proxy(path):
    odb = openOdb(path=path, readOnly=True)
    try:
        sn = 'Step-GlobalTension' if 'Step-GlobalTension' in odb.steps else list(odb.steps.keys())[-1]
        st = odb.steps[sn]
        fr = st.frames[-1]
        inst = odb.rootAssembly.instances['MEMBRANE-1']
        max_label = (inst.nodes[-1].label + 100) if inst.nodes else 100000
        nc = node_coords(fr, max_label)
        try:
            lev = fr.fieldOutputs['LE'].getSubset(position=CENTROID, region=inst)
        except Exception:
            lev = fr.fieldOutputs['E'].getSubset(position=CENTROID, region=inst)
        vals = dict((v.elementLabel, v.data) for v in lev.values)
        m_sqrt = np.zeros(3)
        m_lin = np.zeros(3)
        m_abs = 0.0
        m_zk = np.zeros(3)
        n_elem = 0
        for e in inst.elements:
            conn = [getattr(n, 'label', n) for n in e.connectivity]
            corner = conn[:4] if len(conn) in (4, 8) else (conn[:3] if len(conn) in (3, 6) else None)
            if not corner or e.label not in vals:
                continue
            p = nc[corner]
            if len(corner) == 4:
                cr = np.cross(p[2] - p[0], p[3] - p[1])
            else:
                cr = np.cross(p[1] - p[0], p[2] - p[0])
            area = 0.5 * float(np.linalg.norm(cr))
            if area <= 1e-14:
                continue
            d = np.asarray(vals[e.label], dtype=float)
            if len(d) < 3:
                continue
            e11, e22, e12 = d[0], d[1], d[2]
            ten = np.array([[e11, e12], [e12, e22]], dtype=float)
            w, v = np.linalg.eigh(ten)
            # minimum principal strain = compression; eigenvector is its direction
            n2 = v[:, 0]
            if abs(n2[0]) < 1e-15 and abs(n2[1]) < 1e-15:
                continue
            if n2[0] < 0.0 or (abs(n2[0]) < 1e-15 and n2[1] < 0.0):
                n2 = -n2
            n3 = np.array([n2[0], n2[1], 0.0])
            r = p.mean(axis=0) - R0
            w_abs = abs(float(w[0]))
            m_sqrt += area * math.sqrt(w_abs) * np.cross(r, n3)
            m_lin += area * w_abs * np.cross(r, n3)
            m_zk += area * w_abs * np.cross(r, n3)
            m_abs += area * math.sqrt(w_abs)
            n_elem += 1
        return {'step': sn, 'n_elem': n_elem, 'm1_sqrt': m_sqrt * C_SLOPE * P_SRP,
                'm1_lin': m_lin * C_SLOPE * P_SRP, 'm1_zk': m_zk * C_SLOPE * P_SRP,
                'm0_zk': m_abs * C_SLOPE * P_SRP}
    finally:
        odb.close()


def hf_torque(path):
    odb = openOdb(path=path, readOnly=True)
    try:
        sn, fr, _ = ev._pick_hf_step(odb)
        inst = odb.rootAssembly.instances['MEMBRANE-1']
        conns = [e.connectivity for e in inst.elements]
        max_label = (inst.nodes[-1].label + 100) if inst.nodes else 100000
        # area_ref: use initial frame
        aref = None
        f0 = None
        for s in odb.steps.keys():
            for f in odb.steps[s].frames:
                if abs(f.frameValue) < 1e-12:
                    f0 = f
                    break
            if f0 is not None:
                break
        if f0 is not None:
            try:
                cf = f0.fieldOutputs['COORD'].getSubset(region=inst)
                nc = np.zeros((max_label + 1, 3), dtype=np.float64)
                for b in cf.bulkDataBlocks:
                    if b.nodeLabels is not None:
                        nc[b.nodeLabels] = b.data
                total = 0.0
                for c in conns:
                    p = nc[[int(x) for x in c]]
                    if len(c) == 4:
                        total += 0.5 * np.linalg.norm(np.cross(p[2] - p[0], p[3] - p[1]))
                    elif len(c) == 3:
                        total += 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
                    elif len(c) == 8:
                        for q in ((0, 4, 5, 1), (1, 5, 6, 2), (2, 6, 7, 3), (3, 7, 4, 0)):
                            pp = p[list(q)]
                            total += 0.5 * np.linalg.norm(np.cross(pp[2] - pp[0], pp[3] - pp[1]))
                    elif len(c) == 6:
                        total += 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
                aref = total if total > 0 else None
            except Exception:
                aref = None
        diag = {}
        ev.calc_thrust_loss(fr, inst, conns, max_label, area_ref=aref, diag=diag)
        return {'step': sn, 'tau_x': diag.get('tau_x_nm'), 'tau_y': diag.get('tau_y_nm'),
                'tau_z': diag.get('tau_z_nm'), 'tau_mag': diag.get('torque_nm')}
    finally:
        odb.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='results/lf_proxy/lf_proxy.csv')
    a = ap.parse_args()
    pts = os.path.join(HERE, 'results', 'decision_grid', 'points.csv')
    rows = list(csv.DictReader(io.open(pts, encoding='utf-8', errors='replace')))
    out = []
    for r in rows:
        tag = r['tag']
        lf = os.path.join(HERE, 'aba', 'LF_Analysis_%s.odb' % tag)
        hf = os.path.join(HERE, 'results', 'odb', 'HF_%s_u3free.odb' % tag)
        if not (os.path.exists(lf) and os.path.exists(hf)):
            print('MISSING', tag, os.path.exists(lf), os.path.exists(hf))
            continue
        try:
            lp = lf_proxy(lf)
            ht = hf_torque(hf)
            row = {'tag': tag, 'x_c': r['x_c'], 'd_c': r['d_c']}
            for k in ('m1_sqrt', 'm1_lin', 'm1_zk'):
                v = lp[k]
                row[k + '_x'] = float(v[0]); row[k + '_y'] = float(v[1])
                row[k + '_mag'] = float(np.linalg.norm(v))
            row['m0_zk'] = lp['m0_zk']
            row['hf_tau_x'] = ht['tau_x']; row['hf_tau_y'] = ht['tau_y']
            row['hf_tau_mag'] = ht['tau_mag']
            row['hf_cp_offset'] = None
            out.append(row)
            print(tag, 'proxy', row['m1_sqrt_mag'], 'hf', row['hf_tau_mag'])
        except Exception as e:
            print('ERROR', tag, e)
    if out:
        d = os.path.dirname(a.out)
        if d and not os.path.isdir(d):
            os.makedirs(d)
        fields = list(out[0].keys())
        with io.open(a.out, 'w', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader(); w.writerows(out)
        print('wrote', a.out, len(out))
    return 0


if __name__ == '__main__':
    sys.exit(main())
