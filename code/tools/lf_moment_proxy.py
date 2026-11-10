# -*- coding: utf-8 -*-
"""Scalar loss-density first moment from LF, compared component-wise to HF tau."""
from __future__ import print_function

import argparse
import csv
import io
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from odbAccess import openOdb
from abaqusConstants import CENTROID

R0 = np.array([10.0, 10.0 / 3.0])


def node_coords(frame, max_label):
    nc = np.zeros((max_label + 1, 3), dtype=np.float64)
    for b in frame.fieldOutputs['COORD'].bulkDataBlocks:
        if b.nodeLabels is not None:
            nc[b.nodeLabels] = b.data
    return nc


def lf_moment(path, power=1.0):
    odb = openOdb(path=path, readOnly=True)
    try:
        sn = 'Step-GlobalTension' if 'Step-GlobalTension' in odb.steps else list(odb.steps.keys())[-1]
        fr = odb.steps[sn].frames[-1]
        inst = odb.rootAssembly.instances['MEMBRANE-1']
        max_label = (inst.nodes[-1].label + 100) if inst.nodes else 100000
        nc = node_coords(fr, max_label)
        try:
            ev = fr.fieldOutputs['LE'].getSubset(position=CENTROID, region=inst)
        except Exception:
            ev = fr.fieldOutputs['E'].getSubset(position=CENTROID, region=inst)
        vals = dict((v.elementLabel, v.data) for v in ev.values)
        m0 = mx = my = 0.0
        n = 0
        for e in inst.elements:
            conn = [getattr(nd, 'label', nd) for nd in e.connectivity]
            c = conn[:4] if len(conn) in (4, 8) else (conn[:3] if len(conn) in (3, 6) else None)
            if not c or e.label not in vals:
                continue
            p = nc[c]
            area = (0.5 * np.linalg.norm(np.cross(p[2] - p[0], p[3] - p[1]))
                    if len(c) == 4 else
                    0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0])))
            if area <= 1e-14:
                continue
            d = np.asarray(vals[e.label], dtype=float)
            if len(d) < 3:
                continue
            ten = np.array([[d[0], d[2]], [d[2], d[1]]], dtype=float)
            w = np.linalg.eigvalsh(ten)
            wabs = abs(float(w[0]))
            cen = p.mean(axis=0)[:2] - R0
            v = area * (wabs ** power)
            m0 += v
            mx += v * cen[0]
            my += v * cen[1]
            n += 1
        return {'step': sn, 'n': n, 'M0': m0, 'Mx': mx, 'My': my}
    finally:
        odb.close()


def fnum(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except Exception:
        return None


def ranks(a):
    a = np.asarray(a, float)
    order = np.argsort(a, kind='mergesort')
    r = np.empty(len(a), float)
    r[order] = np.arange(1, len(a) + 1)
    return r


def pear(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 3:
        return None
    x = x - x.mean(); y = y - y.mean()
    den = math.sqrt(float(np.dot(x, x) * np.dot(y, y)))
    return float(np.dot(x, y) / den) if den > 0 else None


def spear(x, y):
    return pear(ranks(x), ranks(y))


def signfrac(a, b):
    n = 0; ok = 0
    for x, y in zip(a, b):
        if x is None or y is None or abs(x) < 1e-15 or abs(y) < 1e-15:
            continue
        n += 1
        if (x > 0) == (y > 0):
            ok += 1
    return (ok / float(n) if n else None), n


def cosine(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='results/lf_proxy/lf_moment.csv')
    ap.add_argument('--power', type=float, default=1.0)
    a = ap.parse_args()
    # HF tau from existing proxy CSV
    hf = {}
    hfp = os.path.join(HERE, 'results', 'lf_proxy', 'lf_proxy.csv')
    with io.open(hfp, encoding='utf-8', errors='replace') as f:
        for r in csv.DictReader(f):
            hf[r['tag']] = (fnum(r['hf_tau_x']), fnum(r['hf_tau_y']), fnum(r['hf_tau_mag']))
    pts = list(csv.DictReader(io.open(os.path.join(HERE, 'results', 'decision_grid', 'points.csv'),
                                      encoding='utf-8', errors='replace')))
    rows = []
    for r in pts:
        tag = r['tag']
        lf = os.path.join(HERE, 'aba', 'LF_Analysis_%s.odb' % tag)
        if not os.path.exists(lf) or tag not in hf:
            continue
        try:
            m = lf_moment(lf, a.power)
            tx, ty, tm = hf[tag]
            rows.append({'tag': tag, 'x_c': r['x_c'], 'd_c': r['d_c'],
                         'M0': m['M0'], 'Mx': m['Mx'], 'My': m['My'],
                         'hf_tau_x': tx, 'hf_tau_y': ty, 'hf_tau_mag': tm})
        except Exception as e:
            print('ERROR', tag, e)
    if not rows:
        raise SystemExit('no rows')
    d = os.path.dirname(a.out)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with io.open(a.out, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    out = {'n': len(rows), 'power': a.power, 'mappings': {}}
    hx = [r['hf_tau_x'] for r in rows]; hy = [r['hf_tau_y'] for r in rows]
    mx = [r['Mx'] for r in rows]; my = [r['My'] for r in rows]
    for name, px, py in (('direct', mx, my), ('cross', [-v for v in my], mx)):
        rec = {'pex_x': pear(px, hx), 'spear_x': spear(px, hx),
               'pear_y': pear(py, hy), 'spear_y': spear(py, hy),
               'sign_x': signfrac(px, hx), 'sign_y': signfrac(py, hy),
               'cos_median': float(np.median([cosine([px[i], py[i]], [hx[i], hy[i]])
                                              for i in range(len(rows))]))}
        out['mappings'][name] = rec
    print(json.dumps(out, ensure_ascii=False, indent=1))
    with io.open(a.out.replace('.csv', '.json'), 'w', encoding='utf-8') as f:
        f.write(json.dumps(out, ensure_ascii=False, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
