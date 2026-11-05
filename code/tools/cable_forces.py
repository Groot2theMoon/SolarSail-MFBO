# -*- coding: utf-8 -*-
"""Extract cable S11 / force from HF ODBs (Abaqus python).

Usage (from code/):
    abaqus python tools/cable_forces.py --all --out results/cable_forces/cable_forces.csv

The HF paper model has five truss instances:
  INST_CABLE_TOP / RIGHT / LEFT  (main, initial S11 = CABLE_PRESTRESS)
  INST_CABLE_CL / CR             (clamp, initial S11 = d_c * CABLE_PRESTRESS)

S11 is read from the last frame of Step-GlobalTension; the initial value from
frame 0.  F = S11 * CABLE_AREA.  The tool also records the cable direction
(from the frame-0 COORD) and the total RF residual so the equilibrium gate can
be checked without reopening the ODB.
"""
import argparse
import csv
import glob
import math
import os

from odbAccess import openOdb

try:
    from core.constants import CABLE_RADIUS, CABLE_PRESTRESS
except Exception:
    CABLE_RADIUS = 5.0e-4
    CABLE_PRESTRESS = 1.71e7

CABLE_AREA = math.pi * (CABLE_RADIUS ** 2)
_MAIN = ('CABLE_TOP', 'CABLE_RIGHT', 'CABLE_LEFT')
_CLAMP = ('CABLE_CL', 'CABLE_CR')


def _tag(path):
    name = os.path.basename(path)
    if name.startswith('HF_'):
        name = name[3:]
    return name.replace('_u3free.odb', '').replace('.odb', '')


def _s11(frame, inst):
    if 'S' not in frame.fieldOutputs:
        return None
    sub = frame.fieldOutputs['S'].getSubset(region=inst)
    vals = [float(v.data[0]) for v in sub.values]
    return vals


def _direction(inst, frame0):
    if 'COORD' not in frame0.fieldOutputs:
        return (0.0, 0.0, 0.0)
    sub = frame0.fieldOutputs['COORD'].getSubset(region=inst)
    vals = list(sub.values)
    if len(vals) < 2:
        return (0.0, 0.0, 0.0)
    coord = {v.nodeLabel: v.data for v in vals if v.nodeLabel is not None}
    conn = list(inst.elements[0].connectivity)
    if len(coord) >= 2 and conn[0] in coord and conn[1] in coord:
        a, b = coord[conn[0]], coord[conn[1]]
    else:
        a, b = vals[0].data, vals[1].data
    d = [float(b[i] - a[i]) for i in range(3)]
    n = math.sqrt(sum(x * x for x in d))
    return tuple(x / n for x in d) if n > 0 else (0.0, 0.0, 0.0)


def _rf_sum(frame):
    if 'RF' not in frame.fieldOutputs:
        return None
    vals = frame.fieldOutputs['RF'].values
    sx = sy = sz = 0.0
    sa = 0.0
    for v in vals:
        d = v.data
        sx += float(d[0]); sy += float(d[1]); sz += float(d[2])
        sa += math.sqrt(float(d[0]) ** 2 + float(d[1]) ** 2 + float(d[2]) ** 2)
    return sx, sy, sz, sa


def extract(path, points):
    odb = openOdb(path=path, readOnly=True)
    row_tag = _tag(path)
    p = points.get(row_tag, {})
    rows = []
    try:
        step_name = 'Step-GlobalTension' if 'Step-GlobalTension' in odb.steps else list(odb.steps.keys())[-1]
        step = odb.steps[step_name]
        if not step.frames:
            raise RuntimeError('no frames')
        frame0 = step.frames[0]
        last = step.frames[-1]
        for name in sorted(odb.rootAssembly.instances.keys()):
            if not name.upper().startswith('INST_CABLE'):
                continue
            inst = odb.rootAssembly.instances[name]
            s_last = _s11(last, inst)
            s_init = _s11(frame0, inst)
            if not s_last:
                continue
            s11 = sum(s_last) / len(s_last)
            s11_i = (sum(s_init) / len(s_init)) if s_init else float('nan')
            f = s11 * CABLE_AREA
            f_i = s11_i * CABLE_AREA
            dx, dy, dz = _direction(inst, frame0)
            rows.append({
                'tag': row_tag, 'x_c': p.get('x_c', ''), 'd_c': p.get('d_c', ''),
                'inst': name, 'S11_Pa': '%.10e' % s11, 'F_N': '%.10e' % f,
                'F_over_init': ('%.10e' % (f / f_i)) if f_i == f_i and f_i != 0.0 else 'nan',
                'S11_init_Pa': '%.10e' % s11_i, 'dx': '%.10e' % dx,
                'dy': '%.10e' % dy, 'dz': '%.10e' % dz,
            })
        rf = _rf_sum(last)
        if rf is not None:
            sx, sy, sz, sa = rf
            for r in rows:
                r['rf_sum_x_N'] = '%.10e' % sx
                r['rf_sum_y_N'] = '%.10e' % sy
                r['rf_sum_z_N'] = '%.10e' % sz
                r['rf_sum_abs_N'] = '%.10e' % sa
                r['rf_residual_frac'] = ('%.10e' % (math.sqrt(sx * sx + sy * sy + sz * sz) / sa)) if sa > 0 else 'nan'
        return rows
    finally:
        odb.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--all', action='store_true', help='all results/odb/HF_*.odb')
    ap.add_argument('--odb', default=None)
    ap.add_argument('--out', default='results/cable_forces/cable_forces.csv')
    ap.add_argument('--points', default='results/decision_grid/points.csv')
    args = ap.parse_args()

    points = {}
    if os.path.exists(args.points):
        with open(args.points, newline='') as f:
            for r in csv.DictReader(f):
                points[r.get('tag', '')] = r

    if args.all:
        paths = sorted(glob.glob('results/odb/HF_*.odb'))
    elif args.odb:
        paths = [args.odb]
    else:
        raise SystemExit('use --all or --odb')

    all_rows = []
    for path in paths:
        try:
            rows = extract(path, points)
            all_rows.extend(rows)
            print('[ok] %s %d cables' % (_tag(path), len(rows)))
        except Exception as ex:
            print('[fail] %s: %s' % (path, ex))

    fields = ['tag', 'x_c', 'd_c', 'inst', 'S11_Pa', 'F_N', 'F_over_init',
              'S11_init_Pa', 'dx', 'dy', 'dz',
              'rf_sum_x_N', 'rf_sum_y_N', 'rf_sum_z_N', 'rf_sum_abs_N',
              'rf_residual_frac']
    out_dir = os.path.dirname(args.out)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    with open(args.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(all_rows)
    print('wrote %s (%d rows)' % (args.out, len(all_rows)))


if __name__ == '__main__':
    main()
