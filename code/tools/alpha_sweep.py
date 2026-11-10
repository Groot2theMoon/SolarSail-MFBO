# -*- coding: utf-8 -*-
"""Alpha/clock sweep over existing HF ODBs (no solves).

For every ODB and (alpha, clock):
  tau(alpha, phi), cp(alpha, phi), hf_alpha = 1 - actual/ideal(alpha).

Run (Abaqus python, from code/):
  abaqus python tools/alpha_sweep.py --out results/alpha_sweep/alpha_sweep.csv
"""
from __future__ import print_function

import argparse
import csv
import glob
import hashlib
import io
import os
import re
import sys

import numpy as np

_HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from odbAccess import openOdb  # noqa: E402
import eval_abaqus as ev  # noqa: E402

DEFAULT_ALPHAS = (0.0, 2.5, 5.0, 7.5, 10.0)
DEFAULT_PHIS = (0.0, 90.0)


def parse_tag(path):
    m = re.search(r'HF_(.+?)_u3free\.odb$', os.path.basename(path))
    return m.group(1) if m else os.path.splitext(os.path.basename(path))[0]


def sidecar(tag):
    p = os.path.join(_HERE, 'results', 'results_%s.txt' % tag)
    d = {}
    if not os.path.exists(p):
        return d
    with io.open(p, encoding='utf-8', errors='replace') as f:
        for line in f:
            if '=' in line:
                k, v = line.strip().split('=', 1)
                d[k] = v
    return d


def area_ref_of(odb, inst, conns):
    """Pristine reference area from frameValue=0, same convention as get_hf()."""
    f0 = None
    for sn in odb.steps.keys():
        for fr in odb.steps[sn].frames:
            if abs(fr.frameValue) < 1e-12:
                f0 = fr
                break
        if f0 is not None:
            break
    if f0 is None or 'COORD' not in f0.fieldOutputs.keys():
        return None
    max_label = (inst.nodes[-1].label + 100) if inst.nodes else 100000
    nc = np.zeros((max_label + 1, 3), dtype=np.float64)
    cf = f0.fieldOutputs['COORD'].getSubset(region=inst)
    for b in cf.bulkDataBlocks:
        if b.nodeLabels is not None:
            nc[b.nodeLabels] = b.data
    total = 0.0
    q = [c for c in conns if len(c) in (4, 8)]
    if q:
        a = np.array([c[:4] for c in q], dtype=np.int32)
        total += 0.5 * np.linalg.norm(
            np.cross(nc[a[:, 2]] - nc[a[:, 0]], nc[a[:, 3]] - nc[a[:, 1]]), axis=1).sum()
    t = [c for c in conns if len(c) in (3, 6)]
    if t:
        a = np.array([c[:3] for c in t], dtype=np.int32)
        total += 0.5 * np.linalg.norm(
            np.cross(nc[a[:, 1]] - nc[a[:, 0]], nc[a[:, 2]] - nc[a[:, 0]]), axis=1).sum()
    return float(total) if total > 0.0 else None


def sha1(path):
    h = hashlib.sha1()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()[:12]


def evaluate(path, alpha, phi):
    odb = openOdb(path=path, readOnly=True)
    try:
        step_name, frame, names = ev._pick_hf_step(odb)
        if step_name is None or frame is None:
            raise RuntimeError('no HF frame: %s' % (names,))
        if frame.frameValue < 0.99:
            raise RuntimeError('incomplete frame=%.4f' % frame.frameValue)
        inst = odb.rootAssembly.instances['MEMBRANE-1']
        conns = [e.connectivity for e in inst.elements]
        max_label = (inst.nodes[-1].label + 100) if inst.nodes else 100000
        aref = area_ref_of(odb, inst, conns)
        us = ev.sun_vector(alpha, phi)
        ideal = ev.ideal_thrust_coeff(alpha)
        extra = {'n_total': 0, 'n_backward': 0, 'area_total': 0.0, 'area_backward': 0.0}
        diag = {}
        ev.calc_thrust_loss(frame, inst, conns, max_label, area_ref=aref,
                            extra_out=extra, diag=diag, u_sun=us, ideal_coeff=ideal)
        nb = (extra['area_backward'] / extra['area_total']
              if extra['area_total'] > 0.0 else float('nan'))
        return {
            'step': step_name, 'area_ref_m2': aref, 'ideal_coeff': ideal,
            'hf_alpha': diag.get('hf_signed', float('nan')),
            'torque_nm': diag.get('torque_nm', float('nan')),
            'tau_x_nm': diag.get('tau_x_nm', float('nan')),
            'tau_y_nm': diag.get('tau_y_nm', float('nan')),
            'tau_z_nm': diag.get('tau_z_nm', float('nan')),
            'cp_x_m': diag.get('cp_x_m', float('nan')),
            'cp_y_m': diag.get('cp_y_m', float('nan')),
            'cp_offset_m': diag.get('cp_offset_m', float('nan')),
            'force_x_n': diag.get('force_x_n', float('nan')),
            'force_y_n': diag.get('force_y_n', float('nan')),
            'force_z_n': diag.get('force_z_n', float('nan')),
            'force_mag_n': diag.get('force_magnitude_n', float('nan')),
            'n_backward_frac': nb,
        }
    finally:
        odb.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='results/alpha_sweep/alpha_sweep.csv')
    ap.add_argument('--odb', default=None,
                    help='single ODB; default = results/odb/HF_*.odb')
    ap.add_argument('--baseline', default=None,
                    help='optional baseline ODB, e.g. results/variants/o1_baseline_clampfree/HF_xc0p2_dc0p0001_u3free.odb')
    ap.add_argument('--alphas', default=','.join('%.1f' % x for x in DEFAULT_ALPHAS))
    ap.add_argument('--phis', default=','.join('%.1f' % x for x in DEFAULT_PHIS))
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--exclude', default='',
                    help='comma-separated tag substrings to drop from results/odb')
    a = ap.parse_args()

    paths = ([a.odb] if a.odb else sorted(glob.glob(os.path.join(_HERE, 'results', 'odb', 'HF_*.odb'))))
    if a.exclude:
        drops = [x.strip() for x in a.exclude.split(',') if x.strip()]
        paths = [p for p in paths if not any(d in parse_tag(p) for d in drops)]
    if a.baseline and os.path.exists(a.baseline):
        paths.append(a.baseline)
    if a.limit > 0:
        paths = paths[:a.limit]
    alphas = [float(x) for x in a.alphas.split(',') if x.strip()]
    phis = [float(x) for x in a.phis.split(',') if x.strip()]
    if not paths:
        raise SystemExit('no ODBs found')

    d = os.path.dirname(a.out)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    fields = ('odb', 'tag', 'case', 'config_sha1', 'lf_model_fp', 'step',
              'alpha_deg', 'clock_deg', 'alpha_rad', 'clock_rad',
              'u_sun_x', 'u_sun_y', 'u_sun_z', 'ideal_coeff', 'area_ref_m2',
              'hf_alpha', 'hf_sensored', 'torque_nm', 'tau_x_nm', 'tau_y_nm', 'tau_z_nm',
              'cp_x_m', 'cp_y_m', 'cp_offset_m',
              'force_x_n', 'force_y_n', 'force_z_n', 'force_mag_n',
              'n_backward_frac', 'odb_sha1')
    rows = []
    for pi, path in enumerate(paths):
        base_tag = parse_tag(path)
        is_base = bool(a.baseline and os.path.abspath(path) == os.path.abspath(a.baseline))
        tag = ('baseline_' + base_tag) if is_base else base_tag
        case = 'baseline' if (is_base or 'baseline' in path.lower() or 'noclamp' in path.lower()) else 'clamp'
        sc = sidecar(base_tag)
        cfg = sc.get('config_sha1', '')
        lfp = sc.get('lf_model_fp', '')
        oh = sha1(path)
        print('[%d/%d] %s' % (pi + 1, len(paths), tag))
        for alpha in alphas:
            for phi in phis:
                try:
                    r = evaluate(path, alpha, phi)
                    row = {'odb': os.path.relpath(path, _HERE).replace('\\', '/'),
                           'tag': tag, 'case': case,
                           'config_sha1': cfg, 'lf_model_fp': lfp,
                           'alpha_deg': alpha, 'clock_deg': phi,
                           'alpha_rad': float(np.deg2rad(alpha)),
                           'clock_rad': float(np.deg2rad(phi)),
                           'u_sun_x': float(ev.sun_vector(alpha, phi)[0]),
                           'u_sun_y': float(ev.sun_vector(alpha, phi)[1]),
                           'u_sun_z': float(ev.sun_vector(alpha, phi)[2]),
                           'hf_sensored': max(r['hf_alpha'], 0.0) if r['hf_alpha'] == r['hf_alpha'] else float('nan'),
                           'odb_sha1': oh}
                    row.update(r)
                    rows.append(row)
                except Exception as e:
                    print('  ERROR alpha=%.1f phi=%.1f: %s' % (alpha, phi, e))
                    rows.append({'odb': os.path.relpath(path, _HERE).replace('\\', '/'),
                                 'tag': tag, 'case': case, 'config_sha1': cfg,
                                 'lf_model_fp': lfp, 'alpha_deg': alpha, 'clock_deg': phi,
                                 'error': str(e), 'odb_sha1': oh})
    with io.open(a.out, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields + ('error',), extrasaction='ignore')
        w.writeheader()
        for row in rows:
            w.writerow(row)
    print('wrote %s : %d rows' % (a.out, len(rows)))


if __name__ == '__main__':
    main()
