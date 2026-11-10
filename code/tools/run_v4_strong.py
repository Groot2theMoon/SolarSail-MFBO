# -*- coding: utf-8 -*-
"""V4 strong clamp: 1st vs 2nd order tau at d_c=0.5, x0.2/x0.8.
Records wall seconds per run.
"""
from __future__ import print_function

import csv
import io
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VAR = os.path.join(HERE, 'results', 'variants', 'v4_strong')
ODB = os.path.join(HERE, 'results', 'odb')
RES = os.path.join(HERE, 'results')
DC = 0.5


def ensure(p):
    if not os.path.isdir(p):
        os.makedirs(p)


def val(path):
    d = {}
    if not os.path.exists(path):
        return d
    with io.open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if '=' in line:
                k, v = line.strip().split('=', 1)
                try:
                    d[k] = float(v)
                except Exception:
                    d[k] = v
    return d


def run(tag, xc, order):
    env = dict(os.environ)
    env['SOLARSAIL_SKIP_SEED_GATE'] = '1'
    if order == 2:
        env['SOLARSAIL_MODE_TABLE'] = os.path.join(
            HERE, 'seeds', 'modes_xc%s_o2.txt' % ('%g' % xc).replace('.', 'p'))
    cmd = 'abaqus cae noGUI=run_abaqus.py -- HF %g %g' % (xc, DC)
    if order == 2:
        cmd += ' --config configs/o2_0p1.json'
    log = os.path.join(VAR, 'run_%s.log' % tag)
    print(tag, cmd, flush=True)
    t0 = time.time()
    with io.open(log, 'w', encoding='utf-8', errors='replace') as lf:
        p = subprocess.Popen(cmd, cwd=HERE, env=env, shell=True,
                             stdout=lf, stderr=subprocess.STDOUT)
        rc = p.wait()
    wall = time.time() - t0
    src_odb = os.path.join(ODB, 'HF_xc%s_dc%s_u3free.odb'
                           % (('%g' % xc).replace('.', 'p'), ('%g' % DC).replace('.', 'p')))
    src_res = os.path.join(RES, 'results_xc%s_dc%s.txt'
                           % (('%g' % xc).replace('.', 'p'), ('%g' % DC).replace('.', 'p')))
    d = val(src_res)
    copy_odb = os.path.join(VAR, '%s.odb' % tag)
    copy_res = os.path.join(VAR, '%s.txt' % tag)
    if os.path.exists(src_odb):
        shutil.copy2(src_odb, copy_odb)
    if os.path.exists(src_res):
        shutil.copy2(src_res, copy_res)
    return {'tag': tag, 'order': order, 'x_c': xc, 'd_c': DC, 'rc': rc,
            'wall_s': wall, 'hf': d.get('hf'), 'torque_nm': d.get('torque_nm'),
            'cp_offset_m': d.get('cp_offset_m'), 'tau_x_nm': d.get('tau_x_nm'),
            'tau_y_nm': d.get('tau_y_nm'), 'odb': copy_odb, 'results': copy_res}


def main():
    ensure(VAR)
    rows = []
    for xc in (0.2, 0.8):
        for order in (1, 2):
            tag = '%s_x%g_d0p5' % ('1st' if order == 1 else '2nd', xc)
            rows.append(run(tag, xc, order))
    with io.open(os.path.join(VAR, 'v4_strong.csv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    for r in rows:
        print(r)
    return 0


if __name__ == '__main__':
    sys.exit(main())
