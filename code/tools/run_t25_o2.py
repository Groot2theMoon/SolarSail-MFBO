# -*- coding: utf-8 -*-
"""Two-point 2.5 um second-order postbuckling with the t25 mode table."""
from __future__ import print_function

import csv
import io
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VAR = os.path.join(HERE, 'results', 'variants', 't25_o2')
ODB = os.path.join(HERE, 'results', 'odb')
RES = os.path.join(HERE, 'results')


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


def run(tag, xc):
    env = dict(os.environ)
    env['SOLARSAIL_MODE_TABLE'] = os.path.join(
        HERE, 'modes_t25_o2_x%s.txt' % ('%g' % xc).replace('.', 'p'))
    env['SOLARSAIL_SKIP_SEED_GATE'] = '1'
    cmd = 'abaqus cae noGUI=run_abaqus.py -- HF %g 0.0001 --config configs/t25_o2_post.json' % xc
    print(tag, cmd, flush=True)
    t0 = time.time()
    log = os.path.join(VAR, 'run_%s.log' % tag)
    with io.open(log, 'w', encoding='utf-8', errors='replace') as lf:
        p = subprocess.Popen(cmd, cwd=HERE, env=env, shell=True,
                             stdout=lf, stderr=subprocess.STDOUT)
        rc = p.wait()
    wall = time.time() - t0
    src_odb = os.path.join(ODB, 'HF_xc%s_dc0p0001_u3free.odb' % ('%g' % xc).replace('.', 'p'))
    src_res = os.path.join(RES, 'results_xc%s_dc0p0001.txt' % ('%g' % xc).replace('.', 'p'))
    d = val(src_res)
    odb = os.path.join(VAR, '%s.odb' % tag)
    res = os.path.join(VAR, '%s.txt' % tag)
    if os.path.exists(src_odb):
        shutil.copy2(src_odb, odb)
    if os.path.exists(src_res):
        shutil.copy2(src_res, res)
    return {'tag': tag, 'x_c': xc, 'd_c': 0.0001, 'rc': rc, 'wall_s': wall,
            'hf': d.get('hf'), 'torque_nm': d.get('torque_nm'),
            'cp_offset_m': d.get('cp_offset_m'), 'odb': odb, 'results': res}


def main():
    ensure(VAR)
    rows = [run('t25_x0p2', 0.2), run('t25_x0p8', 0.8)]
    with io.open(os.path.join(VAR, 't25_o2.csv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    for r in rows:
        print(r)
    return 0


if __name__ == '__main__':
    sys.exit(main())
