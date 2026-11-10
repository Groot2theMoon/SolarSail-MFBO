# -*- coding: utf-8 -*-
"""Run S2: 10 seeds x (1 baseline + 7 low-loss clamp designs) = 80 HF runs.

Run from code/ with normal Python. Existing results/odb 25-point set is backed up
first; each S2 ODB/result is copied to results/variants/uq_s2 and the 25-point
set is restored at the end.
"""
from __future__ import print_function

import csv
import io
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ODB = os.path.join(HERE, 'results', 'odb')
RES = os.path.join(HERE, 'results')
VAR = os.path.join(HERE, 'results', 'variants', 'uq_s2')
PRE = os.path.join(HERE, 'results', 'variants', 'pre_s2_25')
OUT = os.path.join(VAR, 's2_runs.csv')

CLAMP = [
    (0.7, 0.01, 'xc0p7_dc0p01'),
    (0.7, 0.0556, 'xc0p7_dc0p0556'),
    (0.7, 0.0001, 'xc0p7_dc0p0001'),
    (0.7, 0.1, 'xc0p7_dc0p1'),
    (0.8, 0.0001, 'xc0p8_dc0p0001'),
    (0.8, 0.01, 'xc0p8_dc0p01'),
    (0.2, 0.0556, 'xc0p2_dc0p0556'),
]


def ensure_dir(p):
    if not os.path.isdir(p):
        os.makedirs(p)


def copy_if(p, q):
    if os.path.exists(p):
        shutil.copy2(p, q)


def backup25():
    ensure_dir(PRE)
    if not os.path.exists(PRE):
        return
    for fn in os.listdir(ODB) if os.path.isdir(ODB) else []:
        if fn.lower().endswith('.odb'):
            copy_if(os.path.join(ODB, fn), os.path.join(PRE, fn))
    for fn in os.listdir(RES):
        if fn.startswith('results_') and fn.endswith('.txt'):
            copy_if(os.path.join(RES, fn), os.path.join(PRE, fn))
    copy_if(os.path.join(RES, 'all_results.csv'), os.path.join(PRE, 'all_results.csv'))


def restore25():
    for fn in os.listdir(PRE) if os.path.isdir(PRE) else []:
        if fn.lower().endswith('.odb'):
            copy_if(os.path.join(PRE, fn), os.path.join(ODB, fn))
        elif fn.startswith('results_') and fn.endswith('.txt'):
            copy_if(os.path.join(PRE, fn), os.path.join(RES, fn))
    copy_if(os.path.join(PRE, 'all_results.csv'), os.path.join(RES, 'all_results.csv'))


def result_hf(path):
    if not os.path.exists(path):
        return None
    with io.open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('hf='):
                try:
                    return float(line.strip().split('=', 1)[1])
                except Exception:
                    return None
    return None


def run_one(tag, table, xc, dc, baseline):
    env = dict(os.environ)
    env['SOLARSAIL_MODE_TABLE'] = table
    env['SOLARSAIL_SKIP_SEED_GATE'] = '1'
    cmd = 'abaqus cae noGUI=run_abaqus.py -- HF %g %g' % (xc, dc)
    if baseline:
        cmd += ' --config configs/baseline_noclamp.json'
    log = os.path.join(VAR, 'run_%s.log' % tag)
    print('RUN', tag, cmd)
    with io.open(log, 'w', encoding='utf-8', errors='replace') as lf:
        p = subprocess.Popen(cmd, cwd=HERE, env=env, shell=True,
                             stdout=lf, stderr=subprocess.STDOUT)
        rc = p.wait()
    odb = os.path.join(ODB, 'HF_xc%g_dc%g_u3free.odb' % (xc, dc))
    # job name is fixed; use the actual tag-derived result path below
    odb = os.path.join(ODB, 'HF_%s_u3free.odb' % tag.replace('s2_', ''))
    rp = os.path.join(RES, 'results_%s.txt' % tag.replace('s2_', ''))
    # run_abaqus names outputs by x_c/d_c; locate them robustly
    src_odb = os.path.join(ODB, 'HF_xc%s_dc%s_u3free.odb'
                           % (('%g' % xc).replace('.', 'p'), ('%g' % dc).replace('.', 'p')))
    src_res = os.path.join(RES, 'results_xc%s_dc%s.txt'
                           % (('%g' % xc).replace('.', 'p'), ('%g' % dc).replace('.', 'p')))
    hf = result_hf(src_res)
    copy_if(src_odb, os.path.join(VAR, '%s.odb' % tag))
    copy_if(src_res, os.path.join(VAR, '%s.txt' % tag))
    copy_if(log, os.path.join(VAR, '%s.run.log' % tag))
    return rc, hf, os.path.join(VAR, '%s.odb' % tag), os.path.join(VAR, '%s.txt' % tag)


def main():
    ensure_dir(VAR)
    backup25()
    rows = []
    for seed in range(10):
        tag = 's2_base_d%d' % seed
        tab = os.path.join(HERE, 'seeds', 'uq_s2', 's2_base_d%d.txt' % seed)
        rc, hf, odb, res = run_one(tag, tab, 0.2, 0.0001, True)
        rows.append((tag, 'baseline', 0.2, 0.0001, seed, rc, hf, odb, res))
        with io.open(OUT, 'w', encoding='utf-8', newline='') as f:
            w = csv.writer(f)
            w.writerow(['tag', 'kind', 'x_c', 'd_c', 'seed', 'rc', 'hf', 'odb', 'results'])
            w.writerows(rows)
    for xc, dc, base_tag in CLAMP:
        for seed in range(10):
            tag = 's2_%s_d%d' % (base_tag, seed)
            tab = os.path.join(HERE, 'seeds', 'uq_s2', 's2_%s_dc%s_four_d%d.txt'
                               % (base_tag.split('_dc')[0], ('%g' % dc).replace('.', 'p'), seed))
            rc, hf, odb, res = run_one(tag, tab, xc, dc, False)
            rows.append((tag, 'clamp', xc, dc, seed, rc, hf, odb, res))
            with io.open(OUT, 'w', encoding='utf-8', newline='') as f:
                w = csv.writer(f)
                w.writerow(['tag', 'kind', 'x_c', 'd_c', 'seed', 'rc', 'hf', 'odb', 'results'])
                w.writerows(rows)
    restore25()
    print('S2 DONE rows=%d' % len(rows))
    return 0


if __name__ == '__main__':
    sys.exit(main())
