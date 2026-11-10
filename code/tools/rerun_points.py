# -*- coding: utf-8 -*-
"""Re-run the 25-point HF decision grid from points.csv.

Sequential by default (shared all_results.csv and license). Existing ODB+result
pairs are skipped. Run from code/ with: python tools/rerun_points.py
"""
from __future__ import print_function

import csv
import io
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POINTS = os.path.join(HERE, 'results', 'decision_grid', 'points.csv')
ODB_DIR = os.path.join(HERE, 'results', 'odb')
RES_DIR = os.path.join(HERE, 'results')


def has_pair(tag):
    odb = os.path.join(ODB_DIR, 'HF_%s_u3free.odb' % tag)
    res = os.path.join(RES_DIR, 'results_%s.txt' % tag)
    if not (os.path.exists(odb) and os.path.exists(res)):
        return False
    try:
        with io.open(res, encoding='utf-8', errors='replace') as f:
            return any(l.startswith('hf=') for l in f)
    except Exception:
        return False


def main():
    rows = []
    with io.open(POINTS, encoding='utf-8', errors='replace') as f:
        for r in csv.DictReader(f):
            try:
                rows.append((float(r['x_c']), float(r['d_c']), r['tag']))
            except Exception:
                pass
    if not os.path.isdir(ODB_DIR):
        os.makedirs(ODB_DIR)
    env = dict(os.environ)
    env['SOLARSAIL_SKIP_SEED_GATE'] = '1'
    print('points=%d' % len(rows))
    for i, (xc, dc, tag) in enumerate(rows, 1):
        if has_pair(tag):
            print('[%d/%d] %s SKIP (exists)' % (i, len(rows), tag))
            continue
        cmd = 'abaqus cae noGUI=run_abaqus.py -- HF %g %g' % (xc, dc)
        log = os.path.join(RES_DIR, 'rerun_%s.log' % tag)
        print('[%d/%d] %s  %s' % (i, len(rows), tag, cmd))
        with io.open(log, 'w', encoding='utf-8', errors='replace') as lf:
            p = subprocess.Popen(cmd, cwd=HERE, env=env, shell=True,
                                 stdout=lf, stderr=subprocess.STDOUT)
            rc = p.wait()
        print('    rc=%d  log=%s' % (rc, log))
        if rc != 0:
            print('    FAILED - continuing')
    print('done')
    return 0


if __name__ == '__main__':
    sys.exit(main())
