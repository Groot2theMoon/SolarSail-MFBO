# -*- coding: utf-8 -*-
"""Run the S1' 18-case RVE matrix (s0 x mode x seed).

Run with normal Python 3 from code/:
    python tools/subgrid_run.py
"""
import csv
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # code/
S0S = (0.005, 0.02)
MODES = ('biaxial', 'uniaxial', 'shear')
SEEDS = (0, 1, 2)
EPS = 0.005
OUT = os.path.join(ROOT, 'results/subgrid/subgrid_runs.csv')


def job_name(s0, mode, seed):
    return 'RVE_s0_%s_%s_s%d' % (('%g' % s0).replace('.', 'p'), mode, seed)


def done():
    if not os.path.exists(OUT):
        return set()
    with open(OUT, newline='') as f:
        return {r.get('job', '') for r in csv.DictReader(f)}


def main():
    have = done()
    for s0 in S0S:
        for mode in MODES:
            for seed in SEEDS:
                name = job_name(s0, mode, seed)
                if name in have:
                    print('[skip]', name)
                    continue
                print('[gen]', name)
                subprocess.run([sys.executable, 'tools/subgrid_inp.py',
                                '--s0', str(s0), '--mode', mode,
                                '--seed', str(seed), '--eps', str(EPS)],
                               cwd=ROOT, check=True)
                print('[run]', name)
                subprocess.run(['abaqus', 'job=' + name, 'input=' + name + '.inp',
                                'interactive'], cwd=os.path.join(ROOT, 'subgrid'),
                               check=False, shell=True)
                odb = os.path.join(ROOT, 'subgrid', name + '.odb')
                if not os.path.exists(odb):
                    print('[fail] no ODB', name)
                    continue
                print('[post]', name)
                subprocess.run(['abaqus', 'python', 'tools/subgrid_post.py',
                                '--odb', 'subgrid/' + name + '.odb',
                                '--s0', str(s0), '--mode', mode, '--seed', str(seed),
                                '--eps', str(EPS)], cwd=ROOT, check=False, shell=True)


if __name__ == '__main__':
    main()
