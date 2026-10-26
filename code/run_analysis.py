#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""run_analysis.py - run one explicit JSON config through buckle/LF/HF.

Usage:
    python run_analysis.py --config configs/paper_sd632.json
    python run_analysis.py --config configs/paper_sd632.json --dry-run

The config has a metadata-only "RUN" section (ignored by Abaqus scripts):
    x_c, d_c, fidelities, make_seed, force_seed, buckle, buckle_case, clean_lf

No environment variables are used.  Every command is printed and executed with
cwd=code.
"""
from __future__ import print_function

import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import analysis_config as AC
import aba_imperfection as A


def _run(cmd, dry=False):
    print('  $ %s' % cmd)
    if dry:
        return 0
    p = subprocess.run(cmd, shell=True, cwd=HERE)
    if p.returncode:
        print('  !!! exit=%d' % p.returncode)
    return p.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()

    cfg_path = os.path.abspath(args.config)
    cfg = AC.load_config(cfg_path)
    AC.print_config(cfg)
    run = cfg.get('RUN', {}) if isinstance(cfg.get('RUN', {}), dict) else {}
    x_c = float(run.get('x_c', 0.5))
    d_c = float(run.get('d_c', 1.5))
    fids = [str(x).upper() for x in run.get('fidelities', ['LF', 'HF'])]
    print('[RUN] x_c=%.4g d_c=%.4g fidelities=%s config=%s'
          % (x_c, d_c, fids, cfg_path))

    if bool(run.get('make_seed', True)):
        cmd = ('"%s" "%s" --config "%s" --x_c %g --d_c %g'
               % (sys.executable, os.path.join(HERE, 'make_seed.py'),
                  cfg_path, x_c, d_c))
        if bool(run.get('force_seed', True)):
            cmd += ' --force'
        rc = _run(cmd, args.dry_run)
        if rc:
            return rc

    if bool(run.get('clean_lf', False)):
        odb = os.path.join(HERE, 'aba', A.lf_odb_name(x_c, d_c))
        for p in (odb, odb + A.LF_SIDECAR_SUFFIX):
            if os.path.exists(p):
                print('  rm %s' % p)
                if not args.dry_run:
                    try:
                        os.remove(p)
                    except OSError as exc:
                        print('  !!! remove failed: %s' % exc)
                        return 1

    for fid in fids:
        cmd = ('abaqus cae noGUI="%s" -- %s %g %g --config "%s"'
               % (os.path.join(HERE, 'run_abaqus.py'), fid, x_c, d_c, cfg_path))
        rc = _run(cmd, args.dry_run)
        if rc:
            return rc

    if bool(run.get('buckle', False)):
        case = str(run.get('buckle_case', 'paper'))
        cmd = ('abaqus cae noGUI="%s" -- --config "%s" --case %s %g'
               % (os.path.join(HERE, 'run_abaqus_buckle.py'),
                  cfg_path, case, x_c))
        rc = _run(cmd, args.dry_run)
        if rc:
            return rc

    print('[RUN] done')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
