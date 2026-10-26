#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""analysis_config.py - one explicit JSON config for buckle/LF/HF runs.

No environment-variable overrides are used.  The config path and buckle CASE are
passed on the command line:

    abaqus cae noGUI=run_abaqus.py -- HF 0.5 1.5 --config configs/paper_sd632.json
    abaqus cae noGUI=run_abaqus_buckle.py -- 0.5 1e-4 0 \
        --case clamp_only --config configs/paper_sd632.json
    python make_seed.py --config configs/paper_sd632.json ...
    python run_analysis.py --config configs/paper_sd632.json

The JSON is a flat override map.  Keys that match a script global are applied;
unknown keys (for example the "RUN" section consumed by run_analysis.py) are
ignored by Abaqus scripts.  This keeps the first migration step small while
making the config authoritative for every key it contains.
"""
from __future__ import print_function

import hashlib
import json
import os
import sys


def _split_flag(argv, flag):
    rest = []
    value = None
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == flag:
            if i + 1 >= len(argv):
                raise RuntimeError('%s 뒤에 값이 필요합니다.' % flag)
            value = argv[i + 1]
            i += 2
        elif a.startswith(flag + '='):
            value = a.split('=', 1)[1]
            i += 1
        else:
            rest.append(a)
            i += 1
    return rest, value


def split_cli(argv):
    """Return (rest, config_path, case).

    --config and --case can appear anywhere.  The remaining tokens keep the
    Abaqus launcher tokens; run_abaqus.py reads its last three as fidelity/x/d.
    """
    rest, cfg = _split_flag(argv, '--config')
    rest, case = _split_flag(rest, '--case')
    return rest, cfg, case


def load_config(path=None):
    if not path:
        return {}
    p = os.path.abspath(os.path.expanduser(path))
    if not os.path.isfile(p):
        raise RuntimeError('config 파일을 찾을 수 없습니다: %s' % p)
    try:
        with open(p, 'r', encoding='utf-8') as fh:
            cfg = json.load(fh)
    except Exception as exc:
        raise RuntimeError('config JSON 읽기 실패 (%s): %s' % (p, exc))
    if not isinstance(cfg, dict):
        raise RuntimeError('config 최상위는 JSON object 여야 합니다: %s' % p)
    cfg['_CONFIG_PATH'] = p
    return cfg


def config_hash(cfg):
    """Deterministic short hash of the applied config (RUN is metadata only)."""
    skip = {'_CONFIG_PATH', 'RUN', '_comment', '_notes'}
    body = dict((k, v) for k, v in (cfg or {}).items() if k not in skip)
    raw = json.dumps(body, sort_keys=True, separators=(',', ':'),
                     ensure_ascii=True)
    return hashlib.sha1(raw.encode('utf-8')).hexdigest()[:12]


def _as_tuple(v):
    if isinstance(v, (list, tuple)):
        return tuple(v)
    return v


def apply_globals(g, cfg):
    """Apply flat config keys to a script's globals.

    Only keys already present in the script are overridden.  A few derived
    flags are kept consistent when their parents are supplied.
    """
    if not cfg:
        return
    skip = {'_CONFIG_PATH', 'RUN', '_comment', '_notes'}
    tuple_keys = ('IMPERFECTION_MODES', 'STABILIZE_GLOBAL_LADDER',
                  'STABILIZE_POST_LADDER')
    for key, value in cfg.items():
        if key in skip:
            continue
        if key not in g:
            continue
        g[key] = _as_tuple(value) if key in tuple_keys else value

    # derived flags used by run_abaqus.py
    if 'PAPER_INCLUDE_CLAMPS' in cfg and 'NO_CLAMP' in g:
        g['NO_CLAMP'] = not bool(g['PAPER_INCLUDE_CLAMPS'])
    if 'PAPER_MESH_ORDER' in cfg and 'MESH_ORDER' not in cfg and 'MESH_ORDER' in g:
        g['MESH_ORDER'] = (int(g['PAPER_MESH_ORDER'])
                           if g.get('PAPER_MATCH') else 1)
    if 'PAPER_MATCH' in cfg:
        if 'LOAD_MODE' not in cfg and 'LOAD_MODE' in g:
            g['LOAD_MODE'] = 'cable_stress' if g['PAPER_MATCH'] else 'disp'
        if ('PAPER_EDGE_BC_FROM_INITIAL' not in cfg
                and 'PAPER_EDGE_BC_FROM_INITIAL' in g):
            g['PAPER_EDGE_BC_FROM_INITIAL'] = bool(g['PAPER_MATCH'])
    if 'CABLE_RADIUS' in cfg and 'CABLE_AREA' in g:
        try:
            g['CABLE_AREA'] = 3.141592653589793 * (float(g['CABLE_RADIUS']) ** 2)
        except Exception:
            pass
    if 'HF_JOB' in cfg and 'HF_ODB' in g:
        g['HF_ODB'] = str(g['HF_JOB']) + '.odb'


def log_config(emit, cfg, case=None):
    path = (cfg or {}).get('_CONFIG_PATH', '<defaults>')
    emit('[CONFIG] path=%s sha1=%s case=%s'
         % (path, config_hash(cfg), case if case is not None else '-'))


def print_config(cfg):
    print('[CONFIG] path=%s sha1=%s'
          % ((cfg or {}).get('_CONFIG_PATH', '<defaults>'), config_hash(cfg)))


if __name__ == '__main__':
    # Tiny diagnostic: python analysis_config.py configs/paper_sd632.json
    if len(sys.argv) < 2:
        print('usage: python analysis_config.py <config.json>')
        raise SystemExit(2)
    print_config(load_config(sys.argv[1]))
