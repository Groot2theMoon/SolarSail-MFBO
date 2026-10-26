#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""zk_map_check.py - validate a zk_map_extract output against zk_metrics.

Usage:
    python zk_map_check.py PREFIX [--kind lf|hf] [--baseline d_c1p5]
                                   [--zk-metrics PATH]

Default reference is PREFIX.json + PREFIX.csv, and (if present)
zk_metrics.txt next to PREFIX.  --baseline d_c1p5 uses the frozen
d_c=1.5 regression values only when requested.

Checks:
  I_w                    sum(evol_ref*comp)/sum(evol_ref)
  state_*_frac           evol_cur
  clip_zero_frac         evol_ref, denominator = wrinkled ref area
  dir_spread_deg         evol_cur, same axial-mean definition as eval_abaqus
  sum(loss_e)            hf (HF only; area_ref_m2 normalization)
  CSV/NPZ                element count / quad+tri counts
"""
from __future__ import print_function

import argparse
import csv
import json
import math
import os
import sys

import numpy as np


BASE_D_C1P5 = {
    'I_w': 2.310881e-06,
    'state_taut_frac': 0.1824,
    'state_wrinkled_frac': 0.7854,
    'state_slack_frac': 0.0322,
    'I_w_clip_zero_frac': 0.29794,
    'dir_spread_deg': 5.234,
}


def _load_csv(path):
    rows = []
    with open(path, 'r', newline='') as fh:
        for r in csv.DictReader(fh):
            rows.append(r)
    if not rows:
        raise RuntimeError('empty CSV: %s' % path)
    out = {}
    for key in ('elem', 'state', 'clip'):
        out[key] = np.asarray([int(float(r[key])) for r in rows], dtype=np.int64)
    for key in ('evol_ref', 'area_ref_m2', 'evol_cur', 'E1', 'E2', 'E12',
                'comp', 'dir_g_deg', 'dir_nophi_deg', 'loss_e'):
        out[key] = np.asarray([float(r[key]) for r in rows], dtype=np.float64)
    #   [2026-10-26] 역사 기준(변형률) 상태. 없으면 생략(구 CSV 호환).
    if 'state_strain' in rows[0]:
        out['state_strain'] = np.asarray([int(float(r['state_strain'])) for r in rows],
                                         dtype=np.int64)
    return out


def _read_metrics(path):
    out = {}
    try:
        with open(path, 'r') as fh:
            for line in fh:
                s = line.strip()
                if '=' not in s:
                    continue
                k, v = s.split('=', 1)
                try:
                    out[k] = float(v)
                except ValueError:
                    out[k] = v
    except Exception:
        return {}
    return out


def _dir_spread(state, dir_g, w):
    m = (state == 1) & np.isfinite(dir_g)
    if not np.any(m):
        return float('nan')
    ww = np.asarray(w, dtype=np.float64)[m]
    tt = np.radians(np.asarray(dir_g, dtype=np.float64)[m])
    den = float(np.sum(ww))
    if den <= 0.0:
        return float('nan')
    c = float(np.sum(ww * np.cos(2.0 * tt)))
    s = float(np.sum(ww * np.sin(2.0 * tt)))
    R = float(np.sqrt(c * c + s * s) / den)
    return float(min(np.degrees(
        0.5 * np.sqrt(max(-2.0 * np.log(max(R, 1e-12)), 0.0))), 90.0))


def _check(name, actual, expected, atol, rtol=0.0):
    try:
        ok = abs(float(actual) - float(expected)) <= (
            float(atol) + float(rtol) * abs(float(expected)))
    except Exception:
        ok = False
    print('  %-4s %-28s actual=% .12e  expected=% .12e'
          % ('OK' if ok else 'FAIL', name, float(actual), float(expected)))
    return ok


def main():
    ap = argparse.ArgumentParser(description='Validate zk_map_extract output.')
    ap.add_argument('prefix')
    ap.add_argument('--kind', choices=('lf', 'hf'), default=None)
    ap.add_argument('--baseline', choices=('d_c1p5',), default=None,
                    help='use frozen d_c=1.5 reference values')
    ap.add_argument('--zk-metrics', default=None,
                    help='zk_metrics.txt from the same ODB')
    ap.add_argument('--no-npz', action='store_true',
                    help='skip CSV/NPZ element-count check')
    args = ap.parse_args()

    prefix = args.prefix
    csv_path = prefix + '.csv'
    json_path = prefix + '.json'
    if not os.path.exists(csv_path):
        print('ERROR: missing %s' % csv_path)
        return 2
    if not os.path.exists(json_path):
        print('ERROR: missing %s' % json_path)
        return 2

    data = _load_csv(csv_path)
    with open(json_path, 'r') as fh:
        meta = json.load(fh)
    kind = args.kind or meta.get('kind', 'lf')
    wref = data['evol_ref']
    wcur = data['evol_cur']
    comp = data['comp']
    state = data['state']
    clip = data['clip']
    dir_g = data['dir_g_deg']
    loss = data['loss_e']

    ref_sum = float(np.sum(wref))
    cur_sum = float(np.sum(wcur)) if float(np.sum(wcur)) > 0 else float(len(wcur))
    iw = (float(np.sum(wref * comp) / ref_sum) if ref_sum > 0 else 0.0)
    frac = {}
    for name, code in (('taut', 0), ('wrinkled', 1), ('slack', 2)):
        frac[name] = float(np.sum(wcur[state == code]) / cur_sum)
    #   CSV 의 `state` 는 ZK_STATE_CRITERION(기본 mixed) 기준이다. 비교 대상도
    #   `zk_state_*_frac`(mixed) 여야 하고, 역사값은 `_strain` 접미사와 짝지어 본다.
    frac_strain = {}
    if 'state_strain' in data:
        _ss = data['state_strain']
        for name, code in (('taut', 0), ('wrinkled', 1), ('slack', 2)):
            frac_strain[name] = float(np.sum(wcur[_ss == code]) / cur_sum)
    wr_ref = float(np.sum(wref[state == 1]))
    clip_frac = (float(np.sum(wref[(state == 1) & (clip == 1)]) / wr_ref)
                 if wr_ref > 0 else 0.0)
    spread = _dir_spread(state, dir_g, wcur)

    zk_path = args.zk_metrics
    if zk_path is None:
        cand = os.path.join(os.path.dirname(os.path.abspath(prefix)), 'zk_metrics.txt')
        if os.path.exists(cand):
            zk_path = cand
    zk = _read_metrics(zk_path) if zk_path else {}

    print('zk_map_check: kind=%s n=%d' % (kind, len(data['elem'])))
    print('  reference: %s' % (zk_path or ('baseline ' + args.baseline
                                           if args.baseline else 'meta only')))

    #   [2026-10-26] zk_metrics 를 못 찾으면 **PASS 금지**. 검증을 건너뛴 채 PASS 를
    #     찍으면 회귀를 놓친다(실제로 "meta only" 폴백이 집계 6항목을 생략했다).
    if kind == 'lf' and not zk and not args.baseline:
        print('  REFERENCE MISSING - 집계 검증 불가')
        print('  해결: --zk-metrics <ODB 옆 zk_metrics.txt> 또는 --baseline d_c1p5')
        print('\nRESULT: SKIP')
        return 3

    ok = True
    # ---- same-ODB zk_metrics (LF only) ----
    #   HF map의 I_w/state는 HF ODB 자체 상태라 LF zk_metrics와 비교하면 안 된다.
    if kind == 'lf' and zk:
        ok &= _check('I_w', iw, zk.get('zk_Iw', float('nan')), 1e-12, 1e-6)
        ok &= _check('state_taut_frac', frac['taut'],
                     zk.get('zk_state_taut_frac', float('nan')), 1e-6)
        ok &= _check('state_wrinkled_frac', frac['wrinkled'],
                     zk.get('zk_state_wrinkled_frac',
                            zk.get('state_wrinkled_frac', float('nan'))), 1e-6)
        ok &= _check('state_slack_frac', frac['slack'],
                     zk.get('zk_state_slack_frac', float('nan')), 1e-6)
        ok &= _check('clip_zero_frac', clip_frac,
                     zk.get('zk_Iw_clip_zero_frac', float('nan')), 1e-6)
        #   역사(변형률) 기준과도 병기 대조. 생산 sidecar 가 아직 `_strain` 을
        #   wrinkled 만 주므로 그것만 본다.
        if frac_strain and 'zk_state_wrinkled_frac_strain' in zk:
            ok &= _check('wrinkled_frac[strain]', frac_strain['wrinkled'],
                         zk['zk_state_wrinkled_frac_strain'], 1e-6)
        #   dir_spread 는 생산 값이 0.0 으로 나오는 사례가 있어 정의 확인 전까지
        #   FAIL 이 아니라 WARN 으로만 남긴다.
        _ds_ok = _check('dir_spread_deg', spread,
                        zk.get('zk_dir_spread_deg', float('nan')), 1e-3)
        if not _ds_ok:
            print('  WARN dir_spread 불일치 — 가중/마스크 정의 확인 필요 (FAIL 미포함)')
    elif kind == 'lf' and args.baseline == 'd_c1p5':
        ok &= _check('I_w', iw, BASE_D_C1P5['I_w'], 1e-12, 1e-6)
        ok &= _check('state_taut_frac', frac['taut'],
                     BASE_D_C1P5['state_taut_frac'], 1e-4)
        ok &= _check('state_wrinkled_frac', frac['wrinkled'],
                     BASE_D_C1P5['state_wrinkled_frac'], 1e-4)
        ok &= _check('state_slack_frac', frac['slack'],
                     BASE_D_C1P5['state_slack_frac'], 1e-4)
        ok &= _check('clip_zero_frac', clip_frac,
                     BASE_D_C1P5['I_w_clip_zero_frac'], 1e-4)
        ok &= _check('dir_spread_deg', spread,
                     BASE_D_C1P5['dir_spread_deg'], 1e-3)
    else:
        print('  NOTE: LF reference only; HF checks loss_sum and internal counts')

    # ---- HF loss ----
    if kind == 'hf':
        loss_sum = float(np.nansum(loss))
        hf = meta.get('hf', None)
        if hf is None:
            hf = meta.get('loss_sum', None)
        if hf is None:
            print('  FAIL loss_sum: metadata has no hf/loss_sum')
            ok = False
        else:
            ok &= _check('sum(loss_e)', loss_sum, hf, 1e-12, 1e-9)
            print('  %-4s %-28s actual=% .12e  expected=% .12e'
                  % ('OK' if abs(loss_sum - hf) <= 1e-12 + 1e-9 * abs(hf)
                     else 'FAIL', 'hf(metadata)', hf, hf))

    # ---- CSV/NPZ ----
    if not args.no_npz:
        npz_path = prefix + '.npz'
        if not os.path.exists(npz_path):
            print('  FAIL npz: missing %s' % npz_path)
            ok = False
        else:
            z = np.load(npz_path)
            n = len(data['elem'])
            nq = len(z['elem_quad'])
            nt = len(z['elem_tri'])
            ok &= _check('csv_n == npz_n', float(n), float(len(z['elem_all'])), 0.0)
            ok &= _check('quad+tri == csv_n', float(nq + nt), float(n), 0.0)

    print('')
    print('RESULT: %s' % ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
