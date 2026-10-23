#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""zk_map_extract.py - ODB -> element-wise ZK/I_w/loss map CSV+NPZ+JSON.

Abaqus Python (Abaqus 2023+, Python 3.10):
    abaqus python zk_map_extract.py <odb> --kind {lf|hf} [--out PREFIX]

Outputs (PREFIX defaults to zk_map for LF, hf_map for HF):
    PREFIX.csv    one row per element
    PREFIX.npz    exact polygon geometry (quad/tri separated)
    PREFIX.json   aggregate metrics / provenance

The extractor deliberately reuses hf_by_state.py / eval_abaqus.py:
  - element strains, LE->E principal conversion, I_w integration
  - per-element SRP loss (sum(loss_e) == hf)
This avoids a second, drifting implementation of the metrics.
"""
from __future__ import print_function

import argparse
import csv
import json
import os
import sys

import numpy as np
from odbAccess import openOdb

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from eval_abaqus import (LF_STEP_PREFERENCE, HF_STEP_PREFERENCE, ZK_NU,
                         ZK_DIR_R_MIN, _zk_reference_evol,
                         _zk_element_local_angles, _zk_weighted_percentile)
from hf_by_state import (pick_step, get_instance, element_connectivity,
                         node_coords, elem_areas_normals, elem_strains,
                         principal_and_states, iw_metrics, element_losses)


def _evol_map(frame, inst):
    out = {}
    try:
        if 'EVOL' in frame.fieldOutputs.keys():
            fld = frame.fieldOutputs['EVOL'].getSubset(region=inst)
            for val in fld.values:
                try:
                    d = val.data
                    v = float(d[0] if hasattr(d, '__len__') else d)
                    if np.isfinite(v) and v > 0.0:
                        out[int(val.elementLabel)] = v
                except Exception:
                    pass
    except Exception:
        out = {}
    return out


def _find_ref_frame(odb):
    """first frameValue=0 frame with COORD."""
    for sn in odb.steps.keys():
        for fr in odb.steps[sn].frames:
            try:
                if abs(float(fr.frameValue)) < 1e-12 and 'COORD' in fr.fieldOutputs.keys():
                    return fr
            except Exception:
                continue
    return None


def _wrap_deg(x):
    return ((x + 90.0) % 180.0) - 90.0


def _jf(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else None
    except Exception:
        return None


def _element_centers(node_xyz, conns):
    out = np.zeros((len(conns), 3), dtype=np.float64)
    for i, c in enumerate(conns):
        out[i] = np.mean(node_xyz[np.asarray(c, dtype=np.int64)], axis=0)
    return out


def _dir_metrics(inst, labels, e11, e22, e12, mask_wr, wcur):
    """per-element global/local n1 + axial mean spread (same as eval_abaqus)."""
    nan = float('nan')
    phi = _zk_element_local_angles(inst, labels)
    if phi is None:
        phi = np.full(len(labels), nan, dtype=np.float64)
    phi = np.asarray(phi, dtype=np.float64)
    theta_local = 0.5 * np.arctan2(2.0 * e12, e11 - e22)
    theta_g = theta_local + phi
    theta_g = ((theta_g + 0.5 * np.pi) % np.pi) - 0.5 * np.pi
    valid = mask_wr & np.isfinite(theta_g)
    dir_g = np.full(len(labels), nan, dtype=np.float64)
    dir_nophi = np.full(len(labels), nan, dtype=np.float64)
    dir_g[mask_wr & np.isfinite(theta_g)] = np.degrees(theta_g[mask_wr & np.isfinite(theta_g)])
    dir_nophi[mask_wr & np.isfinite(theta_local)] = np.degrees(
        theta_local[mask_wr & np.isfinite(theta_local)])
    dir_g = _wrap_deg(dir_g)
    dir_nophi = _wrap_deg(dir_nophi)

    spread = nan
    if np.any(valid):
        ww = np.asarray(wcur, dtype=np.float64)[valid]
        tt = theta_g[valid]
        den = float(np.sum(ww))
        if den > 0.0:
            c = float(np.sum(ww * np.cos(2.0 * tt)))
            s = float(np.sum(ww * np.sin(2.0 * tt)))
            R = float(np.sqrt(c * c + s * s) / den)
            spread = float(min(np.degrees(
                0.5 * np.sqrt(max(-2.0 * np.log(max(R, 1e-12)), 0.0))), 90.0))
    return dir_g, dir_nophi, phi, spread


def main():
    ap = argparse.ArgumentParser(description='Extract ZK/I_w/HF map from an ODB.')
    ap.add_argument('odb')
    ap.add_argument('--kind', choices=('lf', 'hf'), required=True)
    ap.add_argument('--out', default=None, help='output prefix')
    args = ap.parse_args()

    kind = args.kind
    prefix = args.out or ('zk_map' if kind == 'lf' else 'hf_map')
    outdir = os.path.dirname(os.path.abspath(prefix))
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir)

    odb = None
    try:
        odb = openOdb(args.odb, readOnly=True)
        inst = get_instance(odb)
        if inst is None:
            print('[zk_map] ERROR: no instance')
            return 1
        pref = LF_STEP_PREFERENCE if kind == 'lf' else HF_STEP_PREFERENCE
        step_name, step = pick_step(odb, pref)
        if step is None:
            print('[zk_map] ERROR: no step')
            return 1
        frame = step.frames[-1]
        frame_value = float(frame.frameValue)
        labels, conns = element_connectivity(inst)
        if len(labels) == 0:
            print('[zk_map] ERROR: no elements')
            return 1
        order = np.argsort(labels)
        labels = np.asarray(labels, dtype=np.int64)[order]
        conns = [conns[int(i)] for i in order]
        n = len(labels)

        maxlab = 0
        for c in conns:
            if len(c):
                maxlab = max(maxlab, int(max(c)))
        node_xyz_def = node_coords(frame, maxlab)
        ref_frame = _find_ref_frame(odb)
        if ref_frame is None:
            print('[zk_map] ERROR: frameValue=0 COORD not found')
            return 1
        node_xyz_ref = node_coords(ref_frame, maxlab)

        areas_def, normals = elem_areas_normals(frame, conns, maxlab)
        areas_ref, _ = elem_areas_normals(ref_frame, conns, maxlab)
        area_ref_total = float(np.sum(areas_ref))

        evol_ref_map = _zk_reference_evol(odb, inst)
        evol_cur_map = _evol_map(frame, inst)
        wref = np.asarray([evol_ref_map.get(int(l), 0.0) for l in labels],
                          dtype=np.float64)
        wcur = np.asarray([evol_cur_map.get(int(l), 0.0) for l in labels],
                          dtype=np.float64)
        wref_sum = float(np.sum(wref))
        wcur_sum = float(np.sum(wcur))
        if wcur_sum <= 0.0:
            wcur = np.ones(n, dtype=np.float64)
            wcur_sum = float(n)

        # ---- strains / state / I_w ----
        e11 = np.full(n, np.nan)
        e22 = np.full(n, np.nan)
        e12 = np.full(n, np.nan)
        strain_field = None
        got = elem_strains(frame, inst)
        if got[0] is not None:
            lab_e, (a, b, c), strain_field = got
            pos = dict((int(l), i) for i, l in enumerate(lab_e))
            for i, l in enumerate(labels):
                j = pos.get(int(l))
                if j is not None:
                    e11[i], e22[i], e12[i] = float(a[j]), float(b[j]), float(c[j])
        elif kind == 'lf':
            print('[zk_map] ERROR: LF strain field missing (%s)' % (got[2],))
            return 1

        e1, e2, m_taut, m_wr, m_sl = principal_and_states(
            e11, e22, e12, strain_field or 'E')
        state = np.full(n, -1, dtype=np.int64)
        state[m_taut] = 0
        state[m_wr] = 1
        state[m_sl] = 2

        wr = np.maximum(np.abs(e2) - ZK_NU * e1, 0.0)
        sl = np.abs(e2) + ZK_NU * np.abs(e1)
        comp = np.zeros(n, dtype=np.float64)
        comp[m_wr] = wr[m_wr]
        comp[m_sl] = sl[m_sl]
        comp[~np.isfinite(comp)] = 0.0
        clip = ((state == 1) & (comp <= 0.0)).astype(np.int64)

        iw_val = float(np.sum(wref * comp) / wref_sum) if wref_sum > 0 else 0.0
        pos = comp > 0.0
        iw_p95 = (_zk_weighted_percentile(comp[pos], wref[pos], 95.0)
                  if np.any(pos) else 0.0)
        iw_num = float(np.sum(wref * comp))
        wr_num = float(np.sum(wref[m_wr] * comp[m_wr]))
        sl_num = float(np.sum(wref[m_sl] * comp[m_sl]))
        iw_wr_share = (wr_num / iw_num) if iw_num > 0 else 0.0
        iw_sl_share = (sl_num / iw_num) if iw_num > 0 else 0.0
        wr_area = float(np.sum(wref[m_wr]))
        clip_frac = (float(np.sum(wref[m_wr & (comp <= 0.0)])) / wr_area
                     if wr_area > 0 else 0.0)

        frac_taut = float(np.sum(wcur[m_taut]) / wcur_sum)
        frac_wr = float(np.sum(wcur[m_wr]) / wcur_sum)
        frac_sl = float(np.sum(wcur[m_sl]) / wcur_sum)

        dir_g, dir_nophi, phi, dir_spread = _dir_metrics(
            inst, labels, e11, e22, e12, m_wr, wcur)

        # ---- HF loss ----
        loss_e = np.full(n, np.nan)
        hf_total = None
        if kind == 'hf':
            loss_e = element_losses(areas_def, normals, area_ref_total)
            hf_total = float(np.sum(loss_e))

        centers = _element_centers(node_xyz_def, conns)

        # ---- CSV ----
        cols = ['elem', 'cx', 'cy', 'cz', 'evol_ref', 'area_ref_m2',
                'evol_cur', 'state', 'E1', 'E2', 'E12', 'comp', 'clip',
                'dir_g_deg', 'dir_nophi_deg', 'loss_e']
        csv_path = prefix + '.csv'
        with open(csv_path, 'w', newline='') as fh:
            wrtr = csv.writer(fh)
            wrtr.writerow(cols)
            for i in range(n):
                wrtr.writerow([
                    int(labels[i]),
                    '%0.12e' % centers[i, 0], '%0.12e' % centers[i, 1],
                    '%0.12e' % centers[i, 2],
                    '%0.12e' % wref[i], '%0.12e' % areas_ref[i],
                    '%0.12e' % wcur[i], int(state[i]),
                    '%0.12e' % e1[i], '%0.12e' % e2[i], '%0.12e' % e12[i],
                    '%0.12e' % comp[i], int(clip[i]),
                    '%0.12e' % dir_g[i], '%0.12e' % dir_nophi[i],
                    ('%0.12e' % loss_e[i]) if np.isfinite(loss_e[i]) else 'nan'])

        # ---- NPZ ----
        qi = [i for i, c in enumerate(conns) if len(c) == 4]
        ti = [i for i, c in enumerate(conns) if len(c) == 3]
        npz_path = prefix + '.npz'
        np.savez(
            npz_path,
            node_xyz_ref=node_xyz_ref,
            node_xyz_def=node_xyz_def,
            conn_quad=(np.asarray([conns[i] for i in qi], dtype=np.int64)
                       if qi else np.zeros((0, 4), dtype=np.int64)),
            elem_quad=labels[qi] if qi else np.zeros((0,), dtype=np.int64),
            conn_tri=(np.asarray([conns[i] for i in ti], dtype=np.int64)
                      if ti else np.zeros((0, 3), dtype=np.int64)),
            elem_tri=labels[ti] if ti else np.zeros((0,), dtype=np.int64),
            elem_all=labels,
        )

        # ---- JSON ----
        conversion = ('green_lagrange_from_log' if strain_field == 'LE'
                      else ('none' if strain_field == 'E' else None))
        meta = {
            'kind': kind, 'odb': os.path.abspath(args.odb),
            'step': step_name, 'frame': int(len(step.frames)),
            'frame_value': frame_value, 'n_elem': int(n),
            'n_wrinkled': int(np.count_nonzero(m_wr)),
            'strain_field': strain_field, 'strain_measure': 'E',
            'strain_conversion': conversion, 'phi_applied': True,
            'zk_nu': float(ZK_NU),
            'area_ref_m2_total': area_ref_total,
            'evol_ref_total': wref_sum, 'evol_cur_total': wcur_sum,
            'I_w': _jf(iw_val), 'I_w_p95': _jf(iw_p95),
            'I_w_area_total': _jf(wref_sum), 'I_w_weights': 'ref_evol',
            'I_w_wr_share': _jf(iw_wr_share),
            'I_w_slack_share': _jf(iw_sl_share),
            'I_w_clip_zero_frac': _jf(clip_frac),
            'state_taut_frac': _jf(frac_taut),
            'state_wrinkled_frac': _jf(frac_wr),
            'state_slack_frac': _jf(frac_sl),
            'dir_spread_deg': _jf(dir_spread),
        }
        if kind == 'hf':
            meta['hf'] = _jf(hf_total)
            meta['loss_sum'] = _jf(hf_total)
        with open(prefix + '.json', 'w') as fh:
            json.dump(meta, fh, indent=2, sort_keys=True)

        print('[zk_map] %s -> %s.csv / .npz / .json' % (kind, prefix))
        print('[zk_map] step=%s frame_value=%s n=%d field=%s' %
              (step_name, frame_value, n, strain_field))
        print('[zk_map] I_w=%.6e state=%.4f/%.4f/%.4f clip=%.4f dir_spread=%.3f' %
              (iw_val, frac_taut, frac_wr, frac_sl, clip_frac, dir_spread))
        if kind == 'hf':
            print('[zk_map] hf=sum(loss_e)=%.10e' % hf_total)
        return 0
    except Exception as exc:
        print('[zk_map] ERROR: %s: %s' % (type(exc).__name__, exc))
        return 1
    finally:
        if odb is not None:
            try:
                odb.close()
            except Exception:
                pass


if __name__ == '__main__':
    sys.exit(main())
