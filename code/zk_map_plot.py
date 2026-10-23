#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""zk_map_plot.py - render zk_map_extract outputs (normal Python, matplotlib).

Usage:
    python zk_map_plot.py --lf zk_map [--hf hf_map] [--outdir figures]

Figures:
    fig1_state.png          LF taut/wrinkled/slack
    fig2_excess.png         LF I_w integrand (clip hatched)
    fig3_dir.png            LF n1 direction, phi-applied and no-phi panels
    fig4_lf_hf.png          LF potential vs HF loss, reference coordinates
    fig4b_hf_deformed.png   HF loss on deformed coordinates (optional)

LF maps are explicitly labelled as "potential" (u3=0), not measured amplitude.
The direction panel is diagnostic until the shell local-frame question is settled.
"""
from __future__ import print_function

import argparse
import csv
import json
import math
import os

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.colors import Normalize, BoundaryNorm, ListedColormap
from matplotlib.patches import Patch


STATE_COLORS = ListedColormap(('0.72', '#d62728', '#1f77b4'))
STATE_NORM = BoundaryNorm((-0.5, 0.5, 1.5, 2.5), STATE_COLORS.N)


def _load(prefix):
    with open(prefix + '.json', 'r') as fh:
        meta = json.load(fh)
    rows = []
    with open(prefix + '.csv', 'r', newline='') as fh:
        for r in csv.DictReader(fh):
            rows.append(r)
    data = {}
    for key in ('elem', 'state', 'clip'):
        data[key] = np.asarray([int(float(r[key])) for r in rows], dtype=np.int64)
    for key in ('cx', 'cy', 'cz', 'evol_ref', 'area_ref_m2', 'evol_cur',
                'E1', 'E2', 'E12', 'comp', 'dir_g_deg', 'dir_nophi_deg',
                'loss_e'):
        data[key] = np.asarray([float(r[key]) for r in rows], dtype=np.float64)
    z = np.load(prefix + '.npz')
    return meta, data, z


def _polygons(z, data, use_def=False):
    nc = z['node_xyz_def' if use_def else 'node_xyz_ref']
    verts = []
    labels = []
    by_label = {}
    for i, label in enumerate(z['elem_quad']):
        by_label[int(label)] = np.asarray(z['conn_quad'][i], dtype=np.int64)
    for i, label in enumerate(z['elem_tri']):
        by_label[int(label)] = np.asarray(z['conn_tri'][i], dtype=np.int64)
    for label in z['elem_all']:
        c = by_label.get(int(label))
        if c is None:
            continue
        verts.append(nc[c][:, :2])
        labels.append(int(label))
    order = np.argsort(np.asarray(labels, dtype=np.int64))
    verts = [verts[i] for i in order]
    labels = np.asarray(labels, dtype=np.int64)[order]
    # data rows are sorted by elem; align by label just in case
    idx = dict((int(l), i) for i, l in enumerate(data['elem']))
    out = dict((k, np.asarray([data[k][idx[int(l)]] for l in labels]))
               for k in ('state', 'clip', 'evol_ref', 'evol_cur', 'comp',
                         'dir_g_deg', 'dir_nophi_deg', 'loss_e'))
    out['labels'] = labels
    return nc, verts, out


def _axes(ax, nc, title):
    xy = np.asarray(nc[:, :2])
    xy = xy[np.any(np.abs(xy) > 0, axis=1)]
    if len(xy) == 0:
        return
    xmin, ymin = np.min(xy, axis=0)
    xmax, ymax = np.max(xy, axis=0)
    pad = 0.03 * max(xmax - xmin, ymax - ymin)
    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_aspect('equal')
    #   논문용 좌표축: x/y[m] 와 눈금을 명시한다(격자 계단형으로 보이는 문제 방지).
    ax.set_xlabel('x [m]', fontsize=9)
    ax.set_ylabel('y [m]', fontsize=9)
    ax.tick_params(axis='both', which='major', labelsize=8,
                   direction='out', length=3, width=0.6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title(title, fontsize=10)


def _collection(ax, verts, vals, cmap, norm, edge='none', lw=0.0, zorder=2):
    pc = PolyCollection(verts, cmap=cmap, norm=norm, edgecolors=edge,
                        linewidths=lw, zorder=zorder)
    pc.set_array(np.asarray(vals, dtype=np.float64))
    ax.add_collection(pc)
    return pc


def _save(fig, outdir, name, dpi):
    path = os.path.join(outdir, name)
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print('[plot] %s' % path)


def _fig1(meta, data, z, outdir, dpi):
    nc, verts, d = _polygons(z, data)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    pc = _collection(ax, verts, d['state'], STATE_COLORS, STATE_NORM)
    _axes(ax, nc, 'LF potential (u3=0) - state')
    cb = fig.colorbar(pc, ax=ax, ticks=(0, 1, 2), fraction=0.035, pad=0.02)
    cb.ax.set_yticklabels(('taut', 'wrinkled', 'slack'))
    ax.legend(handles=[Patch(color='0.72', label='taut'),
                       Patch(color='#d62728', label='wrinkled'),
                       Patch(color='#1f77b4', label='slack')],
              loc='lower right', fontsize=8, frameon=True)
    _save(fig, outdir, 'fig1_state.png', dpi)


def _fig2(meta, data, z, outdir, dpi):
    nc, verts, d = _polygons(z, data)
    comp = d['comp']
    wr = d['state'] == 1
    vmax = float(np.nanpercentile(comp[wr], 95)) if np.any(wr) else 0.0
    if not np.isfinite(vmax) or vmax <= 0.0:
        vmax = float(np.nanmax(comp)) if np.any(np.isfinite(comp)) else 1.0
    norm = Normalize(0.0, vmax)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    pc = _collection(ax, verts, comp, 'viridis', norm, edge='none', lw=0.0)
    _axes(ax, nc, 'LF potential (u3=0) - I_w integrand')
    cb = fig.colorbar(pc, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label('comp = max(|E2|-nu*E1, 0) or slack term')
    cv = [verts[i] for i in range(len(verts)) if int(d['clip'][i]) == 1]
    if cv:
        hatch = PolyCollection(cv, facecolors='none', edgecolors='k',
                               linewidths=0.15, hatch='///', zorder=4)
        ax.add_collection(hatch)
        ax.legend(handles=[Patch(facecolor='none', edgecolor='k', hatch='///',
                                 label='wrinkled, comp<=0 (clip)')],
                  loc='lower right', fontsize=8, frameon=True)
    _save(fig, outdir, 'fig2_excess.png', dpi)


def _fig3(meta, data, z, outdir, dpi):
    nc, verts, d = _polygons(z, data)
    centers = np.asarray([np.mean(v, axis=0) for v in verts])
    fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for ax, key, title in ((axs[0], 'dir_g_deg', 'n1 + phi (diagnostic)'),
                           (axs[1], 'dir_nophi_deg', 'n1 no phi (local)')):
        m = (d['state'] == 1) & np.isfinite(d[key]) & (d['comp'] > 0.0)
        _axes(ax, nc, 'LF potential (u3=0) - ' + title)
        if np.any(m):
            p95 = float(np.nanpercentile(d['comp'][m], 95))
            p95 = p95 if p95 > 0 else 1.0
            ell = 0.8 * np.sqrt(np.clip(d['comp'][m] / p95, 0.0, 1.0))
            ang = np.radians(d[key][m])
            xy = centers[m]
            dx = 0.5 * ell * np.cos(ang)
            dy = 0.5 * ell * np.sin(ang)
            segs = [((x - dx[i], y - dy[i]), (x + dx[i], y + dy[i]))
                    for i, (x, y) in enumerate(xy)]
            ax.add_collection(LineCollection(segs, colors='k',
                                             linewidths=0.45, alpha=0.75,
                                             zorder=3))
        ax.text(0.02, 0.02, 'length ~ sqrt(comp), normalized',
                transform=ax.transAxes, fontsize=7)
    fig.suptitle('LF potential (u3=0) - direction is diagnostic; phi frame provisional',
                 fontsize=10)
    _save(fig, outdir, 'fig3_dir.png', dpi)


def _fig4(meta_lf, data_lf, z_lf, meta_hf, data_hf, z_hf, outdir, dpi):
    nc_lf, verts_lf, d_lf = _polygons(z_lf, data_lf, use_def=False)
    nc_hf, verts_hf, d_hf = _polygons(z_hf, data_hf, use_def=False)
    fig, axs = plt.subplots(1, 2, figsize=(12.0, 4.6))
    comp = d_lf['comp']
    wr = d_lf['state'] == 1
    vmax_lf = float(np.nanpercentile(comp[wr], 95)) if np.any(wr) else 1.0
    vmax_lf = vmax_lf if vmax_lf > 0 else 1.0
    pcl = _collection(axs[0], verts_lf, comp, 'viridis',
                      Normalize(0.0, vmax_lf))
    _axes(axs[0], nc_lf, 'LF state (LF ODB, u3=0) - I_w integrand')
    cbl = fig.colorbar(pcl, ax=axs[0], fraction=0.035, pad=0.02)
    cbl.set_label('comp (strain)')
    loss = d_hf['loss_e']
    vmax_hf = float(np.nanpercentile(loss, 95)) if np.any(np.isfinite(loss)) else 1.0
    vmax_hf = vmax_hf if vmax_hf > 0 else 1.0
    pch = _collection(axs[1], verts_hf, loss, 'magma',
                      Normalize(0.0, vmax_hf))
    _axes(axs[1], nc_hf, 'HF loss (HF ODB state) - reference layout')
    cbh = fig.colorbar(pch, ax=axs[1], fraction=0.035, pad=0.02)
    cbh.set_label('loss_e (relative to hf denominator)')
    fig.suptitle('LF potential vs HF loss - reference element layout', fontsize=11)
    fig.text(0.5, 0.01,
             'State/field in each panel is computed from its own ODB; '
             'element labels align but LF and HF states can differ.',
             ha='center', fontsize=8)
    _save(fig, outdir, 'fig4_lf_hf.png', dpi)


def _fig4b(meta_hf, data_hf, z_hf, outdir, dpi):
    nc, verts, d = _polygons(z_hf, data_hf, use_def=True)
    loss = d['loss_e']
    vmax = float(np.nanpercentile(loss, 95)) if np.any(np.isfinite(loss)) else 1.0
    vmax = vmax if vmax > 0 else 1.0
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    pc = _collection(ax, verts, loss, 'magma', Normalize(0.0, vmax))
    _axes(ax, nc, 'HF loss (HF ODB state, deformed)')
    cb = fig.colorbar(pc, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label('loss_e')
    _save(fig, outdir, 'fig4b_hf_deformed.png', dpi)


def main():
    ap = argparse.ArgumentParser(description='Plot zk_map_extract outputs.')
    ap.add_argument('--lf', required=True, help='LF output prefix')
    ap.add_argument('--hf', default=None, help='HF output prefix')
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--dpi', type=int, default=200)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    meta_lf, data_lf, z_lf = _load(args.lf)
    _fig1(meta_lf, data_lf, z_lf, args.outdir, args.dpi)
    _fig2(meta_lf, data_lf, z_lf, args.outdir, args.dpi)
    _fig3(meta_lf, data_lf, z_lf, args.outdir, args.dpi)
    if args.hf:
        meta_hf, data_hf, z_hf = _load(args.hf)
        _fig4(meta_lf, data_lf, z_lf, meta_hf, data_hf, z_hf,
              args.outdir, args.dpi)
        _fig4b(meta_hf, data_hf, z_hf, args.outdir, args.dpi)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
