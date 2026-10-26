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
    fig4b_hf_deformed.png   HF loss on 3D (x,y,u3) deformed surface (optional)
    fig5_loss_accum.png     HF loss accumulation curve (optional)

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
from matplotlib.patches import Circle, Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


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


def _polygons(z, data, use_def=False, three_d=False):
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
        verts.append(nc[c] if three_d else nc[c][:, :2])
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


def _vmax(vals, q):
    v = np.asarray(vals, dtype=np.float64).ravel()
    v = v[np.isfinite(v)]
    if v.size == 0:
        return 1.0
    out = float(np.percentile(v, q))
    if not np.isfinite(out) or out <= 0.0:
        out = float(np.max(v))
    return out if out > 0.0 else 1.0


def _clamp_overlay(ax, xc):
    """Reference-coordinate clamp/cable markers (x_c=0.5 -> (5,5),(15,5))."""
    for sx in (-1.0, 1.0):
        x = 10.0 + sx * 10.0 * float(xc)
        y = 10.0 - 10.0 * float(xc)
        ax.add_patch(Circle((x, y), 0.21, fill=False, edgecolor='cyan',
                            linewidth=1.2, zorder=6))
        ax.plot([x], [y], 'c+', markersize=5, zorder=7)
        dx = sx * 0.5 / np.sqrt(2.0)
        dy = 0.5 / np.sqrt(2.0)
        ax.plot([x, x + dx], [y, y + dy], color='cyan', lw=1.2, zorder=6)
    for x, y in ((10.0, 10.0), (20.0, 0.0), (0.0, 0.0)):
        ax.plot([x], [y], 'k^', markersize=4, zorder=6)
    ax.text(10.0, 0.55, 'clamp patches + clamp cables (x_c=%.2f)' % float(xc),
            ha='center', fontsize=7, color='c', zorder=7)


def _fig1(meta, data, z, outdir, dpi):
    nc, verts, d = _polygons(z, data)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    _collection(ax, verts, d['state'], STATE_COLORS, STATE_NORM)
    _axes(ax, nc, 'LF potential (u3=0) - state')
    #   범주형 정보는 범례 하나로만 표시한다(범례는 데이터 밖).
    ax.legend(handles=[Patch(color='0.72', label='taut'),
                       Patch(color='#d62728', label='wrinkled'),
                       Patch(color='#1f77b4', label='slack')],
              loc='upper left', bbox_to_anchor=(1.01, 1.0),
              borderaxespad=0.0, fontsize=8, frameon=False)
    fig.subplots_adjust(right=0.78)
    _save(fig, outdir, 'fig1_state.png', dpi)


def _fig2(meta, data, z, outdir, dpi, clip_q):
    nc, verts, d = _polygons(z, data)
    comp = d['comp']
    wr = d['state'] == 1
    vmax = _vmax(comp[wr] if np.any(wr) else comp, clip_q)
    norm = Normalize(0.0, vmax)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    pc = _collection(ax, verts, comp, 'viridis', norm, edge='none', lw=0.0)
    _axes(ax, nc, 'LF potential (u3=0) - I_w integrand')
    cb = fig.colorbar(pc, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label('comp = max(|E2|-nu*E1, 0) or slack term (p%.1f)' % clip_q)
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
    state = d['state']
    delta = (d['dir_g_deg'] - d['dir_nophi_deg'] + 90.0) % 180.0 - 90.0
    _dm = (state == 1) & np.isfinite(delta)
    max_delta = float(np.nanmax(np.abs(delta[_dm]))) if np.any(_dm) else 0.0
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.4))
    for ax, key, title in ((axs[0], 'dir_g_deg', 'n1 + phi (diagnostic)'),
                           (axs[1], 'dir_nophi_deg', 'n1 no phi (local)')):
        _axes(ax, nc, 'LF potential (u3=0) - ' + title)
        m = (state == 1) & np.isfinite(d[key]) & (d['comp'] > 0.0)
        idx = np.flatnonzero(m)
        if idx.size:
            p95 = float(np.nanpercentile(d['comp'][m], 95))
            p95 = p95 if p95 > 0 else 1.0
            #   subsample to avoid a black blob at 10k elements
            stride = max(1, int(math.ceil(idx.size / 600.0)))
            sel = idx[::stride]
            ell = 0.8 * np.sqrt(np.clip(d['comp'][sel] / p95, 0.0, 1.0))
            ang = np.radians(d[key][sel])
            xy = centers[sel]
            u = ell * np.cos(ang)
            v = ell * np.sin(ang)
            ax.quiver(xy[:, 0], xy[:, 1], u, v, angles='xy',
                      scale_units='xy', scale=1.0, pivot='mid',
                      width=0.0022, headwidth=3.0, headlength=5.0,
                      headaxislength=4.5, color='k', alpha=0.8, zorder=4)
            #   reference arrow: full length 0.8 m corresponds to comp=p95
            x0, x1 = ax.get_xlim()
            y0, y1 = ax.get_ylim()
            ax.quiver(x0 + 0.05 * (x1 - x0), y0 + 0.09 * (y1 - y0),
                      0.8, 0.0, angles='xy', scale_units='xy', scale=1.0,
                      pivot='tail', width=0.003, headwidth=3.0,
                      headlength=5.0, headaxislength=4.5,
                      color='#d62728', zorder=6)
            ax.text(x0 + 0.05 * (x1 - x0), y0 + 0.13 * (y1 - y0),
                    'L_ref=0.8 m (comp=p95); every %dth' % stride,
                    fontsize=7, color='#d62728')
    fig.suptitle('LF potential (u3=0) - direction diagnostic; '
                 'phi provisional; max|dtheta|=%.2f deg (state=wrinkled)'
                 % max_delta, fontsize=10)
    _save(fig, outdir, 'fig3_dir.png', dpi)


def _fig4(meta_lf, data_lf, z_lf, meta_hf, data_hf, z_hf, outdir, dpi,
          clip_q, xc):
    nc_lf, verts_lf, d_lf = _polygons(z_lf, data_lf, use_def=False)
    nc_hf, verts_hf, d_hf = _polygons(z_hf, data_hf, use_def=False)
    fig, axs = plt.subplots(1, 2, figsize=(12.0, 4.8))
    comp = d_lf['comp']
    wr = d_lf['state'] == 1
    vmax_lf = _vmax(comp[wr] if np.any(wr) else comp, clip_q)
    pcl = _collection(axs[0], verts_lf, comp, 'viridis',
                      Normalize(0.0, vmax_lf))
    _axes(axs[0], nc_lf, 'LF state (LF ODB, u3=0) - I_w integrand')
    cbl = fig.colorbar(pcl, ax=axs[0], fraction=0.035, pad=0.02)
    cbl.set_label('comp (strain, p%.1f)' % clip_q)
    loss = d_hf['loss_e']
    vmax_hf = _vmax(loss, clip_q)
    pch = _collection(axs[1], verts_hf, loss, 'magma',
                      Normalize(0.0, vmax_hf))
    _axes(axs[1], nc_hf, 'HF loss (HF ODB state) - reference layout')
    cbh = fig.colorbar(pch, ax=axs[1], fraction=0.035, pad=0.02)
    cbh.set_label('loss_e (relative, p%.1f)' % clip_q)
    _clamp_overlay(axs[0], xc)
    _clamp_overlay(axs[1], xc)
    fig.suptitle('LF potential vs HF loss - reference element layout', fontsize=11)
    fig.subplots_adjust(left=0.05, right=0.93, bottom=0.20, top=0.88,
                        wspace=0.28)
    fig.text(0.5, 0.035,
             'State/field in each panel is computed from its own ODB; '
             'element labels align but LF and HF states can differ. '
             'Cyan circles/lines: clamp patches/cables; black triangles: vertices.',
             ha='center', fontsize=8)
    _save(fig, outdir, 'fig4_lf_hf.png', dpi)


def _fig4b(meta_hf, data_hf, z_hf, outdir, dpi, clip_q, z_exag):
    nc, verts, d = _polygons(z_hf, data_hf, use_def=True, three_d=True)
    loss = d['loss_e']
    vmax = _vmax(loss, clip_q)
    verts3 = []
    for v in verts:
        w = np.array(v, copy=True, dtype=np.float64)
        w[:, 2] *= float(z_exag)
        verts3.append(w)
    fig = plt.figure(figsize=(8.0, 5.0))
    ax = fig.add_subplot(111, projection='3d')
    pc = Poly3DCollection(verts3, cmap='magma', norm=Normalize(0.0, vmax),
                          edgecolors='none', linewidths=0.0)
    pc.set_array(loss)
    ax.add_collection3d(pc)
    xr = [float(np.nanmin(nc[:, 0])), float(np.nanmax(nc[:, 0]))]
    yr = [float(np.nanmin(nc[:, 1])), float(np.nanmax(nc[:, 1]))]
    zr = [float(np.nanmin(nc[:, 2])) * z_exag,
          float(np.nanmax(nc[:, 2])) * z_exag]
    if zr[1] - zr[0] <= 1e-12:
        zr = [zr[0] - 0.5, zr[0] + 0.5]
    ax.set_xlim(*xr)
    ax.set_ylim(*yr)
    ax.set_zlim(*zr)
    ax.set_box_aspect((xr[1] - xr[0], yr[1] - yr[0],
                       max(zr[1] - zr[0], 0.4)))
    ax.view_init(elev=28, azim=-60)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('u3 x %g [m]' % z_exag)
    ax.set_title('HF loss (HF ODB state, deformed)')
    cb = fig.colorbar(pc, ax=ax, fraction=0.030, pad=0.10)
    cb.set_label('loss_e (relative, p%.1f)' % clip_q)
    _save(fig, outdir, 'fig4b_hf_deformed.png', dpi)


def _fig5(meta, data, z, outdir, dpi):
    loss = np.asarray(data['loss_e'], dtype=np.float64)
    m = np.isfinite(loss) & (loss >= 0.0)
    if not np.any(m):
        print('[plot] fig5 skipped: no finite loss_e')
        return
    #   [2026-10-26] x축을 **면적 기준**으로 바꾼다. 요소 개수 기준은 메쉬를 정련하면
    #     같은 fold 띠가 4배 많은 요소를 차지해 비중이 자동으로 희석된다 -> 메쉬 간
    #     비교가 불가능하다. 면적(evol_ref) 기준이 메쉬 불변 지표다.
    idx = np.flatnonzero(m)
    order = idx[np.argsort(-loss[idx])]
    s = loss[order]
    total = float(np.sum(s))
    if total <= 0.0:
        print('[plot] fig5 skipped: zero total loss')
        return
    cum = np.cumsum(s) / total * 100.0
    w = np.asarray(data['evol_ref'], dtype=np.float64)[order]
    w = np.maximum(w, 0.0)
    if float(np.sum(w)) > 0.0:
        frac = np.cumsum(w) / float(np.sum(w)) * 100.0
        xlabel = 'top X% of AREA (sorted by loss_e)'
    else:                                   # 면적 가중이 없으면 요소 기준 폴백
        frac = np.arange(1, s.size + 1, dtype=np.float64) / s.size * 100.0
        xlabel = 'top X% elements sorted by loss_e'
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(frac, cum, 'k-', lw=1.6, label='area-weighted (mesh-invariant)')
    _el = np.arange(1, s.size + 1, dtype=np.float64) / s.size * 100.0
    ax.plot(_el, cum, ':', color='0.55', lw=1.0, label='element count (mesh-dependent)')
    for q in (1.0, 5.0, 20.0):
        i = int(np.searchsorted(frac, q))
        i = max(0, min(s.size - 1, i))
        ax.plot([frac[i]], [cum[i]], 'o', color='#d62728', ms=4)
        ax.annotate('top %.0f%% area -> %.1f%% of hf' % (q, cum[i]),
                    xy=(frac[i], cum[i]), xytext=(8, -12),
                    textcoords='offset points', fontsize=8)
    ax.legend(loc='lower right', fontsize=7, frameon=False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('cumulative loss [% of hf]')
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)
    ax.set_title('HF loss accumulation (localization)')
    _save(fig, outdir, 'fig5_loss_accum.png', dpi)


def main():
    ap = argparse.ArgumentParser(description='Plot zk_map_extract outputs.')
    ap.add_argument('--lf', required=True, help='LF output prefix')
    ap.add_argument('--hf', default=None, help='HF output prefix')
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--dpi', type=int, default=200)
    ap.add_argument('--xc', type=float, default=0.5,
                    help='x_c for clamp-patch overlay in fig4')
    ap.add_argument('--clip-percentile', type=float, default=99.9,
                    help='colorbar upper percentile for comp/loss_e')
    ap.add_argument('--z-exag', type=float, default=1000.0,
                    help='u3 vertical exaggeration for fig4b')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    meta_lf, data_lf, z_lf = _load(args.lf)
    _fig1(meta_lf, data_lf, z_lf, args.outdir, args.dpi)
    _fig2(meta_lf, data_lf, z_lf, args.outdir, args.dpi,
          args.clip_percentile)
    _fig3(meta_lf, data_lf, z_lf, args.outdir, args.dpi)
    if args.hf:
        meta_hf, data_hf, z_hf = _load(args.hf)
        _fig4(meta_lf, data_lf, z_lf, meta_hf, data_hf, z_hf,
              args.outdir, args.dpi, args.clip_percentile, args.xc)
        _fig4b(meta_hf, data_hf, z_hf, args.outdir, args.dpi,
               args.clip_percentile, args.z_exag)
        _fig5(meta_hf, data_hf, z_hf, args.outdir, args.dpi)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
