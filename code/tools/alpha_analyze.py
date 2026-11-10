# -*- coding: utf-8 -*-
"""Analyze the alpha/clock sweep and draw summary plots (no solves)."""
from __future__ import print_function

import argparse
import csv
import io
import json
import math
import os
import sys

import numpy as np


def fnum(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except Exception:
        return None


def rankdata(a):
    a = np.asarray(a, dtype=float)
    order = np.argsort(a, kind='mergesort')
    ranks = np.empty(len(a), dtype=float)
    ranks[order] = np.arange(1, len(a) + 1, dtype=float)
    # average ties
    vals, inv, cnt = np.unique(a, return_inverse=True, return_counts=True)
    if np.any(cnt > 1):
        sums = np.zeros(len(vals), dtype=float)
        for i, r in enumerate(ranks):
            sums[inv[i]] += r
        for j in range(len(vals)):
            ranks[inv == j] = sums[j] / cnt[j]
    return ranks


def spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3:
        return None
    rx, ry = rankdata(x), rankdata(y)
    rx -= rx.mean()
    ry -= ry.mean()
    den = math.sqrt(float(np.dot(rx, rx) * np.dot(ry, ry)))
    return float(np.dot(rx, ry) / den) if den > 0 else None


def load(path):
    rows = []
    with io.open(path, encoding='utf-8', errors='replace') as f:
        for r in csv.DictReader(f):
            if r.get('error'):
                continue
            vals = {k: fnum(r.get(k)) for k in
                    ('alpha_deg', 'clock_deg', 'torque_nm', 'tau_x_nm', 'tau_y_nm',
                     'hf_alpha', 'hf_sensored', 'cp_x_m', 'cp_y_m', 'cp_offset_m',
                     'n_backward_frac')}
            if vals['alpha_deg'] is None or vals['clock_deg'] is None or vals['torque_nm'] is None:
                continue
            vals.update({'tag': r['tag'], 'case': r.get('case', ''), 'odb': r['odb']})
            rows.append(vals)
    return rows


def pct_spread(v):
    v = np.asarray([x for x in v if x is not None and math.isfinite(x)], dtype=float)
    if len(v) < 3:
        return None, None, len(v)
    p50 = np.percentile(v, 50)
    p90 = np.percentile(v, 90)
    q1, q3 = np.percentile(v, 25), np.percentile(v, 75)
    r1 = p90 / p50 if p50 != 0 else None
    r2 = (q3 - q1) / p50 if p50 != 0 else None
    return r1, r2, len(v)


def series_by(rows, alpha, phi, key='torque_nm'):
    return {r['tag']: r[key] for r in rows
            if abs(r['alpha_deg'] - alpha) < 1e-9 and abs(r['clock_deg'] - phi) < 1e-9
            and r.get(key) is not None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', default='results/alpha_sweep/alpha_sweep.csv')
    ap.add_argument('--outdir', default='results/alpha_sweep')
    a = ap.parse_args()
    rows = load(a.csv)
    if not rows:
        raise SystemExit('no valid rows')
    if not os.path.isdir(a.outdir):
        os.makedirs(a.outdir)
    alphas = sorted(set(r['alpha_deg'] for r in rows))
    phis = sorted(set(r['clock_deg'] for r in rows))
    tags = sorted(set(r['tag'] for r in rows))
    summary = {'n_rows': len(rows), 'n_designs': len(tags),
               'alphas': alphas, 'phis': phis, 'spread': {}, 'rank': {},
               'clock': {}, 'monotonic': {}, 'hf': {}}

    spread_rows = []
    for phi in phis:
        for alpha in alphas:
            vals = list(series_by(rows, alpha, phi).values())
            r1, r2, n = pct_spread(vals)
            summary['spread'].setdefault(str(phi), {})[str(alpha)] = {
                'p90_over_p50': r1, 'iqr_over_p50': r2, 'n': n,
                'p25_tau': float(np.percentile(vals, 25)) if vals else None,
                'median_tau': float(np.median(vals)) if vals else None,
                'p75_tau': float(np.percentile(vals, 75)) if vals else None,
                'p90_tau': float(np.percentile(vals, 90)) if vals else None}
            spread_rows.append((phi, alpha, r1, r2, n))
        s0 = summary['spread'][str(phi)].get(str(alphas[0]), {})
        s10 = summary['spread'][str(phi)].get(str(alphas[-1]), {})
        summary['spread'][str(phi)]['ratio_10_over_0_p90p50'] = (
            (s10.get('p90_over_p50') / s0.get('p90_over_p50'))
            if s0.get('p90_over_p50') and s10.get('p90_over_p50') else None)

    # rankings vs alpha=0
    for phi in phis:
        base = series_by(rows, alphas[0], phi)
        for alpha in alphas:
            cur = series_by(rows, alpha, phi)
            common = sorted(set(base) & set(cur))
            rho = spearman([base[t] for t in common], [cur[t] for t in common]) if len(common) >= 3 else None
            summary['rank'].setdefault(str(phi), {})[str(alpha)] = {
                'spearman_vs_alpha0': rho, 'n': len(common)}

    # clock dependence at each alpha
    for alpha in alphas:
        s0 = series_by(rows, alpha, phis[0])
        s90 = series_by(rows, alpha, phis[-1]) if len(phis) > 1 else {}
        common = sorted(set(s0) & set(s90))
        rho = spearman([s0[t] for t in common], [s90[t] for t in common]) if len(common) >= 3 else None
        ratio = None
        if common:
            ratio = float(np.median([s90[t] / s0[t] for t in common if s0[t] != 0]))
        summary['clock'][str(alpha)] = {
            'spearman_phi0_vs_phi90': rho,
            'median_tau_ratio_90_over_0': ratio,
            'n': len(common)}

    # monotonic alpha per design/phi
    for phi in phis:
        by_tag = {}
        for alpha in alphas:
            for t, v in series_by(rows, alpha, phi).items():
                by_tag.setdefault(t, {})[alpha] = v
        pos = neg = non = 0
        for t, d in by_tag.items():
            xs = [a for a in alphas if a in d]
            if len(xs) < 3:
                continue
            rho = spearman(xs, [d[a] for a in xs])
            if rho is None:
                non += 1
            elif rho > 0.999:
                pos += 1
            elif rho < -0.999:
                neg += 1
            else:
                non += 1
        summary['monotonic'][str(phi)] = {'increasing': pos, 'decreasing': neg,
                                          'non_monotonic': non,
                                          'n_designs': len(by_tag)}

    # hf spread
    for phi in phis:
        for alpha in alphas:
            vals = [r['hf_sensored'] for r in rows
                    if abs(r['alpha_deg'] - alpha) < 1e-9
                    and abs(r['clock_deg'] - phi) < 1e-9
                    and r.get('hf_sensored') is not None]
            if vals:
                summary['hf'].setdefault(str(phi), {})[str(alpha)] = {
                    'median': float(np.median(vals)),
                    'p90': float(np.percentile(vals, 90)),
                    'max': float(np.max(vals))}

    with io.open(os.path.join(a.outdir, 'alpha_summary.json'), 'w', encoding='utf-8') as f:
        f.write(json.dumps(summary, ensure_ascii=False, indent=1))
    with io.open(os.path.join(a.outdir, 'alpha_spread.csv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f)
        w.writerow(['clock_deg', 'alpha_deg', 'p90_over_p50', 'iqr_over_p50', 'n_designs'])
        w.writerows(spread_rows)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception as e:
        print('matplotlib unavailable: %s' % e)
        print(json.dumps(summary, ensure_ascii=False, indent=1))
        return 0

    colors = {phis[0]: '#1f77b4', phis[-1]: '#d62728'}
    # 1) spread
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
    for phi in phis:
        vals = [summary['spread'][str(phi)][str(a)] for a in alphas]
        axs[0].plot(alphas, [v['p90_over_p50'] for v in vals], 'o-', label='phi=%g' % phi,
                    color=colors.get(phi))
        p25 = [v['p25_tau'] for v in vals]
        p50 = [v['median_tau'] for v in vals]
        p75 = [v['p75_tau'] for v in vals]
        p90 = [v['p90_tau'] for v in vals]
        axs[1].fill_between(alphas, p25, p75, color=colors.get(phi), alpha=0.15)
        axs[1].plot(alphas, p50, 'o-', label='phi=%g median' % phi, color=colors.get(phi))
        axs[1].plot(alphas, p25, ':', color=colors.get(phi), lw=1)
        axs[1].plot(alphas, p90, '--', color=colors.get(phi), lw=1)
    axs[0].set_xlabel('alpha [deg]'); axs[0].set_ylabel('P90/P50 tau')
    axs[1].set_xlabel('alpha [deg]'); axs[1].set_ylabel('|tau| [N m] (P25/P50/P90)')
    for ax in axs:
        ax.grid(alpha=0.3)
        ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, 'alpha_spread.svg'))
    fig.savefig(os.path.join(a.outdir, 'alpha_spread.png'), dpi=160)
    plt.close(fig)

    # 2) rank correlation
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    for phi in phis:
        rhos = [summary['rank'][str(phi)][str(a)]['spearman_vs_alpha0'] for a in alphas]
        ax.plot(alphas, rhos, 'o-', label='phi=%g' % phi, color=colors.get(phi))
    ax.axhline(0.9, color='gray', ls='--', lw=1, label='0.9')
    ax.axhline(0.7, color='gray', ls=':', lw=1, label='0.7')
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel('alpha [deg]'); ax.set_ylabel('Spearman vs alpha=0')
    ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, 'alpha_rank.svg'))
    fig.savefig(os.path.join(a.outdir, 'alpha_rank.png'), dpi=160)
    plt.close(fig)

    # 3) clock scatter at largest alpha
    if len(phis) > 1:
        s0 = series_by(rows, alphas[-1], phis[0])
        s90 = series_by(rows, alphas[-1], phis[-1])
        common = sorted(set(s0) & set(s90))
        fig, ax = plt.subplots(figsize=(4.8, 4.4))
        ax.scatter([s0[t] for t in common], [s90[t] for t in common], s=24)
        lo = min([s0[t] for t in common] + [s90[t] for t in common])
        hi = max([s0[t] for t in common] + [s90[t] for t in common])
        ax.plot([lo, hi], [lo, hi], 'k--', lw=1)
        ax.set_xlabel('tau phi=0 [N m]'); ax.set_ylabel('tau phi=90 [N m]')
        ax.set_title('alpha=%g deg' % alphas[-1])
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, 'alpha_clock.svg'))
        fig.savefig(os.path.join(a.outdir, 'alpha_clock.png'), dpi=160)
        plt.close(fig)

    # 4) tau vs alpha
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    for phi in phis:
        for tag in tags:
            ys = []
            for alpha in alphas:
                d = series_by(rows, alpha, phi)
                ys.append(d.get(tag, np.nan))
            ax.plot(alphas, ys, '-', lw=0.7, color=colors.get(phi), alpha=0.45)
    ax.set_xlabel('alpha [deg]'); ax.set_ylabel('|tau| [N m]')
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, 'alpha_tau.svg'))
    fig.savefig(os.path.join(a.outdir, 'alpha_tau.png'), dpi=160)
    plt.close(fig)

    print(json.dumps(summary, ensure_ascii=False, indent=1))
    print('wrote', a.outdir)
    return 0


if __name__ == '__main__':
    sys.exit(main())
