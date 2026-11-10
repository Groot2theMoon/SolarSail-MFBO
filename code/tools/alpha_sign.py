# -*- coding: utf-8 -*-
"""Sign stability of torque components across alpha/clock, from alpha_sweep.csv.

For each design and clock:
  - sign(tau_x), sign(tau_y) at alpha=0
  - retention at alpha=10 and across all alphas
  - flips: tag/component/alpha

Run with normal Python 3:
  python tools/alpha_sign.py --csv results/alpha_sweep/alpha_sweep.csv
"""
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


def sign(x, eps=1e-12):
    if x is None or not math.isfinite(x):
        return None
    if abs(x) <= eps:
        return 0
    return 1 if x > 0 else -1


def load(path):
    out = []
    with io.open(path, encoding='utf-8', errors='replace') as f:
        for r in csv.DictReader(f):
            if r.get('error'):
                continue
            a = fnum(r.get('alpha_deg'))
            p = fnum(r.get('clock_deg'))
            tx = fnum(r.get('tau_x_nm'))
            ty = fnum(r.get('tau_y_nm'))
            if None in (a, p, tx, ty):
                continue
            out.append({'tag': r['tag'], 'case': r.get('case', ''),
                        'alpha': a, 'phi': p, 'tx': tx, 'ty': ty})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', default='results/alpha_sweep/alpha_sweep.csv')
    ap.add_argument('--out', default='results/alpha_sweep/alpha_sign.json')
    ap.add_argument('--eps', type=float, default=1e-12)
    ap.add_argument('--plot', default='results/alpha_sweep/alpha_sign.svg')
    a = ap.parse_args()
    rows = load(a.csv)
    if not rows:
        raise SystemExit('no rows')
    phis = sorted(set(r['phi'] for r in rows))
    alphas = sorted(set(r['alpha'] for r in rows))
    tags = sorted(set(r['tag'] for r in rows))
    out = {'eps': a.eps, 'alphas': alphas, 'phis': phis, 'n_designs': len(tags),
           'retention_0_to_10': {}, 'retention_all_alpha': {}, 'flips': []}

    def get(phi, alpha):
        return {r['tag']: r for r in rows
                if abs(r['phi'] - phi) < 1e-9 and abs(r['alpha'] - alpha) < 1e-9}

    for phi in phis:
        d0 = get(phi, alphas[0])
        d10 = get(phi, alphas[-1])
        common = sorted(set(d0) & set(d10))
        rec = {}
        for comp in ('tx', 'ty'):
            same = both = 0
            flips = []
            for t in common:
                s0 = sign(d0[t][comp], a.eps)
                s1 = sign(d10[t][comp], a.eps)
                if s0 == s1:
                    same += 1
                else:
                    flips.append({'tag': t, 'component': comp,
                                  'alpha_deg': alphas[-1],
                                  'sign_at_0': s0, 'sign_at_10': s1,
                                  'value_at_0': d0[t][comp], 'value_at_10': d10[t][comp]})
                if s0 == s1 and s0 != 0:
                    both += 1
            rec[comp] = {'same': same, 'n': len(common),
                         'fraction': same / float(len(common)) if common else None,
                         'both_same_nonzero': both,
                         'flips': flips}
            out['flips'].extend([dict(f, phi=phi) for f in flips])
        # both components at once
        both_same = 0
        for t in common:
            if (sign(d0[t]['tx'], a.eps) == sign(d10[t]['tx'], a.eps)
                    and sign(d0[t]['ty'], a.eps) == sign(d10[t]['ty'], a.eps)):
                both_same += 1
        rec['both_components'] = {'same': both_same, 'n': len(common),
                                  'fraction': both_same / float(len(common)) if common else None}
        out['retention_0_to_10'][str(phi)] = rec

        # all-alpha consistency vs alpha=0
        all_common = sorted(set.intersection(*[set(get(phi, al)) for al in alphas]))
        rec_all = {}
        for comp in ('tx', 'ty'):
            ok = 0
            flips = []
            for t in all_common:
                s0 = sign(get(phi, alphas[0])[t][comp], a.eps)
                bad = []
                for al in alphas[1:]:
                    s = sign(get(phi, al)[t][comp], a.eps)
                    if s != s0:
                        bad.append({'alpha_deg': al, 'sign': s})
                if bad:
                    flips.append({'tag': t, 'component': comp, 'sign_at_0': s0,
                                  'changes': bad})
                else:
                    ok += 1
            rec_all[comp] = {'same_all': ok, 'n': len(all_common),
                              'fraction': ok / float(len(all_common)) if all_common else None,
                              'flips': flips}
            out['flips'].extend([dict(f, phi=phi) for f in flips])
        out['retention_all_alpha'][str(phi)] = rec_all

    print(json.dumps(out, ensure_ascii=False, indent=1))
    d = os.path.dirname(a.out)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with io.open(a.out, 'w', encoding='utf-8') as f:
        f.write(json.dumps(out, ensure_ascii=False, indent=1))

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap
        cmap = ListedColormap(['#d62728', '#eeeeee', '#1f77b4'])  # -1,0,+1
        fig, axs = plt.subplots(2, len(phis), figsize=(4.0 * len(phis), 6.0),
                                squeeze=False)
        for j, phi in enumerate(phis):
            for i, comp in enumerate(('tx', 'ty')):
                mat = np.zeros((len(tags), len(alphas)))
                for r in range(len(tags)):
                    for c, al in enumerate(alphas):
                        mat[r, c] = sign(get(phi, al).get(tags[r], {}).get(comp), a.eps)
                ax = axs[i][j]
                ax.imshow(mat, aspect='auto', cmap=cmap, vmin=-1, vmax=1)
                ax.set_xticks(range(len(alphas)))
                ax.set_xticklabels(['%g' % x for x in alphas])
                ax.set_yticks(range(len(tags)))
                ax.set_yticklabels(tags if j == 0 else [''] * len(tags), fontsize=5)
                ax.set_title('%s  phi=%g' % (comp, phi))
                ax.set_xlabel('alpha [deg]')
        fig.tight_layout()
        fig.savefig(a.plot)
        fig.savefig(a.plot.replace('.svg', '.png'), dpi=160)
        plt.close(fig)
        print('wrote', a.plot)
    except Exception as e:
        print('plot skipped: %s' % e)
    return 0


if __name__ == '__main__':
    sys.exit(main())
