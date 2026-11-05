# -*- coding: utf-8 -*-
"""Gate and analyze results/cable_forces/cable_forces.csv.

G1 initial S11: main = CABLE_PRESTRESS, clamp = d_c * CABLE_PRESTRESS (1%).
G2 symmetry: RIGHT/LEFT and CL/CR final force ratios in [0.95, 1.05].
G3 equilibrium: total RF residual < 1e-6 (global); cable-only vector
   residual is reported separately (the membrane/BC carries the rest).
"""
import csv
import math
import os
import statistics
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CABLE = os.path.join(ROOT, 'results/cable_forces/cable_forces.csv')
POINTS = os.path.join(ROOT, 'results/decision_grid/points.csv')
ALL = os.path.join(ROOT, 'results/all_results.csv')
PRESTRESS = 1.71e7
MAIN = {'INST_CABLE_TOP', 'INST_CABLE_RIGHT', 'INST_CABLE_LEFT'}


def _f(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None


def _rank(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    ranks = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        r = 0.5 * (i + j) + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = r
        i = j + 1
    return ranks


def spearman(x, y):
    pairs = [(a, b) for a, b in zip(x, y) if a is not None and b is not None]
    if len(pairs) < 3:
        return float('nan')
    a = [p[0] for p in pairs]; b = [p[1] for p in pairs]
    ra, rb = _rank(a), _rank(b)
    ma, mb = statistics.mean(ra), statistics.mean(rb)
    num = sum((u - ma) * (v - mb) for u, v in zip(ra, rb))
    da = math.sqrt(sum((u - ma) ** 2 for u in ra))
    db = math.sqrt(sum((v - mb) ** 2 for v in rb))
    return num / (da * db) if da * db else float('nan')


def main():
    with open(CABLE, newline='') as f:
        rows = list(csv.DictReader(f))
    with open(POINTS, newline='') as f:
        points = {r['tag']: r for r in csv.DictReader(f)}
    lf = {}
    if os.path.exists(ALL):
        with open(ALL, newline='') as f:
            for r in csv.DictReader(f):
                if r.get('fidelity') != 'LF':
                    continue
                x, d = _f(r.get('x_c')), _f(r.get('d_c'))
                l1, l3 = _f(r.get('lf1')), _f(r.get('lf3'))
                if x is not None and d is not None and l1 is not None:
                    lf[(round(x, 8), round(d, 8))] = (l1, l3)

    by_tag = defaultdict(dict)
    for r in rows:
        by_tag[r['tag']][r['inst']] = r

    g1_bad, g1_max = [], 0.0
    for tag, insts in by_tag.items():
        p = points.get(tag, {})
        dc = _f(p.get('d_c')) or 0.0
        for inst, r in insts.items():
            s0 = _f(r.get('S11_init_Pa'))
            exp = PRESTRESS if inst in MAIN else dc * PRESTRESS
            if s0 is None or exp <= 0:
                g1_bad.append((tag, inst, s0, exp)); continue
            err = abs(s0 - exp) / exp
            g1_max = max(g1_max, err)
            if err > 0.01:
                g1_bad.append((tag, inst, s0, exp, err))
    print('G1 initial S11: max_rel_err=%.3e bad=%d' % (g1_max, len(g1_bad)))
    for b in g1_bad[:5]:
        print('  bad', b)

    ratios = {'RIGHT/LEFT': [], 'CL/CR': []}
    vec_res = []
    rf_res = []
    for tag, insts in by_tag.items():
        def fv(k):
            r = insts.get(k)
            return _f(r.get('F_N')) if r else None
        rl, rr = fv('INST_CABLE_LEFT'), fv('INST_CABLE_RIGHT')
        cl, cr = fv('INST_CABLE_CL'), fv('INST_CABLE_CR')
        if rl and rr:
            ratios['RIGHT/LEFT'].append(rr / rl)
        if cl and cr:
            ratios['CL/CR'].append(cr / cl)
        vx = vy = vz = 0.0
        sa = 0.0
        for r in insts.values():
            f = _f(r.get('F_N')) or 0.0
            vx += f * (_f(r.get('dx')) or 0.0)
            vy += f * (_f(r.get('dy')) or 0.0)
            vz += f * (_f(r.get('dz')) or 0.0)
            sa += abs(f)
        vec_res.append(math.sqrt(vx * vx + vy * vy + vz * vz) / sa if sa else float('nan'))
        rr = _f(next(iter(insts.values())).get('rf_residual_frac'))
        if rr is not None:
            rf_res.append(rr)

    print('G2 symmetry:')
    for k, v in ratios.items():
        print('  %-12s min=%.6f max=%.6f pass=%s'
              % (k, min(v), max(v), all(0.95 <= x <= 1.05 for x in v)))
    print('G3 global RF residual: max=%.3e pass=%s'
          % (max(rf_res), max(rf_res) < 1e-6))
    print('G3 cable-only vector residual: min=%.4f median=%.4f max=%.4f'
          % (min(vec_res), statistics.median(vec_res), max(vec_res)))

    print('\nP3 distribution of F/F_init:')
    inst_order = ['INST_CABLE_TOP', 'INST_CABLE_LEFT', 'INST_CABLE_RIGHT',
                  'INST_CABLE_CL', 'INST_CABLE_CR']
    for inst in inst_order:
        vals = [_f(by_tag[t][inst].get('F_over_init')) for t in by_tag if inst in by_tag[t]]
        vals = [v for v in vals if v is not None]
        if not vals:
            continue
        mean = statistics.mean(vals)
        cv = statistics.pstdev(vals) / abs(mean) if mean else float('nan')
        print('  %-18s min=%.6g median=%.6g max=%.6g CV=%.4f'
              % (inst, min(vals), statistics.median(vals), max(vals), cv))

    print('\nP3 Spearman(F/F_init, ...):')
    print('  %-18s %8s %8s %8s %8s %8s' % ('inst', 'x_c', 'd_c', 'lf1', 'lf3', 'zk'))
    for inst in inst_order:
        data = []
        for tag, insts in by_tag.items():
            if inst not in insts:
                continue
            p = points.get(tag, {})
            x, d = _f(p.get('x_c')), _f(p.get('d_c'))
            val = _f(insts[inst].get('F_over_init'))
            if x is None or d is None or val is None:
                continue
            pair = lf.get((round(x, 8), round(d, 8)))
            data.append((val, x, d,
                         pair[0] if pair else None,
                         pair[1] if pair else None,
                         _f(p.get('zk_state_wrinkled_frac'))))
        if not data:
            continue
        cols = list(zip(*data))
        print('  %-18s %8.3f %8.3f %8.3f %8.3f %8.3f'
              % (inst, spearman(cols[0], cols[1]), spearman(cols[0], cols[2]),
                 spearman(cols[0], cols[3]), spearman(cols[0], cols[4]),
                 spearman(cols[0], cols[5])))


if __name__ == '__main__':
    main()
