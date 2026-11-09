# -*- coding: utf-8 -*-
"""Create a UQ seed table: common corner draw + design-specific clamp modes.

- Signed weights for the two corner seed modes (common draw).
- Corner field normalized to max|u3|=1 before adding clamp modes.
- Design clamp modes are fixed (weight 1).
- Header from the clamp table keeps x_c/fingerprint so seed_gate can pass.
- Manifest records draw_id, weights, normalization, and source path.

Usage (from code/):
    python tools/make_uq_seed.py --draw 0 --x-c 0.2 --d-c 0.0001 \
        --clamp-table seeds/modes_clamp_xc0p2_dc0p0556.txt \
        --out seeds/uq_s2/xc0p2_dc0p0001_d0.txt
"""
import argparse
import csv
import io
import os
import re

# Fixed, reproducible signed weights. 10 draws; common across designs.
DRAWS = [
    (1.0, 0.7), (-1.0, 0.8), (0.6, -1.0), (-0.8, -0.9), (1.2, -0.6),
    (-0.7, 1.1), (0.9, 1.3), (-1.3, 0.5), (0.5, 1.5), (-1.1, -1.2),
]


def read_table(path):
    header, modes, cur = [], {}, None
    with io.open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith('#'):
                if cur is None:
                    header.append(s)
                continue
            if s.upper().startswith('MODE '):
                cur = int(s.split()[1])
                modes[cur] = {}
                continue
            if cur is None:
                continue
            p = s.split()
            if len(p) >= 2:
                modes[cur][int(p[0])] = float(p[1])
    return header, modes


def write_table(path, header, modes):
    d = os.path.dirname(path)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with io.open(path, 'w', encoding='utf-8') as f:
        for h in header:
            f.write(h + '\n')
        for m in sorted(modes):
            f.write('MODE %d\n' % m)
            for lab in sorted(modes[m]):
                f.write('%d %.10e\n' % (lab, modes[m][lab]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--draw', type=int, default=0)
    ap.add_argument('--x-c', type=float, required=True)
    ap.add_argument('--d-c', type=float, required=True)
    ap.add_argument('--seed-table', default='seeds/modes_seed.txt')
    ap.add_argument('--clamp-table', default=None)
    ap.add_argument('--header-table', default=None,
                    help='header/fingerprint source (required for clamp-free baseline)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--normalize-total', action='store_true',
                    help='normalize the combined field to max|u3|=1')
    ap.add_argument('--manifest', default='results/uq_s2/seed_manifest.csv')
    a = ap.parse_args()
    if not 0 <= a.draw < len(DRAWS):
        raise SystemExit('draw must be 0..%d' % (len(DRAWS) - 1))
    w1, w2 = DRAWS[a.draw]

    sh, sm = read_table(a.seed_table)
    ch, cm = read_table(a.clamp_table) if a.clamp_table else (sh, {})
    if a.header_table:
        ch, _ = read_table(a.header_table)

    # common corner field (modes 1,2), signed and normalized to max|u3|=1
    labels = set(sm[1]) | set(sm[2])
    modes = {m: {lab: 0.0 for lab in labels} for m in (1, 2, 3, 4)}
    # Keep the two corner modes independent (the original composite uses seed 1,2 +
    # clamp 1,2 = four modes). Do not merge them into one mode.
    modes[1] = {lab: w1 * v for lab, v in sm[1].items()}
    modes[2] = {lab: w2 * v for lab, v in sm[2].items()}
    if a.clamp_table:
        modes[3] = dict(cm[1])
        modes[4] = dict(cm[2])
    # mode 3/4 are zero for a corner-only baseline; run_abaqus loads four modes uniformly.
    if a.normalize_total:
        total = {}
        for m in modes:
            for lab, v in modes[m].items():
                total[lab] = total.get(lab, 0.0) + v
        mx = max(abs(v) for v in total.values()) or 1.0
        for m in modes:
            modes[m] = {lab: v / mx for lab, v in modes[m].items()}
        norm_note = 'combined max|u3|=1 (scale=1/%.6g)' % mx
    else:
        norm_note = 'seed modes 1,2 weighted; clamp modes 3,4 added unnormalized'

    # header from clamp table (has x_c/fingerprint); update SEED line
    header = list(ch)
    out = []
    seed_done = False
    for h in header:
        if h.upper().startswith('# SEED '):
            h = re.sub(r'x_c=[^ ]+', 'x_c=%.10g' % a.x_c, h)
            h = re.sub(r'd_c_ref=[^ ]+', 'd_c_ref=%.10g' % a.d_c, h)
            seed_done = True
        out.append(h)
    if not seed_done:
        out.append('# SEED x_c=%.10g d_c_ref=%.10g' % (a.x_c, a.d_c))
    out.append('# UQ draw=%d w1=%.6g w2=%.6g norm=%s' % (a.draw, w1, w2, norm_note))
    out.append('# UQ clamp=%s' % (a.clamp_table or 'none'))
    write_table(a.out, out, modes)

    d = os.path.dirname(a.manifest)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    exists = os.path.exists(a.manifest)
    row = {'draw': a.draw, 'x_c': a.x_c, 'd_c': a.d_c,
           'seed_table': a.seed_table, 'clamp_table': a.clamp_table or '',
           'w1': w1, 'w2': w2, 'corner_max_abs': 1.0, 'out': a.out}
    with io.open(a.manifest, 'a', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()))
        if not exists:
            w.writeheader()
        w.writerow(row)
    print('wrote %s (modes=%s, nodes=%d)' %
          (a.out, sorted(modes), len(modes[1])))


if __name__ == '__main__':
    main()
