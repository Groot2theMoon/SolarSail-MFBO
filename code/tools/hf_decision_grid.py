#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""hf 결정변수 검정 - 5x5 (x_c x d_c) HF 격자를 돌리고 (A) hf 순위 vs (B) 비-hf 구조지표 순위를 비교한다.

    python tools/hf_decision_grid.py --points 0.1:1e-4,0.1:0.01 --resume
    python tools/hf_decision_grid.py --analyze-only

설계
  A = hf 오름차순(작을수록 좋다).
  B = 구조지표 3종의 **순위 평균**: slope_backward_area_pct, zk_state_wrinkled_frac,
      slope_mean_neg_frac (모두 클수록 나쁘다). hf/slope_f*_hf_* 는 hf 파생이므로 제외한다.
  비교: Spearman rho, top-k(기본 5) 선택 겹침, 평균 |Δrank|, Δrho^2 = 1-rho^2.
  무작위 기준: 같은 n 에서 순열 1000회의 평균 겹침/rho.

판정
  PASS(hf 가 결정적)  : top-k 겹침 <= k-2 이고 rho < 0.7
  FAIL(hf 무의미)     : top-k 겹침 >= k-1 이고 rho >= 0.8
  AMBIGUOUS           : 그 사이. 수치를 그대로 기록한다.

각 점은 `abaqus cae noGUI=run_abaqus.py -- LF/HF x d` 로 실행하고
`results/results_xc<tag>_dc<tag>.txt` 의 key=value 를 읽는다. `--resume` 은 hf 키가 있는 점을 건너뛴다.
"""
from __future__ import print_function

import argparse
import csv
import io
import json
import os
import random
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import aba_imperfection as A                                                # noqa: E402

KV = re.compile(r'^([A-Za-z_][\w.\-]*)\s*=\s*(.*?)\s*$')
X_C = [0.1, 0.2, 0.5, 0.7, 0.8]
D_C = [1e-4, 0.01, 0.0556, 0.1, 0.5]
STRUCT_KEYS = ('slope_backward_area_pct', 'zk_state_wrinkled_frac', 'slope_mean_neg_frac')
OUTDIR = os.path.join(HERE, 'results', 'decision_grid')


def point_paths(x_c, d_c):
    xt, dt = A.value_tag(x_c), A.value_tag(d_c)
    return (os.path.join(HERE, 'results',
                         'results_xc%s_dc%s.txt' % (xt, dt)),
            'xc%s_dc%s' % (xt, dt))


def read_kv(path):
    out = {}
    if not os.path.exists(path):
        return out
    with io.open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            m = KV.match(line.strip())
            if m:
                out[m.group(1)] = m.group(2)
    return out


def run_abq(fid, x_c, d_c):
    cmd = 'abaqus cae noGUI=run_abaqus.py -- %s %g %g' % (fid, x_c, d_c)
    print('  $ ' + cmd)
    t0 = time.time()
    p = subprocess.run(cmd, shell=True, cwd=HERE, capture_output=True,
                       text=True, encoding='utf-8', errors='replace')
    dt = time.time() - t0
    tail = (p.stdout or '')[-800:] + (p.stderr or '')[-400:]
    if p.returncode:
        print('  !!! %s rc=%d (%.1f s)' % (fid, p.returncode, dt))
        print(tail)
    return p.returncode, dt


def collect_point(x_c, d_c, resume):
    path, tag = point_paths(x_c, d_c)
    if resume:
        cur = read_kv(path)
        if cur.get('hf') not in (None, '', 'nan'):
            print('[GRID] skip %s (hf=%s)' % (tag, cur.get('hf')))
            return {'x_c': x_c, 'd_c': d_c, 'tag': tag, 'skipped': True, **cur}
    print('\n[GRID] %s  (x_c=%g d_c=%g)' % (tag, x_c, d_c))
    rc_lf, t_lf = run_abq('LF', x_c, d_c)
    if rc_lf:
        return {'x_c': x_c, 'd_c': d_c, 'tag': tag, 'error': 'LF rc=%d' % rc_lf}
    rc_hf, t_hf = run_abq('HF', x_c, d_c)
    kv = read_kv(path)
    if rc_hf:
        kv['error'] = 'HF rc=%d' % rc_hf
    kv.update({'x_c': x_c, 'd_c': d_c, 'tag': tag, 't_lf_s': round(t_lf, 1),
               't_hf_s': round(t_hf, 1)})
    print('  -> hf=%s  t=%.1f+%.1f s' % (kv.get('hf'), t_lf, t_hf))
    return kv


def _rank(values, reverse=False):
    """값 -> 순위(1=최소). 동률은 평균순위."""
    order = sorted(range(len(values)), key=lambda i: values[i], reverse=reverse)
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def spearman(a, b):
    n = len(a)
    if n < 3:
        return float('nan')
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = sum((x - ma) ** 2 for x in a) ** 0.5
    db = sum((y - mb) ** 2 for y in b) ** 0.5
    return num / (da * db) if da and db else float('nan')


def analyze(rows, k=5):
    pts = [r for r in rows if r.get('hf') not in (None, '', 'nan')
           and all(r.get(m) not in (None, '', 'nan') for m in STRUCT_KEYS)]
    if len(pts) < 3:
        print('[ANALYZE] 유효 점 부족: %d' % len(pts))
        return None
    hf = [float(r['hf']) for r in pts]
    ra = _rank(hf)                       # 작을수록 좋다
    metric_ranks = [_rank([float(r[m]) for r in pts]) for m in STRUCT_KEYS]
    rb = [sum(mr[i] for mr in metric_ranks) / len(metric_ranks) for i in range(len(pts))]
    rho = spearman(ra, rb)
    top_a = set(sorted(range(len(pts)), key=lambda i: ra[i])[:k])
    top_b = set(sorted(range(len(pts)), key=lambda i: rb[i])[:k])
    overlap = len(top_a & top_b)
    mad = sum(abs(ra[i] - rb[i]) for i in range(len(pts))) / len(pts)
    rng = random.Random(12345)
    sims = []
    for _ in range(1000):
        perm = list(range(len(pts)))
        rng.shuffle(perm)
        rp = [0.0] * len(pts)
        for pos, idx in enumerate(perm):
            rp[idx] = pos + 1.0
        sims.append((len(top_a & set(sorted(range(len(pts)), key=lambda i: rp[i])[:k])),
                     spearman(ra, rp)))
    rand_overlap = sum(s[0] for s in sims) / len(sims)
    rand_rho = sum(s[1] for s in sims) / len(sims)
    if overlap <= max(1, k - 2) and rho < 0.7:
        verdict = 'PASS (hf is decisive: selection changes)'
    elif overlap >= k - 1 and rho >= 0.8:
        verdict = 'FAIL (hf is not decisive at this scale: selection is stable)'
    else:
        verdict = 'AMBIGUOUS'
    out = {'n': len(pts), 'k': k, 'rho': rho, 'rho2': rho * rho,
           'delta_rho2': 1 - rho * rho, 'top_overlap': overlap,
           'mean_abs_drank': mad, 'random_overlap': rand_overlap,
           'random_rho': rand_rho, 'verdict': verdict,
           'topA': [pts[i]['tag'] for i in sorted(top_a)],
           'topB': [pts[i]['tag'] for i in sorted(top_b)]}
    print('\n[ANALYZE] n=%d k=%d' % (len(pts), k))
    print('  rho(A,B)=%.3f  delta_rho2=%.3f  top-overlap=%d/%d  mean|drank|=%.2f'
          % (rho, 1 - rho * rho, overlap, k, mad))
    print('  random baseline: overlap=%.2f  rho=%.3f' % (rand_overlap, rand_rho))
    print('  topA=%s' % out['topA'])
    print('  topB=%s' % out['topB'])
    print('  => %s' % verdict)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--points', default='',
                    help='x:d;x:d ... (생략하면 5x5 전체)')
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--analyze-only', action='store_true')
    ap.add_argument('--k', type=int, default=5)
    a = ap.parse_args(argv)

    os.makedirs(OUTDIR, exist_ok=True)
    if a.points:
        pairs = []
        for part in re.split(r'[;,]', a.points):
            if not part.strip():
                continue
            xs, ds = part.split(':')
            pairs.append((float(xs), float(ds)))
    else:
        pairs = [(x, d) for x in X_C for d in D_C]

    rows = []
    if not a.analyze_only:
        for x, d in pairs:
            rows.append(collect_point(x, d, a.resume))
        with io.open(os.path.join(OUTDIR, 'points.json'), 'w', encoding='utf-8') as f:
            f.write(json.dumps(rows, ensure_ascii=False, indent=1))
        cols = ['x_c', 'd_c', 'tag', 'hf', 'slope_backward_area_pct',
                'zk_state_wrinkled_frac', 'slope_mean_neg_frac', 'dn_deg',
                'energy_verdict', 'stabilize_used', 't_lf_s', 't_hf_s', 'error']
        with io.open(os.path.join(OUTDIR, 'points.csv'), 'w', encoding='utf-8',
                     newline='') as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
            w.writeheader()
            for r in rows:
                w.writerow(r)

    pj = os.path.join(OUTDIR, 'points.json')
    if os.path.exists(pj) and (a.analyze_only or not rows):
        with io.open(pj, encoding='utf-8') as f:
            rows = json.load(f)
    res = analyze(rows, a.k)
    if res:
        with io.open(os.path.join(OUTDIR, 'analysis.json'), 'w', encoding='utf-8') as f:
            f.write(json.dumps(res, ensure_ascii=False, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
