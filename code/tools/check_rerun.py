# -*- coding: utf-8 -*-
"""Compare re-run 25-point results against the archived points.csv."""
from __future__ import print_function

import csv
import io
import math
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def read_result(path):
    d = {}
    if not os.path.exists(path):
        return d
    with io.open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if '=' in line:
                k, v = line.strip().split('=', 1)
                try:
                    d[k] = float(v)
                except Exception:
                    d[k] = v
    return d


def main():
    pts = os.path.join(HERE, 'results', 'decision_grid', 'points.csv')
    rows = list(csv.DictReader(io.open(pts, encoding='utf-8', errors='replace')))
    diffs = []
    missing = []
    for r in rows:
        tag = r['tag']
        res = read_result(os.path.join(HERE, 'results', 'results_%s.txt' % tag))
        if 'hf' not in res:
            missing.append(tag)
            continue
        old = float(r['hf'])
        new = float(res['hf'])
        rel = abs(new - old) / abs(old) if abs(old) > 1e-30 else abs(new - old)
        diffs.append((tag, old, new, rel))
    diffs.sort(key=lambda x: x[3], reverse=True)
    print('n=%d missing=%d' % (len(rows), len(missing)))
    if missing:
        print('missing:', missing)
    for tag, old, new, rel in diffs[:10]:
        print('%-22s old=%.10e new=%.10e rel=%.3e' % (tag, old, new, rel))
    if diffs:
        print('max rel=%.3e  median rel=%.3e' % (diffs[0][3],
              sorted(x[3] for x in diffs)[len(diffs) // 2]))
    return 0


if __name__ == '__main__':
    sys.exit(main())
