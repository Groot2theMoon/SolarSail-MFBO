# -*- coding: utf-8 -*-
import csv, io, os, shutil
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VAR = os.path.join(HERE, 'results', 'variants', 'uq_s2_strong')
tag = 's2strong_xc0p2_dc1p0_d2'
src_odb = os.path.join(HERE, 'results', 'odb', 'HF_xc0p2_dc1_u3free.odb')
src_res = os.path.join(HERE, 'results', 'results_xc0p2_dc1.txt')
shutil.copy2(src_odb, os.path.join(VAR, tag + '.odb'))
shutil.copy2(src_res, os.path.join(VAR, tag + '.txt'))
hf = None
for line in io.open(src_res, encoding='utf-8', errors='replace'):
    if line.startswith('hf='):
        hf = float(line.strip().split('=', 1)[1])
p = os.path.join(VAR, 's2_strong_runs.csv')
rows = list(csv.DictReader(io.open(p, encoding='utf-8', errors='replace')))
fields = list(rows[0].keys())
for r in rows:
    if r['tag'] == tag:
        r['hf'] = '%.16e' % hf
        r['results'] = os.path.join(VAR, tag + '.txt')
        r['odb'] = os.path.join(VAR, tag + '.odb')
with io.open(p, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader(); w.writerows(rows)
print('updated', tag, hf)
