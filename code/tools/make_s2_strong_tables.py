# -*- coding: utf-8 -*-
"""Generate seed tables for the strong-clamp S2 extension.

d_c=0.5: x_c = 0.1, 0.2, 0.5, 0.7, 0.8
d_c=1.0: x_c = 0.2, 0.5
"""
import os
import subprocess

py = r'C:\Users\user\anaconda3\python.exe'
designs = [
    (0.1, 0.5, 'xc0p1'), (0.2, 0.5, 'xc0p2'), (0.5, 0.5, 'xc0p5'),
    (0.7, 0.5, 'xc0p7'), (0.8, 0.5, 'xc0p8'),
    (0.2, 1.0, 'xc0p2'), (0.5, 1.0, 'xc0p5'),
]
outdir = os.path.join('seeds', 'uq_s2')
os.makedirs(outdir, exist_ok=True)
for xc, dc, tag in designs:
    for seed in range(10):
        out = os.path.join(outdir, 's2strong_%s_dc%s_d%d.txt'
                           % (tag, str(dc).replace('.', 'p'), seed))
        cmd = [py, 'tools/make_uq_seed.py', '--draw', str(seed),
               '--x-c', str(xc), '--d-c', str(dc),
               '--seed-table', 'seeds/modes_seed.txt',
               '--clamp-table', 'seeds/modes_clamp_%s_dc0p0556.txt' % tag,
               '--header-table', 'seeds/modes_%s.txt' % tag,
               '--normalize-total', '--out', out]
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL)
print('generated strong tables', len(designs) * 10)
