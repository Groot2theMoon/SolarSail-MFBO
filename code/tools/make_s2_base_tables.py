# -*- coding: utf-8 -*-
"""Generate the 10 baseline S2 seed tables (corner-only mode source)."""
import os
import subprocess

py = r'C:\Users\user\anaconda3\python.exe'
outdir = os.path.join('seeds', 'uq_s2')
os.makedirs(outdir, exist_ok=True)
for seed in range(10):
    out = os.path.join(outdir, 's2_base_d%d.txt' % seed)
    cmd = [py, 'tools/make_uq_seed.py', '--draw', str(seed),
           '--x-c', '0.2', '--d-c', '0.0001',
           '--seed-table', 'seeds/modes_seed.txt',
           '--clamp-table', 'seeds/modes_clamp_xc0p2_dc0p0556.txt',
           '--header-table', 'seeds/modes_xc0p2.txt',
           '--normalize-total', '--out', out]
    subprocess.check_call(cmd, stdout=subprocess.DEVNULL)
print('generated baseline tables 10')
