# -*- coding: utf-8 -*-
"""Generate a subgrid RVE INP (S4, 2 mm, wavy initial geometry with s0).

Run with normal Python 3 (not Abaqus):
    python tools/subgrid_inp.py --s0 0.01 --mode biaxial --seed 0
"""
import argparse
import math
import os

L = 0.5
H = 0.002
NX = int(round(L / H))
NY = NX
LAM = L / 5.0
K = 2.0 * math.pi / LAM
THICK = 5.0e-6
E = 2.5e9
NU = 0.34
RHO = 1420.0


def amp(s0):
    return math.sqrt(2.0 * s0) / K


def phase(seed):
    return {0: (0.0, 0.0), 1: (0.0, math.pi / 2.0), 2: (math.pi / 2.0, 0.0)}.get(seed % 3, (0.0, 0.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--s0', type=float, default=0.01)
    ap.add_argument('--mode', default='biaxial', choices=('biaxial', 'uniaxial', 'shear'))
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--eps', type=float, default=0.005)
    ap.add_argument('--outdir', default='subgrid')
    a = ap.parse_args()
    if not os.path.isdir(a.outdir):
        os.makedirs(a.outdir)
    name = 'RVE_s0_%s_%s_s%d' % (('%g' % a.s0).replace('.', 'p'), a.mode, a.seed)
    path = os.path.join(a.outdir, name + '.inp')
    A = amp(a.s0)
    ph = phase(a.seed)

    with open(path, 'w', encoding='ascii') as f:
        f.write('*HEADING\nRVE s0=%g mode=%s seed=%d\n' % (a.s0, a.mode, a.seed))
        f.write('*NODE\n')
        for j in range(NY + 1):
            y = j * H
            for i in range(NX + 1):
                x = i * H
                z = A * (math.sin(K * x + ph[0]) + math.sin(K * y + ph[1]))
                f.write('%d, %.12e, %.12e, %.12e\n' % (j * (NX + 1) + i + 1, x, y, z))
        f.write('*ELEMENT, TYPE=S4, ELSET=ALL\n')
        for j in range(NY):
            for i in range(NX):
                n0 = j * (NX + 1) + i + 1
                f.write('%d, %d, %d, %d, %d\n' % (j * NX + i + 1, n0, n0 + 1,
                                                  n0 + 1 + NX + 1, n0 + NX + 1))
        for nm, pred in (
            ('LEFT', lambda x, y: abs(x) < 1e-9),
            ('RIGHT', lambda x, y: abs(x - L) < 1e-9),
            ('BOTTOM', lambda x, y: abs(y) < 1e-9),
            ('TOP', lambda x, y: abs(y - L) < 1e-9),
        ):
            ids = []
            for j in range(NY + 1):
                for i in range(NX + 1):
                    x, y = i * H, j * H
                    if pred(x, y):
                        ids.append(j * (NX + 1) + i + 1)
            f.write('*NSET, NSET=%s\n' % nm)
            for q in range(0, len(ids), 20):
                f.write(', '.join(str(v) for v in ids[q:q + 20]) + '\n')
        f.write('*MATERIAL, NAME=KAPTON\n*ELASTIC\n%g, %g\n*DENSITY\n%g\n' % (E, NU, RHO))
        f.write('*SHELL SECTION, ELSET=ALL, MATERIAL=KAPTON, OFFSET=0.0\n%g\n' % THICK)
        f.write('*BOUNDARY\n')
        if a.mode in ('biaxial', 'uniaxial'):
            f.write('LEFT, 1, 1, 0.0\nRIGHT, 1, 1, 0.0\nBOTTOM, 2, 2, 0.0\n')
            if a.mode == 'biaxial':
                f.write('TOP, 2, 2, 0.0\n')
        else:
            f.write('BOTTOM, 1, 1, 0.0\nBOTTOM, 2, 2, 0.0\nLEFT, 2, 2, 0.0\nRIGHT, 2, 2, 0.0\n')
        f.write('*STEP, NAME=Step-Tension, NLGEOM=YES\n*STATIC\n0.02, 1.0, 1e-8, 0.05\n')
        f.write('*BOUNDARY, OP=MOD\n')
        if a.mode == 'biaxial':
            f.write('RIGHT, 1, 1, %.12e\nTOP, 2, 2, %.12e\n' % (a.eps * L, a.eps * L))
        elif a.mode == 'uniaxial':
            f.write('RIGHT, 1, 1, %.12e\n' % (a.eps * L))
        else:
            f.write('TOP, 1, 1, %.12e\n' % (a.eps * L))
        f.write('*OUTPUT, FIELD, NUMBER INTERVAL=1\n')
        f.write('*NODE OUTPUT\nCOORD, U\n*ELEMENT OUTPUT\nS, LE\n')
        f.write('*END STEP\n')
    print('wrote', path, 'A=%.6e' % A)


if __name__ == '__main__':
    main()
