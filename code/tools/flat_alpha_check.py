# -*- coding: utf-8 -*-
"""Check the alpha normalization on an undeformed (flat) ODB frame.

For a perfect flat plate, hf(alpha) must be ~0 for every alpha.
Run with abaqus python.
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from odbAccess import openOdb
import eval_abaqus as ev

path = sys.argv[1] if len(sys.argv) > 1 else 'results/odb/HF_xc0p8_dc0p0001_u3free.odb'
odb = openOdb(path=path, readOnly=True)
try:
    f0 = None
    for sn in odb.steps.keys():
        for fr in odb.steps[sn].frames:
            if abs(fr.frameValue) < 1e-12:
                f0 = fr
                break
        if f0 is not None:
            break
    if f0 is None:
        raise RuntimeError('no frameValue=0 frame')
    inst = odb.rootAssembly.instances['MEMBRANE-1']
    conns = [e.connectivity for e in inst.elements]
    max_label = (inst.nodes[-1].label + 100) if inst.nodes else 100000
    print('step frame=0  elements=%d' % len(conns))
    for alpha in (0.0, 2.5, 5.0, 7.5, 10.0):
        us = ev.sun_vector(alpha, 0.0)
        ideal = ev.ideal_thrust_coeff(alpha)
        diag = {}
        extra = {'n_total': 0, 'n_backward': 0, 'area_total': 0.0, 'area_backward': 0.0}
        ev.calc_thrust_loss(f0, inst, conns, max_label, extra_out=extra, diag=diag,
                            u_sun=us, ideal_coeff=ideal)
        print('alpha=%4.1f  ideal=%.12f  hf=%.3e  Fz=%.6e  n_backward=%.3e'
              % (alpha, ideal, diag.get('hf_signed', float('nan')),
                 diag.get('force_z_n', float('nan')),
                 (extra['area_backward'] / extra['area_total']
                  if extra['area_total'] else float('nan'))))
finally:
    odb.close()
