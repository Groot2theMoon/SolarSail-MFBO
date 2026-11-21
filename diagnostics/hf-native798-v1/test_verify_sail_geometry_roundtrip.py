import json
from pathlib import Path
import tempfile
from types import SimpleNamespace as Item
import unittest

import numpy as np

from verify_sail_geometry_roundtrip import compare_arrays, energies, identity, membrane, native_field


class RoundTripTests(unittest.TestCase):
    def test_label_join_keeps_ordered_midsides_signs_and_zero_bound(self):
        nodes = {n: (float(n), 0., 0.) for n in range(1, 9)}
        topo = {2: ('S8R5', (1,2,3,4,5,6,7,8)), 3: ('STRI65', (1,2,3,5,6,7))}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            labels = sorted(nodes, reverse=True)
            values = dict(node_labels=np.array(labels, dtype=np.int64),
                          node_xyz_m=np.array([nodes[n] for n in labels]),
                          element_labels_S8R5=np.array([2], dtype=np.int64),
                          connectivity_S8R5=np.array([topo[2][1]], dtype=np.int64),
                          face_sign_S8R5=np.array([1], dtype=np.int8),
                          element_labels_STRI65=np.array([3], dtype=np.int64),
                          connectivity_STRI65=np.array([topo[3][1]], dtype=np.int64),
                          face_sign_STRI65=np.array([-1], dtype=np.int8))
            def check(mutation):
                np.savez(root/'geometry.npz', **dict(values, **mutation))
                npz = identity(root/'geometry.npz')
                npz['path'] = 'geometry.npz'
                (root/'geometry.json').write_text(json.dumps(dict(npz=npz, node_count=8, element_counts={'S8R5':1,'STRI65':1})))
                return compare_arrays(root, nodes, topo, {2:1,3:-1})
            self.assertTrue(check({})['bitwise_equal_including_signed_zero'])
            displaced = values['node_xyz_m'].copy()
            displaced[0,2] = 1.e-30
            with self.assertRaisesRegex(ValueError, 'mapping roundoff'):
                check({'node_xyz_m': displaced})
            with self.assertRaisesRegex(ValueError, 'topology'):
                check({'connectivity_S8R5': np.array([[1,2,3,4,6,5,7,8]], dtype=np.int64)})
            with self.assertRaisesRegex(ValueError, 'face sign'):
                check({'face_sign_STRI65': np.array([1], dtype=np.int8)})

    def test_native_position_owner_duplicate_and_both_precision_attributes(self):
        def v(owner, p, n, precision='SINGLE_PRECISION', **kw):
            return Item(instance=Item(name=owner),position=p,nodeLabel=n,precision=precision,**kw)
        vals = [v('other','NODAL',1,data=(7.,7.,7.)),v('mem','INTEGRATION_POINT',2,data=(9.,9.,9.)),
                v('mem','NODAL',1,data=(1.,2.,3.)),v('mem','NODAL',2,'DOUBLE_PRECISION',dataDouble=(4.,5.,6.))]
        self.assertEqual(native_field(Item(values=vals),'mem',{1,2}),
                         ({1:(1.,2.,3.),2:(4.,5.,6.)},{'SINGLE_PRECISION':1,'DOUBLE_PRECISION':1}))
        with self.assertRaisesRegex(ValueError,'coverage'):
            native_field(Item(values=vals[:-1]),'mem',{1,2})
        with self.assertRaisesRegex(ValueError,'duplicate'):
            native_field(Item(values=vals+[vals[-1]]),'mem',{1,2})

    def test_energy_ratios_require_a_unique_finite_matching_time_and_nonzero_denominator(self):
        h = Item(description='whole model', historyOutputs={
            'ALLIE':Item(data=[(0.,0.),(.5,4.),(.7,9.),(1.,10.)]),
            'ALLSD':Item(data=[(0.,0.),(.5,.2),(.7,2.),(.7,3.)]),
            'ALLWK':Item(data=[(0.,0.),(1.,8.)])})
        actual = energies(Item(historyRegions={'whole':h}))
        self.assertEqual(actual['missing'], ['ALLAE'])
        a,b,c,d = actual['ALLSD_over_ALLIE']['whole']
        self.assertEqual((a['ratio'],a['reason']), (None,'ALLIE is zero'))
        self.assertEqual((b['time'],b['ratio']),(.5,.05))
        self.assertIsNone(c['ratio'])
        self.assertIsNone(d['ratio'])

    def test_only_empty_instance_placement_can_claim_identity(self):
        deck = ('*Part, name=Membrane\n*Node\n' + ''.join(f'{i}, {i}, 0, 0\n' for i in range(1,7)) +
                '*Element, type=STRI65\n1,1,2,3,4,5,6\n*End Part\n*Instance, name=MEMBRANE-1, part=MEMBRANE\n')
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'candidate.inp'
            p.write_text(deck+'*End Instance\n')
            self.assertEqual(membrane(p,'MEMBRANE-1')[1],{1:('STRI65',(1,2,3,4,5,6))})
            p.write_text(deck+'0., 0., -1.\n*End Instance\n')
            with self.assertRaisesRegex(ValueError,'untransformed'):
                membrane(p,'MEMBRANE-1')


if __name__ == '__main__':
    unittest.main()
