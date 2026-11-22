"""Fabricated ordinary-Python fixtures for the one shared seed recipe."""

import math
from pathlib import Path
import tempfile
import unittest

from prepare_paired_seedb import edit, invert, pin_all, reflection, rows, sha, verify_coverage, verify_edit


DECK = (b'*Part, name=MEMBRANE\r\n*Node\r\n9, 0., 0., 2e-9\r\n'
        b'5, 2., 0., -7e-10\r\n3, 1., 0.5, 1.2e-6\r\n7, 1., 1., -6e-7\r\n'
        b'*End Part\r\n*Initial Conditions, type=stress\r\n'
        b'Inst_cable_CL.All_Elem, 171000., 0., 0.,\r\n'
        b'Inst_cable_Left.All_Elem, 1.71e+07, 0., 0.,\r\n*Step, name=Step-GlobalTension\r\n')


class SeedTests(unittest.TestCase):
    def test_reverse_selected_contribution_then_match_original_rms_once(self):
        A = {9:2e-9, 5:-7e-10, 3:3., 7:4.}
        # Contribution includes the recorded sign. The raw field is (-1, 6).
        p2 = {9:1., 5:2., 3:2., 7:-1.}
        B, labels, k, start, raw = invert(A,p2,{5,9})
        self.assertEqual(labels,[3,7]); self.assertEqual(k,5/math.sqrt(37))
        self.assertAlmostEqual(B[3],-5/math.sqrt(37)); self.assertAlmostEqual(B[7],30/math.sqrt(37))
        self.assertEqual(B[5],A[5]); self.assertEqual(B[9],A[9]); self.assertAlmostEqual(start,k*raw)
        self.assertNotAlmostEqual(B[3],A[3]-5/math.sqrt(37))

    def test_exact_edit_roundtrip_pair_and_edges(self):
        z={5:-7e-10,9:2e-9,3:math.nextafter(1e-9,2e-9),7:-2e-7}; edge={5,9}
        raw=edit(DECK,rows(DECK),z,edge)
        a,changes=verify_edit(DECK,raw,z,edge)
        self.assertEqual(changes,[3,7]);self.assertEqual(a[3][1][2].hex(),z[3].hex())
        for lab in edge:
            self.assertEqual(rows(DECK)[lab][0],a[lab][0]);self.assertIn(a[lab][0],raw)
        other=DECK.replace(b'171000.',b'0.')
        b,_=verify_edit(other,edit(other,rows(other),z,edge),z,edge)
        self.assertEqual({i:v[0] for i,v in a.items()},{i:v[0] for i,v in b.items()})

    def test_independent_verifier_rejects_center_axis_stress_and_solver_changes(self):
        z={5:-7e-10,9:2e-9,3:1e-8,7:-2e-7}; edge={5,9}
        raw=edit(DECK,rows(DECK),z,edge)
        for a,b in ((b'2e-9',b'0.'),(b'5, 2.',b'5, 3.'),(b'3, 1.,',b'3, 1.01,'),
                    (b'1.71e+07',b'1.7e+07'),(b'171000.',b'8550000.'),(b'name=Step-GlobalTension',b'name=changed')):
            with self.subTest(a=a), self.assertRaisesRegex(ValueError,'changed|mismatch'):
                verify_edit(DECK,raw.replace(a,b),z,edge)
        with self.assertRaisesRegex(ValueError,'geometric edge'):
            verify_edit(DECK,edit(DECK,rows(DECK),{**z,9:0.},{5}),z,edge)

    def test_nonfinite_zero_rms_duplicate_missing_and_include_rejected(self):
        for A,p in (({1:1.},{2:0.}),({1:math.nan},{1:1.}),({1:0.},{1:0.}),({1:2.},{1:1.})):
            with self.assertRaises(ValueError): invert(A,p,set())
        for deck in (DECK+b'*Include, input=other.inp\r\n',DECK.replace(b'5, 2.',b'9, 2.')):
            with self.assertRaises(ValueError):rows(deck)
        with self.assertRaisesRegex(ValueError,'ordering'):
            verify_edit(DECK,DECK.replace(b'5, 2., 0., -7e-10\r\n',b''),{9:2e-9,5:-7e-10,3:0.,7:0.},{5,9})

    def test_independent_midside_and_reflection_require_complete_labels(self):
        m={4:[1,2],5:[2,3],6:[1,3]}; conn={11:('STRI65',(1,2,3,4,5,6))}
        verify_coverage([dict.fromkeys(range(1,7),1.)],{1,2,3},m,conn)
        with self.assertRaisesRegex(ValueError,'incomplete'):
            verify_coverage([dict.fromkeys(range(1,6),1.)],{1,2,3},m,conn)
        with self.assertRaisesRegex(ValueError,'map'):
            verify_coverage([dict.fromkeys(range(1,7),1.)],{1,2,3},{**m,6:[2,3]},conn)
        planar=DECK.replace(b'0., 2e-9',b'0., 0.').replace(b'0., -7e-10',b'0., 0.')
        ref=rows(planar);self.assertEqual(reflection(ref,{5,9}),{9:5,5:9,3:3,7:7})
        with self.assertRaisesRegex(ValueError,'reflection'):reflection(ref,{5})

    def test_wrong_source_digest_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'deck';p.write_bytes(DECK)
            ref={'A':{'path':str(p),'sha256':sha(DECK),'bytes':len(DECK)}}
            self.assertEqual(pin_all(ref)['A'],p)
            p.write_bytes(DECK.replace(b'171000.',b'900000.'))
            with self.assertRaisesRegex(ValueError,'input identity'):pin_all(ref)
