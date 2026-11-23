"""Fabricated inputs; never call a vendor executable or read project evidence."""
from decimal import Decimal as D
import math
import unittest

from relocate_clamp_position import endpoint, layout, points, relocate, select, verify


REF = {i:tuple(map(D,xyz)) for i,xyz in enumerate([
    (0,0,0),(20,0,0),(10,10,0),(3,3,0),(17,3,0),(4,4,0),(16,4,0),
    ('3.1',3,0),('16.9',3,0),('4.1',4,0),('15.9',4,0),
],1)}


def deck():
    rows=['*Heading','*Part, name=MEMBRANE','*Node','1, 0., 0., 1e-15','*End Part','*Assembly, name=A',
          '*Instance, name=MEMBRANE-1, part=MEMBRANE','*End Instance',
          '*Node','1, 10., 10., 0.']
    for side,x,rp,angle in [('CL',3,4,135),('CR',17,5,45)]:
        rows += [f'*Instance, name=Inst_cable_{side}, part=Cable_{side}',f'{x}., 3., 0.',
                 f'{x}., 3., 0., {x}., 3., 1., {angle}.','*End Instance',
                 '*Node',f'{rp}, {x}., 3., 0.']
        for name in [f'{side}_RP_Set',f'RP_{side}_Set']:
            rows.extend([f'*Nset, nset={name}',f'{rp},'])
        for name,node in [('START',1),('END',2)]:
            rows.extend([f'*Nset, nset=CABLE_{side}_{name}SET, instance=Inst_cable_{side}',f'{node},'])
    for set_,members in [('All_Edges_NoClamp','1, 2, 3, 6, 7'),('DIAG_EDGES','3, 6, 7'),('BASE_EDGE','1, 2'),
                         ('CL_Nodes','4, 8'),('CR_Nodes','5, 9'),('Left_Nodes','1'),('Right_Nodes','2'),('Top_Nodes','3')]:
        rows.extend([f'*Nset, nset={set_}, instance=MEMBRANE-1',members])
    for side in ('CL','CR'):
        for name in [f'{side}_Nodes',f'RP_{side}_Set',f'CABLE_{side}_STARTSET']:
            rows.extend([f'*Surface, type=NODE, name={name}_CNS_, internal',name+', 1.'])
        rows.extend([f'*Coupling, constraint name={side}_Coupling, ref node={side}_RP_Set, surface={side}_Nodes_CNS_',
                     '*Distributing, weighting method=UNIFORM',f'*Tie, name=Tie_{side}, adjust=yes',
                     f'CABLE_{side}_STARTSET_CNS_, RP_{side}_Set_CNS_'])
    rows+=['*End Assembly','*Boundary','All_Edges_NoClamp, 3, 3']
    for side in ('CL','CR'):
        rows += [f'CABLE_{side}_ENDSET, {i}, {i}' for i in (1,2,3)]
    rows+=['*Initial Conditions, type=STRESS','Inst_cable_CL.All_Elem, 171000., 0., 0.,',
           '*Step, name=Step-GlobalTension','*Static, stabilize=0.0002, allsdtol=0','0.001, 1., 1e-08, 0.005','*End Step','']
    return '\r\n'.join(rows).encode('ascii')


class RelocationTests(unittest.TestCase):
    def setUp(self):
        self.src=deck(); self.old=select(REF,D('.7')); self.new=select(REF,D('.6'))
        self.dst=relocate(self.src,self.old,self.new)

    def test_identity_and_moved_roles(self):
        self.assertEqual(relocate(self.src,self.old,self.old),self.src)
        links=verify(self.src,self.dst,self.new)
        self.assertEqual(links['CL']['RP_and_start_m'],[4.,4.,0.]); self.assertEqual(links['CR']['RP_and_start_m'],[16.,4.,0.])
        p=layout(self.dst)
        self.assertEqual(p['sets']['CL_Nodes'][1],{6,10}); self.assertEqual(p['sets']['CR_Nodes'][1],{7,11})
        self.assertEqual(p['sets']['All_Edges_NoClamp'][1],{1,2,3,4,5})
        self.assertEqual(p['sets']['DIAG_EDGES'][1],{3,4,5})
        for k,sign in [('CL',-1),('CR',1)]:
            x,y,z=links[k]['end_m']; a,b,_=links[k]['RP_and_start_m']
            self.assertAlmostEqual(x-a,sign*math.sqrt(.125)); self.assertAlmostEqual(y-b,math.sqrt(.125)); self.assertEqual(z,0.)
        self.assertNotIn(b'\n',self.dst.replace(b'\r\n',b''))

    def test_wrong_exact_center_is_not_snapped(self):
        q={**REF,6:(D('4.00001'),D('4.00001'),D(0))}
        with self.assertRaisesRegex(ValueError,'exact radius-zero center'):
            select(q,D('.6'))

    def test_bad_or_stale_membership_blocks(self):
        cl=self.dst.split(b'*Nset, nset=CL_Nodes, instance=MEMBRANE-1\r\n')[1].split(b'\r\n')[0]
        start=b'*Nset, nset=CL_Nodes, instance=MEMBRANE-1\r\n'
        for row in (b'4, 8',b'6',b'6, 6',b'7, 10'):
            with self.subTest(row=row), self.assertRaisesRegex(ValueError,'stale|duplicate'):
                verify(self.src,self.dst.replace(start+cl,start+row),self.new)
        with self.assertRaisesRegex(ValueError,'alias'):
            verify(self.src,self.dst.replace(b'*Nset, nset=RP_CL_Set\r\n4,',b'*Nset, nset=RP_CL_Set\r\n5,'),self.new)

    def test_protect_membrane_corner_stress_and_controls(self):
        for a,b in [(b'1e-15',b'0.'),(b'10., 10., 0.',b'11., 10., 0.'),(b'171000.',b'0.'),(b'allsdtol=0',b'allsdtol=0.05')]:
            with self.subTest(original=a), self.assertRaisesRegex(ValueError,'unapproved bytes'):
                verify(self.src,self.dst.replace(a,b),self.new)

    def test_nonidentity_rotation_moves_both_global_axis_points(self):
        # Noncoincident translation and rotation centre: translation before rotation.
        self.assertEqual(endpoint(tuple(map(D,(3,2,0))),tuple(map(D,(1,2,0,1,2,1,90))),D(0)),[1.0000000000000002,4.,0.])
        bad=self.dst.replace(b'4., 4., 0., 4., 4., 1., 135.',b'4., 4., 0., 3., 3., 1., 135.')
        with self.assertRaisesRegex(ValueError,'axis'):
            verify(self.src,bad,self.new)


if __name__=='__main__':
    unittest.main()
