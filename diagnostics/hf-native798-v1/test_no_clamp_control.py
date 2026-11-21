import unittest

from make_no_clamp_control import transform, verify_intervention


def fixture():
    source, control, drops, masks = [], [], [], []
    def keep(text):
        rows = text.splitlines(keepends=True)
        source.extend(rows)
        control.extend(rows)
    def drop(text):
        at = len(source) + 1
        source.extend(text.splitlines(keepends=True))
        drops.append((at,len(source)))
    def make_set(name, old, new):
        keep(f'*Nset, nset={name}, instance=MEMBRANE-1\n')
        at = len(source)+1
        source.append(', '.join(map(str,old))+'\n')
        control.append(', '.join(map(str,new))+'\n')
        masks.append((at,at,new))

    keep('*Heading\n*Preprint, model=YES\n*Part, name=Cable_Top\n*Node\n1, 0, 0, 0\n*End Part\n')
    for x in ('CL','CR'):
        drop(f'*Part, name=Cable_{x}\n*Node\n1, 0, 0, 0\n*Element, type=T3D2\n1, 1, 2\n*End Part\n')
    keep('*Part, name=MEMBRANE\n*Node\n1, 0, 0, 0\n2, 1, 0, 0\n3, 3, 3, 1.9E-10\n4, 17, 3, 8.1E-11\n'
         '*Element, type=STRI65\n3, 1, 2, 3, 5, 6, 7\n*End Part\n*Assembly, name=Assembly\n'
         '*Instance, name=MEMBRANE-1, part=MEMBRANE\n*End Instance\n*Instance, name=Inst_Cable_Top, part=Cable_Top\n*End Instance\n')
    for x in ('CL','CR'):
        drop(f'*Instance, name=Inst_cable_{x}, part=Cable_{x}\n0, 0, 0\n*End Instance\n')
    keep('*Node\n1, 10, 10, 0\n')
    for x, rp in (('CL',4),('CR',5)):
        drop(f'*Node\n{rp}, 3, 3, 0\n')
    keep('*Nset, nset=BASE_EDGE, instance=MEMBRANE-1\n1\n')
    make_set('DIAG_EDGES',[2],[2,3,4])
    make_set('All_Edges_NoClamp',[1,2],[1,2,3,4])
    for x,rp,node in (('CL',4,3),('CR',5,4)):
        drop(f'*Nset, nset=CABLE_{x}_STARTSET, instance=Inst_cable_{x}\n1,\n')
        drop(f'*Nset, nset=CABLE_{x}_ENDSET, instance=Inst_cable_{x}\n2,\n')
        drop(f'*Nset, nset={x}_Nodes, instance=MEMBRANE-1\n{node}\n')
        drop(f'*Nset, nset=RP_{x}_Set\n{rp}\n')
        drop(f'*Nset, nset={x}_RP_Set\n{rp}\n')
        for name in (f'{x}_Nodes',f'RP_{x}_Set',f'CABLE_{x}_STARTSET'):
            drop(f'*Surface, type=NODE, name={name}_CNS_, internal\n{name}, 1.\n')
        drop(f'** Constraint: {x}_Coupling\n*Coupling, constraint name={x}_Coupling, ref node={x}_RP_Set, surface={x}_Nodes_CNS_\n*Distributing, weighting method=UNIFORM\n')
        drop(f'** Constraint: Tie_{x}\n*Tie, name=Tie_{x}, adjust=yes\nCABLE_{x}_STARTSET_CNS_, RP_{x}_Set_CNS_\n')
    keep('*End Assembly\n*Material, name=Kevlar\n*Elastic\n6e10, 0.36\n*Boundary\nAll_Edges_NoClamp, 3, 3\n')
    for x in ('CL','CR'):
        drop(f'** Name: BC_{x}\n*Boundary\nCABLE_{x}_ENDSET, 3, 3\n')
        drop(f'** Name: Disp_Control_{x}\n*Boundary\nCABLE_{x}_ENDSET, 1, 1\nCABLE_{x}_ENDSET, 2, 2\n')
        drop(f'** Name: Prestress_{x}\n*Initial Conditions, type=STRESS\nInst_cable_{x}.All_Elem, 171000., 0., 0.,\n')
    keep('*Initial Conditions, type=STRESS\nInst_Cable_Top.All_Elem, 1.71e7, 0, 0\n'
         '*Step, name=Step-GlobalTension\n*Static, stabilize=0.0002, allsdtol=0\n0.001, 1., 1e-8, 0.005\n*End Step\n')
    return (''.join(source).encode(), ''.join(control).encode(), drops, masks,
            {'unseeded':{'edges':[1,2],'excluded_edges':[3,4],
                         'clouds':{'clamp_left':[3],'clamp_right':[4]}},'intended_edge':[1,2,3,4]})


class TestControl(unittest.TestCase):
    def test_exact_intervention_matches_independently_authored_deck(self):
        source, control, drops, masks, membership = fixture()
        actual, receipt = transform(source,membership)
        self.assertEqual(actual,control)
        self.assertEqual(receipt['restored_u3'],[3,4])
        verify_intervention(source,actual,drops,masks)

    def test_verifier_rejects_unapproved_mutations(self):
        a,b,d,m,_ = fixture()
        for old,new in ((b'*Instance, name=Inst_Cable_Top, part=Cable_Top\n*End Instance\n',b''),
                        (b'*End Assembly',b'*Instance, name=Inst_cable_CL, part=Cable_CL\n*End Instance\n*End Assembly'),
                        (b'1.9E-10',b'0.0'), (b'3, 3, 3, 1.9E-10\n4, 17, 3, 8.1E-11',b'4, 17, 3, 8.1E-11\n3, 3, 3, 1.9E-10'),
                        (b'1, 2, 3, 4\n',b'1, 2, 3, 5\n'), (b'6e10',b'7e10'), (b'1.71e7',b'0'), (b'0.0002',b'0.002'),
                        (b'3, 1, 2, 3, 5, 6, 7',b'3, 1, 2, 3, 6, 5, 7')):
            with self.subTest(old=old):
                with self.assertRaises(ValueError):
                    verify_intervention(a,b.replace(old,new,1),d,m)

    def test_unknown_dependency_and_retained_shared_alias_are_not_silently_deleted(self):
        a,_,_,_,f = fixture()
        for marker, reason in ((b'*Equation\n2\n', 'unsupported keyword'),
                               (b'*Nset, nset=Mixed\nCABLE_CL_STARTSET\n', 'unresolved dependent'),
                               (b'*Nset, nset=Mixed\n4, 1\n', 'shares the clamp RP')):
            with self.subTest(marker=marker):
                with self.assertRaisesRegex(ValueError,reason):
                    transform(a.replace(b'*End Assembly',marker+b'*End Assembly'),f)


if __name__ == '__main__':
    unittest.main()
