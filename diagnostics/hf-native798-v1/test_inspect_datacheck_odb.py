import json
from pathlib import Path
import tempfile
from types import SimpleNamespace as NS
import unittest

import inspect_datacheck_odb as inspection


class Repository:
    """Expose only indexing and keys, as used for native repositories."""
    def __init__(self, **items):
        self.items = items

    def keys(self):
        return list(self.items)

    def __getitem__(self, key):
        return self.items[key]


def node(n):
    return NS(label=n)


def region(nodes, names=None, internal=True):
    return NS(nodes=nodes, instanceNames=names, isInternal=internal)


def instance(**sets):
    return NS(nodeSets=Repository(**sets), surfaces=Repository(), elementSets=Repository(),
              elements=[NS(type='T3D2', connectivity=(1,2))])


def input_excerpt():
    result = []
    for n,k in enumerate(inspection.CABLES,1):
        pre = f'CABLE_{k}'
        for name, value, scope in ((pre+'_STARTSET',1,f', instance=I_{k}'),
                                   (pre+'_ENDSET',2,f', instance=I_{k}'), (f'RP_{k}_SET',n,'')):
            result.extend([f'*Nset, nset={name}{scope}', f'{value},'])
            result.extend([f'*Surface, type=NODE, name={name}_CNS_, internal', f'{name}, 1.'])
        result.extend([f'*Tie, name=Tie_{k}, adjust=yes',f'{pre}_STARTSET_CNS_, RP_{k}_SET_CNS_'])
    return '\n'.join(result)+'\n'


class NativeSelectionTests(unittest.TestCase):
    def test_flat_and_grouped_arrays_do_not_merge_same_labels(self):
        grouped = region(([node(2)], [node(2)], [node(2)]), ('I_LEFT','I_TOP',None))
        self.assertEqual(inspection.set_nodes(grouped), {('I_LEFT',2),('I_TOP',2),('ASSEMBLY',2)})
        self.assertEqual(inspection.set_nodes(region([node(8)], ('I_TOP',))), {('I_TOP',8)})
        self.assertEqual(inspection.set_nodes(region([node(8)]), 'I_LEFT'), {('I_LEFT',8)})
        self.assertEqual(inspection.set_nodes(region([])), set())
        with self.assertRaisesRegex(ValueError,'namespaces'):
            inspection.set_nodes(region([node(8)]))

    def test_inspection_distinguishes_warning_mismatch_and_tie_coverage(self):
        expected_scopes = ('I_TOP','I_LEFT','I_RIGHT')
        sets = {'warn4':region(([node(2)],)*3, expected_scopes),
                'WARN5':region(([node(2)],)*3, ('I_TOP','I_LEFT',None)),
                'INTERNAL_TIED_SECONDARY':region(([node(1)],)*2, ('I_TOP','I_RIGHT')),
                'INTERNAL_UNTIED_SECONDARY':region(([node(1)],)*2, ('I_LEFT','I_RIGHT')),
                'USER_TIED_SECONDARY':region([node(1)], ('I_CL',), internal=False),
                'WarnNodeTieCandidate':region([node(1)], ('I_CR',))}
        assembly = NS(nodeSets=Repository(**sets), elementSets=Repository(),
                      surfaces=Repository(RP_TOP_SET_CNS_=region([node(1)], (None,))), elements=[],
                      instances=Repository(**{f'I_{n}':instance() for n in inspection.CABLES}))
        frames = [NS(incrementNumber=0,frameValue=0.,fieldOutputs=Repository(U=object()))]
        odb = NS(rootAssembly=assembly,steps=Repository(Static=NS(frames=frames)))
        dat = '\n'.join(f'Boundary conditions are specified on inactive dof {i} at 3 nodes.\n'
                        f'The nodes have been identified in node set WARN{i}.' for i in range(4,7))
        dat += ('\nFOR *TIE PAIR (ASSEMBLY_CABLE_TOP_STARTSET_CNS_-ASSEMBLY_RP_TOP_SET_CNS_),\n'
                'ADJUSTMENT WAS SPECIFIED BUT NO NODE WAS ADJUSTED MORE THAN THE\nADJUSTMENT DISTANCE = 2.22000E-16.')
        report = inspection.inspect(odb,input_excerpt(),dat)
        warns = {r['dof']:r for r in report['inactive_dofs']}
        self.assertEqual(warns[4]['status'],'match')
        self.assertEqual(warns[4]['nodes'], [dict(node=(n,2),incident_element_types=['T3D2']) for n in sorted(expected_scopes)])
        self.assertEqual((warns[5]['status'],warns[5]['unexpected']),('mismatch',[('ASSEMBLY',2)]))
        self.assertEqual(warns[6]['status'],'not_exposed_or_unresolved')
        ties = {r['tie']:r for r in report['ties']}
        self.assertEqual({k:v['coverage'] for k,v in ties.items()},
                         dict(TIE_TOP='confirmed_tied',TIE_LEFT='confirmed_untied',TIE_RIGHT='not_exposed_or_unresolved',
                              TIE_CL='not_exposed_or_unresolved',TIE_CR='not_exposed_or_unresolved'))
        self.assertTrue(ties['TIE_TOP']['odb_surfaces']['main']['match'])
        self.assertEqual(ties['TIE_TOP']['adjustment_text']['threshold'],'2.22000E-16')
        self.assertIsNone(ties['TIE_TOP']['effective_position_tolerance'])
        self.assertEqual(report['steps']['Static'][0]['keys'],['U'])
        self.assertEqual(report['adjustment_fields'],[])

    def test_wrong_input_tie_roles_are_not_accepted(self):
        with self.assertRaisesRegex(ValueError,'wrong secondary/main roles'):
            inspection.deck_targets(input_excerpt().replace('CABLE_TOP_STARTSET_CNS_, RP_TOP_SET_CNS_',
                                                             'RP_TOP_SET_CNS_, CABLE_TOP_STARTSET_CNS_',1))

    def test_attribution_fault_does_not_select_another_same_named_case(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            records = root/'run-records'
            records.mkdir()
            want = 'a'*64
            args = 'job=DC input=approved.inp datacheck cpus=4 interactive'
            base = dict(cwd=str(root),args=args)
            control = dict(**base,job='DC',attempts=1,closed=True,timeout=False,remaining=[],
                           input_sha256_before=want,input_sha256_after=want)
            for name, data in (('control-outcome.json',control),('datacheck.launch.json',base),
                               ('datacheck.completion.json',dict(status='launcher_exited'))):
                (records/name).write_text(json.dumps(data), encoding='utf-8')
            (root/'DC.dat').write_text('ANALYSIS DATACHECK COMPLETE\n',encoding='utf-8')
            self.assertEqual(inspection.attribution(root/'DC.odb',root/'approved.inp',want)[0],'DC')
            with self.assertRaisesRegex(ValueError,'attribution'):
                inspection.attribution(root/'other.odb',root/'approved.inp',want)
            with self.assertRaisesRegex(ValueError,'attribution'):
                inspection.attribution(root/'DC.odb',root/'approved.inp','0'*64)


if __name__ == '__main__':
    unittest.main()
