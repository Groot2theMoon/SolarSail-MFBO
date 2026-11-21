"""Portable regression for the guarded CAE attachment audit."""
import ast
from pathlib import Path
import sys
from types import SimpleNamespace as NS
import unittest
from unittest.mock import patch

CODE = Path(__file__).resolve().parent / 'package' / 'staged-code'
sys.path.insert(0, str(CODE))
import build_only_guard as guard


class Feature:
    def __init__(self, feature_id):
        self.id = feature_id


class ReferencePoint:
    __slots__ = ('coordinates',)

    def __init__(self, coordinates):
        self.coordinates = coordinates


class Assembly:
    def __init__(self):
        self.sets = {}

    def getCoordinates(self, *, entity):
        if not isinstance(entity, ReferencePoint):
            raise TypeError('Expected a ReferencePoint, not a Feature')
        return entity.coordinates


def rp_set(*points):
    return NS(referencePoints=points)


class AttachmentRPTests(unittest.TestCase):
    def test_driver_records_creation_provenance_for_all_five_attachments(self):
        driver = CODE / 'run_abaqus.py'
        tree = ast.parse(driver.read_text(encoding='utf-8'))
        (body,) = (n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)
                   and n.name == '_record_final_assembly')
        assembly = Assembly()
        scope = dict(a=assembly, _build_only=guard, aba_grid_mesh=None,
                     _membership_before={}, _state2_seeded={}, INSTANCE_NAME='MEMBRANE-1',
                     HF_JOB=guard.JOB, CABLE_PRESTRESS=1, CLAMP_PRESTRESS=2)
        cases = [
            ('Top', 'Cable_Top', 'V1', 'rp1_obj', 'top', 101),
            ('Right', 'Cable_Right', 'V2', 'rp2_obj', 'right', 203),
            ('Left', 'Cable_Left', 'V3', 'rp3_obj', 'left', 307),
            ('CL', 'cable_CL', 'V_CL', 'rp_cl_obj', 'clamp_left', 409),
            ('CR', 'cable_CR', 'V_CR', 'rp_cr_obj', 'clamp_right', 503),
        ]
        wanted = {}
        for i, (name, cable, point_name, feature_name, key, feature_id) in enumerate(cases):
            point = (float(i), float(2*i), 0.)
            scope[point_name], scope[feature_name] = point, Feature(feature_id)
            assembly.sets[name + '_RP_Set'] = rp_set(ReferencePoint(point))
            assembly.sets['RP_' + name + '_Set'] = rp_set(ReferencePoint(point))
            assembly.sets[name + '_Nodes'] = NS()
            for label, position in enumerate(('START', 'END'), 1):
                node = NS(instanceName='Inst_' + cable, label=label, coordinates=point)
                assembly.sets[cable.upper() + '_' + position + 'SET'] = NS(nodes=(node,))
            wanted[key] = dict(instance=None, label=None, feature_id=feature_id,
                               xyz_m=point, alias_identity='verify_written_node')
        # Only stub the earlier membrane checks; exercise the actual driver RP/cable audit.
        with patch.object(guard, 'final_sets', return_value={}), \
             patch.object(guard, 'require_xy_mirror'), patch.object(guard, 'record') as record:
            exec(compile(ast.Module(body=[body], type_ignores=[]), str(driver), 'exec'), scope)
            scope['_record_final_assembly']()
        (attachments,) = [call.args[1] for call in record.call_args_list
                          if call.args[0] == 'attachments.json']
        self.assertEqual({row['key']: row['rp'] for row in attachments['rows']}, wanted)

    def test_rejects_missing_empty_multiple_or_misplaced_in_either_set(self):
        expected = (1., 2., 0.)
        for name in ('control', 'alias'):
            for bad, reason in [(None, 'Missing RP set'), ((), 'Not a singleton RP set'),
                                ((expected, expected), 'Not a singleton RP set'),
                                (((1., 2., 1.),), 'Unexpected RP placement')]:
                with self.subTest(name=name, reason=reason):
                    a = Assembly()
                    a.sets = {key: rp_set(ReferencePoint(expected)) for key in ('control', 'alias')}
                    if bad is None:
                        del a.sets[name]
                    else:
                        a.sets[name] = rp_set(*(ReferencePoint(p) for p in bad))
                    with self.assertRaisesRegex(guard.GuardViolation, reason + ': ' + name):
                        guard.attachment_rp(a, 'control', 'alias', Feature(101), expected)


if __name__ == '__main__':
    unittest.main()
