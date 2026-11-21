"""Portable optics regressions and fabricated, no-ODB export checks."""

from copy import deepcopy
import hashlib
import io
import json
import math
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch
import zipfile
import zlib

import numpy as np

import export_sail_geometry as exporter
import sail_optics as optics


_SCRIPT = Path(optics.__file__).resolve()
_NORMAL = (0., 0., -1.)
_OBLIQUE = (-3. / 5., 0., -4. / 5.)
_TRI = ((0, 0, 0), (2, 0, 0), (0, 1, 0))
_TILT = ((0, 0, 0), (1, 0, 0), (0, 1, -1))
_PANEL = ((0, 0, 0), (20, 0, 0), (10, 10, 0))
_CONTROL = ((0, 0, 0), (10, 0, 0), (5, 5, 0))
_SMALL = ((0, 0, 0), (8, 0, 0), (4, 4, 0))
_UNIT_PLACEMENT = ((np.eye(3), np.zeros(3)),)


def _json(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    return path


def _read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def _digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _ref(base, path):
    return {'path': path.relative_to(base).as_posix(), 'sha256': _digest(path)}


def _arrays(vertices, *, kind='S3'):
    # Nonconsecutive labels, and coordinate rows in the opposite order to connectivity.
    labels = 11 + 7 * np.arange(len(vertices), dtype=np.int64)
    return {'node_labels': labels[::-1], 'node_xyz_m': np.array(vertices[::-1], dtype=np.float64),
            f'element_labels_{kind}': np.array([101], dtype=np.int64),
            f'connectivity_{kind}': labels.reshape(1, -1),
            f'face_sign_{kind}': np.ones(1, dtype=np.int8)}


def _bundle(root, case_id, vertices=_PANEL, *, kind='S3', array_changes=None, **meta_changes):
    directory = root / case_id
    directory.mkdir()
    arrays = _arrays(vertices, kind=kind)
    for key, value in (array_changes or {}).items():
        if value is None:
            del arrays[key]
        else:
            arrays[key] = value
    npz = directory / 'geometry.npz'
    np.savez(npz, **arrays)
    path = _json(directory / 'geometry.json', {
        'schema_version': 1, 'source_kind': 'synthetic', 'case_id': case_id,
        'frame': 'panel_local_m', 'npz': _ref(directory, npz), 'node_count': len(vertices),
        'element_counts': {kind: 1}, 'export': None, **meta_changes})
    return {'role': 'clamped', 'geometry': _ref(root, path)}


def _initial_case(role, path):
    return {'role': role, 'geometry': None, 'odb_path': path, 'attempt_id': None,
            'records': {}, 'run': None, 'producer': None, 'geometry_basis': None,
            'settings': dict.fromkeys(('membrane', 'load', 'solver', 'supports', 'clamp', 'perturbation')),
            'diagnostics': {'status': 'unresolved', 'observations': [], 'reason': 'fabricated; no job'}}


def _study(cases, contrasts=None, **fields):
    return {'schema_version': 1, 'source_kind': 'synthetic', 'layout_id': 'nominal_square_v1',
            'optics': 'default', 'cases': cases, 'contrasts': contrasts or {}, **fields}


def write_verification_inputs(directory):
    """Write only labelled synthetic bundles and explicitly nonexistent ODB paths."""
    root = Path(directory)
    synthetic = root / 'mixed-synthetic'
    synthetic.mkdir(parents=True)
    cases = {'control': _bundle(synthetic, 'control'),
             'usable': _bundle(synthetic, 'usable', _CONTROL),
             'bad': _bundle(synthetic, 'bad', element_counts={'S3': 1, 'S4': 1})}
    cases['control']['role'] = 'control'
    pairs = {name: {'case_id': name, 'control_id': 'control'} for name in ('usable', 'bad')}
    odb = root / 'missing-odb'
    odb.mkdir()
    candidates = {key: _initial_case(case['role'], f'never-created/{key}.odb')
                  for key, case in cases.items()}
    return (_json(synthetic / 'study.json', _study(cases, pairs, optics='reflecting_only', pressure_pa=0)),
            _json(odb / 'study.initial.json', _study(candidates, pairs, source_kind='odb')))


def write_malformed_input(directory, *, corrupt_deflate):
    """Create a mixed study whose malformed member passes all file-digest checks."""
    root = Path(directory)
    root.mkdir(parents=True)
    cases = {'usable': _bundle(root, 'usable'), 'bad': _bundle(root, 'bad')}
    bad_dir = root / 'bad'
    npz = bad_dir / 'geometry.npz'
    member = 'node_labels.npy'
    stream = io.BytesIO()
    with zipfile.ZipFile(npz) as original, zipfile.ZipFile(stream, 'w') as rewritten:
        for name in original.namelist():
            payload, compression = original.read(name), zipfile.ZIP_STORED
            if name == member:
                if corrupt_deflate:
                    compression = zipfile.ZIP_DEFLATED
                else:
                    payload = b'not an NPY array'
            rewritten.writestr(name, payload, compress_type=compression)
        header = rewritten.getinfo(member).header_offset
    data = bytearray(stream.getvalue())
    if corrupt_deflate:
        name_len, extra_len = struct.unpack_from('<HH', data, header + 26)
        # BTYPE=3 is reserved; fail DEFLATE itself before the ZIP CRC check.
        data[header + 30 + name_len + extra_len] = 0x07
    npz.write_bytes(data)
    meta = bad_dir / 'geometry.json'
    _json(meta, {**_read(meta), 'npz': _ref(bad_dir, npz)})
    cases['bad']['geometry'] = _ref(root, meta)
    return _json(root / 'study.json', _study(cases, optics='reflecting_only'))


def _cli(command, manifest, out):
    return subprocess.run([sys.executable, str(_SCRIPT), command, '--manifest', str(manifest),
                           '--out-dir', str(out)], capture_output=True, text=True, check=False)


def _facet(vertices, sign=1):
    return optics.Facets(np.array([vertices], dtype=np.float64), np.array([sign], dtype=np.int8))


class AnalyticalTests(unittest.TestCase):
    def assert_load(self, result, force, moment, *, area=1., length=math.sqrt(5.)):
        for actual, expected, floor in zip(result, (force, moment), (1e-12 * area, 1e-12 * area * length)):
            expected = np.array(expected)
            bound = floor + 1e-10 * np.maximum(np.abs(actual), np.abs(expected))
            self.assertTrue(np.isfinite(actual).all() and np.all(np.abs(actual - expected) <= bound),
                            (actual, expected, bound))

    def test_flat_absorption_and_reflection_at_normal_and_oblique(self):
        for constants, normal, oblique in (
            ((1., 0.), (0, 0, -1), (-12/25, 0, -16/25)),
            ((0., 1.), (0, 0, -2), (0, 0, -32/25)),
        ):
            for ray, force in ((_NORMAL, normal), (_OBLIQUE, oblique)):
                for direct in (False, True):
                    with self.subTest(optics=constants, ray=ray, direct=direct):
                        result = optics._route(_facet(_TRI), np.array(ray), constants, _UNIT_PLACEMENT,
                                               direct=direct, origin=np.array([2/3, 1/3, 0.]))[0]
                        self.assert_load(result, force, (0, 0, 0))

    def test_flat_square_force_and_cancellation(self):
        for direct in (False, True):
            result = optics._route(_facet(_PANEL), np.array(_NORMAL), (0., 1.),
                                   optics._placements(), direct=direct).sum(axis=0)
            self.assert_load(result, (0, 0, -800), (0, 0, 0), area=400, length=20)

    def test_moment_reference_with_signed_components(self):
        for direct in (False, True):
            f, ma = optics._route(_facet(_TRI), np.array(_NORMAL), (1., 0.), _UNIT_PLACEMENT,
                                  direct=direct)[0]
            fp, mb = optics._route(_facet(_TRI), np.array(_NORMAL), (1., 0.), _UNIT_PLACEMENT,
                                   direct=direct, origin=np.array([1., 0, 0]))[0]
            self.assert_load((f, ma), (0, 0, -1), (-1/3, 2/3, 0))
            self.assert_load((fp, mb), (0, 0, -1), (-1/3, -1/3, 0))
            np.testing.assert_allclose(mb - ma, (0, -1, 0), rtol=0, atol=1e-12)

    def test_ordered_base_endpoints_and_apex_in_all_placements(self):
        expected = [((-10, -10, 0), (10, -10, 0), (0, 0, 0)),
                    ((10, -10, 0), (10, 10, 0), (0, 0, 0)),
                    ((10, 10, 0), (-10, 10, 0), (0, 0, 0)),
                    ((-10, 10, 0), (-10, -10, 0), (0, 0, 0))]
        for (r, t), vertices in zip(optics._placements(), expected, strict=True):
            np.testing.assert_array_equal(np.asarray(_PANEL) @ r.T + t, vertices)

    def test_tilt_winding_cancellation_and_retained_axial_moment(self):
        placements = tuple((r, np.zeros(3)) for r, _ in optics._placements())
        reverse = (_TILT[0], _TILT[2], _TILT[1])
        for direct in (False, True):
            for shape in (_facet(_TILT), _facet(reverse, -1)):
                result = optics._route(shape, np.array(_NORMAL), (0., 1.), placements, direct=direct)
                self.assert_load(result[0], (0, -1/2, -1/2), (-1/3, 1/6, -1/6),
                                 area=0.5, length=math.sqrt(2.))
                self.assert_load(result.sum(axis=0), (0, 0, -2), (0, 0, -2/3),
                                 area=2, length=math.sqrt(2.))

    def test_s8r5_warp_and_stri65_use_exact_ordered_facets(self):
        square = ((0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0), (1, 0, 1),
                  (2, 1, 0), (1, 2, 0), (0, 1, 0))
        centre = (1, 1, 0)
        a, b, c, d, ab, bc, cd, da = square
        expected = [(a, ab, centre), (a, centre, da), (ab, b, bc), (ab, bc, centre),
                    (centre, bc, c), (centre, c, cd), (da, centre, cd), (da, cd, d)]
        np.testing.assert_array_equal(optics._split(np.array([square], dtype=float), 'S8R5'), expected)
        a, b, c = _TRI
        ab, bc, ca = (1, 0, .1), (1, .5, .2), (0, .5, .3)
        pts = [(a, b, c, ab, bc, ca)]
        expected = [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]
        np.testing.assert_array_equal(optics._split(np.array(pts), 'STRI65'), expected)
        np.testing.assert_array_equal(optics._split(np.array([square[:4]]), 'S4'),
                                      [(square[0], square[1], square[2]), (square[0], square[2], square[3])])

    def test_oblique_rotated_translated_panel(self):
        ray = np.array([1/3, 2/3, -2/3])
        for direct in (False, True):
            result = optics._route(_facet(_TILT), ray, (0., 1.), optics._placements(), direct=direct)
            self.assert_load(result[1], (1/2, 0, -1/2), (29/6, 14/3, 29/6), area=.5, length=20)
            block = optics._blocks(*result[1], ray)
            self.assertAlmostEqual(block['ray_projection'], 1/2, delta=1e-12)
            np.testing.assert_allclose(block['transverse_force'], (1/3, -1/3, -1/6), rtol=0, atol=1e-12)

    def test_back_face_sign_only_is_exactly_dark(self):
        for direct in (False, True):
            result = optics._route(_facet(_TRI, -1), np.array(_NORMAL), (0.073744, 0.926256),
                                   _UNIT_PLACEMENT, direct=direct)[0]
            np.testing.assert_array_equal(result, np.zeros((2, 3)))

    def test_nonunit_and_nonfinite_rays_are_rejected(self):
        for ray in ((0, 0, 0), (1, 1, 1), (np.nan, 0, 0), (0, 0, np.inf)):
            with self.subTest(ray=ray), self.assertRaises(optics.InputError):
                optics._unit_ray(ray)
        np.testing.assert_array_equal(optics._unit_ray((0, 0, 1 + 5e-13)), (0, 0, 1))


class NumericalPublicationTests(unittest.TestCase):
    def _load(self, pressure):
        return optics._load(_facet(_PANEL), np.array(_NORMAL), (0., 1.), optics._placements(), pressure)

    def _row(self, force, pressure=1.):
        normalized = optics._blocks(np.array(force), np.zeros(3), np.array(_NORMAL))
        physical = {key: (np.asarray(value) * pressure).tolist() for key, value in normalized.items()}
        return {'validation': 'passed', 'normalized': normalized, 'physical': physical}

    def test_both_direction_operands_must_strictly_clear_floor(self):
        floor = 4e-10
        for case in (floor / 2, floor, floor * 2):
            for ref in (floor / 2, floor, floor * 2):
                with self.subTest(case=case, ref=ref):
                    row = optics._comparison(self._row([case, 0, 0]), self._row([0, ref, 0]),
                                             'flat', 1., flat=True)
                    if case > floor and ref > floor:
                        self.assertAlmostEqual(row['angle_normalized_deg'], 90.)
                        self.assertAlmostEqual(row['angle_physical_deg'], 90.)
                    else:
                        for key in ('angle_normalized_deg', 'angle_physical_deg'):
                            self.assertIsNone(row[key])
                            self.assertIn({'code': 'force_at_or_below_floor', 'target': key}, row['reasons'])

    def test_pressure_scaling_zero_omission_and_overflow(self):
        omitted = self._load(None)
        positive = self._load(2.)
        zero = self._load(0.)
        huge = self._load(1e307)
        self.assertIsNone(omitted['physical'])
        self.assertEqual(positive['physical']['force'], [0., 0., -1600.])
        self.assertEqual(positive['physical']['ray_projection'], 1600.)
        self.assertEqual(zero['physical']['force'], [0., 0., 0.])
        self.assertIsNone(huge['physical'])
        self.assertTrue(any(item['code'] == 'physical_overflow' for item in huge['reasons']))
        for load in (positive, zero, huge):
            self.assertEqual(load['normalized'], omitted['normalized'])
            self.assertEqual(load['validation'], 'passed')
        for load, press, code in ((omitted, None, 'pressure_omitted'), (zero, 0., 'zero_pressure'),
                                  (huge, 1e307, 'physical_overflow')):
            comp = optics._comparison(load, load, 'flat', press, flat=True)
            self.assertEqual(comp['angle_normalized_deg'], 0.)
            self.assertIsNone(comp['angle_physical_deg'])
            self.assertIn({'code': code, 'target': 'angle_physical_deg'}, comp['reasons'])

    def test_signed_deltas_for_both_reference_types(self):
        case, ref = self._row([3., 0., 0.], 2.), self._row([0., 4., 0.], 2.)
        for ref_id in ('flat', 'control'):
            comp = optics._comparison(case, ref, ref_id, 2., flat=ref_id == 'flat')
            self.assertEqual(comp['normalized']['force'], [3., -4., 0.])
            self.assertEqual(comp['normalized']['force_norm'], -1.)
            self.assertEqual(comp['physical']['force'], [6., -8., 0.])
            self.assertEqual(comp['physical']['force_norm'], -2.)

    def test_comparison_overflow_does_not_invalidate_operands(self):
        case, ref = self._row([1e308, 0., 0.]), self._row([-1e308, 0., 0.])
        comp = optics._comparison(case, ref, 'flat', 1., flat=True)
        self.assertIsNone(comp['normalized'])
        self.assertIsNone(comp['physical'])
        self.assertIn({'code': 'nonfinite', 'target': 'normalized.force'}, comp['reasons'])
        self.assertEqual(comp['angle_normalized_deg'], 180.)
        self.assertEqual(case['normalized']['force'][0], 1e308)
        self.assertEqual(ref['normalized']['force'][0], -1e308)

    def test_nonfinite_route_or_derived_load_withholds_even_at_zero_pressure(self):
        vectors = np.zeros((4, 2, 3))
        for xyz in ((np.inf, 0, 0), (1.5e308, 1.5e308, 0)):
            with self.subTest(xyz=xyz):
                vectors[0, 0] = xyz
                with patch.object(optics, '_route', return_value=vectors):
                    row = self._load(0.)
                self.assertEqual(row['validation'], 'failed')
                self.assertIsNone(row['normalized'])
                self.assertIsNone(row['physical'])
                self.assertTrue(any(n['code'] == 'nonfinite' for n in row['reasons']))

    def test_route_mismatch_withholds_only_dependents_at_that_attitude(self):
        formulas = optics._route
        surfaces = {'control': _facet(_CONTROL), 'case': _facet(_SMALL)}
        study = _study(dict.fromkeys(surfaces), pressure_pa=0)
        evidence = {'study': {}, 'contrasts': {'pair': {'case_id': 'case', 'control_id': 'control',
                                                       'status': 'eligible'}}}
        for fault, width in (('flat', 20), ('control', 10), ('case', 8)):
            def one_wrong_route(facets, ray, *args, direct, **kwargs):
                values = formulas(facets, ray, *args, direct=direct, **kwargs)
                if direct and facets.xyz[0, 1, 0] == width and tuple(ray) == _NORMAL:
                    values[0, 0, 0] += 1.
                return values
            with self.subTest(fault=fault), patch.object(optics, '_route', side_effect=one_wrong_route):
                result = optics._evaluate(study, evidence, surfaces, {})
                failed = result['loads'][fault]['normal']
                self.assertIsNone(failed['normalized'])
                self.assertIsNone(failed['physical'])
                self.assertIn({'code': 'route_mismatch', 'target': 'normalized.force.0'}, failed['reasons'])
                for other in ({'flat', 'control', 'case'} - {fault}):
                    self.assertEqual(result['loads'][other]['normal']['validation'], 'passed')
                rows = result['comparisons']['case']
                self.assertEqual(rows['normal']['flat']['normalized'] is None, fault != 'control')
                self.assertEqual(rows['normal']['control']['normalized'] is None, fault != 'flat')
                for entry in rows['alpha5_phi45'].values():
                    self.assertIsNotNone(entry['normalized'])
                self.assertEqual(evidence['contrasts']['pair']['status'], 'eligible')


class CommandTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)

    def _run(self, command, path, name, expected=0):
        out = self.root / name
        run = _cli(command, path, out)
        self.assertEqual(run.returncode, expected, run.stderr)
        if expected:
            self.assertIn('sail_optics', run.stderr)
            self.assertFalse(out.exists())
        return out

    def test_all_element_types_produce_load_from_versioned_bundles(self):
        square = ((0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0))
        cases = {'tri': _bundle(self.root, 'tri', _TRI),
                 'quad': _bundle(self.root, 'quad', square, kind='S4'),
                 'tri6': _bundle(self.root, 'tri6', _TRI + ((1, 0, 0), (1, .5, 0), (0, .5, 0)),
                                 kind='STRI65'),
                 'quad8': _bundle(self.root, 'quad8', square + ((1, 0, 0), (2, 1, 0), (1, 2, 0), (0, 1, 0)),
                                  kind='S8R5')}
        path = _json(self.root / 'types.json', _study(cases, optics='absorbing_only'))
        result = _read(self._run('evaluate', path, 'out') / 'optical_results.json')
        for key, area in (('tri', 1), ('quad', 4), ('tri6', 1), ('quad8', 4)):
            with self.subTest(kind=key):
                force = result['loads'][key]['normal']['normalized']['force']
                np.testing.assert_allclose(force, [0, 0, -4 * area], atol=1e-12, rtol=0)

    def test_invalid_geometry_and_digests_are_per_case(self):
        cases = {'good': _bundle(self.root, 'good'),
                 'missing': {'role': 'clamped', 'geometry': {'path': 'missing.json', 'sha256': '0' * 64}},
                 'missing_npz': _bundle(self.root, 'missing_npz'),
                 'dtype': _bundle(self.root, 'dtype', array_changes={'node_xyz_m': np.zeros((3, 3), np.float32)}),
                 'unknown_node': _bundle(self.root, 'unknown_node', array_changes={
                     'connectivity_S3': np.array([[11, 18, 999]], dtype=np.int64)}),
                 'duplicate_node': _bundle(self.root, 'duplicate_node', array_changes={
                     'node_labels': np.array([1, 1, 2], dtype=np.int64)}),
                 'nonfinite': _bundle(self.root, 'nonfinite', ((0, 0, 0), (1, 0, 0), (1, 1, np.inf))),
                 'bad_sign': _bundle(self.root, 'bad_sign', array_changes={
                     'face_sign_S3': np.array([0], dtype=np.int8)}),
                 'empty_type': _bundle(self.root, 'empty_type', element_counts={'S3': 0}),
                 'bad_type': _bundle(self.root, 'bad_type', element_counts={'S7': 1}),
                 'sign_row_missing': _bundle(self.root, 'sign_row_missing', array_changes={'face_sign_S3': None}),
                 'zero_area': _bundle(self.root, 'zero_area', ((0, 0, 0), (1, 0, 0), (2, 0, 0))),
                 'swapped': _bundle(self.root, 'swapped'),
                 'wrong_kind': _bundle(self.root, 'wrong_kind', source_kind='odb'),
                 'bad_hash': _bundle(self.root, 'bad_hash'),
                 'bad_npz_hash': _bundle(self.root, 'bad_npz_hash')}
        # The edited manifest has a valid digest but belongs to another case.
        swapped = self.root / 'swapped' / 'geometry.json'
        _json(swapped, {**_read(swapped), 'case_id': 'different-case'})
        cases['swapped']['geometry'] = _ref(self.root, swapped)
        (self.root / 'unreadable').mkdir()
        cases['unreadable'] = {'role': 'clamped', 'geometry': {'path': 'unreadable', 'sha256': '0' * 64}}
        cases['bad_hash']['geometry']['sha256'] = 'f' * 64
        (self.root / 'bad_npz_hash' / 'geometry.npz').write_bytes(b'not the hashed bytes')
        (self.root / 'missing_npz' / 'geometry.npz').unlink()
        path = _json(self.root / 'cases.json', _study(cases))
        out = self._run('evaluate', path, 'out')
        evidence, result = _read(out / 'evidence.json'), _read(out / 'optical_results.json')
        for case_id in cases:
            status = {'good': 'qualified', 'missing': 'missing', 'missing_npz': 'missing'}.get(case_id, 'unverified')
            self.assertEqual(evidence['cases'][case_id]['status'], status, case_id)
            actual = result['loads'][case_id]['normal']['normalized']
            if case_id == 'good':
                self.assertAlmostEqual(actual['force'][2], -400 * (0.073744 + 2 * 0.926256))
            else:
                self.assertIsNone(actual, case_id)

    def _assert_malformed_archive_isolated(self, *, corrupt_deflate):
        path = write_malformed_input(self.root / 'input', corrupt_deflate=corrupt_deflate)
        npz = path.parent / 'bad' / 'geometry.npz'
        with zipfile.ZipFile(npz) as archive:
            if corrupt_deflate:
                with self.assertRaises(zlib.error):
                    archive.read('node_labels.npy')
            else:
                self.assertEqual(archive.read('node_labels.npy'), b'not an NPY array')
        inv = self._run('inventory', path, 'inv')
        ev = self._run('evaluate', path, 'ev')
        evidence = _read(inv / 'evidence.json')
        self.assertEqual(evidence, _read(ev / 'evidence.json'))
        self.assertEqual({key: row['status'] for key, row in evidence['cases'].items()},
                         {'usable': 'qualified', 'bad': 'unverified'})
        bad = evidence['cases']['bad']
        study_ref = _read(path)['cases']['bad']['geometry']
        geometry = path.parent / study_ref['path']
        self.assertEqual(bad['checked_files'], [
            {'path': str(geometry.resolve()), 'sha256': study_ref['sha256']},
            {'path': str(npz.resolve()), 'sha256': _read(geometry)['npz']['sha256']},
        ])
        self.assertEqual(bad['reasons'], [{'code': 'export_invalid', 'target': 'bad.geometry.npz'}])
        result = _read(ev / 'optical_results.json')
        for attitude, usable in result['loads']['usable'].items():
            self.assertEqual(usable['validation'], 'passed')
            self.assertEqual(result['loads']['bad'][attitude]['validation'], 'not_attempted')
            self.assertIsNone(result['loads']['bad'][attitude]['normalized'])
        np.testing.assert_allclose(result['loads']['usable']['normal']['normalized']['force'],
                                   (0, 0, -800), atol=4e-10, rtol=1e-10)

    def test_named_npy_with_non_array_bytes_does_not_abort_study(self):
        self._assert_malformed_archive_isolated(corrupt_deflate=False)

    def test_invalid_deflate_payload_does_not_abort_study(self):
        self._assert_malformed_archive_isolated(corrupt_deflate=True)

    def test_valid_mixed_synthetic_study_from_independent_commands(self):
        synthetic, _ = write_verification_inputs(self.root / 'inputs')
        inv = self._run('inventory', synthetic, 'inv')
        ev = self._run('evaluate', synthetic, 'ev')
        evidence = _read(ev / 'evidence.json')
        self.assertEqual(_read(inv / 'evidence.json'), evidence)
        self.assertFalse((inv / 'optical_results.json').exists())
        self.assertEqual({key: row['status'] for key, row in evidence['cases'].items()},
                         {'usable': 'qualified', 'bad': 'unverified', 'control': 'qualified'})
        self.assertEqual(evidence['contrasts']['usable']['status'], 'eligible')
        self.assertEqual(evidence['contrasts']['bad']['status'], 'blocked')
        self.assertEqual(evidence['screen']['status'], 'not_applicable')
        result = _read(ev / 'optical_results.json')
        self.assertEqual(result['study']['sha256'], _digest(synthetic))
        self.assertEqual(result['evidence']['sha256'], _digest(ev / 'evidence.json'))
        self.assertEqual(result['evaluator_version'], _digest(_SCRIPT))
        self.assertEqual(result['source_kind'], 'synthetic')
        self.assertEqual(result['layout_id'], 'nominal_square_v1')
        self.assertEqual(result['replication'], 'one_rotated_panel_shape')
        self.assertEqual(result['faceting_method'], 'piecewise_planar_triangles')
        self.assertEqual(result['faceting_version'], 1)
        self.assertEqual(result['moment_origin_m'], [0, 0, 0])
        self.assertEqual(result['units'], {'angle': 'deg', 'normalized_force': 'm2', 'normalized_moment': 'm3',
                                           'physical_force': 'N', 'physical_moment': 'N_m'})
        sine, cosine = math.sin(math.pi / 36), math.cos(math.pi / 36)
        expected_rays = {'normal': (0, 0, -1), 'alpha5_phi0': (-sine, 0, -cosine),
                         'alpha5_phi45': (-sine / math.sqrt(2), -sine / math.sqrt(2), -cosine)}
        self.assertEqual(result['rays'].keys(), expected_rays.keys())
        for key, expected in expected_rays.items():
            np.testing.assert_allclose(result['rays'][key], expected, rtol=0, atol=1e-12)
            load = result['loads']['usable'][key]
            np.testing.assert_allclose(load['normalized']['force'], (0, 0, -200 * expected[2]**2),
                                       rtol=1e-10, atol=4e-10)
            np.testing.assert_allclose(load['normalized']['moment'], (0, 0, 0), rtol=0, atol=8e-9)
            self.assertEqual(load['physical']['force'], [0, 0, 0])
            self.assertIsNone(result['loads']['bad'][key]['normalized'])
            for ref in ('flat', 'control'):
                pair = result['comparisons']['usable'][key][ref]
                self.assertAlmostEqual(pair['normalized']['force_norm'], -600 * expected[2]**2, delta=4e-10)
                self.assertIsNone(result['comparisons']['bad'][key][ref]['normalized'])
            flat = result['comparisons']['usable'][key]['flat']
            self.assertAlmostEqual(flat['angle_normalized_deg'], 0., delta=1e-12)
            self.assertIsNone(flat['angle_physical_deg'])
            self.assertIn({'code': 'zero_pressure', 'target': 'angle_physical_deg'}, flat['reasons'])
        self.assertEqual(result['completion']['screen']['status'], 'not_applicable')

    def test_all_missing_odb_paths_no_optical_output(self):
        _, initial = write_verification_inputs(self.root / 'inputs')
        inv = self._run('inventory', initial, 'inv')
        ev = self._run('evaluate', initial, 'ev')
        evidence = _read(inv / 'evidence.json')
        self.assertEqual(evidence, _read(ev / 'evidence.json'))
        self.assertEqual(evidence['source_kind'], 'odb')
        for row in evidence['cases'].values():
            self.assertEqual(row['status'], 'missing')
            self.assertEqual(row['checked_files'], [])
        self.assertTrue(all(row['status'] == 'blocked' for row in evidence['contrasts'].values()))
        result = _read(ev / 'optical_results.json')
        self.assertIn('flat', result['loads'])
        for by_attitude in result['loads'].values():
            for row in by_attitude.values():
                self.assertEqual(row['validation'], 'not_attempted')
                self.assertIsNone(row['normalized'])
                self.assertIsNone(row['physical'])
        self.assertEqual(result['completion']['screen']['status'], 'not_complete')

    def test_present_initial_odb_is_hashed_but_cannot_qualify_without_evidence(self):
        _, initial = write_verification_inputs(self.root / 'inputs')
        manifest = _read(initial)
        fake = initial.parent / 'fake.odb'
        fake.write_bytes(b'fabricated sentinel, not an ODB')
        manifest['cases']['usable']['odb_path'] = 'fake.odb'
        manifest['cases']['bad']['odb_path'] = '.'
        _json(initial, manifest)
        inv = self._run('inventory', initial, 'inv')
        ev = self._run('evaluate', initial, 'ev')
        result = _read(inv / 'evidence.json')
        self.assertEqual(result, _read(ev / 'evidence.json'))
        usable, bad = (result['cases'][k] for k in ('usable', 'bad'))
        self.assertEqual(usable['status'], 'unverified')
        self.assertEqual(usable['checked_files'], [{'path': str(fake.resolve()), 'sha256': _digest(fake)}])
        self.assertIn({'code': 'export_missing', 'target': 'usable.geometry'}, usable['reasons'])
        self.assertEqual(bad['status'], 'unverified')
        self.assertEqual(bad['checked_files'], [])
        self.assertEqual(bad['reasons'], [{'code': 'source_unreadable', 'target': 'bad.odb_path'}])

    def test_invalid_studies_exit_two_without_artifacts(self):
        path = self.root / 'study.json'
        for index, settings in enumerate(({'pressure_pa': -1}, {'pressure_pa': True}, {'pressure_pa': float('inf')},
                                          {'optics': 'custom'}, {'layout_id': 'wrong'}, {'extra': 1},
                                          {'schema_version': True})):
            with self.subTest(settings=settings):
                _json(path, _study({}, **settings))
                self._run('evaluate', path, f'out{index}', expected=2)
        path.write_text('{"schema_version":1,"cases":{"x":{},"x":{}}}', encoding='utf-8')
        self._run('inventory', path, 'out-duplicate', expected=2)
        _, initial = write_verification_inputs(self.root / 'inputs')
        for key, value in (('parent', {'path': 'initial.json', 'sha256': 'f' * 64}), ('optics', 'absorbing_only')):
            _json(path, {**_read(initial), key: value})
            self._run('inventory', path, 'out-' + key, expected=2)

    def test_cannot_overwrite_and_fatal_write_leaves_no_partial_results(self):
        synthetic, _ = write_verification_inputs(self.root / 'inputs')
        out = self._run('inventory', synthetic, 'out')
        digest = _digest(out / 'evidence.json')
        run = _cli('evaluate', synthetic, out)
        self.assertEqual(run.returncode, 2)
        self.assertEqual(digest, _digest(out / 'evidence.json'))
        self.assertFalse((out / 'optical_results.json').exists())
        broken = self.root / 'failed-publish'
        real_open = Path.open

        def reject_second(path, *args, **kwargs):
            if path == broken / 'optical_results.json':
                raise PermissionError('fabricated write error')
            return real_open(path, *args, **kwargs)

        with patch.object(Path, 'open', reject_second), self.assertRaises(PermissionError):
            optics._publish(broken, {'evidence.json': b'first', 'optical_results.json': b'second'})
        self.assertFalse(broken.exists())


class ExportTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.odb_path = self.root / 'fabricated.odb'
        self.odb_path.write_bytes(b'fabricated test bytes; never given to Abaqus')
        records = {'deck': ['submitted deck'], 'completed': ['completed attempt'],
                   'unseeded': ['reference anchors'], 'coating': ['reflective sides by element']}
        files = {key: _ref(self.root, _json(self.root / (key + '.json'), value))
                 for key, value in records.items()}
        self.run = dict(job_id='fake-job', attempt_id='attempt-1', driver_revision='rev1',
                        modified_source=None, submitted_deck='deck', completed=True,
                        instance='MEMBRANE-1', terminal_step='Step-GlobalTension',
                        step_duration=1., frame_index=1)
        self.basis = {'units': 'm', 'R': [[0, -1, 0], [1, 0, 0], [0, 0, 1]], 't_m': [4, -5, 1],
                      'anchors': [
                          {'node_label': 11, 'odb_xyz_m': [5, 4, -1], 'panel_xyz_m': [0, 0, 0]},
                          {'node_label': 33, 'odb_xyz_m': [5, -16, -1], 'panel_xyz_m': [20, 0, 0]},
                          {'node_label': 55, 'odb_xyz_m': [15, -6, -1], 'panel_xyz_m': [10, 10, 0]}],
                      'sides': {'55': 'SPOS', '41': 'SNEG', '10': 'SPOS', '3': 'SPOS', '2': 'SNEG'}}
        self.case = _initial_case('clamped', 'fabricated.odb')
        self.case.update(attempt_id='attempt-1', records=files,
                         run={'value': self.run, 'refs': [self.ref('deck'), self.ref('completed')]},
                         geometry_basis={'value': self.basis,
                                         'refs': [self.ref('unseeded'), self.ref('coating')]})
        self.manifest = _study(
            {'low': self.case, 'high': _initial_case('clamped', None), 'control': _initial_case('control', None)},
            {key: {'case_id': key, 'control_id': 'control'} for key in ('low', 'high')}, source_kind='odb')
        self.coords = {11: (0., 0., 0.), 33: (20., 0., 0.), 55: (10., 10., .125),
                       12: (0., 0., 0.), 13: (2., 0., 0.), 14: (2., 2., 0.), 15: (0., 2., 0.),
                       16: (1., 0., 0.), 17: (2., 1., 0.), 18: (1., 2., 0.), 19: (0., 1., 0.),
                       20: (1., 1., 0.)}
        self.install_geometry()

    @staticmethod
    def ref(record_id):
        return {'record_id': record_id, 'locator': 'fixture line 1'}

    def install_geometry(self):
        node = {label: SimpleNamespace(label=label) for label in self.coords}
        self.instance = SimpleNamespace(
            name='MEMBRANE-1', nodes=list(node.values()), elements=[
                SimpleNamespace(label=55, type='S3', connectivity=(node[11], 33, node[55])),
                SimpleNamespace(label=10, type='S4', connectivity=(12, 13, 14, 15)),
                SimpleNamespace(label=41, type='S3', connectivity=(11, 55, 33)),
                SimpleNamespace(label=3, type='STRI65', connectivity=(12, 13, 14, 16, 17, 20)),
                SimpleNamespace(label=2, type='S8R5', connectivity=(12, 13, 14, 15, 16, 17, 18, 19))])
        self.values = []
        self.nodal = object()
        for label, (x, y, z) in reversed(list(self.coords.items())):
            precision, prop = ('DOUBLE_PRECISION', 'dataDouble') if label == 55 else ('SINGLE_PRECISION', 'data')
            self.values.append(SimpleNamespace(nodeLabel=label, position='NODAL', instance=self.instance,
                                               precision=precision, **{prop: (y + 5, 4 - x, z - 1)}))
        integration_point = SimpleNamespace(position='INTEGRATION_POINT', instance=self.instance,
                                            precision='SINGLE_PRECISION', data=(4., 5., 6.))

        def subset(*, region, position=None):
            if region is not self.instance:
                raise AssertionError('wrong COORD region')
            if position is self.nodal:
                return SimpleNamespace(values=[v for v in self.values if v.position == 'NODAL'])
            if position is None:
                return SimpleNamespace(values=[integration_point, *self.values])
            raise AssertionError('expected the vendor NODAL symbol')

        self.integration_point = integration_point
        self.field = SimpleNamespace(getSubset=Mock(side_effect=subset))
        self.frame = SimpleNamespace(frameValue=.99, fieldOutputs={'COORD': self.field})
        self.step = SimpleNamespace(domain='TIME', timePeriod=1., frames=[
            SimpleNamespace(frameValue=0., fieldOutputs={'COORD': object()}), self.frame])
        self.odb = SimpleNamespace(
            steps={'Step-GlobalTension': self.step}, close=Mock(),
            rootAssembly=SimpleNamespace(instances={'MEMBRANE-1': self.instance,
                                                     'OTHER': SimpleNamespace(name='OTHER')}))
        self.open = Mock(return_value=self.odb)

    def invoke(self, *, code=0, out='out', manifest=None, supports_vendor=True):
        manifest = self.manifest if manifest is None else manifest
        path = _json(self.root / 'study.initial.json', manifest)
        destination = self.root / out
        runtime = {'abaqus': 'fabricated', 'python': 'fabricated', 'numpy': 'fabricated'}
        vendor = (patch.object(exporter, '_vendor', return_value=(self.open, self.nodal, runtime, ()))
                  if supports_vendor else patch.object(exporter, '_vendor', side_effect=ImportError('no vendor')))
        with vendor, patch('sys.stderr', new_callable=io.StringIO) as stderr:
            self.assertEqual(exporter.main(['--manifest', str(path), '--case', 'low',
                                            '--out-dir', str(destination)]), code, stderr.getvalue())
        if code == 1:
            self.assertFalse(destination.exists())
        return destination, stderr.getvalue()

    def test_fabricated_export_reads_final_instance_scoped_coords_and_exact_topology(self):
        before = _digest(self.odb_path)
        destination, _ = self.invoke()
        self.open.assert_called_once_with(path=str(self.odb_path.resolve()), readOnly=True)
        self.field.getSubset.assert_called_once_with(region=self.instance, position=self.nodal)
        self.odb.close.assert_called_once_with()
        self.assertEqual(_digest(self.odb_path), before)
        path = destination / 'geometry.json'
        meta = optics.read_export_geometry(path)
        self.assertEqual(meta['case_id'], 'low')
        self.assertEqual(meta['source_kind'], 'odb')
        self.assertEqual(meta['element_counts'], {'S3': 2, 'S4': 1, 'STRI65': 1, 'S8R5': 1})
        self.assertEqual(meta['export']['odb'], {'path': str(self.odb_path.resolve()), 'sha256': before})
        self.assertEqual((meta['export']['instance'], meta['export']['step'], meta['export']['step_duration'],
                          meta['export']['frame_index'], meta['export']['frame_value'],
                          meta['export']['coord_precision']),
                         ('MEMBRANE-1', 'Step-GlobalTension', 1, 1, .99, 'mixed'))
        canonical = json.dumps(self.basis, separators=(',', ':'), sort_keys=True, ensure_ascii=True)
        self.assertEqual(meta['export']['basis_sha256'], hashlib.sha256(canonical.encode('utf-8')).hexdigest())
        npz = destination / 'geometry.npz'
        with np.load(npz, allow_pickle=False) as bundle:
            np.testing.assert_array_equal(bundle['node_labels'], list(self.coords))
            np.testing.assert_array_equal(bundle['node_xyz_m'], list(self.coords.values()))
            np.testing.assert_array_equal(bundle['element_labels_S3'], [55, 41])
            np.testing.assert_array_equal(bundle['connectivity_S3'], [[11, 33, 55], [11, 55, 33]])
            np.testing.assert_array_equal(bundle['face_sign_S3'], [1, -1])
            np.testing.assert_array_equal(bundle['connectivity_STRI65'], [[12, 13, 14, 16, 17, 20]])
            np.testing.assert_array_equal(bundle['connectivity_S8R5'], [[12, 13, 14, 15, 16, 17, 18, 19]])
            np.testing.assert_array_equal(bundle['face_sign_S8R5'], [-1])
        npz.write_bytes(b'changed fabricated bundle')
        with self.assertRaisesRegex(optics.InputError, 'NPZ digest mismatch'):
            optics.read_export_geometry(path)

    def test_identity_mapping_and_both_uniform_precision_paths(self):
        self.basis.update(R=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], t_m=[0, 0, 0])
        for anchor in self.basis['anchors']:
            anchor['odb_xyz_m'] = anchor['panel_xyz_m']
        for name, symbol, field in (('single', 'SINGLE_PRECISION', 'data'),
                                    ('double', 'DOUBLE_PRECISION', 'dataDouble')):
            with self.subTest(precision=name):
                self.values[:] = [SimpleNamespace(nodeLabel=label, position='NODAL', instance=self.instance,
                                                  precision=symbol, **{field: xyz})
                                  for label, xyz in self.coords.items()]
                out, _ = self.invoke(out=name)
                self.assertEqual(optics.read_export_geometry(out / 'geometry.json')['export']['coord_precision'],
                                 name)
                with np.load(out / 'geometry.npz', allow_pickle=False) as mesh:
                    np.testing.assert_array_equal(mesh['node_xyz_m'], list(self.coords.values()))

    def test_filtered_coord_does_not_fill_a_missing_node_from_integration_points(self):
        self.integration_point.nodeLabel = self.values.pop().nodeLabel
        _, message = self.invoke(code=1)
        self.assertIn('COORD does not cover every membrane node', message)
        self.field.getSubset.assert_called_once_with(region=self.instance, position=self.nodal)
        self.odb.close.assert_called_once_with()

    def test_defensive_position_check_rejects_invalid_filtered_result(self):
        self.field.getSubset.side_effect = None
        self.field.getSubset.return_value = SimpleNamespace(values=[self.integration_point, *self.values])
        _, message = self.invoke(code=1)
        self.assertIn('COORD value is not nodal in the requested membrane instance', message)
        self.field.getSubset.assert_called_once_with(region=self.instance, position=self.nodal)
        self.odb.close.assert_called_once_with()

    def test_no_export_for_wrong_step_instance_or_frame(self):
        for index, mutate in enumerate((
            lambda: self.odb.steps.update(Later=self.step),
            lambda: setattr(self.step, 'domain', 'FREQUENCY'),
            lambda: setattr(self.step, 'timePeriod', 2.),
            lambda: setattr(self.step, 'frames', []),
            lambda: setattr(self.frame, 'frameValue', .989999),
            lambda: setattr(self.frame, 'frameValue', math.inf),
            lambda: self.frame.fieldOutputs.clear(),
            lambda: self.odb.rootAssembly.instances.pop('MEMBRANE-1'),
            lambda: self.run.update(frame_index=0),
        )):
            with self.subTest(index=index):
                self.run['frame_index'] = 1
                self.install_geometry()
                mutate()
                self.invoke(code=1, out=f'rejected-{index}')
                self.odb.close.assert_called_once_with()

    def test_does_not_fall_back_to_previous_coord_frame(self):
        self.step.frames[0].fieldOutputs['COORD'] = self.field
        self.frame.fieldOutputs.clear()
        _, message = self.invoke(code=1)
        self.assertIn('final frame 1 does not contain COORD', message)
        self.field.getSubset.assert_not_called()

    def test_unresolved_or_contradictory_inputs_fail_before_vendor_access(self):
        corrections = (
            (self.run, {'step_duration': 2.}),
            (self.run, {'completed': False}),
            (self.run, {'attempt_id': 'unselected'}),
            (self.basis, {'R': [[0, -1, 0], [1, 0, 0], [0, 0, -1]]}),
            (self.basis, {'t_m': [4, -5, 1.01]}),
            (self.basis, {'units': 'mm'}),
            (self.case, {'geometry_basis': None}),
            (self.case['records'], {'deck': {'path': 'deck.json', 'sha256': 'f' * 64}}),
        )
        for index, (obj, overrides) in enumerate(corrections):
            with self.subTest(index=index), patch.dict(obj, overrides):
                self.invoke(code=1, out=f'rejected-{index}')
                self.open.assert_not_called()

    def test_invalid_instance_geometry_is_not_partially_exported(self):
        for index, mutate in enumerate((
            lambda: self.instance.nodes.clear(),
            lambda: self.instance.elements.clear(),
            lambda: self.values.pop(),
            lambda: self.values.append(self.values[0]),
            lambda: setattr(self.values[0], 'instance', SimpleNamespace(name='OTHER')),
            lambda: setattr(self.values[0], 'precision', 'UNKNOWN_PRECISION'),
            lambda: setattr(self.values[0], 'data', (0, 0, math.nan)),
            lambda: setattr(self.instance.elements[0], 'type', 'S9'),
            lambda: setattr(self.instance.elements[0], 'connectivity', (11, 33, 33)),
            lambda: self.basis['sides'].pop('10'),
        )):
            with self.subTest(index=index):
                self.install_geometry()
                mutate()
                self.invoke(code=1, out=f'rejected-{index}')
                self.odb.close.assert_called_once_with()

    def test_changed_odb_identity_is_rejected_after_close(self):
        self.odb.close.side_effect = lambda: self.odb_path.write_bytes(b'external fabricated edit')
        _, message = self.invoke(code=1)
        self.assertIn('ODB SHA-256 changed', message)
        self.odb.close.assert_called_once_with()

    def test_output_errors_and_vendor_unavailability_never_publish_a_bundle(self):
        out = self.root / 'existing'
        out.mkdir()
        sentinel = out / 'geometry.json'
        sentinel.write_bytes(b'leave me alone')
        self.invoke(code=2, out=out.name)
        self.assertEqual(sentinel.read_bytes(), b'leave me alone')
        self.invoke(code=2, out='synthetic', manifest=_study({}))
        _, message = self.invoke(code=1, out='no-vendor', supports_vendor=False)
        self.assertIn('no vendor', message)
        self.open.assert_not_called()


def _review(value, *records):
    return {'value': value, 'refs': [{'record_id': record, 'locator': 'fabricated test facts'} for record in records]}


def _basis_digest(basis):
    return hashlib.sha256(json.dumps(basis, sort_keys=True, separators=(',', ':')).encode('ascii')).hexdigest()


def _support(node, dof=3, kind='membrane_node'):
    return {'target_kind': kind, 'target': node, 'step': 'Initial', 'dof': dof, 'law': 'fixed',
            'parameters': {'value_m': {'value': '0', 'unit': 'm'}}}


def _fabricated_real(root, *, copy_high=False):
    """Exercise the ODB contract with explicitly fabricated bytes; no vendor API."""
    root.mkdir()
    refmesh = root / 'fabricated-reference.txt'
    refmesh.write_bytes(b'fabricated reference mesh with labels 11, 18 and 25')
    mesh = _digest(refmesh)
    cases = {}
    for key, vertices, ratio, stress in (('control', _PANEL, None, None), ('low', _CONTROL, '0.01', '1.71'),
                                        ('high', _SMALL, '0.02', '3.42')):
        odb = root / (key + '.odb')
        odb.write_bytes(b'fabricated opaque product ' + key.encode())
        deck = root / (key + '.inp')
        deck.write_bytes(b'fabricated submitted deck ' + key.encode())
        record = root / (key + '-run.txt')
        record.write_bytes(b'fabricated job, attempt, field, and copy facts ' + key.encode())
        case = _initial_case('control' if ratio is None else 'clamped', odb.name)
        case['attempt_id'] = key + '-1'
        run = {'job_id': key + '-job', 'attempt_id': case['attempt_id'], 'driver_revision': 'fixture-rev',
               'modified_source': None, 'submitted_deck': 'deck', 'completed': True, 'instance': 'M-1',
               'terminal_step': 'Finish', 'step_duration': 1, 'frame_index': 2}
        case['records'] = {'deck': _ref(root, deck), 'job': _ref(root, record), 'reference': _ref(root, refmesh)}
        case['run'] = _review(run, 'deck', 'job')
        current_file = _ref(root, odb)
        product = {'job_id': run['job_id'], 'attempt_id': run['attempt_id'], 'product_id': 'product',
                   'product_sha256': current_file['sha256'], 'current_id': 'product',
                   'current_file': current_file, 'copies': []}
        if key == 'high' and copy_high:
            product.update(product_sha256=None, current_id='final', copies=[
                {'from_id':'product','to_id':'backup','sha256':None,'refs':_review(None,'job')['refs']},
                {'from_id':'backup','to_id':'final','sha256':current_file['sha256'],'refs':_review(None,'job')['refs']}])
        case['producer'] = _review(product, 'job')
        basis = {'units': 'm', 'R': np.eye(3, dtype=int).tolist(), 't_m': [0,0,0],
                 'anchors': [{'node_label': n, 'odb_xyz_m': list(p), 'panel_xyz_m': list(p)}
                             for n,p in zip((11,18,25),_PANEL)], 'sides': {'101': 'SPOS'}}
        exp = {'odb': {**current_file, 'path': '../'+odb.name}, 'job_id': run['job_id'],
               'attempt_id': run['attempt_id'], 'instance': run['instance'], 'step': run['terminal_step'],
               'step_duration': 1, 'frame_index': 2, 'frame_value': 1, 'coord_precision': 'single',
               'runtime': {'abaqus':'fixture','python':'fixture','numpy':'fixture'},
               'adapter_sha256': 'c'*64, 'basis_sha256': _basis_digest(basis)}
        case['geometry'] = _bundle(root, key, vertices, source_kind='odb', export=exp)['geometry']
        receipt = root / key / 'export_verification.json'
        checks = {name: {'status': 'passed', 'evidence': [{'file': {**_ref(root, record), 'path': '../'+record.name},
                 'locator':'fabricated verification proof'}]} for name in ('terminal_frame','anchors','coordinates',
                 'topology','face_signs','portable_readback','odb_unchanged','negative_inputs')}
        _json(receipt, {'schema_version':1, 'geometry':{**case['geometry'],'path':'geometry.json'},
                        'mapping_kind':'identity','node_count':3,'element_counts':{'S3':1},
                        'status':'passed','unexercised':[], 'checks':checks,
                        **{k:v for k,v in exp.items() if k not in ('job_id','attempt_id','basis_sha256')}})
        case['records']['roundtrip'] = _ref(root, receipt)
        case['geometry_basis'] = _review(basis, 'reference', 'roundtrip')
        clamp = 'absent' if ratio is None else {'x_c':'0.7','d_c':ratio,'corner_initial_pa':'171.00',
                 'clamp_initial_pa':stress,'hardware':{'left':{'cable':'cl'},'right':{'cable':'cr'}},
                 'connections':{'left':{'reference_point':'CL'},'right':{'reference_point':'CR'}}}
        supports = {'common': [_support(18), _support(25)] + ([_support(11)] if ratio is None else []),
                    'replaced_common': [] if ratio is None else [_support(11)],
                    'attachment': [] if ratio is None else [_support('CL', kind='reference_point')],
                    'attachment_patch': [] if ratio is None else [11]}
        settings = {'membrane': {'reference_mesh_sha256':mesh,'node_labels':[11,18,25],
                                 'material':{'id':'fixture'},'thickness_m':'0.000005'},
                    'clamp': clamp, 'supports': supports,
                    'load': {'mode':'cable_stress','corner_initial_pa':dict.fromkeys(('apex','first_base','second_base'),'171'),
                             'shared_steps':[{'name':'Finish','controls':{'period':{'value':'1.0','unit':'s'}}}]},
                    'solver': {'shared_options':{'nlgeom':True}, 'applied_stabilization_by_step':{'Finish':{'value':'0.0002','unit':'1'}}},
                    'perturbation': {'reference_mesh_sha256':mesh,'frame':'panel_local_m','node_map':[[n,n] for n in (11,18,25)],
                                     'delta_m':[[11,'0','0','0.00010'],[18,'0','-0','0'],[25,'0','0','0']],
                                     'precision':{'encoding':'exact_decimal','basis':'fabricated exact input coordinates'}}}
        case['settings'] = {k:_review(v,'job', 'reference') for k,v in settings.items()}
        cases[key] = case
    return _study(cases, {k:{'case_id':k,'control_id':'control'} for k in ('low','high')}, source_kind='odb')


def _freeze_real(root, study):
    initial = deepcopy(study)
    for case in initial['cases'].values():
        case['geometry'] = None
    parent = _json(root/'study.initial.json', initial)
    return _json(root/'study.linked.json', {**study, 'parent': _ref(root, parent)})


class RealEvidenceTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)

    def commands(self, path):
        for command in ('inventory', 'evaluate'):
            result = _cli(command, path, path.parent / command)
            self.assertEqual(result.returncode, 0, result.stderr)
        evidence = _read(path.parent/'inventory'/'evidence.json')
        self.assertEqual(evidence, _read(path.parent/'evaluate'/'evidence.json'))
        return evidence, _read(path.parent/'evaluate'/'optical_results.json')

    def test_real_contract_direct_and_documented_copy_chain(self):
        root = self.root/'valid'
        study = _fabricated_real(root, copy_high=True)
        study['cases']['high']['settings']['perturbation']['value']['delta_m'].reverse()
        study['cases']['high']['settings']['load']['value']['corner_initial_pa']['apex'] = '1.71e2'
        evidence, result = self.commands(_freeze_real(root, study))
        self.assertEqual([r['status'] for r in evidence['cases'].values()], ['qualified'] * 3)
        self.assertEqual(evidence['screen']['status'], 'eligible')
        self.assertEqual(result['completion']['screen']['status'], 'complete')
        self.assertEqual(result['source_kind'], 'odb')
        self.assertIsNone(result['pressure_pa'])
        load = result['loads']['low']['normal']
        np.testing.assert_allclose(load['normalized']['force'], [0,0,-100*(.073744+2*.926256)], atol=4e-10, rtol=0)
        self.assertIsNone(load['physical'])
        self.assertEqual(evidence['cases']['low']['diagnostics']['status'], 'unresolved')

    def test_broken_producer_chain_isolated_even_when_current_digest_matches(self):
        root = self.root/'broken'
        study = _fabricated_real(root, copy_high=True)
        study['cases']['high']['producer']['value']['copies'][1]['from_id'] = 'wrong-intermediate'
        (root/'low-run.txt').write_bytes(b'changed supporting record')
        e,r = self.commands(_freeze_real(root, study))
        self.assertEqual(e['cases']['control']['status'], 'qualified')
        self.assertEqual(e['cases']['high']['status'], 'unverified')
        self.assertIn({'code':'producer_unresolved','target':'high.producer'}, e['cases']['high']['reasons'])
        self.assertIn({'code':'record_conflict','target':'low.run.job'}, e['cases']['low']['reasons'])
        self.assertEqual(r['loads']['high']['normal']['validation'], 'not_attempted')

    def test_final_digest_without_origin_and_altered_odb_do_not_qualify(self):
        root = self.root/'no-origin'
        study = _fabricated_real(root)
        high = study['cases']['high']
        high['records']['final_only'] = high['producer']['value']['current_file']
        high['producer']['refs'] = _review(None,'final_only')['refs']
        (root/'low.odb').write_bytes(b'other current opaque bytes')
        e,_ = self.commands(_freeze_real(root,study))
        self.assertIn({'code':'producer_unresolved','target':'high.producer'}, e['cases']['high']['reasons'])
        self.assertIn({'code':'record_conflict','target':'low.producer'}, e['cases']['low']['reasons'])

    def test_real_export_rejects_conflicting_frame_attempt_basis_and_face_rows(self):
        def change_export(root, case, **values):
            p = root / case['geometry']['path']
            doc = _read(p)
            doc['export'].update(values)
            _json(p,doc)
            case['geometry'] = _ref(root,p)
        for field,value in [('frame_index',1),('frame_value',.9),('attempt_id','other'),('basis_sha256','1'*64)]:
            with self.subTest(field=field):
                root = self.root/field
                study = _fabricated_real(root)
                change_export(root,study['cases']['high'],**{field:value})
                low = study['cases']['low']
                basis = low['geometry_basis']['value']
                basis['sides']['101'] = 'SNEG'
                change_export(root,low,basis_sha256=_basis_digest(basis))
                e,r = self.commands(_freeze_real(root,study))
                self.assertEqual(e['cases']['control']['status'], 'qualified')
                self.assertEqual(e['cases']['high']['status'], 'unverified')
                self.assertIn({'code':'record_conflict','target':'low.geometry_basis.sides'},e['cases']['low']['reasons'])
                self.assertIsNone(r['loads']['high']['normal']['normalized'])

    def test_invalid_mapping_and_uncovered_native_precision_withhold_loads(self):
        root = self.root/'coverage'
        study = _fabricated_real(root)
        study['cases']['high']['geometry_basis']['value']['R'][2][2] = -1
        low = study['cases']['low']
        p = root/low['records']['roundtrip']['path']
        _json(p,{**_read(p),'coord_precision':'double'})
        low['records']['roundtrip'] = _ref(root,p)
        e,_ = self.commands(_freeze_real(root,study))
        self.assertIn({'code':'record_conflict','target':'low.export_verification'}, e['cases']['low']['reasons'])
        self.assertIn({'code':'record_conflict','target':'high.geometry_basis'}, e['cases']['high']['reasons'])

    def test_unknown_field_or_restraint_and_known_difference_preserve_other_pair(self):
        for mode in ('unknown','field','support','outside_patch'):
            with self.subTest(mode=mode):
                root = self.root/mode
                study = _fabricated_real(root)
                high = study['cases']['high']['settings']
                if mode == 'field':
                    high['perturbation']['value']['delta_m'][0][-1] = '0.00010000000000000000001'
                    high['solver'] = None
                elif mode == 'support':
                    high['supports']['value']['common'].pop()
                    high['perturbation'] = None
                elif mode == 'outside_patch':
                    high['supports']['value']['replaced_common'] = [_support(25)]
                    high['supports']['value']['common'].pop()
                else:
                    high['perturbation']['value']['node_map'].pop()
                e,r = self.commands(_freeze_real(root,study))
                self.assertEqual(e['cases']['high']['status'], 'qualified')
                self.assertEqual(e['contrasts']['high']['status'], 'unmatched' if mode in ('field','support') else 'blocked')
                self.assertEqual(e['contrasts']['low']['status'], 'eligible')
                self.assertEqual(r['completion']['contrasts']['low']['status'], 'complete')
                self.assertEqual(r['completion']['screen']['status'], 'partial')
                for row in r['comparisons']['high'].values():
                    self.assertIsNone(row['control']['normalized'])
                    self.assertIsNotNone(row['flat']['normalized'])

    def test_two_complete_pairs_do_not_imply_eligible_screen(self):
        root = self.root/'different-clamp'
        study = _fabricated_real(root)
        settings = study['cases']['high']['settings']
        settings['clamp']['value']['hardware']['left']['cable'] = 'other'
        settings['supports']['value']['attachment'][0]['dof'] = 2
        e,r = self.commands(_freeze_real(root,study))
        self.assertEqual(e['screen']['status'], 'unmatched')
        self.assertEqual({n['target'] for n in e['screen']['reasons']},
                         {'screen.clamp.hardware','screen.attachment'})
        self.assertTrue(all(x['status']=='complete' for x in r['completion']['contrasts'].values()))
        self.assertEqual(r['completion']['screen']['status'], 'partial')

    def test_eligible_pair_requires_passing_operands_at_every_attitude(self):
        root = self.root/'failed-load'
        path = _freeze_real(root,_fabricated_real(root))
        inv, ev = root/'inventory', root/'evaluate'
        self.assertEqual(_cli('inventory',path,inv).returncode,0)
        route = optics._route
        def fail_one(facets, ray, *args, direct):
            load = route(facets,ray,*args,direct=direct)
            if direct and np.max(facets.xyz[:,:,0]) == 10 and ray[0] < 0 and ray[1] == 0:
                load[0,0,2] += 1
            return load
        with patch.object(optics,'_route',side_effect=fail_one):
            self.assertEqual(optics.main(['evaluate','--manifest',str(path),'--out-dir',str(ev)]),0)
        self.assertEqual(_read(inv/'evidence.json'),_read(ev/'evidence.json'))
        r = _read(ev/'optical_results.json')
        self.assertEqual(r['loads']['low']['alpha5_phi0']['validation'], 'failed')
        self.assertIsNone(r['comparisons']['low']['alpha5_phi0']['control']['normalized'])
        self.assertIsNotNone(r['comparisons']['low']['alpha5_phi45']['control']['normalized'])
        self.assertEqual(r['completion']['contrasts']['high']['status'], 'complete')
        self.assertEqual(r['completion']['contrasts']['low']['status'], 'not_complete')
        self.assertIn('route_mismatch',{n['code'] for n in r['completion']['contrasts']['low']['reasons']})

    def test_initial_parent_digest_and_resolved_selection_must_agree(self):
        root = self.root/'parent'
        p = _freeze_real(root,_fabricated_real(root))
        linked = _read(p)
        linked['cases']['high']['producer']['value']['current_file']['path'] = str(root/'high.odb')
        _json(p,linked)
        e,_ = self.commands(p)
        self.assertEqual(e['screen']['status'], 'eligible')
        linked['cases']['high']['attempt_id'] = 'substituted'
        _json(p,linked)
        out = root/'changed-selection'
        actual = _cli('inventory',p,out)
        self.assertEqual(actual.returncode,2)
        self.assertIn('frozen initial selection',actual.stderr)
        self.assertFalse(out.exists())
        (root/'study.initial.json').write_bytes(b'not the originally approved bytes')
        out = root/'bad-parent-digest'
        actual = _cli('evaluate',p,out)
        self.assertEqual(actual.returncode,2)
        self.assertIn('parent digest mismatch',actual.stderr)
        self.assertFalse(out.exists())


if __name__ == '__main__':
    unittest.main()
