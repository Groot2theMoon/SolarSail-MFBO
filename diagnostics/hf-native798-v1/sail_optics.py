"""Portable whole-sail optics and path-only initial ODB inventory (Part A)."""

import argparse
from dataclasses import dataclass
import hashlib
import json
import lzma
import math
from pathlib import Path
import re
import stat
import sys
import zipfile
import zlib

import numpy as np


_LAYOUT = 'nominal_square_v1'
_OPTICS = {'default': (0.073744, 0.926256), 'absorbing_only': (1., 0.),
           'reflecting_only': (0., 1.)}
_WIDTHS = {'S3': 3, 'S4': 4, 'STRI65': 6, 'S8R5': 8}
_FACETS = {'S3': ((0, 1, 2),), 'S4': ((0, 1, 2), (0, 2, 3)),
           'STRI65': ((0, 3, 5), (3, 1, 4), (5, 4, 2), (3, 4, 5)),
           'S8R5': ((0, 4, 8), (0, 8, 7), (4, 1, 5), (4, 5, 8),
                    (8, 5, 2), (8, 2, 6), (7, 8, 6), (7, 6, 3))}
_FLAT = np.array([[[0., 0., 0.], [20., 0., 0.], [10., 10., 0.]]])
_ORIGIN = np.zeros(3)
_LOCAL_ORIGIN = np.array([10., 10. / 3., 0.])
_DELTA_FIELDS = ('force', 'moment', 'force_norm', 'moment_norm',
                 'ray_projection', 'transverse_force')
_SOLVER_ROOTS = tuple(Path(__file__).resolve().parents[2] / rel
                      for rel in ('results', 'code/results'))


class InputError(ValueError):
    """An invalid contract, not a failed numerical route."""


@dataclass(frozen=True)
class Facets:
    xyz: np.ndarray  # (facet, vertex, component)
    signs: np.ndarray


def _keys(value, required, optional=()):
    if not isinstance(value, dict):
        raise InputError('expected a JSON object')
    missing = set(required) - value.keys()
    extra = value.keys() - set(required) - set(optional)
    if missing or extra:
        raise InputError(f'missing fields {sorted(missing)}; unknown fields {sorted(extra)}')
    return value


def _object_pairs(pairs):
    obj = {}
    for key, value in pairs:
        if key in obj:
            raise InputError(f'duplicate JSON key: {key}')
        obj[key] = value
    return obj


def _json_float(text):
    value = float(text)
    if not math.isfinite(value):
        raise InputError(f'nonfinite JSON number: {text}')
    return value


def _decode(raw):
    return json.loads(raw.decode('utf-8'), object_pairs_hook=_object_pairs,
                      parse_float=_json_float, parse_constant=_json_float)


def _version(obj):
    if type(obj['schema_version']) is not int or obj['schema_version'] != 1:
        raise InputError('expected schema_version 1')


def _text(value, name):
    if not isinstance(value, str) or not value or '\0' in value:
        raise InputError(f'{name} must be a nonempty string')
    return value


def _file_ref(ref):
    _keys(ref, ('path', 'sha256'))
    _text(ref['path'], 'path')
    if not isinstance(ref['sha256'], str) or not re.fullmatch('[0-9a-f]{64}', ref['sha256']):
        raise InputError('sha256 must contain 64 lowercase hexadecimal digits')


def _identity(path, digest):
    return {'path': str(path), 'sha256': digest}


def _stream_sha256(stream):
    digest = hashlib.sha256()
    for chunk in iter(lambda: stream.read(1024 * 1024), b''):
        digest.update(chunk)
    return digest.hexdigest()


def _referent(base, ref):
    return (base / ref['path']).resolve()


def _enum(value, allowed, name):
    if not isinstance(value, str) or value not in allowed:
        raise InputError(f'unsupported {name}: {value!r}')
    return value


def _initial_candidate(case):
    # Only validate the intake envelope. B/C fact payloads are not interpreted or opened.
    if case['geometry'] is not None:
        raise InputError('Part A requires null geometry for initial ODB candidates')
    for key in ('odb_path', 'attempt_id'):
        if case[key] is not None:
            _text(case[key], key)
    if not isinstance(case['records'], dict):
        raise InputError('records must be an object')
    for record_id, ref in case['records'].items():
        _text(record_id, 'record ID')
        _file_ref(ref)
    settings = _keys(case['settings'], ('membrane', 'load', 'solver', 'supports', 'clamp',
                                       'perturbation'))
    for item in (case['run'], case['producer'], case['geometry_basis'], *settings.values()):
        if item is not None:
            _keys(item, ('value', 'refs'))
            if not isinstance(item['refs'], list) or not item['refs']:
                raise InputError('reviewed values require source references')
            for ref in item['refs']:
                _keys(ref, ('record_id', 'locator'))
                for key in ref:
                    _text(ref[key], key)
    diagnostic = _keys(case['diagnostics'], ('status', 'observations', 'reason'))
    _enum(diagnostic['status'], ('recorded', 'unresolved'), 'diagnostics status')
    if not isinstance(diagnostic['observations'], list):
        raise InputError('diagnostic observations must be an array')
    if diagnostic['reason'] is not None or diagnostic['status'] == 'unresolved':
        _text(diagnostic['reason'], 'diagnostics reason')


def _study(raw):
    data = _keys(_decode(raw), ('schema_version', 'source_kind', 'layout_id', 'optics',
                                'cases', 'contrasts'), ('pressure_pa', 'parent'))
    _version(data)
    kind = _enum(data['source_kind'], ('synthetic', 'odb'), 'source kind')
    if 'parent' in data:
        raise InputError('linked ODB studies are unavailable in Part A')
    _enum(data['layout_id'], (_LAYOUT,), 'layout')
    _enum(data['optics'], _OPTICS if kind == 'synthetic' else ('default',), 'optics')
    if 'pressure_pa' in data:
        pressure = data['pressure_pa']
        if type(pressure) not in (int, float) or not math.isfinite(pressure) or pressure < 0:
            raise InputError('pressure_pa must be a finite nonnegative number')
    cases, contrasts = data['cases'], data['contrasts']
    if not isinstance(cases, dict) or not isinstance(contrasts, dict):
        raise InputError('cases and contrasts must be objects')
    for case_id, case in cases.items():
        _text(case_id, 'case ID')
        if case_id == 'flat':
            raise InputError('flat is a reserved case ID')
        fields = ('role', 'geometry')
        if kind == 'odb':
            fields += ('odb_path', 'attempt_id', 'records', 'run', 'producer', 'settings',
                       'geometry_basis', 'diagnostics')
        _keys(case, fields)
        _enum(case['role'], ('control', 'clamped'), 'case role')
        if kind == 'synthetic':
            _file_ref(case['geometry'])
        else:
            _initial_candidate(case)
    paired = set()
    for contrast_id, pair in contrasts.items():
        _text(contrast_id, 'contrast ID')
        _keys(pair, ('case_id', 'control_id'))
        for key, role in (('case_id', 'clamped'), ('control_id', 'control')):
            member = _text(pair[key], key)
            if member not in cases or cases[member]['role'] != role:
                raise InputError(f'{contrast_id}: {key} must name a {role} case')
        if pair['case_id'] in paired:
            raise InputError(f"duplicate contrast for {pair['case_id']}")
        paired.add(pair['case_id'])
    if kind == 'odb' and (len(cases) != 3 or len(paired) != 2):
        raise InputError('initial ODB intake requires one control, two clamped cases, and both contrasts')
    return data


def _split(points, kind):
    if kind == 'S8R5':
        centre = (points[:, :4] * 0.25).sum(axis=1, keepdims=True)
        points = np.concatenate((points, centre), axis=1)
    return points[:, _FACETS[kind], :].reshape(-1, 3, 3)


def _array(bundle, key, dtype, shape):
    value = bundle[key]
    if not isinstance(value, np.ndarray):
        raise InputError(f'{key}: expected an NPY array')
    if value.dtype != np.dtype(dtype) or value.shape != shape:
        raise InputError(f'{key}: expected {dtype} array of shape {shape}')
    return value


def _crosses(xyz):
    return np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])


def _geometry(bundle, n, counts, *, keep_facets=True):
    _mesh_counts(n, counts)
    expected = {'node_labels', 'node_xyz_m'}
    for kind in counts:
        expected.update(f'{name}_{kind}' for name in ('element_labels', 'connectivity', 'face_sign'))
    keys = list(bundle.keys())
    if len(keys) != len(expected) or set(keys) != expected:
        raise InputError('NPZ groups do not exactly match the declared element types')
    labels = _array(bundle, 'node_labels', 'int64', (n,))
    coords = _array(bundle, 'node_xyz_m', 'float64', (n, 3))
    order = np.argsort(labels)
    sorted_labels = labels[order]
    if np.any(sorted_labels[1:] == sorted_labels[:-1]) or not np.isfinite(coords).all():
        raise InputError('duplicate node labels or nonfinite node coordinates')
    xyz, signs, element_labels = [], [], []
    for kind, count in counts.items():
        element_labels.append(_array(bundle, f'element_labels_{kind}', 'int64', (count,)))
        conn = _array(bundle, f'connectivity_{kind}', 'int64', (count, _WIDTHS[kind]))
        sign = _array(bundle, f'face_sign_{kind}', 'int8', (count,))
        if not np.all((sign == 1) | (sign == -1)):
            raise InputError(f'{kind}: invalid face sign')
        idx = sorted_labels.searchsorted(conn)
        if np.any(idx >= n) or not np.array_equal(sorted_labels[idx], conn):
            raise InputError(f'{kind}: connectivity references an unknown label')
        with np.errstate(over='ignore', invalid='ignore'):
            triangles = _split(coords[order[idx]], kind)
            if np.any(np.all(_crosses(triangles) == 0, axis=1)):
                raise InputError(f'{kind}: zero-area facet')
        if keep_facets:
            xyz.append(triangles)
            signs.append(np.repeat(sign, len(_FACETS[kind])))
    all_elements = np.concatenate(element_labels)
    if np.unique(all_elements).size != all_elements.size:
        raise InputError('duplicate element labels')
    if keep_facets:
        return Facets(np.concatenate(xyz), np.concatenate(signs))


def _note(code, target):
    return {'code': code, 'target': target}


def _sorted_notes(notes):
    return [_note(code, target) for code, target in sorted({(n['code'], n['target']) for n in notes})]


def _mesh_counts(n, counts):
    if type(n) is not int or n <= 0 or not isinstance(counts, dict) or not counts:
        raise InputError('geometry must declare a nonempty mesh')
    if any(kind not in _WIDTHS or type(count) is not int or count <= 0
           for kind, count in counts.items()):
        raise InputError('invalid element types or declared counts')


def _geometry_manifest(raw):
    meta = _keys(_decode(raw), ('schema_version', 'source_kind', 'case_id', 'frame',
                               'npz', 'node_count', 'element_counts', 'export'))
    _version(meta)
    _text(meta['case_id'], 'case ID')
    _enum(meta['frame'], ('panel_local_m',), 'geometry frame')
    _mesh_counts(meta['node_count'], meta['element_counts'])
    return meta


def _read_npz(meta, base, checked, *, keep_facets=True):
    ref = meta['npz']
    _file_ref(ref)
    path = _referent(base, ref)
    with path.open('rb') as stream:
        digest = _stream_sha256(stream)
        checked.append(_identity(path, digest))
        if digest != ref['sha256']:
            raise InputError('NPZ digest mismatch')
        stream.seek(0)
        with zipfile.ZipFile(stream) as archive:
            names = archive.namelist()
            if len(set(names)) != len(names) or any(not name.endswith('.npy') for name in names):
                raise InputError('invalid NPZ members')
        stream.seek(0)
        with np.load(stream, allow_pickle=False) as bundle:
            return _geometry(bundle, meta['node_count'], meta['element_counts'],
                             keep_facets=keep_facets)


def read_export_geometry(path):
    """Read back a v1 adapter bundle, without opening its ODB or qualifying it."""
    path = Path(path).resolve()
    meta = _geometry_manifest(path.read_bytes())
    _enum(meta['source_kind'], ('odb',), 'export source kind')
    export = _keys(meta['export'], ('odb', 'job_id', 'attempt_id', 'instance', 'step',
                                   'step_duration', 'frame_index', 'frame_value', 'coord_precision',
                                   'basis_sha256', 'runtime', 'adapter_sha256'))
    _file_ref(export['odb'])
    for key in ('job_id', 'attempt_id', 'instance', 'step'):
        _text(export[key], key)
    for key in ('basis_sha256', 'adapter_sha256'):
        if not isinstance(export[key], str) or not re.fullmatch('[0-9a-f]{64}', export[key]):
            raise InputError(f'invalid {key}')
    runtime = _keys(export['runtime'], ('abaqus', 'python', 'numpy'))
    for key, version in runtime.items():
        _text(version, f'{key} runtime version')
    _enum(export['coord_precision'], ('single', 'double', 'mixed'), 'coordinate precision')
    if (type(export['step_duration']) not in (float, int) or export['step_duration'] != 1.0 or
            type(export['frame_index']) is not int or export['frame_index'] < 0 or
            type(export['frame_value']) not in (float, int) or
            not math.isfinite(export['frame_value']) or export['frame_value'] < 0.99):
        raise InputError('invalid export step duration or saved frame')
    _read_npz(meta, path.parent, [], keep_facets=False)
    return meta


def _bundle(case_id, base, ref, checked):
    target = f'{case_id}.geometry'
    geometry_hash = None
    try:
        path = _referent(base, ref)
        raw = path.read_bytes()
        geometry_hash = hashlib.sha256(raw).hexdigest()
        checked.append(_identity(path, geometry_hash))
        if geometry_hash != ref['sha256']:
            raise InputError('geometry manifest digest mismatch')
        meta = _geometry_manifest(raw)
        if meta['source_kind'] != 'synthetic' or meta['case_id'] != case_id or meta['export'] is not None:
            raise InputError('case ID, source kind or export disagrees with synthetic study')
        target = f'{case_id}.geometry.npz'
        facets = _read_npz(meta, path.parent, checked)
    except FileNotFoundError:
        return 'missing', [_note('source_missing', target)], geometry_hash, None
    except OSError:
        return 'unverified', [_note('source_unreadable', target)], geometry_hash, None
    except (ValueError, TypeError, KeyError, OverflowError, EOFError, RuntimeError,
            zipfile.BadZipFile, zlib.error, lzma.LZMAError):
        return 'unverified', [_note('export_invalid', target)], geometry_hash, None
    return 'qualified', [], geometry_hash, facets


def _odb_intake(base, case_id, value):
    target = f'{case_id}.odb_path'
    if value is None:
        return 'missing', [_note('source_missing', target)]
    try:
        if not stat.S_ISREG((base / value).stat().st_mode):
            return 'unverified', [_note('source_unreadable', target)]
    except FileNotFoundError:
        return 'missing', [_note('source_missing', target)]
    except OSError:
        return 'unverified', [_note('source_unreadable', target)]
    return 'unverified', [_note('phase_c_not_available', target)]


def _inventory(study, path, identity, *, keep_surfaces):
    synthetic = study['source_kind'] == 'synthetic'
    cases, surfaces = {}, {}
    for case_id, case in study['cases'].items():
        checked, digest = [], None
        if synthetic:
            status, notes, digest, facets = _bundle(case_id, path.parent, case['geometry'], checked)
            if keep_surfaces and facets is not None:
                surfaces[case_id] = facets
        else:
            status, notes = _odb_intake(path.parent, case_id, case['odb_path'])
        cases[case_id] = {'role': case['role'], 'status': status, 'reasons': notes,
                          'checked_files': checked, 'diagnostics': None, 'geometry_sha256': digest}
    contrasts = {}
    for pair_id, pair in study['contrasts'].items():
        notes = [_note('case_unqualified', member) for member in pair.values()
                 if cases[member]['status'] != 'qualified']
        contrasts[pair_id] = {**pair, 'status': 'blocked' if notes else 'eligible',
                               'reasons': _sorted_notes(notes)}
    screen_notes = [_note('match_unresolved', pair_id) for pair_id in contrasts]
    return {'schema_version': 1, 'source_kind': study['source_kind'], 'study': identity,
            'cases': cases, 'contrasts': contrasts, 'screen': {
                'status': 'not_applicable' if synthetic else 'blocked',
                'reasons': [] if synthetic else _sorted_notes(screen_notes)}}, surfaces


def _unit_ray(ray):
    ray = np.asarray(ray, dtype=np.float64)
    if ray.shape != (3,) or not np.isfinite(ray).all():
        raise InputError('expected a finite unit ray in sail axes')
    length = math.hypot(*ray)
    if abs(length - 1.) > 1e-12:
        raise InputError('expected a finite unit ray in sail axes')
    return ray / length


def _rays():
    alpha = math.radians(5)
    sin, cos = math.sin(alpha), math.cos(alpha)
    diagonal = -sin / math.sqrt(2.)
    return {key: _unit_ray(ray) for key, ray in (
        ('normal', (0., 0., -1.)), ('alpha5_phi0', (-sin, 0., -cos)),
        ('alpha5_phi45', (diagonal, diagonal, -cos)))}


def _placements():
    rotations = np.array([[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
                          [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
                          [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]], dtype=np.float64)
    return tuple((r, r @ [-10., -10., 0.]) for r in rotations)


def _integrate(xyz, signs, ray, optics, origin):
    cross = _crosses(xyz)
    twice_area = np.hypot.reduce(cross, axis=1)
    normal = (cross / twice_area[:, None]) * signs[:, None]
    c = np.maximum(0., -(normal @ ray))
    a0, r0 = optics
    forces = (0.5 * twice_area * c)[:, None] * (a0 * ray - (2. * r0 * c)[:, None] * normal)
    arm = (xyz / 3.).sum(axis=1) - origin
    return np.stack((forces.sum(axis=0), np.cross(arm, forces).sum(axis=0)))


def _route(facets, ray, optics, placements, *, direct, origin=_ORIGIN, local_origin=_LOCAL_ORIGIN):
    loads = []
    with np.errstate(over='ignore', invalid='ignore', divide='ignore', under='ignore'):
        for r, t in placements:
            if direct:
                loads.append(_integrate(facets.xyz @ r.T + t, facets.signs, ray, optics, origin))
            else:
                f, m = _integrate(facets.xyz, facets.signs, r.T @ ray, optics, local_origin)
                force = r @ f
                loads.append(np.stack((force, r @ m + np.cross(r @ local_origin + t - origin, force))))
    return np.stack(loads)


def _allowances(area=400., length=20.):
    return np.array([[1e-12 * area], [1e-12 * area * length]])


def _blocks(f, m, ray):
    along = float(f @ ray)
    return {'force': f.tolist(), 'moment': m.tolist(), 'force_norm': math.hypot(*f),
            'moment_norm': math.hypot(*m), 'ray_projection': along,
            'transverse_force': (f - along * ray).tolist()}


def _finite_fields(block):
    return [key for key, value in block.items() if not np.isfinite(value).all()]


def _scaled_load(normalized, pressure):
    block = {key: ([x * pressure for x in value] if isinstance(value, list) else value * pressure)
             for key, value in normalized.items() if key != 'panels'}
    panels = {key: {name: [x * pressure for x in value] for name, value in panel.items()}
              for key, panel in normalized['panels'].items()}
    bad = _finite_fields(block) + [f'panels.{key}.{field}' for key, panel in panels.items()
                                  for field in _finite_fields(panel)]
    return (None, [_note('physical_overflow', f'physical.{key}') for key in bad]) if bad else (
        {**block, 'panels': panels}, [])


def _load(facets, ray, optics, placements, pressure):
    row = {'validation': 'not_attempted', 'normalized': None, 'physical': None, 'reasons': []}
    if facets is None:
        row['reasons'] = [_note('case_unavailable', 'load')]
        return row
    with np.errstate(over='ignore', invalid='ignore', under='ignore'):
        local = _route(facets, ray, optics, placements, direct=False)
        direct = _route(facets, ray, optics, placements, direct=True)
        actual, expected = local.sum(axis=0), direct.sum(axis=0)
        diff = np.abs(actual - expected)
        allow = _allowances() + 1e-10 * np.maximum(np.abs(actual), np.abs(expected))
        notes = []
        for name, vectors in (('local', local), ('direct', direct)):
            for panel, kind, axis in np.argwhere(~np.isfinite(vectors)):
                notes.append(_note('nonfinite', f'{name}.panels.{panel}.{("force", "moment")[kind]}.{axis}'))
        for kind, axis in np.ndindex((2, 3)):
            field = f'{("force", "moment")[kind]}.{axis}'
            if not math.isfinite(diff[kind, axis]) or not math.isfinite(allow[kind, axis]):
                notes.append(_note('nonfinite', f'normalized.{field}'))
            elif diff[kind, axis] > allow[kind, axis]:
                notes.append(_note('route_mismatch', f'normalized.{field}'))
        normalized = _blocks(*actual, ray) if not notes else None
        if normalized is not None:
            notes.extend(_note('nonfinite', f'normalized.{field}') for field in _finite_fields(normalized))
    if notes:
        row.update(validation='failed', reasons=_sorted_notes(notes))
        return row
    normalized['panels'] = {str(i): {'force': f.tolist(), 'moment': m.tolist()}
                            for i, (f, m) in enumerate(local)}
    if pressure is None:
        physical, notes = None, [_note('pressure_omitted', 'physical')]
    else:
        physical, notes = _scaled_load(normalized, pressure)
    return {'validation': 'passed', 'normalized': normalized, 'physical': physical,
            'reasons': _sorted_notes(notes)}


def _delta(case, ref):
    return {key: ([x - y for x, y in zip(case[key], ref[key])] if isinstance(case[key], list)
                  else case[key] - ref[key]) for key in _DELTA_FIELDS}


def _angle(case, ref, floor):
    if not (case['force_norm'] > floor and ref['force_norm'] > floor):
        return None
    f = np.asarray(case['force']) / case['force_norm']
    g = np.asarray(ref['force']) / ref['force_norm']
    return math.degrees(math.atan2(math.hypot(*np.cross(f, g)), float(f @ g)))


def _comparison(case, ref, ref_id, pressure, *, eligible=True, flat=False):
    row = {'reference_id': ref_id, 'normalized': None, 'physical': None, 'reasons': []}
    if flat:
        row.update(angle_normalized_deg=None, angle_physical_deg=None)

    def withhold(code, *fields):
        row['reasons'].extend(_note(code, field) for field in fields)

    missing = [(value, code) for value, code in ((case, 'case_unavailable'),
                                                (ref, 'reference_unavailable'))
               if value['validation'] != 'passed']
    if missing or not eligible:
        for _, code in missing or [(None, 'ineligible')]:
            withhold(code, 'normalized', 'physical')
            if flat:
                withhold(code, 'angle_normalized_deg', 'angle_physical_deg')
        row['reasons'] = _sorted_notes(row['reasons'])
        return row
    for key in ('normalized', 'physical'):
        x, y = case[key], ref[key]
        if x is None or y is None:
            withhold('pressure_omitted' if pressure is None else 'physical_overflow', key)
            continue
        delta = _delta(x, y)
        bad = _finite_fields(delta)
        if bad:
            withhold('nonfinite' if key == 'normalized' else 'physical_overflow',
                     *(f'{key}.{field}' for field in bad))
        else:
            row[key] = delta
    if flat:
        floor = float(_allowances()[0, 0])
        angle = _angle(case['normalized'], ref['normalized'], floor)
        row['angle_normalized_deg'] = angle
        if angle is None:
            withhold('force_at_or_below_floor', 'angle_normalized_deg')
        if pressure is None:
            code = 'pressure_omitted'
        elif pressure == 0:
            code = 'zero_pressure'
        elif case['physical'] is None or ref['physical'] is None:
            code = 'physical_overflow'
        elif (angle is None or case['physical']['force_norm'] <= floor * pressure or
              ref['physical']['force_norm'] <= floor * pressure):
            code = 'force_at_or_below_floor'
        else:
            code = None
            row['angle_physical_deg'] = angle
        if code is not None:
            withhold(code, 'angle_physical_deg')
    row['reasons'] = _sorted_notes(row['reasons'])
    return row


def _evaluate(study, evidence, surfaces, evidence_identity):
    synthetic = study['source_kind'] == 'synthetic'
    optics, pressure = _OPTICS[study['optics']], study.get('pressure_pa')
    rays, placements = _rays(), _placements()
    flat = Facets(_FLAT, np.array([1], dtype=np.int8)) if synthetic else None
    loads = {key: {attitude: _load(facets, ray, optics, placements, pressure)
                   for attitude, ray in rays.items()}
             for key, facets in [('flat', flat)] + [(key, surfaces.get(key)) for key in study['cases']]}
    controls = {pair['case_id']: pair for pair in evidence['contrasts'].values()}
    comparisons = {}
    for key in study['cases']:
        rows = {}
        for attitude, load in loads[key].items():
            entries = {'flat': _comparison(load, loads['flat'][attitude], 'flat', pressure, flat=True)}
            if key in controls:
                pair = controls[key]
                ref = pair['control_id']
                entries['control'] = _comparison(load, loads[ref][attitude], ref, pressure,
                                                  eligible=pair['status'] == 'eligible')
            rows[attitude] = entries
        comparisons[key] = rows
    completion = {'contrasts': {pair_id: {
        'status': 'not_applicable' if synthetic else 'not_complete',
        'reasons': [] if synthetic else [_note('ineligible', pair_id)]}
        for pair_id in evidence['contrasts']},
        'screen': {'status': 'not_applicable' if synthetic else 'not_complete',
                   'reasons': [] if synthetic else [_note('ineligible', 'screen')]}}
    return {'schema_version': 1, 'source_kind': study['source_kind'], 'study': evidence['study'],
            'evidence': evidence_identity, 'evaluator_version': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'layout_id': _LAYOUT, 'transforms': [{'R': r.tolist(), 't_m': t.tolist()} for r, t in placements],
            'moment_origin_m': _ORIGIN.tolist(), 'replication': 'one_rotated_panel_shape',
            'faceting_method': 'piecewise_planar_triangles', 'faceting_version': 1,
            'A0': optics[0], 'R0': optics[1], 'pressure_pa': pressure,
            'units': {'normalized_force': 'm2', 'normalized_moment': 'm3', 'physical_force': 'N',
                      'physical_moment': 'N_m', 'angle': 'deg'},
            'rays': {key: ray.tolist() for key, ray in rays.items()}, 'loads': loads,
            'comparisons': comparisons, 'completion': completion}


def _encoded(obj):
    return (json.dumps(obj, ensure_ascii=True, allow_nan=False, indent=2, sort_keys=True) + '\n').encode('utf-8')


def _output_path(path):
    out = path.resolve()
    if any(out.is_relative_to(root.resolve()) for root in _SOLVER_ROOTS):
        raise InputError('choose an output directory outside solver results')
    return out


def _publish(out, items):
    out.mkdir()  # Exclusive reservation; never reuse an earlier output directory.
    created = []
    try:
        for name, content in items.items():
            path = out / name
            with path.open('xb') as stream:
                created.append(path)
                stream.write(content)
    except OSError:
        for path in created:
            path.unlink()
        out.rmdir()
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('inventory', 'evaluate'))
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--out-dir', required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        path, out = args.manifest.resolve(), _output_path(args.out_dir)
        raw = path.read_bytes()
        study = _study(raw)
        evidence, surfaces = _inventory(study, path, _identity(path, hashlib.sha256(raw).hexdigest()),
                                        keep_surfaces=args.command == 'evaluate')
        encoded = _encoded(evidence)
        items = {'evidence.json': encoded}
        if args.command == 'evaluate':
            evidence_id = _identity(out / 'evidence.json', hashlib.sha256(encoded).hexdigest())
            items['optical_results.json'] = _encoded(_evaluate(study, evidence, surfaces, evidence_id))
        _publish(out, items)
    except (OSError, ValueError, OverflowError, RuntimeError) as error:
        print(f'sail_optics: {error}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
