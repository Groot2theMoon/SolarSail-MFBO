"""Portable whole-sail optics, with reviewed real-export evidence and matching."""

import argparse
from dataclasses import dataclass
from decimal import Decimal, DecimalException, localcontext
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


def _initial_candidate(case, *, linked=False):
    # Preserve the adapter's initial-manifest envelope. Only portable C accepts links.
    if not linked and case['geometry'] is not None:
        raise InputError('initial ODB candidates require null geometry')
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


def _study(raw, *, allow_linked=False):
    data = _keys(_decode(raw), ('schema_version', 'source_kind', 'layout_id', 'optics',
                                'cases', 'contrasts'), ('pressure_pa', 'parent'))
    _version(data)
    kind = _enum(data['source_kind'], ('synthetic', 'odb'), 'source kind')
    linked = 'parent' in data
    if linked and (not allow_linked or kind != 'odb'):
        raise InputError('only a portable linked ODB study may specify a parent')
    if linked:
        _file_ref(data['parent'])
    _enum(data['optics'], _OPTICS if kind == 'synthetic' else ('default',), 'optics')
    _enum(data['layout_id'], (_LAYOUT,), 'layout')
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
            _initial_candidate(case, linked=linked)
    paired = set()
    for contrast_id, pair in contrasts.items():
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


def _selection(study, base):
    """Compare the frozen selection; rebasing file paths does not change it."""
    def file_ref(ref):
        return {**ref, 'path': str(_referent(base, ref))}

    comparable = {k: v for k, v in study.items() if k != 'parent'}
    cases = {}
    for case_id, item in study['cases'].items():
        case = {k: v for k, v in item.items() if k != 'geometry'}
        if case['odb_path'] is not None:
            case['odb_path'] = str((base / case['odb_path']).resolve())
        case['records'] = {k: file_ref(v) for k, v in case['records'].items()}
        producer = case['producer']
        if producer is not None and isinstance(producer['value'], dict):
            current = producer['value'].get('current_file')
            if (isinstance(current, dict) and isinstance(current.get('path'), str) and
                    current['path'] and '\0' not in current['path']):
                case['producer'] = {**producer, 'value': {**producer['value'], 'current_file': file_ref(current)}}
        cases[case_id] = case
    return _encoded({**comparable, 'cases': cases})


def _check_parent(study, path):
    if 'parent' not in study:
        return
    ref = study['parent']
    parent = _referent(path.parent, ref)
    raw = parent.read_bytes()
    if hashlib.sha256(raw).hexdigest() != ref['sha256']:
        raise InputError('linked study initial parent digest mismatch')
    initial = _study(raw)
    if initial['source_kind'] != 'odb' or _selection(study, path.parent) != _selection(initial, parent.parent):
        raise InputError('linked study changed the frozen initial selection or facts')


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


class _EvidenceError(InputError):
    def __init__(self, code, target):
        super().__init__(f'{code}: {target}')
        self.note = _note(code, target)


class _CaseFiles:
    def __init__(self, case_id, base, records, cache):
        self.case_id, self.base, self.records, self.cache = case_id, base, records, cache
        self.checked = {}

    def digest(self, path):
        if not stat.S_ISREG(path.stat().st_mode):
            raise InputError('expected an ordinary file')
        if path not in self.cache:
            with path.open('rb') as stream:
                self.cache[path] = _stream_sha256(stream)
        value = self.cache[path]
        self.checked[path] = _identity(path, value)
        return value

    def read(self, ref, *, base=None, code='record_conflict', target):
        _file_ref(ref)
        path = _referent(self.base if base is None else base, ref)
        try:
            actual = self.digest(path)
        except (OSError, InputError):
            raise _EvidenceError('record_unavailable', target) from None
        if actual != ref['sha256']:
            raise _EvidenceError(code, target)
        return path

    def refs(self, refs, field):
        if not isinstance(refs, list) or not refs:
            raise _EvidenceError('record_unavailable', f'{self.case_id}.{field}')
        paths = {}
        for ref in refs:
            _keys(ref, ('record_id', 'locator'))
            key = _text(ref['record_id'], 'record ID')
            target = f'{self.case_id}.{field}.{key}'
            if key not in self.records or not _text(ref['locator'], 'source locator').strip():
                raise _EvidenceError('record_unavailable', target)
            paths[key] = self.read(self.records[key], target=target)
        return paths


def _finite(value):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise InputError('expected a finite JSON number')
    return value


def _decimal(value):
    if not isinstance(value, str) or not re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?', value, re.ASCII):
        raise InputError('expected an exact decimal string')
    return Decimal(value)


def _sha(value):
    _file_ref({'path': '.', 'sha256': value})
    return value


def _list(value):
    if not isinstance(value, list):
        raise InputError('expected a list')
    return value


def _labels(values):
    rows = _list(values)
    if any(type(x) is not int for x in rows) or len(set(rows)) != len(rows):
        raise InputError('node labels must be unique integers')
    return frozenset(rows)


def _typed(value):
    """Normalize typed property leaves; retain order only for arrays."""
    if isinstance(value, dict):
        if set(value) == {'value', 'unit'}:
            return ('quantity', _decimal(value['value']), _text(value['unit'], 'unit'))
        return ('map', tuple(sorted((k, _typed(v)) for k, v in value.items())))
    if isinstance(value, list):
        return ('list', tuple(map(_typed, value)))
    if type(value) in (bool, int, str):
        return (type(value).__name__, value)
    raise InputError('untyped property')


def _run_facts(run, paths, case):
    _keys(run, ('job_id', 'attempt_id', 'driver_revision', 'modified_source', 'submitted_deck',
                'completed', 'instance', 'terminal_step', 'step_duration', 'frame_index'))
    for key in ('job_id', 'attempt_id', 'driver_revision', 'submitted_deck', 'instance', 'terminal_step'):
        _text(run[key], key)
    index, modified = run['frame_index'], run['modified_source']
    if (run['submitted_deck'] not in paths or len(set(paths.values())) < 2 or
            (modified is not None and _text(modified, 'modified_source') not in paths) or
            run['attempt_id'] != case['attempt_id'] or run['completed'] is not True or
            _finite(run['step_duration']) != 1 or (index is not None and (type(index) is not int or index < 0))):
        raise InputError('invalid completed run, source binding or terminal selection')
    return run


def _producer_facts(item, paths, run, odb, files):
    _keys(item, ('job_id', 'attempt_id', 'product_id', 'product_sha256', 'copies', 'current_id', 'current_file'))
    for key in ('job_id', 'attempt_id', 'product_id', 'current_id'):
        _text(item[key], key)
    file = item['current_file']
    _file_ref(file)
    current = _referent(files.base, file)
    if run is None or not set(paths.values()) - {current, files.deck}:
        raise InputError('no run-linked producer record beyond the current ODB and deck')
    if (file['sha256'] != odb['sha256'] or str(current) != odb['path'] or
            any(item[k] != run[k] for k in ('job_id', 'attempt_id'))):
        raise _EvidenceError('record_conflict', f'{files.case_id}.producer')
    copies = _list(item['copies'])
    digest = item['product_sha256']
    if digest is not None and _sha(digest) != odb['sha256']:
        raise _EvidenceError('record_conflict', f'{files.case_id}.producer.product_sha256')
    if not copies and digest is None:
        raise InputError('direct product digest missing')
    node = item['product_id']
    seen = {node}
    for row in copies:
        _keys(row, ('from_id', 'to_id', 'sha256', 'refs'))
        to = _text(row['to_id'], 'copy ID')
        if row['from_id'] != node or to in seen or (row['sha256'] is not None and _sha(row['sha256']) != file['sha256']):
            raise InputError('discontinuous or conflicting producer copy chain')
        if not set(files.refs(row['refs'], f'producer.copies.{to}').values()) - {current}:
            raise InputError('copy cites only the current ODB')
        node = to
        seen.add(to)
    if node != item['current_id']:
        raise InputError('copy chain does not end at current ID')
    return item


def _geometry_basis(basis):
    def matrix(value, shape):
        try:
            def check(x):
                return [check(v) for v in x] if isinstance(x, list) else _finite(x)
            array = np.asarray(check(value), dtype=np.float64)
        except (TypeError, ValueError, OverflowError) as error:
            raise InputError('invalid geometry mapping') from error
        if array.shape != shape:
            raise InputError('invalid geometry mapping shape')
        return array

    _keys(basis, ('units', 'R', 't_m', 'anchors', 'sides'))
    _enum(basis['units'], ('m',), 'reference units')
    r, t = matrix(basis['R'], (3, 3)), matrix(basis['t_m'], (3,))
    with np.errstate(over='ignore', invalid='ignore'):
        if (not np.allclose(r.T @ r, np.eye(3), rtol=0, atol=1e-12) or
                not math.isclose(np.linalg.det(r), 1., abs_tol=1e-12, rel_tol=0)):
            raise InputError('mapping is not a proper rotation')
    rows = _list(basis['anchors'])
    if len(rows) != 3:
        raise InputError('expected three ordered reference anchors')
    labels = []
    for anchor, expected in zip(rows, _FLAT[0]):
        _keys(anchor, ('node_label', 'odb_xyz_m', 'panel_xyz_m'))
        labels.append(anchor['node_label'])
        ref = matrix(anchor['odb_xyz_m'], (3,))
        if (not np.array_equal(matrix(anchor['panel_xyz_m'], (3,)), expected) or
                not np.allclose(r @ ref + t, expected, atol=1e-9, rtol=0)):
            raise InputError('reference anchor does not map to its required panel position')
    signs = {}
    if not isinstance(basis['sides'], dict) or not basis['sides']:
        raise InputError('missing element side evidence')
    for label, side in basis['sides'].items():
        if not re.fullmatch(r'(?:0|-?[1-9][0-9]*)', label):
            raise InputError('invalid element label in side evidence')
        signs[int(label)] = 1 if _enum(side, ('SPOS', 'SNEG'), 'face side') == 'SPOS' else -1
    kind = 'identity' if np.array_equal(r, np.eye(3)) and not t.any() else 'rigid'
    raw = json.dumps(basis, ensure_ascii=True, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('utf-8')
    return _labels(labels), signs, kind, hashlib.sha256(raw).hexdigest()


def _clamp(value, role):
    if role == 'control':
        return _enum(value, ('absent',), 'control clamp evidence')
    _keys(value, ('x_c', 'd_c', 'corner_initial_pa', 'clamp_initial_pa', 'hardware', 'connections'))
    x, d, corner, stress = (_decimal(value[k]) for k in ('x_c', 'd_c', 'corner_initial_pa', 'clamp_initial_pa'))
    with localcontext() as ctx:
        ctx.prec = max(1, len(d.as_tuple().digits) + len(corner.as_tuple().digits))
        if corner <= 0 or stress != d * corner:
            raise InputError('clamp preload does not match the applied corner-to-clamp ratio')
    for key in ('hardware', 'connections'):
        sides = _keys(value[key], ('left', 'right'))
        if any(not isinstance(v, dict) or not v for v in sides.values()):
            raise InputError('attachment evidence is incomplete')
    rp_names = frozenset(_text(v['reference_point'], 'clamp reference point')
                         for v in value['connections'].values() if 'reference_point' in v)
    return {'x_c': x, 'd_c': d, 'corner': corner, 'stress': stress, 'hardware': _typed(value['hardware']),
            'connections': _typed(value['connections']), 'rp_names': rp_names}


def _membrane(value):
    _keys(value, ('reference_mesh_sha256', 'node_labels', 'material', 'thickness_m'))
    labels = _labels(value['node_labels'])
    thickness = _decimal(value['thickness_m'])
    if not labels or thickness <= 0 or not isinstance(value['material'], dict) or not value['material']:
        raise InputError('incomplete membrane evidence')
    return _sha(value['reference_mesh_sha256']), labels, _typed(value['material']), thickness


def _load_facts(value):
    _keys(value, ('mode', 'corner_initial_pa', 'shared_steps'))
    _enum(value['mode'], ('cable_stress',), 'preload mode')
    corners = _keys(value['corner_initial_pa'], ('apex', 'first_base', 'second_base'))
    steps, names = [], set()
    for row in _list(value['shared_steps']):
        _keys(row, ('name', 'controls'))
        name = _text(row['name'], 'load step')
        if name in names or not isinstance(row['controls'], dict):
            raise InputError('invalid shared step')
        names.add(name)
        steps.append((name, _typed(row['controls'])))
    if not steps:
        raise InputError('shared steps are missing')
    return tuple(_decimal(corners[k]) for k in ('apex', 'first_base', 'second_base')), tuple(steps)


def _solver(value):
    _keys(value, ('shared_options', 'applied_stabilization_by_step'))
    if any(not isinstance(v, dict) for v in value.values()):
        raise InputError('expected solver property maps')
    return _typed(value)


def _supports(value, role, clamp):
    _keys(value, ('common', 'attachment', 'replaced_common', 'attachment_patch'))
    patch = _labels(value['attachment_patch'])
    sets = {}
    for key in ('common', 'attachment', 'replaced_common'):
        rows = set()
        for row in _list(value[key]):
            _keys(row, ('target_kind', 'target', 'step', 'dof', 'law', 'parameters'))
            kind = _enum(row['target_kind'], ('membrane_node', 'reference_point'), 'support target kind')
            target = row['target']
            valid = type(target) is int if kind == 'membrane_node' else bool(_text(target, 'reference point'))
            dof = row['dof']
            if not valid or type(dof) is not int or not 1 <= dof <= 6 or not isinstance(row['parameters'], dict):
                raise InputError('invalid support row')
            support = (kind, target, _text(row['step'], 'step'), dof, _text(row['law'], 'law'), _typed(row['parameters']))
            if support in rows:
                raise InputError('duplicate support row')
            rows.add(support)
        sets[key] = rows
    if role == 'control' and (patch or sets['attachment'] or sets['replaced_common']):
        raise InputError('control has attachment-linked restraints')
    if sets['common'] & sets['replaced_common'] or any(
            k != 'membrane_node' or n not in patch for k, n, *_ in sets['replaced_common']):
        raise InputError('common replacements are not confined to the declared attachment patch')
    if any(not (n in patch if k == 'membrane_node' else clamp is not None and n in clamp['rp_names'])
           for k, n, *_ in sets['attachment']):
        raise InputError('attachment restraints are not on the declared clamp patch or reference points')
    return {**sets, 'patch': patch}


def _perturbation(value, mesh):
    _keys(value, ('reference_mesh_sha256', 'frame', 'node_map', 'delta_m', 'precision'))
    _enum(value['frame'], ('panel_local_m',), 'applied-field frame')
    precision = _keys(value['precision'], ('encoding', 'basis'))
    _enum(precision['encoding'], ('exact_decimal', 'binary64_exact'), 'applied-field precision')
    _text(precision['basis'], 'applied-field precision evidence')
    if mesh is None or _sha(value['reference_mesh_sha256']) != mesh[0]:
        raise InputError('applied field has no matching reviewed reference mesh')
    mapping = _list(value['node_map'])
    if any(not isinstance(row, list) or len(row) != 2 for row in mapping):
        raise InputError('invalid applied-field node map')
    _labels([row[0] for row in mapping])
    if _labels([row[1] for row in mapping]) != mesh[1]:
        raise InputError('applied-field mapping is incomplete')
    offsets = _list(value['delta_m'])
    if any(not isinstance(row, list) or len(row) != 4 for row in offsets) or _labels([row[0] for row in offsets]) != mesh[1]:
        raise InputError('applied-field offsets are incomplete')
    return {row[0]: tuple(_decimal(v) for v in row[1:]) for row in offsets}


def _export_data(meta):
    exp = _keys(meta['export'], ('odb', 'job_id', 'attempt_id', 'instance', 'step', 'step_duration',
                                 'frame_index', 'frame_value', 'coord_precision', 'basis_sha256',
                                 'runtime', 'adapter_sha256'))
    _file_ref(exp['odb'])
    for key in ('job_id', 'attempt_id', 'instance', 'step'):
        _text(exp[key], key)
    _sha(exp['adapter_sha256'])
    _sha(exp['basis_sha256'])
    for version in _keys(exp['runtime'], ('abaqus', 'python', 'numpy')).values():
        _text(version, 'runtime version')
    _enum(exp['coord_precision'], ('single', 'double', 'mixed'), 'coordinate precision')
    if (_finite(exp['step_duration']) != 1 or _finite(exp['frame_value']) < 0.99 or
            type(exp['frame_index']) is not int or exp['frame_index'] < 0):
        raise InputError('export does not identify a covered terminal frame')
    return exp


def _round_trip(paths, mapping, meta, exp, bundle_ref, files, face_paths):
    target = f'{files.case_id}.export_verification'
    receipts = {p for p in paths.values() if p.name == 'export_verification.json'}
    if len(receipts) != 1:
        raise _EvidenceError('record_unavailable', target)
    path, = receipts
    rec = _keys(_decode(path.read_bytes()),
                ('schema_version', 'geometry', 'odb', 'adapter_sha256', 'runtime', 'instance', 'step',
                 'step_duration', 'frame_index', 'frame_value', 'node_count', 'element_counts',
                 'coord_precision', 'mapping_kind', 'checks', 'unexercised', 'status'))
    _version(rec)
    _file_ref(rec['odb'])
    scope = ('adapter_sha256', 'runtime', 'coord_precision')
    _mesh_counts(rec['node_count'], rec['element_counts'])
    if (rec['status'] != 'passed' or rec['mapping_kind'] != mapping or
            any(rec[k] != exp[k] for k in scope) or not meta['element_counts'].keys() <= rec['element_counts'].keys()):
        raise _EvidenceError('record_conflict', target)
    for item in _list(rec['unexercised']):
        _text(item, 'unexercised path')
    checks = _keys(rec['checks'], ('terminal_frame', 'anchors', 'coordinates', 'topology', 'face_signs',
                                 'portable_readback', 'odb_unchanged', 'negative_inputs'))
    for name, check in checks.items():
        _keys(check, ('status', 'evidence'))
        if check['status'] != 'passed' or not _list(check['evidence']):
            raise _EvidenceError('record_conflict', f'{target}.{name}')
        for row in check['evidence']:
            _keys(row, ('file', 'locator'))
            if not _text(row['locator'], 'verification locator').strip():
                raise _EvidenceError('record_unavailable', f'{target}.{name}')
            files.read(row['file'], base=path.parent, target=f'{target}.{name}')
    verified = files.read(rec['geometry'], base=path.parent, target=target)
    if rec['geometry']['sha256'] == bundle_ref['sha256']:
        keys = ('instance', 'step', 'step_duration', 'frame_index', 'frame_value')
        if (any(rec[k] != exp[k] for k in keys) or rec['odb']['sha256'] != exp['odb']['sha256'] or
                _referent(path.parent, rec['odb']) != Path(files.odb['path']) or
                any(rec[k] != meta[k] for k in ('node_count', 'element_counts'))):
            raise _EvidenceError('record_conflict', target)
    else:
        # Check coverage from the previously verified geometry without opening its ODB.
        sample = _geometry_manifest(verified.read_bytes())
        sample_exp = _export_data(sample)
        keys = (*scope, 'instance', 'step', 'step_duration', 'frame_index', 'frame_value')
        if (sample['source_kind'] != 'odb' or any(rec[k] != sample_exp[k] for k in keys) or
                any(rec[k] != sample[k] for k in ('node_count', 'element_counts')) or
                rec['odb']['sha256'] != sample_exp['odb']['sha256'] or
                _referent(path.parent, rec['odb']) != _referent(verified.parent, sample_exp['odb'])):
            raise _EvidenceError('record_conflict', target)
        sample_npz = files.read(sample['npz'], base=verified.parent, target=target)
        with np.load(sample_npz, allow_pickle=False) as arrays:
            for kind, signs in face_paths.items():
                tested = _array(arrays, f'face_sign_{kind}', 'int8', (rec['element_counts'][kind],))
                if not signs <= set(tested):
                    raise _EvidenceError('record_conflict', f'{target}.face_signs.{kind}')


_BUNDLE_ERRORS = (ValueError, TypeError, KeyError, OverflowError, EOFError, RuntimeError,
                  zipfile.BadZipFile, zlib.error, lzma.LZMAError)


def _real_bundle(ref, run, basis, files, *, keep_facets):
    target = f'{files.case_id}.geometry'
    _file_ref(ref)
    path = _referent(files.base, ref)
    try:
        if not stat.S_ISREG(path.stat().st_mode):
            raise InputError('expected an ordinary geometry manifest')
        raw = path.read_bytes()
    except FileNotFoundError:
        raise _EvidenceError('export_missing', target) from None
    digest = hashlib.sha256(raw).hexdigest()
    files.checked[path] = _identity(path, digest)
    if digest != ref['sha256']:
        raise _EvidenceError('export_invalid', target)
    meta = _geometry_manifest(raw)
    exp = _export_data(meta)
    if meta['source_kind'] != 'odb' or meta['case_id'] != files.case_id:
        raise InputError('case ID or source kind conflicts with real manifest')
    if (_referent(path.parent, exp['odb']) != Path(files.odb['path']) or
            exp['odb']['sha256'] != files.odb['sha256'] or
            any(exp[k] != run[k] for k in ('job_id', 'attempt_id', 'instance', 'step_duration')) or
            exp['step'] != run['terminal_step'] or
            (run['frame_index'] is not None and exp['frame_index'] != run['frame_index']) or
            exp['basis_sha256'] != basis[3]):
        raise _EvidenceError('record_conflict', f'{target}.export')
    npz_ref = meta['npz']
    _file_ref(npz_ref)
    target += '.npz'
    npz = _referent(path.parent, npz_ref)
    face_paths = {}
    try:
        if not stat.S_ISREG(npz.stat().st_mode):
            raise InputError('expected a regular NPZ')
        with npz.open('rb') as stream:
            actual = _stream_sha256(stream)
            files.checked[npz] = _identity(npz, actual)
            if actual != npz_ref['sha256']:
                raise InputError('NPZ digest mismatch')
            stream.seek(0)
            with zipfile.ZipFile(stream) as archive:
                names = archive.namelist()
                if len(names) != len(set(names)) or any(not name.endswith('.npy') for name in names):
                    raise InputError('invalid NPZ members')
            stream.seek(0)
            with np.load(stream, allow_pickle=False) as bundle:
                arrays = {k: bundle[k] for k in bundle.files}
        counts = meta['element_counts']
        facets = _geometry(arrays, meta['node_count'], counts, keep_facets=keep_facets)
        labels = set(arrays['node_labels'].tolist())
        if not basis[0] <= labels:
            raise _EvidenceError('record_conflict', f'{files.case_id}.geometry_basis.anchors')
        if sum(counts.values()) != len(basis[1]):
            raise _EvidenceError('record_conflict', f'{files.case_id}.geometry_basis.sides')
        for kind in counts:
            signs = arrays[f'face_sign_{kind}']
            ids = arrays[f'element_labels_{kind}']
            if any(basis[1].get(label) != sign for label, sign in zip(ids.tolist(), signs.tolist())):
                raise _EvidenceError('record_conflict', f'{files.case_id}.geometry_basis.sides')
            face_paths[kind] = set(signs)
    except _EvidenceError:
        raise
    except FileNotFoundError:
        raise _EvidenceError('export_missing', target) from None
    except (OSError, *_BUNDLE_ERRORS):
        raise _EvidenceError('export_invalid', target) from None
    return digest, facets, meta, exp, face_paths, labels


def _real_case(base, case_id, case, cache, *, keep_surfaces):
    files = _CaseFiles(case_id, base, case['records'], cache)
    notes, settings = [], {}

    def fact(item, field, validate, *, required=False, failure='record_conflict'):
        if item is None:
            if required:
                notes.append(_note('producer_unresolved' if field == 'producer' else 'record_unavailable',
                                   f'{case_id}.{field}'))
            return None
        try:
            paths = files.refs(item['refs'], field)
            return validate(item['value'], paths)
        except _EvidenceError as error:
            notes.append(error.note)
        except (OSError, ValueError, OverflowError, TypeError, KeyError, DecimalException):
            notes.append(_note(failure, f'{case_id}.{field}'))
        return None

    diagnostics = case['diagnostics']
    def result(status, digest=None, facets=None):
        return ({'role': case['role'], 'status': status, 'reasons': _sorted_notes(notes),
                 'checked_files': list(files.checked.values()), 'diagnostics': diagnostics,
                 'geometry_sha256': digest}, facets, settings)

    def source_fault(status, code):
        notes.append(_note(code, f'{case_id}.odb_path'))
        return result(status)

    if case['odb_path'] is None:
        return source_fault('missing', 'source_missing')
    try:
        source = (base / case['odb_path']).resolve()
        files.odb = _identity(source, files.digest(source))
    except FileNotFoundError:
        return source_fault('missing', 'source_missing')
    except (OSError, InputError, RuntimeError):
        return source_fault('unverified', 'source_unreadable')

    run = fact(case['run'], 'run', lambda v, p: _run_facts(v, p, case), required=True)
    basis = fact(case['geometry_basis'], 'geometry_basis', lambda v, _: _geometry_basis(v), required=True)
    files.deck = _referent(base, case['records'][run['submitted_deck']]) if run is not None else None
    producer = fact(case['producer'], 'producer', lambda v, p: _producer_facts(v, p, run, files.odb, files),
                    required=True, failure='producer_unresolved')
    raw = case['settings']
    for key, decode in (('membrane', _membrane), ('load', _load_facts), ('solver', _solver),
                        ('clamp', lambda v: _clamp(v, case['role'])),
                        ('supports', lambda v: _supports(v, case['role'], settings['clamp'])),
                        ('perturbation', lambda v: _perturbation(v, settings['membrane']))):
        settings[key] = fact(raw[key], f'settings.{key}', lambda v, _: decode(v), required=key == 'clamp')
    clamp, load = settings['clamp'], settings['load']
    if isinstance(clamp, dict) and load is not None and any(v != clamp['corner'] for v in load[0]):
        notes.append(_note('record_conflict', f'{case_id}.settings.clamp.corner_initial_pa'))
        settings['clamp'] = None
    mesh, supports = settings['membrane'], settings['supports']
    if mesh is not None and supports is not None:
        nodes = supports['patch'] | {node for key in ('common', 'attachment', 'replaced_common')
                                    for kind, node, *_ in supports[key] if kind == 'membrane_node'}
        if not nodes <= mesh[1]:
            notes.append(_note('record_conflict', f'{case_id}.settings.supports.labels'))
            settings['supports'] = None
    observations = []
    for row in diagnostics['observations']:
        try:
            _keys(row, ('step', 'metric', 'value', 'units', 'refs'))
            for key in ('step', 'metric', 'units'):
                _text(row[key], key)
            _finite(row['value'])
            files.refs(row['refs'], 'diagnostics')
            observations.append(row)
        except _EvidenceError as error:
            notes.append(error.note)
        except (OSError, ValueError, OverflowError, TypeError, KeyError):
            notes.append(_note('record_conflict', f'{case_id}.diagnostics'))
    if len(observations) != len(diagnostics['observations']):
        diagnostics = {'status': 'unresolved', 'observations': observations,
                       'reason': 'one or more reviewed observations lack valid current evidence'}
    if case['geometry'] is None:
        notes.append(_note('export_missing', f'{case_id}.geometry'))
    if run is None or producer is None or basis is None or settings['clamp'] is None or case['geometry'] is None:
        return result('unverified')
    digest = None
    try:
        digest, facets, meta, exp, face_paths, labels = _real_bundle(case['geometry'], run, basis, files,
                                                                   keep_facets=keep_surfaces)
        if settings['membrane'] is not None and settings['membrane'][1] != labels:
            raise _EvidenceError('record_conflict', f'{case_id}.settings.membrane.node_labels')
        basis_records = files.refs(case['geometry_basis']['refs'], 'geometry_basis')
        _round_trip(basis_records, basis[2], meta, exp, case['geometry'], files, face_paths)
    except _EvidenceError as error:
        notes.append(error.note)
        return result('unverified', digest)
    except (OSError, *_BUNDLE_ERRORS):
        notes.append(_note('export_invalid', f'{case_id}.geometry'))
        return result('unverified', digest)
    return result('qualified', digest, facets)


def _matching(pair_id, left, right):
    notes = []
    for field in ('membrane', 'load', 'solver', 'supports', 'perturbation'):
        a, b = left[field], right[field]
        if a is None or b is None:
            notes.append(_note('match_unresolved', f'{pair_id}.{field}'))
        elif (a['common'] | a['replaced_common'] != b['common'] if field == 'supports' else a != b):
            notes.append(_note('match_difference', f'{pair_id}.{field}'))
    return notes


def _match_status(notes):
    if any(n['code'] == 'match_difference' for n in notes):
        return 'unmatched'
    return 'blocked' if notes else 'eligible'


def _screen(contrasts, cases, settings):
    notes = [reason for pair in contrasts.values() if pair['status'] != 'eligible' for reason in pair['reasons']]
    clamped = [k for k, case in cases.items() if case['role'] == 'clamped']
    if all(cases[k]['status'] == 'qualified' for k in clamped):
        a, b = (settings[k] for k in clamped)
        for field in ('x_c', 'hardware', 'connections'):
            if a['clamp'][field] != b['clamp'][field]:
                notes.append(_note('match_difference', f'screen.clamp.{field}'))
        if a['clamp']['d_c'] == b['clamp']['d_c']:
            notes.append(_note('match_difference', 'screen.clamp.distinct_preloads'))
        if a['supports'] is None or b['supports'] is None:
            notes.append(_note('match_unresolved', 'screen.attachment'))
        elif any(a['supports'][k] != b['supports'][k] for k in ('attachment', 'replaced_common', 'patch')):
            notes.append(_note('match_difference', 'screen.attachment'))
    return {'status': _match_status(notes), 'reasons': _sorted_notes(notes)}


def _inventory(study, path, identity, *, keep_surfaces):
    synthetic = study['source_kind'] == 'synthetic'
    cases, surfaces, settings, cache = {}, {}, {}, {}
    for case_id, case in study['cases'].items():
        if synthetic:
            checked = []
            status, notes, digest, facets = _bundle(case_id, path.parent, case['geometry'], checked)
            cases[case_id] = {'role': case['role'], 'status': status, 'reasons': notes,
                              'checked_files': checked, 'diagnostics': None, 'geometry_sha256': digest}
        else:
            cases[case_id], facets, settings[case_id] = _real_case(path.parent, case_id, case, cache,
                                                                 keep_surfaces=keep_surfaces)
        if keep_surfaces and facets is not None:
            surfaces[case_id] = facets
    contrasts = {}
    for pair_id, pair in study['contrasts'].items():
        notes = [_note('case_unqualified', member) for member in pair.values()
                 if cases[member]['status'] != 'qualified']
        if not notes and not synthetic:
            notes = _matching(pair_id, settings[pair['case_id']], settings[pair['control_id']])
        contrasts[pair_id] = {**pair, 'status': _match_status(notes), 'reasons': _sorted_notes(notes)}
    screen = {'status': 'not_applicable', 'reasons': []} if synthetic else _screen(contrasts, cases, settings)
    return {'schema_version': 1, 'source_kind': study['source_kind'], 'study': identity,
            'cases': cases, 'contrasts': contrasts, 'screen': screen}, surfaces


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


def _completion(evidence, loads, comparisons, *, synthetic):
    if synthetic:
        return {'contrasts': {key: {'status': 'not_applicable', 'reasons': []} for key in evidence['contrasts']},
                'screen': {'status': 'not_applicable', 'reasons': []}}
    decisions = {}
    for key, pair in evidence['contrasts'].items():
        case, ref = pair['case_id'], pair['control_id']
        notes = [] if pair['status'] == 'eligible' else [_note('ineligible', key)]
        for attitude, case_load in loads[case].items():
            prefix = f'{key}.{attitude}'
            for member, load in ((case, case_load), (ref, loads[ref][attitude])):
                if load['validation'] != 'passed' or load['normalized'] is None:
                    why = load['reasons'] or [_note('case_unavailable', 'normalized')]
                    notes.extend(_note(n['code'], f'{prefix}.{member}.{n["target"]}') for n in why)
            delta = comparisons[case][attitude]['control']
            if delta['normalized'] is None:
                notes.extend(_note(n['code'], f'{prefix}.{n["target"]}') for n in delta['reasons']
                             if n['code'] == 'nonfinite' and n['target'].startswith('normalized'))
        decisions[key] = {'status': 'not_complete' if notes else 'complete', 'reasons': _sorted_notes(notes)}
    n_complete = sum(d['status'] == 'complete' for d in decisions.values())
    eligible = evidence['screen']['status'] == 'eligible'
    notes = [] if eligible else [_note('ineligible', 'screen')]
    notes.extend(n for d in decisions.values() for n in d['reasons'])
    status = ('complete' if eligible and n_complete == len(decisions) else
              'partial' if n_complete else 'not_complete')
    return {'contrasts': decisions, 'screen': {'status': status, 'reasons': _sorted_notes(notes)}}


def _evaluate(study, evidence, surfaces, evidence_identity):
    synthetic = study['source_kind'] == 'synthetic'
    optics, pressure = _OPTICS[study['optics']], study.get('pressure_pa')
    rays, placements = _rays(), _placements()
    flat = Facets(_FLAT, np.array([1], dtype=np.int8)) if synthetic or 'parent' in study else None
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
    completion = _completion(evidence, loads, comparisons, synthetic=synthetic)
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
        study = _study(raw, allow_linked=True)
        _check_parent(study, path)
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
