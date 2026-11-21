"""Export one named terminal membrane frame to the portable v1 geometry format."""

import argparse
import hashlib
import io
import json
import math
import operator
from pathlib import Path
import re
import sys

try:
    import numpy as np
    import sail_optics as optics
except (ImportError, OSError) as error:
    if __name__ != '__main__':
        raise
    print(f'export_sail_geometry: portable dependency unavailable: {error}', file=sys.stderr)
    raise SystemExit(1) from error


_ANCHORS = np.array([[0., 0., 0.], [20., 0., 0.], [10., 10., 0.]])
_PRECISION = {'SINGLE_PRECISION': ('single', 'data'), 'DOUBLE_PRECISION': ('double', 'dataDouble')}
_SIDE_SIGNS = {'SPOS': 1, 'SNEG': -1}


def _sha256(path):
    with path.open('rb') as stream:
        return optics._stream_sha256(stream)


def _reviewed(case, name):
    record = case[name]
    if record is None:
        raise optics.InputError(f'{name}: unresolved evidence')
    return record


def _run(case):
    run = optics._keys(_reviewed(case, 'run')['value'],
                       ('job_id', 'attempt_id', 'driver_revision', 'modified_source',
                        'submitted_deck', 'completed', 'instance', 'terminal_step',
                        'step_duration', 'frame_index'))
    for key in ('job_id', 'attempt_id', 'driver_revision', 'submitted_deck', 'instance', 'terminal_step'):
        optics._text(run[key], key)
    if run['modified_source'] is not None:
        optics._text(run['modified_source'], 'modified_source record ID')
    if run['attempt_id'] != case['attempt_id'] or run['completed'] is not True:
        raise optics.InputError('run must identify the selected completed attempt')
    if type(run['step_duration']) not in (int, float) or run['step_duration'] != 1.0:
        raise optics.InputError('terminal coverage unverified: documented step duration must be 1.0')
    index = run['frame_index']
    if index is not None and (type(index) is not int or index < 0):
        raise optics.InputError('expected frame_index must be a nonnegative integer or null')
    return run


def _check_references(case, run, base):
    run_refs = {ref['record_id'] for ref in case['run']['refs']}
    if run['submitted_deck'] not in run_refs:
        raise optics.InputError('reviewed run must cite the submitted deck')
    record_ids = run_refs | {ref['record_id'] for ref in case['geometry_basis']['refs']}
    if run['modified_source'] is not None:
        record_ids.add(run['modified_source'])
    checked = {}
    for record_id in sorted(record_ids):
        ref = case['records'].get(record_id)
        if ref is None:
            raise optics.InputError(f'unknown evidence record: {record_id}')
        path = optics._referent(base, ref)
        if path not in checked:
            checked[path] = _sha256(path)
        if checked[path] != ref['sha256']:
            raise optics.InputError(f'evidence record digest mismatch: {record_id}')


def _finite_array(value, shape, name):
    def numeric(item):
        return all(numeric(part) for part in item) if isinstance(item, list) else type(item) in (int, float)

    if not numeric(value):
        raise optics.InputError(f'{name} must use numeric JSON components')
    array = np.asarray(value, dtype=np.float64)
    if array.shape != shape or not np.isfinite(array).all():
        raise optics.InputError(f'{name} must be finite with shape {shape}')
    return array


def _basis(value):
    basis = optics._keys(value, ('units', 'R', 't_m', 'anchors', 'sides'))
    optics._enum(basis['units'], ('m',), 'geometry units')
    r = _finite_array(basis['R'], (3, 3), 'R')
    t = _finite_array(basis['t_m'], (3,), 't_m')
    with np.errstate(over='ignore', invalid='ignore'):
        if (not np.allclose(r.T @ r, np.eye(3), rtol=0, atol=1e-12) or
                not math.isclose(np.linalg.det(r), 1., rel_tol=0, abs_tol=1e-12)):
            raise optics.InputError('R is not a proper rotation')
    anchors = basis['anchors']
    if not isinstance(anchors, list) or len(anchors) != 3:
        raise optics.InputError('basis must name the two ordered base endpoints and apex')
    labels = set()
    for anchor, expected in zip(anchors, _ANCHORS):
        optics._keys(anchor, ('node_label', 'odb_xyz_m', 'panel_xyz_m'))
        label = anchor['node_label']
        if type(label) is not int or label in labels:
            raise optics.InputError('anchor labels must be distinct integers')
        labels.add(label)
        local = _finite_array(anchor['panel_xyz_m'], (3,), 'anchor panel_xyz_m')
        ref = _finite_array(anchor['odb_xyz_m'], (3,), 'anchor odb_xyz_m')
        with np.errstate(over='ignore', invalid='ignore'):
            if not np.array_equal(local, expected) or not np.allclose(r @ ref + t, expected, rtol=0, atol=1e-9):
                raise optics.InputError(f'anchor {label} does not map to its canonical reference position')
    sides = basis['sides']
    if not isinstance(sides, dict) or not sides:
        raise optics.InputError('basis must document the reflective side of each element')
    signs = {}
    for encoded, side in sides.items():
        if not re.fullmatch(r'(?:0|-?[1-9][0-9]*)', encoded):
            raise optics.InputError(f'invalid Abaqus element label in sides: {encoded!r}')
        signs[int(encoded)] = _SIDE_SIGNS[optics._enum(side, _SIDE_SIGNS, 'reflective side')]
    digest = hashlib.sha256(json.dumps(basis, ensure_ascii=True, sort_keys=True,
                                       separators=(',', ':'), allow_nan=False).encode('utf-8')).hexdigest()
    return r, t, labels, signs, digest


def _terminal(odb, run):
    name = run['terminal_step']
    names = list(odb.steps.keys())
    if not names or names[-1] != name:
        raise optics.InputError(f'ODB does not end at the documented terminal step {name!r}')
    step = odb.steps[name]
    if str(step.domain) != 'TIME' or float(step.timePeriod) != run['step_duration']:
        raise optics.InputError(f'{name}: ODB time period does not match the documented duration 1.0')
    index = len(step.frames) - 1
    if index < 0 or (run['frame_index'] is not None and run['frame_index'] != index):
        raise optics.InputError(f'{name}: expected index must identify the final saved frame')
    frame = step.frames[index]
    value = float(frame.frameValue)
    if not math.isfinite(value) or value < 0.99:
        raise optics.InputError(f'{name}: final frame does not reach frameValue >= 0.99')
    if 'COORD' not in frame.fieldOutputs.keys():
        raise optics.InputError(f'{name}: final frame {index} does not contain COORD')
    return frame.fieldOutputs['COORD'], index, value


def _coordinates(field, instance, labels, r, t, nodal):
    rows = {label: i for i, label in enumerate(labels)}
    coords = np.empty((len(rows), 3), dtype=np.float64)
    seen = set()
    precisions = set()
    for value in field.getSubset(region=instance, position=nodal).values:
        if str(value.position) != 'NODAL' or value.instance.name != instance.name:
            raise optics.InputError('COORD value is not nodal in the requested membrane instance')
        label = operator.index(value.nodeLabel)
        if label not in rows or label in seen:
            raise optics.InputError(f'COORD label {label} is unknown or duplicated')
        spec = _PRECISION.get(str(value.precision))
        if spec is None:
            raise optics.InputError(f'COORD label {label} has unknown precision')
        precision, attr = spec
        point = np.asarray(getattr(value, attr), dtype=np.float64)
        if point.shape != (3,) or not np.isfinite(point).all():
            raise optics.InputError(f'COORD label {label} does not have three finite components')
        coords[rows[label]] = point
        seen.add(label)
        precisions.add(precision)
    if len(seen) != len(rows):
        raise optics.InputError('COORD does not cover every membrane node')
    with np.errstate(over='ignore', invalid='ignore'):
        return coords @ r.T + t, next(iter(precisions)) if len(precisions) == 1 else 'mixed'


def _topology(instance, signs):
    groups, seen = {}, set()
    for element in instance.elements:
        kind, label = element.type, operator.index(element.label)
        if kind not in optics._WIDTHS:
            raise optics.InputError(f'unsupported element type: {kind}')
        if label not in signs or label in seen:
            raise optics.InputError(f'element {label}: missing side evidence or duplicate label')
        connectivity = [operator.index(getattr(node, 'label', node)) for node in element.connectivity]
        if len(connectivity) != optics._WIDTHS[kind]:
            raise optics.InputError(f'element {label}: incorrect connectivity width for {kind}')
        groups.setdefault(kind, []).append((label, connectivity, signs[label]))
        seen.add(label)
    if not groups:
        raise optics.InputError('empty membrane element set')
    if seen != signs.keys():
        raise optics.InputError('reflective-side evidence names elements outside the membrane')
    arrays, counts = {}, {}
    for kind, typed_rows in groups.items():
        labels, connectivity, face_signs = zip(*typed_rows)
        arrays[f'element_labels_{kind}'] = np.asarray(labels, dtype=np.int64)
        arrays[f'connectivity_{kind}'] = np.asarray(connectivity, dtype=np.int64)
        arrays[f'face_sign_{kind}'] = np.asarray(face_signs, dtype=np.int8)
        counts[kind] = len(typed_rows)
    return arrays, counts


def _geometry(odb, run, r, t, anchors, signs, nodal):
    instances = odb.rootAssembly.instances
    name = run['instance']
    if name not in instances.keys():
        raise optics.InputError(f'membrane instance {name!r} not found')
    instance = instances[name]
    if instance.name != name:
        raise optics.InputError(f'membrane instance does not match {name!r}')
    field, index, frame_value = _terminal(odb, run)
    labels = np.fromiter((operator.index(node.label) for node in instance.nodes), dtype=np.int64)
    unique_labels = set(labels)
    if not unique_labels or len(unique_labels) != len(labels) or not anchors <= unique_labels:
        raise optics.InputError('membrane nodes are empty, duplicated, or missing a reference anchor')
    xyz, precision = _coordinates(field, instance, labels, r, t, nodal)
    arrays, counts = _topology(instance, signs)
    arrays.update(node_labels=labels, node_xyz_m=xyz)
    optics._geometry(arrays, len(labels), counts, keep_facets=False)
    return arrays, counts, index, frame_value, precision


def _vendor():
    import odbAccess
    from abaqusConstants import NODAL

    # Bind the release label to the running vendor interpreter, not PATH or ODB producer metadata.
    parts = Path(sys.executable).resolve().parts
    if (len(parts) < 6 or parts[-6].lower() != 'estproducts' or
            tuple(p.lower() for p in parts[-4:]) != ('win_b64', 'code', 'bin', 'smapython.exe') or
            not re.fullmatch(r'20[0-9]{2}', parts[-5])):
        raise optics.InputError(f'cannot identify the Abaqus release from {sys.executable!r}')
    runtime = {'abaqus': parts[-5], 'python': sys.version.split()[0], 'numpy': np.__version__}
    vendor_errors = tuple(error for error in (getattr(odbAccess, 'OdbError', None),)
                          if isinstance(error, type) and issubclass(error, Exception))
    return odbAccess.openOdb, NODAL, runtime, vendor_errors


def _export(study, manifest, case_id):
    case = study['cases'].get(case_id)
    if case is None or case['odb_path'] is None:
        raise optics.InputError(f'{case_id!r}: no explicitly named ODB candidate')
    run = _run(case)
    r, t, anchors, signs, basis_digest = _basis(_reviewed(case, 'geometry_basis')['value'])
    _check_references(case, run, manifest.parent)
    path = (manifest.parent / case['odb_path']).resolve()
    digest = _sha256(path)
    open_odb, nodal, runtime, vendor_errors = _vendor()
    try:
        odb = open_odb(path=str(path), readOnly=True)
        try:
            arrays, counts, index, frame_value, precision = _geometry(odb, run, r, t, anchors, signs, nodal)
        finally:
            odb.close()
    except vendor_errors as error:
        raise optics.InputError(f'ODB access failed: {error}') from error
    if _sha256(path) != digest:
        raise optics.InputError('ODB SHA-256 changed during extraction')
    buffer = io.BytesIO()
    np.savez(buffer, **arrays)
    npz = buffer.getbuffer()
    meta = {'schema_version': 1, 'source_kind': 'odb', 'case_id': case_id, 'frame': 'panel_local_m',
            'npz': optics._identity('geometry.npz', hashlib.sha256(npz).hexdigest()),
            'node_count': len(arrays['node_labels']), 'element_counts': counts,
            'export': {'odb': optics._identity(path, digest), 'job_id': run['job_id'],
                       'attempt_id': run['attempt_id'], 'instance': run['instance'], 'step': run['terminal_step'],
                       'step_duration': 1., 'frame_index': index, 'frame_value': frame_value,
                       'coord_precision': precision, 'basis_sha256': basis_digest, 'runtime': runtime,
                       'adapter_sha256': _sha256(Path(__file__).resolve())}}
    return {'geometry.npz': npz, 'geometry.json': optics._encoded(meta)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--case', required=True)
    parser.add_argument('--out-dir', required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        manifest = args.manifest.resolve()
        out = optics._output_path(args.out_dir)
        if args.out_dir.is_symlink() or out.exists() or not out.parent.is_dir():
            raise optics.InputError('export requires a new output directory under an existing parent')
        study = optics._study(manifest.read_bytes())
        if study['source_kind'] != 'odb':
            raise optics.InputError('export requires an initial ODB manifest')
        optics._text(args.case, 'case ID')
    except (OSError, ValueError, RuntimeError, OverflowError) as error:
        print(f'export_sail_geometry: {error}', file=sys.stderr)
        return 2
    try:
        items = _export(study, manifest, args.case)
    except (OSError, ValueError, RuntimeError, OverflowError, ImportError, KeyError, AttributeError, TypeError) as error:
        print(f'export_sail_geometry: {error}', file=sys.stderr)
        return 1
    try:
        optics._publish(out, items)
    except OSError as error:
        print(f'export_sail_geometry: {error}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
