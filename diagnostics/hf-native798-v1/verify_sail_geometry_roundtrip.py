"""Independently inspect one named ODB and verify an identity-mapped geometry export."""

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import numpy as np

_WIDTH = {'S8R5': (8, 4), 'STRI65': (6, 3)}
_ENERGY = ('ALLIE', 'ALLSD', 'ALLAE', 'ALLWK')
_ANCHORS = ((1, (0., 0., 0.)), (400, (20., 0., 0.)), (10201, (10., 10., 0.)))


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def identity(path):
    path = Path(path).resolve(strict=True)
    with path.open('rb') as source:
        h = hashlib.sha256()
        for chunk in iter(lambda: source.read(1 << 20), b''):
            h.update(chunk)
    return {'path': str(path), 'sha256': h.hexdigest()}


def unique(pairs):
    out = {}
    for key, value in pairs:
        require(key not in out, f'duplicate key/label: {key}')
        out[key] = value
    return out


def load(path):
    return json.loads(path.read_text(encoding='utf-8'), object_pairs_hook=unique)


def membrane(path, instance=None):
    rows = [text.strip() for text in path.read_text(encoding='ascii').splitlines()
            if text.strip() and not text.lstrip().startswith('**')]
    def block(card, end):
        starts = [i for i, text in enumerate(rows) if re.fullmatch(card, text, re.I)]
        require(len(starts) == 1, f'{path.name}: missing or ambiguous {card}')
        start = starts[0]
        stop = next(i for i in range(start + 1, len(rows)) if rows[i].casefold() == end)
        return rows[start + 1:stop]

    if instance is not None:
        placement = block(r'\*instance,\s*name=' + re.escape(instance) + r',\s*part=membrane', '*end instance')
        require(not placement, 'submitted membrane instance is not an untransformed placement')
    nodes, elements = [], []
    kind = None
    for row in block(r'\*part,\s*name=membrane', '*end part'):
        if row.startswith('*'):
            if row.lower() == '*node':
                kind = 'node'
            elif row.lower().startswith('*element,'):
                match = re.fullmatch(r'\*element,\s*type=(\w+)', row, re.I)
                require(match and match[1].upper() in _WIDTH, 'unsupported membrane element card')
                kind = match[1].upper()
            else:
                kind = None
        elif kind:
            label, *fields = row.rstrip(',').split(',')
            if kind == 'node':
                coords = tuple(map(float, fields))
                require(len(coords) == 3 and all(map(math.isfinite, coords)), f'bad coordinates for {label}')
                nodes.append((int(label), coords))
            else:
                slots = tuple(map(int, fields))
                require(len(slots) == _WIDTH[kind][0] and len(set(slots)) == len(slots), f'bad element {label}')
                elements.append((int(label), (kind, slots)))
    nodes, elements = unique(nodes), unique(elements)
    require(nodes and elements and all(n in nodes for _, e in elements.values() for n in e), 'incomplete membrane')
    return nodes, elements


def reference(recipe, paths):
    points, elements = membrane(paths['reference_mesh'])
    submitted, ordered = membrane(paths['input'], recipe['instance'])
    ref, sides = load(paths['reference_identity']), load(paths['coating_sides'])
    require(ref['reference_mesh_sha256'] == recipe['files']['reference_mesh']['sha256'] and
            ref['coating_sides_sha256'] == recipe['files']['coating_sides']['sha256'], 'unbound reference receipt')
    require(ref['units'] == 'm' and ref['R'] == np.eye(3).tolist() and ref['t_m'] == [0, 0, 0] and
            ref['reflective_side'].startswith('+z in the unseeded panel frame'), 'unsupported frame or diagnostic side')
    anchors = [{'node_label': label, 'odb_xyz_m': list(x), 'panel_xyz_m': list(x)} for label, x in _ANCHORS]
    require(ref['anchors'] == anchors and len(points) == ref['node_count'] and
            {k: n for k, n in Counter(k for k, _ in elements.values()).items()} == ref['element_counts'],
            'frozen anchors or counts disagree')
    require(points.keys() == submitted.keys() and elements == ordered, 'input and reference label/topology mismatch')
    require(all(p[2] == 0. for p in points.values()), 'reference is not unseeded and planar')
    for label, expected in _ANCHORS:
        require(points[label] == expected and submitted[label][:2] == expected[:2], f'anchor {label} disagrees')
    labels = sorted(points)
    # Reference-grid decimals and submitted X/Y are compared at native mesh storage only.
    # This check is not used on exported coordinates.
    x = np.asarray([points[n][:2] for n in labels], dtype=np.float32)
    y = np.asarray([submitted[n][:2] for n in labels], dtype=np.float32)
    require(np.array_equal(x, y), 'submitted X/Y do not preserve the unseeded reference at mesh storage precision')
    signs = {}
    for label, (kind, slots) in elements.items():
        p = [points[n] for n in slots[:_WIDTH[kind][1]]]
        cross = math.fsum(a[0]*b[1] - a[1]*b[0] for a, b in zip(p, p[1:] + p[:1]))
        require(cross != 0. and sides.get(str(label)) == ('SPOS' if cross > 0. else 'SNEG'),
                f'reference winding/coating mismatch for element {label}')
        signs[label] = 1 if cross > 0. else -1
    require(len(sides) == len(signs), 'unexpected coating-side rows')
    facts = {'anchors': anchors, 'units': 'm', 'mapping_kind': 'identity', 'placement': 'empty membrane *INSTANCE placement',
             'reference_plane': 'unseeded z=0; submitted z may be seeded', 'matched_reference_labels': len(points),
             'matched_reference_ordered_elements': len(elements), 'unchanged_xy_at_mesh_storage': True,
             'sign_counts': dict(Counter(signs.values()))}
    return elements, points, signs, facts


def native_field(field, instance, labels):
    points, precisions = [], Counter()
    # Independently filter the native field by its reported position and owner.
    # No exporter subdivision, mapping, sign or field-subset helper is used.
    for v in field.values:
        if str(v.position) != 'NODAL' or v.instance is None or v.instance.name != instance:
            continue
        precision = str(v.precision)
        require(precision in ('SINGLE_PRECISION', 'DOUBLE_PRECISION'), 'unexpected precision')
        xyz = tuple(float(x) for x in (v.data if precision == 'SINGLE_PRECISION' else v.dataDouble))
        require(len(xyz) == 3 and all(map(math.isfinite, xyz)), f'invalid native vector at {v.nodeLabel}')
        points.append((v.nodeLabel, xyz))
        precisions[precision] += 1
    data = unique(points)
    require(data.keys() == labels, 'instance-scoped nodal field coverage differs from the reference')
    return data, dict(precisions)


def max_z(data):
    label = min(data, key=lambda n: (-abs(data[n][2]), n))
    z = data[label][2]
    return {'node_label': label, 'signed_m': z, 'absolute_m': abs(z)}


def energies(step):
    histories, ratios = {}, {}
    def finite(value):
        number = float(value)
        return number if math.isfinite(number) else None

    for name in step.historyRegions.keys():
        region = step.historyRegions[name]
        series = {}
        for key in _ENERGY:
            if key in region.historyOutputs.keys():
                series[key] = [[finite(t), finite(y)] for t, y in region.historyOutputs[key].data]
        if not series:
            continue
        histories[name] = {'description': str(region.description), 'series': series}
        times = {}
        for key in ('ALLSD', 'ALLIE'):
            counts = Counter(t for t, _ in series.get(key, []))
            times[key] = {t: y if counts[t] == 1 else None for t, y in series.get(key, []) if t is not None}
        ratio = []
        for t in sorted(times['ALLSD'].keys() | times['ALLIE'].keys()):
            s, e = (times[k].get(t) for k in ('ALLSD', 'ALLIE'))
            reason = ('no unique finite matched-time pair' if s is None or e is None else
                      'ALLIE is zero' if e == 0. else None)
            value = None if reason else finite(s/e)
            ratio.append({'time': t, 'ratio': value, 'reason': reason or ('nonfinite quotient' if value is None else None)})
        ratios[name] = ratio
    return {'time_units': 's (analysis step time)', 'energy_units': 'J (SI model units)', 'regions': histories,
            'missing': [k for k in _ENERGY if not any(k in r['series'] for r in histories.values())], 'ALLSD_over_ALLIE': ratios}


def compare_arrays(directory, nodes, topo, signs):
    meta = load(directory/'geometry.json')
    path = directory/'geometry.npz'
    require(meta['npz'] == {'path': 'geometry.npz', 'sha256': identity(path)['sha256']}, 'geometry.npz digest mismatch')
    counts = dict(Counter(k for k, _ in topo.values()))
    require(meta['node_count'] == len(nodes) and meta['element_counts'] == counts, 'geometry meta counts disagree with ODB')
    with np.load(path, allow_pickle=False) as data:
        keys = ['node_labels', 'node_xyz_m'] + [f'{f}_{k}' for k in counts for f in ('element_labels','connectivity','face_sign')]
        require(sorted(data.files) == sorted(keys), 'NPZ member names are not the expected unique groups')
        def array(key, dtype, shape):
            a = data[key]
            require(isinstance(a, np.ndarray) and a.dtype == np.dtype(dtype) and a.shape == shape, f'bad array {key}')
            return a
        labels = array('node_labels', 'int64', (len(nodes),))
        xyz = array('node_xyz_m', 'float64', (len(nodes),3))
        require(len(set(labels)) == len(labels) and set(labels) == nodes.keys(), 'export node label coverage disagrees')
        expected = np.array([nodes[n] for n in labels], dtype=np.float64)
        # Identity R, zero t: the independent R*x+t is exactly the promoted native value.
        error = np.abs(xyz - expected)
        bound = (16 * 2.**-52) * np.abs(expected)
        require(np.isfinite(xyz).all() and np.all(error <= bound), 'export coordinate exceeds mapping roundoff')
        known = set()
        for kind, n in counts.items():
            ids = array('element_labels_'+kind, 'int64', (n,))
            links = array('connectivity_'+kind, 'int64', (n,_WIDTH[kind][0]))
            faces = array('face_sign_'+kind, 'int8', (n,))
            for label, row, sign in zip(ids, links, faces):
                require(label not in known and topo.get(label) == (kind, tuple(row)), f'ordered topology mismatch: {label}')
                require(sign == signs[label], f'row face sign mismatch: {label}')
                known.add(label)
        require(known == topo.keys(), 'topology coverage disagrees')
    return {'matched_nodes': len(nodes), 'matched_components': 3 * len(nodes), 'matched_element_counts': counts,
            'matched_connectivity_slots': sum(len(slots) for _, slots in topo.values()),
            'matched_face_signs': len(known), 'max_component_error_m': error.max(axis=0).tolist(),
            'max_component_bound_m': bound.max(axis=0).tolist(),
            'bitwise_equal_including_signed_zero': bool(np.array_equal(xyz.view(np.uint64), expected.view(np.uint64)))}


def inspect(recipe, paths, bundle):
    import odbAccess
    import sail_optics  # Runtime import check only. It is not an independent geometry oracle.
    exe = Path(sys.executable).resolve()
    require(tuple(x.lower() for x in exe.parts[-4:]) == ('win_b64','code','bin','smapython.exe') and
            exe.parts[-6].lower() == 'estproducts' and re.fullmatch(r'20[0-9]{2}', exe.parts[-5]), 'not installed Abaqus Python')
    runtime = dict(abaqus=exe.parts[-5], python=sys.version.split()[0], numpy=np.__version__)
    origins = {module.__name__: identity(module.__file__) for module in (odbAccess, np, sail_optics)}
    require(origins['sail_optics']['sha256'] == recipe['portable_sha256'], 'unverified portable import')
    expected, reference_nodes, signs, frame_basis = reference(recipe, paths)
    odb = odbAccess.openOdb(path=str(paths['odb']), readOnly=True)
    try:
        names = list(odb.steps.keys())
        require(names == [recipe['terminal_step']], 'unexpected ODB step inventory')
        step = odb.steps[names[-1]]
        require(str(step.domain) == 'TIME' and float(step.timePeriod) == 1., 'not a unit-duration TIME step')
        frames = [{'index': i, 'time': float(f.frameValue), 'increment': int(f.incrementNumber),
                   'fields': sorted(f.fieldOutputs.keys())} for i, f in enumerate(step.frames)]
        terminal = frames[-1]
        require(terminal['time'] == 1. and terminal['increment'] == recipe['terminal_increment'] and
                'COORD' in terminal['fields'], 'last frame does not match the completed run/required COORD')
        frame = step.frames[-1]
        inst = odb.rootAssembly.instances[recipe['instance']]
        require(inst.name == recipe['instance'], 'unexpected membrane owner')
        node_labels = [n.label for n in inst.nodes]
        require(len(set(node_labels)) == len(node_labels) and set(node_labels) == reference_nodes.keys(), 'ODB node labels differ')
        topo = unique((e.label, (e.type, tuple(e.connectivity))) for e in inst.elements)
        require(topo == expected, 'ODB typed ordered connectivity differs from submitted/unseeded reference')
        native_xyz, precisions = native_field(frame.fieldOutputs['COORD'], inst.name, reference_nodes.keys())
        coord_precision = 'mixed' if len(precisions) != 1 else next(iter(precisions)).split('_')[0].lower()
        u3 = {'reason':'U field not present in the terminal frame'}
        if 'U' in terminal['fields']:
            try:
                u, u_precision = native_field(frame.fieldOutputs['U'], inst.name, reference_nodes.keys())
                u3 = {'max':max_z(u), 'precision_counts':u_precision}
            except ValueError as ex:
                u3 = {'reason': str(ex)}
        observation = {'instance':inst.name, 'step':step.name, 'step_duration':float(step.timePeriod), 'frames':frames,
                       'terminal':terminal, 'node_count':len(node_labels), 'element_counts':dict(Counter(e[0] for e in topo.values())),
                       'coord_precision':coord_precision, 'coord_precision_counts':precisions,
                       'max_abs_coord_z':max_z(native_xyz), 'terminal_u3':u3, 'step_histories':energies(step),
                       'job_data':{k:str(getattr(odb.jobData,k,None)) for k in ('analysisCode','version','precision')}}
        if bundle is not None:
            meta = load(bundle/'geometry.json')
            exp = meta['export']
            manifest = load(paths['manifest'])
            case = manifest['cases'][recipe['case_id']]
            basis = case['geometry_basis']['value']
            require(basis == {'units':'m', 'R':np.eye(3).tolist(), 't_m':[0,0,0], 'anchors':frame_basis['anchors'],
                              'sides':load(paths['coating_sides'])}, 'manifest basis differs from the independent reference')
            basis_bytes = json.dumps(basis, ensure_ascii=True, sort_keys=True, separators=(',',':'), allow_nan=False).encode('ascii')
            require(meta['schema_version'] == 1 and meta['source_kind'] == 'odb' and meta['case_id'] == recipe['case_id'] and
                    meta['frame'] == 'panel_local_m' and exp == dict(odb=recipe['files']['odb'],job_id=recipe['job'],
                    attempt_id=recipe['attempt'], instance=inst.name, step=step.name, step_duration=1.,
                    frame_index=terminal['index'],frame_value=terminal['time'],coord_precision=coord_precision,
                    basis_sha256=hashlib.sha256(basis_bytes).hexdigest(),runtime=runtime,adapter_sha256=recipe['adapter_sha256']),
                    'bundle identity, runtime or terminal metadata disagree')
            observation['comparison'] = compare_arrays(bundle, native_xyz, topo, signs)
    finally:
        odb.close()
    return {'runtime':runtime, 'executable':str(exe), 'import_origins':origins,
            'helper':identity(__file__), 'reference':frame_basis, 'actual':observation}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--request', required=True, type=Path)
    p.add_argument('--bundle', type=Path)
    p.add_argument('--out', required=True, type=Path)
    args = p.parse_args()
    require(not args.out.exists() and not args.out.is_symlink(), 'use a fresh report path')
    recipe = load(args.request)
    paths = {k: (args.request.resolve().parent/v['path']).resolve(strict=True) for k,v in recipe['files'].items()}
    def snapshot():
        return {k:identity(path) for k,path in paths.items()}
    before = snapshot()
    require(all(before[k]['sha256'] == ref['sha256'] for k,ref in recipe['files'].items()), 'source file digest mismatch')
    result = inspect(recipe, paths, args.bundle)
    after = snapshot()
    require(before == after, 'source identity changed during read-only access')
    with args.out.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(dict(result, before=before, after=after, read_only=True), stream, indent=2, allow_nan=False)
        stream.write('\n')
    print('Identity-mapped ODB comparison passed.' if args.bundle else 'Terminal native inspection complete.')


if __name__ == '__main__':
    main()
