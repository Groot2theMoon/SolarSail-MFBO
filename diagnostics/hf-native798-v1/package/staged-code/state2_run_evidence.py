"""Diagnostic-only evidence and fail-closed checks for these generated decks."""
import hashlib
import json
import math
import os
from pathlib import Path


TYPES = {'S8R5': 8, 'STRI65': 6}
ANCHORS = ((0., 0., 0.), (20., 0., 0.), (10., 10., 0.))


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, payload):
    with open(path, 'x', encoding='utf-8', newline='\n') as fh:
        json.dump(payload, fh, sort_keys=True, indent=2, allow_nan=False)
        fh.write('\n')


def header(line):
    key, *options = (part.strip().upper() for part in line.split(','))
    return key, dict((k.strip(), v.strip()) for opt in options if '=' in opt
                     for k, v in (opt.split('=', 1),))


def membrane(path):
    """Read the project's generated inline Membrane part, not general Abaqus decks."""
    nodes, elements = {}, {}
    active, kind, found = False, None, 0
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        line = line.strip()
        if not line or line.startswith('**'):
            continue
        if line.startswith('*'):
            key, options = header(line)
            if key == '*PART':
                active = options.get('NAME') == 'MEMBRANE'
                found += int(active)
            elif key == '*END PART':
                active = False
            elif key == '*INCLUDE':
                raise ValueError('include not supported in this diagnostic deck')
            kind = 'node' if key == '*NODE' else options.get('TYPE') if key == '*ELEMENT' else None
            if active and key == '*ELEMENT':
                require(kind in TYPES, 'unexpected membrane element type: %r' % kind)
            continue
        if active and kind:
            label, *fields = (field.strip() for field in line.rstrip(',').split(','))
            label = int(label)
            if kind == 'node':
                xyz = tuple(float(x) for x in fields)
                require(len(xyz) == 3 and all(map(math.isfinite, xyz)) and label not in nodes,
                        'invalid/duplicate membrane node %s' % label)
                nodes[label] = xyz
            else:
                conn = tuple(map(int, fields))
                require(len(conn) == TYPES[kind] and label not in elements,
                        'invalid/duplicate membrane element %s' % label)
                elements[label] = (kind, conn)
    require(found == 1 and nodes and elements, 'expected one nonempty Membrane part')
    return nodes, elements


def capture_nodes(part):
    return {int(n.label): tuple(map(float, n.coordinates)) for n in part.nodes}


def close(a, b):
    # Allow CAE coordinate storage and the deck's decimal serialization, not a shape tolerance.
    return math.isclose(a, b, rel_tol=2e-7, abs_tol=5e-15)


def write_applied(before, after, requested):
    ref, _ = membrane(os.environ['STATE2_REFERENCE_MESH'])
    require(before.keys() == after.keys() == requested.keys() == ref.keys(),
            'applied perturbation does not cover the exact reference node labels')
    with open('applied-nodes.tsv', 'x', encoding='ascii', newline='\n') as fh:
        fh.write('label\tbefore_x\tbefore_y\tbefore_z\tafter_x\tafter_y\tafter_z\tdx\tdy\tdz\trequested_dz\n')
        for label in sorted(ref):
            p, q = before[label], after[label]
            d = tuple(y - x for x, y in zip(p, q))
            dz = float(requested[label])
            require(all(close(a, b) for a, b in zip(p, ref[label]))
                    and d[:2] == (0., 0.) and close(d[2], dz),
                    'applied perturbation mismatch at node %d' % label)
            fh.write('%d\t%s\n' % (label, '\t'.join(x.hex() for x in (*p, *q, *d, dz))))


def check_deck(path, seeded):
    ref_nodes, ref_elements = membrane(os.environ['STATE2_REFERENCE_MESH'])
    nodes, elements = membrane(path)
    require(ref_nodes.keys() == nodes.keys() and ref_elements == elements,
            'submitted membrane label/type/ordered-connectivity mismatch')
    if seeded:
        with open('applied-nodes.tsv', encoding='ascii') as fh:
            next(fh)
            expected = {int(row[0]): tuple(map(float.fromhex, row[4:7]))
                        for row in (line.split() for line in fh)}
    else:
        expected = ref_nodes
    require(expected.keys() == nodes.keys() and all(
        close(a, b) for label in nodes for a, b in zip(nodes[label], expected[label])),
        'submitted membrane coordinates differ from intended or applied coordinates')

    step = position = field = None
    steps, has_output, instance = [], set(), False
    gt_static, controls = [], set()
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        line = line.strip()
        if not line or line.startswith('**'):
            continue
        if line.startswith('*'):
            key, options = header(line)
            if key == '*INSTANCE' and options.get('NAME') == 'MEMBRANE-1':
                require(options.get('PART') == 'MEMBRANE' and not instance,
                        'unexpected membrane instance binding')
                instance = True
                position = 'placement'
            else:
                position = key
            if key == '*STEP':
                step = options['NAME']
                steps.append(step)
                field = None
            elif key == '*STATIC' and step == 'STEP-GLOBALTENSION':
                gt_static.append(float(options.get('STABILIZE', 'nan')))
            elif key == '*CONTROLS' and step == 'STEP-GLOBALTENSION':
                controls.add(options.get('PARAMETERS'))
            elif key == '*OUTPUT':
                field = 'FIELD' in (opt.strip().upper() for opt in line.split(',')[1:])
        elif position == 'placement':
            raise ValueError('unexpected membrane instance translation or rotation')
        elif position == '*NODE OUTPUT' and field:
            has_output.update((step, k.strip().upper()) for k in line.split(','))
    require(instance, 'named membrane instance absent')
    if seeded:
        require(steps == ['STEP-GLOBALTENSION'] and ('STEP-GLOBALTENSION', 'COORD') in has_output,
                'HF terminal step or nodal COORD request differs from the diagnostic plan')
        require(gt_static == [0.0002] and
                {'LINE SEARCH', 'TIME INCREMENTATION', 'FIELD'} <= controls,
                'HF first-rung stabilization or legacy controls missing')
    else:
        require(steps == ['STEP-GLOBALTENSION', 'STEP-BUCKLE']
                and any((step_name, 'U') in has_output for step_name in steps),
                'buckle step sequence or nodal U request differs from the diagnostic plan')
        require(gt_static == [0.0005], 'unexpected buckle base-state stabilization')


def fresh_job(name):
    require(not any(Path('.').glob(name + '.*')), 'refusing existing job files: %s' % name)


def record_settings(tag, namespace, keys, override=None):
    data = {key: namespace[key] for key in keys}
    data.update(override or {})
    write_json('effective-%s.json' % tag, data)


def freeze_reference(out):
    import aba_grid_mesh
    out = Path(out)
    path = out / 'reference-mesh.inp'
    require(not path.exists(), 'refusing existing reference mesh')
    aba_grid_mesh.write_inp(str(path), base=20., height=10., seed_div=200., elem_order=2)
    nodes, elements = membrane(path)
    require(sorted(nodes) == list(range(1, 30502)) and sorted(elements) == list(range(1, 10101)),
            'unexpected reference label range')
    counts = {kind: 0 for kind in TYPES}
    area = 0.
    winding, edge_label = {}, {}
    for label, (kind, conn) in sorted(elements.items()):
        counts[kind] += 1
        n = TYPES[kind] // 2
        corner, mids = conn[:n], conn[n:]
        twice_area = 0.
        for i, a in enumerate(corner):
            b, m = corner[(i + 1) % n], mids[i]
            edge = tuple(sorted((a, b)))
            require(edge_label.setdefault(edge, m) == m, 'shared mid-side label differs')
            require(all(close(nodes[m][j], (nodes[a][j] + nodes[b][j]) / 2.) for j in range(3)),
                    'mid-side coordinates differ from specified slot')
            twice_area += nodes[a][0] * nodes[b][1] - nodes[b][0] * nodes[a][1]
        require(twice_area != 0., 'zero-area reference element')
        area += abs(twice_area) / 2.
        winding[str(label)] = 'SPOS' if twice_area > 0. else 'SNEG'
    require(counts == {'S8R5': 9900, 'STRI65': 200} and math.isclose(area, 100., abs_tol=1e-9),
            'wrong reference type counts or area')
    mesh_rows = [[label, *map(float.hex, nodes[label])] for label in sorted(nodes)]
    conn_rows = [[label, elements[label][0], *elements[label][1]] for label in sorted(elements)]
    def row_hash(rows):
        return hashlib.sha256(json.dumps(rows, separators=(',', ':')).encode('ascii')).hexdigest()
    anchors = [{'node_label': next(label for label, xyz in nodes.items() if xyz == target),
                'odb_xyz_m': target, 'panel_xyz_m': target} for target in ANCHORS]
    write_json(out / 'coating-sides.json', winding)
    write_json(out / 'reference-identity.json', {
        'reference_mesh_sha256': digest(path), 'nodes_sha256': row_hash(mesh_rows),
        'ordered_connectivity_sha256': row_hash(conn_rows), 'node_count': len(nodes),
        'element_counts': counts, 'area_m2': area, 'anchors': anchors, 'units': 'm',
        'R': [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 't_m': [0, 0, 0],
        'reflective_side': '+z in the unseeded panel frame (diagnostic model convention)',
        'coating_sides_sha256': digest(out / 'coating-sides.json'),
        'face_counts': {s: sum(x == s for x in winding.values()) for s in ('SPOS', 'SNEG')},
    })


if __name__ == '__main__':
    import sys
    freeze_reference(sys.argv[1])
