"""Read one receipt-bound Abaqus ODB, without changing it or opening any other ODB."""
import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import traceback

import numpy as np
from inspect_sail_coordinate_storage import identity, with_odb

ENERGIES = ('ALLSD', 'ALLIE', 'ALLSE')


def save(path, value):
    temporary = path.with_suffix('.tmp')
    with temporary.open('w', encoding='utf-8') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')
    temporary.replace(path)


def data(v):
    return v.dataDouble if str(v.precision) == 'DOUBLE_PRECISION' else v.data


def scalar(x):
    a = float(x)
    return a if math.isfinite(a) else None


def describe(x, ids=None):
    x = np.asarray(x, dtype=np.float64)
    ok = np.isfinite(x)
    result = {'count': len(x), 'finite': int(ok.sum())}
    if not ok.any():
        return result
    safe = x[ok]
    lo, hi = int(safe.argmin()), int(safe.argmax())
    result.update(min=float(safe[lo]), max=float(safe[hi]), rms=float(np.sqrt(np.mean(safe**2))),
                  mean=float(safe.mean()))
    if ids is not None:
        labels = np.asarray(ids)[ok]
        result.update(min_label=int(labels[lo]), max_label=int(labels[hi]))
    return result


def node_values(field, instance, kind):
    f = field.getSubset(region=instance)
    stem = 'COOR' if kind == 'COORD' and 'COOR3' in f.componentLabels else kind
    index = tuple(f.componentLabels).index(stem + '3')
    rows = sorted((v.nodeLabel, float(data(v)[index])) for v in f.values if str(v.position) == 'NODAL')
    labels = np.asarray([i for i, _ in rows], dtype=np.int64)
    if len(labels) != len(np.unique(labels)):
        raise ValueError('duplicate nodal labels in ' + kind)
    return labels, np.asarray([x for _, x in rows], dtype=np.float64)


def energy_region(region):
    values = {name: [[scalar(t), scalar(v)] for t, v in region.historyOutputs[name].data]
              for name in ENERGIES if name in region.historyOutputs}
    ie = {time: value for time, value in values.get('ALLIE', []) if time is not None}
    ratios = []
    for t, sd in values.get('ALLSD', []):
        denom = ie.get(t)
        if sd is not None and denom is not None and denom > 0:
            ratio = scalar(sd / denom)
            if ratio is not None:
                ratios.append([t, ratio])
    return {'history': values, 'ratio_exact_matched_times': ratios,
            'max_ratio': max(ratios, key=lambda row: row[1]) if ratios else None}


def stress_summary(field, membrane):
    f = field.getSubset(region=membrane)
    columns = tuple(f.componentLabels)
    axes = [columns.index(k) for k in ('S11', 'S22', 'S12')]
    groups = defaultdict(list)
    for v in f.values:
        if str(v.position) != 'INTEGRATION_POINT':
            continue
        sp = v.sectionPoint
        x = data(v)
        a, b, s = (float(x[i]) for i in axes)
        center = (a + b) / 2
        radius = math.hypot((a - b) / 2, s)
        groups[(sp.number if sp else None, sp.description if sp else None)].append(
            (v.elementLabel, center - radius, center + radius))
    result = []
    for (number, description), rows in groups.items():
        labels, lo, hi = np.asarray(rows).T
        finite = np.isfinite(lo)
        result.append(dict(number=number, description=description, minimum_in_plane=describe(lo, labels),
                           maximum_in_plane=describe(hi, labels),
                           fraction_min_compressive=float(np.mean(lo[finite] < 0)) if finite.any() else None))
    return result


def node_field(field):
    labels = tuple(field.componentLabels)
    rows = defaultdict(dict)
    instances = Counter()
    for v in field.values:
        if str(v.position) != 'NODAL':
            continue
        inst = v.instance.name if v.instance is not None else None
        instances[inst] += 1
        if v.nodeLabel in rows[inst]:
            raise ValueError('duplicate nodal vector: ' + str((inst, v.nodeLabel)))
        rows[inst][v.nodeLabel] = [scalar(x) for x in data(v)]
    return labels, rows, instances


def force_comparison(vf, structural):
    # Each element-nodal structural contribution is first assembled at its
    # own node. Norms are then compared on common support, never net resultants.
    vl, vv, _ = node_field(vf)
    sl = tuple(structural.componentLabels)
    vi = [vl.index('VF' + str(i)) for i in (1, 2, 3)]
    si = [sl.index('NFORC' + str(i)) for i in (1, 2, 3)]
    sf = defaultdict(lambda: np.zeros(3))
    for v in structural.values:
        if str(v.position) not in ('ELEMENT_NODAL', 'NODAL'):
            continue
        val = data(v)
        sf[(v.instance.name if v.instance is not None else None, v.nodeLabel)] += [val[i] for i in si]
    result = []
    for inst, group in vv.items():
        points = [(np.asarray([value[i] for i in vi], dtype=float), sf[(inst, node)])
                  for node, value in group.items() if (inst, node) in sf]
        if not points:
            result.append(dict(instance=inst, vf_nodes=len(group), common_nodes=0))
            continue
        axes = np.asarray(points)
        vl2, sl2 = (scalar(v) for v in np.linalg.norm(axes, axis=(0, 2)))
        result.append(dict(instance=inst, vf_nodes=len(group), common_nodes=len(points),
                           vf_l2=vl2, structural_l2=sl2, ratio=vl2/sl2 if sl2 and vl2 is not None else None))
    return result


def inventory(odb):
    return {'instances': [{'name': name, 'nodes': len(inst.nodes), 'elements': len(inst.elements)}
                          for name, inst in odb.rootAssembly.instances.items()],
            'steps': [{'name': name, 'procedure': str(step.procedure), 'domain': str(step.domain),
                       'time_period': scalar(step.timePeriod), 'prior_total_time': scalar(step.totalTime),
                       'history': [{'name': key, 'description': reg.description, 'position': str(reg.position),
                                    'outputs': {k: len(v.data) for k, v in reg.historyOutputs.items()}}
                                   for key, reg in step.historyRegions.items()],
                       'frames': [{'number': j, 'increment': fr.incrementNumber, 'time': scalar(fr.frameValue),
                                   'fields': {k: {'components': list(v.componentLabels),
                                                  'positions': sorted({str(loc.position) for loc in v.locations})}
                                              for k, v in fr.fieldOutputs.items()}}
                                  for j, fr in enumerate(step.frames)]} for name, step in odb.steps.items()]}


def extract(odb, request, out):
    inv = inventory(odb)
    save(out / 'inventory.json', inv)
    instance = odb.rootAssembly.instances[request['instance']]
    terminal = odb.steps[request['step']].frames[request['frame']]
    if request['frame'] != len(odb.steps[request['step']].frames) - 1 or terminal.incrementNumber != request['increment']:
        raise ValueError('terminal frame does not match bound receipt')
    histories, frames, reactions, cables, errors, arrays = {}, [], {}, {}, [], {}
    initial = {}

    def optional(label, fn):
        try:
            return fn()
        except Exception:
            errors.append({'section': label, 'traceback': traceback.format_exc()})
            return None

    for step_name, step in odb.steps.items():
        histories[step_name] = {name: energy_region(region) for name, region in step.historyRegions.items()
                                if set(ENERGIES).intersection(region.historyOutputs.keys())}
        for j, fr in enumerate(step.frames):
            end = step_name == request['step'] and j == request['frame']
            row = dict(step=step_name, frame=j, increment=fr.incrementNumber, time=scalar(fr.frameValue))
            for kind in ('U', 'COORD'):
                if kind not in fr.fieldOutputs:
                    continue
                pair = optional((step_name, j, kind), lambda: node_values(fr.fieldOutputs[kind], instance, kind))
                if pair is None:
                    continue
                labels, a = pair
                row[kind] = describe(a, labels)
                row[kind]['full_instance'] = len(labels) == len(instance.nodes)
                if len(labels):
                    first = initial.setdefault(kind, (labels, a))
                    if np.array_equal(labels, first[0]):
                        row[kind]['change_from_first_recorded'] = describe(a - first[1], labels)
                    if end:
                        arrays[kind + '_label'] = labels
                        arrays[kind + '_terminal'] = a
                        if np.array_equal(labels, first[0]):
                            arrays[kind + '_initial'] = first[1]
            if 'VF' in fr.fieldOutputs and 'NFORC' in fr.fieldOutputs:
                row['force_comparison'] = optional((step_name, j, 'VF/NFORC'),
                                                   lambda: force_comparison(fr.fieldOutputs['VF'], fr.fieldOutputs['NFORC']))
            frames.append(row)
    for kind in ('RF', 'RM', 'VF'):
        if kind not in terminal.fieldOutputs:
            continue
        result = optional(kind, lambda: node_field(terminal.fieldOutputs[kind]))
        if result is not None:
            components, values, count = result
            reactions[kind] = {'components': components, 'instances': [
                {'instance': name, 'node_count': count[name],
                 'nonzero_or_nonfinite': {k: a for k, a in rows.items() if any(v is None or v != 0 for v in a)}}
                for name, rows in values.items()]}
    for name, inst in odb.rootAssembly.instances.items():
        if 'CABLE' not in name.upper():
            continue
        detail = {'reference_nodes': [[int(n.label), [scalar(v) for v in n.coordinates]] for n in inst.nodes],
                  'elements': [[int(e.label), str(e.type), list(map(int, e.connectivity))] for e in inst.elements],
                  'fields': {}}
        for kind in ('S', 'SF', 'NFORC', 'EVOL', 'COORD'):
            if kind not in terminal.fieldOutputs:
                continue

            def cable_field():
                f = terminal.fieldOutputs[kind].getSubset(region=inst)
                rows = []
                for v in f.values:
                    pos = str(v.position)
                    a = data(v)
                    rows.append(dict(position=pos, label=v.nodeLabel if pos == 'NODAL' else v.elementLabel,
                                     values=[scalar(x) for x in np.atleast_1d(a)]))
                return dict(components=list(f.componentLabels), records=rows)

            result = optional((name, kind), cable_field)
            if result is not None:
                detail['fields'][kind] = result
        cables[name] = detail
    stress = optional('membrane in-plane S', lambda: stress_summary(terminal.fieldOutputs['S'], instance)) \
        if 'S' in terminal.fieldOutputs else None
    if arrays:
        np.savez_compressed(out / 'membrane-endpoints.npz', **arrays)
    return dict(case=request['case'], histories=histories, frames=frames, terminal_reactions=reactions,
                terminal_cables=cables, terminal_membrane_stress=stress, optional_errors=errors)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--request', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    request = json.loads(args.request.read_text(encoding='utf-8-sig'))
    expected = {'control-A': ('control-ledger-seed.inp', 'SailStorageA_Control_FULL01.odb'),
                'clamp-A': ('tieprint-ledger-z.inp', 'SailStorageA_Clamp_FULL01.odb')}
    names = expected[request['case']]
    pins = {k: identity(request[k]['path']) for k in ('input', 'odb')}
    for (key, pin), name in zip(pins.items(), names):
        p = Path(pin['path'])
        if (pin['sha256'] != request[key]['sha256'] or p.name != name or p.parent.name != request['case'] or
                p.parent.parent.name != 'sail-a-pair-output-precision1'):
            raise ValueError('run receipt mismatch: ' + str(pin))
    if Path(pins['odb']['path']).with_suffix('.lck').exists():
        raise ValueError('ODB is locked before extraction')
    args.out.mkdir(parents=False, exist_ok=False)
    save(args.out / 'before.json', pins)
    from odbAccess import openOdb
    try:
        result = with_odb(openOdb, pins['odb']['path'], lambda odb: extract(odb, request, args.out))
        save(args.out / 'extracted.json', result)
    finally:
        after = {k: identity(pin['path']) for k, pin in pins.items()}
        save(args.out / 'after.json', after)
        if after != pins or Path(pins['odb']['path']).with_suffix('.lck').exists():
            raise ValueError('source changed or a lock appeared during extraction')


if __name__ == '__main__':
    main()
