"""Summarize available saved terminal fields of one pinned Standard ODB."""
import argparse
from collections import defaultdict
import json
from pathlib import Path

from extract_nonlinear_pair import data, describe, energy_region, save
from inspect_sail_coordinate_storage import identity, with_odb


def components(field, instance):
    """Keep native section points and integration-point components distinct."""
    subset = field.getSubset(region=instance)
    labels = tuple(subset.componentLabels)
    rows = defaultdict(list)
    for value in subset.values:
        if str(value.position) != 'INTEGRATION_POINT':
            continue
        if value.instance.name != instance.name:
            raise ValueError('unexpected instance in stored field subset')
        point = value.sectionPoint
        key = (point.number, point.description) if point is not None else (None, None)
        vals = tuple(float(x) for x in data(value))
        if not labels or len(vals) != len(labels):
            raise ValueError('native component width mismatch')
        rows[key].append((value.elementLabel, vals))
    if not rows:
        return {'status': 'unavailable', 'reason': 'no native integration-point records'}
    result = []
    for (number, description), values in rows.items():
        ids, vals = zip(*values)
        result.append({'section_point': number, 'description': description,
                       'components': {name: describe(values, ids) for name, values in zip(labels, zip(*vals))}})
    return {'status': 'recorded', 'position': 'INTEGRATION_POINT', 'section_points': result}


def summarize(odb, request):
    instance = odb.rootAssembly.instances[request['instance']]
    step = odb.steps[request['step']]
    index = request['frame']
    frame = step.frames[index]
    if index != len(step.frames)-1 or frame.incrementNumber != request['increment'] or frame.frameValue != 1.:
        raise ValueError('wrong terminal frame')
    cables = {n:i for n,i in odb.rootAssembly.instances.items() if i.elements and 'CABLE' in n.upper()}
    if len(cables) != 5:
        raise ValueError('did not identify five cable instances')
    fields = frame.fieldOutputs

    def recorded(name, inst):
        if name not in fields:
            return {'status': 'unavailable', 'reason': 'field not saved'}
        return components(fields[name], inst)

    return {'terminal': {'step': request['step'], 'index': index, 'increment': frame.incrementNumber,
                         'time': frame.frameValue, 'field_names': sorted(fields.keys())},
            'membrane': {name: recorded(name, instance) for name in ('S','E','LE')},
            'cable_axial_stress': {n: recorded('S', inst) for n, inst in cables.items()},
            'energy': {n: energy_region(r) for n,r in step.historyRegions.items()
                       if 'ALLSD' in r.historyOutputs or 'ALLIE' in r.historyOutputs},
            'limits': 'Stored components; native material/section coordinates; no stress or strain allowables assumed.'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--request', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    q = json.loads(a.request.read_text(encoding='utf-8-sig'))
    before = {k: identity(q[k]['path']) for k in ('input','odb')}
    if any(q[k]['sha256'] != before[k]['sha256'] for k in before):
        raise ValueError('identity mismatch')
    lock = Path(before['odb']['path']).with_suffix('.lck')
    if lock.exists() or a.out.exists():
        raise ValueError('writer present or output already exists')
    from odbAccess import openOdb
    result = with_odb(openOdb, before['odb']['path'], lambda odb: summarize(odb,q))
    after = {k: identity(before[k]['path']) for k in before}
    if before != after or lock.exists():
        raise ValueError('source changed or lock appeared')
    save(a.out, dict(result, before=before, after=after, readOnly=True, handle_closed=True))


if __name__ == '__main__':
    main()
