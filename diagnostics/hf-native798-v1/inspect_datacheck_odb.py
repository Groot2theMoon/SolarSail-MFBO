"""Inspect a closed Abaqus/Standard datacheck ODB, without modifying its input."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import re
import sys


CABLES = ('TOP', 'LEFT', 'RIGHT', 'CL', 'CR')
ASSEMBLY = 'ASSEMBLY'
WARNING = re.compile(r'Boundary conditions are specified on inactive dof ([456]) at (\d+) nodes?\.\s*'
                     r'The nodes have been identified in node set ([\w.-]+)\.', re.I)
PAIR = re.compile(r'FOR \*TIE PAIR\s+\((\S+?)-(\S+?)\),\s+'
                  r'ADJUSTMENT WAS SPECIFIED BUT NO NODE WAS ADJUSTED MORE THAN THE\s+'
                  r'ADJUSTMENT DISTANCE\s*=\s*([\d.E+-]+)\.', re.I)


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def digest(path):
    before = path.stat()
    with path.open('rb') as stream:
        sha = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            sha.update(block)
    after = path.stat()
    require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns),
            f'file changed while hashing: {path.name}')
    return dict(path=str(path), bytes=after.st_size, mtime_ns=after.st_mtime_ns, sha256=sha.hexdigest())


def read_text(path):
    return path.read_text(encoding='utf-8-sig', errors='replace')


def attribution(odb, deck, want):
    root = odb.parent
    controls = [root/'run-records'/n for n in ('control-outcome.json', 'datacheck.launch.json', 'datacheck.completion.json')]
    control, launch, done = [json.loads(read_text(p)) for p in controls]
    job = control['job']
    data = root/(job+'.dat')
    text = read_text(data)
    args = f'job={job} input={deck.name} datacheck cpus=4 interactive'
    require(deck.parent == root == Path(control['cwd']).resolve() == Path(launch['cwd']).resolve()
            and odb.name == job+'.odb' and control['args'] == launch['args'] == args
            and control['attempts'] == 1 and control['closed'] is True and control['timeout'] is False
            and control['remaining'] == [] and done['status'] == 'launcher_exited'
            and control['input_sha256_before'] == control['input_sha256_after'] == want
            and re.search(r'^\s*ANALYSIS DATACHECK COMPLETE\s*$', text, re.M)
            and not (root/(job+'.lck')).exists(), 'datacheck attribution or closure does not match')
    return job, [digest(p) for p in [*controls, data]], text


def deck_targets(text):
    """Read only the five authored singleton cable/RP chains, not general input syntax."""
    cables = {f'CABLE_{k}_{end}SET' for k in CABLES for end in ('START', 'END')}
    rps = {f'RP_{k}_SET' for k in CABLES}
    sets, surfaces, ties = {}, {}, {}
    lines = text.splitlines()
    for i, line in enumerate(lines):
        key, _, tail = line.strip().upper().partition(',')
        if key not in ('*NSET', '*SURFACE', '*TIE'):
            continue
        opts = dict(x.strip().split('=', 1) for x in tail.split(',') if '=' in x)
        name = opts.get('NSET' if key == '*NSET' else 'NAME')
        if name not in (cables | rps if key == '*NSET' else
                        {f'{n}_CNS_' for n in cables | rps} if key == '*SURFACE' else
                        {f'TIE_{k}' for k in CABLES}):
            continue
        values = tuple(x.strip().upper() for x in lines[i+1].split(',') if x.strip())
        require(i+2 == len(lines) or lines[i+2].lstrip().startswith('*'), f'non-single-row input definition: {name}')
        if key == '*NSET':
            require(name not in sets and len(values) == 1 and values[0].isdigit(), f'ambiguous singleton {name}')
            namespace = opts.get('INSTANCE', ASSEMBLY)
            require((namespace == ASSEMBLY) == (name in rps), f'wrong input namespace for {name}')
            sets[name] = (namespace, int(values[0]))
        elif key == '*SURFACE':
            require(name not in surfaces and opts.get('TYPE') == 'NODE' and len(values) == 2,
                    f'ambiguous node surface {name}')
            surfaces[name] = values[0]
        else:
            require(name not in ties and len(values) == 2, f'ambiguous tie {name}')
            ties[name] = (values, i+1)
    result = {}
    for k in CABLES:
        inn, out, rp = f'CABLE_{k}_STARTSET', f'CABLE_{k}_ENDSET', f'RP_{k}_SET'
        pair, line = ties['TIE_'+k]
        require(tuple(surfaces[s] for s in pair) == (inn, rp), f'wrong secondary/main roles for TIE_{k}')
        require(sets[inn][0] == sets[out][0] and sets[inn] != sets[out], f'wrong cable endpoints for {k}')
        result[k] = dict(inner=sets[inn], outer=sets[out], rp=sets[rp], secondary_surface=pair[0], main_surface=pair[1], input_line=line)
    require(len({r['inner'] for r in result.values()}) == len(CABLES), 'tie status cannot be attributed to unique secondary nodes')
    return result


def namespace(name):
    return ASSEMBLY if name is None or name == '' else str(name).upper()


def set_nodes(region, owner=None):
    """Documented OdbSet node arrays: flat for one namespace, grouped for many."""
    nodes, names = region.nodes, region.instanceNames
    if len(nodes) == 0:
        return set()
    scopes = (owner,) if owner is not None else tuple(names or ())
    groups = (nodes,) if hasattr(nodes[0], 'label') else nodes
    require(len(scopes) == len(groups), 'node arrays lack unambiguous instance namespaces')
    found = set()
    for scope, group in zip(scopes, groups):
        for node in group:
            found.add((namespace(scope), int(node.label)))
    return found


def sets_inventory(assembly):
    scoped = [(None, assembly)] + [(n, assembly.instances[n]) for n in assembly.instances.keys()]
    repositories = []
    rows = []
    for owner, part in scoped:
        keys = sorted(part.nodeSets.keys())
        surfaces = sorted(part.surfaces.keys())
        repositories.append(dict(owner=owner or ASSEMBLY, node_sets=keys, surfaces=surfaces,
                                 element_sets=sorted(part.elementSets.keys())))
        rows.extend((owner, name, part.nodeSets[name]) for name in keys)
    return repositories, rows


def set_record(owner, name, region):
    row = dict(scope=owner or ASSEMBLY, name=name, locator=f'{owner or "rootAssembly"}.nodeSets[{name!r}]',
               internal=bool(region.isInternal))
    try:
        row['nodes'] = sorted(set_nodes(region, owner))
    except ValueError as e:
        row['nodes'] = None
        row['reason'] = str(e)
    return row


def incident_types(assembly, targets):
    types = {n: set() for n in targets}
    for owner in [None, *assembly.instances.keys()]:
        part = assembly if owner is None else assembly.instances[owner]
        for element in part.elements:
            labels = element.connectivity
            for scope, label in zip(element.instanceNames if owner is None else (owner,)*len(labels), labels):
                key = (namespace(scope), int(label))
                if key in types:
                    types[key].add(element.type)
    return types


def warning_rows(text, rows, expected, assembly):
    warnings = {}
    selected = set()
    for m in WARNING.finditer(text):
        dof, count, name = m.groups()
        require(dof not in warnings, f'duplicate inactive-DOF warning: {dof}')
        sets = [set_record(o,n,s) for o,n,s in rows if n.upper() == name.upper()]
        actual = (set().union(*(set(x['nodes']) for x in sets))
                  if sets and all(x['nodes'] is not None for x in sets) else None)
        warnings[dof] = dict(dof=int(dof), reported_count=int(count), warning_set=name, source_sets=sets,
                             dat_line=text.count('\n', 0, m.start())+1, actual=actual)
        if actual is not None:
            selected.update(actual)
    types = incident_types(assembly, selected | expected)
    for dof in '456':
        row = warnings.get(dof, dict(dof=int(dof), source_sets=[], actual=None))
        actual = row.pop('actual')
        row.update(status=('not_exposed_or_unresolved' if actual is None else
                           'match' if actual == expected and len(actual) == row['reported_count']
                           and all(types[n] == {'T3D2'} for n in actual) else 'mismatch'),
                   expected=sorted(expected), unexpected=None if actual is None else sorted(actual - expected),
                   missing=None if actual is None else sorted(expected - actual),
                   nodes=None if actual is None else [dict(node=n, incident_element_types=sorted(types[n])) for n in sorted(actual)])
        warnings[dof] = row
    return list(warnings.values())


def tie_marker(name):
    words = re.sub(r'(?<=[a-z0-9])(?=[A-Z])', ' ', name).upper()
    words = re.sub(r'\b(?:UN[ _-]?TIED|NOT[ _-]?TIED)\b', 'UNTIED', words.replace('_', ' '))
    tokens = set(re.findall('[A-Z]+', words)) & {'TIED', 'UNTIED'}
    return next(iter(tokens)) if len(tokens) == 1 else None


def ties(assembly, text, rows, targets):
    evidence = [dict(set_record(o,n,s), explicit_status=tie_marker(n)) for o,n,s in rows
                if s.isInternal and re.search(r'TIE',n,re.I)]
    pairs = {(m[1].upper(), m[2].upper()): dict(threshold=m[3], dat_line=text.count('\n',0,m.start())+1,
             finding='specified adjustment; no node adjusted more than the reported distance') for m in PAIR.finditer(text)}
    result = []
    keys = set(assembly.surfaces.keys())
    for k, r in targets.items():
        corroboration = {}
        for role, node in [('secondary',r['inner']), ('main',r['rp'])]:
            name = r[role+'_surface']
            matches = [n for n in keys if n.upper() in (name, ASSEMBLY+'_'+name)]
            try:
                observed = sorted(set_nodes(assembly.surfaces[matches[0]])) if len(matches) == 1 else None
            except ValueError:
                observed = None
            corroboration[role] = dict(names=matches,nodes=observed,match=None if observed is None else observed == [node])
        hits = [x for x in evidence if x['explicit_status'] and x['nodes'] is not None and r['inner'] in x['nodes']]
        statuses = {x['explicit_status'] for x in hits}
        status = ('confirmed_'+next(iter(statuses)).lower() if len(statuses) == 1
                  and not any(x['match'] is False for x in corroboration.values()) else 'not_exposed_or_unresolved')
        adjustment = pairs.get(tuple(ASSEMBLY+'_'+r[role+'_surface'] for role in ('secondary','main')))
        result.append(dict(tie='TIE_'+k, **r, coverage=status, status_sets=hits, odb_surfaces=corroboration,
                           adjustment_text=adjustment, effective_position_tolerance=None))
    return result, evidence


def inspect(odb, deck_text, dat):
    assembly = odb.rootAssembly
    repositories, rows = sets_inventory(assembly)
    frames, adjustment = {}, []
    targets = deck_targets(deck_text)
    inner = {r['inner'] for r in targets.values()}
    for name in odb.steps.keys():
        steps = odb.steps[name]
        frames[name] = []
        for i, frame in enumerate(steps.frames):
            keys = sorted(frame.fieldOutputs.keys())
            frames[name].append(dict(index=i, increment=frame.incrementNumber, value=frame.frameValue, keys=keys))
            for key in ('TIEADJUST', 'TIEDSTATUS'):
                if key not in keys:
                    continue
                for x in frame.fieldOutputs[key].values:
                    if str(x.position) != 'NODAL' or x.instance is None:
                        continue
                    node = (namespace(x.instance.name), int(x.nodeLabel))
                    if node not in inner:
                        continue
                    data = x.dataDouble if str(x.precision) == 'DOUBLE_PRECISION' else x.data
                    data = tuple(data) if hasattr(data, '__iter__') else (data,)
                    require(all(math.isfinite(v) for v in data), f'nonfinite native {key}')
                    adjustment.append(dict(node=node, step=name, frame=i, field=key, precision=str(x.precision), data=data))
    coverage, evidence = ties(assembly, dat, rows, targets)
    return dict(repositories=repositories, steps=frames,
                inactive_dofs=warning_rows(dat, rows, {targets[k]['outer'] for k in CABLES[:3]}, assembly),
                ties=coverage, internal_tie_candidates=evidence, adjustment_fields=adjustment,
                dat_tie_table_markers=[dict(line=i,text=line.strip()) for i,line in enumerate(dat.splitlines(),1)
                                       if re.search(r'TIED\s+TO\s+(?:THE\s+)?(?:MAIN|MASTER)|TIED\s+NODES|UNTIED\s+NODES|SECONDARY\s+NODE\s+',line,re.I)],
                limits=['no position tolerance or exact zero adjustment inferred from missing output',
                        'Explicit-only fields are optional in this Standard datacheck; model coordinates are not adjustment evidence'])


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    for flag in ('odb', 'input', 'out-dir'):
        cli.add_argument('--'+flag, required=True, type=Path)
    cli.add_argument('--expected-input-sha256', required=True)
    opt = cli.parse_args()
    require(not opt.out_dir.exists() and not opt.out_dir.is_symlink(), 'out-dir must be new')
    odb_file, deck = opt.odb.resolve(strict=True), opt.input.resolve(strict=True)
    want = opt.expected_input_sha256
    require(re.fullmatch('[0-9a-f]{64}',want) and not opt.odb.is_symlink() and not opt.input.is_symlink()
            and odb_file.suffix.lower() == '.odb', 'invalid input identity or path')
    job, records, dat = attribution(odb_file, deck, want)
    before = dict(odb=digest(odb_file), input=digest(deck))
    require(before['odb']['bytes'] > 0 and before['input']['sha256'] == want, 'empty ODB or unexpected input digest')
    opt.out_dir.mkdir()
    report = dict(job=job, before=before, source_records=records, helper=digest(Path(__file__).resolve()),
                  runtime=dict(python=sys.version, executable=sys.executable, host=platform.node(), native_header=dat.splitlines()[3:9]),
                  read_only=True)
    try:
        import odbAccess
        odb = odbAccess.openOdb(path=str(odb_file), readOnly=True)
        try:
            report['runtime']['odb_job_data'] = {k: str(getattr(odb.jobData,k,None)) for k in ('version','analysisCode','precision')}
            report['result'] = inspect(odb, read_text(deck), dat)
        finally:
            odb.close()
    finally:
        report['after'] = dict(odb=digest(odb_file), input=digest(deck))
        report['identities_unchanged'] = before == report['after']
        with (opt.out_dir/'inspection.json').open('x', encoding='utf-8', newline='\n') as out:
            json.dump(report,out,indent=2,allow_nan=False)
            out.write('\n')
    require(report['identities_unchanged'], 'input or ODB changed during inspection')
    print('Read-only datacheck ODB inspection complete.')


if __name__ == '__main__':
    main()
