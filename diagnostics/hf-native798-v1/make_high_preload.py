"""Derive and independently check one higher-preload input from its frozen baseline."""

import argparse
import difflib
import hashlib
import json
from pathlib import Path
import re

from make_no_clamp_control import cards, require
from verify_sail_geometry_roundtrip import membrane

CABLES = {
    'inst_cable_cl': ('cable_cl', '171000.'),
    'inst_cable_cr': ('cable_cr', '171000.'),
    'inst_cable_left': ('cable_botl', '1.71e+07'),
    'inst_cable_right': ('cable_botr', '1.71e+07'),
    'inst_cable_top': ('cable_top', '1.71e+07'),
}
CLAMPS = ('inst_cable_cl', 'inst_cable_cr')
DEFAULT_CLAMP_PA = 8_550_000
CORNER_PA = 17_100_000


def digest(data):
    return hashlib.sha256(data).hexdigest()


def _clamp_token(clamp_pa):
    require(type(clamp_pa) is int and clamp_pa > 0 and clamp_pa != 171_000,
            'replacement clamp stress must be a positive changed integer in Pa')
    return f'{clamp_pa}.'.encode('ascii')


def derive(source, *, clamp_pa=DEFAULT_CLAMP_PA):
    replacement = _clamp_token(clamp_pa)
    lines = source.splitlines(keepends=True)
    parsed = list(cards(lines))
    parts, instances = {}, {}
    for i, card in enumerate(parsed):
        if card.keyword == 'part':
            name = card.options['name']
            require(name not in parts, f'duplicate part: {name}')
            end = next(j for j in range(i + 1, len(parsed)) if parsed[j].keyword == 'end part')
            parts[name] = [p for p in parsed[i + 1:end] if p.keyword == 'elset' and p.options.get('elset') == 'all_elem']
        elif card.keyword == 'instance' and card.options['name'] in CABLES:
            instance = card.options['name']
            require(instance not in instances, f'duplicate instance: {instance}')
            instances[instance] = card.options['part']
    require(instances == {k: p for k, (p, _) in CABLES.items()}, 'cable instance/part binding differs')
    require(all(len(parts.get(part, ())) == 1 and parts[part][0].rows for part in instances.values()),
            'missing, duplicate or empty cable All_Elem set')

    assignments = set()
    for card in parsed:
        if card.keyword != 'initial conditions' or card.options.get('type') != 'stress':
            continue
        require(not card.scope and card.options == {'type': 'stress'}, 'unexpected initial stress scope/options')
        for i, (target, *values) in card.rows:
            require(target not in assignments, f'duplicate stress assignment: {target}')
            assignments.add(target)
            if target == 'membrane-1.all_elem':
                require(values == ['500.', '500.', '0.'], 'unexpected membrane stress')
                continue
            instance, _, elset = target.partition('.')
            require(instance in CABLES and elset == 'all_elem', f'unexpected stress target: {target}')
            require(values == [CABLES[instance][1], '0.', '0.'], f'unexpected stress value: {target}')
            if instance in CLAMPS:
                fields = lines[i].split(b',')
                require(fields[1] == b' 171000.', f'unexpected clamp token: {target}')
                fields[1] = b' ' + replacement
                lines[i] = b','.join(fields)
    require(assignments == {'membrane-1.all_elem', *(x + '.all_elem' for x in CABLES)}, 'missing initial stress assignment')
    return b''.join(lines)


def verify_delta(source, candidate, *, clamp_pa=DEFAULT_CLAMP_PA):
    """Compare every byte outside two independently specified authored rows."""
    replacement = _clamp_token(clamp_pa)
    old = source.splitlines(keepends=True)
    new = candidate.splitlines(keepends=True)
    require(len(old) == len(new), 'changed line count')
    pattern = re.compile(rb'(Inst_cable_C[LR]\.All_Elem, )171000\.(, 0\., 0\.,[ \t]*(?:\r?\n)?)')
    seen = {}
    for i, (a, b) in enumerate(zip(old, new), 1):
        if a == b:
            continue
        row = pattern.fullmatch(a)
        require(row is not None and b == row[1] + replacement + row[2], f'unapproved byte change at line {i}')
        name = row[1].decode('ascii').split('.')[0]
        require(name not in seen and i > 1 and old[i - 2].strip().lower() == b'*initial conditions, type=stress',
                f'duplicate or out-of-card stress edit at line {i}')
        seen[name] = i
    require(seen.keys() == {'Inst_cable_CL', 'Inst_cable_CR'}, 'expected two and only two clamp stress edits')
    return seen


def raw_node_rows(data):
    """Capture the submitted membrane node tokens, including signed/tiny Z."""
    rows = re.findall(rb'(?mi)^\*part, name=membrane\r?\n\*node\r?\n([^*]+)', data)
    require(len(rows) == 1 and rows[0], 'missing or ambiguous membrane node card')
    return rows[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'control', 'out', 'diff', 'receipt'):
        parser.add_argument('--' + name, type=Path, required=True)
    for name in ('source', 'control'):
        parser.add_argument('--' + name + '-sha256', required=True)
    parser.add_argument('--clamp-pa', type=int, default=DEFAULT_CLAMP_PA)
    a = parser.parse_args()
    outputs = (a.out, a.diff, a.receipt)
    require(len({p.resolve() for p in outputs}) == len(outputs) and not any(p.exists() for p in outputs),
            'all outputs must be fresh and distinct')
    source, control = a.source.read_bytes(), a.control.read_bytes()
    require(digest(source) == a.source_sha256 and digest(control) == a.control_sha256, 'pinned input digest differs')
    candidate = derive(source, clamp_pa=a.clamp_pa)
    changes = verify_delta(source, candidate, clamp_pa=a.clamp_pa)
    baseline_nodes, baseline_elements = membrane(a.source, 'MEMBRANE-1')
    control_nodes, control_elements = membrane(a.control, 'MEMBRANE-1')
    nodes = raw_node_rows(source)
    require(len(baseline_nodes) == 30501 and baseline_nodes == control_nodes and baseline_elements == control_elements and
            nodes == raw_node_rows(candidate) == raw_node_rows(control), 'shared initial membrane field or topology differs')
    with a.out.open('xb') as out:
        out.write(candidate)
    with a.diff.open('xb') as out:
        out.writelines(difflib.diff_bytes(difflib.unified_diff, source.splitlines(keepends=True),
                      candidate.splitlines(keepends=True), b'DIRECT01', f'clamp-{a.clamp_pa}-Pa'.encode('ascii')))
    with a.receipt.open('x', encoding='utf-8') as out:
        json.dump({'source_sha256': digest(source), 'control_sha256': digest(control), 'candidate_sha256': digest(candidate),
                   'candidate_bytes': len(candidate), 'changed_source_lines': changes, 'other_bytes_unchanged': True,
                   'shared_applied_membrane': {'node_count': len(baseline_nodes), 'element_count': len(baseline_elements),
                                               'ordered_node_rows_sha256': digest(nodes), 'byte_identical': True},
                   'x_c': 0.7, 'd_c_lower': 0.01, 'd_c_higher': a.clamp_pa / CORNER_PA,
                   'axial_stress_pa': {'clamp_original': 171000, 'clamp_applied': a.clamp_pa, 'corner_retained': CORNER_PA}},
                  out, indent=2)
        out.write('\n')


if __name__ == '__main__':
    main()
