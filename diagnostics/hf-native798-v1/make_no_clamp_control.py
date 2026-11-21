"""Derive a hardware-free control from this diagnostic ABAQUS deck family only."""

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re


KEYWORDS = set('heading;preprint;part;node;element;nset;elset;solid section;end part;shell section;'
               'assembly;instance;end instance;surface;coupling;distributing;tie;end assembly;material;'
               'density;elastic;boundary;initial conditions;step;static;controls;restart;output;'
               'node output;element output;end step'.split(';'))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


@dataclass(frozen=True)
class Card:
    start: int
    stop: int
    keyword: str
    options: dict
    rows: tuple
    scope: str


def cards(lines):
    starts = [i for i, line in enumerate(lines) if line.startswith(b'*') and not line.startswith(b'**')]
    scope = ''
    for start, stop in zip(starts, starts[1:] + [len(lines)]):
        head, *args = lines[start].decode('ascii').strip().lower()[1:].split(',')
        require(head in KEYWORDS, f'unsupported keyword: {head}')
        opt = {}
        for arg in args:
            k, _, v = arg.strip().partition('=')
            require(k and k not in opt, f'ambiguous card at line {start+1}')
            opt[k] = v.strip()
        rows = tuple((i, tuple(x.strip().lower() for x in lines[i].decode('ascii').strip().rstrip(',').split(',')))
                     for i in range(start+1, stop) if lines[i].strip() and not lines[i].startswith(b'**'))
        yield Card(start, stop, head, opt, rows, scope)
        if head in ('part','assembly','instance','step'):
            require(not scope or (head == 'instance' and scope == 'assembly'), 'unsupported scope nesting')
            scope = head
        elif head in ('end part','end assembly','end instance','end step'):
            require(scope == head[4:], 'unbalanced deck scope')
            scope = 'assembly' if scope == 'instance' else ''
    require(not scope, 'unterminated deck scope')


def transform(raw, membership):
    lines = raw.splitlines(keepends=True)
    parsed = list(cards(lines))
    native = membership['unseeded']
    edges, extra, full = (set(seq) for seq in (native['edges'], native['excluded_edges'], membership['intended_edge']))
    require(len(extra) == 2 and not edges & extra and edges | extra == full, 'invalid frozen edge partition')
    removed, replaced = {}, {}

    def select(card, tokens, *, scope=None):
        result = [i for i, c in enumerate(parsed) if c.keyword == card and
                  all(c.options.get(k) == v for k, v in tokens.items()) and (scope is None or c.scope == scope)]
        require(len(result) == 1, f'expected one {card} {tokens}')
        return result[0]

    def delete(i, tag, *, remark=None, until=None):
        first, last = parsed[i].start, parsed[i if until is None else until]
        stop = last.rows[-1][0]+1 if last.rows else last.start+1
        if remark:
            require(lines[first-1].strip().lower().startswith(remark.encode()), f'missing owning comment: {tag}')
            first -= 1
        require(not set(range(first, stop)) & removed.keys(), f'overlapping deletion: {tag}')
        removed.update((n, tag) for n in range(first, stop))

    def fields(i):
        return [value for _, row in parsed[i].rows for value in row]

    dead = set()
    for short, side in (('cl','clamp_left'), ('cr','clamp_right')):
        part, inst, inner, outer = f'cable_{short}', f'inst_cable_{short}', f'cable_{short}_startset', f'cable_{short}_endset'
        rp, alias, cloud, couple, tie = f'rp_{short}_set', f'{short}_rp_set', f'{short}_nodes', f'{short}_coupling', f'tie_{short}'
        for block, name, opts in (('part',part,{}), ('instance',inst,{'part':part})):
            i = select(block, {'name':name, **opts})
            stop = next(k for k in range(i+1,len(parsed)) if parsed[k].keyword == 'end '+block)
            delete(i, name, until=stop)
        for name, target, expected in ((inner,inst,{1}), (outer,inst,{2}), (cloud,'membrane-1',set(native['clouds'][side]))):
            i = select('nset', {'nset':name,'instance':target}, scope='assembly')
            require(set(map(int, fields(i))) == expected, f'membership differs: {name}')
            delete(i,name)
        rp_labels = []
        for name in (rp, alias):
            i = select('nset', {'nset':name}, scope='assembly')
            values = fields(i)
            require(len(values) == 1 and values[0].isdigit() and not parsed[i].options.get('instance'), 'not an RP singleton')
            rp_labels.append(values[0]); delete(i,name)
        require(len(set(rp_labels)) == 1, 'RP alias differs')
        definitions = [i for i,c in enumerate(parsed) if c.keyword == 'node' and c.scope == 'assembly' and
                       any(row[0] == rp_labels[0] for _,row in c.rows)]
        require(len(definitions) == 1 and len(parsed[definitions[0]].rows) == 1, 'RP is not an exclusive definition')
        for c in parsed:
            if c.keyword == 'nset' and c.scope == 'assembly' and 'instance' not in c.options and c.options['nset'] not in {rp, alias}:
                require(all(rp_labels[0] not in row for _,row in c.rows), 'another assembly set shares the clamp RP')
        delete(definitions[0],rp+' node')
        for name in (cloud,rp,inner):
            surf = name+'_cns_'
            i = select('surface',{'name':surf},scope='assembly')
            require(fields(i) == [name,'1.'], f'unexpected surface: {surf}')
            delete(i,surf)
            dead.add(surf)
        i = select('coupling', {'constraint name':couple,'ref node':alias,'surface':cloud+'_cns_'},scope='assembly')
        require(parsed[i+1].keyword == 'distributing' and not parsed[i+1].rows, 'unexpected clamp coupling law')
        delete(i,couple,remark='** constraint: '+couple,until=i+1)
        i = select('tie',{'name':tie},scope='assembly')
        require(fields(i) == [inner+'_cns_',rp+'_cns_'], 'tie ownership differs')
        delete(i,tie,remark='** constraint: '+tie)
        bc = [i for i,c in enumerate(parsed) if c.keyword == 'boundary' and any(row[0]==outer for _,row in c.rows)]
        require(len(bc) == 2 and sorted(tuple(row[1:] for _,row in parsed[i].rows) for i in bc) ==
                [(('1','1'),('2','2')), (('3','3'),)] and all(not parsed[i].scope and all(row[0]==outer for _,row in parsed[i].rows) for i in bc),
                'unsupported clamp BC ownership')
        for i in bc:
            delete(i,outer+' boundary',remark='** name:')
        stress = [i for i,c in enumerate(parsed) if c.keyword == 'initial conditions' and
                  any(row[0]==inst+'.all_elem' for _,row in c.rows)]
        require(len(stress)==1 and len(parsed[stress[0]].rows)==1 and
                parsed[stress[0]].options=={'type':'stress'} and not parsed[stress[0]].scope, 'unexpected clamp initial stress')
        delete(stress[0],inst+' stress',remark='** name:')
        dead.update((part,inst,inner,outer,rp,alias,cloud,couple,tie))

    bas = select('nset', {'nset':'base_edge','instance':'membrane-1'}, scope='assembly')
    diag = select('nset', {'nset':'diag_edges','instance':'membrane-1'}, scope='assembly')
    edge = select('nset', {'nset':'all_edges_noclamp','instance':'membrane-1'}, scope='assembly')
    base = set(map(int,fields(bas)))
    diagonal = set(map(int,fields(diag)))
    require(set(map(int,fields(edge))) == edges and base | diagonal == edges and not base & diagonal, 'unexpected edge groups')
    require(sum(len(fields(i)) for i in (bas,diag)) == len(edges), 'duplicate edge labels')
    require(sum(c.keyword=='boundary' and not c.scope and c.rows and c.rows[0][1]==('all_edges_noclamp','3','3')
                for c in parsed) == 1, 'missing homogeneous Initial edge U3')
    for i, values in ((edge,full),(diag,diagonal|extra)):
        first, last = parsed[i].rows[0][0], parsed[i].rows[-1][0]
        require(tuple(n for n,_ in parsed[i].rows) == tuple(range(first,last+1)), 'noncontiguous edge set data')
        newline = b'\r\n' if lines[first].endswith(b'\r\n') else b'\n'
        seq = sorted(values)
        replaced[first] = (last, b''.join((', '.join(map(str,seq[k:k+16]))).encode()+newline for k in range(0,len(seq),16)))
    for c in parsed:
        if c.start not in removed:
            tokens = set(re.findall(r'[a-z_][a-z_0-9-]*', ' '.join(c.options.values()) + ' ' + ' '.join(fields_for(c))))
            require(not dead & tokens, f'unresolved dependent definition at line {c.start+1}: {sorted(dead & tokens)}')
    out = []
    i = 0
    while i < len(lines):
        if i in replaced:
            last, content = replaced[i]; out.append(content); i = last+1
        else:
            if i not in removed:
                out.append(lines[i])
            i += 1
    return b''.join(out), {'removed_source_lines': [[i+1, removed[i]] for i in sorted(removed)],
                           'replaced_source_ranges': [[i+1, v[0]+1] for i,v in sorted(replaced.items())],
                           'restored_u3': sorted(extra), 'before_edge_count':len(edges), 'after_edge_count':len(full)}


def fields_for(card):
    return (value for _, row in card.rows for value in row)


def verify_intervention(source, candidate, deleted_ranges, changed_sets):
    """Independent byte comparison against reviewed source lines, not the edit plan emitted by transform."""
    src = source.splitlines(keepends=True)
    got = candidate.splitlines(keepends=True)
    drops = {i for a,b in deleted_ranges for i in range(a-1,b)}
    masks = {a-1:(b, members) for a,b,members in changed_sets}
    i = j = 0
    while i < len(src):
        if i in drops:
            i += 1
        elif i in masks:
            end, expected = masks[i]
            ids = []
            while j < len(got) and not got[j].startswith(b'*'):
                row = got[j].strip().rstrip(b',').split(b',')
                require(all(v.strip().isdigit() for v in row), 'not a bounded explicit edge set')
                ids.extend(int(v) for v in row)
                j += 1
            require(sorted(ids) == sorted(expected), f'wrong edge membership near source line {i+1}')
            i = end
        else:
            require(j < len(got) and src[i] == got[j], f'unapproved edit near source line {i+1}')
            i += 1
            j += 1
    require(j == len(got), 'unexpected trailing lines')




def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('source','membership','out','receipt'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--source-sha256',required=True)
    p.add_argument('--membership-sha256',required=True)
    a = p.parse_args()
    source, member = a.source.read_bytes(), a.membership.read_bytes()
    require(sha(source)==a.source_sha256 and sha(member)==a.membership_sha256, 'pinned source identity differs')
    require(not a.out.exists() and not a.receipt.exists() and a.out != a.receipt, 'outputs must be new')
    candidate, evidence = transform(source, json.loads(member))
    with a.out.open('xb') as out:
        out.write(candidate)
    with a.receipt.open('x',encoding='utf-8') as receipt:
        json.dump(dict(evidence,source_sha256=sha(source),membership_sha256=sha(member),candidate_sha256=sha(candidate)),receipt,indent=2)
        receipt.write('\n')


if __name__ == '__main__':
    main()
