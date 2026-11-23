"""Move the paired HF clamp attachments in the qualified inline deck dialect."""
import argparse
from decimal import Decimal as D
import difflib
import json
import math
from pathlib import Path
import re

from make_no_clamp_control import cards, require, sha


SIDES = {'CL': (4, 135), 'CR': (5, 45)}
SET_NAMES = ('All_Edges_NoClamp', 'DIAG_EDGES', 'CL_Nodes', 'CR_Nodes')
RADIUS = D('0.21')
BASE = D(20)
BINDS = {'CL': 'clamp_left', 'CR': 'clamp_right'}


def points(xc):
    require(xc.is_finite() and D(0) < xc < D(1), 'position must lie inside the slanted edges')
    h = 10 - 10*xc
    return {'CL': (h, h, D(0)), 'CR': (20-h, h, D(0))}


def reference_nodes(raw):
    block, = re.findall(rb'(?im)^\*part, name=membrane\r?\n\*node\r?\n([^*]+)', raw)
    out = {}
    for line in block.splitlines():
        label, *xyz = line.split(b',')
        label = int(label)
        require(len(xyz) == 3 and label not in out, 'ambiguous membrane node block')
        out[label] = tuple(q.strip().decode('ascii') for q in xyz)
    return out


def select(ref, xc):
    geom = points(xc)
    edge = {n for n, (x,y,_) in ref.items() if y == 0 or x == y or x+y == BASE}
    tolerant_edge = {n for n,(x,y,_) in ref.items() if abs(y)<D('0.0001') or abs(x-y)<D('0.00015') or abs(x+y-BASE)<D('0.00015')}
    require(edge == tolerant_edge, 'exact and recorded edge selector disagree')
    base = {n for n in edge if ref[n][1] == 0}
    patches, centers, margins = {}, {}, {}
    for side, p in geom.items():
        matches = [n for n in edge if ref[n] == p]
        require(len(matches) == 1, 'missing/ambiguous exact radius-zero center on edge: '+side)
        centers[side] = matches[0]
        def d2(x):
            return sum((a-b)**2 for a,b in zip(x,p))
        dist = {n:d2(x) for n,x in ref.items()}
        patches[side] = {n for n,d in dist.items() if d <= RADIUS**2}
        margins[side] = min(abs(math.sqrt(d)-float(RADIUS)) for d in dist.values())
        require(margins[side] > 0, 'a node falls exactly on the sphere boundary; native outcome unresolved')
    require(not patches['CL'] & patches['CR'] and len(set(centers.values())) == 2 and
            {(BASE-ref[n][0], *ref[n][1:]) for n in patches['CL']} == {ref[n] for n in patches['CR']},
            'clamp clouds overlap or are not XY mirrors')
    e = edge - set(centers.values())
    require(base <= e, 'clamp exclusion intersects the base edge')
    return {'points':geom, 'center_labels':centers, 'G':edge, 'base':base, 'sphere_boundary_margin_m':margins,
            'sets':{'All_Edges_NoClamp':e, 'DIAG_EDGES':e-base, **{s+'_Nodes':p for s,p in patches.items()}}}


def one(items, message):
    items = list(items)
    require(len(items) == 1, message)
    return items[0]


def layout(raw):
    lines = raw.splitlines(keepends=True)
    parsed = list(cards(lines))
    def member(name, owner=None):
        expected = {'nset':name.lower()}
        if owner:
            expected['instance'] = owner
        q = one((q for q in parsed if q.keyword == 'nset' and q.scope == 'assembly' and q.options == expected),
                'missing/duplicate/wrong-owner set '+name)
        indices = [i for i,_ in q.rows]
        require(indices == list(range(q.start+1, q.start+1+len(indices))) and indices, 'unsupported member body '+name)
        values = [int(v) for _,row in q.rows for v in row]
        require(len(values) == len(set(values)), 'duplicate member in '+name)
        return (indices[0],indices[-1]+1), set(values)
    groups = {name:member(name,'membrane-1') for name in (*SET_NAMES,'BASE_EDGE','Left_Nodes','Right_Nodes','Top_Nodes')}
    instances, rps, slots = {}, {}, {}
    for side,(label,angle) in SIDES.items():
        inst='inst_cable_'+side.lower()
        q = one(((i,q) for i,q in enumerate(parsed) if q.keyword == 'instance' and q.scope == 'assembly' and
                 q.options.get('name') == inst), 'ambiguous clamp instance '+side)
        j,q=q
        require(q.options == {'name':inst,'part':'cable_'+side.lower()} and len(q.rows)==2 and
                parsed[j+1].keyword=='end instance', 'unsupported clamp rigid placement '+side)
        (ti,t),(ai,axis)=q.rows
        require(len(t)==3 and len(axis)==7, 'unexpected translation/rotation arity '+side)
        instances[side] = tuple(map(D,t)),tuple(map(D,axis))
        slots[side+'_instance'] = (ti,ai)
        require(member(side+'_RP_Set')[1] == member('RP_'+side+'_Set')[1] == {label}, 'wrong tie/coupling RP alias '+side)
        (pi,fields) = one(((i,row) for q in parsed if q.keyword=='node' and q.scope=='assembly' and not q.options
                           for i,row in q.rows if int(row[0])==label), 'missing/ambiguous assembly clamp RP '+side)
        require(len(fields)==4, 'invalid assembly RP row')
        rps[side] = tuple(map(D,fields[1:]));slots[side+'_rp'] = pi
        prefix='cable_'+side.lower()
        require(member(prefix+'_STARTSET',inst)[1]=={1} and member(prefix+'_ENDSET',inst)[1]=={2}, 'wrong clamp endpoint label '+side)
        expected={'constraint name':side.lower()+'_coupling','ref node':side.lower()+'_rp_set','surface':side.lower()+'_nodes_cns_'}
        where=one((j for j,q in enumerate(parsed) if q.keyword=='coupling' and q.scope=='assembly' and q.options==expected),
                  'missing/wrong clamp coupling '+side)
        co=parsed[where+1]
        require(co.keyword=='distributing' and co.options=={'weighting method':'uniform'} and not co.rows,'changed distributing law')
        for surf in (side.lower()+'_nodes','rp_'+side.lower()+'_set',prefix+'_startset'):
            surf_card=one((q for q in parsed if q.keyword=='surface' and q.scope=='assembly' and
                           q.options=={'type':'node','name':surf+'_cns_','internal':''}), 'wrong clamp node surface '+surf)
            require([row for _,row in surf_card.rows]==[(surf,'1.')], 'wrong clamp surface membership '+surf)
        tie=one((q for q in parsed if q.keyword=='tie' and q.scope=='assembly' and
                 q.options=={'name':'tie_'+side.lower(),'adjust':'yes'}), 'wrong clamp tie options')
        require([row for _,row in tie.rows]==[(prefix+'_startset_cns_','rp_'+side.lower()+'_set_cns_')], 'incorrect clamp tie roles')
        outer = prefix+'_endset'
        restraint = [(q.scope,q.options,row[1:]) for q in parsed if q.keyword=='boundary' for _,row in q.rows if row[0]==outer]
        require(sorted((scope,tuple(opts),v) for scope,opts,v in restraint)==[('',(),(str(i),str(i))) for i in (1,2,3)],
                'clamp endpoint restraints changed')
        require(instances[side][1][6]==angle, 'changed clamp cable angle')
    return {'lines':lines, 'sets':groups, 'instances':instances, 'rps':rps, 'slots':slots}


def shift(line, edits):
    fields = line.split(b',')
    for k,delta in edits.items():
        if delta:
            raw=fields[k].strip()
            changed=D(raw.decode('ascii'))+delta
            require(changed==changed.to_integral_value(), 'unsupported noninteger position encoding')
            fields[k]=fields[k].replace(raw, f'{int(changed)}.'.encode('ascii'), 1)
    return b','.join(fields)


def relocate(raw, old, new):
    p = layout(raw)
    lines=p['lines']; edits=[]
    for name in SET_NAMES:
        (a,b),members=p['sets'][name]
        if members != new['sets'][name]:
            ids=sorted(new['sets'][name]); sep=b'\r\n' if lines[a].endswith(b'\r\n') else b'\n'
            rows=[b', '.join(f'{n:5d}'.encode('ascii') for n in ids[i:i+16])+sep for i in range(0,len(ids),16)]
            edits.append((a,b,rows))
    for s in SIDES:
        dx,dy=(new['points'][s][j]-old['points'][s][j] for j in (0,1))
        ti,ai=p['slots'][s+'_instance']; pi=p['slots'][s+'_rp']
        lines[ti]=shift(lines[ti],{0:dx,1:dy})
        lines[ai]=shift(lines[ai],{0:dx,1:dy,3:dx,4:dy})
        lines[pi]=shift(lines[pi],{1:dx,2:dy})
    for a,b,rows in sorted(edits, reverse=True):
        lines[a:b]=rows
    return b''.join(lines)


def protected(p):
    lines=p['lines'][:]
    drops=[p['sets'][n][0] for n in SET_NAMES]
    for k,v in p['slots'].items():
        drops.extend((i,i+1) for i in (v if k.endswith('_instance') else (v,)))
    for a,b in sorted(drops, reverse=True):
        lines[a:b]=[b'<allowed body>\n']
    return b''.join(lines)


def endpoint(translation, axis, local_x):
    # Abaqus applies translation, then rotation about the supplied global line.
    x,y,z=(translation[0]+local_x-axis[0],translation[1]-axis[1],translation[2]-axis[2])
    cosine=math.cos(math.radians(axis[-1])); sine=math.sin(math.radians(axis[-1]))
    return [float(axis[0])+float(x)*cosine-float(y)*sine, float(axis[1])+float(x)*sine+float(y)*cosine, float(axis[2]+z)]


def verify(source,candidate,expected):
    orig, current = layout(source),layout(candidate)
    require(protected(orig)==protected(current), 'unapproved bytes changed outside the eight permitted bodies')
    require(all(current['sets'][k][1]==v for k,v in expected['sets'].items()) and
            current['sets']['BASE_EDGE'][1]==expected['base'], 'stale/wrong/incomplete clamp or edge membership')
    out={}
    for side,(label,angle) in SIDES.items():
        p=expected['points'][side]
        t,a=current['instances'][side]
        require(current['rps'][side]==t==p and a[:3]==p and a[3:6]==(p[0],p[1],p[2]+1) and a[6]==angle,
                'wrong RP/cable translation or nonidentity rotation axis '+side)
        x=endpoint(t,a,D(0)); y=endpoint(t,a,D('0.5'))
        require(x==list(map(float,p)), 'cable start and assembly RP do not coincide '+side)
        out[side]={'RP_label':label,'xy_center_label':expected['center_labels'][side], 'RP_and_start_m':x,
                   'end_m':y,'cable_nodes':[1,2],'axis_angle_deg':angle,'length_m':0.5,
                   'cloud_labels':sorted(expected['sets'][side+'_Nodes'])}
    return out


def bind(original,spec):
    obj={name:Path(spec[name]['path']).read_bytes() for name in ('reference','native','policy','config')}
    for name,b in obj.items():
        require(sha(b)==spec[name]['sha256'] and len(b)==spec[name]['bytes'], 'input pin mismatch: '+name)
    cfg=json.loads(obj['config']); nat=json.loads(obj['native']); init=json.loads(obj['policy'])
    require(cfg['PATCH_R_CLAMP']==float(RADIUS) and cfg['CLAMP_EXCL_R']==0 and cfg['BASE']==20 and cfg['HEIGHT']==10 and
            cfg['CABLE_PRESTRESS']==17_100_000 and cfg['LOAD_MODE']=='cable_stress' and not cfg['RELEASE_EDGE_U3'], 'frozen case rule mismatch')
    ref={n:tuple(map(D,xyz)) for n,xyz in reference_nodes(obj['reference']).items()}
    seeded=reference_nodes(original)
    require(len(seeded)==len(ref)==30501 and seeded.keys()==ref.keys() and all(z==0 for _,_,z in ref.values()), 'unseeded mesh coverage differs')
    old=select(ref,D('.7')); recorded=nat['unseeded']; masks=init['labels']
    require(old['G']==set(nat['intended_edge'])==set(masks['G']) and old['sets']['All_Edges_NoClamp']==set(recorded['edges'])==set(masks['E']) and
            set(old['center_labels'].values())==set(recorded['excluded_edges'])==set(masks['X']) and
            old['base']==set(masks['base']) and old['sets']['DIAG_EDGES']==set(masks['diagonal']), 'recomputed baseline edge cells disagree')
    p=layout(original)
    for name,key in [('Left_Nodes','left'),('Right_Nodes','right'),('Top_Nodes','top')]:
        require(p['sets'][name][1]==set(recorded['clouds'][key]), 'corner membership mismatch')
    for side,src in BINDS.items():
        require(old['sets'][side+'_Nodes']==set(recorded['clouds'][src])==set(masks['clouds'][src]), 'recomputed baseline sphere disagrees')
    require(relocate(original,old,old)==original, 'identity relocation changed bytes')
    verify(original,original,old)
    return ref,seeded,old


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--spec', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args=ap.parse_args(); spec=json.loads(args.spec.read_bytes()); baseline=Path(spec['source']['path']).read_bytes()
    require(sha(baseline)==spec['source']['sha256'] and len(baseline)==spec['source']['bytes'], 'submitted deck bytes differ')
    ref,seeded,old=bind(baseline,spec); entries={}
    prior=old['center_labels']; old_x=set(prior.values())
    for design in spec['designs']:
        case=design['id'];xc=D(design['x_c']); require(re.fullmatch('X[0-9]{3}',case) and xc!=D('.7'), 'unexpected design identity')
        new=select(ref,xc); x=set(new['center_labels'].values()); members=new['sets']['All_Edges_NoClamp']
        require(len(new['G'])==800 and len(members)==798 and old_x<=members and not x & members, 'wrong restored or excluded center scope')
        candidate=relocate(baseline,old,new); links=verify(baseline,candidate,new)
        require(reference_nodes(candidate)==seeded, 'membrane coordinate tokens changed')
        root=args.out/case;root.mkdir(); inp=root/(case+'.inp'); diff=root/'input.patch'
        with inp.open('xb') as f: f.write(candidate)
        with diff.open('xb') as f:
            f.writelines(difflib.diff_bytes(difflib.unified_diff,baseline.splitlines(keepends=True),candidate.splitlines(keepends=True),b'DIRECT01',inp.name.encode('ascii')))
        entries[case]={'x_c':str(xc),'input_path':str(inp),'sha256':sha(candidate),'attachment':links,'edge_count':len(new['G']),
                       'direct_edge_restraint_count':len(members),'restored_center_labels':sorted(old_x),'excluded_center_labels':sorted(x),
                       'seed_Z_tokens_at_old_and_new_centers':{n:seeded[n][2] for n in sorted(old_x|x)},
                       'sphere_boundary_margin_m':new['sphere_boundary_margin_m']}
    with (args.out/'attachment-ledger.json').open('x',encoding='utf-8') as f:
        json.dump({'identity_0p7_byte_exact':True,'identity_sha256':sha(baseline),'cases':entries},f,indent=2); f.write('\n')
    print('Identity and position-only dependency gates passed.')


if __name__=='__main__':
    main()
