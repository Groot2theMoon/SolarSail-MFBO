"""Derive one shared interior-normal field, never launch a solver."""

import argparse
import csv
from decimal import Decimal
import difflib
import hashlib
import json
import math
from pathlib import Path
import re
import struct

from make_no_clamp_control import require, verify_intervention
from verify_sail_geometry_roundtrip import membrane

FAMILIES = ('paper','clamp_only')
CENTRES = {5161,5301}
MEMBRANE = re.compile(rb'(?mi)^\*part, name=membrane\r?\n\*node\r?\n([^*]+)')
FLOAT = rb'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
ROW = re.compile(rb'([ \t]*(\d+)[ \t]*,[ \t]*('+FLOAT+rb')[ \t]*,[ \t]*('+FLOAT+
                 rb')[ \t]*,[ \t]*)('+FLOAT+rb')([ \t]*(?:\r?\n)?)')


def load(path):
    return json.loads(path.read_bytes())


def sha(x):
    return hashlib.sha256(x).hexdigest()


def pin_all(spec):
    out = {}
    for key, record in spec.items():
        p = Path(record['path'])
        require(p.is_file() and p.stat().st_size==record['bytes'] and sha(p.read_bytes())==record['sha256'],
                f'input identity mismatch: {key}')
        out[key] = p
    return out


def rows(raw):
    require(not re.search(rb'(?mi)^[ \t]*\*(?:include|imperfection)\b',raw), 'unsupported include or imperfection card')
    require(len(re.findall(rb'(?mi)^\*part, name=membrane\s*$',raw))==1, 'ambiguous membrane part')
    sets = list(MEMBRANE.finditer(raw))
    require(len(sets)==1, 'missing or ambiguous membrane NODE card')
    start = sets[0].start(1)
    result = {}
    for line in sets[0][1].splitlines(keepends=True):
        m = ROW.fullmatch(line)
        require(m is not None, 'unexpected membrane coordinate row')
        label = int(m[2]); coord = tuple(float(m[i]) for i in (3,4,5))
        require(label not in result and all(map(math.isfinite,coord)), 'duplicate label or nonfinite coordinate')
        result[label] = (line,coord,start+m.start(5),start+m.end(5))
        start += len(line)
    return result


def table(path):
    ans = {}
    for text in path.read_text(encoding='utf-8').splitlines():
        v = text.strip()
        if not v or v.startswith('#'):
            continue
        name,value = v.split()
        if name=='MODE':
            key = int(value)
            require(key==len(ans)+1, 'duplicate or unordered mode')
            ans[key] = {}
            continue
        label = int(name); a = float(value)
        require(ans and label not in ans[key] and math.isfinite(a), 'missing block, duplicate label or nonfinite mode')
        ans[key][label] = value
    return ans


def verify_coverage(selected, corners, mids, topology):
    seen = {}
    for kind,c in topology.values():
        n = 4 if kind=='S8R5' else 3
        for i,l in enumerate(c[n:]):
            parents = sorted((c[i],c[(i+1)%n]))
            require(seen.setdefault(l,parents)==parents, 'inconsistent shared midside parents')
    require(seen==mids and not (corners & seen.keys()), 'midside map/topology conflict')
    require(all(m.keys()==corners|mids.keys() for m in selected), 'incomplete lifted modes/midsides')


def check_modes(p, spec, labels, connectivity):
    cfg, transfer, verified, review, native = map(load,(p['config'],p['transfer_manifest'],p['transfer_review'],
                                                       p['source_review'],p['native_modes']))
    rec = transfer['recipe']
    require(cfg['IMPERFECTION_MODES']==[1,2,3,4] and cfg['THICKNESS']==5e-6 and cfg['IMPERFECTION_AMPL_T']==.1 and
            rec['mode_order']==['paper:1','paper:2','clamp_only:1','clamp_only:2'] and rec['coefficients']==[1]*4 and
            rec['component']=='U3' and verified['numerical_transfer']=='verified' and native['readOnly'] and
            native['manifest_sha256']==spec['transfer_manifest']['sha256'], 'unverified family, recipe or amplitude')
    amp = cfg['THICKNESS']*cfg['IMPERFECTION_AMPL_T']
    require(verified['amplitude']['per_mode_coefficient']==amp and review['native_extraction']['sha256']==spec['native_modes']['sha256'],
            'amplitude or native record identity mismatch')
    all_lifted = table(p['table'])
    require(set(all_lifted)=={1,2,3,4},'expected exactly four lifted blocks')
    selected = list(all_lifted.values())
    mids = {int(k):sorted(v['edge_endpoints']) for k,v in load(p['midmap'])['midsides'].items()}
    corners = labels-mids.keys()
    verify_coverage(selected,corners,mids,connectivity)
    require(len(corners)==10201 and len(mids)==20300 and all(set(m)<=corners for m in mids.values()), 'wrong corner/midside coverage')
    bound_max=0.;mid_max=0.;ident=[]
    for i,family in enumerate(FAMILIES):
        base = table(p[family]); source=transfer['source_candidates'][family]; n=native['sources'][family]
        qual=review['sources'][family]
        require(qual['solver']['converged_count']>=2 and n['before']['sha256']==source['odb']['sha256']==n['after']['sha256'] and
                source['source_table']['sha256']==spec[family]['sha256'] and n['instance']=='MEMBRANE-1',
                'selected family does not bind to the source record')
        for k in (1,2):
            f = n['modes'][k-1]; proof=qual['selected_modes'][k-1]; block=selected[2*i+k-1]
            require(f['native_mode']==f['table_mode']==k and proof['quantized_rows_match'] and
                    f['frame_index']==proof['frame_index']==source['picks'][k-1]['frame_index'] and
                    f['eigenvalue_literal']==proof['eigenvalue_literal']==source['picks'][k-1]['eigenvalue_literal'] and
                    f['sign_node']==proof['sign_node']==source['picks'][k-1]['sign_node'], 'mode, frame, eigenvalue or sign mismatch')
            raw = {u[0]:u[3] for u in f['native_u']}
            require(raw.keys()==base[k].keys()==corners and len(f['native_u'])==len(corners) and all(map(math.isfinite,raw.values())),
                    'incomplete or nonfinite native source field')
            peak=abs(raw[f['sign_node']]);require(peak>0 and peak==max(map(abs,raw.values())), 'native sign node is not peak')
            sign=math.copysign(1,raw[f['sign_node']]); orig=base[k]
            require(all(format(sign*val/peak,'.9e')==orig[j] and float(orig[j])==float(block[j]) for j,val in raw.items()),
                    'source normalization or lifted corner field mismatch')
            for j,(a,b) in mids.items():
                value=(float(orig[a])+float(orig[b]))/2; text=block[j]; got=float(text)
                e=int(text.lower().split('e')[1]); bound=.5*10.**(e-10)+4*max(math.ulp(value),math.ulp(got))
                error=abs(got-value)
                require(error<=bound, 'midside exceeds the qualified serialization bound')
                bound_max=max(bound_max,bound);mid_max=max(mid_max,error)
            stored_header=p[family].read_text(encoding='utf-8').splitlines()[f['table_metadata_line']-1]
            ident.append({'family':family,'native_mode':k,'observed_frame_index':f['frame_index'],
                          'eigenvalue':f['eigenvalue_literal'],'stored_header':stored_header,'source_sign_multiplier':sign,
                          'source_peak':peak,'coefficient':1,'composite_block':2*i+k,'source_odb_sha256':n['before']['sha256']})
    return [{j:float(x) for j,x in b.items()} for b in selected],amp,ident,mid_max,bound_max


def reflection(ref, edge):
    # Exact decimal labelled reference, not the (rounded/seeded/deformed) submitted X/Y.
    coord = {i:tuple(Decimal(x.decode('ascii')) for x in ROW.fullmatch(v[0]).group(3,4,5)) for i,v in ref.items()}
    by_xyz = {x:i for i,x in coord.items()}
    require(len(by_xyz)==len(coord),'reference labels do not uniquely determine reflection')
    lo=min(x[0] for x in coord.values()); hi=max(x[0] for x in coord.values())
    out={i:by_xyz.get((lo+hi-x,y,z)) for i,(x,y,z) in coord.items()}
    require(all(out.get(j)==i for i,j in out.items()) and {out[i] for i in edge}==edge,'invalid reflection label map')
    return out


def rms(field):
    r=math.sqrt(math.fsum(x*x for x in field)/len(field))
    require(math.isfinite(r) and r>0,'nonfinite or zero RMS')
    return r


def invert(A,p2,edge):
    require(A.keys()==p2.keys() and edge<=A.keys(),'incomplete contribution or invalid edge mask')
    inner=sorted(A.keys()-edge)
    require(inner and all(map(math.isfinite,(*A.values(),*p2.values()))),'empty interior or nonfinite input field')
    a=[A[i] for i in inner]; signs={i:A[i]-2*p2[i] for i in inner}
    ar,br=rms(a),rms(list(signs.values())); scale=ar/br
    require(math.isfinite(scale) and scale>0,'invalid amplitude scale')
    return {i:A[i] if i in edge else scale*signs[i] for i in A},inner,scale,ar,br


def edit(raw, parsed, new, edge):
    p=[];cursor=0
    for label,(_,_,s,e) in parsed.items():
        if label in edge:
            continue
        p.extend((raw[cursor:s],format(new[label],'.17g').encode('ascii')))
        cursor=e
    return b''.join((*p,raw[cursor:]))


def verify_edit(before,after,new,edge):
    # Independent byte-by-byte reverse audit; do not call edit() to construct its oracle.
    a,b=rows(before),rows(after)
    require(list(a)==list(b) and a.keys()==new.keys() and edge<=a.keys(),'label or ordering mismatch')
    fixed=bytearray(); pos=0;changed=[]
    for j,(line,coord,s,e) in b.items():
        old,oldxyz,bs,be=a[j]
        if j in edge:
            require(line==old,'changed geometric edge row')
            continue
        require(coord[:2]==oldxyz[:2] and coord[2].hex()==new[j].hex(),'Z ledger or in-plane coordinate mismatch')
        fixed.extend(after[pos:s]); fixed.extend(before[bs:be]); pos=e
        if after[s:e]!=before[bs:be]:
            changed.append(j)
    fixed.extend(after[pos:]);require(fixed==before,'bytes outside permitted interior Z tokens changed')
    return b,changed


def get_stats(v, I):
    x=[v[i] for i in I]; n=rms(x); j=max(I,key=lambda i:abs(v[i]))
    return {'rms_m':n,'peak_abs_m':abs(v[j]),'peak_label':j,'peak_over_rms':abs(v[j])/n,'mean_m':math.fsum(x)/len(x)}


def f32(a):
    return struct.unpack('<f',struct.pack('<f',a))[0]


def floats_on_native_ledger(path,base,ref,requested):
    done=set(); mx=0.; bound=0.
    with path.open(newline='',encoding='ascii') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            j=int(r['label']);require(j in base and j not in done,'missing or duplicate ledger label');done.add(j)
            xyz=tuple(float.fromhex(r['after_'+s]) for s in 'xyz')
            pre=tuple(float.fromhex(r['before_'+s]) for s in 'xyz')
            w=float.fromhex(r['requested_dz'])
            require(w==requested[j] and pre[2]==ref[j][1][2]==0 and xyz[:2]==pre[:2] and
                    all(f32(t)==s==f32(u) for t,s,u in zip(base[j][1][:2],xyz[:2],ref[j][1][:2])) and
                    xyz[2].hex()==base[j][1][2].hex() and xyz[2]==f32(w), 'A ledger not bound to actual applied geometry')
            mx=max(mx,abs(w-xyz[2])); b=2.**(math.frexp(xyz[2])[1]-25) if xyz[2] else 2.**-150
            bound=max(bound,b)
    require(done==base.keys(),'incomplete A ledger')
    return mx,bound


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pins',required=True,type=Path);p.add_argument('--out',required=True,type=Path)
    a=p.parse_args(); out=a.out; spec=load(a.pins); f=pin_all(spec)
    with (out/'seedb-recipe.json').open('x',encoding='utf-8') as file:
        dat=prepare(f,spec,out);json.dump(dat,file,indent=2,allow_nan=False);file.write('\n')


def prepare(f,spec,out):
    direct,control=(f[k].read_bytes() for k in ('clamp','control'))
    r=load(f['reviewed_ranges']);verify_intervention(direct,control,r['deleted'],r['replace_set_membership'])
    original, ref, ctrl=rows(direct),rows(f['reference'].read_bytes()),rows(control)
    require(original.keys()==ref.keys()==ctrl.keys() and len(ref)==30501 and all(original[i][0]==ctrl[i][0] for i in ref),
            'A coordinate fields not byte-identical or label coverage differs')
    pts,conn=membrane(f['reference']);require(len(pts)==30501 and all(p[2]==0. for p in pts.values()),'not a planar unseeded reference')
    for k in ('clamp','control'):
        require(membrane(f[k],'MEMBRANE-1')[1]==conn,'applied ordered quadratic topology mismatch')
    m=load(f['membership']); G=set(m['intended_edge']); E=set(m['unseeded']['edges']); X=set(m['unseeded']['excluded_edges'])
    policy=load(f['boundary_manifest'])['labels']
    require(len(G)==800 and len(E)==798 and X==CENTRES and not E&X and G==E|X<=ref.keys() and
            all(set(policy[k])==s for k,s in (('G',G),('E',E),('X',X))), 'wrong authoritative G/E or exclusion labels')
    require(set(r['replace_set_membership'][0][2])==G,'control not bound to full geometric edge')
    modes, amp, lineage,mid,bound=check_modes(f,spec,ref.keys(),conn)
    requested={}
    with f['A_offset_csv'].open(newline='',encoding='ascii') as z:
        for row in csv.DictReader(z):
            i=int(row['label']); require(i in ref and i not in requested,'bad qualified field label')
            v=0.
            for mode in modes:
                v+=amp*mode[i]
            v=0. if i in E else v
            require(tuple(float(row[k]) for k in ('dx_m','dy_m','dz_m'))==(0.,0.,v),'qualified offset does not recombine')
            requested[i]=v
    require(requested.keys()==ref.keys(),'missing qualified applied-field labels')
    aerr,abound=floats_on_native_ledger(f['build_ledger'],original,ref,requested)
    A={j:r[1][2]-ref[j][1][2] for j,r in original.items()}; p2={i:amp*x for i,x in modes[1].items()}
    require(all(A[i]==0 for i in E) and all(A[i]!=0 for i in X),'A edge values contradict verified policy')
    B,I,k,ra,rb=invert(A,p2,G); map_lr=reflection(ref,G)
    # Numerical resolution witness, not a branch distance: bounded displayed storage/transfer quantization.
    ea=abound; er=ea+2*amp*bound; require(rb>er and ra>ea,'recipe scale unresolved under recorded precision')
    kl,kh=(ra-ea)/(rb+er),(ra+ea)/(rb-er)
    amp_error=ea+kh*er+max(abs(A[i]-2*p2[i]) for i in I)*max(abs(k-kl),abs(kh-k))+4*max(math.ulp(B[i]) for i in I)
    distinct={}
    for name,refl,sign in (('A',False,1),('-A',False,-1),('reflection_A',True,1),('-reflection_A',True,-1)):
        d=[B[i]-sign*A[map_lr[i] if refl else i] for i in I]
        nr=math.sqrt(math.fsum(x*x for x in d)/len(d))
        distinct[name]={'rms_residual_m':nr,'max_abs_residual_m':max(map(abs,d)),'numerical_precision_scale_m':amp_error}
        require(nr>amp_error,'mixture not distinguishable from '+name)
    decks={};results={}
    for name,raw,r0 in (('control',control,ctrl),('clamp',direct,original)):
        changed=edit(raw,r0,B,G); rb0,labs=verify_edit(raw,changed,B,G)
        read_rms=rms([rb0[j][1][2]-ref[j][1][2] for j in I]);require(abs(read_rms-ra)<=1e-12*ra,'written RMS mismatch')
        decks[name]=(changed,rb0)
        results[name]={'sha256':sha(changed),'bytes':len(changed),'changed_token_count':len(labs),'parsed_interior_rms_m':read_rms,
                       'initial_edge_rows_unchanged':800,'only_interior_Z_tokens_changed':True}
    require(all(decks['control'][1][j][0]==decks['clamp'][1][j][0] for j in ref),'B pair membrane bytes differ')
    verify_intervention(decks['clamp'][0],decks['control'][0],r['deleted'],r['replace_set_membership'])
    for name,(data,_) in decks.items():
        raw=control if name=='control' else direct
        with (out/(name+'-seedb.inp')).open('xb') as fp: fp.write(data)
        with (out/(name+'-seedb.diff')).open('xb') as fp:
            fp.writelines(difflib.diff_bytes(difflib.unified_diff,raw.splitlines(keepends=True),data.splitlines(keepends=True),
                                            (name+'-A').encode(),(name+'-B').encode()))
    with (out/'seedb-ledger.csv').open('x',encoding='ascii',newline='') as fp:
        wr=csv.writer(fp);wr.writerow(('label','geometric_edge','reflection_label','z_reference_hex','applied_z_A_hex','p2_contribution_hex','z_B_hex'))
        wr.writerows((i,int(i in G),map_lr[i],ref[i][1][2].hex(),A[i].hex(),p2[i].hex(),B[i].hex()) for i in sorted(ref))
    result={'decks':results,'pair_membrane_rows_byte_identical':True,'original_structural_intervention_verified':True,
            'unchanged_full_G':sorted(G),'E':sorted(E),'edge_z_tokens':{j:ROW.fullmatch(original[j][0])[5].decode('ascii') for j in sorted(G)},
            'seedb_ledger_sha256':sha((out/'seedb-ledger.csv').read_bytes()),'shared_geometry_rows':len(ref)}
    with (out/'paired-input-verification.json').open('x',encoding='utf-8') as fp: json.dump(result,fp,indent=2,allow_nan=False)
    pin_all(spec)
    return {'inputs':spec,'mode_lineage':lineage,'geometry':'one prescribed interior-only normal field, not a native mode',
            'formula':'rawB=wA-2*P2; B=k*rawB on I; same X/Y and all geometric edge tokens',
            'interior_labels':len(I),'original_per_mode_amplitude_m':amp,'scale':k,'rawB_rms_m':rb,
            'precision':{'A_requested_vs_stored_max_m':aerr,'A_binary32_rounding_scale_m':abound,
                         'max_midside_transfer_residual':mid,'max_transfer_text_rounding_scale':bound,
                         'emission':'binary64 .17g round trip, no second application'},
            'distinctness_vs_A_symmetries':distinct,'interior':{'A':get_stats(A,I),'B':get_stats(B,I)},
            'whole_field':{'A':get_stats(A,list(ref)),'B':get_stats(B,list(ref))},'edges_frozen':True}


if __name__=='__main__':
    main()
