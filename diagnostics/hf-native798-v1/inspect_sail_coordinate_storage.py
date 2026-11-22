"""Read-only coordinate/translation receipt and independent native readback."""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re
import struct
import sys

import numpy as np

KINDS = {'S8R5':8,'STRI65':6}
PREC = {'SINGLE_PRECISION': (32,'data','localCoordSystem'), 'DOUBLE_PRECISION':(64,'dataDouble','localCoordSystemDouble')}
FIELDS = ('COORD','U','UT')


def require(ok,why):
    if not ok:
        raise ValueError(why)


def identity(p):
    p=Path(p).resolve(strict=True)
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    s=p.stat()
    return {'path':str(p),'sha256':h.hexdigest(),'bytes':s.st_size,'mtime_ns':s.st_mtime_ns}


def duplicate_free(items):
    v={}
    for k,x in items:
        require(k not in v,'duplicate label/key: '+str(k));v[k]=x
    return v


def columns(name,labels):
    labels=tuple(labels)
    require(len(labels)==len(set(labels)),'duplicate component name')
    if name=='COORD':
        for stem in ('COOR','COORD'):
            if set(labels)=={stem+str(i) for i in (1,2,3)}:
                return tuple(labels.index(stem+str(i)) for i in (1,2,3))
        raise ValueError('ambiguous COORD direction labels '+str(labels))
    expected={name+str(i) for i in (1,2,3)}
    require(expected<=set(labels) and set(labels)-expected<={'UR1','UR2','UR3'},
            'ambiguous translation components '+str(labels))
    return tuple(labels.index(name+str(i)) for i in (1,2,3))


def get_system(v,accessor):
    """Distinguish no metadata from an API that lacks the metadata attribute."""
    missing=object()
    local=getattr(v,accessor,missing)
    if local is missing:
        return 0, None
    if local is None or isinstance(local, (tuple,list)) and not local:
        return 1, None
    a=np.asarray(local,dtype='float64')
    require(a.shape==(3,3) and np.isfinite(a).all(),'unknown local-system metadata')
    return 2, a


class _SinglePrecisionMismatch(ValueError):
    def __init__(self, detail):
        self.detail=detail
        super().__init__(f"{detail['field']}: value contradicts native single precision")


def _binary32_exact(bits64):
    """Test IEEE binary32 representability from the bits of a Python float."""
    magnitude = bits64 & ((1 << 63) - 1)
    if not magnitude:
        return True  # Both signs of zero.
    exponent = (magnitude >> 52) - 1023
    if not -149 <= exponent <= 127:
        return False  # Also excludes binary64 subnormals and nonfinite values.
    # Binary32 normals have 23 fraction bits. Subnormals instead lie on the
    # fixed 2**-149 lattice; one binary64 significand bit weighs 2**(exponent-52).
    discarded = 29 if exponent >= -126 else -97 - exponent
    return not magnitude & ((1 << discarded) - 1)


def _check_single(raw, picked, use, components, *, field, instance, position, label,
                  precision, accessor):
    for slot, scalar in enumerate(picked):
        bits64, = struct.unpack('>Q', struct.pack('>d', scalar))
        if _binary32_exact(bits64):
            continue
        i = use[slot]
        value = raw[i]
        raise _SinglePrecisionMismatch({
            'instance':instance,'field':field,'requested_position':'NODAL','returned_position':position,
            'node_label':label,'component':components[i],'accessor_component_index':i,
            'native_precision':precision,'accessor':accessor,
            'accessor_value_repr':repr(value),
            'accessor_value_type':f'{type(value).__module__}.{type(value).__qualname__}',
            'checked_python_float_hex':scalar.hex(),
            'checked_binary64_bits':f'0x{bits64:016x}',
            'operation_order':['read native accessor result','convert scalar to Python float',
                               'select component by label','test exact binary32 representability with integer bit operations'],
            'coordinate_transform_before_check':False,
        })


def extract_field(frame,name,instance,labels,nodal):
    return _extract_field(frame,name,instance,labels,nodal,keep=True)


def _extract_field(frame,name,instance,labels,nodal,*,keep):
    if name not in frame.fieldOutputs:
        return {}, {'status':'missing','reason':'field not in this frame'}
    field=frame.fieldOutputs[name]
    comp=tuple(field.componentLabels); use=columns(name,comp)
    got={}
    for v in field.getSubset(region=instance,position=nodal).values:
        require(str(v.position)=='NODAL' and v.instance is not None and v.instance.name==instance.name,
                f'{name}: returned field position/instance differs')
        key=int(v.nodeLabel); precision=str(v.precision)
        require(key not in got and key in labels and precision in PREC,f'{name}: duplicate/foreign label or unknown precision')
        bits,accessor,coord_access=PREC[precision]
        raw=tuple(getattr(v,accessor))
        full=tuple(float(q) for q in raw)
        require(len(full)==len(comp) and all(map(math.isfinite,full)),f'{name}: invalid components')
        picked=tuple(full[i] for i in use)
        if bits==32:
            _check_single(raw,picked,use,comp,field=name,instance=instance.name,
                          position=str(v.position),label=key,precision=precision,accessor=accessor)
        sc,sm=get_system(v,coord_access)
        got[key]=(picked,bits,sc,sm) if keep else None
    require(got.keys()==set(labels),f'{name}: incomplete nodal coverage ({len(got)} vs {len(labels)})')
    if not keep:
        return {}, {'status':'complete','checked_nodes':len(got)}
    rows=[got[i] for i in labels]
    systems=[k for k,r in enumerate(rows) if r[2]==2]
    output={name+'_xyz':np.asarray([r[0] for r in rows],dtype='float64'),
            name+'_bits':np.asarray([r[1] for r in rows],dtype='uint8'),
            name+'_local_code':np.asarray([r[2] for r in rows],dtype='uint8'),
            name+'_local_rows':np.asarray(systems,dtype='int64'),
            name+'_local_axes':np.asarray([rows[k][3] for k in systems],dtype='float64').reshape(-1,3,3)}
    return output,{'status':'complete','components':list(comp),'selected_columns':list(use),'location':'NODAL','instance':instance.name,
         'precision_counts':dict(Counter(str(r[1]) for r in rows)), 'local_system_counts':dict(Counter(str(r[2]) for r in rows))}


def terminal(odb,c):
    require(list(odb.steps.keys())==[c['step']],'unexpected step inventory')
    step=odb.steps[c['step']]
    require(str(step.domain)=='TIME' and float(step.timePeriod)==1.,'unit TIME step required')
    require(len(step.frames)==c['frame']+1 and c['frame']>0,'wrong recorded terminal frame index')
    t=step.frames[c['frame']]
    require(float(t.frameValue)==1. and int(t.incrementNumber)==c['increment'],'wrong terminal time or increment')
    require(c['instance'] in odb.rootAssembly.instances,'named instance missing')
    I=odb.rootAssembly.instances[c['instance']]
    require(I.name==c['instance'],'instance identity mismatch')
    return step,I,t


def _frame_fields(step, terminal_frame, case, instance, labels, nodal, *, keep):
    fields={}
    out={}
    for which,index,frame in (('terminal',case['frame'],terminal_frame),('initial_frame',0,step.frames[0])):
        if which=='initial_frame' and (frame.frameValue!=0 or frame.incrementNumber!=0):
            fields[which]={'status':'unavailable','reason':'first saved frame is not zero-time increment zero'}
            continue
        meta={}
        for name in FIELDS:
            try:
                data,info=_extract_field(frame,name,instance,labels,nodal,keep=keep)
            except _SinglePrecisionMismatch as ex:
                ex.detail.update(step=case['step'],frame_role=which,frame_index=index,
                                 frame_time=float(frame.frameValue),increment=int(frame.incrementNumber))
                raise
            if keep:
                out.update({which+'_'+k:v for k,v in data.items()})
            meta[name]=info
        fields[which]={'all_field_names':list(frame.fieldOutputs.keys()),'time':float(frame.frameValue),
                       'increment':int(frame.incrementNumber),'fields':meta}
    return out,fields


def _probe_fields(odb,case,nodal):
    s,instance,t=terminal(odb,case)
    labels={int(n.label) for n in instance.nodes}
    require(len(labels)==len(instance.nodes)==30501,'invalid membrane node coverage')
    require('COORD' in t.fieldOutputs,'terminal nodal COORD missing; cannot replay precision check')
    _,fields=_frame_fields(s,t,case,instance,labels,nodal,keep=False)
    return fields


def extract(odb,c,nodal):
    s,I,t=terminal(odb,c)
    nodes=duplicate_free((int(n.label),(tuple(map(float,n.coordinates)),str(getattr(n,'precision',None)))) for n in I.nodes)
    labels=sorted(nodes)
    xyz=np.asarray([nodes[i][0] for i in labels],dtype='float64')
    require(len(labels)==30501 and xyz.shape==(len(labels),3) and np.isfinite(xyz).all(),'invalid membrane node inventory')
    out={'node_labels':np.asarray(labels,dtype='int64'),'mesh_reference_xyz':xyz}
    topo=duplicate_free((int(e.label),(e.type,tuple(map(int,e.connectivity)))) for e in I.elements)
    for typ,count in (('S8R5',9900),('STRI65',200)):
        v=[(k,topo[k][1]) for k in sorted(topo) if topo[k][0]==typ]
        require(len(v)==count and all(len(n)==KINDS[typ] and set(n)<=nodes.keys() for _,n in v),'wrong element topology')
        out['element_labels_'+typ]=np.asarray([k for k,_ in v],dtype='int64')
        out['connectivity_'+typ]=np.asarray([n for _,n in v],dtype='int64')
    require(len(topo)==10100,'unexpected membrane element type')
    arrays,fields=_frame_fields(s,t,c,I,labels,nodal,keep=True)
    out.update(arrays)
    return out,{'fields':fields,'mesh_api_precision':dict(Counter(v[1] for v in nodes.values())),
         'job_data':{k:str(getattr(odb.jobData,k,None)) for k in ('analysisCode','version','precision')}}


def pack3(v):
    return struct.pack('>ddd',*v)


def recheck(odb,c,arrays,meta):
    """Independent traversal; no subset, extraction columns, mapper or exporter."""
    S=odb.steps[c['step']];T=S.frames[-1];I=odb.rootAssembly.instances[c['instance']]
    require(len(S.frames)-1==c['frame'] and T.frameValue==1 and T.incrementNumber==c['increment'] and I.name==c['instance'],
            'readback terminal identity differs')
    labels={int(n):i for i,n in enumerate(arrays['node_labels'])}
    known_nodes=set()
    for n in I.nodes:
        j=int(n.label);require(j in labels and j not in known_nodes,'readback mesh labels differ')
        require(pack3(n.coordinates)==pack3(arrays['mesh_reference_xyz'][labels[j]]),'mesh reference readback differs')
        known_nodes.add(j)
    require(known_nodes==labels.keys(),'readback mesh coverage differs')
    elements={}
    for kind in KINDS:
        for i,row in zip(arrays['element_labels_'+kind],arrays['connectivity_'+kind]):
            require(int(i) not in elements,'readback repeated elements')
            elements[int(i)]=(kind,tuple(map(int,row)))
    for elem in I.elements:
        require(elements.pop(elem.label,None)==(elem.type,tuple(elem.connectivity)),'readback ordered topology differs')
    require(not elements,'readback element coverage differs')
    verified={}
    for scope,frame in (('terminal',T),('initial_frame',S.frames[0])):
        recorded=meta['fields'][scope]
        if 'fields' not in recorded:
            continue
        for name in FIELDS:
            key=scope+'_'+name
            if key+'_xyz' not in arrays:
                require(name not in frame.fieldOutputs,'missing-field status differs')
                continue
            field=frame.fieldOutputs[name]
            comp=tuple(field.componentLabels)
            require(list(comp)==recorded['fields'][name]['components'],'component labels changed')
            suffixes={label[-1]:i for i,label in enumerate(comp) if (re.fullmatch(r'COORD?[123]',label) if name=='COORD'
                            else re.fullmatch(name+r'[123]',label))}
            require(set(suffixes)=={'1','2','3'},'readback translation components missing')
            seen=set(); matrices=dict(zip(arrays[key+'_local_rows'],arrays[key+'_local_axes']))
            for entry in field.values:
                if str(entry.position)!='NODAL' or entry.instance is None or entry.instance.name!=c['instance']:
                    continue
                lab=int(entry.nodeLabel)
                require(lab in labels and lab not in seen,'readback missing/duplicate node label')
                seen.add(lab); j=labels[lab]
                if entry.precision.__str__()=='DOUBLE_PRECISION':
                    bits,native,loc=64,entry.dataDouble,'localCoordSystemDouble'
                else:
                    require(str(entry.precision)=='SINGLE_PRECISION','readback unknown precision')
                    bits,native,loc=32,entry.data,'localCoordSystem'
                expect=tuple(native[suffixes[i]] for i in '123')
                require(bits==int(arrays[key+'_bits'][j]) and pack3(expect)==pack3(arrays[key+'_xyz'][j]),
                        f'{key} node {lab}: independent value/precision or signed-zero bits disagree')
                absent=object();source=getattr(entry,loc,absent)
                code=0 if source is absent else 1 if source is None or len(source)==0 else 2
                require(code==arrays[key+'_local_code'][j] and (code!=2 or np.array_equal(matrices[j],source)),
                        'readback local metadata differs')
            require(seen==labels.keys(),key+': readback coverage differs')
            verified[key]={'node_count':len(seen),'components':3*len(seen),'promoted_numeric_and_signed_zero_bits_equal':True}
    return {'node_labels':len(labels),'ordered_elements':len(I.elements),'fields':verified,'mesh_API_values_bitwise_equal':True}


def with_odb(open_odb,path,fn):
    obj=open_odb(path=str(path),readOnly=True)
    try:
        return fn(obj)
    finally:
        obj.close()


def _write_native_fields(path, arrays):
    with path.open('xb') as f:
        np.savez_compressed(f, **arrays)


def _save_probe(out, case, before, run, result):
    after={k:identity(case[k]['path']) for k in before}
    out.mkdir()
    receipt={**run,**result,'before':before,'after':after,'inputs_unchanged':before==after,
             'readOnly':True,'handle_closed':True,'geometry_or_raw_bundle_written':False}
    with (out/'precision-probe.json').open('x',encoding='utf-8',newline='\n') as f:
        json.dump(receipt,f,indent=2,allow_nan=False)
        f.write('\n')
    require(before==after,'input or ODB identity changed')


def _run_probe(open_odb,nodal,out,case,before,run):
    try:
        fields=with_odb(open_odb,before['odb']['path'],lambda odb:_probe_fields(odb,case,nodal))
    except _SinglePrecisionMismatch as ex:
        _save_probe(out,case,before,run,{'status':'precision_mismatch','first_mismatch':ex.detail})
        print(ex,file=sys.stderr)
        return 1
    _save_probe(out,case,before,run,{'status':'not_reproduced','inspected_fields':fields})
    return 0


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--request',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--first-mismatch-only',action='store_true',
                   help='read one pinned ODB once, stop at the first precision mismatch, and write only a diagnostic receipt')
    a=p.parse_args();c=json.loads(a.request.read_bytes(),object_pairs_hook=duplicate_free)
    require(not a.out.exists() and not a.out.is_symlink(),'new output required')
    before={k:identity(c[k]['path']) for k in ('input','odb')}
    require(all(before[k]['sha256']==c[k]['sha256'] for k in before),'input or ODB pin differs')
    import odbAccess
    from abaqusConstants import NODAL
    ex=Path(sys.executable).resolve();parts=ex.parts
    require(tuple(p.lower() for p in parts[-4:])==('win_b64','code','bin','smapython.exe') and parts[-6].lower()=='estproducts',
            'not installed Abaqus Python')
    rt={'abaqus':parts[-5],'python':sys.version.split()[0],'numpy':np.__version__}
    require(rt==c['runtime'],'unexpected runtime')
    run={'case':c['case'],'request':identity(a.request),'helper':identity(__file__),'runtime':rt,
         'odbAccess':identity(odbAccess.__file__)}
    if a.first_mismatch_only:
        return _run_probe(odbAccess.openOdb,NODAL,a.out,c,before,run)
    try:
        arrays,metadata=with_odb(odbAccess.openOdb,before['odb']['path'],lambda x:extract(x,c,NODAL))
    except _SinglePrecisionMismatch as ex:
        _save_probe(a.out,c,before,run,{'status':'precision_mismatch','first_mismatch':ex.detail})
        raise
    a.out.mkdir()
    raw=a.out/'native-fields.npz'
    _write_native_fields(raw, arrays)
    del arrays
    with np.load(raw,allow_pickle=False) as b:
        verified=with_odb(odbAccess.openOdb,before['odb']['path'],lambda x:recheck(x,c,{k:b[k] for k in b.files},metadata))
    after={k:identity(c[k]['path']) for k in before}
    require(before==after,'input or ODB identity changed')
    result={**run,
            'instance':c['instance'],'step':c['step'],'frame_index':c['frame'],'terminal_increment':c['increment'],
            'before':before,'after':after,'readOnly':True,'closed_both_handles':True,'native_arrays':identity(raw),
            **metadata,'independent_readback':verified}
    with (a.out/'precision.json').open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print('Closed read-only extraction and independent readback:',c['case'])


if __name__=='__main__':
    sys.exit(main())
