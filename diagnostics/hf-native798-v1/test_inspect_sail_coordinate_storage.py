import contextlib
import io
import json
import math
from pathlib import Path
from types import SimpleNamespace as NS
import struct
import sys
import tempfile
import unittest

import numpy as np

import inspect_sail_coordinate_storage as n


class Field:
    def __init__(self,columns,values):
        self.componentLabels=columns; self.values=values

    def getSubset(self,*,region,position):
        return NS(values=[s for s in self.values if s.instance is region and s.position==position])


class Value:
    def __init__(self,instance,label,bits,xyz,location='NODAL'):
        self.instance=instance;self.nodeLabel=label;self.precision='SINGLE_PRECISION' if bits==32 else 'DOUBLE_PRECISION'
        self.val=tuple(xyz);self.position=location
        self.localCoordSystem=self.localCoordSystemDouble=None

    @property
    def data(self):
        if self.precision=='DOUBLE_PRECISION':raise AssertionError('data on DOUBLE')
        return self.val

    @property
    def dataDouble(self):
        if self.precision=='SINGLE_PRECISION':raise AssertionError('dataDouble on SINGLE')
        return self.val


# Fixed IEEE encodings, independent of narrowing or the production predicate.
SINGLE_FIXTURES = (
    (0x00000001, '0x1p-149'), (0x80000001, '-0x1p-149'), (0x007fffff, '0x1.fffffcp-127'),
    (0x807fffff, '-0x1.fffffcp-127'), (0x00800000, '0x1p-126'), (0x80800000, '-0x1p-126'),
    (0x00800001, '0x1.000002p-126'), (0x80800001, '-0x1.000002p-126'), (0x7f7fffff, '0x1.fffffep+127'),
    (0xff7fffff, '-0x1.fffffep+127'), (0x00000000, '0x0p+0'), (0x80000000, '-0x0p+0'),
)


class AuditTests(unittest.TestCase):
    def test_native_precision_width_mixed_and_rotation_slots(self):
        i=NS(name='MEM');foreign=NS(name='OTHER'); val=1.+2.**-50
        u=Field(('UR1','U3','U1','U2'),[Value(i,2,64,(8.,val,-0.,.5)),Value(i,1,32,(8.,.25,4.,-1.)),
                  Value(i,1,32,(8.,3.,7.,9.),'INTEGRATION_POINT'),Value(foreign,1,64,(8.,3.,7.,9.))])
        frame=NS(fieldOutputs={'U':u})
        rows,meta=n.extract_field(frame,'U',i,(1,2),'NODAL')
        self.assertEqual(meta['precision_counts'],{'32':1,'64':1})
        np.testing.assert_array_equal(rows['U_xyz'],[[4.,-1.,.25],[-0.,.5,val]])
        self.assertTrue(np.signbit(rows['U_xyz'][1,0]));self.assertNotEqual(float(np.float32(val)),val)
        with self.assertRaises(AssertionError): _=u.values[0].data
        with self.assertRaises(AssertionError): _=u.values[1].dataDouble
        with self.assertRaises(ValueError): n.columns('U',('U1','U2','U3','U4'))
        with self.assertRaises(ValueError): n.columns('U',('U1','U1','U3'))

    def test_missing_duplicate_position_bad_components_and_bad_frame(self):
        i=NS(name='MEM');j=NS(name='OTHER'); nodal=Value(i,1,32,(1.,2.,3.));
        def parse(vary):
            field=Field(('U1','U2','U3'),vary)
            return n.extract_field(NS(fieldOutputs={'U':field}),'U',i,(1,),'NODAL')
        self.assertEqual(n.extract_field(NS(fieldOutputs={}), 'U',i,(1,),'NODAL')[1]['status'],'missing')
        for v in ([nodal,nodal],[Value(j,1,32,(1.,2.,3.))],[Value(i,1,32,(1.,2.,3.),'INTEGRATION_POINT')],
                  [Value(i,1,64,(1.,math.inf,3.))]):
            with self.assertRaises(ValueError):parse(v)
        bad=NS(componentLabels=('U1','U2','U3'),getSubset=lambda **_:NS(values=[Value(i,1,32,(1.,2.,3.),'ELEMENT_NODAL')]))
        with self.assertRaisesRegex(ValueError,'position'):
            n.extract_field(NS(fieldOutputs={'U':bad}),'U',i,(1,),'NODAL')
        f=NS(frameValue=1.,incrementNumber=3)
        odb=NS(steps={'S':NS(domain='TIME',timePeriod=1,frames=[None,f])},rootAssembly=NS(instances={'MEM':i}))
        c={'step':'S','frame':1,'increment':3,'instance':'MEM'}
        self.assertIs(n.terminal(odb,c)[2],f)
        for k,v in (('frame',0),('increment',4),('instance','OTHER')):
            with self.assertRaises(ValueError):n.terminal(odb,{**c,k:v})

    def test_read_only_call_closes_even_on_error(self):
        log=[]
        def opened(*,path,readOnly):
            log.append((path,readOnly));return NS(close=lambda:log.append('close'))
        def fail(_):raise ValueError('fixture read error')
        with self.assertRaisesRegex(ValueError,'fixture read error'):
            n.with_odb(opened,'fabricated.odb',fail)
        self.assertEqual(log,[('fabricated.odb',True),'close'])

    def test_full_readback_catches_signed_zero_and_wrong_topology(self):
        i=NS(name='MEM',nodes=[NS(label=7,coordinates=(0.,1.,2.))],elements=[NS(label=12,type='STRI65',connectivity=(7,)*6)])
        v=Value(i,7,64,(-0.,1.+2.**-50,2.));frame=NS(frameValue=1.,incrementNumber=7,fieldOutputs={'COORD':Field(('COOR1','COOR2','COOR3'),[v])})
        ini=NS(fieldOutputs={});odb=NS(steps={'S':NS(frames=[ini,frame])},rootAssembly=NS(instances={'MEM':i}))
        fields,info=n.extract_field(frame,'COORD',i,[7],'NODAL')
        arrays={'node_labels':np.array([7]),'mesh_reference_xyz':np.array([[0.,1.,2.]]),
                'element_labels_S8R5':np.array([]),'connectivity_S8R5':np.zeros((0,8)),
                'element_labels_STRI65':np.array([12]),'connectivity_STRI65':np.array([[7]*6]),
                **{'terminal_'+k:x for k,x in fields.items()}}
        fact={'fields':{'terminal':{'fields':{'COORD':info}},'initial_frame':{'fields':{}}}}
        case={'step':'S','instance':'MEM','frame':1,'increment':7}
        self.assertEqual(n.recheck(odb,case,arrays,fact)['fields']['terminal_COORD']['components'],3)
        with self.assertRaisesRegex(ValueError,'topology'):
            n.recheck(odb,case,{**arrays,'connectivity_STRI65':np.array([[8]*6])},fact)
        with self.assertRaisesRegex(ValueError,'signed-zero'):
            n.recheck(odb,case,{**arrays,'terminal_COORD_xyz':np.abs(arrays['terminal_COORD_xyz'])},fact)

    def test_native_accessor_array_npz_and_independent_readback_preserve_ieee_bits(self):
        codes, hexfloats = zip(*SINGLE_FIXTURES)
        backing = np.array(codes, dtype=np.uint32)
        accessors = backing.view(np.float32).reshape(-1, 3)
        expected = [struct.unpack('>Q', struct.pack('>d', float.fromhex(h)))[0] for h in hexfloats]
        labels = [7, 11, 20, 25]
        i = NS(name='MEM', nodes=[NS(label=k, coordinates=(0.,1.,2.)) for k in labels], elements=[])
        values = [Value(i,k,32,row) for k,row in zip(labels,accessors)]
        raw = [scalar for v in values for scalar in v.data]
        self.assertTrue(all(isinstance(s, np.float32) for s in raw))
        self.assertEqual([struct.pack('>d', float(s)).hex() for s in raw], [f'{b:016x}' for b in expected])
        terminal = NS(frameValue=1.,incrementNumber=7,
                      fieldOutputs={'COORD':Field(('COOR1','COOR2','COOR3'),list(reversed(values)))})
        initial = NS(frameValue=0.,incrementNumber=0,fieldOutputs={})
        step = NS(frames=[initial,terminal])
        odb = NS(steps={'S':step},rootAssembly=NS(instances={'MEM':i}))
        case = {'step':'S','instance':'MEM','frame':1,'increment':7}
        fields, fact = n._frame_fields(step,terminal,case,i,labels,'NODAL',keep=True)
        arrays = {'node_labels':np.array(labels),'mesh_reference_xyz':np.array([(0.,1.,2.)]*len(labels)),
                  'element_labels_S8R5':np.array([],dtype='int64'),'connectivity_S8R5':np.empty((0,8),dtype='int64'),
                  'element_labels_STRI65':np.array([],dtype='int64'),'connectivity_STRI65':np.empty((0,6),dtype='int64'),
                  **fields}
        key = 'terminal_COORD_xyz'
        self.assertEqual(arrays[key].dtype,np.dtype('float64'))
        self.assertEqual(arrays[key].view(np.uint64).ravel().tolist(),expected)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'native-fields.npz'
            n._write_native_fields(path,arrays)
            with np.load(path,allow_pickle=False) as archive:
                loaded = dict(archive)
            self.assertEqual(loaded[key].dtype,np.dtype('float64'))
            actual = loaded[key].view(np.uint64).ravel().tolist()
            self.assertEqual(actual,expected)
            readback = n.recheck(odb,case,loaded,{'fields':fact})
            self.assertEqual(readback['fields']['terminal_COORD']['components'],len(codes))
            corrupted = loaded[key].copy()
            corrupted.view(np.uint64)[0,0] = 0
            with self.assertRaisesRegex(ValueError,'readback disagree|bits disagree'):
                n.recheck(odb,case,{**loaded,key:corrupted},{'fields':fact})
        self.assertEqual(backing.tolist(),list(codes))
        print('ACCESSOR_BINARY32_TO_RELOADED_BINARY64 '+json.dumps(
            [[f'{src:08x}',f'{dst:016x}'] for src,dst in zip(codes,actual)]),flush=True)

    def test_binary64_only_values_rejected_at_single_precision_gate(self):
        # Between binary32 values, below the subnormal lattice, and above max finite.
        only_double = ('0x1.0000000000001p+0','-0x1.0000000000001p+0',
                       '0x1p-150','-0x1.000002p-149','0x1.fffffep-127',
                       '0x1.fffffe0000001p+127','0x1p+128','0x0.0000000000001p-1022')
        i = NS(name='MEM')
        for h in only_double:
            with self.subTest(value=h):
                value = (struct.unpack('>d',bytes.fromhex('0000000000000001'))[0]
                         if h.startswith('0x0.') else float.fromhex(h))
                def field(precision):
                    return NS(fieldOutputs={'COORD':Field(('COOR1','COOR2','COOR3'),
                                                          [Value(i,7,precision,(0.,0.,value))])})
                with self.assertRaises(n._SinglePrecisionMismatch) as caught:
                    n.extract_field(field(32),'COORD',i,[7],'NODAL')
                self.assertEqual(caught.exception.detail['checked_python_float_hex'], value.hex())
                rows,_ = n.extract_field(field(64),'COORD',i,[7],'NODAL')
                self.assertEqual(rows['COORD_xyz'].view(np.uint64)[0,2],
                                 struct.unpack('>Q',struct.pack('>d',value))[0])

    def test_first_mismatch_identifies_frame_component_and_unmodified_accessor(self):
        class RawFloat(float):
            def __repr__(self):
                return 'accessor(' + float(self).hex() + ')'

        i=NS(name='MEM')
        bad=RawFloat(1.+2.**-25)
        for which in ('terminal','initial_frame'):
            with self.subTest(which=which):
                def frame(role,time,inc):
                    value=Value(i,7,32,(4.,bad if role==which else 1.,-0.))
                    field=Field(('COOR3','COOR1','COOR2'),[value])
                    return NS(frameValue=time,incrementNumber=inc,fieldOutputs={'COORD':field})
                initial=frame('initial_frame',0.,0);terminal=frame('terminal',1.,12)
                with self.assertRaises(n._SinglePrecisionMismatch) as caught:
                    n._frame_fields(NS(frames=[initial,terminal]),terminal,{'step':'S','frame':9},
                                    i,{7},'NODAL',keep=True)
                d=caught.exception.detail
                self.assertEqual((d['step'],d['frame_role'],d['frame_index'],d['frame_time'],d['increment']),
                                 ('S',which,9,1.,12) if which=='terminal' else ('S',which,0,0.,0))
                self.assertEqual((d['instance'],d['field'],d['requested_position'],d['returned_position'],d['node_label']),
                                 ('MEM','COORD','NODAL','NODAL',7))
                self.assertEqual((d['component'],d['accessor_component_index']),('COOR1',1))
                self.assertEqual((d['native_precision'],d['accessor']),('SINGLE_PRECISION','data'))
                self.assertEqual(d['accessor_value_repr'],repr(bad))
                self.assertTrue(d['accessor_value_type'].endswith('.RawFloat'))
                self.assertEqual(d['checked_python_float_hex'],float(bad).hex())
                self.assertEqual(d['checked_binary64_bits'],'0x3ff0000008000000')
                self.assertFalse(d['coordinate_transform_before_check'])

    def test_single_precision_gate_checks_accessor_not_a_transformed_coordinate(self):
        i=NS(name='MEM');value=Value(i,7,32,(1.,0.,0.))
        u=math.sqrt(0.5)
        rotation=np.array([[u,-u,0.],[u,u,0.],[0.,0.,1.]])
        value.localCoordSystem=tuple(map(tuple,rotation))
        f=NS(fieldOutputs={'COORD':Field(('COOR1','COOR2','COOR3'),[value])})
        out,_=n.extract_field(f,'COORD',i,{7},'NODAL')
        np.testing.assert_array_equal(out['COORD_xyz'][0],[1.,0.,0.])
        np.testing.assert_array_equal(out['COORD_local_axes'][0],rotation)
        transformed=rotation@out['COORD_xyz'][0]+[8.,0.,0.]
        self.assertFalse(np.array_equal(transformed,transformed.astype('float32').astype('float64')))

    def test_stop_only_probe_writes_first_rejection_after_closure_without_readback(self):
        self.check_probe(mismatch=True)

    def test_stop_only_probe_reports_nonreproduction_without_export(self):
        self.check_probe(mismatch=False)

    def check_probe(self, *, mismatch):
        count=30501
        i=NS(name='MEM',nodes=[NS(label=k) for k in range(1,count+1)])
        if mismatch:
            values=[Value(i,7,32,(1.,1.+2.**-25,2.)),Value(i,8,32,(1.+2.**-26,1.,0.))]
            never=NS(getSubset=lambda **_:self.fail('probe accessed the next field after rejection'),
                     componentLabels=('U1','U2','U3'))
            fields={'COORD':Field(('COOR1','COOR2','COOR3'),values),'U':never}
        else:
            fields={'COORD':Field(('COOR1','COOR2','COOR3'),
                                 [Value(i,1,32,(1.,2.,np.array([1],dtype=np.uint32).view(np.float32)[0]))]+
                                 [Value(i,k,32,(1.,2.,-0.)) for k in range(2,count+1)])}
        terminal=NS(frameValue=1.,incrementNumber=5,fieldOutputs=fields)
        initial=NS(frameValue=0.,incrementNumber=0,fieldOutputs={})
        step=NS(domain='TIME',timePeriod=1.,frames=[initial,terminal])
        calls=[]
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'receipt';c={'step':'S','frame':1,'increment':5,'instance':'MEM'}
            for key in ('input','odb'):
                source=root/key
                source.write_bytes(key.encode('ascii'))
                c[key]=n.identity(source)
            before={k:c[k] for k in ('odb','input')}
            def close():
                self.assertFalse(out.exists())
                calls.append('closed')
            def open_odb(**kw):
                calls.append(kw)
                return NS(steps={'S':step},rootAssembly=NS(instances={'MEM':i}),close=close)
            with contextlib.redirect_stderr(io.StringIO()):
                code=n._run_probe(open_odb,'NODAL',out,c,before,{'case':'fixture'})
            self.assertEqual(code,1 if mismatch else 0)
            self.assertEqual(calls,[{'path':before['odb']['path'],'readOnly':True},'closed'])
            self.assertEqual([p.name for p in out.iterdir()],['precision-probe.json'])
            d=json.loads((out/'precision-probe.json').read_bytes())
            self.assertTrue(d['inputs_unchanged'])
            self.assertEqual(d['before'],d['after'])
            self.assertTrue(d['handle_closed'])
            self.assertFalse(d['geometry_or_raw_bundle_written'])
            if mismatch:
                self.assertEqual(d['status'],'precision_mismatch')
                self.assertEqual(d['first_mismatch']['node_label'],7)
                self.assertEqual(d['first_mismatch']['component'],'COOR2')
            else:
                self.assertEqual(d['status'],'not_reproduced')
                self.assertEqual(d['inspected_fields']['terminal']['fields']['COORD']['checked_nodes'],count)
                self.assertNotIn('first_mismatch',d)

if __name__=='__main__':
    if sys.argv[1:] == ['--no-odb']:
        import odbAccess

        attempts = []
        def forbid_open(*args, **kwargs):
            attempts.append(1)
            raise AssertionError('native openOdb is forbidden during fabricated tests')

        odbAccess.openOdb = forbid_open
        if odbAccess.openOdb is not forbid_open:
            raise RuntimeError('could not guard native openOdb')
        print('NATIVE_TEST_RUNTIME '+json.dumps(
            {'executable':sys.executable,'python':sys.version.split()[0],'numpy':np.__version__,
             'odbAccess':odbAccess.__file__}),flush=True)
        result = unittest.main(argv=[sys.argv[0]],verbosity=2,exit=False).result
        passed = result.wasSuccessful() and not attempts
        print(f'NO_ODB_TESTS {"PASS" if passed else "FAIL"} tests={result.testsRun} native_open_attempts={len(attempts)}',
              flush=True)
        sys.exit(0 if passed else 1)
    unittest.main()
