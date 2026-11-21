"""Isolated CAE write-input gate. No solver-facing Job API is exposed."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import sys

ATTEMPT = 'hf-edge-native798-v1'
JOB = 'SailState2_HF_Native798_Q2_01'
INSTANCE = 'MEMBRANE-1'
FRAME = {'R': [[1,0,0],[0,1,0],[0,0,1]], 't_m': [0,0,0], 'units': 'm'}
CHANGED = [5018, 5160, 5302, 5440, *range(30357, 30365)]
CENTERS = {5161: 1.6910175327000004e-10, 5301: 9.3843595115e-11}
CLOUD_NAMES = {'Top':'top','Right':'right','Left':'left','CL':'clamp_left','CR':'clamp_right'}
LABEL_KEYS = {'G','X','E','base','diagonal','clouds','changed'}
SOURCES = {
    'membership': 'baseline-native-memberships.json',
    'reference_mesh': 'reference-mesh.inp',
    'reference_identity': 'reference-identity.json',
    'coating_sides': 'coating-sides.json',
    'base_config': 'hf-q2.json',
    'mode_table': 'verified-transferred-table.txt',
}
PRODUCTS = {'config':'hf-build-config.json','offsets':'derived-offsets.csv',
            'coordinates':'derived-seeded-coordinates.csv'}
BASELINE_SHA = '686e1b9649915c736a4bd7f1a23e2c9dbb263cd73f489b1be55b8ffa3f91fd2b'


class GuardViolation(BaseException):
    """Never caught as a recoverable driver/probe error."""


class BuildComplete(BaseException):
    """Only the completed, write-only driver can signal success."""

    def __init__(self, job):
        require(isinstance(job, WriteOnlyJob) and job.successful_writes == 1,
                'Driver has no successful native write')
        self.count = job.successful_writes


def require(yes,reason):
    if not yes: raise GuardViolation(reason)


def process_hook(event,args):
    if event=='subprocess.Popen' or event.startswith(('os.exec','os.spawn','os.posix_spawn','os.system','os.startfile','os.fork','_winapi.CreateProcess')):
        raise GuardViolation('Subprocess forbidden during CAE build: '+event)


def install_process_guard():
    sys.addaudithook(process_hook)
    try:sys.audit('subprocess.Popen',None,None,None,None)
    except GuardViolation:return
    raise GuardViolation('Subprocess audit hook did not arm')


class WriteOnlyJob:
    def __init__(self, job):
        self._job = job
        self._attempted = False
        self.successful_writes = 0

    def writeInput(self, **kwargs):
        require(not self._attempted, 'Only one input write authorized')
        self._attempted = True
        result = self._job.writeInput(**kwargs)
        self.successful_writes = 1
        return result
    def submit(self,*args,**kwargs):raise GuardViolation('Job.submit is forbidden')
    def waitForCompletion(self,*args,**kwargs):raise GuardViolation('Job.waitForCompletion is forbidden')


class WriteOnlyMdb:
    def __init__(self,mdb):self._mdb=mdb;self._allocated=False
    @property
    def models(self):return self._mdb.models
    @property
    def jobs(self):return self._mdb.jobs.keys()
    def Model(self,**kwargs):return self._mdb.Model(**kwargs)
    def Job(self,**kwargs):
        require(kwargs.get('name')==JOB and not self._allocated and not self.jobs,'Unauthorized or repeated job creation')
        require(kwargs.get('numCpus')==kwargs.get('numDomains')==4,'Configuration changed')
        self._allocated=True
        return WriteOnlyJob(self._mdb.Job(**kwargs))
    def JobFromInputFile(self,*args,**kwargs):raise GuardViolation('JobFromInputFile is forbidden')


def arm_mdb(mdb):
    require(os.environ.get('SAIL_BUILD_ONLY')=='1','Not entered through guarded build wrapper')
    return WriteOnlyMdb(mdb)


def strict_pairs(pairs):
    value = {}
    for key, item in pairs:
        require(key not in value, f'Duplicate JSON key: {key}')
        value[key] = item
    return value


def invalid_constant(value):
    raise GuardViolation(f'Nonfinite JSON value: {value}')


def audit_state(path):
    return json.loads(Path(path).read_text(encoding='utf-8'),
                      object_pairs_hook=strict_pairs, parse_constant=invalid_constant)


def identity(root, path):
    path = path.resolve(strict=True)
    data = path.read_bytes()
    return dict(path=path.relative_to(root).as_posix(), sha256=hashlib.sha256(data).hexdigest(), bytes=len(data))


def label_array(array, labels, name):
    require(isinstance(array,list) and array and all(type(i) is int and i in labels for i in array)
            and array == sorted(set(array)), f'Invalid label array: {name}')
    return array


def check_policy(package, manifest):
    require(set(manifest) == {'schema_version','policy','attempt','job','frame','baseline','parents','products',
                              'labels','lineage','peak'} and type(manifest['schema_version']) is int and
            manifest['schema_version'] == 1 and
            manifest['policy'] == 'flat_reference_native_798' and manifest['attempt'] == ATTEMPT and
            manifest['job'] == JOB and manifest['frame'] == FRAME, 'Unexpected policy/frame/attempt')
    for key, expected in (('baseline',SOURCES),('products',PRODUCTS)):
        require(set(manifest[key]) == set(expected), f'Manifest {key} paths differ')
        for name, leaf in expected.items():
            require(manifest[key][name] == identity(package,package/'inputs'/leaf),
                    f'Manifest {key} hash/path mismatch: {name}')
    require(manifest['baseline']['membership']['sha256'] == BASELINE_SHA,
            'Unexpected historical membership identity')
    labels = manifest['labels']
    require(set(labels) == LABEL_KEYS and set(labels['clouds']) == set(CLOUD_NAMES.values()),
            'Bad label contract')
    universe = range(1,30502)
    for name, array in labels.items():
        if name != 'clouds':
            label_array(array,universe,name)
    for name, array in labels['clouds'].items():
        label_array(array,universe,name)
    base = audit_state(package/'inputs'/SOURCES['membership'])
    expected = base['unseeded']
    require(labels['G'] == base['intended_edge'] and labels['X'] == expected['excluded_edges'] == list(CENTERS)
            and labels['E'] == expected['edges'] and labels['clouds'] == expected['clouds']
            and len(labels['G']) == 800 and len(labels['E']) == 798 and
            set(labels['E']) == set(labels['G']) - set(labels['X']) and labels['changed'] == CHANGED and
            set(labels['base']).isdisjoint(labels['diagonal']) and
            sorted(labels['base']+labels['diagonal']) == labels['E'] and len(labels['base']) == 401,
            'Manifest baseline mask/cloud/split drift')
    require([r['label'] for r in manifest['lineage']] == CHANGED and
            all(set(r) == {'label','source_labels','parent_dz_m','derived_dz_m'} and
                r['derived_dz_m'] == 0 and math.isfinite(r['parent_dz_m']) and r['parent_dz_m'] != 0
                for r in manifest['lineage']), 'Incomplete or invalid lineage')
    before = audit_state(package/'inputs'/SOURCES['base_config'])
    after = audit_state(package/'inputs'/PRODUCTS['config'])
    require(after == dict(before,HF_JOB=JOB), 'Effective config differs beyond job name')
    return labels


def record(name,d):
    with Path(name).open('x',encoding='utf-8') as f:json.dump(d,f,indent=2,allow_nan=False)


def coordinates(obj):
    value = tuple(map(float,obj.coordinates))
    require(len(value) == 3 and all(map(math.isfinite,value)), 'Invalid native coordinate')
    return value


def snapshot(nodes, instance=None):
    table = {}
    for node in nodes:
        label = int(node.label)
        require(label not in table and label > 0, 'Duplicate or invalid native node label')
        if instance is not None:
            require(node.instanceName == instance, 'Wrong native node instance')
        table[label] = coordinates(node)
    require(bool(table),'Empty native node array')
    return table


def storage_matches(value, expected):
    return value == expected or value == struct.unpack('<f',struct.pack('<f',expected))[0]


def native_mesh(part, instance, reference, elements, imported_grid):
    require(imported_grid.read_bytes() == Path(os.environ['STATE2_REFERENCE_MESH']).read_bytes(),
            'Generated import grid differs from frozen bytes')
    before = snapshot(part.nodes)
    require(instance.name == INSTANCE and snapshot(instance.nodes,INSTANCE) == before,
            'Unexpected membrane instance placement or labels')
    require(before.keys() == reference.keys() and all(p[2] == 0. and
            all(storage_matches(v,w) for v,w in zip(p,reference[i])) for i,p in before.items()),
            'Imported coordinate or label drift')
    for mesh in (part,instance):
        typed = {}
        for element in mesh.elements:
            label = int(element.label)
            require(label not in typed, 'Duplicate native element label')
            # MeshElement.connectivity is an array of internal indices, not labels.
            typed[label] = (str(element.type),tuple(int(mesh.nodes[i].label) for i in element.connectivity))
        require(typed == elements, 'Native element type or ordered connectivity mismatch')
    return before


def reference_lineage(elements, label):
    parents = set()
    for kind, conn in elements.values():
        corners = 4 if kind == 'S8R5' else 3
        if label in conn[corners:]:
            slot = conn.index(label)-corners
            parents.add(tuple(sorted((conn[slot],conn[(slot+1)%corners]))))
    require(len(parents) <= 1, f'Ambiguous reference midside {label}')
    return list(next(iter(parents))) if parents else [label]


def projection(before, combined, manifest, elements, path):
    labels = manifest['labels']
    require(not Path('applied-nodes.tsv').exists(), 'Geometry already edited')
    require(before.keys() == combined.keys() and all(p[2] == 0 for p in before.values())
            and all(map(math.isfinite,combined.values())), 'Bad planar reference or modal coverage')
    with path.open(encoding='ascii',newline='') as stream:
        reader = csv.reader(stream)
        require(next(reader) == ['label','dx_m','dy_m','dz_m'], 'Offset format mismatch')
        offsets, previous = {}, 0
        for label,x,y,z in reader:
            i, dz = int(label), float(z)
            require(i > previous and float(x) == float(y) == 0 and math.isfinite(dz),
                    'Invalid projected offset row')
            offsets[i], previous = dz, i
    edge = set(labels['E'])
    changed = [i for i in labels['E'] if combined[i] != 0.]
    require(changed == CHANGED == labels['changed'] and offsets.keys() == before.keys() and
            all(offsets[i] == (0. if i in edge else combined[i]) for i in offsets) and
            all(offsets[i] == value for i,value in CENTERS.items()), 'Native projection differs from candidate')
    require([row['label'] for row in manifest['lineage']] == changed, 'Incomplete native lineage')
    for row in manifest['lineage']:
        i = row['label']
        require(row['source_labels'] == reference_lineage(elements,i) and
                row['parent_dz_m'] == combined[i] and row['derived_dz_m'] == offsets[i],
                f'Native lineage mismatch: {i}')
    for name, values in (('parent',combined),('derived',offsets)):
        peak = max(sorted(values), key=lambda i: abs(values[i]))
        require(manifest['peak'][name] == dict(label=peak,abs_dz_m=abs(values[peak])),
                f'Native peak mismatch: {name}')
    return offsets


def resolve_nodes(container, expected, instance=None):
    nodes = container.sequenceFromLabels(labels=tuple(expected))
    require(set(snapshot(nodes,instance)) == set(expected) and len(nodes) == len(expected),
            'Native sequenceFromLabels did not resolve exact labels')
    return nodes


class ApplyOnce:
    def __init__(self):
        self._attempted = False

    def apply(self, part, before, offsets):
        require(not self._attempted and offsets.keys() == before.keys() and
                all(p[2] == 0 for p in before.values()), 'Imperfection already applied or bad reference')
        nodes = resolve_nodes(part.nodes, sorted(before))
        require(snapshot(nodes) == before and all(map(math.isfinite,offsets.values())),
                'Pre-edit coordinates or offsets drifted')
        edited = tuple((x,y,z+offsets[node.label]) for node in nodes for x,y,z in (before[node.label],))
        self._attempted = True
        part.editNode(nodes=nodes,coordinates=edited)
        after = snapshot(part.nodes)
        require(after.keys() == before.keys() and all(after[i][:2] == p[:2] and
                storage_matches(after[i][2],p[2]+offsets[i]) for i,p in before.items()),
                'Native edit changed unintended components')
        return after


def membership(instance,grid,base,height,clamp_centers,clamp_r,centers,expected):
    xyz = snapshot(instance.nodes,INSTANCE)
    x = sorted(grid.clamp_exclude_labels(instance,clamp_centers,clamp_r))
    edge = sorted(grid.boundary_node_labels(instance,base=base,height=height,exclude_labels=x))
    clouds = {k:sorted(snapshot(instance.nodes.getByBoundingSphere(center=coord,radius=r),INSTANCE))
              for k,(coord,r) in centers.items()}
    base_edge = [i for i in edge if abs(xyz[i][1]) < 1e-9]
    split = set(base_edge)
    result = dict(X=x,E=edge,base=base_edge,diagonal=[i for i in edge if i not in split],clouds=clouds)
    require(all(result[k] == expected[k] for k in result), 'Native unseeded membership differs from baseline')
    return result


def expected_sets(membership):
    return {'All_Edges_NoClamp':membership['E'],'BASE_EDGE':membership['base'],
            'DIAG_EDGES':membership['diagonal'],
            **{name+'_Nodes': membership['clouds'][key] for name,key in CLOUD_NAMES.items()}}


def make_label_set(assembly, instance, name, expected):
    assembly.Set(name=name,nodes=resolve_nodes(instance.nodes,expected,INSTANCE))
    result = assembly.sets[name]
    require(sorted(snapshot(result.nodes,INSTANCE)) == expected, f'Wrong constructed set: {name}')
    return result


def final_sets(assembly, membership, after):
    instance = assembly.instances[INSTANCE]
    require(snapshot(instance.nodes,INSTANCE) == after, 'Dependent instance differs from part ledger')
    recorded = {}
    for name, expected in expected_sets(membership).items():
        require(name in assembly.sets, f'Missing final set: {name}')
        actual = sorted(snapshot(assembly.sets[name].nodes,INSTANCE))
        require(actual == expected, f'Wrong final set: {name}')
        recorded[name] = actual
    return recorded


def require_xy_mirror(grid, left, right):
    xy = [grid.patch_coords_from_set(side) for side in (left,right)]
    require(all(array and all(len(p) == 2 and all(map(math.isfinite,p)) for p in array)
                for array in xy), 'Unavailable or invalid clamp XY')
    ok, reason = grid.patch_nodes_mirror_ok(*xy,axis=10.,tol=1e-6)
    require(ok, f'Clamp XY mirror check failed: {reason}')


def prepare_patch(path):
    """Check the raw semantics before the legacy patcher; use LF for its working copy."""
    path = Path(path)
    data = path.read_bytes()
    text = data.decode('utf-8').replace('\r\n','\n')
    require('\r' not in text,'Unsupported raw line terminators')
    for line in text.splitlines():
        require(line.partition(',')[0].strip().casefold() != '*controls', 'Unexpected raw controls')
    normalized = text.encode('utf-8')
    if normalized != data:
        path.write_bytes(normalized)


def seal(root, completion, package_digest):
    require(type(completion) is BuildComplete and completion.count == 1 and
            len(package_digest) == 64 and all(c in '0123456789abcdef' for c in package_digest),
            'Missing successful driver completion or approved inventory digest')
    out = root/'staged-code/aba'
    required = ('applied-nodes.tsv','effective-hf.json','native-memberships.json','grid_Membrane.inp',
                'attachments.json',JOB+'.raw.inp',JOB+'.inp')
    optional = ('hf_run_log.txt',)
    names = {p.name for p in out.iterdir()}
    require(set(required) <= names <= set(required+optional) and
            all((out/n).is_file() and not (out/n).is_symlink() for n in names),
            'Missing or unexpected build artifacts')
    record(out/'build-seal.json',dict(schema_version=1,attempt=ATTEMPT,job=JOB,
           package_sha256=package_digest,native_write_count=completion.count,
           raw=identity(root,out/(JOB+'.raw.inp')),sealed=identity(root,out/(JOB+'.inp'))))
