"""Abaqus/CAE noGUI entry: one write-only model, no solver submission."""

from importlib.util import find_spec
import os
from pathlib import Path
import runpy
import sys

# CAE execfile may not supply __file__; the launcher pins this cwd.
sys.dont_write_bytecode = True
CODE = Path.cwd().resolve()
sys.path.insert(0,str(CODE))
from build_only_guard import ATTEMPT, BuildComplete, install_process_guard, require, seal

install_process_guard()
digest = os.environ.get('SAIL_APPROVED_PACKAGE_SHA256','')
require(len(digest) == 64 and all(c in '0123456789abcdef' for c in digest) and
        CODE.parent.name == 'sail-'+ATTEMPT, 'Missing package approval or wrong root')
import hashlib
require(hashlib.sha256((CODE.parent/'build-files.sha256').read_bytes()).hexdigest() == digest,
        'Approved inventory digest mismatch')
for name in ('aba_grid_mesh','aba_imperfection','analysis_config','build_only_guard','core',
             'ctx','post','solver','state2_run_evidence'):
    spec = find_spec(name)
    require(spec is not None and spec.origin and Path(spec.origin).resolve().is_relative_to(CODE),
            'Unexpected staged import origin: '+name)
require(find_spec('coalescence_check') is None, 'Unpackaged optional module is visible')
inputs=CODE.parent/'inputs'
require(CODE.name=='staged-code' and not (CODE/'aba').exists() and (inputs/'boundary-policy-manifest.json').is_file(),
        'Wrong or preexisting execution workspace')
for k in ('MFBO_CODE_DIR','SOLARSAIL_MODE_ODB','SOLARSAIL_MODE_DAT','SOLARSAIL_BASE_ODB','SOLARSAIL_SKIP_SEED_GATE'):
    os.environ.pop(k,None)
os.environ.update(SAIL_BUILD_ONLY='1',SOLARSAIL_MODE_TABLE=str(inputs/'verified-transferred-table.txt'),
    STATE2_REFERENCE_MESH=str(inputs/'reference-mesh.inp'),SAIL_BOUNDARY_MANIFEST=str(inputs/'boundary-policy-manifest.json'),
    SAIL_DERIVED_OFFSETS=str(inputs/'derived-offsets.csv'))
sys.argv=[str(CODE/'run_abaqus.py'),'--config',str(inputs/'hf-build-config.json'),'HF','0.7','0.01']
print('BUILD-ONLY GUARD ARMED before driver',sys.executable,sys.version)
try:
    runpy.run_path(str(CODE/'run_abaqus.py'),run_name='__main__')
except BuildComplete as completed:
    seal(CODE.parent,completed,digest)
else:
    require(False,'Driver returned without write-only completion')
