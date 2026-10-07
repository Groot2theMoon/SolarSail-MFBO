# -*- coding: utf-8 -*-
"""buckle_stub_harness.py — Abaqus 없이 run_abaqus_buckle.py 의 build_model() 을 실행한다.

무엇을 하는가 (배선 검증. 물리 검증이 아니다)
    `from abaqus import *` 를 가짜 모듈로 갈아끼우고, 실행 블록(제출/대기) **앞까지**만 잘라
    로컬 파이썬으로 build_model() 을 호출한 뒤, '어떤 BC/CLOAD 가 어느 스텝에 어떤 값으로
    만들어졌는가'를 기록한다.

무엇을 잡는가 (실측으로 잡힌 것)
    1) 없는 BC 를 참조하는 KeyError — cload 라우트에 Disp_Control_Right 가 없는데 섭동 블록이
       참조해 모델 생성 자체가 죽었다.
    2) 같은 이름을 두 createStepName 으로 만들어 **나중 것이 0 으로 덮어쓴** 사고
       (CF_Corner_* 가 Initial/GlobalTension 양쪽에 생성됨).
    둘 다 py_compile 로는 안 잡히고, Abaqus 로는 라이선스 1회 + 수십 분 뒤에야 드러난다.

무엇을 못 보는가
    lambda, 수렴, 변위장, 모드 형상 — 물리 판정은 실측(라이선스)으로만 나온다.

사용
    python3 buckle_stub_harness.py <case> [x_c] [--out PATH]
        case : CASE 이름 (seed|control_none|clamp_lf|paper_s1)
    python3 buckle_stub_harness.py --legacy <combo> [x_c]
        combo: 리팩토링 **전** 소스용 CLAMP_MODE (none|passive|driven|fixed|cload)

    케이스 상수는 소스 텍스트에서 치환한 뒤 실행한다(= '상수를 고치고 다시 돌리는' 워크플로의 모의).
    치환한 값은 실행 전에 그대로 찍으므로, 무엇으로 돌았는지 로그가 스스로 증명한다.

종료코드: 0 정상 / 2 사용법 / 3 소스 잘림 실패 / 4 head 실행 실패 / 5 build_model 예외
"""
import io
import json
import os
import re
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'run_abaqus_buckle.py')
CUT_MARKER = 'print("=" * 78)'          # 실행 블록의 시작(소스에 두 번 나온다 -> 첫 번째에서 자른다)

# 케이스가 요구하는 상수 — run_abaqus_buckle.py 의 _require_constant 와 같은 계약을 여기서도 만족시킨다.
# (하네스가 하나를 빠뜨리면 스크립트 자신의 계약 검사가 예외를 올려 드러난다)
# ⚠️ 값은 **파이썬 리터럴 그대로** 쓴다(문자열은 따옴표 포함) — 치환이 값을 그대로 끼워 넣는다.
#    실측: 따옴표를 빼면 `PRETENSION_MODE = corner2` 가 되어 NameError 로 죽는다.
CASE_CONSTANTS = {
    'seed':         {'PRETENSION_MODE': "'corner2'", 'CORNER_ANGLE_DEG': '28.6'},
    'control_none': {'PRETENSION_MODE': "'corner2'", 'CORNER_ANGLE_DEG': '28.6'},
    'clamp_lf':     {'PRETENSION_MODE': "'corner2'", 'CORNER_ANGLE_DEG': '28.6'},
    'paper_s1':     {'PRETENSION_MODE': "'paper3'",  'CORNER_ANGLE_DEG': '45.0'},
}

# ---------------------------------------------------------------------------
# 목(mock) — CAE 객체 표면이 넓으므로 범용 기록 목 + 필요한 곳만 실제 동작
# ---------------------------------------------------------------------------


def _brief(v):
    if isinstance(v, (str, int, float, bool)) or v is None:
        return v
    if isinstance(v, Obj):
        return '<obj %s>' % v._name.rsplit('.', 1)[-1]
    try:
        return '|'.join(str(x) for x in v)
    except TypeError:
        return str(v)


class Log(object):
    def __init__(self):
        self.rows = []

    def add(self, obj, call, args, kwargs):
        self.rows.append({'obj': obj, 'call': call,
                          'args': [_brief(a) for a in args],
                          'kwargs': dict((k, _brief(v)) for k, v in kwargs.items())})


class Obj(object):
    """무엇이든 기록하는 목. 속성 접근은 (캐시된) 자식 Obj, 호출은 기록 후 자식 Obj 반환.

    주의: 첫 인자를 `name` 으로 두면 BC(name=...) 처럼 `name=` 키워드를 쓰는 호출과 충돌한다
    (실측: TypeError: Obj.__init__() got multiple values for argument 'name').
    """

    def __init__(self, obj_name, log, **kw):
        object.__setattr__(self, '_name', obj_name)
        object.__setattr__(self, '_log', log)
        for k, v in kw.items():
            object.__setattr__(self, k, v)

    def __getattr__(self, k):
        log = object.__getattribute__(self, '_log')
        child = Obj('%s.%s' % (object.__getattribute__(self, '_name'), k), log)
        object.__setattr__(self, k, child)      # 반복 접근이 같은 객체를 주도록 캐시
        return child

    def __call__(self, *a, **kw):
        log = object.__getattribute__(self, '_log')
        full = object.__getattribute__(self, '_name')
        head, _, leaf = full.rpartition('.')
        log.add(head or full, leaf, a, kw)
        return Obj(full + '()', log)


class Node(object):
    __slots__ = ('label', 'coordinates')

    def __init__(self, label, coords):
        self.label = label
        self.coordinates = tuple(coords)


class NodeArray(object):
    """instance.nodes 목 — getByBoundingSphere 는 **반경을 실제로 존중**한다.

    반경을 무시하면 두 방향으로 오진한다(스킬 기록):
      - 전부 돌려주면 All_Edges_NoClamp 가 부당하게 채워져 '빈 셋' 결함이 숨는다
      - 하나도 안 돌려주면 getClosest 폴백으로 가고 make_set 이 RuntimeError 를 올린다
    """

    def __init__(self, nodes):
        self._nodes = list(nodes)
        self._by_label = dict((n.label, n) for n in self._nodes)

    def __len__(self):
        return len(self._nodes)

    def __iter__(self):
        return iter(self._nodes)

    def getByBoundingSphere(self, center=None, radius=None):
        c, r2 = center, float(radius) ** 2
        out = []
        for n in self._nodes:
            d2 = ((n.coordinates[0] - c[0]) ** 2 + (n.coordinates[1] - c[1]) ** 2
                  + (n.coordinates[2] - (c[2] if len(c) > 2 else 0.0)) ** 2)
            if d2 <= r2:
                out.append(n)
        return NodeArray(out)

    def getClosest(self, coordinates=None):
        c = coordinates
        best, bd = None, None
        for n in self._nodes:
            d = ((n.coordinates[0] - c[0]) ** 2 + (n.coordinates[1] - c[1]) ** 2
                 + (n.coordinates[2] - (c[2] if len(c) > 2 else 0.0)) ** 2)
            if bd is None or d < bd:
                best, bd = n, d
        return [(best, bd ** 0.5)]

    def sequenceFromLabels(self, labels):
        return NodeArray([self._by_label[l] for l in labels])


def _lattice_nodes(base=20.0, height=10.0, h=0.5):
    """삼각형(y>=0, y<=x, y<=base-x) 안의 격자 노드. 메쉬 재현이 아니라 **배선 검증용** 골격이다.

    경계 3변 위에 노드를 정확히 놓아야 aba_grid_mesh.boundary_node_labels 가 빈 셋을 돌려주지 않는다.
    h=0.5 면 각도가 45도인 사선이 격자 대각선과 정확히 겹친다.
    """
    out, lab = [], 0
    n = int(round(base / h))
    m = int(round(height / h))
    for j in range(m + 1):
        for i in range(n + 1):
            x, y = i * h, j * h
            if y > height + 1e-9:
                continue
            if y < -1e-9 or y > x + 1e-9 or y > base - x + 1e-9:
                continue
            lab += 1
            out.append(Node(lab, (x, y, 0.0)))
    return out


class Step(Obj):
    pass


class BC(Obj):
    """경계조건/집중하중 목 — 스텝별 값 이력을 남긴다(0 덮어쓰기 사고를 잡기 위해)."""

    def __init__(self, name, log, created_in, **kw):
        Obj.__init__(self, 'bc:%s' % name, log, name=name)
        log.add('model', 'CREATE', (name,), dict(kw, createStepName=created_in))
        object.__setattr__(self, 'history', [{'step': created_in, 'values': dict(kw)}])

    def setValuesInStep(self, stepName=None, **kw):
        self._log.add(self._name, 'setValuesInStep', (), dict(kw, stepName=stepName))
        self.history.append({'step': stepName, 'values': dict(kw)})

    def deactivate(self, stepName=None):
        self._log.add(self._name, 'deactivate', (stepName,), {})
        self.history.append({'step': stepName, 'values': {'deactivated': True}})


class Part(Obj):
    def __init__(self, name, log):
        Obj.__init__(self, name, log)

    def Set(self, **kw):
        self._log.add(self._name, 'Set', (), kw)
        return Obj('%s.Set(%s)' % (self._name, kw.get('name')), self._log)


class Instance(Obj):
    def __init__(self, name, log, nodes=None):
        Obj.__init__(self, name, log)
        object.__setattr__(self, 'nodes', nodes if nodes is not None else NodeArray([]))
        object.__setattr__(self, 'sets', {})
        self.sets['All'] = Obj('%s.sets[All]' % name, log)


class Assembly(Obj):
    def __init__(self, name, log):
        Obj.__init__(self, name, log)
        object.__setattr__(self, 'sets', {})
        object.__setattr__(self, 'referencePoints', {})
        object.__setattr__(self, '_rp', 5000)

    def Set(self, **kw):
        self._log.add(self._name, 'Set', (), kw)
        s = Obj('%s.Set(%s)' % (self._name, kw.get('name')), self._log)
        if kw.get('name') is not None:
            self.sets[kw['name']] = s
        return s

    def ReferencePoint(self, point=None):
        self._rp += 1
        self._log.add(self._name, 'ReferencePoint', (point,), {})
        rp = Obj('RP%d' % self._rp, self._log, id=self._rp, point=point)
        self.referencePoints[self._rp] = rp
        return rp

    def Instance(self, name=None, part=None, dependent=None):
        self._log.add(self._name, 'Instance', (), {'name': name, 'dependent': dependent})
        inst = Instance(name, self._log, nodes=NODES)
        return inst

    def rotate(self, **kw):
        self._log.add(self._name, 'rotate', (), kw)
        return Obj('%s.rotate' % self._name, self._log)

    def translate(self, **kw):
        self._log.add(self._name, 'translate', (), kw)
        return Obj('%s.translate' % self._name, self._log)

    def regenerate(self, **kw):
        self._log.add(self._name, 'regenerate', (), kw)
        return Obj('%s.regenerate' % self._name, self._log)

    def DatumCsysByDefault(self, *a, **kw):
        self._log.add(self._name, 'DatumCsysByDefault', a, kw)
        return Obj('%s.DatumCsysByDefault' % self._name, self._log)


class Model(Obj):
    def __init__(self, name, log):
        Obj.__init__(self, name, log)
        object.__setattr__(self, 'steps', {})
        object.__setattr__(self, 'boundaryConditions', {})
        object.__setattr__(self, 'concentratedForces', {})
        object.__setattr__(self, 'parts', {})
        object.__setattr__(self, 'rootAssembly', Assembly('%s.assembly' % name, log))

    def Part(self, name=None, dimensionality=None, type=None):
        self._log.add(self._name, 'Part', (), {'name': name})
        p = Part('%s.Part(%s)' % (self._name, name), self._log)
        self.parts[name] = p
        return p

    def StaticStep(self, name=None, **kw):
        self._log.add(self._name, 'StaticStep', (), dict(kw, name=name))
        self.steps[name] = Step('%s.steps[%s]' % (self._name, name), self._log, **kw)
        return self.steps[name]

    def BuckleStep(self, name=None, **kw):
        self._log.add(self._name, 'BuckleStep', (), dict(kw, name=name))
        self.steps[name] = Step('%s.steps[%s]' % (self._name, name), self._log, **kw)
        return self.steps[name]

    def DisplacementBC(self, name=None, **kw):
        self.boundaryConditions[name] = BC(name, self._log, kw.pop('createStepName', None), **kw)
        return self.boundaryConditions[name]

    def _reject_initial_load(self, kind, name, kw):
        """[법칙] Abaqus 는 **하중**을 Initial 스텝에 만들 수 없다.

        ValueError: The specified step either does not exist or is the Initial step.
        BC 는 Initial 에 만들어도 된다 — 그래서 이 검사는 하중 생성기에만 건다.
        (2026-10-07 실측: clamp_lf 배선이 Initial 에 CF 를 만들어 모델 생성이 죽었는데,
         관대한 목이 통과시켜 Abaqus 왕복 1회를 잃었다. 그래서 목에 법칙을 넣는다 —
         스텁 하네스의 존재 이유가 바로 이런 '라이선스 없이 잡히는 실수' 이다.)
        """
        if kw.get('createStepName') == 'Initial':
            raise ValueError(
                "[harness] %s %r 을 Initial 스텝에 만들 수 없습니다 (Abaqus: 'The specified "
                "step either does not exist or is the Initial step'). DEAD 값을 만들 때 주고 "
                "Step-Buckle 에서 setValuesInStep 으로 LIVE 를 주세요." % (kind, name))

    def ConcentratedForce(self, name=None, **kw):
        self._reject_initial_load('ConcentratedForce', name, kw)
        self.concentratedForces[name] = BC(name, self._log, kw.pop('createStepName', None), **kw)
        return self.concentratedForces[name]


class Mdb(Obj):
    def __init__(self, log):
        Obj.__init__(self, 'mdb', log)
        object.__setattr__(self, 'models', {})
        object.__setattr__(self, 'jobs', {})

    def Model(self, name=None, **kw):
        self._log.add('mdb', 'Model', (), {'name': name})
        m = Model(name, self._log)
        self.models[name] = m
        return m


# ---------------------------------------------------------------------------
# 스텁 모듈 + 실행
# ---------------------------------------------------------------------------


class _AnyCallable(object):
    """혼합대소문자 이름(`from mesh import ElemType` 등)을 위한 관대한 목.

    배선 검증에는 '호출되었는가'만 필요하므로 어떤 인자든 받아 기록용 Obj 를 돌려준다.
    """

    def __init__(self, name):
        self._name = name

    def __call__(self, *a, **kw):
        try:
            LOG.add('stub:%s' % self._name, 'call', a, kw)
        except NameError:
            pass
        return Obj('stub:%s()' % self._name, LOG)


def _shim_numpy():
    """numpy 가 없는 환경(에이전트 로컬)에서 스크립트가 쓰는 부분만 제공한다.

    run_abaqus_buckle.py 의 np 사용은 sin/cos/deg2rad/degrees/arctan2/pi/linalg.norm/array 뿐이다.
    사용자 머신(Abaqus python)에는 numpy 가 있으므로 이 심은 로컬 전용이고,
    이 하네스가 표준 라이브러리만으로 돌게 만든다(다른 도구들과 같은 방침).
    """
    import math
    m = types.ModuleType('numpy')

    class _Vec(list):
        def __truediv__(self, s):
            return _Vec([v / s for v in self])
        __div__ = __truediv__

    m.pi = math.pi
    m.sin = math.sin
    m.cos = math.cos
    m.deg2rad = math.radians
    m.degrees = math.degrees
    m.arctan2 = math.atan2
    m.array = lambda x: _Vec(x)
    m.linalg = types.ModuleType('numpy.linalg')
    m.linalg.norm = lambda v: math.sqrt(sum(float(x) * float(x) for x in v))
    return m


def _stub_modules(head):
    """`from abaqus import *` 등이 대문자 심볼을 전부 가져가므로 소스에서 이름을 긁어 채운다."""
    consts = sorted(set(re.findall(r'\b([A-Z][A-Z0-9_]{1,})\b', head)))
    mods = {}
    for mod in ('abaqus', 'abaqusConstants', 'step', 'mesh', 'regionToolset', 'interaction'):
        m = types.ModuleType(mod)
        m.__all__ = consts
        def _ga(k, _c=tuple(consts)):
            if k in _c:
                return k            # SUBSPACE / ON / DISTRIBUTING ... = 자기 이름 문자열
            return _AnyCallable(k)  # ElemType 등 혼합대소문자 이름
        m.__getattr__ = _ga
        mods[mod] = m
    # regionToolset.Region(...) 이 호출되므로 호출 가능한 목을 준다
    mods['regionToolset'].Region = lambda *a, **kw: Obj('region', LOG)
    return mods


def _substitute(head, case, legacy):
    """케이스 상수를 소스 텍스트에서 치환한다(= '상수 고치고 다시 돌리기' 의 모의). 치환 내역을 돌려준다."""
    subs = []
    if legacy:
        for name, val in (('CLAMP_MODE', case),):
            pat = re.compile(r'^%s\s*=\s*[\'"][^\'"]*[\'"]' % name, re.M)
            head, n = pat.subn("%s = '%s'" % (name, val), head, count=1)
            subs.append((name, val, n))
        pat = re.compile(r"^LOAD_MODE\s*=\s*['\"][^'\"]*['\"]", re.M)
        head, n = pat.subn("LOAD_MODE = 'cable'", head, count=1)
        subs.append(('LOAD_MODE', 'cable', n))
        return head, subs
    pat = re.compile(r"^CASE\s*=\s*['\"][^'\"]*['\"]", re.M)
    head, n = pat.subn("CASE = '%s'" % case, head, count=1)
    subs.append(('CASE', case, n))
    for name, val in sorted(CASE_CONSTANTS.get(case, {}).items()):
        pat = re.compile(r'^%s\s*=\s*[\'"]?[^\'"\s]+[\'"]?' % name, re.M)
        head, n2 = pat.subn('%s = %s' % (name, val), head, count=1)
        subs.append((name, val, n2))
    return head, subs


def main(argv):
    args = [a for a in argv[1:] if not a.startswith('--')]
    out_path = None
    for a in argv[1:]:
        if a.startswith('--out='):
            out_path = a.split('=', 1)[1]
    legacy = '--legacy' in argv
    if not args:
        sys.stderr.write(__doc__)
        return 2
    case = args[0]
    x_c = float(args[1]) if len(args) > 1 else 0.5
    if not legacy and case not in CASE_CONSTANTS:
        sys.stderr.write('알 수 없는 케이스 %r (가능: %s)\n' % (case, '|'.join(sorted(CASE_CONSTANTS))))
        return 2

    src = io.open(SRC, encoding='utf-8').read()
    head, sep, _ = src.partition(CUT_MARKER)
    if not sep:
        sys.stderr.write('FAIL 실행 블록 마커(%r)를 찾지 못했습니다 — 파일 구조가 바뀌었습니까?\n' % CUT_MARKER)
        return 3
    head, subs = _substitute(head, case, legacy)

    print('=' * 78)
    print('[harness] src   = %s' % SRC)
    print('[harness] mode  = %s   case=%s   x_c=%s' % ('legacy' if legacy else 'CASE', case, x_c))
    for name, val, n in subs:
        print('[harness] subst %-18s = %-10s (치환 %d회)%s'
              % (name, val, n, '' if n == 1 else '   <<< 치환 실패 — 상수 이름을 확인하라'))
    print('[harness] 주의: lambda/수렴/모드 형상은 이 도구가 볼 수 없다(배선만 본다).')
    print('=' * 78)

    global LOG, NODES
    LOG = Log()
    NODES = NodeArray(_lattice_nodes())
    print('[harness] 골격 메쉬 노드 %d개 (경계 위 노드 포함)' % len(NODES))

    g = {'__name__': '__stub__', '__file__': SRC}
    saved_argv, saved_cwd = list(sys.argv), os.getcwd()
    for name, m in _stub_modules(head).items():
        sys.modules[name] = m
    sys.path.insert(0, HERE)
    try:
        import numpy            # noqa: F401  (있으면 진짜 numpy 를 쓴다 — 사용자 머신과 동일 조건)
    except ImportError:
        sys.modules['numpy'] = _shim_numpy()
        print('[harness] numpy 없음 -> 내장 심(shim) 사용 (표준 라이브러리만)')
    import aba_grid_mesh as _gm            # 실제 모듈(순수 파이썬) — 좌표 판정은 Abaqus 불필요
    sys.modules['aba_grid_mesh'] = _gm
    g['mdb'] = Mdb(LOG)
    sys.argv = [SRC, str(x_c)]             # 스크립트의 parse_args 에 하네스 인자가 새지 않게
    try:
        exec(compile(head, SRC, 'exec'), g)
    except SystemExit:
        raise
    except Exception as e:
        sys.stderr.write('FAIL head 실행 실패(%s): %r\n' % (type(e).__name__, e))
        import traceback
        traceback.print_exc()
        return 4
    finally:
        sys.argv = saved_argv
        os.chdir(saved_cwd)

    build_model = g.get('build_model')
    if not callable(build_model):
        sys.stderr.write('FAIL build_model 을 찾지 못했습니다\n')
        return 5
    try:
        model_name = build_model(g['DISP'])
    except Exception as e:
        sys.stderr.write('FAIL build_model 예외(%s): %r\n' % (type(e).__name__, e))
        import traceback
        traceback.print_exc()
        _dump(LOG, out_path, case, x_c, ok=False)
        return 5

    print('[harness] build_model(%r) 반환 = %r  — 배선 기록 %d건'
          % (g['DISP'], model_name, len(LOG.rows)))
    _report(g, LOG)
    _dump(LOG, out_path, case, x_c, ok=True)
    return 0


def _report(g, log):
    """사람이 읽는 배선 요약 — 이게 하네스의 실제 산출물이다."""
    print('-' * 78)
    print('[배선] 경계조건 / 집중하중  (스텝별 값 이력)')
    for name, bc in sorted(g['mdb'].models[g['MODEL_PREFIX']].boundaryConditions.items()):
        for h in bc.history:
            print('  BC    %-20s %-20s %s'
                  % (name, h['step'], dict((k, _brief(v)) for k, v in h['values'].items())))
    for name, cf in sorted(g['mdb'].models[g['MODEL_PREFIX']].concentratedForces.items()):
        for h in cf.history:
            print('  CLOAD %-20s %-20s %s'
                  % (name, h['step'], dict((k, _brief(v)) for k, v in h['values'].items())))
    print('[배선] 스텝 정의 순서: %s'
          % ' -> '.join('%s(%s)' % (r['kwargs'].get('name'), r['call'])
                        for r in log.rows if r['call'] in ('StaticStep', 'BuckleStep')))
    print('[배선] 노드셋 생성: %s'
          % ', '.join(r['kwargs'].get('name') for r in log.rows
                      if r['call'] == 'Set' and r['kwargs'].get('name')))
    print('-' * 78)


def _dump(log, out_path, case, x_c, ok):
    if not out_path:
        return
    payload = {'case': case, 'x_c': x_c, 'build_model_ok': ok,
               'creates': [r for r in log.rows if r['call'] in ('CREATE', 'DisplacementBC',
                                                                'ConcentratedForce', 'Set', 'Tie',
                                                                'Coupling')],
               'all_calls': len(log.rows)}
    with io.open(out_path, 'w', encoding='utf-8') as f:
        f.write(json.dumps(payload, ensure_ascii=False, indent=1))
    print('[harness] 기록 저장: %s' % out_path)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
