#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""[§19.107] run_abaqus_buckle.py PULL_TOP 하네스 - Abaqus 없이 **실제 실행**한다.

무엇을 검증하나 (실해석이 아니라 '배선(BC/스텝 호출)' 검증이다)
    run_abaqus_buckle.py 의 구속조건 구역(`# ---- 구속조건 ----` ~ `return model_name`)을
    **소스에서 떼어** 가짜 my_model / 가짜 어셈블리로 그대로 실행하고, DisplacementBC 생성과
    setValuesInStep 호출을 기록해 판정한다. 텍스트 grep 이 아니라 실행이라서
    "분기 안에 넣었는데 실제로는 안 걸린다" 류의 오류가 잡힌다.

검증 항목
  (i)   회귀: PULL_TOP=False 의 BC 호출 추적 + 콘솔 로그가 **PULL_TOP 도입 이전 소스(git)** 와
        문자열 동일 (두 케이스: seed / clamp_only). Disp_Control_Top 없음, PULL_TOP 로그 0줄.
  (ii)  PULL_TOP=True:
        - BC_Anchor_Top 에 u1/u2 **없고** u3/ur1~ur3 만 (region=A1 그대로)
        - Disp_Control_Top: Initial(u1=u2=SET, region=A1) -> GlobalTension(0, +DISP)
        - Step-Buckle: seed 는 (0, +PERTURBATION)(세 꼭짓점 동일 크기), clamp_only 는 (0,0)(코너와 동일)
        - 로그 정확히 1줄, 실제 적용값(프리텐션) 포함, cp949 인코딩 가능
  (iii) 상수 단일 경로: 소스에 `SOLARSAIL_PULL_TOP` 문자열 **0회**(env 오버라이드 제거 확인),
        상수 블록을 os/environ 없이 **실제 exec** 해서 기본값 False / 로그 0줄 / 헬퍼변수 0,
        그리고 `PULL_TOP` 줄만 True 로 바꿔치기한 소스가 (ii) 의 True 동작과 **같은 trace** 를 내는지.
  (iv)  잡이름: False 는 종전 이름과 동일, True 는 '_3t' 가 붙고 80자 이하(두 케이스).
  (v)   방향: TOP_PULL_COS/SIN 이 소스의 Cable_Top 방향 (0,1) 과 같은 규약인지.
  (vi)  금지문자(em/en dash, 화살표, smart quote, NBSP, bullet, ellipsis) 0 - cp949 로그 안전.

사용
    python3 buckle_pull_top_stub.py        # PASS/FAIL, 실패 시 EXIT=1
Abaqus/라이선스 불필요(표준 라이브러리만).
종료코드: 0 OK / 1 FAIL / 2 구역 추출 실패
"""
from __future__ import print_function

import io
import math
import os
import re
import subprocess
import sys
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC_PATH = os.path.join(HERE, 'run_abaqus_buckle.py')
REL = 'code/run_abaqus_buckle.py'

MARK_A = u'    # ---- 구속조건 ----'
MARK_B = u'    return model_name'

FORBIDDEN_RE = re.compile(u'[\u2010\u2011\u2012\u2013\u2014\u2015'
                          u'\u2190\u2191\u2192\u2193\u21d2'
                          u'\u2018\u2019\u201c\u201d'
                          u'\u00a0\u2022\u2026]')

FAIL = []
N = [0]


def check(name, cond, detail=None):
    N[0] += 1
    ok = bool(cond)
    print("  %-4s %s%s" % ('OK' if ok else 'FAIL', name,
                           ('   <- ' + str(detail)) if (detail is not None and not ok) else ''))
    if not ok:
        FAIL.append(name)
    return ok


def read_source(path):
    with io.open(path, encoding='utf-8') as f:
        return f.read()


def git(args):
    return subprocess.check_output(['git', '-C', ROOT] + args).decode('utf-8', 'replace')


def baseline_source(cur):
    """PULL_TOP 도입 이전 소스. (text, ref) 반환."""
    head = git(['show', 'HEAD:%s' % REL])
    if 'PULL_TOP' not in head:
        return head, 'HEAD'
    try:
        hs = git(['log', '--reverse', '--format=%H', '-S', 'PULL_TOP', '--', REL]).split()
        if hs:
            parent = hs[0] + '^'
            return git(['show', '%s:%s' % (parent, REL)]), '%s^ (도입 커밋의 부모)' % hs[0][:8]
    except Exception:
        pass
    return git(['show', 'HEAD~1:%s' % REL]), 'HEAD~1'


# ------------------------------------------------------------------ 구역 떼기
def region_bc(src):
    a, b = src.find(MARK_A), src.find(MARK_B)
    if a < 0 or b < 0 or b <= a:
        return None
    return textwrap.dedent(src[a:b])


# ------------------------------------------------------------------ 가짜 모델
class FakeBC(object):
    def __init__(self, name, rec):
        self.name, self.rec = name, rec

    def setValuesInStep(self, **kw):
        self.rec.append(('setValuesInStep', self.name, list(kw.items())))

    def deactivate(self, step):
        self.rec.append(('deactivate', self.name, step))


class FakeModel(object):
    def __init__(self):
        self.rec = []
        self.boundaryConditions = _BCMap(self.rec)

    def DisplacementBC(self, **kw):
        self.rec.append(('DisplacementBC', list(kw.items())))
        nm = kw.get('name')
        if nm is not None:
            self.boundaryConditions[nm] = FakeBC(nm, self.rec)

    def trace(self):
        return '\n'.join(repr(r) for r in self.rec)


class _BCMap(dict):
    def __init__(self, rec):
        dict.__init__(self)
        self.rec = rec

    def __missing__(self, key):
        bc = FakeBC(key, self.rec)
        self[key] = bc
        return bc


class FakeSets(dict):
    """셋 이름 조회 -> 그 이름 문자열. 없는 키도 이름으로 만들어 준다."""

    def __missing__(self, key):
        self[key] = 'SET:%s' % key
        return self[key]


class FakeAssembly(object):
    def __init__(self):
        self.sets = FakeSets()


class FakeInstance(object):
    def __init__(self):
        self.nodes = list(range(1000))
        self.sets = FakeSets()


class FakeGridMesh(object):
    """aba_grid_mesh 의존 3함수만. 호출은 기록한다(구역 실행에 필요한 전부)."""

    def __init__(self):
        self.calls = []

    def clamp_exclude_labels(self, instance, centers, radius):
        self.calls.append(('clamp_exclude_labels', centers, radius))
        return set()

    def boundary_node_labels(self, instance, base=20.0, height=10.0, exclude_labels=None):
        self.calls.append(('boundary_node_labels', base, height, tuple(sorted(exclude_labels or ()))))
        return set(range(0, 100, 5))

    def make_set(self, assembly, name, instance, labels):
        self.calls.append(('make_set', name, len(labels)))
        assembly.sets[name] = 'SET:%s' % name


NS_BASE = {
    # abaqusConstants 의 SET 은 실제로는 심볼 -> 스텁에서는 문자열로 충분(값 비교만 한다).
    'SET': 'SET', 'TAG': '[run_abaqus_buckle]',
    'BASE': 20.0, 'HEIGHT': 10.0,
    'CLEAN_PERTURBATION': True,
    'DISP_TOP_OVER_CORNER': 1.4142135623,
    'PRETENSION_MODE': 'corner2',
    'PATTERN_SIGN': 1.0,
}

#   [§19.111] 자유변 면외 구속 스위치 3줄 - **소스에서 실행해** 얻는다(하네스가 값을 만들지 않는다).
#     False 경로는 RELEASE_EDGE_U3 만 보지만, 스프링 상수도 주입해 두어야 소스가 True 로 남은 채
#     실행해도 NameError 대신 'trace 불일치'로 정직하게 실패한다.
EDGE_LINE = re.compile(r'(?m)^(RELEASE_EDGE_U3|EDGE_SPRING_ALPHA|EDGE_SPRING_K) = ([^#\n]+)')


def edge_consts(src):
    pairs = [(m.group(1), m.group(2).strip()) for m in EDGE_LINE.finditer(src)]
    if not pairs:
        return {}
    ns = {}
    exec(compile('\n'.join('%s = %s' % (k, v) for k, v in pairs), '<edge_consts>', 'exec'), ns)
    return dict((k, ns[k]) for k, _ in pairs)



def make_ns(src, pt, case, logs, grid=None, scale=None):
    ang = math.radians(28.6)
    ns = dict(NS_BASE)
    ns.update({
        'my_model': FakeModel(),
        'a': FakeAssembly(),
        'inst_memb': FakeInstance(),
        'aba_grid_mesh': grid or FakeGridMesh(),
        'A1': 'A1', 'A2': 'A2', 'A3': 'A3',
        'end_cl': 'end_cl', 'end_cr': 'end_cr',
        'V_CL': (5.0, 5.0), 'V_CR': (15.0, 5.0),
        'CLAMP_EXCL_R': 0.21, 'CLAMP_PULL': 1.0e-6,
        'cos_val': math.cos(ang), 'sin_val': math.sin(ang),
        'DIR_CL': (-1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0)),
        'DIR_CR': (1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0)),
        '_clamp_target': lambda side: 'CLAMP:%s' % side,
        'PULL_TOP': pt,
        'TOP_PULL_COS': field(src, 'TOP_PULL_COS'),
        'TOP_PULL_SIN': field(src, 'TOP_PULL_SIN'),
        #   [§19.109] 소스의 값을 그대로 주입한다(없으면 None -> 구역이 참조하면 NameError 로 FAIL).
        'TOP_PULL_SCALE': (field(src, 'TOP_PULL_SCALE') if scale is None else float(scale)),
    })
    if case == 'seed':
        ns.update({'CASE': 'seed', 'HAS_CABLES': True, 'HAS_CLAMPS': False,
                   'HAS_CLAMP_CABLES': False, 'CLAMP_LOAD': 'none', 'LIVE_DRIVER': 'corner'})
    else:
        ns.update({'CASE': 'clamp_only', 'HAS_CABLES': False, 'HAS_CLAMPS': True,
                   'HAS_CLAMP_CABLES': True, 'CLAMP_LOAD': 'base', 'LIVE_DRIVER': 'clamp'})
        ns.update({'A1': 'RP_Top_Set', 'A2': 'RP_Right_Set', 'A3': 'RP_Left_Set'})
    ns['DISP'] = field(src, 'DISP_GLOBAL')
    ns['disp'] = ns['DISP']          # build_model(disp) 의 인자
    ns['PERTURBATION'] = field(src, 'PERTURBATION')
    #   [§19.111] 자유변 면외 구속 스위치(소스 값 그대로) + Abaqus 심볼.
    ns.update(edge_consts(src))
    ns.setdefault('ON', 'ON')
    ns.setdefault('OFF', 'OFF')

    def _print(*args):
        logs.append(' '.join(str(x) for x in args))

    ns['print'] = _print
    return ns


def field(src, name):
    m = re.search(r'(?m)^%s = ([^#\n]+)' % re.escape(name), src)
    if not m:
        return None
    try:
        return float(eval(m.group(1).strip(), {'__builtins__': {}}, {}))
    except Exception:
        return None


def run(region, src, pt, case, scale=None):
    logs = []
    ns = make_ns(src, pt, case, logs, scale=scale)
    exec(compile(region, '<buckle_bc:%s>' % case, 'exec'), ns)
    return ns['my_model'], logs, ns


def find_bc(rec, name):
    for r in rec:
        if r[0] == 'DisplacementBC':
            d = dict(r[1])
            if d.get('name') == name:
                return d
    return None


def find_steps(rec, name):
    out = {}
    for r in rec:
        if r[0] == 'setValuesInStep' and r[1] == name:
            d = dict(r[2])
            out[d.get('stepName')] = d
    return out


# ------------------------------------------------------------------ 상수 블록 exec
CONST_A = u'PULL_TOP = False'
CONST_B = u'TOP_PULL_SIN = -1.0'
PULL_TOP_LINE = re.compile(r'(?m)^PULL_TOP = (True|False)[ \t]*(?:#.*)?$')
#   [§19.109] TOP_PULL_SCALE 줄 - 사용자가 하는 편집과 같은 1줄 바꿔치기용.
SCALE_LINE = re.compile(r'(?m)^TOP_PULL_SCALE = ([^#\n]+)')
#   소스에 남아 있으면 안 되는 env 이름. 하네스는 이 값을 **읽지 않는다**(os.environ 미참조) --
#   여기서는 '소스에 0회' 를 세기 위한 needle 로만 쓴다.
ENV_NEEDLE = u'SOLARSAIL' u'_PULL_TOP'


def const_block(src):
    a, b = src.find(CONST_A), src.find(CONST_B)
    if a < 0 or b < 0 or b <= a:
        return None
    return src[a:b + len(CONST_B)]


def const_value(src):
    """소스의 `PULL_TOP = True|False` 줄을 **읽는다**(하네스가 값을 주입하지 않는다)."""
    m = PULL_TOP_LINE.search(src)
    return None if m is None else (m.group(1) == 'True')


def flip_constant(src, value):
    """사용자가 하는 편집과 같은 조작: 상수 줄 1줄만 바꿔치기한 소스를 돌려준다."""
    return PULL_TOP_LINE.sub('PULL_TOP = %s' % ('True' if value else 'False'), src, count=1)


def flip_scale(src, value):
    """[§19.109] TOP_PULL_SCALE 줄 1줄만 바꿔치기(사용자 편집과 같은 조작)."""
    return SCALE_LINE.sub('TOP_PULL_SCALE = %r' % float(value), src, count=1)


def scale_value(src):
    m = SCALE_LINE.search(src)
    if m is None:
        return None
    try:
        return float(m.group(1).strip())
    except Exception:
        return None


def exec_const(src):
    """상수 블록을 os/environ 없이 실제 실행 -> (PULL_TOP, logs, ns).

    블록이 아직 env 를 본다면 'os' 가 없어 NameError 로 **즉시 FAIL** 한다.
    """
    block = const_block(src)
    if block is None:
        return None, [], {}
    logs = []

    def _print(*a):
        logs.append(' '.join(str(x) for x in a))

    ns = {'TAG': '[run_abaqus_buckle]', 'print': _print}
    exec(compile(block, '<pull_top_consts>', 'exec'), ns)
    return ns['PULL_TOP'], logs, ns


# ------------------------------------------------------------------ 잡이름
def job_name(src, pt, case):
    """JOB_NAME 대입식을 실제로 eval 한다(포맷 인자 개수/태그 위치 검증 포함)."""
    a = src.find('JOB_NAME = ')
    b = src.find('if len(JOB_NAME) > 80', a)
    if a < 0 or b < 0:
        return None
    expr = src[a + len('JOB_NAME = '):b].strip()
    expr = expr.rstrip('\n').rstrip()
    ns = {'_value_tag': lambda v: ('%g' % float(v)).replace('+', '').replace('-', 'm')
          .replace('.', 'p'),
          'x_c': 0.5, 'DISP': field(src, 'DISP_GLOBAL'), '_dc_tag': '0p0556',
          '_f_tag': 100, 'ELEM_TAG': 's4', 'MESH_TAG': 'grid', 'CASE': case,
          'CLAMP_DIR': 'normal', '_excl_tag': '_ex210', '_rot_tag': '',
          '_live_tag': '_plv01000', '_pt_tag': ('_3t' if pt else '')}
    try:
        return eval(expr, {'__builtins__': {'int': int, 'round': round, 'float': float}}, ns)
    except Exception as e:
        return 'EVAL_FAIL:%r' % (e,)


def cable_top_dir(src):
    m = re.search(r"connect_cable\([^)]*'Cable_Top'[^)]*\)", src)
    if not m:
        return None
    t = re.findall(r'\((-?[0-9.]+),\s*(-?[0-9.]+),\s*-?[0-9.]+\)', m.group(0))
    return (float(t[-1][0]), float(t[-1][1])) if t else None


def _enc_ok(s):
    try:
        s.encode('cp949')
        return True
    except Exception:
        return False


# ------------------------------------------------------------------ 본체
def main():
    src = read_source(SRC_PATH)
    base, ref = baseline_source(src)
    print("[harness] PULL_TOP 하네스 - 소스 %s / 기저 %s" % (os.path.basename(SRC_PATH), ref))

    r_new, r_old = region_bc(src), region_bc(base)
    if r_new is None or r_old is None:
        print("[harness] FAIL: 구역 추출 실패 (마커: %r ~ %r)" % (MARK_A.strip(), MARK_B.strip()))
        return 2
    check('구역 추출: 현재/기저 (줄 수 %d / %d)' % (r_new.count('\n'), r_old.count('\n')),
          r_new.count('\n') > 20 and r_old.count('\n') > 20)

    DG = field(src, 'DISP_GLOBAL')
    PB = field(src, 'PERTURBATION')
    tc, ts = field(src, 'TOP_PULL_COS'), field(src, 'TOP_PULL_SIN')
    cdir = cable_top_dir(src)
    print("[harness] DISP_GLOBAL=%.6g PERTURBATION=%.6g  TOP_PULL_COS=%.6g TOP_PULL_SIN=%.6g"
          % (DG, PB, tc, ts))
    print("[harness] Cable_Top 방향(소스) = %s" % (cdir,))
    check('방향: (TOP_PULL_COS, -TOP_PULL_SIN) == Cable_Top 방향',
          cdir is not None and abs(tc - cdir[0]) < 1e-12 and abs(-ts - cdir[1]) < 1e-12,
          '%s vs %s' % ((tc, -ts), cdir))
    check('방향: Top 은 +y (u1=0, u2=+1)', abs(tc) < 1e-12 and abs(-ts - 1.0) < 1e-12, (tc, -ts))

    # ---------------------------------------------------------------- (i) 회귀
    print("\n[(i)] PULL_TOP=False 회귀 (기저 %s 대비)" % ref)
    TRACE_OLD = {}
    for case in ('seed', 'clamp_only'):
        mN, lgN, _ = run(r_new, src, False, case)
        mO, lgO, _ = run(r_old, base, False, case)
        TRACE_OLD[case] = mO.trace()
        check('%s: BC 호출 추적 문자열 동일 (%d 호출)' % (case, len(mN.rec)),
              mN.trace() == mO.trace(),
              '\n--- new ---\n' + mN.trace() + '\n--- base ---\n' + mO.trace())
        check('%s: 콘솔 로그 문자열 동일' % case, lgN == lgO, (lgN, lgO))
        check('%s: Disp_Control_Top 없음' % case,
              'Disp_Control_Top' not in mN.boundaryConditions, sorted(mN.boundaryConditions))
        check('%s: PULL_TOP 로그 0줄' % case, not any('PULL_TOP' in s for s in lgN), lgN)
        a = find_bc(mN.rec, 'BC_Anchor_Top')
        check('%s: BC_Anchor_Top 은 u1~ur3 완전 고정' % case,
              a and all(k in a and a[k] == 'SET'
                        for k in ('u1', 'u2', 'u3', 'ur1', 'ur2', 'ur3')), a)

    # ---------------------------------------------------------------- (ii) True
    print("\n[(ii)] PULL_TOP=True")
    TRACE_NEW, LOG_NEW = {}, {}
    for case in ('seed', 'clamp_only'):
        mT, lgT, ns = run(r_new, src, True, case)
        TRACE_NEW[case], LOG_NEW[case] = mT.trace(), lgT
        a = find_bc(mT.rec, 'BC_Anchor_Top')
        print("         [%s] BC_Anchor_Top = %s" % (case, list(a.items())))
        check('%s: BC_Anchor_Top 에 u1/u2 없음' % case, 'u1' not in a and 'u2' not in a, a)
        check('%s: BC_Anchor_Top 은 u3+회전 3개 유지' % case,
              all(k in a and a[k] == 'SET' for k in ('u3', 'ur1', 'ur2', 'ur3')), a)
        check('%s: BC_Anchor_Top region 은 A1 그대로' % case, a.get('region') == ns['A1'], a)

        t = find_bc(mT.rec, 'Disp_Control_Top')
        check('%s: Disp_Control_Top 생성 1회' % case, t is not None, None)
        if t:
            print("         [%s] Disp_Control_Top = %s" % (case, list(t.items())))
            check('%s: Initial u1=u2=SET, region=A1' % case,
                  t.get('u1') == 'SET' and t.get('u2') == 'SET'
                  and t.get('createStepName') == 'Initial' and t.get('region') == ns['A1'], t)
        steps = find_steps(mT.rec, 'Disp_Control_Top')
        for st in sorted(steps):
            print("         [%s] %-18s u1=%.9g u2=%.9g"
                  % (case, st, steps[st]['u1'], steps[st]['u2']))
        gt = steps.get('Step-GlobalTension', {})
        check('%s: GlobalTension = (0, +DISP)' % case,
              gt and abs(gt.get('u1', 9e9)) < 1e-18 and abs(gt.get('u2', 9e9) - DG) < 1e-18, gt)
        bk = steps.get('Step-Buckle', {})
        if case == 'seed':
            check('seed: Step-Buckle = (0, +PERTURBATION) (세 꼭짓점 동일 크기)',
                  bk and abs(bk.get('u1', 9e9)) < 1e-18
                  and abs(bk.get('u2', 9e9) - PB) < 1e-18, bk)
            rs = find_steps(mT.rec, 'Disp_Control_Right').get('Step-Buckle', {})
            check('seed: Right 와 Top 의 Buckle 변위 크기 동일',
                  abs(math.hypot(rs.get('u1', 0), rs.get('u2', 0))
                      - math.hypot(bk.get('u1', 0), bk.get('u2', 0))) < 1e-18, (rs, bk))
        else:
            check('clamp_only: Step-Buckle = (0,0) (코너와 같은 0 명시)',
                  bk and abs(bk.get('u1', 9e9)) < 1e-18 and abs(bk.get('u2', 9e9)) < 1e-18, bk)
            rs = find_steps(mT.rec, 'Disp_Control_Right').get('Step-Buckle', {})
            check('clamp_only: Right 도 (0,0) - Top 과 동일 처리',
                  abs(rs.get('u1', 9)) < 1e-18 and abs(rs.get('u2', 9)) < 1e-18, rs)
        pl = [s for s in lgT if 'PULL_TOP=True' in s]
        check('%s: PULL_TOP 로그 정확히 1줄' % case, len(pl) == 1, lgT)
        if pl:
            print("         [%s] LOG| %s" % (case, pl[0]))
            check('%s: 로그에 실제 적용값(프리텐션) 포함' % case,
                  ('%.4e' % DG) in pl[0], pl[0])
        check('%s: 로그 cp949 인코딩 가능(§19.95)' % case, all(_enc_ok(s) for s in lgT), lgT)

    # ---------------------------------------------------------------- (iii) 상수 단일 경로
    print("\n[(iii)] 상수 단일 경로 (env 오버라이드 없음)")
    check('소스에 %s 문자열 0회' % ENV_NEEDLE, src.count(ENV_NEEDLE) == 0, src.count(ENV_NEEDLE))

    blk = const_block(src)
    check('상수 블록 추출 (PULL_TOP = ... ~ TOP_PULL_SIN)', blk is not None)
    check('상수 블록에 os/environ 참조 0',
          blk is not None and 'os.' not in blk and 'environ' not in blk,
          None if blk is None else [l for l in blk.split('\n') if 'os.' in l or 'environ' in l])
    check('PULL_TOP 대입 정확히 1줄', len(PULL_TOP_LINE.findall(src)) == 1,
          PULL_TOP_LINE.findall(src))

    pt_c, lg_c, ns_c = exec_const(src)      # os 없이 실행: env 를 보면 NameError 로 즉시 FAIL
    check('기본 상수 False (os 없이 exec)', pt_c is False, (pt_c, lg_c))
    check('상수 해석은 로그 0줄', lg_c == [], lg_c)
    check('상수 해석이 헬퍼변수를 남기지 않음(_PULL_TOP* 0)',
          not [k for k in ns_c if k.startswith('_PULL_TOP')], sorted(ns_c))

    srcT = flip_constant(src, True)
    check('상수 줄 바꿔치기: True 로 1줄만 변경', const_value(srcT) is True and srcT != src)
    print("         True 로 바꾼 줄: %s" % [l for l in srcT.split('\n') if l.startswith('PULL_TOP = ')])
    for case in ('seed', 'clamp_only'):
        mC0, _, _ = run(r_new, src, const_value(src), case)
        check('%s: 상수 False 실행이 기저 trace 와 동일' % case,
              mC0.trace() == TRACE_OLD[case], mC0.trace())
        mC1, lgC1, _ = run(region_bc(srcT), srcT, const_value(srcT), case)
        check('%s: 상수만 True 로 바꿔치기 -> (ii) trace 와 동일' % case,
              mC1.trace() == TRACE_NEW[case],
              '\n--- flip ---\n' + mC1.trace() + '\n--- (ii) ---\n' + TRACE_NEW[case])
        check('%s: 상수 True 로그가 (ii) 와 동일' % case, lgC1 == LOG_NEW[case],
              (lgC1, LOG_NEW[case]))
        nT = job_name(srcT, const_value(srcT), case)
        check('%s: 상수 True 잡이름이 _3t' % case, bool(nT) and nT.endswith('_3t'), nT)

    # ---------------------------------------------------------------- (iv) 잡이름
    print("\n[(iv)] 잡이름 (두 런이 공존해야 비교가 된다)")
    for case in ('seed', 'clamp_only'):
        nF, nT = job_name(src, False, case), job_name(src, True, case)
        nO = job_name(base, False, case)
        print("         %-10s False=%s (%d자)" % (case, nF, len(nF)))
        print("         %-10s True =%s (%d자)" % (case, nT, len(nT)))
        check('%s: False 이름이 기저와 동일' % case, nF == nO, (nF, nO))
        check('%s: True 이름이 _3t 로 갈린다' % case, nT.endswith('_3t') and nT != nF, (nF, nT))
        check('%s: True 이름 80자 이하' % case, len(nT) <= 80, len(nT))

    # ---------------------------------------------------------------- (vii) 0.389
    print("\n[(vii)] TOP_PULL_SCALE = 0.389 (3-tendon 힘 평형) - Top 만 스케일, Right/Left 불변")
    _ang = field(src, 'CORNER_ANGLE_DEG')
    _lt, _lb = field(src, 'LEN_TOP'), field(src, 'LEN_BOT')
    eq_formula = 2.0 * math.sin(math.radians(_ang)) * (_lt / _lb)
    print("         유도: 2*sin(%.3f deg)*(LEN_TOP/LEN_BOT) = 2*%.6f*(%.3f/%.3f) = %.9f"
          % (_ang, math.sin(math.radians(_ang)), _lt, _lb, eq_formula))
    check('유도값이 0.389 와 일치 (|x-0.389| < 5e-4)', abs(eq_formula - 0.389) < 5e-4, eq_formula)
    check('기본 TOP_PULL_SCALE = 1.0 (평형 보정 없음)', scale_value(src) == 1.0, scale_value(src))
    check('TOP_PULL_SCALE 대입 정확히 1줄', len(SCALE_LINE.findall(src)) == 1,
          SCALE_LINE.findall(src))

    SQ = 0.389
    srcS = flip_scale(src, SQ)
    check('상수 줄 바꿔치기: 0.389 로 1줄만 변경',
          srcS != src and scale_value(srcS) == SQ and scale_value(src) == 1.0)
    print("         0.389 로 바꾼 줄: %s" % [l for l in srcS.split('\n')
                                              if l.startswith('TOP_PULL_SCALE = ')])
    for case in ('seed', 'clamp_only'):
        m1, lg1, _ = run(r_new, src, True, case)                 # 기본 1.0
        mS, lgS, _ = run(region_bc(srcS), srcS, True, case)      # 0.389
        s1 = find_steps(m1.rec, 'Disp_Control_Top')
        sS = find_steps(mS.rec, 'Disp_Control_Top')
        print("         [%s] GlobalTension/Buckle Top (1.0 -> 0.389)" % case)
        for st in ('Step-GlobalTension', 'Step-Buckle'):
            print("             %-18s u1 %.6e -> %.6e / u2 %.6e -> %.6e"
                  % (st, s1[st]['u1'], sS[st]['u1'], s1[st]['u2'], sS[st]['u2']))
        gt = sS.get('Step-GlobalTension', {})
        check('%s: GlobalTension Top u2 == 0.389 x (1.0 값) (절대오차 1e-12)' % case,
              abs(gt.get('u2', 9e9) - s1['Step-GlobalTension']['u2'] * SQ) < 1e-12,
              (gt.get('u2'), s1['Step-GlobalTension']['u2'] * SQ))
        check('%s: GlobalTension Top u1 == 0 (방향 불변)' % case, abs(gt.get('u1', 9e9)) < 1e-18, gt)
        bk = sS.get('Step-Buckle', {})
        if case == 'seed':
            check('seed: Buckle Top u2 == 0.389 x (1.0 값) (절대오차 1e-12)',
                  abs(bk.get('u2', 9e9) - s1['Step-Buckle']['u2'] * SQ) < 1e-12,
                  (bk.get('u2'), s1['Step-Buckle']['u2'] * SQ))
        else:
            check('clamp_only: Buckle Top 은 (0,0) 유지 (스케일과 무관)',
                  abs(bk.get('u1', 9)) < 1e-18 and abs(bk.get('u2', 9)) < 1e-18, bk)
        # Right/Left/클램프: 계수를 곱하지 않는다 -> 호출열이 1.0 런과 완전 동일
        _oth = ['Disp_Control_Right', 'Disp_Control_Left']
        if case == 'clamp_only':
            _oth += ['Disp_Control_CL', 'Disp_Control_CR', 'BC_Clamp_CL', 'BC_Clamp_CR']
        for nm in _oth:
            a = [r for r in m1.rec if r[0] in ('setValuesInStep', 'DisplacementBC')
                 and (r[1] == nm if r[0] == 'setValuesInStep' else dict(r[1]).get('name') == nm)]
            b = [r for r in mS.rec if r[0] in ('setValuesInStep', 'DisplacementBC')
                 and (r[1] == nm if r[0] == 'setValuesInStep' else dict(r[1]).get('name') == nm)]
            check('%s: %s 호출 불변 (%d회)' % (case, nm, len(a)), a == b and a, (a, b))
        # 로그: 1.0 은 추가 줄 0(바이트 동일), 0.389 는 정확히 1줄 + 실제 적용값
        check('%s: 1.0 로그에 TOP_PULL_SCALE 줄 0' % case,
              not any('TOP_PULL_SCALE' in s for s in lg1), lg1)
        _sc = [s for s in lgS if 'TOP_PULL_SCALE' in s]
        check('%s: 0.389 로그 정확히 1줄' % case, len(_sc) == 1, lgS)
        if _sc:
            print("         [%s] LOG| %s" % (case, _sc[0]))
            check('%s: 로그에 실제 적용값(%.6e) 포함' % (case, gt.get('u2', 0.0)),
                  ('%.6e' % gt.get('u2', 0.0)) in _sc[0], _sc[0])
            check('%s: 로그 cp949 인코딩 가능' % case, _enc_ok(_sc[0]), _sc[0])

    # ---------------------------------------------------------------- (vi) 금지문자
    print("\n[(vi)] 새 문자열 문자 검사")
    added = [l for l in src.split('\n')
             if 'PULL_TOP' in l or 'Disp_Control_Top' in l or 'TOP_PULL' in l or '_pt_tag' in l]
    bad = [l for l in added if FORBIDDEN_RE.search(l)]
    check('새 줄(%d개) 금지문자 0' % len(added), not bad, bad)

    print("\nRESULT: %s (%d checks, %d fail)" % ('PASS' if not FAIL else 'FAIL', N[0], len(FAIL)))
    if FAIL:
        for f in FAIL:
            print("  FAILED: %s" % f)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
