#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""[§19.111] RELEASE_EDGE_U3 하네스 - Abaqus 없이 **실제 실행**한다.

무엇을 검증하나 (실해석이 아니라 '배선(BC/스프링 호출)' 검증이다)
    run_abaqus.py / run_abaqus_buckle.py 에서 자유변 면외 구속 구역을 **소스에서 떼어**
    가짜 my_model / 가짜 어셈블리로 그대로 실행하고, DisplacementBC / SpringDashpotToGround
    호출을 기록해 판정한다. 텍스트 grep 이 아니라 실행이라서
    "분기 안에 넣었는데 실제로는 안 걸린다" 류의 오류가 잡힌다.

    구역
      H1  HF  `# 면외 z 구속용 노드셋(All_Edges_NoClamp)` ~ `#   [2026-10-12 §19.66]`
      H2  HF  `# Step-Buckle 삭제 시 BC_Edges_Only_Z` ~ `# 초기 상태에서 최종 장력` (Postbuckle 재생성)
      B1  buckle `# ---- Buckle 스텝: 논문 (c) 3변 u3=0 유지` ~ `if HAS_CABLES and LIVE_DRIVER`

    검증 항목
      (i)   기본값: 두 스크립트 모두 RELEASE_EDGE_U3 = False / ALPHA = 1e-3 / K = ALPHA*0.035.
      (ii)  False 회귀: 세 구역의 BC 호출 추적이 **도입 이전 소스(git)** 와 문자열 동일,
            스프링 호출 0회, RELEASE 로그 0줄.
      (iii) True 동작: BC_Edges_Only_Z **없음**, SpringDashpotToGround **1회**,
            dof=3 / springStiffness == EDGE_SPRING_K / region == All_Edges_NoClamp,
            로그 정확히 1줄 + n == 노드수.
      (iv)  면외 강체모드: HF 는 3점 RP(Top/Right/Left) u3=0 이 새로 걸리고,
            좌굴은 BC_Right_Z/BC_Left_Z(u3=SET) + BC_Anchor_Top(u3=SET) 로 이미 충분 -> 추가 0.
      (v)   alpha 스윕 1e-3 / 1e-1 / 1.0: springStiffness 가 정확히 비례(절대오차 1e-15).
      (vi)  HF Postbuckle 구역은 True 에서 BC_Edges_Only_Z 재생성 0회(스프링은 model-level).
      (vii) 새 문자열 금지문자(em/en dash, 화살표, smart quote, NBSP, bullet, ellipsis) 0,
            로그 cp949 인코딩 가능.

사용
    python3 release_edge_u3_stub.py        # PASS/FAIL, 실패 시 EXIT=1
Abaqus/라이선스 불필요(표준 라이브러리만).
종료코드: 0 OK / 1 FAIL / 2 구역 추출 실패
"""
from __future__ import print_function

import io
import os
import re
import subprocess
import sys
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
HF_PATH = os.path.join(HERE, 'run_abaqus.py')
BK_PATH = os.path.join(HERE, 'run_abaqus_buckle.py')
HF_REL = 'code/run_abaqus.py'
BK_REL = 'code/run_abaqus_buckle.py'
NEEDLE = 'RELEASE_EDGE_U3'

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


def baseline_source(rel):
    """RELEASE_EDGE_U3 도입 이전 소스. (text, ref) 반환. 커밋 전/후 모두 동작한다."""
    head = git(['show', 'HEAD:%s' % rel])
    if NEEDLE not in head:
        return head, 'HEAD'
    try:
        hs = git(['log', '--reverse', '--format=%H', '-S', NEEDLE, '--', rel]).split()
        if hs:
            return (git(['show', '%s^:%s' % (hs[0], rel)]),
                    '%s^ (도입 커밋의 부모)' % hs[0][:8])
    except Exception:
        pass
    return git(['show', 'HEAD~1:%s' % rel]), 'HEAD~1'


# ------------------------------------------------------------------ 가짜 모델 (호출 기록)
class FakeBC(object):
    def __init__(self, name, rec):
        self.name, self.rec = name, rec

    def setValuesInStep(self, **kw):
        self.rec.append(('setValuesInStep', self.name, list(kw.items())))

    def deactivate(self, step):
        self.rec.append(('deactivate', self.name, step))


class _BCMap(dict):
    def __init__(self, rec):
        dict.__init__(self)
        self.rec = rec

    def __missing__(self, key):
        bc = FakeBC(key, self.rec)
        self[key] = bc
        return bc


class FakeEngineeringFeatures(object):
    """Abaqus/CAE: rootAssembly.engineeringFeatures.SpringDashpotToGround(**kw).

    덱에서는 *ELEMENT, TYPE=SPRINGA(노드-지면 1개씩) + *SPRING, DOF=3 으로 나간다.
    """

    def __init__(self, rec):
        self.rec = rec

    def SpringDashpotToGround(self, **kw):
        self.rec.append(('SpringDashpotToGround', list(kw.items())))


class FakeAssembly(object):
    def __init__(self, rec=None):
        self.rec = rec if rec is not None else []
        self.sets = _SetMap()
        self.engineeringFeatures = FakeEngineeringFeatures(self.rec)


class _SetMap(dict):
    def __missing__(self, key):
        self[key] = 'SET:%s' % key
        return self[key]


class FakeModel(object):
    def __init__(self):
        self.rec = []
        self.boundaryConditions = _BCMap(self.rec)
        self.rootAssembly = FakeAssembly(self.rec)

    def DisplacementBC(self, **kw):
        self.rec.append(('DisplacementBC', list(kw.items())))
        nm = kw.get('name')
        if nm is not None:
            self.boundaryConditions[nm] = FakeBC(nm, self.rec)

    def trace(self):
        return '\n'.join(repr(r) for r in self.rec)


class FakeGridMesh(object):
    """aba_grid_mesh 의존 3함수만. 경계 노드 라벨은 고정(=스프링 대상 노드수 검증용)."""

    LABELS = tuple(range(1000, 1158))     # 158개 - 격자 경계 노드 규모

    def __init__(self):
        self.calls = []

    def clamp_exclude_labels(self, instance, centers, radius):
        self.calls.append(('clamp_exclude_labels', centers, radius))
        return set()

    def boundary_node_labels(self, instance, base=20.0, height=10.0, exclude_labels=None):
        self.calls.append(('boundary_node_labels', base, height, tuple(sorted(exclude_labels or ()))))
        return set(self.LABELS)

    def make_set(self, assembly, name, instance, labels):
        self.calls.append(('make_set', name, len(labels)))
        assembly.sets[name] = 'SET:%s' % name


class FakeInstance(object):
    def __init__(self):
        self.nodes = list(range(10000))
        self.sets = _SetMap()


# ------------------------------------------------------------------ 구역 떼기
HF_A = u'# 면외 z 구속용 노드셋(All_Edges_NoClamp)'
HF_B = u'#   [2026-10-12 §19.66]'
HF_PA = u'# Step-Buckle 삭제 시 BC_Edges_Only_Z'
HF_PB = u'# 초기 상태에서 최종 장력'
BK_A = u'# ---- Buckle 스텝: 논문 (c) 3변 u3=0 유지'
BK_B = u'if HAS_CABLES and LIVE_DRIVER !='


def region(src, mark_a, mark_b, dedent=False):
    a, b = src.find(mark_a), src.find(mark_b)
    if a < 0 or b < 0 or b <= a:
        return None
    if dedent:
        #   마커 문자열은 들여쓰기를 포함하지 않는다 -> 줄 시작으로 되돌린 뒤 dedent 해야
        #   첫 줄만 0칸으로 남아 dedent 가 무효화되는 일이 없다.
        a = src.rfind('\n', 0, a) + 1
        return textwrap.dedent(src[a:b])
    return src[a:b]


# ------------------------------------------------------------------ 상수 (소스에서 실행)
EDGE_LINE = re.compile(r'(?m)^(RELEASE_EDGE_U3|EDGE_SPRING_ALPHA|EDGE_SPRING_K) = ([^#\n]+)')


def edge_consts(src):
    pairs = [(m.group(1), m.group(2).strip()) for m in EDGE_LINE.finditer(src)]
    if len(pairs) != 3:
        return {}
    ns = {}
    exec(compile('\n'.join('%s = %s' % (k, v) for k, v in pairs), '<edge_consts>', 'exec'), ns)
    return dict((k, ns[k]) for k, _ in pairs)


def flip_flag(src, value):
    return re.sub(r'(?m)^RELEASE_EDGE_U3 = (?:True|False)', 'RELEASE_EDGE_U3 = %s' % value,
                  src, count=1)


def flip_alpha(src, value):
    return re.sub(r'(?m)^EDGE_SPRING_ALPHA = [^#\n]+', 'EDGE_SPRING_ALPHA = %r' % float(value),
                  src, count=1)


def field(src, name):
    m = re.search(r'(?m)^%s = ([^#\n]+)' % re.escape(name), src)
    if not m:
        return None
    try:
        return float(eval(m.group(1).strip(), {'__builtins__': {}}, {}))
    except Exception:
        return None


# ------------------------------------------------------------------ 구역 실행
def hf_ns(src):
    logs = []
    model = FakeModel()

    def _emit(*a):
        logs.append(' '.join(str(x) for x in a))

    ns = {
        '__name__': 'run_abaqus_bc', 'my_model': model, 'a': model.rootAssembly,
        'inst_memb': FakeInstance(), 'aba_grid_mesh': FakeGridMesh(),
        'BASE': field(src, 'BASE'), 'HEIGHT': field(src, 'HEIGHT'),
        'CLAMP_EXCL_R': field(src, 'CLAMP_EXCL_R'), 'NO_CLAMP': True,
        'V_CL': (10.0, 10.0), 'V_CR': (10.0, 10.0),
        'ON': 'ON', 'OFF': 'OFF', 'emit': _emit,
    }
    ns.update(edge_consts(src))
    return ns, logs


def bk_ns(src):
    logs = []
    model = FakeModel()

    def _print(*a):
        logs.append(' '.join(str(x) for x in a))

    ns = {
        '__name__': 'run_abaqus_buckle_bc', 'my_model': model, 'a': model.rootAssembly,
        'inst_memb': FakeInstance(), 'aba_grid_mesh': FakeGridMesh(),
        'BASE': field(src, 'BASE'), 'HEIGHT': field(src, 'HEIGHT'),
        'CLAMP_EXCL_R': field(src, 'CLAMP_EXCL_R'), 'HAS_CLAMPS': True,
        'V_CL': (10.0, 10.0), 'V_CR': (10.0, 10.0),
        'CLEAN_PERTURBATION': True, 'TAG': '[run_abaqus_buckle]',
        'ON': 'ON', 'OFF': 'OFF', 'print': _print,
    }
    ns.update(edge_consts(src))
    return ns, logs


def find_calls(rec, kind, name=None):
    out = []
    for r in rec:
        if r[0] != kind:
            continue
        d = dict(r[1])
        if kind == 'DisplacementBC' and name is not None and d.get('name') != name:
            continue
        if kind == 'SpringDashpotToGround' and name is not None and d.get('name') != name:
            continue
        out.append(d)
    return out


def main():
    hf = read_source(HF_PATH)
    bk = read_source(BK_PATH)
    print("[harness] RELEASE_EDGE_U3 하네스 - %s / %s"
          % (os.path.basename(HF_PATH), os.path.basename(BK_PATH)))

    r_h1, r_h2 = region(hf, HF_A, HF_B), region(hf, HF_PA, HF_PB, dedent=True)
    r_b1 = region(bk, BK_A, BK_B, dedent=True)
    check('구역 추출: HF 노드셋/BC', r_h1 is not None)
    check('구역 추출: HF Postbuckle 재생성', r_h2 is not None)
    check('구역 추출: buckle BC', r_b1 is not None)
    if FAIL:
        return 2
    hf0, hf_ref = baseline_source(HF_REL)
    bk0, bk_ref = baseline_source(BK_REL)
    r_h1_0, r_h2_0 = region(hf0, HF_A, HF_B), region(hf0, HF_PA, HF_PB, dedent=True)
    r_b1_0 = region(bk0, BK_A, BK_B, dedent=True)
    check('구역 추출: 기저(HF %s / buckle %s)' % (hf_ref, bk_ref),
          r_h1_0 is not None and r_h2_0 is not None and r_b1_0 is not None)

    # ---------------------------------------------------------------- (i) 기본값
    print("\n[(i)] 기본값 (상수 단일 경로)")
    ec, ecb = edge_consts(hf), edge_consts(bk)
    check('HF: 세 상수 대입 각 1줄', len(EDGE_LINE.findall(hf)) == 3, EDGE_LINE.findall(hf))
    check('buckle: 세 상수 대입 각 1줄', len(EDGE_LINE.findall(bk)) == 3, EDGE_LINE.findall(bk))
    check('HF: RELEASE_EDGE_U3 = False (기본)', ec.get('RELEASE_EDGE_U3') is False, ec)
    check('buckle: RELEASE_EDGE_U3 = False (기본)', ecb.get('RELEASE_EDGE_U3') is False, ecb)
    check('HF: EDGE_SPRING_ALPHA = 1e-3', ec.get('EDGE_SPRING_ALPHA') == 1.0e-3, ec)
    check('buckle: EDGE_SPRING_ALPHA = 1e-3', ecb.get('EDGE_SPRING_ALPHA') == 1.0e-3, ecb)
    check('HF: EDGE_SPRING_K == ALPHA * 0.035 (실행값)',
          ec.get('EDGE_SPRING_K') == ec.get('EDGE_SPRING_ALPHA') * 0.035, ec)
    check('buckle: EDGE_SPRING_K == ALPHA * 0.035 (실행값)',
          ecb.get('EDGE_SPRING_K') == ecb.get('EDGE_SPRING_ALPHA') * 0.035, ecb)
    check('상수 블록에 os/environ 참조 0',
          not any(k in s for s in (r_h1 or '', r_b1 or '') for k in ('os.environ', 'getenv')))

    # ---------------------------------------------------------------- (ii) False 회귀
    print("\n[(ii)] False 회귀 (기저 %s 대비)" % hf_ref)
    ns, lg = hf_ns(hf)
    exec(compile(r_h1, '<hf_edge>', 'exec'), ns)
    exec(compile(r_h2, '<hf_post>', 'exec'), ns)
    tr_hf = ns['my_model'].trace()
    ns0, lg0 = hf_ns(hf0)
    exec(compile(r_h1_0, '<hf_edge_base>', 'exec'), ns0)
    exec(compile(r_h2_0, '<hf_post_base>', 'exec'), ns0)
    tr_hf0 = ns0['my_model'].trace()
    check('HF: BC 호출 추적이 기저와 문자열 동일 (%d 호출)' % len(ns['my_model'].rec),
          tr_hf == tr_hf0, '\n--- new ---\n' + tr_hf + '\n--- base ---\n' + tr_hf0)
    rec = ns['my_model'].rec
    check('HF(False): BC_Edges_Only_Z 2회 (Step-Buckle + Postbuckle 재생성)',
          [d.get('createStepName') for d in find_calls(rec, 'DisplacementBC', 'BC_Edges_Only_Z')]
          == ['Step-Buckle', 'Step-Postbuckle'],
          find_calls(rec, 'DisplacementBC', 'BC_Edges_Only_Z'))
    check('HF(False): 스프링 호출 0회', not find_calls(rec, 'SpringDashpotToGround'))
    check('HF(False): RELEASE 로그 0줄', not any(NEEDLE in s for s in lg), lg)

    nsb, lgb = bk_ns(bk)
    exec(compile(r_b1, '<bk_edge>', 'exec'), nsb)
    tr_bk = nsb['my_model'].trace()
    nsb0, lgb0 = bk_ns(bk0)
    exec(compile(r_b1_0, '<bk_edge_base>', 'exec'), nsb0)
    tr_bk0 = nsb0['my_model'].trace()
    check('buckle: BC 호출 추적이 기저와 문자열 동일 (%d 호출)' % len(nsb['my_model'].rec),
          tr_bk == tr_bk0, '\n--- new ---\n' + tr_bk + '\n--- base ---\n' + tr_bk0)
    check('buckle(False): BC_Edges_Only_Z 1회 (Initial)',
          [d.get('createStepName')
           for d in find_calls(nsb['my_model'].rec, 'DisplacementBC', 'BC_Edges_Only_Z')]
          == ['Initial'], find_calls(nsb['my_model'].rec, 'DisplacementBC', 'BC_Edges_Only_Z'))
    check('buckle(False): 스프링 0회', not find_calls(nsb['my_model'].rec, 'SpringDashpotToGround'))
    check('buckle(False): RELEASE 로그 0줄', not any(NEEDLE in s for s in lgb), lgb)
    print("         False HF trace   | %s" % tr_hf.replace('\n', ' || '))
    print("         False buckle trace| %s" % tr_bk.replace('\n', ' || '))

    # ---------------------------------------------------------------- (iii) True 동작
    print("\n[(iii)] True 동작")
    hfT = flip_flag(hf, 'True')
    check('상수 줄 바꿔치기: HF True 로 1줄만 변경', hfT != hf and edge_consts(hfT).get('RELEASE_EDGE_U3') is True)
    nsT, lgT = hf_ns(hfT)
    exec(compile(region(hfT, HF_A, HF_B), '<hf_edge_true>', 'exec'), nsT)
    exec(compile(region(hfT, HF_PA, HF_PB, dedent=True), '<hf_post_true>', 'exec'), nsT)
    recT = nsT['my_model'].rec
    spring = find_calls(recT, 'SpringDashpotToGround')
    print("         True HF trace    | %s" % nsT['my_model'].trace().replace('\n', ' || '))
    check('HF(True): BC_Edges_Only_Z 0회', not find_calls(recT, 'DisplacementBC', 'BC_Edges_Only_Z'))
    check('HF(True): SpringDashpotToGround 1회', len(spring) == 1, spring)
    if spring:
        s = spring[0]
        print("         spring = %s" % list(s.items()))
        check('HF(True): name=EdgeZ_Spring', s.get('name') == 'EdgeZ_Spring', s.get('name'))
        check('HF(True): dof=3', s.get('dof') == 3, s.get('dof'))
        check('HF(True): springStiffness == EDGE_SPRING_K (정확히)',
              s.get('springStiffness') == edge_consts(hfT).get('EDGE_SPRING_K'), s.get('springStiffness'))
        check('HF(True): region == All_Edges_NoClamp',
              s.get('region') == 'SET:All_Edges_NoClamp', s.get('region'))
        check('HF(True): orientation=None / springBehavior=ON / dashpotBehavior=OFF',
              s.get('orientation') is None and s.get('springBehavior') == 'ON'
              and s.get('dashpotBehavior') == 'OFF', s)
    n_labels = len(FakeGridMesh.LABELS)
    re_ok = [s for s in lgT if NEEDLE in s]
    check('HF(True): RELEASE 로그 정확히 1줄', len(re_ok) == 1, lgT)
    if re_ok:
        print("         LOG| %s" % re_ok[0])
        check('HF(True): 로그 n == 스프링 대상 노드수 (%d)' % n_labels,
              ('n=%d' % n_labels) in re_ok[0], re_ok[0])
        check('HF(True): 로그에 k/alpha 값 포함',
              ('%.3e' % edge_consts(hfT)['EDGE_SPRING_K']) in re_ok[0]
              and ('%.1e' % edge_consts(hfT)['EDGE_SPRING_ALPHA']) in re_ok[0], re_ok[0])
        check('HF(True): 로그 cp949 인코딩 가능', _enc_ok(re_ok[0]), re_ok[0])
    check('HF(True): Postbuckle 구역이 BC_Edges_Only_Z 재생성 0회 (스프링은 model-level)',
          not find_calls(recT, 'DisplacementBC', 'BC_Edges_Only_Z')
          and len(find_calls(recT, 'SpringDashpotToGround')) == 1)

    # 3점 u3 (HF): 코너 RP 에 u3=0
    corner = dict((d.get('name'), d) for d in find_calls(recT, 'DisplacementBC')
                  if str(d.get('name', '')).endswith('_Corner_Z'))
    print("         3점 u3 BC = %s" % sorted(corner))
    for nm, rp in (('BC_Top_Corner_Z', 'RP_Top_Set'),
                   ('BC_Right_Corner_Z', 'RP_Right_Set'),
                   ('BC_Left_Corner_Z', 'RP_Left_Set')):
        d = corner.get(nm)
        check('HF(True): %s region=%s u3=0 Initial' % (nm, rp),
              d is not None and d.get('u3') == 0 and d.get('createStepName') == 'Initial'
              and d.get('region') == 'SET:%s' % rp, d)
    check('HF(False): 3점 u3 BC 0회 (종전 배선 불변)',
          not [k for k in ns['my_model'].boundaryConditions if str(k).endswith('_Corner_Z')])

    # 좌굴: 이미 있는 코너 u3 (추가 0)
    bkT = flip_flag(bk, 'True')
    nsbT, lgbT = bk_ns(bkT)
    exec(compile(region(bkT, BK_A, BK_B, dedent=True), '<bk_edge_true>', 'exec'), nsbT)
    recbT = nsbT['my_model'].rec
    print("         True buckle trace| %s" % nsbT['my_model'].trace().replace('\n', ' || '))
    springb = find_calls(recbT, 'SpringDashpotToGround')
    check('buckle(True): BC_Edges_Only_Z 0회',
          not find_calls(recbT, 'DisplacementBC', 'BC_Edges_Only_Z'))
    check('buckle(True): SpringDashpotToGround 1회', len(springb) == 1, springb)
    if springb:
        s = springb[0]
        check('buckle(True): dof=3 / k == EDGE_SPRING_K / region == All_Edges_NoClamp',
              s.get('dof') == 3 and s.get('springStiffness') == edge_consts(bkT)['EDGE_SPRING_K']
              and s.get('region') == 'SET:All_Edges_NoClamp', s)
    re_b = [s for s in lgbT if NEEDLE in s]
    check('buckle(True): RELEASE 로그 정확히 1줄', len(re_b) == 1, lgbT)
    if re_b:
        print("         LOG| %s" % re_b[0])
        check('buckle(True): 로그 n == 스프링 대상 노드수 (%d)' % n_labels,
              ('n=%d' % n_labels) in re_b[0], re_b[0])
    check('buckle(True): 코너 u3 추가 구속 0회 (BC_Right_Z/BC_Left_Z 가 이미 u3=SET)',
          not [k for k in nsbT['my_model'].boundaryConditions
               if str(k).endswith('_Corner_Z')])

    # ---------------------------------------------------------------- (iv) 기존 u3 근거
    print("\n[(iv)] Right/Left u3 구속 근거 (소스 문자열)")
    check('buckle: BC_Right_Z/BC_Left_Z 에 u3=SET',
          len(re.findall(r"name='BC_Right_Z'.*?u3=SET", bk, re.S)) >= 1
          and len(re.findall(r"name='BC_Left_Z'.*?u3=SET", bk, re.S)) >= 1)
    check('buckle: BC_Anchor_Top 에 u3=SET',
          len(re.findall(r"name='BC_Anchor_Top'.*?u3=SET", bk, re.S)) >= 1)
    check('HF: BC_Right/BC_Left(u3=0)는 케이블 바깥 끝(end_c2/end_c3)에 걸린다',
          len(re.findall(r"name='BC_Right'.*?region=end_c2", hf, re.S)) >= 1
          and len(re.findall(r"name='BC_Left'.*?region=end_c3", hf, re.S)) >= 1)
    check('HF: 코너 RP u3=0 은 기저 소스에 없음 (True 에서 새로 추가되는 구속)',
          'RP_Right_Set' in hf and 'RP_Right_Set' not in
          ''.join(l for l in hf0.split('\n') if 'DisplacementBC' in l))

    # ---------------------------------------------------------------- (v) alpha 스윕
    print("\n[(v)] alpha 스윕 (1e-3 / 1e-1 / 1.0) - springStiffness 비례성")
    ks = []
    for alpha in (1.0e-3, 1.0e-1, 1.0):
        srcA = flip_flag(flip_alpha(hf, alpha), 'True')
        ca = edge_consts(srcA)
        check('alpha=%g: K == alpha*0.035 (소스 실행값)' % alpha,
              ca.get('EDGE_SPRING_K') == alpha * 0.035, ca)
        nsA, _lgA = hf_ns(srcA)
        exec(compile(region(srcA, HF_A, HF_B), '<hf_edge_alpha>', 'exec'), nsA)
        sp = find_calls(nsA['my_model'].rec, 'SpringDashpotToGround')
        got = sp[0].get('springStiffness') if sp else None
        print("         alpha=%-6g -> k=%.12e (호출값 %.12e)" % (alpha, alpha * 0.035, got))
        check('alpha=%g: springStiffness == EDGE_SPRING_K (절대오차 1e-15)'
              % alpha, got is not None and abs(got - alpha * 0.035) <= 1e-15, got)
        ks.append(got)
    if len(ks) == 3 and all(k for k in ks):
        check('스윕 비례: k(1e-1) == 100 x k(1e-3) (절대오차 1e-15)',
              abs(ks[1] - 100.0 * ks[0]) <= 1e-15, (ks[1], 100.0 * ks[0]))
        check('스윕 비례: k(1.0) == 1000 x k(1e-3) (절대오차 1e-15)',
              abs(ks[2] - 1000.0 * ks[0]) <= 1e-15, (ks[2], 1000.0 * ks[0]))
        check('스윕 단조 증가 (k1e-3 < k1e-1 < k1.0)', ks[0] < ks[1] < ks[2], ks)
    else:
        check('스윕 값 수집', False, ks)

    # ---------------------------------------------------------------- (vii) 문자 검사
    print("\n[(vii)] 새 문자열 문자 검사")
    added = [l for l in (hf.split('\n') + bk.split('\n'))
             if NEEDLE in l or 'EDGE_SPRING' in l or 'EdgeZ_Spring' in l or '_Corner_Z' in l]
    bad = [l for l in added if FORBIDDEN_RE.search(l)]
    check('새 줄(%d개) 금지문자 0' % len(added), not bad, bad)
    check('로그 문구가 두 스크립트에서 동일',
          ('[BC] RELEASE_EDGE_U3=True -> BC_Edges_Only_Z 제거, EdgeZ_Spring ' in hf)
          and ('[BC] RELEASE_EDGE_U3=True -> BC_Edges_Only_Z 제거, EdgeZ_Spring ' in bk))

    print("\nRESULT: %s (%d checks, %d fail)" % ('PASS' if not FAIL else 'FAIL', N[0], len(FAIL)))
    if FAIL:
        for f in FAIL:
            print("  FAILED: %s" % f)
        return 1
    return 0


def _enc_ok(s):
    try:
        s.encode('cp949')
        return True
    except Exception:
        return False


if __name__ == '__main__':
    sys.exit(main())
