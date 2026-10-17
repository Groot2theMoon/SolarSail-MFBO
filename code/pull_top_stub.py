#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""[§19.105] PULL_TOP 하네스 - Abaqus 없이 **실제 실행**한다.

무엇을 검증하나 (실제 해석이 아니라 '배선(BC/스텝 호출)' 검증이다)
    run_abaqus.py 에서 세 구역을 **소스에서 떼어** 가짜 my_model 로 그대로 실행하고,
    BoundaryCondition 생성/스텝별 setValuesInStep 호출을 기록해 판정한다:

      R1  구속조건 구역 (`# 구속조건 생성` ~ `# Step Buckle 에서는`)
      R2  최종하중 구역 1 (`Step-HighTension` - LF 경로)
      R3  최종하중 구역 2 (`Step-Postbuckle` - HF 경로)

    검증 항목:
      (i)   PULL_TOP=False: BC 호출 추적이 **PULL_TOP 도입 이전 소스(git)** 와 문자열 동일,
            `Disp_Control_Top` 없음, 새 로그 0줄(회귀).
      (ii)  PULL_TOP=True: BC_Anchor 에 u1/u2 **없고** u3/ur1~ur3 만 남는다,
            `Disp_Control_Top` 이 Initial u1=u2=0 -> GlobalTension DISP_GLOBAL ->
            Buckle PERTURBATION -> HighTension/Postbuckle GLOBAL_FINAL 로 구동,
            로그 1줄(+방향 1줄) 존재, cp949 인코딩 가능.
      (iii) 네 조합(PULL_TOP x NO_CLAMP) 모두 예외 없이 통과.
      (iv)  Top 방향(top_cos, -top_sin) 이 소스의 Cable_Top 케이블 방향과 일치.
      (v)   새 로그에 금지문자(em-dash/en-dash/화살표/smart quote/NBSP/bullet/ellipsis) 0.

사용
    python3 pull_top_stub.py        # PASS/FAIL, 실패 시 EXIT=1
Abaqus/라이선스 불필요(표준 라이브러리만).
종료코드: 0 OK / 1 FAIL
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
RA = os.path.join(HERE, 'run_abaqus.py')
REL = 'code/run_abaqus.py'

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


# ------------------------------------------------------------------ 가짜 모델 (호출 기록)
class FakeBC(object):
    def __init__(self, name, rec):
        self.name = name
        self.rec = rec

    def setValuesInStep(self, **kw):
        self.rec.append(('setValuesInStep', self.name, list(kw.items())))

    def deactivate(self, step):
        self.rec.append(('deactivate', self.name, step))


class FakeModel(object):
    """DisplacementBC(...)/boundaryConditions[name].setValuesInStep(...) 만 기록한다."""

    def __init__(self):
        self.rec = []
        self.boundaryConditions = _BCMap(self.rec)
        self.fieldOutputRequests = {'F-Output-1': FakeBC('F-Output-1', self.rec)}

    def DisplacementBC(self, **kw):
        self.rec.append(('DisplacementBC', list(kw.items())))
        nm = kw.get('name')
        if nm is not None:
            self.boundaryConditions[nm] = FakeBC(nm, self.rec)

    def StaticStep(self, **kw):
        self.rec.append(('StaticStep', list(kw.items())))

    def trace(self):
        return '\n'.join(repr(r) for r in self.rec)


class _BCMap(dict):
    """없는 이름을 조회하면 즉시 등록한다(최종하중 구역만 떼어 실행할 때 필요)."""

    def __init__(self, rec):
        dict.__init__(self)
        self.rec = rec

    def __missing__(self, key):
        bc = FakeBC(key, self.rec)
        self[key] = bc
        return bc


# ------------------------------------------------------------------ 소스 구역 떼기
MARK_BC_A = u'# 구속조건 생성'
MARK_BC_B = u'# Step Buckle 에서는 전체 z 구속에서'


def region_bc(src):
    a = src.find(MARK_BC_A)
    b = src.find(MARK_BC_B)
    if a < 0 or b < 0 or b <= a:
        return None
    return src[a:b]


def region_drive(src, idx):
    """i번째 최종하중 구역: 'for name, sign in [(Right...' 줄 ~ 다음 fieldOutputRequests 줄 직전."""
    hits = [m.start() for m in re.finditer(r'(?m)^\s*for name, sign in \[\('
                                            r"'Disp_Control_Right', 1\)", src)]
    if len(hits) <= idx:
        return None
    a = hits[idx]
    f = src.find("my_model.fieldOutputRequests", a)
    if f < 0:
        return None
    return textwrap.dedent(src[a:f])


def field(src, name):
    """모듈 수준 상수 값을 소스에서 뽑아 float 으로 (주석 제거)."""
    m = re.search(r'(?m)^%s = ([^#\n]+)' % re.escape(name), src)
    if not m:
        return None
    try:
        return float(eval(m.group(1).strip(), {'__builtins__': {}}, {}))
    except Exception:
        return None


def cable_top_dir(src):
    """소스의 Cable_Top connect_cable 방향 튜플 (u1, u2)."""
    m = re.search(r"connect_cable\(\s*'Cable_Top'[^)]*\)", src)
    if not m:
        return None
    t = re.findall(r'\((-?[0-9.]+),\s*(-?[0-9.]+),\s*-?[0-9.]+\)', m.group(0))
    if not t:
        return None
    return (float(t[-1][0]), float(t[-1][1]))


# ------------------------------------------------------------------ 구역 실행
def make_ns(src, pt, nc):
    logs = []

    def _emit(*a):
        logs.append(' '.join(str(x) for x in a))

    ang = math.radians(28.6)
    ns = {
        '__name__': 'run_abaqus_bc', 'math': math,
        'my_model': FakeModel(),
        'NO_CLAMP': nc, 'PULL_TOP': pt,
        'end_c1': 'end_c1', 'end_c2': 'end_c2', 'end_c3': 'end_c3',
        'end_cl': 'end_cl', 'end_cr': 'end_cr',
        'DISP_GLOBAL': field(src, 'DISP_GLOBAL'),
        'GLOBAL_FINAL': field(src, 'GLOBAL_FINAL'),
        'PERTURBATION': field(src, 'PERTURBATION'),
        'top_cos': field(src, 'top_cos'), 'top_sin': field(src, 'top_sin'),
        'cos_val': math.cos(ang), 'sin_val': math.sin(ang),
        'emit': _emit,
    }
    # 소스에 있는 CLAMP_* 를 그대로 (없으면 계산)
    for k in ('CLAMP_PULL', 'CLAMP_FINAL', 'CLAMP_PERT'):
        v = field(src, k)
        if v is not None:
            ns[k] = v
    ns.setdefault('CLAMP_PULL', ns['DISP_GLOBAL'] * 0.0556)
    ns.setdefault('CLAMP_FINAL', ns['GLOBAL_FINAL'] * 0.0556)
    ns.setdefault('CLAMP_PERT', ns['PERTURBATION'] * 0.0556)
    ns['sqrt2'] = math.sqrt(2.0)
    return ns, logs


def run(region, src, pt, nc):
    ns, logs = make_ns(src, pt, nc)
    exec(compile(region, '<run_abaqus_bc>', 'exec'), ns)
    return ns['my_model'], logs


# ------------------------------------------------------------------ 기저(baseline) 소스
def git(args):
    return subprocess.check_output(['git', '-C', ROOT] + args).decode('utf-8', 'replace')


def baseline_source(cur):
    """PULL_TOP 도입 이전 소스. (text, ref) 반환."""
    head = git(['show', 'HEAD:%s' % REL])
    if 'PULL_TOP' not in head:
        return head, 'HEAD'
    try:
        # -S 는 추가/삭제 커밋을 모두 낸다 -> 가장 오래된 것이 도입 커밋
        hs = git(['log', '--reverse', '--format=%H', '-S', 'PULL_TOP', '--', REL]).split()
        if hs:
            parent = hs[0] + '^'
            return (git(['show', '%s:%s' % (parent, REL)]),
                    '%s^ (PULL_TOP 도입 커밋의 부모)' % hs[0][:8])
    except Exception:
        pass
    return git(['show', 'HEAD~1:%s' % REL]), 'HEAD~1'


def find(rec, kind, name=None):
    out = []
    for r in rec:
        if r[0] != kind:
            continue
        if kind == 'DisplacementBC':
            d = dict(r[1])
            if d.get('name') == name or name is None:
                out.append(d)
        else:
            # setValuesInStep: ('setValuesInStep', bc_name, [(k, v), ...]) -> 키는 stepName
            d = dict(r[2])
            if r[1] == name or name is None:
                out.append((d.get('stepName'), d))
    return out


# ------------------------------------------------------------------ 본체
def main():
    src = read_source(RA)
    base, ref = baseline_source(src)
    print("[harness] PULL_TOP 하네스 - 소스 %s / 기저 %s" % (os.path.basename(RA), ref))

    r_bc, r_bc0 = region_bc(src), region_bc(base)
    check('구역 추출: BC 구역(현재/기저)', r_bc is not None and r_bc0 is not None)
    r_fin = region_drive(src, 0)
    r_post = region_drive(src, 1)
    r_fin0 = region_drive(base, 0)
    r_post0 = region_drive(base, 1)
    check('구역 추출: HighTension/Postbuckle(현재)', r_fin is not None and r_post is not None)
    check('구역 추출: HighTension/Postbuckle(기저)', r_fin0 is not None and r_post0 is not None)
    if FAIL:
        return 1

    DG, GF, PB = field(src, 'DISP_GLOBAL'), field(src, 'GLOBAL_FINAL'), field(src, 'PERTURBATION')
    tc, ts = field(src, 'top_cos'), field(src, 'top_sin')
    cdir = cable_top_dir(src)
    print("[harness] DISP_GLOBAL=%.6g GLOBAL_FINAL=%.6g PERTURBATION=%.6g" % (DG, GF, PB))
    print("[harness] top_cos=%.6f top_sin=%.6f  (u1=v*cos, u2=-v*sin) -> (u1,u2)=(%.3f,%.3f)"
          % (tc, ts, tc, -ts))
    print("[harness] Cable_Top 방향(소스) = %s" % (cdir,))

    # (iv) 방향 코사인이 케이블 방향과 일치
    check('방향: (top_cos, -top_sin) == Cable_Top 방향',
          cdir is not None and abs(tc - cdir[0]) < 1e-12 and abs(ts + cdir[1]) < 1e-12,
          '%s vs %s' % ((tc, -ts), cdir))
    check('방향: Top 은 +y (u1=0, u2=+1)', abs(tc) < 1e-12 and abs(-ts - 1.0) < 1e-12,
          (tc, -ts))

    # ---------------------------------------------------------------- (i) 회귀
    print("\n[(i)] PULL_TOP=False 회귀 (기저 %s 대비)" % ref)
    mF, lgF = run(r_bc, src, False, False)
    mB, lgB = run(r_bc0, base, False, False)
    check('BC_Anchor/Right/Left/CL/CR 추적 바이트 동일', mF.trace() == mB.trace(),
          '\n--- new ---\n' + mF.trace() + '\n--- base ---\n' + mB.trace())
    check('Disp_Control_Top 없음', 'Disp_Control_Top' not in mF.boundaryConditions,
          list(mF.boundaryConditions))
    check('새 로그 0줄 (PULL_TOP 로그 없음)',
          lgF == lgB and not any('PULL_TOP' in s for s in lgF), lgF)
    aF = find(mF.rec, 'DisplacementBC', 'BC_Anchor')[0]
    print("         BC_Anchor(False) = %s" % (list(aF.items()),))
    check('BC_Anchor(False) 는 u1~ur3 완전 고정 (u1/u2/u3/ur1..3=0)',
          all(k in aF and aF[k] == 0 for k in ('u1', 'u2', 'u3', 'ur1', 'ur2', 'ur3')), aF)

    # ---------------------------------------------------------------- (ii) PULL_TOP=True
    print("\n[(ii)] PULL_TOP=True")
    mT, lgT = run(r_bc, src, True, False)
    aT = find(mT.rec, 'DisplacementBC', 'BC_Anchor')[0]
    print("         BC_Anchor(True)  = %s" % (list(aT.items()),))
    check('BC_Anchor: u1/u2 없음', 'u1' not in aT and 'u2' not in aT, aT)
    check('BC_Anchor: u3=0 + 회전 3개 유지',
          all(k in aT and aT[k] == 0 for k in ('u3', 'ur1', 'ur2', 'ur3')), aT)
    check('BC_Anchor: region 은 end_c1 그대로', aT.get('region') == 'end_c1', aT)
    top = find(mT.rec, 'DisplacementBC', 'Disp_Control_Top')
    check('Disp_Control_Top 생성 1회', len(top) == 1, len(top))
    if top:
        print("         Disp_Control_Top(Initial) = %s" % (list(top[0].items()),))
        check('Disp_Control_Top Initial: u1=0,u2=0, step=Initial, region=end_c1',
              top[0].get('u1') == 0 and top[0].get('u2') == 0
              and top[0].get('createStepName') == 'Initial'
              and top[0].get('region') == 'end_c1', top[0])
    steps = dict((s, d) for s, d in find(mT.rec, 'setValuesInStep', 'Disp_Control_Top'))
    for s in sorted(steps):
        print("         Disp_Control_Top %-18s u1=%.9g u2=%.9g" % (s, steps[s]['u1'], steps[s]['u2']))
    gt = steps.get('Step-GlobalTension', {})
    bk = steps.get('Step-Buckle', {})
    check('GlobalTension: u1=DISP_GLOBAL*top_cos, u2=-DISP_GLOBAL*top_sin',
          gt and abs(gt.get('u1', 9e9) - DG * tc) < 1e-18
          and abs(gt.get('u2', 9e9) - (-DG * ts)) < 1e-18, gt)
    check('GlobalTension: 실제 값 (0, +165um)', gt and abs(gt.get('u1', 9)) < 1e-18
          and abs(gt.get('u2', 0) - DG) < 1e-18, gt)
    check('Step-Buckle: PERTURBATION 도 세 꼭짓점 동일',
          bk and abs(bk.get('u2', 0) - PB) < 1e-18, bk)
    check('새 로그 1줄: [BC] PULL_TOP=True ...', sum(1 for s in lgT if s.startswith('[BC] PULL_TOP=True')) == 1, lgT)
    check('방향 로그 1줄 (Top/Right/Left 동일 형식)',
          sum(1 for s in lgT if 'cos=' in s and 'sin=' in s) == 1, lgT)
    check('로그 cp949 인코딩 가능(§19.95)', all(_enc_ok(s) for s in lgT), lgT)
    for s in lgT:
        print("         LOG| %s" % s)

    # Right/Left 와 '같은 형식'인지: 같은 스텝 집합을 받는가
    rsteps = set(s for s, _ in find(mT.rec, 'setValuesInStep', 'Disp_Control_Right'))
    tsteps = set(steps)
    check('Right 와 Top 의 구동 스텝 집합 동일', rsteps == tsteps, (sorted(rsteps), sorted(tsteps)))

    # ---------------------------------------------------------------- 최종하중 구역
    print("\n[(ii-b)] 최종 하중 구역 (HighTension / Postbuckle)")
    for tag, rnew, rold, step in (('HighTension', r_fin, r_fin0, 'Step-HighTension'),
                                  ('Postbuckle', r_post, r_post0, 'Step-Postbuckle')):
        mF2, _ = run(rnew, src, False, False)
        mB2, _ = run(rold, base, False, False)
        check('%s: False 추적 기저와 동일' % tag, mF2.trace() == mB2.trace() and mF2.rec,
              '\n--- new ---\n' + mF2.trace() + '\n--- base ---\n' + mB2.trace())
        check('%s: False 에서 Disp_Control_Top 없음' % tag,
              not find(mF2.rec, 'setValuesInStep', 'Disp_Control_Top'))
        mT2, _ = run(rnew, src, True, False)
        d = dict(find(mT2.rec, 'setValuesInStep', 'Disp_Control_Top'))
        print("         %-14s Top u1=%.9g u2=%.9g" % (tag, d[step]['u1'], d[step]['u2']))
        check('%s: Top = (0, +GLOBAL_FINAL)' % tag,
              abs(d[step]['u1']) < 1e-18 and abs(d[step]['u2'] - GF) < 1e-18, d)
        # Right/Left 와 같은 스텝/같은 크기 (크기는 벡터 크기 sqrt(u1^2+u2^2) 기준)
        rgt = dict(find(mT2.rec, 'setValuesInStep', 'Disp_Control_Right'))
        lft = dict(find(mT2.rec, 'setValuesInStep', 'Disp_Control_Left'))

        def _mag(d):
            return math.hypot(d['u1'], d['u2'])

        check('%s: 세 꼭짓점 변위 크기 동일 (|u|=GLOBAL_FINAL)' % tag,
              abs(_mag(d[step]) - GF) < 1e-18
              and abs(_mag(rgt[step]) - GF) < 1e-18
              and abs(_mag(lft[step]) - GF) < 1e-18,
              (_mag(d[step]), _mag(rgt[step]), _mag(lft[step])))

    # ---------------------------------------------------------------- (iii) 네 조합
    print("\n[(iii)] 네 조합 (PULL_TOP x NO_CLAMP)")
    for pt in (False, True):
        for nc in (False, True):
            try:
                mm, ll = run(r_bc, src, pt, nc)
                ok = True
                err = ''
            except Exception as e:
                ok, err = False, repr(e)
            print("         PULL_TOP=%-5s NO_CLAMP=%-5s -> %s" % (pt, nc, 'OK' if ok else 'FAIL ' + err))
            check('조합 PULL_TOP=%s NO_CLAMP=%s 예외 없음' % (pt, nc), ok, err)
            if ok:
                _want = 5 + (4 if not nc else 0) + (1 if pt else 0)
                check('  조합 BC 개수: %d (Top 포함)' % _want,
                      len(mm.boundaryConditions) == _want,
                      sorted(mm.boundaryConditions))

    # ---------------------------------------------------------------- (v) 금지문자
    print("\n[(v)] 새 문자열 문자 검사")
    added = [l for l in src.split('\n') if 'PULL_TOP' in l or 'Disp_Control_Top' in l
             or 'top_cos' in l or 'top_sin' in l]
    bad = [l for l in added if FORBIDDEN_RE.search(l)]
    check('새 줄(%d개) 금지문자 0' % len(added), not bad, bad)

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
