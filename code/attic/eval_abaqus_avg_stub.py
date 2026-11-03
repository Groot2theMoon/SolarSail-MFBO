# -*- coding: utf-8 -*-
"""eval_abaqus_avg_stub.py - Abaqus 없이 eval_abaqus_avg.main() 을 **실제 실행**한다 (§19.99)

무엇을 검증하나 (배선/로직 검증. 물리 검증이 아니다)
    가짜 odbAccess + 가짜 ODB(프레임/인스턴스/COORD bulkDataBlocks)를 끼우고
    eval_abaqus_avg.py 를 실제로 돌린다. 검증 항목:
      (1) 수직입사 자기검증 assert 통과 (hf_normal_self == eval_abaqus.calc_thrust_loss)
          + assert 가 **비어있지 않다**는 음성통제(1e-8 드리프트는 잡고 1e-10 은 통과)
      (2) Murphey 2002 앵커: 플랫 돛 phi=5도 total = 1.1083% +- 1e-4 (해석식 + 자체 코드 경로)
      (3) MC 수렴 + SE ~ 1/sqrt(N) 스케일
      (4) disk(전력가중) != hemi(입체각 균일), 그리고 완전평판 돛에서 wrinkle 정의는 정확히 0
          (total 정의는 해석값으로 수렴: hemi 0.7404291 / disk 0.5897910)
      (5) 메모리 블록 분할이 결과를 바꾸지 않는다(블록=1 강제 == 기본)
      (6) 스텝/프레임 선택(기본=Step-Postbuckle, 없으면 마지막 스텝), --frame 1 = 미변형 프레임
      (7) --flip/--no-flip 이 실제로 뒤집힌 요소에 작용한다(뒤집힌 면적비 > 0)
      (8) JSON 저장(기본 경로/지정 경로), odb.close() 호출(핸들 누수 없음)
      (9) 실패 경로: 없는 ODB / 없는 스텝 / 범위 밖 프레임 / 사용법 오류 -> rc != 0
      (10) 금지문자(em/en-dash, 화살표, smart quote, NBSP, bullet, ellipsis) 0 + 출력 cp949 인코딩

무엇을 못 보나
    실제 ODB 파서, 실제 메쉬, 실제 좌굴 형상, Abaqus numpy 버전. 값의 물리적 타당성은 실측만 판정한다.

사용 (numpy 가 있는 python 으로. Abaqus/라이선스 불필요)
    python3 eval_abaqus_avg_stub.py
종료코드: 0 OK / 1 FAIL / 2 numpy 없음
"""
import contextlib
import io
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import time
import types

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

if sys.version_info[0] < 3:                       # 저장소 스크립트는 py2/3 겸용이지만 하네스는 py3 전용
    print('  [stub] py3 필요'); sys.exit(2)
try:
    import numpy as np
except ImportError:
    print('  [stub] numpy 가 필요하다(Abaqus 대신 numpy 있는 python 으로 실행).')
    sys.exit(2)

CODE = HERE
SCRIPT = os.path.join(CODE, 'eval_abaqus_avg.py')

# 금지 문자(§19.95 교훈: cp949 콘솔에서 UnicodeEncodeError 로 도구 자체가 죽는다)
FORBIDDEN_RE = re.compile(u'[\u2010\u2011\u2012\u2013\u2014\u2015'
                          u'\u2190\u2191\u2192\u2193\u21d2'
                          u'\u2018\u2019\u201c\u201d'
                          u'\u00a0\u2022\u2026]')

FAIL = []
N = [0]


def check(name, cond, detail=None):
    N[0] += 1
    ok = bool(cond)
    print('  %-4s %s%s' % ('OK' if ok else 'FAIL', name,
                           ('   <- ' + str(detail)) if (detail is not None and not ok) else ''))
    if not ok:
        FAIL.append(name)
    return ok


# ------------------------------------------------------------------ 가짜 ODB
_CUR = {'odb': None}


class Blk(object):
    def __init__(self, nodeLabels, data):
        self.nodeLabels = nodeLabels
        self.data = data


class FO(object):
    def __init__(self, blocks):
        self.bulkDataBlocks = blocks

    def getSubset(self, region=None, position=None):
        return self


class Frame(object):
    def __init__(self, t, coords):
        self.frameValue = t
        self.fieldOutputs = {'COORD': FO([Blk(np.array(sorted(coords.keys())),
                                              np.array([coords[k] for k in sorted(coords.keys())]))])}


class Step(object):
    def __init__(self, frames):
        self.frames = frames


class Node(object):
    def __init__(self, label, xyz):
        self.label = label
        self.coordinates = xyz


class Elem(object):
    def __init__(self, label, conn):
        self.label = label
        self.connectivity = conn


class Inst(object):
    def __init__(self, name, nodes, elems):
        self.name = name
        self.nodes = nodes
        self.elements = elems


class Odb(object):
    def __init__(self, steps, instances):
        self.steps = steps
        self.rootAssembly = types.SimpleNamespace(instances=instances)
        self.closed = False

    def close(self):
        self.closed = True


def grid_coords(nx, ny, cell, amp, phase=1.0):
    """(i,j) 노드 좌표. amp>0 이면 z = amp*sin(2pi i/nx)*sin(2pi j/ny) 주름(높이장)."""
    co = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            z = 0.0
            if amp:
                z = (amp * phase * math.sin(2.0 * math.pi * i / float(nx))
                     * math.sin(2.0 * math.pi * j / float(ny)))
            co[lab(i, j, ny)] = (i * cell, j * cell, z)
    return co


def lab(i, j, ny):
    return i * (ny + 1) + j + 1


def elems_of(nx, ny, folded=(), hybrid=True):
    """hybrid=True -> 왼쪽 절반 quad, 오른쪽 절반 tri (실제 메쉬가 혼합이므로).
    folded = 뒤집을 요소 순번(연결성 순서를 뒤집어 n_z<0 을 만든다)."""
    out = []
    for i in range(nx):
        for j in range(ny):
            if hybrid and i < nx // 2:
                out.append((lab(i, j, ny), lab(i + 1, j, ny),
                            lab(i + 1, j + 1, ny), lab(i, j + 1, ny)))
            else:
                out.append((lab(i, j, ny), lab(i + 1, j, ny), lab(i + 1, j + 1, ny)))
                out.append((lab(i, j, ny), lab(i + 1, j + 1, ny), lab(i, j + 1, ny)))
    fixed = []
    for k, c in enumerate(out):
        fixed.append(tuple(reversed(c)) if k in folded else c)
    return fixed


def build_odb(nx=4, ny=3, cell=0.5, amp=0.30, folded=(), steps=('GT', 'PB'), flat=False):
    """가짜 ODB. GT = Step-GlobalTension, PB = Step-Postbuckle.
    frameValue=0 프레임(미변형) + 변형 프레임을 둔다. flat=True 면 전부 평평."""
    conns = elems_of(nx, ny, folded=folded)
    nodes = [Node(k, v) for k, v in sorted(grid_coords(nx, ny, cell, 0.0).items())]
    elems = [Elem(k + 1, conns[k]) for k in range(len(conns))]
    inst = Inst('MEMBRANE-1', nodes, elems)
    co0 = grid_coords(nx, ny, cell, 0.0)
    coF = grid_coords(nx, ny, cell, 0.0 if flat else amp)
    coM = grid_coords(nx, ny, cell, 0.0 if flat else amp, phase=0.6)
    st = {}
    if 'GT' in steps:
        st['Step-GlobalTension'] = Step([Frame(0.0, co0), Frame(1.0, coM)])
    if 'PB' in steps:
        st['Step-Postbuckle'] = Step([Frame(0.0, co0), Frame(0.5, coM), Frame(1.0, coF)])
    if 'CT' in steps:
        st['Step-ClampTension'] = Step([Frame(0.0, co0), Frame(1.0, coF)])
    return Odb(st, {'MEMBRANE-1': inst})


def install_fake_abaqus():
    odb_mod = types.ModuleType('odbAccess')
    odb_mod.openOdb = lambda path, readOnly=True: _CUR['odb']
    abc = types.ModuleType('abaqusConstants')
    abc.CENTROID = 'CENTROID'
    abc.MAX_INPLANE_PRINCIPAL = 'MAX_INPLANE_PRINCIPAL'
    abc.MIN_INPLANE_PRINCIPAL = 'MIN_INPLANE_PRINCIPAL'
    sys.modules['odbAccess'] = odb_mod
    sys.modules['abaqusConstants'] = abc


def run_main(argv, odb, cwd=None):
    """eval_abaqus_avg.main() 을 실제로 돌리고 (rc, stdout) 을 돌려준다."""
    import eval_abaqus_avg as M
    _CUR['odb'] = odb
    cwd0 = os.getcwd()
    if cwd:
        os.chdir(cwd)
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            rc = M.main(list(argv))
    finally:
        if cwd:
            os.chdir(cwd0)
    return rc, buf.getvalue()


def kv_of(out):
    for ln in out.splitlines():
        if ln.startswith('[eval_avg:KV]'):
            toks = ln.split('KV]', 1)[1].split()
            d = {}
            for t in toks:
                if '=' in t:
                    k, v = t.split('=', 1)
                    d[k] = v
            return d
    return {}


def gof(kv, key):
    try:
        return float(kv[key])
    except Exception:
        return float('nan')


def main():
    tmpd = tempfile.mkdtemp(prefix='avgstub_')
    fake = os.path.join(tmpd, 'HF_Postbuckle.odb')
    with io.open(fake, 'wb') as fh:
        fh.write(b'x' * 4096)
    install_fake_abaqus()
    import eval_abaqus_avg as M                                          # noqa: E402
    import eval_abaqus as E                                              # noqa: E402

    print('--- 0) 컴파일 / 상수 ---')
    p = subprocess.run([sys.executable, '-m', 'py_compile', SCRIPT],
                       capture_output=True, text=True, cwd=CODE)
    check('py_compile', p.returncode == 0, p.stderr[-300:])
    check('상수 import (A0/R0/IDEAL/U_SUN) 가 eval_abaqus 와 동일',
          (M.A0 == E.A0 and M.R0 == E.R0 and M.IDEAL_THRUST_COEFF == E.IDEAL_THRUST_COEFF
           and abs(M.IDEAL_THRUST_COEFF - (2 * E.R0 + E.A0)) < 1e-15),
          (M.A0, M.R0, M.IDEAL_THRUST_COEFF))

    # ---- 해석 기준값(완전평판): total 정의의 방향평균
    A0, R0, I = M.A0, M.R0, M.IDEAL_THRUST_COEFF
    hemi_flat = 1.0 - (A0 * 0.5 + 2.0 * R0 * (1.0 / 4.0)) / I            # c ~ U[0,1)
    disk_flat = 1.0 - (A0 * (2.0 / 3.0) + 2.0 * R0 * (2.0 / 5.0)) / I    # (x,y) 단위원 균일
    print('    해석(완전평판 total): hemi=%.7f  disk=%.7f' % (hemi_flat, disk_flat))

    # ================= (1)+(2) 주름진 가짜 돛으로 본 실행 =================
    print('--- 1) 주름 돛 실행: 자기검증 / 앵커 / MC / 스캔 / JSON ---')
    odb = build_odb(amp=0.35, folded=(3, 7))
    outdir = os.path.join(tmpd, 'run1')
    os.makedirs(outdir)
    rc, out = run_main([fake, '--n-mc', '3000', '--seed', '12345'], odb, cwd=outdir)
    check('rc=0', rc == 0, rc)
    kv = kv_of(out)
    check('KV 줄 존재 + hf_ref>0', bool(kv) and gof(kv, 'hf_ref') > 0.0, kv)
    check('hf_normal_self == hf_ref (KV)',
          abs(gof(kv, 'hf_normal_self') - gof(kv, 'hf_ref')) <= 1e-12 * abs(gof(kv, 'hf_ref')),
          (kv.get('hf_ref'), kv.get('hf_normal_self')))
    lines = out.splitlines()
    sc = [l for l in lines if '자기검증(ideal=' in l]
    check('자기검증 줄 2개(total/wrinkle) + rel<1e-9', len(sc) == 2
          and all('rel=' in l and l.rstrip().endswith('OK') for l in sc), sc)
    check('뒤집기 편향 줄에 뒤방향 면적비 > 0 (folded 요소 존재)',
          any('뒤방향(n_z<0) 면적비' in l for l in lines)
          and float([l for l in lines if '뒤방향(n_z<0) 면적비' in l][0]
                    .split('면적비')[1].split('%')[0]) > 0.0, lines[:3])
    check('MC 4조합 출력', all(('MC %s' % k) in out
                              for k in ('hemi_total', 'disk_total',
                                        'hemi_wrinkle', 'disk_wrinkle')))
    ht, dt = gof(kv, 'hemi_total'), gof(kv, 'disk_total')
    check('hemi_total != disk_total (샘플러 가중이 실제로 다르다)', abs(ht - dt) > 1e-3,
          (ht, dt))
    check('SE>0 이고 SE << mean (통계 유효)',
          gof(kv, 'hemi_total_se') > 0 and gof(kv, 'hemi_total_se') < 0.5 * ht,
          (kv.get('hemi_total_se'), ht))
    check('경사 스캔표 5행 + phi=5 줄 존재',
          sum(1 for l in lines if re.match(r'^\[eval_avg\]\s+[0-9]+\.[0-9]\s', l)) == 5
          and any(re.match(r'^\[eval_avg\]\s+5\.0\s', l) for l in lines))
    jpath = os.path.join(outdir, 'results', 'avg_HF_Postbuckle_Step-Postbuckle.json')
    check('JSON 기본 경로 저장', os.path.exists(jpath), jpath)
    if os.path.exists(jpath):
        rec = json.load(io.open(jpath, encoding='utf-8'))
        check('JSON 필수 키', all(k in rec for k in ('odb', 'step', 'hf_ref',
                                                     'hf_normal_self', 'mc', 'tilt_scan',
                                                     'n_mc', 'seed', 'murphey_ok')))
        check('JSON hf_normal_self == hf_ref (1e-15)',
              abs(rec['hf_normal_self']['total'] - rec['hf_ref'])
              <= 1e-15 * max(abs(rec['hf_ref']), 1e-30))
        check('JSON odb 신원(절대경로/mtime/나이) 기록',
              rec['odb']['path'] == os.path.abspath(fake) and rec['odb']['mtime']
              and rec['odb']['age_human'])
        check('JSON tilt_scan 5행 + flat_total phi=5 = 1.1083%',
              len(rec['tilt_scan']) == 5
              and abs(rec['tilt_scan'][2]['flat_total'] - 0.0110828196) < 1e-9)
        check('JSON mc 4조합', len(rec['mc']) == 4)
    check('odb.close() 호출(핸들 누수 없음)', odb.closed is True)

    # ================= (2) Murphey 앵커 =================
    print('--- 2) Murphey 2002 5도 앵커 (항등식) ---')
    rc, out = run_main([fake, '--check'], build_odb(amp=0.35, folded=(3, 7)), cwd=outdir)
    check('--check rc=0', rc == 0, (rc, out.splitlines()[-1:]))
    check('--check 가 앵커 1.1083% 를 찍는다', '1.1083' in out and 'Murphey' in out)
    check('--check 3개 항목 모두 OK', out.count('OK') >= 3, out.count('OK'))
    old_anchor = M.MURPHEY_5DEG
    M.MURPHEY_5DEG = 0.5
    rc_bad, out_bad = run_main([fake, '--check'], build_odb(amp=0.35, folded=(3, 7)),
                               cwd=outdir)
    M.MURPHEY_5DEG = old_anchor
    check('앵커 상수를 틀리게 두면 --check rc=1 (네거티브 통제)',
          rc_bad == 1 and 'FAIL' in out_bad, rc_bad)

    # ================= (1) 자기검증 assert 가 비어있지 않다 =================
    print('--- 3) 자기검증 음성통제 (assert 실효성) ---')
    g = None
    buf = io.StringIO()
    _CUR['odb'] = build_odb(amp=0.35, folded=(3, 7))
    with contextlib.redirect_stdout(buf):
        g = M.read_geometry(fake, None, None, True)
    check('read_geometry 성공', g is not None)
    nrm, ar, aref, asum, href = g['normals'], g['areas'], g['area_ref'], g['area_sum'], g['hf_ref']
    v, rel, smode = M.vertical_self_check(nrm, ar, aref, asum, href, 'total')
    check('자기검증 rel < 1e-9 (실측)', rel < 1e-9, rel)
    print('    hf_ref=%.10e  hf_normal_self=%.10e  %s=%.3e' % (href, v, smode, rel))
    bad = False
    try:
        M.vertical_self_check(nrm, ar, aref, asum, href * (1.0 + 1e-8), 'total')
    except AssertionError:
        bad = True
    check('1e-8 드리프트는 잡는다 (assert 작동)', bad)
    ok10 = True
    try:
        M.vertical_self_check(nrm, ar, aref, asum, href * (1.0 + 1e-10), 'total')
    except AssertionError:
        ok10 = False
    check('1e-10 드리프트는 통과(과도하게 조이지 않음)', ok10)
    n2 = nrm.copy()
    n2[0] = n2[0] + 1e-3
    caught = False
    try:
        M.vertical_self_check(n2, ar, aref, asum, href, 'total')
    except AssertionError:
        caught = True
    check('법선 1개 0.001 틀어짐도 잡는다', caught)

    # ================= (3)(4) 완전평판 돛: 해석값 수렴 / wrinkle=0 / 샘플러 차이 =================
    print('--- 4) 완전평판 돛: 방향평균 수렴 / wrinkle=0 / 샘플러 차이 ---')
    flat_odb = build_odb(amp=0.0, flat=True)
    rcd, outd = run_main([fake, '--n-mc', '20000', '--seed', '7'], flat_odb, cwd=outdir)
    kvd = kv_of(outd)
    check('평평한 돛 rc=0 + hf_ref~0', rcd == 0 and abs(gof(kvd, 'hf_ref')) < 1e-12,
          (rcd, kvd.get('hf_ref')))
    check('평판 hemi_total = 해석 %.6f (tol 5e-3)' % hemi_flat,
          abs(gof(kvd, 'hemi_total') - hemi_flat) < 5e-3,
          (gof(kvd, 'hemi_total'), hemi_flat))
    check('평판 disk_total = 해석 %.6f (tol 5e-3)' % disk_flat,
          abs(gof(kvd, 'disk_total') - disk_flat) < 5e-3,
          (gof(kvd, 'disk_total'), disk_flat))
    check('평판 hemi != disk (균일 vs 전력가중)',
          abs(gof(kvd, 'hemi_total') - gof(kvd, 'disk_total')) > 0.1,
          (gof(kvd, 'hemi_total'), gof(kvd, 'disk_total')))
    check('평판 wrinkle 정의는 0 (둘 다, tol 1e-10)',
          abs(gof(kvd, 'hemi_wrinkle')) < 1e-10 and abs(gof(kvd, 'disk_wrinkle')) < 1e-10,
          (kvd.get('hemi_wrinkle'), kvd.get('disk_wrinkle')))
    print('    평판 MC(seed=7,N=20000): hemi=%.6f disk=%.6f | wrinkle hemi=%.2e disk=%.2e '
          '| hf_ref=%.2e'
          % (gof(kvd, 'hemi_total'), gof(kvd, 'disk_total'), gof(kvd, 'hemi_wrinkle'),
             gof(kvd, 'disk_wrinkle'), gof(kvd, 'hf_ref')))

    # ================= (3) SE ~ 1/sqrt(N) =================
    print('--- 5) MC 수렴 / SE ~ 1/sqrt(N) 스케일 ---')
    buf = io.StringIO()
    _CUR['odb'] = flat_odb                     # 완전평판 -> 해석값과 직접 비교
    with contextlib.redirect_stdout(buf):
        gf = M.read_geometry(fake, None, None, True)
    est = {}
    for nn in (500, 2000, 8000):
        s = M.sample_dirs('disk', nn, seed=7)
        hf, _raw = M.hf_over_dirs(gf['normals'], gf['areas'], s, 'total',
                                  gf['area_ref'], gf['area_sum'])
        est[nn] = M.stats_of(hf)
        print('    N=%-5d mean=%.6f dev=%.2e SE=%.2e dev/SE=%.2f'
              % (nn, est[nn]['mean'], est[nn]['mean'] - disk_flat, est[nn]['se'],
                 (est[nn]['mean'] - disk_flat) / est[nn]['se']))
        check('N=%d |dev| < 3*SE' % nn,
              abs(est[nn]['mean'] - disk_flat) < 3.0 * est[nn]['se'] + 1e-6,
              (est[nn]['mean'], disk_flat, est[nn]['se']))
    check('SE 단조 감소 (N 증가)', est[8000]['se'] < est[2000]['se'] < est[500]['se'],
          [est[k]['se'] for k in (500, 2000, 8000)])
    ratio = est[500]['se'] / est[8000]['se']
    check('SE(500)/SE(8000) ~ sqrt(16)=4 (2.5..6.5)', 2.5 <= ratio <= 6.5, ratio)

    # ================= (5) 블록 분할 불변 =================
    print('--- 6) 메모리 블록 분할이 결과를 바꾸지 않는다 ---')
    s = M.sample_dirs('disk', 300, seed=3)
    h1, r1 = M.hf_over_dirs(gf['normals'], gf['areas'], s, 'total', gf['area_ref'],
                            gf['area_sum'])
    old_budget = M.MEM_BUDGET
    M.MEM_BUDGET = 4.0                        # 요소 40개 -> 블록 1개(방향 1개씩)
    h2, r2 = M.hf_over_dirs(gf['normals'], gf['areas'], s, 'total', gf['area_ref'],
                            gf['area_sum'])
    M.MEM_BUDGET = old_budget
    check('블록=1 강제 == 기본 (tol 1e-12)', np.allclose(h1, h2, rtol=1e-12, atol=1e-14)
          and np.allclose(r1, r2, rtol=1e-12, atol=1e-14),
          float(np.max(np.abs(h1 - h2))))
    check('블록 크기 상한 계산 (요소 10k -> 200, 요소 100 -> 500)',
          min(M.DEFAULT_BLOCK, int(M.MEM_BUDGET / 10000)) == 200
          and min(M.DEFAULT_BLOCK, int(M.MEM_BUDGET / 100)) == 500)

    # ================= (6) 스텝 / 프레임 선택 =================
    print('--- 7) 스텝 / 프레임 선택 ---')
    ct_odb = build_odb(amp=0.3, steps=('GT', 'CT'))
    rc, out = run_main([fake, '--n-mc', '200'], ct_odb, cwd=outdir)
    check('Step-Postbuckle 없으면 마지막 스텝(Step-ClampTension)',
          rc == 0 and kv_of(out).get('step') == 'Step-ClampTension'
          and '기본(마지막 step)' in out, kv_of(out).get('step'))
    rc, out = run_main([fake, '--step', 'Step-GlobalTension', '--n-mc', '200'], ct_odb,
                       cwd=outdir)
    check('--step 인자 존중 + [인자] 표기',
          rc == 0 and kv_of(out).get('step') == 'Step-GlobalTension' and '[인자]' in out)
    rc, out = run_main([fake, '--frame', '1', '--n-mc', '100'],
                       build_odb(amp=0.35, folded=(3, 7)), cwd=outdir)
    check('--frame 1 = 미변형 프레임(frameValue=0) -> hf_ref~0',
          rc == 0 and abs(gof(kv_of(out), 'hf_ref')) < 1e-12, kv_of(out).get('hf_ref'))
    rc, out = run_main([fake, '--frame', '99'], build_odb(amp=0.35), cwd=outdir)
    check('범위 밖 --frame -> rc=1 + 이유', rc == 1 and '범위 밖' in out, rc)
    rc, out = run_main([fake, '--step', 'Step-NoSuchStep'], build_odb(amp=0.35), cwd=outdir)
    check('없는 --step -> rc=1 + 있는 step 목록', rc == 1 and 'Step-Postbuckle' in out, rc)

    # ================= (7) flip =================
    print('--- 8) --flip / --no-flip ---')
    fold_odb = build_odb(amp=0.35, folded=(3, 7))
    rcA, outA = run_main([fake, '--n-mc', '200'], fold_odb, cwd=outdir)
    rcB, outB = run_main([fake, '--n-mc', '200', '--no-flip'], fold_odb, cwd=outdir)
    hfA, hfB = gof(kv_of(outA), 'hf_ref'), gof(kv_of(outB), 'hf_ref')
    check('flip ON/OFF rc=0', rcA == 0 and rcB == 0)
    check('뒤집힌 요소 때문에 hf(no-flip) > hf(flip)', hfB > hfA,
          (hfA, hfB))
    check('뒤집기 OFF 로그 표기', '뒤집기=OFF' in outB and '뒤집기=ON' in outA)
    rc, out = run_main([fake, '--no-flip', '--n-mc', '200', '--check'], fold_odb, cwd=outdir)
    check('--no-flip 에서도 자기검증 통과(rc=0)', rc == 0, rc)
    conns = elems_of(4, 3)
    nq = sum(1 for c in conns if len(c) == 4)
    nt = sum(1 for c in conns if len(c) == 3)
    check('topology/기하 카운트: quad+tri 혼합 인식 (%d/%d)' % (nq, nt),
          ('quad %d / tri %d' % (nq, nt)) in outA and '요소 %d개' % (nq + nt) in outA,
          [l for l in outA.splitlines() if '기하:' in l])

    # ================= (8) --ideal / --json =================
    print('--- 9) --ideal / --json ---')
    rc, out = run_main([fake, '--ideal', 'total', '--n-mc', '200'], fold_odb, cwd=outdir)
    check('--ideal total -> 2조합', rc == 0 and 'MC hemi_wrinkle' not in out
          and 'MC hemi_total' in out and 'MC disk_total' in out)
    rc, out = run_main([fake, '--ideal', 'wrinkle', '--n-mc', '200'], fold_odb, cwd=outdir)
    check('--ideal wrinkle -> 2조합', rc == 0 and 'MC hemi_total' not in out
          and 'MC hemi_wrinkle' in out)
    jp = os.path.join(tmpd, 'custom.json')
    rc, out = run_main([fake, '--json', jp, '--n-mc', '100', '--ideal', 'wrinkle'],
                       fold_odb, cwd=outdir)
    check('--json 지정 경로 저장 + mc 키 2개',
          rc == 0 and os.path.exists(jp)
          and len(json.load(io.open(jp, encoding='utf-8'))['mc']) == 2)

    # ================= (9) 실패 경로 / 사용법 =================
    print('--- 10) 실패 경로 / 사용법 ---')
    rc, out = run_main([os.path.join(tmpd, 'NOPE.odb')], fold_odb, cwd=outdir)
    check('없는 ODB -> rc=1 + 명시 메시지',
          rc == 1 and 'ODB 파일이 없다' in out and '!!!' in out, rc)
    for argv, why in (([], '인자 없음'), (['--nope'], '모르는 옵션'),
                      ([fake, '--ideal', 'bogus'], '잘못된 ideal'),
                      ([fake, '--n-mc', '0'], 'n-mc=0'),
                      ([fake, '--frame', '0'], 'frame=0'),
                      ([fake, '--n-mc', 'x'], 'n-mc 비정수')):
        rc, out = run_main(argv, fold_odb, cwd=outdir)
        check('사용법 오류(%s) -> rc=2 + usage' % why,
              rc == 2 and 'usage:' in out, (rc, out.splitlines()[:1]))
    rc, out = run_main(['--help'], fold_odb, cwd=outdir)
    check('--help rc=0 + usage', rc == 0 and 'usage:' in out, rc)

    # ================= (10) 금지문자 / 인코딩 =================
    print('--- 11) 금지문자 / cp949 인코딩 ---')
    for path in (SCRIPT, os.path.join(CODE, 'eval_abaqus_avg_stub.py')):
        src = io.open(path, encoding='utf-8').read()
        bad = FORBIDDEN_RE.findall(src)
        check('%s 금지문자 0' % os.path.basename(path), len(bad) == 0, set(bad))
    all_out = outA + outB + outd
    check('실행 출력에 금지문자 0', len(FORBIDDEN_RE.findall(all_out)) == 0,
          set(FORBIDDEN_RE.findall(all_out)))
    try:
        all_out.encode('cp949')
        enc = 'OK'
    except Exception as e:
        enc = 'FAIL(%s)' % e
    check('출력 전체 cp949 인코딩 가능(Windows 콘솔)', enc == 'OK', enc)
    check('새 스크립트가 eval_abaqus.py 를 수정하지 않았다',
          subprocess.run(['git', 'diff', '--quiet', '--', 'eval_abaqus.py'],
                         cwd=CODE).returncode == 0)

    print('[stub-test] eval_abaqus_avg (19.99) %s  (%d/%d)'
          % ('OK' if not FAIL else 'FAIL', N[0] - len(FAIL), N[0]))
    if FAIL:
        print('  실패 항목: %s' % FAIL)
    return 0 if not FAIL else 1


if __name__ == '__main__':
    sys.exit(main())
