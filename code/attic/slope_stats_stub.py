# -*- coding: utf-8 -*-
"""slope_stats_stub.py - Abaqus 없이 slope_stats.main() 을 **실제 실행**한다 (§19.101)

무엇을 검증하나 (배선/로직/해석앵커 검증. 실제 ODB 물리 검증이 아니다)
    가짜 odbAccess + 가짜 ODB(프레임/인스턴스/COORD bulkDataBlocks/S)를 끼우고
    slope_stats.py 를 실제로 돌린다. 검증 항목:

      (i)   평평한 프레임: 모든 theta=0, hf_actual~0, hf_pred~0 -> 비율 nan(0 으로 나눠 inf 없음)
      (ii)  해석 앵커: 모든 요소 법선을 z_hat 에서 **0.298도** 균일 기울인 fixture
              -> hf_pred_abs ~ 3.94e-5 (실측 (0.5,0.3) 값) 이고 k*theta^2 와 정확식이 일치
              -> theta_fluct = 0 (법선 전부 동일 -> 국소평균 = 자기 자신)
      (iii) 단일 파장 사인파 fixture(주름) -> theta_fluct rms < theta_abs rms (거시/주름 분리)
      (iv)  패치별 rms 10x10 / 최악 패치 위치 / 방위각 18빈 히스토그램 출력 OK
            + 사인파 fixture 에서 배향성비 ~ 1 (정렬), 가우스 bump fixture 에서 최악 패치가
              bump 근처에 온다
      (v)   금지문자(em/en-dash, 화살표, smart quote, NBSP, bullet, ellipsis) 0 + 출력 cp949 인코딩
      (vi)  기존 4개 스크립트(eval_abaqus.py/eval_abaqus_avg.py/base_state_probe.py/
            run_abaqus.py)가 수정되지 않았다 (git diff --quiet)

    추가: 스텝/프레임 선택, 실패 경로(rc=1/2), JSON(기본/지정 경로), odb.close(),
          k_slope == k_theory, 국소평균 이웃수, 순수함수(area_stats/patch_grid/azimuth_hist/
          compression_summary) 단위 검증.

사용 (numpy 가 있는 python 으로. Abaqus/라이선스 불필요)
    python3 slope_stats_stub.py
종료코드: 0 OK / 1 FAIL / 2 numpy 없음
"""
import contextlib
import io
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import types

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

if sys.version_info[0] < 3:
    print('  [stub] py3 필요')
    sys.exit(2)
try:
    import numpy as np
except ImportError:
    print('  [stub] numpy 가 필요하다(Abaqus 대신 numpy 있는 python 으로 실행).')
    sys.exit(2)

CODE = HERE
SCRIPT = os.path.join(CODE, 'slope_stats.py')

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


class COORDfo(object):
    def __init__(self, blocks):
        self.bulkDataBlocks = blocks

    def getSubset(self, region=None, position=None):
        return self


class SecPt(object):
    def __init__(self, n):
        self.number = n


class SVAL(object):
    def __init__(self, el, data, inst='MEMBRANE-1', secpt=None):
        self.elementLabel = el
        self.data = data
        self.instance = types.SimpleNamespace(name=inst) if inst else None
        self.sectionPoint = SecPt(secpt)


class Sfo(object):
    def __init__(self, vals):
        self.values = vals

    def getSubset(self, position=None, region=None):
        return self


class Frame(object):
    def __init__(self, t, coords, svals=None):
        self.frameValue = t
        labs = np.array(sorted(coords.keys()), dtype=np.int64)
        self.fieldOutputs = {
            'COORD': COORDfo([Blk(labs, np.array([coords[k] for k in sorted(coords.keys())],
                                                 dtype=np.float64))]),
            'S': Sfo(svals or [])}


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


# ------------------------------------------------------------------ 메쉬 생성기
def lab(i, j, ny):
    return i * (ny + 1) + j + 1


def grid_coords(nx, ny, cell, shape=None):
    co = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            x = i * cell
            y = j * cell
            z = 0.0 if shape is None else shape(x, y)
            co[lab(i, j, ny)] = (x, y, z)
    return co


def elems_of(nx, ny, folded=(), hybrid=True):
    """hybrid=True -> 왼쪽 절반 quad, 오른쪽 절반 tri (실제 메쉬가 혼합이므로)."""
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


def mk_svals(conns):
    """요소 k: s22 = -1000(짝수) / +1000(홀수), s11 = +500, s12 = 0."""
    vals = []
    for k, c in enumerate(conns):
        s22 = -1000.0 if (k % 2 == 0) else 1000.0
        vals.append(SVAL(k + 1, [500.0, s22, 0.0]))
    return vals


def build_odb(nx=8, ny=8, cell=0.5, shape=None, folded=(), steps=('GT', 'PB'),
              with_stress=True):
    """가짜 ODB. GT/CT/PB 스텝. 첫 프레임(frameValue=0) = 미변형(평평)."""
    conns = elems_of(nx, ny, folded=folded)
    nodes = [Node(k, v) for k, v in sorted(grid_coords(nx, ny, cell).items())]
    elems = [Elem(k + 1, conns[k]) for k in range(len(conns))]
    inst = Inst('MEMBRANE-1', nodes, elems)
    co0 = grid_coords(nx, ny, cell)
    coF = grid_coords(nx, ny, cell, shape)
    sv = mk_svals(conns) if with_stress else None
    st = {}
    if 'GT' in steps:
        st['Step-GlobalTension'] = Step([Frame(0.0, co0, sv), Frame(1.0, coF, sv)])
    if 'PB' in steps:
        st['Step-Postbuckle'] = Step([Frame(0.0, co0, sv), Frame(0.5, coF, sv),
                                      Frame(1.0, coF, sv)])
    if 'CT' in steps:
        st['Step-ClampTension'] = Step([Frame(0.0, co0, sv), Frame(1.0, coF, sv)])
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
    import slope_stats as M
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
        if ln.startswith('[slope:KV]'):
            d = {}
            for t in ln.split('KV]', 1)[1].split():
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


# ------------------------------------------------------------------ fixtures
TILT_DEG = 0.298


def shape_tilt(x, y):
    """z = tan(0.298 deg) * x -> 모든 요소 법선이 z_hat 에서 정확히 0.298도."""
    return math.tan(math.radians(TILT_DEG)) * x


def shape_sine(x, y):
    """1D 사인파 주름: z = A*sin(2*pi*x/lam). 법선 수평성분이 +-x 방향 -> 배향성 ~ 1."""
    return 0.05 * math.sin(2.0 * math.pi * x / 4.0)


def shape_bump(x, y):
    """가우스 bump: 최악 패치가 bump 근처(0.8, 3.2)에 오는지 검증용."""
    return 0.1 * math.exp(-((x - 0.8) ** 2 + (y - 3.2) ** 2) / (2.0 * 0.35 ** 2))


# ------------------------------------------------------------------ 순수함수 단위검증
def pure_tests():
    import slope_stats as M
    print('--- 순수함수 단위검증 ---')
    # area_stats: 두 값 0.0 과 2.0, 면적 1:3 -> mean=1.5, rms=sqrt(3)=1.732
    th = np.array([0.0, 2.0])
    w = np.array([1.0, 3.0])
    st = M.area_stats(th, w)
    check('area_stats 면적가중 mean (0,2 | 1:3 -> 1.5)', abs(st['mean_rad'] - 1.5) < 1e-12,
          st['mean_rad'])
    check('area_stats 면적가중 rms (-> sqrt(3))',
          abs(st['rms_rad'] - math.sqrt(3.0)) < 1e-12, st['rms_rad'])
    check('area_stats deg 병기 (rms deg = deg(rms rad))',
          abs(st['rms_deg'] - math.degrees(math.sqrt(3.0))) < 1e-12, st['rms_deg'])
    check('area_stats p50 면적중앙값 = 2.0 (면적 3 이 0.0 쪽 1 보다 큼)',
          abs(st['p50_rad'] - 2.0) < 1e-12, st['p50_rad'])
    check('area_stats max=2.0 / 빈입력 nan', abs(st['max_rad'] - 2.0) < 1e-12
          and M.area_stats(np.zeros(0), np.zeros(0))['rms_rad'] != M.area_stats(np.zeros(0), np.zeros(0))['rms_rad'])

    # patch_grid: (0.1,0.1) 에 큰 값 -> 최악 패치가 (0,0)
    cx = np.array([0.1, 5.0, 9.9])
    cy = np.array([0.1, 5.0, 9.9])
    t = np.array([1.0, 0.1, 0.2])
    a = np.array([1.0, 1.0, 1.0])
    grid, worst, bounds = M.patch_grid(cx, cy, t, a)
    check('patch_grid 최악 = (0,0) rms=1.0', worst is not None and worst[0] == 0
          and worst[1] == 0 and abs(worst[2] - 1.0) < 1e-12, worst)
    check('patch_grid 빈 칸 nan + 격자 10x10',
          len(grid) == 10 and all(len(r) == 10 for r in grid)
          and sum(1 for r in grid for v in r if v == v) == 3, None)

    # azimuth_hist: (1,0) -> phi=0 -> 0번빈, 2배각 비=1
    nrm = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    arc = np.array([1.0, 1.0, 1.0])
    frac, orient, dom = M.azimuth_hist(nrm, arc)
    check('azimuth_hist x축 정렬 -> 배향성비=1, 빈 0/9',
          abs(orient - 1.0) < 1e-12 and frac[0] > 0 and frac[9] > 0, (orient, frac[0], frac[9]))
    # 등방(4방향, 2배각이 정확히 상쇄되게 45도 간격) -> 2배각 비 ~ 0
    n2 = np.array([[math.cos(a), math.sin(a), 0.0] for a in
                   (0.0, math.pi / 4.0, math.pi / 2.0, 3.0 * math.pi / 4.0)])
    f2, o2, d2 = M.azimuth_hist(n2, np.ones(4))
    check('azimuth_hist 4방향 45도 균등 -> 배향성비 ~ 0', o2 < 1e-12, o2)

    # compression_summary: 손계산 대조
    areas = [1.0, 1.0, 2.0]
    labels = [1, 2, 3]
    stress = {1: (500.0, -1000.0, 0.0),      # minP=-1000<0, mean=-250<0
              2: (500.0, 1000.0, 0.0),       # minP= +500>0, mean= 750>0
              3: (0.0, -100.0, 0.0)}         # minP= -100<0, mean= -50<0
    c = M.compression_summary(areas, labels, stress)
    check('compression_summary minP<0 면적비 = (1+2)/4 = 0.75',
          abs(c['area_w_minp_neg'] - 0.75) < 1e-12, c['area_w_minp_neg'])
    check('compression_summary mean<0 면적비 = (1+2)/4 = 0.75',
          abs(c['area_w_mean_neg'] - 0.75) < 1e-12, c['area_w_mean_neg'])
    check('compression_summary minP 면적평균 = (1*-1000 + 1*500 + 2*-100)/4 = -175',
          abs(c['minp_area_mean_pa'] - (-175.0)) < 1e-12, c['minp_area_mean_pa'])

    #   _json_safe: 비표준 NaN/Infinity 토큰을 만들지 않는다
    js = M._json_safe({'a': float('nan'), 'b': float('inf'), 'c': [1.0, float('nan')],
                       'd': np.float64('2.5'), 'e': 3, 'f': 'x', 'g': True})
    check('_json_safe: NaN/Inf -> None, 나머지 보존',
          js['a'] is None and js['b'] is None and js['c'] == [1.0, None]
          and abs(js['d'] - 2.5) < 1e-12 and js['e'] == 3 and js['f'] == 'x' and js['g'] is True,
          js)

    # local_mean_normals: 평면(z=0)에서는 국소평균 == 자기 법선 -> 요동 0
    #   (face_geometry 는 coords 를 **행 인덱스**로 참조한다 - 실 ODB 에서는 라벨 인덱스)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0]])
    quads = [(0, 1, 4, 3), (1, 2, 5, 4)]
    a, n, ctr = M.face_geometry(coords, quads, [], True)
    mn, cnt, h, R = M.local_mean_normals(ctr, n, a)
    dev = M.angle_between(n, mn)
    check('local_mean_normals 평면 -> 요동 0', float(np.max(dev)) < 1e-12, float(np.max(dev)))
    check('local_mean_normals 이웃수 >= 1 (자기 포함)', bool(np.all(cnt >= 1)), cnt.tolist())


# ------------------------------------------------------------------ §19.102 검증 헬퍼
def _unsetenv(k):
    if k in os.environ:
        del os.environ[k]


def _enc_ok(s):
    try:
        s.encode('cp949')
        return True
    except Exception:
        return False


def _isnum(v):
    try:
        float(v)
        return True
    except Exception:
        return False


def _csv_cols_block(src):
    m = re.search(r"_cols = \[(.*?)\]", src, re.S)
    return m.group(0) if m else ''


def _new_run_abaqus_lines(src):
    """§19.102 로 새로 쓴 run_abaqus.py 줄만 골라낸다(마커 기반 - 커밋 전후 무관하게 결정적)."""
    marks = ('[SLOPE]', '[ODBCOPY]', '§19.102', '_archive_slope',
             'SOLARSAIL_NO_SLOPE', 'SOLARSAIL_NO_ODB_COPY')
    return [ln for ln in src.splitlines() if any(m in ln for m in marks)]


def extract_wiring(src_path):
    """run_abaqus.py 에서 _archive_slope 함수 소스만 떼어낸다(다음 def _r13_ct_keys 전까지)."""
    src = io.open(src_path, encoding='utf-8').read()
    m = re.search(r'\ndef _archive_slope\(.*?\n(?=def _r13_ct_keys\()', src, re.S)
    return (m.group(0) if m else None), src


# ------------------------------------------------------------------ §19.102 (i) compute_stats API
def stats_api_tests(fake, M, np):
    print('--- §19.102 (i) compute_stats dict + 목적함수 후보 f1~f4 ---')
    _CUR['odb'] = build_odb(nx=10, ny=8, cell=0.5, shape=shape_sine)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        S = M.compute_stats(fake)
    check('compute_stats -> dict', isinstance(S, dict), type(S))
    if not isinstance(S, dict):
        return
    need = ('theta_abs_rms_deg', 'theta_abs_p95_deg', 'theta_abs_p99_deg', 'theta_abs_max_deg',
            'worst_patch_rms_deg', 'worst_patch_xy', 'orient_ratio', 'minp_frac_neg',
            'area_ratio', 'hf_actual', 'ratio_abs', 'hf_pred_exact',
            'f1_hf_pred', 'f2_p95', 'f3_worst_patch', 'f4_cvar90')
    check('필수 키 전부 존재 (목적함수 + 기존 통계)',
          all(k in S for k in need), [k for k in need if k not in S])
    check('f1_hf_pred == hf_pred_abs (항등)',
          abs(S['f1_hf_pred'] - S['hf_pred_abs']) < 1e-18,
          (S['f1_hf_pred'], S['hf_pred_abs']))
    check('f2_p95 >= f1_hf_pred', S['f2_p95'] >= S['f1_hf_pred'] - 1e-18,
          (S['f2_p95'], S['f1_hf_pred']))
    check('f4_cvar90 >= f2_p95', S['f4_cvar90'] >= S['f2_p95'] - 1e-18,
          (S['f4_cvar90'], S['f2_p95']))
    check('f3_worst_patch >= f1_hf_pred', S['f3_worst_patch'] >= S['f1_hf_pred'] - 1e-18,
          (S['f3_worst_patch'], S['f1_hf_pred']))
    check('theta2_q95 >= theta2_mean (면적가중 분위수 정의)',
          S['theta2_q95'] >= S['theta2_mean'] - 1e-18, None)
    check('theta2_cvar90_mean >= theta2_q95', S['theta2_cvar90_mean'] >= S['theta2_q95'] - 1e-18,
          (S['theta2_cvar90_mean'], S['theta2_q95']))
    print('    sine fixture: f1=%.6e f2=%.6e f3=%.6e f4=%.6e'
          % (S['f1_hf_pred'], S['f2_p95'], S['f3_worst_patch'], S['f4_cvar90']))
    items = M.scalar_items(S, prefix='slope_')
    check('scalar_items: 전부 slope_ 접두 + key=value',
          bool(items) and all(t.startswith('slope_') and '=' in t for t in items), items[:2])
    check('scalar_items: f1~f4 포함',
          any(t.startswith('slope_f1_hf_pred=') for t in items)
          and any(t.startswith('slope_f4_cvar90=') for t in items), None)
    check('scalar_items: 중첩(dict/list) 키 제외 (patches/theta_abs/grid/note 없음)',
          not any(t.startswith(('slope_patches=', 'slope_theta_abs=', 'slope_grid_abs_rad=',
                                'slope_note=')) for t in items), None)
    _CUR['odb'] = build_odb(nx=8, ny=8, cell=0.5, shape=None)
    with contextlib.redirect_stdout(io.StringIO()):
        S0 = M.compute_stats(fake)
    check('평평 프레임: f1=f2=f3=f4 ~ 0',
          all(abs(S0[k]) < 1e-30 for k in ('f1_hf_pred', 'f2_p95', 'f3_worst_patch', 'f4_cvar90')),
          [S0[k] for k in ('f1_hf_pred', 'f2_p95', 'f3_worst_patch', 'f4_cvar90')])


# ------------------------------------------------------------------ §19.102 (ii) run_abaqus 배선
def wiring_tests(fake, tmpd, M):
    print('--- §19.102 (ii) run_abaqus 배선: ODB 태깅 + slope_ 아카이빙 + 환경변수 스위치 ---')
    rasrc_path = os.path.join(CODE, 'run_abaqus.py')
    fn_src, rasrc = extract_wiring(rasrc_path)
    check('run_abaqus.py 에서 _archive_slope 소스 추출', bool(fn_src))
    if not fn_src:
        return
    cdir = os.path.join(tmpd, 'code')
    rdir = os.path.join(cdir, 'aba')
    os.makedirs(rdir)
    with io.open(os.path.join(rdir, 'HF_Postbuckle.odb'), 'wb') as fh:
        fh.write(b'o' * 2048)
    os.makedirs(os.path.join(cdir, 'results'))
    ptxt = os.path.join(cdir, 'results', 'results_xc0p5_dc0p3.txt')
    base = 'x_c=0.5\nd_c=0.3\nfidelity=HF\nhf=3.9e-05\n'
    with io.open(ptxt, 'w', encoding='utf-8') as fh:
        fh.write(base)

    logs = []

    def _emit(*a):
        logs.append(' '.join(str(x) for x in a))

    from aba_imperfection import value_tag
    ns = {'os': os, 'sys': sys, 'io': io, 'shutil': shutil, 'math': math, 'time': time,
          'value_tag': value_tag, 'emit': _emit, 'slope_stats': M,
          '_RUN': rdir, '_HERE': cdir, 'HF_ODB': 'HF_Postbuckle.odb'}
    exec(compile(fn_src, '<run_abaqus:_archive_slope>', 'exec'), ns)
    fn = ns['_archive_slope']

    dst = os.path.join(cdir, 'results', 'odb', 'HF_xc0p5_dc0p3.odb')
    jp = os.path.join(cdir, 'results', 'slope_HF_xc0p5_dc0p3_Step-Postbuckle.json')

    # (a) 정상 경로 (스위치 미설정)
    _unsetenv('SOLARSAIL_NO_SLOPE')
    _unsetenv('SOLARSAIL_NO_ODB_COPY')
    _CUR['odb'] = build_odb(nx=8, ny=8, cell=0.5, shape=shape_tilt, folded=(3,))
    with contextlib.redirect_stdout(io.StringIO()):
        fn(0.5, 0.3)
    txt = io.open(ptxt, encoding='utf-8').read()
    check('정상: ODB 태깅 복사 results/odb/HF_xc0p5_dc0p3.odb', os.path.exists(dst))
    check('정상: 점별 txt 에 slope_ 접두 key=value append',
          'slope_f1_hf_pred=' in txt and 'slope_worst_patch_rms_deg=' in txt
          and 'slope_theta_abs_rms_deg=' in txt, None)
    check('정상: 기존 RESULTS 계약 필드 유지 (fidelity/hf)',
          'fidelity=HF' in txt and 'hf=3.9e-05' in txt, None)
    check('정상: JSON 저장 (slope_HF_xc0p5_dc0p3_Step-Postbuckle.json)', os.path.exists(jp), jp)
    check('정상: [ODBCOPY] 로그에 대상 경로/크기', any(l.startswith('[ODBCOPY]')
          and 'HF_xc0p5_dc0p3.odb' in l and '2048 B' in l for l in logs), logs[-3:])
    check('정상: [SLOPE] 로그', any(l.startswith('[SLOPE]') for l in logs), logs[-3:])
    #   (iii) txt 파서: slope_ 키 개수 + 수치 파싱
    parsed = {}
    for ln in txt.splitlines():
        if ln.startswith('slope_') and '=' in ln:
            k, v = ln.split('=', 1)
            parsed[k] = v
    check('(iii) txt 파서: slope_ 키 >= 20 개', len(parsed) >= 20, len(parsed))
    check('(iii) txt 파서: f1~f4 + 통계 키 수치 파싱',
          all(_isnum(parsed.get(k, 'x')) for k in
              ('slope_f1_hf_pred', 'slope_f2_p95', 'slope_f3_worst_patch', 'slope_f4_cvar90',
               'slope_theta_abs_rms_deg', 'slope_worst_patch_rms_deg', 'slope_orient_ratio',
               'slope_minp_frac_neg', 'slope_area_ratio', 'slope_hf_actual',
               'slope_hf_pred_exact')),
          {k: parsed.get(k) for k in ('slope_f1_hf_pred', 'slope_worst_patch_rms_deg')})
    if os.path.exists(jp):
        rec = json.load(io.open(jp, encoding='utf-8'))
        check('JSON objectives 블록에 f1~f4',
              all(k in rec.get('objectives', {})
                  for k in ('f1_hf_pred', 'f2_p95', 'f3_worst_patch', 'f4_cvar90')),
              sorted(rec.get('objectives', {}).keys()))

    # (b) SOLARSAIL_NO_SLOPE=1 -> 배선 전체 생략
    for f in (dst, jp):
        if os.path.exists(f):
            os.remove(f)
    with io.open(ptxt, 'w', encoding='utf-8') as fh:
        fh.write(base)
    logs[:] = []
    os.environ['SOLARSAIL_NO_SLOPE'] = '1'
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            fn(0.5, 0.3)
    finally:
        _unsetenv('SOLARSAIL_NO_SLOPE')
    txt = io.open(ptxt, encoding='utf-8').read()
    check('NO_SLOPE=1: ODB 복사 안 함', not os.path.exists(dst))
    check('NO_SLOPE=1: slope_ 줄 없음', 'slope_f1_hf_pred=' not in txt)
    check('NO_SLOPE=1: SKIP 로그 1줄', any('SOLARSAIL_NO_SLOPE=1' in l and 'SKIP' in l
                                          for l in logs), logs)

    # (c) SOLARSAIL_NO_ODB_COPY=1 -> 복사만 생략, slope_ 는 남음 (경로 분기)
    logs[:] = []
    with io.open(ptxt, 'w', encoding='utf-8') as fh:
        fh.write(base)
    os.environ['SOLARSAIL_NO_ODB_COPY'] = '1'
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            fn(0.5, 0.3)
    finally:
        _unsetenv('SOLARSAIL_NO_ODB_COPY')
    txt = io.open(ptxt, encoding='utf-8').read()
    check('NO_ODB_COPY=1: ODB 복사 안 함', not os.path.exists(dst))
    check('NO_ODB_COPY=1: slope_ 는 남음 + SKIP 로그',
          'slope_f1_hf_pred=' in txt and any('[ODBCOPY] SKIP' in l for l in logs), None)

    # (iv) ASCII/금지문자/cp949 + CSV 열 스키마 불변
    check('(iv) 추출 배선 소스 금지문자 0', len(FORBIDDEN_RE.findall(fn_src)) == 0,
          sorted(set(FORBIDDEN_RE.findall(fn_src))))
    _add = '\n'.join(_new_run_abaqus_lines(rasrc))
    check('(iv) run_abaqus.py **추가분** 금지문자 0 (기존 줄은 불변)',
          len(FORBIDDEN_RE.findall(_add)) == 0 and len(_add) > 0,
          sorted(set(FORBIDDEN_RE.findall(_add))))
    check('(iv) 점별 txt cp949 인코딩 가능 (Windows 콘솔/파일)', _enc_ok(txt), None)
    check('(iv) all_results.csv 열 스키마 불변 (slope 미포함, hf 마지막)',
          "'lf1', 'lf2', 'lf3', 'hf']" in _csv_cols_block(rasrc)
          and 'slope' not in _csv_cols_block(rasrc), None)


# ------------------------------------------------------------------ main
def main():
    tmpd = tempfile.mkdtemp(prefix='slopestub_')
    fake = os.path.join(tmpd, 'HF_Postbuckle.odb')
    with io.open(fake, 'wb') as fh:
        fh.write(b'x' * 4096)
    install_fake_abaqus()
    import slope_stats as M                                             # noqa: E402

    print('--- 0) 컴파일 / 상수 ---')
    p = subprocess.run([sys.executable, '-m', 'py_compile', SCRIPT],
                       capture_output=True, text=True, cwd=CODE)
    check('py_compile', p.returncode == 0, p.stderr[-300:])
    check('k_slope == k_theory (A0/2+3R0)/IDEAL',
          abs(M.K_SLOPE - M.K_SLOPE_THEORY) < 1e-5,
          (M.K_SLOPE, M.K_SLOPE_THEORY))
    print('    k(사용)=%.6f  k_theory=%.6f  U_SUN=%s'
          % (M.K_SLOPE, M.K_SLOPE_THEORY, np.asarray(M.U_SUN).tolist()))

    pure_tests()

    outdir = os.path.join(tmpd, 'run')
    os.makedirs(outdir)

    # ================= §19.102: compute_stats API + run_abaqus 배선 =================
    stats_api_tests(fake, M, np)
    wiring_tests(fake, tmpd, M)

    # ================= (i) 평평한 프레임 =================
    print('--- (i) 평평한 프레임: theta=0 / hf=0 / 비율 nan ---')
    flat = build_odb(nx=8, ny=8, cell=0.5, shape=None)
    rc, out = run_main([fake], flat, cwd=outdir)
    kv = kv_of(out)
    check('평평 프레임 rc=0', rc == 0, rc)
    check('hf_actual ~ 0 (< 1e-12)', abs(gof(kv, 'hf_actual')) < 1e-12, kv.get('hf_actual'))
    check('theta_abs rms = 0 (deg)', abs(gof(kv, 'theta_abs_rms_deg')) < 1e-9,
          kv.get('theta_abs_rms_deg'))
    check('theta_fluct rms = 0 (deg)', abs(gof(kv, 'theta_fluct_rms_deg')) < 1e-9,
          kv.get('theta_fluct_rms_deg'))
    check('hf_pred_abs = 0 이고 비율은 nan (inf/0-division 없음)',
          abs(gof(kv, 'hf_pred_abs')) < 1e-12 and math.isnan(gof(kv, 'ratio_abs')),
          (kv.get('hf_pred_abs'), kv.get('ratio_abs')))
    check('모든 비율 nan (abs/fluct/exact)',
          math.isnan(gof(kv, 'ratio_abs')) and math.isnan(gof(kv, 'ratio_fluct'))
          and math.isnan(gof(kv, 'ratio_exact')),
          (kv.get('ratio_abs'), kv.get('ratio_fluct'), kv.get('ratio_exact')))
    check('비율 미정의 로그 1줄', '비율은 정의되지 않는다' in out)
    check('patch/azimuth 압축장 줄이 평평 프레임에서도 나온다',
          any('최악 패치' in l for l in out.splitlines())
          and any('배향성비' in l for l in out.splitlines())
          and any('압축장 요약(R-13 정의)' in l for l in out.splitlines()), None)
    print('    %s' % [l.strip() for l in out.splitlines() if 'hf_actual=' in l][0])

    # ================= (ii) 해석 앵커: 0.298도 균일 기울기 =================
    print('--- (ii) 해석 앵커: z_hat 에서 0.298도 균일 기울기 ---')
    tilt = build_odb(nx=8, ny=8, cell=0.5, shape=shape_tilt, folded=(3,))
    rc, out = run_main([fake], tilt, cwd=outdir)
    kv = kv_of(out)
    check('균일기울기 rc=0', rc == 0, rc)
    check('theta_abs rms = 0.298도 (tol 1e-9)',
          abs(gof(kv, 'theta_abs_rms_deg') - TILT_DEG) < 1e-9, kv.get('theta_abs_rms_deg'))
    check('theta_fluct rms ~ 0 (법선 전부 동일 -> 국소평균 = 자기)',
          abs(gof(kv, 'theta_fluct_rms_deg')) < 1e-9, kv.get('theta_fluct_rms_deg'))
    hpa = gof(kv, 'hf_pred_abs')
    check('hf_pred_abs ~ 3.94e-05 (실측 (0.5,0.3), tol 1%)',
          abs(hpa - 3.9408139330e-05) < 0.01 * 3.9408139330e-05, hpa)
    check('hf_pred_abs == k*<theta^2> (k_theory 와 tol 1e-6 상대)',
          abs(hpa - M.K_SLOPE_THEORY * math.radians(TILT_DEG) ** 2)
          < 1e-6 * hpa, (hpa, M.K_SLOPE_THEORY * math.radians(TILT_DEG) ** 2))
    check('hf_pred_exact_norm (= 1-<f>) == hf_pred_abs (tol 1e-3 상대: 작은 각 근사 검증)',
          abs(gof(kv, 'hf_pred_exact_norm') - hpa) < 1e-3 * hpa,
          (gof(kv, 'hf_pred_exact_norm'), hpa))
    check('hf_pred_exact(요구식) == 1 - area_ratio*(1-hf_pred_exact_norm) (항등식)',
          abs(gof(kv, 'hf_pred_exact')
              - (1.0 - gof(kv, 'area_ratio') * (1.0 - gof(kv, 'hf_pred_exact_norm'))))
          < 1e-12, None)
    check('비율 hf_pred_exact_norm / hf_actual = 1 (tol 1e-4): 정확 SRP 식이 참조구현과 일치',
          abs(gof(kv, 'ratio_exact_norm') - 1.0) < 1e-4, kv.get('ratio_exact_norm'))
    check('면적비 Sum(A_e)/area_ref = 1/cos(0.298도) (tol 1e-9)',
          abs(gof(kv, 'area_ratio') - 1.0 / math.cos(math.radians(TILT_DEG))) < 1e-9,
          kv.get('area_ratio'))
    check('요구식 비율(=0.658)은 면적비 항 때문에 1 이 아니다 (로그에 명시)',
          '면적비 항 때문에 1 이 아님' in out and abs(gof(kv, 'ratio_exact') - 1.0) > 0.3,
          kv.get('ratio_exact'))
    check('비율 hf_pred_abs / hf_actual ~ 1 (tol 1%)',
          abs(gof(kv, 'ratio_abs') - 1.0) < 0.01, kv.get('ratio_abs'))
    check('비율 hf_pred_fluct / hf_actual ~ 0 (요동 0)',
          abs(gof(kv, 'ratio_fluct')) < 1e-6, kv.get('ratio_fluct'))
    check('뒤집힌 요소(folded) 있어도 flip 으로 복구(뒤방향 면적비는 0이 아님)',
          '뒤방향(n_z<0) 면적비' in out, None)
    print('    %s' % [l.strip() for l in out.splitlines() if 'hf_pred_abs' in l][0])
    print('    %s' % [l.strip() for l in out.splitlines() if 'hf_pred_exact ' in l][0])
    print('    %s' % [l.strip() for l in out.splitlines() if 'hf_pred_exact_norm' in l][0])
    print('    %s' % [l.strip() for l in out.splitlines()
                     if 'hf_pred_exact_norm / hf_actual' in l][0])
    print('    %s' % [l.strip() for l in out.splitlines() if 'hf_actual=' in l][0])

    # ================= (iii) 단일 파장 사인파: 요동 < 절대 =================
    print('--- (iii) 사인파 주름: theta_fluct rms < theta_abs rms ---')
    sine = build_odb(nx=40, ny=20, cell=0.1, shape=shape_sine)
    rc, out = run_main([fake], sine, cwd=outdir)
    kv = kv_of(out)
    check('사인파 rc=0', rc == 0, rc)
    ta, tf = gof(kv, 'theta_abs_rms_deg'), gof(kv, 'theta_fluct_rms_deg')
    check('theta_abs rms > 0.5도 (주름 신호 있음)', ta > 0.5, ta)
    check('theta_fluct rms < 0.7 * theta_abs rms (거시/주름 분리 작동)', tf < 0.7 * ta,
          (tf, ta))
    print('    theta_abs rms=%.4f deg  theta_fluct rms=%.4f deg  비=%.3f'
          % (ta, tf, (tf / ta) if ta else float('nan')))
    print('    hf_actual=%.6e  hf_pred_abs=%.6e  ratio_abs=%.4f  ratio_fluct=%.4f'
          % (gof(kv, 'hf_actual'), gof(kv, 'hf_pred_abs'), gof(kv, 'ratio_abs'),
             gof(kv, 'ratio_fluct')))
    check('사인파는 hf_pred_abs 가 hf_actual 과 같은 자릿수(비율 0.5~2)',
          0.5 < gof(kv, 'ratio_abs') < 2.0, kv.get('ratio_abs'))
    check('사인파 배향성비 ~ 1 (법선 수평성분이 +-x 정렬)',
          gof(kv, 'orient_ratio') > 0.999, kv.get('orient_ratio'))

    # ================= (iv) 패치 / 최악패치 / 방위각 히스토그램 =================
    print('--- (iv) 패치 10x10 / 최악 패치 위치 / 방위각 18빈 ---')
    patch_rows = [l for l in out.splitlines() if re.match(r'^\[slope\]\s+y\d+ ', l)]
    check('패치맵 헤더 + 10행 출력', any('패치별 rms [deg]' in l for l in out.splitlines())
          and len(patch_rows) == 10, len(patch_rows))
    check('패치맵 행마다 10칸(접두 2토큰 + 10값)', all(len(r.split()) == 12 for r in patch_rows),
          patch_rows[:1])
    check('최악 패치(위치 포함) 줄 있음',
          any('최악 패치: (ix=' in l and '중심=' in l for l in out.splitlines()), None)
    azline = [l for l in out.splitlines() if '배향성비(2배각)' in l]
    check('방위각 헤더 1줄 + 18빈 목록 1줄',
          len(azline) == 1 and len([l for l in out.splitlines() if '0-20:0.' in l]) == 1, None)
    bump = build_odb(nx=20, ny=20, cell=0.2, shape=shape_bump)
    rcb, outb = run_main([fake], bump, cwd=outdir)
    kvb = kv_of(outb)
    wl = [l for l in outb.splitlines() if '최악 패치: (ix=' in l]
    bx = None
    by = None
    if wl:
        mm = re.search(r'중심=\(([-\d.]+), ([-\d.]+)\)', wl[0])
        if mm:
            bx, by = float(mm.group(1)), float(mm.group(2))
    check('bump fixture rc=0', rcb == 0, rcb)
    near = (bx is not None and by is not None
            and abs(bx - 0.8) < 1.0 and abs(by - 3.2) < 1.0)
    check('최악 패치가 bump(0.8, 3.2) 근처 (tol 1.0 m)', near, (bx, by))
    check('bump fixture 도 패치맵 10행', len([l for l in outb.splitlines()
                                             if re.match(r'^\[slope\]\s+y\d+ ', l)]) == 10, None)
    if wl:
        print('    %s' % wl[0].strip())

    # ================= (v) 스텝/프레임/json/실패경로 =================
    print('--- (v) 스텝/프레임 선택 / JSON / 실패 경로 ---')
    ct = build_odb(nx=6, ny=6, cell=0.5, shape=shape_tilt, steps=('GT', 'CT'))
    rc, out = run_main([fake], ct, cwd=outdir)
    check('Step-Postbuckle 없으면 마지막 스텝(Step-ClampTension)',
          rc == 0 and kv_of(out).get('step') == 'Step-ClampTension'
          and '기본(마지막 step)' in out, kv_of(out).get('step'))
    rc, out = run_main([fake, '--step', 'Step-GlobalTension'], ct, cwd=outdir)
    check('--step 인자 존중 + [인자] 표기',
          rc == 0 and kv_of(out).get('step') == 'Step-GlobalTension' and '[인자]' in out, None)
    rc, out = run_main([fake, '--frame', '1'], tilt, cwd=outdir)
    kf = kv_of(out)
    check('--frame 1 = 미변형 프레임 -> theta=0, rc=0',
          rc == 0 and abs(gof(kf, 'theta_abs_rms_deg')) < 1e-9, kf.get('theta_abs_rms_deg'))
    rc, out = run_main([fake, '--frame', '99'], tilt, cwd=outdir)
    check('범위 밖 --frame -> rc=1 + 이유', rc == 1 and '범위 밖' in out, rc)
    rc, out = run_main([fake, '--step', 'Step-NoSuchStep'], tilt, cwd=outdir)
    check('없는 --step -> rc=1 + 있는 step 목록',
          rc == 1 and 'Step-Postbuckle' in out and 'Step-GlobalTension' in out, rc)
    rc, out = run_main([os.path.join(tmpd, 'NOPE.odb')], tilt, cwd=outdir)
    check('없는 ODB -> rc=1 + 명시 메시지', rc == 1 and 'ODB 파일이 없다' in out, rc)
    for argv, why in (([], '인자 없음'), (['--nope'], '모르는 옵션'),
                      ([fake, '--frame', '0'], 'frame=0'),
                      ([fake, '--frame', 'x'], 'frame 비정수')):
        rc, out = run_main(argv, tilt, cwd=outdir)
        check('사용법 오류(%s) -> rc=2 + usage' % why, rc == 2 and 'usage:' in out, rc)
    rc, out = run_main(['--help'], tilt, cwd=outdir)
    check('--help rc=0 + usage', rc == 0 and 'usage:' in out, rc)
    #   JSON 내용 검사는 **기울기 fixture**(비율이 정의되는 프레임)로 다시 써서 본다.
    rcv, outv = run_main([fake], tilt, cwd=outdir)
    check('JSON 재생성용 기울기 fixture rc=0', rcv == 0, rcv)
    jpath = os.path.join(outdir, 'results', 'slope_HF_Postbuckle_Step-Postbuckle.json')
    check('JSON 기본 경로 저장', os.path.exists(jpath), jpath)
    if os.path.exists(jpath):
        rec = json.load(io.open(jpath, encoding='utf-8'))
        check('JSON 필수 키',
              all(k in rec for k in ('odb', 'step', 'theta_abs', 'theta_fluct', 'hf',
                                     'patches', 'azimuth', 'compression', 'k_slope',
                                     'area_ref', 'n_elem')),
              sorted(rec.keys()))
        check('JSON hf.ratio_exact_norm == 1 (tol 1e-4)',
              abs(rec['hf']['ratio_exact_norm'] - 1.0) < 1e-4, rec['hf']['ratio_exact_norm'])
        check('JSON 패치 격자 10x10 + worst', len(rec['patches']['rms_deg_abs']) == 10
              and rec['patches']['worst_abs'] is not None, None)
        check('JSON 방위각 18빈', len(rec['azimuth']['area_frac']) == 18, None)
        check('JSON odb 신원(절대경로/mtime)', rec['odb']['path'] == os.path.abspath(fake)
              and bool(rec['odb']['mtime']), rec['odb'])
        check('JSON NaN 없음(유한)', 'NaN' not in io.open(jpath, encoding='utf-8').read(), None)
    #   평평 프레임 JSON: 비율이 null 로 나가야 한다(비표준 NaN 토큰 금지).
    rc, outf = run_main([fake, '--frame', '1'], tilt, cwd=outdir)
    txt = io.open(jpath, encoding='utf-8').read()
    recf = json.loads(txt)
    check('평평 프레임 JSON: NaN/Infinity 토큰 0 + ratio null',
          'NaN' not in txt and 'Infinity' not in txt
          and recf['hf']['ratio_abs'] is None and recf['hf']['ratio_exact_norm'] is None,
          recf['hf'].get('ratio_abs'))
    jp2 = os.path.join(tmpd, 'custom.json')
    rc, out = run_main([fake, '--json', jp2], tilt, cwd=outdir)
    check('--json 지정 경로 저장', rc == 0 and os.path.exists(jp2), rc)
    rc, out = run_main([fake], tilt, cwd=outdir)
    check('odb.close() 호출(핸들 누수 없음)', _CUR['odb'].closed is True, None)

    # ================= (vi) 금지문자 / cp949 / 기존 스크립트 무수정 =================
    print('--- (vi) 금지문자 / cp949 / 기존 스크립트 무수정 ---')
    for path in (SCRIPT, os.path.join(CODE, 'slope_stats_stub.py')):
        src = io.open(path, encoding='utf-8').read()
        bad = FORBIDDEN_RE.findall(src)
        check('%s 금지문자 0' % os.path.basename(path), len(bad) == 0, sorted(set(bad)))
    all_out = out + outb + outv + outf
    check('실행 출력에 금지문자 0', len(FORBIDDEN_RE.findall(all_out)) == 0,
          sorted(set(FORBIDDEN_RE.findall(all_out))))
    try:
        all_out.encode('cp949')
        enc = 'OK'
    except Exception as e:
        enc = 'FAIL(%s)' % e
    check('출력 전체 cp949 인코딩 가능(Windows 콘솔)', enc == 'OK', enc)
    for f in ('eval_abaqus.py', 'eval_abaqus_avg.py', 'base_state_probe.py',
              'aba_imperfection.py'):
        #   run_abaqus.py 는 §19.102 의 **수정 대상**이므로 여기서 제외한다(배선 검증은 (ii)).
        r = subprocess.run(['git', 'diff', '--quiet', '--', f], cwd=CODE)
        check('기존 스크립트 무수정: %s' % f, r.returncode == 0, r.returncode)

    print('[stub-test] slope_stats (19.101) %s  (%d/%d)'
          % ('OK' if not FAIL else 'FAIL', N[0] - len(FAIL), N[0]))
    if FAIL:
        print('  실패 항목: %s' % FAIL)
    return 0 if not FAIL else 1


if __name__ == '__main__':
    sys.exit(main())
