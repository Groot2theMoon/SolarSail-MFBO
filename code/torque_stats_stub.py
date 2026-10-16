# -*- coding: utf-8 -*-
"""torque_stats_stub.py - Abaqus 없이 torque_stats.main() 을 **실제 실행**한다 (§19.103)

무엇을 검증하나 (배선/로직/해석앵커 검증. 실제 ODB 물리 검증이 아니다)
    가짜 odbAccess + 가짜 ODB(프레임/인스턴스/COORD bulkDataBlocks)를 끼우고
    torque_stats.py 를 실제로 돌린다. 검증 항목:

      (i)   평평 + 대칭 돛 -> 도심 기준 M ~ 0, tau_rel ~ 0 (가장 강한 sanity)
      (ii)  균일 기울기 fixture -> |M| 이 해석값(지렛대 x 힘 총합)과 일치.
            독립 루프(손계산)로 요소별 f_e 와 M 을 다시 만들어 대조한다.
            균일 법선이면 도심 기준 M = 0 항등식도 함께 확인(힘 방향이 하나면 지렛대가 안 남는다).
      (iii) 비대칭 bump fixture -> tau_rel > 0 이고 면외(M_x,M_y) 성분 지배
      (iv)  대칭 사인파 주름(z(x) 홀함수, y 대칭 메쉬) -> tau_rel ~ 0 (대칭이면 토크도 없다)
      (v)   run_abaqus.py 배선을 **소스에서 떼어** 실제 실행: 정상 경로 / SOLARSAIL_NO_TORQUE=1 경로,
            tau_ 키 append, _archive_result -> _archive_slope -> _archive_torque 순서
      (vi)  금지문자(em/en-dash, 화살표, smart quote, NBSP, bullet, ellipsis) 0 + 출력 cp949 인코딩
            + 기존 5개 스크립트(eval_abaqus/eval_abaqus_avg/base_state_probe/aba_imperfection/
            slope_stats) git diff 무수정

    추가: CLI(--step/--frame/--n-mc/--json), 실패 경로(rc=1/2), JSON(기본/지정 경로),
          odb.close(), 순수함수(dist_stats/tau_rel_of/moments_over_dirs) 단위 검증, 재현성(같은 시드).

사용 (numpy 가 있는 python 으로. Abaqus/라이선스 불필요)
    python3 torque_stats_stub.py
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
SCRIPT = os.path.join(CODE, 'torque_stats.py')

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


class Frame(object):
    def __init__(self, t, coords):
        self.frameValue = t
        labs = np.array(sorted(coords.keys()), dtype=np.int64)
        self.fieldOutputs = {
            'COORD': COORDfo([Blk(labs, np.array([coords[k] for k in sorted(coords.keys())],
                                                 dtype=np.float64))])}


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


def grid_coords(nx, ny, cell, shape=None, shift=(0.0, 0.0)):
    """**원점 중심** 격자(shift 만큼 평행이동). x = (i - nx/2)*cell + shift[0], y = ... + shift[1].

    (실제 돛은 삼각형이라 원점이 도심이 아니다(V3=(0,0,0) 이 꼭짓점) - 그 상황은 shift 로
     재현한다. 대칭 앵커를 만들려면 shift=(0,0) 인 원점 대칭 메쉬가 필요하다.)
    """
    co = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            x = (i - nx / 2.0) * cell + shift[0]
            y = (j - ny / 2.0) * cell + shift[1]
            z = 0.0 if shape is None else shape(x, y)
            co[lab(i, j, ny)] = (x, y, z)
    return co


def elems_of(nx, ny, hybrid=False, folded=()):
    """hybrid=False -> 전부 quad(대칭 보존) / True -> 왼쪽 절반 quad, 오른쪽 절반 tri."""
    out = []
    for i in range(nx):
        for j in range(ny):
            if (not hybrid) or i < nx // 2:
                out.append((lab(i, j, ny), lab(i + 1, j, ny),
                            lab(i + 1, j + 1, ny), lab(i, j + 1, ny)))
            else:
                out.append((lab(i, j, ny), lab(i + 1, j, ny), lab(i + 1, j + 1, ny)))
                out.append((lab(i, j, ny), lab(i + 1, j + 1, ny), lab(i, j + 1, ny)))
    fixed = []
    for k, c in enumerate(out):
        fixed.append(tuple(reversed(c)) if k in folded else c)
    return fixed


def build_odb(nx=8, ny=8, cell=0.5, shape=None, hybrid=False, folded=(),
              steps=('GT', 'PB'), flat_zero=True, shift=(0.0, 0.0)):
    """가짜 ODB. 첫 프레임(frameValue=0) = 미변형(평평). flat_zero=False 면 그 프레임을 뺀다."""
    conns = elems_of(nx, ny, hybrid=hybrid, folded=folded)
    nodes = [Node(k, v) for k, v in sorted(grid_coords(nx, ny, cell).items())]
    elems = [Elem(k + 1, conns[k]) for k in range(len(conns))]
    inst = Inst('MEMBRANE-1', nodes, elems)
    co0 = grid_coords(nx, ny, cell, shift=shift)
    coF = grid_coords(nx, ny, cell, shape, shift=shift)
    st = {}
    if 'GT' in steps:
        st['Step-GlobalTension'] = Step([Frame(0.0, co0), Frame(1.0, coF)])
    if 'PB' in steps:
        fr = [Frame(0.5, coF), Frame(1.0, coF)]
        if flat_zero:
            fr.insert(0, Frame(0.0, co0))
        st['Step-Postbuckle'] = Step(fr)
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
    import torque_stats as M
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
        if ln.startswith('[torque:KV]'):
            d = {}
            for t in ln.split('KV]', 1)[1].split():
                if '=' in t:
                    k, v = t.split('=', 1)
                    d[k] = v
            return d
    return {}


def gof(kv, key):
    """KV 값 float. 접두 없는 이름도 받는다(점별 txt/KV 는 tau_ 접두를 붙인다)."""
    v = kv.get(key)
    if v is None:
        v = kv.get('tau_' + key)
    try:
        return float(v)
    except Exception:
        return float('nan')


# ------------------------------------------------------------------ fixtures
def shape_tilt(x, y):
    """z = tan(5도) x -> 모든 요소 법선이 z_hat 에서 정확히 5도 기울어진 하나의 평면."""
    return math.tan(math.radians(5.0)) * x


def shape_bump(x, y):
    """가우스 bump(비대칭 위치) -> 토크가 0 이 아니고 면외(M_x,M_y) 지배."""
    return 0.10 * math.exp(-((x + 0.80) ** 2 + (y - 0.60) ** 2) / (2.0 * 0.40 ** 2))


def shape_sine(x, y):
    """x 에만 의존하는 홀함수(z(-x) = -z(x)) 주름. y 대칭 메쉬면 M = 0 이어야 한다."""
    return 0.05 * math.sin(2.0 * math.pi * x / 4.0)


# ------------------------------------------------------------------ 독립 손계산
def hand_moment(nx, ny, cell, shape, x_ref_kind='origin', shift=(0.0, 0.0)):
    """요소별 a/n/중심 -> f -> M 을 **독립 루프**로 다시 만든다(해석 대조용).

    x_ref_kind='origin' 이면 원점, 'centroid' 면 같은 루프로 구한 면적가중 도심.
    """
    from torque_stats import A0, R0
    s = np.array([0.0, 0.0, -1.0], dtype=np.float64)
    co = grid_coords(nx, ny, cell, shape, shift=shift)
    conns = elems_of(nx, ny, hybrid=False)
    A = 0.0
    m1 = np.zeros(3)
    rows = []
    for c in conns:
        v = [np.asarray(co[l], dtype=np.float64) for l in c]
        cr = np.cross(v[2] - v[0], v[3] - v[1])
        nn = float(np.linalg.norm(cr))
        a = 0.5 * nn
        n = (cr / nn) if nn > 0.0 else np.zeros(3)
        if n[2] < 0.0:
            n = -n
        cc = max(0.0, min(1.0, float(np.dot(n, -s))))
        f = a * (A0 * cc * s + 2.0 * R0 * cc * cc * (-n))
        xc = 0.25 * (v[0] + v[1] + v[2] + v[3])
        rows.append((a, n, xc, f))
        A += a
        m1 += a * xc
    x_ref = np.zeros(3) if x_ref_kind == 'origin' else (m1 / A)
    M = np.zeros(3)
    sf = 0.0
    for a, n, xc, f in rows:
        M += np.cross(xc - x_ref, f)
        sf += float(np.linalg.norm(f))
    return M, sf, A, x_ref


# ------------------------------------------------------------------ 순수함수 단위검증
def pure_tests():
    import torque_stats as M
    print('--- 순수함수 단위검증 ---')
    # tau_rel_of: |M|=2, sum_f=4, L=1.5 -> 2/(4*1.5)=1/3
    tr, den = M.tau_rel_of(np.array([2.0]), np.array([4.0]), 1.5, 1e-30)
    check('tau_rel_of: |M|/(sum|f|*L_ref) = 2/(4*1.5)', abs(float(tr[0]) - 1.0 / 3.0) < 1e-15, tr)
    check('tau_rel_of: 분모 = sum|f|*L_ref', abs(float(den[0]) - 6.0) < 1e-15, den)
    tr0, _d = M.tau_rel_of(np.array([0.0]), np.array([0.0]), 10.0, 1e-30)
    check('tau_rel_of: 분모 0 -> nan (inf/0-division 금지)', math.isnan(float(tr0[0])), tr0)

    # dist_stats: nan 제거 + 통계
    st = M.dist_stats(np.array([1.0, 3.0, float('nan')]))
    check('dist_stats: nan 제거 n=2, mean=2, min/max', st['n'] == 2 and st['n_nan'] == 1
          and abs(st['mean'] - 2.0) < 1e-15 and st['min'] == 1.0 and st['max'] == 3.0, st)
    check('dist_stats: SE = std/sqrt(n) (ddof=1)', abs(st['se'] - st['std'] / math.sqrt(2.0)) < 1e-15, st)

    # 모멘트: 원점 기준에 요소 하나만 있으면 M = x x f (정의 그대로)
    a = np.array([2.0])
    n = np.array([[0.0, 0.0, 1.0]])
    ce = np.array([[1.0, 2.0, 0.0]])
    s = np.array([[0.0, 0.0, -1.0]])
    r = M.moments_over_dirs(a, n, ce, s, np.zeros(3))
    k = M.IDEAL_THRUST_COEFF
    check('moments_over_dirs: 단일요소 M = x x f (손계산)',
          np.allclose(r['M'][0], np.cross(np.array([1.0, 2.0, 0.0]),
                                          np.array([0.0, 0.0, -2.0 * k])), atol=1e-15),
          r['M'][0])
    check('moments_over_dirs: sum|f| = A*(A0+2R0)', abs(float(r['sum_f'][0]) - 2.0 * k) < 1e-15,
          r['sum_f'][0])
    # 도심 기준으로 옮기면 그 요소에서 M = 0 (지렛대가 없으므로)
    r2 = M.moments_over_dirs(a, n, ce, s, np.array([1.0, 2.0, 0.0]))
    check('moments_over_dirs: x_ref = 요소중심 -> M = 0', abs(float(r2['absM'][0])) < 1e-15,
          r2['absM'][0])

    # dir_at_tilt / dirs_to_phi_psi 왕복
    d = M.dir_at_tilt(30.0, 90.0)
    ph, ps = M.dirs_to_phi_psi(d.reshape(1, 3))
    check('dir_at_tilt -> dirs_to_phi_psi 왕복 (30도, 90도)',
          abs(float(ph[0]) - 30.0) < 1e-9 and abs(float(ps[0]) - 90.0) < 1e-9, (ph, ps))
    s0 = M.dir_at_tilt(0.0, 210.0)
    check('phi=0 -> s_hat = U_SUN(0,0,-1) (수직입사)',
          np.allclose(s0, np.array([0.0, 0.0, -1.0]), atol=1e-15), s0.tolist())
    # 디스크 MC: 단위벡터 + 시드 재현성 + 전력가중(cos 분포)
    s1 = M.sample_dirs_disk(50, seed=7)
    s2 = M.sample_dirs_disk(50, seed=7)
    check('sample_dirs_disk: 단위벡터', np.allclose(np.linalg.norm(s1, axis=1), 1.0, atol=1e-12), None)
    check('sample_dirs_disk: 같은 시드 -> 같은 방향', np.allclose(s1, s2), None)
    check('sample_dirs_disk: z 성분이 모두 <= 0 (반구 전방향)', bool(np.all(s1[:, 2] <= 0)), None)
    grid, gphi, gpsi = M.grid_dirs()
    check('grid_dirs: %d phi x %d psi' % (len(M.PHI_GRID_DEG), len(M.PSI_GRID_DEG)),
          grid.shape == (len(M.PHI_GRID_DEG) * len(M.PSI_GRID_DEG), 3), grid.shape)
    check('grid_dirs: 첫 방향이 phi=0(수직입사) = U_SUN',
          np.allclose(grid[0], np.array([0.0, 0.0, -1.0]), atol=1e-15), grid[0].tolist())
    check('grid_dirs: psi 상한 330/60도 (Heaton 60도 포함)',
          max(gphi) == 60.0 and max(gpsi) == 330.0, (max(gphi), max(gpsi)))

    # scalar_items: 중첩 키 제외 + tau_ 접두
    items = M.scalar_items({'tau_a': 1.0, 'tau_b': float('nan'), 'note': 'x',
                            'grid': {'k': 1}, 'tau_s': 'g', 'tau_flag': True, 'tau_n': 3})
    check('scalar_items: tau_ 접두 + nan 문자열 + 중첩/note 제외',
          items == ['tau_a=1', 'tau_b=nan', 'tau_flag=True', 'tau_n=3', 'tau_s=g'], items)
    js = M._json_safe({'a': float('nan'), 'b': float('inf'), 'c': [1.0, float('nan')],
                       'd': np.float64('2.5'), 'e': 3, 'f': 'x', 'g': True})
    check('_json_safe: NaN/Inf -> None, 나머지 보존',
          js['a'] is None and js['b'] is None and js['c'] == [1.0, None]
          and abs(js['d'] - 2.5) < 1e-12 and js['e'] == 3 and js['f'] == 'x' and js['g'] is True,
          js)


# ------------------------------------------------------------------ §19.103 (v) run_abaqus 배선
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
    """§19.103 로 새로 쓴 run_abaqus.py 줄만 골라낸다(마커 기반 - 커밋 전후 무관하게 결정적)."""
    marks = ('[TORQUE]', '§19.103', '_archive_torque', 'SOLARSAIL_NO_TORQUE', 'tau_')
    return [ln for ln in src.splitlines() if any(m in ln for m in marks)]


def extract_wiring(src_path):
    """run_abaqus.py 에서 _archive_torque 함수 소스만 떼어낸다(다음 def _r13_ct_keys 전까지)."""
    src = io.open(src_path, encoding='utf-8').read()
    m = re.search(r'\ndef _archive_torque\(.*?\n(?=def _r13_ct_keys\()', src, re.S)
    return (m.group(0) if m else None), src


def wiring_tests(tmpd, M):
    print('--- (v) run_abaqus 배선: tau_ 아카이빙 + 환경변수 스위치 + 순서 ---')
    rasrc_path = os.path.join(CODE, 'run_abaqus.py')
    fn_src, rasrc = extract_wiring(rasrc_path)
    check('run_abaqus.py 에서 _archive_torque 소스 추출', bool(fn_src))
    if not fn_src:
        return
    #   호출 순서는 **소스 인덱스**로 확인한다(_archive_result -> _archive_slope -> _archive_torque).
    i_res = rasrc.find("_archive_result(x_c, d_c, 'HF', _rc, _extra)")
    i_slp = rasrc.find('_archive_slope(x_c, d_c)')
    i_tor = rasrc.find('_archive_torque(x_c, d_c)')
    check('소스 순서: _archive_result < _archive_slope < _archive_torque',
          -1 < i_res < i_slp < i_tor, (i_res, i_slp, i_tor))
    check('호출이 _archive_slope **뒤** (점별 txt 재작성 뒤 append 성립)',
          i_slp != -1 and i_tor > i_slp, (i_slp, i_tor))

    cdir = os.path.join(tmpd, 'code')
    rdir = os.path.join(cdir, 'aba')
    os.makedirs(rdir)
    with io.open(os.path.join(rdir, 'HF_Postbuckle.odb'), 'wb') as fh:
        fh.write(b'o' * 2048)
    os.makedirs(os.path.join(cdir, 'results'))
    ptxt = os.path.join(cdir, 'results', 'results_xc0p5_dc0p3.txt')
    #   _archive_result('w') -> _archive_slope append 를 흉내낸 상태
    base = 'x_c=0.5\nd_c=0.3\nfidelity=HF\nhf=3.9e-05\nslope_f1_hf_pred=4.2e-05\n'
    with io.open(ptxt, 'w', encoding='utf-8') as fh:
        fh.write(base)

    logs = []

    def _emit(*a):
        logs.append(' '.join(str(x) for x in a))

    from aba_imperfection import value_tag
    ns = {'os': os, 'sys': sys, 'io': io, 'shutil': shutil, 'math': math, 'time': time,
          'value_tag': value_tag, 'emit': _emit, 'torque_stats': M,
          '_RUN': rdir, '_HERE': cdir, 'HF_ODB': 'HF_Postbuckle.odb'}
    exec(compile(fn_src, '<run_abaqus:_archive_torque>', 'exec'), ns)
    fn = ns['_archive_torque']

    jp = os.path.join(cdir, 'results', 'torque_HF_Postbuckle_Step-Postbuckle.json')

    # (a) 정상 경로 (스위치 미설정)
    _unsetenv('SOLARSAIL_NO_TORQUE')
    _CUR['odb'] = build_odb(nx=8, ny=8, cell=0.5, shape=shape_tilt, hybrid=True, folded=(3,))
    with contextlib.redirect_stdout(io.StringIO()):
        fn(0.5, 0.3)
    txt = io.open(ptxt, encoding='utf-8').read()
    check('정상: 점별 txt 에 tau_ 접두 key=value append',
          'tau_rel_plane=' in txt and 'tau_abs=' in txt and 'tau_mx=' in txt
          and 'tau_worst_phi_deg=' in txt, None)
    check('정상: 기존 RESULTS 계약 필드 + slope_ 줄 보존',
          'fidelity=HF' in txt and 'hf=3.9e-05' in txt and 'slope_f1_hf_pred=' in txt, None)
    check('정상: tau_ 줄이 slope_ 줄 **뒤**에 붙는다(append 순서)',
          txt.find('slope_f1_hf_pred=') < txt.find('tau_rel_plane='), None)
    check('정상: JSON 저장 (torque_HF_Postbuckle_Step-Postbuckle.json)', os.path.exists(jp), jp)
    check('정상: [TORQUE] 로그에 key 개수와 tau_rel_plane/cosw/max',
          any(l.startswith('[TORQUE]') and 'tau_ 접두' in l and 'tau_rel_plane=' in l
              for l in logs), logs[-3:])
    parsed = {}
    for ln in txt.splitlines():
        if ln.startswith('tau_') and '=' in ln:
            k, v = ln.split('=', 1)
            parsed[k] = v
    check('txt 파서: tau_ 키 >= 25 개', len(parsed) >= 25, len(parsed))
    check('txt 파서: 필수 키 9종 존재',
          all(k in parsed for k in ('tau_rel_plane', 'tau_rel_cosw', 'tau_rel_max', 'tau_abs',
                                    'tau_mx', 'tau_my', 'tau_mz', 'tau_worst_phi_deg',
                                    'tau_flat_rel')),
          sorted(parsed.keys())[:6])
    check('txt 파서: 수치 키가 수치로 파싱',
          all(_isnum(parsed.get(k, 'x')) for k in
              ('tau_rel_plane', 'tau_rel_cosw', 'tau_rel_max', 'tau_abs', 'tau_mx',
               'tau_my', 'tau_mz', 'tau_worst_phi_deg', 'tau_flat_rel', 'l_ref',
               'area_ref', 'n_elem')),
          {k: parsed.get(k) for k in ('tau_rel_plane', 'tau_flat_rel')})
    #   (ii) 해석 앵커를 **배선 경로에서도** 확인: 균일 기울기 = 하나의 법선 -> 도심 기준 M = 0
    check('배선 경로(균일 기울기): tau_abs ~ 0 (도심 기준)',
          abs(float(parsed.get('tau_abs', 'nan'))) < 1e-9, parsed.get('tau_abs'))
    if os.path.exists(jp):
        rec = json.load(io.open(jp, encoding='utf-8'))
        check('JSON section=19.103 + tau 블록 + sweep(grid/disk) + flat_anchor',
              rec.get('section') == '19.103' and 'tau' in rec and 'grid' in rec.get('sweep', {})
              and 'disk_mc' in rec.get('sweep', {}) and 'flat_anchor' in rec,
              sorted(rec.keys()))

    # (b) SOLARSAIL_NO_TORQUE=1 -> 배선 전체 생략
    if os.path.exists(jp):
        os.remove(jp)
    with io.open(ptxt, 'w', encoding='utf-8') as fh:
        fh.write(base)
    logs[:] = []
    os.environ['SOLARSAIL_NO_TORQUE'] = '1'
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            fn(0.5, 0.3)
    finally:
        _unsetenv('SOLARSAIL_NO_TORQUE')
    txt = io.open(ptxt, encoding='utf-8').read()
    check('NO_TORQUE=1: tau_ 줄 없음', 'tau_rel_plane=' not in txt)
    check('NO_TORQUE=1: JSON 안 씀', not os.path.exists(jp))
    check('NO_TORQUE=1: SKIP 로그 1줄',
          any('SOLARSAIL_NO_TORQUE=1' in l and 'SKIP' in l for l in logs), logs)

    # (c) ODB 가 없으면 경고 1줄(런 사망 금지)
    logs[:] = []
    os.remove(os.path.join(rdir, 'HF_Postbuckle.odb'))
    with io.open(ptxt, 'w', encoding='utf-8') as fh:
        fh.write(base)
    with contextlib.redirect_stdout(io.StringIO()):
        fn(0.5, 0.3)
    check('ODB 없음: [TORQUE] 경고 1줄 + tau_ 줄 없음',
          any(l.startswith('[TORQUE] !!!') for l in logs)
          and 'tau_rel_plane=' not in io.open(ptxt, encoding='utf-8').read(), logs)
    with io.open(os.path.join(rdir, 'HF_Postbuckle.odb'), 'wb') as fh:
        fh.write(b'o' * 2048)

    # (iv) 소스 문자열 검사
    check('(iv) 추출 배선 소스 금지문자 0', len(FORBIDDEN_RE.findall(fn_src)) == 0,
          sorted(set(FORBIDDEN_RE.findall(fn_src))))
    _add = '\n'.join(_new_run_abaqus_lines(rasrc))
    check('(iv) run_abaqus.py **추가분** 금지문자 0 (기존 줄은 불변)',
          len(FORBIDDEN_RE.findall(_add)) == 0 and len(_add) > 0,
          sorted(set(FORBIDDEN_RE.findall(_add))))
    check('(iv) 점별 txt cp949 인코딩 가능 (Windows 콘솔/파일)', _enc_ok(txt), None)
    check('(iv) all_results.csv 열 스키마 불변 (tau 미포함, hf 마지막)',
          "'lf1', 'lf2', 'lf3', 'hf']" in _csv_cols_block(rasrc)
          and 'tau' not in _csv_cols_block(rasrc), None)
    check('(iv) 실패(HF_FAIL) 경로에 _archive_torque 호출 없음',
          '_archive_torque' not in rasrc.split('HF_FAIL')[1][:400], None)


# ------------------------------------------------------------------ main
def main():
    tmpd = tempfile.mkdtemp(prefix='torquestub_')
    fake = os.path.join(tmpd, 'HF_Postbuckle.odb')
    with io.open(fake, 'wb') as fh:
        fh.write(b'x' * 4096)
    install_fake_abaqus()
    import torque_stats as M                                             # noqa: E402

    print('--- 0) 컴파일 / 상수 ---')
    p = subprocess.run([sys.executable, '-m', 'py_compile', SCRIPT],
                       capture_output=True, text=True, cwd=CODE)
    check('py_compile', p.returncode == 0, p.stderr[-300:])
    check('평판 이상계수 = A0 + 2*R0 (참조구현과 동일 상수)',
          abs(M.IDEAL_THRUST_COEFF - (M.A0 + 2.0 * M.R0)) < 1e-15, M.IDEAL_THRUST_COEFF)
    print('    A0=%.6f R0=%.6f IDEAL=%.6f U_SUN=%s'
          % (M.A0, M.R0, M.IDEAL_THRUST_COEFF, np.asarray(M.U_SUN).tolist()))

    pure_tests()

    outdir = os.path.join(tmpd, 'run')
    os.makedirs(outdir)

    # ================= (i) 평평 + 대칭 돛 =================
    print('--- (i) 평평 + 대칭 돛: M ~ 0 / tau_rel ~ 0 ---')
    flat = build_odb(nx=8, ny=8, cell=0.5, shape=None, hybrid=True, folded=(3,))
    rc, out = run_main([fake], flat, cwd=outdir)
    kv = kv_of(out)
    check('평평 돛 rc=0', rc == 0, rc)
    check('도심 기준 |M| < 1e-9 (모멘트 상쇄)', abs(gof(kv, 'tau_abs')) < 1e-9, kv.get('tau_abs'))
    check('원점 기준 |M| < 1e-9 (원점 대칭 메쉬 -> 지렛대도 상쇄)',
          'tau_abs_origin' in kv and abs(gof(kv, 'tau_abs_origin')) < 1e-9, kv.get('tau_abs_origin'))
    check('tau_rel_plane < 1e-12', abs(gof(kv, 'tau_rel_plane')) < 1e-12, kv.get('tau_rel_plane'))
    check('tau_rel_max < 1e-12 (평판은 모든 방향에서 토크 0)', abs(gof(kv, 'tau_rel_max')) < 1e-12,
          kv.get('tau_rel_max'))
    check('tau_mx/tau_my/tau_mz ~ 0',
          max(abs(gof(kv, 'tau_mx')), abs(gof(kv, 'tau_my')), abs(gof(kv, 'tau_mz'))) < 1e-9,
          (kv.get('tau_mx'), kv.get('tau_my'), kv.get('tau_mz')))
    check('tau_flat_rel < 1e-12 (frame 0 앵커, 도심 기준)', abs(gof(kv, 'tau_flat_rel')) < 1e-12,
          kv.get('tau_flat_rel'))
    check('tau_flat_rel(원점) < 1e-12 (원점 대칭)', abs(gof(kv, 'tau_flat_rel_origin')) < 1e-12,
          kv.get('tau_flat_rel_origin'))
    print('    %s' % [l.strip() for l in out.splitlines() if '수직입사(도심기준): |M|' in l][0])
    print('    %s' % [l.strip() for l in out.splitlines() if 'tau_rel = |M|' in l][0])

    # ================= (ii) 균일 기울기 = 해석값 =================
    print('--- (ii) 균일 기울기 fixture: |M| = 지렛대 x 힘 총합(해석) ---')
    nx, ny, cell = 8, 8, 0.5
    #   실제 돛은 삼각형(V3=(0,0,0) 이 꼭짓점)이라 원점이 도심이 아니다 -> shift 로 그 상황을
    #   재현한다. 평행이동해 두면 '원점 기준' 지렛대 항이 0 이 아니게 되어 항등식이 살아난다.
    SH = (10.0, 3.0)
    tilt = build_odb(nx=nx, ny=ny, cell=cell, shape=shape_tilt, shift=SH)
    rc, out = run_main([fake], tilt, cwd=outdir)
    kv = kv_of(out)
    check('균일 기울기(평행이동) rc=0', rc == 0, rc)
    Mh, sfh, ah, xrefh = hand_moment(nx, ny, cell, shape_tilt, 'origin', shift=SH)
    hs = float(np.linalg.norm(Mh))
    check('원점 기준 |M| == 손계산 |M| = |sum x_e x f_e| (tol 1e-12 상대)',
          hs > 1.0 and abs(gof(kv, 'tau_abs_origin') - hs) <= 1e-12 * hs,
          (gof(kv, 'tau_abs_origin'), hs))
    check('sum|f_e| == 손계산 (tol 1e-11 상대: KV 는 %.12g 로 인쇄된다)',
          abs(gof(kv, 'sum_f_plane') - sfh) < 1e-11 * sfh,
          (gof(kv, 'sum_f_plane'), sfh))
    check('L_ref == sqrt(area_ref) (tol 1e-15)',
          abs(gof(kv, 'l_ref') - math.sqrt(gof(kv, 'area_ref'))) < 1e-15 * gof(kv, 'area_ref'),
          (gof(kv, 'l_ref'), gof(kv, 'area_ref')))
    check('tau_rel(원점) == |M_origin|/(sum|f|*L_ref) (정의 항등식, tol 1e-9 상대)',
          abs(gof(kv, 'tau_rel_origin_plane')
              - gof(kv, 'tau_abs_origin') / (gof(kv, 'sum_f_plane') * gof(kv, 'l_ref')))
          < 1e-9 * gof(kv, 'tau_rel_origin_plane'),
          kv.get('tau_rel_origin_plane'))
    #   인쇄 정밀도(%.12g)에 가리지 않는 **전정밀** 대조: compute_stats 반환값으로 다시 본다.
    Mc, sfc, ac, xrefc = hand_moment(nx, ny, cell, shape_tilt, 'centroid', shift=SH)
    _CUR['odb'] = tilt
    with contextlib.redirect_stdout(io.StringIO()):
        S = M.compute_stats(fake)
    check('compute_stats -> dict (CLI 없이 직접 호출)', isinstance(S, dict), type(S))

    def _rel(a, b, tol=1e-13):
        return abs(float(a) - float(b)) <= tol * max(abs(float(b)), 1.0)

    check('전정밀: |M_origin| == 손계산 (tol 1e-13 상대)', _rel(S['tau_abs_origin'], hs),
          (S['tau_abs_origin'], hs))
    check('전정밀: sum|f| == 손계산 (tol 1e-13 상대)', _rel(S['sum_f_plane'], sfh),
          (S['sum_f_plane'], sfh))
    check('전정밀: tau_rel(원점) == |M|/(sum|f|*L_ref) (tol 1e-13 상대)',
          _rel(S['tau_rel_origin_plane'],
               S['tau_abs_origin'] / (S['sum_f_plane'] * S['l_ref'])), None)
    check('전정밀: l_ref == sqrt(area_ref)',
          abs(S['l_ref'] - math.sqrt(S['area_ref'])) <= 1e-15 * S['area_ref'], None)
    check('전정밀: 도심 == 손계산 도심 (tol 1e-13 상대)',
          _rel(S['x_ref_x'], float(xrefc[0])) and _rel(S['x_ref_y'], float(xrefc[1])),
          (S['x_ref_x'], S['x_ref_y'], list(xrefc)))
    check('전정밀: tau_rel_max >= tau_rel_plane (격자에 수직입사가 포함)',
          S['tau_rel_max'] >= S['tau_rel_plane'] - 1e-15, None)
    check('전정밀: min <= p50 <= p95 <= max (격자 분포 단조)',
          0.0 <= S['tau_rel_grid_min'] <= S['tau_rel_grid_p50'] <= S['tau_rel_grid_p95']
          <= S['tau_rel_grid_max'], None)
    #   균일 법선(힘 방향 하나) -> **도심 기준** M = 0 항등식. 평행이동 뒤에도 성립해야 한다.
    check('손계산(도심 기준) M ~ 0 - 균일 법선이면 지렛대가 남지 않는다',
          float(np.linalg.norm(Mc)) < 1e-9, float(np.linalg.norm(Mc)))
    check('구현(도심 기준) |M| ~ 0 (동일 항등식, 평행이동 무관)',
          abs(gof(kv, 'tau_abs')) < 1e-9, kv.get('tau_abs'))
    check('메쉬 도심 == 손계산 도심 (tol 1e-12)',
          abs(gof(kv, 'x_ref_x') - float(xrefc[0])) <= 1e-12 * max(abs(float(xrefc[0])), 1.0)
          and abs(gof(kv, 'x_ref_y') - float(xrefc[1])) <= 1e-12 * max(abs(float(xrefc[1])), 1.0),
          (gof(kv, 'x_ref_x'), gof(kv, 'x_ref_y'), list(xrefc)))
    print('    손계산 |M_origin|=%.10e N m  sum|f|=%.10e N  A=%.6f m^2  도심=%s'
          % (hs, sfh, ah, ['%.6f' % v for v in xrefc]))
    print('    %s' % [l.strip() for l in out.splitlines() if '수직입사(원점기준)' in l][0])

    # 원점 대칭 배치(shift 없음)에서는 지렛대 항도 상쇄된다 - 기준점 선택의 효과를 나란히 본다.
    rc, out = run_main([fake], build_odb(nx=nx, ny=ny, cell=cell, shape=shape_tilt), cwd=outdir)
    kv0 = kv_of(out)
    check('원점 중심 동일 기울기 -> tau_rel(원점) ~ 0 (지렛대 상쇄)',
          rc == 0 and abs(gof(kv0, 'tau_rel_origin_plane')) < 1e-12,
          kv0.get('tau_rel_origin_plane'))
    check('평행이동하면 tau_rel(원점) 이 커진다(> 1) - 기준점이 만드는 항',
          gof(kv, 'tau_rel_origin_plane') > 1.0, kv.get('tau_rel_origin_plane'))
    check('두 경우 모두 tau_rel(도심) ~ 0 (형상 몫만 남는다)',
          abs(gof(kv0, 'tau_rel_plane')) < 1e-12 and abs(gof(kv, 'tau_rel_plane')) < 1e-12,
          (kv0.get('tau_rel_plane'), kv.get('tau_rel_plane')))

    # ================= (iii) 비대칭 bump =================
    print('--- (iii) 비대칭 bump: tau_rel > 0 + 면외 지배 ---')
    bump = build_odb(nx=20, ny=20, cell=0.2, shape=shape_bump)
    rc, out = run_main([fake], bump, cwd=outdir)
    kvb = kv_of(out)
    check('bump rc=0', rc == 0, rc)
    check('tau_rel_plane > 1e-4 (형상이 토크를 만든다)', gof(kvb, 'tau_rel_plane') > 1e-4,
          kvb.get('tau_rel_plane'))
    check('면외 비중 > 0.5 (M_x,M_y 지배)', gof(kvb, 'tau_out_share') > 0.5, kvb.get('tau_out_share'))
    check('면내 비중 < 0.1 (M_z 는 거의 없다)', gof(kvb, 'tau_in_share') < 0.1, kvb.get('tau_in_share'))
    check('|M_z| < |M_xy| (면내가 더 작다)',
          abs(gof(kvb, 'tau_mz')) < abs(gof(kvb, 'tau_mout')), (kvb.get('tau_mz'), kvb.get('tau_mout')))
    check('tau_rel_max >= tau_rel_plane (격자가 수직입사를 포함)', gof(kvb, 'tau_rel_max')
          >= gof(kvb, 'tau_rel_plane') - 1e-15,
          (kvb.get('tau_rel_max'), kvb.get('tau_rel_plane')))
    check('tau_rel_cosw > 0 (전력가중 MC 평균도 0 이 아니다)', gof(kvb, 'tau_rel_cosw') > 0.0,
          kvb.get('tau_rel_cosw'))
    check('MC SE <= std (정의)', gof(kvb, 'tau_rel_cosw_se') <= gof(kvb, 'tau_rel_cosw_std') + 1e-30,
          (kvb.get('tau_rel_cosw_se'), kvb.get('tau_rel_cosw_std')))
    check('격자 nan 0 (0~60도에서는 분모가 살아 있다)', gof(kvb, 'tau_grid_n_nan') == 0,
          kvb.get('tau_grid_n_nan'))
    print('    %s' % [l.strip() for l in out.splitlines() if '면외(M_x,M_y)' in l][0])
    print('    %s' % [l.strip() for l in out.splitlines() if '방향 격자(' in l][0])
    print('    %s' % [l.strip() for l in out.splitlines() if '디스크 MC(' in l][0])

    # ================= (iv) 대칭 사인파 =================
    print('--- (iv) 대칭 사인파 주름: tau_rel ~ 0 ---')
    sine = build_odb(nx=40, ny=20, cell=0.1, shape=shape_sine)
    rc, out = run_main([fake], sine, cwd=outdir)
    kvs = kv_of(out)
    check('사인파 rc=0', rc == 0, rc)
    check('tau_abs < 1e-9 (z 홀함수 + y 대칭 -> M = 0)', abs(gof(kvs, 'tau_abs')) < 1e-9,
          kvs.get('tau_abs'))
    check('tau_rel_plane < 1e-12', abs(gof(kvs, 'tau_rel_plane')) < 1e-12, kvs.get('tau_rel_plane'))
    check('tau_rel_max < 1e-12 (모든 방향에서 대칭이 유지된다)',
          abs(gof(kvs, 'tau_rel_max')) < 1e-12, kvs.get('tau_rel_max'))
    check('tau_rel_cosw < 1e-12', abs(gof(kvs, 'tau_rel_cosw')) < 1e-12, kvs.get('tau_rel_cosw'))
    check('대칭 주름: tau_rel_cosw ~ 0 (모든 방향 평균도 0)',
          abs(gof(kvs, 'tau_rel_cosw')) < 1e-12, kvs.get('tau_rel_cosw'))
    #   두 축이 다르다는 대조: 같은 대칭 주름 ODB 로 hf(참조 구현)를 재면 0 이 아니다.
    import slope_stats as SS
    _CUR['odb'] = sine
    with contextlib.redirect_stdout(io.StringIO()):
        SS_ = SS.compute_stats(fake)
    check('대칭 주름: hf_actual > 1e-6 (추력손실은 있다) 인데 tau_rel ~ 0 (두 축이 다르다)',
          SS_ is not None and float(SS_['hf_actual']) > 1e-6
          and abs(gof(kvs, 'tau_rel_plane')) < 1e-12,
          (None if SS_ is None else SS_['hf_actual'], kvs.get('tau_rel_plane')))
    print('    %s' % [l.strip() for l in out.splitlines() if '수직입사(도심기준): |M|' in l][0])
    print('    %s' % [l.strip() for l in out.splitlines() if 'tau_rel_cosw' in l and 'append' not in l][-1])

    # ================= (v) 스텝/프레임/n-mc/json/실패경로 =================
    print('--- (v) 스텝/프레임 --n-mc / JSON / 실패 경로 ---')
    ct = build_odb(nx=6, ny=6, cell=0.5, shape=shape_tilt, steps=('GT', 'CT'))
    rc, out = run_main([fake], ct, cwd=outdir)
    check('Step-Postbuckle 없으면 마지막 스텝(Step-ClampTension)',
          rc == 0 and kv_of(out).get('step') == 'Step-ClampTension'
          and '기본(마지막 step)' in out, kv_of(out).get('step'))
    rc, out = run_main([fake, '--step', 'Step-GlobalTension'], ct, cwd=outdir)
    check('--step 인자 존중 + [인자] 표기',
          rc == 0 and kv_of(out).get('step') == 'Step-GlobalTension' and '[인자]' in out, None)
    rc, out = run_main([fake, '--frame', '1'], tilt, cwd=outdir)
    check('--frame 1 = 미변형 프레임 -> tau_rel 0, rc=0',
          rc == 0 and abs(gof(kv_of(out), 'tau_rel_plane')) < 1e-12,
          kv_of(out).get('tau_rel_plane'))
    rc, out = run_main([fake, '--frame', '99'], tilt, cwd=outdir)
    check('범위 밖 --frame -> rc=1 + 이유', rc == 1 and '범위 밖' in out, rc)
    rc, out = run_main([fake, '--step', 'Step-NoSuchStep'], tilt, cwd=outdir)
    check('없는 --step -> rc=1 + 있는 step 목록',
          rc == 1 and 'Step-Postbuckle' in out and 'Step-GlobalTension' in out, rc)
    rc, out = run_main([os.path.join(tmpd, 'NOPE.odb')], tilt, cwd=outdir)
    check('없는 ODB -> rc=1 + 명시 메시지', rc == 1 and 'ODB 파일이 없다' in out, rc)
    rc, out = run_main([fake, '--n-mc', '150'], bump, cwd=outdir)
    check('--n-mc 150 -> tau_n_mc=150 / tau_disk_n=150', rc == 0
          and kv_of(out).get('tau_n_mc') == '150', kv_of(out).get('tau_n_mc'))
    #   같은 시드 재현성
    rc, o1 = run_main([fake, '--n-mc', '300'], bump, cwd=outdir)
    rc, o2 = run_main([fake, '--n-mc', '300'], bump, cwd=outdir)
    check('같은 시드 -> tau_rel_cosw 재현(결정적)',
          kv_of(o1).get('tau_rel_cosw') == kv_of(o2).get('tau_rel_cosw'),
          (kv_of(o1).get('tau_rel_cosw'), kv_of(o2).get('tau_rel_cosw')))
    for argv, why in (([], '인자 없음'), (['--nope'], '모르는 옵션'),
                      ([fake, '--frame', '0'], 'frame=0'),
                      ([fake, '--frame', 'x'], 'frame 비정수'),
                      ([fake, '--n-mc', '0'], 'n-mc=0'),
                      ([fake, '--n-mc', 'x'], 'n-mc 비정수')):
        rc, out = run_main(argv, tilt, cwd=outdir)
        check('사용법 오류(%s) -> rc=2 + usage' % why, rc == 2 and 'usage:' in out, rc)
    rc, out = run_main(['--help'], tilt, cwd=outdir)
    check('--help rc=0 + usage', rc == 0 and 'usage:' in out, rc)
    jpath = os.path.join(outdir, 'results', 'torque_HF_Postbuckle_Step-Postbuckle.json')
    rcv, outv = run_main([fake], bump, cwd=outdir)
    check('JSON 기본 경로 저장', rcv == 0 and os.path.exists(jpath), jpath)
    if os.path.exists(jpath):
        rec = json.load(io.open(jpath, encoding='utf-8'))
        need = ('odb', 'step', 'frame_index', 'n_elem', 'area_ref', 'l_ref', 'x_ref_centroid',
                'tau', 'flat_anchor', 'sweep', 'moment_plane_centroid', 'moment_plane_origin',
                'note', 'section')
        check('JSON 필수 키', all(k in rec for k in need), sorted(rec.keys()))
        check('JSON tau 블록에 성분/비율/최악 방향',
              all(k in rec['tau'] for k in ('abs_plane_centroid', 'rel_plane_centroid', 'mx', 'my',
                                            'mz', 'out_share', 'in_share', 'rel_max_all_dirs',
                                            'worst_phi_deg', 'worst_psi_deg', 'rel_mean_disk')),
              sorted(rec['tau'].keys()))
        check('JSON sweep.grid.tau_rel 길이 = phi*psi (%d x %d)'
              % (len(M.PHI_GRID_DEG), len(M.PSI_GRID_DEG)),
              len(rec['sweep']['grid']['tau_rel']) == len(M.PHI_GRID_DEG) * len(M.PSI_GRID_DEG),
              len(rec['sweep']['grid']['tau_rel']))
        check('JSON sweep.disk_mc n/seed + stats',
              rec['sweep']['disk_mc']['n'] == M.DEFAULT_N_MC
              and rec['sweep']['disk_mc']['seed'] == M.TORQUE_SEED
              and 'mean' in rec['sweep']['disk_mc']['stats'], rec['sweep']['disk_mc']['n'])
        check('JSON sweep.sampler 에 phi 상한/블록/바닥값',
              all(k in rec['sweep']['sampler'] for k in ('phi_grid_deg', 'psi_grid_deg', 'n_dirs',
                                                         'n_mc', 'seed', 'block',
                                                         'sum_f_floor_n', 'disk_mode')),
              sorted(rec['sweep']['sampler'].keys()))
        check('JSON odb 신원(절대경로/mtime)', rec['odb']['path'] == os.path.abspath(fake)
              and bool(rec['odb']['mtime']), rec['odb'])
        check('JSON NaN/Infinity 비표준 토큰 0',
              'NaN' not in io.open(jpath, encoding='utf-8').read()
              and 'Infinity' not in io.open(jpath, encoding='utf-8').read(), None)
    jp2 = os.path.join(tmpd, 'custom_torque.json')
    rc, out = run_main([fake, '--json', jp2], bump, cwd=outdir)
    check('--json 지정 경로 저장', rc == 0 and os.path.exists(jp2), rc)
    #   frameValue=0 프레임이 없는 ODB -> tau_flat 은 nan (0 으로 채우지 않는다)
    noz = build_odb(nx=6, ny=6, cell=0.5, shape=shape_tilt, flat_zero=False)
    for _sn in list(noz.steps.keys()):
        noz.steps[_sn].frames = [f for f in noz.steps[_sn].frames
                                 if abs(float(f.frameValue)) > 1e-12]
    rc, out = run_main([fake], noz, cwd=outdir)
    kzn = kv_of(out)
    check('frameValue=0 없음 -> rc=0 + tau_flat_rel=nan + 경고 로그',
          rc == 0 and math.isnan(gof(kzn, 'tau_flat_rel')) and 'tau_flat 은 nan' in out,
          (rc, kzn.get('tau_flat_rel')))
    rc, out = run_main([fake], bump, cwd=outdir)
    check('odb.close() 호출(핸들 누수 없음)', _CUR['odb'].closed is True, None)

    # ================= (vi) 금지문자 / cp949 / 기존 스크립트 무수정 =================
    print('--- (vi) 금지문자 / cp949 / 기존 스크립트 무수정 ---')
    for path in (SCRIPT, os.path.join(CODE, 'torque_stats_stub.py')):
        src = io.open(path, encoding='utf-8').read()
        bad = FORBIDDEN_RE.findall(src)
        check('%s 금지문자 0' % os.path.basename(path), len(bad) == 0, sorted(set(bad)))
    all_out = out + outv + o1 + o2
    check('실행 출력에 금지문자 0', len(FORBIDDEN_RE.findall(all_out)) == 0,
          sorted(set(FORBIDDEN_RE.findall(all_out))))
    try:
        all_out.encode('cp949')
        enc = 'OK'
    except Exception as e:
        enc = 'FAIL(%s)' % e
    check('출력 전체 cp949 인코딩 가능(Windows 콘솔)', enc == 'OK', enc)
    for f in ('eval_abaqus.py', 'eval_abaqus_avg.py', 'base_state_probe.py',
              'aba_imperfection.py', 'slope_stats.py'):
        #   run_abaqus.py 는 §19.103 의 **수정 대상**이므로 여기서 제외한다(배선 검증은 (v)).
        r = subprocess.run(['git', 'diff', '--quiet', '--', f], cwd=CODE)
        check('기존 스크립트 무수정: %s' % f, r.returncode == 0, r.returncode)

    print('[stub-test] torque_stats (19.103) %s  (%d/%d)'
          % ('OK' if not FAIL else 'FAIL', N[0] - len(FAIL), N[0]))
    if FAIL:
        print('  실패 항목: %s' % FAIL)
    return 0 if not FAIL else 1


if __name__ == '__main__':
    sys.exit(main())
