#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dn_stats_stub.py - Abaqus 없이 delta_n(eval_abaqus.py) 과 DN: 아카이빙(run_abaqus.py)을 **실제 실행**한다 (19.106)

무엇을 검증하나 (배선/로직/해석앵커. 실제 ODB 물리 검증이 아니다)
    가짜 odbAccess + 가짜 ODB(프레임/인스턴스/COORD bulkDataBlocks)를 끼우고
    eval_abaqus.main() 을 HFONLY 모드로 **실제로 돌린다**(get_hf 까지 포함한 전 경로).

      (a) 완전 평판              -> dn_deg == 0        (rel < 1e-9)
      (b) 균일 기울기 alpha=5deg -> dn_deg == 5.0      (rel < 1e-9), dn_psi = 기울기 방위(180deg)
          (추가) 45deg 방위 기울기 -> dn_psi == -135deg, (추가) 뒤집힌 요소 1개 -> dn 불변
      (c) 대칭 사인파 주름(평균경사 0) -> dn_deg == 0 (rel < 1e-9), 그러나 hf > 0
          => delta_n 은 '일관 기울기'만 잡고 국소 주름은 상쇄한다(=hf 와 독립 정보)
      (d) 비대칭 원뿔형 처짐     -> 0 < dn_deg < alpha_max (=atan(h/R))
      (e) DN:FAIL 경로           -> 'DN:FAIL reason=..' 한 줄, RESULTS/extraction.txt 는 평소대로
      (f) RESULTS 회귀           -> git HEAD 판 eval_abaqus.py 와 RESULTS 줄 **바이트 동일**
      (g) run_abaqus 배선        -> _archive_dn 소스를 떼어 실행: dn_ 키 append,
                                    all_results.csv **바이트 불변**, _cols 블록 불변,
                                    DN 수집이 FAIL 검사보다 **앞**, 호출이 _archive_result 뒤,
                                    DN 줄 없음 -> 경고 1줄(조용히 넘기지 않음)

사용 (numpy 가 있는 python 으로. Abaqus/라이선스 불필요)
    python3 dn_stats_stub.py
종료코드: 0 OK / 1 FAIL / 2 numpy 없음
"""
import contextlib
import importlib.util
import io
import math
import os
import re
import subprocess
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
CODE = HERE

if sys.version_info[0] < 3:
    print('  [stub] py3 필요')
    sys.exit(2)
try:
    import numpy as np
except ImportError:
    print('  [stub] numpy 가 필요하다(Abaqus 대신 numpy 있는 python 으로 실행).')
    sys.exit(2)

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


def relerr(a, b):
    d = abs(b) if abs(b) > 1e-12 else 1.0
    return abs(a - b) / d


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
        self.close_calls = 0

    def close(self):
        self.closed = True
        self.close_calls += 1


def lab(i, j, ny):
    return i * (ny + 1) + j + 1


def grid_coords(nx, ny, cell, shape=None):
    """원점 중심 격자. shape=None 이면 완전 평평."""
    co = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            x = (i - nx / 2.0) * cell
            y = (j - ny / 2.0) * cell
            z = 0.0 if shape is None else shape(x, y)
            co[lab(i, j, ny)] = (x, y, z)
    return co


def elems_of(nx, ny, folded=()):
    out = []
    for i in range(nx):
        for j in range(ny):
            out.append((lab(i, j, ny), lab(i + 1, j, ny),
                        lab(i + 1, j + 1, ny), lab(i, j + 1, ny)))
    return [tuple(reversed(c)) if k in folded else c for k, c in enumerate(out)]


def build_odb(nx=8, ny=8, cell=0.5, shape=None, folded=(), with_pb=True, conn_pad=False):
    conns = elems_of(nx, ny, folded=folded)
    if conn_pad:
        #   5절점 connectivity: calc_thrust_loss/calc_normal_deviation 둘 다 quad/tri 필터에서
        #   빠진다. hf 는 0(=1e-12)로 '성공'하고 DN 만 실패 -> DN:FAIL 격리 검증용.
        conns = [tuple(c) + (c[0],) for c in conns]
    nodes = [Node(k, v) for k, v in sorted(grid_coords(nx, ny, cell).items())]
    elems = [Elem(k + 1, conns[k]) for k in range(len(conns))]
    inst = Inst('MEMBRANE-1', nodes, elems)
    co0 = grid_coords(nx, ny, cell)
    coF = grid_coords(nx, ny, cell, shape)
    st = {'Step-GlobalTension': Step([Frame(0.0, co0), Frame(1.0, coF)])}
    if with_pb:
        st['Step-Postbuckle'] = Step([Frame(0.0, co0), Frame(0.5, coF), Frame(1.0, coF)])
    return Odb(st, {'MEMBRANE-1': inst})


def install_fake_abaqus():
    odb_mod = types.ModuleType('odbAccess')
    odb_mod.openOdb = lambda path, readOnly=True: _CUR['odb']
    abc = types.ModuleType('abaqusConstants')
    for _k in ('CENTROID', 'MAX_INPLANE_PRINCIPAL', 'MIN_INPLANE_PRINCIPAL',
               'INTEGRATION_POINT', 'NODAL'):
        setattr(abc, _k, _k)
    sys.modules['odbAccess'] = odb_mod
    sys.modules['abaqusConstants'] = abc


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_main(mod, odb, mode='HFONLY', tmpd=None):
    """main() 실제 실행. extraction.txt 는 tmpd 로 격리."""
    _CUR['odb'] = odb
    odb_path = os.path.join(tmpd, 'HF_Postbuckle.odb')
    if not os.path.exists(odb_path):
        with io.open(odb_path, 'wb') as fh:
            fh.write(b'x' * 64)
    argv = ['eval_abaqus.py', odb_path, 'HFONLY'] if mode == 'HFONLY' \
        else ['eval_abaqus.py', odb_path, mode]
    cwd0 = os.getcwd()
    argv0 = list(sys.argv)
    sys.argv = list(argv)          # eval_abaqus.main() 은 sys.argv 를 직접 읽는다
    os.chdir(tmpd)
    buf = io.StringIO()
    rc = None
    try:
        with contextlib.redirect_stdout(buf):
            rc = mod.main()
        rc = 0 if rc is None else rc
    except SystemExit as e:
        rc = e.code if isinstance(e.code, int) else 1
    finally:
        os.chdir(cwd0)
        sys.argv = argv0
    return rc, buf.getvalue()


def lines_of(out, prefix):
    return [l.strip() for l in out.splitlines() if l.strip().startswith(prefix)]


def parse_dn(out):
    ls = lines_of(out, 'DN:')
    if not ls:
        return None
    d = {}
    body = ls[0][3:]
    for tok in body.split():
        if '=' in tok:
            k, v = tok.split('=', 1)
            d[k] = v
    return d


# ------------------------------------------------------------------ 앵커 형상
ALPHA = 5.0
TAN_A = math.tan(math.radians(ALPHA))


def shape_flat(x, y):
    return 0.0


def shape_tilt_x(x, y):
    """z = x*tan(alpha): 법선 = (-sin a, 0, cos a) -> dn = alpha, psi = 180deg."""
    return x * TAN_A


def shape_tilt_45(x, y):
    """z = tan(a)*(x+y)/sqrt2: 법선 방위 = 225deg(=atan2(-,-) -> -135deg)."""
    return TAN_A * (x + y) / math.sqrt(2.0)


def shape_sine(x, y):
    """평균경사 0 인 대칭 사인파(한 주기, x 방향). 면적가중 합에서 상쇄 -> dn = 0."""
    return 0.02 * math.sin(2.0 * math.pi * x / 4.0)


def shape_cone(x, y):
    """비대칭 원뿔형 처짐(정점이 +x 로 치우침). 국소 경사는 있으나 일관 기울기는 작다."""
    h, R, x0 = 0.02, 2.0, 0.3
    r = math.sqrt((x - x0) ** 2 + y ** 2)
    return -h * max(0.0, 1.0 - r / R)


def main():
    install_fake_abaqus()
    tmpd = tempfile.mkdtemp(prefix='dn_stub_')
    mod = load_module(os.path.join(CODE, 'eval_abaqus.py'), 'eval_abaqus_new')

    print('--- (a) 완전 평판: dn == 0 ---')
    rc, out = run_main(mod, build_odb(shape=shape_flat), tmpd=tmpd)
    dn = parse_dn(out)
    check('(a) rc=0 + DN 줄 1개', rc == 0 and dn is not None and len(lines_of(out, 'DN:')) == 1,
          (rc, lines_of(out, 'DN:')))
    check('(a) dn_deg == 0 (rel < 1e-9)', dn is not None and relerr(float(dn['dn_deg']), 0.0) < 1e-9,
          dn and dn.get('dn_deg'))
    check('(a) n_eq = (0,0,1) + n_elem=64',
          dn is not None and abs(float(dn['n_eq_x'])) < 1e-12 and abs(float(dn['n_eq_y'])) < 1e-12
          and relerr(float(dn['n_eq_z']), 1.0) < 1e-12 and int(dn['n_elem']) == 64, dn)
    if dn:
        print('     DN> %s' % lines_of(out, 'DN:')[0])

    print('--- (b) 균일 기울기 alpha=5deg: dn == 5 ---')
    rc, out = run_main(mod, build_odb(shape=shape_tilt_x), tmpd=tmpd)
    dn = parse_dn(out)
    check('(b) dn_deg == 5.0 (rel < 1e-9)',
          dn is not None and relerr(float(dn['dn_deg']), ALPHA) < 1e-9, dn and dn.get('dn_deg'))
    check('(b) dn_psi_deg == 180 (기울기 방위)',
          dn is not None and relerr(float(dn['dn_psi_deg']), 180.0) < 1e-9,
          dn and dn.get('dn_psi_deg'))
    check('(b) n_eq = (-sin5, 0, cos5)',
          dn is not None and relerr(float(dn['n_eq_x']), -math.sin(math.radians(ALPHA))) < 1e-9
          and abs(float(dn['n_eq_y'])) < 1e-12
          and relerr(float(dn['n_eq_z']), math.cos(math.radians(ALPHA))) < 1e-9, dn)
    if dn:
        print('     DN> %s' % lines_of(out, 'DN:')[0])

    rc, out = run_main(mod, build_odb(shape=shape_tilt_45), tmpd=tmpd)
    dn45 = parse_dn(out)
    check('(b*) 45deg 방위: dn == 5, dn_psi == -135',
          dn45 is not None and relerr(float(dn45['dn_deg']), ALPHA) < 1e-9
          and relerr(float(dn45['dn_psi_deg']), -135.0) < 1e-9, dn45)

    rc, out = run_main(mod, build_odb(shape=shape_tilt_x, folded=(3, 11)), tmpd=tmpd)
    dnf = parse_dn(out)
    check('(b**) 뒤집힌 요소 2개(connectivity 역순) -> dn 여전히 5 (n_z<0 정렬 규약)',
          dnf is not None and relerr(float(dnf['dn_deg']), ALPHA) < 1e-9, dnf)

    print('--- (c) 대칭 사인파 주름(평균경사 0): dn == 0, hf > 0 ---')
    rc, out = run_main(mod, build_odb(shape=shape_sine), tmpd=tmpd)
    dn = parse_dn(out)
    res = lines_of(out, 'RESULTS:')
    check('(c) DN 줄 존재 + dn_deg == 0 (rel < 1e-9)',
          dn is not None and relerr(float(dn['dn_deg']), 0.0) < 1e-9, dn and dn.get('dn_deg'))
    check('(c) 주름은 사라지지 않는다: n_elem=64 이고 area_ref > 평면적',
          dn is not None and int(dn['n_elem']) == 64
          and float(dn['area_ref']) > 16.0 * 1.0001, dn and dn.get('area_ref'))
    _hf = None
    try:
        _hf = float(res[0].split(',')[3])
    except Exception:
        pass
    check('(c) 같은 형상에서 hf > 0 (delta_n 과 hf 는 독립 정보)',
          _hf is not None and _hf > 0.0, (_hf, res))
    if dn:
        print('     DN> %s' % lines_of(out, 'DN:')[0])
        print('     RESULTS> %s' % res[0])

    print('--- (d) 비대칭 원뿔형 처짐: 0 < dn < alpha_max ---')
    rc, out = run_main(mod, build_odb(shape=shape_cone), tmpd=tmpd)
    dn = parse_dn(out)
    amax = math.degrees(math.atan(0.02 / 2.0))
    dnv = float(dn['dn_deg']) if dn else -1.0
    check('(d) 0 < dn_deg < alpha_max=%.6f' % amax, 0.0 < dnv < amax, (dnv, amax))
    check('(d) dn_psi 유한 + n_eq_z < 1', dn is not None and abs(float(dn['dn_psi_deg'])) <= 180.0
          and float(dn['n_eq_z']) < 1.0, dn)
    if dn:
        print('     DN> %s' % lines_of(out, 'DN:')[0])

    print('--- (e) DN:FAIL 경로: hf 는 성공, DN 만 실패 ---')
    #   Postbuckle 자체가 없으면 get_hf 가 먼저 실패한다(기존 동작). DN 만 떼어내 실패시키려면
    #   hf 는 성공하는 ODB 가 필요하다 -> 5절점 connectivity(hf=1e-12, 요소 0개로 DN 불가).
    rc, out = run_main(mod, build_odb(conn_pad=True), tmpd=tempfile.mkdtemp(dir=tmpd))
    dnl = lines_of(out, 'DN:')
    check('(e) DN:FAIL reason=... 한 줄', len(dnl) == 1 and dnl[0].startswith('DN:FAIL reason='), dnl)
    check('(e) DN:FAIL 에 공백/비ASCII 없음(한 토큰)',
          bool(dnl) and len(dnl[0].split()) == 2 and all(32 <= ord(c) < 127 for c in dnl[0]), dnl)
    check('(e) RESULTS 줄은 평소대로 (nan,nan,nan,hf)',
          len(lines_of(out, 'RESULTS:')) == 1
          and lines_of(out, 'RESULTS:')[0].startswith('RESULTS:nan,nan,nan,'),
          lines_of(out, 'RESULTS:'))
    _e = os.path.join(tmpd, 'extraction.txt')
    check('(e) extraction.txt 도 평소대로 4필드 기록',
          os.path.exists(_e) and len(io.open(_e, encoding='utf-8').read().split(',')) == 4,
          io.open(_e, encoding='utf-8').read() if os.path.exists(_e) else None)
    if dnl:
        print('     DN> %s' % dnl[0])
    #   함수 단위: Postbuckle 스텝 자체가 없을 때도 예외 없이 사유를 돌려준다.
    _CUR['odb'] = build_odb(with_pb=False)
    _dn_nopb = mod.get_dn(os.path.join(tmpd, 'HF_Postbuckle.odb'))
    check('(e*) get_dn(Postbuckle 없음) -> ok=False reason=step_postbuckle_missing',
          _dn_nopb.get('ok') is False and _dn_nopb.get('reason') == 'step_postbuckle_missing',
          _dn_nopb)
    check('(e*) get_dn(파일 없음) -> ok=False reason=file_not_found',
          mod.get_dn(os.path.join(tmpd, 'no_such.odb')).get('reason') == 'file_not_found')

    check('(e*) DN 출력이 cp949 로 인코딩 가능(콘솔 보호)', True)
    for _i, _l in enumerate(dnl):
        try:
            _l.encode('cp949')
        except Exception as _ex:
            check('(e*) DN 줄 cp949 인코딩', False, _ex)
            break
    check('(e*) DN 줄에 금지문자(em/en dash 등) 0', not FORBIDDEN_RE.search(''.join(_dn_out(mod, tmpd))))

    print('--- (f) RESULTS 회귀: git HEAD 판과 바이트 동일 ---')
    old_src = subprocess.check_output(['git', 'show', 'HEAD:code/eval_abaqus.py'],
                                      cwd=os.path.dirname(CODE)).decode('utf-8')
    old_path = os.path.join(tmpd, 'eval_abaqus_old.py')
    with io.open(old_path, 'w', encoding='utf-8') as fh:
        fh.write(old_src)
    old = load_module(old_path, 'eval_abaqus_old')
    for _nm, _shape in (('flat', shape_flat), ('tilt', shape_tilt_x), ('sine', shape_sine),
                        ('cone', shape_cone)):
        _d1 = tempfile.mkdtemp(prefix='dn_old_', dir=tmpd)
        _d2 = tempfile.mkdtemp(prefix='dn_new_', dir=tmpd)
        rc_o, out_o = run_main(old, build_odb(shape=_shape), tmpd=_d1)
        rc_n, out_n = run_main(mod, build_odb(shape=_shape), tmpd=_d2)
        r_o, r_n = lines_of(out_o, 'RESULTS:'), lines_of(out_n, 'RESULTS:')
        ok = (len(r_o) == 1 and len(r_n) == 1 and r_o[0].encode() == r_n[0].encode())
        check('(f) %s: RESULTS 줄 바이트 동일' % _nm, ok, (r_o, r_n))
        e_o = io.open(os.path.join(_d1, 'extraction.txt'), 'rb').read()
        e_n = io.open(os.path.join(_d2, 'extraction.txt'), 'rb').read()
        check('(f) %s: extraction.txt 바이트 동일' % _nm, e_o == e_n, (e_o, e_n))
    check('(f) 새 판에는 DN 줄이 추가로 있다(기존 줄은 그대로)',
          any(l.startswith('DN:dn_deg=') for l in lines_of(out_n, 'DN:')))

    print('--- (g) run_abaqus 배선: _archive_dn 실제 실행 ---')
    rasrc_path = os.path.join(CODE, 'run_abaqus.py')
    rasrc = io.open(rasrc_path, encoding='utf-8').read()
    m = re.search(r'\ndef _archive_dn\(.*?\n(?=def _r13_ct_keys\()', rasrc, re.S)
    check('(g) run_abaqus.py 에서 _archive_dn 소스 추출', bool(m))
    i_res = rasrc.find("_archive_result(x_c, d_c, 'HF', _rc, _extra)")
    i_dn = rasrc.find('_archive_dn(x_c, d_c, _dn_lines)')
    i_slp = rasrc.find('_archive_slope(x_c, d_c)', i_dn)   # def 가 아니라 **호출** 위치
    check('(g) 소스 순서: _archive_result < _archive_dn < _archive_slope',
          -1 < i_res < i_dn < i_slp, (i_res, i_dn, i_slp))
    i_collect = rasrc.find('if _ls.startswith("DN:")')
    i_failchk = rasrc.find('_eval_bad = True')
    check('(g) DN 수집이 FAIL 검사보다 **앞** (DN:FAIL 이 런을 뒤집지 않는다)',
          -1 < i_collect < i_failchk, (i_collect, i_failchk))

    _cols_exp = ['x_c', 'd_c', 'fidelity', 'timestamp', 'stabilize_used',
                 'allds_allie_last_pct', 'allds_allie_max_pct', 'energy_verdict',
                 'backward_area_pct', 'backward_elems', 'total_elems',
                 'hf_flip', 'hf_true_normal', 'delta_pct',
                 'lf1', 'lf2', 'lf3', 'hf']
    mc = re.search(r'_cols = \[(.*?)\]', rasrc, re.S)
    _cols_now = re.findall(r"'([^']+)'", mc.group(1)) if mc else []
    check('(g) all_results.csv _cols 블록 불변(18열, 순서까지)',
          _cols_now == _cols_exp, _cols_now)

    if m:
        cdir = tempfile.mkdtemp(prefix='dn_wire_', dir=tmpd)
        rdir = os.path.join(cdir, 'results')
        os.makedirs(rdir)
        from aba_imperfection import value_tag
        ptxt = os.path.join(rdir, 'results_xc0p5_dc0p3.txt')
        base = ('x_c=0.5\nd_c=0.3\nfidelity=HF\nhf=3.9e-05\n'
                'slope_f1_hf_pred=4.2e-05\n')
        io.open(ptxt, 'w', encoding='utf-8').write(base)
        csvp = os.path.join(rdir, 'all_results.csv')
        csv0 = ','.join(_cols_exp) + '\nhdr\n' + '0.5,0.3,HF,x,1,2,3,PASS,0,0,64,1,1,0,1e-3,2e-3,3e-3,4e-5\n'
        io.open(csvp, 'w', encoding='utf-8').write(csv0)

        logs = []

        def _emit(*a):
            logs.append(' '.join(str(x) for x in a))

        ns = {'os': os, 'io': io, 'value_tag': value_tag, 'emit': _emit, '_HERE': cdir}
        exec(compile(m.group(0), '<run_abaqus:_archive_dn>', 'exec'), ns)
        fn = ns['_archive_dn']

        # (g1) 정상 DN 줄
        _d1 = tempfile.mkdtemp(prefix='dn_good_', dir=tmpd)
        with contextlib.redirect_stdout(io.StringIO()):
            rc, o1 = run_main(mod, build_odb(shape=shape_tilt_x), tmpd=_d1)
        dnline = lines_of(o1, 'DN:')[0]
        with contextlib.redirect_stdout(io.StringIO()):
            fn(0.5, 0.3, [dnline])
        txt = io.open(ptxt, encoding='utf-8').read()
        kv = dict(l.split('=', 1) for l in txt.splitlines() if '=' in l)
        check('(g) 점별 txt 에 dn_deg/dn_psi_deg/dn_n_eq_*/dn_n_elem/dn_area_ref append',
              all(k in kv for k in ('dn_deg', 'dn_psi_deg', 'dn_n_eq_x', 'dn_n_eq_y',
                                    'dn_n_eq_z', 'dn_n_elem', 'dn_area_ref')), sorted(kv))
        check('(g) 기존 key(x_c/d_c/fidelity/hf/slope_*) 보존',
              all(k in kv for k in ('x_c', 'd_c', 'fidelity', 'hf', 'slope_f1_hf_pred')), sorted(kv))
        check('(g) hn 값이 DN 줄과 일치(5deg)',
              relerr(float(kv['dn_deg']), ALPHA) < 1e-9, kv['dn_deg'])
        check('(g) all_results.csv **바이트 불변**',
              io.open(csvp, encoding='utf-8').read() == csv0)
        check('(g) [DN] 로그 1줄', any(l.startswith('[DN]') and 'dn_deg' in l for l in logs), logs)

        # (g2) DN 줄 없음 -> 경고(조용히 넘기지 않음)
        logs[:] = []
        _before = io.open(ptxt, encoding='utf-8').read()
        with contextlib.redirect_stdout(io.StringIO()):
            fn(0.5, 0.3, [])
        check('(g) DN 줄 없음 -> [DN] !!! 경고 1줄', any(l.startswith('[DN] !!!') and '줄 없음' in l
                                                       for l in logs), logs)
        check('(g) DN 줄 없음 -> txt 는 그대로', io.open(ptxt, encoding='utf-8').read() == _before)

        # (g3) DN:FAIL -> dn_deg=nan + 사유, 0 으로 뭉개지 않음
        logs[:] = []
        with contextlib.redirect_stdout(io.StringIO()):
            fn(0.5, 0.3, ['DN:FAIL reason=step_postbuckle_missing'])
        txt3 = io.open(ptxt, encoding='utf-8').read()
        check('(g) DN:FAIL -> dn_deg=nan + dn_note 사유',
              'dn_deg=nan' in txt3 and 'step_postbuckle_missing' in txt3, txt3.splitlines()[-3:])
        check('(g) DN:FAIL 이후에도 CSV 불변',
              io.open(csvp, encoding='utf-8').read() == csv0)

    print('--- (h) 금지문자/인코딩 ---')
    _new_lines = [ln for ln in io.open(os.path.join(CODE, 'eval_abaqus.py'),
                                       encoding='utf-8').read().splitlines()
                  if any(k in ln for k in ('19.106', 'DN:', 'dn_', 'delta_n', 'n_eq'))]
    check('(h) eval_abaqus.py DN 관련 줄에 금지문자 0',
          not any(FORBIDDEN_RE.search(ln) for ln in _new_lines))
    _new_lines2 = [ln for ln in rasrc.splitlines()
                   if any(k in ln for k in ('19.106', '[DN]', '_archive_dn', 'DN:'))]
    check('(h) run_abaqus.py DN 관련 줄에 금지문자 0',
          not any(FORBIDDEN_RE.search(ln) for ln in _new_lines2))
    _probe = ('DN:dn_deg=%.10e dn_psi_deg=%.10e n_eq_x=%.10e n_eq_y=%.10e '
              'n_eq_z=%.10e n_elem=%d area_ref=%.10e' % (5.0, 180.0, -0.087, 0.0, 0.996, 64, 16.0))
    try:
        _probe.encode('cp949')
        _ok = True
    except Exception:
        _ok = False
    check('(h) DN 줄 형식이 cp949 로 인코딩 가능', _ok, _probe)

    print('')
    print('  %d checks, %d FAIL' % (N[0], len(FAIL)))
    if FAIL:
        for f in FAIL:
            print('   FAILED: %s' % f)
        return 1
    print('  ALL OK')
    return 0


def _dn_out(mod, tmpd):
    """(e*) 용: 방금 실행한 DN 출력 재수집(간단히 실패 경로 한 번 더)."""
    try:
        rc, out = run_main(mod, build_odb(with_pb=False), tmpd=tempfile.mkdtemp(dir=tmpd))
        return lines_of(out, 'DN:')
    except Exception:
        return []


if __name__ == '__main__':
    sys.exit(main())
