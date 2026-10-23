#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eval_lf_stub.py - Abaqus 없이 eval_abaqus.get_lf/get_zk/main 을 검증한다.

가짜 ODB 를 만들어 다음을 확인한다(라이선스 0):

  1) paper path (Step-HighTension 없음)
       - get_lf 가 Step-GlobalTension 을 선택한다.
       - lf1/lf2/lf3 계산값.
       - ZK (state_wrinkled_frac, intensity mean/p95, global direction).
       - zk_Iw (reference weighting, slack, positive-part) sidecar.
       - RESULTS: 3필드 / extraction.txt / zk_metrics.txt.

  2) old disp path (Step-HighTension 있음)
       - get_lf 가 Step-HighTension 을 선택한다(구 기준선 보존).
       - HighTension 프레임이 0 이면 RESULTS 도 0 이다.

  3) E 필드 없음
       - ZK:FAIL reason=E_field_missing.
       - 그래도 RESULTS 3필드 / extraction.txt 는 정상이다.

  4) I_w positive-part clip
       - E1>0, E2<0 이지만 |E2|-nu*E1 <= 0 인 요소는 I_w 에서 0.
       - slack 항만 남는지, clip fraction / share 가 맞는지 확인.

  5) LE fallback
       - Abaqus NLGEOM=ON 처럼 E 대신 LE 만 있는 ODB.
       - LE 주성분을 E=0.5(exp(2*LE)-1) 로 변환해 같은 zk_Iw/ZK 가 나오는지 확인.

사용:
    python3 eval_lf_stub.py
종료코드: 0 PASS / 1 FAIL
"""
import contextlib
import io
import math
import os
import shutil
import sys
import tempfile
import types

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

# ---------------------------------------------------------------------------
# fake abaqusConstants / odbAccess: eval_abaqus import **전에** 주입해야 한다.
# ---------------------------------------------------------------------------
_CONST = types.ModuleType('abaqusConstants')
_CONST.CENTROID = 'CENTROID'
_CONST.INTEGRATION_POINT = 'INTEGRATION_POINT'
_CONST.MAX_INPLANE_PRINCIPAL = 'MAX_INPLANE_PRINCIPAL'
_CONST.MIN_INPLANE_PRINCIPAL = 'MIN_INPLANE_PRINCIPAL'
sys.modules['abaqusConstants'] = _CONST

_CURRENT = {'odb': None}


def _fake_openOdb(path, readOnly=True):
    return _CURRENT['odb']


_ODB_MOD = types.ModuleType('odbAccess')
_ODB_MOD.openOdb = _fake_openOdb
sys.modules['odbAccess'] = _ODB_MOD

import eval_abaqus as E   # noqa: E402  (fake modules must be installed first)


# ---------------------------------------------------------------------------
# fake ODB objects (only the API surface eval_abaqus actually touches)
# ---------------------------------------------------------------------------
class _Val(object):
    def __init__(self, elementLabel=None, data=None, nodeLabel=None):
        self.elementLabel = elementLabel
        self.data = data
        self.nodeLabel = nodeLabel


class _Field(object):
    def __init__(self, values):
        self.values = list(values)

    def getSubset(self, region=None, position=None):
        return self

    def getScalarField(self, invariant=None):
        out = []
        for v in self.values:
            d = v.data
            s11, s22, s12 = float(d[0]), float(d[1]), float(d[2])
            m = 0.5 * (s11 + s22)
            r = (0.25 * (s11 - s22) ** 2 + s12 * s12) ** 0.5
            val = (m + r) if invariant == 'MAX_INPLANE_PRINCIPAL' else (m - r)
            out.append(_Val(v.elementLabel, val))
        return _Field(out)


class _Frame(object):
    def __init__(self, fields, frameValue=1.0):
        self.fieldOutputs = fields
        self.frameValue = frameValue


class _Step(object):
    def __init__(self, frames):
        self.frames = list(frames)


class _Node(object):
    def __init__(self, label, coordinates):
        self.label = label
        self.coordinates = coordinates


class _Elem(object):
    def __init__(self, label, connectivity):
        self.label = label
        self.connectivity = connectivity


class _Instance(object):
    def __init__(self, name, nodes, elements):
        self.name = name
        self.nodes = list(nodes)
        self.elements = list(elements)


class _Root(object):
    def __init__(self, instance):
        self.instances = {instance.name: instance}


class _Odb(object):
    def __init__(self, steps, instance):
        self.steps = steps
        self.rootAssembly = _Root(instance)

    def close(self):
        pass


# ---------------------------------------------------------------------------
# synthetic model
#   elem 1: taut       S=(100, 50)   E=(1e-3,  5e-4)  vol=1  local1=+x
#   elem 2: wrinkled   S=(100,-50)   E=(3e-4, -2e-4)  vol=1  local1=30 deg
#   elem 3: slack      S=(-100,-50)  E=(-1e-4,-3e-4)  vol=2  local1=+x
#
#   reference EVOL (frameValue=0): elem1=1, elem2=2, elem3=2  (A_ref=5)
#   (public ZK 는 현재 EVOL 1,1,2 를 그대로 쓴다 -> state_wrinkled_frac=1/4)
#
# expected lf1 = (1+2)/4 = 0.75
# expected lf2 = -(0.34/2.5e9)*100*(-50) = 6.8e-7   (elem 2 only)
# expected lf3 = mean(50/100, 50/100) = 0.5          (elems 1,2)
# expected ZK  : wrinkled frac = 1/4 = 0.25
#                intensity mean = p95 = 2e-4
#                dir = 30 deg (elem 2 local1 rotated by 30 deg)
# expected I_w = [2*max(2e-4 - 0.34*3e-4, 0)
#                 + 2*(3e-4 + 0.34*1e-4)] / 5 = 1.728e-4
#                I_w_p95 = 3.34e-4
#                slack share = 6.68e-4 / 8.64e-4 = 0.773148...
#                wr share    = 1.96e-4 / 8.64e-4 = 0.226852...
#                clip zero frac = 0
# ---------------------------------------------------------------------------
ELEMENTS = [
    (1, (1, 2, 3, 4), (100.0, 50.0, 0.0),
     (1e-3, 5e-4, 0.0, 0.0, 0.0, 0.0), 1.0),
    (2, (20, 21, 22, 23), (100.0, -50.0, 0.0),
     (3e-4, -2e-4, 0.0, 0.0, 0.0, 0.0), 1.0),
    (3, (30, 31, 32, 33), (-100.0, -50.0, 0.0),
     (-1e-4, -3e-4, 0.0, 0.0, 0.0, 0.0), 2.0),
]

#   clip test 는 같은 모델에서 elem 2 만 E1=1e-3 으로 바꾼다 -> |E2|-nu*E1 < 0.
REF_VOL = {1: 1.0, 2: 2.0, 3: 2.0}
CLIP_ELEMENTS = [
    (1, ELEMENTS[0][1], ELEMENTS[0][2], ELEMENTS[0][3], ELEMENTS[0][4]),
    (2, ELEMENTS[1][1], ELEMENTS[1][2],
     (1e-3, -2e-4, 0.0, 0.0, 0.0, 0.0), ELEMENTS[1][4]),
    (3, ELEMENTS[2][1], ELEMENTS[2][2], ELEMENTS[2][3], ELEMENTS[2][4]),
]


def _green_to_log_le(e11, e22):
    """대각 Green-Lagrange E1/E2 -> logarithmic LE1/LE2."""
    m = 0.5 * (e11 + e22)
    r = (0.25 * (e11 - e22) ** 2) ** 0.5
    e1 = m + r
    e2 = m - r
    return 0.5 * math.log1p(2.0 * e1), 0.5 * math.log1p(2.0 * e2)


def _le_elements(elements):
    """같은 E 상태를 LE 필드로 표현한 합성 요소 목록(테스트용)."""
    out = []
    for (lbl, conn, s, e, vol) in elements:
        le1, le2 = _green_to_log_le(e[0], e[1])
        out.append((lbl, conn, s, (le1, le2, 0.0, 0.0, 0.0, 0.0), vol))
    return out


LE_ELEMENTS = _le_elements(ELEMENTS)

NODES = [
    (1, (0.0, 0.0, 0.0)), (2, (1.0, 0.0, 0.0)),
    (3, (1.0, 1.0, 0.0)), (4, (0.0, 1.0, 0.0)),
    # elem 2: p2-p1 = (cos30, sin30)
    (20, (10.0, 0.0, 0.0)),
    (21, (10.866025403784439, 0.5, 0.0)),
    (22, (10.366025403784439, 1.366025403784439, 0.0)),
    (23, (9.5, 0.866025403784439, 0.0)),
    (30, (20.0, 0.0, 0.0)), (31, (21.0, 0.0, 0.0)),
    (32, (21.0, 1.0, 0.0)), (33, (20.0, 1.0, 0.0)),
]


def _frame(elements, with_E=True, strain_key='E'):
    fields = {
        'S': _Field([_Val(lbl, s) for (lbl, conn, s, e, vol) in elements]),
        'EVOL': _Field([_Val(lbl, vol) for (lbl, conn, s, e, vol) in elements]),
    }
    if with_E:
        fields[strain_key] = _Field([_Val(lbl, e)
                                     for (lbl, conn, s, e, vol) in elements])
    return _Frame(fields)


def _ref_frame(elements):
    """frameValue=0 의 기준 EVOL 만 가진 프레임(I_w 전용)."""
    return _Frame(
        {'EVOL': _Field([_Val(lbl, REF_VOL.get(lbl, vol))
                         for (lbl, conn, s, e, vol) in elements])},
        frameValue=0.0)


def _make_odb(include_high=False, with_E=True, high_zero=True, elements=None,
              strain_key='E'):
    elements = elements or ELEMENTS
    instance = _Instance(
        'MEMBRANE-1',
        [_Node(lbl, c) for (lbl, c) in NODES],
        [_Elem(lbl, conn) for (lbl, conn, s, e, vol) in elements],
    )
    steps = {'Step-GlobalTension': _Step([_ref_frame(elements),
                                          _frame(elements, with_E=with_E,
                                                 strain_key=strain_key)])}
    if include_high:
        if high_zero:
            elems = [(lbl, conn, (0.0, 0.0, 0.0),
                      (0.0, 0.0, 0.0, 0.0, 0.0, 0.0), vol)
                     for (lbl, conn, s, e, vol) in elements]
        else:
            elems = elements
        steps['Step-HighTension'] = _Step([_ref_frame(elems),
                                           _frame(elems, with_E=with_E,
                                                  strain_key=strain_key)])
    return _Odb(steps, instance)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _run_main(odb, odb_path, workdir, mode='LF'):
    _CURRENT['odb'] = odb
    old_argv = list(sys.argv)
    old_cwd = os.getcwd()
    sys.argv = ['eval_abaqus.py', odb_path, mode]
    os.chdir(workdir)
    buf = io.StringIO()
    rc = 0
    try:
        with contextlib.redirect_stdout(buf):
            try:
                E.main()
            except SystemExit as ex:
                rc = ex.code if ex.code is not None else 0
    finally:
        sys.argv = old_argv
        os.chdir(old_cwd)
    return rc, buf.getvalue()


def _results(out):
    for line in out.splitlines():
        if line.startswith('RESULTS:'):
            body = line.split(':', 1)[1].strip()
            if body.upper() == 'FAIL':
                return None
            try:
                return [float(x) for x in body.split(',')]
            except ValueError:
                return None
    return None


def _zk_line(out):
    for line in out.splitlines():
        if line.startswith('ZK:'):
            return line
    return None


def _kv_line(line):
    out = {}
    for tok in line.split(':', 1)[1].split():
        if '=' in tok:
            k, v = tok.split('=', 1)
            out[k] = v
    return out


def _read(path):
    with open(path, 'r') as f:
        return f.read()


def _read_kv(path):
    out = {}
    for line in _read(path).splitlines():
        s = line.strip()
        if not s or '=' not in s:
            continue
        k, v = s.split('=', 1)
        out[k] = v
    return out


def _num(kv, key):
    try:
        return float(kv.get(key, 'nan'))
    except (TypeError, ValueError):
        return float('nan')


def _close(actual, expected, atol=1e-15, rtol=0.0):
    try:
        return abs(float(actual) - float(expected)) <= (
            atol + rtol * abs(float(expected)))
    except (TypeError, ValueError):
        return False


def _check(name, cond, detail=''):
    if detail is None:
        detail = ''
    print("  %-4s %s%s" % ('OK' if cond else 'FAIL', name,
                           (' - ' + str(detail)) if detail else ''))
    return bool(cond)


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------
def main():
    ok = True
    root = tempfile.mkdtemp(prefix='eval_lf_stub_')
    try:
        # ---------------- 1) paper path ----------------
        d1 = os.path.join(root, 'paper')
        os.makedirs(d1)
        p1 = os.path.join(d1, 'LF_paper.odb')
        with open(p1, 'wb') as f:
            f.write(b'')
        rc, out = _run_main(_make_odb(include_high=False, with_E=True), p1, d1)
        res = _results(out)
        zk = _zk_line(out)

        ok &= _check('paper: exit 0', rc == 0, 'rc=%r' % rc)
        ok &= _check('paper: RESULTS has 3 fields',
                     res is not None and len(res) == 3, repr(res))
        if res is not None and len(res) == 3:
            ok &= _check('paper: lf1 == 0.75', abs(res[0] - 0.75) < 1e-12,
                         '%.12g' % res[0])
            ok &= _check('paper: lf2 == 6.8e-7', abs(res[1] - 6.8e-7) < 1e-15,
                         '%.12g' % res[1])
            ok &= _check('paper: lf3 == 0.5', abs(res[2] - 0.5) < 1e-12,
                         '%.12g' % res[2])
        ok &= _check('paper: ZK line exists',
                     zk is not None and zk.startswith('ZK:'), repr(zk))
        ok &= _check('paper: ZK is OK (not FAIL)',
                     zk is not None and zk.startswith('ZK:') and 'FAIL' not in zk, repr(zk))
        if zk is not None and zk.startswith('ZK:') and 'FAIL' not in zk:
            kv = _kv_line(zk)
            ok &= _check('paper: ZK wrinkled frac == 0.25',
                         abs(float(kv.get('state_wrinkled_frac', 'nan')) - 0.25) < 1e-12,
                         kv.get('state_wrinkled_frac'))
            ok &= _check('paper: ZK intensity mean == 2e-4',
                         abs(float(kv.get('zk_intensity_mean', 'nan')) - 2e-4) < 1e-15,
                         kv.get('zk_intensity_mean'))
            ok &= _check('paper: ZK intensity p95 == 2e-4',
                         abs(float(kv.get('zk_intensity_p95', 'nan')) - 2e-4) < 1e-15,
                         kv.get('zk_intensity_p95'))
            ok &= _check('paper: ZK dir == 30 deg',
                         abs(float(kv.get('zk_dir_deg', 'nan')) - 30.0) < 1e-9,
                         kv.get('zk_dir_deg'))
        ext = os.path.join(d1, 'extraction.txt')
        ok &= _check('paper: extraction.txt exists', os.path.exists(ext))
        if os.path.exists(ext):
            parts = [p for p in _read(ext).strip().split(',') if p != '']
            ok &= _check('paper: extraction.txt has 3 fields', len(parts) == 3, repr(parts))
        zkm = os.path.join(d1, 'zk_metrics.txt')
        ok &= _check('paper: zk_metrics.txt exists', os.path.exists(zkm))
        if os.path.exists(zkm):
            txt = _read(zkm)
            ok &= _check('paper: zk status=OK', 'status=OK' in txt,
                         txt.replace('\n', '; '))
            ok &= _check('paper: zk step=Step-GlobalTension',
                         'step=Step-GlobalTension' in txt)
            ok &= _check('paper: zk file wrinkled frac=0.25',
                         'state_wrinkled_frac=0.25' in txt)
            zkv = _read_kv(zkm)
            ok &= _check('paper: public zk_area_total == 4 (current EVOL)',
                         _close(_num(zkv, 'zk_area_total'), 4.0),
                         zkv.get('zk_area_total'))
            ok &= _check('paper: zk_Iw == 1.728e-4',
                         _close(_num(zkv, 'zk_Iw'), 1.728e-4),
                         zkv.get('zk_Iw'))
            ok &= _check('paper: zk_Iw_p95 == 3.34e-4',
                         _close(_num(zkv, 'zk_Iw_p95'), 3.34e-4),
                         zkv.get('zk_Iw_p95'))
            ok &= _check('paper: zk_Iw reference area == 5',
                         _close(_num(zkv, 'zk_Iw_area_total'), 5.0),
                         zkv.get('zk_Iw_area_total'))
            ok &= _check('paper: zk_Iw weights == ref_evol',
                         zkv.get('zk_Iw_weights') == 'ref_evol',
                         zkv.get('zk_Iw_weights'))
            ok &= _check('paper: zk_Iw slack share == 0.773148',
                         _close(_num(zkv, 'zk_Iw_slack_share'),
                                0.7731481481481481, atol=1e-12),
                         zkv.get('zk_Iw_slack_share'))
            ok &= _check('paper: zk_Iw wr share == 0.226852',
                         _close(_num(zkv, 'zk_Iw_wr_share'),
                                0.22685185185185186, atol=1e-12),
                         zkv.get('zk_Iw_wr_share'))
            ok &= _check('paper: zk_Iw no clipped wrinkled',
                         _close(_num(zkv, 'zk_Iw_clip_zero_frac'), 0.0),
                         zkv.get('zk_Iw_clip_zero_frac'))
            ok &= _check('paper: zk_Iw ref_missing == 0',
                         _close(_num(zkv, 'zk_Iw_ref_missing'), 0.0),
                         zkv.get('zk_Iw_ref_missing'))
            ok &= _check('paper: strain_field == E',
                         zkv.get('strain_field') == 'E',
                         zkv.get('strain_field'))

        # ---------------- 2) old disp path ----------------
        d2 = os.path.join(root, 'disp')
        os.makedirs(d2)
        p2 = os.path.join(d2, 'LF_disp.odb')
        with open(p2, 'wb') as f:
            f.write(b'')
        rc, out = _run_main(_make_odb(include_high=True, with_E=True, high_zero=True),
                            p2, d2)
        res = _results(out)
        zk = _zk_line(out)

        ok &= _check('disp: exit 0', rc == 0, 'rc=%r' % rc)
        ok &= _check('disp: RESULTS has 3 fields',
                     res is not None and len(res) == 3, repr(res))
        ok &= _check('disp: ZK is OK (not FAIL)',
                     zk is not None and zk.startswith('ZK:') and 'FAIL' not in zk, repr(zk))
        if res is not None and len(res) == 3:
            ok &= _check('disp: HighTension selected (all lf == 0)',
                         max(abs(x) for x in res) < 1e-15, repr(res))
        if zk is not None and zk.startswith('ZK:') and 'FAIL' not in zk:
            kv = _kv_line(zk)
            ok &= _check('disp: ZK wrinkled frac == 0',
                         abs(float(kv.get('state_wrinkled_frac', 'nan'))) < 1e-15,
                         kv.get('state_wrinkled_frac'))
        zkm = os.path.join(d2, 'zk_metrics.txt')
        if os.path.exists(zkm):
            txt = _read(zkm)
            ok &= _check('disp: zk step=Step-HighTension',
                         'step=Step-HighTension' in txt)

        # ---------------- 3) E missing ----------------
        d3 = os.path.join(root, 'noE')
        os.makedirs(d3)
        p3 = os.path.join(d3, 'LF_noE.odb')
        with open(p3, 'wb') as f:
            f.write(b'')
        rc, out = _run_main(_make_odb(include_high=False, with_E=False), p3, d3)
        res = _results(out)
        zk = _zk_line(out)

        ok &= _check('noE: exit 0', rc == 0, 'rc=%r' % rc)
        ok &= _check('noE: RESULTS has 3 fields',
                     res is not None and len(res) == 3, repr(res))
        ok &= _check('noE: ZK:FAIL reason=E_field_missing',
                     zk == 'ZK:FAIL reason=E_field_missing', repr(zk))
        ext = os.path.join(d3, 'extraction.txt')
        ok &= _check('noE: extraction.txt exists', os.path.exists(ext))
        zkm = os.path.join(d3, 'zk_metrics.txt')
        ok &= _check('noE: zk_metrics.txt exists', os.path.exists(zkm))
        if os.path.exists(zkm):
            txt = _read(zkm)
            ok &= _check('noE: zk status=FAIL', 'status=FAIL' in txt,
                         txt.replace('\n', '; '))
            ok &= _check('noE: zk reason=E_field_missing',
                         'reason=E_field_missing' in txt)

        # ---------------- 4) I_w positive-part clip ----------------
        d4 = os.path.join(root, 'clip')
        os.makedirs(d4)
        p4 = os.path.join(d4, 'LF_clip.odb')
        with open(p4, 'wb') as f:
            f.write(b'')
        rc, out = _run_main(_make_odb(elements=CLIP_ELEMENTS), p4, d4)
        res = _results(out)
        ok &= _check('clip: exit 0', rc == 0, 'rc=%r' % rc)
        ok &= _check('clip: RESULTS has 3 fields',
                     res is not None and len(res) == 3, repr(res))
        zkm = os.path.join(d4, 'zk_metrics.txt')
        ok &= _check('clip: zk_metrics.txt exists', os.path.exists(zkm))
        if os.path.exists(zkm):
            zkv = _read_kv(zkm)
            ok &= _check('clip: public zk_intensity_mean == 2e-4 (unchanged)',
                         _close(_num(zkv, 'zk_intensity_mean'), 2e-4),
                         zkv.get('zk_intensity_mean'))
            ok &= _check('clip: zk_Iw == 1.336e-4 (slack only)',
                         _close(_num(zkv, 'zk_Iw'), 1.336e-4),
                         zkv.get('zk_Iw'))
            ok &= _check('clip: zk_Iw_p95 == 3.34e-4',
                         _close(_num(zkv, 'zk_Iw_p95'), 3.34e-4),
                         zkv.get('zk_Iw_p95'))
            ok &= _check('clip: zk_Iw slack share == 1.0',
                         _close(_num(zkv, 'zk_Iw_slack_share'), 1.0),
                         zkv.get('zk_Iw_slack_share'))
            ok &= _check('clip: zk_Iw wr share == 0.0',
                         _close(_num(zkv, 'zk_Iw_wr_share'), 0.0),
                         zkv.get('zk_Iw_wr_share'))
            ok &= _check('clip: zk_Iw clip_zero_frac == 1.0',
                         _close(_num(zkv, 'zk_Iw_clip_zero_frac'), 1.0),
                         zkv.get('zk_Iw_clip_zero_frac'))
            ok &= _check('clip: zk_Iw reference area == 5',
                         _close(_num(zkv, 'zk_Iw_area_total'), 5.0),
                         zkv.get('zk_Iw_area_total'))
            ok &= _check('clip: zk_Iw weights == ref_evol',
                         zkv.get('zk_Iw_weights') == 'ref_evol',
                         zkv.get('zk_Iw_weights'))

        # ---------------- 5) LE fallback (Abaqus NLGEOM=ON) ----------------
        d5 = os.path.join(root, 'le')
        os.makedirs(d5)
        p5 = os.path.join(d5, 'LF_le.odb')
        with open(p5, 'wb') as f:
            f.write(b'')
        rc, out = _run_main(
            _make_odb(elements=LE_ELEMENTS, strain_key='LE'), p5, d5)
        res = _results(out)
        zk = _zk_line(out)
        ok &= _check('LE: exit 0', rc == 0, 'rc=%r' % rc)
        ok &= _check('LE: RESULTS has 3 fields',
                     res is not None and len(res) == 3, repr(res))
        ok &= _check('LE: ZK is OK (not FAIL)',
                     zk is not None and zk.startswith('ZK:') and 'FAIL' not in zk,
                     repr(zk))
        zkm = os.path.join(d5, 'zk_metrics.txt')
        ok &= _check('LE: zk_metrics.txt exists', os.path.exists(zkm))
        if os.path.exists(zkm):
            zkv = _read_kv(zkm)
            ok &= _check('LE: strain_field == LE',
                         zkv.get('strain_field') == 'LE',
                         zkv.get('strain_field'))
            ok &= _check('LE: strain_measure == E',
                         zkv.get('strain_measure') == 'E',
                         zkv.get('strain_measure'))
            ok &= _check('LE: conversion == green_lagrange_from_log',
                         zkv.get('strain_conversion') == 'green_lagrange_from_log',
                         zkv.get('strain_conversion'))
            ok &= _check('LE: zk_Iw == 1.728e-4 (same as E)',
                         _close(_num(zkv, 'zk_Iw'), 1.728e-4),
                         zkv.get('zk_Iw'))
            ok &= _check('LE: zk_Iw_p95 == 3.34e-4',
                         _close(_num(zkv, 'zk_Iw_p95'), 3.34e-4),
                         zkv.get('zk_Iw_p95'))
            ok &= _check('LE: public zk_intensity_mean == 2e-4',
                         _close(_num(zkv, 'zk_intensity_mean'), 2e-4),
                         zkv.get('zk_intensity_mean'))
            ok &= _check('LE: zk_dir == 30 deg',
                         _close(_num(zkv, 'zk_dir_deg'), 30.0, atol=1e-9),
                         zkv.get('zk_dir_deg'))
            ok &= _check('LE: zk_Iw weights == ref_evol',
                         zkv.get('zk_Iw_weights') == 'ref_evol',
                         zkv.get('zk_Iw_weights'))
    finally:
        shutil.rmtree(root, ignore_errors=True)

    print('')
    print('RESULT: %s' % ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
