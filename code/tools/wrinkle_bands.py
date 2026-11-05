#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""wrinkle_bands.py - HF ODB 의 u3 를 파장 대역으로 분해한다 (Plan 1).

    python tools/wrinkle_bands.py --synthetic          # G1 + taper study (라이선스 0)
    abaqus python tools/wrinkle_bands.py --odb <HF.odb> --tag <tag> --out results/bands

정의
    u3(x,y) 참조좌표 기준. 삼각형 경계에서 0으로 가는 smooth taper 를 곱한 뒤
    201x101 사각 격자에서 fft2.  hf_b = K_SLOPE * (1/N^2) * sum_{k in band} k^2 |U_k|^2.
    (Parseval: mean|grad u|^2 = (1/N^2) sum k^2 |U_k|^2. K_SLOPE=1.461716, slope_stats 와 동일.)

대역 (h=0.1 m)
    B1 lambda > 5     B2 2-5     B3 1-2     B4 0.5-1
    B5 0.2-0.5  -> 2-5 elements, **Nyquist 한계: 정량 보고 금지**. 'unresolved' 로만 기록.
    lambda < 0.2 -> aliasing, 보고 금지.

게이트
    G1  합성 단일 대역장(lambda=6): B1 >= 95 %, 인접(B2) <= 2 %.
    G2a 합성장에서 sum_b hf_b / hf_taper = 1.000 +- 0.01 (직교성).
    G2b 실 ODB 에서 hf_taper / hf_odb = 0.85-1.15 (taper/격자/스텐실 오차).

주의
    - u3 는 **참조좌표**에서 읽는다(변형 후 좌표를 쓰면 파장이 왜곡된다).
    - detrend 하지 않는다(B1 포함). hf 는 raw field 기울기에서 나오기 때문이다.
"""
from __future__ import print_function

import argparse
import io
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
#   slope_stats.py 는 import 시 odbAccess 를 요구한다(Abaqus 밖에서 실패). 상수만 복제한다.
#   값/의미는 slope_stats.py:115 와 동일하다.
K_SLOPE = 1.461716

BASE, HEIGHT, H = 20.0, 10.0, 0.1
NX, NY = int(round(BASE / H)) + 1, int(round(HEIGHT / H)) + 1   # 201 x 101
AREA = 0.5 * BASE * HEIGHT
BANDS = [
    ('B1', 5.0, float('inf')),
    ('B2', 2.0, 5.0),
    ('B3', 1.0, 2.0),
    ('B4', 0.5, 1.0),
]
B5 = ('B5', 0.2, 0.5)


def triangle_mask():
    x = np.arange(NX) * H
    y = np.arange(NY) * H
    X, Y = np.meshgrid(x, y)
    inside = (Y <= X + 1e-9) & (Y <= (BASE - X) + 1e-9) & (Y >= 0)
    return X, Y, inside


def taper_weights(X, Y, width):
    """삼각형 세 변까지의 거리로 만든 smoothstep 곱. width=0 이면 mask(0/1)."""
    if width <= 0:
        _, _, m = triangle_mask()
        return m.astype(float)
    d = np.minimum(np.minimum(Y, X - Y), BASE - X - Y)     # 세 변까지 거리
    t = np.clip(d / float(width), 0.0, 1.0)
    w = 0.5 * (1.0 - np.cos(np.pi * t))
    _, _, inside = triangle_mask()
    return np.where(inside, w, 0.0)


def k_grid():
    kx = 2.0 * math.pi * np.fft.fftfreq(NX, d=H)
    ky = 2.0 * math.pi * np.fft.fftfreq(NY, d=H)
    KX, KY = np.meshgrid(kx, ky)
    return np.sqrt(KX ** 2 + KY ** 2)


def decompose(u3, width):
    """u3 (NY,NX) -> dict: bands hf_b, fractions, total, taper ratio inputs."""
    w = taper_weights(*triangle_mask()[:2], width=width)
    ut = np.asarray(u3, dtype=np.float64) * w
    U = np.fft.fft2(ut)
    K = k_grid()
    #   Parseval: (1/N^2) sum_k k^2|U_k|^2 = mean over the **rectangular** grid.
    #   실제 돛은 삼각형(노드 n_in)이므로 사각형 평균을 n_in 기준으로 환산한다.
    _, _, inside = triangle_mask()
    n_in = float(np.count_nonzero(inside))
    scale = float(NX * NY) / n_in
    power = (K ** 2) * (np.abs(U) ** 2) / float(NX * NY) ** 2
    hf_taper = K_SLOPE * float(np.sum(power)) * scale
    out = {'hf_taper': hf_taper, 'bands': {}, 'unresolved': 0.0}
    lam = np.full_like(K, np.inf)
    nz = K > 0
    lam[nz] = 2.0 * math.pi / K[nz]
    covered = np.zeros_like(K, dtype=bool)
    for name, lo, hi in BANDS:
        sel = (K > 0) & (lam > lo) & (lam <= hi)
        out['bands'][name] = K_SLOPE * float(np.sum(power[sel])) * scale
        covered |= sel
    sel5 = (lam > B5[1]) & (lam <= B5[2])
    out['unresolved_b5'] = K_SLOPE * float(np.sum(power[sel5])) * scale
    covered |= sel5
    out['residual'] = K_SLOPE * float(np.sum(power[(K > 0) & ~covered])) * scale
    return out, ut


def synthetic_field(modes):
    """modes: [(lambda_m, amplitude)] -> u3 on triangle (no taper)."""
    X, Y, inside = triangle_mask()
    u = np.zeros_like(X)
    for lam, amp in modes:
        u += amp * np.sin(2.0 * math.pi * X / lam)
    return np.where(inside, u, 0.0)


def synthetic_report(width):
    rows = []
    tests = [
        ('single_B1', [(6.0, 1.0)], 'B1', None),
        ('single_B3', [(1.5, 1.0)], 'B3', None),
        ('multi', [(6.0, 1.0), (1.5, 0.5), (0.4, 0.25)], None, None),
    ]
    for name, modes, main, _ in tests:
        u = synthetic_field(modes)
        d, _ = decompose(u, width)
        total = sum(d['bands'].values()) + d['unresolved_b5'] + d['residual']
        frac = dict((k, (v / total if total else 0.0)) for k, v in d['bands'].items())
        frac['B5'] = d['unresolved_b5'] / total if total else 0.0
        frac['residual'] = d['residual'] / total if total else 0.0
        row = {'test': name, 'width': width, 'total_over_taper': total / d['hf_taper'],
               'fractions': frac}
        if main:
            others = sorted(((k, v) for k, v in frac.items() if k != main),
                            key=lambda kv: -kv[1])
            row['main_frac'] = frac[main]
            row['adjacent_frac'] = others[0][1]
        rows.append(row)
    return rows


def taper_study(widths):
    rows = []
    for w in widths:
        r = synthetic_report(w)
        s = [x for x in r if x['test'] == 'single_B1'][0]
        rows.append({'width': w, 'B1_frac': s['main_frac'],
                     'leak_max': s['adjacent_frac'],
                     'G2a_err': abs(s['total_over_taper'] - 1.0)})
    return rows


# --------------------------------------------------------------------------- ODB
def read_odb(path):
    """(x_ref, y_ref, u3) — 참조좌표와 면외변위. Abaqus python 전용."""
    from odbAccess import openOdb
    odb = openOdb(path=path, readOnly=True)
    try:
        step = None
        for nm in ('Step-GlobalTension', 'Step-Postbuckle'):
            if nm in odb.steps:
                step = odb.steps[nm]
                break
        if step is None:
            step = odb.steps[list(odb.steps.keys())[-1]]
        frame = step.frames[-1]
        inst = odb.rootAssembly.instances['MEMBRANE-1']
        coord = frame.fieldOutputs['COORD']
        disp = frame.fieldOutputs['U']
        c = dict((v.nodeLabel, v.data) for v in coord.values)
        d = dict((v.nodeLabel, v.data) for v in disp.values)
        xs, ys, zs = [], [], []
        for lab in inst.nodes:
            n = lab.label
            if n not in c or n not in d:
                continue
            xs.append(c[n][0] - d[n][0])
            ys.append(c[n][1] - d[n][1])
            zs.append(d[n][2])
        return np.array(xs), np.array(ys), np.array(zs), step.name, frame.frameValue
    finally:
        odb.close()


def read_element_slopes(path):
    """요소 법선 기울기 s=(n_x/|n_z|, n_y/|n_z|), 참조 중심, 면적, theta².

    slope_stats 의 theta_abs 와 같은 정의(요소 법선 vs z축)를 쓴다. 노드 u3 미분보다
    hf_pred_abs 와 직접 대응한다. Abaqus python 전용.
    """
    from odbAccess import openOdb
    odb = openOdb(path=path, readOnly=True)
    try:
        step = None
        for nm in ('Step-GlobalTension', 'Step-Postbuckle'):
            if nm in odb.steps:
                step = odb.steps[nm]
                break
        if step is None:
            step = odb.steps[list(odb.steps.keys())[-1]]
        frame = step.frames[-1]
        inst = odb.rootAssembly.instances['MEMBRANE-1']
        c = dict((v.nodeLabel, v.data) for v in frame.fieldOutputs['COORD'].values)
        d = dict((v.nodeLabel, v.data) for v in frame.fieldOutputs['U'].values)
        sx, sy, cx, cy, ar, th2 = [], [], [], [], [], []
        for e in inst.elements:
            conn = list(e.connectivity)
            if len(conn) not in (3, 4):
                continue
            if any(n not in c or n not in d for n in conn):
                continue
            ref = np.array([c[n] - d[n] for n in conn])       # 참조좌표
            cur = np.array([c[n] for n in conn])              # 변형좌표
            if len(conn) == 4:
                cr = np.cross(cur[2] - cur[0], cur[3] - cur[1])
                area = 0.5 * float(np.linalg.norm(cr))
            else:
                cr = np.cross(cur[1] - cur[0], cur[2] - cur[0])
                area = 0.5 * float(np.linalg.norm(cr))
            nrm = float(np.linalg.norm(cr))
            if area <= 1e-14 or nrm <= 0.0:
                continue
            n = cr / nrm
            if n[2] < 0.0:                                    # 면외 반전 제거(hf 와 동일)
                n = -n
            nz = abs(float(n[2]))
            if nz < 1e-9:
                continue
            sx.append(float(n[0] / nz))
            sy.append(float(n[1] / nz))
            cen = ref.mean(axis=0)
            cx.append(float(cen[0]))
            cy.append(float(cen[1]))
            ar.append(area)
            th2.append(math.atan2(math.hypot(n[0], n[1]), nz) ** 2)
        return (np.array(sx), np.array(sy), np.array(cx), np.array(cy),
                np.array(ar), np.array(th2), step.name, frame.frameValue)
    finally:
        odb.close()


def grid_slopes(sx, sy, cx, cy, ar):
    """요소 중심 -> 200x100 셀. 셀별 면적가중 평균 sx,sy, 면적합."""
    SX = np.zeros((NY - 1, NX - 1))
    SY = np.zeros((NY - 1, NX - 1))
    AA = np.zeros((NY - 1, NX - 1))
    for a, b, x, y, w in zip(sx, sy, cx, cy, ar):
        ix = min(int(x / H), NX - 2)
        iy = min(int(y / H), NY - 2)
        SX[iy, ix] += a * w
        SY[iy, ix] += b * w
        AA[iy, ix] += w
    ok = AA > 0.0
    SX[ok] /= AA[ok]
    SY[ok] /= AA[ok]
    return SX, SY, AA, ok


def decompose_slopes(SX, SY, AA, width):
    """slope 벡터장 -> 대역별 hf_b. sum_b hf_b = K*mean(|s|^2) (면적가중)."""
    ny, nx = SX.shape
    m = np.zeros((ny, nx))
    valid = AA > 0.0
    m[valid] = 1.0
    if width > 0:
        X, Y = np.meshgrid((np.arange(nx) + 0.5) * H, (np.arange(ny) + 0.5) * H)
        #   셀 격자(200x100)용 테이퍼: 노드 격자(201x101)와 shape 가 다르다.
        d = np.minimum(np.minimum(Y, X - Y), BASE - X - Y)
        t = np.clip(d / float(width), 0.0, 1.0)
        m = (0.5 * (1.0 - np.cos(np.pi * t))) * m
    #   셀 면적가중: 삼각형 경계 셀은 면적이 h^2 의 절반이다. sqrt(AA/h^2) 를 곱해야
    #   Parseval 합이 요소 면적가중 평균과 정확히 일치한다(균일 가중은 저손실 점에서 19% 오차).
    wa = AA / (H * H)
    sw = np.sqrt(np.maximum(wa, 0.0))
    Ux = np.fft.fft2(SX * sw * m)
    Uy = np.fft.fft2(SY * sw * m)
    kx = 2 * math.pi * np.fft.fftfreq(nx, d=H)
    ky = 2 * math.pi * np.fft.fftfreq(ny, d=H)
    KX, KY = np.meshgrid(kx, ky)
    K = np.sqrt(KX ** 2 + KY ** 2)
    power = (np.abs(Ux) ** 2 + np.abs(Uy) ** 2) / float(nx * ny) ** 2
    _sum_a = float(np.sum(AA))
    scale = (H * H * float(nx * ny) / _sum_a) if _sum_a > 0 else 1.0
    hf_t = K_SLOPE * float(np.sum(power)) * scale
    lam = np.full_like(K, np.inf)
    nz = K > 0
    lam[nz] = 2 * math.pi / K[nz]
    out = {'hf_taper': hf_t, 'bands': {}, 'unresolved': 0.0}
    covered = np.zeros_like(K, dtype=bool)
    for name, lo, hi in BANDS:
        sel = (K > 0) & (lam > lo) & (lam <= hi)
        out['bands'][name] = K_SLOPE * float(np.sum(power[sel])) * scale
        covered |= sel
    sel5 = (lam > B5[1]) & (lam <= B5[2])
    out['unresolved_b5'] = K_SLOPE * float(np.sum(power[sel5])) * scale
    covered |= sel5
    out['residual'] = K_SLOPE * float(np.sum(power[(K > 0) & ~covered])) * scale
    #   taper 없는 참조값(면적가중 mean |s|²) — G2b 비교용
    wsum = float(np.sum(AA[valid]))
    out['hf_ref'] = (K_SLOPE * float(np.sum(AA[valid] * (SX[valid] ** 2 + SY[valid] ** 2)))
                     / wsum) if wsum > 0 else 0.0
    return out


def run_odb_slope(path, tag, width, hf_odb):
    sx, sy, cx, cy, ar, th2, step, fv = read_element_slopes(path)
    SX, SY, AA, ok = grid_slopes(sx, sy, cx, cy, ar)
    d = decompose_slopes(SX, SY, AA, width)
    total = sum(d['bands'].values()) + d['unresolved_b5'] + d['residual']
    res = {'tag': tag, 'odb': os.path.basename(path), 'step': step, 'frame': fv,
           'mode': 'element_slope', 'taper_m': width, 'hf_taper': d['hf_taper'],
           'hf_ref': d['hf_ref'], 'hf_odb_reported': hf_odb,
           'g2b_ratio': (d['hf_taper'] / hf_odb) if hf_odb else None,
           'bands': d['bands'], 'b5_unresolved': d['unresolved_b5'],
           'residual': d['residual'], 'sum_bands': total,
           'fractions': dict((k, (v / total if total else 0.0)) for k, v in d['bands'].items()),
           'fraction_b5_unresolved': (d['unresolved_b5'] / total if total else 0.0)}
    res['fractions']['residual'] = d['residual'] / total if total else 0.0
    return res


def grid_from_nodes(xs, ys, u3):
    """노드를 201x101 격자에 매핑. 격자 밖/중복은 오류로 드러낸다."""
    ix = np.rint(xs / H).astype(int)
    iy = np.rint(ys / H).astype(int)
    ok = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
    if not np.all(ok):
        raise RuntimeError('격자 밖 노드 %d개' % int(np.sum(~ok)))
    u = np.zeros((NY, NX), dtype=np.float64)
    seen = np.zeros((NY, NX), dtype=bool)
    for i, j, z in zip(iy, ix, u3):
        if seen[i, j]:
            raise RuntimeError('중복 노드 매핑 (%d,%d)' % (i, j))
        seen[i, j] = True
        u[i, j] = z
    _, _, inside = triangle_mask()
    if int(np.count_nonzero(seen)) != int(np.count_nonzero(inside)):
        raise RuntimeError('노드 수 %d != 삼각 격자 %d'
                           % (int(np.count_nonzero(seen)), int(np.count_nonzero(inside))))
    return u


def run_odb(path, tag, width, hf_odb):
    xs, ys, u3, step, fv = read_odb(path)
    u = grid_from_nodes(xs, ys, u3)
    d, _ = decompose(u, width)
    total = sum(d['bands'].values()) + d['unresolved_b5'] + d['residual']
    res = {'tag': tag, 'odb': os.path.basename(path), 'step': step, 'frame': fv,
           'taper_m': width, 'hf_taper': d['hf_taper'],
           'hf_odb_reported': hf_odb,
           'g2b_ratio': (d['hf_taper'] / hf_odb) if hf_odb else None,
           'bands': d['bands'], 'b5_unresolved': d['unresolved_b5'],
           'residual': d['residual'], 'sum_bands': total,
           'fractions': dict((k, (v / total if total else 0.0))
                             for k, v in d['bands'].items()),
           'fraction_b5_unresolved': (d['unresolved_b5'] / total if total else 0.0)}
    res['fractions']['residual'] = d['residual'] / total if total else 0.0
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--synthetic', action='store_true')
    ap.add_argument('--taper', type=float, default=1.0, help='taper 폭 [m]')
    ap.add_argument('--taper-study', default='0.5,1.0,1.5,2.0')
    ap.add_argument('--odb')
    ap.add_argument('--mode', choices=('element_slope', 'u3'), default='element_slope')
    ap.add_argument('--tag', default='?')
    ap.add_argument('--hf-odb', type=float, default=0.0,
                    help='results 의 slope_f1_hf_pred (G2b 기준)')
    ap.add_argument('--all', action='store_true', help='results/odb 의 HF_*.odb 전부')
    ap.add_argument('--out')
    a = ap.parse_args(argv)

    if a.synthetic:
        widths = [float(x) for x in a.taper_study.split(',') if x.strip()]
        study = taper_study(widths)
        print('width  B1_frac  leak_max  G2a_err')
        for r in study:
            print('  %.2f   %.4f   %.4f    %.2e'
                  % (r['width'], r['B1_frac'], r['leak_max'], r['G2a_err']))
        chosen = [r for r in study if r['B1_frac'] >= 0.95 and r['leak_max'] <= 0.02]
        print('chosen widths (G1 pass):', [r['width'] for r in chosen] or 'NONE')
        if a.out:
            with io.open(a.out, 'w', encoding='utf-8') as f:
                f.write(json.dumps({'study': study, 'detail': synthetic_report(a.taper)},
                                   ensure_ascii=False, indent=1))
        return 0

    if a.all:
        import glob
        rows = []
        for odb in sorted(glob.glob(os.path.join(HERE, 'results', 'odb', 'HF_*.odb'))):
            tag = os.path.basename(odb)[3:-len('_u3free.odb')]
            rp = os.path.join(HERE, 'results', 'results_%s.txt' % tag)
            hf_pred = 0.0
            if os.path.exists(rp):
                d = {}
                with io.open(rp, encoding='utf-8', errors='replace') as f:
                    for line in f:
                        if '=' in line:
                            k, v = line.strip().split('=', 1)
                            d[k] = v
                try:
                    hf_pred = float(d.get('slope_f1_hf_pred', 0.0))
                except Exception:
                    hf_pred = 0.0
            r = (run_odb_slope if a.mode == 'element_slope' else run_odb)(
                odb, tag, a.taper, hf_pred)
            rows.append(r)
            print('%-20s g2b=%.3f  B1=%.2e B2=%.2e B3=%.2e B4=%.3f  B5(unres)=%.3f'
                  % (tag, r['g2b_ratio'] or float('nan'), r['fractions']['B1'],
                     r['fractions']['B2'], r['fractions']['B3'], r['fractions']['B4'],
                     r['fraction_b5_unresolved']))
        if a.out:
            with io.open(a.out, 'w', encoding='utf-8') as f:
                f.write(json.dumps(rows, ensure_ascii=False, indent=1))
        return 0

    if not a.odb:
        raise SystemExit('--odb 또는 --synthetic 또는 --all 이 필요하다')
    res = (run_odb_slope if a.mode == 'element_slope' else run_odb)(
        a.odb, a.tag, a.taper, a.hf_odb)
    print(json.dumps(res, ensure_ascii=False, indent=1))
    if a.out:
        with io.open(a.out, 'w', encoding='utf-8') as f:
            f.write(json.dumps(res, ensure_ascii=False, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
