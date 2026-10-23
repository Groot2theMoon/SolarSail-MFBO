# -*- coding: utf-8 -*-
"""hf_by_state.py — hf 를 LF 상태(taut/wrinkled/slack) 별로 분해한다.

왜
    hf 는 요소별 손실의 면적가중 평균이라 상태별로 **정확히 가법**이다.
    taut 영역은 압축이 없어 주름이 생길 수 없으므로, 그 영역의 손실은
    전부 '전역 곡면' 성분이다. 즉 hf_taut 는 전역 곡면 기저의 **실측값**이다.

        hf_total = Σ_s  f_s * hf_s        (s = taut, wrinkled, slack)

    판정
        hf_taut/hf_total ≈ f_taut(0.18)  -> 손실이 상태와 무관 = 전역 곡면 지배
        hf_taut/hf_total << f_taut        -> 주름 영역에 집중 = 주름 지배

사용법
    abaqus python hf_by_state.py <LF.odb> <HF.odb> [--json out.json]

주의
    - 두 ODB 는 같은 메쉬여야 한다(요소 라벨 1:1 대응).
    - 상태는 LF(평평, u3=0) 기준이다. 'LF 예측 상태'로 읽어야 한다.
    - hf 정의/정규화는 eval_abaqus.calc_thrust_loss 와 동일하게 맞춘다
      (분자=변형 면적, 분모=미변형 기준 면적, IDEAL=A0+2R0).
"""
import sys
import json

import numpy as np
from odbAccess import openOdb

from eval_abaqus import (A0, R0, IDEAL_THRUST_COEFF, U_SUN, ZK_NU,
                         LF_STEP_PREFERENCE, HF_STEP_PREFERENCE,
                         _zk_reference_evol, _zk_get_e_field)


_P = '[state]'


# ------------------------------------------------------------------ 유틸
def pick_step(odb, pref):
    """선호 목록 -> (stepName, step). 프레임 있는 첫 스텝. 없으면 (None, None)."""
    for name in pref:
        if name in odb.steps.keys() and len(odb.steps[name].frames) > 0:
            return name, odb.steps[name]
    for name in odb.steps.keys():
        if len(odb.steps[name].frames) > 0:
            return name, odb.steps[name]
    return None, None


def get_instance(odb):
    for name in odb.rootAssembly.instances.keys():
        inst = odb.rootAssembly.instances[name]
        if len(inst.elements) > 0:
            return inst
    return None


def element_connectivity(inst):
    """(labels(N,), conns: list[list[int]]) — 요소 라벨 순서와 conns 순서 대응."""
    labs, conns = [], []
    for el in inst.elements:
        labs.append(el.label)
        conns.append([nd.label for nd in el.connectivity])
    return np.asarray(labs, dtype=np.int64), conns


def node_coords(frame, maxlab):
    nc = np.zeros((maxlab + 1, 3), dtype=np.float64)
    for block in frame.fieldOutputs['COORD'].bulkDataBlocks:
        if block.nodeLabels is not None:
            nc[block.nodeLabels] = block.data
    return nc


def elem_areas_normals(frame, conns, maxlab):
    """calc_thrust_loss 와 같은 식: quad 는 (v2-v0)x(v3-v1), tri 는 (v1-v0)x(v2-v0)."""
    nc = node_coords(frame, maxlab)
    n = len(conns)
    areas = np.zeros(n, dtype=np.float64)
    normals = np.zeros((n, 3), dtype=np.float64)
    for i, c in enumerate(conns):
        v = nc[np.asarray(c, dtype=np.int64)]
        if len(c) == 4:
            cr = np.cross(v[2] - v[0], v[3] - v[1])
        else:
            cr = np.cross(v[1] - v[0], v[2] - v[0])
        nrm = float(np.linalg.norm(cr))
        areas[i] = 0.5 * nrm
        if nrm > 1e-12:
            normals[i] = cr / nrm
    normals[normals[:, 2] < 0.0] *= -1.0          # 뒤집기 보정(생산 코드와 동일)
    return areas, normals


def elem_strains(frame, inst):
    """요소별 평균 (E11,E22,E12) + 필드명. integration point -> 요소 평균."""
    fld, fname, why = _zk_get_e_field(frame, inst)
    if fld is None:
        return None, None, why
    acc = {}
    for val in fld.values:
        lab = val.elementLabel
        d = val.data
        a = acc.get(lab)
        if a is None:
            acc[lab] = [float(d[0]), float(d[1]), float(d[3]), 1.0]
        else:
            a[0] += float(d[0])
            a[1] += float(d[1])
            a[2] += float(d[3])
            a[3] += 1.0
    labels = np.asarray(sorted(acc.keys()), dtype=np.int64)
    arr = np.asarray([acc[l] for l in labels], dtype=np.float64)
    cnt = arr[:, 3]
    return labels, (arr[:, 0] / cnt, arr[:, 1] / cnt, arr[:, 2] / cnt), fname


def principal_and_states(e11, e22, e12, field_name):
    """주성분 -> (LE 면) Green-Lagrange 변환 -> taut/wrinkled/slack 마스크."""
    mean_e = 0.5 * (e11 + e22)
    rad = np.sqrt((0.5 * (e11 - e22)) ** 2 + e12 ** 2)
    e1 = mean_e + rad
    e2 = mean_e - rad
    #   고유벡터는 LE/E 에서 동일하므로 주성분만 변환하면 정확하다.
    if field_name == 'LE':
        e1 = 0.5 * np.expm1(2.0 * e1)
        e2 = 0.5 * np.expm1(2.0 * e2)
    m_taut = (e2 > 0.0)
    m_wr = (e1 > 0.0) & (e2 <= 0.0)
    m_sl = (e1 <= 0.0)
    return e1, e2, m_taut, m_wr, m_sl


# ------------------------------------------------------------------ 메인
def main(argv):
    if len(argv) < 2:
        print(_P + ' 사용법: abaqus python hf_by_state.py <LF.odb> <HF.odb> [--json out.json]')
        return 2
    lf_path, hf_path = argv[0], argv[1]
    js_path = None
    if '--json' in argv:
        js_path = argv[argv.index('--json') + 1]

    # ---------- LF: 상태
    lf = openOdb(lf_path, readOnly=True)
    try:
        lf_inst = get_instance(lf)
        lf_sname, lf_step = pick_step(lf, LF_STEP_PREFERENCE)
        if lf_step is None:
            print(_P + ' !!! LF 스텝 없음')
            return 1
        lf_frame = lf_step.frames[-1]
        got = elem_strains(lf_frame, lf_inst)
        if got[0] is None:
            print(_P + ' !!! LF 변형률 없음: %s' % got[2])
            return 1
        lf_labels, (e11, e22, e12), lf_fname = got
        e1, e2, m_taut, m_wr, m_sl = principal_and_states(e11, e22, e12, lf_fname)
        wref = _zk_reference_evol(lf, lf_inst)
        print(_P + ' LF: %s / step=%s frame=%d / field=%s / 요소 %d'
              % (lf_path, lf_sname, len(lf_step.frames), lf_fname, len(lf_labels)))
        if not wref:
            print(_P + ' !!! LF frameValue=0 EVOL 없음 -> 기준면적 가중 불가')
            return 1
        print(_P + ' LF 상태: taut=%.4f wrinkled=%.4f slack=%.4f'
              % (float(np.sum(m_taut)) / len(m_taut),
                 float(np.sum(m_wr)) / len(m_wr),
                 float(np.sum(m_sl)) / len(m_sl)))
    finally:
        lf.close()

    state_by_label = {}
    for i, lab in enumerate(lf_labels):
        if m_wr[i]:
            state_by_label[int(lab)] = 'wrinkled'
        elif m_taut[i]:
            state_by_label[int(lab)] = 'taut'
        else:
            state_by_label[int(lab)] = 'slack'

    # ---------- HF: 요소별 손실
    hf = openOdb(hf_path, readOnly=True)
    try:
        hf_inst = get_instance(hf)
        hf_sname, hf_step = pick_step(hf, HF_STEP_PREFERENCE)
        if hf_step is None:
            print(_P + ' !!! HF 스텝 없음')
            return 1
        hf_frame = hf_step.frames[-1]
        hf_labels, conns = element_connectivity(hf_inst)
        maxlab = 0
        for c in conns:
            maxlab = max(maxlab, max(c))
        areas, normals = elem_areas_normals(hf_frame, conns, maxlab)

        #   미변형 기준 면적: COORD frameValue=0
        f0 = None
        for sn in hf.steps.keys():
            for fr in hf.steps[sn].frames:
                if abs(float(fr.frameValue)) < 1e-12 and 'COORD' in fr.fieldOutputs.keys():
                    f0 = fr
                    break
            if f0 is not None:
                break
        if f0 is None:
            print(_P + ' !!! HF frameValue=0 COORD 없음 -> hf 정규화 불가')
            return 1
        areas0, _ = elem_areas_normals(f0, conns, maxlab)
        area_ref_total = float(np.sum(areas0))
        print(_P + ' HF: %s / step=%s frame=%d / 요소 %d / 미변형 면적=%.6f m2'
              % (hf_path, hf_sname, len(hf_step.frames), len(hf_labels), area_ref_total))

        cos_t = np.clip(normals @ (-np.asarray(U_SUN, dtype=np.float64)), 0.0, 1.0)
        #   요소 손실(면적×):  A_e*(IDEAL - A0*c - 2R0*c^3)
        loss_e = areas * (IDEAL_THRUST_COEFF - A0 * cos_t - 2.0 * R0 * cos_t ** 3)
        hf_total = float(np.sum(loss_e) / (area_ref_total * IDEAL_THRUST_COEFF))
        print(_P + ' HF hf(재계산) = %.10e' % hf_total)
    finally:
        hf.close()

    # ---------- 조인 + 집계
    wref_sum = float(sum(wref.values()))
    agg = {}
    missing = 0
    for i, lab in enumerate(hf_labels):
        st = state_by_label.get(int(lab))
        if st is None:
            missing += 1
            continue
        a = agg.get(st)
        if a is None:
            a = agg[st] = {'loss': 0.0, 'aref': 0.0, 'adeform': 0.0,
                           'n': 0, 'sx': 0.0, 'sy': 0.0}
        a['loss'] += float(loss_e[i])
        a['adeform'] += float(areas[i])
        a['n'] += 1
        a['aref'] += float(wref.get(int(lab), 0.0))

    print()
    print(_P + ' %-10s %8s %8s %14s %8s %10s' %
          ('state', 'f_s', 'n', 'hf_s', 'w_s', 'share'))
    print(_P + ' ' + '-' * 62)
    denom = area_ref_total * IDEAL_THRUST_COEFF
    tot_w = 0.0
    rows = {}
    for st in ('taut', 'wrinkled', 'slack'):
        a = agg.get(st)
        if a is None or a['aref'] <= 0.0:
            print(_P + ' %-10s %8s %8s %14s %8s %8s' % (st, '-', '-', '-', '-', '-'))
            continue
        f_s = a['aref'] / wref_sum
        w_s = a['loss'] / denom
        hf_s = w_s / f_s if f_s > 0.0 else float('nan')
        share = w_s / hf_total if hf_total > 0.0 else float('nan')
        tot_w += w_s
        rows[st] = {'f_s': f_s, 'n': a['n'], 'hf_s': hf_s,
                    'w_s': w_s, 'share': share,
                    'area_deformed': a['adeform']}
        print(_P + ' %-10s %8.4f %8d %14.6e %8.3f %9.2f%%'
              % (st, f_s, a['n'], hf_s, w_s, 100.0 * share))
    print(_P + ' ' + '-' * 62)
    print(_P + ' %-10s %8.4f %8s %14.6e %8.3f' %
          ('합계', 1.0, '', hf_total, tot_w))

    if missing:
        print(_P + ' (경고) LF 에 없는 HF 요소 라벨 %d개 — 집계에서 제외' % missing)

    # ---------- 판정
    t = rows.get('taut')
    if t:
        r = t['share']
        f = t['f_s']
        print()
        print(_P + ' 판정: hf_taut 비중=%.2f%% vs taut 면적비=%.2f%%' % (100.0 * r, 100.0 * f))
        if r > 0.7 * f:
            print(_P + '   -> 손실이 상태와 무관하게 분포 = **전역 곡면 지배**')
        elif r < 0.3 * f:
            print(_P + '   -> 손실이 wrinkled/slack 에 집중 = **주름 지배**')
        else:
            print(_P + '   -> 혼합. wrinkled 영역 hf_s 와 I_w 를 따로 비교해야 한다')

    if js_path:
        with open(js_path, 'w') as fh:
            json.dump({'lf': lf_path, 'hf': hf_path,
                       'hf_total': hf_total, 'hf_total_sum': tot_w,
                       'area_ref_total': area_ref_total,
                       'by_state': rows}, fh, indent=1)
        print(_P + ' JSON 저장: %s' % js_path)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
