# -*- coding: utf-8 -*-
"""slope_stats.py - HF ODB 의 요소별 기울기 분포 진단 + hf 사슬 검증 (§19.101)

무엇을 하는가
    HF(후좌굴) ODB 의 변형 프레임에서 요소별 법선 n_e 를 참조구현
    (eval_abaqus.calc_thrust_loss) 과 **문자 그대로 같은 식**으로 계산하고, 두 가지
    기울기를 **구분해서** 분포를 낸다.

      절대 기울기 theta_abs,e = angle(n_e, z_hat)
          설계 평면(z_hat) 기준. 추력손실 hf 와 직접 연결되는 양이다.
          거시 형상(전막 휨/비틀림)과 주름이 **함께** 들어간다.
      요동 기울기 theta_fluct,e = angle(n_e, n_local_mean)
          주변 국소 평균 법선 기준. 국소 평균은 요소 중심 기준 k-최근접
          (반경 3*요소크기) 이웃 법선의 면적가중 평균이다. 거시 형상을 빼고
          **순수 주름(요동)** 몫만 남긴다.

    그리고 이 스크립트의 목적은 마지막 '사슬 검증'이다.

      hf_actual     = eval_abaqus.calc_thrust_loss(...)   <- 참조구현 그대로
      hf_pred_abs   = k * <theta_abs^2>,   k = 1.461716 (면적가중)
      hf_pred_fluct = k * <theta_fluct^2>
      hf_pred_exact      = 1 - Sum(A_e*f_e) / (area_ref*IDEAL),  f_e = A0*c_e + 2*R0*c_e^3
                           c_e = clip(n_e . (-U_SUN), 0, 1)          (요구된 식 그대로)
      hf_pred_exact_norm = 1 - Sum(A_e*f_e) / (Sum(A_e)*IDEAL) = 1 - <f_e>_A  (면적 정규화판)

    [주의 - 왜 두 개인가] 요구식 hf_pred_exact 는 분모를 **미변형** area_ref 로 고정하는데 분자는
      **변형** A_e 를 쓴다. 참조구현은 loss = Sum(A_e)*IDEAL - Sum(A_e*f_e) 이므로
          hf_actual = (Sum(A_e)/area_ref) * (1 - <f_e>_A)
      이고, 요구식과는 앞의 (Sum(A_e)/area_ref) 항만큼 다르다. 실측 (0.5,0.3) 에서 이 비는
      1.0000247 인데 hf(=3.94e-5)와 **같은 자릿수**라 비율이 0.374 로 나와 구현 드리프트로
      오독된다. 그래서 면적 정규화판(분모 = Sum(A_e))을 함께 내고, 그 비율을 드리프트 지표로 쓴다
      (1 에 가까우면 이 스크립트의 SRP 식이 참조구현과 같다는 뜻). 면적비도 함께 찍는다.

    k 의 출처(유도): 작은 각 전개에서
        1 - (A0*cos + 2*R0*cos^3)/(A0+2*R0)
          = [ (A0/2 + 3*R0) / IDEAL_THRUST_COEFF ] * theta^2 + O(theta^4)
      이므로 k_theory = (A0/2 + 3*R0)/IDEAL_THRUST_COEFF = 1.461716.
      이 스크립트는 주어진 k = 1.461716 을 쓰되 k_theory 도 함께 찍어 대조한다.

기대(역산으로 확립된 관계)
    hf_pred_abs / hf_actual   는 약 1   (거시 형상 포함 -> hf 와 같은 양)
    hf_pred_fluct / hf_actual 는 더 작다(순수 주름 몫만)

무엇을 더 내는가
    - 면적가중 분포 통계(mean/rms/p50/p95/p99/max, 라디안 + 도)
    - 10x10 격자 패치별 rms + 최악 패치(위치 포함)
    - 방위각 히스토그램 18빈(법선 수평성분) + 배향성 비(0=등방, 1=정렬)
    - 압축장 요약 1줄(R-13 정의 재사용: minP<0 면적비, 면내평균<0 면적비, minP 면적가중 평균)

읽는 기하 (eval_abaqus.get_hf / eval_abaqus_avg.py 와 같은 방식)
    odbAccess.openOdb(readOnly) -> odb.rootAssembly.instances['MEMBRANE-1'] ->
    instance.elements 의 connectivity, frame.fieldOutputs['COORD'].bulkDataBlocks 좌표,
    frameValue=0 프레임에서 구한 미변형(pristine) 기준면적 area_ref.
    quad: cross(v2-v0, v3-v1) / tri: cross(v1-v0, v2-v0), 면적 = 0.5*|cross|,
    법선 = cross/|cross|, n_z<0 이면 뒤집기(참조구현 기본과 동일).

메모리
    (요소 x 요소) 행렬을 만들지 않는다. 국소 평균은 2D 격자 셀 버킷 + 3x3 이웃 조회로
    구한다. brute force 10100x10100 = 1.02e8 거리계산(float64 = 2.4 GB)을 피한다.

사용 (code/ 에서)
    abaqus python slope_stats.py <odb> [--step NAME] [--frame N] [--json PATH]

재사용 API (§19.102)
    compute_stats(odb_path, step=None, frame=None) -> dict | None
        main() 의 계산부를 그대로 감싼 것. 면적가중 통계 + 사슬검증 + 목적함수 후보 f1~f4 를
        dict 로 돌려준다. run_abaqus.py 가 성공한 HF 런에서 **인프로세스로** 호출해 점별
        results txt 에 'slope_' 접두로 아카이빙한다(scalar_items). 서브프로세스 금지.
    목적함수 후보(전부 면적가중 theta_abs^2 기반, 단위 = hf 와 동일):
        f1_hf_pred     = k*mean(theta_abs^2)              (= 기존 hf_pred_abs)
        f2_p95         = k*q95(theta_abs^2)
        f3_worst_patch = k*max_patch(mean(theta_abs^2))   (10x10 격자 최악 패치)
        f4_cvar90      = k*mean(theta_abs^2 | theta_abs^2 >= q90(theta_abs^2))
    순서 관계(theta^2 >= 0 이라 항상): f1 <= f2 <= f4,  f1 <= f3.

종료코드
    0 정상 / 1 실패(ODB 열기/스텝/프레임/기하/JSON) / 2 사용법

한계 (정직하게)
    - 기하만 읽는다. 하중/재료/경계조건은 재실행하지 않는다.
    - k 는 작은 각 근사다. 각이 큰 요소(수십 도)에서는 k*theta^2 보다 정확식이 맞다.
      둘 다 출력하므로 차이를 직접 볼 수 있다.
    - 국소 평균 반경은 3*요소크기(고정 상수)다. 파장/스펙트럼을 **가정하지 않는다**.
      대신 반경을 바꾸면 theta_fluct 가 달라진다 - 로그에 쓴 반경을 함께 찍는다.
"""
from __future__ import print_function

import io
import json
import math
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

# eval_abaqus.py 를 **고치지 않고** 상수/참조구현을 그대로 가져온다.
try:
    from eval_abaqus import (A0, R0, IDEAL_THRUST_COEFF, U_SUN,
                             calc_thrust_loss)                              # noqa: F401
    from odbAccess import openOdb                                          # noqa: F401
except ImportError as _e:
    print("[slope] odbAccess/eval_abaqus 를 못 읽었습니다 (%s)" % _e)
    print("[slope] 'abaqus python slope_stats.py <odb>' 로 실행해야 합니다.")
    sys.exit(2)

_PFX = '[slope]'

# ------------------------------------------------------------------ 상수
#   k: hf ~= k * <theta^2>. 역산으로 확립된 값(작은 각 전개와 일치).
K_SLOPE = 1.461716
#   k_theory = (A0/2 + 3*R0)/(A0+2*R0). 위 값을 상수로 쓰되 이론식도 찍어 대조한다.
K_SLOPE_THEORY = (0.5 * A0 + 3.0 * R0) / IDEAL_THRUST_COEFF

INSTANCE_NAME = 'MEMBRANE-1'
DEFAULT_STEP_PREF = ('Step-Postbuckle',)   # 있으면 이것, 없으면 마지막 스텝
RING_FACTOR = 3.0                          # 국소 평균 반경 = RING_FACTOR * 요소크기
N_PATCH = 10                               # 10x10 격자 패치
AZ_BINS = 18                               # 방위각 18빈 (bin = 20도)
STALE_WARN_S = 21600.0                     # 6시간 (§19.94 관례)
#   hf_actual 이 이 값 이하면 비율이 정의되지 않는다 -> nan 으로 둔다(0 으로 나눠 inf 를 만들지 않는다).
#   평평한 프레임(미변형/완전평판)에서는 hf_actual 이 정확히 0 이 아니라 +-1e-19 수준의 부동소수
#   찌꺼기가 될 수 있다(eval_abaqus_avg.py 의 SELF_CHECK_FLOOR 와 같은 문제).
HF_RATIO_FLOOR = 1e-12

USAGE = """usage: abaqus python slope_stats.py <odb> [options]

  --step NAME    읽을 스텝 (기본: Step-Postbuckle 이 있으면 그것, 없으면 마지막 스텝)
  --frame N      스텝 안 프레임 번호, 1-based (기본: 마지막 프레임)
  --json PATH    결과 JSON 경로 (기본 results/slope_<odb stem>_<step>.json)
  -h, --help     이 도움말
"""


# ------------------------------------------------------------------ 인자
def parse_args(argv):
    """수동 파싱. (args, err) 반환. err 가 있으면 사용법 오류(rc=2).

    argparse 를 쓰지 않는다: Abaqus 내장 python 이 구형일 수 있고, 이 저장소의 다른
    프로브(base_state_probe.py)도 sys.argv 를 직접 읽는 관례를 따른다.
    """
    a = {'odb': None, 'step': None, 'frame': None, 'json': None, 'help': False}
    i = 0
    while i < len(argv):
        t = argv[i]
        if t in ('-h', '--help'):
            a['help'] = True
            i += 1
            continue
        if t in ('--step', '--frame', '--json'):
            if i + 1 >= len(argv):
                return None, '%s 뒤에 값이 필요하다' % t
            v = argv[i + 1]
            if t == '--step':
                a['step'] = v
            elif t == '--json':
                a['json'] = v
            else:
                try:
                    n = int(v)
                except ValueError:
                    return None, '--frame 은 정수여야 한다: %r' % v
                if n < 1:
                    return None, '--frame 은 1-based(1 이상)여야 한다: %d' % n
                a['frame'] = n
            i += 2
            continue
        if t.startswith('-'):
            return None, '알 수 없는 옵션: %r' % t
        if a['odb'] is None:
            a['odb'] = t
        else:
            return None, 'ODB 인자는 하나만 받는다(추가: %r)' % t
        i += 1
    return a, None


# ------------------------------------------------------------------ ODB 신원
def odb_identity(odb_path, t_start):
    """[§19.94] ODB 신원 로그 - stale ODB 를 로그만으로 판정 가능하게 한다."""
    rec = {'path': os.path.abspath(odb_path), 'size_bytes': None, 'size_mb': None,
           'mtime': None, 'age_s': None, 'age_human': None, 'stale': None}
    try:
        st = os.stat(rec['path'])
        age = t_start - st.st_mtime
        if age < 600.0:
            age_txt = '%.1f 초' % age
        elif age < 86400.0:
            age_txt = '%.2f 시간' % (age / 3600.0)
        else:
            age_txt = '%.2f 일' % (age / 86400.0)
        rec['size_bytes'] = int(st.st_size)
        rec['size_mb'] = st.st_size / 1048576.0
        rec['mtime'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(st.st_mtime))
        rec['age_s'] = age
        rec['age_human'] = age_txt
        rec['stale'] = bool(age > STALE_WARN_S)
        print(_PFX + " ODB 신원: %s | 크기 %.2f MB (%d B) | mtime %s | 프로브 실행 %s | 나이 %s"
              % (rec['path'], rec['size_mb'], rec['size_bytes'], rec['mtime'],
                 time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(t_start)), age_txt))
        if rec['stale']:
            print(_PFX + " !!! ODB 나이 %s (>6시간) - 이 런에서 만든 파일이 아닐 가능성이 높다"
                  "(stale). 이 숫자를 이 설계점의 값으로 읽지 마라." % age_txt)
    except Exception as e:
        print(_PFX + " ODB 신원 확인 실패(무시): %s" % e)
    return rec


# ------------------------------------------------------------------ 기하
def _node_coords(coord_field, max_node_label):
    """bulkDataBlocks -> (max_label+1, 3) 배열. calc_thrust_loss 와 같은 방식."""
    nc = np.zeros((max_node_label + 1, 3), dtype=np.float64)
    for block in coord_field.bulkDataBlocks:
        if block.nodeLabels is not None:
            nc[block.nodeLabels] = block.data
    return nc


def face_geometry(coords, quads, tris, flip_normals):
    """요소별 면적/법선/중심. **eval_abaqus.calc_thrust_loss 와 문자 그대로 같은 식**.

    quad: cross(v2-v0, v3-v1) / tri: cross(v1-v0, v2-v0), 면적=0.5*|cross|,
    법선=cross/|cross|. n_z<0 이면 뒤집기. 퇴화요소(|cross|~0)는 법선 0, 면적 0.
    중심은 요소 절점 좌표의 평균(= 면적가중 중심과 동일하지 않지만 k-NN 격자용으로 동일).
    """
    areas = []
    normals = []
    centers = []
    if len(quads):
        q = np.asarray(quads, dtype=np.int32)
        v0, v1, v2, v3 = coords[q[:, 0]], coords[q[:, 1]], coords[q[:, 2]], coords[q[:, 3]]
        cross = np.cross(v2 - v0, v3 - v1)
        norms = np.linalg.norm(cross, axis=1)
        areas.append(0.5 * norms)
        nrm = np.zeros_like(cross)
        ok = norms > 1e-12
        nrm[ok] = cross[ok] / norms[ok, np.newaxis]
        normals.append(nrm)
        centers.append(0.25 * (v0 + v1 + v2 + v3))
    if len(tris):
        t = np.asarray(tris, dtype=np.int32)
        v0, v1, v2 = coords[t[:, 0]], coords[t[:, 1]], coords[t[:, 2]]
        cross = np.cross(v1 - v0, v2 - v0)
        norms = np.linalg.norm(cross, axis=1)
        areas.append(0.5 * norms)
        nrm = np.zeros_like(cross)
        ok = norms > 1e-12
        nrm[ok] = cross[ok] / norms[ok, np.newaxis]
        normals.append(nrm)
        centers.append((v0 + v1 + v2) / 3.0)
    if not areas:
        return (np.zeros(0), np.zeros((0, 3)), np.zeros((0, 3)))
    areas = np.concatenate(areas)
    normals = np.concatenate(normals, axis=0)
    centers = np.concatenate(centers, axis=0)
    if flip_normals:
        normals[normals[:, 2] < 0.0] *= -1.0
    return areas, normals, centers


def _pristine_area_ref(odb, instance, connectivity, max_label):
    """[§19.89] 미변형(pristine) 기준면적: frameValue=0 프레임에서 계산.

    eval_abaqus.get_hf 의 같은 블록을 그대로 옮겼다(참조구현과 동일 정의).
    """
    try:
        f0 = None
        for sn in odb.steps.keys():
            for fr in odb.steps[sn].frames:
                if abs(float(fr.frameValue)) < 1e-12:
                    f0 = fr
                    break
            if f0 is not None:
                break
        if f0 is None or 'COORD' not in f0.fieldOutputs.keys():
            return None
        cf = f0.fieldOutputs['COORD'].getSubset(region=instance)
        nc = _node_coords(cf, max_label)
        s = 0.0
        q = [c for c in connectivity if len(c) == 4]
        if q:
            a = np.array(q, dtype=np.int32)
            s += 0.5 * np.linalg.norm(
                np.cross(nc[a[:, 2]] - nc[a[:, 0]], nc[a[:, 3]] - nc[a[:, 1]]), axis=1).sum()
        t = [c for c in connectivity if len(c) == 3]
        if t:
            a = np.array(t, dtype=np.int32)
            s += 0.5 * np.linalg.norm(
                np.cross(nc[a[:, 1]] - nc[a[:, 0]], nc[a[:, 2]] - nc[a[:, 0]]), axis=1).sum()
        return float(s) if s > 0.0 else None
    except Exception as e:
        print(_PFX + " area_ref 판독 실패(%s) -> 변형 면적 폴백" % e)
        return None


def read_slope_geometry(odb_path, step_arg, frame_arg):
    """ODB 를 한 번 열어 좌표/법선/면적/중심/area_ref/hf_actual/응력을 뽑는다. 실패면 None."""
    odb = None
    try:
        odb = openOdb(path=odb_path, readOnly=True)
    except Exception as e:
        print(_PFX + " !!! ODB 열기 실패: %s" % e)
        return None
    try:
        names = list(odb.steps.keys())
        if not names:
            print(_PFX + " !!! ODB 에 step 이 없다.")
            return None
        if step_arg is not None:
            if step_arg not in names:
                print(_PFX + " !!! step '%s' 가 이 ODB 에 없다. 있는 step: %s"
                      % (step_arg, names))
                return None
            sname, how = step_arg, '인자'
        else:
            sname, how = None, None
            for cand in DEFAULT_STEP_PREF:
                if cand in names:
                    sname, how = cand, '기본(있음)'
                    break
            if sname is None:
                sname, how = names[-1], '기본(마지막 step)'
        step = odb.steps[sname]
        frames = step.frames
        nfr = len(frames)
        if nfr == 0:
            print(_PFX + " !!! step '%s' 에 프레임이 없다." % sname)
            return None
        if frame_arg is None:
            fi = nfr
        else:
            fi = frame_arg
            if fi > nfr:
                print(_PFX + " !!! --frame %d 가 범위 밖이다(step '%s' 프레임 %d개, 1-based)."
                      % (fi, sname, nfr))
                return None
        frame = frames[fi - 1]
        fval = float(getattr(frame, 'frameValue', float('nan')))
        print(_PFX + " step=%s [%s] | frame=%d/%d | frameValue=%s"
              % (sname, how, fi, nfr, ('%.6g' % fval) if fval == fval else 'nan'))
        if fval == fval and fval < 0.99:
            print(_PFX + " (주의) 이 프레임 frameValue=%.4f < 0.99 - 미완료 런일 수 있다." % fval)

        try:
            instance = odb.rootAssembly.instances[INSTANCE_NAME]
        except Exception:
            try:
                avail = list(odb.rootAssembly.instances.keys())
            except Exception:
                avail = '?'
            print(_PFX + " !!! instance '%s' 가 없다. 있는 instance: %s" % (INSTANCE_NAME, avail))
            return None

        conns = []
        el_labels = []
        for e in instance.elements:
            conns.append(e.connectivity)
            el_labels.append(e.label)
        quads, labels_q = [], []
        tris, labels_t = [], []
        for i, c in enumerate(conns):
            if len(c) == 4:
                quads.append(c)
                labels_q.append(el_labels[i])
            elif len(c) == 3:
                tris.append(c)
                labels_t.append(el_labels[i])
        if not quads and not tris:
            print(_PFX + " !!! 쓸 수 있는 요소(4절점/3절점)가 없다: 요소 %d개 (연결성 길이 %s)"
                  % (len(conns), sorted(set(len(c) for c in conns)) if conns else []))
            return None
        connectivity = quads + tris

        if instance.nodes:
            max_label = instance.nodes[-1].label + 100
        else:
            max_label = 100000

        area_ref = _pristine_area_ref(odb, instance, connectivity, max_label)

        flip_st = {'n_total': 0, 'n_backward': 0, 'area_total': 0.0, 'area_backward': 0.0}
        try:
            cf = frame.fieldOutputs['COORD'].getSubset(region=instance)
            coords = _node_coords(cf, max_label)
            areas, normals, centers = face_geometry(coords, quads, tris, True)
        except Exception as e:
            print(_PFX + " !!! 변형 프레임 좌표/법선 계산 실패: %s" % e)
            return None
        areas_sum = float(np.sum(areas))
        if areas_sum <= 0.0:
            print(_PFX + " !!! 변형 프레임 총면적이 0 이다. 계산 불가.")
            return None
        if area_ref is None:
            area_ref = areas_sum
            print(_PFX + " (주의) area_ref 를 못 구해 변형 면적 %.6f m^2 로 대체한다." % area_ref)
        else:
            print(_PFX + " 미변형(pristine) 기준면적 area_ref = %.6f m^2 (frameValue=0)" % area_ref)
        print(_PFX + " 기하: 요소 %d개 (quad %d / tri %d) | 변형 총면적 %.6f m^2 | 면적비 변형/미변형 %.6f"
              % (len(areas), len(quads), len(tris), areas_sum, areas_sum / area_ref))

        # 참조구현 그대로 hf 산출 (get_hf 와 같은 호출: flip 기본 ON, 같은 area_ref).
        try:
            hf_actual = float(calc_thrust_loss(frame, instance, connectivity, max_label,
                                               area_ref=area_ref, extra_out=flip_st))
            hf_true = float(calc_thrust_loss(frame, instance, connectivity, max_label,
                                             area_ref=area_ref, flip_normals=False))
        except Exception as e:
            print(_PFX + " !!! eval_abaqus.calc_thrust_loss 호출 실패: %s" % e)
            return None
        fb = (100.0 * flip_st['area_backward'] / flip_st['area_total']
              if flip_st['area_total'] > 0.0 else float('nan'))
        print(_PFX + " hf_actual=%.10e (참조구현 calc_thrust_loss) | hf_true(no-flip)=%.10e | "
              "뒤방향(n_z<0) 면적비 %.4f%% (%d/%d)"
              % (hf_actual, hf_true, fb, flip_st['n_backward'], flip_st['n_total']))

        # 압축장: R-13 정의 재사용용 (base_state_probe.py 의 S(CENTROID) 읽기 방식 참고).
        stress = {}
        try:
            from abaqusConstants import CENTROID
            try:
                sf = frame.fieldOutputs['S'].getSubset(position=CENTROID)
            except Exception:
                sf = frame.fieldOutputs['S']
            acc = {}
            for v in sf.values:
                d = v.data
                if len(d) < 3:
                    continue
                acc.setdefault(v.elementLabel, []).append((float(d[0]), float(d[1]), float(d[2])))
            for lab, lst in acc.items():
                m = len(lst)
                stress[lab] = (sum(x[0] for x in lst) / m,
                               sum(x[1] for x in lst) / m,
                               sum(x[2] for x in lst) / m)
            print(_PFX + " 응력(S, CENTROID): 요소 %d개 확보 (압축장 요약용)" % len(stress))
        except Exception as e:
            stress = {}
            print(_PFX + " (주의) S 필드를 못 읽어 압축장 요약을 생략한다: %s" % e)

        return {'step': sname, 'step_how': how, 'frame_index': fi, 'n_frames': nfr,
                'frame_value': fval, 'n_quad': len(quads), 'n_tri': len(tris),
                'labels': labels_q + labels_t, 'areas': areas, 'normals': normals,
                'centers': centers, 'area_ref': float(area_ref), 'area_sum': areas_sum,
                'hf_actual': hf_actual, 'hf_true_no_flip': hf_true,
                'backward_area_pct': fb, 'stress': stress}
    finally:
        try:
            odb.close()
        except Exception:
            pass


# ------------------------------------------------------------------ 국소 평균
def local_mean_normals(centers, normals, areas):
    """요소 중심 기준 반경 R = 3*요소크기 이웃 법선의 면적가중 평균(단위화).

    2D 격자 셀 버킷(셀 크기 = R) + 3x3 이웃 조회. (요소x요소) 거리행렬을 만들지 않는다.
    거리는 설계평면(x,y) 기준 2D 다 - z 요동이 큰 주름에서도 이웃 정의가 흔들리지 않게.
    자기 자신(d=0)은 항상 포함된다.
    """
    n = int(centers.shape[0])
    out = normals.copy()
    cnt = np.zeros(n, dtype=np.int64)
    if n == 0:
        return out, cnt, float('nan'), float('nan')
    pos = areas > 0.0
    h = math.sqrt(float(np.mean(areas[pos]))) if np.any(pos) else float('nan')
    if not (h == h) or h <= 0.0:
        return out, cnt, h, float('nan')
    R = RING_FACTOR * h
    R2 = R * R
    x0, y0 = float(np.min(centers[:, 0])), float(np.min(centers[:, 1]))
    ix = np.floor((centers[:, 0] - x0) / R).astype(np.int64)
    iy = np.floor((centers[:, 1] - y0) / R).astype(np.int64)
    buckets = {}
    for k in range(n):
        buckets.setdefault((int(ix[k]), int(iy[k])), []).append(k)
    grid = {}
    for key, lst in buckets.items():
        grid[key] = np.asarray(lst, dtype=np.int64)
    for k in range(n):
        acc = np.zeros(3, dtype=np.float64)
        cx, cy = centers[k, 0], centers[k, 1]
        kx, ky = int(ix[k]), int(iy[k])
        m = 0
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                lst = grid.get((kx + dx, ky + dy))
                if lst is None:
                    continue
                d2 = (centers[lst, 0] - cx) ** 2 + (centers[lst, 1] - cy) ** 2
                sel = lst[d2 <= R2]
                if sel.size:
                    w = areas[sel]
                    acc += np.sum(w[:, np.newaxis] * normals[sel], axis=0)
                    m += int(sel.size)
        nrm = math.sqrt(float(np.dot(acc, acc)))
        out[k] = acc / nrm if nrm > 1e-30 else normals[k]
        cnt[k] = m
    return out, cnt, h, R


def tilt_from_z(normals):
    """theta_abs = angle(n, z_hat). 단위 법선 기준."""
    return np.arccos(np.clip(normals[:, 2], -1.0, 1.0))


def angle_between(n1, n2):
    """두 단위벡터 사이 각(라디안)."""
    return np.arccos(np.clip(np.sum(n1 * n2, axis=1), -1.0, 1.0))


# ------------------------------------------------------------------ 통계
def _wquantile(theta, w, p, W):
    order = np.argsort(theta, kind='mergesort')
    th = theta[order]
    cw = np.cumsum(w[order])
    idx = int(np.searchsorted(cw, p * W))
    if idx > th.size - 1:
        idx = th.size - 1
    if idx < 0:
        idx = 0
    return float(th[idx])


def area_stats(theta, w):
    """면적가중 통계. theta [rad], w 면적가중. rad/deg 병기 dict."""
    out = {'n': int(theta.size), 'mean_rad': float('nan'), 'rms_rad': float('nan'),
           'p50_rad': float('nan'), 'p95_rad': float('nan'), 'p99_rad': float('nan'),
           'max_rad': float('nan'), 'min_rad': float('nan'),
           'mean_deg': float('nan'), 'rms_deg': float('nan'), 'p50_deg': float('nan'),
           'p95_deg': float('nan'), 'p99_deg': float('nan'), 'max_deg': float('nan'),
           'min_deg': float('nan')}
    W = float(np.sum(w))
    if theta.size == 0 or W <= 0.0:
        return out
    out['mean_rad'] = float(np.sum(w * theta) / W)
    out['rms_rad'] = math.sqrt(float(np.sum(w * theta * theta) / W))
    out['p50_rad'] = _wquantile(theta, w, 0.50, W)
    out['p95_rad'] = _wquantile(theta, w, 0.95, W)
    out['p99_rad'] = _wquantile(theta, w, 0.99, W)
    out['max_rad'] = float(np.max(theta))
    out['min_rad'] = float(np.min(theta))
    for k in ('mean', 'rms', 'p50', 'p95', 'p99', 'max', 'min'):
        out[k + '_deg'] = math.degrees(out[k + '_rad'])
    return out


def patch_grid(cx, cy, theta, w, n=N_PATCH):
    """10x10 패치별 면적가중 rms [rad]. (grid, worst, bounds).

    grid: list of list, grid[iy][ix] = rms_rad (요소 없으면 nan).
    worst: (iy, ix, rms_rad, center_x, center_y) 또는 None.
    """
    if len(cx) == 0:
        return ([[float('nan')] * n for _ in range(n)], None, None)
    x0, x1 = float(np.min(cx)), float(np.max(cx))
    y0, y1 = float(np.min(cy)), float(np.max(cy))
    dx = (x1 - x0) if x1 > x0 else 1.0
    dy = (y1 - y0) if y1 > y0 else 1.0
    num = [[0.0] * n for _ in range(n)]
    den = [[0.0] * n for _ in range(n)]
    for k in range(len(cx)):
        ix = int((cx[k] - x0) / dx * n)
        iy = int((cy[k] - y0) / dy * n)
        if ix < 0:
            ix = 0
        if ix >= n:
            ix = n - 1
        if iy < 0:
            iy = 0
        if iy >= n:
            iy = n - 1
        num[iy][ix] += w[k] * theta[k] * theta[k]
        den[iy][ix] += w[k]
    grid = []
    worst = None
    for iy in range(n):
        row = []
        for ix in range(n):
            if den[iy][ix] > 0.0:
                r = math.sqrt(num[iy][ix] / den[iy][ix])
                row.append(r)
                if worst is None or r > worst[2]:
                    worst = (iy, ix, r, x0 + (ix + 0.5) * dx / n, y0 + (iy + 0.5) * dy / n)
            else:
                row.append(float('nan'))
        grid.append(row)
    return grid, worst, (x0, x1, y0, y1)


def azimuth_hist(normals, areas, n_bins=AZ_BINS):
    """법선 수평성분 (n_x, n_y) 방위각 18빈. 면적가중 비율.

    (n_x, n_y) 는 높이장 z=u3 에 대해 -grad(u3) 방향이다(수직에서 멀지 않을 때).
    배향성 비 = |Sum A exp(i*2*phi)| / Sum A (2배각: 0=등방, 1=완전 정렬).
    """
    hxy = normals[:, :2]
    phi = np.mod(np.arctan2(hxy[:, 1], hxy[:, 0]), 2.0 * math.pi)
    W = float(np.sum(areas))
    frac = np.zeros(n_bins, dtype=np.float64)
    if W <= 0.0:
        return frac, float('nan'), 0
    idx = np.floor(phi / (2.0 * math.pi) * n_bins).astype(np.int64)
    idx = np.clip(idx, 0, n_bins - 1)
    for k in range(n_bins):
        frac[k] = float(np.sum(areas[idx == k])) / W
    c = float(np.sum(areas * np.cos(2.0 * phi))) / W
    s = float(np.sum(areas * np.sin(2.0 * phi))) / W
    return frac, math.sqrt(c * c + s * s), int(np.argmax(frac))


def _principal(s11, s22, s12):
    """base_state_probe.py 의 _principal 과 같은 식(최대/최소 면내 주응력)."""
    mean = 0.5 * (s11 + s22)
    dev = math.sqrt((0.5 * (s11 - s22)) ** 2 + s12 * s12)
    return mean + dev, mean - dev


def compression_summary(areas, labels, stress):
    """R-13 정의 재사용: 면적가중 minP<0 비율 / 면내평균<0 비율 / minP 면적가중 평균."""
    out = {'area_w_minp_neg': float('nan'), 'area_w_mean_neg': float('nan'),
           'minp_area_mean_pa': float('nan'), 'n': 0}
    if not stress:
        return out
    wm = 0.0
    w_minp = 0.0
    w_mean = 0.0
    acc = 0.0
    n = 0
    for lab, a in zip(labels, areas):
        s = stress.get(lab)
        if s is None or a <= 0.0:
            continue
        _s1, s2 = _principal(s[0], s[1], s[2])
        m = 0.5 * (s[0] + s[1])
        wm += a
        if s2 < 0.0:
            w_minp += a
        if m < 0.0:
            w_mean += a
        acc += a * s2
        n += 1
    if wm > 0.0:
        out['area_w_minp_neg'] = w_minp / wm
        out['area_w_mean_neg'] = w_mean / wm
        out['minp_area_mean_pa'] = acc / wm
    out['n'] = n
    return out


# ------------------------------------------------------------------ 표시
def print_stats(tag, st):
    print(_PFX + "   %-16s mean=%.6g (%.4f deg)  rms=%.6g (%.4f deg)"
          % (tag, st['mean_rad'], st['mean_deg'], st['rms_rad'], st['rms_deg']))
    print(_PFX + "   %-16s p50=%.4f  p95=%.4f  p99=%.4f  max=%.4f  (deg)"
          % ('', st['p50_deg'], st['p95_deg'], st['p99_deg'], st['max_deg']))


def print_patch_map(tag, grid, worst, bounds):
    print(_PFX + " === 패치별 rms [deg] (%s, %dx%d, 면적가중) ===" % (tag, N_PATCH, N_PATCH))
    if bounds is None:
        print(_PFX + "   (요소 없음)")
        return
    x0, x1, y0, y1 = bounds
    print(_PFX + "   x 범위 [%.3f, %.3f] m / y 범위 [%.3f, %.3f] m" % (x0, x1, y0, y1))
    for iy in range(N_PATCH - 1, -1, -1):
        cells = []
        for ix in range(N_PATCH):
            v = grid[iy][ix]
            cells.append('   .' if v != v else '%4.3f' % math.degrees(v))
        print(_PFX + "   y%d %s" % (iy + 1, ' '.join(cells)))
    if worst is not None:
        iy, ix, r, wcx, wcy = worst
        print(_PFX + "   최악 패치: (ix=%d, iy=%d) 중심=(%.3f, %.3f) m  rms=%.4f deg"
              % (ix, iy, wcx, wcy, math.degrees(r)))


def _json_safe(obj):
    """NaN/Inf 를 None 으로 바꾼다.

    JSON 표준에는 NaN/Infinity 가 없다(json.dumps 는 비표준 'NaN' 토큰을 쓴다). 엄격한
    파서(strict JSON)가 그 파일에서 죽는다. 읽지 못한 값/빈 패치는 null 로 남긴다
    (0 으로 채우지 않는다 - 조용한 0 금지).
    """
    if isinstance(obj, dict):
        return dict((k, _json_safe(v)) for k, v in obj.items())
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, bool):
        return obj
    if isinstance(obj, (int, float)):
        x = float(obj)
        return x if math.isfinite(x) else None
    if isinstance(obj, np.generic):
        try:
            return _json_safe(obj.item())
        except Exception:
            return obj
    return obj


# ------------------------------------------------------------------ 재사용 API + 목적함수 후보 (§19.102)
def _odb_ident_quiet(odb_path):
    """ODB 신원 dict (출력 없음). build_json_record 가 ident 를 못 받았을 때 쓴다."""
    rec = {'path': os.path.abspath(odb_path), 'size_bytes': None, 'size_mb': None,
           'mtime': None, 'age_s': None, 'age_human': None, 'stale': None}
    try:
        st = os.stat(rec['path'])
        rec['size_bytes'] = int(st.st_size)
        rec['size_mb'] = st.st_size / 1048576.0
        rec['mtime'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(st.st_mtime))
    except Exception:
        pass
    return rec


def compute_stats(odb_path, step=None, frame=None):
    """HF ODB -> 면적가중 통계 + 사슬검증 + 목적함수 후보 f1~f4 dict. 실패하면 None (§19.102).

    main() 의 계산부를 그대로 감싼 **재사용** 함수다. run_abaqus.py 가 성공한 HF 런에서
    **인프로세스로** 호출한다(서브프로세스 금지 - 이미 Abaqus 안이다). [slope] 접두 로그가
    기하/국소평균/목적함수 후보 몇 줄 나간다.

    반환 dict 의 **스칼라 키**는 run_abaqus.py 가 점별 results txt 에 'slope_' 접두로 아카이빙한다
    (scalar_items). 배열/중첩(dict/list/tuple) 키(grid_abs_rad, patches, theta_abs, note ...)는 제외.
    """
    g = read_slope_geometry(odb_path, step, frame)
    if g is None:
        return None

    areas = g['areas']
    normals = g['normals']
    centers = g['centers']
    valid = np.isfinite(areas) & (areas > 0.0)
    if not np.any(valid):
        print(_PFX + " !!! 유효 요소(면적>0)가 0개다.")
        return None
    if int(np.count_nonzero(valid)) != int(areas.size):
        print(_PFX + " (주의) 퇴화 요소 %d개 제외(면적<=0)."
              % (int(areas.size) - int(np.count_nonzero(valid))))
    v_areas = areas[valid]
    v_norm = normals[valid]

    # ---- (1) 국소 평균 법선 -> 요동 기울기
    mn, cnt, h_elem, R_loc = local_mean_normals(centers, normals, areas)
    print(_PFX + " 국소평균: 요소크기 h=sqrt(mean A)=%.6f m | 반경 R=%.1f*h=%.6f m | "
          "평균 이웃 %.1f개 (격자 셀 버킷, (요소x요소) 행렬 없음)"
          % (h_elem, RING_FACTOR, R_loc, float(np.mean(cnt[valid])) if R_loc == R_loc else float('nan')))

    th_abs = tilt_from_z(normals)[valid]
    th_fluct = angle_between(normals, mn)[valid]

    st_abs = area_stats(th_abs, v_areas)
    st_fluct = area_stats(th_fluct, v_areas)

    sum_a = float(np.sum(v_areas))
    cos_w = float(np.sum(v_areas * np.cos(th_abs)) / sum_a)

    # ---- (2) 사슬 검증 수치 (참조구현 정의 재사용)
    hf_actual = float(g['hf_actual'])
    m2_abs = float(np.sum(v_areas * th_abs * th_abs) / sum_a)
    m2_fl = float(np.sum(v_areas * th_fluct * th_fluct) / sum_a)
    hf_pred_abs = K_SLOPE * m2_abs
    hf_pred_fluct = K_SLOPE * m2_fl

    c_exact = np.clip(np.dot(v_norm, -np.asarray(U_SUN, dtype=np.float64)), 0.0, 1.0)
    f_e = A0 * c_exact + 2.0 * R0 * c_exact ** 3
    area_ratio = sum_a / g['area_ref']
    num_exact = float(np.sum(v_areas * f_e))
    hf_pred_exact = 1.0 - num_exact / (g['area_ref'] * IDEAL_THRUST_COEFF)
    hf_pred_exact_norm = 1.0 - num_exact / (sum_a * IDEAL_THRUST_COEFF)

    def _ratio(num, den):
        return (num / den) if (den is not None and den > HF_RATIO_FLOOR) else float('nan')

    r_abs = _ratio(hf_pred_abs, hf_actual)
    r_fl = _ratio(hf_pred_fluct, hf_actual)
    r_ex = _ratio(hf_pred_exact, hf_actual)
    r_exn = _ratio(hf_pred_exact_norm, hf_actual)

    # ---- (3) 패치별 rms + 최악 패치
    vcx = centers[valid, 0]
    vcy = centers[valid, 1]
    grid_abs, worst_abs, bounds = patch_grid(vcx, vcy, th_abs, v_areas)
    grid_fl, worst_fl, _b2 = patch_grid(vcx, vcy, th_fluct, v_areas)

    # ---- (4) 방위각 히스토그램
    az_frac, orient, dom = azimuth_hist(v_norm, v_areas)

    # ---- (5) 압축장 요약
    comp = compression_summary(areas, g['labels'], g['stress'])

    # ---- (6) 목적함수 후보 f1~f4 (§19.102). 전부 면적가중, 단위 = hf 와 동일.
    #   theta_abs^2 의 면적가중 분포에서 (i) 평균 (ii) 95백분위 (iii) 최악 10x10 패치 평균 (iv) CVaR90.
    #   면적가중 분위수는 _wquantile(정렬 + 누적면적) 로 낸다.
    th2 = th_abs * th_abs
    W = sum_a
    q90 = _wquantile(th2, v_areas, 0.90, W)
    q95 = _wquantile(th2, v_areas, 0.95, W)
    f1_hf_pred = K_SLOPE * m2_abs
    f2_p95 = K_SLOPE * q95
    worst_rms_rad = float(worst_abs[2]) if worst_abs is not None else float('nan')
    f3_worst_patch = K_SLOPE * (worst_rms_rad * worst_rms_rad)
    _sel = th2 >= q90
    _den = float(np.sum(v_areas[_sel]))
    cvar90_mean = (float(np.sum(v_areas[_sel] * th2[_sel]) / _den)
                   if _den > 0.0 else float('nan'))
    f4_cvar90 = K_SLOPE * cvar90_mean

    print(_PFX + " --- 목적함수 후보 f1~f4 (k=%.6f, 면적가중, hf 와 같은 단위) ---" % K_SLOPE)
    print(_PFX + " f1_hf_pred     = %.10e  (k*mean(theta^2) = hf_pred_abs)" % f1_hf_pred)
    print(_PFX + " f2_p95         = %.10e  (k*q95(theta^2), q95=%.6e rad^2)" % (f2_p95, q95))
    print(_PFX + " f3_worst_patch = %.10e  (k*max_patch mean(theta^2), 최악 rms=%.4f deg)"
          % (f3_worst_patch, math.degrees(worst_rms_rad)))
    print(_PFX + " f4_cvar90      = %.10e  (k*mean(theta^2 | theta^2>=q90), q90=%.6e rad^2)"
          % (f4_cvar90, q90))
    if f1_hf_pred > 0.0:
        print(_PFX + " 비 f3/f1 = %.4f   f4/f1 = %.4f   f2/f1 = %.4f   (기대: f1 <= f2 <= f4)"
              % (f3_worst_patch / f1_hf_pred, f4_cvar90 / f1_hf_pred, f2_p95 / f1_hf_pred))

    return {
        'step': g['step'], 'step_how': g['step_how'], 'frame': g['frame_index'],
        'n_frames': g['n_frames'], 'frame_value': g['frame_value'],
        'n_elem': int(areas.size), 'n_quad': g['n_quad'], 'n_tri': g['n_tri'],
        'n_valid': int(np.count_nonzero(valid)),
        'n_degenerate_removed': int(areas.size) - int(np.count_nonzero(valid)),
        'area_ref': float(g['area_ref']), 'area_sum': float(g['area_sum']),
        'area_ratio': float(area_ratio), 'backward_area_pct': float(g['backward_area_pct']),
        'elem_size_m': float(h_elem), 'r_local_m': float(R_loc),
        'mean_neighbors': (float(np.mean(cnt[valid])) if R_loc == R_loc else float('nan')),
        'ring_factor': RING_FACTOR,
        'k_slope': K_SLOPE, 'k_slope_theory': K_SLOPE_THEORY,
        'theta_abs_mean_deg': st_abs['mean_deg'], 'theta_abs_rms_deg': st_abs['rms_deg'],
        'theta_abs_rms_rad': st_abs['rms_rad'],
        'theta_abs_p50_deg': st_abs['p50_deg'], 'theta_abs_p95_deg': st_abs['p95_deg'],
        'theta_abs_p99_deg': st_abs['p99_deg'], 'theta_abs_max_deg': st_abs['max_deg'],
        'theta_abs_mean_rad': st_abs['mean_rad'], 'theta_abs_p95_rad': st_abs['p95_rad'],
        'theta_fluct_rms_deg': st_fluct['rms_deg'], 'theta_fluct_p95_deg': st_fluct['p95_deg'],
        'theta_fluct_p99_deg': st_fluct['p99_deg'], 'theta_fluct_max_deg': st_fluct['max_deg'],
        'cos_theta_abs_area_mean': cos_w,
        'hf_actual': hf_actual, 'hf_true_no_flip': float(g['hf_true_no_flip']),
        'hf_pred_abs': hf_pred_abs, 'hf_pred_fluct': hf_pred_fluct,
        'hf_pred_exact': hf_pred_exact, 'hf_pred_exact_norm': hf_pred_exact_norm,
        'm2_abs': m2_abs, 'm2_fluct': m2_fl,
        'ratio_abs': r_abs, 'ratio_fluct': r_fl,
        'ratio_exact': r_ex, 'ratio_exact_norm': r_exn,
        'worst_patch_rms_deg': (float('nan') if worst_abs is None else math.degrees(worst_abs[2])),
        'worst_patch_xy': ('' if worst_abs is None else '%.6g,%.6g' % (worst_abs[3], worst_abs[4])),
        'worst_patch_ix': (-1 if worst_abs is None else int(worst_abs[1])),
        'worst_patch_iy': (-1 if worst_abs is None else int(worst_abs[0])),
        'worst_patch_fluct_rms_deg': (float('nan') if worst_fl is None else math.degrees(worst_fl[2])),
        'orient_ratio': orient, 'az_dominant_bin': int(dom),
        'minp_frac_neg': comp['area_w_minp_neg'],
        'mean_neg_frac': comp['area_w_mean_neg'],
        'minp_area_mean_pa': comp['minp_area_mean_pa'], 'compression_n': comp['n'],
        'theta2_mean': m2_abs, 'theta2_q90': q90, 'theta2_q95': q95,
        'theta2_cvar90_mean': cvar90_mean, 'theta2_max': float(np.max(th2)),
        'f1_hf_pred': f1_hf_pred, 'f2_p95': f2_p95,
        'f3_worst_patch': f3_worst_patch, 'f4_cvar90': f4_cvar90,
        # ---- 중첩/배열 (출력 + JSON 재사용; scalar_items 는 건너뛴다)
        'theta_abs': st_abs, 'theta_fluct': st_fluct,
        'grid_abs_rad': grid_abs, 'grid_fluct_rad': grid_fl,
        'worst_abs_tuple': worst_abs, 'worst_fluct_tuple': worst_fl,
        'patch_bounds': bounds, 'az_frac': az_frac,
        'hf': {'actual': hf_actual, 'actual_no_flip': float(g['hf_true_no_flip']),
               'pred_abs': hf_pred_abs, 'pred_fluct': hf_pred_fluct,
               'pred_exact': hf_pred_exact, 'pred_exact_norm': hf_pred_exact_norm,
               'm2_abs': m2_abs, 'm2_fluct': m2_fl,
               'area_ratio_deformed_over_ref': area_ratio,
               'ratio_abs': r_abs, 'ratio_fluct': r_fl, 'ratio_exact': r_ex,
               'ratio_exact_norm': r_exn},
        'patches': {'n': N_PATCH,
                    'rms_deg_abs': [[(float('nan') if v != v else math.degrees(v))
                                     for v in row] for row in grid_abs],
                    'rms_deg_fluct': [[(float('nan') if v != v else math.degrees(v))
                                       for v in row] for row in grid_fl],
                    'worst_abs': (None if worst_abs is None else
                                  {'iy': worst_abs[0], 'ix': worst_abs[1],
                                   'rms_deg': math.degrees(worst_abs[2]),
                                   'center_x': worst_abs[3], 'center_y': worst_abs[4]}),
                    'worst_fluct': (None if worst_fl is None else
                                    {'iy': worst_fl[0], 'ix': worst_fl[1],
                                     'rms_deg': math.degrees(worst_fl[2]),
                                     'center_x': worst_fl[3], 'center_y': worst_fl[4]})},
        'azimuth': {'n_bins': AZ_BINS, 'bin_deg': 360.0 / AZ_BINS,
                    'area_frac': [float(x) for x in az_frac],
                    'orient_ratio_double_angle': orient, 'dominant_bin': int(dom)},
        'compression': comp,
        'note': ('theta_abs: n_e vs z_hat (macro+wrinkle, links to hf). '
                 'theta_fluct: n_e vs local area-weighted mean normal within '
                 'R=%.1f*elem_size (wrinkle only). k=(A0/2+3*R0)/IDEAL from small-angle '
                 'expansion of 1-(A0*cos+2*R0*cos^3)/(A0+2*R0).' % RING_FACTOR),
    }


def default_json_path(odb_path, step):
    """slope_stats JSON 기본 경로: results/slope_<odb stem>_<step>.json (기존 관례)."""
    stem = os.path.splitext(os.path.basename(odb_path))[0]
    return os.path.join('results', 'slope_%s_%s.json' % (stem, step))


def write_json_record(rec, out_path):
    """JSON 레코드 -> 파일. 디렉터리 자동 생성. 실패면 예외를 올린다(호출자가 처리)."""
    d = os.path.dirname(os.path.abspath(out_path))
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with io.open(out_path, 'w', encoding='utf-8') as f:
        f.write(json.dumps(_json_safe(rec), ensure_ascii=True, indent=2, sort_keys=True))
    return os.path.abspath(out_path)


def build_json_record(S, odb_path, ident=None, t_start=None):
    """compute_stats 결과(S) -> slope_stats JSON 레코드. 기존 스키마 그대로 + objectives(§19.102)."""
    if ident is None:
        ident = _odb_ident_quiet(odb_path)
    if t_start is None:
        t_start = time.time()
    return {
        'script': 'slope_stats.py', 'section': '19.101',
        'generated': time.strftime('%Y-%m-%d %H:%M:%S'),
        'probe_run_epoch': t_start, 'odb': ident,
        'step': S['step'], 'step_how': S['step_how'], 'frame_index': S['frame'],
        'n_frames': S['n_frames'], 'frame_value': S['frame_value'],
        'n_elem': S['n_elem'], 'n_quad': S['n_quad'], 'n_tri': S['n_tri'],
        'n_degenerate_removed': S['n_degenerate_removed'],
        'area_ref': S['area_ref'], 'area_deformed_sum': S['area_sum'],
        'area_ratio_deformed_pristine': ((S['area_sum'] / S['area_ref'])
                                         if S['area_ref'] else float('nan')),
        'backward_area_pct': S['backward_area_pct'],
        'k_slope': S['k_slope'], 'k_slope_theory': S['k_slope_theory'],
        'ring_factor': S['ring_factor'], 'elem_size_m': S['elem_size_m'],
        'r_local_m': S['r_local_m'], 'mean_neighbors': S['mean_neighbors'],
        'theta_abs': S['theta_abs'], 'theta_fluct': S['theta_fluct'],
        'cos_theta_abs_area_mean': S['cos_theta_abs_area_mean'],
        'hf': S['hf'], 'patches': S['patches'], 'azimuth': S['azimuth'],
        'compression': S['compression'],
        'objectives': {'section': '19.102', 'k_slope': S['k_slope'],
                       'f1_hf_pred': S['f1_hf_pred'], 'f2_p95': S['f2_p95'],
                       'f3_worst_patch': S['f3_worst_patch'], 'f4_cvar90': S['f4_cvar90'],
                       'theta2_mean': S['theta2_mean'], 'theta2_q90': S['theta2_q90'],
                       'theta2_q95': S['theta2_q95'],
                       'theta2_cvar90_mean': S['theta2_cvar90_mean'],
                       'note': 'area-weighted theta_abs^2 statistics: f1=k*mean, '
                               'f2=k*q95, f3=k*worst-10x10-patch-mean, f4=k*CVaR90. hf units.'},
        'note': S['note'],
    }


_SCALAR_SKIP = ('note',)


def scalar_items(S, prefix='slope_', skip=_SCALAR_SKIP):
    """dict -> 'prefix+key=value' 줄 목록(키 정렬). **스칼라만**(int/float/str/bool/np 스칼라).

    NaN/Inf -> 'nan'. 배열/중첩(dict/list/tuple) 키와 skip 키는 건너뛴다(§19.102 점별 txt 아카이빙용).
    """
    out = []
    if not S:
        return out
    for k in sorted(S):
        if k in skip:
            continue
        v = S[k]
        if isinstance(v, bool):
            s = 'True' if v else 'False'
        elif isinstance(v, np.generic):
            try:
                f = float(v)
            except Exception:
                continue
            s = 'nan' if not math.isfinite(f) else ('%.12g' % f)
        elif isinstance(v, float):
            s = 'nan' if not math.isfinite(v) else ('%.12g' % v)
        elif isinstance(v, int):
            s = '%d' % v
        elif isinstance(v, str):
            s = v
        else:
            continue
        out.append('%s%s=%s' % (prefix, k, s))
    return out


# ------------------------------------------------------------------ main
def main(argv=None):
    t_start = time.time()
    if argv is None:
        argv = sys.argv[1:]
    print(_PFX + " build: slope_stats.py (19.101) | numpy %s | U_SUN=%s | A0=%.6f R0=%.6f IDEAL=%.6f"
          % (np.__version__, np.asarray(U_SUN, dtype=float).tolist(),
             A0, R0, IDEAL_THRUST_COEFF))
    print(_PFX + " build+: §19.102 목적함수 후보 f1~f4 + compute_stats() 재사용 API + "
          "run_abaqus.py 자동배선 지원")
    print(_PFX + " k(사용)=%.6f  k_theory=(A0/2+3*R0)/IDEAL=%.6f  차이=%.3e"
          % (K_SLOPE, K_SLOPE_THEORY, abs(K_SLOPE - K_SLOPE_THEORY)))
    args, err = parse_args(argv)
    if err:
        print(_PFX + " 사용법 오류: %s" % err)
        print(USAGE)
        return 2
    if args['help']:
        print(USAGE)
        return 0
    if args['odb'] is None:
        print(_PFX + " 사용법 오류: ODB 경로가 필요하다.")
        print(USAGE)
        return 2
    if not os.path.exists(args['odb']):
        print(_PFX + " !!! ODB 파일이 없다: %s" % os.path.abspath(args['odb']))
        return 1

    ident = odb_identity(args['odb'], t_start)
    S = compute_stats(args['odb'], args['step'], args['frame'])
    if S is None:
        return 1

    #   §19.102: 계산부는 compute_stats() 로 옮겼다. 여기서는 출력/JSON/KV 가 쓰는 이름으로 푼다.
    st_abs = S['theta_abs']
    st_fluct = S['theta_fluct']
    hf_actual = S['hf_actual']
    m2_abs = S['m2_abs']
    m2_fl = S['m2_fluct']
    hf_pred_abs = S['hf_pred_abs']
    hf_pred_fluct = S['hf_pred_fluct']
    hf_pred_exact = S['hf_pred_exact']
    hf_pred_exact_norm = S['hf_pred_exact_norm']
    area_ratio = S['area_ratio']
    r_abs = S['ratio_abs']
    r_fl = S['ratio_fluct']
    r_ex = S['ratio_exact']
    r_exn = S['ratio_exact_norm']
    grid_abs = S['grid_abs_rad']
    grid_fl = S['grid_fluct_rad']
    worst_abs = S['worst_abs_tuple']
    worst_fl = S['worst_fluct_tuple']
    bounds = S['patch_bounds']
    az_frac = S['az_frac']
    orient = S['orient_ratio']
    dom = S['az_dominant_bin']
    comp = S['compression']
    cos_w = S['cos_theta_abs_area_mean']

    print(_PFX + " --- 절대 기울기 theta_abs = angle(n_e, z_hat) [설계평면 기준, hf 와 직접 연결] ---")
    print_stats('theta_abs', st_abs)
    print(_PFX + " --- 요동 기울기 theta_fluct = angle(n_e, n_local_mean) [거시형상 제거, 순수 주름] ---")
    print_stats('theta_fluct', st_fluct)

    print(_PFX + " 면적가중 <cos theta_abs> = %.8f" % cos_w)

    # ---- (2) 사슬 검증 수치 (hf_actual/m2_*/hf_pred_*/ratio_* 는 compute_stats 가 dict 로 반환)

    print(_PFX + " ===== 사슬 검증: hf = k*<theta^2>, k=%.6f =====" % K_SLOPE)
    print(_PFX + " hf_actual      = %.10e   (eval_abaqus.calc_thrust_loss, 참조구현 그대로)"
          % hf_actual)
    print(_PFX + " <theta_abs^2>  = %.10e rad^2  (면적가중)  -> hf_pred_abs   = %.10e"
          % (m2_abs, hf_pred_abs))
    print(_PFX + " <theta_fl^2>   = %.10e rad^2  (면적가중)  -> hf_pred_fluct = %.10e"
          % (m2_fl, hf_pred_fluct))
    print(_PFX + " hf_pred_exact  = %.10e   (요구식 1 - Sum(A_e*f_e)/(area_ref*IDEAL))"
          % hf_pred_exact)
    print(_PFX + " 면적비 Sum(A_e)/area_ref = %.8f  (변형/미변형)"
          % area_ratio)
    print(_PFX + " hf_pred_exact_norm = %.10e   (면적 정규화판 1 - <f_e>_A = 참조구현의 1-<f> 항)"
          % hf_pred_exact_norm)
    if hf_actual > HF_RATIO_FLOOR:
        print(_PFX + " 비율 hf_pred_abs        / hf_actual = %.6f   (기대: 약 1, 거시형상 포함)"
              % r_abs)
        print(_PFX + " 비율 hf_pred_fluct      / hf_actual = %.6f   (기대: < 1, 순수 주름 몫)"
              % r_fl)
        print(_PFX + " 비율 hf_pred_exact      / hf_actual = %.6f   (요구식. 면적비 항 때문에 1 이 아님)"
              % r_ex)
        print(_PFX + " 비율 hf_pred_exact_norm / hf_actual = %.6f   (드리프트 지표: ~면적비이면 정합)"
              % r_exn)
    else:
        print(_PFX + " hf_actual=%.3e <= %.0e 이다(평평/미변형 프레임). 비율은 정의되지 않는다"
              " -> nan 로 둔다(0 으로 나누지 않는다)." % (hf_actual, HF_RATIO_FLOOR))

    # ---- (3) 패치별 rms + 최악 패치 (grid/worst/bounds 는 compute_stats 결과)
    print_patch_map('theta_abs', grid_abs, worst_abs, bounds)
    if worst_fl is not None:
        print(_PFX + " 최악 패치(theta_fluct): (ix=%d, iy=%d) 중심=(%.3f, %.3f) m  rms=%.4f deg"
              % (worst_fl[1], worst_fl[0], worst_fl[3], worst_fl[4],
                 math.degrees(worst_fl[2])))

    # ---- (4) 방위각 히스토그램 (az_frac/orient/dom 는 compute_stats 결과)
    print(_PFX + " 방위각 히스토그램(법선 수평성분, %d빈=%.0f도, 면적가중): 배향성비(2배각)=%.4f "
          "(0=등방, 1=정렬)  최대빈=%d..%d도 (면적비 %.4f)"
          % (AZ_BINS, 360.0 / AZ_BINS, orient, dom * int(360.0 / AZ_BINS),
             (dom + 1) * int(360.0 / AZ_BINS), az_frac[dom]))
    print(_PFX + "   " + ' '.join('%d-%.0f:%.4f' % (int(k * 360.0 / AZ_BINS),
                                                    (k + 1) * 360.0 / AZ_BINS, az_frac[k])
                                  for k in range(AZ_BINS)))
    print(_PFX + "   (수평성분 (n_x,n_y) 는 높이장 z=u3 의 -grad(u3) 방향이다. 2배각 비는 "
          "선(line) 배향성 -> 축 방향 정렬만 본다.)")

    # ---- (5) 압축장 요약 (R-13 정의 재사용; comp 는 compute_stats 결과)
    if comp['n'] > 0:
        print(_PFX + " 압축장 요약(R-13 정의): 면적가중 minP<0 비율=%.4f | 면내평균<0 비율=%.4f | "
              "minP 면적가중 평균=%+.4g Pa (n=%d)"
              % (comp['area_w_minp_neg'], comp['area_w_mean_neg'],
                 comp['minp_area_mean_pa'], comp['n']))
        print(_PFX + "   (출처: base_state_probe.py 의 R-13 정의를 최소 복사 - S(CENTROID) 요소평균, "
              "면적 가중)")
    else:
        print(_PFX + " 압축장 요약: S 값이 없어 nan (정의: minP<0 면적비 / 면내평균<0 면적비 / minP 평균)")

    # ---- (6) JSON (기존 스키마 그대로 + objectives 블록 추가, §19.102)
    out_path = args['json']
    if out_path is None:
        out_path = default_json_path(args['odb'], S['step'])
    rec = build_json_record(S, args['odb'], ident=ident, t_start=t_start)
    try:
        print(_PFX + " JSON 저장: %s" % write_json_record(rec, out_path))
    except Exception as e:
        print(_PFX + " !!! JSON 저장 실패: %s" % e)
        return 1

    # ---- (7) 기계 파싱용 KV 1줄
    kv = {'hf_actual': hf_actual, 'hf_pred_abs': hf_pred_abs,
          'hf_pred_fluct': hf_pred_fluct, 'hf_pred_exact': hf_pred_exact,
          'hf_pred_exact_norm': hf_pred_exact_norm,
          'ratio_abs': r_abs, 'ratio_fluct': r_fl, 'ratio_exact': r_ex,
          'ratio_exact_norm': r_exn, 'area_ratio': area_ratio,
          'theta_abs_rms_deg': st_abs['rms_deg'], 'theta_fluct_rms_deg': st_fluct['rms_deg'],
          'theta_abs_p95_deg': st_abs['p95_deg'],
          'theta_abs_p99_deg': st_abs['p99_deg'], 'theta_abs_max_deg': st_abs['max_deg'],
          'worst_patch_rms_deg': S['worst_patch_rms_deg'],
          'orient_ratio': orient, 'minp_frac_neg': comp['area_w_minp_neg'],
          'area_ref': S['area_ref'], 'n_elem': S['n_elem'],
          'f1_hf_pred': S['f1_hf_pred'], 'f2_p95': S['f2_p95'],
          'f3_worst_patch': S['f3_worst_patch'], 'f4_cvar90': S['f4_cvar90'],
          'theta2_q90': S['theta2_q90'], 'theta2_q95': S['theta2_q95'],
          'step': S['step'], 'frame': S['frame']}
    print('[slope:KV] ' + ' '.join('%s=%s' % (k, kv[k]) for k in sorted(kv)))
    print(_PFX + " 완료 (rc=0)")
    return 0


if __name__ == '__main__':
    sys.exit(main())
