# -*- coding: utf-8 -*-
"""torque_stats.py - HF ODB 의 요소별 SRP 힘벡터 모멘트(토크) 지표 tau_rel (§19.103)

무엇을 하는가
    HF(후좌굴) ODB 의 변형 프레임에서 요소별 법선 n_e 를 참조구현
    (eval_abaqus.calc_thrust_loss) 과 **문자 그대로 같은 식**으로 계산하고,
    그 힘벡터에 **지렛대를 곱해서** 모멘트를 낸다. 참조구현은 요소 힘을 전부 합해
    스칼라(추력손실 hf)만 남기고 **위치 정보를 버린다** - 위치를 살린 것이 이 스크립트다.

        f_e = A_e * ( A0*c_e*U_SUN + 2*R0*c_e^2*(-n_e) ),   c_e = clip(n_e . (-U_SUN), 0, 1)
        M   = sum_e ( x_e - x_ref ) x f_e,                  x_e = 요소 중심(변형 좌표)

    기준점 x_ref 는 두 가지를 **모두** 낸다:
        (a) 면적가중 도심(변형 좌표에서 계산) - 돛이 **자기 무게중심 주위로** 도는 성분.
            균일 법선(평판/균일 기울기)이면 M = 0 항등식이다(tau 가 '형상 대비'를 본다).
        (b) 원점(0,0,0) - 실제 돛은 삼각형(V3=(0,0,0), V2=(20,0,0), V1=(10,10,0))이라
            **원점이 도심이 아니다**. 평판에서도 |M_origin| = |x_c x F| != 0 이다
            (순수 지렛대 항). 두 값을 같이 봐야 '형상이 만든 것'과 '기준점이 만든 것'이 갈린다.

    지표
        tau_abs  = |M|                                  [N m]
        tau_rel  = |M| / ( sum_e |f_e| * L_ref ),   L_ref = sqrt(area_ref)
                   무차원. 힘 총량으로 정규화했으므로 '형상이 만드는 회전 비율'이다.
        면외/면내: M_x, M_y(면외, 돛을 기울이는 축) vs M_z(면내, 돛면 회전축).
                   자세 교란은 면외 지배일 것으로 본다 - 가정하지 않고 매번 비를 찍는다.
        tau_flat = 미변형·평평 프레임(frameValue=0)에 같은 계산.
                   평판 + 대칭 메쉬면 **도심 기준 M ~ 0** 이 sanity anchor 다.

    방향 의존(두 방식 모두)
        (i)  격자: phi(0..60도) x psi(0..330도, 30도 간격) -> tau_rel(s_hat) 분포 + 최악 방향
        (ii) 단위 디스크 MC: 밀도가 cos(phi) 에 비례하는 **전력가중** 무작위 샘플
             (eval_abaqus_avg.sample_dirs('disk') 와 같은 식) -> 평균/SE/분위수

문헌 근거(왜 hf 와 별개 지표인가)
    - Zeiders 2005 (NTRS 20050209955): 평균제곱 기울기는 추력 *크기*의 지표이고
      "...it might not guarantee flat sail performance for **the moments** when there are
      large-scale disturbances..." - 즉 RMSE 하나로 **모멘트**는 환원되지 않는다.
    - Heaton, Miller, Ahmad 2017 (NTRS 20170001488, NEA Scout Thrust and Torque Model):
      0도에서 설계형상과 평탄형상의 **추력이 722 uN 으로 동일**한데, 60도에서 **토크는 2.07배**.
      "sail thrust magnitude is not affected much by shape; sail torque is strongly affected
      by shape." -> 형상이 실제로 때리는 축은 토크다.
    - 그래서 hf(크기 축)와 tau_rel(회전 축)을 **둘 다** 남긴다. 어느 쪽이 최적화 신호인지는
      실측 ODB 의 x_c 스캔으로 판정한다(§19.103 사전 판정 기준: 대칭 x_c 에서 tau_rel 이
      작고 x_c 가 비대칭으로 갈수록 크게 자라면 토크가 목적함수로 승격될 후보다).

읽는 기하 (eval_abaqus.get_hf / slope_stats.py 와 같은 방식 - 그대로 복제)
    odbAccess.openOdb(readOnly) -> odb.rootAssembly.instances['MEMBRANE-1'] ->
    instance.elements 의 connectivity, frame.fieldOutputs['COORD'].bulkDataBlocks 좌표,
    frameValue=0 프레임에서 구한 미변형(pristine) 기준면적 area_ref.
    quad: cross(v2-v0, v3-v1) / tri: cross(v1-v0, v2-v0), 면적 = 0.5*|cross|,
    법선 = cross/|cross|, n_z<0 이면 뒤집기(참조구현 기본과 동일).
    [주의] slope_stats.py 는 수정 금지 대상이라 import 대신 **복제**했다. 기하 판독이
      두 곳에 있으므로 한쪽을 고치면 다른 쪽도 함께 봐야 한다(의도적 중복).

메모리
    (요소 x 요소) 행렬을 만들지 않는다. 방향도 블록으로 나눈다(블록당 요소수 x 방향수 상한
    MEM_BUDGET). 10100 요소 x 4000 방향 x 3성분 float64 = 970 MB 를 한 번에 만들지 않는다.

사용 (code/ 에서)
    abaqus python torque_stats.py <odb> [--step NAME] [--frame N] [--n-mc 4000] [--json PATH]

재사용 API (§19.103 배선)
    compute_stats(odb_path, step=None, frame=None, n_mc=4000) -> dict | None
        run_abaqus.py 가 성공한 HF 런에서 **인프로세스로** 호출해 점별 results txt 에
        'tau_' 접두로 아카이빙한다(scalar_items). 서브프로세스 금지 - 이미 Abaqus 안이다.

종료코드
    0 정상 / 1 실패(ODB 열기/스텝/프레임/기하/JSON) / 2 사용법

한계 (정직하게)
    - 기하만 읽는다. 하중/재료/경계조건은 재실행하지 않는다.
    - 힘에 **변형 면적** A_e 를 쓴다(참조구현과 동일). 기준면적 area_ref 는 L_ref 정규화에만
      쓰므로, 면적이 1% 늘면 tau_rel 은 그 반대만큼 움직인다(면적비를 함께 찍는다).
    - tau_rel 에 문헌 허용치는 **없다**. 우리가 정의한 무차원 지표이고, 절대값이 아니라
      설계점 간 비교(x_c 스캔)로 읽어야 한다.
    - phi 격자는 0~60도다. 스칠 때(phi -> 90도) sum|f_e| -> 0 이라 tau_rel 의 분모가 0 으로
      가는데, 그 영역은 nan 으로 두고 개수를 찍는다(0 으로 나눠 inf 를 만들지 않는다).
    - 평판 frameValue=0 프레임이 없으면 tau_flat 은 nan 이다(없는 값을 0 으로 채우지 않는다).
    - 실제 ODB 로는 아직 돌려보지 않았다(로컬에 ODB 가 없다). 스텁 하네스는 배선/로직/
      해석앵커만 검증한다 - tau_rel 의 실제 크기는 실측 ODB 로만 판정된다.
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

# eval_abaqus.py 를 **고치지 않고** 상수/참조식을 그대로 가져온다.
try:
    from eval_abaqus import (A0, R0, IDEAL_THRUST_COEFF, U_SUN)             # noqa: F401
    from odbAccess import openOdb                                          # noqa: F401
except ImportError as _e:
    print("[torque] odbAccess/eval_abaqus 를 못 읽었습니다 (%s)" % _e)
    print("[torque] 'abaqus python torque_stats.py <odb>' 로 실행해야 합니다.")
    sys.exit(2)

_PFX = '[torque]'

# ------------------------------------------------------------------ 상수
INSTANCE_NAME = 'MEMBRANE-1'
DEFAULT_STEP_PREF = ('Step-Postbuckle',)   # 있으면 이것, 없으면 마지막 스텝
DEFAULT_N_MC = 4000                        # 디스크 MC 방향 샘플 수(기본)
TORQUE_SEED = 12345                        # eval_abaqus_avg.DEFAULT_SEED 와 같은 관례(재현 가능)
#   방향 격자: phi = 태양 방향의 법선 기준 기울기(도), psi = 방위각(도).
#   phi 상한 60도는 Heaton 2017 이 토크 비(2.07배)를 보고한 각도다. 90도 근처는 sum|f_e| -> 0
#   이라 tau_rel 분모가 죽는다(그 영역은 별도로 nan 개수를 찍는다).
PHI_GRID_DEG = (0.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 45.0, 60.0)
PSI_GRID_DEG = (0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0)
DEFAULT_BLOCK = 200                        # 블록당 방향 수 상한
MEM_BUDGET = 2.0e6                         # 블록당 (요소수 x 방향수) 상한(eval_abaqus_avg 관례)
#   tau_rel 분모(sum|f_e| * L_ref)의 하한. 이보다 작으면 tau_rel 은 정의되지 않는다 -> nan.
#   상대 하한(area_ref * IDEAL 대비)으로 잡아 메쉬 크기에 무관하게 만든다.
SUM_F_FLOOR_REL = 1e-12
STALE_WARN_S = 21600.0                     # 6시간 (§19.94 관례)

USAGE = """usage: abaqus python torque_stats.py <odb> [options]

  --step NAME    읽을 스텝 (기본: Step-Postbuckle 이 있으면 그것, 없으면 마지막 스텝)
  --frame N      스텝 안 프레임 번호, 1-based (기본: 마지막 프레임)
  --n-mc N       단위 디스크 MC 방향 샘플 수 (기본 4000)
  --json PATH    결과 JSON 경로 (기본 results/torque_<odb stem>_<step>.json)
  -h, --help     이 도움말
"""


# ------------------------------------------------------------------ 인자
def parse_args(argv):
    """수동 파싱. (args, err) 반환. err 가 있으면 사용법 오류(rc=2).

    argparse 를 쓰지 않는다: Abaqus 내장 python 이 구형일 수 있고, 이 저장소의 다른
    프로브(slope_stats.py/base_state_probe.py)도 sys.argv 를 직접 읽는 관례를 따른다.
    """
    a = {'odb': None, 'step': None, 'frame': None, 'n_mc': DEFAULT_N_MC,
         'json': None, 'help': False}
    i = 0
    while i < len(argv):
        t = argv[i]
        if t in ('-h', '--help'):
            a['help'] = True
            i += 1
            continue
        if t in ('--step', '--frame', '--n-mc', '--json'):
            if i + 1 >= len(argv):
                return None, '%s 뒤에 값이 필요하다' % t
            v = argv[i + 1]
            if t == '--step':
                a['step'] = v
            elif t == '--json':
                a['json'] = v
            elif t == '--n-mc':
                try:
                    n = int(v)
                except ValueError:
                    return None, '--n-mc 는 정수여야 한다: %r' % v
                if n < 1:
                    return None, '--n-mc 는 1 이상이어야 한다: %d' % n
                a['n_mc'] = n
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


# ------------------------------------------------------------------ 기하 (slope_stats.py 복제)
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
    중심은 요소 절점 좌표의 평균 - 면적가중 중심과 같지 않지만 참조구현/slope_stats 와 동일하다
    (면적은 요소별로 따로 쓰므로 모멘트 합에는 영향이 없다).
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


def _frame_zero(odb):
    """frameValue ~ 0 (미변형/평판) 프레임. eval_abaqus.get_hf 와 같은 탐색 순서."""
    try:
        for sn in odb.steps.keys():
            for fr in odb.steps[sn].frames:
                if abs(float(fr.frameValue)) < 1e-12:
                    return fr
    except Exception:
        pass
    return None


def _geom_of_frame(frame, instance, quads, tris, max_label):
    """프레임 하나의 (areas, normals, centers). COORD 가 없으면 None."""
    try:
        cf = frame.fieldOutputs['COORD'].getSubset(region=instance)
        coords = _node_coords(cf, max_label)
    except Exception as e:
        print(_PFX + " !!! COORD 판독 실패: %s" % e)
        return None
    return face_geometry(coords, quads, tris, True)


def read_torque_geometry(odb_path, step_arg, frame_arg):
    """ODB 를 한 번 열어 좌표/법선/면적/중심/area_ref/평판 프레임 기하를 뽑는다. 실패면 None."""
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

        if instance.nodes:
            max_label = instance.nodes[-1].label + 100
        else:
            max_label = 100000

        # ---- 미변형(pristine) 프레임: area_ref(=L_ref 정규화) + tau_flat 계산용 기하
        f0 = _frame_zero(odb)
        area_ref = None
        flat = None
        if f0 is None:
            print(_PFX + " (경고) frameValue=0 프레임이 없다 -> area_ref/tau_flat 을 못 낸다.")
        else:
            flat = _geom_of_frame(f0, instance, quads, tris, max_label)
            if flat is None:
                print(_PFX + " (경고) 미변형 프레임 기하 판독 실패 -> tau_flat 은 nan.")
            else:
                area_ref = float(np.sum(flat[0]))
                if area_ref <= 0.0:
                    area_ref, flat = None, None
                    print(_PFX + " (경고) 미변형 총면적이 0 이다 -> area_ref/tau_flat 은 nan.")

        # ---- 변형 프레임 기하
        g = _geom_of_frame(frame, instance, quads, tris, max_label)
        if g is None:
            return None
        areas, normals, centers = g
        areas_sum = float(np.sum(areas))
        if areas_sum <= 0.0:
            print(_PFX + " !!! 변형 프레임 총면적이 0 이다. 계산 불가.")
            return None

        if area_ref is None:
            area_ref = areas_sum
            print(_PFX + " (주의) area_ref 를 못 구해 변형 면적 %.6f m^2 로 대체한다(-> L_ref 도 변형 기준)."
                  % area_ref)
        else:
            print(_PFX + " 미변형(pristine) 기준면적 area_ref = %.6f m^2 (frameValue=0)" % area_ref)
        print(_PFX + " 기하: 요소 %d개 (quad %d / tri %d) | 변형 총면적 %.6f m^2 | 면적비 변형/미변형 %.6f"
              % (areas.size, len(quads), len(tris), areas_sum, areas_sum / area_ref))

        return {'step': sname, 'step_how': how, 'frame_index': fi, 'n_frames': nfr,
                'frame_value': fval, 'n_quad': len(quads), 'n_tri': len(tris),
                'labels': labels_q + labels_t, 'areas': areas, 'normals': normals,
                'centers': centers, 'area_ref': float(area_ref), 'area_sum': areas_sum,
                'flat': flat}
    finally:
        try:
            odb.close()
        except Exception:
            pass


# ------------------------------------------------------------------ 방향 샘플러
def dir_at_tilt(phi_deg, psi_deg=0.0):
    """고정 기울기 방향 1개. s_hat = -(sin(p)cos(t), sin(p)sin(t), cos(p)).

    eval_abaqus_avg.dir_at_tilt 와 **같은 식**이다(수정 금지 대상이라 복제). s_hat 은
    광자의 진행방향이고 참조구현의 U_SUN 과 같은 역할이다(phi=0 -> U_SUN).
    """
    p = math.radians(phi_deg)
    t = math.radians(psi_deg)
    return np.array([-(math.sin(p) * math.cos(t)), -(math.sin(p) * math.sin(t)),
                     -math.cos(p)], dtype=np.float64)


def grid_dirs():
    """phi x psi 격자 방향. (s_dirs(K,3), phi(K), psi(K)) 반환."""
    phis, psis = [], []
    for p in PHI_GRID_DEG:
        for t in PSI_GRID_DEG:
            phis.append(float(p))
            psis.append(float(t))
    d = np.vstack([dir_at_tilt(p, t) for p, t in zip(phis, psis)])
    return d, np.asarray(phis, dtype=np.float64), np.asarray(psis, dtype=np.float64)


def sample_dirs_disk(n, seed=TORQUE_SEED):
    """단위 디스크 균일 (x,y) -> 밀도가 cos(phi) 에 비례하는 **전력가중** 샘플.

    eval_abaqus_avg.sample_dirs('disk', ...) 와 **같은 식**이다(수정 금지 대상이라 복제).
    """
    rng = np.random.RandomState(seed)
    r = np.sqrt(rng.rand(int(n)))
    psi = 2.0 * math.pi * rng.rand(int(n))
    x = r * np.cos(psi)
    y = r * np.sin(psi)
    z = np.sqrt(np.maximum(1.0 - x * x - y * y, 0.0))
    return -np.stack((x, y, z), axis=1)


def dirs_to_phi_psi(s_dirs):
    """s_hat -> (phi, psi) 도. dir_at_tilt 의 역변환."""
    s = np.asarray(s_dirs, dtype=np.float64)
    cz = np.clip(-s[:, 2], -1.0, 1.0)
    phi = np.degrees(np.arccos(cz))
    psi = np.degrees(np.mod(np.arctan2(-s[:, 1], -s[:, 0]), 2.0 * math.pi))
    return phi, psi


# ------------------------------------------------------------------ 모멘트
def _block_size(n_elem, n_dirs, n_block=DEFAULT_BLOCK):
    """블록당 방향 수 = min(n_block, MEM_BUDGET/요소수). (요소x방향) 행렬을 한 번에 만들지 않는다."""
    return max(1, min(int(n_block), int(MEM_BUDGET / float(max(int(n_elem), 1))), int(n_dirs)))


def moments_over_dirs(areas, normals, centers, s_dirs, x_ref, n_block=DEFAULT_BLOCK):
    """방향 K개에 대해 M(s_hat) = sum_e (x_e - x_ref) x f_e 와 sum_e |f_e| 를 낸다.

    f_e = a_e*(A0*c_e*s_hat + 2*R0*c_e^2*(-n_e)),  c_e = clip(n_e . (-s_hat), 0, 1)
    참조구현의 f_abs/f_refl 과 같은 식이고, U_SUN 자리에 방향 s_hat 을 넣은 것뿐이다.

    반환: dict(M=(K,3), absM=(K,), sum_f=(K,), n_nan=int, block=int)
    """
    a = np.asarray(areas, dtype=np.float64)
    n = np.asarray(normals, dtype=np.float64)
    ce = np.asarray(centers, dtype=np.float64)
    s = np.asarray(s_dirs, dtype=np.float64)
    xr = np.asarray(x_ref, dtype=np.float64).reshape(3)
    K = int(s.shape[0])
    M = np.zeros((K, 3), dtype=np.float64)
    SF = np.zeros(K, dtype=np.float64)
    if a.size == 0 or K == 0:
        return {'M': M, 'absM': np.zeros(K), 'sum_f': SF, 'n_nan': K,
                'block': 0, 'n_elem': int(a.size)}
    blk = _block_size(a.size, K, n_block)
    nx, ny, nz = n[:, 0], n[:, 1], n[:, 2]
    rx, ry, rz = ce[:, 0] - xr[0], ce[:, 1] - xr[1], ce[:, 2] - xr[2]
    for i in range(0, K, blk):
        j = min(i + blk, K)
        sb = s[i:j]
        c = np.clip(-np.dot(n, sb.T), 0.0, 1.0)              # (n_elem, K_blk)
        c2 = c * c
        #   f = a*(A0*c*s_hat - 2*R0*c^2*n)  (성분별로 만든다 - (n,K,3) 을 한 번에 만들지 않는다)
        fx = a[:, None] * (A0 * c * sb[None, :, 0] - 2.0 * R0 * c2 * nx[:, None])
        fy = a[:, None] * (A0 * c * sb[None, :, 1] - 2.0 * R0 * c2 * ny[:, None])
        fz = a[:, None] * (A0 * c * sb[None, :, 2] - 2.0 * R0 * c2 * nz[:, None])
        #   r x f = (ry*fz - rz*fy, rz*fx - rx*fz, rx*fy - ry*fx).  ry @ fz 는 (K,) (BLAS gemv).
        M[i:j, 0] = ry @ fz - rz @ fy
        M[i:j, 1] = rz @ fx - rx @ fz
        M[i:j, 2] = rx @ fy - ry @ fx
        SF[i:j] = np.sum(np.sqrt(fx * fx + fy * fy + fz * fz), axis=0)
    absM = np.sqrt(M[:, 0] ** 2 + M[:, 1] ** 2 + M[:, 2] ** 2)
    return {'M': M, 'absM': absM, 'sum_f': SF, 'n_nan': 0, 'block': int(blk),
            'n_elem': int(a.size)}


def tau_rel_of(absM, sum_f, l_ref, floor_abs):
    """tau_rel = |M| / (sum|f| * L_ref). 분모가 하한 이하면 nan(0 으로 나누지 않는다)."""
    absM = np.asarray(absM, dtype=np.float64)
    den = np.asarray(sum_f, dtype=np.float64) * float(l_ref)
    out = np.full(absM.shape, float('nan'), dtype=np.float64)
    ok = np.isfinite(den) & (den > float(floor_abs))
    out[ok] = absM[ok] / den[ok]
    return out, den


def dist_stats(v):
    """nan 제거 후 min/p50/p95/max/평균/SE/표준편차. (eval_abaqus_avg.stats_of 관례)"""
    v = np.asarray(v, dtype=np.float64)
    good = np.isfinite(v)
    w = v[good]
    out = {'n': int(w.size), 'n_total': int(v.size), 'n_nan': int(v.size) - int(w.size),
           'mean': float('nan'), 'se': float('nan'), 'std': float('nan'),
           'min': float('nan'), 'p50': float('nan'), 'p95': float('nan'),
           'max': float('nan')}
    if w.size == 0:
        return out
    out['mean'] = float(np.mean(w))
    sd = float(np.std(w, ddof=1)) if w.size > 1 else 0.0
    out['std'] = sd
    out['se'] = sd / math.sqrt(w.size)
    out['min'] = float(np.min(w))
    out['max'] = float(np.max(w))
    q = np.percentile(w, [50.0, 95.0])
    out['p50'], out['p95'] = float(q[0]), float(q[1])
    return out


def _argmax_finite(v):
    """nan 을 무시한 argmax. 유한값이 없으면 -1."""
    v = np.asarray(v, dtype=np.float64)
    if not np.any(np.isfinite(v)):
        return -1
    return int(np.nanargmax(v))


# ------------------------------------------------------------------ 재사용 API (§19.103)
def compute_stats(odb_path, step=None, frame=None, n_mc=DEFAULT_N_MC):
    """HF ODB -> 토크 지표 dict. 실패하면 None (§19.103).

    main() 의 계산부를 그대로 감싼 **재사용** 함수다. run_abaqus.py 가 성공한 HF 런에서
    **인프로세스로** 호출한다(서브프로세스 금지 - 이미 Abaqus 안이다). [torque] 접두 로그가
    기하/모멘트/방향분포 몇 줄 나간다.

    반환 dict 의 **스칼라 키**는 run_abaqus.py 가 점별 results txt 에 'tau_' 접두로 아카이빙한다
    (scalar_items). 중첩 키(sweep/flat/moment_*/grid/disk)는 JSON 전용이고 scalar_items 는 건너뛴다.
    """
    g = read_torque_geometry(odb_path, step, frame)
    if g is None:
        return None

    areas = g['areas']
    normals = g['normals']
    centers = g['centers']
    valid = np.isfinite(areas) & (areas > 0.0)
    if not np.any(valid):
        print(_PFX + " !!! 유효 요소(면적>0)가 0개다.")
        return None
    n_degen = int(areas.size) - int(np.count_nonzero(valid))
    if n_degen:
        print(_PFX + " (주의) 퇴화 요소 %d개 제외(면적<=0)." % n_degen)
    v_areas = areas[valid]
    v_norm = normals[valid]
    v_ctr = centers[valid]
    sum_a = float(np.sum(v_areas))
    area_ref = float(g['area_ref'])
    l_ref = math.sqrt(area_ref)
    floor_abs = SUM_F_FLOOR_REL * area_ref * IDEAL_THRUST_COEFF * l_ref

    # ---- (1) 기준점 두 개: 면적가중 도심(변형 좌표) + 원점
    x_cen = np.array([float(np.sum(v_areas * v_ctr[:, 0]) / sum_a),
                      float(np.sum(v_areas * v_ctr[:, 1]) / sum_a),
                      float(np.sum(v_areas * v_ctr[:, 2]) / sum_a)], dtype=np.float64)
    x_org = np.zeros(3, dtype=np.float64)
    d_cen = float(np.linalg.norm(x_cen))
    print(_PFX + " 기준점: 면적가중 도심 x_ref=(%.6f, %.6f, %.6f) m (원점에서 %.6f m) / 원점 (0,0,0)"
          % (x_cen[0], x_cen[1], x_cen[2], d_cen))
    print(_PFX + " L_ref = sqrt(area_ref=%.6f) = %.6f m | sum|f_e| 하한 = %.3e N"
          % (area_ref, l_ref, floor_abs))

    # ---- (2) 수직입사(s_hat = U_SUN) 모멘트: 도심 / 원점
    s_plane = np.asarray(U_SUN, dtype=np.float64).reshape(1, 3)
    r_c = moments_over_dirs(v_areas, v_norm, v_ctr, s_plane, x_cen)
    r_o = moments_over_dirs(v_areas, v_norm, v_ctr, s_plane, x_org)
    tau_plane_arr, den_plane = tau_rel_of(r_c['absM'], r_c['sum_f'], l_ref, floor_abs)
    tau_org_arr, _den_o = tau_rel_of(r_o['absM'], r_o['sum_f'], l_ref, floor_abs)
    tau_abs = float(r_c['absM'][0])
    tau_rel_plane = float(tau_plane_arr[0])
    tau_abs_origin = float(r_o['absM'][0])
    tau_rel_origin_plane = float(tau_org_arr[0])
    mx, my, mz = float(r_c['M'][0, 0]), float(r_c['M'][0, 1]), float(r_c['M'][0, 2])
    m_out = math.sqrt(mx * mx + my * my)
    m_z_abs = abs(mz)
    out_share = (m_out / tau_abs) if tau_abs > 0.0 else float('nan')
    in_share = (m_z_abs / tau_abs) if tau_abs > 0.0 else float('nan')
    out_over_in = (m_out / m_z_abs) if m_z_abs > 0.0 else float('nan')
    sum_f_plane = float(r_c['sum_f'][0])

    print(_PFX + " 수직입사(도심기준): |M|=%.10e N m | M=(%.6e, %.6e, %.6e) | sum|f_e|=%.6e N"
          % (tau_abs, mx, my, mz, sum_f_plane))
    print(_PFX + " 수직입사(도심기준): tau_rel = |M|/(sum|f|*L_ref) = %.10e"
          % (tau_rel_plane if tau_rel_plane == tau_rel_plane else float('nan')))
    print(_PFX + " 수직입사(원점기준): |M|=%.10e N m | tau_rel = %.10e  (평판이면 지렛대 항만 남는다)"
          % (tau_abs_origin,
             tau_rel_origin_plane if tau_rel_origin_plane == tau_rel_origin_plane else float('nan')))
    print(_PFX + " 면외(M_x,M_y)=%.6e N m (비중 %.4f) | 면내 M_z=%.6e N m (비중 %.4f) | 면외/면내=%.4f"
          % (m_out, out_share, mz, in_share, out_over_in))
    print(_PFX + " (출처: Zeiders 2005 NTRS 20050209955 - RMS slope 는 모멘트를 보장하지 않음 / "
          "Heaton 2017 NTRS 20170001488 - 0도 추력 동일(722 uN), 60도 토크 2.07배)")

    # ---- (3) 평판(미변형) 프레임 sanity anchor
    flat_abs = flat_rel = flat_rel_org = float('nan')
    flat_ok = False
    if g['flat'] is not None:
        a0, n0, c0 = g['flat']
        v0 = np.isfinite(a0) & (a0 > 0.0)
        if np.any(v0):
            s0 = float(np.sum(a0[v0]))
            xc0 = np.array([float(np.sum(a0[v0] * c0[v0, 0]) / s0),
                            float(np.sum(a0[v0] * c0[v0, 1]) / s0),
                            float(np.sum(a0[v0] * c0[v0, 2]) / s0)], dtype=np.float64)
            rf = moments_over_dirs(a0[v0], n0[v0], c0[v0], s_plane, xc0)
            tau_f, _d = tau_rel_of(rf['absM'], rf['sum_f'], l_ref, floor_abs)
            flat_abs = float(rf['absM'][0])
            flat_rel = float(tau_f[0])
            ro = moments_over_dirs(a0[v0], n0[v0], c0[v0], s_plane, x_org)
            tau_fo, _d2 = tau_rel_of(ro['absM'], ro['sum_f'], l_ref, floor_abs)
            flat_rel_org = float(tau_fo[0])
            flat_ok = True
    print(_PFX + " 평판(미변형 frame 0): tau_flat_rel(도심)=%.10e | tau_flat_rel(원점)=%.10e"
          % (flat_rel if flat_rel == flat_rel else float('nan'),
             flat_rel_org if flat_rel_org == flat_rel_org else float('nan')))
    if flat_ok:
        print(_PFX + "   (앵커: 평판 + 대칭 메쉬면 도심 기준 M ~ 0 이어야 한다. 원점 기준은 "
              "돛이 원점에 대해 비대칭(삼각형)이면 지렛대 항 때문에 0 이 아니다.)")
    else:
        print(_PFX + "   (경고: 미변형 프레임을 못 읽어 tau_flat 은 nan - 0 으로 채우지 않는다.)")

    # ---- (4) 방향 격자
    g_dirs, g_phi, g_psi = grid_dirs()
    rg = moments_over_dirs(v_areas, v_norm, v_ctr, g_dirs, x_cen)
    tau_g, _dg = tau_rel_of(rg['absM'], rg['sum_f'], l_ref, floor_abs)
    sg = dist_stats(tau_g)
    ig = _argmax_finite(tau_g)
    g_worst_phi = float(g_phi[ig]) if ig >= 0 else float('nan')
    g_worst_psi = float(g_psi[ig]) if ig >= 0 else float('nan')
    g_worst_rel = float(tau_g[ig]) if ig >= 0 else float('nan')
    print(_PFX + " 방향 격자(%d phi x %d psi = %d방향, %d블록): tau_rel min=%.6e p50=%.6e p95=%.6e max=%.6e"
          " (nan %d)"
          % (len(PHI_GRID_DEG), len(PSI_GRID_DEG), g_dirs.shape[0], rg['block'],
             sg['min'], sg['p50'], sg['p95'], sg['max'], sg['n_nan']))
    print(_PFX + "   최악 격자 방향: phi=%.1f도 psi=%.1f도 tau_rel=%.10e"
          % (g_worst_phi, g_worst_psi, g_worst_rel))

    # ---- (5) 디스크 MC(전력가중)
    m_dirs = sample_dirs_disk(int(n_mc), seed=TORQUE_SEED)
    rm = moments_over_dirs(v_areas, v_norm, v_ctr, m_dirs, x_cen)
    tau_m, _dm = tau_rel_of(rm['absM'], rm['sum_f'], l_ref, floor_abs)
    sm = dist_stats(tau_m)
    im = _argmax_finite(tau_m)
    m_phi, m_psi = dirs_to_phi_psi(m_dirs)
    m_worst_phi = float(m_phi[im]) if im >= 0 else float('nan')
    m_worst_psi = float(m_psi[im]) if im >= 0 else float('nan')
    m_worst_rel = float(tau_m[im]) if im >= 0 else float('nan')
    print(_PFX + " 디스크 MC(전력가중, seed=%d, n=%d): tau_rel mean=%.6e SE=%.3e std=%.6e min=%.6e "
          "p50=%.6e p95=%.6e max=%.6e (nan %d)"
          % (TORQUE_SEED, int(n_mc), sm['mean'], sm['se'], sm['std'], sm['min'],
             sm['p50'], sm['p95'], sm['max'], sm['n_nan']))
    print(_PFX + "   최악 MC 방향: phi=%.1f도 psi=%.1f도 tau_rel=%.10e" % (m_worst_phi, m_worst_psi, m_worst_rel))

    # ---- (6) 전체 최악(격자 vs MC)
    if ig < 0 and im < 0:
        tau_rel_max, worst_phi, worst_psi, worst_src = float('nan'), float('nan'), float('nan'), 'none'
    elif im < 0 or (ig >= 0 and g_worst_rel >= m_worst_rel):
        tau_rel_max, worst_phi, worst_psi, worst_src = g_worst_rel, g_worst_phi, g_worst_psi, 'grid'
    else:
        tau_rel_max, worst_phi, worst_psi, worst_src = m_worst_rel, m_worst_phi, m_worst_psi, 'disk_mc'
    print(_PFX + " 최악 방향(전체): phi=%.1f도 psi=%.1f도 tau_rel_max=%.10e (출처 %s)"
          % (worst_phi, worst_psi, tau_rel_max, worst_src))

    return {
        # ---- 스텝/프레임 신원
        'step': g['step'], 'step_how': g['step_how'], 'frame': g['frame_index'],
        'n_frames': g['n_frames'], 'frame_value': g['frame_value'],
        # ---- 기하
        'n_elem': int(areas.size), 'n_quad': g['n_quad'], 'n_tri': g['n_tri'],
        'n_valid': int(np.count_nonzero(valid)), 'n_degenerate_removed': n_degen,
        'area_ref': area_ref, 'area_sum': sum_a, 'l_ref': l_ref,
        'area_ratio': sum_a / area_ref,
        # ---- 기준점
        'x_ref_x': float(x_cen[0]), 'x_ref_y': float(x_cen[1]), 'x_ref_z': float(x_cen[2]),
        'x_ref_from_origin': d_cen,
        # ---- 수직입사 모멘트(도심) + 성분
        'tau_abs': tau_abs, 'tau_rel_plane': tau_rel_plane,
        'tau_mx': mx, 'tau_my': my, 'tau_mz': mz,
        'tau_mout': m_out, 'tau_out_share': out_share, 'tau_in_share': in_share,
        'tau_out_over_in': out_over_in, 'sum_f_plane': sum_f_plane,
        # ---- 수직입사 모멘트(원점)
        'tau_abs_origin': tau_abs_origin, 'tau_rel_origin_plane': tau_rel_origin_plane,
        # ---- 평판 앵커
        'tau_flat_abs': flat_abs, 'tau_flat_rel': flat_rel,
        'tau_flat_rel_origin': flat_rel_org, 'tau_flat_ok': bool(flat_ok),
        # ---- 방향 격자
        'tau_rel_grid_mean': sg['mean'], 'tau_rel_grid_se': sg['se'], 'tau_rel_grid_min': sg['min'],
        'tau_rel_grid_p50': sg['p50'], 'tau_rel_grid_p95': sg['p95'],
        'tau_rel_grid_max': sg['max'], 'tau_grid_n': sg['n'], 'tau_grid_n_nan': sg['n_nan'],
        'tau_grid_worst_phi_deg': g_worst_phi, 'tau_grid_worst_psi_deg': g_worst_psi,
        # ---- 디스크 MC(전력가중)
        'tau_rel_cosw': sm['mean'], 'tau_rel_cosw_se': sm['se'], 'tau_rel_cosw_std': sm['std'],
        'tau_rel_cosw_min': sm['min'], 'tau_rel_cosw_p50': sm['p50'],
        'tau_rel_cosw_p95': sm['p95'], 'tau_rel_cosw_max': sm['max'],
        'tau_disk_n': sm['n'], 'tau_disk_n_nan': sm['n_nan'],
        'tau_disk_worst_phi_deg': m_worst_phi, 'tau_disk_worst_psi_deg': m_worst_psi,
        'tau_disk_worst_rel': m_worst_rel,
        # ---- 전체 최악 / 샘플러 설정
        'tau_rel_max': tau_rel_max, 'tau_worst_phi_deg': worst_phi,
        'tau_worst_psi_deg': worst_psi, 'tau_worst_source': worst_src,
        'tau_rel_p95': sg['p95'], 'tau_rel_mean': sm['mean'],
        'tau_n_mc': int(n_mc), 'tau_seed': int(TORQUE_SEED),
        'tau_n_phi': len(PHI_GRID_DEG), 'tau_n_psi': len(PSI_GRID_DEG),
        'tau_n_dirs': int(g_dirs.shape[0]), 'tau_block': int(rg['block']),
        'tau_sum_f_floor_n': float(floor_abs),
        # ---- 중첩/배열 (출력 + JSON 재사용; scalar_items 는 건너뛴다)
        'grid': {'phi_deg': [float(x) for x in PHI_GRID_DEG],
                 'psi_deg': [float(x) for x in PSI_GRID_DEG],
                 'n_dirs': int(g_dirs.shape[0]),
                 'tau_rel': [float(x) for x in tau_g],
                 'tau_rel_phi_deg': [float(x) for x in g_phi],
                 'tau_rel_psi_deg': [float(x) for x in g_psi],
                 'stats': dict(sg),
                 'worst': {'phi_deg': g_worst_phi, 'psi_deg': g_worst_psi,
                           'tau_rel': g_worst_rel}},
        'disk': {'n': int(n_mc), 'seed': int(TORQUE_SEED), 'stats': dict(sm),
                 'worst': {'phi_deg': m_worst_phi, 'psi_deg': m_worst_psi,
                           'tau_rel': m_worst_rel}},
        'moment_plane_centroid': {'x_ref': [float(x) for x in x_cen],
                                  'M': [mx, my, mz], 'absM': tau_abs,
                                  'sum_f': sum_f_plane, 'l_ref': l_ref,
                                  'tau_rel': tau_rel_plane},
        'moment_plane_origin': {'x_ref': [0.0, 0.0, 0.0],
                                'M': [float(r_o['M'][0, 0]), float(r_o['M'][0, 1]),
                                      float(r_o['M'][0, 2])],
                                'absM': tau_abs_origin, 'sum_f': float(r_o['sum_f'][0]),
                                'tau_rel': tau_rel_origin_plane},
        'flat': {'ok': bool(flat_ok), 'absM': flat_abs, 'tau_rel': flat_rel,
                 'tau_rel_origin': flat_rel_org},
        'note': ('f_e = A_e*(A0*c_e*s_hat + 2*R0*c_e^2*(-n_e)), c_e = clip(n_e.(-s_hat),0,1); '
                 'M = sum_e (x_e - x_ref) x f_e. tau_rel = |M| / (sum|f_e| * L_ref), '
                 'L_ref = sqrt(area_ref). x_ref = area-weighted centroid (deformed coords) '
                 'and origin, both reported. Literature: Zeiders 2005 (NTRS 20050209955) RMS '
                 'slope does not guarantee moments; Heaton 2017 (NTRS 20170001488) shape leaves '
                 'thrust unchanged (722 uN) but changes torque 2.07x at 60 deg.'),
    }


def default_json_path(odb_path, step):
    """torque_stats JSON 기본 경로: results/torque_<odb stem>_<step>.json (기존 관례)."""
    stem = os.path.splitext(os.path.basename(odb_path))[0]
    return os.path.join('results', 'torque_%s_%s.json' % (stem, step))


def _json_safe(obj):
    """NaN/Inf 를 None 으로 바꾼다(slope_stats.py 와 같은 규칙).

    JSON 표준에는 NaN/Infinity 가 없다(json.dumps 는 비표준 'NaN' 토큰을 쓴다). 엄격한
    파서가 그 파일에서 죽는다. 읽지 못한 값/빈 통계는 null 로 남긴다(0 으로 채우지 않는다).
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


def write_json_record(rec, out_path):
    """JSON 레코드 -> 파일. 디렉터리 자동 생성. 실패면 예외를 올린다(호출자가 처리)."""
    d = os.path.dirname(os.path.abspath(out_path))
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with io.open(out_path, 'w', encoding='utf-8') as f:
        f.write(json.dumps(_json_safe(rec), ensure_ascii=True, indent=2, sort_keys=True))
    return os.path.abspath(out_path)


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


def build_json_record(S, odb_path, ident=None, t_start=None):
    """compute_stats 결과(S) -> torque_stats JSON 레코드(분포/성분/샘플러 설정 포함)."""
    if ident is None:
        ident = _odb_ident_quiet(odb_path)
    if t_start is None:
        t_start = time.time()
    return {
        'script': 'torque_stats.py', 'section': '19.103',
        'generated': time.strftime('%Y-%m-%d %H:%M:%S'),
        'probe_run_epoch': t_start, 'odb': ident,
        'step': S['step'], 'step_how': S['step_how'], 'frame_index': S['frame'],
        'n_frames': S['n_frames'], 'frame_value': S['frame_value'],
        'n_elem': S['n_elem'], 'n_quad': S['n_quad'], 'n_tri': S['n_tri'],
        'n_valid': S['n_valid'], 'n_degenerate_removed': S['n_degenerate_removed'],
        'area_ref': S['area_ref'], 'area_deformed_sum': S['area_sum'],
        'area_ratio_deformed_pristine': S['area_ratio'], 'l_ref': S['l_ref'],
        'x_ref_centroid': [S['x_ref_x'], S['x_ref_y'], S['x_ref_z']],
        'x_ref_from_origin': S['x_ref_from_origin'],
        'tau': {'abs_plane_centroid': S['tau_abs'], 'rel_plane_centroid': S['tau_rel_plane'],
                'abs_plane_origin': S['tau_abs_origin'],
                'rel_plane_origin': S['tau_rel_origin_plane'],
                'mx': S['tau_mx'], 'my': S['tau_my'], 'mz': S['tau_mz'],
                'm_out_of_plane': S['tau_mout'], 'out_share': S['tau_out_share'],
                'in_share': S['tau_in_share'], 'out_over_in': S['tau_out_over_in'],
                'sum_f_plane': S['sum_f_plane'], 'l_ref': S['l_ref'],
                'rel_max_all_dirs': S['tau_rel_max'],
                'worst_phi_deg': S['tau_worst_phi_deg'],
                'worst_psi_deg': S['tau_worst_psi_deg'],
                'worst_source': S['tau_worst_source'],
                'rel_p95_grid': S['tau_rel_p95'], 'rel_mean_disk': S['tau_rel_mean']},
        'flat_anchor': S['flat'],
        'sweep': {'grid': S['grid'], 'disk_mc': S['disk'],
                  'sampler': {'phi_grid_deg': S['tau_n_phi'], 'psi_grid_deg': S['tau_n_psi'],
                              'n_dirs': S['tau_n_dirs'], 'n_mc': S['tau_n_mc'],
                              'seed': S['tau_seed'], 'block': S['tau_block'],
                              'sum_f_floor_n': S['tau_sum_f_floor_n'],
                              'disk_mode': 'unit disk uniform (x,y) -> cos(phi) power weighted',
                              'grid_phi_max_deg': max(PHI_GRID_DEG)}},
        'moment_plane_centroid': S['moment_plane_centroid'],
        'moment_plane_origin': S['moment_plane_origin'],
        'note': S['note'],
    }


_SCALAR_SKIP = ('note',)


def scalar_items(S, prefix='tau_', skip=_SCALAR_SKIP):
    """dict -> 'prefix+key=value' 줄 목록(키 정렬). **스칼라만**(int/float/str/bool/np 스칼라).

    접두는 키마다 **한 번만** 붙인다 - 이 스크립트의 dict 키는 이미 'tau_' 를 달고 있으므로
    (예: tau_rel_plane) 두 번 붙이면 tau_tau_... 가 된다(§19.102 의 slope_stats 는 키가 접두
    없는 이름이라 항상 붙였지만, 여기서는 이름 자체가 접두를 포함한다).
    NaN/Inf -> 'nan'. 중첩(dict/list/tuple) 키와 skip 키는 건너뛴다(§19.103 점별 txt 아카이빙용).
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
        kk = k if (prefix and k.startswith(prefix)) else (prefix + k)
        out.append('%s=%s' % (kk, s))
    return out


# ------------------------------------------------------------------ main
def main(argv=None):
    t_start = time.time()
    if argv is None:
        argv = sys.argv[1:]
    print(_PFX + " build: torque_stats.py (19.103) | numpy %s | U_SUN=%s | A0=%.6f R0=%.6f IDEAL=%.6f"
          % (np.__version__, np.asarray(U_SUN, dtype=float).tolist(),
             A0, R0, IDEAL_THRUST_COEFF))
    print(_PFX + " 정의: f_e=A_e*(A0*c_e*s_hat + 2*R0*c_e^2*(-n_e)), c_e=clip(n_e.(-s_hat),0,1);"
          " M=sum (x_e-x_ref) x f_e; tau_rel=|M|/(sum|f_e|*L_ref), L_ref=sqrt(area_ref)")
    print(_PFX + " 문헌: Zeiders 2005(NTRS 20050209955) RMS slope 는 모멘트를 보장하지 않음 | "
          "Heaton 2017(NTRS 20170001488) 0도 추력 동일(722 uN) - 60도 토크 2.07배")
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
    S = compute_stats(args['odb'], args['step'], args['frame'], n_mc=args['n_mc'])
    if S is None:
        return 1

    out_path = args['json']
    if out_path is None:
        out_path = default_json_path(args['odb'], S['step'])
    rec = build_json_record(S, args['odb'], ident=ident, t_start=t_start)
    try:
        print(_PFX + " JSON 저장: %s" % write_json_record(rec, out_path))
    except Exception as e:
        print(_PFX + " !!! JSON 저장 실패: %s" % e)
        return 1

    # 기계 파싱용 KV 1줄 (점별 txt 아카이빙과 **같은 키 집합** + 접두 없는 step/frame)
    kv = {}
    for _it in scalar_items(S):
        _k, _v = _it.split('=', 1)
        kv[_k] = _v
    kv['step'] = S['step']
    kv['frame'] = S['frame']
    print('[torque:KV] ' + ' '.join('%s=%s' % (k, kv[k]) for k in sorted(kv)))
    print(_PFX + " 완료 (rc=0)")
    return 0


if __name__ == '__main__':
    sys.exit(main())
