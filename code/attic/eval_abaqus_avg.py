# -*- coding: utf-8 -*-
"""eval_abaqus_avg.py - 태양돛 HF ODB 의 **반구 전방향 평균** 추력손실 (§19.99)

eval_abaqus.py 는 **수정하지 않는다**. 상수(A0/R0/IDEAL_THRUST_COEFF/U_SUN)와 참조구현
(calc_thrust_loss/get_hf)을 그대로 import 해서 쓴다.

무엇을 하는가
    기존 eval_abaqus.calc_thrust_loss 는 태양이 수직입사(U_SUN = (0,0,-1))한다고 가정하고
    변형된 돛의 추력손실 hf 를 **연직 1방향**에 대해 낸다. 이 스크립트는 **FE 재실행 없이**
    같은 ODB 기하를 한 번 읽고, 태양 방향(s_hat)만 반구 전방향에 대해 몬테카를로로 샘플링해
    방향평균 추력손실을 낸다. (방향마다 요소 수만큼의 내적/거듭제곱만 다시 계산하면 된다.)

이상(분모) 정의 두 가지 - 반드시 구분해서 읽어라 (--ideal)
    ideal='total'   : 분모 = area_ref * IDEAL_THRUST_COEFF = area_ref * (A0 + 2*R0)
                      '수직입사 완전평판'을 기준으로 한 손실이다. 방향 s_hat 이 기울면
                      완전평판도 이미 손실을 가지므로, 이 값은 **기울기 자체의 손실까지 포함**한다.
                      따라서 문헌(Murphey 2002: 5도 기울기에서 1.108%)과 직접 비교 가능하다.
    ideal='wrinkle' : 분모 = area_ref * (A0*cosf + 2*R0*cosf^3),  cosf = clip(-s_hat.z, 0, 1)
                      **같은 방향의 완전평판**을 기준으로 한 손실이다. 기울기 효과가 분자와
                      분모에서 같은 형태로 상쇄되므로 '주름 때문에 추가로 잃은 몫'만 격리된다.
                      완전평판 돛이면 어떤 방향에서도 정확히 0 이다.
    s_hat = U_SUN (수직입사)이면 두 정의가 모두 eval_abaqus.calc_thrust_loss 와 정확히 같아야
    하고, 그 사실을 실행 중 assert 로 확인한다(구현 드리프트 차단). assert 하나라도 깨지면 rc=1.

방향 샘플러 (둘 다 출력한다. --n-mc, --seed)
    mode='hemi' : 반구 위 **입체각 균일** 샘플. s_hat = -(sin(p)*cos(t), sin(p)*sin(t), cos(p)),
                  cos(p) ~ U[0,1) 이므로 입체각 균일. 단순 균일평균.
    mode='disk' : **전력(cos) 가중** 평균. 단위원 (x,y) 균일 -> s_hat = -(x, y, sqrt(1-x^2-y^2)).
                  복사압은 cos(입사각)에 비례해 실리므로 물리적으로 맞는 가중이다.
                  (단위원 균일밀도 = 입체각밀도 x cos(p). 실제로 그렇다.)
    두 샘플러의 평균은 서로 다르다(기울기 가중이 다르므로). 둘 다 출력해 차이를 보게 한다.

읽는 기하 (eval_abaqus.py 와 **동일한 방식** - get_hf()/calc_thrust_loss() 참조)
    odbAccess.openOdb(readOnly) -> odb.rootAssembly.instances['MEMBRANE-1'] ->
    instance.elements 의 connectivity, frame.fieldOutputs['COORD'].bulkDataBlocks 좌표,
    frameValue=0 프레임에서 구한 미변형(pristine) 기준면적 area_ref.
    사각형은 cross(v2-v0, v3-v1), 삼각형은 cross(v1-v0, v2-v0), 면적 = 0.5*|cross|,
    법선 = cross/|cross|, n_z<0 이면 뒤집기(= --flip 기본). 전부 참조구현과 같은 식이다.

메모리
    (요소수 x 방향수) 행렬을 한 번에 만들지 않는다. 방향을 블록으로 나눠 처리하고
    블록 크기 = min(--block 기본 500, MEM_BUDGET/요소수) 로 상한을 둔다.

출력
    기본: hf_ref(수직입사 참조), hf_normal_self(자기 기하 수직입사), MC 4조합 통계
    (평균, SE=std/sqrt(N), 표준편차, min, p1, p50, p99, max), 경사 스캔표(0/2/5/10/20도),
    ODB 신원(절대경로, 크기, mtime, 프로브 실행시각, 나이).

사용 (code/ 에서)
    abaqus python eval_abaqus_avg.py <odb> [--step NAME] [--frame N] [--n-mc 2000]
                                     [--seed 12345] [--ideal both|total|wrinkle]
                                     [--json PATH] [--flip|--no-flip] [--check]

종료코드
    0 정상 / 1 실패(ODB 열기/스텝/기하/자기검증/JSON) / 2 사용법

한계 (정직하게)
    - **기하만** 다시 쓴다. 하중/재료/경계조건은 재실행하지 않는다. 즉 '이 변형형상이 반구
      전방향 태양광을 받으면 평균적으로 얼마를 잃는가'를 묻는 것이지, '다른 태양 방향에 대해
      다시 해석하면 형상이 어떻게 달라지는가'를 묻는 것이 아니다.
    - 'wrinkle' 정의는 분모가 같은-방향 평판이라 매끄럽다. 그러나 완전평판 기준은 방향이
      스칠 때(s_hat.z -> 0-) 0 에 가까워지므로, 뒤집힌 면이 많은 실제 돛에서는 비율이 커질 수
      있다. 그래서 p99/max 를 함께 찍는다(꼬리를 숨기지 않는다).
    - hemi(입체각 균일)는 스칠 때 1/cos 형태의 꼬리를 가질 수 있어 평균 수렴이 disk 보다 느리다.
      SE 를 함께 찍고 N 을 늘려 확인하라. 물리적으로 맞는 가중은 disk(전력가중)다.
    - hf 는 참조구현과 같이 `max(loss, 0)` 으로 자른다(§19.82 정의). 자른 개수를 clip 으로 찍는다.
"""
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
#   get_hf 는 상태확인용으로 import 만 한다(여기서는 프레임/스텝을 인자로 받으므로
#   Step-Postbuckle/마지막프레임이 하드코딩된 get_hf 를 직접 부르지 않는다).
try:
    from eval_abaqus import (A0, R0, IDEAL_THRUST_COEFF, U_SUN,
                             calc_thrust_loss, get_hf)                       # noqa: F401
    from odbAccess import openOdb                                        # noqa: F401
except ImportError as _e:
    print("[eval_avg] odbAccess/abaqusConstants 를 못 읽었습니다 (%s)" % _e)
    print("[eval_avg] 'abaqus python eval_abaqus_avg.py <odb>' 로 실행해야 합니다.")
    sys.exit(2)

_PFX = '[eval_avg]'

INSTANCE_NAME = 'MEMBRANE-1'
DEFAULT_STEP_PREF = ('Step-Postbuckle',)      # 있으면 이것, 없으면 마지막 스텝
DEFAULT_N_MC = 2000
DEFAULT_SEED = 12345
DEFAULT_BLOCK = 500                           # 방향 블록 크기 상한
MEM_BUDGET = 2.0e6                            # 블록당 (요소수 x 방향수) 상한
STALE_WARN_S = 21600.0                        # 6시간 (§19.94 관례)
SELF_CHECK_TOL = 1e-9
# hf_ref 가 이보다 작으면 상대비교가 정의되지 않는다(분모 0). 실측 사고: --frame 1 처럼
#   **미변형(평평한) 프레임**을 읽으면 hf_ref 가 정확히 0 이 아니라 +-1e-19 수준의 부동소수
#   찌꺼기가 될 수 있고, 그때 상대비교는 0/0 이 되어 구현이 멀쩡해도 assert 가 터진다.
#   그 영역에서는 절대비교(둘 다 ~1e-13 이하)로 전환한다.
SELF_CHECK_FLOOR = 1e-12
# Murphey 2002: 완전평판 5도 기울기의 추력손실 = 1.108%.
#   (A0*cos5 + 2*R0*cos^3 5)/(A0+2*R0) = 0.9889172 -> 손실 1.1083% (아래 --check 에서 항등식 확인)
MURPHEY_5DEG = 0.0110828196
MURPHEY_TOL = 1e-4
TILT_SCAN_DEG = (0.0, 2.0, 5.0, 10.0, 20.0)

USAGE = """usage: abaqus python eval_abaqus_avg.py <odb> [options]

  --step NAME    읽을 스텝 (기본: Step-Postbuckle 이 있으면 그것, 없으면 마지막 스텝)
  --frame N      스텝 안 프레임 번호, 1-based (기본: 마지막 프레임)
  --n-mc N       몬테카를로 방향 샘플 수 (기본 2000)
  --seed S       난수 시드 (기본 12345)
  --ideal K      both | total | wrinkle (기본 both)
  --json PATH    결과 JSON 경로 (기본 results/avg_<odb stem>_<step>.json)
  --flip         n_z<0 요소 법선을 뒤집는다 (기본)
  --no-flip      뒤집지 않는다(원시 법선)
  --check        검증 모드: 수직 자기검증 + Murphey 5도 앵커 항등식만 확인하고 rc 반환
  -h, --help     이 도움말
"""


# ------------------------------------------------------------------ 인자
def parse_args(argv):
    """수동 파싱. (args, err) 를 돌려준다. err 가 있으면 사용법 오류(rc=2).

    argparse 를 쓰지 않는다: Abaqus 내장 python 이 구형일 수 있고, 이 저장소의 다른
    프로브(base_state_probe.py)도 sys.argv 를 직접 읽는 관례를 따른다.
    """
    a = {'odb': None, 'step': None, 'frame': None, 'n_mc': DEFAULT_N_MC,
         'seed': DEFAULT_SEED, 'ideal': 'both', 'json': None, 'flip': True,
         'check': False, 'help': False}
    i = 0
    while i < len(argv):
        t = argv[i]
        if t in ('-h', '--help'):
            a['help'] = True
            i += 1
            continue
        if t in ('--step', '--frame', '--n-mc', '--seed', '--ideal', '--json'):
            if i + 1 >= len(argv):
                return None, '%s 뒤에 값이 필요하다' % t
            v = argv[i + 1]
            if t == '--step':
                a['step'] = v
            elif t == '--json':
                a['json'] = v
            elif t == '--ideal':
                if v not in ('both', 'total', 'wrinkle'):
                    return None, "--ideal 은 both|total|wrinkle 중 하나여야 한다: %r" % v
                a['ideal'] = v
            else:
                try:
                    n = int(v)
                except ValueError:
                    return None, '%s 는 정수여야 한다: %r' % (t, v)
                if t == '--n-mc':
                    if n < 1:
                        return None, '--n-mc 는 1 이상이어야 한다: %d' % n
                    a['n_mc'] = n
                elif t == '--seed':
                    a['seed'] = n
                else:
                    if n < 1:
                        return None, '--frame 은 1-based(1 이상)여야 한다: %d' % n
                    a['frame'] = n
            i += 2
            continue
        if t == '--flip':
            a['flip'] = True
        elif t == '--no-flip':
            a['flip'] = False
        elif t == '--check':
            a['check'] = True
        elif t.startswith('-'):
            return None, '알 수 없는 옵션: %r' % t
        else:
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
    """요소별 면적과 법선. **eval_abaqus.calc_thrust_loss 와 문자 그대로 같은 식**을 쓴다.

    quad: cross(v2-v0, v3-v1) / tri: cross(v1-v0, v2-v0), 면적=0.5*|cross|, 법선=cross/|cross|.
    n_z<0 이면 뒤집기(flip_normals). 퇴화요소(|cross|~0)는 법선 0, 면적 0.
    """
    areas = []
    normals = []
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
    if not areas:
        return (np.zeros(0), np.zeros((0, 3)))
    areas = np.concatenate(areas)
    normals = np.concatenate(normals, axis=0)
    if flip_normals:
        normals[normals[:, 2] < 0.0] *= -1.0
    return areas, normals


def read_geometry(odb_path, step_arg, frame_arg, flip):
    """ODB 를 한 번만 열어 좌표/법선/면적/area_ref/hf_ref 를 뽑는다. 실패면 None."""
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
            print(_PFX + " (주의) 이 프레임 frameValue=%.4f < 0.99 - 미완료 런일 수 있다."
                  % fval)

        try:
            instance = odb.rootAssembly.instances[INSTANCE_NAME]
        except Exception:
            try:
                avail = list(odb.rootAssembly.instances.keys())
            except Exception:
                avail = '?'
            print(_PFX + " !!! instance '%s' 가 없다. 있는 instance: %s"
                  % (INSTANCE_NAME, avail))
            return None

        connectivity = [e.connectivity for e in instance.elements]
        quads = [c for c in connectivity if len(c) == 4]
        tris = [c for c in connectivity if len(c) == 3]
        if not quads and not tris:
            print(_PFX + " !!! 쓸 수 있는 요소(4절점/3절점)가 없다: 요소 %d개 (연결성 길이 %s)"
                  % (len(connectivity),
                     sorted(set(len(c) for c in connectivity)) if connectivity else []))
            return None

        # 배열 크기 할당용 최대 노드 라벨 (get_hf 와 동일)
        if instance.nodes:
            max_label = instance.nodes[-1].label + 100
        else:
            max_label = 100000

        # [§19.89] 미변형(pristine) 기준면적: frameValue=0 프레임에서 계산 (get_hf 와 동일)
        area_ref = None
        try:
            f0 = None
            for sn in odb.steps.keys():
                for fr in odb.steps[sn].frames:
                    if abs(float(fr.frameValue)) < 1e-12:
                        f0 = fr
                        break
                if f0 is not None:
                    break
            if f0 is not None and 'COORD' in f0.fieldOutputs.keys():
                cf = f0.fieldOutputs['COORD'].getSubset(region=instance)
                nc0 = _node_coords(cf, max_label)
                a0, _n0 = face_geometry(nc0, quads, tris, False)
                s0 = float(np.sum(a0))
                if s0 > 0.0:
                    area_ref = s0
                    print(_PFX + " 미변형(pristine) 기준면적 area_ref = %.6f m^2 (frameValue=0)"
                          % area_ref)
        except Exception as e:
            print(_PFX + " area_ref 판독 실패(%s) -> 변형 면적 폴백" % e)

        try:
            coord_field = frame.fieldOutputs['COORD']
            coords = _node_coords(coord_field, max_label)
            areas, normals = face_geometry(coords, quads, tris, flip)
        except Exception as e:
            print(_PFX + " !!! 변형 프레임 좌표/법선 계산 실패: %s" % e)
            return None
        area_sum = float(np.sum(areas))
        if area_sum <= 0.0:
            print(_PFX + " !!! 변형 프레임 총면적이 0 이다. 계산 불가.")
            return None
        if area_ref is None:
            area_ref = area_sum
            print(_PFX + " (주의) area_ref 를 못 구해 변형 면적 %.6f m^2 로 대체한다." % area_ref)

        # 참조구현으로 수직입사 hf + 뒤집기 편향 실측 (get_hf 와 같은 호출)
        st = {'n_total': 0, 'n_backward': 0, 'area_total': 0.0, 'area_backward': 0.0}
        try:
            hf_ref = float(calc_thrust_loss(frame, instance, connectivity, max_label,
                                            area_ref=area_ref, flip_normals=flip,
                                            extra_out=st))
            hf_true = float(calc_thrust_loss(frame, instance, connectivity, max_label,
                                             area_ref=area_ref, flip_normals=False))
        except Exception as e:
            print(_PFX + " !!! eval_abaqus.calc_thrust_loss 호출 실패: %s" % e)
            return None
        fb = (100.0 * st['area_backward'] / st['area_total']
              if st['area_total'] > 0.0 else float('nan'))
        dl = (100.0 * (hf_ref - hf_true) / hf_true) if hf_true else float('nan')
        print(_PFX + " 기하: 요소 %d개 (quad %d / tri %d) | 변형 총면적 %.6f m^2 | "
              "면적비 변형/미변형 %.6f" % (len(areas), len(quads), len(tris), area_sum,
                                          area_sum / area_ref if area_ref else float('nan')))
        print(_PFX + " 뒤집기=%s | 뒤방향(n_z<0) 면적비 %.4f%% (%d/%d) | hf_flip=%.6e "
              "hf_true=%.6e delta=%.3f%%"
              % ('ON' if flip else 'OFF', fb, st['n_backward'], st['n_total'],
                 hf_ref, hf_true, dl))
        g = {'step': sname, 'step_how': how, 'frame_index': fi, 'n_frames': nfr,
             'frame_value': fval, 'n_quad': len(quads), 'n_tri': len(tris),
             'n_elem': len(areas), 'areas': areas, 'normals': normals,
             'area_ref': float(area_ref), 'area_sum': area_sum, 'hf_ref': hf_ref,
             'hf_true': hf_true, 'backward_area_pct': fb,
             'n_backward': st['n_backward'], 'n_total': st['n_total']}
        return g
    finally:
        try:
            odb.close()
        except Exception:
            pass


# ------------------------------------------------------------------ 목적함수
def _ideal_coeff(ideal, s_dirs):
    """방향별 '같은-방향 완전평판' 추력계수. ideal='total' 이면 방향과 무관한 상수."""
    if ideal == 'total':
        return np.full(len(s_dirs), IDEAL_THRUST_COEFF, dtype=np.float64)
    cosf = np.clip(-s_dirs[:, 2], 0.0, 1.0)
    return A0 * cosf + 2.0 * R0 * cosf ** 3


def hf_block(normals, areas, s_dirs, ideal, area_ref, area_sum):
    """방향 K개에 대한 hf 를 **한 블록**으로 계산한다. (hf_clamped(K), hf_raw(K)) 반환.

    변형된 요소 i 의 광선방향 s_hat(=광자의 진행방향, 참조구현의 U_SUN 역할)에 대한
    cos(theta_i) = clip(-n_i . s_hat, 0, 1) 이고, 그 요소가 광선방향으로 내는 추력은
    a_i*(A0*c_i + 2*R0*c_i^3) 이다(수직입사 c=1 이면 참조구현과 동일한 식).
    hf = (Sigma a_i * coeff(s) - actual) / (area_ref * coeff(s)),  coeff = 이상 계수.
    분모가 0 이면 nan(스칠 때 wrinkle 정의의 특이점).
    """
    c = np.clip(-np.dot(normals, s_dirs.T), 0.0, 1.0)                 # (n_elem, K)
    actual = np.sum(areas[:, None] * (A0 * c + 2.0 * R0 * c ** 3), axis=0)
    coeff = _ideal_coeff(ideal, s_dirs)
    loss = area_sum * coeff - actual
    den = area_ref * coeff
    raw = np.full(len(s_dirs), float('nan'), dtype=np.float64)
    ok = den > 0.0
    raw[ok] = loss[ok] / den[ok]
    hf = np.where(np.isfinite(raw), np.maximum(raw, 0.0), float('nan'))
    return hf, raw


def hf_over_dirs(normals, areas, s_dirs, ideal, area_ref, area_sum,
                 n_block=DEFAULT_BLOCK):
    """방향별 hf 를 **블록으로 나눠** 계산한다(메모리 상한).

    블록 크기 = min(n_block, MEM_BUDGET/요소수). 10k 요소 x 2000 방향 float64 = 160MB 를
    한 번에 만들지 않기 위한 것. 결과는 방향 순서 그대로 이어붙인다.
    """
    n = int(normals.shape[0])
    blk = max(1, min(int(n_block), int(MEM_BUDGET / float(max(n, 1)))))
    k = len(s_dirs)
    hf = np.empty(k, dtype=np.float64)
    raw = np.empty(k, dtype=np.float64)
    for i in range(0, k, blk):
        j = min(i + blk, k)
        _h, _r = hf_block(normals, areas, s_dirs[i:j], ideal, area_ref, area_sum)
        hf[i:j] = _h
        raw[i:j] = _r
    return hf, raw


def sample_dirs(mode, n, seed=None):
    """반구 전방향 방향 샘플러. s_hat 은 **광자의 진행방향**(참조구현 U_SUN 의 역할)이다.

    hemi : 입체각 균일 (cos(phi) ~ U[0,1)) -> 단순 균일평균
    disk : 단위 디스크 균일 (x,y) -> 밀도가 cos(phi) 에 비례하는 전력가중 평균
    """
    rng = np.random.RandomState(seed)
    if mode == 'hemi':
        cz = rng.rand(n)
        psi = 2.0 * math.pi * rng.rand(n)
        st = np.sqrt(np.maximum(1.0 - cz * cz, 0.0))
        return -np.stack((st * np.cos(psi), st * np.sin(psi), cz), axis=1)
    if mode == 'disk':
        r = np.sqrt(rng.rand(n))
        psi = 2.0 * math.pi * rng.rand(n)
        x = r * np.cos(psi)
        y = r * np.sin(psi)
        z = np.sqrt(np.maximum(1.0 - x * x - y * y, 0.0))
        return -np.stack((x, y, z), axis=1)
    raise ValueError("mode 는 'hemi' 또는 'disk' 여야 한다: %r" % (mode,))


def dir_at_tilt(phi_deg, psi_deg=0.0):
    """고정 기울기 방향 1개. s_hat = -(sin(p)cos(t), sin(p)sin(t), cos(p))."""
    p = math.radians(phi_deg)
    t = math.radians(psi_deg)
    return np.array([-(math.sin(p) * math.cos(t)), -(math.sin(p) * math.sin(t)),
                     -math.cos(p)], dtype=np.float64)


def _flat_sail_hf(phi_deg, ideal):
    """**완전평판** 돛을 기울기 phi 로 비췄을 때의 hf (해석식, 기하 없음).

    actual = a*(A0*cos + 2*R0*cos^3) 이므로
      total   : 1 - (A0 cos + 2 R0 cos^3)/(A0+2R0)
      wrinkle : 0 (분자가 분모와 같은 식 -> 정확히 0)
    phi=5도 total = 0.0110828 = Murphey 2002 의 1.108% (아래 --check 앵커).
    """
    if ideal == 'wrinkle':
        return 0.0
    c = math.cos(math.radians(phi_deg))
    return 1.0 - (A0 * c + 2.0 * R0 * c ** 3) / IDEAL_THRUST_COEFF


def stats_of(v):
    """nan 제거 후 평균/SE/표준편차/분위수."""
    v = np.asarray(v, dtype=np.float64)
    good = np.isfinite(v)
    w = v[good]
    nn = int(w.size)
    out = {'n': nn, 'n_total': int(v.size), 'n_nan': int(v.size) - nn,
           'mean': float('nan'), 'se': float('nan'), 'std': float('nan'),
           'min': float('nan'), 'p1': float('nan'), 'p50': float('nan'),
           'p99': float('nan'), 'max': float('nan')}
    if nn == 0:
        return out
    out['mean'] = float(np.mean(w))
    sd = float(np.std(w, ddof=1)) if nn > 1 else 0.0
    out['std'] = sd
    out['se'] = sd / math.sqrt(nn)
    out['min'] = float(np.min(w))
    out['max'] = float(np.max(w))
    q = np.percentile(w, [1.0, 50.0, 99.0])
    out['p1'], out['p50'], out['p99'] = float(q[0]), float(q[1]), float(q[2])
    return out


# ------------------------------------------------------------------ 자기검증
def vertical_self_check(normals, areas, area_ref, area_sum, hf_ref, ideal):
    """수직입사(s_hat = U_SUN)에서 자기 기하 결과가 참조구현과 같은지 assert.

    두 정의 모두 cosf=1 이 되어 참조구현과 같은 식으로 퇴화하므로, 어느 쪽이든 일치해야 한다.
    hf_ref 가 SELF_CHECK_FLOOR 이하(미변형/완전평판 프레임: 비율이 정의되지 않음)면
    절대비교로 전환한다. 두 값 모두 그 영역에서는 ~1e-13 이하이다.
    """
    s = np.asarray(U_SUN, dtype=np.float64).reshape(1, 3)
    hf, _raw = hf_block(normals, areas, s, ideal, area_ref, area_sum)
    val = float(hf[0])
    mode = 'rel'
    if abs(hf_ref) > SELF_CHECK_FLOOR:
        rel = abs(val - hf_ref) / abs(hf_ref)
    else:
        rel = abs(val - hf_ref)
        mode = 'abs(hf_ref~0)'
    assert rel < SELF_CHECK_TOL, (
        "자기검증 실패(ideal=%s): hf_normal_self=%.17g vs hf_ref=%.17g %s=%.3e >= %.0e"
        % (ideal, val, hf_ref, mode, rel, SELF_CHECK_TOL))
    return val, rel, mode


def anchor_check(verbose=True):
    """Murphey 2002 5도 앵커를 **항등식**으로 확인한다(플랫 돛 공식 + 자체 코드 경로)."""
    ok = True
    f = _flat_sail_hf(5.0, 'total')
    d = abs(f - MURPHEY_5DEG)
    ok &= (d < MURPHEY_TOL)
    if verbose:
        print(_PFX + ":check flat formula   phi=5 total = %.10f (%.4f%%) vs Murphey "
              "%.4f%% | d=%.3e < %.0e %s"
              % (f, 100.0 * f, 100.0 * MURPHEY_5DEG, d, MURPHEY_TOL, 'OK' if d < MURPHEY_TOL else 'FAIL'))
    # 같은 값을 **자체 hf_block 코드 경로**로도 계산 (합성 완전평판 요소 1개)
    fl_n = np.array([[0.0, 0.0, 1.0]])
    fl_a = np.array([1.0])
    hf, _raw = hf_block(fl_n, fl_a, dir_at_tilt(5.0).reshape(1, 3), 'total', 1.0, 1.0)
    d2 = abs(float(hf[0]) - MURPHEY_5DEG)
    ok &= (d2 < MURPHEY_TOL)
    if verbose:
        print(_PFX + ":check code path      phi=5 total = %.10f | d=%.3e < %.0e %s"
              % (float(hf[0]), d2, MURPHEY_TOL, 'OK' if d2 < MURPHEY_TOL else 'FAIL'))
    # 완전평판 wrinkle 정의는 모든 방향에서 0 (항등식)
    d3 = 0.0
    for phi in (0.0, 5.0, 20.0, 45.0, 80.0):
        hw, _r = hf_block(fl_n, fl_a, dir_at_tilt(phi).reshape(1, 3), 'wrinkle', 1.0, 1.0)
        d3 = max(d3, abs(float(hw[0])))
    ok &= (d3 < 1e-15)
    if verbose:
        print(_PFX + ":check wrinkle=0      완전평판 wrinkle 최대 |hf| = %.3e < 1e-15 %s"
              % (d3, 'OK' if d3 < 1e-15 else 'FAIL'))
    return ok


# ------------------------------------------------------------------ 스캔표
def tilt_scan(normals, areas, area_ref, area_sum):
    """고정 기울기 스캔: 플랫 돛 공식 vs 이 ODB 기하. 교차항(=기울기 몫)을 눈에 보이게."""
    rows = []
    for phi in TILT_SCAN_DEG:
        s = dir_at_tilt(phi).reshape(1, 3)
        ht, _ = hf_block(normals, areas, s, 'total', area_ref, area_sum)
        hw, _ = hf_block(normals, areas, s, 'wrinkle', area_ref, area_sum)
        rows.append({'phi_deg': phi,
                     'flat_total': _flat_sail_hf(phi, 'total'),
                     'flat_wrinkle': _flat_sail_hf(phi, 'wrinkle'),
                     'sail_total': float(ht[0]),
                     'sail_wrinkle': float(hw[0]),
                     'sail_extra_after_tilt': float(ht[0]) - _flat_sail_hf(phi, 'total')})
    return rows


# ------------------------------------------------------------------ 표시
def print_stats(tag, st, n_clip):
    print(_PFX + " MC %-18s N=%d mean=%.6e (%.4f %%) SE=%.3e std=%.3e "
          "min=%.4e p1=%.4e p50=%.4e p99=%.4e max=%.4e nan=%d clip=%d"
          % (tag, st['n'], st['mean'], 100.0 * st['mean'], st['se'], st['std'],
             st['min'], st['p1'], st['p50'], st['p99'], st['max'], st['n_nan'], n_clip))


def print_scan(rows):
    print(_PFX + " 경사 스캔 (psi=0 고정, 플랫 돛 공식 vs 이 ODB 기하)")
    print(_PFX + "   %-8s %-12s %-13s %-13s %-13s %-13s"
          % ('phi_deg', 'flat_total', 'flat_wrinkle', 'sail_total', 'sail_wrinkle',
             'sail-flat'))
    for r in rows:
        print(_PFX + "   %-8.1f %-12.6f %-13.6f %-13.6f %-13.6f %-13.6f"
              % (r['phi_deg'], r['flat_total'], r['flat_wrinkle'], r['sail_total'],
                 r['sail_wrinkle'], r['sail_extra_after_tilt']))
    print(_PFX + "   (flat_total 은 완전평판 해석식, sail_* 는 이 ODB 변형기하. "
          "sail-flat 은 같은 기울기에서 주름 때문에 추가로 잃은 몫.)")


# ------------------------------------------------------------------ main
def main(argv=None):
    t_start = time.time()
    if argv is None:
        argv = sys.argv[1:]
    print(_PFX + " build: eval_abaqus_avg.py (19.99) | numpy %s | U_SUN=%s | "
          "A0=%.6f R0=%.6f IDEAL=%.6f"
          % (np.__version__, np.asarray(U_SUN, dtype=float).tolist(),
             A0, R0, IDEAL_THRUST_COEFF))
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
    g = read_geometry(args['odb'], args['step'], args['frame'], args['flip'])
    if g is None:
        return 1

    ideals = ['total', 'wrinkle'] if args['ideal'] == 'both' else [args['ideal']]

    # ---- (1) 수직입사 자기검증: 자체 기하 x 자체 목적함수 == eval_abaqus.calc_thrust_loss
    self_vals = {}
    for ideal in ideals:
        try:
            v, rel, smode = vertical_self_check(g['normals'], g['areas'], g['area_ref'],
                                                g['area_sum'], g['hf_ref'], ideal)
        except AssertionError as e:
            print(_PFX + " !!! %s" % e)
            print(_PFX + " !!! 구현 드리프트: 이 스크립트의 목적함수가 eval_abaqus.py 와 다르다.")
            return 1
        self_vals[ideal] = v
        print(_PFX + " 자기검증(ideal=%s): hf_normal_self=%.10e vs hf_ref=%.10e "
              "(calc_thrust_loss) | %s=%.3e < %.0e OK"
              % (ideal, v, g['hf_ref'], smode, rel, SELF_CHECK_TOL))

    # ---- (2) 앵커 항등식 (Murphey 2002: 5도 기울기 1.108%)
    a_ok = anchor_check()
    if not a_ok:
        print(_PFX + " !!! Murphey 5도 앵커 항등식 실패.")
        if args['check']:
            return 1
    if args['check']:
        print(_PFX + ":check 결과: 자기검증 %d개 OK, 앵커 %s -> rc=%d"
              % (len(self_vals), 'OK' if a_ok else 'FAIL', 0 if a_ok else 1))
        return 0 if a_ok else 1

    # ---- (3) MC 4(또는 2)조합
    s_dirs = {}
    for mode in ('hemi', 'disk'):
        s_dirs[mode] = sample_dirs(mode, args['n_mc'], args['seed'])
    mc = {}
    for mode in ('hemi', 'disk'):
        for ideal in ideals:
            hf, raw = hf_over_dirs(g['normals'], g['areas'], s_dirs[mode], ideal,
                                   g['area_ref'], g['area_sum'])
            key = '%s_%s' % (mode, ideal)
            st = stats_of(hf)
            n_clip = int(np.count_nonzero(np.isfinite(raw) & (raw < 0.0)))
            st['n_clipped'] = n_clip
            st['mean_raw'] = float(np.mean(raw[np.isfinite(raw)])) if st['n'] else float('nan')
            mc[key] = st
            print_stats(key, st, n_clip)

    # ---- (4) 경사 스캔표
    rows = tilt_scan(g['normals'], g['areas'], g['area_ref'], g['area_sum'])
    print_scan(rows)

    # ---- (5) JSON
    out_path = args['json']
    if out_path is None:
        stem = os.path.splitext(os.path.basename(args['odb']))[0]
        out_path = os.path.join('results', 'avg_%s_%s.json' % (stem, g['step']))
    rec = {'script': 'eval_abaqus_avg.py', 'section': '19.99',
           'generated': time.strftime('%Y-%m-%d %H:%M:%S'),
           'probe_run_epoch': t_start,
           'odb': ident,
           'step': g['step'], 'step_how': g['step_how'], 'frame_index': g['frame_index'],
           'n_frames': g['n_frames'], 'frame_value': g['frame_value'],
           'flip_normals': bool(args['flip']),
           'n_elem': g['n_elem'], 'n_quad': g['n_quad'], 'n_tri': g['n_tri'],
           'area_ref': g['area_ref'], 'area_deformed_sum': g['area_sum'],
           'backward_area_pct': g['backward_area_pct'],
           'hf_ref': g['hf_ref'], 'hf_true_no_flip': g['hf_true'],
           'hf_normal_self': self_vals,
           'self_check_tol': SELF_CHECK_TOL,
           'murphey_5deg': MURPHEY_5DEG, 'murphey_ok': bool(a_ok),
           'n_mc': args['n_mc'], 'seed': args['seed'],
           'mc': mc, 'tilt_scan': rows,
           'note': 'ideal=total: flat normal-incidence reference (Murphey-comparable). '
                   'ideal=wrinkle: same-direction flat reference (isolates wrinkle). '
                   'hemi=uniform solid angle, disk=power(cos)-weighted.'}
    try:
        d = os.path.dirname(os.path.abspath(out_path))
        if d and not os.path.isdir(d):
            os.makedirs(d)
        with io.open(out_path, 'w', encoding='utf-8') as f:
            f.write(json.dumps(rec, ensure_ascii=True, indent=2, sort_keys=True))
        print(_PFX + " JSON 저장: %s" % os.path.abspath(out_path))
    except Exception as e:
        print(_PFX + " !!! JSON 저장 실패: %s" % e)
        return 1

    # ---- (6) 기계 파싱용 KV (RESULTS 계약과 무관한 별도 이름공간)
    kv = {'hf_ref': g['hf_ref'], 'hf_normal_self': self_vals[ideals[0]],
          'area_ref': g['area_ref'], 'n_elem': g['n_elem'],
          'n_mc': args['n_mc'], 'seed': args['seed'], 'step': g['step'],
          'frame': g['frame_index']}
    for k in ('hemi_total', 'disk_total', 'hemi_wrinkle', 'disk_wrinkle'):
        kv[k] = mc[k]['mean'] if k in mc else float('nan')
    for k in ('hemi_total', 'disk_total', 'hemi_wrinkle', 'disk_wrinkle'):
        kv[k + '_se'] = mc[k]['se'] if k in mc else float('nan')
    print('[eval_avg:KV] ' + ' '.join('%s=%s' % (k, kv[k]) for k in sorted(kv)))
    print(_PFX + " 완료 (rc=0)")
    return 0


if __name__ == '__main__':
    sys.exit(main())
