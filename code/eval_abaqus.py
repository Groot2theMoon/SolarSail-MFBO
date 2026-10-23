"""
아바쿠스 해석이 종료된 후, 그 결과 파일인 .odb 를 읽어서 필요한 lf, hf값을 계산하는 코드
run_abaqus.py 코드 내에서 호출되어 사용되며, 
"abaqus python eval_abaqus.py (lf_odb파일) (hf_odb파일) HF" 형태로 호출

Return Format (stdout):
    "RESULTS:val1,val2,..." 형식으로 출력하여 상위 프로세스가 파싱할 수 있게 함.
    - LF Mode: "RESULTS:lf1,lf2,lf3"
    - HF Mode: "RESULTS:lf1,lf2,lf3,hf" 
    - 별도 줄: "ZK:state_wrinkled_frac=... zk_intensity_mean=... zk_intensity_p95=... zk_dir_deg=..."
      (Zhang & Kiendl 변형률 기반 주름 상태. RESULTS 계약과 무관하다.)
      E 필드가 없으면 "ZK:FAIL reason=..." 한 줄만 내고 RESULTS 는 평소대로 낸다.

LF : 선형/프리스트레스 해석에서 얻을 수 있는 지표들
    - lf1 (Compressed Volume Ratio): 압축 주응력이 발생하는 요소의 부피 비율.
    - lf2 (Wrinkling Energy Density): 인장-압축 혼합 응력 상태에서의 변형 에너지 밀도.
    - lf3 (Stress Anisotropy): 주응력 간의 비대칭성 비율.
    - ZK (Zhang & Kiendl): Green-Lagrange E 의 주성분 E1/E2 로 taut/wrinkled/slack 을
      판정하고, wrinkled 영역의 |E2| 평균/p95 와 n1 대표 방향각을 낸다.
      Abaqus NLGEOM=ON 이 `E` 대신 `LE`(logarithmic)를 출력하면 주성분을
      `E = 0.5(exp(2*LE)-1)` 로 변환해 같은 Green-Lagrange 기준으로 계산한다.
      `zk_metrics.txt` 에는 추가로 기준면적 기반 압축 변형률 `zk_Iw` 를 남긴다
      (public `ZK:` 4필드는 불변).

HF :
    - hf (Thrust Loss): 이상적인 평판 대비, 주름진 형상으로 인해 손실된 추력의 양.
    - Solar Radiation Pressure (SRP) 모델을 적용하여 변형된 요소들의 법선 벡터 변화를 적분.
    - DN (delta_n): 면적가중 등가법선의 이상법선 대비 편차각. 별도 `DN:` 줄로만 낸다.

LF 상태 스텝 (2026-10-21):
    - paper 경로(PAPER_MATCH=True): LF 모델의 최종(그리고 유일한) 해석 스텝은
      Step-GlobalTension 이다. 케이블/클램프 초기응력 + u3=0 프리스트레스 상태.
    - 구 disp 경로(PAPER_MATCH=False): Step-HighTension 이 최종 상태다.
    - 그래서 선택 우선순위는 HighTension -> GlobalTension -> ClampTension -> 마지막 스텝.
      (GlobalTension 을 먼저 고르면 구 disp 기준선이 중간 상태를 읽어 값이 바뀐다.)
"""
from odbAccess import openOdb
from abaqusConstants import *
import numpy as np
import sys
import os

# Solar Radiation Pressure (SRP) Constants
R0 = 0.926256
A0 = 0.073744
U_SUN = np.array([0.0, 0.0, -1.0])  # 태양광 입사 방향
# 이상적인 평판의 추력 계수 (2*R0 + A0)
IDEAL_THRUST_COEFF = 2.0 * R0 + A0 

# ---------------------------------------------------------------------------
# LF 상태 스텝 선택 (2026-10-21)
#   paper 경로(PAPER_MATCH=True): LF 모델에 Step-HighTension 이 없다 ->
#     Step-GlobalTension 이 **최종 프리스트레스** 상태다.
#   구 disp 경로(PAPER_MATCH=False): Step-HighTension 이 최종 상태다.
#   우선순위 = HighTension -> GlobalTension -> ClampTension -> 마지막 스텝.
#   (GlobalTension 을 먼저 고르면 구 disp 기준선이 중간 상태를 읽어 lf1~lf3 가 바뀐다.)
# ---------------------------------------------------------------------------
LF_STEP_PREFERENCE = ('Step-HighTension', 'Step-GlobalTension', 'Step-ClampTension')


def _pick_lf_step(odb):
    """(step_name, frame, available_steps) 반환. 프레임이 있는 첫 후보/마지막 스텝."""
    try:
        names = list(odb.steps.keys())
    except Exception:
        names = []
    for name in LF_STEP_PREFERENCE:
        if name in names:
            try:
                frames = odb.steps[name].frames
                if len(frames) > 0:
                    return name, frames[-1], names
            except Exception:
                pass
    last = None
    for name in names:
        try:
            frames = odb.steps[name].frames
            if len(frames) > 0:
                last = (name, frames[-1])
        except Exception:
            pass
    if last is not None:
        return last[0], last[1], names
    return None, None, names


def get_lf(odb_path):
    """
    LF ODB 의 **최종 상태**에서 주름 발생 가능성을 나타내는 대리 지표(Proxy Metrics) 3종을 추출.
    
    Returns:
        list: [lf1, lf2, lf3]
    """
    try:
        if not os.path.exists(odb_path):
            print("LF Error: File not found (%s)" % odb_path)
            return [None, None, None]

        odb = openOdb(path=odb_path, readOnly=True)
        
        step_name, frame, _names = _pick_lf_step(odb)
        if step_name is None or frame is None:
            print("LF Error: no step with frames. available=%s" % (_names,))
            odb.close()
            return [None, None, None]
        print("[LF] state step = %s (available=%s)" % (step_name, _names))
        instance = odb.rootAssembly.instances['MEMBRANE-1']
        
        if 'EVOL' not in frame.fieldOutputs.keys():
            print("LF Error: EVOL field not found.")
            odb.close()
            return [None, None, None]
            
        evol_field = frame.fieldOutputs['EVOL'].getSubset(region=instance)
        
        vol_map = {}
        for val in evol_field.values:
            vol_map[val.elementLabel] = val.data
            
        total_vol = sum(vol_map.values())
        if total_vol == 0:
            odb.close()
            return [None, None, None]
        
        stress_field = frame.fieldOutputs['S'].getSubset(position=CENTROID, region=instance)
        s_max_field = stress_field.getScalarField(invariant=MAX_INPLANE_PRINCIPAL)
        s_min_field = stress_field.getScalarField(invariant=MIN_INPLANE_PRINCIPAL)
        
        s1_temp = {}
        s2_temp = {}
        
        for val in s_max_field.values:
            lbl = val.elementLabel
            if lbl not in s1_temp: s1_temp[lbl] = []
            s1_temp[lbl].append(val.data)
            
        for val in s_min_field.values:
            lbl = val.elementLabel
            if lbl not in s2_temp: s2_temp[lbl] = []
            s2_temp[lbl].append(val.data)
            
        sorted_labels = sorted(vol_map.keys())
        
        vols = []
        sig_1 = []
        sig_2 = []
        
        for lbl in sorted_labels:
            # 부피 (EVOL은 요소당 1개 값)
            vols.append(vol_map[lbl])
            
            # 응력 (여러 적분점의 평균값 사용)
            if lbl in s1_temp:
                sig_1.append(np.mean(s1_temp[lbl]))
            else:
                sig_1.append(0.0)
                
            if lbl in s2_temp:
                sig_2.append(np.mean(s2_temp[lbl]))
            else:
                sig_2.append(0.0)

        vols = np.array(vols)
        sig_1 = np.array(sig_1)
        sig_2 = np.array(sig_2)
        
        #LF1: 압축 응력(음수)이 발생하는 요소의 부피 비율
        # 압축 응력이 발생하는 영역이 넓을수록, 주름이 발생할 수 있는 영역이 넓어짐
        mask_comp = (sig_2 < -1e-7)
        lf1_area_ratio = np.sum(vols[mask_comp]) / total_vol
        
        #LF2: Wrinkling Energy Density
        # 인장(sig1>0)과 압축(sig2<0)이 공존하는 영역의 변형 에너지
        mask_tc = (sig_1 > 1e-7) & (sig_2 < -1e-7)
        if np.any(mask_tc):
            nu, E = 0.34, 2.5e9
            # Wrinkling Energy Formula
            e_density = -(nu / E) * sig_1[mask_tc] * sig_2[mask_tc]
            lf2_energy_density = np.sum(e_density * vols[mask_tc]) / np.sum(vols[mask_tc])
        else:
            lf2_energy_density = 0.0
            
        #LF3: Stress Anisotropy Ratio
        max_c = np.max(np.abs(sig_2)) if len(sig_2) > 0 else 0
        if max_c > 1e-9:
            # 유의미한 압축 응력이 있는 곳에서만 계산 (Peak의 5% 이상)
            mask_aniso = (np.abs(sig_2) > 0.05 * max_c) & (sig_1 > 1e-7)
            if np.any(mask_aniso):
                lf3_anisotropy = np.mean(np.abs(sig_2[mask_aniso] / sig_1[mask_aniso]))
            else:
                lf3_anisotropy = 0.0
        else:
            lf3_anisotropy = 0.0
            
        odb.close()
        return [float(lf1_area_ratio), float(lf2_energy_density), float(lf3_anisotropy)]

    except Exception as e:
        print("LF Error: " + str(e))
        try: odb.close()      # P1-1: 열린 odb 핸들이 남으면 다음 실행이 잠긴다
        except Exception: pass
        return [None, None, None]


# ---------------------------------------------------------------------------
# [paper-match] ZK (Zhang & Kiendl) 변형률 기반 주름 지표 (2026-10-21)
#
#   LF ODB 의 최종 상태(paper: Step-GlobalTension, disp: Step-HighTension)에서
#   Green-Lagrange 변형률 E 를 읽어 주성분 E1 >= E2 와 주방향 n1 을 계산한다.
#   Abaqus/Standard NLGEOM=ON 은 ODB 에 `E` 대신 `LE` 를 쓰는 경우가 있다. 그때는
#   주성분을 E = 0.5*(exp(2*LE)-1) 로 변환한다(고유벡터는 동일).
#
#   상태(요소별, 면적가중):
#       taut      : E2 > 0
#       wrinkled  : E1 > 0 and E2 <= 0
#       slack     : E1 <= 0
#   intensity   : wrinkled 영역의 |E2| 면적가중 평균 / p95 (public ZK)
#   I_w         : 기준면적(frameValue=0 EVOL) 가중 **전체** 압축 변형률 (sidecar)
#                 I_w = (1/A_ref)[Σ_wr w_ref·max(|E2|-νE1,0)
#                                  + Σ_slack w_ref·(|E2|+ν|E1|)]
#                 public `zk_intensity_mean` 과 달리 taut=0, slack 포함, ν 보정, 기준면적.
#   direction   : wrinkled 영역 n1 의 면적가중 axial mean (global x-y, mod 180 deg)
#
#   Abaqus shell E 성분은 요소 local frame(1축 = node1->node2)이다. E1/E2 는 회전
#   불변이라 그대로 쓰지만, n1 각도는 요소 local 1축 각도를 더해 global x-y 로
#   회전해야 한다(격자 S4 는 local1=+x, S3 경계는 회전 필요).
#
#   출력 계약:
#       "ZK:state_wrinkled_frac=.. zk_intensity_mean=.. zk_intensity_p95=.. zk_dir_deg=.."
#       계산 불가: "ZK:FAIL reason=.." (ASCII 토큰). RESULTS 3/4필드는 절대 건드리지 않는다.
# ---------------------------------------------------------------------------
ZK_DIR_R_MIN = 1.0e-3    # axial mean 이 의미를 가지려면 R >= 이 값 (등방/직교면 nan)
ZK_NU = 0.34             # Kapton 푸아송비 — `zk_Iw` 의 압축 변형률 보정에만 사용


def _zk_fail(reason):
    try:
        _s = str(reason)
    except Exception:
        _s = 'unknown'
    _s = ''.join(ch if (33 <= ord(ch) < 127) else '_' for ch in _s)
    return {'ok': False, 'reason': (_s or 'unknown')}


def _zk_float(x):
    """EVOL/스칼라 값 -> float. tuple/list/ndarray 면 첫 성분. 실패는 None."""
    try:
        if hasattr(x, '__len__') and not isinstance(x, str):
            return float(x[0])
        return float(x)
    except Exception:
        return None


def _zk_field_len(field):
    try:
        return len(field.values)
    except Exception:
        return 0


def _zk_reference_evol(odb, instance):
    """frameValue=0 프레임의 EVOL 을 {elementLabel: volume} 로 반환. 없으면 {}.

    `I_w` 의 기준면적 가중이다. frameValue=0 의 EVOL 은 미변형(초기) 부피이므로,
    두께가 일정한 막에서는 기준면적 가중과 같다. `hf` 의 `area_ref` 와 같은
    '미변형 기준'을 쓰되, `I_w` 는 요소별 가중이 필요해 EVOL 을 그대로 쓴다.
    """
    ref = {}
    try:
        f0 = None
        for sn in odb.steps.keys():
            for fr in odb.steps[sn].frames:
                try:
                    if (abs(float(fr.frameValue)) < 1e-12
                            and 'EVOL' in fr.fieldOutputs.keys()):
                        f0 = fr
                        break
                except Exception:
                    continue
            if f0 is not None:
                break
        if f0 is None:
            return ref
        ef = f0.fieldOutputs['EVOL'].getSubset(region=instance)
        for val in ef.values:
            v = _zk_float(val.data)
            if v is not None and np.isfinite(v) and v > 0.0:
                ref[val.elementLabel] = v
    except Exception:
        return {}
    return ref


def _zk_get_e_field(frame, instance):
    """E(우선) 또는 LE 를 integration point -> centroid -> raw 순으로 시도한다.

    반환: (field, field_name, 실패사유). field_name 은 'E' 또는 'LE'.
    Abaqus NLGEOM=ON 은 E 대신 LE 를 저장할 수 있다(get_zk 에서 E 로 변환).
    """
    keys = frame.fieldOutputs.keys()
    for _name in ('E', 'LE'):
        if _name not in keys:
            continue
        raw = frame.fieldOutputs[_name]
        for _pos in (INTEGRATION_POINT, CENTROID):
            try:
                cand = raw.getSubset(position=_pos, region=instance)
            except Exception:
                continue
            if _zk_field_len(cand) > 0:
                return cand, _name, None
        try:
            cand = raw.getSubset(region=instance)
        except Exception:
            cand = raw
        if _zk_field_len(cand) > 0:
            return cand, _name, None
    if ('E' not in keys) and ('LE' not in keys):
        return None, None, 'E_field_missing'
    return None, None, 'E_field_empty'


def _zk_element_local_angles(instance, labels):
    """요소 local 1축(node1->node2)의 global x-y 각도[rad]. 불가능하면 None."""
    try:
        node_xy = {}
        for nd in instance.nodes:
            c = nd.coordinates
            node_xy[nd.label] = (float(c[0]), float(c[1]))
        conn = {}
        for el in instance.elements:
            conn[el.label] = el.connectivity
    except Exception:
        return None
    out = []
    for lbl in labels:
        c = conn.get(lbl)
        if not c or len(c) < 2:
            out.append(float('nan'))
            continue
        p1 = node_xy.get(c[0])
        p2 = node_xy.get(c[1])
        if p1 is None or p2 is None:
            out.append(float('nan'))
            continue
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        #   S4 는 Abaqus local 1축이 1-2 와 4-3 edge 의 평균 방향이다(격자 사각형은 +x).
        #   S3 는 1-2 edge 방향이다.
        if len(c) >= 4:
            p3 = node_xy.get(c[2])
            p4 = node_xy.get(c[3])
            if p3 is not None and p4 is not None:
                dx = 0.5 * (dx + (p3[0] - p4[0]))
                dy = 0.5 * (dy + (p3[1] - p4[1]))
        if (abs(dx) + abs(dy)) < 1e-15:
            out.append(float('nan'))
            continue
        out.append(float(np.arctan2(dy, dx)))
    return np.array(out, dtype=np.float64)


def _zk_weighted_percentile(values, weights, q):
    v = np.asarray(values, dtype=np.float64).ravel()
    w = np.asarray(weights, dtype=np.float64).ravel()
    if v.size == 0:
        return 0.0
    order = np.argsort(v)
    v = v[order]
    w = w[order]
    cw = np.cumsum(w)
    if cw[-1] <= 0.0:
        return float(np.percentile(v, q))
    target = (float(q) / 100.0) * cw[-1]
    i = int(np.searchsorted(cw, target, side='left'))
    if i >= v.size:
        i = v.size - 1
    return float(v[i])


def get_zk(odb_path):
    """LF ODB 의 최종 상태에서 Zhang & Kiendl 지표를 계산한다. 예외를 던지지 않는다."""
    odb = None
    try:
        if not odb_path or not os.path.exists(odb_path):
            return _zk_fail('file_not_found')
        odb = openOdb(path=odb_path, readOnly=True)
        step_name, frame, names = _pick_lf_step(odb)
        if step_name is None or frame is None:
            return _zk_fail('lf_state_step_missing')
        try:
            instance = odb.rootAssembly.instances['MEMBRANE-1']
        except Exception:
            return _zk_fail('membrane_instance_missing')

        e_field, strain_field, why = _zk_get_e_field(frame, instance)
        if e_field is None:
            _fail = _zk_fail(why or 'E_field_missing')
            #   진단: 선택 프레임에 실제로 어떤 field 가 있는지 sidecar 에 남긴다.
            try:
                _fail['fields'] = ','.join(sorted(frame.fieldOutputs.keys()))
                _fail['frame_value'] = frame.frameValue
            except Exception:
                pass
            return _fail

        # ---- E 성분을 요소별 평균 (section point/integration point 평균) ----
        e_acc = {}
        for val in e_field.values:
            lbl = getattr(val, 'elementLabel', None)
            if lbl is None:
                continue
            d = val.data
            try:
                if len(d) >= 6:
                    e11, e22, e12 = float(d[0]), float(d[1]), float(d[3])
                elif len(d) >= 3:
                    e11, e22, e12 = float(d[0]), float(d[1]), float(d[2])
                else:
                    continue
            except Exception:
                continue
            if not (np.isfinite(e11) and np.isfinite(e22) and np.isfinite(e12)):
                continue
            e_acc.setdefault(lbl, []).append((e11, e22, e12))
        if not e_acc:
            return _zk_fail('E_field_empty')

        # ---- EVOL (면적가중). 없으면 균등가중 폴백 ----
        vol_map = {}
        try:
            if 'EVOL' in frame.fieldOutputs.keys():
                evol_field = frame.fieldOutputs['EVOL'].getSubset(region=instance)
                for val in evol_field.values:
                    v = _zk_float(val.data)
                    if v is not None and np.isfinite(v) and v > 0.0:
                        vol_map[val.elementLabel] = v
        except Exception:
            vol_map = {}
        weight_mode = 'evol' if vol_map else 'uniform'
        #   `I_w` 전용 기준면적 가중: frameValue=0 EVOL. 없으면 현재 EVOL/균등으로 폴백.
        ref_vol_map = _zk_reference_evol(odb, instance)

        labels_all = sorted(e_acc.keys())
        labels, e11_l, e22_l, e12_l, w_l, wref_l = [], [], [], [], [], []
        ref_missing = 0
        for lbl in labels_all:
            arr = np.asarray(e_acc[lbl], dtype=np.float64)
            m = np.mean(arr, axis=0)
            w = vol_map.get(lbl, 1.0 if not vol_map else 0.0)
            if w <= 0.0:
                continue
            wr = ref_vol_map.get(lbl)
            if wr is None or not np.isfinite(wr) or wr <= 0.0:
                ref_missing += 1
                wr = w
            labels.append(lbl)
            e11_l.append(float(m[0]))
            e22_l.append(float(m[1]))
            e12_l.append(float(m[2]))
            w_l.append(float(w))
            wref_l.append(float(wr))
        if not e11_l:
            return _zk_fail('EVOL_zero_for_all_E_elements')

        e11 = np.asarray(e11_l, dtype=np.float64)
        e22 = np.asarray(e22_l, dtype=np.float64)
        e12 = np.asarray(e12_l, dtype=np.float64)
        w = np.asarray(w_l, dtype=np.float64)
        wref = np.asarray(wref_l, dtype=np.float64)
        wsum = float(np.sum(w))
        if not (wsum > 0.0):
            return _zk_fail('EVOL_zero_total')
        if ref_vol_map and ref_missing == 0:
            iw_weights = 'ref_evol'
        elif ref_vol_map and ref_missing < len(labels):
            iw_weights = 'ref_evol_partial'
        elif vol_map:
            iw_weights = 'evol_current'
        else:
            iw_weights = 'uniform'

        # ---- 주변형률 (local 성분; 고유값은 회전 불변) ----
        mean_e = 0.5 * (e11 + e22)
        rad = np.sqrt((0.5 * (e11 - e22)) ** 2 + e12 ** 2)
        e1 = mean_e + rad
        e2 = mean_e - rad
        #   Abaqus NLGEOM=ON 은 LE(대수 변형률)를 쓴다. 주성분은
        #   E = 0.5*(exp(2*LE)-1) 로 Green-Lagrange 로 정확히 변환된다(고유벡터 동일).
        strain_measure = 'E'
        strain_conversion = 'none'
        if strain_field == 'LE':
            e1 = 0.5 * np.expm1(2.0 * e1)
            e2 = 0.5 * np.expm1(2.0 * e2)
            strain_conversion = 'green_lagrange_from_log'

        mask_taut = (e2 > 0.0)
        mask_wrinkled = (e1 > 0.0) & (e2 <= 0.0)
        mask_slack = (e1 <= 0.0)
        wt = w / wsum
        frac_wr = float(np.sum(wt[mask_wrinkled]))
        frac_taut = float(np.sum(wt[mask_taut]))
        frac_slack = float(np.sum(wt[mask_slack]))

        if np.any(mask_wrinkled):
            vals = np.abs(e2[mask_wrinkled])
            wts = w[mask_wrinkled]
            denom = float(np.sum(wts))
            int_mean = float(np.sum(vals * wts) / denom) if denom > 0.0 else 0.0
            int_p95 = _zk_weighted_percentile(vals, wts, 95.0)
        else:
            int_mean = 0.0
            int_p95 = 0.0

        # ---- I_w (sidecar): 기준면적 가중 전체 압축 변형률 ----
        #   wrinkled: max(|E2|-nu*E1, 0)  (positive part; 준-taut 은 0)
        #   slack   : |E2| + nu*|E1|     (E1<=0, E2<=0)
        #   taut    : 0
        #   분모는 public EVOL 이 아니라 frameValue=0 EVOL 의 합(A_ref)이다.
        wr_comp = np.maximum(np.abs(e2) - ZK_NU * e1, 0.0)
        slack_comp = np.abs(e2) + ZK_NU * np.abs(e1)
        comp = np.zeros_like(e1)
        comp[mask_wrinkled] = wr_comp[mask_wrinkled]
        comp[mask_slack] = slack_comp[mask_slack]
        iw_den = float(np.sum(wref))
        iw_num = float(np.sum(wref * comp))
        iw_mean = (iw_num / iw_den) if iw_den > 0.0 else 0.0
        _pos = comp > 0.0
        iw_p95 = (_zk_weighted_percentile(comp[_pos], wref[_pos], 95.0)
                  if np.any(_pos) else 0.0)
        wr_num = float(np.sum(wref[mask_wrinkled] * wr_comp[mask_wrinkled]))
        slack_num = float(np.sum(wref[mask_slack] * slack_comp[mask_slack]))
        iw_slack_share = (slack_num / iw_num) if iw_num > 0.0 else 0.0
        iw_wr_share = (wr_num / iw_num) if iw_num > 0.0 else 0.0
        wr_area = float(np.sum(wref[mask_wrinkled]))
        clip_area = float(np.sum(wref[mask_wrinkled & (wr_comp <= 0.0)]))
        iw_clip_zero_frac = (clip_area / wr_area) if wr_area > 0.0 else 0.0

        # ---- 방향: local 1축 회전 후 axial mean (2*theta) ----
        dir_deg = float('nan')
        dir_spread = float('nan')
        dir_note = ''
        if np.any(mask_wrinkled):
            phi = _zk_element_local_angles(instance, labels)
            if phi is None:
                dir_note = 'local_axis_unavailable'
            else:
                theta_local = 0.5 * np.arctan2(2.0 * e12, e11 - e22)
                theta_g = theta_local + phi
                theta_g = ((theta_g + 0.5 * np.pi) % np.pi) - 0.5 * np.pi
                valid = mask_wrinkled & np.isfinite(theta_g)
                if np.any(valid):
                    ww = w[valid]
                    tt = theta_g[valid]
                    denom = float(np.sum(ww))
                    if denom > 0.0:
                        c = float(np.sum(ww * np.cos(2.0 * tt)))
                        s = float(np.sum(ww * np.sin(2.0 * tt)))
                        R = float(np.sqrt(c * c + s * s) / denom)
                        if R >= ZK_DIR_R_MIN:
                            dir_deg = float(np.degrees(0.5 * np.arctan2(s, c)))
                            dir_deg = ((dir_deg + 90.0) % 180.0) - 90.0
                        else:
                            dir_note = 'direction_degenerate_R=%.3e' % R
                        dir_spread = float(min(np.degrees(
                            0.5 * np.sqrt(max(-2.0 * np.log(max(R, 1e-12)), 0.0))), 90.0))
                else:
                    dir_note = 'local_axis_invalid'

        try:
            _nfr = int(len(odb.steps[step_name].frames))
        except Exception:
            _nfr = -1
        return {'ok': True, 'step': step_name,
                'frame': _nfr,
                'state_wrinkled_frac': frac_wr,
                'zk_intensity_mean': int_mean,
                'zk_intensity_p95': int_p95,
                'zk_dir_deg': dir_deg,
                'zk_state_taut_frac': frac_taut,
                'zk_state_slack_frac': frac_slack,
                'zk_dir_spread_deg': dir_spread,
                'zk_n_elem': int(len(labels)),
                'zk_n_wrinkled': int(np.count_nonzero(mask_wrinkled)),
                'zk_area_total': wsum,
                'weights': weight_mode,
                'strain_field': strain_field,
                'strain_measure': strain_measure,
                'strain_conversion': strain_conversion,
                'zk_Iw': iw_mean,
                'zk_Iw_p95': iw_p95,
                'zk_Iw_area_total': iw_den,
                'zk_Iw_weights': iw_weights,
                'zk_Iw_slack_share': iw_slack_share,
                'zk_Iw_wr_share': iw_wr_share,
                'zk_Iw_clip_zero_frac': iw_clip_zero_frac,
                'zk_Iw_ref_missing': int(ref_missing),
                'zk_Iw_nu': float(ZK_NU),
                'note': dir_note}
    except Exception as e:
        return _zk_fail('exception_%s' % type(e).__name__)
    finally:
        if odb is not None:
            try:
                odb.close()
            except Exception:
                pass


def _zk_fmt(x, fmt='%.10e'):
    try:
        xf = float(x)
        if not np.isfinite(xf):
            return 'nan'
        return fmt % xf
    except Exception:
        return 'nan'


def format_zk_line(zk):
    """ZK: 한 줄. RESULTS 3/4필드와 완전히 분리된 별도 계약이다."""
    if not zk or not zk.get('ok'):
        return 'ZK:FAIL reason=%s' % ((zk or {}).get('reason', 'unknown'))
    return ('ZK:state_wrinkled_frac=%s zk_intensity_mean=%s '
            'zk_intensity_p95=%s zk_dir_deg=%s'
            % (_zk_fmt(zk.get('state_wrinkled_frac'), '%.10e'),
               _zk_fmt(zk.get('zk_intensity_mean'), '%.10e'),
               _zk_fmt(zk.get('zk_intensity_p95'), '%.10e'),
               _zk_fmt(zk.get('zk_dir_deg'), '%.10e')))


def write_zk_metrics(zk, odb_path):
    """zk_metrics.txt (cwd) 기록. 실패해도 런을 죽이지 않는다."""
    try:
        lines = ['odb=%s' % os.path.basename(odb_path)]
        if not zk or not zk.get('ok'):
            lines.append('status=FAIL')
            lines.append('reason=%s' % ((zk or {}).get('reason', 'unknown')))
            for _k in ('fields', 'frame_value'):
                if _k in (zk or {}):
                    lines.append('%s=%s' % (_k, zk[_k]))
        else:
            lines.append('status=OK')
            lines.append('step=%s' % zk.get('step', '?'))
            lines.append('frame=%s' % zk.get('frame', '?'))
            for k in ('state_wrinkled_frac', 'zk_intensity_mean', 'zk_intensity_p95',
                      'zk_dir_deg', 'zk_state_taut_frac', 'zk_state_slack_frac',
                      'zk_dir_spread_deg', 'zk_n_elem', 'zk_n_wrinkled',
                      'zk_area_total', 'weights', 'note',
                      'strain_field', 'strain_measure', 'strain_conversion',
                      'zk_Iw', 'zk_Iw_p95', 'zk_Iw_area_total',
                      'zk_Iw_weights', 'zk_Iw_slack_share', 'zk_Iw_wr_share',
                      'zk_Iw_clip_zero_frac', 'zk_Iw_ref_missing', 'zk_Iw_nu'):
                if k in zk:
                    lines.append('%s=%s' % (k, zk[k]))
        with open('zk_metrics.txt', 'w') as f:
            f.write('\n'.join(lines) + '\n')
    except Exception as e:
        print('[ZK] zk_metrics.txt write failed: %s' % e)


def calc_thrust_loss(frame, instance, connectivity, max_node_label, area_ref=None,
                     flip_normals=True, extra_out=None):
    """
    메쉬의 각 요소 법선 벡터를 계산하고,
    SRP 모델을 적용하여 실제 추력 벡터를 적분, 이상적 추력과의 차이를 계산함.
    
    run_abaqus.py 에서 사각형+삼각형 메쉬의 조합으로 생성하기 때문에, 메쉬 형태에 따라 계산을 분기했음.

    Abaqus python api에서 모든 메쉬 전체에 대해 루프를 돌리면 계산시간이 메쉬의 수에 따라 크게 영향받으므로
    bulkDataBlocks를 통해 좌표 데이터를 통으로 가져온 뒤, Numpy 벡터연산으로 계산속도를 개선함
    """
    # 절점 좌표 추출 (Vectorized)
    coord_field = frame.fieldOutputs['COORD']

    node_coords = np.zeros((max_node_label + 1, 3), dtype=np.float64)
    
    # bulkDataBlocks를 순회하며 좌표 채우기
    for block in coord_field.bulkDataBlocks:
        if block.nodeLabels is not None:
            node_coords[block.nodeLabels] = block.data
    
    loss_sum = 0.0
    ideal_area_sum = 0.0

    # 4절점 요소 (Quad) 필터링
    quad_conns = [c for c in connectivity if len(c) == 4]
    if quad_conns:
        quads = np.array(quad_conns, dtype=np.int32)
        v0 = node_coords[quads[:, 0]]
        v1 = node_coords[quads[:, 1]]
        v2 = node_coords[quads[:, 2]]
        v3 = node_coords[quads[:, 3]]
        
        cross_prod = np.cross(v2 - v0, v3 - v1)
        norms = np.linalg.norm(cross_prod, axis=1)
        areas = 0.5 * norms
        
        normals = np.zeros_like(cross_prod)
        valid = norms > 1e-12
        normals[valid] = cross_prod[valid] / norms[valid, np.newaxis]
        
        #   [§19.90] 뒤집기 **전에** 진짜 역방향(n_z<0) 요소를 계수한다 — 뒤집기 편향 진단용.
        if extra_out is not None:
            _bw = normals[:, 2] < 0.0
            extra_out['n_total'] += int(len(areas))
            extra_out['n_backward'] += int(np.count_nonzero(_bw))
            extra_out['area_total'] += float(np.sum(areas))
            extra_out['area_backward'] += float(np.sum(areas[_bw]))
        if flip_normals:
            normals[normals[:, 2] < 0.0] *= -1.0
        
        # SRP 계산
        cos_theta = np.clip(np.dot(normals, -U_SUN), 0.0, 1.0)
        f_abs = A0 * cos_theta[:, np.newaxis] * areas[:, np.newaxis] * U_SUN
        f_refl = 2.0 * R0 * (cos_theta**2)[:, np.newaxis] * areas[:, np.newaxis] * (-normals)
        actual_thrust = np.dot(np.sum(f_abs + f_refl, axis=0), U_SUN)
        
        loss_sum += (np.sum(areas) * IDEAL_THRUST_COEFF - actual_thrust)
        ideal_area_sum += np.sum(areas)

    # 3절점 요소 (Tri) 필터링
    tri_conns = [c for c in connectivity if len(c) == 3]
    if tri_conns:
        tris = np.array(tri_conns, dtype=np.int32)
        v0 = node_coords[tris[:, 0]]
        v1 = node_coords[tris[:, 1]]
        v2 = node_coords[tris[:, 2]]
        
        cross_prod = np.cross(v1 - v0, v2 - v0)
        norms = np.linalg.norm(cross_prod, axis=1)
        areas = 0.5 * norms
        
        normals = np.zeros_like(cross_prod)
        valid = norms > 1e-12
        normals[valid] = cross_prod[valid] / norms[valid, np.newaxis]
        #   [§19.90] 뒤집기 **전에** 진짜 역방향(n_z<0) 요소를 계수한다 — 뒤집기 편향 진단용.
        if extra_out is not None:
            _bw = normals[:, 2] < 0.0
            extra_out['n_total'] += int(len(areas))
            extra_out['n_backward'] += int(np.count_nonzero(_bw))
            extra_out['area_total'] += float(np.sum(areas))
            extra_out['area_backward'] += float(np.sum(areas[_bw]))
        if flip_normals:
            normals[normals[:, 2] < 0.0] *= -1.0
        
        cos_theta = np.clip(np.dot(normals, -U_SUN), 0.0, 1.0)
        f_abs = A0 * cos_theta[:, np.newaxis] * areas[:, np.newaxis] * U_SUN
        f_refl = 2.0 * R0 * (cos_theta**2)[:, np.newaxis] * areas[:, np.newaxis] * (-normals)
        actual_thrust = np.dot(np.sum(f_abs + f_refl, axis=0), U_SUN)
        
        loss_sum += (np.sum(areas) * IDEAL_THRUST_COEFF - actual_thrust)
        ideal_area_sum += np.sum(areas)

    #   [§19.82] **상대** 손실로 반환한다.
    #     종전 반환은 `loss_sum` = Σ(면적×IDEAL_THRUST_COEFF) − 실제추력 = **절대량**이었고,
    #     정규화용으로 누적한 `ideal_area_sum` 을 쓰지 않았다(의도는 상대화였는데 빠져 있었다).
    #     절대량을 비율로 오독하면 "0.00755 = 0.7%" 처럼 100배 작게 보인다.
    #     상대 = 손실 / (총면적 × IDEAL_THRUST_COEFF) = 평판 이상추력 대비 손실률.
    #   [§19.89] 분모는 **미변형(pristine) 기준 면적**을 쓴다(사용자 지시).
    #     종전 `ideal_area_sum` 은 *변형 후* 총면적이라 설계점마다 기준이 미세하게 달라진다.
    #     원형 평면적을 기준으로 고정하면 `hf` 가 "원형 평면 대비 성능 손실률"로 정의되어
    #     설계 간 비교가 정합해진다. area_ref 를 못 구하면 종전 동작으로 폴백한다.
    _den_area = area_ref if (area_ref and area_ref > 0.0) else ideal_area_sum
    if _den_area > 0.0 and IDEAL_THRUST_COEFF > 0.0:
        return float(max(loss_sum, 0.0) / (_den_area * IDEAL_THRUST_COEFF))
    return float(max(loss_sum, 1e-12))

def get_hf(odb_path):
    """
    Step-Postbuckle의 마지막 프레임(최종 목표 변위 도달 시점)에서
    추력 손실(Thrust Loss)을 계산하여 반환함. (포스트버클링이 목표 변위에 도달하지 못하는 경우를 대비)
    """
    try:
        odb = openOdb(path=odb_path, readOnly=True)
        
        if 'Step-Postbuckle' not in odb.steps:
            print("Error: 'Step-Postbuckle' not found in ODB.")
            odb.close()
            return None
            
        step = odb.steps['Step-Postbuckle']
        
        last_frame = step.frames[-1]
        current_time = last_frame.frameValue
        
        if current_time < 0.99:
            print("Warning: Job did not complete (Time = {:.4f}).".format(current_time))
            odb.close()       # P1-1
            return None
            # 수렴 실패 시 페널티 값을 리턴하거나, 현재 상태라도 계산할지 결정
        
        instance = odb.rootAssembly.instances['MEMBRANE-1']
        
        # Element Connectivity 캐싱
        connectivity = [e.connectivity for e in instance.elements]
        
        # 배열 크기를 할당하기 위해 최대 노드 라벨 탐색
        if instance.nodes:
            max_label = instance.nodes[-1].label + 100
        else:
            max_label = 100000

        #   [§19.89] 미변형 기준 면적: 첫 스텝의 frameValue=0 프레임(원형 메쉬 기하)에서 계산.
        _area_ref = None
        try:
            _f0 = None
            for _sn in odb.steps.keys():
                for _fr in odb.steps[_sn].frames:
                    if abs(_fr.frameValue) < 1e-12:
                        _f0 = _fr
                        break
                if _f0 is not None:
                    break
            if _f0 is not None and 'COORD' in _f0.fieldOutputs.keys():
                _cf = _f0.fieldOutputs['COORD'].getSubset(region=instance)
                _nc = np.zeros((max_label + 1, 3), dtype=np.float64)
                for _b in _cf.bulkDataBlocks:
                    if _b.nodeLabels is not None:
                        _nc[_b.nodeLabels] = _b.data
                _sum = 0.0
                _q = [c for c in connectivity if len(c) == 4]
                if _q:
                    _a = np.array(_q, dtype=np.int32)
                    _sum += 0.5 * np.linalg.norm(
                        np.cross(_nc[_a[:, 2]] - _nc[_a[:, 0]],
                                 _nc[_a[:, 3]] - _nc[_a[:, 1]]), axis=1).sum()
                _t = [c for c in connectivity if len(c) == 3]
                if _t:
                    _a = np.array(_t, dtype=np.int32)
                    _sum += 0.5 * np.linalg.norm(
                        np.cross(_nc[_a[:, 1]] - _nc[_a[:, 0]],
                                 _nc[_a[:, 2]] - _nc[_a[:, 0]]), axis=1).sum()
                if _sum > 0.0:
                    _area_ref = float(_sum)
                    print("[HF] undeformed area (pristine) = %.4f m^2" % _area_ref)
        except Exception as _are:
            print("[HF] undeformed area probe failed (%s) -> deformed area fallback" % _are)
        #   [§19.90] 뒤집기 편향 실측: 같은 프레임을 뒤집기 ON/OFF 두 번 적분한다.
        #     진짜 n_z<0 요소의 면적비가 편향을 결정한다 → 가정하지 않고 매 런 기록.
        _st = {'n_total': 0, 'n_backward': 0, 'area_total': 0.0, 'area_backward': 0.0}
        hf_metric = calc_thrust_loss(last_frame, instance, connectivity, max_label,
                                     area_ref=_area_ref, extra_out=_st)
        hf_true = calc_thrust_loss(last_frame, instance, connectivity, max_label,
                                   area_ref=_area_ref, flip_normals=False)
        _fb = (100.0 * _st['area_backward'] / _st['area_total']
               if _st['area_total'] > 0.0 else float('nan'))
        _dl = (100.0 * (hf_metric - hf_true) / hf_true) if hf_true else float('nan')
        print("[HF] flip bias: backward area=%.4f%% (%d/%d elem) hf_flip=%.6e "
              "hf_true=%.6e delta=%.3f%%"
              % (_fb, _st['n_backward'], _st['n_total'], hf_metric, hf_true, _dl))
        try:
            with open('flip_bias.txt', 'w') as _bf:
                _bf.write("backward_area_pct=%.6f\n" % _fb)
                _bf.write("backward_elems=%d\n" % _st['n_backward'])
                _bf.write("total_elems=%d\n" % _st['n_total'])
                _bf.write("hf_flip=%.10e\n" % hf_metric)
                _bf.write("hf_true_normal=%.10e\n" % hf_true)
                _bf.write("delta_pct=%.6f\n" % _dl)
        except Exception as _be:
            print("[HF] flip_bias.txt write failed: %s" % _be)
        
        odb.close()
        return float(hf_metric)

    except Exception as e:
        print("HF Extraction Error: " + str(e))
        try: odb.close()      # P1-1
        except Exception: pass
        return None

# ---------------------------------------------------------------------------
#   [19.106] 형상 지표: 면적가중 등가법선 편차 delta_n (및 방위 delta_psi)
#
#   왜: 참조 문헌(Pino 2019 ASR 63(11):3675, Vulpetti 2021 ASR 67(9):2669)은 형상 영향을
#       '등가법선 편차 delta_n' 하나로 정량화한다(6 deg / 7.774 deg). 우리 지표(lf1..lf3, hf)
#       중 이 값만 없었다.
#
#   정의(이상 법선 = (0,0,1), U_SUN=(0,0,-1) 기준 면):
#       n_eq   = normalize( sum_e A_e * n_e )      (A_e 요소 면적, n_e 요소 단위법선)
#       delta_n   [deg] = acos( clip(n_eq . (0,0,1), -1, 1) )
#       delta_psi [deg] = atan2(n_eq_y, n_eq_x)
#
#   성질: **일관 기울기(whole-surface bias)만** 잡고 국소 주름은 면적가중에서 상쇄된다.
#         (앵커 (c): 평균경사 0 인 대칭 사인파 -> delta_n = 0. hf 는 0 이 아니다.)
#         즉 delta_n 은 hf 와 **독립**인 정보(지향 오차)를 준다.
#
#   출력 계약: RESULTS:lf1,lf2,lf3,hf 4필드는 **절대 건드리지 않는다**(아카이빙/MFBO 무효화 방지).
#             별도 한 줄  "DN:dn_deg=.. dn_psi_deg=.. n_eq_x=.. n_eq_y=.. n_eq_z=.. n_elem=.. area_ref=.."
#             계산 불가면 "DN:FAIL reason=.." 한 줄만 내고 나머지 경로는 평소대로 진행한다.
# ---------------------------------------------------------------------------
DN_IDEAL_NORMAL = np.array([0.0, 0.0, 1.0], dtype=np.float64)


def _dn_fail(reason):
    """DN 실패 사유를 ASCII 토큰 하나로 정규화(공백 금지 -> key=value 파싱 보존)."""
    try:
        _s = str(reason)
    except Exception:
        _s = 'unknown'
    _s = ''.join(ch if (33 <= ord(ch) < 127) else '_' for ch in _s)
    return {'ok': False, 'reason': (_s or 'unknown')}


def _dn_node_coords(frame, instance, max_node_label):
    """프레임의 COORD 를 (label, 3) 배열로 전개한다(calc_thrust_loss 와 같은 방식)."""
    coord_field = frame.fieldOutputs['COORD']
    try:
        coord_field = coord_field.getSubset(region=instance)
    except Exception:
        pass
    node_coords = np.zeros((max_node_label + 1, 3), dtype=np.float64)
    for block in coord_field.bulkDataBlocks:
        if block.nodeLabels is not None:
            node_coords[block.nodeLabels] = block.data
    return node_coords


def calc_normal_deviation(frame, instance, connectivity, max_node_label):
    """[19.106] 면적가중 등가법선 n_eq 와 편차 delta_n / 방위 delta_psi 를 계산한다.

    기존 calc_thrust_loss 와 **같은** 요소 면적/법선 규약을 쓴다(사각형은 대각선 외적,
    삼각형은 두 변 외적, n_z<0 은 뒤집어 +z 면을 본다). 새 함수이므로 기존 계산은 불변.

    Returns dict:
        ok=True  -> dn_deg, dn_psi_deg, n_eq_x/y/z, n_elem, area_ref
        ok=False -> reason
    """
    try:
        nc = _dn_node_coords(frame, instance, max_node_label)

        areas_l = []
        normals_l = []

        quad_conns = [c for c in connectivity if len(c) == 4]
        if quad_conns:
            q = np.array(quad_conns, dtype=np.int32)
            cp = np.cross(nc[q[:, 2]] - nc[q[:, 0]], nc[q[:, 3]] - nc[q[:, 1]])
            nrm = np.linalg.norm(cp, axis=1)
            ok = nrm > 1e-12
            areas_l.append(0.5 * nrm[ok])
            normals_l.append(cp[ok] / nrm[ok, np.newaxis])

        tri_conns = [c for c in connectivity if len(c) == 3]
        if tri_conns:
            t = np.array(tri_conns, dtype=np.int32)
            cp = np.cross(nc[t[:, 1]] - nc[t[:, 0]], nc[t[:, 2]] - nc[t[:, 0]])
            nrm = np.linalg.norm(cp, axis=1)
            ok = nrm > 1e-12
            areas_l.append(0.5 * nrm[ok])
            normals_l.append(cp[ok] / nrm[ok, np.newaxis])

        if not areas_l:
            return _dn_fail('no_supported_element')
        areas = np.concatenate(areas_l)
        normals = np.concatenate(normals_l)
        if areas.size == 0:
            return _dn_fail('zero_area_all_elements')

        area_total = float(np.sum(areas))
        if not (area_total > 0.0):
            return _dn_fail('zero_area_total')

        #   +z 면을 보도록 정렬(hf 계산의 flip_normals=True 와 동일 규약).
        flip = normals[:, 2] < 0.0
        if np.any(flip):
            normals[flip] = -normals[flip]

        acc = np.sum(areas[:, np.newaxis] * normals, axis=0)
        acc_norm = float(np.linalg.norm(acc))
        if not (acc_norm > 0.0):
            return _dn_fail('degenerate_normal_sum')
        n_eq = acc / acc_norm

        cos_dn = float(np.clip(float(np.dot(n_eq, DN_IDEAL_NORMAL)), -1.0, 1.0))
        dn_deg = float(np.degrees(np.arccos(cos_dn)))
        dn_psi_deg = float(np.degrees(np.arctan2(float(n_eq[1]), float(n_eq[0]))))

        return {'ok': True,
                'dn_deg': dn_deg,
                'dn_psi_deg': dn_psi_deg,
                'n_eq_x': float(n_eq[0]),
                'n_eq_y': float(n_eq[1]),
                'n_eq_z': float(n_eq[2]),
                'n_elem': int(areas.size),
                'area_ref': area_total}
    except Exception as e:
        return _dn_fail('exception_' + type(e).__name__)


def get_dn(odb_path):
    """[19.106] HF ODB 의 최종 Postbuckle 프레임에서 delta_n 을 계산한다(get_hf 와 같은 프레임).

    실패해도 예외를 던지지 않는다(dict(ok=False, reason=..) 반환) - 호출부가 DN:FAIL 한 줄만 낸다.
    """
    odb = None
    try:
        if not odb_path or not os.path.exists(odb_path):
            return _dn_fail('file_not_found')
        odb = openOdb(path=odb_path, readOnly=True)
        if 'Step-Postbuckle' not in odb.steps:
            return _dn_fail('step_postbuckle_missing')
        step = odb.steps['Step-Postbuckle']
        if not step.frames:
            return _dn_fail('no_frame')
        instance = odb.rootAssembly.instances['MEMBRANE-1']
        connectivity = [e.connectivity for e in instance.elements]
        if not connectivity:
            return _dn_fail('no_element')
        max_label = (instance.nodes[-1].label + 100) if instance.nodes else 100000
        return calc_normal_deviation(step.frames[-1], instance, connectivity, max_label)
    except Exception as e:
        return _dn_fail('exception_' + type(e).__name__)
    finally:
        if odb is not None:
            try: odb.close()
            except Exception: pass


def format_dn_line(dn):
    """[19.106] DN: 한 줄 포맷. RESULTS: 4필드 줄과 **완전히 분리**된 별도 줄이다."""
    if not dn or not dn.get('ok'):
        _r = (dn or {}).get('reason', 'unknown')
        return 'DN:FAIL reason=%s' % _r
    return ('DN:dn_deg=%.10e dn_psi_deg=%.10e n_eq_x=%.10e n_eq_y=%.10e '
            'n_eq_z=%.10e n_elem=%d area_ref=%.10e'
            % (dn['dn_deg'], dn['dn_psi_deg'], dn['n_eq_x'], dn['n_eq_y'],
               dn['n_eq_z'], dn['n_elem'], dn['area_ref']))


def main():
    # 인자 파싱 루틴
    args = sys.argv
    # sys.argv 예시: ['eval_abaqus.py', 'LF.odb', 'LF'] 또는 ['eval_abaqus.py', 'LF.odb', 'HF.odb', 'HF']
    
    try:
        mode = args[-1].upper()
        if mode == "LF":
            lf_metrics = get_lf(args[-2])
            if any(v is None for v in lf_metrics):
                print("RESULTS:FAIL"); sys.exit(1)   # E-7: 0 은 실패를 숨긴다
            results = list(lf_metrics[:3])
            #   [paper-match] ZK 는 별도 줄/파일로만 낸다(RESULTS 3필드 계약 불변).
            #     E 가 없어도 ZK:FAIL 만 내고 LF 결과는 정상적으로 낸다.
            _zk = get_zk(args[-2])
            print(format_zk_line(_zk))
            write_zk_metrics(_zk, args[-2])
        elif mode == "HFONLY":
            #   [§19.87] HF 단독 평가: 짝 LF ODB 가 없어도 **hf 지표만** 낸다.
            #     get_hf 는 HF ODB 의 Postbuckle 프레임만 읽으므로 LF 없이 성립한다.
            #     LF 자리에는 nan 을 채운다(LF 지표는 그 점에서 정의되지 않음).
            hf_metric = get_hf(args[-2])
            if hf_metric is None:
                print("HF Error: hf metric extraction failed")
                return 1
            #   [19.106] delta_n 은 별도 줄로만 낸다(RESULTS 4필드 계약 불변).
            print(format_dn_line(get_dn(args[-2])))
            with open('extraction.txt', 'w') as _f:
                _f.write("nan,nan,nan,%s" % hf_metric)
            print("RESULTS:nan,nan,nan,%s" % hf_metric)
            return 0
        elif mode == "HF":
            hf_metric = get_hf(args[-2])
            lf_metrics = get_lf(args[-3])
            if hf_metric is None or any(v is None for v in lf_metrics):
                print("RESULTS:FAIL"); sys.exit(1)   # E-7: 0 은 실패를 숨긴다
            #   [paper-match] ZK 는 **짝 LF ODB** 의 최종 상태에서 계산한다.
            #     HF 의 hf/DN 은 그대로이고, ZK 는 별도 줄/파일로만 나간다.
            _zk = get_zk(args[-3])
            print(format_zk_line(_zk))
            write_zk_metrics(_zk, args[-3])
            #   [19.106] delta_n 은 별도 줄로만 낸다(RESULTS 4필드 계약 불변).
            print(format_dn_line(get_dn(args[-2])))
            results = list(lf_metrics[:3]) + [hf_metric]
        else:
            # [2026-10-05] 인자 형식 오류를 **조용히 FAIL 로 뭉개지 않는다.**
            #   실측: `eval_abaqus.py A.odb B.odb` 로 두 ODB 를 비교하려다 'HF_X025I.ODB'
            #   가 mode 가 되어 여기로 떨어졌고, 사용자는 원인을 알 수 없었다.
            #   이 스크립트는 두 ODB 비교기가 아니라 **단일 ODB 지표 추출기**다.
            print("ERROR: 마지막 인자가 HF 또는 LF 여야 합니다. 받은 값: %r" % args[-1])
            print("사용법:")
            print("  abaqus python eval_abaqus.py <LF.odb> LF")
            print("  abaqus python eval_abaqus.py <LF.odb> <HF.odb> HF")
            print("  (두 ODB 를 비교하려면 compare_odb_fields.py 를 쓰세요.)")
            print("RESULTS:FAIL"); sys.exit(1)
    except Exception as e:
        print("eval error:", e)
        print("RESULTS:FAIL"); sys.exit(1)   # E-7: 0 은 실패를 숨긴다

    output_str = ",".join(map(str, results))
    print("RESULTS:" + output_str)
    
    with open('extraction.txt', 'w') as f:
        # LF : lf1, lf2, lf3
        # HF : lf1, lf2, lf3, hf
        full_data = list(lf_metrics)
        if mode == "HF": full_data.append(hf_metric)
        f.write(",".join(map(str, full_data)))

if __name__ == "__main__":
    main()
    