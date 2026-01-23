"""
아바쿠스 해석이 종료된 후, 그 결과 파일인 .odb 를 읽어서 필요한 lf, hf값을 계산하는 코드
run_abaqus.py 코드 내에서 호출되어 사용되며, 
"abaqus python eval_abaqus.py (lf_odb파일) (hf_odb파일) HF" 형태로 호출

Return Format (stdout):
    "RESULTS:val1,val2,..." 형식으로 출력하여 상위 프로세스가 파싱할 수 있게 함.
    - LF Mode: "RESULTS:lf1,lf2,lf3"
    - HF Mode: "RESULTS:lf1,lf2,lf3,hf" 

LF : 선형 해석에서 얻을 수 있는 지표들
    - lf1 (Compressed Volume Ratio): 압축 주응력이 발생하는 요소의 부피 비율.
    - lf2 (Wrinkling Energy Density): 인장-압축 혼합 응력 상태에서의 변형 에너지 밀도.
    - lf3 (Stress Anisotropy): 주응력 간의 비대칭성 비율.

HF :
    - hf (Thrust Loss): 이상적인 평판 대비, 주름진 형상으로 인해 손실된 추력의 양.
    - Solar Radiation Pressure (SRP) 모델을 적용하여 변형된 요소들의 법선 벡터 변화를 적분.
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

def get_lf(odb_path):
    """
    선형 해석 결과(LF ODB)에서 주름 발생 가능성을 나타내는 대리 지표(Proxy Metrics) 3종을 추출.
    
    Returns:
        list: [lf1, lf2, lf3]
    """
    try:
        if not os.path.exists(odb_path):
            print("LF Error: File not found (%s)" % odb_path)
            return [0.0, 0.0, 0.0]

        odb = openOdb(path=odb_path, readOnly=True)
        
        if 'Step-HighTension' in odb.steps:
            step_name = 'Step-HighTension'
        else:
            print("LF Error: Target step not found.")
            return [0.0, 0.0, 0.0]

        frame = odb.steps[step_name].frames[-1]
        instance = odb.rootAssembly.instances['MEMBRANE-1'] 
        
        if 'EVOL' not in frame.fieldOutputs.keys():
            print("LF Error: EVOL field not found.")
            return [0.0, 0.0, 0.0]
            
        evol_field = frame.fieldOutputs['EVOL'].getSubset(region=instance)
        
        vol_map = {}
        for val in evol_field.values:
            vol_map[val.elementLabel] = val.data
            
        total_vol = sum(vol_map.values())
        if total_vol == 0: return [0.0, 0.0, 0.0]
        
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
        return [1e6, -1e6, -1e6]

def calc_thrust_loss(frame, instance, connectivity, max_node_label):
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
        normals[normals[:, 2] < 0.0] *= -1.0
        
        cos_theta = np.clip(np.dot(normals, -U_SUN), 0.0, 1.0)
        f_abs = A0 * cos_theta[:, np.newaxis] * areas[:, np.newaxis] * U_SUN
        f_refl = 2.0 * R0 * (cos_theta**2)[:, np.newaxis] * areas[:, np.newaxis] * (-normals)
        actual_thrust = np.dot(np.sum(f_abs + f_refl, axis=0), U_SUN)
        
        loss_sum += (np.sum(areas) * IDEAL_THRUST_COEFF - actual_thrust)
        ideal_area_sum += np.sum(areas)

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
            return 1e6
            
        step = odb.steps['Step-Postbuckle']
        
        last_frame = step.frames[-1]
        current_time = last_frame.frameValue
        
        if current_time < 0.99:
            print("Warning: Job did not complete (Time = {:.4f}).".format(current_time))
            # 수렴 실패 시 페널티 값을 리턴하거나, 현재 상태라도 계산할지 결정
            # 일단 계산은 하되 경고를 띄우는 방식
        
        instance = odb.rootAssembly.instances['MEMBRANE-1']
        
        # Element Connectivity 캐싱
        connectivity = [e.connectivity for e in instance.elements]
        
        # 배열 크기를 할당하기 위해 최대 노드 라벨 탐색
        if instance.nodes:
            max_label = instance.nodes[-1].label + 100
        else:
            max_label = 100000

        hf_metric = calc_thrust_loss(last_frame, instance, connectivity, max_label)
        
        odb.close()
        return float(hf_metric)

    except Exception as e:
        print("HF Extraction Error: " + str(e))
        return 1e6 # 에러 발생시 아주 큰 값 < 최소화 문제에서 페널티

def main():
    # 인자 파싱 루틴
    args = sys.argv
    # sys.argv 예시: ['eval_abaqus.py', 'LF.odb', 'LF'] 또는 ['eval_abaqus.py', 'LF.odb', 'HF.odb', 'HF']
    
    try:
        mode = args[-1].upper()
        if mode == "LF":
            path_lf = args[-2]
            # 모든 LF 지표(lf1, lf2, lf3)를 계산함
            lf_metrics = get_lf(path_lf)
            results = [lf_metrics[0], lf_metrics[1], lf_metrics[2]] 
            
        elif mode == "HF":
            path_hf = args[-2]
            path_lf = args[-3]
            
            lf_metrics = get_lf(path_lf)
            hf_metric = get_hf(path_hf) 
            # Trace-aware를 위해 lf, hf를 모두 반환
            results = [lf_metrics[0], lf_metrics[1], lf_metrics[2], hf_metric]
        else:
            sys.exit(1)
    except:
        sys.exit(1)

    output_str = ",".join(map(str, results))
    print("RESULTS:" + output_str)
    
    with open('extraction_full.txt', 'w') as f:
        # LF : lf1, lf2, lf3
        # HF : lf1, lf2, lf3, hf
        full_data = list(lf_metrics)
        if mode == "HF": full_data.append(hf_metric)
        f.write(",".join(map(str, full_data)))

if __name__ == "__main__":
    main()
    