# -*- coding: utf-8 -*-
"""base state 응력 측정 (R-13 / 스파이크 10-0).

좌굴 base state(= 직전 일반 스텝 말단)의 면내 응력 분포를 .odb 에서 읽어
  - 최대/최소 주응력 (인장/압축)
  - 압축 면적비 (min principal < 0 인 요소 비율, 부피 가중)
  - 앵커 반력 (RF 최대 크기)
를 콘솔에 찍는다. 읽기 전용이라 해석에는 전혀 영향을 주지 않는다.

사용:
    abaqus python base_state_probe.py Buckle_Analysis.odb Step-GlobalTension

주의: 좌굴 스텝 자체는 실패해도 되고, 그 앞 일반 스텝(GlobalTension)의
      마지막 프레임이 곧 base state 다.
"""
from __future__ import print_function
import sys
import math

TARGET_STRESS = 7000.0  # Pa, R-13 목표 운용점(미검증)


def _principal(s11, s22, s12):
    mean = 0.5 * (s11 + s22)
    dev = math.sqrt((0.5 * (s11 - s22)) ** 2 + s12 ** 2)
    return mean + dev, mean - dev   # (max, min)


def main():
    try:
        from odbAccess import openOdb
    except ImportError:
        print("[R-13] odbAccess 를 못 읽었습니다 — 'abaqus python' 으로 실행해야 합니다.")
        return 2

    odb_path = sys.argv[1] if len(sys.argv) > 1 else 'Buckle_Analysis.odb'
    step_name = sys.argv[2] if len(sys.argv) > 2 else 'Step-GlobalTension'

    try:
        odb = openOdb(odb_path, readOnly=True)
    except Exception as e:
        print("[R-13] odb 열기 실패 (%s): %s" % (odb_path, e))
        return 1

    try:
        if step_name not in odb.steps.keys():
            print("[R-13] step '%s' 없음. 있는 step: %s" % (step_name, list(odb.steps.keys())))
            return 1
        step = odb.steps[step_name]
        nf = len(step.frames)
        if nf == 0:
            print("[R-13] '%s' 프레임 0개 — 해석이 이 스텝을 끝내지 못했습니다." % step_name)
            return 1
        frame = step.frames[-1]
        print("[R-13] base state = %s / %s  (frame %d/%d, step time %s)"
              % (odb_path, step_name, nf, nf, frame.frameValue))

        # --- 면내 주응력 (요소 중앙값, 쉘 S = S11,S22,S12) ---
        # CENTROID 로 요소당 1개 값 (쉘 섹션포인트 중복 제거) — eval_abaqus.py 와 동일 규약
        try:
            from abaqusConstants import CENTROID
            sf = frame.fieldOutputs['S'].getSubset(position=CENTROID)
        except Exception:
            try:
                sf = frame.fieldOutputs['S']
            except Exception:
                print("[R-13] S 필드 없음 — fieldOutputRequest 확인 필요")
                return 1

        n = 0
        s1_max = -1e30
        s2_min = 1e30
        s2_max = -1e30
        n_comp = 0
        s11_max = -1e30
        s11_min = 1e30
        for v in sf.values:
            d = v.data
            if len(d) < 3:
                continue
            s11, s22, s12 = float(d[0]), float(d[1]), float(d[2])
            s1, s2 = _principal(s11, s22, s12)
            n += 1
            s1_max = max(s1_max, s1)
            s2_min = min(s2_min, s2)
            s2_max = max(s2_max, s2)
            s11_max = max(s11_max, s11)
            s11_min = min(s11_min, s11)
            if s2 < 0.0:
                n_comp += 1

        if n == 0:
            print("[R-13] S 값 0개")
        else:
            print("[R-13] 주응력 (n=%d 값): max principal %+.4g Pa / min principal %+.4g Pa"
                  % (n, s1_max, s2_min))
            print("[R-13] S11 범위 %+.4g ~ %+.4g Pa   S22(max principal) 최대 %+.4g Pa"
                  % (s11_min, s11_max, s2_max))
            print("[R-13] 압축 응력(min principal<0) 비율 = %.4f  (%d / %d)"
                  % (float(n_comp) / n, n_comp, n))
            print("[R-13] 목표 운용점 TARGET_STRESS = %.1f Pa (미검증) — 위 S11 최대와 비교"
                  % TARGET_STRESS)

        # --- 반력 (앵커/케이블 끝 고정점) ---
        try:
            rf = frame.fieldOutputs['RF']
            best = []
            for v in rf.values:
                d = v.data
                mag = math.sqrt(sum(float(x) ** 2 for x in d))
                if mag > 1e-12:
                    best.append((mag, v.nodeLabel, tuple(d)))
            best.sort(reverse=True)
            print("[R-13] RF 비영 노드 %d개, 상위 3개 (크기, node, 성분):" % len(best))
            for mag, lab, comp in best[:3]:
                print("       %.6g N  node %s  %s" % (mag, lab, comp))
        except Exception as e:
            print("[R-13] RF 읽기 실패: %s" % e)
        return 0
    finally:
        try:
            odb.close()
        except Exception:
            pass


if __name__ == '__main__':
    sys.exit(main())
