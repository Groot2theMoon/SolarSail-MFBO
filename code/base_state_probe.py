# -*- coding: utf-8 -*-
"""base state 응력 측정 v2 (R-13 / 스파이크 10-0).

좌굴 base state(= 직전 일반 스텝 말단)의 면내 응력 분포를 .odb 에서 읽어

  1) 단면점(section point) 별 통계 — 압축이 '진짜 면내 압축'인지 굽힘(표면) artifact 인지 분리
  2) 면적 가중(요소 EVOL) 통계        — 요소 개수 비율이 아니라 실제 면적 비율
  3) 면내 평균응력 (S11+S22)/2 < 0 비율 — 단면점에 무관한 압축 지표
  4) 앵커 반력 (RF, instance 이름 포함 = 케이블 장력)

읽기 전용. 해석에는 전혀 영향을 주지 않는다.

사용:
    abaqus python base_state_probe.py Buckle_Analysis.odb Step-GlobalTension

판독:
  * 면내 평균응력<0 면적비 ~0 : base state 인장 지배 -> 좌굴 추출 실패는 '수치' 문제
  * 0.2 이상                  : base state 가 slack/주름 상태 ->
                                Abaqus 문서상 '예압이 좌굴하중 위' 케이스
                                (subspace 는 수렴 실패, Lanczos 는 사용 금지)
"""
from __future__ import print_function
import sys
import math

TARGET_STRESS = 7000.0   # Pa, R-13 목표 운용점(미검증)


def _principal(s11, s22, s12):
    mean = 0.5 * (s11 + s22)
    dev = math.sqrt((0.5 * (s11 - s22)) ** 2 + s12 ** 2)
    return mean + dev, mean - dev      # (max, min)


def _aggregate(rows):
    """rows: list of (s11,s22,s12,w). 면적 가중 통계를 돌려준다."""
    wsum = sum(r[3] for r in rows)
    if wsum <= 0.0:                       # EVOL 이 없으면 균등 가중으로 대체
        rows = [(r[0], r[1], r[2], 1.0) for r in rows]
        wsum = float(len(rows))
    w_comp_mean = 0.0     # (S11+S22)/2 < 0
    w_comp_min = 0.0      # min principal < 0
    s1mx = -1e30
    s2mn = 1e30
    mean_acc = 0.0
    for s11, s22, s12, w in rows:
        s1, s2 = _principal(s11, s22, s12)
        m = 0.5 * (s11 + s22)
        if m < 0.0:
            w_comp_mean += w
        if s2 < 0.0:
            w_comp_min += w
        s1mx = max(s1mx, s1)
        s2mn = min(s2mn, s2)
        mean_acc += m * w
    return (len(rows), w_comp_mean / wsum, w_comp_min / wsum, mean_acc / wsum, s1mx, s2mn)


def _fmt(lab, n, cm, cmin, mavg, s1mx, s2mn):
    return ("  %-20s n=%6d  면적비(mean<0)=%6.4f  면적비(minP<0)=%6.4f  "
            "평균면내=%+10.4g  maxP=%+10.4g  minP=%+10.4g"
            % (lab, n, cm, cmin, mavg, s1mx, s2mn))


def main():
    try:
        from odbAccess import openOdb
    except ImportError:
        print("[R-13] odbAccess 를 못 읽었습니다 - 'abaqus python' 으로 실행해야 합니다.")
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
            print("[R-13] '%s' 프레임 0개 - 해석이 이 스텝을 끝내지 못했습니다." % step_name)
            return 1
        frame = step.frames[-1]
        print("[R-13] base state = %s / %s  (frame %d/%d, step time %s)"
              % (odb_path, step_name, nf, nf, frame.frameValue))

        # 요소 체적 (면적 가중용; 쉘은 두께 균일 -> 면적비와 동일)
        vols = {}
        try:
            for v in frame.fieldOutputs['EVOL'].values:
                d = v.data
                vols[v.elementLabel] = float(d[0]) if isinstance(d, (tuple, list)) else float(d)
            print("[R-13] EVOL 요소 %d개 확보 (면적 가중 사용)" % len(vols))
        except Exception as e:
            print("[R-13] EVOL 없음 -> 균등 가중으로 대체 (%s)" % e)

        # CENTROID 로 요소당 1개 값 (쉘은 단면점마다 1개)
        try:
            from abaqusConstants import CENTROID
            sf = frame.fieldOutputs['S'].getSubset(position=CENTROID)
        except Exception:
            sf = frame.fieldOutputs['S']

        groups = {}      # section number -> rows
        for v in sf.values:
            d = v.data
            if len(d) < 3:
                continue
            try:
                num = v.sectionPoint.number
            except Exception:
                num = None
            w = vols.get(v.elementLabel, 0.0)
            groups.setdefault(num, []).append(
                (float(d[0]), float(d[1]), float(d[2]), w if w > 0 else 0.0))

        if not groups:
            print("[R-13] S 값 0개")
            return 1

        key = lambda x: (x is None, x)
        print("[R-13] 단면점 %d개: %s" % (len(groups), sorted(groups.keys(), key=key)))
        for num in sorted(groups.keys(), key=key):
            n, cm, cmin, mavg, s1mx, s2mn = _aggregate(groups[num])
            print(_fmt('sec %s' % ('-' if num is None else num), n, cm, cmin, mavg, s1mx, s2mn))

        allrows = [r for rows in groups.values() for r in rows]
        n, cm, cmin, mavg, s1mx, s2mn = _aggregate(allrows)
        print(_fmt('ALL', n, cm, cmin, mavg, s1mx, s2mn))

        nums = sorted([k for k in groups.keys() if k is not None])
        mid = nums[len(nums) // 2] if nums else None
        if mid is not None and len(nums) > 1:
            n, cm, cmin, mavg, s1mx, s2mn = _aggregate(groups[mid])
            print(_fmt('MIDSECTION(%s)' % mid, n, cm, cmin, mavg, s1mx, s2mn))
            print("[R-13] >>> 헤드라인(중앙단면 면적가중): 면내평균<0 비율=%.4f  "
                  "min principal<0 비율=%.4f  평균 면내응력=%+.4g Pa" % (cm, cmin, mavg))
            if TARGET_STRESS:
                print("[R-13] >>> 목표 TARGET_STRESS=%.0f Pa 대비 평균 면내응력 배율 = %.3f"
                      % (TARGET_STRESS, mavg / TARGET_STRESS))
            print("[R-13] >>> 판정: 면내평균<0 비율 ~0 이면 base state 인장지배(수치문제) / "
                  "0.2 이상이면 slack·주름 상태(문서상 '예압>좌굴하중' 케이스)")
        elif mid is not None:
            print("[R-13] >>> 헤드라인: 면내평균<0 비율=%.4f  minP<0 비율=%.4f" % (cm, cmin))

        # 앵커 반력 (instance 이름 포함)
        try:
            rf = frame.fieldOutputs['RF']
            best = []
            for v in rf.values:
                d = v.data
                mag = math.sqrt(sum(float(x) ** 2 for x in d))
                if mag > 1e-12:
                    inst = ''
                    try:
                        inst = v.instance.name
                    except Exception:
                        pass
                    best.append((mag, v.nodeLabel, inst,
                                 tuple(round(float(x), 8) for x in d)))
            best.sort(reverse=True)
            print("[R-13] RF 비영 노드 %d개, 상위 8개:" % len(best))
            for mag, lab, inst, comp in best[:8]:
                print("       %.6g N  node %s  inst=%s  %s" % (mag, lab, inst, comp))
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
