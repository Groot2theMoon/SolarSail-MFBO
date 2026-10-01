# -*- coding: utf-8 -*-
"""base state 응력 측정 v2 (R-13 / 스파이크 10-0).

좌굴 base state(= 직전 일반 스텝 말단)의 면내 응력 분포를 .odb 에서 읽어

  1) 단면점(section point) 별 통계 — 압축이 '진짜 면내 압축'인지 굽힘(표면) artifact 인지 분리
  2) 면적 가중(요소 EVOL) 통계        — 요소 개수 비율이 아니라 실제 면적 비율
  3) 면내 평균응력 (S11+S22)/2 < 0 비율 — 단면점에 무관한 압축 지표
  4) 앵커 반력 (RF, instance 이름 포함 = 케이블 장력)
  5) 면외변위 |u3| 통계      — 주름이 실제로 발생했는지 (B안 검증용)

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
import os

TARGET_STRESS = 7000.0   # Pa, R-13 목표 운용점(미검증)

# 중앙 영역 반경 [m]. 논문의 7000 Pa 는 '사분면 중앙(the centre of each solar sail quadrant)'의
# 막 응력이다. 우리 모델은 사분면 하나(삼각형)이므로 그 중앙은 무게중심 = ((10+20+0)/3, (10+0+0)/3)
# = (10.0, 3.333) 이다. 반경 안의 요소만 모아 면적가중 통계를 낸다.
#   - MIDSECTION 은 '두께 방향 단면점' 통계이지 공간적 중앙이 아니다(혼동 주의).
#   - 좌표를 못 읽어도 기존 출력은 그대로 나온다(조용히 건너뛴다).
R_CENTRE = 1.0
CENTRE_XY = (10.0, 10.0 / 3.0)


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


BUILD = '2026-10-01a (centre+radius)'
#   개정 이력: (a) 중앙 영역(CENTRE) 통계 추가, (b) 4번째 인자로 중앙 반경 덮어쓰기.
#   사용자 실행형 도구이므로 빌드를 첫 줄에 찍는다. 출력에 이 줄이 없거나 반경이 반영되지
#   않으면 옛 리비전이다(실제로 2026-10-01 에 이 함정에 한 번 걸렸다).


def main():
    print('[R-13] build %s' % BUILD)
    try:
        from odbAccess import openOdb
    except ImportError:
        print("[R-13] odbAccess 를 못 읽었습니다 - 'abaqus python' 으로 실행해야 합니다.")
        return 2

    odb_path = sys.argv[1] if len(sys.argv) > 1 else 'Buckle_Analysis.odb'
    step_name = sys.argv[2] if len(sys.argv) > 2 else 'Step-GlobalTension'
    # [2026-09-29] 선택 인자 3: 목표 step time. 두 런을 '같은 하중 수준'에서 대조할 때 쓴다
    #   (예: C_n^a 완화 런과 기본 기준 런을 step time 0.05 에서 비교). 생략하면 마지막 프레임.
    # [2026-10-01] 선택 인자 4: 중앙 영역 반경 [m]. 논문의 'centre of each quadrant' 정의에
    #   따라 중앙 응력이 달라지므로(반경 1 m 에서 3958 Pa, 전막 평균 6634 Pa) 민감도를 본다.
    global R_CENTRE
    if len(sys.argv) > 4:
        try:
            R_CENTRE = float(sys.argv[4])
            print("[R-13] 중앙 반경을 인자로 덮어씀: R_CENTRE = %.4g m (기본 1.0)" % R_CENTRE)
        except ValueError:
            print("[R-13] 4번째 인자는 숫자(중앙 반경, m)여야 합니다: %r" % sys.argv[4])
            return 2

    t_target = None
    if len(sys.argv) > 3:
        try:
            t_target = float(sys.argv[3])
        except ValueError:
            print("[R-13] 3번째 인자는 숫자(step time)여야 합니다: %r" % sys.argv[3])
            return 2

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
        sel_idx = nf
        if t_target is not None:
            # frames 는 repository 일 수도 배열일 수도 있으므로 .index() 를 쓰지 않고 직접 센다.
            best_d, best_i = None, 0
            for i, f in enumerate(step.frames):
                d = abs(f.frameValue - t_target)
                if best_d is None or d < best_d:
                    best_d, best_i, frame = d, i, f
            sel_idx = best_i + 1
            print("[R-13] 목표 step time %.6g -> 가장 가까운 프레임 (frame %d/%d, step time %s)"
                  % (t_target, sel_idx, nf, frame.frameValue))
        print("[R-13] base state = %s / %s  (frame %d/%d, step time %s)"
              % (odb_path, step_name, sel_idx, nf, frame.frameValue))

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

        # 요소 라벨 -> 중심 좌표. 실패해도 치명적이지 않다(중앙 통계만 생략).
        cent = {}
        try:
            for inst in odb.rootAssembly.instances.values():
                nl = {}
                for nd in inst.nodes:
                    nl[nd.label] = nd.coordinates
                for el in inst.elements:
                    cs = [nl[c] for c in el.connectivity if c in nl]
                    if cs:
                        cent[el.label] = (sum(c[0] for c in cs) / len(cs),
                                          sum(c[1] for c in cs) / len(cs))
            print("[R-13] 요소 중심 좌표 %d개 확보 (중앙 반경 %.3g m)" % (len(cent), R_CENTRE))
        except Exception as e:
            print("[R-13] 요소 좌표를 못 읽음 -> 중앙 영역 통계 생략 (%s)" % e)

        groups = {}      # section number -> rows
        centrerows = []  # 중앙 반경 안의 rows
        for v in sf.values:
            d = v.data
            if len(d) < 3:
                continue
            try:
                num = v.sectionPoint.number
            except Exception:
                num = None
            w = vols.get(v.elementLabel, 0.0)
            row = (float(d[0]), float(d[1]), float(d[2]), w if w > 0 else 0.0)
            groups.setdefault(num, []).append(row)
            xy = cent.get(v.elementLabel)
            if xy is not None:
                dx = xy[0] - CENTRE_XY[0]
                dy = xy[1] - CENTRE_XY[1]
                if dx * dx + dy * dy <= R_CENTRE * R_CENTRE:
                    centrerows.append((row, num))

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
                print("[R-13] >>> (참고용) 전막 평균 대비 TARGET_STRESS 배율 = %.3f  -- 논문 대조 대상 아님(논문은 전막 평균을 언급하지 않음). 대조는 아래 CENTRE"
                      % (mavg / TARGET_STRESS))
            print("[R-13] >>> 판정: 면내평균<0 비율 ~0 이면 base state 인장지배(수치문제) / "
                  "0.2 이상이면 slack·주름 상태(문서상 '예압>좌굴하중' 케이스)")
        # 논문과 직접 비교되는 지표: 사분면 중앙 영역(반경 R_CENTRE)의 응력
        if centrerows:
            for lab, sel in (('CENTRE(mid)', [r for r, nm in centrerows if nm == mid]),
                             ('CENTRE(all)', [r for r, nm in centrerows])):
                if not sel:
                    continue
                n, cm, cmin, mavg, s1mx, s2mn = _aggregate(sel)
                print(_fmt('%s r=%.3g' % (lab, R_CENTRE), n, cm, cmin, mavg, s1mx, s2mn))
            sel_all = [r for r, nm in centrerows]
            mavg_c = _aggregate(sel_all)[3]
            print("[R-13] >>> CENTRE 헤드라인: 중심 (%.3g, %.3g) 반경 %.3g m / n=%d / "
                  "평균 면내응력=%+.4g Pa" % (CENTRE_XY[0], CENTRE_XY[1], R_CENTRE,
                                              len(sel_all), mavg_c))
            print("[R-13] >>> 논문 대조용: 7000 Pa 대비 중앙 평균 면내응력 배율 = %.3f"
                  % (mavg_c / TARGET_STRESS))
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

        # 면외변위 |u3| 통계 — B안(trigger + 자연 주름)에서 '주름이 실제로 발생했는지'를
        # 같은 로그로 확인하기 위한 것. A안/프리텐션 단계에서는 u3=0 고정이라 ~0 이 나온다.
        try:
            T_MEMB = 5.0e-6                        # 막 두께 [m] (코드 상수)
            uf = frame.fieldOutputs['U']
            u3 = []
            inst_u3 = {}
            for v in uf.values:
                d = v.data
                if len(d) < 3:
                    continue
                z = abs(float(d[2]))
                u3.append(z)
                try:
                    inst = v.instance.name
                except Exception:
                    inst = ''
                inst_u3.setdefault(inst, []).append(z)
            if u3:
                u3s = sorted(u3)
                n3 = len(u3s)
                f2 = sum(1 for x in u3s if x > 2.0 * T_MEMB) / float(n3)
                f20 = sum(1 for x in u3s if x > 20.0 * T_MEMB) / float(n3)
                print("[R-13] |u3| 통계 (n=%d): max=%.3e m  p99=%.3e  median=%.3e  "
                      "|u3|>2t 비율=%.4f  |u3|>20t 비율=%.4f"
                      % (n3, u3s[-1], u3s[int(0.99 * (n3 - 1))], u3s[n3 // 2], f2, f20))
                _u3ave = sum(u3s) / float(n3)
                print("[R-13] >>> 논문 대조용 진폭: u_z,max=%.4g m (=%.1f t)  u_z,ave(|u3| mean)=%.4g m (=%.1f t)"
                % (u3s[-1], u3s[-1] / T_MEMB, _u3ave, _u3ave / T_MEMB))
                print("[R-13] >>> 참조: Galhofo2022 Table A.1 STRI65+S8R5 u_z,max=2.284e-04 m (45.7 t), 2x12 주름")
                print("[R-13] >>> (논문 Tables 3/4 는 진폭을 u_z,ave 로 보고한다) 1차 요소는 주름 수 비교에 부적합")
                for inst in sorted(inst_u3.keys()):
                    arr = sorted(inst_u3[inst])
                    print("       inst=%-18s n=%-6d max|u3|=%.3e m" % (inst, len(arr), arr[-1]))
                print("[R-13] >>> 판정(주름): max|u3| 가 막 두께(%.1e m)의 수백배 이상이면 "
                      "주름 발달 / ~0 이면 평탄(주름 미발생)" % T_MEMB)
        except Exception as e:
            print("[R-13] U(u3) 읽기 실패: %s" % e)
        return 0
    finally:
        try:
            odb.close()
        except Exception:
            pass


if __name__ == '__main__':
    sys.exit(main())
