"""
[compare-states] 두 ODB 를 여러 step time 에서 나란히 비교한다 (읽기 전용)

용도
  같은 모델의 두 런(예: 기본 C_n^a 로 돌린 대조군 vs C_n^a=1 완화 런)이
  **공통으로 도달한 구간**에서 같은 해인지 판정한다. 한 점 비교는 우연일 수 있으므로
  여러 step time 을 한 번에 훑는다.
  PowerShell 루프/파싱이 필요 없다 - 양쪽 ODB 를 직접 읽어 표를 만든다.

사용법 (인자는 전부 위치 인자 - 아래 주의 참조)
  abaqus python compare_states.py <odbA> <odbB> <step> [times] [tol] [node]
    times : 쉼표 구분 step time 목록. 생략하면 두 런의 공통 구간을 6등분
    tol   : 판정 허용 상대차 %% (기본 2.0)
    node  : |RF| 를 읽을 절점 번호 (기본 2)
  예) abaqus python compare_states.py HF_ls.odb HF_rc.odb Step-Postbuckle 0.02,0.04,0.06,0.08,0.10,0.118

주의 (실측 2026-09-30)
  `abaqus python` 실행기는 `--key=value` 형식 인자를 **거부**한다:
    ABAQUS Error: Argument "--times=..." is not a valid argument.
                   It has both equality sign and prepended dash.
  따라서 이 스크립트는 **위치 인자만** 쓴다. 대시가 붙은 인자는 방어적으로 무시한다.

판정
  스칼라 지표(평균 면내응력 / 압축 면적비 / maxP / minP / |RF| / |u3|)에 대해
  모든 step time 에서 상대차가 --tol % 이내면 "같은 해" 로 본다.
  국소량(max|u3|)은 증분 이력 차이로 흔들리므로 tol 을 별도로 보라.

주의
  - 면내응력 계산은 base_state_probe.py 의 _aggregate 를 그대로 import 해 쓴다(같은 양 보장).
  - ODB 의 instance.nodes / elements 는 repository 가 아니라 배열이다(keys() 없음).
  - abaqusConstants.CENTROID 를 못 import 하면 전체 단면점으로 대체한다(base_state_probe 와 동일).
"""
import io
import math
import sys


def _norm_times(s):
    out = []
    for part in s.split(','):
        part = part.strip()
        if part:
            out.append(float(part))
    return out


def _nearest_frame(step, t):
    best_d, best_i, best_f = None, 0, None
    for i, f in enumerate(step.frames):
        d = abs(f.frameValue - t)
        if best_d is None or d < best_d:
            best_d, best_i, best_f = d, i, f
    return best_i + 1, best_f


def _thrust_loss(frame, instance, conn, max_label):
    """eval_abaqus.calc_thrust_loss 를 그대로 쓴다(같은 목적함수를 보장). 실패하면 None."""
    try:
        from eval_abaqus import calc_thrust_loss
    except Exception:
        return None
    try:
        return float(calc_thrust_loss(frame, instance, conn, max_label))
    except Exception:
        return None


def _metrics(frame, odb, node, instances, tk=None):
    """프레임 하나에서 비교용 지표를 뽑는다."""
    from base_state_probe import _aggregate

    # --- 면적 가중 (EVOL) ---
    vols = {}
    try:
        for v in frame.fieldOutputs['EVOL'].values:
            d = v.data
            vols[v.elementLabel] = float(d[0]) if isinstance(d, (tuple, list)) else float(d)
    except Exception:
        pass

    try:
        from abaqusConstants import CENTROID
        sf = frame.fieldOutputs['S'].getSubset(position=CENTROID)
    except Exception:
        sf = frame.fieldOutputs['S']

    groups = {}
    for v in sf.values:
        d = v.data
        if len(d) < 3:
            continue
        try:
            num = v.sectionPoint.number
        except Exception:
            num = None
        w = vols.get(v.elementLabel, 0.0)
        groups.setdefault(num, []).append((float(d[0]), float(d[1]), float(d[2]), w if w > 0 else 0.0))
    rows = groups.get(5) or groups.get(None) or (list(groups.values())[0] if groups else [])
    n, sh_mean, sh_min, mean_s, s1mx, s2mn = _aggregate(rows)

    # --- out-of-plane ---
    u3 = []
    try:
        for v in frame.fieldOutputs['U'].values:
            if len(v.data) >= 3:
                u3.append(v.data[2])
    except Exception:
        pass
    a3 = sorted(abs(x) for x in u3)
    u3_max = a3[-1] if a3 else float('nan')
    u3_med = a3[len(a3) // 2] if a3 else float('nan')

    # --- 케이블 루트 반력 크기 ---
    rf = {}
    try:
        rf_all = frame.fieldOutputs['RF']
        for inst_name in instances:
            inst = odb.rootAssembly.instances.get(inst_name) if hasattr(odb.rootAssembly.instances, 'get') \
                else (odb.rootAssembly.instances[inst_name] if inst_name in odb.rootAssembly.instances else None)
            if inst is None:
                continue
            d = None
            for v in rf_all.getSubset(region=inst).values:
                if getattr(v, 'nodeLabel', None) == node:
                    d = v.data
                    break
            if d is not None:
                rf[inst_name] = math.sqrt(sum(x * x for x in d))
    except Exception:
        pass

    hl = None
    if tk is not None:
        hl = _thrust_loss(frame, tk['inst'], tk['conn'], tk['max_label'])
    return dict(n=n, mean=mean_s, sh_mean=sh_mean, sh_min=sh_min, maxP=s1mx, minP=s2mn,
                u3max=u3_max, u3med=u3_med, rf=rf, hl=hl)


def _pct(a, b):
    if a == 0 or not (a == a) or not (b == b):
        return float('nan')
    return (b - a) / abs(a) * 100.0


def main():
    try:
        from odbAccess import openOdb
    except ImportError:
        print("[compare-states] odbAccess 를 못 읽었습니다 - 'abaqus python' 으로 실행해야 합니다.")
        return 2

    # `--key=value` 는 abaqus python 실행기가 거부하므로 위치 인자만 쓴다.
    argv = [a for a in sys.argv[1:] if not a.startswith('-')]
    if len(argv) < 3:
        print(__doc__)
        return 2
    odb_a, odb_b, step_name = argv[0], argv[1], argv[2]
    # 7번째 위치 인자 = 시각 불일치 허용 % (기본 2.0). 이보다 크면 그 지점은 판정에서 제외한다.
    gap_tol = float(argv[6]) if len(argv) > 6 else 2.0
    times = _norm_times(argv[3]) if len(argv) > 3 else None
    tol = float(argv[4]) if len(argv) > 4 else 2.0
    node = int(argv[5]) if len(argv) > 5 else 2
    instances = ["INST_CABLE_RIGHT", "INST_CABLE_LEFT", "INST_CABLE_CL", "INST_CABLE_CR", "INST_CABLE_TOP"]

    odbs, steps = [], []
    for path in (odb_a, odb_b):
        try:
            odb = openOdb(path, readOnly=True)
        except Exception as e:
            print("[compare-states] odb 열기 실패 (%s): %s" % (path, e))
            return 1
        if step_name not in odb.steps.keys():
            print("[compare-states] %s 에 step '%s' 없음. 있는 step: %s"
                  % (path, step_name, list(odb.steps.keys())))
            return 1
        odbs.append(odb)
        steps.append(odb.steps[step_name])

    # 추력손실용 재료 (Odb 당 한 번)
    tks = []
    for odb in odbs:
        try:
            inst = odb.rootAssembly.instances['MEMBRANE-1']
            conn = [e.connectivity for e in inst.elements]
            max_label = (inst.nodes[-1].label + 100) if inst.nodes else 100000
            tks.append(dict(inst=inst, conn=conn, max_label=max_label))
        except Exception:
            tks.append(None)

    if times is None:
        # 두 런의 공통 구간을 자동으로 6 등분
        end = min(steps[0].frames[-1].frameValue, steps[1].frames[-1].frameValue)
        times = [round(end * k / 6.0, 6) for k in range(1, 7)]

    lbl = [odb_a, odb_b]
    print("=" * 96)
    print("[compare-states] %s  vs  %s   step=%s   tol=%.1f%%" % (odb_a, odb_b, step_name, tol))
    print("[compare-states] 공통 종료 step time = %.6g / 비교 지점 %d개"
          % (min(steps[0].frames[-1].frameValue, steps[1].frames[-1].frameValue), len(times)))
    print("=" * 96)

    worst = {}
    for t in times:
        res = []
        for k in range(2):
            idx, fr = _nearest_frame(steps[k], t)
            if fr is None:
                print("[compare-states] 프레임이 없습니다: %s" % lbl[k])
                return 1
            res.append((idx, fr.frameValue, _metrics(fr, odbs[k], node, instances, tks[k])))
        print()
        gap = abs(res[1][1] - res[0][1]) / max(res[0][1], 1e-30) * 100.0
        # 앞선 검증: 지표 차이가 이 시각 차이에 그대로 비례한다(측정 비 1.03).
        #   시각이 어긋난 지점은 같은 상태 비교가 아니므로 verdict 에서 제외한다.
        usable = gap <= gap_tol
        print()
        print("### 목표 step time %.6g    (A frame %d @ %.6g / B frame %d @ %.6g)"
              % (t, res[0][0], res[0][1], res[1][0], res[1][1]))
        print("    시각 불일치 %+.2f%%  ->  %s"
              % (gap, "판정 대상" if usable else "** 판정 제외 (프레임 격자 차이) **"))
        print("  %-22s %14s %14s %9s" % ("지표", "A", "B", "차이"))
        for key, name, kind in [("mean", "평균 면내응력 [Pa]", "int"),
                                ("sh_min", "압축면적비 minP<0", "int"),
                                ("sh_mean", "압축면적비 mean<0", "int"),
                                ("maxP", "maxP [Pa]", "int"),
                                ("minP", "minP [Pa]", "int"),
                                ("u3max", "max|u3| [m]", "loc"),
                                ("u3med", "median|u3| [m]", "loc")]:
            va, vb = res[0][2][key], res[1][2][key]
            d = _pct(va, vb)
            if usable:
                worst.setdefault(key, []).append(abs(d))
            if d != d:
                mark = "  (기준값 0 - 비교 불가)"
            elif abs(va) < 1e-9:
                mark = "  (A 가 0 에 가까움 - 비교 불가)"
            else:
                mark = "  <= tol" if abs(d) <= tol else "  ** 초과 **"
            print("  %-22s %14.6g %14.6g %+8.2f%%%s" % (name, va, vb, d, mark))
        ha, hb = res[0][2].get('hl'), res[1][2].get('hl')
        if ha is not None and hb is not None:
            d = _pct(ha, hb)
            if usable:
                worst.setdefault("thrust_loss", []).append(abs(d))
            mark = "" if abs(d) <= tol else "  ** 초과 **"
            print("  %-22s %14.6g %14.6g %+8.2f%%%s" % ("추력손실(목적함수)", ha, hb, d, mark))
        for inst_name in instances:
            va = res[0][2]['rf'].get(inst_name)
            vb = res[1][2]['rf'].get(inst_name)
            if va is None or vb is None:
                continue
            d = _pct(va, vb)
            key = "rf:" + inst_name
            if usable:
                worst.setdefault(key, []).append(abs(d))
            if d != d:
                mark = "  (기준값 0 - 비교 불가)"
            elif abs(va) < 1e-9:
                mark = "  (A 가 0 에 가까움 - 비교 불가)"
            else:
                mark = "  <= tol" if abs(d) <= tol else "  ** 초과 **"
            print("  %-22s %14.6g %14.6g %+8.2f%%%s" % ("|RF| " + inst_name.replace('INST_CABLE_', ''),
                                                       va, vb, d, mark))

    print()
    print("=" * 96)
    print("[compare-states] 요약: 각 지표의 지점별 상대차 최대값")
    ok_all = True
    skipped = []
    for key in sorted(worst):
        mx = max(worst[key])
        if mx != mx:                      # nan = 기준값이 0 이라 상대차가 정의되지 않음
            skipped.append(key)
            continue
        flag = "OK" if mx <= tol else "초과"
        if mx > tol and key not in ("u3max", "u3med"):
            ok_all = False
        print("   %-30s 최대 %7.2f%%   [%s]" % (key, mx, flag))
    if skipped:
        print("   비교 제외(기준값이 0 또는 0 근처라 상대차 무의미): %s" % ", ".join(sorted(skipped)))
    print()
    if ok_all:
        print("[compare-states] >>> 적분량 지표가 모두 tol(%.1f%%) 이내 -> 같은 해로 판정" % tol)
    else:
        print("[compare-states] >>> tol 을 넘는 지표가 있습니다 -> 같은 해로 볼 수 없음")
    print("[compare-states] 참고: max|u3|/median|u3| 는 국소량이라 증분 이력 차이로 흔들린다(tol 별도 적용).")
    print("=" * 96)
    return 0


if __name__ == '__main__':
    sys.exit(main())
