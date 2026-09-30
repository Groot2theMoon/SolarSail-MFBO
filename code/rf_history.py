"""
[RF-history] 프레임별 케이블 반력 이력 (읽기 전용)

**무엇을 재는가**: 각 케이블 루트 절점의 반력 **크기 |RF| = sqrt(Rx^2+Ry^2+Rz^2)** 의 이력.
  절대 성분을 더하면 안 된다 - 인스턴스마다 부호가 반대라 서로 상쇄되어
  물리적 의미가 없는 스칼라가 나온다(2026-09-30 실측 오류).

**무엇을 알 수 있는가**: 장력 이력.
  - 단조 증가 -> 슬랙/주름 이벤트 없이 팽팽해지는 경로
  - 감소 구간 -> 국소 슬랙 또는 주름 이벤트 (진단용)

**무엇을 알 수 없는가 (중요)**: 이 도구로 Riks 필요성을 판정할 수 없다.
  변위 제어 해석에서 반력이 줄어드는 구간(음의 강성)은 그대로 따라갈 수 있다.
  Riks(아크길이 제어)가 정말 필요한 경우는 **델타가 되돌아가는 스냅백/폴드** 뿐이고,
  변위 제어 정적 런이 완주했다면 그 경로에 스냅백은 없었다는 증거다
  (있었다면 변위 제어로 완주하지 못한다). 2026-09-30 정정.

용도

사용법
  abaqus python rf_history.py <odb> <step> [--nodes 2] [--dof 1] [--instances A,B,C]
  예: abaqus python rf_history.py HF_rc.odb Step-Postbuckle

주의
  - ODB 의 instance.nodes / elements 는 repository 가 아니라 배열이다.
    (keys() 가 없다 - 2026-09-29 확립) 이 스크립트는 keys() 를 쓰지 않는다.
  - RF 는 절점 필드이므로 instance 별 getSubset(region=...) 로 읽는다.
  - 노드 번호 기본값 2 는 케이블 인스턴스의 루트 절점(base_state_probe 의 RF 출력과 동일).
"""
import math
import sys

DEFAULT_INSTANCES = ["INST_CABLE_RIGHT", "INST_CABLE_LEFT", "INST_CABLE_CL",
                     "INST_CABLE_CR", "INST_CABLE_TOP"]
DEFAULT_NODE = 2


def _vals(field, node_label):
    """필드에서 특정 절점의 벡터를 찾아 돌려준다. 없으면 None."""
    for v in field.values:
        if getattr(v, 'nodeLabel', None) == node_label:
            return v.data
    return None


def main():
    try:
        from odbAccess import openOdb
    except ImportError:
        print("[RF-history] odbAccess 를 못 읽었습니다 - 'abaqus python' 으로 실행해야 합니다.")
        return 2

    argv = sys.argv[1:]
    if len(argv) < 2:
        print(__doc__)
        return 2
    odb_path, step_name = argv[0], argv[1]
    node = DEFAULT_NODE
    comp = None
    instances = list(DEFAULT_INSTANCES)
    for a in argv[2:]:
        if a.startswith('--nodes='):
            node = int(a.split('=', 1)[1])
        elif a.startswith('--dof='):
            comp = int(a.split('=', 1)[1]) - 1
        elif a.startswith('--instances='):
            instances = a.split('=', 1)[1].split(',')

    try:
        odb = openOdb(odb_path, readOnly=True)
    except Exception as e:
        print("[RF-history] odb 열기 실패 (%s): %s" % (odb_path, e))
        return 1

    try:
        if step_name not in odb.steps.keys():
            print("[RF-history] step '%s' 없음. 있는 step: %s" % (step_name, list(odb.steps.keys())))
            return 1
        step = odb.steps[step_name]
        nf = len(step.frames)
        if nf == 0:
            print("[RF-history] 프레임 0개")
            return 1

        print("=" * 78)
        print("[RF-history] %s / %s / node %d%s / 프레임 %d개"
              % (odb_path, step_name, node,
                 ("" if comp is None else " dof %d" % (comp + 1)), nf))
        print("=" * 78)

        hist = []                      # (frameValue, {inst: value}, total)
        for i, fr in enumerate(step.frames):
            try:
                rf = fr.fieldOutputs['RF']
            except Exception:
                continue
            row, total, missing = {}, 0.0, []
            for inst_name in instances:
                inst = odb.rootAssembly.instances[inst_name] if inst_name in odb.rootAssembly.instances else None
                if inst is None:
                    missing.append(inst_name)
                    continue
                sub = rf.getSubset(region=inst)
                d = _vals(sub, node)
                if d is None:
                    missing.append(inst_name)
                    continue
                val = math.sqrt(sum(x * x for x in d)) if comp is None else d[comp]
                row[inst_name] = val
                total += val
            if missing and i == 0:
                print("[RF-history] 주의: 다음 인스턴스/절점에서 RF 를 못 읽었습니다 -> %s" % ", ".join(missing))
            hist.append((fr.frameValue, row, total))

        if not hist:
            print("[RF-history] RF 를 읽은 프레임이 없습니다.")
            return 1

        hdr = "  %-4s %-12s %12s" % ("fr", "stepTime", "|RF|sum")
        for n in instances:
            hdr += " %12s" % n.replace("INST_CABLE_", "")
        print(hdr)
        show = list(range(min(6, len(hist)))) + (["..."] if len(hist) > 12 else []) + \
               list(range(max(6, len(hist) - 6), len(hist)))
        for k in show:
            if k == "...":
                print("   ...")
                continue
            t, row, total = hist[k]
            line = "  %-4d %-12.6g %12.6g" % (k, t, total)
            for n in instances:
                line += " %12.6g" % row.get(n, float('nan'))
            print(line)

        # --- 단조성 판정 (총 반력 기준) ---
        print()
        drops = []
        for k in range(1, len(hist)):
            if hist[k][2] < hist[k - 1][2]:
                drops.append((k, hist[k][0], hist[k - 1][2], hist[k][2]))
        t_end = hist[-1][0]
        print("[RF-history] 프레임 %d개, step time %.6g 까지 (양: 각 케이블 루트 |RF|)" % (len(hist), t_end))
        if not drops:
            print("[RF-history] >>> 장력 이력: 단조 증가 (감소 0회)")
            print("[RF-history]     => 슬랙/주름 이벤트 없이 팽팽해지는 경로.")
        else:
            worst = min(drops, key=lambda d: (d[3] - d[2]) / d[2] if d[2] else 0.0)
            rel = (worst[3] - worst[2]) / worst[2] * 100 if worst[2] else 0.0
            print("[RF-history] >>> 장력 이력: 감소 %d회 (전체 %d 프레임)" % (len(drops), len(hist)))
            print("[RF-history]     최대 감소 %.2f%%  (frame %d, step time %.6g, %.6g -> %.6g)"
                  % (rel, worst[0], worst[1], worst[2], worst[3]))
            print("[RF-history]     감소 프레임: %s%s" % ([d[0] for d in drops][:12],
                                                       " ..." if len(drops) > 12 else ""))
            print("[RF-history]     해석: 국소 슬랙 또는 주름 이벤트. 감소폭이 반복적으로 0.1%% 미만이면")
            print("[RF-history]           수치잡음, 한 번의 큰 감소(수 %%)면 물리적 슬랙 이벤트입니다.")
        # --- 인스턴스별 판정 (총합이 단조여도 개별 케이블은 이완할 수 있다) ---
        print()
        print("[RF-history] 인스턴스별 (각 케이블 장력의 거동):")
        for inst_name in instances:
            seq = [row.get(inst_name) for _, row, _ in hist if row.get(inst_name) is not None]
            if len(seq) < 2:
                continue
            label = inst_name.replace("INST_CABLE_", "")
            ups = sum(1 for k in range(1, len(seq)) if seq[k] > seq[k - 1])
            dns = sum(1 for k in range(1, len(seq)) if seq[k] < seq[k - 1])
            chg = (seq[-1] - seq[0]) / seq[0] * 100 if seq[0] else 0.0
            if dns == 0:
                tag = "단조 증가"
            elif ups == 0:
                tag = "단조 감소 = 지속적 이완"
            else:
                tag = "비단조 (증가 %d / 감소 %d)" % (ups, dns)
            print("   %-6s %-28s %8.4g -> %8.4g N  (%+7.2f%%)"
                  % (label, tag, seq[0], seq[-1], chg))

        print()
        print("[RF-history] 주의: 이 결과로 Riks 필요성을 판정할 수 없습니다. 변위 제어 해석은")
        print("[RF-history]       반력이 줄어드는(음의 강성) 구간을 그대로 따라갑니다. Riks 가 필요한")
        print("[RF-history]       것은 delta 가 되돌아가는 스냅백/폴드뿐이고, 변위 제어 런이 완주했다면")
        print("[RF-history]       그 경로에 스냅백은 없었다는 증거입니다.")
        print("=" * 78)
        return 0
    finally:
        try:
            odb.close()
        except Exception:
            pass


if __name__ == '__main__':
    sys.exit(main())
