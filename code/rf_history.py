"""
[RF-history] 프레임별 케이블 반력 이력 -> 리밋포인트/스냅백 신호 판정 (읽기 전용)

용도
  "이 문제가 리밋포인트를 지나는가?" 를 기존 ODB 만으로 판정한다.
  지정 변위 제어 정적 해석에서 스텝 타임은 강제로 단조 증가하므로, 하중 자체의
  단조성은 .sta 로는 알 수 없다. 대신 **반력**을 보면 된다:
    - 반력이 단조 증가  -> 스냅백 신호 없음 -> Riks 의 경로 추적 이점 없음
    - 반력이 감소하는 구간 존재 -> 구조가 하중을 내려놓는 구간 = 리밋포인트/스냅백
      -> 하중·변위 제어 뉴턴은 원리적으로 그 구간을 넘지 못한다 (Riks 필수)

사용법
  abaqus python rf_history.py <odb> <step> [--nodes 2] [--dof 1] [--instances A,B,C]
  예: abaqus python rf_history.py HF_rc.odb Step-Postbuckle

주의
  - ODB 의 instance.nodes / elements 는 repository 가 아니라 배열이다.
    (keys() 가 없다 - 2026-09-29 확립) 이 스크립트는 keys() 를 쓰지 않는다.
  - RF 는 절점 필드이므로 instance 별 getSubset(region=...) 로 읽는다.
  - 노드 번호 기본값 2 는 케이블 인스턴스의 루트 절점(base_state_probe 의 RF 출력과 동일).
"""
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
                val = sum(d) if comp is None else d[comp]
                row[inst_name] = val
                total += val
            if missing and i == 0:
                print("[RF-history] 주의: 다음 인스턴스/절점에서 RF 를 못 읽었습니다 -> %s" % ", ".join(missing))
            hist.append((fr.frameValue, row, total))

        if not hist:
            print("[RF-history] RF 를 읽은 프레임이 없습니다.")
            return 1

        hdr = "  %-4s %-12s %12s" % ("fr", "stepTime", "sum")
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
        print("[RF-history] 프레임 %d개, step time %.6g 까지" % (len(hist), t_end))
        if not drops:
            print("[RF-history] >>> 판정: 총 반력이 단조 증가 (감소 0회)")
            print("[RF-history]     => 스냅백/리밋포인트 신호 없음.")
            print("[RF-history]     => Riks 의 '경로 추적' 이점이 이 하중 경로에는 없습니다.")
            print("[RF-history]        (남는 이점은 컷백 오버헤드와 lambda_end 종료뿐)")
        else:
            worst = min(drops, key=lambda d: (d[3] - d[2]) / d[2] if d[2] else 0.0)
            rel = (worst[3] - worst[2]) / worst[2] * 100 if worst[2] else 0.0
            print("[RF-history] >>> 판정: 총 반력 감소 %d회 (전체 %d 프레임)" % (len(drops), len(hist)))
            print("[RF-history]     최대 감소 %.2f%%  (frame %d, step time %.6g, %.6g -> %.6g)"
                  % (rel, worst[0], worst[1], worst[2], worst[3]))
            idx = [d[0] for d in drops][:12]
            print("[RF-history]     감소 프레임: %s%s" % (idx, " ..." if len(drops) > 12 else ""))
            print("[RF-history]     => 리밋포인트/스냅백 신호. 하중·변위 제어 뉴턴은 이 구간을")
            print("[RF-history]        원리적으로 넘지 못하므로 Riks(아크길이 제어)가 정당화됩니다.")
            print("[RF-history]        단, 감소폭이 수치잡음 수준인지 먼저 확인하세요 (아래 참고).")
            print("[RF-history]        참고: 감소폭이 0.1%% 미만이고 반복적으로 나타나면 수치잡음,")
            print("[RF-history]              한 번의 큰 감소(수 %%)면 물리적 스냅백입니다.")
        print("=" * 78)
        return 0
    finally:
        try:
            odb.close()
        except Exception:
            pass


if __name__ == '__main__':
    sys.exit(main())
