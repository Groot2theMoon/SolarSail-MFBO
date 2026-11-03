# -*- coding: utf-8 -*-
"""특정 노드/자유도의 국소 진단 (읽기 전용).

왜: HF 런이 increment 926 근처에서 한 노드의 면외 자유도(u3)가 부호를 바꾸며 커지는
    진동을 보이고 `DISP. CORRECTION TOO LARGE` -> `APPEARS TO BE DIVERGING` 으로
    1/4 컷백을 반복한다 (관측 2026-09-29: node 970 DOF 3). 그 노드가
      (a) 경계/부착부(클램프 패치, 케이블, 모서리 구속)에 붙어 있는지,
      (b) 내부 주름 봉우리인지
    에 따라 처방이 완전히 달라진다: (a)면 하중경로/결합 문제, (b)면 국소 모드(요소/감쇠).

사용:
    abaqus python node_probe.py HF_Postbuckle.odb Step-Postbuckle 970
    abaqus python node_probe.py HF_Postbuckle.odb Step-Postbuckle 970 3

출력:
  1) 좌표
  2) 그 노드가 속한 모든 노드셋 (어셈블리 + 인스턴스)  -> 부착부 여부 판정
  3) 인접 요소와 그 요소들의 마지막 프레임 u3          -> 이웃 대비 국소 이상 판정
  4) 프레임별 u1,u2,u3 이력                            -> 진동/발산 확인

해석에는 영향이 없다(읽기 전용).
"""
from __future__ import print_function
import sys


def _find_node(inst, label):
    """인스턴스의 nodes 를 라벨로 찾는다.

    실측(2026-09-29): 실제 ODB 에서 `inst.nodes` 는 repository 가 아니라
    OdbMeshNodeArray(keys() 없는 시퀀스)다 ->
      AttributeError: 'OdbMeshNodeArray' object has no attribute 'keys'
    그래서 (a) keys() 가 있으면 repository 로, (b) 없으면 시퀀스로 순회한다.
    알 수 없는 원소 타입은 조용히 건너뛴다(도구가 트레이스백으로 죽지 않게).
    """
    nodes = getattr(inst, 'nodes', None)
    if nodes is None:
        return None
    if hasattr(nodes, 'keys'):
        try:
            return nodes[label]
        except Exception:
            return None
    try:
        for nd in nodes:
            if int(getattr(nd, 'label', -1)) == label:
                return nd
    except TypeError:
        return None
    return None


def _labels_of(nodeset):
    """어셈블리 셋은 (instanceName, label) 튜플, 인스턴스 셋은 노드 객체를 준다."""
    out = []
    for n in nodeset.nodes:
        if isinstance(n, (tuple, list)):
            out.append(int(n[1]))
        else:
            lab = getattr(n, 'label', None)
            if lab is not None:
                out.append(int(lab))
    return out


def main():
    try:
        from odbAccess import openOdb
    except ImportError:
        print("[node-probe] odbAccess 를 못 읽었습니다 - 'abaqus python' 으로 실행하세요.")
        return 2
    if len(sys.argv) < 4:
        print("사용: abaqus python node_probe.py <odb> <step> <nodeLabel> [dof=3]")
        return 2
    odb_path = sys.argv[1]
    step_name = sys.argv[2]
    label = int(sys.argv[3])
    dof = int(sys.argv[4]) if len(sys.argv) > 4 else 3

    odb = openOdb(odb_path, readOnly=True)
    try:
        if step_name not in odb.steps.keys():
            print("[node-probe] 스텝이 없습니다: %s (있는 스텝: %s)" % (step_name, list(odb.steps.keys())))
            return 1
        step = odb.steps[step_name]
        frames = len(step.frames)
        print("=" * 70)
        print("[node-probe] %s / %s / 노드 %d (dof %d) / 프레임 %d개"
              % (odb_path, step_name, label, dof, frames))
        print("=" * 70)
        if frames == 0:
            print("  프레임이 없습니다 (수렴 프레임 미기록).")
            return 0

        ra = odb.rootAssembly
        # 1) 소속 인스턴스 + 좌표
        target = None
        for iname in ra.instances.keys():
            inst = ra.instances[iname]
            nd = _find_node(inst, label)
            if nd is not None:
                target = inst
                print("[1] 인스턴스: %s" % iname)
                print("    좌표 (x,y,z) = %s" % (tuple(nd.coordinates),))
                break
        if target is None:
            print("[1] 노드 %d 를 어셈블리에서 찾지 못했습니다." % label)
            return 1

        # 2) 소속 노드셋
        print("[2] 소속 노드셋:")
        found = 0
        for holder, tag in ((ra, "assembly"), (target, "instance")):
            try:
                names = holder.nodeSets.keys()
            except Exception:
                continue
            for nm in names:
                try:
                    labs = _labels_of(holder.nodeSets[nm])
                except Exception:
                    continue
                if label in labs:
                    print("    %-8s %s  (구성 노드 %d개)" % (tag, nm, len(labs)))
                    found += 1
        if found == 0:
            print("    없음 -> 어느 구속/부착 셋에도 속하지 않는 자유 노드")

        # 3) 인접 요소 + 이웃의 마지막 프레임 u3
        last = step.frames[-1]
        u3 = {}
        try:
            for v in last.fieldOutputs['U'].values:
                try:
                    u3[int(v.nodeLabel)] = float(v.data[2])
                except Exception:
                    pass
        except Exception as e:
            print("[3] U 읽기 실패: %s" % e)

        adj = []
        if target.elements is not None:
            for el in target.elements:
                try:
                    conn = [int(c) for c in el.connectivity]
                except Exception:
                    continue
                if label in conn:
                    adj.append((el.label, conn, getattr(el, 'type', '?')))
        print("[3] 인접 요소 %d개 (각 요소 노드들의 마지막 프레임 u3):" % len(adj))
        for elabel, conn, etype in adj:
            vals = []
            for c in conn:
                if c in u3:
                    vals.append("%d:%+.3e" % (c, u3[c]))
            print("    el %-8d (type %s)" % (elabel, etype))
            print("        " + "  ".join(vals))
        if label in u3:
            print("    -> 이 노드 u3 = %+.6e m" % u3[label])

        # 4) 이력
        print("[4] 프레임별 이력 (앞 5 + 뒤 8):")
        idxs = list(range(min(5, frames))) + list(range(max(0, frames - 8), frames))
        seen = set()
        for i in idxs:
            if i in seen:
                continue
            seen.add(i)
            f = step.frames[i]
            val = None
            try:
                for v in f.fieldOutputs['U'].values:
                    if int(v.nodeLabel) == label:
                        val = v.data
                        break
            except Exception:
                pass
            if val is None:
                continue
            print("    frame %-4d stepTime %-12.6g u1=%+.6e u2=%+.6e u3=%+.6e"
                  % (i, f.frameValue, val[0], val[1], val[2]))
        print("=" * 70)
        return 0
    finally:
        try:
            odb.close()
        except Exception:
            pass


if __name__ == '__main__':
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception:
        import traceback
        print("[node-probe] 실패 - 아래 추적을 그대로 보내주세요:")
        traceback.print_exc()
        sys.exit(3)
