#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""미완주(stall/컷백) HF 잡의 **현재 상태**를 읽는다 (2026-10-12).

    abaqus python probe_partial_hf.py <job_or_odb> [GLOBAL_FINAL] [CLAMP_FINAL]

무엇을 답하나
  ① 어디까지 갔나: 스텝별 프레임 수 + 마지막 프레임의 스텝 시간(=진행률 tau)
  ② 그때 **규정된(prescribed) 케이블/클램프 변위** = tau x 상수 (기본 1e-3 m = GLOBAL_FINAL)
     — StaticStep 기본 RAMP 이므로 스텝 시간에 선형이다.
  ③ 실제 변형: 마지막 프레임의 U 요약(max|U1|,|U2|,|U3| = 주름 진폭 u_z,max) + 노드 위치
  ④ 잡음/반력: RF 가 있으면 최대/합, 그리고 케이블·RP 노드셋이 보이면 그 값도 찍는다.

주의: 미완주 상태를 MFBO 목적함수로 쓰면 안 된다(수렴 안 한 값). 진단·경향 판단용이다.
"""
from __future__ import print_function
import os
import sys

from odbAccess import openOdb


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    p = argv[1]
    if not p.lower().endswith('.odb'):
        p = p + '.odb'
    if not os.path.exists(p):
        print("!!! ODB 없음: %s" % p)
        return 1
    gfin = float(argv[2]) if len(argv) > 2 else 1.0e-3
    cfin = float(argv[3]) if len(argv) > 3 else gfin     # 클램프는 d_c 배 — 모르면 같은 값으로 둔다

    odb = openOdb(path=p, readOnly=True)
    print("=" * 72)
    print("[partial] %s" % os.path.basename(p))
    print("  규정변위 스케일 가정: GLOBAL_FINAL=%.6g m / CLAMP_FINAL=%.6g m" % (gfin, cfin))

    order = []
    for name in odb.steps.keys():
        st = odb.steps[name]
        fr = st.frames
        t_last = fr[-1].frameValue if fr else -1.0
        order.append(name)
        print("  STEP %-18s frames=%-5d last_step_time=%.6e  -> 코너 %.4e m / 클램프 %.4e m"
              % (name, len(fr), t_last, t_last * gfin, t_last * cfin))

    # 마지막 스텝의 마지막 프레임에서 필드 요약
    last_step = odb.steps[order[-1]]
    inst = None
    for iname in odb.rootAssembly.instances.keys():
        if 'MEMBRANE' in iname.upper():
            inst = odb.rootAssembly.instances[iname]
    print("  instance: %s" % (inst.name if inst else '(없음)'))
    if last_step.frames and inst is not None and 'U' in last_step.frames[-1].fieldOutputs.keys():
        fr = last_step.frames[-1]
        U = fr.fieldOutputs['U'].getSubset(region=inst).values
        mx = [0.0, 0.0, 0.0]
        arg = [0, 0, 0]
        for v in U:
            for k in range(3):
                if abs(v.data[k]) > abs(mx[k]):
                    mx[k], arg[k] = v.data[k], v.nodeLabel
        print("  U: max|U1|=%.4e (node %d) / max|U2|=%.4e (node %d) / max|U3|=%.4e (node %d)"
              % (abs(mx[0]), arg[0], abs(mx[1]), arg[1], abs(mx[2]), arg[2]))
        print("     (max|U3| = 주름 진폭 u_z,max — 두께 5e-06 m 대비 %.1f t)"
              % (abs(mx[2]) / 5.0e-6))
    if last_step.frames and 'RF' in last_step.frames[-1].fieldOutputs.keys():
        rf = last_step.frames[-1].fieldOutputs['RF'].values
        s = [0.0, 0.0, 0.0]
        mx = [0.0, 0.0, 0.0]
        for v in rf:
            for k in range(3):
                s[k] += v.data[k]
                if abs(v.data[k]) > abs(mx[k]):
                    mx[k] = v.data[k]
        print("  RF: 합(%.4e, %.4e, %.4e) / max(%.4e, %.4e, %.4e)" % tuple(s + mx))

    # 케이블/RP 노드셋이 보이면 그 값
    cand = []
    for src, nm in ((odb.rootAssembly, 'assembly'),):
        for k in src.nodeSets.keys():
            if any(t in k.upper() for t in ('CABLE', 'RP_', 'CLAMP', 'END')):
                cand.append((nm, k))
    for k in list(getattr(inst, 'nodeSets', {}).keys()) if inst is not None else []:
        if any(t in k.upper() for t in ('CABLE', 'RP_', 'CLAMP', 'END')):
            cand.append((inst.name, k))
    print("  후보 노드셋 %d개" % len(cand))
    for owner, k in cand:
        try:
            ns = odb.rootAssembly.nodeSets[k] if owner == 'assembly' else inst.nodeSets[k]
            vals = last_step.frames[-1].fieldOutputs['U'].getSubset(region=ns).values
            if vals:
                v0 = vals[0]
                print("    %-28s node %-6d U=(%.4e, %.4e, %.4e)"
                      % (k, v0.nodeLabel, v0.data[0], v0.data[1], v0.data[2]))
        except Exception as e:
            print("    %-28s (읽기 실패: %s)" % (k, e))
    odb.close()
    print("=" * 72)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
