# -*- coding: utf-8 -*-
"""mesh_quality_probe.py -- 왜곡 요소가 어디에 몰려 있는지 ODB 에서 직접 확인한다.

배경: .dat 에 "***WARNING: 1711 elements are distorted ... element set WarnElemDistorted" 가
  뜨는데, 1711 은 x_c=0.5 / 0.25 두 덱에서 **완전히 동일**했다. 즉 임퍼펙션·하중과 무관한
  메쉬 생성의 결정론적 결과다. 원인 후보는 (a) apex 폭 0.28 m 에 요소 0.10 m 가 2.8개만
  들어가는 좁은 영역, (b) technique=FREE + algorithm=MEDIAL_AXIS, (c) deviationFactor=0.1.

사용법: abaqus python mesh_quality_probe.py <odb>
출력  : 왜곡 요소 개수, 요소 중심 (x,y) 분포(높이 구간별), 정상 요소와의 비교
"""
from __future__ import print_function
import sys

BUILD = '2026-10-02a (mesh-quality)'


def main():
    print('[MQ] build %s' % BUILD)
    if len(sys.argv) < 2:
        print('[MQ] usage: abaqus python mesh_quality_probe.py <odb>')
        return 1
    odb_path = sys.argv[1]
    try:
        from odbAccess import openOdb
    except Exception:
        print('[MQ] odbAccess 를 못 읽었습니다 - abaqus python 으로 실행해야 합니다.')
        return 1
    odb = openOdb(odb_path, readOnly=True)

    inst_name = None
    for k in odb.rootAssembly.instances.keys():
        if 'MEMBRANE' in k.upper():
            inst_name = k
            break
    if inst_name is None:
        print('[MQ] MEMBRANE 인스턴스를 못 찾음. 있는 것: %s' % list(odb.rootAssembly.instances.keys()))
        return 1
    inst = odb.rootAssembly.instances[inst_name]
    print('[MQ] 인스턴스 = %s / 요소 %d개 / 노드 %d개' % (inst_name, len(inst.elements), len(inst.nodes)))

    # 왜곡 요소 집합을 이름 변형으로 찾는다
    warn = None
    cand = []
    for key in odb.rootAssembly.elementSets.keys():
        cand.append(key)
        if 'DISTORT' in key.upper() or 'WARN' in key.upper():
            warn = odb.rootAssembly.elementSets[key]
            print('[MQ] 왜곡 집합 발견(어셈블리): %s' % key)
            break
    if warn is None:
        for key in inst.elementSets.keys():
            cand.append('inst:' + key)
            if 'DISTORT' in key.upper() or 'WARN' in key.upper():
                warn = inst.elementSets[key]
                print('[MQ] 왜곡 집합 발견(인스턴스): %s' % key)
                break
    if warn is None:
        print('[MQ] 왜곡 집합을 못 찾음. 있는 집합: %s' % cand[:20])
        print('[MQ] -> .dat 의 WarnElemDistorted 가 ODB 에 저장되지 않은 경우입니다.')
        return 1

    # 노드 라벨 -> 좌표 (한 번만 훑는다. getNodeFromLabel 의존 제거)
    nlab = {}
    for nd in inst.nodes:
        nlab[nd.label] = nd.coordinates
    print('[MQ] 노드 좌표 %d개 확보' % len(nlab))

    # 요소 라벨 -> 중심 좌표
    lab2c = {}
    missing = 0
    for el in inst.elements:
        cs = [nlab.get(n) for n in el.connectivity]
        cs = [c for c in cs if c is not None]
        if not cs:
            missing += 1
            continue
        cx = sum(c[0] for c in cs) / float(len(cs))
        cy = sum(c[1] for c in cs) / float(len(cs))
        lab2c[el.label] = (cx, cy)
    print('[MQ] 요소 중심 좌표 %d개 확보 (좌표 못 찾은 요소 %d개)' % (len(lab2c), missing))

    bad_labels = [e.label for e in warn.elements]
    print('[MQ] 왜곡 요소 수 = %d / 전체 %d (%.2f%%)'
          % (len(bad_labels), len(inst.elements), 100.0 * len(bad_labels) / max(1, len(inst.elements))))

    # 높이 구간별 분포 (apex 는 y=10, 밑변은 y=0)
    bins = [(0, 2), (2, 4), (4, 6), (6, 8), (8, 9), (9, 9.5), (9.5, 10.01)]
    print('[MQ] 높이(y) 구간별 왜곡 요소 수:')
    pts = [lab2c[l] for l in bad_labels if l in lab2c]
    tot = [lab2c[l] for l in lab2c.keys()]
    for lo, hi in bins:
        nb = sum(1 for (x, y) in pts if lo <= y < hi)
        nt = sum(1 for (x, y) in tot if lo <= y < hi)
        pct = (100.0 * nb / nt) if nt else 0.0
        print('    y=[%.1f,%.1f)  왜곡 %5d / 전체 %5d  = %6.2f %%' % (lo, hi, nb, nt, pct))
    if pts:
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        print('[MQ] 왜곡 요소 좌표 범위: x=[%.3f, %.3f]  y=[%.3f, %.3f]'
              % (min(xs), max(xs), min(ys), max(ys)))
        # 밑변 근처(모서리 45도) 비율
        near_base_corner = sum(1 for (x, y) in pts if y < 1.0)
        print('[MQ] 밑변 1 m 이내(y<1): %d개 (%.1f%%)'
              % (near_base_corner, 100.0 * near_base_corner / len(pts)))
    odb.close()
    return 0


if __name__ == '__main__':
    sys.exit(main())
