# -*- coding: utf-8 -*-
"""probe_grid_import.py — CAE 에서 `.inp` 격자 import 경로만 **30초**에 검증한다.

왜 필요한가
    `Part.addNodes` / `addElements` / `deleteMesh` 는 **`odb.Part` 전용**이라
    `mdb.models[].Part` 에는 없다(실측 2026-10-06/10-07: AttributeError).
    그래서 CAE 스크립트로 격자를 모델에 넣는 경로는 `.inp` import 뿐이다:
        mdb.models[<name>].PartFromInputFile(inputFileName=<abs path to .inp>)
    그 경로가 **이 기계·이 Abaqus 버전에서** 실제로 되는지, 파트 이름·개수가 기대와
    맞는지를 큰 런(수십 분) 전에 확인한다.

사용:
    abaqus cae noGUI=probe_grid_import.py

기대 출력:
    [grid] PartFromInputFile grid_probe.inp -> part 'Membrane' (노드 10201, S4 9900, S3 200)
    [probe] 파트 노드 10201 / 요소 10100   (기대와 일치)
    [probe] x 범위 0.0000 .. 20.0000
    [probe] OK — 이 경로가 이 기계에서 작동한다
"""
from __future__ import print_function
import os

from abaqus import mdb                                     # noqa: F401
import aba_grid_mesh


def main():
    m = mdb.Model(name='ProbeGrid')
    inp = os.path.join(os.getcwd(), 'grid_probe.inp')
    part, n_nod, n_s4, n_s3 = aba_grid_mesh.import_grid_part(
        m, part_name='Membrane', base=20.0, height=10.0, seed_div=200.0, inp_path=inp)

    try:
        print('[probe] import 된 파트 이름 = %r  (이름 규칙을 여기서 확인한다)' % (part.name,))
    except Exception as e:
        print('[probe] 파트 이름 조회 생략 (%s)' % type(e).__name__)
    n_parts = len(part.elements)
    print('[probe] 파트 노드 %d / 요소 %d   (기대 노드 10201, 요소 10100)'
          % (len(part.nodes), n_parts))
    if n_nod != 10201 or n_s4 != 9900 or n_s3 != 200 or n_parts != 10100:
        print('[probe] !! 기대와 다르다 — .inp 또는 import 경로를 확인해야 한다')
        return 1

    xs = [nd.coordinates[0] for nd in part.nodes]
    print('[probe] x 범위 %.4f .. %.4f  (삼각형 폭 0..20)' % (min(xs), max(xs)))

    #   요소 타입: .inp 의 *Element, type=S4/S3 이 그대로 들어왔는지 (있으면 확인, 없으면 건너뜀)
    try:
        et = part.elements[0].type
        print('[probe] 첫 요소 타입 = %s' % (et,))
    except Exception as e:
        print('[probe] 요소 타입 조회 생략 (%s)' % type(e).__name__)

    #   좌우대칭: x -> 20-x 짝이 노드 집합에 있어야 한다(§19.13 의 전제)
    key = set((round(nd.coordinates[0], 6), round(nd.coordinates[1], 6)) for nd in part.nodes)
    miss = sum(1 for (x, y) in key if (round(20.0 - x, 6), y) not in key)
    print('[probe] 좌우대칭: 미러짝 누락 %d / %d  (0 이어야 한다)' % (miss, len(key)))

    print('[probe] OK — 이 경로가 이 기계에서 작동한다. 이제 본 런을 돌려도 좋다.')
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
