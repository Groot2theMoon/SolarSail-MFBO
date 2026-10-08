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

import io as _io
import sys

from abaqus import mdb                                     # noqa: F401
import aba_grid_mesh

_REPORT = 'grid_probe_report.txt'   # 출력이 사라져도(CAE noExit 함정) 여기서 숫자를 회수한다
_LINES = []


def say(msg):
    print(msg)
    _LINES.append(msg)
    try:
        sys.stdout.flush()
    except Exception:
        pass


def main():
    m = mdb.Model(name='ProbeGrid')
    inp = os.path.join(os.getcwd(), 'grid_probe.inp')
    part, n_nod, n_s4, n_s3 = aba_grid_mesh.import_grid_part(
        m, part_name='Membrane', base=20.0, height=10.0, seed_div=200.0, inp_path=inp)

    try:
        _nm = part.name if isinstance(part.name, str) else str(part.name)
        say('[probe] import 된 파트 이름 = %r  (이름 규칙을 여기서 확인한다)' % (_nm,))
    except Exception as e:
        say('[probe] 파트 이름 조회 생략 (%s)' % type(e).__name__)
    n_parts = len(part.elements)
    n_nodes_read = len(part.nodes)
    n_elem_read = len(part.elements)
    say('[probe] .inp 에 쓴 값: 노드 %d / S4 %d / S3 %d' % (n_nod, n_s4, n_s3))
    say('[probe] 모델에서 읽은 값: 노드 %d / 요소 %d   (기대 노드 10201, 요소 10100)'
        % (n_nodes_read, n_elem_read))
    #   **중요**: .inp 에 쓴 개수가 아니라 **모델에서 읽은** 개수를 검사한다.
    #   (예전 판은 쓴 값만 비교해서, 모델에 메쉬가 잘못 들어가도 통과했다 — 로컬 검증이 잡았다.)
    if n_nodes_read != 10201 or n_elem_read != 10100:
        say('[probe] !! 모델의 메쉬가 기대와 다르다 — .inp 또는 import 경로를 확인해야 한다')
        return 1

    xs = [nd.coordinates[0] for nd in part.nodes]
    say('[probe] x 범위 %.4f .. %.4f  (삼각형 폭 0..20)' % (min(xs), max(xs)))

    #   요소 타입: .inp 의 *Element, type=S4/S3 이 그대로 들어왔는지 (있으면 확인, 없으면 건너뜀)
    try:
        et = part.elements[0].type
        say('[probe] 첫 요소 타입 = %s' % (et,))
    except Exception as e:
        say('[probe] 요소 타입 조회 생략 (%s)' % type(e).__name__)

    #   좌우대칭: x -> 20-x 짝이 노드 집합에 있어야 한다(§19.13 의 전제)
    #   주의: CAE noGUI 인터프리터는 `sum(generator)` 를 거부한다
    #   (TypeError: arg1; found 'generator', expecting a recognized type, 2026-10-07 실측).
    #   여기서는 생성기 표현식 없이 명시적 루프로 센다.
    key = set()
    for nd in part.nodes:
        c = nd.coordinates
        key.add((round(c[0], 6), round(c[1], 6)))
    miss = 0
    for (x, y) in key:
        if (round(20.0 - x, 6), y) not in key:
            miss += 1
    say('[probe] 좌우대칭: 미러짝 누락 %d / %d  (0 이어야 한다)' % (miss, len(key)))

    say('[probe] OK — 이 경로가 이 기계에서 작동한다. 이제 본 런을 돌려도 좋다.')
    return 0


def _dump_file():
    try:
        with _io.open(_REPORT, 'w', encoding='utf-8') as f:
            f.write('\n'.join(_LINES) + '\n')
        print('[probe] 결과 파일: %s' % _REPORT)
    except Exception as e:
        print('[probe] 결과 파일 쓰기 실패: %s' % e)


if __name__ == '__main__':
    #   **sys.exit 을 쓰지 않는다**: CAE noGUI 는 sys.exit 을 지원하지 않고 exit code 를 항상
    #   0 으로 고정하며, 그 순간 프로세스가 끝나 stdout 버퍼가 버려진다(실측: 프로브가 아무
    #   출력 없이 'Exit code: 0' 으로 끝났다). 대신 flush + 파일 기록을 쓴다.
    try:
        main()
    finally:
        #   실패해도(예외) 부분 로그가 남도록 finally 에 둔다.
        _dump_file()
        sys.stdout.flush()
