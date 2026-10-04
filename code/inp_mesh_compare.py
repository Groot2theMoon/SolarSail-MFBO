#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""두 Abaqus .inp 의 메쉬(노드/요소) 정합을 비교한다.  build 2026-10-04a (mesh-compare)

왜 필요한가 (실측 2026-10-04)
    *IMPERFECTION 은 **노드 라벨**로 주입된다. 모드 소스와 HF 의 메쉬가 다르면 엉뚱한 노드에
    섭동이 조용히 들어간다. 요소 타입을 바꾸면 CAE 메셔가 **요소 구성**을 바꾸는 관찰이 있었다
    (셸 덱 = S4 단일 / 막 덱 = M3D3+M3D4, 총 요소 수는 둘 다 18,200).
    그래서 '요소코드가 같다'가 아니라 **노드 집합과 좌표가 같은가**를 직접 본다.
    노드만 같으면 모드 소스=셸 / HF=막 하이브리드가 성립한다.

사용 (Abaqus 불필요, 표준 라이브러리만):
    python inp_mesh_compare.py <A.inp> <B.inp>
"""
from __future__ import print_function
import sys, os, re, hashlib

BUILD = '2026-10-04a (mesh-compare)'
_RE_NODE = re.compile(r'^\*Node\b')
_RE_ELEM = re.compile(r'^\*Element\s*,\s*type\s*=\s*([^,\s]+)', re.I)
_RE_KW = re.compile(r'^\*')
_RE_DATA = re.compile(r'^\s*([0-9]+)\s*,(.*)$')


def parse_inp(path):
    """-> (nodes{label:(x,y,z)}, elem_counts{type:int}, node_blocks:int, elem_blocks:int)"""
    nodes, elems = {}, {}
    nb = eb = 0
    mode, etype = None, None
    with open(path, 'r', errors='replace') as f:
        for line in f:
            m = _RE_NODE.match(line)
            if m:
                mode, etype, nb = 'node', None, nb + 1
                continue
            m = _RE_ELEM.match(line)
            if m:
                etype = m.group(1).strip()
                mode = 'elem'
                eb += 1
                elems.setdefault(etype, 0)
                continue
            if _RE_KW.match(line):
                mode, etype = None, None
                continue
            if mode is None:
                continue
            if mode == 'node':
                m = _RE_DATA.match(line)
                if not m:
                    continue
                lab = int(m.group(1))
                rest = m.group(2)
                try:
                    xyz = tuple(float(t) for t in rest.split(',')[:3])
                    if len(xyz) == 3:
                        nodes[lab] = xyz
                except ValueError:
                    pass                      # nset 등 좌표 없는 줄
            elif mode == 'elem':
                if _RE_DATA.match(line):
                    elems[etype] += 1
    return nodes, elems, nb, eb


def mesh_digest(nodes):
    """노드 라벨+좌표(1e-9 m 반올림)로 만든 지문. 좌표까지 같아야 같은 메쉬다."""
    h = hashlib.md5()
    for lab in sorted(nodes):
        h.update(('%d|%.9f|%.9f|%.9f\n' % ((lab,) + nodes[lab])).encode('utf-8'))
    return h.hexdigest()


def _show(path, nodes, elems, nb, eb):
    print('  파일        : %s' % path)
    print('  노드        : %d개 (블록 %d개, 라벨 %d~%d)'
          % (len(nodes), nb, min(nodes) if nodes else 0, max(nodes) if nodes else 0))
    print('  노드 지문   : %s' % mesh_digest(nodes))
    print('  요소        : %s (블록 %d개)'
          % (', '.join('%s=%d' % (k, v) for k, v in sorted(elems.items())), eb))


def main():
    print('=' * 78)
    print('메쉬 정합 비교  build %s' % BUILD)
    print('=' * 78)
    if len(sys.argv) < 3:
        print('사용: python inp_mesh_compare.py <A.inp> <B.inp>')
        return 2
    pa, pb = sys.argv[1], sys.argv[2]
    for q in (pa, pb):
        if not os.path.exists(q):
            print('!!! 파일 없음: %s' % q)
            return 2
    nA, eA, nbA, ebA = parse_inp(pa)
    nB, eB, nbB, ebB = parse_inp(pb)

    print('\n[A]'); _show(pa, nA, eA, nbA, ebA)
    print('\n[B]'); _show(pb, nB, eB, nbB, ebB)

    dA, dB = mesh_digest(nA), mesh_digest(nB)
    same_nodes = (dA == dB)
    only_A = sorted(set(nA) - set(nB))
    only_B = sorted(set(nB) - set(nA))
    same_elems = (eA == eB)

    print('\n' + '-' * 78)
    print('[판정]')
    print('  노드 집합+좌표 : %s' % ('**동일**' if same_nodes else '**다름**'))
    if not same_nodes:
        print('    A 에만 있는 라벨 %d개 %s' % (len(only_A), only_A[:8]))
        print('    B 에만 있는 라벨 %d개 %s' % (len(only_B), only_B[:8]))
    print('  요소 구성      : %s' % ('동일' if same_elems else '다름'))
    if not same_elems:
        print('    A: %s' % dict(sorted(eA.items())))
        print('    B: %s' % dict(sorted(eB.items())))

    print('\n[결론]')
    if same_nodes and same_elems:
        print('  RESULT:IDENTICAL — 메쉬가 완전히 같다. 모드 주입 안전.')
    elif same_nodes:
        print('  RESULT:NODES_ONLY — **노드는 같고 요소 구성만 다르다.**')
        print('    -> *IMPERFECTION 은 노드 라벨로 주입되므로 **모드 주입은 안전하다.**')
        print('    -> 모드 소스=셸 / HF=막 하이브리드가 성립한다.')
    else:
        print('  RESULT:DIFFERENT — 노드가 다르다. 모드 주입이 엉뚱한 노드로 간다.')
        print('    -> 두 모델의 요소코드를 같게 맞춰야 한다.')
    print('=' * 78)
    return 0 if same_nodes else 1


if __name__ == '__main__':
    sys.exit(main())
