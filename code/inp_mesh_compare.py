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
    # [2026-10-05] 라벨 -> 좌표 **집합**. 왜 집합인가: 어셈블리 .inp 에는 파트 레벨 노드 블록과
    #   인스턴스 레벨 노드 블록이 함께 있어 **같은 라벨이 여러 번** 나온다(실측: 모드 7블록 /
    #   HF 13블록). dict 로 덮어쓰면 '어느 블록이 마지막이냐'에 따라 지문이 달라져
    #   **정합하는 두 메쉬를 '다름'으로 오탐**했다(실측 2026-10-05).
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
                        nodes.setdefault(lab, set()).add(xyz)
                except ValueError:
                    pass                      # nset 등 좌표 없는 줄
            elif mode == 'elem':
                if _RE_DATA.match(line):
                    elems[etype] += 1
    return nodes, elems, nb, eb


def mesh_digest(nodes):
    """라벨별 좌표 집합을 정렬해 만든 지문. 어느 블록이 먼저/나중이든 결과가 같다."""
    h = hashlib.md5()
    for lab in sorted(nodes):
        for c in sorted(nodes[lab]):
            h.update(('%d|%.9f|%.9f|%.9f\n' % ((lab,) + c)).encode('utf-8'))
    return h.hexdigest()


def _struct(elems):
    """구조 요소(셸/막)만. 케이블(T3D2) 등은 모델마다 파트 수가 달라도 메쉬 정합과 무관하다."""
    return {k: v for k, v in elems.items()
            if k.startswith('S') or k.startswith('M3D') or k.startswith('STRI')}


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

    labA, labB = set(nA), set(nB)
    only_A = sorted(labA - labB)
    only_B = sorted(labB - labA)
    common = labA & labB
    # 겹치는 라벨에 대해 '좌표 집합의 교집합이 비어있지 않은' 라벨의 비율
    n_coord_ok = sum(1 for l in common if nA[l] & nB[l])
    coord_frac = (n_coord_ok / len(common)) if common else 0.0
    sA, sB = _struct(eA), _struct(eB)

    print('\n' + '-' * 78)
    print('[판정]')
    print('  노드 라벨      : A %d개 / B %d개' % (len(labA), len(labB)))
    print('    A 에만 %d개 %s' % (len(only_A), only_A[:6]))
    print('    B 에만 %d개 %s' % (len(only_B), only_B[:6]))
    print('  좌표 일치      : 겹치는 라벨 중 %d/%d = %.2f%%'
          % (n_coord_ok, len(common), 100.0 * coord_frac))
    print('  구조 요소      : A %s / B %s' % (dict(sorted(sA.items())), dict(sorted(sB.items()))))
    print('  전체 요소      : A %s / B %s' % (dict(sorted(eA.items())), dict(sorted(eB.items()))))

    labels_ok = (not only_A) and (not only_B)
    struct_ok = (sA == sB)
    coords_ok = (coord_frac >= 1.0)
    full_ok = labels_ok and struct_ok and coords_ok and (eA == eB)

    print('\n[결론]')
    if full_ok:
        print('  RESULT:IDENTICAL — 라벨/좌표/요소가 모두 같다.')
    elif labels_ok and struct_ok and coords_ok:
        print('  RESULT:NODES_MATCH — **노드 라벨과 좌표, 셸/막 요소 수가 같다.**')
        print('    -> *IMPERFECTION 은 노드 라벨로 주입되므로 **주입 안전하다.**')
        print('    -> 전체 요소의 차이는 케이블 파트 수 같은 구조적 차이다(메쉬 무관).')
    elif labels_ok and coords_ok:
        print('  RESULT:STRUCT_DIFFER — 노드는 같지만 셸/막 요소 수가 다르다: %s vs %s'
              % (dict(sorted(sA.items())), dict(sorted(sB.items()))))
    else:
        print('  RESULT:DIFFERENT — 노드가 다르다. 모드 주입이 엉뚱한 노드로 간다.')
        print('    -> 두 모델의 요소코드/메쉬 설정을 같게 맞춰야 한다.')
    print('=' * 78)
    return 0 if (labels_ok and struct_ok and coords_ok) else 1


if __name__ == '__main__':
    sys.exit(main())
