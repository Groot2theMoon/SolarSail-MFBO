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

BUILD = '2026-10-05a (mesh-count)'
_RE_NODE = re.compile(r'^\*Node\b')
_RE_ELEM = re.compile(r'^\*Element\s*,\s*type\s*=\s*([^,\s]+)', re.I)
_RE_KW = re.compile(r'^\*')
_RE_DATA = re.compile(r'^\s*([0-9]+)\s*,(.*)$')


def parse_inp(path):
    """-> (nodes{label:set(coords)}, elems{type:int}, node_blocks, elem_blocks, blocks[(i,n)])"""
    # [2026-10-05] 라벨 -> 좌표 **집합**.
    #   함정 1: 어셈블리 .inp 에는 파트 레벨 + 인스턴스 레벨 노드 블록이 함께 있어 같은 라벨이
    #           여러 번 나온다(실측: 모드 7블록 / HF 13블록). dict 로 덮어쓰면 마지막 블록 값이
    #           남아 정합하는 두 메쉬를 오탐한다.
    #   함정 2(더 근본적): **각 파트는 독립적인 노드 라벨 공간을 갖는다**(멤브레인도 1번부터,
    #           케이블도 1번부터). 라벨로 합치면 다른 파트의 노드가 섞여 좌표 비교가 무의미해진다
    #           (실측: 정합하는 두 덱에서 "좌표 일치 4.78%"). => 좌표는 판정에 쓰지 않는다.
    nodes, elems = {}, {}
    nb = eb = 0
    blocks = []              # [(블록번호, 노드수)] — 부품별 구조를 눈으로 보기 위함
    cur_n = 0
    mode, etype = None, None
    with open(path, 'r', errors='replace') as f:
        for line in f:
            m = _RE_NODE.match(line)
            if m:
                if mode == 'node':
                    blocks.append((nb, cur_n))
                mode, etype, nb, cur_n = 'node', None, nb + 1, 0
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
                        cur_n += 1
                except ValueError:
                    pass                      # nset 등 좌표 없는 줄
            elif mode == 'elem':
                if _RE_DATA.match(line):
                    elems[etype] += 1
    if mode == 'node':
        blocks.append((nb, cur_n))
    return nodes, elems, nb, eb, blocks


def _struct(elems):
    """셸/막 요소만. 케이블(T3D2) 등은 모델마다 파트 수가 달라도 메쉬 정합과 무관하다."""
    return {k: v for k, v in elems.items()
            if k.startswith('S') or k.startswith('M3D') or k.startswith('STRI')}


def _show(path, nodes, elems, nb, eb, blocks):
    print('  파일        : %s' % path)
    print('  노드        : %d개 (블록 %d개, 라벨 %d~%d)'
          % (len(nodes), nb, min(nodes) if nodes else 0, max(nodes) if nodes else 0))
    print('  노드 블록별 : %s' % (', '.join('%d개' % c for _, c in blocks) or '-'))
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
    nA, eA, nbA, ebA, blkA = parse_inp(pa)
    nB, eB, nbB, ebB, blkB = parse_inp(pb)
    print('\n[A]'); _show(pa, nA, eA, nbA, ebA, blkA)
    print('\n[B]'); _show(pb, nB, eB, nbB, ebB, blkB)

    labA, labB = set(nA), set(nB)
    only_A = sorted(labA - labB)
    only_B = sorted(labB - labA)
    sA, sB = _struct(eA), _struct(eB)

    # [2026-10-05 교정] 좌표 비교를 판정에서 뺐다. .inp 에서 **각 파트는 독립적인 노드 라벨 공간**을
    #   가진다(멤브레인도 1번부터, 케이블도 1번부터). 라벨로 합치면 서로 다른 파트의 노드가 섞여
    #   "좌표 일치 4.78%" 같은 무의미한 값이 나온다(실측). 신뢰할 수 있는 지표만 쓴다:
    #     (1) 총 노드 수  (2) 셸/막 요소 수  (3) A 에만 / B 에만 라벨 수
    #   **최종 기능 판정은 HF 자신의 로그다**:
    #       [IMPERFECTION] 기하 섭동 적용: <찾은>/<전체> 노드
    #   가 N/N 이면 모드 소스의 라벨을 전부 찾았다는 뜻 = 주입 매핑이 성립한다.
    print('\n' + '-' * 78)
    print('[판정]  (좌표 비교는 다부품 .inp 에서 성립하지 않아 제외한다 — 위 주석 참조)')
    print('  노드 총수      : A %d / B %d  %s'
          % (len(labA), len(labB), '일치' if len(labA) == len(labB) else '**다름**'))
    print('    A 에만 %d개 %s' % (len(only_A), only_A[:6]))
    print('    B 에만 %d개 %s' % (len(only_B), only_B[:6]))
    print('  셸/막 요소     : A %s / B %s  %s'
          % (dict(sorted(sA.items())), dict(sorted(sB.items())), '일치' if sA == sB else '**다름**'))
    print('  전체 요소      : A %s / B %s  (케이블 파트 수 차이는 무관)'
          % (dict(sorted(eA.items())), dict(sorted(eB.items()))))
    counts_ok = (len(labA) == len(labB)) and (sA == sB)
    print('\n[결론]')
    if counts_ok and (eA == eB):
        print('  RESULT:IDENTICAL_COUNTS — 노드 총수와 전체 요소가 같다.')
    elif counts_ok:
        print('  RESULT:COUNT_MATCH — **노드 총수와 셸/막 요소 수가 같다.**')
        print('    -> 임퍼펙션 주입의 필요조건을 만족한다. 전체 요소 차이는 케이블 파트 수 차이다.')
    else:
        print('  RESULT:COUNT_DIFFER — 노드 총수 또는 셸/막 요소 수가 다르다.')
        print('    -> 두 모델의 메쉬 설정/덱 생성 시점을 맞춰야 한다.')
    print('\n  ※ 최종 기능 판정은 HF 로그의 아래 줄이다(이 도구보다 강한 증거):')
    print('       [IMPERFECTION] 기하 섭동 적용: <찾은>/<전체> 노드')
    print('     가 <찾은> == <전체> 이면 모드 소스의 라벨을 전부 찾은 것 = 주입 매핑 성립.')
    print('=' * 78)
    return 0 if counts_ok else 1


if __name__ == '__main__':
    sys.exit(main())
