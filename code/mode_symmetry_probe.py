# -*- coding: utf-8 -*-
"""mode_symmetry_probe.py — 좌굴모드의 **좌우대칭성** 진단과 대칭/반대칭 정규화.

문제(2026-10-07 사용자 지적): 수렴에 성공한 모드들이 모두 클램프 근처에만 몰려 있거나,
수렴에 실패하거나, 양쪽에 대칭적으로 생기지 않고 **한쪽에만** 생긴다.

원인 후보 셋:
  (1) **모델 메쉬가 x=10 대칭이 아니다.** 현행 모델은 자유 메쉬다
      (`seedPart(BASE/SEED_DIV, deviationFactor=0.1)` + `QUAD_DOMINATED/FREE/MEDIAL_AXIS`,
       약 1.82만 요소 -> 인스턴스 노드 18442개). 자유 메쉬는 대칭 영역에서도 좌우대칭이 아니다.
      대조: `aba_grid_mesh.generate_grid` 의 격자 메쉬는 노드 10201 / 요소 10100 에서
      미러짝 누락 **0개**로 완전 대칭임을 확인했지만, 그 격자는 §10 이후 **미사용**이다.
  (2) **축퇴 고유공간의 기저 선택.** 좌/우 클램프 국소모드가 (거의) 축퇴하면 솔버는 그 2차원
      공간의 아무 기저나 돌려준다(보통 한쪽에 몰린 국소 벡터). 실측 x_c=0.35 에서
      λ1=λ2 가 정확히 같고 두 모드의 무게중심이 (4.352,2.440)/(4.352,2.439) 로 **같은 쪽**이었다.
  (3) **base state 자체의 대칭 깨짐.** 주름 개시 근처(DEAD_FRAC=0.5 가 개시 43.7~61 % 대)에서는
      비선형 해가 한쪽으로 분기할 수 있다.

이 도구는 (1)을 **정량화**하고 (2)(3)에 **정규화된 답**을 준다:
  - 노드 미러쌍 커버율 + 편차(최대/평균) [m]      <- 메쉬 비대칭 지표
  - 각 모드의 좌/우 에너지 분율(sidedness)        <- "한쪽에만 생겼다"의 정량 판정
  - 각 모드의 대칭조합 u+ 와 반대칭조합 u-          <- 양쪽이 모두 있는 표준 기저
      모델이 정확히 대칭이면 u+ 와 u- 는 **같은 λ 의 고유벡터**다(미러 연산자가 해밀토니안과
      교환하므로). 대칭이 조금 깨져도 섭동(imperfection) 형상으로는 더 정규적이고,
      x_c 를 바꿔가며 비교할 때 **기저 선택의 임의성**이 사라진다.

사용:
    abaqus python mode_symmetry_probe.py <odb> <step> [n_modes] [instance]
    abaqus python mode_symmetry_probe.py selftest
읽기 전용(라이선스 0). 순수 함수부(mirror_partner/symmetrize/sidedness/centroid_x)는
Abaqus 없이 selftest 로 검증된다.
"""
from __future__ import print_function
import math
import os
import sys
from collections import defaultdict

MIRROR_AXIS_X = 10.0     # 모델 좌우대칭축 (BASE/2). 미러는 x -> 2*axis - x, y 는 그대로.
MIRROR_TOL_M = 5.0e-3    # 미러쌍 매칭 허용오차 [m]. 자유 메쉬는 정확히 맞지 않으므로 넉넉히.
U3_EPS = 1.0e-30         # sidedness 분모 보호


def mirror_partner(coords, axis=MIRROR_AXIS_X, tol=MIRROR_TOL_M):
    """coords=[(x,y,z),...] -> (partner, stats).

    partner[i] = 미러상 (2*axis-x, y) 에 가장 가까운 노드 인덱스(tol 이내), 없으면 -1.
    stats = dict(n, miss, cover, dev_max, dev_mean)  — 메쉬 비대칭의 직접 지표.
    순수 파이썬(격자 해시로 후보 축소). numpy 불필요.
    """
    n = len(coords)
    partner = [-1] * n
    cell = max(tol, 1.0e-6) * 4.0
    grid = defaultdict(list)
    for i, p in enumerate(coords):
        grid[(int(math.floor(p[0] / cell)), int(math.floor(p[1] / cell)))].append(i)
    devs = []
    for i, p in enumerate(coords):
        mx, my = 2.0 * axis - p[0], p[1]
        cx, cy = int(math.floor(mx / cell)), int(math.floor(my / cell))
        best, bd = -1, None
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for j in grid.get((cx + dx, cy + dy), ()):
                    q = coords[j]
                    d = math.hypot(q[0] - mx, q[1] - my)
                    if bd is None or d < bd:
                        best, bd = j, d
        if bd is not None and bd <= tol:
            partner[i] = best
            devs.append(bd)
    miss = sum(1 for j in partner if j < 0)
    stats = dict(n=n, miss=miss,
                 cover=(float(n - miss) / n) if n else 0.0,
                 dev_max=(max(devs) if devs else None),
                 dev_mean=((sum(devs) / len(devs)) if devs else None))
    return partner, stats


def symmetrize(u3, partner):
    """u3 -> (sym, anti). 짝이 없으면(partner<0) 자기 값을 양쪽에 쓴다.

    sym  = (u + M u)/2   — 양쪽에 **같은 부호**로 나타나는 성분
    anti = (u - M u)/2   — 양쪽에 **반대 부호**로 나타나는 성분
    sym + anti == u  (짝이 있는 노드에서 정확히 성립)
    """
    n = len(u3)
    sym = [0.0] * n
    anti = [0.0] * n
    for i in range(n):
        j = partner[i]
        if j < 0:
            sym[i] = u3[i]
            anti[i] = u3[i]
        else:
            sym[i] = 0.5 * (u3[i] + u3[j])
            anti[i] = 0.5 * (u3[i] - u3[j])
    return sym, anti


def sidedness(u3, coords, axis=MIRROR_AXIS_X):
    """|u3| 에너지의 (좌분율, 우분율). 한쪽에만 몰리면 (1,0) 또는 (0,1) 에 가깝다."""
    L = R = 0.0
    for v, p in zip(u3, coords):
        w = v * v
        if p[0] < axis:
            L += w
        elif p[0] > axis:
            R += w
    t = L + R
    if t <= U3_EPS:
        return None, None
    return L / t, R / t


def centroid_x(u3, coords):
    """|u3| 가중 무게중심 x. 대칭조합이면 axis 에 와야 한다."""
    wsum = sum(abs(v) for v in u3)
    if wsum <= U3_EPS:
        return None
    return sum(abs(v) * p[0] for v, p in zip(u3, coords)) / wsum


def _fmt(v, nd=3):
    return '  -  ' if v is None else ('%.*f' % (nd, v))


def report(odb_path, step_name, n_modes=4, instance=None, axis=MIRROR_AXIS_X):
    """ODB 를 읽어 모드별 sidedness + 대칭/반대칭 조합 지표를 찍는다."""
    from buckle_mode_report import read_odb_modes, DEFAULT_INSTANCE
    if instance is None:
        instance = DEFAULT_INSTANCE
    nodes, disps, metrics, labels = read_odb_modes(
        odb_path, step_name, n_modes=n_modes, instance=instance,
        dat_hint=os.path.splitext(odb_path)[0] + '.dat')
    coords = [(p[0], p[1], p[2] if len(p) > 2 else 0.0) for p in nodes]
    partner, st = mirror_partner(coords, axis=axis)
    print('=' * 78)
    print('[SYM] %s / %s   노드 %d개 (인스턴스 %s)' % (os.path.basename(odb_path), step_name,
                                                       len(coords), instance))
    print('[SYM] 대칭축 x=%.4g  미러쌍 커버 %.4f (%d/%d)  편차 평균 %s m / 최대 %s m'
          % (axis, st['cover'], st['n'] - st['miss'], st['n'],
             _fmt(st['dev_mean'], 6), _fmt(st['dev_max'], 6)))
    if st['cover'] < 0.999:
        print('[SYM] >>> 메쉬가 x=%.4g 대칭이 아니다. 이 자체가 "한쪽에만 모드가 생기는" 원인 후보다.'
              % axis)
    else:
        print('[SYM] >>> 메쉬는 대칭이다. 한쪽 편향이 보이면 원인은 축퇴 기저/베이스 상태 쪽이다.')
    print('[SYM] 판정 기준: sidedness (L/R) 이 0.9/0.1 을 넘으면 **한쪽 편향**')
    print('')
    for k in range(len(disps)):
        u3 = [float(v) for v in disps[k]]
        L, R = sidedness(u3, coords, axis=axis)
        sy, an = symmetrize(u3, partner)
        tag = ''
        if L is not None:
            tag = '한쪽 편향' if max(L, R) >= 0.9 else ('대칭' if abs(L - R) < 0.2 else '중간')
        print('  M%-3d sidedness L/R = %s/%s  %s' % (k + 1, _fmt(L), _fmt(R), tag))
        print('        원래모드 centroid_x=%s  |  대칭조합 u+ centroid_x=%s  '
              '반대칭조합 u- centroid_x=%s'
              % (_fmt(centroid_x(u3, coords)), _fmt(centroid_x(sy, coords)),
                 _fmt(centroid_x(an, coords))))
    print('')
    print('[SYM] 권고: 섭동(imperfection)으로는 u+ / u- 를 쓰면 양쪽이 모두 살아 있고,')
    print('      x_c 를 바꿔가며 비교할 때 기저 선택의 임의성이 사라진다.')
    return 0


def selftest():
    """순수 함수만 시험한다(Abaqus 불필요)."""
    cnt = [0, 0]      # [ok, tot]  — 중첩 함수에서 py2 호환으로 리스트를 쓴다(global 금지)
    # --- 완전 대칭인 장난감 메쉬 (축 x=0) ---
    coords = [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0),
              (1.0, 1.0, 0.0), (-1.0, 1.0, 0.0),
              (0.0, 0.5, 0.0)]
    partner, st = mirror_partner(coords, axis=0.0, tol=1e-6)

    def chk(name, cond):
        cnt[1] += 1
        if cond:
            cnt[0] += 1
            print('  OK   %s' % name)
        else:
            print('  FAIL %s' % name)

    chk('대칭 메쉬: 미러쌍 누락 0', st['miss'] == 0 and st['cover'] == 1.0)
    chk('대칭축 위 노드는 자기 자신과 짝', partner[4] == 4)
    chk('거울 노드끼리 짝', partner[0] == 1 and partner[2] == 3)

    # --- 한쪽(+x)에만 있는 모드 ---
    u3 = [1.0, 0.0, 1.0, 0.0, 0.0]
    L, R = sidedness(u3, coords, axis=0.0)
    chk('한쪽 모드 sidedness = (0,1)', (L, R) == (0.0, 1.0))

    sy, an = symmetrize(u3, partner)
    chk('sym+anti == u', all(abs(sy[i] + an[i] - u3[i]) < 1e-12 for i in range(5)))
    Ls, Rs = sidedness(sy, coords, axis=0.0)
    chk('대칭조합은 좌우 반반 (0.5,0.5)', abs(Ls - 0.5) < 1e-12 and abs(Rs - 0.5) < 1e-12)
    chk('대칭조합 centroid_x == 0', abs(centroid_x(sy, coords) - 0.0) < 1e-12)
    chk('원래 모드 centroid_x == 1 (한쪽)', abs(centroid_x(u3, coords) - 1.0) < 1e-12)

    # --- 비대칭 메쉬는 누락으로 보고되어야 한다 ---
    bad = [(1.0, 0.0, 0.0), (-1.02, 0.0, 0.0)]      # 미러가 0.02 m 어긋남
    _, st2 = mirror_partner(bad, axis=0.0, tol=1e-3)
    chk('비대칭 메쉬는 누락 검출', st2['miss'] == 2 and st2['cover'] == 0.0)
    _, st3 = mirror_partner(bad, axis=0.0, tol=5e-2)
    chk('허용오차 내면 짝 성립', st3['miss'] == 0)

    # --- 기하 예외: u3 가 항등 0 ---
    Lz, Rz = sidedness([0.0] * 5, coords, axis=0.0)
    chk('u3==0 이면 sidedness None', Lz is None and Rz is None)
    chk('u3==0 이면 centroid None', centroid_x([0.0] * 5, coords) is None)

    print('')
    print('selftest: %d/%d 통과' % (cnt[0], cnt[1]))
    return 0 if cnt[0] == cnt[1] else 1


def main(argv):
    if len(argv) >= 2 and argv[1] == 'selftest':
        return selftest()
    if len(argv) < 3:
        print('사용: abaqus python mode_symmetry_probe.py <odb> <step> [n_modes] [instance]')
        print('      abaqus python mode_symmetry_probe.py selftest')
        return 2
    odb = argv[1]
    step = argv[2]
    nm = int(argv[3]) if len(argv) > 3 else 4
    inst = argv[4] if len(argv) > 4 else None
    try:
        return report(odb, step, n_modes=nm, instance=inst)
    except Exception as e:
        print('[SYM] 실패: %s: %s' % (type(e).__name__, e))
        print('[SYM] odbAccess 가 필요하므로 `abaqus python` 으로 실행해야 한다.')
        return 4


if __name__ == '__main__':
    sys.exit(main(sys.argv))
