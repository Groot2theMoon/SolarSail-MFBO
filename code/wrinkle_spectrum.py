"""wrinkle_spectrum.py — 기존 ODB 에서 주름 파장/개수를 측정한다 (라이선스 소모 없음).

사용법 (code\\ 디렉터리에서):
    abaqus python wrinkle_spectrum.py aba\\HF_Postbuckle.odb
    abaqus python wrinkle_spectrum.py aba\\HF_Postbuckle.odb Step-ClampTension

목적: "포스트버클링이 어느 상태에서, 어떤 주름 패턴으로 막히는가"를 정량화한다.
  - 진행도(step time) + |u3|/t 진폭 분포
  - 변형 전 x-y 좌표로 bin-average 한 u3 맵을 ASCII 로 출력 -> 주름 개수/파장을 눈으로 확인
  - 중앙 가로줄의 부호 변화 횟수 = 웨이브 수 (논문의 'number of wrinkles' 와 같은 종류의 지표)

주의: 이 스크립트는 저자(에이전트)가 로컬에 Abaqus ODB 가 없어 실행 검증을 못 했다.
      그래서 모든 단계에서 진행 상황을 출력하고, 실패하면 어느 단계인지 바로 보이게 했다.
      ODB 객체 속성 접근은 최소화했다 (CAE Set.name 사고 재발 방지).
"""

import sys

T = 5.0e-6          # 막 두께 [m]
NC, NR = 64, 28     # ASCII 맵 격자 (가로 x, 세로 y)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'aba\\HF_Postbuckle.odb'
    print('=' * 78)
    print('주름 스펙트럼 판독: %s' % path)
    print('  두께 t = %.3e m / 격자 %dx%d' % (T, NC, NR))
    print('=' * 78)

    from odbAccess import openOdb
    odb = openOdb(path, readOnly=True)

    # --- 1) 스텝/프레임 나열 ---
    step_names = list(odb.steps.keys())
    print('[1] 스텝 목록:')
    for sn in step_names:
        print('    %-26s frames=%d' % (sn, len(odb.steps[sn].frames)))

    # --- 2) 판독 대상 스텝: 인자로 주면 그 스텝, 아니면 '프레임이 있는 마지막 스텝' ---
    if len(sys.argv) > 2:
        step_name = sys.argv[2]
    else:
        step_name = None
        for sn in step_names:
            if len(odb.steps[sn].frames) > 0:
                step_name = sn
    if step_name is None or step_name not in step_names:
        raise RuntimeError('판독할 스텝을 정하지 못했습니다. 인자로 지정하세요: %s' % step_names)
    frames = odb.steps[step_name].frames
    if len(frames) == 0:
        raise RuntimeError("스텝 '%s' 에 프레임이 없습니다 (number frames: 0)." % step_name)
    fr = frames[-1]
    print('[2] 판독 대상: %s / 마지막 프레임  step time = %s  (frame %d/%d)'
          % (step_name, fr.frameValue, len(frames), len(frames)))

    # --- 3) MEMBRANE 인스턴스 찾기 (이름은 키에서 가져온다: .name 접근 회피) ---
    inst_key = None
    for nm in list(odb.rootAssembly.instances.keys()):
        if 'MEMBRANE' in nm.upper():
            inst_key = nm
            break
    if inst_key is None:
        raise RuntimeError('MEMBRANE 인스턴스가 없습니다: %s'
                           % list(odb.rootAssembly.instances.keys()))
    inst = odb.rootAssembly.instances[inst_key]
    print('[3] 인스턴스: %s  nodes=%d' % (inst_key, len(inst.nodes)))

    # --- 4) 변형 전 좌표 (binning 기준) ---
    coords = {}
    xmin = ymin = 1.0e30
    xmax = ymax = -1.0e30
    for n in inst.nodes:
        x = float(n.coordinates[0])
        y = float(n.coordinates[1])
        coords[n.label] = (x, y)
        if x < xmin:
            xmin = x
        if x > xmax:
            xmax = x
        if y < ymin:
            ymin = y
        if y > ymax:
            ymax = y
    sx = (xmax - xmin) or 1.0
    sy = (ymax - ymin) or 1.0
    print('[4] 형상 범위(변형 전): x [%.4g, %.4g]  y [%.4g, %.4g]'
          % (xmin, xmax, ymin, ymax))

    # --- 5) u3 추출 ---
    if 'U' not in fr.fieldOutputs.keys():
        raise RuntimeError("프레임에 'U' 필드가 없습니다: %s" % list(fr.fieldOutputs.keys()))
    fo = fr.fieldOutputs['U']
    try:
        sub = fo.getSubset(region=inst, position=2)   # position=NODAL
    except Exception:
        sub = fo.getSubset(region=inst)
    u3 = {}
    for v in sub.values:
        d = v.data
        if len(d) < 3:
            continue
        u3[v.nodeLabel] = float(d[2])
    print('[5] u3 값 개수 = %d  (좌표 보유 노드 %d)' % (len(u3), len(coords)))
    if not u3:
        raise RuntimeError('u3 값이 하나도 없습니다.')

    # --- 6) 통계 ---
    vals = sorted(u3.values(), key=lambda z: abs(z))
    n = len(vals)

    def q(p):
        return vals[min(n - 1, int(p * n))]

    print('[6] |u3| 통계 (n=%d): max=%.6e (%.1f t)  p99=%.3e  p50=%.3e'
          % (n, vals[-1], vals[-1] / T, q(0.99), q(0.50)))
    for k in (2, 10, 50, 200):
        c = sum(1 for z in vals if z > k * T)
        print('       |u3| > %4dt : %7d  (%6.2f%%)' % (k, c, 100.0 * c / n))

    # --- 7) ASCII 맵 (bin 평균) ---
    acc = {}
    for lab, xy in coords.items():
        if lab not in u3:
            continue
        i = min(NC - 1, int((xy[0] - xmin) / sx * NC))
        j = min(NR - 1, int((xy[1] - ymin) / sy * NR))
        a = acc.get((i, j))
        if a is None:
            acc[(i, j)] = [u3[lab], 1]
        else:
            a[0] += u3[lab]
            a[1] += 1

    print('[7] u3 맵 (bin 평균, 세로=+y, 가로=+x)')
    print('    범례: . |u|<2t   - / +  2~10t   x / o  10~50t   @ / #  >=50t')
    print('    ' + '-' * NC)
    for j in range(NR - 1, -1, -1):
        row = []
        for i in range(NC):
            a = acc.get((i, j))
            if a is None or a[1] == 0:
                row.append(' ')
                continue
            m = a[0] / a[1]
            r = abs(m) / T
            if r < 2.0:
                row.append('.')
            elif r < 10.0:
                row.append('+' if m > 0 else '-')
            elif r < 50.0:
                row.append('o' if m > 0 else 'x')
            else:
                row.append('#' if m > 0 else '@')
        print('    ' + ''.join(row))
    print('    ' + '-' * NC)

    # --- 8) 중앙 가로줄/세로줄의 부호 변화 = 웨이브 수 ---
    jm = NR // 2
    im = NC // 2
    row_x = [acc[(i, jm)][0] / acc[(i, jm)][1] for i in range(NC) if (i, jm) in acc]
    col_y = [acc[(im, j)][0] / acc[(im, j)][1] for j in range(NR) if (im, j) in acc]
    wx = sum(1 for k in range(1, len(row_x)) if row_x[k] * row_x[k - 1] < 0)
    wy = sum(1 for k in range(1, len(col_y)) if col_y[k] * col_y[k - 1] < 0)
    print('[8] 부호 변화(웨이브 추정):')
    print('       가로 중앙줄: 유효 bin=%2d  부호변화=%2d' % (len(row_x), wx))
    print('       세로 중앙줄: 유효 bin=%2d  부호변화=%2d' % (len(col_y), wy))
    if wx > 0:
        lam = 2.0 * sx / (wx + 1)
        print('       가로 추정 파장 ~ %.4g m = %.0f t' % (lam, lam / T))
    if wy > 0:
        lam = 2.0 * sy / (wy + 1)
        print('       세로 추정 파장 ~ %.4g m = %.0f t' % (lam, lam / T))

    print()
    print('[해석 기준] 웨이브 수 1~2 = 긴 파장 소수 주름(국소 폴드 위험) /')
    print('            수십 개 = 다중 웨이브(논문 모드 seed 상태, 분산된 slack).')
    odb.close()


main()
