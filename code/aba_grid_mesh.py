"""aba_grid_mesh.py — 삼각 막을 '내부 균일 정사각형 + 경계 삼각형' 격자로 채운다.

목적
    Galhofo(2022) 의 참조 메쉬를 재현한다. 논문은 `10100 quadratic thin-shell
    elements (STRI65+S8R5)` 를 쓰지만 Table A.1 에서 **같은 10,100 요소** 로
    `S3+S4`(1차 셸)도 보고했다(u_z,max 2.003e-4 m / 2x3 wrinkles). 즉 논문의
    요소 수는 1차/2차가 같은 메쉬 토폴로지다.
    이 모듈이 만드는 격자는 **S4 9900 + S3 200 = 10100** 으로 그 메쉬와 일치한다.

    현재 자유 메쉬(약 18,200 S4, deviationFactor 0.05 + ADVANCING_FRONT)는
    요소 수가 논문의 1.8배이고 조성도 다르다.

기하
    삼각형 정점 = (0,0), (BASE,0), (BASE/2,HEIGHT).   BASE=20, HEIGHT=10.
    세 변의 방정식:  y = 0  /  y = x  /  y = BASE - x
    클램프는 두 사선변(y=x, y=BASE-x) 위에 붙는다:  clamp_coord_L/R 참조.

방법
    h x h 셀 격자를 만들고, 각 셀을 삼각형의 3개 반평면으로 **정확히 클리핑**한다
    (Sutherland-Hodgman). 결과가
      - 4각형 -> S4  (내부 균일 정사각형)
      - 3각형 -> S3  (경계에서 잘린 직각/사선 삼각형)
      - 5각형 이상 -> 부채꼴 분할로 S3 여러 개  (정점 근처에서만 발생)
    이 방식은 면적을 **기계 정밀도로 보존**한다. 무게중심 판정 방식은 경계
    셀을 통째로 넣거나 빼서 면적 오차가 0.5 % 났다(실측, 폐기).

검증 (Abaqus 불필요)
    generate_grid(...) 는 순수 파이썬이며 (nodes, s4, s3) 를 돌려준다.
    - 면적 합 == BASE*HEIGHT/2 (오차 < 1e-12)
    - 요소 수 == 논문과 비교
    - 경계 노드가 y=0 / y=x / y=BASE-x 위 (tol=1e-4) 에 정확히 온다
      -> All_Edges_NoClamp 의 좌표 판정과 일관

주의 (Abaqus 측)
    CAE 에서는 seedPart/setMeshControls/generateMesh 경로를 쓰지 않고
    part.deleteMesh() -> part.addNodes(nodeData=...) -> part.addElements(...)
    로 직접 채운다. addNodes/addElements 는 **튜플**을 받는다(Set(nodes=) 와 달리
    MeshNodeArray 를 요구하지 않는다). 라벨은 1부터 연속으로 준다.
"""

TOL = 1e-9
BOUNDARY_TOL = 1.0e-4   # run_abaqus.py 의 All_Edges_NoClamp 와 같은 값

# 삼각형의 3개 반평면:  s * f(p) >= 0 이 내부
_HALFPLANES = (
    (lambda p: p[1],                 1.0),   # y >= 0
    (lambda p: p[1] - p[0],         -1.0),   # y - x <= 0        -> y <= x
    (lambda p: p[1] - (20.0 - p[0]),-1.0),   # y + x - 20 <= 0   -> y <= 20 - x
)


def _halfplanes(base):
    return (
        (lambda p: p[1],                    1.0),
        (lambda p: p[1] - p[0],            -1.0),
        (lambda p: p[1] - (base - p[0]),   -1.0),
    )


def _clip(poly, f, s):
    """반평면 s*f(p) >= 0 으로 볼록 다각형 클리핑 (Sutherland-Hodgman)."""
    out = []
    n = len(poly)
    for k in range(n):
        a, b = poly[k], poly[(k + 1) % n]
        fa, fb = s * f(a), s * f(b)
        if fa >= -TOL:
            out.append(a)
        if (fa > TOL and fb < -TOL) or (fa < -TOL and fb > TOL):
            t = fa / (fa - fb)
            out.append((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
    return out


def _area(poly):
    a = 0.0
    for k in range(len(poly)):
        x1, y1 = poly[k]
        x2, y2 = poly[(k + 1) % len(poly)]
        a += x1 * y2 - x2 * y1
    return abs(a) / 2.0


def generate_grid(base=20.0, height=10.0, seed_div=200.0):
    """순수 파이썬 격자 생성. (nodes, s4, s3) 반환.

    nodes : [(x, y), ...]  사용된 노드의 좌표 (요소에 등장하는 순서)
    s4    : [(i0,i1,i2,i3), ...]  nodes 인덱스
    s3    : [(i0,i1,i2), ...]     nodes 인덱스
    """
    h = float(base) / float(seed_div)
    nx = int(round(base / h))
    ny = int(round(height / h))
    hp = _halfplanes(base)

    idx = {}
    nodes = []

    def nid(p):
        key = (round(p[0], 10), round(p[1], 10))
        i = idx.get(key)
        if i is None:
            i = len(nodes)
            idx[key] = i
            nodes.append((key[0], key[1]))
        return i

    s4, s3 = [], []
    for j in range(ny):
        for i in range(nx):
            poly = [(i * h, j * h), ((i + 1) * h, j * h),
                    ((i + 1) * h, (j + 1) * h), (i * h, (j + 1) * h)]
            for f, s in hp:
                poly = _clip(poly, f, s)
                if len(poly) < 3:
                    break
            if len(poly) < 3 or _area(poly) < 1e-14:
                continue
            if len(poly) == 4:
                s4.append(tuple(nid(p) for p in poly))
            elif len(poly) == 3:
                s3.append(tuple(nid(p) for p in poly))
            else:                       # 5각형 이상 -> 부채꼴 분할
                a0 = nid(poly[0])
                for k in range(1, len(poly) - 1):
                    s3.append((a0, nid(poly[k]), nid(poly[k + 1])))
    return nodes, s4, s3


def grid_report(base=20.0, height=10.0, seed_div=200.0):
    """격자 통계 + 검증 요약을 문자열로 돌려준다 (로그/검증용)."""
    nodes, s4, s3 = generate_grid(base, height, seed_div)
    h = float(base) / float(seed_div)
    hp = _halfplanes(base)
    # 면적
    def parea(poly):
        return _area(poly)
    area = 0.0
    for t in s4:
        p = [nodes[k] for k in t]
        area += abs((p[1][0]-p[0][0])*(p[3][1]-p[0][1]) - (p[3][0]-p[0][0])*(p[1][1]-p[0][1]))
    for t in s3:
        p = [nodes[k] for k in t]
        area += _area(p)
    target = base * height / 2.0
    # 경계 노드
    nb = sum(1 for x, y in nodes
             if abs(y) < BOUNDARY_TOL
             or abs(y - x) < BOUNDARY_TOL * 1.5
             or abs(y - (base - x)) < BOUNDARY_TOL * 1.5)
    lines = [
        "[grid] h=%.6g  nx=%d ny=%d  nodes=%d  S4=%d  S3=%d  total=%d"
        % (h, int(round(base / h)), int(round(height / h)), len(nodes), len(s4), len(s3), len(s4) + len(s3)),
        "[grid] area=%.9f m^2  (target %.9f)  err=%.3e m^2" % (area, target, area - target),
        "[grid] boundary nodes (tol=%.1e): %d   <- y=0 / y=x / y=BASE-x 위"
        % (BOUNDARY_TOL, nb),
    ]
    return "\n".join(lines), dict(nodes=len(nodes), s4=len(s4), s3=len(s3),
                                  total=len(s4) + len(s3), area=area, target=target,
                                  boundary_nodes=nb)


def fill_part(part, base=20.0, height=10.0, seed_div=200.0, verbose=True):
    """Abaqus CAE Part 를 격자로 채운다 (seedPart 경로 대체).

      part.deleteMesh()
      part.addNodes(nodeData=((label,(x,y,z)), ...))
      part.addElements(elementData=((label,(n1,n2,n3,n4)), ...), type=S4)
      part.addElements(elementData=((label,(n1,n2,n3)),   ...), type=S3)

    addNodes/addElements 는 **튜플**을 받는다 (a.Set(nodes=) 와 달리 MeshNodeArray 를
    요구하지 않는다). 라벨은 1부터 연속으로 준다.
    """
    from abaqusConstants import S4, S3   # noqa: F401  (심볼)

    nodes, s4, s3 = generate_grid(base, height, seed_div)
    part.deleteMesh()
    part.addNodes(nodeData=tuple(
        (i + 1, (float(x), float(y), 0.0)) for i, (x, y) in enumerate(nodes)))
    if s4:
        part.addElements(elementData=tuple(
            (k + 1, tuple(n + 1 for n in t)) for k, t in enumerate(s4)), type=S4)
    if s3:
        part.addElements(elementData=tuple(
            (k + 1 + len(s4), tuple(n + 1 for n in t)) for k, t in enumerate(s3)), type=S3)
    if verbose:
        rep, _ = grid_report(base, height, seed_div)
        print(rep)
    return len(s4), len(s3), len(nodes)


if __name__ == '__main__':
    rep, st = grid_report()
    print(rep)
    print("\n[grid] 논문 대조: Galhofo Table A.1 '10100 elements' 의 S3+S4 행은"
          "\n       u_z,max 2.003e-4 m / 2x3 wrinkles. 격자가 10100 이면 같은 토폴로지다.")
