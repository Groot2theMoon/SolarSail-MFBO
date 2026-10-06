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


def boundary_node_labels(instance, base=20.0, height=10.0,
                         tol=BOUNDARY_TOL, exclude_labels=()):
    """인스턴스에서 삼각형 3변(y=0 / y=x / y=BASE-x) 위의 노드 라벨을 고른다.

    edge->노드 API 를 쓰지 않고 좌표로 판정한다. 이유 두 가지:
      (1) CAE 버전에 따라 edge->노드 API 가 다르다.
      (2) fill_part 의 orphan mesh 에서는 인스턴스의 기하 edge 가 비어 있다.
    exclude_labels 에 든 라벨은 뺀다(클램프 부착 구간 제외용).

    **세 스크립트(HF / mode / buckle)가 반드시 이 함수를 써야** 같은 노드 집합이
    나온다. 한 곳에서 손으로 판정하면 tolerance 나 조건이 어긋나 모드 이식이 깨진다.

    판정은 격자에서 정확하다: generate_grid 는 경계 노드를 사선 위에 정확히 놓고,
    격자 h=0.1 일 때 tol=1e-4 로 400개(y=0 201 / y=x 101 / y=BASE-x 101, 꼭짓점
    3개 중복)가 잡힌다(실측).
    """
    ex = set(exclude_labels)
    out = []
    for n in instance.nodes:
        if n.label in ex:
            continue
        x, y = n.coordinates[0], n.coordinates[1]
        if (abs(y) < tol) or (abs(y - x) < tol * 1.5) or (abs(y - (base - x)) < tol * 1.5):
            out.append(n.label)
    return out


def clamp_exclude_labels(instance, centers, radius):
    """클램프/패치 중심들에서 반경 radius 내 노드 라벨(면외 구속에서 제외할 대상)."""
    ex = set()
    for c in centers:
        for n in instance.nodes.getByBoundingSphere(center=c, radius=radius):
            ex.add(n.label)
    return ex


def make_set(assembly, name, instance, labels):
    """라벨 리스트로 어셈블리 노드셋을 만든다.

    a.Set(nodes=...) 는 MeshNodeArray 를 요구한다 — tuple/list 를 넘기면
    'Feature creation failed.' 로 죽는다(실측 2026-10-06). sequenceFromLabels 로 만든다.
    빈 라벨이면 RuntimeError 를 올린다(조용히 빈 셋을 만들면 BC 가 아무 데도
    안 걸려 물리적으로 다른 문제가 된다).
    """
    if not labels:
        raise RuntimeError("%s: 노드를 하나도 찾지 못했다. 좌표/tol 을 확인하라." % name)
    assembly.Set(name=name, nodes=instance.nodes.sequenceFromLabels(tuple(labels)))
    return assembly.sets[name]


def fill_part(part, base=20.0, height=10.0, seed_div=200.0,
              elem_quad=None, elem_tri=None, verbose=True):
    """Abaqus CAE Part 를 격자로 채운다 (seedPart 경로 대체).

      part.deleteMesh()
      part.addNodes(nodeData=((label,(x,y,z)), ...))
      part.addElements(elementData=((label,(n1,n2,n3,n4)), ...), type=elem_quad)
      part.addElements(elementData=((label,(n1,n2,n3)),   ...), type=elem_tri)

    addNodes/addElements 는 **튜플**을 받는다 (a.Set(nodes=) 와 달리 MeshNodeArray 를
    요구하지 않는다). 라벨은 1부터 연속으로 준다.

    elem_quad / elem_tri 는 호출 스크립트의 ELEM_CODE_QUAD / ELEM_CODE_TRI 를
    넘긴다. 여기서 S4/S3 를 하드코딩하면 스크립트의 상수와 어긋날 수 있어
    인자로 받는다(스크립트마다 요소 코드를 바꾸는 실험이 있었다).

    반환: (n_S4, n_S3, n_nodes)

    주의: 이 함수는 **orphan mesh** 를 만든다. 즉 요소가 파트 기하(face)에 붙지
    않는다. 따라서 호출 측은 섹션 할당을 반드시 요소 기반으로 해야 한다:
        p.SectionAssignment(region=p.Set(elements=p.elements, name='All'), ...)
    faces 기반(p.Set(faces=p.faces, ...))으로 두면 요소가 섹션을 못 받아
    'N elements have missing property definitions' 로 입력 단계에서 죽는다.
    같은 이유로 어셈블리 인스턴스의 .edges 도 비게 되므로 geometry 기반 셋
    (All_Edges)은 쓰면 안 된다.
    """
    if elem_quad is None or elem_tri is None:
        from abaqusConstants import S4, S3
        elem_quad = elem_quad if elem_quad is not None else S4
        elem_tri = elem_tri if elem_tri is not None else S3

    nodes, s4, s3 = generate_grid(base, height, seed_div)
    part.deleteMesh()
    part.addNodes(nodeData=tuple(
        (i + 1, (float(x), float(y), 0.0)) for i, (x, y) in enumerate(nodes)))
    if s4:
        part.addElements(elementData=tuple(
            (k + 1, tuple(n + 1 for n in t)) for k, t in enumerate(s4)), type=elem_quad)
    if s3:
        part.addElements(elementData=tuple(
            (k + 1 + len(s4), tuple(n + 1 for n in t)) for k, t in enumerate(s3)), type=elem_tri)
    if verbose:
        rep, _ = grid_report(base, height, seed_div)
        print(rep)
    return len(s4), len(s3), len(nodes)


def write_inp(path, base=20.0, height=10.0, seed_div=200.0):
    """격자를 최소 `.inp` 로 쓴다 (*Node + *Element 만, 셋/섹션 없음).

    용도: mdb 파트에 격자를 넣는 경로를 찾기 위한 재료다. `Part.addNodes` 는 존재하지
    않으므로(odb.Part 전용) CAE 는 `.inp` import 로 orphan mesh part 를 만드는 쪽이
    후보다. 이 파일을 CAE GUI 의 File > Import > Part 에 넣어 보면 된다.
    **정확한 스크립트 API 이름은 아직 확정하지 못했다** — CAE 에서
    File > Macro > Record 로 그 import 를 녹화하면 그 자리에서 확인된다.

    반환: (노드 수, S4 수, S3 수)
    """
    nodes, s4, s3 = generate_grid(base, height, seed_div)
    with io.open(path, 'w', encoding='utf-8') as f:
        f.write("*Heading\n")
        f.write("** uniform grid: interior quads + edge right triangles"
                " (Galhofo 10100 topology)\n")
        f.write("*Node\n")
        for i, (x, y) in enumerate(nodes):
            f.write("%d, %.10g, %.10g, 0.0\n" % (i + 1, x, y))
        f.write("*Element, type=S4\n")
        for k, t in enumerate(s4):
            f.write("%d, %s\n" % (k + 1, ", ".join(str(n + 1) for n in t)))
        f.write("*Element, type=S3\n")
        for k, t in enumerate(s3):
            f.write("%d, %s\n" % (k + 1 + len(s4), ", ".join(str(n + 1) for n in t)))
    return len(nodes), len(s4), len(s3)


if __name__ == '__main__':
    import sys as _sys
    if '--write' in _sys.argv:
        _out = 'grid_10100.inp'
        _n, _q, _t = write_inp(_out)
        print("[grid] wrote %s  nodes=%d S4=%d S3=%d" % (_out, _n, _q, _t))
    rep, st = grid_report()
    print(rep)
    print("\n[grid] 논문 대조: Galhofo Table A.1 '10100 elements' 의 S3+S4 행은"
          "\n       u_z,max 2.003e-4 m / 2x3 wrinkles. 격자가 10100 이면 같은 토폴로지다.")
    print("\n[grid] CAE import 시험용 .inp 를 쓰려면: python aba_grid_mesh.py --write")
