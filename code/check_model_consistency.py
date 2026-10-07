"""check_model_consistency.py — 좌굴 모델과 HF 모델의 '모델 정의' 일치를 검사한다.

왜 필요한가
    run_abaqus_buckle.py 는 HF(run_abaqus.py)에서 케이블만 제외한 **별개 스크립트**다.
    스크립트를 분리하면 모델 정의(메쉬/형상/패치/재질/프리텐션)가 따로 놀 수 있다.
    그런데 좌굴 모드는 **HF 와 같은 노드**에 정의되어야 *IMPERFECTION 으로 이식된다.
    -> 한쪽만 수정하면 모드가 조용히 어긋나고, 그 사실은 런(라이선스)을 낭비한 뒤에야 드러난다.
    그래서 이 검사를 **제출 전에** 돌린다.

검사 방식 (주석을 코드로 오인하지 않기 위해)
    두 파일을 ast 로 파싱해 **docstring 을 비우고 주석을 제거한 뒤**(= 살아있는 코드만),
    공백을 모두 없앤 정규형으로 바꿔 항목별 포함 여부를 비교한다.
    -> "# radius=0.2 처럼 주석에만 있는 문구" 가 통과시키는 false OK 를 막고,
    -> "Step-Trigger -> Step-Buckle" 같은 설명문이 금지어 검사에 걸리는 false FAIL 도 막는다.

사용법
    python3 code/check_model_consistency.py

종료 코드
    0 = 일치 / 1 = 불일치 (어느 줄이 다른지 출력)
"""

import ast
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# [2026-09-28 정정] 검사 대상은 **실제로 제출하는 HF** 다.
#   예전에는 run_abaqus_new.py 를 봤다 — 그 파일은 legacy(1차 요소 S3/S4)이고 아무도 제출하지 않는다.
#   그래서 run_abaqus.py 의 모델 정의(요소/메쉬/결합)를 바꿔도 검사기가 PASS 를 냈다(가짜 초록불).
NEW = os.path.join(HERE, 'run_abaqus.py')
BUCKLE = os.path.join(HERE, 'run_abaqus_buckle.py')
MODE = os.path.join(HERE, 'run_abaqus_mode.py')

# 모델 정의에 해당하는 줄들. HF 쪽과 좌굴 쪽에서 '문자 그대로' 같아야 한다.
# (제외: HF/LF 전용 하중·감쇠·스텝·추출 관련. 그것들은 좌굴 모델에 없어야 정상이다.)
SHARED = [
    # ---- 형상 ----
    "BASE = 20.0",
    "HEIGHT = 10.0",
    "THICKNESS = 5.0e-6",
    "V1 = (BASE/2.0, HEIGHT, 0.0)",
    "V2 = (BASE, 0.0, 0.0)",
    "V3 = (0.0, 0.0, 0.0)",
    "def clamp_coord_L(x): return (10-10*x, 10-10*x, 0)",
    "def clamp_coord_R(x): return (10+10*x, 10-10*x, 0)",
    # ---- 재질/단면 ----
    "mat.Density(table=((1420.0,),))",
    "mat.Elastic(table=((2.5e9, 0.34),))",
    "HomogeneousShellSection(name='Section-Membrane', material='Kapton', thickness=THICKNESS)",
    # ---- 파트/스케치 ----
    "s.Line(point1=V3[:2], point2=V2[:2])",
    "s.Line(point1=V2[:2], point2=V1[:2])",
    "s.Line(point1=V1[:2], point2=V3[:2])",
    "p.BaseShell(sketch=s)",
    # ---- 메쉬 (모드 노드 일치의 핵심) ----
    #   균일 격자(aba_grid_mesh.fill_part, orphan mesh) 전환은 되돌렸다:
    #   'Part' 객체에 addNodes/addElements 가 없다(그 둘은 odb.Part 전용). 격자는
    #   .inp import 경로로 다시 시도한다. 지금은 자유 메쉬가 기준선이다.
    "p.seedPart(size=BASE/SEED_DIV, deviationFactor=0.1)",
    "p.setMeshControls(regions=p.faces, elemShape=QUAD_DOMINATED, technique=FREE, algorithm=MEDIAL_AXIS)",
    "elemTypeQuad = ElemType(elemCode=ELEM_CODE_QUAD, elemLibrary=STANDARD)",
    "elemTypeTri = ElemType(elemCode=ELEM_CODE_TRI, elemLibrary=STANDARD)",
    "ELEM_CODE_QUAD = S4",
    "ELEM_CODE_TRI = S3",
    "SEED_DIV = 200.0",
    "p.setElementType(regions=(p.faces,), elemTypes=(elemTypeQuad, elemTypeTri))",
    # ---- 클램프/정점 강체 패치 ----
    "radius=0.2",
    "u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON",
    "u3=SET",
    # ---- 프리텐션 정의 (좌굴의 alpha 가 곱해지는 기준값) ----
    # [2026-09-28] DISP_GLOBAL 은 여기서 검사하지 않는다: 운용점은 **의도적으로 분기**한 값이다
    #   (HF 165 um = 논문 정합 / 모드 소스 1.8e-5 = 클램프-프리 안정 영역). 대신
    #   check_operating_point() 가 두 값을 직접 읽어 배율을 로그로 남긴다.
    "angle_deg = 28.6",
    # ---- 스텝 GlobalTension (좌굴의 base state 를 만드는 스텝) ----
    "initialInc=0.0001, minInc=1e-8, maxNumInc=1000",
    # ---- 초기응력 ----
    "sigma11=SIGMA0, sigma22=SIGMA0, sigma33=0.0,",
]

# 모드 소스(run_abaqus_mode.py)는 클램프가 없는 모델이므로 clamp_coord_* 정의만 제외한다.
# 나머지는 전부 같아야 한다: *IMPERFECTION 은 노드 라벨로 매핑되기 때문.
SHARED_MODE = [s for s in SHARED if 'clamp_coord' not in s]

# 모드 소스와 HF 가 '선언하면' 달라져도 되는 항목 = base state/스텝 물리 파라미터.
#   이유: 모드 소스는 좌굴모드를 얻기 위해 base state 를 안정 영역에 둘 필요가 있다(사용자 튜닝).
#   메쉬/형상/요소/재료(SHARED 의 앞부분)는 절대 달라지면 안 된다 — 노드 라벨 매핑이 깨진다.
#   선언은 모드 소스 안의 `# CHECKER-DIVERGENCE: <이름,...>` 주석으로만 인정한다(논문 공개 의무).
# ============================================================================
# 좌굴 스크립트 <-> 모드 소스 (계약, 2026-10-06)
# ============================================================================
#   왜: buckle 이 실패한 원인을 '클램프'라고 단정할 수 없었다. 클램프 외에
#   DISP/SIGMA0/PERTURBATION/안정화/프리텐션모드/케이블 이 전부 달랐기 때문이다.
#   그래서 buckle 을 mode.py(ClampFree_Buckle)에 맞추고 클램프만 차이로 남겼다.
#   => 이 목록이 그 계약이다. 여기가 깨지면 '클램프만 다르다'는 전제가 무너진다.
SHARED_BUCKLE_MODE = [
    "DISP_GLOBAL = 1.8e-5",
    "SIGMA0 = 700.0",
    "PERTURBATION = 0.01",
    "MODE_STABILIZATION = 0.0005",
    "PRETENSION_MODE = 'corner2'",
    "DISP_TOP_OVER_CORNER = 1.4142135623",
    "PATTERN_SIGN = 1.0",
    "CABLE_RADIUS = 5.0e-4",
    "CABLE_AREA = np.pi * (CABLE_RADIUS**2)",
    "LEN_TOP = 0.280",
    "LEN_BOT = 0.689",
    #   [2026-10-07] 200/500 -> 10/20. 실측: vectors=500 을 요청해도 subspace 가 14 로
    #   줄어들어 200 개는 원리적으로 못 찾는다. 요청 수 <= 실제 subspace 차원이 조건이다.
    "N_EIG_BUCKLE = 4",
    "BUCKLE_VECTORS = 250",
]

DECLARABLE = [
    # [2026-09-28] DISP_GLOBAL 은 여기서 검사하지 않는다: 운용점은 **의도적으로 분기**한 값이다
    #   (HF 165 um = 논문 정합 / 모드 소스 1.8e-5 = 클램프-프리 안정 영역). 대신
    #   check_operating_point() 가 두 값을 직접 읽어 배율을 로그로 남긴다.
]


def declared_divergences(path):
    """`# CHECKER-DIVERGENCE:` 주석에 적힌 이름들을 돌려준다(정규화)."""
    with io.open(path, encoding='utf-8') as f:
        txt = f.read()
    hit = set()
    for ln in txt.splitlines():
        if 'CHECKER-DIVERGENCE:' in ln:
            tail = ln.split('CHECKER-DIVERGENCE:', 1)[1]
            for piece in tail.replace(';', ',').split(','):
                if norm(piece):
                    hit.add(norm(piece))
    return hit

# 좌굴 스크립트에 '있으면 안 되는' HF/LF 전용 것들 (역할 분리 검사)
FORBIDDEN_IN_BUCKLE = [
    # [2026-10-07] "Step-ClampTension" 은 여기 그대로 둔다.
    #   버클모드 추출에서는 클램프를 **별도 텐션 스텝으로 분리하지 않고** 글로벌 텐션에
    #   포함시키기로 했다(사용자 결정). 클램프를 Buckle(perturbation)에 걸면 LIVE 로만
    #   작용해 base state 에 압축을 만들지 못하므로, GlobalTension 에 둔다.
    "Step-ClampTension",
    "Step-Postbuckle",
    "Step-HighTension",
    "Step-Trigger",
    # 주의: "eval_abaqus.py" 는 금지어가 아니다. _resolve_here() 가 코드 디렉터리를
    #       찾기 위해 그 파일명을 정당하게 참조한다. 대신 '추출 호출' 자체를 막는다.
    "Calling extraction script",
    "BC_Trig_P",
    "TRIG_ON",
    "extraction.txt",
    "HF_ODB",
]

# 모드 스크립트에 '있으면 안 되는' 것들 (역할 분리: 모드 추출만 한다)
FORBIDDEN_IN_MODE = [
    "*IMPERFECTION",       # 모드를 '소비'하면 안 된다 (소비는 HF 담당)
    "Step-Postbuckle",
    "Step-HighTension",
    "Step-ClampTension",
    "HF_Postbuckle",
    "Step-Trigger",
    "BC_Trig_P",
    "TRIG_ON",
    "Calling extraction script",
]


def code_only(path):
    """주석·docstring 을 제거한 '살아있는 코드'를 공백 없는 정규형으로 돌려준다.

    구현 주의(2026-09-22 실제 오류): 처음에는 ast.unparse 를 썼는데 unparse 가 숫자
    리터럴을 재작성한다(5.0e-6 -> 5e-06, 2.5e9 -> 2500000000.0). 그래서 정상 항목이
    false FAIL 로 나왔다. 원문 리터럴을 보존해야 하므로 tokenize 로 토큰을 모은다.
    """
    import tokenize
    with io.open(path, encoding='utf-8') as f:
        src = f.read()

    # docstring 의 줄 범위만 수집 (그 STRING 토큰만 버린다)
    tree = ast.parse(src)
    doc_ranges = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.ClassDef, ast.AsyncFunctionDef)):
            body = getattr(node, 'body', [])
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                v = body[0].value
                doc_ranges.add((v.lineno, v.end_lineno))

    skip = (tokenize.COMMENT, tokenize.NL, tokenize.NEWLINE,
            tokenize.INDENT, tokenize.DEDENT, tokenize.ENCODING)
    parts = []
    for tok in tokenize.generate_tokens(io.StringIO(src).readline):
        if tok.type in skip:
            continue
        if tok.type == tokenize.STRING and (tok.start[0], tok.end[0]) in doc_ranges:
            continue
        parts.append(tok.string)
    return re.sub(r'\s+', '', ' '.join(parts))


def norm(s):
    return re.sub(r'\s+', '', s)


def main():
    if not (os.path.exists(NEW) and os.path.exists(BUCKLE) and os.path.exists(MODE)):
        print("!!! 파일을 찾을 수 없습니다: %s / %s / %s" % (NEW, BUCKLE, MODE))
        return 1
    _legacy = os.path.join(HERE, 'run_abaqus_new.py')
    if os.path.exists(_legacy):
        print("  [참고] run_abaqus_new.py 는 legacy(1차 요소 S3/S4) — 검사 대상이 아니다.")
    a = code_only(NEW)
    b = code_only(BUCKLE)
    c = code_only(MODE)

    print("=" * 74)
    print("모델 정의 일치 검사 (AST 기반: 주석/docstring 제외)")
    print("  %s   (live HF)" % os.path.basename(NEW))
    print("  %s" % os.path.basename(BUCKLE))
    print("  %s   (모드 소스)" % os.path.basename(MODE))
    print("=" * 74)
    bad = []
    for k in SHARED:
        nk = norm(k)
        ia, ib = nk in a, nk in b
        if not (ia and ib):
            bad.append(k)
        print("  %-4s HF=%-6s buckle=%-6s  %s"
              % ('OK' if (ia and ib) else '!!!', ia, ib, k[:60]))

    print()
    print("--- 모드 소스(run_abaqus_mode.py) 와 HF 의 모델 정의 일치 ---")
    print("    *IMPERFECTION 은 노드 라벨로 매핑되므로 메쉬/형상/프리텐션 정의가 같아야 한다.")
    print("    (제외: clamp_coord_* — 이 모델에는 클램프가 없다)")
    declared = declared_divergences(MODE)
    if declared:
        print("    [선언된 분기] %s" % ', '.join(sorted(declared)))
    bad4 = []
    for k in SHARED_MODE:
        nk = norm(k)
        ia, ic = nk in a, nk in c
        decl = ia and (not ic) and any(dn in nk for dn in declared)
        if not (ia and ic) and not decl:
            bad4.append(k)
        print("  %-5s HF=%-5s mode=%-5s  %s"
              % ('OK' if (ia and ic) else ('DIVRG' if decl else '!!!'), ia, ic, k[:60]))
    if declared and not bad4:
        print("    -> 프리텐션/base state 분기 %d건은 모드 소스에 선언됨(메쉬 정의는 전부 일치)."
              % len([k for k in SHARED_MODE
                     if norm(k) in a and norm(k) not in c
                     and any(dn in norm(k) for dn in declared)]))

    print()
    print()
    print("--- 좌굴 스크립트(run_abaqus_buckle.py) <-> 모드 소스(run_abaqus_mode.py) ---")
    print("    클램프만 다른 상태로 비교하려면 이 값들이 같아야 한다(2026-10-06 계약).")
    bad6 = []
    for k in SHARED_BUCKLE_MODE:
        nk = norm(k)
        ib, ic = nk in b, nk in c
        if not (ib and ic):
            bad6.append(k)
        print("  %-5s buckle=%-6s mode=%-6s  %s"
              % ("OK" if (ib and ic) else "!!!", ib, ic, k[:58]))

    # ------------------------------------------------------------------
    # [2026-10-07] PATTERN_SIGN 은 '선언'만 보면 부족하다 — 실제로 곱해지는지 본다.
    #   mode.py 는 PATTERN_SIGN * PERTURBATION * ... 로 섭동을 주는데, buckle 은
    #   PATTERN_SIGN 을 선언만 하고 곱하지 않았다. 값이 1.0 일 때는 수치가 같아
    #   드러나지 않지만, -1.0 등으로 바꾸면 두 스크립트가 조용히 갈라진다.
    #   [2026-10-07 보강] 케이스마다 LIVE 크기가 다르다(paper_s1 은 논문 크기를 쓴다).
    #   리터럴 1종만 보면 **정당한 케이스 추가를 위반으로 잡는다**(실측: paper45 분기를 넣자
    #   이 규칙이 위반 1건을 냈다). 그래서 서명된 크기 목록을 허용하고, 대신
    #   **서명 없이 크기를 대입하는 줄**이 있으면 위반으로 센다(느슨해지는 게 아니라 촘촘해진다).
    _SIGNED_MAGS = ("PERTURBATION", "PAPER_S1_LIVE_M")
    _pat_use_b = sum(1 for s in _SIGNED_MAGS if norm("PATTERN_SIGN * " + s) in b)
    _pat_use_c = 1 if norm("PATTERN_SIGN * PERTURBATION") in c else 0
    with open(BUCKLE, encoding='utf-8') as _f:
        _rawb_sig = _f.read()
    _unsigned = []
    for _ln in _rawb_sig.splitlines():
        _code = _ln.split('#')[0]
        if '=' not in _code:
            continue
        # 규칙의 범위를 좁힌다: **크기를 그대로 별칭에 넣는 줄**(name = PERTURBATION)만 잡는다.
        #   `CLAMP_PERT = PERTURBATION * CLAMP_DC_EFF` 처럼 '크기 x 계수' 로 중간값을 만드는 줄은
        #   정당하다(PATTERN_SIGN 은 그 중간값을 쓰는 자리에서 곱해진다) — 실측으로 이 오탐 2건을 봤다.
        _rhs = _code.split('=', 1)[1].strip()
        if _rhs in _SIGNED_MAGS and 'PATTERN_SIGN' not in _code:
            _unsigned.append(_code.strip())
    print()
    print("--- PATTERN_SIGN 이 실제로 섭동에 곱해지는지 (선언만으로는 부족) ---")
    print("  %-4s %-22s %s" % ('OK' if _pat_use_b else '!!!', 'run_abaqus_buckle.py',
                               '%d/%d종 서명 사용 %s'
                               % (_pat_use_b, len(_SIGNED_MAGS), _SIGNED_MAGS)))
    print("  %-4s %-22s %s" % ('OK' if _pat_use_c else '!!!', 'run_abaqus_mode.py',
                               '%d/1종 서명 사용 (PATTERN_SIGN * PERTURBATION)' % _pat_use_c))
    for _u in _unsigned:
        print("  !!!  서명 없는 LIVE 크기 대입: %s" % _u)
    if not _pat_use_b or not _pat_use_c:
        bad6.append("PATTERN_SIGN 사용(PATTERN_SIGN * PERTURBATION)")
        print("    -> 한쪽만 곱하면 PATTERN_SIGN != 1.0 에서 두 모델이 갈라진다.")
    if _unsigned:
        bad6.append("서명 없는 LIVE 크기 대입 %d줄" % len(_unsigned))

    # ------------------------------------------------------------------
    # [2026-10-07] 케이스 계약 — '클램프만 다르다'는 정합은 케이블 라우트 + 클램프 구속만인
    #   CASE='seed' / 'control_none' 에서만 성립한다. 'clamp_lf' 는 케이블을 만들지 않는
    #   별도 하중 경로(전부 집중하중)이고, 'paper_s1' 은 논문 baseline 재현용이다.
    #   둘 다 실패가 아니라 '비교 대상이 아님'이므로 경고만 찍고 위반으로 세지 않는다.
    #   b 는 code_only() 로 공백이 정규화된 텍스트라 ^ 앵커가 안 맞는다 — 원본을 읽는다.
    with open(BUCKLE, encoding='utf-8') as _f:
        _raw_b = _f.read()

    def _bconst(pat, default='(없음)'):
        _m = re.search(pat, _raw_b, re.M)
        return _m.group(1) if _m else default

    bad7 = []
    _case = _bconst(r"^CASE\s*=\s*'([^']+)'")
    _cang = _bconst(r"^CORNER_ANGLE_DEG\s*=\s*([0-9.eE+-]+)")
    _dfr = _bconst(r"^DEAD_FRAC\s*=\s*([0-9.eE+-]+)")
    _cf0 = _bconst(r"^CORNER_F0\s*=\s*([0-9.eE+-]+)")
    print()
    print("--- 하중 케이스 (buckle 단독 상수) ---")
    print("  run_abaqus_buckle.py  CASE = %s   CORNER_ANGLE_DEG = %s   DEAD_FRAC = %s   CORNER_F0 = %s"
          % (_case, _cang, _dfr, _cf0))
    _CABLE_CASES = ('seed', 'control_none')
    _KNOWN_CASES = _CABLE_CASES + ('clamp_lf', 'paper_s1')
    if _case == '(없음)':
        # [2026-10-07] 여기서 **실패**시킨다. 경고만 찍으면 CASE 줄을 주석 처리한 옛 사본이
        #   PASS 로 통과한다(규칙 18: 문자열 존재 검사는 주석 처리된 쌍둥이에 만족한다).
        bad7.append('CASE 를 읽지 못함(구버전 사본이거나 주석 처리됨)')
        print("  !!!  CASE 를 읽지 못했습니다 — buckle<->mode 계약의 적용 범위를 판정할 수 없다.")
    elif _case not in _KNOWN_CASES:
        bad7.append('알 수 없는 CASE=%s' % _case)
        print("  !!!  CASE='%s' 는 알려진 케이스가 아니다(%s)."
              % (_case, '|'.join(_KNOWN_CASES)))
    elif _case in _CABLE_CASES:
        print("  OK   CASE='%s' — 위 buckle<->mode 정합 계약이 그대로 성립한다." % _case)
    elif _case == 'clamp_lf':
        print("  WARN CASE='clamp_lf' — 케이블 0개(전부 집중하중). buckle<->mode 계약은 적용 대상이 아니다.")
        print("       판정은 모드 형상(면외 비율 / 압축영역 분포 / x_c 이동)으로 한다.")
        if _dfr == '(없음)':
            bad7.append('clamp_lf 인데 DEAD_FRAC 를 읽지 못함')
            print("  !!!  DEAD_FRAC 를 읽지 못했습니다 — clamp_lf 는 DEAD/LIVE 배분이 필수다.")
    else:
        print("  WARN CASE='paper_s1' — 논문 baseline 재현용. 정합 계약 대상 아님.")
        if _cang == '(없음)':
            bad7.append('paper_s1 인데 CORNER_ANGLE_DEG 를 읽지 못함')

    print("--- 좌굴 스크립트에 HF/LF 전용이 섞여 있지 않은지 (0 이어야 정상) ---")
    bad2 = []
    for k in FORBIDDEN_IN_BUCKLE:
        n = b.count(norm(k))
        if n:
            bad2.append((k, n))
        print("  %-4s %-22s %d회" % ('OK' if n == 0 else '!!!', k, n))

    print()
    print("--- 모드 스크립트에 HF/LF 전용이 섞여 있지 않은지 (0 이어야 정상) ---")
    bad5 = []
    for k in FORBIDDEN_IN_MODE:
        n = c.count(norm(k))
        if n:
            bad5.append((k, n))
        print("  %-4s %-22s %d회" % ('OK' if n == 0 else '!!!', k, n))

    # --- EVOL 필드출력 (면적가중 산출에 필요. 2026-09-22 누락을 실측으로 발견) ---
    print()
    print("--- 필드출력 EVOL (base_state_probe 면적가중) ---")
    for _label, _txt in (('run_abaqus_new.py', a), ('run_abaqus_buckle.py', b),
                         ('run_abaqus_mode.py', c)):
        _has = 'EVOL' in _txt
        print("  %s  %s" % ('OK  ' if _has else '!!! ', _label + ' EVOL=' + str(_has)))
        if not _has:
            bad.append('%s: EVOL 필드출력 누락 -> base_state_probe 면적가중 불가' % _label)

    print()
    print("--- 작업 디렉터리 분리 (산출물 혼입 방지) ---")
    def _run_dir_of(path):
        with io.open(path, encoding='utf-8') as f:
            txt = f.read()
        m = re.search(r'RUN_DIR_NAME\s*=\s*["\']([^"\']+)["\']', txt)
        return m.group(1) if m else None
    ra, rb = _run_dir_of(NEW), _run_dir_of(BUCKLE)
    print("  run_abaqus_new.py     RUN_DIR_NAME = %s" % ra)
    print("  run_abaqus_buckle.py  RUN_DIR_NAME = %s" % rb)
    bad3 = []
    if ra is None or rb is None:
        bad3.append('RUN_DIR_NAME 을 읽지 못함')
        print("  !!! RUN_DIR_NAME 을 찾지 못했습니다")
    elif ra == rb:
        bad3.append('same dir')
        print("  !!! 두 스크립트가 같은 디렉터리를 쓴다 -> 산출물이 섞이고")
        print("      *IMPERFECTION, FILE= 이 stale .fil 을 조용히 읽을 수 있다")
    else:
        print("  OK  분리됨 (buckle 산출물은 code/%s 에 쌓인다)" % rb)

    print()
    print("=" * 74)
    if bad or bad2 or bad3 or bad4 or bad5 or bad6 or bad7:
        print("RESULT: !!! 불일치 (HF↔buckle 위반 %d / buckle 역할위반 %d / 디렉터리 위반 %d"
              " / HF↔모드소스 위반 %d / 모드소스 역할위반 %d / buckle↔모드소스 위반 %d"
              " / CASE 계약 위반 %d)"
              % (len(bad), len(bad2), len(bad3), len(bad4), len(bad5), len(bad6), len(bad7)))
        print("        모드 노드가 어긋나면 *IMPERFECTION 이 조용히 실패한다.")
        print("        한쪽을 고쳤으면 다른 쪽도 같이 고쳐라.")
        return 1
    print("RESULT: PASS — HF↔buckle %d / HF↔모드소스 %d / buckle↔모드소스 %d 항목 일치, "
          "역할 분리 위반 0건" % (len(SHARED), len(SHARED_MODE), len(SHARED_BUCKLE_MODE)))
    print("=" * 74)
    return 0



# ---------------------------------------------------------------------------
# 운용점(DISP_GLOBAL) 일치 검사 (2026-09-28 신설)
#   왜 필요한가: HF 의 `DISP_GLOBAL = 0.000005 * PRETENSION_SCALE` 줄이 주석 처리된 뒤에도
#   위의 문자열 검사는 'HF=True' 로 통과했다 — 즉 운용점 변경을 아무 게이트도 잡지 못했다.
#   모드 소스와 HF 의 운용점이 어긋나면 모드 형상이 다른 프리스트레스 상태의 것이 되어
#   노드 섭동(임퍼펙션)이 물리적으로 틀어진다(2026-09-28 실측: 9.2배 어긋난 상태로 HF 를 돌렸다).
#   여기서는 주석을 제거한 살아있는 코드에서 값을 직접 뽑아 비교한다.
# ---------------------------------------------------------------------------
def _live_disp_global(path):
    live = code_only(path)
    hit = re.search(r'DISP_GLOBAL\s*=\s*([0-9.eE+-]+)', live)
    return float(hit.group(1)) if hit else None


def check_operating_point():
    """두 파일의 운용점을 '살아있는 코드'에서 뽑아 기록한다(FAIL 아님 — 선언된 분기).

    2026-09-28 정정: 처음에는 두 값이 같아야 통과로 만들었으나, 그렇게 맞추면
    모드 소스(run_abaqus_mode.py)의 base state 가 불안정해져 모드 추출이 실패한다
    (`THE EIGENVALUES CANNOT BE FOUND ... INSTABILITIES IN THE BASE STATE`, 고유값 전부 음수).
    모드 소스는 '좌굴모드를 뽑을 수 있는 안정 영역'(사용자 튜닝값)에 두는 것이 설계이며,
    운용점 차이는 DECLARABLE 에 선언된 분기다. 다만 그 비율은 매 검사마다 찍어 판독 가능하게 한다.
    """
    hf = _live_disp_global('run_abaqus.py')
    ms = _live_disp_global('run_abaqus_mode.py')
    print("--- 운용점(DISP_GLOBAL): 선언된 분기 (모드 소스는 안정 영역) ---")
    print("    HF = %s   /   mode = %s" % (hf, ms))
    if hf is None or ms is None:
        print("  FAIL  두 파일 중 DISP_GLOBAL 값을 읽지 못했습니다.")
        return False
    ratio = max(hf, ms) / max(min(hf, ms), 1e-30)
    print("  NOTE  HF/mode = %.2f배 (선언된 분기: 모드 소스는 안정 영역에 둔다)" % ratio)
    return True


if not check_operating_point():
    sys.exit(1)

SYM = "DISTRIBUTING|KINEMATIC|STRUCTURAL|UNIFORM_NDOF"   # CAE 가 받는 couplingType 심볼(문자열 아님)


def check_coupling_type():
    """패치 절점 결합 방식이 두 파일에서 같은지 (2026-09-28 신설, KINEMATIC->DISTRIBUTING).

    SHARED 목록에 넣지 않은 이유: legacy 스크립트(new/buckle)는 여전히 KINEMATIC 이라
    SHARED 에 두면 legacy 가 FAIL 한다. 여기서는 HF 와 모드 소스만 비교한다.
    """
    hf = re.search(r"couplingType=(" + SYM + r"),\s*weightingMethod", code_only('run_abaqus.py'))
    ms = re.search(r"COUPLING_TYPE=(" + SYM + r")", code_only('run_abaqus_mode.py'))
    print("--- 패치 결합 방식(Coupling) 일치: 모드 소스 vs HF ---")
    hv = hf.group(1) if hf else None
    mv = ms.group(1) if ms else None
    print("    HF = %s   /   mode = %s" % (hv, mv))
    if hv is None or mv is None:
        print("  FAIL  결합 방식을 읽지 못했습니다.")
        return False
    if hv != mv:
        print("  FAIL  결합 방식이 다릅니다 -> 모드 재추출이 필요합니다.")
        return False
    print("  OK    결합 방식 일치 (%s)" % hv)
    return True


if not check_coupling_type():
    sys.exit(1)

if __name__ == '__main__':
    sys.exit(main())
