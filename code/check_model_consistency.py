"""check_model_consistency.py — 좌굴 모델과 HF 모델의 '모델 정의' 일치를 검사한다.

왜 필요한가
    run_abaqus_buckle.py 는 run_abaqus_new.py 에서 케이블만 제외한 **별개 스크립트**다.
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
NEW = os.path.join(HERE, 'run_abaqus_new.py')
BUCKLE = os.path.join(HERE, 'run_abaqus_buckle.py')

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
    "p.seedPart(size=BASE/200.0, deviationFactor=0.1)",
    "p.setMeshControls(regions=p.faces, elemShape=QUAD_DOMINATED, technique=FREE, algorithm=MEDIAL_AXIS)",
    "elemTypeQuad = ElemType(elemCode=S4, elemLibrary=STANDARD)",
    "elemTypeTri = ElemType(elemCode=S3, elemLibrary=STANDARD)",
    "p.setElementType(regions=(p.faces,), elemTypes=(elemTypeQuad, elemTypeTri))",
    # ---- 클램프/정점 강체 패치 ----
    "radius=0.2",
    "u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON",
    "influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC",
    "u3=SET",
    # ---- 프리텐션 정의 (좌굴의 alpha 가 곱해지는 기준값) ----
    "PRETENSION_SCALE = 10.0",
    "DISP_GLOBAL = 0.000005 * PRETENSION_SCALE",
    "angle_deg = 28.6",
    # ---- 스텝 GlobalTension (좌굴의 base state 를 만드는 스텝) ----
    "initialInc=0.0001, minInc=1e-8, maxNumInc=1000",
    "stabilizationMagnitude=0.0002,",
    # ---- 초기응력 ----
    "sigma11=SIGMA0, sigma22=SIGMA0, sigma33=0.0,",
    "SIGMA0 = 500.0",
]

# 좌굴 스크립트에 '있으면 안 되는' HF/LF 전용 것들 (역할 분리 검사)
FORBIDDEN_IN_BUCKLE = [
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
    if not (os.path.exists(NEW) and os.path.exists(BUCKLE)):
        print("!!! 파일을 찾을 수 없습니다: %s / %s" % (NEW, BUCKLE))
        return 1
    a = code_only(NEW)
    b = code_only(BUCKLE)

    print("=" * 74)
    print("모델 정의 일치 검사 (AST 기반: 주석/docstring 제외)")
    print("  %s" % os.path.basename(NEW))
    print("  %s" % os.path.basename(BUCKLE))
    print("=" * 74)
    bad = []
    for k in SHARED:
        nk = norm(k)
        ia, ib = nk in a, nk in b
        if not (ia and ib):
            bad.append(k)
        print("  %-4s new=%-5s buckle=%-5s  %s"
              % ('OK' if (ia and ib) else '!!!', ia, ib, k[:60]))

    print()
    print("--- 좌굴 스크립트에 HF/LF 전용이 섞여 있지 않은지 (0 이어야 정상) ---")
    bad2 = []
    for k in FORBIDDEN_IN_BUCKLE:
        n = b.count(norm(k))
        if n:
            bad2.append((k, n))
        print("  %-4s %-22s %d회" % ('OK' if n == 0 else '!!!', k, n))

    print()
    print("=" * 74)
    if bad or bad2:
        print("RESULT: !!! 불일치 %d건 (모델 정의) + %d건 (역할 분리)"
              % (len(bad), len(bad2)))
        print("        모드 노드가 어긋나면 *IMPERFECTION 이 조용히 실패한다.")
        print("        한쪽을 고쳤으면 다른 쪽도 같이 고쳐라.")
        return 1
    print("RESULT: PASS — 모델 정의 %d 항목 일치 / 좌굴 스크립트에 HF 전용 요소 없음"
          % len(SHARED))
    print("=" * 74)
    return 0


if __name__ == '__main__':
    sys.exit(main())
