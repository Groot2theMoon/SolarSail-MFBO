#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""hf_settings_stub.py — run_abaqus.py 세팅 리팩터(§19.110)의 **동작 불변** 검증 하네스.

무엇을 검증하나 (Abaqus/라이선스 불필요, 표준 라이브러리만)
    A. 상수 대조 — 리팩터 전/후 소스를 문장 단위로 exec 해 모듈 수준 상수를 {name: repr} 로 뽑고
       키/값을 비교한다. (Abaqus/numpy 는 스텁, mdb 조립부는 실패해도 건너뛴다.)
    B. 실행문 대조 — 상수 정의문을 제외한 **모든 실행문**을 AST 로 비교한다. 상수 참조는 그 값으로
       인라인하므로 `radius=PATCH_R_CLAMP` 와 `radius=0.21` 은 같은 트리로 접힌다.
       -> 154개 실행문(모델 조립/스텝/BC/제출 포함)이 전부 같아야 통과 = '덱을 만드는 호출열' 동일.
    C. 12조합 회귀 — 3구역(BC / Step-HighTension / Step-Postbuckle) x PULL_TOP(2) x NO_CLAMP(2)
       의 BC 호출 추적 + 콘솔 로그를 전/후 비교한다(pull_top_stub.py 의 구역 추출 재사용).

사용
    python3 hf_settings_stub.py                 # 기저 = HEAD~1 (직전 커밋의 소스)
    python3 hf_settings_stub.py --pre HEAD~3    # 다른 기저
    python3 hf_settings_stub.py --pre path/to/old_run_abaqus.py

종료코드: 0 OK / 1 FAIL / 2 사용법
한계: 라이선스가 없으므로 '실제 .inp 바이트'가 아니라 **덱을 만드는 호출열과 상수**를 비교한다.
"""
from __future__ import print_function

import ast
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NEW = os.path.join(HERE, 'run_abaqus.py')
REL = 'code/run_abaqus.py'

#   리팩터가 SETTINGS 로 옮긴(=본문에서 사라진) 상수 정의 이름.
MOVED = set('''RUN_DIR_NAME NUMCPUS BASE_PROBE EIG_RECORD sqrt2 N_EIG MODE_SOURCE
MODE_SOURCE_FIL MODE_SOURCE_DAT MODE_SOURCE_MSG MODE_SOURCE_STEP IMPERFECTION_NAME
IMPERFECTION_MODES IMPERFECTION_MAX_MODES IMPERFECTION_AMPL_T STABILIZE_LADDER
RUN_SELF_BUCKLE_JOB IMPERFECTION_MODE MODE_SOURCE_ODB MODE_SOURCE_STEP_NAME MODE_INSTANCE
SEED_FAMILY_NAME NO_CLAMP PULL_TOP TOP_PULL_SCALE N_EIG_BUCKLE BUCKLE_VECTORS MODEL_NAME
INSTANCE_NAME BASE HEIGHT ELEM_CODE_QUAD ELEM_CODE_TRI SEED_DIV THICKNESS TARGET_STRESS
CLAMP_EXCL_R CABLE_RADIUS CABLE_AREA LEN_TOP LEN_BOT CLAMP_CABLE_LEN V1 V2 V3
PRETENSION_SCALE DISP_GLOBAL PERTURBATION GLOBAL_FINAL angle_deg angle_rad cos_val sin_val
top_cos top_sin USE_GRID_MESH RUN_MODE SIGMA0 _EIGENSOLVER BUCKLE_ODB HF_ODB
CABLE_RHO CABLE_E CABLE_NU PATCH_R_CLAMP
RELEASE_EDGE_U3 EDGE_SPRING_ALPHA EDGE_SPRING_K'''.split())
NEW_NAMES = set(['CABLE_RHO', 'CABLE_E', 'CABLE_NU', 'PATCH_R_CLAMP',
                 #   [§19.111] 자유변 면외 구속 스위치 3개
                 'RELEASE_EDGE_U3', 'EDGE_SPRING_ALPHA', 'EDGE_SPRING_K'])
#   [§19.111] 새 기능이 만든 실행문을 식별하는 표지. 상수 재배치(§19.110) 감사와 기능 추가를
#     구분하기 위해, B 의 '전부 동일' 판정에서 **이 이름을 포함한 문**만 차이를 허용한다.
NEW_FEATURE = 'RELEASE_EDGE_U3'
D_C = 0.3    # 12조합 회귀에서 CLAMP_* 파생값 계산용(임의 고정값)

FAIL = []
N = [0]


def check(name, cond, detail=None):
    N[0] += 1
    ok = bool(cond)
    print("  %-4s %s%s" % ('OK' if ok else 'FAIL', name,
                           ('   <- ' + str(detail)) if (detail is not None and not ok) else ''))
    if not ok:
        FAIL.append(name)
    return ok


def read(path):
    with io.open(path, encoding='utf-8') as f:
        return f.read()


# ------------------------------------------------------------------ A. 상수 대조
STUBS = {
    'abaqus.py': ("class _Any(object):\n"
                  "    def __getattr__(self, k): return _Any()\n"
                  "    def __call__(self, *a, **k): return _Any()\n"
                  "    def __getitem__(self, k): return _Any()\n"
                  "    def __contains__(self, k): return False\n"
                  "    def __iter__(self): return iter(())\n"
                  "mdb = _Any()\n"),
    'abaqusConstants.py': ("S4 = 'AC.S4'\nS3 = 'AC.S3'\nON = 'AC.ON'\nOFF = 'AC.OFF'\n"
                           "def __getattr__(name):\n    return 'AC.' + name\n"),
    'step.py': "",
    'mesh.py': "class ElemType(object):\n    def __init__(self, **kw): pass\n",
    'regionToolset.py': "",
    'interaction.py': "",
    'numpy.py': ("import math\n"
                 "pi = math.pi\n"
                 "def deg2rad(x): return math.radians(x)\n"
                 "def cos(x): return math.cos(x)\n"
                 "def sin(x): return math.sin(x)\n"
                 "def array(*a, **k): return list(a[0]) if a else []\n"
                 "class linalg(object):\n"
                 "    @staticmethod\n"
                 "    def norm(v): return math.sqrt(sum(float(x) * float(x) for x in v))\n"
                 "def __getattr__(name):\n"
                 "    def _f(*a, **k): return 0.0\n"
                 "    return _f\n"),
}


def snapshot(src_text, workdir, stubs):
    """문장 단위 exec -> 모듈 수준 단순값 상수 dict + 콘솔 로그 줄."""
    path = os.path.join(workdir, 'run_abaqus.py')
    with io.open(path, 'w', encoding='utf-8') as f:
        f.write(src_text)
    for name in ('aba_imperfection.py', 'aba_grid_mesh.py'):
        shutil.copy(os.path.join(HERE, name), os.path.join(workdir, name))
    io.open(os.path.join(workdir, 'eval_abaqus.py'), 'w').write('# _resolve_here marker\n')

    tree = ast.parse(src_text, path)
    ns = {'__name__': '__main__', '__file__': path}
    cwd = os.getcwd()
    old_argv = list(sys.argv)
    buf = io.StringIO()
    sys.path.insert(0, stubs)
    sys.path.insert(0, workdir)
    sys.argv = [path, 'HF', '0.5', '0.3']
    os.chdir(workdir)
    try:
        old_out = sys.stdout
        sys.stdout = buf
        try:
            for st in tree.body:
                mod = ast.Module(body=[st], type_ignores=[])
                ast.fix_missing_locations(mod)
                try:
                    exec(compile(mod, path, 'exec'), ns)
                except SystemExit:
                    pass
                except Exception:
                    pass        # 모델 조립부(mdb) 등은 스텁에 없다 -> 건너뛴다
        finally:
            sys.stdout = old_out
    finally:
        os.chdir(cwd)
        sys.argv = old_argv
        sys.path.remove(workdir)
        sys.path.remove(stubs)

    simple = (bool, int, float, str, type(None), tuple, list, dict)
    consts = {}
    for k, v in ns.items():
        if k.startswith('__'):
            continue
        if isinstance(v, simple):
            consts[k] = repr(v)
    consts.pop('_RUN_T0', None)     # 런 시작 시각(monotonic) — 본질적으로 다르다
    return consts, buf.getvalue().splitlines()


# ------------------------------------------------------------------ B. 실행문 대조
class Inline(ast.NodeTransformer):
    def __init__(self, consts):
        self.consts = consts
        self.n = 0

    def visit_Name(self, node):
        if isinstance(node.ctx, ast.Load) and node.id in self.consts:
            try:
                val = ast.literal_eval(self.consts[node.id])
            except Exception:
                return node
            if isinstance(val, (bool, int, float, str)) or val is None:
                self.n += 1
                return ast.copy_location(ast.Constant(value=val), node)
        return node


def statements(src_text, path, consts):
    tree = ast.parse(src_text, path)
    out, dropped, n_inline = [], [], 0
    for st in tree.body:
        if (isinstance(st, ast.Assign) and len(st.targets) == 1
                and isinstance(st.targets[0], ast.Name) and st.targets[0].id in MOVED):
            dropped.append(st.targets[0].id)
            continue
        tr = Inline(consts)
        st = tr.visit(st)
        ast.fix_missing_locations(st)
        n_inline += tr.n
        out.append(ast.unparse(st))
    return out, dropped, n_inline


# ------------------------------------------------------------------ C. 12조합 회귀
def regions(src_text):
    sys.path.insert(0, HERE)
    import pull_top_stub as P
    return P, P.region_bc(src_text), P.region_drive(src_text, 0), P.region_drive(src_text, 1)


def run_region(P, region, src_text, pt, nc):
    ns, logs = P.make_ns(src_text, pt, nc)
    ns['CLAMP_PULL'] = ns['DISP_GLOBAL'] * D_C
    ns['CLAMP_FINAL'] = ns['GLOBAL_FINAL'] * D_C
    ns['CLAMP_PERT'] = ns['PERTURBATION'] * D_C
    exec(compile(region, '<region>', 'exec'), ns)
    return ns['my_model'].trace(), logs


def git(args):
    return subprocess.check_output(['git', '-C', ROOT] + args).decode('utf-8', 'replace')


def baseline(pre_ref):
    """기저 소스 (text, label). 기본: **직전 커밋** = '이번 변경이 방금 무엇을 바꿨나'.

    [§19.111] 종전 기본값은 §19.110 도입 커밋의 부모(리팩터 전 소스)였다. 그 비교는 상수 재배치
      감사에는 맞지만, 그 뒤에 기능 커밋이 쌓이면 '실행문 전부 동일' 이 **구조적으로** 깨진다
      (기능 추가는 실행문을 바꾸는 것이 정상이다 — 그걸 실패로 세면 하네스가 기능 개발을 막는다).
      상수 재배치 감사는 예전과 똑같이 `--pre <§19.110 부모>` 로 재현된다.
    작업트리가 깨끗하면 HEAD 는 내 소스와 같으므로 한 칸 더 내려가 `HEAD~1` 과 비교한다
      (커밋 전/후 모두 '실제로 달라진 것' 을 비교하게 만드는 규칙).
    """
    if pre_ref is not None:
        if os.path.exists(pre_ref):
            return read(pre_ref), pre_ref
        return git(['show', '%s:%s' % (pre_ref, REL)]), '%s:%s' % (pre_ref, REL)
    dirty = git(['status', '--porcelain', '--', REL]).strip() != ''
    ref = 'HEAD' if dirty else 'HEAD~1'
    return git(['show', '%s:%s' % (ref, REL)]), '%s (%s)' % (
        ref, '작업트리 변경 대비' if dirty else '깨끗한 트리 -> 직전 커밋 대비')


def main(argv):
    pre_ref = None
    if len(argv) > 1:
        if argv[1] == '--pre' and len(argv) > 2:
            pre_ref = argv[2]
        elif argv[1] in ('-h', '--help'):
            print(__doc__)
            return 0
        else:
            print('usage: python3 hf_settings_stub.py [--pre <git-ref|path>]')
            return 2

    try:
        pre_text, pre_label = baseline(pre_ref)
    except Exception as e:
        print('!!! 기저 소스를 읽지 못했습니다(%s): %s' % (pre_ref, e))
        return 2
    new_text = read(NEW)
    print('[harness] 현재 %s / 기저 %s' % (os.path.basename(NEW), pre_label))

    tmp = tempfile.mkdtemp(prefix='hf_settings_')
    stubs = os.path.join(tmp, 'stubs')
    os.makedirs(stubs)
    for name, body in STUBS.items():
        with io.open(os.path.join(stubs, name), 'w', encoding='utf-8') as f:
            f.write(body)
    try:
        # ---- A. 상수 대조 (같은 디렉터리 이름에서 연속 실행 -> _HERE 문자열 동일)
        print('\n[A] 상수 대조')
        work = os.path.join(tmp, 'work')
        os.makedirs(work)
        a, la = snapshot(pre_text, work, stubs)
        shutil.rmtree(work)
        os.makedirs(work)
        b, lb = snapshot(new_text, work, stubs)
        print('    pre %d개 / post %d개' % (len(a), len(b)))
        only_pre = sorted(set(a) - set(b))
        only_post = sorted(set(b) - set(a))
        diffs = [(k, a[k], b[k]) for k in sorted(set(a) & set(b)) if a[k] != b[k]]
        check('기존 상수 키 전부 존재 (only-pre 없음)', not only_pre, only_pre)
        check('신규 키는 선언된 목록 안 (%d개)' % len(only_post), set(only_post) <= NEW_NAMES, only_post)
        check('공통 키 값 diff 0', not diffs, diffs)
        check('모듈 수준 콘솔 로그 동일 (%d줄)' % len(la), la == lb,
              [(x, y) for x, y in zip(la + [''] * len(lb), lb + [''] * len(la)) if x != y][:1])
        for k in sorted(only_post):
            print('      + %-14s = %s' % (k, b[k]))

        # ---- B. 실행문 대조
        print('\n[B] 실행문(AST) 대조 — 상수 인라인 후')
        sa, da, ia = statements(pre_text, '<pre>', b)
        sb, db, ib = statements(new_text, '<new>', b)
        #   [§19.111] 상수 인라인을 **하지 않은** 문장도 뽑는다: 스위치 상수는 bool 이라
        #     인라인되면 `if not RELEASE_EDGE_U3:` 가 `if not False:` 로 접혀 표지가 사라진다.
        sa_raw, _, _ = statements(pre_text, '<pre>', {})
        sb_raw, _, _ = statements(new_text, '<new>', {})
        print('    pre %d문 / post %d문 (상수 정의 제외 %d/%d, 인라인 %d/%d회)'
              % (len(sa), len(sb), len(da), len(db), ia, ib))
        bad = [(i, x, y) for i, (x, y) in enumerate(zip(sa, sb))
               if x != y and NEW_FEATURE not in sa_raw[i] and NEW_FEATURE not in sb_raw[i]]
        n_new = sum(1 for x, y in zip(sa, sb) if x != y) - len(bad)
        print('    차이 %d문 (그중 신규 스위치 블록 %d문 -> 허용)' % (len(bad) + n_new, n_new))
        check('실행문 개수 동일', len(sa) == len(sb), '%d vs %d' % (len(sa), len(sb)))
        check('신규 스위치 블록 외 실행문 전부 동일 (%d문)' % (len(sa) - n_new), not bad,
              ('\n--- pre ---\n' + bad[0][1] + '\n--- post ---\n' + bad[0][2]) if bad else None)
        check('제외한 상수 정의 이름: pre-only 없음 / post 추가는 승격분만',
              not (set(da) - set(db)) and (set(db) - set(da)) <= NEW_NAMES,
              'pre-only %s / post-only %s' % (sorted(set(da) - set(db)), sorted(set(db) - set(da))))

        # ---- C. 12조합 회귀
        print('\n[C] 3구역 x PULL_TOP x NO_CLAMP = 12조합 (BC trace + 로그)')
        P, rp_bc, rp_fin, rp_post = regions(pre_text)
        _, rq_bc, rq_fin, rq_post = regions(new_text)
        check('구역 텍스트 동일 (BC/HighTension/Postbuckle)',
              (rp_bc, rp_fin, rp_post) == (rq_bc, rq_fin, rq_post))
        hashes = []
        for label, rp, rq in (('BC', rp_bc, rq_bc), ('HighTension', rp_fin, rq_fin),
                              ('Postbuckle', rp_post, rq_post)):
            if rp is None or rq is None:
                check('구역 추출: %s' % label, False)
                continue
            for pt in (False, True):
                for nc in (False, True):
                    tp, lp = run_region(P, rp, pre_text, pt, nc)
                    tq, lq = run_region(P, rq, new_text, pt, nc)
                    h = hashlib.sha256(tp.encode('utf-8')).hexdigest()
                    check('%s PULL_TOP=%-5s NO_CLAMP=%-5s trace+로그 동일 (%d자)'
                          % (label, pt, nc, len(tp)), tp == tq and lp == lq, h)
                    hashes.append(h)
        print('    12조합 trace 해시: %s' % hashlib.sha256(
            ''.join(hashes).encode('ascii')).hexdigest()[:32])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print('\nRESULT: %s (%d checks, %d fail)'
          % ('PASS' if not FAIL else 'FAIL', N[0], len(FAIL)))
    return 0 if not FAIL else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv))
