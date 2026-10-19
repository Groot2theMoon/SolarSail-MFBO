#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""buckle_compare_stub.py - buckle_compare.py 를 **가짜 입력 파일**로 실제 실행하는 하네스.

왜 필요한가
    buckle_compare.py 의 MAC 핵심은 selftest 로 검증하지만, "파일을 읽어 - 표를 찍고 - 판정하고 -
    종료코드를 낸다" 는 배선은 CLI 를 실제로 돌려야 검증된다. Abaqus/라이선스는 필요 없다
    (가짜 .dat / 모드표 .txt 를 임시 디렉터리에 쓴다).

검증 항목
  (C1) 같은 런 두 번            -> max MAC=1.0000, 판정 '교체 불필요', exit 0
  (C2) 직교 모드 집합           -> max MAC=0.0000, 판정 '재생성 필요', exit 10
  (C3) 부호 반전 런             -> max MAC=1.0000, exit 0
  (C4) lambda 표/비율이 .dat 값과 일치
  (C5) 모드 순서가 뒤바뀐 두 런 -> 최대 쌍이 대각이 아님(행렬 인덱스)
  (C6) .dat 만(벡터 없음)       -> exit 4 + ODB/모드표 안내
  (C7) n_modes 제한이 실제로 먹는다(3요청 -> 2모드 표에서 2행)
  (C8) JSON 산출물 키/값
  (C9) MAC 행렬 모양(행=A, 열=B) - 비대칭 예제로 확인
  (C10) ODB 를 일반 python3 로 열면 exit 3 (odbAccess 없을 때만)
  (C11) selftest 자체가 exit 0

사용
    python3 buckle_compare_stub.py        # PASS/FAIL, 실패 시 EXIT=1
"""
from __future__ import print_function

import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
CMP = os.path.join(HERE, 'buckle_compare.py')

FAIL = []
N = [0]


def check(name, cond, detail=None):
    N[0] += 1
    ok = bool(cond)
    print('  %-4s %s%s' % ('OK' if ok else 'FAIL', name,
                           ('   <- ' + str(detail)) if (detail is not None and not ok) else ''))
    if not ok:
        FAIL.append(name)
    return ok


def run(args):
    p = subprocess.Popen([sys.executable, CMP] + list(args),
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out = p.communicate()[0].decode('utf-8', 'replace')
    return p.returncode, out


def write(path, text):
    with io.open(path, 'w', encoding='utf-8') as f:
        f.write(text)


def dat_text(lams):
    rows = '\n'.join('      %d      %.6e' % (i + 1, v) for i, v in enumerate(lams))
    return (' B U C K L I N G   F A C T O R   O U T P U T\n'
            ' MODE NO.   BUCKLING FACTOR\n' + rows + '\n')


def table_text(modes, lams=None, three_d=False):
    """aba_mode_from_odb.py 형식 모드표. modes = [{label: (u1,u2,u3)}]"""
    out = ['# source=fake', '# step=Step-Buckle', '# instance=MEMBRANE-1']
    for k, d in enumerate(modes, start=1):
        lam = lams[k - 1] if lams and k <= len(lams) else float(k)
        out.append('# MODE %d frame=%d lambda=%.6e raw_max_u3=1.0 sign_node=1'
                   % (k, k + 1, lam))
        out.append('MODE %d' % k)
        for lab in sorted(d):
            u1, u2, u3 = d[lab]
            if three_d:
                out.append('%d %.9e %.9e %.9e' % (lab, u1, u2, u3))
            else:
                out.append('%d %.9e' % (lab, u3))
    return '\n'.join(out) + '\n'


def mkvec(vals):
    return dict((i + 1, (0.0, 0.0, v)) for i, v in enumerate(vals))


def main():
    tmp = tempfile.mkdtemp(prefix='buckle_cmp_')
    try:
        #   노드 12개(임계 MIN_COMMON_NODES=10 을 넘긴다). u1/u2 는 u3 와 직교하도록 잡는다.
        u1 = [1.0, -2.0] + [0.0] * 10
        u2 = [0.5, 0.5, -1.0, 2.0, 1.5, -0.5, 1.0, -1.5, 0.25, -0.75, 0.5, 2.0]
        u3 = [0.0, 0.0, 1.0] + [0.0] * 9   # u1 과 직교(u1.u3=0)
        mA, mB = mkvec(u1), mkvec(u2)

        pa = os.path.join(tmp, 'A.txt')
        pb = os.path.join(tmp, 'B.txt')
        pc = os.path.join(tmp, 'C_signflip.txt')
        pd = os.path.join(tmp, 'D_orth.txt')
        pe = os.path.join(tmp, 'E_swap.txt')
        datA = os.path.join(tmp, 'A.dat')
        datB = os.path.join(tmp, 'B.dat')
        fA = os.path.join(tmp, 'flat.txt')
        fB = os.path.join(tmp, 'flat2.txt')

        lams_a = [8.196910e-03, 1.121200e-02]
        lams_b = [8.240000e-03, 1.130000e-02]
        write(pa, table_text([mA, mB], lams_a))
        write(pb, table_text([mA, mB], lams_b))
        write(datA, dat_text(lams_a))
        write(datB, dat_text(lams_b))
        write(fA, table_text([mA], lams_a[:1]))
        write(fB, table_text([mkvec(u3)], lams_a[:1]))          # 직교
        flip = dict((k, (0.0, 0.0, -v[2])) for k, v in mA.items())
        write(pc, table_text([flip, mB], lams_a))
        write(pd, table_text([mkvec(u3), mB], lams_a))          # 모드1 이 직교
        write(pe, table_text([mB, mA], lams_b))                 # 순서 뒤바뀜

        print('[harness] buckle_compare.py 가짜 입력 하네스 (%s)' % tmp)

        # ---------------------------------------------------------- (C1)
        rc, out = run([pa, pb, '2'])
        print('\n--- (C1) 같은 런 두 번 (exit=%d) ---\n%s' % (rc, out))
        check('C1 exit 0', rc == 0, rc)
        check('C1 max MAC=1.000000', 'max MAC = 1.000000' in out, out)
        check("C1 판정 'MAC>=0.9 -> 모드 소스 교체 불필요'",
              'MAC>=0.9 -> 모드 소스 교체 불필요' in out, out)
        check('C1 A1<->B1 쌍', '(A1, B1)' in out, out)

        # ---------------------------------------------------------- (C2)
        rc, out = run([fA, fB, '1'])
        print('\n--- (C2) 직교 (exit=%d) ---\n%s' % (rc, out))
        check('C2 exit 10 (재생성 필요)', rc == 10, rc)
        check('C2 max MAC=0.000000', 'max MAC = 0.000000' in out, out)
        check("C2 판정 'MAC<0.9 -> 재생성 필요'", 'MAC<0.9 -> 재생성 필요' in out, out)

        # ---------------------------------------------------------- (C3)
        rc, out = run([pc, pb, '2'])
        print('\n--- (C3) 부호 반전 (exit=%d) ---\n%s' % (rc, out))
        check('C3 exit 0', rc == 0, rc)
        check('C3 max MAC=1.000000 (부호 불변)', 'max MAC = 1.000000' in out, out)
        check('C3 쌍은 (A1, B1)', '(A1, B1)' in out, out)

        # ---------------------------------------------------------- (C4) + (C6)
        rc, out = run([datA, datB, '2'])
        print('\n--- (C4/C6) .dat 만 (exit=%d) ---\n%s' % (rc, out))
        check('C4 .dat 에서 lambda 2개 읽음',
              '8.196910e-03' in out and '1.121200e-02' in out, out)
        check('C4 비율 B/A = 1.0053', re.search(r'8\.240000e-03.*1\.0053', out) is not None, out)
        check('C4 벡터 없음 안내(exit 4)', rc == 4 and '모드 벡터를 읽지 못했다' in out, (rc, out))
        check('C6 안내에 aba_mode_from_odb 언급', 'aba_mode_from_odb.py' in out, out)

        # ---------------------------------------------------------- (C5)
        rc, out = run([pe, pb, '2'])
        print('\n--- (C5) 순서 뒤바뀜 (exit=%d) ---\n%s' % (rc, out))
        check('C5 최대 쌍이 대각이 아님 (A1,B2) 또는 (A2,B1)',
              ('(A1, B2)' in out or '(A2, B1)' in out), out)
        check('C5 exit 0', rc == 0, rc)

        # ---------------------------------------------------------- (C6)
        rc, out = run([datA, datB, '2'])
        check('C6 안내에 aba_mode_from_odb 언급', 'aba_mode_from_odb.py' in out, out)

        # ---------------------------------------------------------- (C7)
        rc, out = run([pa, pb, '9'])
        rows = len(re.findall(r'(?m)^  A\d+\s', out))
        check('C7 n_modes=9 요청 -> 있는 2모드만(행 2개)', rows == 2, (rows, out))

        # ---------------------------------------------------------- (C8)
        js = os.path.join(tmp, 'cmp.json')
        rc, _ = run([pa, pb, '2', '0.9', js])
        check('C8 JSON 파일 생성', os.path.exists(js), js)
        if os.path.exists(js):
            with io.open(js, encoding='utf-8') as f:
                d = json.load(f)
            print('         JSON: %s' % json.dumps(d, ensure_ascii=False))
            check('C8 JSON max_mac=1.0 / verdict=ok / pair=[1,1]',
                  abs(d['max_mac'] - 1.0) < 1e-12 and d['verdict'] == 'ok'
                  and d['pair'] == [1, 1], d)
            check('C8 JSON lambda_b 가 .dat 값', abs(d['lambda_b'][0] - 8.24e-3) < 1e-15, d)

        # ---------------------------------------------------------- (C9)
        #   A 는 1모드, B 는 2모드 -> 행 1개 x 열 2개. (A1,B1)=1, (A1,B2)<1
        rc, out = run([fA, pb, '2'])
        print('\n--- (C9) 비대칭 행렬 (exit=%d) ---\n%s' % (rc, out))
        check('C9 행 1개(A1) x 열 2개(B1,B2)',
              len(re.findall(r'(?m)^  A\d+\s', out)) == 1
              and out.count('B1') >= 1 and out.count('B2') >= 1, out)

        # ---------------------------------------------------------- (C10)
        try:
            import odbAccess  # noqa: F401
            has_odb = True
        except Exception:
            has_odb = False
        if not has_odb:
            fake = os.path.join(tmp, 'fake.odb')
            write(fake, 'not an odb')
            rc, out = run([fake, fake, '1'])
            print('\n--- (C10) .odb without odbAccess (exit=%d) ---\n%s' % (rc, out))
            check('C10 exit 3 + abaqus python 안내', rc == 3 and 'abaqus python' in out, (rc, out))
        else:
            print('\n--- (C10) SKIP: 이 인터프리터에 odbAccess 가 있다 ---')

        # ---------------------------------------------------------- (C11)
        rc, out = run(['selftest'])
        print('\n--- (C11) selftest (exit=%d) ---\n%s' % (rc, out))
        check('C11 selftest exit 0 + PASS', rc == 0 and 'RESULT: PASS' in out, (rc, out))

        print('\nRESULT: %s (%d checks, %d fail)'
              % ('PASS' if not FAIL else 'FAIL', N[0], len(FAIL)))
        if FAIL:
            for f in FAIL:
                print('  FAILED: %s' % f)
            return 1
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    sys.exit(main())
