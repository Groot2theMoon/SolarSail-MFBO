# -*- coding: utf-8 -*-
"""clamp 프리셋 4개가 클램프 좌굴모드를 깨뜨리는지 격리 검증(OFAT).

배경
    run_abaqus_buckle.py 는 이제 프리셋 상수가 **조건식 없이 독립**이다
    (PAPER_MATCH 가 더 이상 이들을 제어하지 않는다). 그래서 하나씩 바꿔 원인을 분리할 수 있다.

설계 (CASE='clamp_only' 고정, 나머지 4개만 변화)

    arm 0 base   : corner2 / 1.4142 / 0.01 / 5e-4   <- 이미 성공한 조합(대조군)
    arm 1 user   : corner2 / 2.0    / 0.1  / 2e-4   <- 요청 3종 동시
    arm 2 top    : corner2 / 2.0    / 0.01 / 5e-4   <- DISP_TOP_OVER_CORNER 만
    arm 3 pert   : corner2 / 1.4142 / 0.1  / 5e-4   <- PERTURBATION 만
    arm 4 stab   : corner2 / 1.4142 / 0.01 / 2e-4   <- MODE_STABILIZATION 만
    arm 5 paper3 : paper3  / 1.4142 / 0.01 / 5e-4   <- PRETENSION_MODE 만(마지막 남은 paper 값)

사용
    python sweep_clamp_preset.py --dry-run      # 상수 이름 존재 확인(파일 무변경)
    python sweep_clamp_preset.py                # 6 arm 실행
    python sweep_clamp_preset.py --arms 0,5     # 일부만
    python sweep_clamp_preset.py --disp 1e-4 --clamp-pull 5.56e-6 --excl 0.21

주의
    - run_abaqus_buckle.py 는 **바이트 단위로 원복**한다(finally).
    - 각 arm 의 .dat/.msg 를 sweep_out/ 에 복사한다(잡 이름이 겹쳐도 보존).
"""
from __future__ import print_function

import argparse
import glob
import io
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
BUCKLE = os.path.join(HERE, 'run_abaqus_buckle.py')
OUTDIR = os.path.join(HERE, 'sweep_out')

ARMS = [
    # name,    PRETENSION_MODE, DISP_TOP_OVER_CORNER, PERTURBATION, MODE_STABILIZATION
    ('base',   'corner2', 1.4142135623, 0.01, 5.0e-4),
    ('user',   'corner2', 2.0,          0.1,  2.0e-4),
    ('top',    'corner2', 2.0,          0.01, 5.0e-4),
    ('pert',   'corner2', 1.4142135623, 0.1,  5.0e-4),
    ('stab',   'corner2', 1.4142135623, 0.01, 2.0e-4),
    ('paper3', 'paper3',  1.4142135623, 0.01, 5.0e-4),
]


def _sub1(pat, rep, s, what):
    new, n = re.subn(pat, rep, s, count=1, flags=re.M)
    if n != 1:
        raise RuntimeError('%s 줄을 1개 찾지 못했습니다(발견 %d)' % (what, n))
    return new


def patch_constants(txt, pret, top, pert, stab):
    """CASE 를 clamp_only 로 고정하고 프리셋 상수 4개만 덮어쓴다."""
    txt = _sub1(r"^CASE\s*=\s*'[^']*'", "CASE = 'clamp_only'", txt, 'CASE')
    txt = _sub1(r"^PRETENSION_MODE\s*=\s*'[^']*'",
                "PRETENSION_MODE = '%s'" % pret, txt, 'PRETENSION_MODE')
    txt = _sub1(r"^DISP_TOP_OVER_CORNER\s*=.*$",
                "DISP_TOP_OVER_CORNER = %.10g" % top, txt, 'DISP_TOP_OVER_CORNER')
    txt = _sub1(r"^PERTURBATION\s*=.*$",
                "PERTURBATION = %.10g" % pert, txt, 'PERTURBATION')
    txt = _sub1(r"^MODE_STABILIZATION\s*=.*$",
                "MODE_STABILIZATION = %.10g" % stab, txt, 'MODE_STABILIZATION')
    return txt


def newest_dat(t0):
    cands = [p for p in glob.glob(os.path.join(HERE, 'buckle', '*.dat'))
             if os.path.getmtime(p) >= t0 - 1.0]
    if not cands:
        return None
    return max(cands, key=os.path.getmtime)


def parse_eigen(dat_path):
    """'MODE NO  EIGENVALUE' 표를 읽는다. 파일 내 마지막 표를 쓴다."""
    if not dat_path or not os.path.exists(dat_path):
        return []
    rows, grab = [], False
    with io.open(dat_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            s = line.strip()
            if 'MODE NO' in s and 'EIGENVALUE' in s:
                rows, grab = [], True
                continue
            if not grab:
                continue
            if not s:
                continue
            m = re.match(r'^(\d+)\s+([-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?)$', s)
            if m:
                rows.append(float(m.group(2).replace('D', 'E').replace('d', 'e')))
            else:
                grab = False
    return rows


def parse_msg(msg_path):
    """음수 고유값 경고 합계와 ERROR 줄 수."""
    neg, err = 0, 0
    if not msg_path or not os.path.exists(msg_path):
        return neg, err
    with io.open(msg_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            if 'NEGATIVE EIGENVALUES' in line:
                m = re.search(r'HAS\s+(\d+)\s+NEGATIVE', line)
                if m:
                    neg += int(m.group(1))
            if line.strip().startswith('***ERROR'):
                err += 1
    return neg, err


def verdict(eigs, neg, err):
    if err:
        return 'FAIL', 'ERROR %d건' % err
    if not eigs:
        return 'FAIL', '고유값 0개'
    n = len(eigs)
    amin = min(abs(x) for x in eigs)
    if neg:
        return 'FAIL', '음수 고유값 경고 %d' % neg
    if n < 2:
        return 'WARN', '모드 %d개' % n
    if amin < 1e-7:
        return 'WARN', '준영모드 |λ|=%.2e' % amin
    if amin < 1e-6:
        return 'WARN', '임계 근접 |λ|=%.2e' % amin
    return 'PASS', 'λ %.4e ~ %.4e' % (min(eigs), max(eigs))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--x_c', type=float, default=0.5)
    ap.add_argument('--disp', type=float, default=1e-4)
    ap.add_argument('--clamp-pull', type=float, default=5.56e-6, dest='clamp_pull',
                    help='CLAMP_PULL [m]. disp*d_c -> 기본 1e-4*0.0556')
    ap.add_argument('--excl', type=float, default=0.21)
    ap.add_argument('--arms', default='', help='예: 0,1,5 (생략하면 전부)')
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args(argv)

    if not os.path.exists(BUCKLE):
        raise SystemExit('[sweep] run_abaqus_buckle.py 가 없습니다: %s' % BUCKLE)

    idx = [int(x) for x in a.arms.split(',') if x.strip()] if a.arms else list(range(len(ARMS)))
    for i in idx:
        if not (0 <= i < len(ARMS)):
            raise SystemExit('[sweep] arm 번호 범위 밖: %d' % i)

    if not a.dry_run:
        os.makedirs(OUTDIR, exist_ok=True)

    orig = io.open(BUCKLE, 'rb').read()
    print('[sweep] disp=%.4g  clamp_pull=%.4e  excl=%.2f  arms=%s'
          % (a.disp, a.clamp_pull, a.excl, idx))
    print('[sweep] 고정: CASE=clamp_only  (프리셋 상수 4개만 변화)')
    results = []

    try:
        for i in idx:
            name, pret, top, pert, stab = ARMS[i]
            print('\n===== arm %d/%d  %s : PRET=%s  DISP_TOP=%.10g  PERT=%.4g  STAB=%.4g ====='
                  % (i, len(ARMS) - 1, name, pret, top, pert, stab))
            txt = orig.decode('utf-8')
            txt = patch_constants(txt, pret, top, pert, stab)
            cmd = 'abaqus cae noGUI=run_abaqus_buckle.py -- %g %g %g %g' % (
                a.x_c, a.disp, a.clamp_pull, a.excl)
            print('  $ ' + cmd)

            eigs, neg, err, datp, msgp = [], 0, 0, None, None
            if not a.dry_run:
                with io.open(BUCKLE, 'w', encoding='utf-8') as f:
                    f.write(txt)
                t0 = time.time()
                rc = subprocess.call(cmd, shell=True, cwd=HERE)
                print('  rc=%s  (%.1f s)' % (rc, time.time() - t0))
                datp = newest_dat(t0)
                if datp:
                    shutil.copy2(datp, os.path.join(OUTDIR, 'arm%d_%s.dat' % (i, name)))
                    msgp = os.path.splitext(datp)[0] + '.msg'
                    if os.path.exists(msgp):
                        shutil.copy2(msgp, os.path.join(OUTDIR, 'arm%d_%s.msg' % (i, name)))
                eigs = parse_eigen(datp)
                neg, err = parse_msg(msgp)
                print('  dat=%s' % (os.path.basename(datp) if datp else '없음'))
                print('  λ = %s' % (['%.5e' % x for x in eigs[:6]] or '없음'))
                v, why = verdict(eigs, neg, err)
                print('  => %s (%s)' % (v, why))
            else:
                v, why = 'DRY', ''
            results.append((name, pret, top, pert, stab, len(eigs), v, why))
    finally:
        with open(BUCKLE, 'wb') as f:
            f.write(orig)
        print('\n[sweep] run_abaqus_buckle.py 원복 완료 (%d bytes)' % len(orig))

    print('\n===== 요약 =====')
    print('%-7s %-8s %-10s %-7s %-8s %4s  %-5s %s' %
          ('arm', 'PRET', 'DISP_TOP', 'PERT', 'STAB', 'n', '판정', '비고'))
    for (name, pret, top, pert, stab, n, v, why) in results:
        print('%-7s %-8s %-10.6g %-7.4g %-8.4g %4d  %-5s %s'
              % (name, pret, top, pert, stab, n, v, why))

    print('\n해석')
    print('  arm1 PASS                        -> 세 값 모두 무해. 원인은 PRETENSION_MODE(paper3)')
    print('  arm1 FAIL + arm2/3/4 중 하나 FAIL -> 그 상수가 원인')
    print('  arm1 FAIL + 2/3/4 모두 PASS       -> 세 값의 조합 효과(상호작용)')
    print('  arm5(paper3) FAIL                 -> paper3 단독으로 깨진다 = 진범')
    print('  원본 .dat/.msg 는 %s 에 보존' % OUTDIR)
    return 0


if __name__ == '__main__':
    sys.exit(main())
