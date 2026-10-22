# -*- coding: utf-8 -*-
"""clamp 프리셋과 CLAMP_PULL 이 클램프 좌굴모드를 깨뜨리는지 격리 검증(OFAT).

배경
    1차 스윕 결과 λ×PERT 가 네 arm 에서 전부 2.0962e-07 m 로 동일했다.
    즉 PERTURBATION 은 스케일일 뿐이고 DISP_TOP / STAB 은 λ에 영향이 없었다.
    원인 후보는 CLAMP_PULL=0 하나뿐이므로 그 arm 을 추가한다.

설계 (CASE='clamp_only' 고정)

    arm 0 base   : corner2 / 1.4142 / 0.01 / 5e-4  pull=기본   <- 대조군
    arm 1 user   : corner2 / 2.0    / 0.1  / 2e-4  pull=기본   <- 요청 3종 동시
    arm 2 top    : corner2 / 2.0    / 0.01 / 5e-4  pull=기본   <- DISP_TOP 만
    arm 3 pert   : corner2 / 1.4142 / 0.1  / 5e-4  pull=기본   <- PERTURBATION 만
    arm 4 stab   : corner2 / 1.4142 / 0.01 / 2e-4  pull=기본   <- MODE_STABILIZATION 만
    arm 5 paper3 : paper3  / 1.4142 / 0.01 / 5e-4  pull=기본   <- PRETENSION_MODE 만
    arm 6 dc0    : corner2 / 1.4142 / 0.01 / 5e-4  pull=0      <- CLAMP_PULL=0 (진범 후보)

사용
    python sweep_clamp_preset.py --dry-run
    python sweep_clamp_preset.py                 # 7 arm
    python sweep_clamp_preset.py --arms 0,5,6    # 실패 arm + dc0 재확인
    python sweep_clamp_preset.py --no-retry      # 자동 재시도 끄기

판정
    PASS : 고유값 >= 2, |λ|min >= 1e-6, SYSTEM MATRIX 음수 고유값 < 100
    WARN : 모드 1개, 또는 임계 근접(1e-7 <= |λ|min < 1e-6)
    FAIL : ERROR, 고유값 0개, 또는 음수 고유값 >= 100 (base state 불안정)

주의
    - run_abaqus_buckle.py 는 **바이트 단위로 원복**한다(finally).
    - 고유값이 0개면 1회 자동 재시도한다(완료 마커 미기록 같은 일시 중단 대비).
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

# name, PRETENSION_MODE, DISP_TOP_OVER_CORNER, PERTURBATION, MODE_STABILIZATION, CLAMP_PULL
ARMS = [
    ('base',   'corner2', 1.4142135623, 0.01, 5.0e-4, None),
    ('user',   'corner2', 2.0,          0.1,  2.0e-4, None),
    ('top',    'corner2', 2.0,          0.01, 5.0e-4, None),
    ('pert',   'corner2', 1.4142135623, 0.1,  5.0e-4, None),
    ('stab',   'corner2', 1.4142135623, 0.01, 2.0e-4, None),
    ('paper3', 'paper3',  1.4142135623, 0.01, 5.0e-4, None),
    ('dc0',    'corner2', 1.4142135623, 0.01, 5.0e-4, 0.0),
]

NEG_FATAL = 100   # 이 이상이면 base state 불안정(실측: 실패 1870 / 정상 3)


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
            if not grab or not s:
                continue
            m = re.match(r'^(\d+)\s+([-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?)$', s)
            if m:
                rows.append(float(m.group(2).replace('D', 'E').replace('d', 'e')))
            else:
                grab = False
    return rows


def parse_msg(msg_path):
    """SYSTEM MATRIX 음수 고유값의 최대치와 ERROR 줄 수(합산하지 않는다)."""
    neg, err = 0, 0
    if not msg_path or not os.path.exists(msg_path):
        return neg, err
    with io.open(msg_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            if 'NEGATIVE EIGENVALUES' in line:
                m = re.search(r'HAS\s+(\d+)\s+NEGATIVE', line)
                if m:
                    neg = max(neg, int(m.group(1)))
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
    if neg >= NEG_FATAL:
        return 'FAIL', 'base 불안정(음수 고유값 %d)' % neg
    if n < 2:
        return 'WARN', '모드 %d개' % n
    if amin < 1e-7:
        return 'WARN', '준영모드 |λ|=%.2e' % amin
    if amin < 1e-6:
        return 'WARN', '임계 근접 |λ|=%.2e' % amin
    tail = (' (neg %d)' % neg) if neg else ''
    return 'PASS', 'λ %.4e~%.4e%s' % (min(eigs), max(eigs), tail)


def run_arm(i, a):
    """한 arm 실행. 반환: (eigs, neg, err, datp, secs)."""
    name, pret, top, pert, stab, pull0 = ARMS[i]
    pull = a.clamp_pull if pull0 is None else pull0
    txt = patch_constants(io.open(BUCKLE, 'r', encoding='utf-8').read(),
                          pret, top, pert, stab)
    with io.open(BUCKLE, 'w', encoding='utf-8') as f:
        f.write(txt)
    cmd = 'abaqus cae noGUI=run_abaqus_buckle.py -- %g %g %g %g' % (
        a.x_c, a.disp, pull, a.excl)
    print('  $ ' + cmd)

    tries = 2 if a.retry else 1
    eigs, neg, err, datp, secs = [], 0, 0, None, 0.0
    for k in range(tries):
        if k:
            print('  (고유값 0개 -> 1회 재시도)')
        t0 = time.time()
        subprocess.call(cmd, shell=True, cwd=HERE)
        secs += time.time() - t0
        datp = newest_dat(t0)
        msgp = (os.path.splitext(datp)[0] + '.msg') if datp else None
        eigs = parse_eigen(datp)
        neg, err = parse_msg(msgp)
        if eigs:
            break
    return eigs, neg, err, datp, secs


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--x_c', type=float, default=0.5)
    ap.add_argument('--disp', type=float, default=1e-4)
    ap.add_argument('--clamp-pull', type=float, default=5.56e-6, dest='clamp_pull',
                    help='기본 CLAMP_PULL [m] (arm 의 pull 이 None 일 때). 1e-4*0.0556')
    ap.add_argument('--excl', type=float, default=0.21)
    ap.add_argument('--arms', default='', help='예: 0,5,6 (생략하면 전부)')
    ap.add_argument('--no-retry', dest='retry', action='store_false', default=True,
                    help='고유값 0개일 때 자동 재시도 끄기')
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
    print('[sweep] x_c=%g  disp=%.4g  pull(기본)=%.4e  excl=%.2f  arms=%s'
          % (a.x_c, a.disp, a.clamp_pull, a.excl, idx))
    print('[sweep] 고정: CASE=clamp_only.  음수 고유값 임계=%d' % NEG_FATAL)
    results = []

    try:
        for i in idx:
            name, pret, top, pert, stab, pull0 = ARMS[i]
            pull = a.clamp_pull if pull0 is None else pull0
            print('\n===== arm %d/%d  %s : PRET=%s  DISP_TOP=%.10g  PERT=%.4g  STAB=%.4g  pull=%.3e ====='
                  % (i, len(ARMS) - 1, name, pret, top, pert, stab, pull))
            if a.dry_run:
                patch_constants(orig.decode('utf-8'), pret, top, pert, stab)
                results.append((name, pret, top, pert, stab, pull, 0, 'DRY', '', 0.0))
                continue
            eigs, neg, err, datp, secs = run_arm(i, a)
            if datp:
                shutil.copy2(datp, os.path.join(OUTDIR, 'arm%d_%s.dat' % (i, name)))
                msgp = os.path.splitext(datp)[0] + '.msg'
                if os.path.exists(msgp):
                    shutil.copy2(msgp, os.path.join(OUTDIR, 'arm%d_%s.msg' % (i, name)))
            print('  %.1f s  dat=%s' % (secs, os.path.basename(datp) if datp else '없음'))
            print('  λ = %s' % (['%.5e' % x for x in eigs[:6]] or '없음'))
            v, why = verdict(eigs, neg, err)
            pcrit = (min(eigs) * pert) if eigs else 0.0
            if pcrit:
                print('  P_crit = λ_min x PERT = %.4e m' % pcrit)
            print('  => %s (%s)' % (v, why))
            results.append((name, pret, top, pert, stab, pull, len(eigs), v, why, pcrit))
    finally:
        with open(BUCKLE, 'wb') as f:
            f.write(orig)
        print('\n[sweep] run_abaqus_buckle.py 원복 완료 (%d bytes)' % len(orig))

    print('\n===== 요약 =====')
    print('%-7s %-8s %-9s %-6s %-8s %-10s %3s  %-5s %-11s %s'
          % ('arm', 'PRET', 'DISP_TOP', 'PERT', 'STAB', 'PULL', 'n', '판정', 'P_crit[m]', '비고'))
    for (name, pret, top, pert, stab, pull, n, v, why, pcrit) in results:
        print('%-7s %-8s %-9.6g %-6.4g %-8.4g %-10.3e %3d  %-5s %-11s %s'
              % (name, pret, top, pert, stab, pull, n, v,
                 ('%.4e' % pcrit) if pcrit else '-', why))

    print('\n해석')
    print('  P_crit 열이 일정하면 PERTURBATION 은 스케일일 뿐 물리에 영향 없음')
    print('  arm6(dc0) FAIL + 음수 고유값 >= 100 -> CLAMP_PULL=0 이 진범')
    print('  arm6 PASS                            -> dc0 는 무해. 다른 변수를 봐야 함')
    print('  원본 .dat/.msg 는 %s 에 보존' % OUTDIR)
    return 0


if __name__ == '__main__':
    sys.exit(main())
