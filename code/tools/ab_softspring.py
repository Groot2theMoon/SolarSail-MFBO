#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""RELEASE_EDGE_U3 A/B 실행기 — 3자유변 u3=0 vs. to-ground 스프링.

    python ab_softspring.py --dry-run          # 파일 편집/복원만 검증 (Abaqus 호출 없음)
    python ab_softspring.py                    # Arm B 실제 실행

왜 스크립트인가: 스위치를 켠 채로 두면 이후 모든 런이 오염된다. 이 실행기는
  git clean 확인 -> A 아카이브 -> 상수 1줄 편집 -> Abaqus 실행 -> **원본 바이트 복원** -> sha1 검증
순서를 강제하고, 실패해도 finally 에서 복원한다.

전제: seed 표는 이미 만들어져 있어야 한다(make_seed 를 부르지 않는다 = 모드 기저 고정).
"""
from __future__ import print_function

import argparse
import hashlib
import io
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # [P1] tools/ -> code/
SRC = os.path.join(HERE, 'run_abaqus.py')
PAT = re.compile(rb'(?m)^(RELEASE_EDGE_U3 = )(True|False)(\r?\n)')

# A/B 비교에서 사람이 먼저 볼 키 (있으면 출력, 없으면 생략)
KEYS = [
    'hf', 'dn_deg', 'dn_psi_deg',
    'centre_mean_pa', 'centre_minp_pa', 'area_frac_neg', 'minp_frac_neg',
    'uz_max_m', 'uz_min_m', 'uz_ave_m',
    'state_wrinkled_frac', 'zk_intensity_mean', 'zk_Iw',
]


def sha1(path):
    h = hashlib.sha1()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 16), b''):
            h.update(chunk)
    return h.hexdigest()


def read_flag(path):
    raw = open(path, 'rb').read()
    m = PAT.search(raw)
    if not m:
        raise SystemExit('[AB] !!! %s 에서 RELEASE_EDGE_U3 줄을 찾지 못했습니다' % path)
    return m.group(2).decode()


def set_flag(path, value):
    """바이트 단위 치환 — CRLF/LF 를 그대로 보존한다. 정확히 1곳만 바뀌어야 한다."""
    raw = open(path, 'rb').read()
    new, n = PAT.subn(rb'\g<1>' + value.encode() + rb'\g<3>', raw)
    if n != 1:
        raise SystemExit('[AB] !!! RELEASE_EDGE_U3 치환 %d곳 (1곳이어야 함)' % n)
    open(path, 'wb').write(new)
    return raw


def git_clean(path):
    try:
        out = subprocess.run(['git', 'status', '--porcelain', '--', path],
                             cwd=HERE, capture_output=True, text=True).stdout.strip()
    except Exception as e:                                    # noqa: BLE001
        print('[AB] git 확인 실패(%s) — 계속하지만 커밋 상태는 직접 확인하세요' % e)
        return True
    if out:
        print('[AB] !!! %s 에 커밋되지 않은 변경이 있습니다:\n%s' % (os.path.basename(path), out))
        return False
    return True


def read_kv(path):
    kv = {}
    if not os.path.exists(path):
        return kv
    for line in io.open(path, encoding='utf-8', errors='replace'):
        s = line.strip()
        if '=' in s and not s.startswith('#'):
            k, v = s.split('=', 1)
            kv[k.strip()] = v.strip()
    return kv


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def compare(a_path, b_path):
    A, B = read_kv(a_path), read_kv(b_path)
    print('\n[AB] ---- Arm A (u3=0) vs Arm B (spring) ----')
    print('%-22s %14s %14s %12s' % ('key', 'A', 'B', 'B/A-1'))
    rows = 0
    for k in KEYS:
        if k not in A and k not in B:
            continue
        va, vb = fnum(A.get(k)), fnum(B.get(k))
        rel = ('%+.1f%%' % (100.0 * (vb / va - 1.0))) if (va not in (None, 0) and vb is not None) else '-'
        print('%-22s %14s %14s %12s' % (k, A.get(k, '-'), B.get(k, '-'), rel))
        rows += 1
    if not rows:
        print('  (비교 가능한 키가 없습니다 — 파일 경로 확인)')
    print('[AB] A=%s\n[AB] B=%s' % (a_path, b_path))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--x_c', type=float, default=0.2)
    ap.add_argument('--d_c', type=float, default=1e-4)
    ap.add_argument('--config', default=os.path.join('configs', 'seed.json'))
    ap.add_argument('--dry-run', action='store_true', help='Abaqus 호출 없이 편집/복원만 검증')
    ap.add_argument('--no-archive', action='store_true', help='A 결과 txt 이동 생략')
    a = ap.parse_args(argv)

    def tag(x):
        return ('%g' % x).replace('.', 'p').replace('-', 'm')

    stem = 'xc%s_dc%s' % (tag(a.x_c), tag(a.d_c))
    live = os.path.join(HERE, 'results', 'results_%s.txt' % stem)
    arch = os.path.join(HERE, 'results', 'AB_A_u3fixed_results_%s.txt' % stem)
    log = os.path.join(HERE, 'results', 'AB_B_u3spring_%s.log' % stem)

    if read_flag(SRC) != 'False':
        raise SystemExit('[AB] 현재 RELEASE_EDGE_U3 가 False 가 아닙니다 — 먼저 되돌리세요.')
    if not git_clean(SRC):
        raise SystemExit('[AB] 중단: run_abaqus.py 가 dirty 합니다.')
    if not os.path.exists(live) and not a.no_archive:
        print('[AB] 경고: A 결과 파일이 없습니다 (%s)' % live)

    # ---- A 아카이브 (B 가 같은 파일에 append 하므로 먼저 옮긴다) ----
    if not a.dry_run and not a.no_archive and os.path.exists(live):
        if os.path.exists(arch):
            print('[AB] 기존 아카이브를 덮어씁니다: %s' % os.path.basename(arch))
        shutil.move(live, arch)
        print('[AB] A 아카이브: %s' % os.path.basename(arch))

    before = sha1(SRC)
    original = None
    rc = 1
    try:
        original = set_flag(SRC, 'True')
        print('[AB] RELEASE_EDGE_U3 -> True (sha1 %s)' % before)
        if a.dry_run:
            print('[AB] dry-run: Abaqus 호출 생략. 편집 확인 =', read_flag(SRC))
            rc = 0
        else:
            cmd = ['abaqus', 'cae', 'noGUI=run_abaqus.py', '--',
                   'HF', '%g' % a.x_c, '%g' % a.d_c, '--config', a.config]
            cmdline = subprocess.list2cmdline(cmd)   # abaqus 는 .bat -> shell=True 필요
            print('[AB] $', cmdline)
            with io.open(log, 'w', encoding='utf-8', errors='replace') as lf:
                p = subprocess.Popen(cmdline, shell=True, cwd=HERE, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, universal_newlines=True)
                assert p.stdout is not None
                for line in p.stdout:
                    lf.write(line)
                    if '[BC]' in line or line.startswith('RESULTS:') or '[DN]' in line:
                        print('   ', line.rstrip())
                rc = p.wait()
            print('[AB] 종료코드 %d / 로그 %s' % (rc, os.path.basename(log)))
    finally:
        # ---- 무조건 원본 바이트 복원 ----
        if original is not None:
            open(SRC, 'wb').write(original)
            after = sha1(SRC)
            ok = (after == before)
            print('[AB] 복원 %s (sha1 %s)' % ('OK' if ok else '!!! 실패', after))
            if not ok:
                print('[AB] !!! run_abaqus.py 를 git 에서 복구하세요: git checkout -- code/run_abaqus.py')

    if a.dry_run:
        return rc
    compare(arch, live)
    return rc


if __name__ == '__main__':
    sys.exit(main())
