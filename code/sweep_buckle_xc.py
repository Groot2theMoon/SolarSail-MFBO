"""
sweep_buckle_xc.py — 버클모드가 성공하는 x_c 범위를 찾는 스윕 (Abaqus 머신에서 실행).

무엇을 하나
    x_c 목록 각각에 대해 run_abaqus_buckle.py 를 실행한다. 실패하면 f(DEAD_FRAC)를
    사다리 [0.50, 0.35, 0.20, 0.10, 0.05] 로 낮춰 가며 '성공하는 최대 f' 를 찾는다.
    각 시도의 결과(성공/실패 + 실패 시그니처)를 buckle/sweep_xc_results.csv 에 쌓는다.

왜 f 를 낮추나
    f=0.5 는 x_c=0.5 에서만 개시(50.4%) 아래다. 개시 하중비는 x_c 마다 다르고 **비단조**다
    (0.35 -> 19.2% / 0.50 -> 50.4% / 0.65 -> 37.5%). 개시 위의 base state 에는 안정한
    선형모드가 없어 'THE EIGENVALUES CANNOT BE FOUND' 로 죽는다 — 창(N_EIG) 문제가 아니므로
    좌굴 창을 키우는 것으로는 못 넘긴다(음수 고유값 234개 > 요청 수).

무엇을 고정하나
    클램프 당김 = DISP x 0.0556 (실측 동작 운용값). d_c=0.5/1.0 실패는 인위 프로브였고,
    운용점과 무관하다 — 모드 추출의 목적은 'guidance' 이므로 운용 클램프 하나로 충분하다.

주의
    아래 DISP / EXCL_R / CLAMP_DC 는 run_abaqus_buckle.py 의 상수와 **값이 같아야 한다**.
    그쪽 상수를 바꾸면 여기도 같이 고친다(다르면 스윕이 다른 운용점을 시험하게 된다).

사용법 (Windows PowerShell, code 디렉터리에서)
    python sweep_buckle_xc.py --dry-run             # 계획만 출력 (Abaqus 불필요)
    python sweep_buckle_xc.py                       # 전체 스윕 (기본 x_c 0.20~0.80, ~10-30분)
    python sweep_buckle_xc.py --xc 0.30 0.35 0.40   # 일부만
    python sweep_buckle_xc.py --ladder 0.5 0.2      # 사다리 지정(내림차순)
    Ctrl+C 로 중단해도 안전하다. 다시 실행하면 CSV 에 이미 있는 (x_c, f) 는 건너뛴다.

출력
    buckle/sweep_xc_results.csv — x_c, f, verdict, negative, seconds, job, note
    판정: ok(**양수 lambda 존재**) / fail:negative-modes(표는 있으나 전부 음수 = 개시 넘음) /
          fail:eigenvalues(표 없음 — 창·조건화) / fail:other(그 외) / fail:no-job(실행 실패)
    note 열에 n과 lambda 값 기록. 개시 하중비 = f + (1-f) x lambda1 (x_c=0.5 검증: 0.5041).
"""
import argparse
import csv
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
BUCKLE = os.path.join(HERE, 'buckle')
OUT_CSV = os.path.join(BUCKLE, 'sweep_xc_results.csv')

DISP = 1.8e-5        # m — 운용점 (DISP_GLOBAL, mode.py 와 동일)
EXCL_R = 0.21        # m — CLAMP_EXCL_R (이중 구속 방지, 격자 경계 정합)
CLAMP_DC = 0.0556    # 운용 클램프 비율 (CLAMP_DC — 실측 동작값)
DEFAULT_LADDER = [0.50, 0.35, 0.20, 0.10, 0.05]   # 내림차순. 0.05 아래는 조건화 악화(0.018 퇴화)


def buckle_cmd(xc, f):
    """run_abaqus_buckle.py 1회 실행 명령. shell=True — Windows abaqus 는 .bat 이다."""
    return ('abaqus cae noGUI=run_abaqus_buckle.py -- %.4g %.4g %.6g %.4g %.4g'
            % (xc, DISP, DISP * CLAMP_DC, EXCL_R, f))


def read_text(path):
    try:
        with open(path, 'r', errors='replace') as fh:
            return fh.read()
    except OSError:
        return ''


def parse_lams(dat):
    """.dat 의 MODE NO/EIGENVALUE 표에서 lambda 열을 읽는다 (없으면 None)."""
    m = re.search(r'MODE NO\s+EIGENVALUE(.*?)(?:THE ANALYSIS|\Z)', dat, re.S)
    if not m:
        return None
    lams = [float(x) for x in re.findall(
        r'^\s*\d+\s+([-+]?[0-9.]+(?:[Ee][-+]?\d+)?)\s*$', m.group(1), re.M)]
    return lams or None


def classify(job, out):
    """(verdict, note) — .dat/.msg 기반. **양수 lambda 가 하나도 없으면 ok 가 아니다**:
        전부 음수면 base 가 이미 개시를 넘은 상태의 unloading 분기다(시드로 쓸 수 없다).
        개시 하중비 = f + (1-f) x lambda1 (x_c=0.5 검증: 0.5+0.5x8.1969e-3 = 0.5041)."""
    dat = read_text(os.path.join(BUCKLE, job + '.dat')) if job else ''
    msg = read_text(os.path.join(BUCKLE, job + '.msg')) if job else ''
    neg = ''
    m = re.search(r'SYSTEM MATRIX HAS (\d+) NEGATIVE', msg)
    if m:
        neg = m.group(1)
    if not job:
        return 'fail:no-job', neg
    if 'CANNOT BE FOUND' in dat or 'CANNOT BE FOUND' in msg:
        return 'fail:eigenvalues', neg
    lams = parse_lams(dat)
    if lams:
        note = 'n=%d lam=%s' % (len(lams), ','.join('%.4g' % v for v in lams[:4]))
        if any(v > 0 for v in lams):
            return 'ok', note
        return 'fail:negative-modes', note
    if 'HAS BEEN COMPLETED' in dat:
        return 'fail:other', neg          # 완주했는데 고유값 표가 없다 — f 로 못 넘긴다
    return 'fail:other', neg


def run_one(xc, f, timeout):
    cmd = buckle_cmd(xc, f)
    t0 = time.time()
    try:
        p = subprocess.run(cmd, cwd=HERE, shell=True, capture_output=True,
                           text=True, timeout=timeout)
        out = (p.stdout or '') + (p.stderr or '')
    except subprocess.TimeoutExpired:
        return 'fail:timeout', '', time.time() - t0
    m = re.search(r'Submitting Job:\s*(\S+)', out)
    job = m.group(1) if m else ''
    verdict, neg = classify(job, out)
    return verdict, (job, neg), time.time() - t0


def load_done():
    """이미 기록된 (x_c, f) 집합과 x_c 별 성공 여부."""
    done, ok_xc = set(), set()
    if os.path.exists(OUT_CSV):
        with open(OUT_CSV, 'r', newline='', encoding='utf-8') as fh:
            for row in csv.DictReader(fh):
                done.add((round(float(row['x_c']), 4), round(float(row['f']), 4)))
                if row['verdict'] == 'ok':
                    ok_xc.add(round(float(row['x_c']), 4))
    return done, ok_xc


def append_row(row):
    new = not os.path.exists(OUT_CSV)
    os.makedirs(BUCKLE, exist_ok=True)
    with open(OUT_CSV, 'a', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(['x_c', 'f', 'verdict', 'negative', 'seconds', 'job', 'note'])
        w.writerow(row)
    print('  CSV <- %s' % (', '.join(str(v) for v in row)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--xc', nargs='*', type=float,
                    default=[round(0.20 + 0.05 * i, 2) for i in range(13)],
                    help='x_c 목록 (기본 0.20~0.80, 0.05 간격)')
    ap.add_argument('--ladder', nargs='*', type=float, default=DEFAULT_LADDER,
                    help='DEAD_FRAC 사다리, 내림차순')
    ap.add_argument('--timeout', type=float, default=1800.0, help='회당 제한 [s]')
    ap.add_argument('--dry-run', action='store_true', help='계획만 출력')
    a = ap.parse_args()

    ladder = sorted(a.ladder, reverse=True)
    done, ok_xc = load_done()

    plan = []
    for xc in a.xc:
        if round(xc, 4) in ok_xc:
            print('x_c=%.2f — 이미 성공 기록, 건너뜀' % xc)
            continue
        for f in ladder:
            if (round(xc, 4), round(f, 4)) in done:
                continue
            plan.append((xc, f))

    print('=' * 72)
    print('스윕: x_c %d개 x 사다리 %s / 실행할 조합 %d개'
          % (len(a.xc), ladder, len(plan)))
    print('회당 30~60초 예상 -> 약 %.0f~%.0f분' % (len(plan) * 0.5, len(plan) * 1.0))
    print('=' * 72)
    for xc, f in plan:
        print('  x_c=%.2f f=%.2f  |  %s' % (xc, f, buckle_cmd(xc, f)))
    if a.dry_run:
        print('(dry-run — 실행하지 않았다)')
        return 0

    for i, (xc, f) in enumerate(plan, 1):
        print('\n[%d/%d] x_c=%.2f f=%.2f' % (i, len(plan), xc, f))
        verdict, info, dt = run_one(xc, f, a.timeout)
        if isinstance(info, tuple):
            job, neg = info
        else:
            job, neg = info, ''
        append_row(['%.2f' % xc, '%.2f' % f, verdict, neg, '%.1f' % dt, job, ''])
        print('  -> %s (%.0fs)%s' % (verdict, dt, ('  negative=%s' % neg) if neg else ''))
        if verdict == 'ok':
            continue
        if verdict in ('fail:eigenvalues', 'fail:negative-modes'):
            continue                       # 개시 위/넘음 — f 를 낮춰 재시도
        print('  !! %s — f 로 못 넘긴다. 이 x_c 의 나머지 사다리 생략.' % verdict)
        for f2 in ladder:
            if f2 < f and (round(xc, 4), round(f2, 4)) not in done:
                append_row(['%.2f' % xc, '%.2f' % f2, 'skipped', '', '0', '',
                            '앞 단계가 %s' % verdict])
        break
    print('\n완료. 결과: %s' % OUT_CSV)
    return 0


if __name__ == '__main__':
    sys.exit(main())
