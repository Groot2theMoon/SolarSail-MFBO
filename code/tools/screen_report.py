# -*- coding: utf-8 -*-
"""스크리닝 런 판정: .msg / .sta / 로그를 읽어 지표 표와 판정을 출력한다.

  python screen_report.py A_c1 A_c4 A_c8 B_c1 B_c4
  python screen_report.py                 (인자 없으면 기본 목록)

읽는 파일 (잡 이름 <job> 기준, 현재 디렉터리):
  <job>.msg   WALLCLOCK TIME / 증분 수 / 컷백 수 / 반복 수
  <job>.sta   마지막 데이터 줄의 8번째 필드 = STEP TIME/LPF (스텝이 0..1 이므로 곧 람다)
  <job>.log   실행 로그(stdout). 토큰 수와 스레드 수

지표:
  lambda        도달 람다
  lam/inc       증분당 람다  (수렴 제어 효율)
  lam/s         초당 람다    (실제 속도)
  speedup       첫 런(=기준)의 lam/s 대비 배수
"""
import io
import os
import re
import sys

DEFAULT_RUNS = ['A_c1', 'A_c4', 'A_c8', 'B_c1', 'B_c4']
BASELINE = 'A_c1'


def _read(path):
    if not os.path.exists(path):
        return None
    return io.open(path, 'r', encoding='utf-8', errors='replace').read()


def parse_run(name):
    """잡 이름 하나를 읽어 지표 dict 를 돌려준다. 없으면 None."""
    msg = _read(name + '.msg')
    if msg is None:
        return None
    rec = {'run': name, 'lambda': None, 'inc': None, 'cut': None,
           'iter': None, 'wall': None, 'tokens': None, 'threads': None,
           'completed': False, 'error': None, 'lambda_msg': None,
           'lambda_sta': None, 'mismatch': None}

    m = re.search(r'WALLCLOCK TIME \(SEC\)\s*=\s*([\d\.Ee+\-]+)', msg)
    if m:
        rec['wall'] = float(m.group(1))
    m = re.search(r'TOTAL OF\s+(\d+)\s+INCREMENTS', msg)
    if m:
        rec['inc'] = int(m.group(1))
    m = re.search(r'(\d+)\s+CUTBACKS IN AUTOMATIC INCREMENTATION', msg)
    if m:
        rec['cut'] = int(m.group(1))
    m = re.search(r'(\d+)\s+ITERATIONS INCLUDING', msg)
    if m:
        rec['iter'] = int(m.group(1))
    errs = re.findall(r'\*\*\*ERROR:[^\n]*', msg)
    if errs:
        rec['error'] = errs[0].strip()[:90]

    # .sta: STEP TIME/LPF 열. 열 위치가 스텝 종류에 따라 달라지므로 헤더에서 읽는다.
    #   (Riks 실측줄: '2 45 2 0 1 1 0.0503 5.077e-06 R' -> 0.0503 이 인덱스 6)
    #   헤더가 없으면 인덱스 6을 쓰고, 쓴 인덱스와 원본 줄을 항상 출력해 검증 가능하게 한다.
    sta = _read(name + '.sta')
    # 1차 출처: .msg 의 마지막 'FRACTION OF STEP COMPLETED' = 물리적 람다.
    #   (.sta 열 위치가 스텝 종류에 따라 달라지는 문제를 피한다)
    fr = re.findall(r'FRACTION OF STEP COMPLETED\s+([\d\.Ee+\-]+)', msg)
    if fr:
        try:
            rec['lambda_msg'] = float(fr[-1])
        except ValueError:
            pass
    rec['completed'] = not (('THE ANALYSIS HAS NOT BEEN COMPLETED' in msg)
                            or ('THE ANALYSIS HAS NOT BEEN COMPLETED' in (sta or '')))
    rec['lambda_col'] = None
    rec['last_line'] = None
    if sta:
        lines = sta.splitlines()
        col = None
        for l in lines[:80]:
            if 'STEP TIME' in l.upper():
                # 헤더는 2칸 이상 공백으로 묶음이 나뉜다. 첫 묶음(STEP INC ATT SEVERE EQUIL TOTAL)의
                # 열 수를 세고, 그 뒤 묶음에서 'STEP TIME' 을 찾아 실제 필드 인덱스를 만든다.
                groups = [g for g in re.split(r'\s{2,}', l.strip()) if g]
                if groups:
                    base = len(groups[0].split())
                    for i, g in enumerate(groups[1:]):
                        if g.upper().startswith('STEP TIME'):
                            col = base + i
                            break
                break
        data = [l for l in lines if re.match(r'^\s*\d+\s+\d+\s', l)]
        if data:
            use = col if col is not None else 6
            rec['lambda_col'] = use
            rec['last_line'] = data[-1].strip()
            f = [x for x in re.split(r'\s+', data[-1].strip()) if x]
            if len(f) > use:
                try:
                    rec['lambda_sta'] = float(f[use])
                except ValueError:
                    pass
    land = rec.get('lambda_msg')
    if land is None:
        land = rec.get('lambda_sta')
    rec['lambda'] = land
    # 교차검증: .sta 열 위치는 스텝 종류에 따라 다르다(TOTAL TIME/FREQ 열이 비면 인덱스가 당겨진다).
    #   그래서 열을 고정하지 않고, 마지막 줄의 숫자 중 .msg 값과 가장 가까운 것을 찾아 비교한다.
    a = rec.get('lambda_msg')
    if a is not None and rec.get('last_line'):
        nums = []
        for x in re.split(r'\s+', rec['last_line']):
            if re.match(r'^[+-]?\d*\.?\d+([eE][+-]?\d+)?$', x):
                try:
                    nums.append(float(x))
                except ValueError:
                    pass
        if nums:
            best = min(nums, key=lambda v: abs(v - a))
            rec['lambda_sta'] = best
            if abs(best - a) / max(abs(a), abs(best), 1e-30) > 0.02:
                rec['mismatch'] = 'msg=%.5g vs sta 최근접=%.5g' % (a, best)

    # 스레드 수: .log 가 없어도 .msg 의 솔버 배너에 나온다(실측 로그로 확인).
    m = re.search(r'(\d+)\s+THREADS? PER RANK', msg)
    if m:
        rec['threads'] = int(m.group(1))
    log = _read(name + '.log')
    if log:
        m = re.search(r'checked out (\d+) tokens', log)
        if m:
            rec['tokens'] = int(m.group(1))
        m = re.search(r'(\d+)\s+THREADS? PER RANK', log)
        if m:
            rec['threads'] = int(m.group(1))
    return rec


def _f(v, spec):
    return '-' if v is None else (spec % v)


def report(recs, baseline=None):
    baseline = baseline or (recs[0]['run'] if recs else BASELINE)
    base = next((r for r in recs if r['run'] == baseline), None)
    base_rate = None
    if base and base['wall'] and base['lambda']:
        base_rate = base['lambda'] / base['wall']

    width = max([len('run')] + [len(r['run']) for r in recs])
    hdr = ('%-*s %5s %9s %5s %5s %6s %8s %10s %9s %8s %6s'
           % (width, 'run', 'cpus', 'lambda', 'inc', 'cut', 'iter',
              'wall_s', 'lam/inc', 'lam/s', 'speedup', 'tokens'))
    print(hdr)
    print('-' * len(hdr))
    for r in recs:
        rate = (r['lambda'] / r['wall']) if (r['wall'] and r['lambda']) else None
        sp = (rate / base_rate) if (rate and base_rate) else None
        per_inc = (r['lambda'] / r['inc']) if (r['inc'] and r['lambda']) else None
        print('%-*s %5s %9s %5s %5s %6s %8s %10s %9s %8s %6s'
              % (width, r['run'],
                 _f(r['threads'], '%d'), _f(r['lambda'], '%.4g'),
                 _f(r['inc'], '%d'), _f(r['cut'], '%d'), _f(r['iter'], '%d'),
                 _f(r['wall'], '%.5g'), _f(per_inc, '%.3e'),
                 _f(rate, '%.3e'), _f(sp, '%.2fx'), _f(r['tokens'], '%d')))
    print()
    for r in recs:
        flags = []
        if not r['completed']:
            flags.append('미완주')
        if r['error']:
            flags.append(r['error'])
        if r.get('mismatch'):
            flags.append('열 불일치 ' + r['mismatch'])
        if r['inc'] and r['inc'] < 100 and r['lambda'] is not None and r['lambda'] < 0.999:
            flags.append('증분 100 미만에서 중단(스텝 조기 종료?)')
        if flags:
            print('  [%s] %s' % (r['run'], ' | '.join(flags)))

    print()
    print('판정')
    print('  기준 = %s : lam/s = %s' % (baseline, _f(base_rate, '%.3e')))
    for r in recs:
        if r['run'] == baseline or not r['wall'] or not r['lambda']:
            continue
        rate = r['lambda'] / r['wall']
        ratio = (rate / base_rate) if base_rate else None
        if ratio is None:
            continue
        verdict = '개선' if ratio > 1.15 else ('역효과' if ratio < 0.85 else '차이 없음')
        print('  %-*s vs %s : %s  (%.2f배)'
              % (width, r['run'], baseline, verdict, ratio))
    print()
    print('검증용 원본 (.sta 마지막 데이터줄 / 교차검증에 쓴 최근접값)')
    for r in recs:
        print('  %-*s sta=%s  %s' % (width, r['run'], _f(r.get('lambda_sta'), '%.5g'),
                                     (r.get('last_line') or '-')))
    print()
    print('  증분당 람다가 기준보다 작으면 수렴 제어가 증분을 못 키우는 것,')
    print('  크면 같은 증분 수로 더 멀리 간 것. lam/s 가 실제 속도 지표다.')


def main():
    names = [a for a in sys.argv[1:] if not a.startswith('-')] or DEFAULT_RUNS
    recs = [r for r in (parse_run(n) for n in names) if r]
    missing = [n for n in names if not os.path.exists(n + '.msg')]
    if missing:
        print('[알림] 다음 잡의 .msg 가 없어 제외했습니다: %s' % ', '.join(missing))
        print()
    if not recs:
        print('판정할 잡이 없습니다. 먼저 screen_runs.ps1 로 런을 돌리세요.')
        return 1
    report(recs, names[0] if names else None)
    return 0


if __name__ == '__main__':
    sys.exit(main())
