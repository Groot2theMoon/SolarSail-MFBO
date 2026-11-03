"""`.sta`(Abaqus 상태 파일) 판독기 — "이 모델이 완주 가능한가"를 10~20분에 판정한다.

사용법:
    python sta_increment_report.py aba\HF_Postbuckle.sta            (Abaqus 불필요, 표준 파이썬)
    python sta_increment_report.py aba\HF_Postbuckle.sta 1.4 1.0    (증분당 초, 스텝 목표시간)

왜 필요한가 (2026-09-28 실측):
    Step-Postbuckle 은 **증분 크기 리밋사이클**에 빠져 있다. 성공 증분이 Abaqus 기본 배율 1.5 로
    자라다가 상한(약 1~2e-4)에서 **실제로 수렴 실패**하고, 정확히 1/4 로 컷백되어 다시 자란다.
    그 결과 평균 증분이 4.5e-5 수준으로 묶여 완주에 18,000+ 증분(7시간+)이 필요하다.
    -> 완주를 기다리지 말고 **초기 100~200 증분의 '증분 상한'** 만 읽어서
       모델 변형(요소 차수 / 메쉬 / 패치)이 증분 제어를 회복시켰는지 판정한다.
       상한이 3e-4 이상으로 회복되면 그 변형이 원인을 제거한 것이다.

읽는 값:
    - 완료 증분 수 / 시도 수 / 컷백 수·비율
    - 성공 증분 최대값(= 도달한 증분 상한), 평균, 최근 구간 평균
    - 실패(컷백 직전) 증분 크기 범위
    - 목표 스텝시간까지 남은 증분 및 예상 시간
"""
from __future__ import print_function

import sys

def parse(path):
    rows = []
    with open(path, 'r') as f:
        for line in f:
            p = line.split()
            if len(p) < 9:
                continue
            try:
                step = int(p[0]); inc = int(p[1])
                att = int(''.join(ch for ch in p[2] if ch.isdigit()))
                sev = int(p[3]); eq = int(p[4]); tot = int(p[5])
                ttot = float(p[6]); stime = float(p[7]); ds = float(p[8])
            except ValueError:
                continue
            rows.append(dict(step=step, inc=inc, att=att, sev=sev, eq=eq,
                             tot=tot, ttot=ttot, stime=stime, ds=ds))
    return rows

def report(rows, sec_per_inc, target):
    steps = []
    for r in rows:
        if r['step'] not in steps:
            steps.append(r['step'])
    for st in steps:
        rs = [r for r in rows if r['step'] == st]
        incs = {}
        for r in rs:
            incs.setdefault(r['inc'], []).append(r)
        done = len(incs)
        tries = len(rs)
        cut = tries - done
        ok, fail = [], []
        for inc in sorted(incs):
            seq = sorted(incs[inc], key=lambda r: r['att'])
            ok.append(seq[-1]['ds'])                  # 마지막 attempt = 성공
            for r in seq[:-1]:
                fail.append(r['ds'])                  # 앞선 attempt = 컷백
        tail = ok[-100:] if len(ok) >= 100 else ok
        avg_tail = sum(tail) / float(len(tail)) if tail else 0.0
        avg = sum(ok) / float(len(ok)) if ok else 0.0
        stime = rs[-1]['stime']
        print("  [스텝 %d] 완료 증분 %d / 시도 %d / 컷백 %d (%.0f%%)"
              % (st, done, tries, cut, (100.0 * cut / tries) if tries else 0.0))
        print("           성공 증분: 최대 %.4e / 평균 %.4e / 최근100 평균 %.4e"
              % (max(ok) if ok else 0.0, avg, avg_tail))
        if fail:
            print("           실패(컷백 직전) 증분: 최소 %.4e / 최대 %.4e  (%d회)"
                  % (min(fail), max(fail), len(fail)))
        ok_tot = [max(incs[i], key=lambda r: r['att'])['tot'] for i in incs]
        print("           마지막 step time %.4f / 성공 attempt 최대 iterations %d"
              % (stime, max(ok_tot) if ok_tot else 0))
        if stime < target and avg_tail > 0:
            left = (target - stime) / avg_tail
            print("           >>> 완주까지 약 %.0f 증분 남음, %.1f 시간 (%.2f s/증분 가정)"
                  % (left, left * sec_per_inc / 3600.0, sec_per_inc))
        # 판정
        ceil_ok = max(ok) if ok else 0.0
        print("           >>> 판정: 증분 상한 %.3e -> %s"
              % (ceil_ok, "회복(>=3e-4)" if ceil_ok >= 3e-4 else "여전히 묶임(<3e-4)"))

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    path = sys.argv[1]
    sec = float(sys.argv[2]) if len(sys.argv) > 2 else 1.4
    tgt = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    rows = parse(path)
    if not rows:
        print("파싱된 줄이 없습니다:", path)
        return 2
    print("=" * 74)
    print(".sta 증분 리포트: %s  (증분당 %.2f s 가정, 목표 step time %.2f)" % (path, sec, tgt))
    print("=" * 74)
    report(rows, sec, tgt)
    print()
    print("[읽는 법] 성공 증분 최대값이 곧 '이 모델이 넘을 수 있는 최대 걸음'이다.")
    print("          컷백 비율 20%대 + 성공 최대 <3e-4 = 리밋사이클 -> 완주에 수만 증분.")
    return 0

if __name__ == '__main__':
    sys.exit(main())
