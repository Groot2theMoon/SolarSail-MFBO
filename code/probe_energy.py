#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""안정화 에너지 점검 — stabilize 로 "완주시킨" 해석이 신뢰 가능한지 판정한다 (2026-10-13).

    abaqus python probe_energy.py <job_or_odb> [step]

배경
    HF 의 Step-Postbuckle 은 `*Static, stabilize=<값>` 로 돈다. 안정화는 인공 감쇠라서,
    스냅스루를 "잘 흘려보내" 완주를 만들어 줄 수 있다 — 즉 완주가 물리의 결과인지
    감쇠의 결과인지 구분해야 한다. 판정 기준(관행): 안정화 소산 에너지가 내부(변형) 에너지의
    **수 % 이하**면 신뢰, 5% 이상이면 결과를 의심한다.

읽는 것
    ODB 의 전체모델 에너지 이력(ALLIE 내부, ALLSD 안정화소산, ALLKE 운동, ETOTAL)에서
    - 각 스텝 마지막 값과 그 비 ALLSD/ALLIE
    - **시간에 대한 최대 비**(worst case) — 이게 실질 판정치다.

에너지 이력이 없으면
    스텝의 출력 요청에 `*Energy Output`(또는 `*Output, history` + 에너지 변수)이 없어서다.
    그때는 patcher 의 스텝 블록에 한 줄 추가하고 재실행해야 한다 — 스크립트가 그렇게 안내한다.
"""
from __future__ import print_function
import os
import sys

from odbAccess import openOdb

WANT = ('ALLIE', 'ALLSD', 'ALLKE', 'ETOTAL', 'ALLSE', 'ALLPD', 'ALLWK', 'ALLVD')


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    p = argv[1]
    if not p.lower().endswith('.odb'):
        p = p + '.odb'
    if not os.path.exists(p):
        print("!!! ODB 없음: %s" % p)
        return 1
    odb = openOdb(path=p, readOnly=True)
    print("=" * 72)
    print("[energy] %s" % os.path.basename(p))
    only = argv[2] if len(argv) > 2 else None
    found_any = False
    for name in odb.steps.keys():
        if only and only not in name:
            continue
        st = odb.steps[name]
        print("-" * 72)
        print("STEP %s" % name)
        for rname in st.historyRegions.keys():
            r = st.historyRegions[rname]
            have = [k for k in r.historyOutputs.keys() if k.upper().split()[0] in WANT]
            if not have:
                continue
            vals = {}
            for k in have:
                d = r.historyOutputs[k].data
                if d:
                    vals[k.upper().split()[0]] = d
            if not vals:
                continue
            found_any = True
            tmax = max(d[-1][0] for d in vals.values())
            print("  region=%s  (t_last=%.4g)" % (rname, tmax))
            for k in sorted(vals):
                print("    %-8s last=%+.6e  max=%+.6e" % (k, vals[k][-1][1],
                                                           max(abs(v) for _, v in vals[k])))
            if 'ALLSD' in vals and 'ALLIE' in vals:
                sd = {round(t, 10): v for t, v in vals['ALLSD']}
                ie = {round(t, 10): v for t, v in vals['ALLIE']}
                common = sorted(set(sd) & set(ie))
                ratios = [(abs(sd[t]) / abs(ie[t]), t) for t in common if abs(ie[t]) > 0.0]
                if ratios:
                    rmax, tmaxr = max(ratios)
                    print("    ** ALLSD/ALLIE  last=%.4g%%  max=%.4g%% (t=%.4g) **"
                          % (100.0 * (abs(sd[common[-1]]) / abs(ie[common[-1]])
                                      if abs(ie[common[-1]]) > 0 else 0.0),
                             100.0 * rmax, tmaxr))
                    print("       판정: max < 1%% 양호 / < 5%% 허용 / >= 5%% 결과 의심"
                          "(안정화가 물리를 대신했을 수 있음)")
    if not found_any:
        print("!!! 전체모델 에너지 이력이 ODB 에 없습니다.")
        print("    -> HF 스텝에 에너지 출력이 요청되지 않았습니다. patcher 의 스텝 블록에")
        print("       *Energy Output  (또는 *Output, history, variable=preselect 중 에너지)")
        print("       한 줄을 추가하고 재실행하세요.")
    odb.close()
    print("=" * 72)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
