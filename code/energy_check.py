"""에너지 이력(ALLIE / ALLSD / ALLKE) 판독기 — 감쇠가 결과를 만들고 있는지 정량화한다.

사용법:
    abaqus python energy_check.py aba\HF_Postbuckle.odb

배경 (Abaqus 문서 7.1.1 / 일반 지침):
    *STATIC, STABILIZE 로 넣는 인공감쇠의 가격은 ALLSD(정적 감쇠 소산에너지) / ALLIE(총내부에너지) 다.
    Abaqus 매뉴얼은 ALLSD 를 ALLIE 의 5% 이하로 유지하라고 권고한다.
    - <= 5%  : 주름이 감쇠 산물이 아니라고 주장 가능
    - >  5%  : 주름이 감쇠로 지지되고 있을 가능성 -> 요소차수/메쉬 변경 등 물리적 해상도 개선 필요

이 스크립트는 해석을 하지 않는다 (읽기 전용). 실패해도 기존 산출물에 영향 없음.
"""
from __future__ import print_function

import os
import sys

from odbAccess import openOdb

# 판정 임계값 (Abaqus 권고)
LIMIT = 0.05


def _last(hist):
    """historyOutput 객체에서 마지막 (time, value) 를 꺼낸다."""
    try:
        d = hist.data
        if not d:
            return None, None
        return d[-1][0], d[-1][1]
    except Exception:
        return None, None


def main():
    if len(sys.argv) < 2:
        print("usage: abaqus python energy_check.py <odb>")
        return 1

    path = sys.argv[1]
    if not os.path.exists(path):
        print("ERROR: odb 없음: %s" % path)
        return 1

    odb = openOdb(path=path, readOnly=True)
    try:
        print("=" * 78)
        print("에너지 이력 판독: %s" % path)
        print("  판정 기준: ALLSD/ALLIE <= %.0f%% (Abaqus 권고)" % (100.0 * LIMIT))
        print("=" * 78)

        for step_name in odb.steps.keys():
            step = odb.steps[step_name]
            regions = step.historyRegions.keys()
            print("\n[%s]" % step_name)

            # 전체모델 에너지는 보통 'Assembly ASSEMBLY' 에 들어간다. 없으면 전부 훑는다.
            cand = []
            for rname in regions:
                ho = step.historyRegions[rname].historyOutputs.keys()
                if 'ALLIE' in ho or 'ALLSD' in ho:
                    cand.append(rname)
            if not cand:
                print("  (에너지 이력 없음)  historyRegions = %s" % list(regions))
                print("  -> H-Energy 출력요청이 모델에 없거나, 이 스텝에서 비활성")
                continue

            for rname in cand:
                ho = step.historyRegions[rname].historyOutputs
                t = {}
                for key in ('ALLIE', 'ALLSD', 'ALLKE', 'ETOTAL', 'ALLSE'):
                    if key in ho:
                        tt, vv = _last(ho[key])
                        if vv is not None:
                            t[key] = vv
                            if tt is not None:
                                t['_t_' + key] = tt

                ie = t.get('ALLIE')
                sd = t.get('ALLSD')
                ke = t.get('ALLKE')
                print("  region=%s" % rname)
                if ie is not None:
                    print("    ALLIE  = %.6e J   (총 내부에너지, t=%.4g)" % (ie, t.get('_t_ALLIE', -1)))
                if sd is not None:
                    print("    ALLSD  = %.6e J   (정적 감쇠 소산, t=%.4g)" % (sd, t.get('_t_ALLSD', -1)))
                if ke is not None:
                    print("    ALLKE  = %.6e J   (운동에너지)" % ke)

                if ie is None or sd is None:
                    print("    -> ALLIE/ALLSD 중 누락: 비율 계산 불가")
                    continue
                if ie <= 0.0:
                    print("    -> ALLIE 가 0 이하: 비율 계산 불가 (거의 무응력 상태?)")
                    continue

                ratio = sd / ie
                verdict = "OK (<=5%%)" if ratio <= LIMIT else "주의 (>5%%)"
                print("    >>> ALLSD/ALLIE = %.5f  (%.2f%%)  -> %s"
                      % (ratio, 100.0 * ratio, verdict))
                if ratio > LIMIT:
                    print("        주름이 인공감쇠로 지지되고 있을 가능성. "
                          "요소차수/메쉬 해상도 개선을 검토할 것.")
                # 운동에너지도 크면 준정적 가정이 깨진 것
                if ke is not None and ke > 0.01 * ie:
                    print("        (경고) ALLKE/ALLIE = %.4f : 준정적 가정 확인 필요"
                          % (ke / ie))

        print("\n" + "=" * 78)
    finally:
        odb.close()
    return 0


if __name__ == '__main__':
    sys.exit(main())
