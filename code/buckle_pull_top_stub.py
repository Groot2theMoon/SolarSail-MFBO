#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""RETIRED 2026-10-21: run_abaqus_buckle.py 의 PULL_TOP 스위치가 제거되었다.

왜 제거했나
    PULL_TOP=True(상단 꼭짓점 케이블 방향 구동)는 CASE='paper' +
    PRETENSION_MODE='paper3' 와 역할이 완전히 겹쳤고, 기본값이 False 라
    프로덕션 경로에서 한 번도 쓰이지 않았다. 유지보수 비용(BC 배선 37줄 +
    job 이름 _3t 태그 + 이 전용 하네스)만 남아 있어 제거했다.

3-tendon(상단 꼭짓점 구동)을 다시 볼 때
    - 모드 소스: CASE='paper' (PRETENSION_MODE='paper3', 상단/코너 2:1)
    - BC 배선 검증: python3 buckle_stub_harness.py paper 0.5
    - 과거 PULL_TOP 구현은 git 이력에서 확인한다.

이 파일은 더 이상 테스트를 수행하지 않는다(기존 PULL_TOP 소스가 없으므로).
"""
from __future__ import print_function

import sys


def main():
    print("[buckle_pull_top_stub] RETIRED: PULL_TOP was removed on 2026-10-21.")
    print("  Use CASE='paper' (PRETENSION_MODE='paper3') for the 3-tendon source.")
    print("  BC trace: python3 buckle_stub_harness.py paper 0.5")
    return 0


if __name__ == '__main__':
    sys.exit(main())
