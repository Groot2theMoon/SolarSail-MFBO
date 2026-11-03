#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ODB step/frame/field 목록 출력 (진단용, Abaqus Python 전용).

사용:
    abaqus python check_odb_fields.py aba\\LF_Analysis_xc0p5_dc2.odb

목적:
  `eval_abaqus.py` 가 `ZK:FAIL reason=E_field_missing` 을 낼 때,
  (a) ODB 어느 스텝/프레임에 E 가 있는지,
  (b) 선택 프레임에 S/EVOL 만 있고 E 가 없는지
  를 바로 확인한다.
"""
from __future__ import print_function

import sys

try:
    from odbAccess import openOdb
except Exception as exc:                      # pragma: no cover
    print("ERROR: odbAccess import failed: %s" % exc)
    print("  abaqus python check_odb_fields.py <odb> 로 실행하세요.")
    sys.exit(1)


def main(argv):
    if len(argv) < 2:
        print("usage: abaqus python check_odb_fields.py <odb>")
        return 2
    path = argv[1]
    try:
        odb = openOdb(path=path, readOnly=True)
    except Exception as exc:
        print("ERROR: openOdb failed: %s" % exc)
        return 1
    try:
        print("ODB: %s" % path)
        for step_name in odb.steps.keys():
            step = odb.steps[step_name]
            print("STEP %s frames=%d" % (step_name, len(step.frames)))
            if not step.frames:
                print("  !!! no frames")
                continue
            for i, frame in enumerate(step.frames):
                try:
                    fields = sorted(frame.fieldOutputs.keys())
                except Exception as exc:
                    fields = ['<error:%s>' % type(exc).__name__]
                print("  frame[%d] value=%s fields=%s"
                      % (i, frame.frameValue, ','.join(fields)))
            _last = step.frames[-1].fieldOutputs.keys()
            if ('E' not in _last) and ('LE' not in _last):
                print("  !!! last frame has neither E nor LE")
            else:
                print("  OK: last frame has %s"
                      % ('E' if 'E' in _last else 'LE'))
    finally:
        try:
            odb.close()
        except Exception:
            pass
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
