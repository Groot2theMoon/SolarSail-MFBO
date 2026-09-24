# -*- coding: utf-8 -*-
"""ODB 에서 고유모드(좌굴/주파수)를 뽑아 '기하 섭동용 모드표'로 저장한다 (2026-09-22)

왜 .fil 이 아니라 ODB 인가
    *BUCKLE 스텝은 .fil 출력이 금지된다(실측: ClampFree_Buckle.dat:7025
    "***WARNING: FILE OUTPUT IS NOT AVAILABLE FOR BUCKLING ANALYSIS").
    반면 .odb 에는 모드 프레임이 항상 들어 있다(CAE/Viewer 로 모드를 볼 수 있는 그 데이터).
    업스트림 원본 run_abaqus.py 의 docstring 도 이 경로를 계획으로 적어 두었다:
      "2차: 1차 해석 결과(.odb)에서 고유모드를 추출하여 초기 결함(Imperfection)으로 주입"

사용
    abaqus python aba_mode_from_odb.py <src.odb> <step_name> <out.txt> [n_modes] [instance]
예
    abaqus python aba_mode_from_odb.py Buckle_Analysis.odb Step-Buckle modes_ClampFree_Buckle.txt 4

출력 형식 (aba_imperfection.load_mode_table 이 읽는다)
    # MODE 1 raw_max_u3=1.234e-03 sign_node=1234
    MODE 1
    <node_label> <u3_normalized>
    ...
    u3 를 max|u3| = 1 로 모드별 정규화한다(모드 형상의 절대 크기는 의미가 없다).
    부호는 max|u3| 노드를 + 로 고정 -> 같은 ODB 에서 항상 같은 표가 나온다(재현성).

종료 코드: 0 정상 / 2 인자 / 3 odbAccess 없음 / 4 odb 없음 / 5 스텝 없음 / 6 인스턴스 없음 / 7 값 없음
"""
import os
import sys


def _usage():
    print("usage: abaqus python aba_mode_from_odb.py <src.odb> <step_name> <out.txt>"
          " [n_modes] [instance]")
    return 2


def main(argv):
    if len(argv) < 4:
        return _usage()
    odb_path, step_name, out_path = argv[1], argv[2], argv[3]
    n_modes = int(argv[4]) if len(argv) > 4 else 4
    inst_name = argv[5] if len(argv) > 5 else 'MEMBRANE-1'

    try:
        from odbAccess import openOdb
    except ImportError:
        print("[MODES] odbAccess 를 못 읽었습니다 - 'abaqus python' 으로 실행해야 합니다.")
        return 3
    if not os.path.exists(odb_path):
        print("[MODES] 파일 없음: %s" % odb_path)
        return 4

    odb = openOdb(path=odb_path, readOnly=True)
    try:
        if step_name not in odb.steps:
            print("[MODES] 스텝 '%s' 없음. 있는 스텝: %s"
                  % (step_name, sorted(odb.steps.keys())))
            return 5
        step = odb.steps[step_name]
        try:
            inst = odb.rootAssembly.instances[inst_name]
        except KeyError:
            print("[MODES] 인스턴스 '%s' 없음. 있는 것: %s"
                  % (inst_name, sorted(odb.rootAssembly.instances.keys())))
            return 6

        frames = step.frames
        print("[MODES] %s / %s : 프레임 %d개 (요청 %d모드) / 인스턴스 %s"
              % (odb_path, step_name, len(frames), n_modes, inst_name))
        if len(frames) < n_modes:
            print("[MODES] 경고: 프레임 %d개 < 요청 %d모드 -> 있는 만큼만 쓴다"
                  % (len(frames), n_modes))
            n_modes = len(frames)
        if n_modes <= 0:
            print("[MODES] 모드 프레임이 없다 -> 이 스텝은 고유치 해석이 아니다")
            return 7

        out = open(out_path, 'w')
        try:
            out.write("# source=%s\n# step=%s\n# instance=%s\n"
                      % (os.path.abspath(odb_path), step_name, inst_name))
            out.write("# u3 은 모드별 max|u3|=1 정규화. 부호는 max|u3| 노드를 + 로 고정.\n")
            total = 0
            for i in range(n_modes):
                try:
                    fo = frames[i].fieldOutputs['U'].getSubset(region=inst)
                except KeyError:
                    print("[MODES] 모드 %d: 'U' 필드가 없다 -> 중단" % (i + 1))
                    return 7
                u3_at = {}
                for v in fo.values:
                    u3_at[v.nodeLabel] = v.data[2]
                if not u3_at:
                    print("[MODES] 모드 %d: U 값이 비어 있다 -> 중단" % (i + 1))
                    return 7
                sign_node, u3_ref = max(u3_at.items(), key=lambda kv: abs(kv[1]))
                umax = abs(u3_ref)
                if umax <= 0.0:
                    print("[MODES] 모드 %d: u3 가 전부 0 (면외 성분 없음) -> 중단" % (i + 1))
                    return 7
                sign = 1.0 if u3_ref >= 0 else -1.0
                out.write("# MODE %d raw_max_u3=%.6e sign_node=%d\n"
                          % (i + 1, u3_ref, sign_node))
                out.write("MODE %d\n" % (i + 1))
                for lab in sorted(u3_at):
                    out.write("%d %.9e\n" % (lab, sign * u3_at[lab] / umax))
                total += len(u3_at)
                print("[MODES] 모드 %d: 노드 %d개, raw max|u3|=%.3e (기준 노드 %d)"
                      % (i + 1, len(u3_at), umax, sign_node))
        finally:
            out.close()
        print("[MODES] 저장: %s (모드 %d개, 총 %d행)" % (os.path.abspath(out_path), n_modes, total))
    finally:
        odb.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
