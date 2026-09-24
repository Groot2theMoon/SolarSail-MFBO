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

DEFAULT_INSTANCE = 'MEMBRANE-1'


def write_mode_table(odb_path, step_name, out_path, n_modes=4, instance=DEFAULT_INSTANCE,
                     verbose=True):
    """ODB 모드 프레임에서 u3 를 뽑아 모드표를 쓴다. 반환: (쓴 모드 수, 메시지 리스트).

    실패를 조용히 넘기지 않고 예외로 올린다(호출측 -- run_abaqus_mode.py / CLI -- 이 판단한다).
    """
    msgs = []

    def _say(s):
        msgs.append(s)
        if verbose:
            print(s)

    if not os.path.exists(odb_path):
        raise IOError("ODB 없음: %s" % odb_path)
    from odbAccess import openOdb

    odb = openOdb(path=odb_path, readOnly=True)
    try:
        if step_name not in odb.steps:
            raise KeyError("스텝 '%s' 없음. 있는 스텝: %s"
                           % (step_name, sorted(odb.steps.keys())))
        step = odb.steps[step_name]
        inst = odb.rootAssembly.instances.get(instance)
        if inst is None:
            raise KeyError("인스턴스 '%s' 없음. 있는 것: %s"
                           % (instance, sorted(odb.rootAssembly.instances.keys())))

        frames = step.frames
        _say("[MODES] %s / %s : 프레임 %d개 (요청 %d모드) / 인스턴스 %s"
             % (os.path.basename(odb_path), step_name, len(frames), n_modes, instance))
        if len(frames) < n_modes:
            _say("[MODES] 경고: 프레임 %d개 < 요청 %d모드 -> 있는 만큼만 쓴다"
                 % (len(frames), n_modes))
            n_modes = len(frames)
        if n_modes <= 0:
            raise ValueError("모드 프레임이 없다 -> 이 스텝이 고유치 결과를 내지 못했다"
                             " (base state 가 분기점을 넘었거나 요청 수가 부족하다)")

        blocks = []
        for i in range(n_modes):
            try:
                fo = frames[i].fieldOutputs['U'].getSubset(region=inst)
            except KeyError:
                raise KeyError("모드 %d: 'U' 필드가 없다. 이 프레임의 필드: %s"
                               % (i + 1, sorted(frames[i].fieldOutputs.keys())))
            u3_at = {}
            for v in fo.values:
                u3_at[v.nodeLabel] = v.data[2]
            if not u3_at:
                raise ValueError("모드 %d: U 값이 비어 있다" % (i + 1))
            sign_node, u3_ref = max(u3_at.items(), key=lambda kv: abs(kv[1]))
            umax = abs(u3_ref)
            if umax <= 0.0:
                raise ValueError("모드 %d: u3 가 전부 0 (면외 성분 없음)" % (i + 1))
            sign = 1.0 if u3_ref >= 0 else -1.0
            lines = ["# MODE %d raw_max_u3=%.6e sign_node=%d" % (i + 1, u3_ref, sign_node),
                     "MODE %d" % (i + 1)]
            for lab in sorted(u3_at):
                lines.append("%d %.9e" % (lab, sign * u3_at[lab] / umax))
            blocks.append("\n".join(lines) + "\n")
            _say("[MODES] 모드 %d: 노드 %d개, raw max|u3|=%.3e (부호기준 노드 %d)"
                 % (i + 1, len(u3_at), umax, sign_node))

        header = ["# source=%s" % os.path.abspath(odb_path),
                  "# step=%s" % step_name,
                  "# instance=%s" % instance,
                  "# u3 은 모드별 max|u3|=1 정규화. 부호는 max|u3| 노드를 + 로 고정."]
        with open(out_path, 'w') as f:
            f.write("\n".join(header) + "\n")
            f.write("".join(blocks))
        _say("[MODES] 저장: %s (모드 %d개)" % (os.path.abspath(out_path), n_modes))
        return n_modes, msgs
    finally:
        odb.close()


def _usage():
    print("usage: abaqus python aba_mode_from_odb.py <src.odb> <step_name> <out.txt>"
          " [n_modes] [instance]")
    return 2


def main(argv):
    if len(argv) < 4:
        return _usage()
    odb_path, step_name, out_path = argv[1], argv[2], argv[3]
    n_modes = int(argv[4]) if len(argv) > 4 else 4
    inst_name = argv[5] if len(argv) > 5 else DEFAULT_INSTANCE
    try:
        write_mode_table(odb_path, step_name, out_path, n_modes, inst_name)
    except ImportError:
        print("[MODES] odbAccess 를 못 읽었습니다 - 'abaqus python' 으로 실행해야 합니다.")
        return 3
    except IOError as e:
        print("[MODES] %s" % e)
        return 4
    except KeyError as e:
        print("[MODES] %s" % e)
        return 5
    except ValueError as e:
        print("[MODES] %s" % e)
        return 7
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
