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
    λ(하중계수)는 프레임의 frameValue 에서 읽으며, **λ > 0 인 모드만** 표에 넣는다
    (요청 수를 늘리면 음수 λ 모드도 ODB 에 들어오는데, 임퍼펙션에 쓸 모드는 양수 λ 쪽이다).

종료 코드: 0 정상 / 2 인자 / 3 odbAccess 없음 / 4 odb 없음 / 5 스텝 없음 / 6 인스턴스 없음 / 7 값 없음
"""
import os
import sys

DEFAULT_INSTANCE = 'MEMBRANE-1'


def pick_mode_frames(lambdas, n_modes):
    """λ > 0 인 프레임을 앞에서부터 n_modes 개 고른다(음수/NaN λ 는 건너뛴다).

    반환: 0-based 프레임 인덱스 리스트(순서 보존, 길이 <= n_modes).
    순수 함수로 분리한 이유: odbAccess 없이 로컬에서 단위검증할 수 있게 하기 위해서.
    """
    pick = []
    for i, v in enumerate(lambdas):
        if v == v and v > 0.0:          # v == v -> NaN 제외
            pick.append(i)
            if len(pick) >= n_modes:
                break
    return pick


def extract_modes(odb_path, step_name, instance=DEFAULT_INSTANCE, n_modes=4, verbose=True,
                  dat_hint=None):
    """ODB 에서 좌굴모드를 뽑아 '파일 없이' 표를 돌려준다 (run_abaqus.py 의 ODB 직접 주입용).

    반환: (table, meta, picks, lam, msgs)
        table = {mode(1부터): {node_label: u3}}   u3 는 모드별 max|u3|=1 정규화 + 부호 고정
        meta  = {mode: {'frame': i+1, 'lam': λ, 'raw_max_u3': ..., 'sign_node': ...}}
        picks = 선택된 0-based 프레임 인덱스, lam = 프레임별 λ
    규칙(write_mode_table 과 동일): λ > 0 인 프레임만 쓴다.
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
        # (2026-09-28 실측) odb.rootAssembly.instances 는 Abaqus 의 Repository 라서 .get() 이 없다.
        #   .get() 을 쓰면 AttributeError: 'Repository' object has no attribute 'get' 로 죽고,
        #   호출측의 try/except 에 삼켜져 '모드표가 조용히 안 만들어지는' 원인이 됐다.
        #   Repository 는 in / [] / keys() 를 지원한다.
        inst = (odb.rootAssembly.instances[instance]
                if instance in odb.rootAssembly.instances else None)
        if inst is None:
            raise KeyError("인스턴스 '%s' 없음. 있는 것: %s"
                           % (instance, sorted(odb.rootAssembly.instances.keys())))

        frames = step.frames
        lam = []
        for fr in frames:
            try:
                lam.append(float(fr.frameValue))
            except Exception:
                lam.append(float('nan'))
        _say("[MODES] %s / %s : 프레임 %d개 (요청 %d모드) / 인스턴스 %s"
             % (os.path.basename(odb_path), step_name, len(frames), n_modes, instance))
        # (2026-09-28 실측) *BUCKLE 스텝 ODB 의 frame.frameValue 는 고유치가 아니라 '모드 번호'(1,2,3,4)다.
        #   좌굴계수 λ 는 .dat 의 MODE NO 표에 있으므로 dat_hint 가 있으면 그쪽을 우선한다.
        #   (λ 는 모드 선택/로그 표기용이다 — 노드 섭동 물리에는 영향이 없다. 그래도 틀린 값을 쓰지 않는다.)
        lam_src = 'odb.frameValue'
        if dat_hint and os.path.exists(dat_hint):
            try:
                from coalescence_check import parse_eigenvalues
                _dl = [float(x) for x in parse_eigenvalues(open(dat_hint).read())]
            except Exception as _e_dat:
                _dl = []
                _say("[MODES] .dat 고유치 파싱 실패(%s: %s) -> ODB frameValue 사용"
                     % (type(_e_dat).__name__, _e_dat))
            if _dl:
                lam = _dl
                lam_src = os.path.basename(dat_hint)
        if lam:
            _say("[MODES] λ 출처=%s (앞 12개): %s"
                 % (lam_src, ', '.join('%.6e' % v for v in lam[:12])))

        picks = pick_mode_frames(lam, n_modes)
        if not picks:
            raise ValueError("양수 λ 모드가 없다 -> 좌굴모드가 수렴하지 않았다. λ: %s"
                             % ', '.join('%.6e' % v for v in lam[:12]))
        if len(picks) < n_modes:
            _say("[MODES] 경고: 쓸 수 있는 양수 λ 모드 %d개 < 요청 %d개 -> 있는 만큼만 쓴다"
                 % (len(picks), n_modes))
        # 방어: λ 를 .dat 에서 읽었으면 .dat 이 ODB 프레임 수보다 많은 모드를 나열할 수 있다.
        #   (예: 100개 요청 중 4개만 수렴 -> .dat 에는 여러 값, ODB 에는 4개 프레임)
        _over = [p for p in picks if p >= len(frames)]
        if _over:
            _say("[MODES] 경고: λ 출처(%s)의 모드 %s 는 ODB 프레임(%d개) 범위 밖 -> 제외"
                 % (lam_src, [p + 1 for p in _over], len(frames)))
            picks = [p for p in picks if p < len(frames)]
            if not picks:
                raise ValueError("λ 는 %d개인데 ODB 프레임이 %d개다 -> 모드 추출 불가"
                                 % (len(lam), len(frames)))
        _say("[MODES] 선택한 프레임 %s (λ %s)"
             % ([p + 1 for p in picks], ['%.6e' % lam[p] for p in picks]))

        table, meta = {}, {}
        for k, i in enumerate(picks):
            try:
                fo = frames[i].fieldOutputs['U'].getSubset(region=inst)
            except KeyError:
                raise KeyError("프레임 %d: 'U' 필드가 없다. 이 프레임의 필드: %s"
                               % (i + 1, sorted(frames[i].fieldOutputs.keys())))
            u3_at = {}
            for v in fo.values:
                u3_at[v.nodeLabel] = v.data[2]
            if not u3_at:
                raise ValueError("프레임 %d: U 값이 비어 있다" % (i + 1))
            sign_node, u3_ref = max(u3_at.items(), key=lambda kv: abs(kv[1]))
            umax = abs(u3_ref)
            if umax <= 0.0:
                raise ValueError("프레임 %d: u3 가 전부 0 (면외 성분 없음)" % (i + 1))
            sign = 1.0 if u3_ref >= 0 else -1.0
            table[k + 1] = dict((lab, sign * val / umax) for lab, val in u3_at.items())
            meta[k + 1] = {'frame': i + 1, 'lam': lam[i], 'raw_max_u3': umax,
                           'sign_node': sign_node}
            _say("[MODES] 모드 %d <- 프레임 %d: 노드 %d개, λ=%.6e, raw max|u3|=%.3e (부호기준 노드 %d)"
                 % (k + 1, i + 1, len(u3_at), lam[i], umax, sign_node))
        return table, meta, picks, lam, msgs
    finally:
        odb.close()


def write_mode_table(odb_path, step_name, out_path, n_modes=4, instance=DEFAULT_INSTANCE,
                     verbose=True, extra_header=None):
    """ODB 모드에서 모드표 txt 를 쓴다. 반환: (쓴 모드 수, 메시지 리스트).

    파일을 거치지 않고 ODB 에서 바로 읽으려면 run_abaqus.py 의 IMPERFECTION_MODE='odb_direct'
    를 쓴다(같은 extract_modes 를 호출하므로 결과가 같다).
    """
    table, meta, picks, lam, msgs = extract_modes(odb_path, step_name, instance, n_modes, verbose)
    header = ["# source=%s" % os.path.abspath(odb_path),
              "# step=%s" % step_name,
              "# instance=%s" % instance,
              "# u3 은 모드별 max|u3|=1 정규화. 부호는 max|u3| 노드를 + 로 고정.",
              "# table_modes=%d / source_frames=%s / lambdas=%s"
              % (len(picks), [p + 1 for p in picks], ['%.6e' % lam[p] for p in picks])]
    header += list(extra_header or [])      # 예: '# FINGERPRINT sha1=...'
    blocks = []
    for k in sorted(table):
        lines = ["# MODE %d frame=%d lambda=%.6e raw_max_u3=%.6e sign_node=%d"
                 % (k, meta[k]['frame'], meta[k]['lam'], meta[k]['raw_max_u3'],
                    meta[k]['sign_node']),
                 "MODE %d" % k]
        for lab in sorted(table[k]):
            lines.append("%d %.9e" % (lab, table[k][lab]))
        blocks.append("\n".join(lines) + "\n")
    with open(out_path, 'w') as f:
        f.write("\n".join(header) + "\n")
        f.write("".join(blocks))
    _say_save = "[MODES] 저장: %s (모드 %d개)" % (os.path.abspath(out_path), len(table))
    msgs.append(_say_save)
    if verbose:
        print(_say_save)
    return len(table), msgs




def _usage():
    print("usage: abaqus python aba_mode_from_odb.py <src.odb> <step_name> <out.txt>"
          " [n_modes] [instance]")
    return 2


def main(argv):
    if len(argv) < 4:
        return _usage()
    odb_path, step_name, out_path = argv[1], argv[2], argv[3]
    try:
        n_modes = int(argv[4]) if len(argv) > 4 else 4
    except ValueError:
        print("[MODES] 인자 4(모드 수)가 숫자가 아닙니다: %r"
              "  <- 명령을 한 줄로 붙여넣어 뒤 인자가 섞였는지 확인하세요." % (argv[4],))
        return 2
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
    except Exception as e:
        # 예상 못 한 예외도 traceback 으로만 끝나지 않게 ASCII 한 줄 + 고유 종료코드로 남긴다.
        print("[MODES] UNEXPECTED %s: %s" % (type(e).__name__, e))
        return 8
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
