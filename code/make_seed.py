#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""설계점 (x_c, d_c) 의 HF 임퍼펙션 시드표를 만든다 (2026-10-12).

    python3 make_seed.py --x_c 0.5 --d_c 0.0556                # 생성(캐시 있으면 건너뜀)
    python3 make_seed.py --x_c 0.25 --d_c 0.0556 --force       # 강제 재생성
    python3 make_seed.py --x_c 0.5 --d_c 0.0556 --dry-run      # 실행할 명령만 출력(라이선스 0)

무엇을 하는가
    ① seed  가족(꼭짓점 케이블만, x_c 무관) 좌굴 4모드  -> seeds/modes_seed.txt      [1회만]
    ② clamp 가족(클램프 케이블만, x_c·d_c 의존) 2모드   -> seeds/modes_clamp_xc<NNN>_dc<NNN>.txt
    ③ 2+2 합성 + 설계점/지문 헤더                       -> seeds/modes_xc<NNN>_dc<NNN>.txt   <- HF 가 읽는 표
    ④ 검증(노드 수·모드 수·peak|u3|) + seeds/<표>.meta.txt

설계 결정 (code/notes/plan-hf-seed-pipeline-2026-10-12.md)
    Q1=a -> §19.69 에서 **부분 철회**: HF d_c 는 표 이름/헤더(짝 규약)에만 쓰고,
          시드 클램프 당김은 **고정 5.56%**(--seed-clamp-dc)로 뽑는다(d_c=0.5 급에서 모드 0개 실측).
    Q2=a  지문 불일치는 HF 가 **제출 전 중단**한다(이 스크립트는 지문을 표 머리에 새긴다)
    Q3=b  CASE 전환은 소스 패치 - run_abaqus_buckle.py 의 `CASE = '...'` 줄을 잠시 바꾸고
          **finally 에서 원문 바이트로 복원**한다(실패해도 남기지 않는다)
"""
from __future__ import print_function

import argparse
import glob
import io
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
#   [§19.68] Windows 콘솔/파이프 인코딩(cp949) 때문에 자식이 출력 중 죽는 것을 막는다.
#     (실측: aba_mode_from_odb.py 가 파이프로 잡힌 stdout 에서 '—' 를 못 써 죽었다.)
os.environ.setdefault('PYTHONIOENCODING', 'utf-8')
os.environ.setdefault('PYTHONUTF8', '1')
_CHILD_ENV = dict(os.environ, PYTHONIOENCODING='utf-8', PYTHONUTF8='1')
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
sys.path.insert(0, HERE)
import aba_imperfection as A                                                     # noqa: E402

BUCKLE_SCRIPT = os.path.join(HERE, A.BUCKLE_SCRIPT_NAME)
SEEDS_DIR = os.path.join(HERE, A.SEED_DIR_NAME)
RUN_DIR = os.path.join(HERE, 'buckle')
BUCKLE_MODES = 4          # seed 가족에서 뽑는 모드 수(상위 2개만 합성에 쓴다)
CLAMP_MODES = 2           # clamp 가족


# --------------------------------------------------------------------- helpers
def sh(cmd, dry=False, cwd=HERE):
    print("  $ %s" % cmd)
    if dry:
        return 0, ''
    p = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True,
                       encoding='utf-8', errors='replace', env=_CHILD_ENV)
    out = (p.stdout or '')
    #   `abaqus python` 런처는 스크립트가 예외로 죽어도 rc=0 을 돌려줄 때가 있다(실측) ->
    #   종료코드를 믿지 않고 **산출물로 검증**한다. 그래서 자식 출력 꼬리는 항상 보여준다.
    tail = [l for l in out.strip().splitlines() if l.strip()][-6:]
    for l in tail:
        print("  | %s" % l)
    if p.returncode != 0:
        print("  !!! exit=%d" % p.returncode)
        print((p.stderr or '')[-800:])
    return p.returncode, out


def _need(path, what, hint=''):
    """산출물 검증 — rc 가 아니라 파일로 판정한다."""
    if not os.path.exists(path):
        raise SystemExit("[make_seed] %s 가 만들어지지 않았습니다 -> %s\n"
                         "  (위에 찍힌 자식 출력의 마지막 줄을 보세요.)%s"
                         % (what, path, ('\n  ' + hint) if hint else ''))
    if os.path.getsize(path) == 0:
        raise SystemExit("[make_seed] %s 가 비어 있습니다 -> %s" % (what, path))
    print("  [ok] %s (%d bytes)" % (what, os.path.getsize(path)))


def buckle_disp():
    """run_abaqus_buckle.py 의 DISP_GLOBAL(또는 MODE_PRESTRESS)을 읽는다."""
    txt = io.open(BUCKLE_SCRIPT, encoding='utf-8', errors='replace').read()
    m = re.search(r'^DISP_GLOBAL\s*=\s*([0-9.eE+-]+)', txt, re.M)
    if m:
        return float(m.group(1))
    m = re.search(r'^MODE_PRESTRESS\s*=\s*([0-9.eE+-]+)', txt, re.M)
    return float(m.group(1)) if m else 1.8e-5


def buckle_mesh_tag():
    """run_abaqus_buckle.py 의 PAPER_MATCH/PAPER_MESH_ORDER 를 읽어 메쉬 태그를 맞춘다.

    2차 격자(PAPER_MESH_ORDER=2)면 잡 이름이 'grid2' 다 — make_seed 가 1차 이름으로
    ODB 를 찾으면 항상 빈손이 되어 잡을 불필요하게 다시 돌린다(또는 stale 을 집는다).
    """
    txt = io.open(BUCKLE_SCRIPT, encoding='utf-8', errors='replace').read()
    _paper = re.search(r'^PAPER_MATCH\s*=\s*True\b', txt, re.M) is not None
    _o2 = re.search(r'^PAPER_MESH_ORDER\s*=\s*2\b', txt, re.M) is not None
    return 'grid2' if (_paper and _o2) else 'grid'


def seed_case():
    """모드 추출은 **seed + clamp_only 합성** (사용자 결정 2026-10-20).

    cable-free 'paper' 케이스는 run_abaqus_buckle.py 에 수동 실행용으로 남겨 둔다
    (모드표의 1가족으로 쓰지 않는다 — 사용자는 두 가족 합성을 선택했다).
    """
    return 'seed'


def buckle_pert_tag():
    """run_abaqus_buckle.py 의 PERTURBATION -> 잡 이름의 _plvNNNNN 태그."""
    txt = io.open(BUCKLE_SCRIPT, encoding='utf-8', errors='replace').read()
    m = re.search(r'^PERTURBATION\s*=\s*([0-9.eE+-]+)', txt, re.M)
    p = float(m.group(1)) if m else 0.01
    return '_plv%05d' % int(round(p * 1.0e5))


def _table_nodes(path):
    """MODE k 블록마다 노드 줄 수를 센다 (표가 현재 메쉬와 같은 라벨 집합인지 판정)."""
    counts, cur = {}, None
    with io.open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            s = line.strip()
            if s.startswith('MODE '):
                cur = s.split()[1]
                counts[cur] = 0
            elif cur and s and not s.startswith('#'):
                counts[cur] += 1
    return counts


def expected_membrane_nodes():
    """PAPER_MATCH + PAPER_MESH_ORDER=2 -> 30,501 (S8R5/STRI65); 그 외 1차 10,201."""
    txt = io.open(BUCKLE_SCRIPT, encoding='utf-8', errors='replace').read()
    _paper = re.search(r'^PAPER_MATCH\s*=\s*True\b', txt, re.M) is not None
    _o2 = re.search(r'^PAPER_MESH_ORDER\s*=\s*2\b', txt, re.M) is not None
    return 30501 if (_paper and _o2) else 10201


def check_table_nodes(path, what):
    exp = expected_membrane_nodes()
    counts = _table_nodes(path)
    bad = dict((k, v) for k, v in counts.items() if v != exp)
    if not counts or bad:
        raise SystemExit(
            "[make_seed] %s 의 노드 수가 기대 %d 와 다릅니다: %s\n"
            "  -> 1차 표를 2차 메쉬에 쓰려 했거나, buckle 잡이 실패해 **이전 표**를 읽었습니다.\n"
            "  -> buckle/<job>.msg·.dat 의 원인을 먼저 고치고 --force 로 다시 실행하세요.\n"
            "  table=%s" % (what, exp, bad or 'MODE 없음', path))


def fresh(path, t0):
    """이번 실행에서 갱신된 파일인가 (Windows mtime 해상도 여유 1초)."""
    try:
        return os.path.exists(path) and os.path.getmtime(path) >= (t0 - 1.0)
    except OSError:
        return False


def expected_odb(x_c, disp, d_c, case, excl=0.21):
    """buckle 잡 이름 규약으로 ODB 경로를 만든다(없으면 glob 폴백)."""
    dc = A.value_tag(d_c) if case == 'clamp_only' else '0'     # [§19.76] 반올림 없음
    ex = '' if excl == 0.0 else '_ex%03d' % int(round(excl * 1000.0))
    job = ('Buckle_xc%s_d%03dum_dc%s_f%03d_s4_%s_%s_normal%s%s'
           % (A.value_tag(x_c), int(round(disp * 1.0e6)), dc, 100, buckle_mesh_tag(), case, ex,
              buckle_pert_tag()))
    p = os.path.join(RUN_DIR, job + '.odb')
    if os.path.exists(p):
        return p, job
    pat = os.path.join(RUN_DIR, 'Buckle_xc%s_*_%s_*.odb' % (A.value_tag(x_c), case))
    hits = sorted(glob.glob(pat), key=os.path.getmtime, reverse=True)
    return (hits[0], os.path.splitext(os.path.basename(hits[0]))[0]) if hits else (p, job)


class CasePatch(object):
    """Q3=b: CASE 줄을 잠시 바꾸고 finally 에서 **원문 바이트로** 복원한다."""

    def __init__(self, case, dry=False):
        self.case, self.dry, self.orig = case, dry, None

    def __enter__(self):
        if self.dry:
            print("  [patch] CASE = '%s' (dry-run: 파일 안 건드림)" % self.case)
            return self
        self.orig = io.open(BUCKLE_SCRIPT, 'rb').read()
        txt = self.orig.decode('utf-8')
        new, n = re.subn(r"^CASE\s*=\s*'[^']*'", "CASE = '%s'" % self.case, txt, count=1, flags=re.M)
        if n != 1:
            raise RuntimeError('CASE 줄을 찾지 못했습니다(1개여야 함): %d' % n)
        if new == txt:
            print("  [patch] CASE 는 이미 '%s' - 패치 생략" % self.case)
            self.orig = None
        else:
            io.open(BUCKLE_SCRIPT, 'wb').write(new.encode('utf-8'))
            print("  [patch] CASE -> '%s'" % self.case)
        return self

    def __exit__(self, *exc):
        if self.dry or self.orig is None:
            return False
        cur = io.open(BUCKLE_SCRIPT, 'rb').read()
        io.open(BUCKLE_SCRIPT, 'wb').write(self.orig)
        back = io.open(BUCKLE_SCRIPT, 'rb').read()
        if back != self.orig:
            raise RuntimeError('CASE 복원 실패 - run_abaqus_buckle.py 를 수동 확인하세요')
        print("  [patch] CASE 원문 복원 완료%s"
              % ('' if cur != self.orig else ' (변경 없음)'))
        return False


# --------------------------------------------------------------------- 단계들
def main(argv=None):
    ap = argparse.ArgumentParser(description='설계점 HF 시드표 생성')
    ap.add_argument('--x_c', type=float, required=True)
    ap.add_argument('--d_c', type=float, required=True,
                    help='클램프 당김 비율(clamp_pull = disp x d_c). HF 의 x2 와 같은 값')
    ap.add_argument('--disp', type=float, default=None,
                    help='buckle 코너 당김 [m]. 생략하면 run_abaqus_buckle.py 의 DISP_GLOBAL')
    ap.add_argument('--excl', type=float, default=0.21)
    ap.add_argument('--seed-clamp-dc', type=float, default=0.0, dest='seed_clamp_dc',
                    help='[§19.69/§19.70] 시드 클램프 당김 비율. 기본 0.0 = §19.52 에서 실제로 100/100 을 '
                         '낸 설정(clamp_only, dc000). HF 의 d_c 와 분리한다 - 맞추면 d_c=0.5 급에서 '
                         '클램프 base 가 퇴화해 선형모드가 0개 나온다(실측).')
    ap.add_argument('--force', action='store_true', help='캐시 무시하고 재생성')
    ap.add_argument('--dry-run', action='store_true', help='명령만 출력(라이선스 0)')
    ap.add_argument('--refresh-seed-family', action='store_true',
                    help='x_c 무관한 seed 가족도 다시 뽑는다')
    a = ap.parse_args(argv)

    disp = a.disp if a.disp is not None else buckle_disp()
    pull = disp * a.d_c                        # HF 와 짝지을 설계점 값(표 이름/헤더에 쓴다)
    pull_seed = disp * a.seed_clamp_dc         # 시드 클램프 당김 - 고정 비율(§19.69)
    seed_tbl = os.path.join(SEEDS_DIR, 'modes_seed.txt')
    clamp_tbl = os.path.join(SEEDS_DIR, A.seed_table_name(a.x_c, a.seed_clamp_dc, prefix='modes_clamp'))
    out_tbl = os.path.join(SEEDS_DIR, A.composite_table_name(a.x_c))   # [§19.77] d_c 무관
    if not a.dry_run:
        os.makedirs(SEEDS_DIR, exist_ok=True)

    print("[make_seed] x_c=%.4g d_c(HF)=%.5g disp=%.4g clamp_pull(HF)=%.4g excl=%.3g" %
          (a.x_c, a.d_c, disp, pull, a.excl))
    print("[make_seed] 시드 클램프 비율=%.4g -> clamp_pull(시드)=%.4g m "
          "(HF d_c 와 분리 - 비율이 크면 클램프 base 가 퇴화해 모드 0개)"
          % (a.seed_clamp_dc, pull_seed))
    print("[make_seed] seed  표 = %s" % seed_tbl)
    print("[make_seed] clamp 표 = %s" % clamp_tbl)
    print("[make_seed] 합성  표 = %s" % out_tbl)

    # ---- 1) seed 가족 (x_c 무관, 1회) ----
    _sc = seed_case()
    if a.force or a.refresh_seed_family or not os.path.exists(seed_tbl):
        print("\n[1/3] seed 가족 (%s 케이스) - x_c 무관" % _sc)
        odb, _job = expected_odb(a.x_c, disp, 0.0, _sc, a.excl)
        if os.path.exists(odb) and not a.force:
            print("      buckle 산출물 있음 - 잡 생략 (%s)" % os.path.basename(odb))
        else:
            _job_t0 = time.time()
            with CasePatch(_sc, a.dry_run):
                rc, _ = sh('abaqus cae noGUI=run_abaqus_buckle.py -- %g %g 0' % (a.x_c, disp),
                           a.dry_run, cwd=HERE)
                if rc:
                    return rc
            odb, _job = expected_odb(a.x_c, disp, 0.0, _sc, a.excl)
            if not a.dry_run and not fresh(odb, _job_t0):
                raise SystemExit(
                    "[make_seed] buckle ODB 가 **이번 실행에서 갱신되지 않았습니다** — CAE/해석이 실패했습니다.\n"
                    "  이전 ODB 를 읽어 모드표를 만들면 안 됩니다. buckle/*.msg·.dat 의 원인을 확인하세요.\n"
                    "  odb=%s" % odb)
        _t0 = time.time()
        rc, _ = sh('abaqus python aba_mode_from_odb.py "%s" Step-Buckle "%s" %d'
                   % (odb, seed_tbl, BUCKLE_MODES), a.dry_run, cwd=HERE)
        if not a.dry_run:
            _need(seed_tbl, 'seed 가족 모드표')
            if not fresh(seed_tbl, _t0):
                raise SystemExit(
                    "[make_seed] seed 모드표가 **이번 실행에서 갱신되지 않았습니다** — buckle 잡이 실패해\n"
                    "  이전 표가 그대로 남아 있습니다. 이 표를 쓰면 1차/낡은 모드가 HF 에 주입됩니다.\n"
                    "  buckle/*.msg·.dat 의 원인을 먼저 고치고 --force 로 다시 실행하세요.\n"
                    "  table=%s" % seed_tbl)
            check_table_nodes(seed_tbl, 'seed 가족 모드표')
    else:
        print("\n[1/3] seed 가족 표 있음 - 건너뜀 (%s)" % os.path.basename(seed_tbl))

    # ---- 2) clamp 가족 (x_c·d_c 의존) ----
    print("\n[2/3] clamp 가족 (클램프 케이블만) - x_c·d_c 의존")
    if a.force or not os.path.exists(clamp_tbl):
        odb, _job = expected_odb(a.x_c, disp, a.seed_clamp_dc, 'clamp_only', a.excl)
        if os.path.exists(odb) and not a.force:
            print("      buckle 산출물 있음 - 잡 생략 (%s)" % os.path.basename(odb))
        else:
            _job_t0 = time.time()
            with CasePatch('clamp_only', a.dry_run):
                rc, _ = sh('abaqus cae noGUI=run_abaqus_buckle.py -- %g %g %g %g'
                           % (a.x_c, disp, pull_seed, a.excl), a.dry_run, cwd=HERE)
                if rc:
                    return rc
            odb, _job = expected_odb(a.x_c, disp, a.seed_clamp_dc, 'clamp_only', a.excl)
            if not a.dry_run and not fresh(odb, _job_t0):
                raise SystemExit(
                    "[make_seed] clamp buckle ODB 가 **이번 실행에서 갱신되지 않았습니다** — CAE/해석이 실패했습니다.\n"
                    "  이전 ODB 를 읽어 모드표를 만들면 안 됩니다. buckle/*.msg·.dat 의 원인을 확인하세요.\n"
                    "  odb=%s" % odb)
        _t0 = time.time()
        rc, _ = sh('abaqus python aba_mode_from_odb.py "%s" Step-Buckle "%s" %d'
                   % (odb, clamp_tbl, CLAMP_MODES), a.dry_run, cwd=HERE)
        if not a.dry_run:
            _need(clamp_tbl, 'clamp 가족 모드표',
                  hint=('자식 출력에 "프레임 1개 / 양수 λ 모드가 없다" 가 보이면 시드 클램프 당김이 '
                        '너무 큽니다. 기본 --seed-clamp-dc 0.0(=§19.52 성공 설정)을 쓰세요. 현재: %.4g'
                        % a.seed_clamp_dc))
            if not fresh(clamp_tbl, _t0):
                raise SystemExit(
                    "[make_seed] clamp 모드표가 **이번 실행에서 갱신되지 않았습니다** — buckle 잡이 실패해\n"
                    "  이전 표가 그대로 남아 있습니다. 이 표를 쓰면 1차/낡은 모드가 HF 에 주입됩니다.\n"
                    "  buckle/*.msg·.dat 의 원인을 먼저 고치고 --force 로 다시 실행하세요.\n"
                    "  table=%s" % clamp_tbl)
            check_table_nodes(clamp_tbl, 'clamp 가족 모드표')
    else:
        print("      clamp 표 있음 - 건너뜀 (%s)" % os.path.basename(clamp_tbl))

    # ---- 3) 합성 + 헤더 ----
    print("\n[3/3] 합성 2+2 + 설계점 헤더")
    _fp, _hdr = A.seed_header_lines(a.x_c, a.d_c, A.buckle_fingerprint(BUCKLE_SCRIPT),
                                    seed_tbl, '1,2', clamp_tbl, '1,2')
    if not a.force and os.path.exists(out_tbl):
        _cur = A.read_seed_header(out_tbl)
        if _cur.get('fingerprint') == _fp:
            print("      합성표 있음 - 건너뜀 (%s, 지문 %s)" % (os.path.basename(out_tbl), _fp))
            print("\n[make_seed] OK — HF 가 %s 를 읽는다 (d_c 무관, §19.77)" % os.path.basename(out_tbl))
            return 0
    if a.dry_run:
        _rel_d = lambda p: os.path.relpath(p, HERE)
        print("  $ abaqus python combine_mode_tables.py \"%s\" 1.0:\"%s\":1,2 1.0:\"%s\":1,2"
              % (out_tbl, _rel_d(seed_tbl), _rel_d(clamp_tbl)))
        print("  [header] %s" % (A.seed_header_lines(
            a.x_c, a.d_c, A.buckle_fingerprint(BUCKLE_SCRIPT),
            seed_tbl, '1,2', clamp_tbl, '1,2')[1],))
        print("\n[make_seed] dry-run 종료 - 파일을 만들지 않았다.")
        return 0
    _rel = lambda p: os.path.relpath(p, HERE)          # '1.0:seeds/x.txt:1,2' - 드라이브 콜론 회피
    rc, _ = sh('abaqus python combine_mode_tables.py "%s" 1.0:"%s":1,2 1.0:"%s":1,2'
               % (out_tbl, _rel(seed_tbl), _rel(clamp_tbl)), False, cwd=HERE)
    _need(out_tbl, '합성 시드표')

    fp, hdr = A.seed_header_lines(a.x_c, a.d_c, A.buckle_fingerprint(BUCKLE_SCRIPT),
                                  seed_tbl, '1,2', clamp_tbl, '1,2')
    body = io.open(out_tbl, encoding='utf-8').read()
    #   combine 이 쓴 헤더/데이터 위에 우리 헤더를 **맨 앞**에 얹는다(리더는 'MODE ' 를 만나면 멈춘다).
    io.open(out_tbl, 'w', encoding='utf-8').write('\n'.join(hdr) + '\n' + body)

    # ---- 4) 검증 + meta ----
    ok, why, info = A.seed_gate(a.x_c, a.d_c, out_tbl, BUCKLE_SCRIPT)
    tbl = A.load_mode_table(out_tbl, (1, 2, 3, 4))
    nodes = sorted(set(len(v) for v in tbl.values()))
    _exp_nodes = expected_membrane_nodes()
    if nodes != [_exp_nodes]:
        raise SystemExit(
            "[make_seed] 합성표의 노드 수 %s 가 기대 [%d] 와 다릅니다 -> 이 표를 HF 에 쓰면 안 됩니다: %s"
            % (nodes, _exp_nodes, out_tbl))
    peaks = [round(max(abs(x) for x in v.values()), 6) for _k, v in sorted(tbl.items())]
    meta = ('# 시드표 검증 (make_seed.py %s)\n'
            'x_c=%g d_c=%g disp=%.6g clamp_pull=%.6g excl=%.3g\n'
            'fingerprint=%s\nmodes=%s\nnodes_per_mode=%s\npeak_abs_u3=%s\n'
            'gate=%s (%s)\ntable=%s\n'
            % (time.strftime('%Y-%m-%d %H:%M:%S'), a.x_c, a.d_c, disp, pull, a.excl,
               fp, sorted(tbl.keys()), nodes, peaks, ok, why, out_tbl))
    meta_p = os.path.splitext(out_tbl)[0] + '.meta.txt'
    io.open(meta_p, 'w', encoding='utf-8').write(meta)
    print(meta)
    if not ok:
        print("!!! 게이트 실패 - 표를 만들었지만 HF 가 중단한다: %s" % why)
        return 1
    if nodes != [10201] or sorted(tbl.keys()) != [1, 2, 3, 4] or any(p != 1.0 for p in peaks):
        print("!!! 검증 실패 - nodes=%s modes=%s peaks=%s" % (nodes, sorted(tbl.keys()), peaks))
        return 1
    print("[make_seed] OK — HF 가 %s 를 읽는다" % os.path.basename(out_tbl))
    return 0


if __name__ == '__main__':
    sys.exit(main())
