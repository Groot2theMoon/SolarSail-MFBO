#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""buckle_compare.py - 두 좌굴 런의 lambda + 모드형상(MAC)을 비교해 '모드 소스 교체 필요?' 를 판정한다.

왜 필요한가
    HF(run_abaqus.py)에는 PULL_TOP 스위치가 있어 3-tendon(Top 도 당김)으로 돌 수 있는데,
    좌굴모드 소스(run_abaqus_buckle.py)는 이제 §19.107 로 같은 스위치를 갖는다. 그런데
    "두 base state 가 다르니 모드도 다시 뽑아야 하는가?" 는 lambda 만 봐서는 답이 안 나온다
    (lambda 는 스케일만 변하고 형상은 같은 경우가 흔하다). 필요한 것은 **형상의 상관**이다.
    MAC(Modal Assurance Criterion) 은 그 표준 지표다:
        MAC(a,b) = (a.b)^2 / ((a.a)(b.b))      0(직교) ~ 1(같은 형상, 부호/스케일 무관)

입력
    A, B  : 다음 중 아무거나(자동 판별)
      *.txt  모드표 (aba_mode_from_odb.py 출력: `MODE k` + `<label> <u3>` 행)
      *.dat  Abaqus .dat - 고유값 표(MODE NO / EIGENVALUE|BUCKLING FACTOR)를 읽는다.
             벡터는 보통 .dat 에 없으므로, 같은 파일에 모드표 블록이 있으면 그쪽도 읽는다.
      *.odb  odbAccess 로 Step-Buckle 프레임의 U(u1,u2,u3)를 직접 읽는다
             -> **abaqus python** 으로만 가능(일반 python3 는 exit 3 + 대안 안내).
    n_modes : 비교할 모드 수 (기본 4)
    tol     : 합격 임계 MAC (기본 0.9)

사용
    python3 buckle_compare.py <A> <B> [n_modes] [tol] [json_out]
    abaqus python buckle_compare.py <A.odb> <B.odb> 4          <- ODB 경로
    python3 buckle_compare.py selftest                          <- 해석적 앵커(아래)

출력
    (a) 모드별 lambda 나란히 + 비율 B/A
    (b) MAC 행렬 (행=A 모드, 열=B 모드)
    (c) 최대 MAC 과 그 쌍
    (d) 판정 문구:  MAC>=tol -> 모드 소스 교체 불필요   /   <tol -> 재생성 필요

종료코드
    0  계산 완료 + MAC>=tol (교체 불필요)
    1  selftest 실패
    2  사용법
    3  odbAccess 없음(일반 python3 로 .odb 를 열었다)
    4  모드 벡터를 읽을 수 없음(=MAC 불가)
    5  노드 라벨 교집합 부족
    10 계산 완료 + MAC<tol (재생성 필요)

selftest 앵커 (Abaqus 없이 같은 함수들을 실제 실행)
    (A1) 같은 벡터 두 번            -> MAC = 1.0
    (A2) 직교 벡터                  -> MAC = 0.0
    (A3) 부호 반전                  -> MAC = 1.0
    (A4) 스케일 2.5배               -> MAC = 1.0 (크기 불변)
    (A5) 모드 순서가 뒤바뀐 두 런   -> 최대 MAC 이 (1,2)/(2,1) 에서 1.0 (행렬 인덱스 검증)
"""
from __future__ import print_function

import io
import json
import math
import os
import re
import sys

DEFAULT_N_MODES = 4
DEFAULT_TOL = 0.9
DEFAULT_STEP = 'Step-Buckle'
DEFAULT_INSTANCE = 'MEMBRANE-1'
MIN_COMMON_NODES = 10

TAG = '[MAC]'


# ---------------------------------------------------------------------------
# 1) MAC 핵심 (순수 함수 - selftest 가 직접 호출한다)
# ---------------------------------------------------------------------------
def _split_components(modes):
    """모드 집합이 실제로 3성분을 담고 있는지 판정한다.

    aba_mode_from_odb.py 의 모드표는 u3 만 담는다(u1=u2=0 으로 채워진다).
    그 표에 3성분 MAC 을 걸면 0 으로 나눠 None 이 되므로, 성분을 **데이터로** 정한다.
    """
    for d in modes:
        for v in d.values():
            if abs(v[0]) > 0.0 or abs(v[1]) > 0.0:
                return (0, 1, 2)
    return (2,)


def mac(va, vb, comps=(2,)):
    """MAC(Modal Assurance Criterion). 공통 노드 라벨만 쓴다.

    va, vb = {label: (u1, u2, u3)}
    반환: 0..1 실수 또는 None(공통노드 없음 / 영벡터)
    성질: 스케일 불변, 부호 불변(제곱), 같은 형상이면 1, 직교면 0.
    """
    common = [l for l in va if l in vb]
    if not common:
        return None
    dot = aa = bb = 0.0
    for l in common:
        a, b = va[l], vb[l]
        for c in comps:
            dot += a[c] * b[c]
            aa += a[c] * a[c]
            bb += b[c] * b[c]
    if aa <= 0.0 or bb <= 0.0:
        return None
    return (dot * dot) / (aa * bb)


def mac_matrix(modes_a, modes_b, comps_a, comps_b):
    """행=A, 열=B 의 MAC 행렬. 성분은 **양쪽 공통** 성분만(3D x u3-only -> u3)."""
    comps = comps_a if len(comps_a) <= len(comps_b) else comps_b
    ka = sorted(modes_a)
    kb = sorted(modes_b)
    return ka, kb, [[mac(modes_a[i], modes_b[j], comps) for j in kb] for i in ka], comps


def best_pair(ka, kb, mat):
    """(mac, i, j) 최대값. None 셀은 건너뛴다."""
    best = (-1.0, None, None)
    for r, i in enumerate(ka):
        for c, j in enumerate(kb):
            v = mat[r][c]
            if v is None:
                continue
            if v > best[0]:
                best = (v, i, j)
    return best


# ---------------------------------------------------------------------------
# 2) 입력 읽기
# ---------------------------------------------------------------------------
def _parse_dat_lambdas(text, limit=200):
    """Abaqus .dat 의 고유값 표 -> [lambda, ...] (모드 순서).

    두 헤더를 모두 받는다(실측 .dat 에 둘 다 나온다):
        MODE NO.   EIGENVALUE            (EIGENVALUE OUTPUT)
        MODE NO.   BUCKLING FACTOR       (BUCKLING FACTOR OUTPUT - *BUCKLE 의 표준)
    parse_buckle_log.dat_lambdas 는 EIGEN 헤더만 보므로 여기서는 쓰지 않는다
    (좌굴 .dat 은 BUCKLING FACTOR 라 그 함수로는 0개가 나온다).
    """
    lines = text.splitlines()

    def _is_hdr(s):
        c = re.sub(r'\s+', '', s).upper()
        return c.startswith('MODENO') or c.startswith('MODE') and 'BUCKLINGFACTOR' in c

    hdr = [i for i, l in enumerate(lines) if _is_hdr(l)]
    if not hdr:
        return []
    out = []
    for l in lines[hdr[-1] + 1:]:
        p = l.split()
        if not p:
            if out:
                break
            continue
        if not (p[0].isdigit() and len(p) >= 2):
            break
        try:
            out.append(float(p[-1].replace('D', 'E').replace('d', 'e')))
        except ValueError:
            break
        if len(out) >= limit:
            break
    return out


def _parse_mode_blocks(text):
    """모드표 텍스트 -> ({k: {label: (u1,u2,u3)}}, {k: lambda}).

    지원 형식(aba_mode_from_odb.py 와 같은 표):
        # MODE 1 frame=.. lambda=8.19e-03 ...      <- lambda 메타(선택)
        MODE 1
        1234 1.000000000e+00                       <- u3 만
        또는
        1234 0.0 0.0 1.000000000e+00               <- 3성분
    """
    modes, lams = {}, {}
    cur = None
    for raw in text.splitlines():
        s = raw.strip()
        if not s:
            continue
        if s.startswith('#'):
            m = re.search(r'#\s*MODE\s+(\d+).*?lambda\s*=\s*([-+0-9.eEdD]+)', s)
            if m:
                try:
                    lams[int(m.group(1))] = float(m.group(2).replace('D', 'E').replace('d', 'e'))
                except ValueError:
                    pass
            continue
        m = re.match(r'^MODE\s+(\d+)$', s, re.I)
        if m:
            cur = int(m.group(1))
            modes.setdefault(cur, {})
            continue
        if cur is None:
            continue
        p = s.split()
        if len(p) < 2:
            continue
        try:
            lab = int(p[0])
            vals = [float(x) for x in p[1:]]
        except ValueError:
            continue
        if len(vals) == 1:
            modes[cur][lab] = (0.0, 0.0, vals[0])
        elif len(vals) >= 3:
            modes[cur][lab] = (vals[0], vals[1], vals[2])
    return modes, lams


def _load_odb(path, n_modes, step, instance):
    """ODB -> (modes, lams, msgs). odbAccess 필요(abaqus python)."""
    from odbAccess import openOdb
    msgs = []
    odb = openOdb(path=path, readOnly=True)
    try:
        if step not in odb.steps:
            raise KeyError("스텝 '%s' 없음. 있는 스텝: %s" % (step, sorted(odb.steps.keys())))
        if instance not in odb.rootAssembly.instances:
            raise KeyError("인스턴스 '%s' 없음. 있는 것: %s"
                           % (instance, sorted(odb.rootAssembly.instances.keys())))
        inst = odb.rootAssembly.instances[instance]
        frames = list(odb.steps[step].frames)
        #   좌굴 스텝의 frame 1 은 기저 상태(frameValue 0)이고 모드가 아니다(buckle_mode_report 와 같은 규칙).
        if len(frames) > 1:
            try:
                if float(getattr(frames[0], 'frameValue', -1.0)) == 0.0:
                    frames = frames[1:]
                    msgs.append('기저 프레임(frame 1, frameValue 0) 건너뜀')
            except Exception:
                pass
        modes = {}
        for k, fr in enumerate(frames[:n_modes], start=1):
            try:
                vals = fr.fieldOutputs['U'].getSubset(region=inst).values
            except KeyError:
                raise KeyError("프레임 %d 에 'U' 필드가 없다: %s"
                               % (k, sorted(fr.fieldOutputs.keys())))
            d = {}
            for v in vals:
                d[v.nodeLabel] = (float(v.data[0]), float(v.data[1]), float(v.data[2]))
            if d:
                modes[k] = d
        msgs.append('프레임 %d개 -> 모드 %d개 (U: u1,u2,u3)' % (len(frames), len(modes)))
    finally:
        odb.close()

    lams = []
    dat = os.path.splitext(path)[0] + '.dat'
    if os.path.exists(dat):
        with io.open(dat, encoding='utf-8', errors='replace') as f:
            lams = _parse_dat_lambdas(f.read(), limit=max(n_modes, 12))
        if lams:
            msgs.append('lambda 출처 = %s' % os.path.basename(dat))
    return modes, lams, msgs


def load_source(path, n_modes, step=DEFAULT_STEP, instance=DEFAULT_INSTANCE):
    """(modes, lams, kind, msgs) 를 돌려준다. kind 는 사람이 읽는 '무엇을 읽었나'."""
    if not os.path.exists(path):
        raise IOError('파일 없음: %s' % path)
    ext = os.path.splitext(path)[1].lower()
    msgs = []
    if ext == '.odb':
        try:
            modes, lams, m = _load_odb(path, n_modes, step, instance)
        except ImportError:
            raise OdbAccessMissing(path)
        msgs.extend(m)
        kind = 'ODB %s/U (u1,u2,u3)' % step
    else:
        with io.open(path, encoding='utf-8', errors='replace') as f:
            text = f.read()
        modes, table_lams = _parse_mode_blocks(text)
        dat_lams = [] if ext == '.txt' else _parse_dat_lambdas(text, limit=max(n_modes, 12))
        lams = dat_lams or [table_lams[k] for k in sorted(table_lams)]
        if modes:
            kind = '모드표(%s)' % os.path.basename(path)
        else:
            kind = '.dat 고유값 표만'
        if dat_lams:
            msgs.append('lambda 출처 = .dat 표 (%d개)' % len(dat_lams))
        elif table_lams:
            msgs.append('lambda 출처 = 모드표 헤더 (%d개)' % len(table_lams))

    #   모드 번호 정렬 + n_modes 제한
    keep = sorted(modes)[:n_modes]
    modes = dict((k, modes[k]) for k in keep)
    lams = lams[:n_modes]
    return modes, lams, kind, msgs


class OdbAccessMissing(Exception):
    def __init__(self, path):
        Exception.__init__(self, path)
        self.path = path


# ---------------------------------------------------------------------------
# 3) 보고
# ---------------------------------------------------------------------------
def _fmt(v, spec='%.6e'):
    return '-' if v is None else (spec % v)


def _print_lambdas(name_a, name_b, lams_a, lams_b):
    print('\n=== (a) 고유값 lambda (%s) ===' % '모드 순서')
    print('  모드   %-14s %-14s %-9s' % ('A[%s]' % name_a, 'B[%s]' % name_b, 'B/A'))
    n = max(len(lams_a), len(lams_b), 1)
    for i in range(n):
        la = lams_a[i] if i < len(lams_a) else None
        lb = lams_b[i] if i < len(lams_b) else None
        ratio = (lb / la) if (la not in (None, 0.0) and lb is not None) else None
        print('  %4d   %-14s %-14s %-9s'
              % (i + 1, _fmt(la), _fmt(lb), _fmt(ratio, '%.4f')))


def report(name_a, name_b, modes_a, modes_b, lams_a, lams_b, tol, json_out=None):
    _print_lambdas(name_a, name_b, lams_a, lams_b)

    if not modes_a or not modes_b:
        #   lambda 만 있는 입력(.dat 단독). (a) 는 여전히 유용하므로 찍고, MAC 은 '불가'로 끝낸다.
        print('\n=== (b) MAC ===')
        print('  불가: 벡터가 없는 입력이다 (A 모드 %d / B 모드 %d)' % (len(modes_a), len(modes_b)))
        return 4

    comps_a, comps_b = _split_components(modes_a.values()), _split_components(modes_b.values())
    ka, kb, mat, comps = mac_matrix(modes_a, modes_b, comps_a, comps_b)

    ax = 'u3' if comps == (2,) else 'u1,u2,u3'
    print('\n=== (b) MAC 행렬 (행=A, 열=B, 성분=%s) ===' % ax)
    print('        ' + ''.join('%9s' % ('B%d' % j) for j in kb))
    for r, i in enumerate(ka):
        print('  A%-4d ' % i + ''.join(('   %6.4f' % mat[r][c]) if mat[r][c] is not None
                                       else '        -' for c in range(len(kb))))

    v, bi, bj = best_pair(ka, kb, mat)
    print('\n=== (c) 최대 MAC ===')
    if bi is None:
        print('  없음 (공통 노드가 있는 모드 쌍이 없다)')
    else:
        print('  max MAC = %.6f   (A%d, B%d)' % (v, bi, bj))
        print('  대응: A%d <-> B%d' % (bi, bj))

    print('\n=== (d) 판정 (임계 MAC=%.3f) ===' % tol)
    ok = (bi is not None and v >= tol)
    if ok:
        print('  MAC>=%.1f -> 모드 소스 교체 불필요' % tol)
    else:
        print('  MAC<%.1f -> 재생성 필요' % tol)

    if json_out:
        row = {'A': name_a, 'B': name_b, 'components': ax, 'tol': tol,
               'lambda_a': lams_a, 'lambda_b': lams_b,
               'modes_a': ka, 'modes_b': kb,
               'mac': [[None if c is None else round(c, 12) for c in r] for r in mat],
               'max_mac': (None if bi is None else v),
               'pair': [bi, bj], 'verdict': ('ok' if ok else 'regenerate')}
        with io.open(json_out, 'w', encoding='utf-8') as f:
            f.write(json.dumps(row, ensure_ascii=False) + u'\n')
        print('  JSON: %s' % json_out)
    return 0 if ok else 10


# ---------------------------------------------------------------------------
# 4) selftest - 해석적 앵커
# ---------------------------------------------------------------------------
def _vec(labels_vals, comps=3):
    out = {}
    for i, l in enumerate(sorted(labels_vals)):
        out[l] = tuple(labels_vals[l]) if len(labels_vals[l]) == 3 \
            else (0.0, 0.0, labels_vals[l][0])
    return out


def selftest():
    fails = []
    n = [0]

    def ck(name, cond, detail=None):
        n[0] += 1
        print('  %-4s %s%s' % ('OK' if cond else 'FAIL', name,
                               '' if cond else ('   <- ' + str(detail))))
        if not cond:
            fails.append(name)

    #   앵커용 벡터: 노드 1..6, u3 = [1,-2,3,0.5,-1,2.25]
    u3 = {1: 1.0, 2: -2.0, 3: 3.0, 4: 0.5, 5: -1.0, 6: 2.25}
    v = dict((k, (0.0, 0.0, x)) for k, x in u3.items())
    #   직교: 노드 1..3 에만 값 (u3 = [1,0,0]) vs (u3 = [0,1,0])
    o1 = {1: (0.0, 0.0, 1.0), 2: (0.0, 0.0, 0.0), 3: (0.0, 0.0, 0.0)}
    o2 = {1: (0.0, 0.0, 0.0), 2: (0.0, 0.0, 1.0), 3: (0.0, 0.0, 0.0)}

    print('[selftest] MAC 해석적 앵커')
    m = mac(v, v, (2,))
    ck('A1 같은 벡터 두 번 -> MAC=1.0 (rel<1e-14)', abs(m - 1.0) < 1e-14, m)

    m = mac(o1, o2, (2,))
    ck('A2 직교 벡터 -> MAC=0.0', abs(m) < 1e-14, m)

    vn = dict((k, (t[0], t[1], -t[2])) for k, t in v.items())
    m = mac(v, vn, (2,))
    ck('A3 부호 반전 -> MAC=1.0 (rel<1e-14)', abs(m - 1.0) < 1e-14, m)

    s = dict((k, (t[0], t[1], 2.5 * t[2])) for k, t in v.items())
    m = mac(v, s, (2,))
    ck('A4 스케일 2.5배 -> MAC=1.0 (크기 불변)', abs(m - 1.0) < 1e-14, m)

    #   A5 모드 순서 뒤바뀜 + 3성분 혼합(면내 성분 포함) -> 최대 MAC 이 (1,2)
    wa = {1: (1.0, 0.0, 0.0), 2: (0.0, 1.0, 0.5), 3: (0.0, 0.0, 1.0)}
    wb = {1: (0.0, 2.0, 1.0), 2: (3.0, 0.0, 0.0), 3: (0.0, 0.0, -4.0)}
    A = {1: wa, 2: wb}
    B = {1: wb, 2: wa}
    ka, kb, mat, comps = mac_matrix(A, B, _split_components(A.values()),
                                    _split_components(B.values()))
    best = best_pair(ka, kb, mat)
    ck('A5 뒤바뀐 두 런 -> 3성분 MAC 사용', comps == (0, 1, 2), comps)
    ck('A5 최대 MAC=1.0 이 대각이 아닌 (1,2)/(2,1)',
       abs(best[0] - 1.0) < 1e-14 and ((best[1], best[2]) in ((1, 2), (2, 1))), best)
    ck('A5 반대편 비대각도 1.0 (A2,B1)',
       abs(mat[ka.index(2)][kb.index(1)] - 1.0) < 1e-14, mat)

    #   A6 u3-only 표에 3성분 MAC 을 걸면 안 된다(성분 자동 선택 검증)
    u3a = {1: (0.0, 0.0, 1.0), 2: (0.0, 0.0, 0.5)}
    u3b = {1: (0.0, 0.0, 2.0), 2: (0.0, 0.0, 1.0)}
    ka, kb, mat, comps = mac_matrix({1: u3a}, {1: u3b}, _split_components([u3a]),
                                    _split_components([u3b]))
    ck('A6 u3-only 표 -> 성분 u3 자동 선택', comps == (2,) and abs(mat[0][0] - 1.0) < 1e-14,
       (comps, mat))

    #   A7 라벨 교집합 없음 -> None (조용한 0 이 아니라 '판정 불가')
    ck('A7 공통 라벨 없음 -> MAC=None', mac({1: (0, 0, 1.0)}, {2: (0, 0, 1.0)}) is None)

    print('\nRESULT: %s (%d checks, %d fail)'
          % ('PASS' if not fails else 'FAIL', n[0], len(fails)))
    return 0 if not fails else 1


# ---------------------------------------------------------------------------
# 5) main
# ---------------------------------------------------------------------------
def _usage():
    print('usage: python3 buckle_compare.py <A> <B> [n_modes] [tol] [json_out]')
    print('       python3 buckle_compare.py selftest')
    print('  A/B: .txt(모드표) | .dat(고유값 표) | .odb(abaqus python 전용)')
    return 2


def main(argv):
    args = [a for a in argv[1:] if a.strip() not in ('--',)]
    if not args:
        return _usage()
    if args[0].lower() == 'selftest':
        return selftest()
    if len(args) < 2:
        return _usage()
    path_a, path_b = args[0], args[1]
    try:
        n_modes = int(args[2]) if len(args) > 2 else DEFAULT_N_MODES
        tol = float(args[3]) if len(args) > 3 else DEFAULT_TOL
    except ValueError:
        print('%s 인자 해석 실패: %r (n_modes 는 정수, tol 은 실수)' % (TAG, args[2:4]))
        return 2
    json_out = args[4] if len(args) > 4 else None

    print('%s A = %s' % (TAG, path_a))
    print('%s B = %s' % (TAG, path_b))
    print('%s n_modes=%d tol=%.3f' % (TAG, n_modes, tol))

    try:
        ma, la, kind_a, msg_a = load_source(path_a, n_modes)
        mb, lb, kind_b, msg_b = load_source(path_b, n_modes)
    except OdbAccessMissing as e:
        print('%s ODB 는 odbAccess 가 필요하다 (일반 python3 로는 못 연다): %s' % (TAG, e.path))
        print('      (a) abaqus python buckle_compare.py A.odb B.odb %d' % n_modes)
        print('      (b) abaqus python aba_mode_from_odb.py A.odb %s A.txt %d'
              % (DEFAULT_STEP, n_modes))
        print('          python3 buckle_compare.py A.txt B.txt %d' % n_modes)
        return 3
    except (IOError, OSError) as e:
        print('%s 파일 오류: %s' % (TAG, e))
        return 2

    for m in msg_a:
        print('%s   A| %s' % (TAG, m))
    for m in msg_b:
        print('%s   B| %s' % (TAG, m))
    print('%s 벡터 출처: A = %s / B = %s' % (TAG, kind_a, kind_b))
    print('%s 모드 수: A = %d / B = %d' % (TAG, len(ma), len(mb)))

    if not ma or not mb:
        print('%s 모드 벡터를 읽지 못했다 -> MAC 불가.' % TAG)
        print('      .dat 에는 보통 고유값만 있다: ODB(abaqus python) 또는 모드표(.txt)를 쓰라.')
        print('      abaqus python aba_mode_from_odb.py <odb> %s <out.txt> %d'
              % (DEFAULT_STEP, n_modes))
        return report(os.path.basename(path_a), os.path.basename(path_b), ma, mb, la, lb, tol,
                      json_out)

    common = len(set(next(iter(ma.values()))) & set(next(iter(mb.values()))))
    print('%s 노드 라벨 교집합 = %d (첫 모드 기준)' % (TAG, common))
    if common < MIN_COMMON_NODES:
        print('%s 교집합이 %d 미만 -> 두 런이 같은 메쉬가 아니다(MAC 판정 불가).'
              % (TAG, MIN_COMMON_NODES))
        return 5

    return report(os.path.basename(path_a), os.path.basename(path_b), ma, mb, la, lb, tol,
                  json_out)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
