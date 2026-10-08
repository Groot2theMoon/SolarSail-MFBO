# -*- coding: utf-8 -*-
"""buckle_mode_report.py — 좌굴모드 **형상**을 읽고 판정한다. (읽기 전용, 라이선스 0)

왜 필요한가
    지금까지 판정은 lambda 의 부호와 개수뿐이었다. 정작 "이 모드가 우리 주름인가" 는 한 번도
    보지 않았다. 실패 유형 두 가지를 형상으로 가른다:
      - 면내 모드(u3 ~ 0): 임퍼펙션은 z 섭동이므로 아무것도 주입하지 못한다.
      - 경계조건 아티팩트: 모드가 사선변(클램프 부착선, u3=0)에 붙어 있고 x_c 를 바꿔도
        움직이지 않는다. 클램프 물리가 아니므로 시드로 쓰면 틀린 초기결함이 된다.
    그래서 x_c 를 바꾼 런들의 모드가 '따라 움직이는가' 를 상관계수와 무게중심 이동으로 본다.

사용 (abaqus python — 위치 인자만. 런처가 --key=value 의 '=' 를 소비한다)
    abaqus python buckle_mode_report.py <odb> <step> [n_modes] [instance]
    abaqus python buckle_mode_report.py cmp <odbA> <odbB> [<odbC>] <step> [n_modes] [instance]
    python3 buckle_mode_report.py selftest          <- Abaqus 없이 순수 함수만 검증

종료코드: 0 정상 / 2 사용법 / 3 odbAccess 없음 / 4 odb/스텝 없음 / 5 selftest 실패
"""
import os
import sys

DEFAULT_INSTANCE = 'MEMBRANE-1'
DEFAULT_THICKNESS = 5.0e-6      # run_abaqus_buckle.py 의 THICKNESS 와 맞춘다

# 판정 임계값 — 실측 분포를 보고 조정한다(계획 §4: 3번 런 결과 뒤에 확정).
CORR_MOVE = 0.80        # 이 이상 상관이면 '같은 모드가 이동한 것' 으로 본다
CORR_SAME = 0.99        # 이 이상이면 사실상 같은 형상
MOVE_FRAC = 0.25        # 무게중심 이동 / 모드 자체 폭(r90) 이 이 이상이면 '이동' (상대 기준)
PINNED_FRAC = 0.05      # 위 비가 이 미만이면 '고정' = 경계조건 아티팩트 의심
MOVE_TOL_M = 0.15       # r90 을 못 구할 때만 쓰는 절대 폴백
U3_FRAC_MIN = 0.10      # 면외 비율이 이 미만이면 '면내 모드' 로 버린다


# ---------------------------------------------------------------------------
# 순수 함수 — 합성 입력으로 Abaqus 없이 검증한다(selftest)
# ---------------------------------------------------------------------------


def mode_metrics(nodes, disp, thickness=DEFAULT_THICKNESS, top_frac=0.1):
    """모드 하나의 형상 지표.

    nodes = {label: (x, y, z)}, disp = {label: (u1, u2, u3)}  (공통 라벨만 쓰인다)

    반환: {'u3_frac', 'cx', 'cy', 'r90', 'peak_xy', 'over_2t', 'n'}
      u3_frac : max|u3| / max(max|u1|, max|u2|, max|u3|)  — 면내 모드면 O(1e-16)
      cx, cy  : |u3| 가중 무게중심 (진폭이 어디 사는가)
      r90     : |u3| 상위 top_frac 노드들의 무게중심으로부터의 RMS 반경 (분포 폭)
      over_2t : |u3| > 2*thickness 인 노드 비율
    """
    common = [l for l in disp if l in nodes]
    if not common:
        return None
    m1 = max(abs(disp[l][0]) for l in common)
    m2 = max(abs(disp[l][1]) for l in common)
    m3 = max(abs(disp[l][2]) for l in common)
    denom = max(m1, m2, m3)
    if denom <= 0.0:
        return None
    s = sum(abs(disp[l][2]) for l in common)
    if s <= 0.0:
        return {'u3_frac': 0.0, 'cx': None, 'cy': None, 'r90': None, 'rw': None,
                'peak_xy': None, 'over_2t': 0.0, 'n': len(common)}
    cx = sum(abs(disp[l][2]) * nodes[l][0] for l in common) / s
    cy = sum(abs(disp[l][2]) * nodes[l][1] for l in common) / s
    # rw = |u3| 가중 RMS 반경(전체 노드) — **모드의 실제 폭** 대리변수.
    #   r90(상위 10% 의 분산 반경)은 피크 영역만 재므로 모드 폭의 대리변수가 되지 못한다
    #   (합성 검증: s=2 -> 16 으로 넓혀도 r90 이 1.12 -> 1.30 으로 거의 불변).
    rw = (sum(abs(disp[l][2]) * ((nodes[l][0] - cx) ** 2 + (nodes[l][1] - cy) ** 2)
              for l in common) / s) ** 0.5
    peak = max(common, key=lambda l: abs(disp[l][2]))
    k = max(1, int(round(len(common) * top_frac)))
    top = sorted(common, key=lambda l: -abs(disp[l][2]))[:k]
    ts = sum(abs(disp[l][2]) for l in top)
    if ts > 0:
        tcx = sum(abs(disp[l][2]) * nodes[l][0] for l in top) / ts
        tcy = sum(abs(disp[l][2]) * nodes[l][1] for l in top) / ts
        r90 = (sum(abs(disp[l][2]) * ((nodes[l][0] - tcx) ** 2 + (nodes[l][1] - tcy) ** 2)
                   for l in top) / ts) ** 0.5
    else:
        r90 = None
    over = sum(1 for l in common if abs(disp[l][2]) > 2.0 * thickness) / float(len(common))
    return {'u3_frac': m3 / denom, 'cx': cx, 'cy': cy, 'r90': r90, 'rw': rw,
            'peak_xy': (nodes[peak][0], nodes[peak][1]), 'over_2t': over, 'n': len(common)}


def corr_u3(disp_a, disp_b, normalize=True):
    """공통 라벨에서 u3 의 피어슨 상관. 규칙 29: 두 해의 동일성은 노드별 상관으로만 판정한다.

    normalize=False 면 원값으로 비교한다(모드별 max|u3|=1 정규화가 되어 있으면 True 로 충분).
    """
    common = [l for l in disp_a if l in disp_b]
    if len(common) < 3:
        return None
    xs = [disp_a[l][2] for l in common]
    ys = [disp_b[l][2] for l in common]
    n = float(len(common))
    mx, my = sum(xs) / n, sum(ys) / n
    vx = sum((v - mx) ** 2 for v in xs)
    vy = sum((v - my) ** 2 for v in ys)
    if vx <= 0.0 or vy <= 0.0:
        return None
    cov = sum((xs[i] - mx) * (ys[i] - my) for i in range(len(common)))
    return cov / ((vx ** 0.5) * (vy ** 0.5))


def verdict_pair(ma, mb, corr, move_frac=MOVE_FRAC, pinned_frac=PINNED_FRAC,
                 corr_move=CORR_MOVE, corr_same=CORR_SAME, move_tol=MOVE_TOL_M):
    """두 런의 같은 순번 모드가 '따라 움직였는가' 판정.

    MOVES      : 상관이 높고 무게중심이 움직였다  -> 클램프 기인 후보(시드로 쓸 수 있다)
    PINNED     : 상관이 매우 높고 무게중심이 고정  -> 경계조건 아티팩트, 시드로 쓰면 안 된다
    ANTI       : 상관이 음수 = 방향이 반대인 다른 분기(규칙 29)
    DIFFERENT  : 그 밖(모드 순서가 바뀌었거나 형상이 다름)
    UNKNOWN    : 비교 불가(무게중심 없음/공통 노드 부족)

    ⚠️ 이동량은 **모드 자체의 폭(r90)에 상대적으로** 본다. 절대 거리로만 보면 국소 모드
    (r90 ~0.3 m)의 작은 이동은 놓치고 전역 모드(r90 ~수 m)의 큰 이동은 과대평가한다.
    실측(합성 검증): 시프트 0.6 m 는 폭 2.0 m 모드에서 d/r90=0.30(이동), 같은 0.6 m 라도
    폭 0.3 m 모드에서는 d/r90=2.0(이동)이지만 폭 8 m 모드에서는 0.075(사실상 고정)이다.
    """
    if corr is None or ma is None or mb is None:
        return 'UNKNOWN'
    if ma.get('cx') is None or mb.get('cx') is None:
        return 'UNKNOWN'
    d = ((ma['cx'] - mb['cx']) ** 2 + (ma['cy'] - mb['cy']) ** 2) ** 0.5
    # 스케일은 **모드의 실제 폭**(rw, |u3| 가중 RMS 반경)을 쓴다. r90 은 피크 영역만 재서
    #   모드 폭의 대리변수가 못 된다(합성 검증). rw 가 없으면 r90, 그것도 없으면 절대 폴백.
    scale = None
    for _m in (ma, mb):
        for _k in ('rw', 'r90'):
            if _m.get(_k):
                scale = max(scale or 0.0, _m[_k])
                break
    d_rel = (d / scale) if scale else None
    if corr <= 0.0:
        return 'ANTI'
    # 절대 폴백(move_tol)은 r90 을 못 구했을 때만 쓴다 — 둘을 OR 로 묶으면 넓은 모드에서
    #   상대 기준이 무력해진다(합성 검증: r90=4 m 인 넓은 모드의 0.6 m 이동이 MOVES 로 오판됐다).
    if corr >= corr_same and ((d_rel < pinned_frac) if d_rel is not None else (d < move_tol)):
        return 'PINNED'
    if corr >= corr_move and ((d_rel >= move_frac) if d_rel is not None else (d >= move_tol)):
        return 'MOVES'
    return 'DIFFERENT'


def judge_modes(metrics, u3_min=U3_FRAC_MIN):
    """런 하나의 모드 목록에서 '시드로 쓸 만한 모드' 를 고른다."""
    ok, drop = [], []
    for i, m in enumerate(metrics):
        if m is None:
            drop.append((i + 1, '지표 계산 불가'))
        elif m['u3_frac'] < u3_min:
            drop.append((i + 1, 'u3_frac=%.2e < %.2g (면내 모드)' % (m['u3_frac'], u3_min)))
        else:
            ok.append(i + 1)
    return ok, drop


# ---------------------------------------------------------------------------
# ODB 읽기 (Abaqus 전용 — selftest 는 이 함수를 부르지 않는다)
# ---------------------------------------------------------------------------


def read_odb_modes(odb_path, step_name, n_modes=8, instance=DEFAULT_INSTANCE, dat_hint=None):
    """(nodes, disps, lams, msgs) 를 돌려준다.

    nodes = {label: (x,y,z)}, disps = [ {label: (u1,u2,u3)} ] (프레임 순), lams = [lambda]
    lambda 는 ODB frameValue(=모드 번호)가 아니라 .dat 의 MODE NO 표에서 온다(규칙 10).
    """
    msgs = []
    if not os.path.exists(odb_path):
        raise IOError('ODB 없음: %s' % odb_path)
    from odbAccess import openOdb
    odb = openOdb(path=odb_path, readOnly=True)
    try:
        if step_name not in odb.steps:
            raise KeyError("스텝 '%s' 없음. 있는 스텝: %s"
                           % (step_name, sorted(odb.steps.keys())))
        inst = (odb.rootAssembly.instances[instance]
                if instance in odb.rootAssembly.instances else None)
        if inst is None:
            raise KeyError("인스턴스 '%s' 없음. 있는 것: %s"
                           % (instance, sorted(odb.rootAssembly.instances.keys())))
        nodes = {}
        for nd in inst.nodes:
            nodes[nd.label] = (float(nd.coordinates[0]), float(nd.coordinates[1]),
                               float(nd.coordinates[2]) if len(nd.coordinates) > 2 else 0.0)
        msgs.append('[REPORT] %s / %s : 노드 %d개 (인스턴스 %s)'
                    % (os.path.basename(odb_path), step_name, len(nodes), instance))
        disps = []
        _frames = list(odb.steps[step_name].frames)
        for _fi, fr in enumerate(_frames):
            #   [2026-10-07] 좌굴 스텝 ODB 의 **frame 1 은 기저 상태**(frameValue 0)이고 모드가 아니다.
            #   포함하면 세 가지가 망가진다:
            #     (a) M1 이 '평탄한 기저'로 잡혀 u3_frac~0 -> 면내 모드처럼 보인다(거짓 판정)
            #     (b) .dat 의 lambda 목록과 **한 칸씩 어긋난다**
            #     (c) 요청 n개 중 실제 모드는 n-1개만 보인다(모드 하나를 통째로 놓친다)
            #   frameValue==0 은 좌굴 스텝에서만 기저를 뜻한다(주파수 스텝은 주파수라 0 이 아니다).
            if _fi == 0 and len(_frames) > 1:
                try:
                    if float(getattr(fr, 'frameValue', -1.0)) == 0.0:
                        msgs.append('[REPORT] 기저 프레임(frame 1, frameValue 0) 건너뜀 — 모드가 아니다')
                        continue
                except Exception:
                    pass
            d = {}
            try:
                for v in fr.fieldOutputs['U'].values:
                    d[v.nodeLabel] = (float(v.data[0]), float(v.data[1]), float(v.data[2]))
            except KeyError:
                msgs.append('[REPORT] 프레임에 U 출력이 없다 — FieldOutputRequest 확인')
                break
            disps.append(d)
        lams = []
        if dat_hint and os.path.exists(dat_hint):
            try:
                import parse_buckle_log as _pbl
                import io as _io
                lams = _pbl.dat_lambdas(_io.open(dat_hint, encoding='utf-8',
                                                 errors='replace').read(), limit=max(n_modes, 12))
            except Exception as _e:
                msgs.append('[REPORT] .dat lambda 파싱 실패(%s) -> frameValue 사용' % _e)
        if not lams:
            lams = [float(getattr(fr, 'frameValue', 0.0)) for fr in odb.steps[step_name].frames]
            msgs.append('[REPORT] 경고: lambda 출처 = ODB frameValue(=모드 번호). .dat 를 주면 정확해진다')
        return nodes, disps, lams, msgs
    finally:
        try:
            odb.close()
        except Exception:
            pass


def fmt_mode_line(i, lam, met, label=None, note_if_none='(지표 계산 불가)'):
    """모드 한 줄 포맷. **순수 함수** — 면내 모드(met['cx'] is None)에서도 크래시하지 않는다.

    실측 크래시(2026-10-07, control_none M1): 그 모드는 u3 가 **항등적으로 0**(면내 모드)이라
    mode_metrics 가 cx=None 을 돌려주는데, 포맷이 '%7.3f' % None 을 시도해
    `TypeError: must be real number, not NoneType` 로 죽었다. 포맷을 순수 함수로 빼서
    그 입력을 로컬 단위검증에 넣는다.
    """
    _lab = label if label is not None else 'M%d(fr%d)' % (i, i)
    _lam = '%.5e' % lam if lam is not None else '?'
    if met is None:
        return '  %-10s lambda=%-13s %s' % (_lab, _lam, note_if_none)
    _c = ('(%7.3f,%7.3f)' % (met['cx'], met['cy'])
          if met.get('cx') is not None else '(   면내   )')
    _rw = ('%.3f' % met['rw']) if met.get('rw') is not None else '  -  '
    _r90 = ('%.3f' % met['r90']) if met.get('r90') is not None else '  -  '
    return ('  %-10s lambda=%-13s u3_frac=%-8.3g centroid=%s rw=%-7s r90=%-7s >2t=%5.1f%%'
            % (_lab, _lam, met['u3_frac'], _c, _rw, _r90, 100.0 * met['over_2t']))


def _report_one(odb_path, step_name, n_modes, instance, dat_hint, limit):
    nodes, disps, lams, msgs = read_odb_modes(odb_path, step_name, n_modes, instance, dat_hint)
    for m in msgs:
        print(m)
    metric_list, label_list = [], []
    for i in range(min(n_modes, len(disps))):
        # 모드별 정규화: max|u3| = 1 (형상 비교는 스케일 무관)
        d = disps[i]
        m3 = max([abs(v[2]) for v in d.values()] or [0.0])
        if m3 > 0:
            d = dict((k, (v[0] / m3, v[1] / m3, v[2] / m3)) for k, v in d.items())
        met = mode_metrics(nodes, d)
        lam = lams[i] if i < len(lams) else None
        label_list.append('M%d(fr%d)' % (i + 1, i + 1))
        metric_list.append(met)
        print(fmt_mode_line(i + 1, lam, met, label_list[-1]))
    ok, drop = judge_modes(metric_list)
    print('  -> 시드 후보 모드: %s' % (ok or '없음'))
    for i, why in drop:
        print('     버림 M%d: %s' % (i, why))
    return nodes, disps, metric_list, label_list


def cmp_verdicts(da, db, ma, mb):
    """두 ODB 의 모드 대응 판정 목록. **순수 함수** — 프레임 수와 메트릭 수가 달라도 안전하다.

    실측 크래시(2026-10-07): ODB 프레임 수가 메트릭 수보다 많아(요청 4 인데 프레임이 더 있다)
    range(min(len(da), len(db))) 로 돌다 ma[k] 에서 IndexError.
    => 루프 상한은 **메트릭 수**로 잡는다. 반환 = [(모드번호, corr, 판정)].
    """
    out = []
    for k in range(min(len(ma), len(mb))):
        c = corr_u3(da[k], db[k]) if (k < len(da) and k < len(db)) else None
        out.append((k + 1, c, verdict_pair(ma[k], mb[k], c)))
    return out


def _interp(v):
    """문자열 -> (kind, value). 'p2 xc000' 같은 태그에서 x_c 를 뽑는다."""
    import re
    m = re.search(r'_xc(\d{3})', v)
    if m:
        return ('x_c', int(m.group(1)) / 100.0)
    return ('raw', v)


def _cmp(argv, step_name, n_modes, instance):
    """여러 ODB 의 같은 순번 모드를 비교해 MOVES/PINNED 를 판정한다."""
    outs = []
    for p in argv:
        print('--- %s ---' % os.path.basename(p))
        # .dat 을 물려 lambda 를 ODB frameValue(=모드 번호) 대신 실제 고유값으로 읽는다.
        nodes, disps, metrics, labels = _report_one(
            p, step_name, n_modes, instance, os.path.splitext(p)[0] + '.dat', 12)
        outs.append((p, nodes, disps, metrics))
        print()
    bad = False
    for i in range(len(outs) - 1):
        pa, na, da, ma = outs[i]
        pb, nb, db, mb = outs[i + 1]
        print('=== %s  vs  %s ===' % (os.path.basename(pa), os.path.basename(pb)))
        for k, c, v in cmp_verdicts(da, db, ma, mb):
            print('  M%-3d corr(u3)=%-8s  -> %s' % (k, ('%.4f' % c) if c is not None else '  -  ', v))
            if v == 'PINNED':
                bad = True
        print()
    if bad:
        print('판정: PINNED 가 있다 — 그 모드는 경계조건 아티팩트일 수 있다. 시드로 쓰기 전에')
        print('      압축영역(클램프 부착선 밖)에 진폭이 있는지 확인하라.')
    else:
        print('판정: PINNED 없음 — 상관/이동으로 판단한 결과는 위 표를 따른다.')
    return 0


def selftest():
    """합성 모드로 순수 함수를 검증한다(부정 포함). Abaqus 불필요.

    픽스처는 **국소 모드**(Gaussian)를 쓴다 — 실제 좌굴모드에 가깝다. 초기 픽스처였던 원뿔(cone)은
    10x10 전역에서 상관이 0.32 밖에 안 나와 MOVES 판정을 시험하지 못했다(그 실패가 이 픽스처를 고쳤다).
    """
    import math
    nodes, base = {}, {}
    for i in range(21):
        for j in range(21):
            l = i * 100 + j
            x, y = i * 0.5, j * 0.5                       # 0..10 m, 간격 0.5 m
            nodes[l] = (x, y, 0.0)
            base[l] = (0.0, 0.0, math.exp(-((x - 5.0) ** 2 + (y - 5.0) ** 2) / 2.0))

    def bump(cx, cy, s=2.0):
        return dict((l, (0.0, 0.0,
                         math.exp(-((nodes[l][0] - cx) ** 2 + (nodes[l][1] - cy) ** 2) / s)))
                    for l in nodes)

    moved = bump(5.6, 5.0)            # 같은 모양, 0.6 m 이동
    wide0 = bump(5.0, 5.0, s=8.0)     # 훨씬 넓은 모드
    wide = bump(5.6, 5.0, s=8.0)      # 같은 0.6 m 이동이지만 상대 이동량은 작다
    anti = dict((l, (0.0, 0.0, -base[l][2])) for l in base)
    inplane = dict((l, (base[l][2], 0.0, 1e-16)) for l in base)
    #  u3 가 **항등적으로 0** 인 픽스처 — 실측 크래시(control_none M1)의 실제 입력이다.
    #  1e-16 이면 early-return 경로를 타지 않으므로(w 합 > 0) 크래시 회귀를 시험하지 못한다.
    inplane0 = dict((l, (base[l][2], 0.0, 0.0)) for l in base)

    mb = mode_metrics(nodes, base)
    mm = mode_metrics(nodes, moved)
    mw0 = mode_metrics(nodes, wide0)
    mw = mode_metrics(nodes, wide)
    mn = mode_metrics(nodes, inplane)
    mn0 = mode_metrics(nodes, inplane0)
    c_move = corr_u3(base, moved)
    c_anti = corr_u3(base, anti)
    c_wide = corr_u3(wide0, wide)
    v_move = verdict_pair(mb, mm, c_move)
    v_wide = verdict_pair(mw0, mw, c_wide)
    checks = [
        ('면외 비율 O(1)', mb['u3_frac'] > 0.9),
        ('무게중심 = (5,5)', abs(mb['cx'] - 5.0) < 0.2 and abs(mb['cy'] - 5.0) < 0.2),
        ('r90 이 모드 폭 규모(~1-2 m)', mb['r90'] is not None and 0.5 < mb['r90'] < 3.0),
        ('over_2t 비율 유한', 0.0 <= mb['over_2t'] <= 1.0),
        ('이동한 모드 -> MOVES', v_move == 'MOVES'),
        ('같은 모드 -> PINNED', verdict_pair(mb, mb, corr_u3(base, base)) == 'PINNED'),
        ('반전 -> ANTI', verdict_pair(mb, mb, c_anti) == 'ANTI'),
        ('반전 상관 = -1', c_anti is not None and c_anti < -0.999),
        ('면내 모드 -> u3_frac < 0.1', mn['u3_frac'] < U3_FRAC_MIN),
        ('면내 모드가 judge 에서 버려짐', judge_modes([mn])[0] == []),
        ('정상 모드가 judge 를 통과', judge_modes([mb])[0] == [1]),
        ('스케일 민감성: 넓은 모드의 같은 0.6 m 이동은 MOVES 가 아니다', v_wide != 'MOVES'),
        #  포맷 회귀 — 실측 크래시(면내 모드 = u3 항등 0 -> cx=None)를 그대로 입력으로 넣는다
        ('면내 픽스처가 early-return 경로를 탄다(cx=None)', mn0 is not None and mn0['cx'] is None),
        ('크래시 회귀: 면내 모드(cx=None) 포맷', '면내' in fmt_mode_line(1, 1.5e-4, mn0)),
        ('크래시 회귀: met=None 포맷', '계산 불가' in fmt_mode_line(1, None, None)),
        ('정상 모드 포맷에 lambda 가 들어감', '1.48328e-04' in fmt_mode_line(1, 1.48328e-4, mb)),
        #  cmp 회귀 — 프레임 수가 메트릭 수보다 많은 실측 상황
        ('cmp 상한: 프레임 수 > 메트릭 수여도 안전',
         len(cmp_verdicts([base] * 6, [moved] * 6, [mb, mb], [mm, mm])) == 2),
        ('cmp 판정: 같은 모드쌍은 PINNED',
         cmp_verdicts([base], [base], [mb], [mb])[0][2] == 'PINNED'),
    ]
    nfail = 0
    print('  [합성] shift=0.6 m, r90=%.3f -> corr=%.4f, d=%.3f, d/r90=%.3f -> %s'
          % (mb['r90'], c_move,
             ((mb['cx'] - mm['cx']) ** 2 + (mb['cy'] - mm['cy']) ** 2) ** 0.5,
             ((mb['cx'] - mm['cx']) ** 2 + (mb['cy'] - mm['cy']) ** 2) ** 0.5 / mb['r90'],
             v_move))
    for name, ok in checks:
        print('  %-4s %s' % ('OK' if ok else 'FAIL', name))
        nfail += 0 if ok else 1
    print('selftest: %d/%d 통과' % (len(checks) - nfail, len(checks)))
    return 0 if nfail == 0 else 5


def main(argv):
    if len(argv) < 2:
        sys.stderr.write(__doc__)
        return 2
    if argv[1] == 'selftest':
        return selftest()
    try:
        from odbAccess import openOdb           # noqa: F401
    except ImportError:
        sys.stderr.write('odbAccess 없음 — `abaqus python` 으로 실행하거나 selftest 를 쓰세요.\n')
        return 3
    args = list(argv[1:])
    if args and args[0] == 'cmp':
        args = args[1:]
        if len(args) < 3:
            sys.stderr.write('usage: cmp <odbA> <odbB> [<odbC>] <step> [n_modes] [instance]\n')
            return 2
        # 뒤에서부터 step / n_modes / instance 를 집는다(위치 인자만)
        inst = args[-1] if not args[-1].replace('.', '').lstrip('-').isdigit() \
            and not args[-1].startswith('Step') else DEFAULT_INSTANCE
        if inst is not DEFAULT_INSTANCE:
            args = args[:-1]
        nmodes = 8
        if args[-1].isdigit():
            nmodes = int(args[-1])
            args = args[:-1]
        step = args[-1]
        odbs = args[:-1]
        return _cmp(odbs, step, nmodes, inst)
    if len(args) < 2:
        sys.stderr.write(__doc__)
        return 2
    odb, step = args[0], args[1]
    nmodes = int(args[2]) if len(args) > 2 else 8
    inst = args[3] if len(args) > 3 else DEFAULT_INSTANCE
    dat_hint = os.path.splitext(odb)[0] + '.dat'
    try:
        _report_one(odb, step, nmodes, inst, dat_hint, 12)
    except (IOError, KeyError) as e:
        sys.stderr.write('FAIL %s\n' % e)
        return 4
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
