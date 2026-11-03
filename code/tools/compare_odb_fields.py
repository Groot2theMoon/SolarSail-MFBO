# -*- coding: utf-8 -*-
"""두 ODB 의 변위장을 노드별로 비교한다 — '같은 해인가 / 거울상인가 / 다른 해인가'.

사용:
    abaqus python compare_odb_fields.py <A.odb> <B.odb> <step_name> [frameA] [frameB]
    (frame 생략 시 각 스텝의 마지막 프레임. instance 는 MEMBRANE-1)

왜 필요한가:
    포스트버클링은 다중 평형해를 갖는다. 수렴 기준(Rαn/Cαn)이나 임퍼펙션 진폭에 따라
    **다른 국소해**로 갈 수 있는데, 스칼라 지표(면적비/응력/RF)만 보면 두 해가 '비슷해
    보여서' 같은 해인지 판정할 수 없다. 노드별 상관계수로 판정한다.

판정:
    SAME       corr(A, B) >= 0.999          같은 해
    MIRROR     corr(A, -B) >= 0.999         거울상 쌍 (z 부호 반전)
    CLOSE      corr(A, B) >= 0.99           유사하나 동일하지 않음
    DIFFERENT  그 외
"""
import sys

INSTANCE = 'MEMBRANE-1'
THK = 5.0e-06


def read_u3(odb_path, step_name, frame_idx=None):
    from odbAccess import openOdb
    odb = openOdb(odb_path, readOnly=True)
    try:
        if step_name not in odb.steps.keys():
            raise SystemExit('ERROR: %s 에 스텝 %r 이 없습니다. 있는 스텝: %s'
                             % (odb_path, step_name, list(odb.steps.keys())))
        st = odb.steps[step_name]
        nf = len(st.frames)
        if nf == 0:
            raise SystemExit('ERROR: %s / %s 에 프레임이 없습니다' % (odb_path, step_name))
        fi = nf - 1 if frame_idx is None else int(frame_idx)
        fr = st.frames[fi]
        out = {}
        for v in fr.fieldOutputs['U'].values:
            inst = getattr(v, 'instance', None)
            if inst is None or inst.name != INSTANCE:
                continue
            d = v.data
            out[int(v.nodeLabel)] = float(d[2])      # u3
        return out, fi + 1, nf, float(fr.frameValue)
    finally:
        try:
            odb.close()
        except Exception:
            pass


def stats(d):
    vals = list(d.values())
    mx = max(vals); mn = min(vals)
    return dict(n=len(vals), uzmax=mx, uzmin=mn, uzave=(abs(mx) + abs(mn)) / 2.0,
                mean_abs=sum(abs(v) for v in vals) / len(vals))


def verdict(a, b):
    """(판정, corr_same, corr_mirror, corr_mag, maxdiff, rms, ncommon)"""
    try:
        import numpy as np
    except ImportError:
        np = None
    labs = sorted(set(a) & set(b))
    if not labs:
        return 'DIFFERENT', None, None, None, None, 0
    if np is not None:
        va = np.array([a[l] for l in labs]); vb = np.array([b[l] for l in labs])
        if va.std() < 1e-30 or vb.std() < 1e-30:
            return 'DIFFERENT', None, None, None, float(np.max(np.abs(va - vb))), None, len(labs)
        cs = float(np.corrcoef(va, vb)[0, 1])
        cm = float(np.corrcoef(va, -vb)[0, 1])
        # |A| vs |B|: 부호를 지운 크기 패턴이 닮았는가. 거울상/반상관 판정의 보조 지표.
        cg = float(np.corrcoef(np.abs(va), np.abs(vb))[0, 1])
        md = float(np.max(np.abs(va - vb)))
        rms = float(np.sqrt(np.mean((va - vb) ** 2)))
    else:
        n = len(labs)
        ma = sum(a[l] for l in labs) / n; mb = sum(b[l] for l in labs) / n
        sa = sum((a[l] - ma) ** 2 for l in labs) ** 0.5
        sb = sum((b[l] - mb) ** 2 for l in labs) ** 0.5
        if sa < 1e-30 or sb < 1e-30:
            return 'DIFFERENT', None, None, None, None, n
        cs = sum((a[l] - ma) * (b[l] - mb) for l in labs) / (sa * sb)
        cm = sum((a[l] - ma) * (-b[l] - (-mb)) for l in labs) / (sa * sb)
        aa = [abs(a[l]) for l in labs]; bb = [abs(b[l]) for l in labs]
        m2a = sum(aa) / n; m2b = sum(bb) / n
        s2a = sum((x - m2a) ** 2 for x in aa) ** 0.5
        s2b = sum((x - m2b) ** 2 for x in bb) ** 0.5
        cg = (sum((x - m2a) * (y - m2b) for x, y in zip(aa, bb)) / (s2a * s2b)
              if s2a > 1e-30 and s2b > 1e-30 else None)
        md = max(abs(a[l] - b[l]) for l in labs)
        rms = (sum((a[l] - b[l]) ** 2 for l in labs) / n) ** 0.5
    if cs >= 0.999:
        v = 'SAME'
    elif cm >= 0.999:
        v = 'MIRROR'
    elif cs >= 0.99:
        v = 'CLOSE'
    elif cs <= -0.5:
        # [2026-10-05] 실측: HF_x025r vs HF_x025i 가 corr=-0.802, corr(A,-B)=+0.802 였다.
        #   즉 '거울상'도 '무상관'도 아니다. 전역 스냅 방향이 반대인 두 분기이면서
        #   국소 패턴은 완전 대칭이 아닌(클램프가 대칭을 깬다) 상태다.
        #   => 이 등급을 따로 둔다. 크기 지표만 보면 '비슷해' 보이는 게 함정이다.
        v = 'ANTI'
    else:
        v = 'DIFFERENT'
    return v, cs, cm, cg, md, rms, len(labs)


def main():
    if len(sys.argv) < 4:
        print('사용: abaqus python compare_odb_fields.py <A.odb> <B.odb> <step> [fA] [fB]')
        return 2
    pa, pb, step = sys.argv[1], sys.argv[2], sys.argv[3]
    fa = int(sys.argv[4]) if len(sys.argv) > 4 else None
    fb = int(sys.argv[5]) if len(sys.argv) > 5 else None

    a, ia, na, ta = read_u3(pa, step, fa)
    b, ib, nb, tb = read_u3(pb, step, fb)
    sa, sb_ = stats(a), stats(b)

    print('=' * 78)
    print('두 ODB 변위장 비교   instance=%s   step=%s' % (INSTANCE, step))
    print('=' * 78)
    for tag, path, st, fi, nf, tv in (('A', pa, sa, ia, na, ta), ('B', pb, sb_, ib, nb, tb)):
        print('[%s] %s' % (tag, path))
        print('     프레임 %d/%d  frameValue=%.6g  노드 %d개' % (fi, nf, tv, st['n']))
        print('     u_z,max=%+.6e  u_z,min=%+.6e  u_z,ave=%.6e  mean|u3|=%.6e'
              % (st['uzmax'], st['uzmin'], st['uzave'], st['mean_abs']))
        print('     t 정규화(%.1e m): max %.1f t  min %.1f t  ave %.1f t  mean|u3| %.1f t'
              % (THK, st['uzmax'] / THK, st['uzmin'] / THK, st['uzave'] / THK, st['mean_abs'] / THK))

    v, cs, cm, cg, md, rms, ncom = verdict(a, b)
    print('-' * 78)
    print('[판정]  공통 노드 %d개' % ncom)
    print('  corr(A,  B) = %s' % ('%.6f' % cs if cs is not None else 'n/a'))
    print('  corr(A, -B) = %s      <- +1 에 가까우면 완전 거울상, -1 에 가까우면 A 와 동일' % ('%.6f' % cm if cm is not None else 'n/a'))
    print('  corr(|A|,|B|) = %s    <- 부호를 지운 크기 패턴이 닮았는가' % ('%.6f' % cg if cg is not None else 'n/a'))
    print('  max|A-B|    = %s m' % ('%.4e' % md if md is not None else 'n/a'))
    print('  RMS(A-B)    = %s m   (mean|u3|: A %.4e / B %.4e)'
          % ('%.4e' % rms if rms is not None else 'n/a', sa['mean_abs'], sb_['mean_abs']))
    print('  RESULT:%s' % v)
    if v == 'SAME':
        print('  -> 같은 해다. 수렴 기준/설정이 달라도 동일한 평형해로 갔다.')
    elif v == 'MIRROR':
        print('  -> **거울상 쌍**이다. 크기는 같고 z 방향이 반대다.')
        print('     임퍼펙션이 어느 쪽을 택할지 정하지 못했다는 뜻이다(진폭 과소).')
    elif v == 'CLOSE':
        print('  -> 유사하나 동일하지 않다. 같은 분기의 서로 다른 수렴점일 수 있다.')
    elif v == 'ANTI':
        print('  -> **반상관**이다. 전역 스냅 방향이 반대인 두 분기로 갈렸고,')
        print('     클램프가 대칭을 깨므로 완전 거울상(-1)까지는 아니다.')
        print('     크기 지표(면적비/응력/최대변위)만 보면 비슷해 보이는 것이 함정이다.')
    else:
        print('  -> 다른 해다. 수렴 경로/임퍼펙션/감쇠 중 무엇이 결정했는지 따져야 한다.')
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
