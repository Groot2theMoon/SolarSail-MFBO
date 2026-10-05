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
    """(판정, corr_same, corr_mirror, maxdiff, rms)"""
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
            return 'DIFFERENT', None, None, float(np.max(np.abs(va - vb))), None, len(labs)
        cs = float(np.corrcoef(va, vb)[0, 1])
        cm = float(np.corrcoef(va, -vb)[0, 1])
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
        md = max(abs(a[l] - b[l]) for l in labs)
        rms = (sum((a[l] - b[l]) ** 2 for l in labs) / n) ** 0.5
    if cs >= 0.999:
        v = 'SAME'
    elif cm >= 0.999:
        v = 'MIRROR'
    elif cs >= 0.99:
        v = 'CLOSE'
    else:
        v = 'DIFFERENT'
    return v, cs, cm, md, rms, len(labs)


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

    v, cs, cm, md, rms, ncom = verdict(a, b)
    print('-' * 78)
    print('[판정]  공통 노드 %d개' % ncom)
    print('  corr(A,  B) = %s' % ('%.6f' % cs if cs is not None else 'n/a'))
    print('  corr(A, -B) = %s      <- 이게 1 에 가까우면 거울상' % ('%.6f' % cm if cm is not None else 'n/a'))
    print('  max|A-B|    = %s m' % ('%.4e' % md if md is not None else 'n/a'))
    print('  RMS(A-B)    = %s m' % ('%.4e' % rms if rms is not None else 'n/a'))
    print('  RESULT:%s' % v)
    if v == 'SAME':
        print('  -> 같은 해다. 수렴 기준/설정이 달라도 동일한 평형해로 갔다.')
    elif v == 'MIRROR':
        print('  -> **거울상 쌍**이다. 크기는 같고 z 방향이 반대다.')
        print('     임퍼펙션이 어느 쪽을 택할지 정하지 못했다는 뜻이다(진폭 과소).')
    elif v == 'CLOSE':
        print('  -> 유사하나 동일하지 않다. 같은 분기의 서로 다른 수렴점일 수 있다.')
    else:
        print('  -> 다른 해다. 수렴 경로/임퍼펙션/감쇠 중 무엇이 결정했는지 따져야 한다.')
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
