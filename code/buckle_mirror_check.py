"""
buckle_mirror_check.py — 좌굴 ODB의 미러(x -> 20-x) 대칭 진단.

무엇을 하나
    좌굴 잡 ODB의 기저 프레임(Step-Buckle frame 0)과 모드 프레임들의 변위장을 읽어,
    x -> 20-x 미러쌍에 대한 대칭 잔차를 계산한다.
    미러 변환: u1' = -u1, u2' = +u2, u3' = +u3 (벡터의 x성분만 뒤집힌다).

판정 원리 (군론)
    기저가 미러대칭이면 축퇴쌍(λ1≈λ2)의 고유공간은 대칭부 ⊕ 반대칭부로 분해되고,
    그 안의 **정규직교 쌍**은 (대칭±반대칭)/√2 의 회전일 뿐이다 — 즉 둘 다 같은 쪽에
    몰리는 것은 불가능하고, 반드시 한 쪽은 반대편이다(둘 다 한쪽 = 서로의 거울상).
    따라서 '둘 다 한쪽'이 관측되면 기저가 미러대칭을 깬 것이다(개시 넘음 -> 국소좌굴).

보조 검사
    모드1 ± 모드2 의 미러 잔차가 작으면 기저 선택 문제(대칭 기저), 크면 실제 비대칭이다.

사용법 (Abaqus 머신, code 디렉터리, job = cwd 상대 경로·확장자 없이)
    abaqus cae noGUI=buckle_mirror_check.py -- <job> [step]
    예) abaqus cae noGUI=buckle_mirror_check.py -- buckle/Buckle_xc020_d018um_dc006_f020_s4_grid_clamp_lf_normal_ex210 Step-Buckle
    ⚠️ [§19.53] CAE 커널은 argv 앞에 exe 경로를 넣는다 -> 인자는 **꼬리 기준**으로 읽는다
       (앞에서 a[0] 을 쓰던 구현이 `ABQcaeK.exe.odb` 로 죽었다).
    ⚠️ buckle 산출물은 RUN_DIR_NAME='buckle' 아래에 쌓이므로 buckle/ 접두사가 필요하다
       (openOdb(job + '.odb') 를 cwd 상대경로로 연다 — 접두사가 없으면 FileNotFound 로 죽는다).

출력 판독
    기저 잔차/크기 ~ 0   -> 기저 대칭 (모드 쌍은 반드시 좌우 짝; 보이는 쪽 문제는 기저 회전)
    기저 잔차/크기 큼    -> 기저가 대칭을 깼다 (개시 위) — 씨드 재료로는 부적합 판정 근거
"""
from odbAccess import openOdb
import sys

INST = 'MEMBRANE-1'
AXIS = 10.0     # 미러축 x=10 (BASE=20 의 중앙)
TOL = 1e-6


def _is_kernel_token(t):
    """Abaqus CAE 커널이 argv 앞에 끼워 넣는 토큰인가."""
    tl = t.lower()
    return (tl.endswith('.py') or tl.endswith('.exe') or 'nogui=' in tl)


def _args(argv=None):
    """CLI 인자를 꼬리 기준으로 뽑는다 -> (job, step).

    ⚠️ [2026-10-12 §19.53] `abaqus cae noGUI=script.py -- <job> <step>` 로 부르면 sys.argv 앞쪽에
    **커널 실행파일 경로**(.../ABQcaeK.exe)와 `noGUI=...py` 가 들어온다. 예전 구현은 앞에서부터
    a[0] 을 job 으로 써서 exe 경로를 열려다 `OdbError: ...ABQcaeK.exe.odb` 로 죽었다.
    (반면 `abaqus python script.py ...` 는 argv[0]='script.py' 라 무사하다 — 두 호출 방식이 다르다.)
    -> 커널 토큰을 걸러내고 **꼬리에서** 읽는다(run_abaqus_buckle.py 의 sys.argv[-3:] 관례와 동일).
    """
    a = [t for t in list(sys.argv if argv is None else argv)
         if t not in ('--',) and not _is_kernel_token(t)]
    if not a:
        return None, None
    if len(a) >= 2 and a[-1].lower().startswith('step'):
        return a[-2], a[-1]
    return a[-1], 'Step-Buckle'


def collect(fr):
    """프레임에서 (라벨->(u1,u2,u3), field max) 를 모은다. 변위가 없으면 None."""
    try:
        u = fr.fieldOutputs['U']
    except Exception:
        return None
    vals, mx = {}, 0.0
    for v in u.values:
        inst = getattr(v, 'instance', None)
        if inst is None or inst.name != INST:
            continue
        vals[v.nodeLabel] = v.data
        for c in v.data:
            mx = max(mx, abs(c))
    return (vals, mx) if vals else None


def mirror_pairs(inst):
    """(라벨, 미러 라벨) 쌍 목록. x < 10 만 사용(축 위 노드는 자기 자신이라 제외)."""
    coord = {}
    for n in inst.nodes:
        coord[(round(n.coordinates[0], 6), round(n.coordinates[1], 6))] = n.label
    pairs = []
    for (x, y), lab in coord.items():
        if x < AXIS - TOL:
            m = coord.get((round(20.0 - x, 6), round(y, 6)))
            if m is not None:
                pairs.append((lab, m))
    return pairs


def resid(vals, pairs):
    """미러 잔차 (max|du1|, max|du2|, max|du3|, 쌍 수)."""
    r = [0.0, 0.0, 0.0]
    n = 0
    for a, b in pairs:
        if a in vals and b in vals:
            va, vb = vals[a], vals[b]
            d = (va[0] + vb[0], va[1] - vb[1], va[2] - vb[2])
            for i in range(3):
                r[i] = max(r[i], abs(d[i]))
            n += 1
    return r[0], r[1], r[2], n


def main():
    job, step = _args()
    if job is None:
        print('usage: abaqus cae noGUI=buckle_mirror_check.py -- <job> [step]')
        return 1
    odb = openOdb(job + '.odb', readOnly=True)
    inst = odb.rootAssembly.instances[INST]
    pairs = mirror_pairs(inst)
    print('[mirror] 미러쌍 %d개 (x < %g 제외)  job=%s  step=%s' % (len(pairs), AXIS, job, step))
    frames = odb.steps[step].frames
    datas = []
    for i, f in enumerate(frames):
        got = collect(f)
        if not got:
            continue
        datas.append((i, got))
        r1, r2, r3, n = resid(got[0], pairs)
        scale = got[1] if got[1] > 0 else 1.0
        tag = '기저 ' if i == 0 else '모드%d' % i
        print('  %s (frame %d): max|du1|=%.3e  max|du2|=%.3e  max|du3|=%.3e  |  field max=%.3e  잔차/크기=%.2f'
              % (tag, i, r1, r2, r3, got[1], max(r1, r2, r3) / scale))
    modes = [g for (i, g) in datas if i >= 1]
    if len(modes) >= 2:
        v1, v2 = modes[0][0], modes[1][0]
        for label, sgn in (('합(+)', 1.0), ('차(-)', -1.0)):
            comb = {}
            for k in v1:
                if k in v2:
                    comb[k] = tuple(v1[k][j] + sgn * v2[k][j] for j in range(3))
            if not comb:
                continue
            r1, r2, r3, n = resid(comb, pairs)
            mxc = max(max(abs(c) for c in v) for v in comb.values()) or 1.0
            print('  모드1 %s 모드2: max|du1|=%.3e  max|du2|=%.3e  max|du3|=%.3e  잔차/크기=%.2f'
                  % (label, r1, r2, r3, max(r1, r2, r3) / mxc))
    print('  판독: 기저 잔차/크기 ~0 -> 기저 대칭(모드가 한쪽으로 보이면 기저 회전 문제)')
    print('        기저 잔차/크기 큼  -> 기저가 대칭을 깼다(개시 넘음) — 씨드 부적합 근거')
    odb.close()
    return 0


sys.exit(main())
