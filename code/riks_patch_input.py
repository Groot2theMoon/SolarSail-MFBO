# -*- coding: utf-8 -*-
r"""
Step-Postbuckle 을 Riks 스텝으로 바꾼 입력파일을 만든다 (원본은 건드리지 않는다).

  python riks_patch_input.py [원본.inp] [출력.inp] [스텝이름]

기본값: aba\HF_Postbuckle.inp -> aba\HF_Riks.inp, 스텝 Step-Postbuckle

바꾸는 것 (딱 3줄):
  1) *Step 줄의 inc=1000 -> inc=10000            (최대 증분 수; 리크스에서는 아크길이 증분 수)
  2) *Static, stabilize=..., allsdtol=..., continue=NO  ->  *Static, riks
     (STABILIZE/FACTOR 는 RIKS 와 함께 쓸 수 없다 - 키워드 문서)
  3) 그 데이터 줄 -> 리크스 8항목 중 앞 5개:
     dl_in, l_period, dl_min, dl_max, lpf_end
       0.05,   1.0,    1e-05, 0.1,    1.0
     lpf_end=1.0 이므로 BC 가 '스텝에서 지정한 최종값'(= 설계 운용점)에 닿으면 스텝이 끝난다.
     노드/dof/값(6~8항목)은 쓰지 않는다.

검증: v7 프로브에서 확인된 사실만 쓴다.
  - lambda 는 prescribed displacement 의 '값'도 스케일한다 (스텝 중 BC 값 변경이 실제로 진행됨)
  - 리크스 데이터 라인 순서는 dl_in, l_period, dl_min, dl_max, lambda_end, node, dof, value
    (lambda_end 를 3번째에 두면 min/max 로 읽혀 입력단계에서 거부된다)
"""
import io
import os
import sys

_argv = [a for a in sys.argv[1:] if not a.startswith('--')]
LPF = None
LS = ('--line-search' in sys.argv[1:])
RC = ('--relax-corr' in sys.argv[1:])
FORCE_RIKS = ('--riks' in sys.argv[1:])
ARC = None                      # --arc=<dl_in>: 아크길이 직접 지정(문서 기본인 '상한 없음' 동반)
ARC_MAX = None                  # --arc-max=<v>: 굳이 상한을 걸고 싶을 때만
TSTOP = None                     # --tstop=<step time>: *Static 의 timePeriod(2번째 항목) 를 줄여 조기 종료
for _a in sys.argv[1:]:
    if _a.startswith('--tstop='):
        try:
            TSTOP = float(_a.split('=', 1)[1])
        except ValueError:
            print('[중단] --tstop= 뒤에는 숫자(step time)가 와야 합니다: %r' % _a)
            sys.exit(2)
if TSTOP is not None and FORCE_RIKS:
    print("[중단] --tstop 과 --riks 는 함께 쓸 수 없습니다 (Riks 는 timePeriod 대신 아크길이를 받습니다).")
    sys.exit(2)
for a in sys.argv[1:]:
    if a.startswith('--lpf='):
        LPF = float(a.split('=', 1)[1])
    elif a.startswith('--arc='):
        ARC = float(a.split('=', 1)[1])
    elif a.startswith('--arc-max='):
        ARC_MAX = float(a.split('=', 1)[1])
SRC = _argv[0] if len(_argv) > 0 else os.path.join('aba', 'HF_Postbuckle.inp')
DST = _argv[1] if len(_argv) > 1 else os.path.join('aba', 'HF_Riks.inp')
STEP = _argv[2] if len(_argv) > 2 else 'Step-Postbuckle'

RIKS_STATIC = '*Static, riks'
RIKS_DATA = '0.05, 1.0, 1e-05, 0.1, %s'   # 마지막 항목 = lpf_end (기본 1.0)


def riks_data(lpf, arc=None, arc_max=None):
    """리크스 데이터 라인: dl_in, l_period, dl_min, dl_max, lambda_end.

    아바쿠스 *STATIC 키워드 문서(2017, Data line for the Riks method) 기준:
      1항 dl_in   : 0/미지정이면 '스텝의 총 아크길이' 가 기본 = 첫 증분에서 스텝 전체.
      2항 lperiod : 기본 1.0. 첫 증분의 lambda 증가는 dlambda_in = dl_in / lperiod.
      3항 dl_min  : 0이면 min(제안 dl_in, 1e-5 x 총 아크길이) 가 기본.
      4항 dl_max  : **미지정이면 상한 없음**(no upper limit is imposed).
      5항 lpf_end : lambda 상한.

    => 그래서 dl_max 를 지정하는 것은 자동 증분이 아크길이를 '늘릴' 여지를 잘라먹는다.
       이전 실패 런들은 dl_max 를 0.1 / 0.025 로 조여놓았고, 더 좁게 조인 쪽(V1)이 더 못 갔다.
       문서 기본(상한 없음)으로 되돌리려면 4항을 비운다.
    """
    if lpf is None:
        lpf = 1.0
    if arc is not None:
        dl_in = arc
        dl_min = 1e-06
        dl_max_s = '' if arc_max is None else ('%.6g' % arc_max)
    else:
        dl_in = min(0.05, max(0.001, lpf / 5.0))
        dl_min = 1e-07
        dl_max_s = '%.6g' % (arc_max if arc_max is not None
                             else min(0.1, max(0.002, lpf / 2.0)))
    return '%.6g, 1.0, %.6g, %s, %.6g' % (dl_in, dl_min, dl_max_s, lpf)
NEW_INC = 10000

if not os.path.exists(SRC):
    sys.exit('ERROR: 원본을 찾지 못했습니다: %s' % SRC)

lines = io.open(SRC, 'r', encoding='utf-8', errors='replace').read().splitlines()

# 1) 대상 *Step 줄
si = None
for i, l in enumerate(lines):
    ls = l.strip().lower()
    if ls.startswith('*step') and ('name=' + STEP).lower() in ls:
        si = i
        break
if si is None:
    sys.exit('ERROR: *Step, name=%s 를 찾지 못했습니다' % STEP)

# 2) 그 스텝 안의 *Static 줄 (+다음 줄이 데이터 줄)
sidx = None
for j in range(si, min(si + 60, len(lines))):
    if lines[j].strip().lower().startswith('*static'):
        sidx = j
        break
if sidx is None or sidx + 1 >= len(lines):
    sys.exit('ERROR: %s 스텝 안에서 *Static + 데이터 줄을 찾지 못했습니다' % STEP)

print('[찾음] line %d: %s' % (si + 1, lines[si]))
print('[찾음] line %d: %s' % (sidx + 1, lines[sidx]))
print('[찾음] line %d: %s' % (sidx + 2, lines[sidx + 1]))

if lines[sidx].strip().lower().replace(' ', '').startswith('*static,riks'):
    print('[안내] 이미 Riks 스텝입니다 - 그대로 복사합니다.')

# 데이터 줄 형식 확인 (예상: initialInc, timePeriod, minInc, maxInc = 4항목)
old_data = [x.strip() for x in lines[sidx + 1].split(',')]
print('[확인] 기존 데이터 항목 수 = %d (%s)' % (len(old_data), lines[sidx + 1].strip()))
if len(old_data) > 4 and old_data[0] != '':
    print('[경고] 항목이 4개를 넘습니다 - 이 스텝이 정말 일반 Static 인지 확인하세요.')

# [2026-09-29] --tstop: *Static 데이터 줄 = initialInc, timePeriod, minInc, maxInc
#   2번째 항목(timePeriod)을 줄여 스텝을 조기에 끝낸다. 기본 진폭(ramp)이 스텝 타임을 따라가므로
#   두 런을 '같은 하중 수준'에서 대조할 수 있다 (C_n^a 완화 타당성 검증용).
if TSTOP is not None:
    if len(old_data) < 2:
        print('[중단] 데이터 줄에 timePeriod 항목이 없습니다: %r' % lines[sidx + 1])
        sys.exit(2)
    # [2026-09-30] 함정 실측: *Static 데이터 줄의 initialInc/minInc/maxInc 는 **스텝타임 단위**다.
    #   timePeriod 만 0.05 로 줄이면 같은 값이 램프 기준으로는 20배 큰 증분이 되어
    #   (초기증분 0.0001 -> 램프의 2e-3 = 2 um, 완주 런은 1e-4 = 0.1 um) 취약한 상태에서
    #   C_n^a 보정검정이 깨지고 증분이 나노미터급으로 붕괴한다(HF_ref 실측 5.06e-8).
    #   -> 세 증분 파라미터를 같은 비율로 스케일해 **물리적 증분 이력을 그대로 보존**한다.
    _old_tp = float(old_data[1]) if len(old_data) > 1 else 1.0
    _ratio = (TSTOP / _old_tp) if _old_tp else 1.0
    for _k in (0, 2, 3):
        if len(old_data) > _k:
            try:
                old_data[_k] = '%g' % (float(old_data[_k]) * _ratio)
            except ValueError:
                pass
    old_data[1] = '%g' % TSTOP
    lines[sidx + 1] = ', '.join(old_data)
    print('[진단] timePeriod %g -> %g, 증분 파라미터 x%g (물리적 증분 이력 보존)'
          % (_old_tp, TSTOP, _ratio))

# 3) 세 줄 교체
new_step = lines[si]
if 'inc=' in new_step.lower():
    import re
    new_step = re.sub(r'inc=\s*\d+', 'inc=%d' % NEW_INC, new_step, flags=re.I)
    print('[교체] *Step: %s' % new_step)

out = list(lines)
out[si] = new_step
if LS:
    # [2026-09-29] 로그 관측: attempt 1 에서 한 노드의 보정이 부호를 바꾸며 커지고
    #   'DISP. CORRECTION TOO LARGE' -> 'APPEARS TO BE DIVERGING'. 이 서명은 line search 가
    #   겨냥하는 상황이다(보정 방향을 감쇠). 정적 해석 + 기존 안정화를 그대로 두고 삽입만 한다.
    #   *Controls 는 *Step 뒤, 절차 키워드(*Static) 앞에 온다.
    # 배치: *Controls 는 'Type: History data / Level: Step' 이므로 절차 키워드(*Static)와
    #   그 데이터 줄 '뒤'에 온다. *Step 과 *Static 사이에 넣으면 입력단계에서 거부된다(실측).
    # 데이터 줄: Nls = line search 최대 반복. 문서상 기본값이 Newton 스텝에서 **0(비활성)**이므로
    #   반드시 값을 준다(권장 Nls=5).
    out.insert(sidx + 2, '5,')
    out.insert(sidx + 2, '*Controls, parameters=line search')
if RC:
    # [2026-09-29] 실패 서명: 잔차는 통과(5e-7 / 평균 1.87e-3 = 2.7e-4 < Rαn 5e-3)인데
    #   'DISP. CORRECTION TOO LARGE COMPARED TO DISP. INCREMENT' 로 계속 실패한다.
    #   원인은 보정/증분 비 기준 Cαn(기본 1e-2) 인데, 증분이 1.09e-6 m(= 0.011 t)까지
    #   작아지면 보정 2.5e-7 이 23% 가 되어 그 비율은 만족될 수 없다(컷백 함정).
    #   문서: 'in cases where the incremental solution is essentially zero' 에서 Cαn 등을
    #   수정해야 할 수 있고, 'To avoid testing the magnitude of the solution correction,
    #   you can set Cαn to 1.'  (Analysis UG, Commonly used control parameters)
    #   잔차 기준 Rαn 은 건드리지 않는다 -> 평형 정확성 근거는 유지된다.
    out.insert(sidx + 2, ', 1.0, ,')
    out.insert(sidx + 2, '*Controls, parameters=field, field=displacement')
# 리크스 전환 판정: 옵션을 주지 않으면 리크스(기존 기본), --riks 를 주면 옵션과 무관하게 리크스.
#   이렇게 해야 '리크스 + C_n^a 완화' 라는 공정한 조합을 만들 수 있다.
DO_RIKS = ('--riks' in sys.argv[1:]) or (not (LS or RC or (TSTOP is not None)))
if DO_RIKS:
    out[sidx] = RIKS_STATIC
    out[sidx + 1] = riks_data(LPF, ARC, ARC_MAX)
if LPF is not None:
    print('[진단] lpf_end = %g -> 램프의 %g%% 지점에서 스텝을 끝낸다' % (LPF, LPF * 100))

_raw = io.open(SRC, 'rb').read()
_nl = '\r\n' if b'\r\n' in _raw else '\n'
# 가드: 리크스 스텝은 데이터 덱의 마지막 스텝이어야 한다.
#   실측 오류(2026-09-29): ***ERROR: IF A RIKS STEP IS SPECIFIED IT MUST BE THE LAST STEP
#   IN A DATA DECK. ADDITIONAL STEPS MAY BE DEFINED VIA THE *RESTART OPTION.
#   원본에 Postbuckle 뒤 스텝이 남아 있으면 그 입력은 입력단계에서 즉시 죽는다 -> 쓰지 않고 중단한다.
_tail = [] if not DO_RIKS else [(sidx + 3 + k, l) for k, l in enumerate(out[sidx + 2:])
         if l.strip().lower().startswith('*step')]
if _tail:
    print('ERROR: %s 뒤에 스텝이 더 있습니다 - 리크스 스텝은 마지막이어야 합니다.' % STEP)
    for ln, txt in _tail[:5]:
        print('   line %d: %s' % (ln, txt.strip()))
    print('   -> 원본이 최신 HF 입력(aba\\HF_Postbuckle.inp)이 맞는지 확인하세요.')
    print('   -> 스텝을 지우거나, 리크스 스텝을 덱의 끝으로 옮긴 뒤 다시 실행하세요. 쓰지 않았습니다.')
    sys.exit(2)

io.open(DST, 'w', encoding='ascii', errors='replace', newline=_nl).write(_nl.join(out) + _nl)
print('[줄바꿈] %s' % ('CRLF' if _nl == '\r\n' else 'LF'))
print('[저장] %s (%d 줄)' % (DST, len(out)))
print('--- 새 스텝 블록 (앞 7줄, *Controls 포함) ---')
for l in out[si:si + 7]:
    print('   ', l)
print('[다음] abaqus job=%s interactive' % os.path.splitext(os.path.basename(DST))[0])
