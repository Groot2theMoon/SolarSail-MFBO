# -*- coding: utf-8 -*-
"""post.py — 설계점 결과 아카이빙 (results_*.txt / all_results.csv / slope / DN / ZK).

run_abaqus.py 에서 **함수만** 옮겼다(2026-11-04, P2). 수치 경로/계약은 바뀌지 않는다:
  - `RESULTS:` 4필드(lf1,lf2,lf3,hf)는 이 모듈이 읽지도 쓰지도 않는다.
  - `all_results.csv` 는 `archive_result` 만 쓰고, 열 스키마는 고정이다(§19.91).
  - 모든 함수는 예외 안전이다 — 실패해도 해석은 계속되고 경고 1줄만 남는다.

경로는 호출자가 넘긴다(모듈이 스크립트의 `_HERE`/`_RUN`/`HF_ODB` 전역을 읽지 않는다).
"""

import io
import os
import shutil
import sys

from aba_imperfection import lf_sidecar_read, value_tag
from ctx import emit


def archive_result(here, cfg_sha1, lf_odb, x_c, d_c, fidelity, results_csv, extra=None):
    """[§19.88] 설계점별 결과 아카이빙.

    왜: HF 잡 이름이 `HF_Postbuckle` 고정이라 다음 런이 **직전 산출물을 지운다**
    (`run_job_safely` 가 odb/dat/msg 를 제거). 그래서 (0.5,0.3)의 ODB 를 이미 잃었다 —
    결과를 설계점 이름의 텍스트 + 누적 CSV 로 남겨 스윕 분석이 한 파일로 가능하게 한다.

    남기는 것: RESULTS 4필드(lf1,lf2,lf3,hf) + 사다리가 쓴 stabilize_used/에너지 비/판정 + 시각.
    """
    try:
        import time
        _dir = os.path.join(here, 'results')
        if not os.path.isdir(_dir):
            os.makedirs(_dir)
        _p = [str(v).strip() for v in (results_csv or '').split(',')]
        while len(_p) < 4:
            _p.append('nan')
        _kv = {'x_c': x_c, 'd_c': d_c, 'fidelity': fidelity,
               'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
               'lf1': _p[0], 'lf2': _p[1], 'lf3': _p[2], 'hf': _p[3]}
        _kv.update(extra or {})
        #   [CONFIG] 이 점이 **어느 config/모델**로 계산됐는지 남긴다. HF 는 짝 LF ODB 를
        #     재사용할 수 있으므로 "이 런의 config"(cfg_sha1)와 "LF ODB sidecar 의 config"를
        #     따로 적는다(둘이 다르면 재사용이 섞였다는 신호다). CSV 열은 늘리지 않는다.
        _kv.setdefault('config_sha1', cfg_sha1)
        try:
            _lfsd = lf_sidecar_read(lf_odb or '')
        except Exception:
            _lfsd = {}
        _kv.setdefault('lf_config_sha1', _lfsd.get('config_sha1', 'nan'))
        _kv.setdefault('lf_model_fp', _lfsd.get('model_fp', 'nan'))
        _tag = 'results_xc%s_dc%s.txt' % (value_tag(x_c), value_tag(d_c))
        with io.open(os.path.join(_dir, _tag), 'w', encoding='utf-8') as _f:
            for _k in sorted(_kv):
                _f.write('%s=%s\n' % (_k, _kv[_k]))
        #   [§19.91] 열 스키마를 **고정**한다 — 진단키(stabilize/에너지/flip_bias)도 CSV 에 실제로
        #     실려야 스윕 분석이 한 파일로 된다(헤더가 1회만 쓰이므로 열 목록은 상수여야 한다).
        _cols = ['x_c', 'd_c', 'fidelity', 'timestamp', 'stabilize_used',
                 'allds_allie_last_pct', 'allds_allie_max_pct', 'energy_verdict',
                 'backward_area_pct', 'backward_elems', 'total_elems',
                 'hf_flip', 'hf_true_normal', 'delta_pct',
                 'lf1', 'lf2', 'lf3', 'hf']
        _csv = os.path.join(_dir, 'all_results.csv')
        _new = not os.path.exists(_csv)
        with io.open(_csv, 'a', encoding='utf-8') as _f:
            if _new:
                _f.write(','.join(_cols) + '\n')
            _f.write(','.join(str(_kv.get(_c, '')) for _c in _cols) + '\n')
        emit("[ARCHIVE] %s + all_results.csv (1행 추가)" % _tag)
    except Exception as _e:
        emit("[ARCHIVE] 실패(무시): %s" % _e)


def archive_slope(here, run, hf_odb, x_c, d_c):
    """[§19.102] **성공한 HF 런에서만** 호출: (1) HF ODB 설계점 태깅 복사,
    (2) slope_stats.compute_stats 를 인프로세스로 돌려 점별 results txt 에 slope_ 접두로 아카이빙.

    모두 예외 안전 - 실패해도 런을 죽이지 않고 경고 1줄만 남긴다. RESULTS 4필드 계약과
    all_results.csv 열 스키마는 건드리지 않는다(점별 txt append + 별도 JSON 만).
    스위치: SOLARSAIL_NO_SLOPE=1 -> 전체 생략, SOLARSAIL_NO_ODB_COPY=1 -> ODB 복사만 생략.
    로그 접두 [SLOPE]/[ODBCOPY] 는 기존 stdout 필터가 버리지 않는 형식이다.
    """
    _tag = 'xc%s_dc%s_u3free' % (value_tag(x_c), value_tag(d_c))
    _src = os.path.join(run, hf_odb)
    _src = _src if os.path.exists(_src) else (hf_odb if os.path.exists(hf_odb) else None)

    if str(os.environ.get('SOLARSAIL_NO_SLOPE', '')).strip() == '1':
        emit("[SLOPE] SKIP (SOLARSAIL_NO_SLOPE=1): 설계점 %s slope_stats 아카이빙 생략" % _tag)
        return

    # ---- (1) ODB 설계점 태깅 복사: results/odb/HF_<tag>.odb (태그 규칙 = value_tag, 기존 관례)
    _odb_used = _src
    if str(os.environ.get('SOLARSAIL_NO_ODB_COPY', '')).strip() == '1':
        emit("[ODBCOPY] SKIP (SOLARSAIL_NO_ODB_COPY=1): 설계점 %s ODB 복사 생략" % _tag)
    else:
        try:
            if _src is None:
                emit("[ODBCOPY] !!! HF ODB 를 찾지 못했다(%s) - 복사 생략" % hf_odb)
            else:
                _odir = os.path.join(here, 'results', 'odb')
                if not os.path.isdir(_odir):
                    os.makedirs(_odir)
                _dst = os.path.join(_odir, 'HF_%s.odb' % _tag)
                shutil.copy2(_src, _dst)
                _odb_used = _dst
                _sz = os.path.getsize(_dst)
                emit("[ODBCOPY] %s -> %s (%.2f MB, %d B)"
                     % (os.path.abspath(_src), os.path.abspath(_dst), _sz / 1048576.0, _sz))
        except Exception as _e:
            emit("[ODBCOPY] 실패(무시): %s" % _e)

    # ---- (2) slope_stats 인프로세스 실행 (Abaqus 안이므로 서브프로세스 금지)
    try:
        if _odb_used is None:
            emit("[SLOPE] !!! 읽을 ODB 가 없다(%s) - 통계 생략" % hf_odb)
            return
        if here not in sys.path:
            sys.path.insert(0, here)
        import slope_stats
        _S = slope_stats.compute_stats(_odb_used)
        if _S is None:
            emit("[SLOPE] !!! compute_stats 실패(None) - 아카이빙 생략")
            return
        # 점별 results txt 에 slope_ 접두 key=value append (RESULTS 4필드 계약은 그대로)
        _lines = slope_stats.scalar_items(_S, prefix='slope_')
        _ptxt = os.path.join(here, 'results',
                             'results_xc%s_dc%s.txt' % (value_tag(x_c), value_tag(d_c)))
        with io.open(_ptxt, 'a', encoding='utf-8') as _f:
            for _l in _lines:
                _f.write(_l + '\n')
        emit("[SLOPE] %d개 key(slope_ 접두)를 %s 에 append | f1=%.10e f2=%.10e f3=%.10e f4=%.10e"
             % (len(_lines), os.path.basename(_ptxt), _S.get('f1_hf_pred'),
                _S.get('f2_p95'), _S.get('f3_worst_patch'), _S.get('f4_cvar90')))
        # JSON 은 slope_stats 기존 관례: results/slope_<stem>_<step>.json
        try:
            _stem = os.path.splitext(os.path.basename(_odb_used))[0]
            _jpath = os.path.join(here, 'results', 'slope_%s_%s.json' % (_stem, _S['step']))
            _rec = slope_stats.build_json_record(_S, _odb_used)
            slope_stats.write_json_record(_rec, _jpath)
            emit("[SLOPE] JSON 저장: %s" % os.path.abspath(_jpath))
        except Exception as _e:
            emit("[SLOPE] !!! JSON 저장 실패(무시): %s" % _e)
    except Exception as _e:
        emit("[SLOPE] 실패(무시, 런은 계속): %s" % _e)


def archive_dn(here, x_c, d_c, dn_lines):
    """[19.106] eval_abaqus.py 가 낸 `DN:` 줄을 점별 results txt 에 `dn_` 접두로 append.

    계약 보존:
      - `RESULTS:` 4필드(lf1,lf2,lf3,hf)는 **읽지도 쓰지도 않는다**.
      - `all_results.csv` 는 이 함수에서 **열지 않는다**(헤더 1회 규칙 - 열 추가 금지).
      - 실패해도 런을 죽이지 않는다(경고 1줄). DN 줄이 아예 없으면 조용히 넘기지 않고 경고한다.
    호출 위치: archive_result('w') **뒤** (그 함수가 txt 를 다시 쓰므로).
    """
    try:
        _ptxt = os.path.join(here, 'results',
                             'results_xc%s_dc%s.txt' % (value_tag(x_c), value_tag(d_c)))
        if not os.path.exists(_ptxt):
            emit('[DN] !!! 점별 txt 없음(%s) - dn_ 아카이빙 생략' % os.path.basename(_ptxt))
            return
        _kv = {}
        _fail = None
        for _ln in (dn_lines or []):
            _s = str(_ln).strip()
            if not _s.startswith('DN:'):
                continue
            _body = _s[3:].strip()
            if _body.startswith('FAIL'):
                _fail = _body
                continue
            for _tok in _body.split():
                if '=' in _tok:
                    _k, _v = _tok.split('=', 1)
                    _kv[_k if _k.startswith('dn_') else ('dn_' + _k)] = _v

        if _fail is not None and not _kv:
            #   계산 불가: 0 으로 뭉개지 않는다(nan + 사유). CSV 는 그대로.
            with io.open(_ptxt, 'a', encoding='utf-8') as _f:
                _f.write('dn_deg=nan\n')
                _f.write('dn_note=%s\n' % _fail)
            emit('[DN] !!! DN 계산 실패 -> dn_deg=nan (%s)' % _fail)
            return

        if not _kv:
            emit('[DN] !!! 줄 없음 (eval_abaqus.py 출력에 DN: 이 없다) - dn_ 아카이빙 생략')
            return

        _order = ['dn_deg', 'dn_psi_deg', 'dn_n_eq_x', 'dn_n_eq_y', 'dn_n_eq_z',
                  'dn_n_elem', 'dn_area_ref']
        _keys = [k for k in _order if k in _kv] + sorted(k for k in _kv if k not in _order)
        with io.open(_ptxt, 'a', encoding='utf-8') as _f:
            for _k in _keys:
                _f.write('%s=%s\n' % (_k, _kv[_k]))
        emit('[DN] %d개 key(dn_ 접두)를 %s 에 append | dn_deg=%s deg dn_psi_deg=%s '
             'n_eq=(%s,%s,%s) n_elem=%s'
             % (len(_keys), os.path.basename(_ptxt), _kv.get('dn_deg', 'nan'),
                _kv.get('dn_psi_deg', 'nan'), _kv.get('dn_n_eq_x', 'nan'),
                _kv.get('dn_n_eq_y', 'nan'), _kv.get('dn_n_eq_z', 'nan'),
                _kv.get('dn_n_elem', 'nan')))
    except Exception as _e:
        emit('[DN] 실패(무시, 런은 계속): %s' % _e)


def archive_zk(here, x_c, d_c, zk_lines, zkiw_lines=None):
    """[paper-match] eval_abaqus.py 가 낸 `ZK:` 줄을 점별 results txt 에 `zk_` 접두로 append.

    `zk_metrics.txt` sidecar 의 추가 키(`zk_Iw`/기준면적/가중모드)도 함께 append 한다.
    [A2] `ZKIW:` 줄이 stdout 의 **정본**이다 -> sidecar 병합이 실패해도 값이 남는다.

    계약 보존:
      - `RESULTS:` 3/4필드는 읽지도 쓰지도 않는다.
      - `all_results.csv` 는 열지 않는다(헤더 1회 규칙 - 열 추가 금지).
    호출 위치: archive_result('w') 뒤 (그 함수가 txt 를 다시 쓰므로).
    """
    try:
        _ptxt = os.path.join(here, 'results',
                             'results_xc%s_dc%s.txt' % (value_tag(x_c), value_tag(d_c)))
        if not os.path.exists(_ptxt):
            emit('[ZK] !!! 점별 txt 없음(%s) - zk_ 아카이빙 생략' % os.path.basename(_ptxt))
            return
        _kv = {}
        _fail = None
        for _ln in (zk_lines or []):
            _s = str(_ln).strip()
            if not _s.startswith('ZK:'):
                continue
            _body = _s[3:].strip()
            if _body.startswith('FAIL'):
                _fail = _body
                continue
            for _tok in _body.split():
                if '=' in _tok:
                    _k, _v = _tok.split('=', 1)
                    if _k == 'state_wrinkled_frac':
                        _k = 'zk_state_wrinkled_frac'
                    elif not _k.startswith('zk_'):
                        _k = 'zk_' + _k
                    _kv[_k] = _v

        #   [A2] ZKIW: 는 stdout 이 정본이다(sidecar 병합이 실패해도 값이 남는다).
        #     ZK: 4키와 겹치지 않으므로 그대로 덮어쓴다.
        for _ln in (zkiw_lines or []):
            _s = str(_ln).strip()
            if not _s.startswith('ZKIW:'):
                continue
            _body2 = _s[5:].strip()
            if _body2.startswith('FAIL'):
                continue
            for _tok in _body2.split():
                if '=' not in _tok:
                    continue
                _k, _v = _tok.split('=', 1)
                if not _k.startswith('zk_'):
                    _k = 'zk_' + _k      #   sidecar 병합과 같은 키 이름(zk_strain_field 등)
                _kv[_k] = _v

        #   `zk_metrics.txt` sidecar 의 추가 키(`zk_Iw` 등)도 점별 txt 에 남긴다.
        #     ZK: 줄은 public 4키만 갖고, I_w/기준면적/가중모드는 sidecar 에만 있다.
        #   [A1] ZK:FAIL 이어도 sidecar 를 병합한다 — 원인 단서(fields/frame_value)가
        #     거기에만 있고, 다음 런이 zk_metrics.txt 를 덮으면 영구 소실된다.
        #     stale 차단은 "FAIL 에서는 진단 키만 통과"로 대신한다(값 키는 병합 금지).
        _FAILDX = ('fields', 'frame_value', 'reason', 'odb')
        if _fail is not None or _kv:
            try:
                _zp = os.path.join(os.getcwd(), 'zk_metrics.txt')
                if os.path.exists(_zp):
                    with io.open(_zp, 'r', encoding='utf-8', errors='replace') as _zf:
                        for _ln in _zf:
                            _s = _ln.strip()
                            if not _s or '=' not in _s:
                                continue
                            _k, _v = _s.split('=', 1)
                            if _k == 'status':
                                continue
                            if _fail is not None and _k not in _FAILDX:
                                continue
                            if _k == 'state_wrinkled_frac':
                                _k = 'zk_state_wrinkled_frac'
                            elif _k in ('step', 'frame', 'odb', 'weights', 'note'):
                                _k = 'zk_' + _k
                            elif not _k.startswith('zk_'):
                                _k = 'zk_' + _k
                            if _k not in _kv:      #   stdout 우선
                                _kv[_k] = _v
            except Exception as _e:
                emit('[ZK] zk_metrics.txt sidecar merge 실패(무시): %s' % _e)

        with io.open(_ptxt, 'a', encoding='utf-8') as _f:
            if _fail is not None:
                #   계산 불가: 0 으로 뭉개지 않는다(nan + 사유). CSV 는 그대로.
                _f.write('zk_state_wrinkled_frac=nan\n')
                _f.write('zk_note=%s\n' % _fail.replace(' ', '_'))
                for _k in sorted(_kv):
                    _f.write('%s=%s\n' % (_k, _kv[_k]))
                emit('[ZK] !!! ZK 계산 실패 -> zk_state_wrinkled_frac=nan (%s) / 진단키 %d개'
                     % (_fail, len(_kv)))
                return
            if not _kv:
                emit('[ZK] !!! 줄 없음 (eval_abaqus.py 출력에 ZK: 이 없다) - zk_ 아카이빙 생략')
                return
            for _k in sorted(_kv):
                _f.write('%s=%s\n' % (_k, _kv[_k]))
        emit('[ZK] %d개 key를 %s 에 append | state_wrinkled_frac=%s dir=%s'
             % (len(_kv), os.path.basename(_ptxt),
                _kv.get('zk_state_wrinkled_frac', 'nan'), _kv.get('zk_dir_deg', 'nan')))
    except Exception as _e:
        emit('[ZK] 실패(무시, 런은 계속): %s' % _e)


def r13_ct_keys(r13_ct):
    """[§19.97] 점별 results_*.txt 에 넣을 R-13 Step-ClampTension 판독 키.

    **CSV(all_results.csv) 열은 건드리지 않는다** - 헤더가 1회만 쓰이므로 열을 추가하면 기존 행과
    어긋난다(§19.91). 그래서 진단값은 점별 txt(key=value)로만 남긴다.
    읽지 못한 항목은 'nan' + r13_ct_note 에 이유를 남긴다(조용한 0 금지).
    """
    _d = {'r13_ct_centre_mean_pa': 'nan', 'r13_ct_centre_minp_pa': 'nan',
          'r13_ct_area_frac_neg': 'nan', 'r13_ct_minp_frac_neg': 'nan',
          'r13_ct_uz_max_m': 'nan', 'r13_ct_uz_ave_m': 'nan',
          'r13_ct_odb': 'nan', 'r13_ct_odb_mtime': 'nan', 'r13_ct_stale': 'nan',
          'r13_ct_note': 'probe 미실행(BASE_PROBE=False, ODB 없음, 또는 비-HF 경로)'}
    _d.update(r13_ct or {})
    return _d
