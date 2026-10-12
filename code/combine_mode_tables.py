#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""모드표 **합성기** — 여러 좌굴 모드표에서 모드를 뽑아 가중 합성표 하나를 만든다.

[왜 필요한가]  씨드(imperfection) = a * Sum_i w_i * phi_i.
  HF(run_abaqus.py)는 **모드표 하나**만 읽고, 그 표의 모드를 **균일 진폭**으로 더한다
  (aba_imperfection.build_perturbation: {label: dz} = Sum_modes amplitude*phi(label)).
  따라서 가족별 가중치 w_i 는 **표에 구워 넣어야** 한다 -> 이 스크립트가 각 모드의 u3 를 w 로
  스케일해 새 표(모드 1..N 재번호)를 쓴다.

[사용]
  python3 combine_mode_tables.py <out.txt> <w>:<table.txt>:<m1,m2,...> [<w>:<table.txt>:<m...> ...]
예:
  python3 combine_mode_tables.py modes_2p2.txt 1.0:modes_seed.txt:1,2  1.0:modes_clamp.txt:1,2

[규약]  각 소스 모드표의 u3 는 "모드별 max|u3| = 1, 부호는 max 노드를 +" 로 정규화돼 있다
  (aba_mode_from_odb.write_mode_table). 같은 w 를 주면 각 모드가 같은 진폭으로 더해진다.
  w 에 음수를 주면 그 가족의 패턴을 뒤집을 수 있다(상대 부호는 물리적 선택).

[검증]  출력 표는 aba_imperfection.load_mode_table 로 **다시 읽어** 값이 w*원본 인지 확인한다
  (같은 파서를 공유하므로 형식 불일치가 원천 차단된다).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from aba_imperfection import load_mode_table, mode_table_info
except Exception:                                    # pragma: no cover - 폴백(독립 실행용)
    load_mode_table = None
    mode_table_info = lambda p: {}

MODE_HDR = 'MODE'
AMPL_T_DEFAULT = 0.10     # HF 의 IMPERFECTION_AMPL_T 권장값(논문 Galhofo). 요약 리포트에 쓴다.


def read_table(path):
    """{mode: {label: u3}} + {mode: meta} 로 읽는다(load_mode_table 과 같은 규칙)."""
    if not os.path.exists(path):
        raise SystemExit("[합성] 소스 표가 없습니다: %s" % path)
    table, meta, cur = {}, {}, None
    with open(path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith('#'):
                body = s.lstrip('#').strip()
                if body.upper().startswith(MODE_HDR + ' ') and cur is None:
                    pass
                elif body.upper().startswith(MODE_HDR + ' '):
                    parts = body.split()
                    if len(parts) > 2 and parts[2].isdigit():
                        meta[int(parts[2])] = body      # '# MODE k frame=... lambda=...'
                continue
            if s.upper().startswith(MODE_HDR):
                parts = s.split()
                if len(parts) < 2 or not parts[1].isdigit():
                    raise SystemExit("[합성] %s: 잘못된 MODE 행 %r" % (path, s))
                cur = int(parts[1])
                table[cur] = {}
                continue
            if cur is None:
                continue
            parts = s.split()
            if len(parts) < 2:
                raise SystemExit("[합성] %s: 잘못된 행 %r" % (path, s))
            table[cur][int(parts[0])] = float(parts[1])
    for m in sorted(table):
        if not table[m]:
            raise SystemExit("[합성] %s: 모드 %d 가 비어 있습니다" % (path, m))
    return table, meta


def parse_spec(spec):
    """'w:table:1,2' -> (w, table, [modes])"""
    try:
        w_s, table, m_s = spec.split(':', 2)
        w = float(w_s)
        modes = [int(x) for x in m_s.replace(' ', '').split(',') if x]
    except Exception:
        raise SystemExit("[합성] 인자 형식 오류: %r  (기대 'w:table:1,2')" % (spec,))
    if not modes:
        raise SystemExit("[합성] 모드 목록이 비었습니다: %r" % (spec,))
    return w, table, modes


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    out_path, specs = argv[1], [parse_spec(a) for a in argv[2:]]

    combined, meta_out, prov, mode_src = {}, {}, [], {}
    _spec_w = [w for (w, _t, _m) in specs]
    k = 0
    for _si, (w, table_path, modes) in enumerate(specs):
        tbl, meta = read_table(table_path)
        info = mode_table_info(table_path) or {}
        print("[합성] %s  w=%g  모드=%s  (헤더: source=%s)"
              % (table_path, w, modes, info.get('source', '?')))
        for m in modes:
            if m not in tbl:
                raise SystemExit("[합성] %s 에 모드 %d 가 없습니다 (있는 것: %s)"
                                 % (table_path, m, sorted(tbl)))
            k += 1
            combined[k] = {lab: w * val for lab, val in tbl[m].items()}
            mtag = meta.get(m, '# MODE %d' % m)
            mtag = mtag.lstrip('# ').strip()
            meta_out[k] = "%s  <-- src=%s w=%g mode=%d" % (mtag, os.path.basename(table_path), w, m)
            mode_src[k] = _si                     # 이 모드가 어느 소스에서 왔는가(요약 리포트용)
        prov.append("# source%d=%g:%s:%s" % (len(prov) + 1, w, os.path.basename(table_path),
                                             ','.join(str(x) for x in modes)))
        if k > 4:
            print("[합성] 경고: 모드 %d개 -> HF 의 IMPERFECTION_MAX_MODES=4 를 넘습니다(상위 4개만 주입됨)" % k)

    with open(out_path, 'w') as f:
        f.write("# combiner=combine_mode_tables.py\n")
        f.write("# combined_modes=%d\n" % len(combined))
        f.write("\n".join(prov) + "\n")
        f.write("# 주의: 각 모드 u3 = w * (소스 표의 정규화된 u3). HF 는 균일 진폭으로 더한다.\n")
        for m in sorted(combined):
            f.write("# %s\n" % meta_out[m])
            f.write("MODE %d\n" % m)
            for lab in sorted(combined[m]):
                f.write("%d %.9e\n" % (lab, combined[m][lab]))
    print("[합성] 저장: %s (모드 %d개)" % (os.path.abspath(out_path), len(combined)))

    # ---- 합성장 요약 (HF 진폭 기준) — 40분 HF 를 쏘기 전에 씨드 크기/구성을 확인한다 ----
    #   t = 막 두께 5 µm(THICKNESS=5e-6). *IMPERFECTION 규약: u3_주입 = AMPL_T * t * u3_표.
    _t_um, _ampl_t = 5.0, AMPL_T_DEFAULT
    total = {}
    for m in sorted(combined):
        for lab, v in combined[m].items():
            total[lab] = total.get(lab, 0.0) + v
    lab_max = max(total, key=lambda l: abs(total[l]))
    per_src = {}
    for m in sorted(combined):
        if lab_max in combined[m]:
            per_src[mode_src[m]] = per_src.get(mode_src[m], 0.0) + combined[m][lab_max]
    _den = sum(abs(v) for v in per_src.values()) or 1.0
    print("[요약] 합성장: 노드 %d개 / 표 기준 max|Sum w*phi| = %.4g (노드 %d)"
          % (len(total), abs(total[lab_max]), lab_max))
    print("       HF 진폭 %.2f t (=%.3g µm) 적용 시 주입 max|Δz| = %.3g µm (%.3g t)"
          % (_ampl_t, _t_um, _ampl_t * _t_um * abs(total[lab_max]),
             _ampl_t * abs(total[lab_max])))
    print("       최대노드 기여: " + ", ".join(
        "source%d %.0f%%" % (k + 1, 100.0 * abs(v) / _den) for k, v in sorted(per_src.items())))
    for m in sorted(combined):
        _lab = max(combined[m], key=lambda l: abs(combined[m][l]))
        print("       MODE %d (source%d, w=%g): max|u3|=%.4g @node %d"
              % (m, mode_src[m] + 1, _spec_w[mode_src[m]], abs(combined[m][_lab]), _lab))

    if load_mode_table is not None:                  # 실제 HF 로더로 재검증
        want = tuple(sorted(combined))
        got = load_mode_table(out_path, want)
        bad = 0
        for m in want:
            for lab, v in combined[m].items():
                if abs(got[m][lab] - v) > 1e-15:
                    bad += 1
        peak = max(abs(v) for m in combined for v in combined[m].values())
        print("[검증] load_mode_table 재독 OK — 모드 %s / 노드 %d개/모드 / 원본 대비 불일치 %d건 / peak|u3|=%.4g"
              % (list(want), len(combined[want[0]]), bad, peak))
        return 1 if bad else 0
    print("[검증] aba_imperfection 임포트 실패 -> 형식 검증 생략")
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
