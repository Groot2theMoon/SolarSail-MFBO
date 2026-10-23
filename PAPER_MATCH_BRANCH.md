# Galhofo 2022 paper-match branch

Branch: `feature/galhofo-paper-match`
Base: `af89de7` (refactor/buckle-case-clamp-lf)

This branch adds a **paper-matching preset** for Galhofo et al., *Thin-Walled Structures* 181 (2022) 110007.
It is a research branch: the MFBO baseline is preserved behind `PAPER_MATCH = False`.

## What changes when `PAPER_MATCH = True`

| Item | Paper | This branch |
|---|---|---|
| Membrane elements | STRI65 + S8R5 (2nd order), 10,100 | **S4/S3 (1st order), 10,100 / 10,201 nodes** (default); 2nd-order generator exists but is experimental |
| Pre-tension load | cable initial stress 2.69e7 Pa (5PCS), cable ends fixed | **3 cables, sigma11 = 2.69e7 Pa**, ends fixed, no displacement ramp |
| Prescribed displacement | static 0.001 / 0.0005 m (2:1), buckle 0.10 m (all five) | buckle **0.01 m** (paper 0.10 m reference; λ×P invariant); static **shape** 2:1, magnitude **0.0001 m** (symmetric 10-mode base); **LF and HF post-buckling: no vertex displacement** (initial stress only) |
| BC | three edges u3 = 0 | **LF: full-membrane u3=0 from Initial; HF paper: three edges only (u3-free); HF old disp: hold -> release** |
| Stabilization | 2e-4 | **LF GlobalTension 1.0; HF GlobalTension 2e-4 ladder (2e-4/1e-3/1e-2); HF Postbuckle 5e-4 ladder (5e-4/2e-3/1e-2)** |
| Clamps | none (5PCS) | paper reproduction: `NO_CLAMP = True`; physical model: `PAPER_INCLUDE_CLAMPS = True` (clamps + vertex cables; clamp prestress = `d_c × 2.69e7 Pa`) |
| Imperfection | first 4 modes, 0.10 t | same; mode table = **seed + clamp_only composite** (2+2) |

## Constants (edit only these)

`code/run_abaqus.py`:
```python
PAPER_MATCH = True
PAPER_MESH_ORDER = 1        # 1차 기본 (2차는 실험용)
CABLE_PRESTRESS = 2.69e7    # Pa; calibrate to 7000 Pa at CENTRE
PAPER_STABILIZE = 2.0e-4
USE_COMPOSITE_MODE_TABLE = True   # seed + clamp_only 합성표
# [2026-10-21] u3-free HF 분기
STABILIZE_GLOBAL_LF = 1.0             # LF: 기존 값 유지 (u3=0)
STABILIZE_GLOBAL_LADDER = (2e-4, 1e-3, 1e-2)   # HF: 내리는 방향
STABILIZE_POST_LADDER = (5e-4, 2e-3, 1e-2)
# HF Step-GlobalTension initialInc=0.001 / maxNumInc=10000
# LF Step-GlobalTension initialInc=0.01 / maxNumInc=1000 (회귀 기준선 보존)
```

`code/run_abaqus_buckle.py`:
```python
PAPER_MATCH = True
PAPER_MESH_ORDER = 1
MODE_PRESTRESS = 0.0001     # symmetric 10-mode base (2026-10-20)
PAPER_STATIC_DISP = 0.0005  # paper reference only (unstable here)
# derived: CASE='paper', paper3, DISP_TOP_OVER_CORNER=2.0,
#          PERTURBATION=0.01, MODE_STABILIZATION=5e-4
# 2026-10-21 refactor: PULL_TOP/TOP_PULL_*, CLAMP_DIR, BUCKLE_SOLVER/LANCZOS,
#   CLEAN_PERTURBATION, 자유메쉬 경로 제거. PULL_TOP 전용 하네스는 retired
#   (3-tendon 소스는 CASE='paper' 하나).
```

`code/run_abaqus_mode.py`:
```python
PAPER_MATCH = True
PAPER_MESH_ORDER = 1
# mode source keeps DISP_GLOBAL = 0.0001 (symmetric 10-mode base) on purpose;
# modes are only the imperfection shape basis.
```

### Why the static step uses 0.0001 m, not 0.0005 m

The paper's 0.0005 m corner displacement puts the **2nd-order** base state above its bifurcation load.
Measured on the 2nd-order path (S8R5/STRI65, 202,050 equations): 592 negative eigenvalues in the system matrix,
166,009 negative diagonal entries, and the subspace iteration ends with
`THE EIGENVALUES CANNOT BE FOUND ... INSTABILITIES IN THE BASE STATE` (0 modes).
The same failure was already recorded at 1e-3 m. Modes are only the imperfection shape basis;
HF sets the amplitude at 0.10 t, so the mode-extraction prestress does not need to be the paper's.
The 2:1 shape is kept; the branch uses `PERTURBATION=0.01` (the paper's 0.10 m is the reference;
λ scales as 1/P and the eigenvectors are unchanged).

On the **1st-order** mesh, 0.0005 m completes but the first eigenvalues are negative (post-buckled).
The chosen base is **0.0001 m** (paper case, symmetric patches): 10 modes in near-degenerate pairs,

```
-1.58666e-4 / -1.58678e-4
-1.62237e-4 / -1.62241e-4
-3.99094e-4 / -3.99220e-4
-4.05422e-4 / -4.05469e-4
-5.56690e-4 / -5.57332e-4
```

Negative λ is accepted because only the **shapes** feed the imperfection (`ALLOW_NEGATIVE = True`
in `aba_mode_from_odb.py`). The pair split is 0.005–0.1 %, i.e. the left/right symmetry is recovered.

### Element order: 1st order is the default (2026-10-20)

The paper's reference row is STRI65 + S8R5 (2nd order), but that path is not usable for mode extraction here.

| Run | Mesh | Prestress | Result |
|---|---|---|---|
| R1 | S4/S3 1st, 10,201 nodes | 1.8e-5 m, corner2 | **4 modes**, λ 2.3285–4.4945e-4 |
| R2 | S4/S3 1st | 1.8e-5 m, paper3 | 4 modes (see R2 log) |
| R3 | S8R5/STRI65 2nd, 30,501 nodes | 1.8e-5 m | 184 negative eigenvalues / 177,421 negative diagonals -> **0 modes** |
| R3b | S8R5/STRI65 2nd | 0.0005 m | 592 negative eigenvalues / 166,009 negative diagonals -> **0 modes** |
| R4 | S4/S3 1st | 0.0005 m, paper | 10 modes but λ negative (post-buckled) |
| R5 | S4/S3 1st, symmetric patch | 0.0001 m, paper | **10 symmetric pairs** (current base) |

The 2nd-order generator remains in `aba_grid_mesh.py` (`--order 2`) for future experiments,
but the default is 1st order. The paper itself reports an **S3+S4, 10,100-element row**
(Table A.1: u_z,max 2.003e-4 m, 2x3 wrinkles), so the 1st-order comparison is a paper-reported
configuration, not an arbitrary deviation. Amplitude comparisons stay within the paper's
element-sensitivity band (80.1–99.0 t).

Set `PAPER_MATCH = False` in all three to return to the old MFBO baseline.

## Mode table (2026-10-20: paper + clamp_only)

The imperfection table is the **composite of the `paper` family (stored in the `seed` slot) and the
`clamp_only` family** (2+2 modes), not the pure paper-family table. In `run_abaqus.py`:

```
USE_COMPOSITE_MODE_TABLE = True   # False = pure paper family (modes_seed.txt)
```

`make_seed.py` runs both cases and combines them with `combine_mode_tables.py`:

- `seed` slot → **`paper` case**: cable-free, 3-vertex `paper3` static (2:1, 0.0001 m base state),
  buckle 0.10 m. Symmetric right/left patches.
- `clamp_only` case: clamp model, clamp pull `DISP * CLAMP_DC`, buckle clamp perturbation 0.10 m.

The cable `seed` case remains available for manual `run_abaqus_buckle.py` runs,
but it is **not** the mode-table source any more.

### Paper reproduction vs physical model (2026-10-20)

`PAPER_INCLUDE_CLAMPS` splits the two models:

| Model | `PAPER_INCLUDE_CLAMPS` | u3 policy (LF / HF) | Clamps | Vertex cables | Use |
|---|---|---|---|---|---|
| Paper 5PCS | `False` | LF: full u3=0 / HF: three edges only (free) | no | 2.69e7 Pa initial stress | Galhofo comparison |
| Physical | `True` | LF: full u3=0 / HF: three edges only (free) | yes, **initial stress** `CLAMP_PRESTRESS = d_c * CABLE_PRESTRESS` | 2.69e7 Pa initial stress | MFBO / post-buckling |

`NO_CLAMP = not PAPER_INCLUDE_CLAMPS`. `_HOLD_U3_ALL` is fidelity-based: LF always holds the
full membrane `u3=0`; HF paper is u3-free (three edges only); HF old disp holds then releases.
Both cable families are initial-stress; the clamp displacement pull is removed
(`CLAMP_PULL = CLAMP_FINAL = 0`). The composite mode table is the right basis for the physical model.

### Paper-faithful post-buckling: initial stress only (2026-10-20)

In the `cable_stress` path `_DISP_PRE = _DISP_FINAL = 0`, so **LF and HF** apply **no vertex
displacement**. The LF uses the same initial-stress prestress path and keeps the full-membrane
`u3=0` (paper path); its final state is `Step-GlobalTension`. The HF paper path has **no
full-membrane `u3=0` from Initial** (three edges only), so the wrinkles form in
`Step-GlobalTension`; `Step-Postbuckle` is a no-load settling step. The old disp HF path keeps
the old full hold and releases `BC_Stabilize_Z` in `Step-Postbuckle`. The load is the
cable/clamp initial stress only.


`PAPER_HOLD_U3_GLOBAL` is no longer a switch; `_HOLD_U3_ALL` is the fidelity-based policy.

- vertex cable initial stress **2.69e7 Pa**, ends fixed;
- clamp cable initial stress **`d_c × 2.69e7 Pa`**, ends fixed.

That is the Galhofo Fig. 11/17 path, with the clamps added as the physical design. The old MFBO
baseline (`LOAD_MODE='disp'`) keeps the 1 mm `GLOBAL_FINAL`; it is not the paper path.

## LF prestress and metrics (2026-10-21)

| Fidelity | Final step for metrics | Loads / BC | Output |
|---|---|---|---|
| LF (paper, `PAPER_INCLUDE_CLAMPS=True`) | `Step-GlobalTension` | cable 2.69e7 Pa + clamp `d_c*2.69e7 Pa`; full-membrane `u3=0`; no displacement | `RESULTS:lf1,lf2,lf3` + `ZK:` |
| LF (old baseline) | `Step-HighTension` | displacement-driven; preserved when `PAPER_MATCH=False` | `RESULTS:lf1,lf2,lf3` (+ `ZK:` diagnostic) |
| HF (paper) | `Step-Postbuckle` (last frame; no-load settling) | no full u3=0 from Initial; three edges only; wrinkles form in GlobalTension | `RESULTS:lf1,lf2,lf3,hf`, `DN:`, `ZK:` from the paired LF ODB |
| HF (old disp) | `Step-Postbuckle` (last frame) | full u3=0 hold -> release; final load ramp | same, archived as `HF_<tag>_u3free.odb` for paper |

### ZK (Zhang & Kiendl) definition

`ZK:` is a **separate line**; it never changes the `RESULTS` contract.

- `E` is the Green-Lagrange strain (`nlgeom=ON`) at the last frame of the LF state.
  Components are averaged over integration/section points per element before the principal
  values (membrane strain).
- Principal strains `E1 >= E2` in the shell local frame; `E1/E2` are rotation-invariant.
- State (area-weighted by `EVOL`): `taut` if `E2 > 0`; `wrinkled` if `E1 > 0 and E2 <= 0`;
  `slack` if `E1 <= 0`.
- `state_wrinkled_frac` = wrinkled area / total area.
- `zk_intensity_mean` / `zk_intensity_p95` = area-weighted mean / weighted p95 of `|E2|`
  over the wrinkled area (dimensionless strain).
- `zk_dir_deg` = area-weighted axial mean of the major principal strain direction `n1`
  in the **global x-y** plane, measured from +x, mod 180 deg. The element local 1-axis
  (S3: node1 -> node2; S4: average of the 1-2 and 4-3 edges) is used to rotate the
  shell-local `E` components. `nan` if there is no wrinkled area or the axial distribution
  is degenerate (`R < 1e-3`); the spread is in `zk_metrics.txt`.
- If `E` is missing: `ZK:FAIL reason=E_field_missing`; the LF/HF result is still emitted.
- Caveat: on a flat `u3=0` membrane, uniaxial tension gives `E2 < 0` by Poisson contraction,
  so the strain criterion can flag `wrinkled` without true in-plane compression. `lf1`
  (`S2 < 0`) is the complementary stress-based indicator.

`zk_metrics.txt` (in `code/aba` when run by `run_abaqus.py`) contains `status`, `step`, and
the diagnostics (`zk_state_taut_frac`, `zk_state_slack_frac`, `zk_dir_spread_deg`,
`zk_n_elem`, `zk_n_wrinkled`, `zk_area_total`, `weights`, `note`). Successful LF/HF runs
also append `zk_` keys to `results/results_xc..._dc....txt` (the CSV schema is unchanged).

#### `zk_Iw` sidecar (not the public `ZK:` line)

`zk_metrics.txt` additionally records the reference-area compressive strain `zk_Iw`:

```
I_w = (1/A_ref) [ Σ_wr    w_ref · max(|E2| - νE1, 0)
                + Σ_slack w_ref · (|E2| + ν|E1|) ]
```

- `w_ref` = element `EVOL` at `frameValue=0` (undeformed reference; uniform thickness
  makes it area-equivalent). If the reference frame is absent, the current `EVOL` or
  uniform weights are used and `zk_Iw_weights` records the fallback.
- `wrinkled` (`E1>0, E2<=0`) uses the **positive part** `max(|E2|-νE1, 0)`; a wrinkled
  element with `|E2| <= νE1` contributes zero (quasi-taut).
- `slack` (`E1<=0`) contributes `|E2| + ν|E1|`; taut contributes zero.
- `ν = 0.34` (`ZK_NU` in `eval_abaqus.py`).
- Sidecar keys: `zk_Iw`, `zk_Iw_p95`, `zk_Iw_area_total`, `zk_Iw_weights`,
  `zk_Iw_slack_share`, `zk_Iw_wr_share`, `zk_Iw_clip_zero_frac`,
  `zk_Iw_ref_missing`, `zk_Iw_nu`. `zk_Iw_slack_share` / `zk_Iw_wr_share` are the
  shares of the `I_w` numerator from slack / wrinkled elements; `zk_Iw_clip_zero_frac`
  is the fraction of the **wrinkled reference area** whose positive part is zero.
- `zk_Iw` is **not** `zk_intensity_mean`: the public intensity is wrinkled-only `|E2|`
  with current-`EVOL` weights; `I_w` covers the whole reference area, includes slack,
  applies the `ν` correction, and clips the positive part.
- `zk_Iw` is not yet an MFBO LF target; `mfbo.py` still uses `lf1`. No fixed
  `hf ≈ 2.9234·I_w` scale is assumed — a multi-fidelity GP learns the LF–HF relation.
- `python eval_lf_stub.py` (no Abaqus license) checks the synthetic `I_w`:
  reference weights `(1,2,2)`, wrinkled `9.8e-5`, slack `3.34e-4`, the
  `max(|E2|-νE1,0)` clip case, and that the public `ZK:` values are unchanged.

### hf / DN / u_z,ave definitions (unchanged)

- `hf`: `get_hf` integrates the SRP thrust over the last `Step-Postbuckle` frame and returns
  the **relative** loss `loss / (area_ref * IDEAL_THRUST_COEFF)`, where `area_ref` is the
  **undeformed** area at `frameValue=0` (not the deformed area). Normals are flipped to +z;
  `flip_bias.txt` records the backward-area percentage and the no-flip value.
- `DN:dn_deg`: area-weighted equivalent normal deviation from +z at the last `Step-Postbuckle`
  frame; `dn_psi_deg` is its azimuth.
- `u_z,ave`: `base_state_probe.py` on the `Step-Postbuckle` last frame,
  `(|u3,max| + |u3,min|)/2` (paper definition).

## Required run order (Windows, Abaqus)

The table is node-label based. The default mesh is 1st order (10,201 membrane nodes).
The composite table's fingerprint includes `MESH_ORDER`, so `run_abaqus.py` aborts before submitting
if the table was built on a different mesh or base state.

```bat
cd code

:: 0) one-time after switching to this branch: drop old displacement-driven LF ODBs
::    (the sidecar guard also detects them, but deleting is cleaner)
del /q aba\LF_Analysis*.odb

:: 1) regenerate seed + clamp_only families and combine (1st order)
python make_seed.py --x_c 0.5 --d_c 1.0 --force

:: 2) confirm the seed family and node count
findstr /C:"s4_grid_paper" seeds\modes_seed.txt
findstr /C:"nodes_per_mode=[10201]" seeds\modes_xc0p5.meta.txt

:: 3) LF run (same prestress path; Step-GlobalTension is the final state)
abaqus cae noGUI=run_abaqus.py -- LF 0.5 0.5
type aba\zk_metrics.txt

:: 4) HF run
abaqus cae noGUI=run_abaqus.py -- HF 0.5 0.5

:: 5) manual base-state / post-buckle probes (R-13 is also live; see below)
abaqus python base_state_probe.py aba\LF_Analysis_xc0p5_dc0p5.odb Step-GlobalTension 1.0 1.0
abaqus python base_state_probe.py aba\HF_Postbuckle_u3free.odb Step-GlobalTension 1.0 1.0
abaqus python base_state_probe.py aba\HF_Postbuckle_u3free.odb Step-Postbuckle 1.0 1.0
```

### R-13 base-state probe (live by default)

`run_abaqus.py` ships with `BASE_PROBE=True`. The automatic R-13/R-13:CT now
uses the **paired LF ODB first** (`lf_odb_is_current` sidecar); the fixed
`Buckle_Analysis.odb` is **never** used as an implicit fallback. If the LF ODB
is missing, R-13 skips instead of reading an untagged cross-point file.

- The MFBO objective/training data are `RESULTS` (`lf1-3`, `hf`) only. R-13 is a
  human diagnostic; it does not change the model, ODB, `hf`, `DN`, or `ZK`.
- LF is the HF prerequisite (`mfbo.py` runs it first); its `mtime < run start` is
  normal and the sidecar tag is the freshness proof.
- Set `BASE_PROBE=False` to disable the automatic probes for a production MFBO sweep.
- The point `results_*.txt` records `r13_ct_*`; the `all_results.csv` schema is unchanged.

Manual commands (same sources):

```bat
:: paired LF prestress state (same x_c,d_c; no imperfection -> flat u_z baseline)
abaqus python base_state_probe.py aba\LF_Analysis_xc0p5_dc0p5.odb Step-GlobalTension 1.0 1.0

:: exact HF base if the point-tagged archive exists (includes imperfection + RF)
abaqus python base_state_probe.py results\odb\HF_xc0p5_dc0p5_u3free.odb Step-GlobalTension 1.0 1.0

:: HF post-buckle: u_z,ave / compression map
abaqus python base_state_probe.py aba\HF_Postbuckle_u3free.odb Step-Postbuckle 1.0 1.0
```

For debugging only, `SOLARSAIL_BASE_ODB=<path>` overrides the R-13 source explicitly.
The fixed `aba\Buckle_Analysis.odb` has **no point tag** (and the self-buckle job is
`RUN_SELF_BUCKLE_JOB=False`), so do not use it as a per-point R-13 source.

If `CENTRE(mid)` is not 7000 Pa, scale:
`CABLE_PRESTRESS <- CABLE_PRESTRESS * 7000 / measured` and rerun (1–2 iterations).

## LF regression baseline (before the u3-free change)

The LF model must not change. Save the two LF points before/after the patch and compare
`extraction.txt` / `zk_metrics.txt` numerically (the ODB bytes change because LF now also
writes `RF`; the solve is unchanged).

```bat
cd code

:: before the patch (previous commit / git stash), for each point:
abaqus cae noGUI=run_abaqus.py -- LF 0.5 1.0
copy aba\extraction.txt results\lf_baseline_xc0p5_dc1.extraction.txt
copy aba\zk_metrics.txt results\lf_baseline_xc0p5_dc1.zk_metrics.txt

abaqus cae noGUI=run_abaqus.py -- LF 0.5 2.0
copy aba\extraction.txt results\lf_baseline_xc0p5_dc2.extraction.txt
copy aba\zk_metrics.txt results\lf_baseline_xc0p5_dc2.zk_metrics.txt

:: after the patch: remove the pre-patch LF ODBs (RF field was added), then rerun
:: the same two LF points and diff the two txt files
del /q aba\LF_Analysis*.odb
```

The LF keeps `Step-GlobalTension` `initialInc=0.01` / `maxNumInc=1000` and the old
GlobalTension stabilization `1.0`; only the `RF` field output was added.

## Imperfection amplitude sensitivity (HF u3-free)

The HF is free from the start, so the imperfection pattern—not just its role as a
tie-breaker—selects the wrinkle mode. Run the same point with
`IMPERFECTION_AMPL_T = 0.1` and `0.5` and compare `hf`, `u_z,ave`, and `DN`:

```bat
:: edit IMPERFECTION_AMPL_T in code/run_abaqus.py, then:
abaqus cae noGUI=run_abaqus.py -- HF 0.5 1.0
abaqus python base_state_probe.py aba\HF_Postbuckle_u3free.odb Step-Postbuckle 1.0 1.0
```

The mode table does not need regeneration (the amplitude is applied in `build_perturbation`).
If the result depends on the amplitude, it is **initial-condition dependence**, not physics.

## Verification done locally (no Abaqus)

> 2026-10-21 변경분 중 `check_model_consistency.py`와 `make_seed.py --dry-run`은 사용자 PC에서
> 확인했다. `buckle_stub_harness.py`는 첫 실행에서 mock 결함(`Assembly.Set` 이 `.nodes` 를
> 보존하지 않음 / `getClosest` 반환형)이 드러나 보정했으므로 **재실행**이 필요하다.
> **u3-free 변경분(LF u3=0 / HF free, GlobalTension 사다리, `HF_Postbuckle_u3free`)은 이
> 환경에서 실행하지 못했다.** 사용자 PC에서 LF baseline 2점 -> patch -> 재실행 비교,
> HF `0.5 1.0`, `IMPERFECTION_AMPL_T=0.1/0.5` 순으로 확인한다.
> 남은 명령:
> `python -m py_compile run_abaqus.py eval_abaqus.py aba_imperfection.py mfbo.py eval_lf_stub.py run_abaqus_buckle.py check_model_consistency.py buckle_stub_harness.py make_seed.py`,
> `python eval_lf_stub.py`, `python buckle_stub_harness.py paper 0.5`,
> `python buckle_stub_harness.py clamp_only 0.5`, `python buckle_stub_harness.py seed 0.5`.

- `py_compile` on all modified scripts: OK (이전 기준선).
- `check_model_consistency.py`: **PASS — HF↔buckle 21 / HF↔mode source 18** (2026-10-21 사용자 런).
  free-mesh SHARED 9항목 제거 후 예상치와 일치한다.
- `aba_grid_mesh.py --write --order 2` (optional 2nd-order generator): 30,501 nodes / 9,900 S8R5 + 200 STRI65;
  all checked mid-side nodes are exact edge midpoints; shared edges used twice, boundary edges once.
- `buckle_stub_harness.py paper|clamp_only|seed 0.5` (1st order): the paper trace is
  Initial (top fixed, corners 0, edges u3=0) → Static (top 0.0002, corners 0.0001) →
  Buckle (top 0.01, corners 0.01). The clamp_only/seed traces are in the harness `[배선]`
  output. Lambda/convergence are NOT verified by the harness.
  Mock fixes (2026-10-21): `Assembly.Set` keeps `.nodes`; `NodeArray.getClosest` returns
  `(node, distance)`. Rerun the harness after these fixes.
  The PULL_TOP harness (`buckle_pull_top_stub.py`) was deleted with the feature.
- `eval_lf_stub.py` (new; run `python3 eval_lf_stub.py`): checks that a fake ODB selects
  `Step-GlobalTension` when `Step-HighTension` is absent and `Step-HighTension` when present;
  checks lf1/lf2/lf3, ZK, `RESULTS:...`, `extraction.txt`, and
  `ZK:FAIL reason=E_field_missing` with a valid RESULTS line.
- `make_seed.py --x_c 0.5 --d_c 1.0 --dry-run`: selects `CASE='paper'` + `CASE='clamp_only'`,
  job tags `s4_grid_paper` / `s4_grid_clamp_only`, `_plv01000`, then combines 2+2.
  **2026-10-21 buckle refactor changed the model fingerprint/JOB_NAME keys** — regenerate
  the seed tables with `--force` before the first HF run (old tables will fail `seed_gate`).

## Known risks (not yet tested in Abaqus)

1. **Cable initial stress API**: `model.Stress(region=instance.sets['All_Elem'])` relies on part-set
   propagation to dependent instances. This is now used by **both LF and HF** (paper path).
   If Abaqus rejects it, fall back to `*INITIAL CONDITIONS, TYPE=STRESS, ELSET=...` in the keyword block.
2. **2nd-order orphan mesh import** (experimental, not used by default): `PartFromInputFile` with
   S8R5/STRI65 is verified at the `.inp` level only. The default is 1st order.
3. **Paper buckle base state**: 0.0005/0.001 m is unstable on the 2nd-order path (`0 modes`);
   on the 1st-order path 0.0005 completes but all λ are negative (post-buckled). The chosen base
   is **0.0001 m**, which gives 10 symmetric pairs (R5). λ is negative there, so
   `aba_mode_from_odb.py` has `ALLOW_NEGATIVE = True`; only the shapes are used.
   Do not change `MODE_PRESTRESS` without re-checking the pair symmetry and the λ sign.
4. **Stabilization ladder**: HF GlobalTension `(2e-4, 1e-3, 1e-2)` and Postbuckle
   `(5e-4, 2e-3, 1e-2)` are paired; `ALLSD/ALLIE` is read on **both** steps. >5% is recorded
   as `SUSPECT` (a fake, damping-supported solution); check `stabilize_used.txt`/CSV.
5. **LF is also 1st order** on this branch (same model).
6. **Thickness is still 5 µm** (paper 2.5 µm). Amplitude comparisons must be thickness-normalised.
7. **Stale LF ODB**: LF ODB names do not encode the model version. The LF reuse is guarded by
   a `.model.txt` sidecar (`tag=globaltension-v1-c1-u1` for the physical default; the tag
   includes the clamp/u3 sub-mode). Delete old `aba/LF_Analysis*.odb` once before switching
   to this branch; otherwise they are detected as stale and resubmitted (the `mfbo.py`
   prerequisite also checks the sidecar). Note: `mfbo.py` expects the physical-default tag;
   with `PAPER_MATCH=False` or a different clamp/u3 sub-mode it will resubmit LF
   (safe, but slower).
8. **ZK Poisson caveat**: the strain-based state can flag uniaxial tension as `wrinkled`
   (Poisson contraction, `E2 < 0`). Use `lf1` as the complementary stress-based indicator.
9. **Buckle seed tables after the refactor**: `JOB_NAME`/fingerprint keys changed
   (`BUCKLE_FINGERPRINT_KEYS` now excludes `USE_GRID_MESH`/`CLAMP_DIR`/`BUCKLE_SOLVER`).
   Regenerate `seeds/modes_*.txt` with `make_seed.py --force` before the first HF run.

## Removed in the 2026-10-21 buckle refactor

- `PULL_TOP`/`TOP_PULL_*` (3-tendon is now `CASE='paper'`/`paper3` only), `CLAMP_DIR`,
  `BUCKLE_SOLVER`/LANCZOS branch, `CLEAN_PERTURBATION`, free-mesh path,
  the unreachable nested `if CLAMP_LOAD=='none'` block (inside `if CLAMP_LOAD=='base'`),
  and the misleading `CLAMP_PERT` log/value.
- `buckle_pull_top_stub.py` is deleted (see below).
  `buckle_stub_harness.py paper|clamp_only|seed` remains the BC-trace harness.

## Deleted legacy files (2026-10-21)

These files had **no callers** and are superseded, switched off, or completed harnesses.
Remove them with the buckle refactor:

```bat
git rm code/run_abaqus_new.py code/solution_space.py code/solution_lf.py code/buckle_pull_top_stub.py code/release_edge_u3_stub.py code/hf_settings_stub.py code/dn_stats_stub.py
```

| File | Why |
|---|---|
| `run_abaqus_new.py` | legacy B-route (natural wrinkles, no buckle modes); `mfbo.py` submits `run_abaqus.py`; the checker now targets `run_abaqus.py` |
| `solution_space.py` / `solution_lf.py` | early exhaustive-sweep prototypes; superseded by `mfbo.py` + `sweep_*` |
| `buckle_pull_top_stub.py` | PULL_TOP removed; `paper3` is the 3-tendon source |
| `release_edge_u3_stub.py` | `RELEASE_EDGE_U3=False`; the switch stays in the model, the harness is retired |
| `hf_settings_stub.py` | §19.110 settings refactor is complete |
| `dn_stats_stub.py` | DN archiving wiring is complete (`DN:` / `dn_` keys verified) |

The other unreferenced files in the inventory are **not** deleted here; they need
case-by-case review because entry points, standalone calculators, and data generators
can legitimately have no imports.

## Not changed

- `run_abaqus_cable.py` (clamp-free oracle).
- Result CSV schema (`all_results.csv`) and the `RESULTS`/`extraction.txt` contracts.
  (The LF archive now records `fidelity=LF`; it was hardcoded `HF` before. The CSV columns are unchanged.)
- `mfbo.py`: only the LF ODB sidecar freshness check was added; the HF path/CSV schema is unchanged.
- Paper's 7000 Pa is a *result* at the quadrant centre, not a held load.
