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
| Prescribed displacement | static 0.001 / 0.0005 m (2:1), buckle 0.10 m (all five) | buckle 0.10 m **same**; static **shape** 2:1, magnitude **0.0001 m** (symmetric 10-mode base) |
| BC | three edges u3 = 0 | **three edges u3 = 0 from Initial** (HF) |
| Stabilization | 2e-4 | **2e-4 fixed** (no 0.005~0.05 ladder) |
| Clamps | none (5PCS) | paper reproduction: `NO_CLAMP = True`; physical model: `PAPER_INCLUDE_CLAMPS = True` (clamps + vertex cables) |
| Imperfection | first 4 modes, 0.10 t | same; mode table = **seed + clamp_only composite** (2+2) |

## Constants (edit only these)

`code/run_abaqus.py`:
```python
PAPER_MATCH = True
PAPER_MESH_ORDER = 1        # 1차 기본 (2차는 실험용)
CABLE_PRESTRESS = 2.69e7    # Pa; calibrate to 7000 Pa at CENTRE
PAPER_STABILIZE = 2.0e-4
USE_COMPOSITE_MODE_TABLE = True   # seed + clamp_only 합성표
```

`code/run_abaqus_buckle.py`:
```python
PAPER_MATCH = True
PAPER_MESH_ORDER = 1
MODE_PRESTRESS = 0.0001     # symmetric 10-mode base (2026-10-20)
PAPER_STATIC_DISP = 0.0005  # paper reference only (unstable here)
# derived: CASE='paper', paper3, DISP_TOP_OVER_CORNER=2.0,
#          PERTURBATION=0.10, MODE_STABILIZATION=2e-4
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
The 2:1 shape and the 0.10 m buckle step are kept.

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

| Model | `PAPER_INCLUDE_CLAMPS` | Clamps | Vertex cables | Use |
|---|---|---|---|---|
| Paper 5PCS | `False` | no | 2.69e7 Pa initial stress | Galhofo comparison |
| Physical | `True` | yes, `CLAMP_PULL = DISP_GLOBAL * d_c` | 2.69e7 Pa initial stress | MFBO / post-buckling |

`NO_CLAMP = not PAPER_INCLUDE_CLAMPS`. The paper-only model uses the 2e-4 stabilization;
the physical model uses the old `(0.005, 0.01, 0.02, 0.05)` ladder because the clamps cause cutbacks.
The composite mode table (paper + clamp_only) is the right basis for the physical model.

## Required run order (Windows, Abaqus)

The table is node-label based. The default mesh is 1st order (10,201 membrane nodes).
The composite table's fingerprint includes `MESH_ORDER`, so `run_abaqus.py` aborts before submitting
if the table was built on a different mesh or base state.

```bat
cd code

:: 1) regenerate seed + clamp_only families and combine (1st order)
python make_seed.py --x_c 0.5 --d_c 1.0 --force

:: 2) confirm the seed family and node count
findstr /C:"s4_grid_paper" seeds\modes_seed.txt
findstr /C:"nodes_per_mode=[10201]" seeds\modes_xc0p5.meta.txt

:: 3) HF run
abaqus cae noGUI=run_abaqus.py -- HF 0.5 1.0

:: 4) probe the pre-tension and post-buckle states
abaqus python base_state_probe.py HF_Postbuckle.odb Step-GlobalTension 1.0 1.0
abaqus python base_state_probe.py HF_Postbuckle.odb Step-Postbuckle 1.0 1.0
```

If `CENTRE(mid)` is not 7000 Pa, scale:
`CABLE_PRESTRESS <- CABLE_PRESTRESS * 7000 / measured` and rerun (1–2 iterations).

## Verification done locally (no Abaqus)

- `py_compile` on all modified scripts: OK.
- `check_model_consistency.py`: PASS (HF↔buckle 30, HF↔mode source 27).
- `aba_grid_mesh.py --write --order 2` (optional 2nd-order generator): 30,501 nodes / 9,900 S8R5 + 200 STRI65;
  all checked mid-side nodes are exact edge midpoints; shared edges used twice, boundary edges once.
- `buckle_stub_harness.py paper 0.5` (1st order): BC trace is
  Initial (top fixed, corners 0, edges u3=0) → Static (top 0.0002, corners 0.0001) →
  Buckle (top 0.10, corners 0.10). Lambda/convergence are NOT verified by the harness.
- `make_seed.py --x_c 0.5 --d_c 1.0 --dry-run`: selects `CASE='paper'` + `CASE='clamp_only'`,
  job tags `s4_grid_paper` / `s4_grid_clamp_only`, `_plv10000`, then combines 2+2.

## Known risks (not yet tested in Abaqus)

1. **Cable initial stress API**: `model.Stress(region=instance.sets['All_Elem'])` relies on part-set
   propagation to dependent instances. If Abaqus rejects it, fall back to
   `*INITIAL CONDITIONS, TYPE=STRESS, ELSET=...` in the keyword block.
2. **2nd-order orphan mesh import** (experimental, not used by default): `PartFromInputFile` with
   S8R5/STRI65 is verified at the `.inp` level only. The default is 1st order.
3. **Paper buckle base state**: 0.0005/0.001 m is unstable on the 2nd-order path (`0 modes`);
   on the 1st-order path 0.0005 completes but all λ are negative (post-buckled). The chosen base
   is **0.0001 m**, which gives 10 symmetric pairs (R5). λ is negative there, so
   `aba_mode_from_odb.py` has `ALLOW_NEGATIVE = True`; only the shapes are used.
   Do not change `MODE_PRESTRESS` without re-checking the pair symmetry and the λ sign.
4. **No stabilization ladder**: `(2e-4,)` only. A failure yields `RESULTS:FAIL` by design;
   extend to `(2e-4, 5e-4, 1e-3)` if needed, and record which level completed.
5. **LF is also 1st order** on this branch (same model).
6. **Thickness is still 5 µm** (paper 2.5 µm). Amplitude comparisons must be thickness-normalised.

## Not changed

- `run_abaqus_cable.py` (clamp-free oracle).
- MFBO driver (`mfbo.py`) and result CSV schema.
- Paper's 7000 Pa is a *result* at the quadrant centre, not a held load.
