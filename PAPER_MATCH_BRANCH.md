# Galhofo 2022 paper-match branch

Branch: `feature/galhofo-paper-match`
Base: `af89de7` (refactor/buckle-case-clamp-lf)

This branch adds a **paper-matching preset** for Galhofo et al., *Thin-Walled Structures* 181 (2022) 110007.
It is a research branch: the MFBO baseline is preserved behind `PAPER_MATCH = False`.

## What changes when `PAPER_MATCH = True`

| Item | Paper | This branch |
|---|---|---|
| Membrane elements | STRI65 + S8R5 (2nd order), 10,100 | **S8R5 9,900 + STRI65 200**, 30,501 nodes |
| Pre-tension load | cable initial stress 2.69e7 Pa (5PCS), cable ends fixed | **3 cables, sigma11 = 2.69e7 Pa**, ends fixed, no displacement ramp |
| Prescribed displacement | static 0.001 / 0.0005 m (2:1), buckle 0.10 m (all five) | buckle 0.10 m **same**; static **shape** 2:1 but magnitude **1.8e-5 m** for mode extraction |
| BC | three edges u3 = 0 | **three edges u3 = 0 from Initial** (HF) |
| Stabilization | 2e-4 | **2e-4 fixed** (no 0.005~0.05 ladder) |
| Clamps | none (5PCS) | `NO_CLAMP = True` enforced |
| Imperfection | first 4 modes, 0.10 t | same; mode table = **seed + clamp_only composite** (2+2) |

## Constants (edit only these)

`code/run_abaqus.py`:
```python
PAPER_MATCH = True
PAPER_MESH_ORDER = 2        # 1 = old S4/S3, 2 = paper S8R5/STRI65
CABLE_PRESTRESS = 2.69e7    # Pa; calibrate to 7000 Pa at CENTRE
PAPER_STABILIZE = 2.0e-4
USE_COMPOSITE_MODE_TABLE = True   # seed + clamp_only 합성표
```

`code/run_abaqus_buckle.py`:
```python
PAPER_MATCH = True
PAPER_MESH_ORDER = 2
MODE_PRESTRESS = 1.8e-5     # actual mode-extraction base state
PAPER_STATIC_DISP = 0.0005  # paper reference only (unstable here)
# derived: CASE='paper', paper3, DISP_TOP_OVER_CORNER=2.0,
#          PERTURBATION=0.10, MODE_STABILIZATION=2e-4
```

`code/run_abaqus_mode.py`:
```python
PAPER_MATCH = True
PAPER_MESH_ORDER = 2
# mode source keeps DISP_GLOBAL = 1.8e-5 (low pre-stress) on purpose;
# modes are only the imperfection shape basis.
```

### Why the static step uses 1.8e-5 m, not 0.0005 m

The paper's 0.0005 m corner displacement puts this model's base state above its bifurcation load.
Measured on this branch (S8R5/STRI65, 202,050 equations): 592 negative eigenvalues in the system matrix,
166,009 negative diagonal entries, and the subspace iteration ends with
`THE EIGENVALUES CANNOT BE FOUND ... INSTABILITIES IN THE BASE STATE` (0 modes).
The same failure was already recorded at 1e-3 m. Modes are only the imperfection shape basis;
HF sets the amplitude at 0.10 t, so the mode-extraction prestress does not need to be the paper's.
The 2:1 shape and the 0.10 m buckle step are kept.

Set `PAPER_MATCH = False` in all three to return to the old MFBO baseline.

## Mode table (user decision 2026-10-20)

The imperfection table is the **composite of the `seed` family and the `clamp_only` family** (2+2 modes),
not the pure seed-family table. In `run_abaqus.py`:

```
USE_COMPOSITE_MODE_TABLE = True   # False = old pure seed family (modes_seed.txt)
```

`make_seed.py` runs both cases and combines them with `combine_mode_tables.py`:

- `seed` case: cable model, 3-vertex `paper3` static (0.001 / 0.0005 m), buckle 0.10 m.
- `clamp_only` case: clamp model, clamp pull `DISP * CLAMP_DC`, buckle clamp perturbation 0.10 m.

The `paper` case (cable-free, clamp-free) remains available for manual `run_abaqus_buckle.py` runs,
but it is **not** one of the two mode-table families.

## Required run order (Windows, Abaqus)

The table is node-label based. A 1st-order table cannot be used with the 2nd-order mesh.
The composite table's fingerprint includes `MESH_ORDER`, so `run_abaqus.py` aborts before submitting
if the table was built on a different mesh or base state.

```bat
cd code

:: 1) regenerate seed + clamp_only families and combine (2nd order)
python make_seed.py --x_c 0.5 --d_c 1.0 --force

:: 2) confirm the seed family came from a 2nd-order run
findstr /C:"grid2" seeds\modes_seed.txt

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
- `aba_grid_mesh.py --write --order 2`: 30,501 nodes / 9,900 S8R5 + 200 STRI65;
  all checked mid-side nodes are exact edge midpoints; shared edges used twice, boundary edges once.
- `buckle_stub_harness.py paper 0.5`: BC trace is
  Initial (top fixed, corners 0, edges u3=0) → Static (top 0.001, corners 0.0005) →
  Buckle (top 0.10, corners 0.10). Lambda/convergence are NOT verified by the harness.
- `make_seed.py --x_c 0.5 --d_c 1.0 --dry-run`: selects `CASE='seed'` + `CASE='clamp_only'`,
  job tags `s4_grid2_seed` / `s4_grid2_clamp_only`, `_plv10000`, then combines 2+2.

## Known risks (not yet tested in Abaqus)

1. **Cable initial stress API**: `model.Stress(region=instance.sets['All_Elem'])` relies on part-set
   propagation to dependent instances. If Abaqus rejects it, fall back to
   `*INITIAL CONDITIONS, TYPE=STRESS, ELSET=...` in the keyword block.
2. **2nd-order orphan mesh import**: `PartFromInputFile` with S8R5/STRI65 is verified at the `.inp`
   level only. Check `buckle/grid_Membrane.inp` if the CAE import fails.
3. **Paper buckle base state**: 0.0005/0.001 m may be unstable for mode extraction (`CONVERGED=0`).
   The mode source deliberately stays at 1.8e-5 m. If the buckle script fails, lower its
   `DISP_GLOBAL` and regenerate the seed table.
4. **No stabilization ladder**: `(2e-4,)` only. A failure yields `RESULTS:FAIL` by design;
   extend to `(2e-4, 5e-4, 1e-3)` if needed, and record which level completed.
5. **LF is also 2nd order** on this branch (same model), so LF/HF correlation changes.
6. **Thickness is still 5 µm** (paper 2.5 µm). Amplitude comparisons must be thickness-normalised.

## Not changed

- `run_abaqus_cable.py` (clamp-free oracle).
- MFBO driver (`mfbo.py`) and result CSV schema.
- Paper's 7000 Pa is a *result* at the quadrant centre, not a held load.
