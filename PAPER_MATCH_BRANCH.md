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
| Prescribed displacement | static 0.001 / 0.0005 m (2:1), buckle 0.10 m (all five) | **same** (`run_abaqus_buckle.py`, `CASE='paper'`) |
| BC | three edges u3 = 0 | **three edges u3 = 0 from Initial** (HF) |
| Stabilization | 2e-4 | **2e-4 fixed** (no 0.005~0.05 ladder) |
| Clamps | none (5PCS) | `NO_CLAMP = True` enforced |
| Imperfection | first 4 modes, 0.10 t | same (`IMPERFECTION_AMPL_T = 0.1`) |

## Constants (edit only these)

`code/run_abaqus.py`:
```python
PAPER_MATCH = True
PAPER_MESH_ORDER = 2        # 1 = old S4/S3, 2 = paper S8R5/STRI65
CABLE_PRESTRESS = 2.69e7    # Pa; calibrate to 7000 Pa at CENTRE
PAPER_STABILIZE = 2.0e-4
```

`code/run_abaqus_buckle.py`:
```python
PAPER_MATCH = True
PAPER_MESH_ORDER = 2
# derived: CASE='paper', paper3, DISP_GLOBAL=0.0005, DISP_TOP_OVER_CORNER=2.0,
#          PERTURBATION=0.10, MODE_STABILIZATION=2e-4
```

`code/run_abaqus_mode.py`:
```python
PAPER_MATCH = True
PAPER_MESH_ORDER = 2
# mode source keeps DISP_GLOBAL = 1.8e-5 (low pre-stress) on purpose;
# modes are only the imperfection shape basis.
```

Set `PAPER_MATCH = False` in all three to return to the old MFBO baseline.

## Required run order (Windows, Abaqus)

The seed table is **node-label based**. A 1st-order table cannot be used with the 2nd-order mesh.
`run_abaqus.py` aborts before submitting if the table's source ODB name does not contain `grid2`.

```bat
cd code

:: 1) regenerate the seed family from the cable-free paper case (2nd order)
python make_seed.py --x_c 0.5 --d_c 1.0 --force

:: 2) confirm the table came from the 2nd-order paper run
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
- `make_seed.py --x_c 0.5 --d_c 1.0 --dry-run`: selects `CASE='paper'`, job tag `s4_grid2_paper`, `_plv10000`.

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
