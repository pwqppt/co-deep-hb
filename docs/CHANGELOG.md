# CHANGELOG

Newest first. Every entry: what changed and why. Bug entries add what was wrong, how it was found, and how it was fixed.

## 2026-10-09 — Step 11: Fig 2a/2c digitized, pads 4/5/7/8/9 added, orientation test R31–R34 registered

- **Data.** `inputs/ayoub2022_fig2ac_digitized.json` contains Fig. 2(a) and 2(c) for all 10 pads, with QA overlays in `docs/history/`.
  - Fig. 2(c) is the change from 30 °C to 400 °C **after the 2 h hold**.
  - Both figures use the Table 1 E values on the x axis.
- **Pipeline.**
  - `PADS` now includes pads 4, 5, 7, 8, 9 (Table 1 hkl, T_sf 277 °C assumption; pad 8 is a placeholder).
  - Registry step 11: R31–R34, the R26 configuration for elastic pads 4, 5, 7, 9.
  - New `src/tools/step11_compare.py` evaluates the pre-registered criterion.
- **Safety.** All R01–R30 deck physics hashes are unchanged. R31–R34 differ from R26 only in the Cu stiffness lines.
- **Why.** This tests whether the pad-10 shortfall is one orientation-independent factor, which is the precondition for Option B. Pre-registered in `docs/stage1_step11_orientation_test.md`.

## 2026-10-09 — Step 10: paper-literal model and lateral-side screen (new capability, new runs R27–R30)

- **What.**
  - `s1_mesh.build_mesh(..., t_liner=0.0)`: optional barrier-liner ring (MAT 5) around the Cu plus a trench-bottom layer.
  - `stage1_pipeline.py`:
    - liner material (Ta handbook values);
    - optional `hist_layers.txt` with ε′zz per height band (gated by `layers=True`);
    - registry step 10: R27–R30.
- **Why.** No cap/sidewall combination reproduces the magnitude of pad 10's ε′zz (best: R26, 55 %). The missing ~1.8× must come from the lateral side; see `docs/stage1_lateral_screen.md` (pre-registered).
- **Safety.** With `t_liner = 0` and `layers` unset, the mesh and deck are unchanged. The physics hashes of all R01–R26 decks were re-generated and are identical.
- Outcome (R27–R30 solved 2026-10-09): see `docs/stage1_lateral_screen.md` §4. The stopping rule was triggered.
  Added `docs/stage1_summary_for_PI_ko.md`, `docs/stage1_optionB_assessment.md`, `results/diagnostics/step10_outcomes_pad10.csv` and the figure.
- Analytic screen output: `results/diagnostics/lateral_screen_analytic.csv`. The pad-10 diagnostic tables are in `results/diagnostics/`.

## 2026-10-09 — Repository created (structure and paths only)

- **What.** The Stage 1 work was moved out of the Claude Science chat workspace into this repository.
  - Layout: `src/` (pipeline, mesh, validation, tools), `inputs/`, `configs/`, `results/<run>/`, `docs/`.
  - All paths are now relative to the repo root. `stage1_pipeline.py` finds the root from its own
    location.
  - Per-PC settings (MAPDL path, `-np`) move to `configs/machine.json`, which is not in git; the
    defaults are unchanged.
  - `run_cu_bar_test.py` is renamed `src/validation/g1_cu_bar_test.py`.
- **Why.** Multi-session and multi-PC work; nothing may depend on a chat workspace path.
- **No model physics changed and no new cases were added.**
  - The decks generated from the new layout are byte-identical to the solved R01, R24, R25 and R26
    decks, and their physics hashes are equal.
  - Regression check: re-run R01 and R26 from the repo layout. Required, exactly: R01 u = 2.083 nm,
    ε′zz(400) = −0.585e-3; R26 u = 6.093 nm, ε′zz(400) = +0.327e-3.
  - Reference values are in `docs/history/regression_reference_pre_restructure.json`; check with
    `src/tools/check_run.py`.
- **Added.**
  - `src/tools/check_run.py`: compare runs against the pre-restructure reference or the committed
    `HEAD`.
  - `README.md`, `.gitignore` (excludes `.rst` and other MAPDL scratch, generated decks, PDFs,
    `.venv`, `configs/machine.json`), `requirements.txt`, `docs/REFERENCES.md`, `docs/GIT_SETUP.md`.

## 2026-10-09 — Bug fix: "sliding" sidewall was partly bonded (affects R10, R25, R26)

- **What was wrong.**
  - The frictionless sidewall variant is supposed to tie only the normal displacement (no
    separation) between Cu and SiO2.
  - `s1_mesh.build_mesh(sliding=True)` duplicated only the *interior* sidewall nodes. It left 100 Cu
    nodes shared with SiO2, and a shared node is bonded in every DOF, including UZ:
    - the 40 nodes of the Cu top perimeter ring;
    - the 60 nodes along the four vertical corner edges.
  - The constraint equations themselves were correct (UX/UY only).
- **How it was found.**
  - The diagnostic R26 (no cap + "sliding") violated equilibrium. With a free top and frictionless
    sides, the vertical force across any Cu section must be zero. Yet R26 gave:
    - top-layer σz = −90 MPa and mid-height σz = −176 MPa at 400 °C (section force −15.8 µN);
    - pad lengthening 4.17 nm, below the free 5.19 nm;
    - ε′zz(400) = −0.330e-3, the wrong sign.
  - Auditing the R26 deck for every Cu node shared with another material found exactly 121 bottom,
    40 top-rim and 60 corner-edge nodes. No D, CP or wall condition touched Cu. The wall
    displacement used α_Si·(T − T_sf)·L with T_sf = TREF = REFT.
- **How it was fixed.**
  - Corner-edge nodes are duplicated and tied in UX and UY; UZ is free.
  - Without a cap, the top rim is also duplicated (normal tie; corners UX + UY).
  - With a bonded cap, the rim stays shared, because the cap physically joins Cu and SiO2 there.
  - The bottom rim stays shared: the pad bottom is bonded by design.
  - The probes and the Cu top-face average follow the Cu-side copies.
  - Pre-run audit: R26 Cu shares only its 121 bottom nodes (704 CEs); R25 has 660 CEs.
- **Verification.** The fixed R26 passes the equilibrium acceptance set before the run: Cu σz
  volume average ≈ 0 (−4e-9 MPa), mid-height ≈ 0, pad lengthening 5.649 nm > 5.19 nm,
  ε′zz(400) = +0.327e-3 > 0. Bonded decks (R01, R24) were unchanged by the fix.
- **Results invalidated.** The pre-fix R25 and R26 results were overwritten by the fixed runs. R10
  (plastic sliding) had not been run.

## 2026-10-09 — Convergence flag, vertical force split, diagnostics R24–R26

- **Bug: R01 was reported "converged False" although it solved.**
  - `*GET,…,ACTIVE,0,SOLU,CNVG` returns 0 for linear (single-iteration) solutions. R01's `.mntr`
    showed every load step completed in 1 substep and 1 iteration, and the `.out` had no errors.
  - Fix: the verdict now comes from `.mntr` (all load steps present) plus failure strings in
    `.out`. The raw flag is kept as `cnvg_flag_raw`.
- **Added `hist_force.txt`:** area-weighted Cu σz in the top, middle and bottom Cu element layers,
  which splits the vertical restraint between cap and sidewalls. Post-processing only.
- **Added diagnostic runs R24–R26** (pad 10 elastic: no cap / sliding / both).
  - Why: R01 gave ε′zz mirrored relative to the measurement (+1.175 / −0.585e-3 at 30 / 400 °C vs
    −1.197 / +0.600), with Cu σz = 2·σx,y.
  - Releasing σz at R01's lateral state reproduced the measurement (−1.167 / +0.581e-3), and the
    frame and sign were verified. The baseline was not changed.
- Added `--recollect`, which re-parses finished runs without solving.

## 2026-10-08 — Stage 1 pipeline (Phase B-1) and inputs

- `stage1_pipeline.py` and `s1_mesh.py`:
  - structured hex mesh, levels L0–L3, all below 128k nodes/elements;
  - Cu as TB,ANEL rotated to the global frame (as in G3 and cp_check) plus MISO from Ludwick;
  - load step 0 (T_sf → 30 °C), then the measured history;
  - batch and PyMAPDL runners;
  - run registry R01–R23.
- Mapping rule on computed E_z with pointwise flow-curve interpolation; G4a/G4b definitions;
  stress-free temperatures. These are decisions from the planning document, implemented per
  `docs/stage1_model_spec.md` §B.
- Inputs:
  - `inputs/ayoub2022_eps_dev_digitized.json` (Figs. 2b, 3a, 3b);
  - `inputs/ayoub2022_fig5_digitized.json`;
  - `inputs/cu_elastic_constants_vs_T.csv` (Chang & Himmel, UCRL-16697, Table I, report p. 14;
    G_Hill(400)/G(30) = 0.8545).
- `*GET` into array elements replaced by a scalar then assignment; the round-off threshold for
  the rotated stiffness was tightened (an E_z check failed at 4e-9 relative).

## 2026-10-08 — Validation G0, G1, G3

- G0: MAPDL cannot be licensed inside the app's sandbox (AppContainer). All runs are native.
- G1: Cu bar thermal expansion in µMKS, exact match (`src/validation/g1_cu_bar_test.py`).
- G3: anisotropic elasticity and orientation, ESYS and pre-rotated routes, 12/12 within 1 %
  (`src/validation/g3_anisotropic_cu_check.py`).
