# HANDOFF — run status (updated 2026-10-09)

**This file records run status only.** Design decisions (criteria, model choices, scope) are owned by
the main planning document. If a decision seems to have changed, flag it to the project lead; do not
write it here.

## Session-start routine (a new chat reads these, in this order)

1. **Main planning document** (all design decisions):
   https://claude.ai/code/artifact/94558db6-f844-4507-a2ad-f1d335461aeb
   If the chat cannot open the link, the lead attaches or pastes the document.
2. **This file**, `docs/HANDOFF.md`, for run status and next commands.
3. **The repository** on the PC in use: `$HOME\co-deep-hb` (this PC: `C:\Users\USER\co-deep-hb`),
   remote `https://github.com/pwqppt/co-deep-hb.git`.
   - In Claude Science, request host access to that folder.
   - Then `git pull` and read `docs/CHANGELOG.md` (latest entries) and `docs/RUNLIST.md`.
4. Execution rule: ANSYS cannot run inside the app.
   - Every run is handed to the lead as **one** paste-ready PowerShell block (with the `cd`),
     together with the expected runtime and the files it creates.
   - The lead runs it on a PC with ANSYS and reports back.

## Run status

| Run / gate | What | Status | Key numbers (end of heating, 400 °C) |
|---|---|---|---|
| G0 | MAPDL connectivity from the app | done — runs only natively | `docs/history/G0_mapdl_connectivity_report.md` |
| G1 | Cu bar thermal expansion, units | PASS | ΔL 6.1875 µm, σ −680.625 MPa (exact) — `results/validation_G1_cu_bar/` |
| G3 | anisotropic Cu + orientation | PASS (12/12) | `results/validation_G3_anisotropy/` |
| cp_check | ANEL + MISO (Ludwick) syntax | done | `results/validation_cp_check/` |
| R01 | pad 10 elastic, cap, bonded | done | u 2.083 nm, ε′zz −0.585e-3, Cu σz −261 MPa |
| R24 | DIAG: no cap, bonded | done | u 2.817 nm, ε′zz −0.492e-3 |
| R25 | DIAG: cap, sliding (fixed build) | done | u 2.539 nm, ε′zz −0.478e-3, Cu σz −227 MPa |
| R26 | DIAG: no cap, sliding (fixed build) | done; equilibrium acceptance **PASS** | u 6.093 nm, ε′zz +0.327e-3, Cu σz ≈ 0, pad lengthening 5.649 nm |
| R27 | LAT paper-literal: circle, no cap, bonded, no liner | done | ε′zz 30/400 +1.143/−0.569e-3 (opposite sign to measured/Ayoub FEM), u 2.572 nm |
| R28 | LAT R26 + circle | done | ε′zz −0.677/+0.337e-3, gap closed 4 % |
| R29 | LAT R26 + Ta liner 10 nm | done | ε′zz −0.700/+0.348e-3, gap closed 8 % |
| R30 | LAT R26 + Ta liner 25 nm | done | ε′zz −0.755/+0.376e-3, gap closed 18 % |
| R31 | ORI pad 4 elastic, R26 config | done | Δε′zz 1.078 vs Fig2c 2.363 (r 0.456, FAIL) |
| R32 | ORI pad 5 elastic, R26 config | done | Δε′zz 0.989 vs Fig2c 1.854 (r 0.534, pass) |
| R33 | ORI pad 7 elastic, R26 config | done | Δε′zz 1.105 vs Fig2c 2.406 (r 0.459, FAIL) |
| R34 | ORI pad 9 elastic, R26 config | done | Δε′zz 1.103 vs Fig2c 1.854 (r 0.595, FAIL) |
| R02–R07, R10–R23 | see `docs/RUNLIST.md` | not run | — |
| R08 | pad 10 plastic = elastic | no solve needed (deck identical to R01) | — |
| R09 | Cu CTE(T) variant | blocked: needs `inputs/cu_cte_vs_T.csv` (NIST) | — |

**Pending analysis (not yet done):** compare R01 (bonded) and the fixed R25 (sliding, cap kept)
against measured pad 10 ε′zz on heating 30–350 °C. Report both; do not pick.

## Next commands

Restructure regression check R01/R26: **ALL MATCH** (2026-10-09); first commit/push: command issued, completion not yet confirmed.
Step 10 done: pre-registered stopping rule TRIGGERED (best combination 33–44 % < 50 %). Lateral-factor search stopped.
Step 11 (orientation test, pre-registered in docs/stage1_step11_orientation_test.md): R31–R34 done; **verdict FAIL** (pads 9, 7, 4) → Option B not justified, stopped. Without pad 9 (secondary): within tolerance.
Open: author question (docs/stage1_summary_for_PI_ko.md); Option B awaiting PI decision (docs/stage1_optionB_assessment.md) — not run.
Previously: the restructure regression check (R01 and R26 from the repo layout must
reproduce the numbers above exactly), then the first commit and push. See `docs/GIT_SETUP.md`.

## Open issues (run level)

- MAPDL `.out` files are ~13 MB, mostly the echo of 100k N/E lines. Wrapping the mesh block in
  `/NOPR` … `/GOPR` would shrink them, but it changes the deck text, so it is deferred to an explicit
  decision.
- The pipeline records `cnvg_flag_raw` from `*GET,…,SOLU,CNVG`, but that flag is 0 for linear solves.
  The verdict uses `.mntr` and the `.out` failure strings instead.
- `.rst` files are not in git; R01 and R25 `.rst` are kept on PC 1 only (in the old chat workspace).
- Modelling-level open items are tracked in `docs/stage1_model_spec.md` §B11 and the planning
  document.

## Session-end routine (when a gate passes or the lead says "wrap up")

1. Update this file: run status, next commands, open issues.
2. Add a `docs/CHANGELOG.md` entry: what changed and why. Bugs get what was wrong, how it was found,
   and how it was fixed.
3. Give the lead one commit + push block:
   ```
   cd "$HOME\co-deep-hb"; git add -A; git status --short; git commit -m "<summary>"; git push
   ```
4. At a model freeze (e.g. the end of Stage 1), tag it:
   ```
   cd "$HOME\co-deep-hb"; git tag -a v1.0 -m "Stage 1 model freeze"; git push origin v1.0
   ```
