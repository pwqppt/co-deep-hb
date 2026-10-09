# Stage 1 · Step 11 — Orientation test of a single lumped factor (pre-registered 2026-10-09, BEFORE R31–R34)

**Question.** Is the pad-10 shortfall of the R26 configuration (no cap, frictionless no-separation sidewalls; 55 % of the measured magnitude) one orientation-independent factor? If yes, a single lumped lateral-constraint parameter (Option B) is justified. If no, Option B is not justified and we stop.

## 1. Data (step 1, no runs)
- `inputs/ayoub2022_fig2ac_digitized.json`: Fig. 2(a) ε′zz at 30 °C and Fig. 2(c) Δε′zz for all 10 pads.
- QA overlays: `docs/history/qa_fig2a.png`, `docs/history/qa_fig2c.png`.
- **Fig. 2(c) is Δ between 30 °C and 400 °C *after the 2 h hold*.** The model quantity is therefore ε′zz(load step 8, end of hold) − ε′zz(load step 0, 30 °C). For elastic pads this equals the heating-only Δ. Check on pad 10: Fig 2c gives 1.806, Fig 3a heating gives 1.797.
- The Fig. 2 x-axis uses the **Table 1 E values**; the model uses computed E_z.
- Pads 4, 5, 7, 8, 9 were added to the pipeline from Table 1 hkl:

| pad | hkl (Table 1) | Table E | computed E_z (GPa) | class |
|---|---|---|---|---|
| 4 | 1 0.9 0.6 | 169 | 169.39 | elastic |
| 5 | 1 0.89 0.87 | 188 | 188.35 | elastic |
| 7 | 1 0.85 0.56 | 164 | 164.64 | elastic |
| 8 | 0.7 1 0.2 | 132 | 123.18 | plastic (Table E misprint) |
| 9 | 1 0.56 1 | 172 | 164.18 | elastic (Table E misprint) |

## 2. Runs
R31 (pad 4), R32 (pad 5), R33 (pad 7), R34 (pad 9) use the R26 configuration: elastic, L2 mesh, T_sf = 277 °C. Pad 10 = R26, already solved.
- The decks differ from R26 **only in the rotated Cu stiffness** (TBDATA lines) and the title. Checked by line diff.
- R01–R30 decks are unchanged: all physics hashes re-checked.

## 3. Pre-registered criterion (primary, decisive)
r_p = Δ_model,p / Δ_meas,p for p ∈ {4, 5, 7, 9, 10}; r̄ = mean of the five.
tol_p = max(10 %, err_p / Δ_meas,p) (whichever is larger: ±10 % or the pad's own Fig 2c error bar).
**PASS iff |r_p / r̄ − 1| ≤ tol_p for all five pads** → single orientation-independent factor → Option B justified. Otherwise **FAIL** → Option B not justified → stop.
Evaluated by `src/tools/step11_compare.py`, which writes `results/diagnostics/step11_orientation_test.csv`.

**Report only (not in the criterion).**
- u = UZ(Cu top centre) − UZ(top surface, cell edge), 30 → 400 °C, per pad, compared with the Fig 5 point and the Fig 5 regression at the Table E. Fig 5 is Ayoub's FEM output, not a measurement.
- Model slope of u vs E_z over the elastic pads, compared with −0.0132 nm/GPa.
- Plastic pads 1, 2, 3, 6, 8: Fig 2c values with the analytic σz = 0 elastic estimate (table below). No runs.

## 4. Pre-registered prediction (analytic σz = 0 model, κ = 0.863 calibrated on R26)
| pad | Table E | E_z | Δ_meas (Fig 2c) | Δ_model (pred.) | r_p | r_p/r̄ − 1 | tol | pred. |
|---|---|---|---|---|---|---|---|---|
| 9 | 172 | 164.2 | 1.854 ± 0.20 | 1.143 | 0.616 | +16.7% | 10.8% | **fail** |
| 7 | 164 | 164.6 | 2.406 ± 0.20 | 1.139 | 0.474 | -10.3% | 10.0% | **fail** |
| 4 | 169 | 169.4 | 2.363 ± 0.20 | 1.107 | 0.468 | -11.3% | 10.0% | **fail** |
| 5 | 188 | 188.3 | 1.854 ± 0.20 | 0.992 | 0.535 | +1.4% | 11.0% | pass |
| 10 | 198 | 189.5 | 1.806 ± 0.20 | 0.986 | 0.546 | +3.4% | 11.0% | pass |

r̄ (pred.) = 0.528. **Predicted verdict: FAIL**, driven by pad 9 and, marginally, pad 4.
- Why: pads 7 and 9 have almost the same computed E_z (164.6 / 164.2 GPa), so the model gives almost the same Δ (1.139 / 1.143).
- The measurement, however, differs by 30 % (2.41 vs 1.85).
- The measured elastic-pad Δ falls into two levels (pads 7, 4 ≈ 2.4; pads 9, 5, 10 ≈ 1.83). The model is smooth in E_z.

**Known caveat, fixed now.** Pad 9's Table 1 hkl gives E_z = 164.2 GPa, while the Table and Fig. 2 put it at 172 GPa.
- A misprinted hkl for pad 9 would weaken the pad-9 comparison.
- Rule: the primary verdict uses all five pads as specified. A secondary result without pad 9 is reported but **does not override** the primary verdict.

Plastic pads, report only (elastic estimate; it should under-predict because plasticity adds deviatoric strain):
| pad | Table E | E_z | Δ_meas | Δ_model elastic (pred.) | r |
|---|---|---|---|---|---|
| 3 | 68.0 | 68.0 | 3.702 ± 0.20 | 2.484 | 0.671 |
| 6 | 102.0 | 102.3 | 2.843 ± 0.20 | 1.780 | 0.626 |
| 8 | 132.0 | 123.2 | 2.359 ± 0.22 | 1.507 | 0.639 |
| 1 | 127.0 | 127.3 | 2.683 ± 0.20 | 1.462 | 0.545 |
| 2 | 130.0 | 130.3 | 2.359 ± 0.20 | 1.430 | 0.606 |

## 5. What happens next
- **PASS:** pre-register step 12 (Option B) before any calibration run. One phenomenological lateral-constraint parameter, calibrated on pad 10 heating 30–350 °C only and labelled outside the SiO2 literature range. Validation without re-tuning:
  - elastic pads 4, 5, 7, 9 (Fig 2c);
  - pad 3 heating with plasticity ②;
  - the u-vs-E_z slope against Fig 5 within ±10 %.
- **FAIL:** stop. Add the result to the PI summary. The question for the authors stays as written.

## 6. Outcome (added after R31–R34; sections 1–5 unchanged)
| pad | run | E_z | Δ_model | Δ_meas | r_p | r_p/r̄ − 1 | tol | result | u model / Fig 5 (nm) |
|---|---|---|---|---|---|---|---|---|---|
| 9 | R34 | 164.2 | 1.103 | 1.854 | 0.595 | +14.9% | 10.8% | **fail** | 6.21 / 7.59 |
| 7 | R33 | 164.6 | 1.105 | 2.406 | 0.459 | -11.3% | 10.0% | **fail** | 6.21 / 7.69 |
| 4 | R31 | 169.4 | 1.078 | 2.363 | 0.456 | -11.9% | 10.0% | **fail** | 6.18 / 7.65 |
| 5 | R32 | 188.3 | 0.989 | 1.854 | 0.534 | +3.0% | 11.1% | pass | 6.10 / 7.34 |
| 10 | R26 | 189.5 | 0.984 | 1.806 | 0.545 | +5.2% | 11.0% | pass | 6.09 / 7.19 |

- r̄ = 0.518. **PRIMARY VERDICT: FAIL** (pads 9, 7, 4), as predicted.
  - Option B is **not justified**, so we stop.
  - FEM and the analytic prediction agree within 2.2 points on every deviation.
- Secondary result (pad 9 excluded; does not override): deviations of −7.9 / −8.5 / +7.0 / +9.3 %, all within tolerance.
  - The failure therefore depends on pad 9, whose Table 1 hkl (E_z 164.2) is inconsistent with its Table/Fig. 2 E (172).
- u (report only): model slope -0.0045 nm/GPa over the elastic pads, against −0.0132 for the Fig 5 fit. That is about one third. Model u = 6.09–6.21 nm vs Fig 5 7.19–7.69 nm.
- Checks: all four runs exit 0 and converged, global frame, Cu σz ≈ 0, pad lengthening 5.65–5.78 nm.
