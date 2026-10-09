# Stage 1 · Step 10 — Paper-literal model and lateral-side screen (pre-registered 2026-10-09)

Status: **pre-registered before R27–R30 were solved.** All runs: pad 10, elastic, L2 mesh, T_sf = 277 °C.
Selection rule (unchanged): pad 10 heating 30–350 °C decides. Pad 3 is the independent check.

## 0. Basis
- For pad 10, ε′zz[1e-3] = −0.004466·(σlat − σz)[MPa], computed from the paper constants. The hydrostatic part gives exactly 0.
- To match the measurement: σlat − σz = **+268 MPa at 30 °C**, with a slope of **−1.043 MPa/K** (30–350 °C).
- R26 (no cap, sliding, σz ≈ 0) gives +147 MPa and −0.597 MPa/K, i.e. 55 % of the target.
- Half-gap threshold for the stopping rule: **σlat ≥ 207.7 MPa at 30 °C** (slope ≥ 0.820 MPa/K).

## 1. Analytic screen (no runs). `results/diagnostics/lateral_screen_analytic.csv`
- **Model.** A Cu column with σz = 0 sits in a hole in the oxide. The oxide's in-plane strain is set by Si:
  σlat = κ·Δε / (1/M_Cu + 1/K_hole).
  - M_Cu = 259.5 GPa: equibiaxial in-plane modulus of pad 10 with σz = 0.
  - K_hole = E_ox/(1+ν_ox) for the bare oxide. With a liner: + E_l/(1−ν_l²)·t/a.
  - Δε = ∫(α_Cu − α_Si) dT, from T_sf.
- **Calibration.** κ = 0.863 makes the formula reproduce R26 exactly; the uncalibrated formula gives 170.8 vs 147.4 MPa. The oxide carries 81 % of the lateral compliance.
- **Rigid-matrix bound** (σz = 0, oxide replaced by a rigid matrix at α_Si): σlat = 891 MPa (769 MPa after κ). This is 3.3× the target, so the target is physically reachable in principle.

| Factor (value) | Range and source | σlat 30 °C (MPa) | gap closed (30 °C / slope) |
|---|---|---|---|
| E_SiO2 = 60 GPa | 48–70 GPa for SiO2 films (literature table in nanoindentation study); 74.8 ± 3.3 GPa for PECVD SiOx (microbridge); 85 ± 4 GPa for CVD-TEOS (Carlotti) | 126.9 | −17 % / −19 % |
| E_SiO2 = 85 GPa | same | 168.2 | 17 % / 19 % |
| α_SiO2 0.24–0.56 ppm/K | thermal oxide 0.24 ± 0.02 ppm/K; bulk silica 0.56 | ≈ 0 (in-plane strain is set by Si) | second order, not run |
| Si α(T) | Okada & Tokumaru, J. Appl. Phys. 56, 314 (1984): 2.59 → 3.95 ppm/K (30 → 400 °C) | 140.2 | −6 % / −8 % |
| Cu α(T) | compilation values 16.8 → 19.4 ppm/K (30 → 400 °C); **to be replaced by Hahn 1970 / NIST SRM 736** | 160.7 | 11 % / 15 % |
| Cu + Si α(T) | — | 153.5 | 5 % / 7 % |
| Ta liner 5 / 10 / 25 nm, E 186 GPa | Ta bulk E 186 GPa, ν 0.34, α 6.3 ppm/K (handbook). Fig 1a shows a TaTaN liner but **no thickness**. 5–25 nm is the usual damascene barrier range, **not yet backed by a retrieved citation** | 160.3 / 172.7 / 207.0 | 11 / 21 / 49 % |
| liner 10 / 25 nm, E 300 GPa (nitride-like) | TiN films 100–450 GPa (analog; no TaN value retrieved) | 187.2 / 238.2 | 33 / 75 % |
| All non-liner factors, upper ends | E 85 + Cu&Si α(T) | 175.2 | 23 % / 26 % |
| liner 10/186 + E 85 + α(T) | — | 199.8 | 43 % / 49 % |
| liner 25/186 + E 85 + α(T) | — | 233.3 | 71 % / 80 % |

**Screen result.**
- Without a liner, no factor or combination reaches half the gap (best 23 %).
- The liner is the only lever that can. It does so only at the thick end of the range (≳ 20 nm at 186 GPa), or combined with the upper oxide modulus.
- Circular vs square cannot be estimated analytically, so it is run.

## 2. Runs (registry step 10)
| Run | Configuration | Purpose | Pre-registered expectation |
|---|---|---|---|
| R27 | circular pad, no cap, **bonded** sidewalls, no liner | paper-literal model (Ayoub §2: "300 nm single pad embedded in SiO2 … 750 µm Si … symmetric BCs on two sides") | ε′zz(30 °C) **positive**, ≈ +0.9…+1.0e-3 (like R24 +0.987). It will **not** reproduce Ayoub's FEM (−1.181e-3). |
| R28 | R26 + circular pad | shape factor | within ±10 % of R26 (−0.657e-3) |
| R29 | R26 + Ta liner 10 nm, 186 GPa | liner, mid range | ε′zz(30) ≈ −0.77e-3 (21 %) |
| R30 | R26 + Ta liner 25 nm, 186 GPa | liner, thick end | ε′zz(30) ≈ −0.92e-3 (49 %) |

**Equivalences used for R27 instead of a literal 750 µm quarter model.**
- **Substrate.** A 750 µm wafer under a ~1.2 µm film has in-plane strain α_Si·ΔT to within ~4σ_f·t_f/(E_s·t_s) ≈ 5e-6. That is 0.2 % of the 3e-3 mismatch, so the 2 µm slab with α_Si-displaced walls is equivalent.
- **Symmetry.** The cell walls are mirror planes: x = 0 and y = 0 have the normal displacement fixed, and the opposite walls get a uniform normal displacement. This matches symmetric BCs.
  - Pad 10's anisotropy breaks exact mirror symmetry, but the shear stresses this causes are small (|σxz| = 4 MPa against σlat = 260 MPa in R01).
  - The 2.9 µm pitch couples neighbouring pads by only ~(a/r)² ≈ 1 %.

**Liner implementation.**
- One element ring of thickness t outside the 0.3 µm Cu (Cu size unchanged), plus one layer of thickness t under the trench bottom.
- Isotropic Ta properties.
- In R29/R30 the sliding interface is Cu/liner. The liner is bonded to the oxide.

**Extra output (R27–R30 only; post-processing, no physics change).**
- `hist_layers.txt` gives ε′zz in the top Cu element layer and in five equal height fifths.
- Purpose: test whether a point- or surface-level quantity, rather than a Cu volume average, could explain Ayoub's FEM curve.
- The decks of R01–R26 are byte-identical to before; all 26 physics hashes were re-checked.

## 3. Stopping rule (fixed before runs)
- If no single factor, and no combination evaluated within the cited ranges (single-factor FEM results plus the analytic non-liner terms), closes ≥ 50 % of the gap (σlat ≥ 207.7 MPa at 30 °C), **stop**.
- Then write the short summary for the professor: what was tested, the numbers, what the paper leaves unstated, and the question for the authors.
- Because the liner thickness has no retrieved citation, a liner-only "pass" is reported as **conditional on the liner thickness**, which must then be confirmed (TEM, or ask the authors).

## 4. Outcome (added after the runs; sections 0–3 are unchanged)
| Run | predicted ε′zz(30) | outcome ε′zz 30 / 400 °C | slope (1e-6/K) | gap closed |
|---|---|---|---|---|
| R27 | +0.9…+1.0 (does not reproduce Ayoub FEM) | **+1.143** / −0.569 | −4.63 | — (wrong sign) |
| R28 | R26 ± 10 % | −0.677 / +0.337 | 2.74 | 4 % |
| R29 | −0.77 (21 %) | −0.700 / +0.348 | 2.83 | 8 % |
| R30 | −0.92 (49 %) | −0.755 / +0.376 | 3.06 | 18 % |

- **R27.** Positive in every height band (top element layer +1.23; fifths +0.62…+1.58e-3). So the evaluation location cannot explain Ayoub's FEM curve.
- **Liner.** The analytic ring model over-predicted the liner effect by about 2.7× (FEM increment / analytic increment = 0.38 for 10 nm and 0.37 for 25 nm).
- **Best combination.** FEM circle (4 %) + FEM liner 25 nm (18 %) + analytic E_SiO2 85 GPa (17 %, or 6 % after the 0.37 correction) + α(T) (5 %) gives **33–44 % < 50 %**.
- **→ Stopping rule triggered; lateral-factor search stopped.**
