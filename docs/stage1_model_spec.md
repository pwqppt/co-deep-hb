> **Status of this document.** Implementation spec for Stage 1, derived from the main planning document (link in `docs/HANDOFF.md`). Design decisions are owned by the planning document; if the two disagree, the planning document wins and the discrepancy is flagged to the project lead, not silently edited here.

# Stage 1 · Phase A — Model specification for reproducing Ayoub et al. (2022)

**Paper.** B. Ayoub, S. Moreau, S. Lhostis, H. Frémont, S. Mermoz, E. Souchier, E. Deloffre,
S. Escoubas, T. W. Cornelius, O. Thomas, "In-situ characterization of thermomechanical behavior of
copper nano-interconnect for 3D integration," *Microelectronic Engineering* **261**, 111809 (2022),
doi:10.1016/j.mee.2022.111809.

**Source used — substitution notice.** No PDF arrived with the request. The extraction below is
from the **authors' accepted version on HAL (hal-03672631v1, deposited 31 Jan 2023)**. Page numbers
"p.N" refer to that 7-page file: p.1 is the HAL cover sheet and the article runs pp.2–7. Content
should match the journal version, but the layout differs. If you have the publisher PDF, the only
values worth re-checking are Table 1 and Fig. 5.

**Status.** Phase A accepted on 2026-10-08 with the decisions in §B; B-1 accepted the same day. Decks and scripts are generated and dry-run. G(T) from Chang & Himmel is in place, so every run except R09 (NIST α(T)) can now be built. The first solves are running on your machine.

---

## 0. Findings that change the plan

1. **The paper's displacement values are simulations, not measurements.** *(Phase B: the measured ε′<sub>zz</sub> is now the primary target, G4a; Fig. 5 is the secondary implementation check, G4b.)* Laue microdiffraction
   only gives the *deviatoric* strain; the hydrostatic part "is not experimentally accessible"
   (§2, p.3). The "Total displacement in Z" of Fig. 5 (p.6) is COMSOL output. So G4 is a
   **model-to-model** check: our MAPDL vs. their COMSOL, both run with largely the same inputs.
   It validates our implementation, not the physics. The only experimental quantity is
   ε′<sub>zz</sub>(T) (Figs. 2–3). §3.3 below proposes a later, optional G4b against it.
2. **The pad is square, not circular.** The text says "diameter of 300 nm" (§2, p.2). Fig. 1(a)
   (p.3) is the only dimensioned drawing, and it shows a **0.3 µm × 0.3 µm × 0.85 µm** square prism
   with a **TaTaN** liner. That drawing is the basis for the geometry here.
3. **The paper's Cu constants differ from those used in G3.** Paper: C11 = 168.4, C12 = 121.4,
   C44 = 75.4 GPa (§2, p.3, Nye). G3 used 169.1 / 122.2 / 75.4. Stage 1 will use the **paper's**
   values. E<sub>z</sub> moves by ≤ 0.2 %: [100] 66.69, [110] 130.34, [111] 191.15 GPa. G3 still
   stands, because it validated the ANEL convention and the rotation method, not the numbers.
4. **Table 1 contradicts itself for 3 of 10 pads.** Recomputing eq. (1) with the paper's own
   constants reproduces pads 1–7 to within 0.6 GPa, but gives
   pad 8: 123.2 (table says 132), pad 9: 164.2 (172), pad 10: **189.5 (198)**.
   198 GPa is impossible: it exceeds the theoretical maximum E[111] = 191.15 GPa. The likely cause
   is rounded or mistyped hkl values. Fig. 5 is plotted against the Table 1 values, so the
   comparison will use the paper's regression line rather than pad 10's point (§3).
5. **Pure [100]/[110]/[111] pads were never measured.** The nearest are:
   - pad 3, (1, 0.1, 0.02), 68 GPa — close to [100];
   - pad 10, (1, 1, 0.9) — close to [111];
   - pad 2, 130 GPa — matches [110] in **modulus only**; its direction (1, 0.76, 0.24) is about 18°
     away from [110].

   This supports comparing at equal E<sub>z</sub>, which is how the paper itself presents Fig. 5.
6. **Fig. 5 is the 30 → 400 °C change; this is inferred, not stated.** The paper gives neither the
   temperature state nor the reference point for Fig. 5. A closed-form bound settles the state
   (Cu α = 16.5 ppm/°C, Si α = 2.6 ppm/°C, h = 850 nm, paper constants):

   | Orientation | Free expansion 30→400 °C | Laterally locked to Si 30→400 °C | Locked 275→400 °C | Fig. 5 (fit) |
   |---|---|---|---|---|
   | [100] | 5.19 nm | 11.49 nm | 3.88 nm | 9.01 nm |
   | [110] | 5.19 nm | 8.98 nm | 3.03 nm | 8.17 nm |
   | [111] | 5.19 nm | 8.38 nm | 2.83 nm | 7.36 nm |

   Fig. 5 falls between the free and locked bounds only for **ΔT = 370 °C (30 → 400 °C)**. From the
   stress-free temperature (275 → 400 °C), even a fully locked pad reaches only 3.9 nm. The values
   must also be relative to the substrate: the 750 µm Si by itself grows about 720 nm in z over
   30 → 400 °C.

---

## 1. Everything the paper states

| # | Item | Value in paper | Source |
|---|---|---|---|
| **Geometry** ||||
| G1 | Pad lateral size | "diameter of 300 nm" (text); 0.3 µm × 0.3 µm square (drawing) | §2 p.2; Fig. 1(a) p.3 |
| G2 | Pad thickness | "~850 nm"; 0.85 µm | §2 p.2; Fig. 1(a) |
| G3 | Pad spacing | "inter-pad spacing of 2.6 µm"; drawn **edge-to-edge** | §2 p.2; Fig. 1(a) |
| G4 | Thickness source | Cu and SiO2 thicknesses "extracted from TEM images" (no numbers) | §2 p.3 |
| G5 | Cap layer | 60 nm SiN deposited after CMP, against Cu oxidation | §2 p.2 |
| G6 | Other layers drawn | TaTaN liner around Cu; SiN between SiO2 and Si | Fig. 1(a) |
| G7 | FEM domain | single 300 nm pad embedded in SiO2, 750 µm Si substrate | §2 p.3 |
| G8 | FEM symmetry | "Symmetric boundary conditions are applied at two sides" | §2 p.3 |
| G9 | FEM code | COMSOL Multiphysics, 3D | §2 p.3 |
| **Materials** ||||
| M1 | Cu elastic constants | C11 = 168.4, C12 = 121.4, C44 = 75.4 GPa (Nye [21]) | §2 p.3 |
| M2 | Cu orientation | rotation matrix from Laue-derived Euler angles per pad (angles not listed) | §2 p.3 |
| M3 | E<sub>z</sub> formula | eq. (1), direction cosines u, v, w | §3.1 p.3 |
| M4 | Pad orientations | hkl along z and E<sub>z</sub> for 10 pads; E<sub>z</sub> 68–198 GPa | Table 1 p.3 |
| M5 | Si, SiO2 | linear isotropic; E, ν, α "taken from literature at ambient temperature" [18–20] (values not given) | §2 p.3 |
| M6 | Cu CTE | not given; "thermal expansion is orientation-independent in cubic materials" | §3.2 p.5 |
| M7 | Plasticity law | Ludwick σ = σ<sub>y</sub> + K ε<sub>p</sub><sup>N</sup>, eq. (3) | §3.2 p.5 |
| M8 | Plasticity search | DoCE (CCD + LHS), 70 runs; σ<sub>y</sub> 100–400 MPa, K 100–3000 MPa, N 0.1–0.9; fitted to 30 °C and 100 °C data | §3.2 p.5 |
| M9 | Fitted values | pad 1 (127 GPa): σ<sub>y</sub> 180–190, K 650–700, N 0.35–0.40 · pad 2 (130): 185–189, 670–720, 0.30–0.37 · pad 3 (68): 140–146, 650–700, 0.65–0.72 · pad 6 (102): 160–170, 700–725, 0.55–0.60 | §3.2 p.5 |
| M10 | Plasticity scope | plastic only for E < 130 GPa, at T < 150 °C and T > 350 °C; none for E > 130 GPa | §3.1 p.4; §3.2 p.5; §4 p.6 |
| M11 | Creep | seen at 400 °C for E < 130 GPa; **not modeled**; Fig. 5 "does not account for creep" | §3.2 p.5 |
| **Reference state and thermal history** ||||
| T1 | Stress-free T (experiment) | 250–300 °C for all pads | §3.1 p.4 |
| T2 | Stress-free T (FEM) | "taken from the experimental curves"; pad 3 stress-free at T ≈ 275 °C | §3.2 p.4–5; Fig. 4 p.5 |
| T3 | Heating steps | RT (30 °C), 100 °C, then 50 °C steps to 400 °C | §2 p.3 |
| T4 | Ramp / dwell | 15 °C/min; ~5 min stabilization per step | §2 p.3 |
| T5 | Peak hold | 400 °C, 3 acquisitions standing in for a 2 h bonding anneal | §2 p.3 |
| T6 | Cooling | to ambient in 100 °C steps, rate not controlled | §2 p.3 |
| **Results** ||||
| R1 | Measured quantity | deviatoric ε′<sub>zz</sub> only; hydrostatic part not accessible | §2 p.3 |
| R2 | ε′<sub>zz</sub> data | vs. E at 30 °C (Fig. 2a); vs. T for pads 3 and 10 (Fig. 2b); Δ30→400 °C (Fig. 2c); exp. vs. FEM (Fig. 3) | pp.4–5 |
| R3 | FEM displacement | Fig. 5: "Total displacement in Z (nm)" vs. E<sub>z</sub>; fit y = 9.887 − 0.0132 E<sub>z</sub>, R² = 0.955 | Fig. 5 p.6 |
| R4 | Orientation spread | ~2 nm between single pads ("4 nm in bonding"); dishing spec should follow [111] | §3.3 p.6 |
| R5 | Stress state | total stress tensile at 30 °C, zero at ~275 °C, compressive above (pad 3, volume average) | §3.2 p.5; Fig. 4 |

**Not stated anywhere:** element type, mesh, outer boundary conditions beyond G8, Cu–dielectric
interface treatment, hardening rule, yield criterion, temperature dependence of any property,
in-plane Euler angle, and which point or average defines "total displacement".

### Fig. 5, digitized

Markers were detected automatically; axes were calibrated on six tick marks. 1 px = 0.0035 nm, so
digitization error is about ±0.01 nm. The digitized points refit to y = 9.886 − 0.0133 E
(R² = 0.961), against the printed 9.887 − 0.0132 E (R² = 0.955). Pad 8 overlaps pad 2 and could not
be separated.

| Pad | E<sub>z</sub> (Table 1, GPa) | u<sub>z</sub> (nm) | Residual vs. paper fit (nm) |
|---|---|---|---|
| 3 | 68 | 8.897 | −0.092 |
| 6 | 102 | 8.395 | −0.146 |
| 1 | 127 | 8.394 | +0.183 |
| 2 | 130 | 8.346 | +0.175 |
| 7 | 164 | 7.693 | −0.029 |
| 4 | 169 | 7.646 | −0.010 |
| 9 | 172* | 7.590 | −0.027 |
| 5 | 188 | 7.340 | −0.065 |
| 10 | 198* | 7.191 | −0.082 |

\* Table 1 E<sub>z</sub> disagrees with eq. (1) (Finding 4). Scatter about the regression:
RMS 0.109 nm, max 0.183 nm.

---

## 2. Values the paper does not give (all marked ASSUMPTION)

| # | Item | Proposed value | Reason | Sensitivity / variant |
|---|---|---|---|---|
| A1 | Pad shape | **square 300 × 300 nm**, sharp corners | Fig. 1(a) is the only dimensioned drawing; "diameter" is read as width | optional: Ø300 nm cylinder |
| A2 | Pad pitch / cell | **2.9 × 2.9 µm** cell, pad centered (2.6 µm edge-to-edge + 0.3 µm) | arrow in Fig. 1(a) runs edge-to-edge; same pitch assumed in x and y | low: neighbor 2.6 µm away |
| A3 | Dielectric stack | SiO2 top **flush with Cu top** (no dishing); **0.30 µm SiO2 below pad** (total 1.15 µm) | post-CMP planar surface; Fig. 1(a) shows a SiO2 layer under the pads, about ⅓ pad height (qualitative) | moderate; check ±0.2 µm |
| A4 | Bottom SiN (SiO2/Si) | **omitted** | not in the paper's FEM description (G7); thin | low |
| A5 | 60 nm SiN cap | **omitted** in baseline | FEM description lists only Cu/SiO2/Si; no cap at bonding | optional run with cap (E 220 GPa, ν 0.27, α 2.3 ppm/°C) — a top cap would directly restrain protrusion **[Phase B: cap **included** by default (G4a); removed in one variant, R11]** |
| A6 | TaTaN liner | **not meshed**; its mechanical role is bracketed by the two interface variants | thickness not given; ~5–20 nm would need very small elements | — |
| A7 | Si substrate | **750 µm, full thickness**, graded mesh in z; isotropic **E = 130 GPa, ν = 0.28, α = 2.6 ppm/°C** | thickness from paper; α is the RT value of ref. [18] (Okada); E is Si's in-plane ⟨100⟩ modulus | E: very low · α: moderate **[Phase B: replaced by a 2 µm Si slab whose cell walls are displaced by α<sub>Si</sub>·(T−T<sub>sf</sub>)·L (infinite-substrate equivalent, B4); R14 checks 4 µm]** |
| A8 | SiO2 | **E = 72 GPa, ν = 0.17, α = 0.5 ppm/°C** | fused-silica RT values, consistent with ref. [19] (Bansal & Doremus) | low–moderate |
| A9 | Cu CTE | **16.5 ppm/°C, constant** | paper used ambient-T constants for the other materials | **largest magnitude lever**: the 30–400 °C mean is ≈ 17.6 ppm/°C, which raises u<sub>z</sub> by ~6–7 %. Optional run **[Phase B: constant 16.5 ppm/K default; NIST α(T) is variant R09, table to be supplied (B7)]** |
| A10 | Cu elasticity vs. T | temperature-independent | paper gives one set (M1) | — |
| A11 | In-plane crystal angle | same frames that passed G3 ([110]: x = [001], y = [1-10]; [111]: x = [1-10], y = [11-2]) | E<sub>z</sub> is independent of it; Euler angles not published | small for [110] (twofold in-plane); optional 45° check |
| A12 | Cu plasticity, [100] | Ludwick, **pad 3 midpoints**: σ<sub>y</sub> = 143 MPa, K = 675 MPa, N = 0.685 | pad 3 (68 GPa) is the measured pad nearest [100] in both E and direction | bracket: range ends from M9 **[Phase B: superseded by the pad-specific mapping rule, B5]** |
| A13 | Cu plasticity, [110] | Ludwick, **pad 2 midpoints**: σ<sub>y</sub> = 187 MPa, K = 695 MPa, N = 0.335 | pad 2's E (130 GPa) equals E[110] = 130.34 GPa | as above **[Phase B: superseded by B5 (pure [110] is not simulated; pads 3 and 10 are)]** |
| A14 | Cu plasticity, [111] | **elastic only** | paper: no plasticity for E > 130 GPa (M10) | — **[Phase B: kept as rule: E<sub>z</sub> > 130 GPa → elastic (pad 10)]** |
| A15 | Yield / hardening | von Mises yield, **isotropic hardening**, Ludwick tabulated as multilinear (TB,PLAS MISO) on a fine ε<sub>p</sub> grid near 0 (N < 1 gives an infinite initial slope); temperature-independent | paper gives neither; fit used one σ<sub>y</sub> per pad, which implies a J2-type yield | Phase B must confirm TB,ANEL + TB,PLAS MISO is accepted by SOLID185/186 in v261 before any production run. Fallback: elastic-only for [100] and [110], reported as a deviation **[Phase B: accepted; ANEL + MISO syntax reused from cp_check lw_*.dat; temperature dependence per B6]** |
| A16 | Stress-free T | **MP,REFT = TREF = 275 °C** | midpoint of 250–300 °C (T1); equals pad 3's stress-free point (T2) | only matters for the plastic orientations; bracket 250/300 °C **[Phase B: per pad: pad 3 = 275 °C, pad 10 = 277 °C (B3)]** |
| A17 | FEM thermal history | uniform T, quasi-static load steps: **275 → 30 → 100 → 150 → … → 400 (50 °C steps) → 300 → 200 → 100 → 30 °C**; hold at 400 °C omitted (no creep) | mirrors the experiment (T3–T6); first step builds the as-processed state at 30 °C | elastic cases are path-independent **[Phase B: accepted, with output at every step (13 load steps, B2)]** |
| A18 | Outer lateral faces | **full model**; each outer face gets a **coupled normal displacement** (CP: the face stays planar and moves as one) — periodic-array approximation | the cell can expand with the Si. Plain rollers would stop Si's lateral expansion and over-restrain the Cu by ≈ 1.03 · α<sub>Si</sub> · ΔT · h ≈ 0.8 nm, half the orientation spread | optional: free outer faces (another reading of G8) **[Phase B: accepted as flat walls with **prescribed** normal displacement (substrate expansion), B4]** |
| A19 | Bottom | Si bottom UZ = 0, plus 3-2-1 point restraints for rigid-body motion | no other loads; substrate stays flat by symmetry | — |
| A20 | Symmetry cut | **none (full model)**. A quarter model is valid for [100] and [110] in the chosen frames (x = 0 and y = 0 are crystal mirror planes), but **not for [111]**, because (11-2) is not a cubic mirror plane | one mesh topology for every orientation, as you asked | — |
| A21 | Sidewall interface | Variant **B** (bonded): shared nodes. Variant **S** (sliding): frictionless standard contact on the 4 sidewalls (CONTA174/TARGE170, augmented Lagrange), separation allowed. **Pad bottom bonded in both** | paper silent; you asked for both | — **[Phase B: sliding implemented as frictionless, **no separation**: normal displacement tied by CE on duplicated sidewall nodes (B7)]** |
| A22 | Elements / mesh | SOLID185 with enhanced strain or SOLID186. Pad ≈ 20 nm elements at G2 level 0, refined until Δu<sub>center</sub> ≤ 0.1 nm; geometric grading into the Si | to be sized in Phase B against the 128k node/element limit | G2 **[Phase B: SOLID185, KEYOPT(2)=3; four structured levels L0–L3 (B8)]** |
| A23 | Displacement definition | **u<sub>center</sub>** = u<sub>z</sub>(Cu top center) − u<sub>z</sub>(Si/dielectric interface), taken as the 30 → 400 °C change on heating. Also reported: u<sub>edge</sub> (top-face edge midpoint), u<sub>avg</sub> (top-face area average), and protrusion relative to far-field oxide surface | paper silent; subtracting the substrate surface removes the ~720 nm Si thickness expansion, which Fig. 5's magnitudes exclude (Finding 6) | center vs. average expected to differ by a few tenths of a nm; both reported **[Phase B: accepted; u and ε′<sub>zz</sub> definitions fixed in B0/B9]** |

---

## 3. Comparison target (superseded by §B2)

Phase A proposed Fig. 5 protrusion as the primary target. Phase B reverses this: the **measured**
deviatoric strain is primary (G4a), and Fig. 5 is a secondary implementation check (G4b). See §B2.

## 4. Run matrix (superseded by §B10 and `docs/RUNLIST.md`)

## 5. Decisions (resolved 2026-10-08)

| # | Question | Decision |
|---|---|---|
| 1 | Square pad, 2.9 µm cell | accepted; circular pad is one variant (R12) |
| 2 | SiN cap / liner | 60 nm cap **included** (default for G4a); no TaTaN liner |
| 3 | Cu CTE | constant 16.5 ppm/K default; NIST α(T) as one variant |
| 4 | Plasticity | TB,ANEL + TB,PLAS MISO from Ludwick (not TB,XTAL), pad-specific mapping (B5) |
| 5 | Outer walls | coupled flat walls (periodic-like) |
| 6 | Comparison basis | measured ε′<sub>zz</sub>(T) primary; Fig. 5 secondary, ±0.3 nm |

---

# B. Phase B specification (decided 2026-10-08)

## B0. Terminology (used identically in code, files and reports)

- **recess** – initial Cu depth after CMP. **Not used in Stage 1.**
- **protrusion u** – vertical rise of the Cu top surface during heating, from FEM (nm).
- **deviatoric strain ε′<sub>zz</sub>** – the measured Laue quantity. Model equivalent:
  ε′<sub>zz</sub> = ε<sub>zz</sub> − (ε<sub>xx</sub> + ε<sub>yy</sub> + ε<sub>zz</sub>)/3, computed from the **elastic**
  strain (EPEL), volume-averaged over all Cu elements, in the global (sample) frame.
  Thermal strain is isotropic in cubic Cu and so has no deviatoric part. Plastic strain does not
  strain the lattice. Taking EPEL is therefore consistent with what Laue measures.

## B1. Execution

ANSYS is never run inside the app. `src/stage1/stage1_pipeline.py` writes the decks; you run them in
your terminal with the repo's `.venv\Scripts\python.exe`, using
MAPDL Student v261 (`ANSYS261.exe -b -smp -np 4`). Results land in `results/`, which I then
read. Default mode is batch (one MAPDL job per run, exit code and `.out` log kept); `--mode auto`
tries PyMAPDL first. The ANEL + MISO syntax is taken unchanged from `src/validation/cp_check.py` (its `lw_*` decks): ANEL
rotated to the global frame with `ESYS,0`, then `TB,PLAS,…,MISO` with `TBPT,DEFI`.

## B2. Validation hierarchy and pass criteria

| Gate | Quantity | Reference | Pass |
|---|---|---|---|
| **G4a (primary)** | ε′<sub>zz</sub>(T) per simulated pad | measured, Figs. 2(b), 3(a), 3(b), digitized (§B9) | model within the measured error bar at every measured T in 30–350 °C on **heating only, for every pad** (decision 2026-10-08, supersedes addendum B). The measured heating/cooling offset is irreversible behaviour the model does not contain. **Cooling is report-only for all pads**: the model-vs-measured cooling offset is reported as a number (×10⁻³) at each cooling point (300, 200, 30 °C). The 400 °C hold is reported, never pass/fail |
| G4b (secondary) | u = Δ[UZ(Cu top centre) − UZ(top surface at the cell edge)], 30 °C → end of heating at 400 °C (before the hold) | digitized Fig. 5 per pad, plus regression u = 9.887 − 0.0132·E<sub>z</sub> | ±0.3 nm. Pads compared against their **own** Fig. 5 point; the regression is evaluated at the **computed** E<sub>z</sub>. Also reported: u at the Cu top edge (same reference) and u referenced to the Cu pad bottom. If the offset to Fig. 5 is uniform across pads while the slope du/dE<sub>z</sub> matches, it is classified as a **reference-definition difference** |
| G2 | u<sub>center</sub>(400 °C) vs mesh | L0 → L1 → L2 → L3 | change ≤ 0.1 nm between the last two levels |
| G1, G3 | — | — | already passed, not repeated |

**Computed vs. listed E<sub>z</sub>** (eq. 1, paper constants): pad 3 68.04 (table 68), pad 10
**189.53** (table 198), pad 1 127.30 (127), pad 2 130.34 (130), pad 6 102.29 (102) GPa. Pads 8/9/10
are misprinted in Table 1. The regression in G4b uses 189.53 for pad 10 (→ 7.385 nm); the digitized
pad-10 point (7.191 nm) is compared directly.

## B3. Pads, orientation and stress-free temperature

- **Loading direction.** The Table 1 continuous Miller indices, sorted |h| ≥ |k| ≥ |l|:
  pad 3 [1, 0.1, 0.02], pad 10 [1, 1, 0.9].
- **In-plane angle.** The paper gives no Euler angles. 0° is defined as global x = the projection of
  crystal [001] (the cube axis most nearly perpendicular to the loading direction) onto the pad
  plane. Elastic runs are done at 0° and 45° (R01/R02, R03/R04), plus a plastic 45° run (R13).
- **Stress-free temperature T<sub>sf</sub>, with MP,REFT = TREF = T<sub>sf</sub> for every material** (asserted
  by the generator):

| Pad | T<sub>sf</sub> | Source |
|---|---|---|
| 3 | **275 °C** | stated in the paper (§3.2, Fig. 4: stress-free at T ≈ 275 °C); the paper's FEM ε′<sub>zz</sub> crosses zero at 275.7 °C. The measured curve crosses at 269.5 °C, but that curve is plastically bent |
| 10 | **277 °C** | zero crossing of a linear fit to the measured heating curve, Fig. 3(a), RMS residual 0.023 × 10⁻³. The paper's FEM curve crosses at 274 °C. **Not stated in the paper; derived** |
| 1, 2, 6 (optional) | 275 °C | **ASSUMPTION**: no per-pad curve exists; midpoint of 250–300 °C |

The whole stack is assumed stress-free at T<sub>sf</sub>. No intrinsic film stress is included.

**Load step 0 (addendum A).** Every run starts stress-free at the pad's T<sub>sf</sub>
(MP,REFT = TREF = T<sub>sf</sub>) and cools to 30 °C (`LS0_cool_Tsf_to_30`, 12 substeps when plastic).
Only then does the measured history follow. With plasticity, this cool-down already yields. All
reported changes (u, Δε′<sub>zz</sub>) are measured from the end of load step 0.

| Load step | 0 | 1–7 | 8 | 9–12 |
|---|---|---|---|---|
| T (°C) | T<sub>sf</sub> → 30 | 100, 150, 200, 250, 300, 350, 400 | 400 (2 h hold; no creep, so the state is unchanged) | 300, 200, 100, 30 |

## B4. Geometry, substrate equivalent, boundary conditions

- Square Cu pad 0.3 × 0.3 × 0.85 µm in a 2.9 × 2.9 µm cell, full model. 0.30 µm SiO2 under the
  pad, 60 nm SiN cap, no liner.
- **Substrate equivalent (replaces the 750 µm Si).** The model has a 2 µm Si slab, bottom UZ = 0.
  The four outer walls stay flat: x = 0 and y = 0 are fixed normal to themselves, and x = L and
  y = L are displaced normal to themselves by α<sub>Si</sub>·(T − T<sub>sf</sub>)·L at every load step, which is
  the free expansion of the thick wafer. Justification:
  - The film's in-plane stiffness relative to 750 µm of Si is
    (E·t)<sub>film</sub>/(E·t)<sub>Si</sub> = (72 × 1.15 + 220 × 0.06)/(130 × 750) ≈ 1.0 × 10⁻³. So the wafer
    strain equals its free expansion to 0.1 %.
  - Flat walls exclude bending, as in a periodic array.
  - The wafer's own thickness expansion (~0.72 µm over 30 → 400 °C) is removed by measuring u
    relative to the Si top surface.
  - R14 repeats a plastic run with a 4 µm slab to confirm the slab is thick enough.
- Materials, isotropic, ambient: SiO2 E 72 GPa, ν 0.17, α 0.5 ppm/K; Si E 130 GPa, ν 0.28,
  α 2.6 ppm/K; SiN E 220 GPa, ν 0.27, α 2.3 ppm/K (**ASSUMPTION**: typical PECVD values; the paper
  gives none, and the cap stiffness bears directly on u).
- Cu: TB,ANEL with C11 = 168.4, C12 = 121.4, C44 = 75.4 GPa (paper), identical in every variant.

## B5. Plasticity model (decided)

**Model.** TB,ANEL + TB,PLAS,…,MISO. The curve is tabulated from Ludwick σ = σ<sub>0</sub> + K·ε<sub>p</sub><sup>N</sup>
on 19 points, ε<sub>p</sub> = 0 … 0.2, with fine spacing near 0 because N < 1. Yield is von Mises with
isotropic hardening.

**Why not TB,XTAL.**
- Crystal plasticity took ~30 s per substep at 64 elements, i.e. more than 50 min per orientation.
- Ayoub's fits are not consistent with a single critical resolved shear stress: back-computed τ<sub>c</sub>
  varies with CV ≈ 12 %, and pads with near-equal Schmid factors differ 1.3× in σ<sub>0</sub>.

**Parameters** (midpoints of the paper's fitted ranges):

| Pad | E (table) | σ<sub>0</sub> (MPa) | K (MPa) | N |
|---|---|---|---|---|
| 3 | 68 | 143 | 675 | 0.685 |
| 6 | 102 | 165 | 712.5 | 0.575 |
| 1 | 127 | 185 | 675 | 0.375 |
| 2 | 130 | 187 | 695 | 0.335 |

**Mapping rule (addendum C) — on COMPUTED E<sub>z</sub>** (eq. 1, paper constants): pad 3 68.04, pad 6
102.29, pad 1 127.30, pad 2 130.34 GPa.
- E<sub>z</sub> > 130.34 GPa (pad 2's computed value) → elastic only. Pad 10, at 189.53 GPa, is elastic.
- 68.04 ≤ E<sub>z</sub> ≤ 130.34 GPa → the flow curves are **interpolated pointwise**: σ is interpolated
  at the same ε<sub>p</sub> grid points between the two neighbouring fitted pads. σ<sub>0</sub>, K and N are
  not interpolated separately, because the paper states the parameter fits are non-unique. Verified:
  E = 115 GPa gives w = 0.508 between pads 6 and 1, matching the pointwise combination exactly.
- 66.6 ≤ E<sub>z</sub> < 68.04 GPa → pad 3 curve.
- A fitted pad lands exactly on its own curve (w = 0 or 1).

**Elastic and plastic runs.** Every pad is run elastic-only and plastic. Pad 10 "plastic" is
elastic by the rule. R08 confirms this without a solve: its deck is byte-identical to R01's apart
from the title line (verified, SHA-256 match).

## B6. Temperature dependence of plasticity

**Variants** (ANEL constants unchanged in all three):
1. **T-independent:** one MISO curve.
2. **DEFAULT, shear-modulus scaled:** σ<sub>0</sub>(T) and K(T) scaled by G(T)/G(30 °C); N fixed, so the
   whole flow curve scales. G is the **Hill-averaged isotropic shear modulus** from Chang & Himmel's
   C11, C12, C44(T), which is consistent with isotropic J2 plasticity. The **C44-only ratio** is
   computed and reported alongside. **G(400)/G(30) is stated for both bases once the data are
   supplied; neither value can be stated now** (see Source status). MISO tables are given at 30,
   100, 150 … 400 °C.
3. **Stronger reduction (ASSUMPTION):** twice the reduction rate of variant 2,
   f₃(T) = 1 − 2·[1 − G(T)/G(30 °C)].

**Source — resolved 2026-10-08.**
- Values are from Y. A. Chang and L. Himmel, *Temperature Dependence of the Elastic Constants of Cu,
  Ag, and Au above Room Temperature*, report **UCRL-16697**, Lawrence Radiation Laboratory (1966),
  published as *J. Appl. Phys.* 37, 3567 (1966), doi:10.1063/1.1708903.
- The file you provided is the report version, so the page reference is to the report: **Table I,
  "Adiabatic Elastic Constants … for Copper", report p. 14 (PDF p. 18)**, 300–800 K in 50 K steps.
- Ledbetter & Naimon 1974 was not used: per your check it gives C<sub>ij</sub>(T) only as figures.
- **Transcription.** Read from the page image, not the OCR. Checked against the table's own
  identities: C′<sub>L</sub> = ½(C11 + C12 + 2C44) to ≤ 0.004, B<sub>s</sub> = (C11 + 2C12)/3 to ≤ 0.003, and
  C′ = ½(C11 − C12) to ≤ 0.002 (× 10¹² dyn/cm²), all within the table's 0.005 rounding of C11 and
  C12. The directly measured C44, C′ and C′<sub>L</sub> are linear in T to ≤ 0.0005, 0.0002 and 0.0022, as
  the authors state.
- **Input file:** `inputs/cu_elastic_constants_vs_T.csv`, with C11, C12, C44 in GPa at
  26.85–526.85 °C and the source on every row.
- **Adiabatic vs isothermal.** Not material here: for cubic crystals C44 and C′ are identical in the
  adiabatic and isothermal states, and the Hill G depends only on C44 and C11 − C12.

| T (°C) | G<sub>Hill</sub>/G<sub>Hill</sub>(30 °C) | C44/C44(30 °C) | variant ③ factor |
|---|---|---|---|
| 100 | 0.9696 | 0.9754 | 0.9392 |
| 200 | 0.9314 | 0.9397 | 0.8628 |
| 300 | 0.8932 | 0.9040 | 0.7864 |
| 350 | 0.8737 | 0.8855 | 0.7474 |
| **400** | **0.8545** | **0.8677** | **0.709** |

- G<sub>Hill</sub>(30 °C) = 47.62 GPa and G<sub>Hill</sub>(400 °C) = 40.69 GPa from Chang & Himmel. Only the
  ratio is used; ANEL stays at the paper's constants.
- If C11 − C12 is taken as 2C′ (directly measured) instead of the printed C11 and C12, the Hill
  ratio at 400 °C becomes 0.8554 (Δ 0.001).
- Pad 3 σ<sub>0</sub> at 400 °C: variant ① 143 MPa, ② 122.2 MPa, ③ 101.4 MPa.

**350–400 °C check (as specified: a check, not a fit; one-sided form accepted 2026-10-08).** Pad 3 has **no 400 °C measurement before
the hold**. The three 400 °C points are the 40 min, 80 min and 2 h acquisitions (1.20, 1.55,
1.85 × 10⁻³). The heating segment therefore has only its 350 °C end (0.499 × 10⁻³). The model's
400 °C end (no creep) can only be checked as an **upper bound**: it must not exceed the 40 min
value (1.20 × 10⁻³). I will report each variant's 350 °C value, its 400 °C value against that
bound, and Δε′<sub>zz</sub> and Δu between variants. **If ①②③ all satisfy the bound and agree at
350 °C within the error bar, I will state plainly that the data cannot rank them.** The sensitivity
(Δε′<sub>zz</sub>, Δu) is then the deliverable; it feeds the Stage 2 Sobol parameter.

## B7. Other variants (one at a time; pad 3, plastic variant 2 unless stated)

| Run | Variant | Status |
|---|---|---|
| R09 | Cu CTE: NIST instantaneous α(T) vs constant 16.5 ppm/K | BLOCKED: the open sources I could reach give only the 283 K value (16.43 ppm/K), not the 30–400 °C table. Needs `inputs/cu_cte_vs_T.csv` (NIST SRM 736, Hahn, *J. Appl. Phys.* 41, 5096 (1970), doi:10.1063/1.1658614) |
| R10 | Sidewall frictionless sliding | implemented as **no separation**: Cu-side sidewall nodes duplicated, normal displacement tied by CE, tangential free; pad top and bottom perimeters (and square corners) stay bonded. A true contact that can open is **not** modelled (open issue O6) |
| R11 | SiN cap removed | |
| R12 | Circular pad, D = 300 nm | O-grid blend; 4 corner elements reach 175° interior angle (below the 179.9° error limit, `SHPP,WARN`) |
| R13 | In-plane 45° (plastic) | |
| R14 | Si slab 4 µm (substrate-equivalent check; variant 1) | READY |

## B8. Mesh (structured hex, generated in Python, SOLID185 KEYOPT(2)=3)

| Level | Pad elements across | Nodes | Elements | Role |
|---|---|---|---|---|
| L0 | 6 (50 nm) | 18,750 | 16,704 | G2 (R15) |
| L1 | 8 (37.5 nm) | 31,958 | 29,008 | G2 (R16) |
| **L2** | 10 (30 nm) | 51,450 | 47,396 | **production** (= G2 level for R05) |
| L3 | 12 (25 nm) | 82,369 | 76,800 | G2 (R17) |

All levels are below the 128,000 limit. The largest variant model is L3 at 82,369 nodes; the
sliding variant at L2 adds 540 duplicated nodes. All element Jacobians are positive (checked).
Production uses L2 before G2 has run. If G2 shows L2 has not converged, the production runs move to
L3 and are repeated.

## B9. Outputs and the digitized measurement

**Per load step** (13 steps: init → 30, heating 100 … 400, end of hold, cooling 300/200/100/30):
- Cu volume-averaged EPEL (6 components) and ε′<sub>zz</sub>; volume-averaged stress (6 components,
  for comparison with Fig. 4).
- Mean and maximum equivalent plastic strain; maximum Cu von Mises stress.
- Sidewall maximum tensile and compressive normal stress, and maximum shear traction, from Cu-side
  nodal stresses.
- **Protrusion** (all are changes since the end of load step 0 at 30 °C; nm):
  - **u<sub>center</sub>** = Δ[UZ(Cu top centre) − UZ(top cap/dielectric surface at the cell edge,
    mid-side node)]. This is the **G4b primary** quantity.
  - u<sub>edge</sub>: the same, with the Cu top edge (mid-side) in place of the centre.
  - u referenced to the Cu pad bottom: Δ[UZ(top centre) − UZ(bottom centre)].
  - Secondary: Cu top centre and top-face average relative to the Si top; surface centre and
    surface edge relative to the Si top.
- **Frame of ε′<sub>zz</sub> (addendum F).**
  - All elements are created with `ESYS,0`; Cu anisotropy is entered as the stiffness rotated to the
    global frame. The element coordinate system is therefore the global (sample) frame.
  - POST1 sets `RSYS,0`. The volume averages use `ETABLE,…,EPEL,X/Y/Z/XY/YZ/XZ` over the Cu elements.
  - Each run writes `*GET,…,ELEM,n,ATTR,ESYS` for the first and last Cu element to `meta.txt`. The
    collector flags `global_frame = True` only if both are 0 and element 1 is Cu.
  - No crystal-frame quantity enters ε′<sub>zz</sub>.
- Per-step CPU time and convergence flag.

**Per run:** node/element counts, wall time, exit code.

**Digitized measured ε′<sub>zz</sub> (× 10⁻³)** — `inputs/ayoub2022_eps_dev_digitized.json`:
- Error bars are read from the bar caps: ±0.08–0.11 × 10⁻³, consistent with the stated 10⁻⁴
  resolution.
- Digitization uncertainty is 0.02–0.03 × 10⁻³. Points hidden behind markers are reconstructed
  from the caps and flagged in the file.

| T (°C) | pad 3 heat | pad 3 cool | pad 10 heat | pad 10 cool |
|---|---|---|---|---|
| 30 | −1.851 ± 0.080 | −0.760 ± 0.078 | −1.197 ± 0.102 | −1.14 ± 0.10 |
| 100 | −1.109 ± 0.103 | — | −0.795 ± 0.100 | — |
| 150 | −0.740 ± 0.104 | | −0.578 ± 0.100 | |
| 200 | −0.39 ± 0.10 | +0.100 ± 0.109 | −0.367 ± 0.102 | −0.118 ± 0.100 |
| 250 | −0.10 ± 0.10 | | −0.144 ± 0.103 | |
| 300 | +0.156 ± 0.10 | +0.995 ± 0.095 | +0.102 ± 0.106 | +0.380 ± 0.104 |
| 350 | +0.499 ± 0.097 | | +0.321 ± 0.101 | |
| 400 | 40/80/120 min: 1.20 / 1.55 / 1.85 | | +0.600 ± 0.100 | |

Cooling was measured at 400 → 300 → 200 → 30 °C only; there is **no 100 °C cooling point**. The
paper's own FEM curves are stored alongside: pad 10 elastic, and pad 3 elastic and plastic.

## B10. Run order

See `docs/RUNLIST.md` (one command per step, expected files, runtime). Steps follow your order:
1. pad 10 elastic
2. pad 3 elastic
3. pad 3 plastic ②
4. ① and ③
5. pad 10 plastic confirmation
6. variants
7. G2
8. optional pads 1/2/6

## B11. Open issues (flagged, not guessed)

- **O1. Cooling scope — resolved (2026-10-08): cooling is report-only for all pads.** Historical note on why: Pad 3's cooling branch
  is reported without pass/fail. Pad 10's cooling **is** pass/fail. The measured pad-10 cooling
  branch differs from its heating branch: +0.380 vs +0.102 × 10⁻³ at 300 °C, and −0.118 vs −0.367 at
  200 °C. An elastic model retraces its heating line. **Pad 10's G4a on cooling is therefore
  expected to fail at 300 °C and 200 °C by ~0.25–0.28 × 10⁻³, about 2.5 error bars, independent of
  model parameters.** This is why cooling was made report-only; the offset is still reported as a
  number.
- **O2. Unknown FEM averaging.** How Ayoub computed FEM ε′<sub>zz</sub> is not stated: volume average,
  beam-footprint weighting (the 500 nm beam is larger than the 300 nm pad), or a point value. Fig. 4
  says "average volume" for stresses only. We use the Cu volume average.
- **O3. No pre-hold 400 °C measurement for pad 3** (§B6). The 350–400 °C check is an inequality
  only.
- **O4. Fig. 3 caption swaps (a) and (b)** relative to the legends; the legends are used.
- **O5. Processing history before the first 30 °C measurement is unknown.** The model cools
  directly T<sub>sf</sub> → 30 °C. For pad 3 this first cooling already yields plastically, which sets
  the 30 °C starting point.
- **O6. Sliding variant cannot open.** At 30 °C the Cu is in lateral tension (Fig. 4), so a real
  frictionless contact could separate; this variant cannot.
- **O7. Assumed cap properties.** The SiN cap properties are assumed and act directly on u.
- **O8. Untested element/material combinations.** Temperature-dependent MISO (multiple TBTEMP) with
  ANEL, and KEYOPT(2)=3 with ANEL + MISO, have not been exercised in v261. cp_check validated ANEL +
  single-temperature MISO with KEYOPT(2)=0. R06 is the first run that tests KEYOPT(2)=3.
- **O10. Hardening rule (addendum G, reported, not changed).** MISO is **isotropic hardening**; the
  paper does not state its hardening rule. With isotropic hardening the yield surface only expands,
  so after forward yielding (load step 0 cooling, and heating > 350 °C) reverse yielding needs a
  stress change of 2·σ<sub>y</sub>(current), which grows with accumulated ε<sub>p</sub>. Kinematic hardening
  (Bauschinger effect) would reverse-yield earlier: during re-heating after the initial cool-down,
  and on cooling after the hold. Isotropic hardening therefore tends to **understate** plastic
  strain in the later half-cycles, i.e. it under-predicts the curvature of ε′<sub>zz</sub>(T) below 150 °C
  on heating and the heating/cooling hysteresis. The effect is largest for pad 3 (lowest σ<sub>0</sub>,
  largest N).
- **O9. Missing property table.** G(T) is resolved (§B6). The NIST α(T) table for R09 is still
  missing (§B7).
