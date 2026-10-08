#!/usr/bin/env python
r"""
Stage 1 - reproduce Ayoub et al., Microelectron. Eng. 261, 111809 (2022) in MAPDL.

3-D single Cu pad (0.3 x 0.3 x 0.85 um) in SiO2 with a 60 nm SiN cap, 2.9 um
periodic cell, anisotropic Cu (TB,ANEL rotated to the global frame, as validated
in G3 / cp_check lw_*.dat) + optional MISO plasticity from a Ludwick curve.
Thermal history of the experiment: T_sf -> 30 C, 30 -> 100 -> 150 ... 400 C,
2 h hold, 400 -> 300 -> 200 -> 100 -> 30 C. Output at every step.

Units uMKS: um, kg, s, degC  ->  uN, MPa.   MP,REFT == TREF == T_sf (asserted).

USAGE (from the repository root, in your own terminal; ANSYS cannot run inside the app)
    .venv\Scripts\python.exe src\stage1\stage1_pipeline.py --list              # run table + status
    .venv\Scripts\python.exe src\stage1\stage1_pipeline.py --case R01          # one run
    .venv\Scripts\python.exe src\stage1\stage1_pipeline.py --case R01,R02      # several runs, in order
    .venv\Scripts\python.exe src\stage1\stage1_pipeline.py --step 1            # one step of the run order
    .venv\Scripts\python.exe src\stage1\stage1_pipeline.py --dry-run --case all

Options: --mode auto|batch|pymapdl (default batch), --np N, --exe PATH, --force, --recollect.
Machine settings (MAPDL path, cores) come from configs/machine.json if present
(copy configs/machine.example.json), else the MAPDL_EXE environment variable, else the
Student v261 default path. All other paths are relative to the repository root.

Outputs: results/<RUN>/  deck (.dat), MAPDL log (.out), hist_*.txt, meta.txt, record.json
         results/results.json, results/results.csv   (all finished runs)
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent              # <repo>/src/stage1
REPO = HERE.parents[1]                              # <repo>
sys.path.insert(0, str(HERE))
import s1_mesh as MESH  # noqa: E402

_DEFAULT_EXE = r"C:\Program Files\ANSYS Inc\ANSYS Student\v261\ansys\bin\winx64\ANSYS261.exe"
_MACHINE_CFG = REPO / "configs" / "machine.json"     # per-PC, not in git
MACHINE = json.loads(_MACHINE_CFG.read_text()) if _MACHINE_CFG.exists() else {}
EXE = Path(os.environ.get("MAPDL_EXE") or MACHINE.get("mapdl_exe") or _DEFAULT_EXE)
DEFAULT_NP = int(MACHINE.get("nproc", 4))
RESULTS = REPO / "results"
INPUTS = REPO / "inputs"
CSV_ELASTIC_T = INPUTS / "cu_elastic_constants_vs_T.csv"   # Chang & Himmel (for variant 2/3)
CSV_CTE_T = INPUTS / "cu_cte_vs_T.csv"                     # NIST (for the CTE variant)

# ---------------------------------------------------------------- materials (MPa, 1/K)
C11, C12, C44 = 168400.0, 121400.0, 75400.0      # paper Sec. 2 (Nye) - ANEL, all variants
CU_ALPHA_CONST = 16.5e-6
OXIDE = dict(EX=72000.0, PRXY=0.17, ALPX=0.5e-6)    # ASSUMPTION A8 (fused silica, RT)
SILICON = dict(EX=130000.0, PRXY=0.28, ALPX=2.6e-6) # ASSUMPTION A7 (RT)
SIN = dict(EX=220000.0, PRXY=0.27, ALPX=2.3e-6)     # ASSUMPTION A5 (PECVD SiN, RT)

# Ludwick midpoints of the paper's fitted ranges (Sec. 3.2) - sigma0, K in MPa
LUDWICK_FITS = {   # pad: (E_table GPa, sigma0, K, N)
    3: (68.0, 143.0, 675.0, 0.685),
    6: (102.0, 165.0, 712.5, 0.575),
    1: (127.0, 185.0, 675.0, 0.375),
    2: (130.0, 187.0, 695.0, 0.335),
}
EPS_P_GRID = [0.0, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 3e-3, 5e-3,
              7.5e-3, 0.01, 0.015, 0.02, 0.03, 0.05, 0.1, 0.2]
MISO_TEMPS = [30.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0]
VARIANT3_FACTOR = 2.0       # variant 3: twice the reduction rate of variant 2 (ASSUMPTION)

# ---------------------------------------------------------------- pads (Table 1)
PADS = {   # Miller indices along z as printed; stress-free T and its source
    1: dict(hkl=(1, 0.6, 0.37), E_table=127, T_sf=275.0, T_sf_src="ASSUMPTION: midpoint of 250-300 C (no per-pad curve)"),
    2: dict(hkl=(1, 0.76, 0.24), E_table=130, T_sf=275.0, T_sf_src="ASSUMPTION: midpoint of 250-300 C (no per-pad curve)"),
    3: dict(hkl=(1, 0.1, 0.02), E_table=68, T_sf=275.0, T_sf_src="paper Sec. 3.2/Fig. 4: stress-free at ~275 C (FEM curve Fig. 3b crosses 0 at 275.7 C)"),
    6: dict(hkl=(1, 0.24, 0.47), E_table=102, T_sf=275.0, T_sf_src="ASSUMPTION: midpoint of 250-300 C (no per-pad curve)"),
    10: dict(hkl=(1, 1, 0.9), E_table=198, T_sf=277.0, T_sf_src="zero crossing of linear fit to measured heating curve, Fig. 3a (rms 0.023e-3); paper FEM curve crosses at 274 C"),
}

# ---------------------------------------------------------------- thermal history
HISTORY = [  # (tag, branch, T_end, dt_s, nsub_plastic)
    ("LS0_cool_Tsf_to_30", "init", 30.0, 1.0, 12),   # load step 0: stress-free T_sf -> 30 C
    ("h100", "heating", 100.0, 70 / 0.25 + 300, 4),
    ("h150", "heating", 150.0, 50 / 0.25 + 300, 4),
    ("h200", "heating", 200.0, 50 / 0.25 + 300, 4),
    ("h250", "heating", 250.0, 50 / 0.25 + 300, 4),
    ("h300", "heating", 300.0, 50 / 0.25 + 300, 4),
    ("h350", "heating", 350.0, 50 / 0.25 + 300, 4),
    ("h400", "heating", 400.0, 50 / 0.25 + 300, 4),
    ("hold400_2h", "hold", 400.0, 7200.0, 1),
    ("c300", "cooling", 300.0, 100 / 0.25, 5),
    ("c200", "cooling", 200.0, 100 / 0.25, 5),
    ("c100", "cooling", 100.0, 100 / 0.25, 5),
    ("c030", "cooling", 30.0, 70 / 0.25, 4),
]


# ================================================================ crystal math
VOIGT = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]   # MAPDL order x y z xy yz xz


def cubic_D(c11=C11, c12=C12, c44=C44):
    D = np.zeros((6, 6)); D[:3, :3] = c12; np.fill_diagonal(D[:3, :3], c11)
    D[3, 3] = D[4, 4] = D[5, 5] = c44
    return D


def _to_tensor(D):
    C = np.zeros((3, 3, 3, 3))
    for I, (i, j) in enumerate(VOIGT):
        for J, (k, l) in enumerate(VOIGT):
            for a, b in ((i, j), (j, i)):
                for c, d in ((k, l), (l, k)):
                    C[a, b, c, d] = D[I, J]
    return C


def _to_voigt(C):
    return np.array([[C[i, j, k, l] for (k, l) in VOIGT] for (i, j) in VOIGT])


def sorted_hkl(hkl):
    return tuple(sorted((abs(v) for v in hkl), reverse=True))


def rotation(hkl, angle_deg=0.0):
    """R with rows = global x, y, z axes in crystal coordinates; R @ v_crys = v_glob.

    z (row 3) = loading direction n = [hkl] (sorted |h|>=|k|>=|l|).
    0 deg in-plane: global x = projection of crystal [001] (the cube axis most
    nearly perpendicular to n) onto the plane normal to n; +angle rotates x about n."""
    n = np.array(sorted_hkl(hkl), float); n /= np.linalg.norm(n)
    ref = np.array([0.0, 0.0, 1.0])
    x = ref - (ref @ n) * n; x /= np.linalg.norm(x)
    a = math.radians(angle_deg)
    x = math.cos(a) * x + math.sin(a) * np.cross(n, x)
    y = np.cross(n, x)
    R = np.vstack([x, y, n])
    assert np.allclose(R @ R.T, np.eye(3)) and np.isclose(np.linalg.det(R), 1.0)
    return R


def rotated_D(hkl, angle_deg=0.0):
    R = rotation(hkl, angle_deg)
    Cg = np.einsum("ip,jq,kr,ls,pqrs->ijkl", R, R, R, R, _to_tensor(cubic_D()), optimize=True)
    D = _to_voigt(Cg)
    D[np.abs(D) < 1e-10 * np.abs(D).max()] = 0.0   # round-off only
    return 0.5 * (D + D.T)


def Ez_closed_form(hkl, c11=C11, c12=C12, c44=C44):
    den = (c11 - c12) * (c11 + 2 * c12)
    s11, s12, s44 = (c11 + c12) / den, -c12 / den, 1 / c44
    l, m, n = np.array(hkl, float) / np.linalg.norm(hkl)
    return 1 / (s11 - 2 * (s11 - s12 - s44 / 2) * (l*l*m*m + m*m*n*n + n*n*l*l))


def Ez_check(hkl, angle_deg=0.0):
    E_fe = 1.0 / np.linalg.inv(rotated_D(hkl, angle_deg))[2, 2]
    E_cf = Ez_closed_form(hkl)
    assert abs(E_fe - E_cf) / E_cf < 1e-8, (hkl, E_fe, E_cf)
    return E_cf


# ================================================================ plasticity
def _fitted_E(pad):
    return Ez_closed_form(sorted_hkl(PADS[pad]["hkl"])) / 1000.0


def _ludwick_curve(s0, K, N):
    return np.array([s0 + K * (e ** N if e > 0 else 0.0) for e in EPS_P_GRID])


def flow_curve_for(pad, Ez_GPa):
    """Mapping rule (spec B5, addendum C) on COMPUTED E_z.

    E_z > E(pad 2) = 130.34 -> elastic (None); E(pad 3) <= E_z <= E(pad 2) -> flow curves
    interpolated POINTWISE (sigma at equal eps_p) between the neighbouring fitted pads;
    66.6 <= E_z < E(pad 3) -> pad 3 curve. Fitted pads get exactly their own curve."""
    fits = sorted(((_fitted_E(p), p, v) for p, v in LUDWICK_FITS.items()))
    E_lo, E_hi = fits[0][0], fits[-1][0]
    if Ez_GPa > E_hi + 1e-9:
        return None
    if Ez_GPa < E_lo - 1e-9:
        if Ez_GPa < 66.6:
            raise ValueError(f"E_z {Ez_GPa:.2f} GPa below the [100] limit")
        _, p, (_, s0, K, N) = fits[0]
        return dict(curve=_ludwick_curve(s0, K, N), rule=f"66.6<=E<{E_lo:.2f}: pad {p} curve", pads=[p])
    for (Ea, pa, va), (Eb, pb, vb) in zip(fits[:-1], fits[1:]):
        if Ea - 1e-9 <= Ez_GPa <= Eb + 1e-9:
            w = 0.0 if Eb == Ea else (Ez_GPa - Ea) / (Eb - Ea)
            if abs(Ez_GPa - Ea) < 1e-6:
                w = 0.0
            if abs(Ez_GPa - Eb) < 1e-6:
                w = 1.0
            ca, cb = _ludwick_curve(*va[1:]), _ludwick_curve(*vb[1:])
            rule = (f"pad {pa} own curve" if w == 0.0 else f"pad {pb} own curve" if w == 1.0 else
                    f"pointwise interpolation, w = {w:.4f} between pad {pa} (E {Ea:.2f}) and pad {pb} (E {Eb:.2f})")
            return dict(curve=(1 - w) * ca + w * cb, rule=rule, pads=[pa, pb], w=w)
    raise AssertionError


E_ELASTIC_ABOVE = None   # filled after PADS is defined (computed E of pad 2)


def _read_csv_rows(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh) if r and not str(next(iter(r.values()), "")).startswith("#")]
    return rows


def shear_ratio_table():
    """G(T)/G(30 C) from the user-supplied Chang & Himmel constants.

    Returns (hill_ratio, c44_ratio, sources): Hill (Voigt-Reuss-Hill) average is used for
    variant 2; the C44-only ratio is reported alongside (addendum D)."""
    if not CSV_ELASTIC_T.exists():
        raise FileNotFoundError(
            f"{CSV_ELASTIC_T.name} missing. Variants 2/3 need Chang & Himmel (1966) C11, C12, C44 vs T "
            f"with a source column; see {CSV_ELASTIC_T.name}.template")
    rows = _read_csv_rows(CSV_ELASTIC_T)
    if not rows or not all(r.get("source", "").strip() for r in rows):
        raise ValueError(f"{CSV_ELASTIC_T.name}: every row needs a non-empty 'source'")
    T = np.array([float(r["T_C"]) for r in rows])
    c11 = np.array([float(r["C11_GPa"]) for r in rows])
    c12 = np.array([float(r["C12_GPa"]) for r in rows])
    c44 = np.array([float(r["C44_GPa"]) for r in rows])
    if T.min() > 30.0 or T.max() < 400.0:
        raise ValueError("G(T) table must span 30-400 C (Chang & Himmel start at 300 K = 26.85 C)")
    Gv = (c11 - c12 + 3 * c44) / 5
    Gr = 5 * (c11 - c12) * c44 / (4 * c44 + 3 * (c11 - c12))
    G = 0.5 * (Gv + Gr)
    hill = {t: float(np.interp(t, T, G) / np.interp(30.0, T, G)) for t in MISO_TEMPS}
    c44r = {t: float(np.interp(t, T, c44) / np.interp(30.0, T, c44)) for t in MISO_TEMPS}
    return hill, c44r, sorted({r["source"] for r in rows})


def miso_table(base_curve, variant):
    """Return ([(T, [(eps_p, sigma), ...]), ...], notes, ratios) for variant 1/2/3.
    sigma0 and K scale by the same factor and N is fixed, so the whole curve scales."""
    if variant == 1:
        return [(30.0, list(zip(EPS_P_GRID, base_curve)))], ["variant 1: temperature independent"], {}
    hill, c44r, src = shear_ratio_table()
    out = []
    for T in MISO_TEMPS:
        f = hill[T] if variant == 2 else 1.0 - VARIANT3_FACTOR * (1.0 - hill[T])
        out.append((T, list(zip(EPS_P_GRID, f * base_curve))))
    ratios = dict(G_hill_400_over_30=hill[400.0], C44_400_over_30=c44r[400.0], hill=hill, c44=c44r)
    return out, [f"G(T)/G(30C), Hill average, from: {s}" for s in src], ratios


def cte_table():
    if not CSV_CTE_T.exists():
        raise FileNotFoundError(
            f"{CSV_CTE_T.name} missing. The CTE variant needs NIST Cu instantaneous CTE vs T "
            f"with a source column; see {CSV_CTE_T.name}.template")
    rows = _read_csv_rows(CSV_CTE_T)
    if not all(r.get("source", "").strip() for r in rows):
        raise ValueError(f"{CSV_CTE_T.name}: every row needs a non-empty 'source'")
    T = [float(r["T_C"]) for r in rows]
    a = [float(r["alpha_inst_per_K"]) for r in rows]
    if min(T) > 30 or max(T) < 400:
        raise ValueError("CTE table must span 30-400 C")
    return list(zip(T, a)), sorted({r["source"] for r in rows})


# ================================================================ run registry
BASE = dict(model="elastic", variant=1, mesh="L2", shape="square", cap=True,
            sliding=False, cte="const", angle=0.0, t_si=MESH.T_SI)


def _run(rid, pad, step, note, **kw):
    d = dict(BASE); d.update(kw); d.update(id=rid, pad=pad, step=step, note=note)
    return d


RUNS = [
    _run("R01", 10, 1, "pad 10 elastic, in-plane 0 deg"),
    _run("R02", 10, 1, "pad 10 elastic, in-plane 45 deg", angle=45.0),
    _run("R03", 3, 2, "pad 3 elastic, in-plane 0 deg"),
    _run("R04", 3, 2, "pad 3 elastic, in-plane 45 deg", angle=45.0),
    _run("R05", 3, 3, "pad 3 plastic, variant 2 (G-scaled) - DEFAULT", model="plastic", variant=2),
    _run("R06", 3, 4, "pad 3 plastic, variant 1 (T-independent)", model="plastic", variant=1),
    _run("R07", 3, 4, "pad 3 plastic, variant 3 (2x reduction)", model="plastic", variant=3),
    _run("R08", 10, 5, "pad 10 plastic: E_z > 130 GPa -> elastic by rule (deck must equal R01)", model="plastic", variant=2),
    _run("R09", 3, 6, "variant: Cu CTE(T) NIST", model="plastic", variant=2, cte="nist"),
    _run("R10", 3, 6, "variant: sidewall frictionless sliding (normal tied)", model="plastic", variant=2, sliding=True),
    _run("R11", 3, 6, "variant: SiN cap removed", model="plastic", variant=2, cap=False),
    _run("R12", 3, 6, "variant: circular pad D = 300 nm", model="plastic", variant=2, shape="circle"),
    _run("R13", 3, 6, "variant: pad 3 plastic in-plane 45 deg", model="plastic", variant=2, angle=45.0),
    _run("R14", 3, 6, "check: Si slab 4 um (substrate equivalent), variant 1", model="plastic", variant=1, t_si=4.0),
    _run("R15", 3, 7, "G2 mesh L0", model="plastic", variant=2, mesh="L0"),
    _run("R16", 3, 7, "G2 mesh L1", model="plastic", variant=2, mesh="L1"),
    #     G2 level L2 is R05
    _run("R17", 3, 7, "G2 mesh L3", model="plastic", variant=2, mesh="L3"),
    _run("R18", 1, 8, "optional pad 1 elastic"),
    _run("R19", 1, 8, "optional pad 1 plastic v2", model="plastic", variant=2),
    _run("R20", 2, 8, "optional pad 2 elastic"),
    _run("R21", 2, 8, "optional pad 2 plastic v2", model="plastic", variant=2),
    _run("R22", 6, 8, "optional pad 6 elastic"),
    _run("R23", 6, 8, "optional pad 6 plastic v2", model="plastic", variant=2),
    # diagnostics for the R01 sign/protrusion finding (do not change the baseline)
    _run("R24", 10, 9, "DIAG pad 10 elastic, SiN cap removed", cap=False),
    _run("R25", 10, 9, "DIAG pad 10 elastic, sidewall sliding (normal tied)", sliding=True),
    _run("R26", 10, 9, "DIAG pad 10 elastic, no cap + sliding", cap=False, sliding=True),
]
RUN_BY_ID = {r["id"]: r for r in RUNS}
E_ELASTIC_ABOVE = _fitted_E(2)       # 130.34 GPa: above -> elastic; also the G4a cooling-scope split


# ================================================================ deck writer
def _f(v):
    return f"{v:.10g}"


def material_block(run, Tsf):
    """Return (lines, info) for all materials."""
    pad = PADS[run["pad"]]
    hkl = sorted_hkl(pad["hkl"])
    Ez = Ez_check(hkl, run["angle"])
    D = rotated_D(hkl, run["angle"])
    L, info = [], dict(hkl_sorted=hkl, Ez_GPa=Ez / 1000.0, sources=[])
    L += ["! ---- MAT 1: Cu, anisotropic elastic (paper constants), rotated to global frame",
          "TB,ANEL,1,1,,0", "TBTEMP,25"]
    vals = [D[i, j] for j in range(6) for i in range(j, 6)]
    for s in range(0, 21, 6):
        L.append(f"TBDATA,{s + 1}," + ",".join(_f(v) for v in vals[s:s + 6]))
    if run["cte"] == "const":
        L += [f"MP,ALPX,1,{_f(CU_ALPHA_CONST)}", f"MP,ALPY,1,{_f(CU_ALPHA_CONST)}", f"MP,ALPZ,1,{_f(CU_ALPHA_CONST)}"]
        info["cte"] = f"constant {CU_ALPHA_CONST:g} 1/K"
    else:
        tab, src = cte_table()
        L.append("MPTEMP")
        for s in range(0, len(tab), 6):
            L.append(f"MPTEMP,{s + 1}," + ",".join(_f(t) for t, _ in tab[s:s + 6]))
        for lab in ("CTEX", "CTEY", "CTEZ"):
            for s in range(0, len(tab), 6):
                L.append(f"MPDATA,{lab},1,{s + 1}," + ",".join(_f(a) for _, a in tab[s:s + 6]))
        L.append("MPTEMP")
        info["cte"] = "instantaneous CTE table (CTEX/Y/Z)"; info["sources"] += src
    L.append(f"MP,REFT,1,{_f(Tsf)}")
    fc = flow_curve_for(run["pad"], Ez / 1000.0) if run["model"] == "plastic" else None
    info["plastic"] = fc is not None
    if run["model"] == "plastic" and fc is None:
        info["plastic_note"] = f"E_z = {Ez / 1000:.2f} GPa > {E_ELASTIC_ABOVE:.2f} GPa: elastic only by mapping rule"
    if fc is not None:
        tab, src, ratios = miso_table(fc["curve"], run["variant"])
        info["flow_rule"] = fc["rule"]; info["sources"] += src; info["shear_ratios"] = ratios
        info["sigma0_MPa"] = float(fc["curve"][0])
        info["miso_scale"] = {T: c[0][1] / tab[0][1][0][1] for T, c in tab}
        L.append(f"TB,PLAS,1,{len(tab)},{len(EPS_P_GRID)},MISO")
        for T, curve in tab:
            L.append(f"TBTEMP,{_f(T)}")
            for e, s in curve:
                L.append(f"TBPT,DEFI,{_f(e)},{_f(s)}")
    for m, props, name in ((2, OXIDE, "SiO2"), (3, SILICON, "Si"), (4, SIN, "SiN")):
        L.append(f"! ---- MAT {m}: {name}, isotropic, ambient properties")
        L += [f"MP,EX,{m},{_f(props['EX'])}", f"MP,PRXY,{m},{_f(props['PRXY'])}",
              f"MP,ALPX,{m},{_f(props['ALPX'])}", f"MP,REFT,{m},{_f(Tsf)}"]
    return L, info


def build_deck(run):
    pad = PADS[run["pad"]]
    Tsf = pad["T_sf"]
    mesh = MESH.build_mesh(run["mesh"], run["shape"], run["cap"], run["sliding"], run["t_si"])
    if mesh.n_nodes > MESH.STUDENT_LIMIT or mesh.n_elems > MESH.STUDENT_LIMIT:
        raise RuntimeError(f"{run['id']}: {mesh.n_nodes} nodes / {mesh.n_elems} elements exceed Student limit")
    mat_lines, info = material_block(run, Tsf)
    plastic = info["plastic"]
    nls = len(HISTORY)
    d = [f"/TITLE,Stage1 {run['id']} pad{run['pad']} {run['model']} v{run['variant']} {run['mesh']} {run['shape']}",
         "/NERR,200,10000000",
         "/PREP7"]
    if run["shape"] == "circle":
        d.append("SHPP,WARN")
    d += ["ET,1,SOLID185", "KEYOPT,1,2,3            ! simplified enhanced strain"]
    d += mat_lines
    d.append("! ---- nodes")
    d += [f"N,{int(r[0])},{_f(r[1])},{_f(r[2])},{_f(r[3])}" for r in mesh.nodes]
    d += ["TYPE,1", "ESYS,0"]
    for m, conn in mesh.elems.items():
        d.append(f"MAT,{m}")
        d += ["E," + ",".join(str(int(v)) for v in e) for e in conn]
    if mesh.ce_pairs:
        d.append("! ---- frictionless sidewall: normal displacement tied (no separation)")
        for q, (ncu, nox, nxv, nyv) in enumerate(mesh.ce_pairs, start=1):
            terms = [(ncu, "UX", nxv), (ncu, "UY", nyv), (nox, "UX", -nxv), (nox, "UY", -nyv)]
            terms = [t for t in terms if abs(t[2]) > 1e-12]
            first = ",".join(f"{n},{lab},{_f(c)}" for n, lab, c in terms[:3])
            d.append(f"CE,{q},0,{first}")
            if len(terms) > 3:
                d.append(f"CE,{q},," + ",".join(f"{n},{lab},{_f(c)}" for n, lab, c in terms[3:]))
    eps = 1e-6
    zbot = -mesh.t_si
    d += ["ALLSEL,ALL",
          f"NSEL,S,LOC,X,{-eps},{eps}", "CM,WXM,NODE",
          f"NSEL,S,LOC,X,{MESH.L_CELL - eps},{MESH.L_CELL + eps}", "CM,WXP,NODE",
          f"NSEL,S,LOC,Y,{-eps},{eps}", "CM,WYM,NODE",
          f"NSEL,S,LOC,Y,{MESH.L_CELL - eps},{MESH.L_CELL + eps}", "CM,WYP,NODE",
          f"NSEL,S,LOC,Z,{zbot - eps},{zbot + eps}", "CM,ZBOT,NODE",
          "ALLSEL,ALL",
          "*GET,NELEM,ELEM,0,COUNT", "*GET,NNODE,NODE,0,COUNT",
          "FINISH", "",
          "/SOLU", "ANTYPE,STATIC", "NLGEOM,OFF", f"TREF,{_f(Tsf)}",
          "OUTRES,ERASE", "OUTRES,ALL,LAST", "AUTOTS,ON", "KBC,0",
          "CMSEL,S,ZBOT", "D,ALL,UZ,0", "CMSEL,S,WXM", "D,ALL,UX,0", "CMSEL,S,WYM", "D,ALL,UY,0", "ALLSEL,ALL",
          f"*DIM,CPUL,ARRAY,{nls}", f"*DIM,CNVL,ARRAY,{nls}", f"*DIM,TLS,ARRAY,{nls}"]
    t = 0.0
    for ls, (tag, branch, T, dt, nsub) in enumerate(HISTORY, start=1):
        t += dt
        ew = SILICON["ALPX"] * (T - Tsf) * MESH.L_CELL     # substrate-equivalent wall displacement
        ns = nsub if plastic else 1
        d += [f"! ---- LS{ls} {tag} ({branch}) -> {T:g} C",
              f"TIME,{_f(t)}", f"TUNIF,{_f(T)}",
              "CMSEL,S,WXP", f"D,ALL,UX,{_f(ew)}", "CMSEL,S,WYP", f"D,ALL,UY,{_f(ew)}", "ALLSEL,ALL",
              f"NSUBST,{ns},{max(ns * 25, 50)},{max(1, ns // 2)}",
              "*GET,TC0_,ACTIVE,0,TIME,CPU", "SOLVE", "*GET,TC1_,ACTIVE,0,TIME,CPU",
              f"CPUL({ls})=TC1_-TC0_", "*GET,CNV_,ACTIVE,0,SOLU,CNVG", f"CNVL({ls})=CNV_", f"TLS({ls})={_f(T)}"]
    d += ["FINISH", ""]
    # ---------------- post-processing
    pr = mesh.probes
    sw = mesh.sidewall
    d += ["/POST1", "RSYS,0",
          f"NSW={len(sw)}", f"*DIM,SWN,ARRAY,{len(sw)},3"]
    for i, (n, nxv, nyv) in enumerate(sw, start=1):
        d.append(f"SWN({i},1)={n}"); d.append(f"SWN({i},2)={_f(nxv)}"); d.append(f"SWN({i},3)={_f(nyv)}")
    d += [f"NCT={len(mesh.cu_top)}", f"*DIM,CUT,ARRAY,{len(mesh.cu_top)}"]
    d += [f"CUT({i})={n}" for i, n in enumerate(mesh.cu_top, start=1)]
    zz, kk = mesh.z, mesh.k
    z_top = (zz[kk["k_pt"] - 1], zz[kk["k_pt"]]); z_bot = (zz[kk["k_pb"]], zz[kk["k_pb"] + 1])
    k_mid = (kk["k_pb"] + kk["k_pt"]) // 2
    z_mid = (zz[k_mid], zz[k_mid + 1])
    d += [f"*DIM,RD,ARRAY,{nls},4",
          f"*DIM,RA,ARRAY,{nls},10", f"*DIM,RB,ARRAY,{nls},11", f"*DIM,RC,ARRAY,{nls},11",
          f"*DO,LS_,1,{nls}",
          "SET,LS_,LAST",
          "ESEL,S,MAT,,1",
          "ETABLE,EVOL,VOLU"]
    for c in ("X", "Y", "Z", "XY", "YZ", "XZ"):
        d += [f"ETABLE,E{c},EPEL,{c}", f"SMULT,WE{c},E{c},EVOL", f"ETABLE,S{c},S,{c}", f"SMULT,WS{c},S{c},EVOL"]
    if plastic:
        d += ["ETABLE,EPQ,EPPL,EQV", "SMULT,WPQ,EPQ,EVOL"]
    d += ["ETABLE,SEQ,S,EQV", "SSUM", "*GET,SVOL_,SSUM,,ITEM,EVOL",
          "RA(LS_,1)=LS_", "RA(LS_,2)=TLS(LS_)"]
    for j, c in enumerate(("X", "Y", "Z", "XY", "YZ", "XZ")):
        d += [f"*GET,TMP_,SSUM,,ITEM,WE{c}", f"RA(LS_,{3 + j})=TMP_/SVOL_",
              f"*GET,TMP_,SSUM,,ITEM,WS{c}", f"RB(LS_,{2 + j})=TMP_/SVOL_"]
    if plastic:
        d += ["*GET,TMP_,SSUM,,ITEM,WPQ", "RA(LS_,9)=TMP_/SVOL_",
              "ESORT,ETAB,EPQ,0", "*GET,TMP_,SORT,,MAX", "RA(LS_,10)=TMP_", "EUSORT"]
    else:
        d += ["RA(LS_,9)=0", "RA(LS_,10)=0"]
    d += ["ESORT,ETAB,SEQ,0", "*GET,TMP_,SORT,,MAX", "RB(LS_,8)=TMP_", "EUSORT",
          "RB(LS_,1)=LS_",
          "RD(LS_,1)=LS_"]
    for col, (za, zb) in ((2, z_top), (3, z_mid), (4, z_bot)):
        d += ["ESEL,S,MAT,,1", f"ESEL,R,CENT,Z,{_f(za)},{_f(zb)}",
              "ETABLE,LVOL,VOLU", "ETABLE,LSZ,S,Z", "SMULT,LWSZ,LSZ,LVOL", "SSUM",
              "*GET,LV_,SSUM,,ITEM,LVOL", "*GET,LS_SZ,SSUM,,ITEM,LWSZ", f"RD(LS_,{col})=LS_SZ/LV_"]
    d += ["ESEL,S,MAT,,1",
          "NSLE,S",
          "SNMX_=-1E30", "SNMN_=1E30", "STMX_=0",
          "*DO,II_,1,NSW",
          "ND_=SWN(II_,1)", "NXV_=SWN(II_,2)", "NYV_=SWN(II_,3)",
          "*GET,QX_,NODE,ND_,S,X", "*GET,QY_,NODE,ND_,S,Y", "*GET,QXY_,NODE,ND_,S,XY",
          "*GET,QYZ_,NODE,ND_,S,YZ", "*GET,QXZ_,NODE,ND_,S,XZ",
          "SNN_=NXV_*NXV_*QX_+NYV_*NYV_*QY_+2*NXV_*NYV_*QXY_",
          "TX_=NXV_*QX_+NYV_*QXY_", "TY_=NXV_*QXY_+NYV_*QY_", "TZ_=NXV_*QXZ_+NYV_*QYZ_",
          "TS2_=TX_*TX_+TY_*TY_+TZ_*TZ_-SNN_*SNN_",
          "*IF,TS2_,LT,0,THEN", "TS2_=0", "*ENDIF",
          "TS_=SQRT(TS2_)",
          "*IF,SNN_,GT,SNMX_,THEN", "SNMX_=SNN_", "*ENDIF",
          "*IF,SNN_,LT,SNMN_,THEN", "SNMN_=SNN_", "*ENDIF",
          "*IF,TS_,GT,STMX_,THEN", "STMX_=TS_", "*ENDIF",
          "*ENDDO",
          "RB(LS_,9)=SNMX_", "RB(LS_,10)=SNMN_", "RB(LS_,11)=STMX_",
          "ALLSEL,ALL",
          "UAV_=0",
          "*DO,II_,1,NCT", "UAV_=UAV_+UZ(CUT(II_))", "*ENDDO",
          "RC(LS_,1)=LS_",
          f"RC(LS_,2)=UZ({pr['cu_top_center']})",
          f"RC(LS_,3)=UZ({pr['cu_top_edge']})",
          "RC(LS_,4)=UAV_/NCT",
          f"RC(LS_,5)=UZ({pr['cu_bottom_center']})",
          f"RC(LS_,6)=UZ({pr['si_top_corner']})",
          f"RC(LS_,7)=UZ({pr['surface_corner']})",
          f"RC(LS_,8)=UZ({pr['surface_center']})",
          "RC(LS_,9)=CPUL(LS_)", "RC(LS_,10)=CNVL(LS_)",
          f"RC(LS_,11)=UZ({pr['surface_edge_mid']})",
          "*ENDDO",
          "*CFOPEN,hist_strain,txt",
          "*VWRITE,RA(1,1),RA(1,2),RA(1,3),RA(1,4),RA(1,5),RA(1,6),RA(1,7),RA(1,8),RA(1,9),RA(1,10)",
          "(10E22.13)", "*CFCLOS",
          "*CFOPEN,hist_stress,txt",
          "*VWRITE,RB(1,1),RB(1,2),RB(1,3),RB(1,4),RB(1,5),RB(1,6),RB(1,7),RB(1,8),RB(1,9),RB(1,10),RB(1,11)",
          "(11E22.13)", "*CFCLOS",
          "*CFOPEN,hist_disp,txt",
          "*VWRITE,RC(1,1),RC(1,2),RC(1,3),RC(1,4),RC(1,5),RC(1,6),RC(1,7),RC(1,8),RC(1,9),RC(1,10),RC(1,11)",
          "(11E22.13)", "*CFCLOS",
          f"*GET,ESY1_,ELEM,1,ATTR,ESYS",
          f"*GET,ESY2_,ELEM,{len(mesh.elems[MESH.MAT_CU])},ATTR,ESYS",
          f"*GET,MAT1_,ELEM,1,ATTR,MAT",
          "*CFOPEN,hist_force,txt",
          "*VWRITE,RD(1,1),RD(1,2),RD(1,3),RD(1,4)",
          "(4E22.13)", "*CFCLOS",
          "*CFOPEN,meta,txt",
          "*VWRITE,NELEM,NNODE,ESY1_,ESY2_,MAT1_",
          "(5E22.13)", "*CFCLOS",
          "FINISH", ""]
    text = "\n".join(d) + "\n"
    # unit/consistency assertions on the generated text
    refts = {float(v) for v in re.findall(r"^MP,REFT,\d+,(\S+)$", text, re.M)}
    trefs = {float(v) for v in re.findall(r"^TREF,(\S+)$", text, re.M)}
    assert refts == trefs == {Tsf}, (refts, trefs, Tsf)
    meta = dict(run=run, T_sf=Tsf, T_sf_src=pad["T_sf_src"], n_nodes=mesh.n_nodes, n_elems=mesh.n_elems,
                probes=mesh.probes, n_sidewall=len(sw), n_ce=len(mesh.ce_pairs), mat=info,
                history=[dict(ls=i + 1, tag=h[0], branch=h[1], T=h[2]) for i, h in enumerate(HISTORY)])
    # identity hash ignores the title line (run id) so equal physics -> equal hash
    body = text.split("\n", 1)[1]
    meta["physics_sha256"] = hashlib.sha256(body.encode()).hexdigest()
    return text, meta


# ================================================================ status / list
def run_status(run):
    try:
        if run["model"] == "plastic" and run["variant"] in (2, 3):
            pad = PADS[run["pad"]]
            Ez = Ez_closed_form(sorted_hkl(pad["hkl"])) / 1000
            if flow_curve_for(run["pad"], Ez) is not None:
                shear_ratio_table()
        if run["cte"] == "nist":
            cte_table()
    except (FileNotFoundError, ValueError) as exc:
        return "BLOCKED", str(exc).split(".")[0]
    done = (RESULTS / run["id"] / "hist_disp.txt").exists()
    return ("DONE" if done else "READY"), ""


def estimate_minutes(run, n_nodes):
    """Rough wall-time estimate (minutes) for -smp -np 4 sparse solver; calibrate with R01."""
    dof = 3 * n_nodes
    t_solve = 9e-6 * dof ** 1.15 / 60            # minutes per factorisation (~8 s at 155k DOF, 4 cores; +-2x)
    pad = PADS[run["pad"]]
    Ez = Ez_closed_form(sorted_hkl(pad["hkl"])) / 1000
    plastic = run["model"] == "plastic" and flow_curve_for(run["pad"], Ez) is not None
    nsolves = sum(h[4] for h in HISTORY) * 2.5 if plastic else len(HISTORY)
    read = n_nodes / 4000 / 60                   # deck parsing
    t = read + nsolves * t_solve + 0.3
    return t


def list_runs(selected=None):
    print(f"{'ID':4s} {'step':>4s} {'pad':>3s} {'model':8s} {'var':>3s} {'mesh':4s} {'shape':6s} "
          f"{'cap':3s} {'iface':6s} {'cte':5s} {'ang':>4s} {'nodes':>7s} {'elems':>7s} {'~min':>6s}  status")
    for r in RUNS:
        if selected and r["id"] not in selected:
            continue
        m = MESH.build_mesh(r["mesh"], r["shape"], r["cap"], r["sliding"], r["t_si"])
        st, why = run_status(r)
        print(f"{r['id']:4s} {r['step']:4d} {r['pad']:3d} {r['model']:8s} {r['variant']:3d} {r['mesh']:4s} "
              f"{r['shape']:6s} {'yes' if r['cap'] else 'no':3s} {'slide' if r['sliding'] else 'bond':6s} "
              f"{r['cte']:5s} {r['angle']:4.0f} {m.n_nodes:7d} {m.n_elems:7d} {estimate_minutes(r, m.n_nodes):6.1f}  "
              f"{st}{(' - ' + why) if why else ''}")


# ================================================================ runners
def _write_case(run, outdir):
    text, meta = build_deck(run)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / f"{run['id']}.dat").write_text(text)
    (outdir / "deck_meta.json").write_text(json.dumps(meta, indent=1, default=float))
    return text, meta


def run_batch(run, outdir, nproc):
    deck = outdir / f"{run['id']}.dat"
    inp = outdir / f"{run['id']}_batch.dat"
    inp.write_text(deck.read_text() + "/EXIT,NOSAVE\n")
    for f in ("hist_strain.txt", "hist_stress.txt", "hist_disp.txt", "meta.txt"):
        (outdir / f).unlink(missing_ok=True)
    cmd = [str(EXE), "-b", "-smp", "-np", str(nproc), "-dir", str(outdir), "-j", run["id"],
           "-i", str(inp), "-o", str(outdir / f"{run['id']}.out")]
    t0 = time.perf_counter()
    p = subprocess.run(cmd, capture_output=True, text=True, errors="replace")
    return dict(method="batch", exit_code=p.returncode, wall_s=time.perf_counter() - t0, cmd=" ".join(cmd))


_MAPDL = None


def run_pymapdl(run, outdir, nproc):
    global _MAPDL
    from ansys.mapdl.core import launch_mapdl
    if _MAPDL is None:
        _MAPDL = launch_mapdl(exec_file=str(EXE), run_location=str(RESULTS / "_pymapdl"), nproc=nproc,
                              additional_switches="-smp", override=True, start_timeout=180)
    m = _MAPDL
    m.clear()
    m.cwd(str(outdir))
    t0 = time.perf_counter()
    err = None
    try:
        m.input(str(outdir / f"{run['id']}.dat"))
    except Exception as exc:  # keep going, record failure
        err = f"{type(exc).__name__}: {exc}"
    return dict(method="pymapdl", exit_code=0 if err is None else 1, wall_s=time.perf_counter() - t0, error=err)


def collect(run, outdir, exec_info):
    meta = json.loads((outdir / "deck_meta.json").read_text())
    rec = dict(id=run["id"], exec=exec_info, deck=meta, ok=False)
    try:
        A = np.loadtxt(outdir / "hist_strain.txt", ndmin=2)
        B = np.loadtxt(outdir / "hist_stress.txt", ndmin=2)
        C = np.loadtxt(outdir / "hist_disp.txt", ndmin=2)
        mt = np.loadtxt(outdir / "meta.txt", ndmin=1)
    except OSError as exc:
        rec["error"] = f"result files missing: {exc}"
        return rec
    rows = []
    def d(i, a, b):          # change since load step 0 end (30 C, heating start), nm
        return 1000.0 * ((C[i, a] - C[i, b]) - (C[0, a] - C[0, b]))
    # columns of hist_disp: 1 cu_top_center, 2 cu_top_edge, 3 cu_top_avg, 4 cu_bottom_center,
    # 5 si_top_corner, 6 surface_corner, 7 surface_center, 8 cpu, 9 cnvg, 10 surface_edge_mid
    for i in range(A.shape[0]):
        ex, ey, ez = A[i, 2], A[i, 3], A[i, 4]
        h = HISTORY[i]
        rows.append(dict(
            ls=int(A[i, 0]) - 1, tag=h[0], branch=h[1], T=A[i, 1],
            epel_x=ex, epel_y=ey, epel_z=ez, epel_xy=A[i, 5], epel_yz=A[i, 6], epel_xz=A[i, 7],
            eps_dev_zz=ez - (ex + ey + ez) / 3.0,
            eppl_eqv_avg=A[i, 8], eppl_eqv_max=A[i, 9],
            s_x=B[i, 1], s_y=B[i, 2], s_z=B[i, 3], s_xy=B[i, 4], s_yz=B[i, 5], s_xz=B[i, 6],
            seqv_max_cu=B[i, 7], sidewall_sn_max=B[i, 8], sidewall_sn_min=B[i, 9], sidewall_tau_max=B[i, 10],
            # G4b primary (addendum E): Cu top centre relative to the top surface at the cell edge
            u_center_nm=d(i, 1, 10),
            u_edge_nm=d(i, 2, 10),
            u_center_rel_padbottom_nm=d(i, 1, 4),
            u_center_rel_si_nm=d(i, 1, 5),
            u_avg_rel_si_nm=d(i, 3, 5),
            u_surface_edge_rel_si_nm=d(i, 10, 5),
            u_surface_center_rel_si_nm=d(i, 7, 5),
            cpu_s=C[i, 8], converged=int(round(C[i, 9]))))
    conv = convergence_from_files(outdir, run["id"])
    for r_ in rows:
        r_["cnvg_flag_raw"] = r_.pop("converged")
        r_["converged"] = conv["per_ls"].get(r_["ls"] + 1, False)
    fz = outdir / "hist_force.txt"
    if fz.exists():
        Fz = np.loadtxt(fz, ndmin=2)
        for i, r_ in enumerate(rows):
            r_.update(sz_cu_top_layer=Fz[i, 1], sz_cu_mid_layer=Fz[i, 2], sz_cu_bottom_layer=Fz[i, 3])
    rec.update(ok=True, n_elems=int(mt[0]), n_nodes=int(mt[1]), rows=rows, convergence=conv,
               frame_check=dict(rsys=0, cu_elem_esys=[int(mt[2]), int(mt[3])], elem1_mat=int(mt[4]),
                                global_frame=bool(mt[2] == 0 and mt[3] == 0 and mt[4] == 1)),
               cpu_total_s=float(sum(r["cpu_s"] for r in rows)),
               all_converged=all(r["converged"] for r in rows))
    return rec


def convergence_from_files(outdir, rid):
    """Per-load-step completion from <job>.mntr; failure strings from the .out log.

    MAPDL's *GET,...,SOLU,CNVG returns 0 for linear (single-iteration) solutions, so it is
    recorded only as cnvg_flag_raw and not used for the verdict."""
    per_ls, bad = {}, []
    mntr = next((p for p in outdir.glob("*.mntr")), None)
    if mntr is not None:
        for line in mntr.read_text(errors="replace").splitlines():
            parts = line.split()
            if len(parts) > 6 and parts[0].isdigit() and parts[1].isdigit():
                per_ls[int(parts[0])] = True
    out = next((p for p in outdir.glob("*.out")), None)
    if out is not None:
        txt = out.read_text(errors="replace")
        for pat in ("SOLUTION NOT CONVERGED", "not converged", "Solution not converged", "*** ERROR ***"):
            if pat in txt:
                bad.append(pat)
    ok = bool(per_ls) and all(per_ls.get(i, False) for i in range(1, len(HISTORY) + 1)) and not bad
    if bad:
        per_ls = {k: False for k in per_ls}
    return dict(per_ls=per_ls, failures=bad, all_load_steps_completed=ok,
                source=str(mntr.name) if mntr else None)


def write_aggregate():
    recs = []
    for d in sorted(RESULTS.glob("R*")):
        f = d / "record.json"
        if f.exists():
            recs.append(json.loads(f.read_text()))
    (RESULTS / "results.json").write_text(json.dumps(recs, indent=1, default=float))
    cols = ["id", "pad", "Ez_GPa", "model", "variant", "plastic_applied", "mesh", "shape", "cap", "sliding",
            "cte", "angle", "T_sf", "n_nodes", "n_elems", "exit_code", "wall_s", "cpu_total_s", "all_converged",
            "global_frame",
            "ls", "tag", "branch", "T", "eps_dev_zz", "epel_x", "epel_y", "epel_z", "eppl_eqv_avg", "eppl_eqv_max",
            "s_x", "s_y", "s_z", "seqv_max_cu", "sidewall_sn_max", "sidewall_sn_min", "sidewall_tau_max",
            "u_center_nm", "u_edge_nm", "u_center_rel_padbottom_nm", "u_center_rel_si_nm", "u_avg_rel_si_nm",
            "u_surface_edge_rel_si_nm", "u_surface_center_rel_si_nm",
            "sz_cu_top_layer", "sz_cu_mid_layer", "sz_cu_bottom_layer", "cnvg_flag_raw"]
    with open(RESULTS / "results.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for r in recs:
            if not r.get("ok"):
                continue
            run = r["deck"]["run"]
            base = [r["id"], run["pad"], r["deck"]["mat"]["Ez_GPa"], run["model"], run["variant"],
                    r["deck"]["mat"]["plastic"], run["mesh"], run["shape"], run["cap"], run["sliding"],
                    run["cte"], run["angle"], r["deck"]["T_sf"], r["n_nodes"], r["n_elems"],
                    r["exec"]["exit_code"], round(r["exec"]["wall_s"], 2), round(r["cpu_total_s"], 2),
                    r["all_converged"], r.get("frame_check", {}).get("global_frame", False)]
            for row in r["rows"]:
                w.writerow(base + [row.get(c, "") for c in cols[len(base):]])
    return len(recs)


def execute(run_ids, mode, nproc, dry, force):
    RESULTS.mkdir(exist_ok=True)
    r01_hash = None
    for rid in run_ids:
        run = RUN_BY_ID[rid]
        st, why = run_status(run)
        outdir = RESULTS / rid
        if st == "BLOCKED":
            print(f"[{rid}] BLOCKED: {why}")
            continue
        if st == "DONE" and not force and not dry:
            print(f"[{rid}] already done (use --force to re-run)")
            continue
        text, meta = _write_case(run, outdir)
        print(f"[{rid}] {run['note']}: {meta['n_nodes']} nodes, {meta['n_elems']} elements, "
              f"E_z = {meta['mat']['Ez_GPa']:.2f} GPa, T_sf = {meta['T_sf']:g} C, plastic = {meta['mat']['plastic']}")
        if rid == "R08":
            ref_text, ref_meta = build_deck(RUN_BY_ID["R01"])
            same = ref_meta["physics_sha256"] == meta["physics_sha256"]
            (outdir / "R08_check.json").write_text(json.dumps(dict(identical_to_R01=same,
                sha_R08=meta["physics_sha256"], sha_R01=ref_meta["physics_sha256"]), indent=1))
            print(f"[R08] deck physics identical to R01: {same} -> {'no solve needed' if same else 'SOLVE REQUIRED'}")
            if same:
                continue
        if dry:
            continue
        if mode in ("auto", "pymapdl"):
            try:
                info = run_pymapdl(run, outdir, nproc)
            except Exception as exc:
                if mode == "pymapdl":
                    raise
                print(f"[{rid}] PyMAPDL unavailable ({type(exc).__name__}: {exc}); falling back to batch")
                info = run_batch(run, outdir, nproc)
        else:
            info = run_batch(run, outdir, nproc)
        rec = collect(run, outdir, info)
        (outdir / "record.json").write_text(json.dumps(rec, indent=1, default=float))
        if rec["ok"]:
            h400 = next(r for r in rec["rows"] if r["tag"] == "h400")
            print(f"[{rid}] exit {info['exit_code']}, wall {info['wall_s']:.1f} s, converged {rec['all_converged']}; "
                  f"u_center(30->400) = {h400['u_center_nm']:.3f} nm, eps'_zz(400) = {1e3 * h400['eps_dev_zz']:+.3f}e-3")
            print(f"[{rid}]   at 400 C: Cu sigma_z vol-avg = {h400['s_z']:.1f} MPa, mid-height sigma_z = "
                  f"{h400.get('sz_cu_mid_layer', float('nan')):.1f} MPa, pad lengthening = {h400['u_center_rel_padbottom_nm']:.3f} nm")
        else:
            print(f"[{rid}] FAILED (exit {info['exit_code']}): {rec.get('error')} - see {outdir / (rid + '.out')}")
    if not dry:
        n = write_aggregate()
        print(f"aggregate: {n} run records -> {RESULTS / 'results.csv'}")
    if _MAPDL is not None:
        try:
            _MAPDL.exit()
        except Exception:
            pass


def main(argv=None):
    ap = argparse.ArgumentParser(description="Stage 1 Ayoub 2022 reproduction")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--case", default=None, help="R01 | R01,R02 | all")
    ap.add_argument("--step", type=int, default=None, help="run-order step 1..8")
    ap.add_argument("--mode", choices=["auto", "batch", "pymapdl"], default="batch")
    ap.add_argument("--np", dest="nproc", type=int, default=DEFAULT_NP)
    ap.add_argument("--exe", default=None)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--recollect", action="store_true", help="re-parse finished runs (no solve)")
    a = ap.parse_args(argv)
    global EXE
    if a.exe:
        EXE = Path(a.exe)
    if a.list:
        list_runs(); return 0
    if a.recollect:
        for dd in sorted(RESULTS.glob("R*")):
            if (dd / "hist_disp.txt").exists() and dd.name in RUN_BY_ID:
                old = json.loads((dd / "record.json").read_text()) if (dd / "record.json").exists() else {}
                rec = collect(RUN_BY_ID[dd.name], dd, old.get("exec", dict(method="?", exit_code=None, wall_s=float("nan"))))
                (dd / "record.json").write_text(json.dumps(rec, indent=1, default=float))
                print(f"[{dd.name}] re-collected: all_converged = {rec.get('all_converged')}")
        print(f"aggregate: {write_aggregate()} run records"); return 0
    if a.step is not None:
        ids = [r["id"] for r in RUNS if r["step"] == a.step]
    elif a.case:
        ids = [r["id"] for r in RUNS] if a.case == "all" else [s.strip() for s in a.case.split(",")]
    else:
        ap.error("give --list, --case or --step")
    unknown = [i for i in ids if i not in RUN_BY_ID]
    if unknown:
        ap.error(f"unknown run ids {unknown}")
    if not a.dry_run and not EXE.exists():
        print(f"MAPDL executable not found: {EXE}"); return 2
    print(f"{_dt.datetime.now():%Y-%m-%d %H:%M:%S}  runs: {ids}  mode={a.mode} np={a.nproc}  exe={EXE}")
    execute(ids, a.mode, a.nproc, a.dry_run, a.force)
    return 0


if __name__ == "__main__":
    sys.exit(main())
