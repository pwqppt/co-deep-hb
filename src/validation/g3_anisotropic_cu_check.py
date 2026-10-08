#!/usr/bin/env python
"""
G3 validation -- anisotropic Cu elasticity and crystal orientation in MAPDL.

A Cu cube with cubic elastic constants (TB,ANEL) is rotated so that the
crystal [100], [110] and [111] directions each align with global z. A uniform
uniaxial stress sigma_zz = S0 is applied, and the effective modulus
E_z = sigma_zz / eps_zz is compared with the closed-form compliance result

    1/E = S11 - 2 (S11 - S12 - S44/2) (l^2 m^2 + m^2 n^2 + n^2 l^2)

PASS criterion (G3): |E_z,FEA - E_z,exact| / E_z,exact <= 1 % for every case.

Two independent orientation routes are solved for each direction:

  esys   : crystal-frame ANEL matrix + element coordinate system (CSKP + ESYS)
           -- this is how the 3D anisotropic pad model will assign orientation.
  prerot : stiffness rotated in Python (4th-order tensor rotation), entered
           directly in the global frame with ESYS,0
           -- validates the tensor math the surrogate/2D pipeline will reuse.

Secondary diagnostics (reported, not part of the G3 verdict): the transverse
strain ratios -eps_xx/eps_zz and -eps_yy/eps_zz, which depend on the in-plane
crystal axes and therefore catch a wrong frame even when E_z happens to match
(Cu loaded along [110] is auxetic towards [1-10]: nu ~ -0.14).

Units: uMKS, consistent -- length um, mass kg, time s, temperature degC
       => force uN, stress MPa. So C11 = 169.1 GPa is entered as 169100.

Boundary conditions: equilibrated tension on z = 0 and z = L, plus a 3-2-1
point restraint (6 rigid-body DOF only). Reactions are zero, so the stress
state is exactly uniaxial for ANY orientation -- no face is forced to stay
planar, so this set-up is also valid for later random-orientation checks.

USAGE (run natively, NOT inside the Claude Science sandbox -- the ANSYS
licensing client cannot load there):

    python g3_anisotropic_cu_check.py                     # auto: PyMAPDL, fall back to batch
    python g3_anisotropic_cu_check.py --mode batch        # batch only
    python g3_anisotropic_cu_check.py --mode both         # run both, compare
    python g3_anisotropic_cu_check.py --dry-run           # write decks + references, no solve
    python g3_anisotropic_cu_check.py --ndiv 6 --out g3_results_n6

Requires: numpy; ansys-mapdl-core (only for the PyMAPDL route).
Outputs (in --out, default ./g3_results):
    results.json   all per-case numbers, references, verdicts, metadata
    results.csv    one flat row per case
    decks/         the exact APDL decks that were solved
    pymapdl_work/, batch_work/   MAPDL working directories (*.out logs)
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import json
import os
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
AWP_ROOT = Path(os.environ.get(
    "AWP_ROOT261", r"C:\Program Files\ANSYS Inc\ANSYS Student\v261"))
MAPDL_EXE = AWP_ROOT / "ansys" / "bin" / "winx64" / "ANSYS261.exe"

# Cu cubic elastic constants (GPa) -> MPa for uMKS
C11_GPA, C12_GPA, C44_GPA = 169.1, 122.2, 75.4
GPA_TO_MPA = 1000.0

L_CUBE = 1.0          # um, cube edge
S0 = 100.0            # MPa, applied uniaxial tension (linear-elastic, small strain)
TOL_PCT = 1.0         # G3 tolerance on E_z
ZERO_TOL = 1e-4       # absolute tolerance for ratios whose exact value is 0

ORIENTATIONS = {      # loading direction in crystal frame + in-plane x axis
    "100": ((1, 0, 0), (0, 1, 0)),
    "110": ((1, 1, 0), (0, 0, 1)),     # x = [001], y = [1-10]  (orthotropic frame)
    "111": ((1, 1, 1), (1, -1, 0)),    # x = [1-10], y = [11-2]
}
ROUTES = ("esys", "prerot")

# MAPDL Voigt order: x, y, z, xy, yz, xz  (engineering shear strain)
VOIGT = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]
VOIGT_LABELS = ["X", "Y", "Z", "XY", "YZ", "XZ"]


# ----------------------------------------------------------------------------
# Analytic elasticity (pure numpy, no ANSYS)
# ----------------------------------------------------------------------------
def cubic_stiffness_voigt(c11, c12, c44):
    """6x6 stiffness (engineering-shear Voigt, MAPDL order) in the crystal frame."""
    D = np.zeros((6, 6))
    D[:3, :3] = c12
    np.fill_diagonal(D[:3, :3], c11)
    D[3, 3] = D[4, 4] = D[5, 5] = c44
    return D


def voigt_to_tensor(D):
    C = np.zeros((3, 3, 3, 3))
    for I, (i, j) in enumerate(VOIGT):
        for J, (k, l) in enumerate(VOIGT):
            for a, b in ((i, j), (j, i)):
                for c, d in ((k, l), (l, k)):
                    C[a, b, c, d] = D[I, J]
    return C


def tensor_to_voigt(C):
    D = np.zeros((6, 6))
    for I, (i, j) in enumerate(VOIGT):
        for J, (k, l) in enumerate(VOIGT):
            D[I, J] = C[i, j, k, l]
    return D


def rotation_crystal_to_global(load_dir, x_dir):
    """R with R @ v_crystal = v_global; maps crystal load_dir onto global +z.

    Rows of R are the global x, y, z axes expressed in crystal coordinates;
    columns of R are the crystal [100], [010], [001] axes in global coordinates.
    """
    r3 = np.asarray(load_dir, float); r3 /= np.linalg.norm(r3)
    r1 = np.asarray(x_dir, float);    r1 /= np.linalg.norm(r1)
    assert abs(r1 @ r3) < 1e-12, "in-plane x axis must be perpendicular to load axis"
    r2 = np.cross(r3, r1)
    R = np.vstack([r1, r2, r3])
    assert np.allclose(R @ R.T, np.eye(3)) and np.isclose(np.linalg.det(R), 1.0)
    return R


def rotate_stiffness(D_crys, R):
    C = voigt_to_tensor(D_crys)
    Cg = np.einsum("ip,jq,kr,ls,pqrs->ijkl", R, R, R, R, C, optimize=True)
    return tensor_to_voigt(Cg)


def closed_form_Ez(c11, c12, c44, direction):
    """Textbook cubic compliance formula (independent of the tensor route)."""
    den = (c11 - c12) * (c11 + 2 * c12)
    s11 = (c11 + c12) / den
    s12 = -c12 / den
    s44 = 1.0 / c44
    l, m, n = np.asarray(direction, float) / np.linalg.norm(direction)
    inv_E = s11 - 2 * (s11 - s12 - s44 / 2) * (l*l*m*m + m*m*n*n + n*n*l*l)
    return 1.0 / inv_E, dict(S11=s11, S12=s12, S44=s44)


def analytic_reference(key):
    """All exact quantities for one orientation (MPa units)."""
    load_dir, x_dir = ORIENTATIONS[key]
    c11, c12, c44 = (v * GPA_TO_MPA for v in (C11_GPA, C12_GPA, C44_GPA))
    D_crys = cubic_stiffness_voigt(c11, c12, c44)
    R = rotation_crystal_to_global(load_dir, x_dir)
    D_glob = rotate_stiffness(D_crys, R)
    S_glob = np.linalg.inv(D_glob)                     # engineering compliance
    E_formula, S_crys = closed_form_Ez(c11, c12, c44, load_dir)
    E_tensor = 1.0 / S_glob[2, 2]
    # the two analytic routes must agree to round-off before we trust either
    assert abs(E_tensor - E_formula) / E_formula < 1e-10, (key, E_tensor, E_formula)
    strain_per_S0 = S_glob[:, 2]                       # strain response to sigma_zz = 1
    return dict(
        R=R, D_crys=D_crys, D_glob=D_glob,
        E_formula_MPa=E_formula, E_tensor_MPa=E_tensor,
        nu_zx=-strain_per_S0[0] / strain_per_S0[2],
        nu_zy=-strain_per_S0[1] / strain_per_S0[2],
        gam_xy_ratio=strain_per_S0[3] / strain_per_S0[2],
        gam_yz_ratio=strain_per_S0[4] / strain_per_S0[2],
        gam_xz_ratio=strain_per_S0[5] / strain_per_S0[2],
        S11=S_crys["S11"], S12=S_crys["S12"], S44=S_crys["S44"],
    )


# ----------------------------------------------------------------------------
# APDL deck generation
# ----------------------------------------------------------------------------
def _f(x):
    return f"{x:.12e}"


def anel_commands(D):
    """TB,ANEL in stiffness form: C1..C21 = column-wise lower triangle."""
    vals = [D[i, j] for j in range(6) for i in range(j, 6)]
    assert len(vals) == 21
    lines = ["TB,ANEL,1,1,,0", "TBTEMP,25"]
    for start in range(0, 21, 6):
        chunk = ",".join(_f(v) for v in vals[start:start + 6])
        lines.append(f"TBDATA,{start + 1},{chunk}")
    return lines


RESULT_KEYS = [
    "NELEM", "NNODE", "NTOP", "NBOT", "CPU_SOLVE",
    "UZ_TOP", "UZ_BOT",
    "SX_AVG", "SY_AVG", "SZ_AVG", "SXY_AVG", "SYZ_AVG", "SXZ_AVG",
    "EX_AVG", "EY_AVG", "EZ_AVG", "EXY_AVG", "EYZ_AVG", "EXZ_AVG",
    "SZ_DEVMAX",
]


def build_case_deck(key, route, ndiv):
    ref = analytic_reference(key)
    tag = f"{route}_{key}"
    L = L_CUBE
    D = ref["D_crys"] if route == "esys" else ref["D_glob"]
    R = ref["R"]
    a1, a2 = R[:, 0], R[:, 1]          # crystal [100], [010] in global coords

    d = []
    d += [f"/TITLE,G3 Cu cube {route} [{key}] || z  (uMKS: um, MPa)",
          "/PREP7",
          "ET,1,SOLID185",
          "MP,DENS,1,8.96e-15     ! kg/um^3 (unused in static, keeps material complete)"]
    d += anel_commands(D)
    d += [f"BLOCK,0,{_f(L)},0,{_f(L)},0,{_f(L)}"]
    if route == "esys":
        d += ["! crystal frame -> local CS 11 from three keypoints",
              "K,1001,0,0,0",
              f"K,1002,{_f(a1[0])},{_f(a1[1])},{_f(a1[2])}    ! crystal [100]",
              f"K,1003,{_f(a2[0])},{_f(a2[1])},{_f(a2[2])}    ! crystal [010] (xy-plane)",
              "CSKP,11,0,1001,1002,1003",
              "CSYS,0",
              "ESYS,11"]
    else:
        d += ["ESYS,0                 ! stiffness already rotated to global"]
    d += ["MAT,1",
          "TYPE,1",
          f"ESIZE,{_f(L / ndiv)}",
          "MSHAPE,0,3D",
          "MSHKEY,1",
          "VMESH,ALL",
          "*GET,NELEM,ELEM,0,COUNT",
          "*GET,NNODE,NODE,0,COUNT",
          "FINISH",
          "",
          "/SOLU",
          "ANTYPE,STATIC",
          "CSYS,0",
          f"S0={_f(S0)}",
          f"NSEL,S,LOC,Z,{_f(L)}",
          "SF,ALL,PRES,-S0        ! negative pressure = tension",
          "NSEL,S,LOC,Z,0",
          "SF,ALL,PRES,-S0",
          "ALLSEL,ALL",
          "! 3-2-1 restraint: removes 6 rigid-body modes, zero reactions",
          "NA=NODE(0,0,0)",
          f"NB=NODE({_f(L)},0,0)",
          f"NC=NODE(0,{_f(L)},0)",
          "D,NA,UX,0,,,,UY,UZ",
          "D,NB,UY,0,,,,UZ",
          "D,NC,UZ,0",
          "*GET,TCPU0,ACTIVE,0,TIME,CPU",
          "SOLVE",
          "*GET,TCPU1,ACTIVE,0,TIME,CPU",
          "CPU_SOLVE=TCPU1-TCPU0",
          "FINISH",
          "",
          "/POST1",
          "SET,LAST",
          "RSYS,0                 ! report in global frame",
          "CSYS,0",
          "! ---- face-averaged UZ ----",
          f"NSEL,S,LOC,Z,{_f(L)}",
          "*GET,NTOP,NODE,0,COUNT",
          "UZ_TOP=0",
          "ND=0",
          "*DO,II,1,NTOP",
          "ND=NDNEXT(ND)",
          "UZ_TOP=UZ_TOP+UZ(ND)",
          "*ENDDO",
          "UZ_TOP=UZ_TOP/NTOP",
          "NSEL,S,LOC,Z,0",
          "*GET,NBOT,NODE,0,COUNT",
          "UZ_BOT=0",
          "ND=0",
          "*DO,II,1,NBOT",
          "ND=NDNEXT(ND)",
          "UZ_BOT=UZ_BOT+UZ(ND)",
          "*ENDDO",
          "UZ_BOT=UZ_BOT/NBOT",
          "ALLSEL,ALL",
          "! ---- nodal-averaged stress and elastic strain (global frame) ----"]
    comps = list(zip(["X", "Y", "Z", "XY", "YZ", "XZ"],
                     ["X", "Y", "Z", "XY", "YZ", "XZ"]))
    for c, _ in comps:
        d += [f"S{c}_AVG=0", f"E{c}_AVG=0"]
    d += ["SZ_DEVMAX=0", "ND=0", "*DO,II,1,NNODE", "ND=NDNEXT(ND)"]
    for c, _ in comps:
        d += [f"*GET,VTMP,NODE,ND,S,{c}", f"S{c}_AVG=S{c}_AVG+VTMP"]
        if c == "Z":
            d += ["DTMP=ABS(VTMP-S0)",
                  "*IF,DTMP,GT,SZ_DEVMAX,THEN",
                  "SZ_DEVMAX=DTMP",
                  "*ENDIF"]
        d += [f"*GET,VTMP,NODE,ND,EPEL,{c}", f"E{c}_AVG=E{c}_AVG+VTMP"]
    d += ["*ENDDO"]
    for c, _ in comps:
        d += [f"S{c}_AVG=S{c}_AVG/NNODE", f"E{c}_AVG=E{c}_AVG/NNODE"]
    d += ["", f"*CFOPEN,g3_{tag},txt"]
    for k in RESULT_KEYS:
        d += [f"*VWRITE,{k}", f"('{k} ',E24.15)"]
    d += ["*CFCLOS", "FINISH", ""]
    return tag, "\n".join(d) + "\n"


def all_cases():
    return [(route, key) for route in ROUTES for key in ORIENTATIONS]


# ----------------------------------------------------------------------------
# Result parsing and evaluation
# ----------------------------------------------------------------------------
def parse_result_file(path):
    vals = {}
    for line in Path(path).read_text(errors="replace").splitlines():
        m = re.match(r"\s*([A-Z_0-9]+)\s+(\S+)", line)
        if m:
            vals[m.group(1)] = float(m.group(2).replace("D", "E"))
    missing = [k for k in RESULT_KEYS if k not in vals]
    if missing:
        raise ValueError(f"{path}: missing keys {missing}")
    return vals


def _rel_or_abs_ok(fea, ref):
    if abs(ref) < 1e-8:
        return abs(fea) <= ZERO_TOL, abs(fea)
    err = abs(fea - ref) / abs(ref) * 100.0
    return err <= TOL_PCT, err


def evaluate(route, key, raw, wall_s, method):
    ref = analytic_reference(key)
    ez = raw["EZ_AVG"]
    E_strain = raw["SZ_AVG"] / ez
    E_disp = S0 * L_CUBE / (raw["UZ_TOP"] - raw["UZ_BOT"])
    E_ref = ref["E_formula_MPa"]
    err_disp = abs(E_disp - E_ref) / E_ref * 100.0
    err_strain = abs(E_strain - E_ref) / E_ref * 100.0
    g3_pass = err_disp <= TOL_PCT and err_strain <= TOL_PCT

    diag = {}
    for name, fea in [("nu_zx", -raw["EX_AVG"] / ez),
                      ("nu_zy", -raw["EY_AVG"] / ez),
                      ("gam_xy_ratio", raw["EXY_AVG"] / ez),
                      ("gam_yz_ratio", raw["EYZ_AVG"] / ez),
                      ("gam_xz_ratio", raw["EXZ_AVG"] / ez)]:
        ok, err = _rel_or_abs_ok(fea, ref[name])
        diag[name] = dict(fea=fea, exact=ref[name], err=err,
                          err_kind="abs" if abs(ref[name]) < 1e-8 else "pct", ok=ok)
    nonuni = max(abs(raw[k]) for k in
                 ("SX_AVG", "SY_AVG", "SXY_AVG", "SYZ_AVG", "SXZ_AVG")) / S0
    diag["stress_nonuniaxiality"] = dict(fea=nonuni, exact=0.0, err=nonuni,
                                         err_kind="abs", ok=nonuni <= ZERO_TOL)
    diag["sz_uniformity"] = dict(fea=raw["SZ_DEVMAX"] / S0, exact=0.0,
                                 err=raw["SZ_DEVMAX"] / S0, err_kind="abs",
                                 ok=raw["SZ_DEVMAX"] / S0 <= ZERO_TOL)
    return dict(
        method=method, route=route, orientation=f"[{key}]",
        E_exact_GPa=E_ref / 1000, E_disp_GPa=E_disp / 1000,
        E_strain_GPa=E_strain / 1000,
        err_disp_pct=err_disp, err_strain_pct=err_strain, g3_pass=g3_pass,
        diagnostics=diag, diagnostics_pass=all(v["ok"] for v in diag.values()),
        n_elem=int(round(raw["NELEM"])), n_node=int(round(raw["NNODE"])),
        cpu_solve_s=raw["CPU_SOLVE"], wall_s=wall_s, raw=raw,
    )


# ----------------------------------------------------------------------------
# Solvers
# ----------------------------------------------------------------------------
def run_pymapdl(decks, workdir):
    print("\n[PyMAPDL] launching ...")
    try:
        from ansys.mapdl.core import launch_mapdl
    except ImportError as exc:
        print(f"  ansys-mapdl-core not importable: {exc}")
        return None, {}
    workdir.mkdir(parents=True, exist_ok=True)
    meta, rows, mapdl = {}, [], None
    try:
        t0 = time.perf_counter()
        mapdl = launch_mapdl(exec_file=str(MAPDL_EXE), run_location=str(workdir),
                             jobname="g3", nproc=1, additional_switches="-smp",
                             override=True, start_timeout=180)
        meta = dict(launch_s=time.perf_counter() - t0, version=str(mapdl.version))
        print(f"  connected, MAPDL {meta['version']} (launch {meta['launch_s']:.1f} s)")
        for (route, key), (tag, text) in decks.items():
            deck = workdir / f"g3_{tag}.dat"
            deck.write_text(text)
            res = workdir / f"g3_{tag}.txt"
            res.unlink(missing_ok=True)
            mapdl.clear()
            t1 = time.perf_counter()
            mapdl.input(str(deck))
            wall = time.perf_counter() - t1
            rows.append(evaluate(route, key, parse_result_file(res), wall, "pymapdl"))
            print(f"  {tag:12s} done in {wall:.2f} s")
        return rows, meta
    except Exception as exc:
        print(f"  FAILED: {type(exc).__name__}: {exc}")
        meta["error"] = f"{type(exc).__name__}: {exc}"
        return None, meta
    finally:
        if mapdl is not None:
            try:
                mapdl.exit()
            except Exception:
                pass


def run_batch(decks, workdir):
    print("\n[batch] running one MAPDL job for all cases ...")
    if not MAPDL_EXE.exists():
        print(f"  solver not found: {MAPDL_EXE}")
        return None, {"error": "solver not found"}
    workdir.mkdir(parents=True, exist_ok=True)
    master = []
    for (route, key), (tag, text) in decks.items():
        (workdir / f"g3_{tag}.txt").unlink(missing_ok=True)
        master += ["/CLEAR,NOSTART", text]
    master.append("/EXIT,NOSAVE\n")
    inp = workdir / "g3_master.dat"
    inp.write_text("\n".join(master))
    out = workdir / "g3.out"
    cmd = [str(MAPDL_EXE), "-b", "-smp", "-np", "1", "-dir", str(workdir),
           "-j", "g3", "-i", str(inp), "-o", str(out)]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, errors="replace")
    wall = time.perf_counter() - t0
    meta = dict(returncode=proc.returncode, total_wall_s=wall)
    banner = re.search(r"Ansys MAPDL\s+(\S+\s+\S+)\s+Build\s+(\S+)",
                       out.read_text(errors="replace")) if out.exists() else None
    if banner:
        meta["version"] = f"{banner.group(1)} build {banner.group(2)}"
    rows = []
    for (route, key), (tag, _) in decks.items():
        res = workdir / f"g3_{tag}.txt"
        if not res.exists():
            tail = out.read_text(errors="replace")[-2000:] if out.exists() else proc.stderr
            print(f"  missing result for {tag} (rc={proc.returncode}); log tail:\n{tail}")
            meta["error"] = f"missing result for {tag}"
            return None, meta
        rows.append(evaluate(route, key, parse_result_file(res), None, "batch"))
    print(f"  done: {len(rows)} cases, total wall {wall:.1f} s (incl. license checkout)")
    return rows, meta


# ----------------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------------
def print_table(rows):
    hdr = (f"{'method':8s} {'route':7s} {'dir':6s} {'E_exact':>9s} {'E_disp':>9s} "
           f"{'err%':>8s} {'E_strain':>9s} {'err%':>8s} {'G3':>5s} {'diag':>5s} "
           f"{'nel':>5s} {'cpu_s':>7s}")
    print("\n" + hdr + "\n" + "-" * len(hdr))
    for r in rows:
        print(f"{r['method']:8s} {r['route']:7s} {r['orientation']:6s} "
              f"{r['E_exact_GPa']:9.3f} {r['E_disp_GPa']:9.3f} {r['err_disp_pct']:8.4f} "
              f"{r['E_strain_GPa']:9.3f} {r['err_strain_pct']:8.4f} "
              f"{'PASS' if r['g3_pass'] else 'FAIL':>5s} "
              f"{'ok' if r['diagnostics_pass'] else 'CHECK':>5s} "
              f"{r['n_elem']:5d} {r['cpu_solve_s']:7.3f}")


def write_outputs(out, rows, meta):
    out.mkdir(parents=True, exist_ok=True)

    def _clean(o):
        if isinstance(o, dict):
            return {k: _clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [_clean(v) for v in o]
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.floating, np.integer, np.bool_)):
            return o.item()
        return o

    (out / "results.json").write_text(json.dumps(_clean(dict(meta=meta, cases=rows)), indent=2))
    with open(out / "results.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["method", "route", "orientation", "E_exact_GPa", "E_disp_GPa",
                    "err_disp_pct", "E_strain_GPa", "err_strain_pct", "g3_pass",
                    "nu_zx_fea", "nu_zx_exact", "nu_zy_fea", "nu_zy_exact",
                    "diagnostics_pass", "n_elem", "n_node", "cpu_solve_s", "wall_s"])
        for r in rows:
            dg = r["diagnostics"]
            w.writerow([r["method"], r["route"], r["orientation"],
                        f"{r['E_exact_GPa']:.6f}", f"{r['E_disp_GPa']:.6f}",
                        f"{r['err_disp_pct']:.6f}", f"{r['E_strain_GPa']:.6f}",
                        f"{r['err_strain_pct']:.6f}", r["g3_pass"],
                        f"{dg['nu_zx']['fea']:.6f}", f"{dg['nu_zx']['exact']:.6f}",
                        f"{dg['nu_zy']['fea']:.6f}", f"{dg['nu_zy']['exact']:.6f}",
                        r["diagnostics_pass"], r["n_elem"], r["n_node"],
                        f"{r['cpu_solve_s']:.4f}",
                        "" if r["wall_s"] is None else f"{r['wall_s']:.3f}"])


def reference_summary():
    out = {}
    for key in ORIENTATIONS:
        ref = analytic_reference(key)
        out[f"[{key}]"] = dict(
            E_exact_GPa=ref["E_formula_MPa"] / 1000,
            nu_zx=ref["nu_zx"], nu_zy=ref["nu_zy"],
            R_crystal_to_global=ref["R"].tolist())
    return out


# ----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--mode", choices=["auto", "pymapdl", "batch", "both"], default="auto")
    ap.add_argument("--ndiv", type=int, default=4, help="elements per cube edge")
    ap.add_argument("--out", default="g3_results")
    ap.add_argument("--dry-run", action="store_true",
                    help="write decks and analytic references only")
    a = ap.parse_args(argv)

    out = Path(a.out).resolve()
    decks = {(r, k): build_case_deck(k, r, a.ndiv) for r, k in all_cases()}
    (out / "decks").mkdir(parents=True, exist_ok=True)
    for tag, text in decks.values():
        (out / "decks" / f"g3_{tag}.dat").write_text(text)

    meta = dict(
        test="G3 anisotropic Cu elasticity + orientation",
        timestamp=_dt.datetime.now().isoformat(timespec="seconds"),
        host=platform.node(), python=sys.version.split()[0],
        mapdl_exe=str(MAPDL_EXE), mode=a.mode, ndiv=a.ndiv,
        units="uMKS (um, kg, s, degC -> uN, MPa)",
        constants_GPa=dict(C11=C11_GPA, C12=C12_GPA, C44=C44_GPA),
        S0_MPa=S0, L_um=L_CUBE, tol_pct=TOL_PCT, zero_tol=ZERO_TOL,
        references=reference_summary())

    print("Analytic references (closed form == tensor rotation, asserted):")
    for k, v in meta["references"].items():
        print(f"  {k}: E_z = {v['E_exact_GPa']:.4f} GPa   nu_zx = {v['nu_zx']:+.4f}"
              f"   nu_zy = {v['nu_zy']:+.4f}")
    if a.dry_run:
        (out / "references.json").write_text(json.dumps(meta, indent=2))
        print(f"\nDry run: decks + references written to {out}")
        return 0

    rows = []
    if a.mode in ("auto", "pymapdl", "both"):
        r, m = run_pymapdl(decks, out / "pymapdl_work")
        meta["pymapdl"] = m
        rows += r or []
    if a.mode == "batch" or a.mode == "both" or (a.mode == "auto" and not rows):
        r, m = run_batch(decks, out / "batch_work")
        meta["batch"] = m
        rows += r or []

    if not rows:
        meta["verdict"] = "NO_RESULTS"
        write_outputs(out, rows, meta)
        print("\nNo method produced results; see results.json and *.out logs.")
        return 2

    meta["verdict"] = "PASS" if all(r["g3_pass"] for r in rows) else "FAIL"
    meta["diagnostics_verdict"] = ("PASS" if all(r["diagnostics_pass"] for r in rows)
                                   else "CHECK")
    write_outputs(out, rows, meta)
    print_table(rows)
    print(f"\nG3 verdict: {meta['verdict']}   (diagnostics: {meta['diagnostics_verdict']})")
    print(f"Results written to {out}")
    return 0 if meta["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
