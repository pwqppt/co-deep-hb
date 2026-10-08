"""
Crystal plasticity (TB,XTAL) check for MAPDL Student 2026 R1 (v261), single-crystal Cu.

Stages
  elastic : resolve two documentation ambiguities by comparing E_z against the
            closed form at 0.02 % strain (below yield):
              - TB,ELAS (OELM) C4-C6: G or 2G ("2x the orthotropic shear modulus")
              - TB,XTAL ORIE Euler angles: crystal->global (Bunge) or inverse
  cp      : Step 2/3 - displacement-controlled tension to 1 % along z, Schmid check
  ludwick : Step 4  - TB,ANEL (rotated) + isotropic MISO from a Ludwick curve

Run natively (NOT inside the Claude Science sandbox):
    python cp_check.py                 # all stages
    python cp_check.py --stage cp      # one stage

Units uMKS: um, kg, s, degC -> uN, MPa, pJ.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import time
from pathlib import Path

import numpy as np

EXE = r"C:\Program Files\ANSYS Inc\ANSYS Student\v261\ansys\bin\winx64\ANSYS261.exe"
HERE = Path(__file__).resolve().parent
RUN = HERE / "run"

# ---------------------------------------------------------------- inputs
C11, C12, C44 = 168400.0, 121400.0, 75400.0      # MPa
TAU_C = 80.0                                     # MPa
L = 10.0                                         # um cube edge
NDIV = 4                                         # 4x4x4 SOLID185
NSUB_CP = 100                                    # >= 50 substeps
EPS_END = 0.01
OFFSET = 0.002
TOL_PCT = 5.0

# CP flow / hardening (FCC). Thermal part of the slip resistance is 1 % of the
# total, so the flow stress is bracketed in [g_a, g] = [79.2, 80] MPa regardless
# of the activation-energy units; hardening saturates at 81 MPa (near-perfect).
KB_SI = 1.380649e-23
T_ABS = 298.15
FLFCC = dict(gdot0=1.0, p=1.0, q=1.0, ratio=0.01, dF=20.0 * KB_SI * T_ABS, cbya=1.0)
HARD = dict(g0=TAU_C, h0=10.0, gsat=81.0, r=1.0, n=0.0, qcross=1.4)

LUDWICK = dict(sigma0=143.0, K=675.0, N=0.68)

ORIENTS = {
    "[100]":         (np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])),
    "[110]":         (np.array([1.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0])),
    "[111]":         (np.array([1.0, 1.0, 1.0]), np.array([1.0, -1.0, 0.0])),
    "[1 0.1 0.02]":  (np.array([1.0, 0.1, 0.02]), np.array([0.0, 0.0, 1.0])),
}
EXPECTED_RATIO = {"[100]": 2.449, "[110]": 2.449, "[111]": 3.674, "[1 0.1 0.02]": 2.247}


# ---------------------------------------------------------------- crystal math
def compliance():
    S = np.linalg.inv(cubic_voigt())
    return S


def cubic_voigt():
    D = np.zeros((6, 6))
    D[:3, :3] = C12
    np.fill_diagonal(D[:3, :3], C11)
    D[3, 3] = D[4, 4] = D[5, 5] = C44
    return D


def E_dir(d):
    l, m, n = d / np.linalg.norm(d)
    S = compliance()
    s11, s12, s44 = S[0, 0], S[0, 1], S[3, 3]
    return 1.0 / (s11 - 2 * (s11 - s12 - s44 / 2) * (l*l*m*m + m*m*n*n + n*n*l*l))


def schmid_max(d):
    d = d / np.linalg.norm(d)
    normals = [np.array(v, float) for v in ([1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1])]
    best = 0.0
    for nrm in normals:
        nrm = nrm / np.linalg.norm(nrm)
        for b in itertools.product([-1, 0, 1], repeat=3):
            b = np.array(b, float)
            if np.count_nonzero(b) != 2 or abs(b @ nrm) > 1e-12:
                continue
            b /= np.linalg.norm(b)
            best = max(best, abs((nrm @ d) * (b @ d)))
    return best


def frame(load_dir, x_hint):
    """Rows = global x, y, z expressed in crystal coordinates (crystal->global map)."""
    z = load_dir / np.linalg.norm(load_dir)
    x = x_hint - (x_hint @ z) * z
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    return np.vstack([x, y, z])


def euler_zxz(M):
    """Angles (a, b, c) in degrees with M = Rz(a) Rx(b) Rz(c)."""
    b = math.degrees(math.acos(max(-1.0, min(1.0, M[2, 2]))))
    a = math.degrees(math.atan2(M[0, 2], -M[1, 2]))
    c = math.degrees(math.atan2(M[2, 0], M[2, 1]))
    return a, b, c


def rot_zxz(a, b, c):
    a, b, c = map(math.radians, (a, b, c))
    Rz = lambda t: np.array([[math.cos(t), -math.sin(t), 0], [math.sin(t), math.cos(t), 0], [0, 0, 1]])
    Rx = lambda t: np.array([[1, 0, 0], [0, math.cos(t), -math.sin(t)], [0, math.sin(t), math.cos(t)]])
    return Rz(a) @ Rx(b) @ Rz(c)


def euler_for(key, convention):
    """Euler angles for TB,XTAL ORIE.
    bunge  : columns of Rz Rx Rz = crystal axes in the global frame
    inverse: Rz Rx Rz maps global -> crystal instead
    """
    d, xh = ORIENTS[key]
    Q = frame(d, xh)               # v_global = Q v_crystal
    A = Q if convention == "bunge" else Q.T
    ang = euler_zxz(A)
    assert np.allclose(rot_zxz(*ang), A, atol=1e-10)
    return ang


VOIGT = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]   # MAPDL order x y z xy yz xz


def rotate_voigt(D, Q):
    C = np.zeros((3, 3, 3, 3))
    for I, (i, j) in enumerate(VOIGT):
        for J, (k, l) in enumerate(VOIGT):
            for (a, b) in {(i, j), (j, i)}:
                for (c, e) in {(k, l), (l, k)}:
                    C[a, b, c, e] = D[I, J]
    Cg = np.einsum("ia,jb,kc,ld,abcd->ijkl", Q, Q, Q, Q, C)
    return np.array([[Cg[i, j, k, l] for (k, l) in VOIGT] for (i, j) in VOIGT])


# ---------------------------------------------------------------- APDL decks
def f(x):
    return f"{x:.10g}"


def mat_cp(key, convention, gfac):
    S = compliance()
    E, nu, G = 1 / S[0, 0], -S[0, 1] / S[0, 0], C44
    a, b, c = euler_for(key, convention)
    h, fl = HARD, FLFCC
    return [
        "TB,ELAS,1,,,OELM",
        f"TBDATA,1,{f(E)},{f(E)},{f(E)},{f(gfac*G)},{f(gfac*G)},{f(gfac*G)}",
        f"TBDATA,7,{f(nu)},{f(nu)},{f(nu)}",
        "TB,XTAL,1,,,ORIE",   f"TBDATA,1,{f(a)},{f(b)},{f(c)}",
        "TB,XTAL,1,,,NSLFAM", "TBDATA,1,1",
        "TB,XTAL,1,,,FORM",   "TBDATA,1,1",
        "TB,XTAL,1,,,XPARAM", "TBDATA,1,1,0,1,12,0,1",
        "TB,XTAL,1,,,HARD",
        f"TBDATA,1,{f(h['g0'])},{f(h['h0'])},{f(h['gsat'])},{f(h['r'])},{f(h['n'])},{f(h['qcross'])}",
        "TB,XTAL,1,,,FLFCC",
        f"TBDATA,1,{f(fl['gdot0'])},{f(fl['p'])},{f(fl['q'])},{f(fl['ratio'])},{f(fl['dF'])},{f(fl['cbya'])}",
    ]


def ludwick_curve():
    lw = LUDWICK
    eps = [0.0, 0.0002, 0.0005, 0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.05, 0.1]
    return [(e, lw["sigma0"] + lw["K"] * e ** lw["N"]) for e in eps]


def mat_ludwick(key):
    d, xh = ORIENTS[key]
    Dg = rotate_voigt(cubic_voigt(), frame(d, xh))
    vals = [Dg[i, j] for j in range(6) for i in range(j, 6)]       # column-wise lower triangle
    lines = ["TB,ANEL,1,1,,0", "TBTEMP,25"]
    for s in range(0, 21, 6):
        lines.append(f"TBDATA,{s+1}," + ",".join(f(v) for v in vals[s:s+6]))
    pts = ludwick_curve()
    lines.append(f"TB,PLAS,1,1,{len(pts)},MISO")
    lines.append("TBTEMP,25")
    lines += [f"TBPT,DEFI,{f(e)},{f(s)}" for e, s in pts]
    return lines


def deck(case, mat_lines, eps_end, nsub):
    dz = eps_end * L
    return ["/PREP7", "ET,1,SOLID185", *mat_lines,
        f"BLOCK,0,{f(L)},0,{f(L)},0,{f(L)}", f"ESIZE,{f(L/NDIV)}", "VMESH,ALL", "FINISH",
        "/SOLU", "ANTYPE,STATIC", "NLGEOM,ON", "TOFFST,273.15", "TREF,25", "TUNIF,25",
        "NSEL,S,LOC,Z,0", "D,ALL,UZ,0",
        f"NSEL,S,LOC,Z,{f(L)}", f"D,ALL,UZ,{f(dz)}", "NSEL,ALL",
        "D,NODE(0,0,0),UX,0", "D,NODE(0,0,0),UY,0", f"D,NODE({f(L)},0,0),UY,0",
        "TIME,1", f"NSUBST,{nsub},{nsub*10},{nsub}", "OUTRES,ALL,ALL",
        "*GET,cpu0,ACTIVE,0,TIME,CPU", "SOLVE", "*GET,cpu1,ACTIVE,0,TIME,CPU",
        "*GET,cnv,ACTIVE,0,SOLU,CNVG", "dcpu=cpu1-cpu0", "FINISH",
        "/POST1", "SET,LAST", "*GET,tlast,ACTIVE,0,SET,TIME", "*GET,nset,ACTIVE,0,SET,SBST",
        "ALLSEL", "*GET,nel,ELEM,0,COUNT", "*GET,nnd,NODE,0,COUNT",
        "*DEL,res,,NOPR", "*DIM,res,ARRAY,nset,4",
        "*DO,i,1,nset",
        "SET,1,i", "*GET,t,ACTIVE,0,SET,TIME",
        "ALLSEL", "NSEL,S,LOC,Z,0", "*GET,nn,NODE,0,COUNT", "nd=0", "fz=0",
        "*DO,j,1,nn", "nd=NDNEXT(nd)", "*GET,r,NODE,nd,RF,FZ", "fz=fz+r", "*ENDDO",
        "ALLSEL", "ETABLE,szz,S,Z", "SSUM", "*GET,ssz,SSUM,0,ITEM,SZZ",
        "res(i,1)=t", f"res(i,2)={f(dz)}*t", "res(i,3)=-fz", "res(i,4)=ssz/nel",
        "*ENDDO",
        f"*CFOPEN,{case},crv", "*VWRITE,res(1,1),res(1,2),res(1,3),res(1,4)", "(4E20.11)", "*CFCLOS",
        f"*CFOPEN,{case},met", "*VWRITE,cnv,dcpu,nel,nnd,nset,tlast", "(6E20.11)", "*CFCLOS",
        "FINISH"]


# ---------------------------------------------------------------- run + post
class Runner:
    def __init__(self):
        from ansys.mapdl.core import launch_mapdl
        RUN.mkdir(exist_ok=True)
        t0 = time.perf_counter()
        self.m = launch_mapdl(exec_file=EXE, run_location=str(RUN), jobname="cp",
                              nproc=1, additional_switches="-smp", override=True,
                              start_timeout=180)
        print(f"MAPDL {self.m.version} launched in {time.perf_counter()-t0:.1f} s")

    def solve(self, case, lines):
        (RUN / f"{case}.dat").write_text("\n".join(lines) + "\n")
        for ext in ("crv", "met"):
            (RUN / f"{case}.{ext}").unlink(missing_ok=True)
        self.m.clear()
        t0 = time.perf_counter()
        err = None
        try:
            out = self.m.input(str(RUN / f"{case}.dat"))
        except Exception as e:                       # keep going, record verbatim
            out, err = f"{type(e).__name__}: {e}", f"{type(e).__name__}: {e}"
        wall = time.perf_counter() - t0
        (RUN / f"{case}.out").write_text(str(out), errors="replace")
        crv, met = RUN / f"{case}.crv", RUN / f"{case}.met"
        if not (crv.exists() and met.exists()):
            return dict(case=case, ok=False, error=err or "no result files", wall_s=wall)
        c = np.loadtxt(crv, ndmin=2)
        cnv, dcpu, nel, nnd, nset, tlast = np.loadtxt(met)
        warn = [ln.strip() for ln in str(out).splitlines() if "*** WARNING" in ln]
        return dict(case=case, ok=True, error=err, wall_s=wall, cpu_s=float(dcpu),
                    converged=bool(cnv == 1 and abs(tlast - 1.0) < 1e-6),
                    nel=int(nel), nnodes=int(nnd), substeps=int(nset),
                    strain=(c[:, 1] / L).tolist(), sig_rf=(c[:, 2] / L**2).tolist(),
                    sig_avg=c[:, 3].tolist(), n_warnings=len(warn))

    def close(self):
        self.m.exit()


def yield_offset(strain, sig):
    eps, s = np.array([0.0, *strain]), np.array([0.0, *sig])
    E0 = s[1] / eps[1]
    g = s - E0 * (eps - OFFSET)
    for i in range(1, len(g)):
        if g[i] <= 0 < g[i-1]:
            w = g[i-1] / (g[i-1] - g[i])
            return float(s[i-1] + w * (s[i] - s[i-1])), float(E0)
    return float("nan"), float(E0)


# ---------------------------------------------------------------- stages
def stage_elastic(R):
    rows = []
    for conv, gfac in itertools.product(("bunge", "inverse"), (1, 2)):
        for key in ORIENTS:
            case = f"el_{conv[:3]}_g{gfac}_{list(ORIENTS).index(key)}"
            r = R.solve(case, deck(case, mat_cp(key, conv, gfac), 0.0002, 2))
            Ez = r["sig_rf"][-1] / r["strain"][-1] if r["ok"] else float("nan")
            Eref = E_dir(ORIENTS[key][0])
            err = (Ez - Eref) / Eref * 100
            rows.append(dict(convention=conv, shear_input=f"{gfac}*G", orientation=key,
                             Ez_fea_GPa=Ez/1000, Ez_exact_GPa=Eref/1000, err_pct=err,
                             ok=r["ok"], error=r.get("error")))
            print(f"  elastic {conv:7s} C4={gfac}G {key:13s} E_z={Ez/1000:9.3f}  exact={Eref/1000:9.3f}  err={err:+.3f}%")
    best = {}
    for conv, gfac in itertools.product(("bunge", "inverse"), (1, 2)):
        sel = [r for r in rows if r["convention"] == conv and r["shear_input"] == f"{gfac}*G"]
        best[(conv, gfac)] = max(abs(r["err_pct"]) for r in sel)
    (conv, gfac), worst = min(best.items(), key=lambda kv: kv[1])
    print(f"  -> selected convention={conv}, C4-C6 = {gfac}*G (max |err| {worst:.4f} %)")
    return rows, dict(convention=conv, gfac=gfac, max_err_pct=worst)


def stage_cp(R, conv, gfac):
    rows, curves = [], {}
    for idx, key in enumerate(ORIENTS):
        case = f"cp_{idx}"
        r = R.solve(case, deck(case, mat_cp(key, conv, gfac), EPS_END, NSUB_CP))
        row = dict(model="crystal plasticity (TB,XTAL)", orientation=key, ok=r["ok"],
                   error=r.get("error"))
        m = schmid_max(ORIENTS[key][0])
        if r["ok"]:
            sy, E0 = yield_offset(r["strain"], r["sig_rf"])
            sy_avg, _ = yield_offset(r["strain"], r["sig_avg"])
            exp_ratio = EXPECTED_RATIO[key]
            ratio = sy / TAU_C
            err = (ratio - exp_ratio) / exp_ratio * 100
            row.update(sim_yield_MPa=sy, sim_yield_avgSZ_MPa=sy_avg, expected_MPa=exp_ratio*TAU_C,
                       sim_ratio=ratio, expected_ratio=exp_ratio, schmid_check_ratio=1/m,
                       err_pct=err, pass_=abs(err) <= TOL_PCT, E0_GPa=E0/1000,
                       E_exact_GPa=E_dir(ORIENTS[key][0])/1000,
                       sig_at_1pct_MPa=r["sig_rf"][-1], converged=r["converged"],
                       substeps=r["substeps"], cpu_s=r["cpu_s"], wall_s=r["wall_s"],
                       elements=r["nel"], warnings=r["n_warnings"])
            curves[key] = dict(strain=r["strain"], sigma_rf=r["sig_rf"], sigma_avgSZ=r["sig_avg"])
            print(f"  CP {key:13s} sy={sy:8.2f} MPa  ratio={ratio:.3f} (exp {exp_ratio})  err={err:+.2f}%  "
                  f"conv={r['converged']} nsub={r['substeps']} cpu={r['cpu_s']:.2f}s")
        else:
            print(f"  CP {key}: FAILED {r.get('error')}")
        rows.append(row)
    return rows, curves


def stage_ludwick(R):
    rows, curves = [], {}
    lw = LUDWICK
    expected = lw["sigma0"] + lw["K"] * OFFSET ** lw["N"]
    for idx, key in enumerate(ORIENTS):
        case = f"lw_{idx}"
        r = R.solve(case, deck(case, mat_ludwick(key), EPS_END, NSUB_CP))
        row = dict(model="ANEL + MISO (Ludwick)", orientation=key, ok=r["ok"], error=r.get("error"))
        if r["ok"]:
            sy, E0 = yield_offset(r["strain"], r["sig_rf"])
            err = (sy - expected) / expected * 100
            row.update(sim_yield_MPa=sy, expected_MPa=expected, err_pct=err,
                       pass_=abs(err) <= TOL_PCT, E0_GPa=E0/1000,
                       E_exact_GPa=E_dir(ORIENTS[key][0])/1000,
                       sig_at_1pct_MPa=r["sig_rf"][-1], converged=r["converged"],
                       substeps=r["substeps"], cpu_s=r["cpu_s"], wall_s=r["wall_s"],
                       elements=r["nel"], warnings=r["n_warnings"])
            curves[key] = dict(strain=r["strain"], sigma_rf=r["sig_rf"], sigma_avgSZ=r["sig_avg"])
            print(f"  Ludwick {key:13s} sy={sy:8.2f} MPa (exp {expected:.2f})  err={err:+.2f}%  "
                  f"E0={E0/1000:.2f} GPa conv={r['converged']} cpu={r['cpu_s']:.2f}s")
        else:
            print(f"  Ludwick {key}: FAILED {r.get('error')}")
        rows.append(row)
    return rows, curves


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all", choices=["all", "elastic", "cp", "ludwick"])
    ap.add_argument("--convention", default=None, choices=["bunge", "inverse"])
    ap.add_argument("--gfac", type=int, default=None)
    ap.add_argument("--ndiv", type=int, default=NDIV,
                    help="elements per edge (field is homogeneous; TB,XTAL costs ~0.5 s CPU per element per substep)")
    a = ap.parse_args()
    globals()["NDIV"] = a.ndiv

    out = {"inputs": dict(C11=C11, C12=C12, C44=C44, tau_c=TAU_C, L_um=L, ndiv=a.ndiv,
                          nsub=NSUB_CP, eps_end=EPS_END, flfcc=FLFCC, hard=HARD, ludwick=LUDWICK),
           "schmid_factor": {k: schmid_max(v[0]) for k, v in ORIENTS.items()}}
    prev = HERE / "results.json"
    if prev.exists():
        old = json.loads(prev.read_text())
        out.update({k: v for k, v in old.items() if k not in out})

    R = Runner()
    try:
        if a.stage in ("all", "elastic"):
            print("\n[elastic calibration]")
            out["elastic_rows"], out["calibration"] = stage_elastic(R)
        cal = out.get("calibration", {})
        conv = a.convention or cal.get("convention", "bunge")
        gfac = a.gfac or cal.get("gfac", 1)
        if a.stage in ("all", "cp"):
            print(f"\n[crystal plasticity]  convention={conv}  C4-C6={gfac}*G")
            out["cp_rows"], out["cp_curves"] = stage_cp(R, conv, gfac)
            out["cp_material_commands"] = mat_cp("[1 0.1 0.02]", conv, gfac)
        if a.stage in ("all", "ludwick"):
            print("\n[Ludwick fallback]")
            out["ludwick_rows"], out["ludwick_curves"] = stage_ludwick(R)
            out["ludwick_material_commands"] = mat_ludwick("[1 0.1 0.02]")
    finally:
        R.close()

    prev.write_text(json.dumps(out, indent=1, default=float))
    rows = out.get("cp_rows", []) + out.get("ludwick_rows", [])
    keys = ["model", "orientation", "sim_yield_MPa", "expected_MPa", "err_pct", "pass_",
            "converged", "substeps", "cpu_s", "wall_s", "elements", "E0_GPa", "E_exact_GPa",
            "sig_at_1pct_MPa", "warnings", "error"]
    with open(HERE / "results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwritten: {HERE / 'results.json'}, {HERE / 'results.csv'}")


if __name__ == "__main__":
    main()
