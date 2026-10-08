"""
Cu bar 1-D thermal expansion benchmark -- connection + unit-consistency check
for ANSYS Mechanical APDL Student 2026 R1 (v261) on Windows.

Tries PyMAPDL first, then falls back to an MAPDL batch solve, and compares
both against the closed-form answer.

Run natively (NOT inside the Claude Science sandbox -- the ANSYS licensing
client cannot load there):

    python run_cu_bar_test.py

Analytic reference, uMKS units (um, kg, s, degC -> uN, MPa):
    dL  = alpha*dT*L0 =  6.187500 um
    Sxx = -E*alpha*dT = -680.625000 MPa
"""

from __future__ import annotations
import os
import re
import subprocess
import sys
import time
from pathlib import Path

# ----------------------------------------------------------------- constants
AWP_ROOT = Path(os.environ.get("AWP_ROOT261",
                               r"C:\Program Files\ANSYS Inc\ANSYS Student\v261"))
MAPDL_EXE = AWP_ROOT / "ansys" / "bin" / "winx64" / "ANSYS261.exe"
DECK = Path(__file__).with_name("cu_bar_thermal_test.dat")

L0, ALPHA, E_MPA, TREF, TEND = 1000.0, 16.5e-6, 110000.0, 25.0, 400.0
DT = TEND - TREF
UX_REF = ALPHA * DT * L0        #   6.1875 um
SIG_REF = -E_MPA * ALPHA * DT   # -680.625 MPa
TOL_PCT = 1.0


def _verdict(name, got, ref):
    err = abs(got - ref) / abs(ref) * 100.0
    ok = err <= TOL_PCT
    print(f"  {name:28s} FEA={got: .6f}   exact={ref: .6f}   "
          f"err={err:.4f}%   {'PASS' if ok else 'FAIL'}")
    return ok


def report(ux, sig, method, solve_s):
    print(f"\n--- {method} results (tolerance {TOL_PCT}%) ---")
    ok = _verdict("free expansion UX [um]", ux, UX_REF)
    ok &= _verdict("constrained Sxx [MPa]", sig, SIG_REF)
    print(f"  solve wall-clock: {solve_s:.2f} s")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    return ok


# ------------------------------------------------------------------ PyMAPDL
def try_pymapdl():
    print("\n[1] PyMAPDL ...")
    try:
        from ansys.mapdl.core import launch_mapdl
    except ImportError as exc:
        print(f"    ansys-mapdl-core not importable: {exc}")
        return None

    mapdl = None
    try:
        t0 = time.perf_counter()
        mapdl = launch_mapdl(
            exec_file=str(MAPDL_EXE),
            run_location=str(Path.cwd() / "pymapdl_run"),
            jobname="cubar",
            nproc=1,
            additional_switches="-smp",   # avoid Intel MPI bootstrap
            override=True,
            start_timeout=120,
        )
        print(f"    connected: {mapdl.version} (launch {time.perf_counter()-t0:.1f} s)")
        t0 = time.perf_counter()
        mapdl.input(str(DECK))
        solve_s = time.perf_counter() - t0
        ux = float(mapdl.parameters["UX2"])
        sig = float(mapdl.parameters["SIG"])
        return report(ux, sig, "PyMAPDL", solve_s)
    except Exception as exc:
        print(f"    FAILED: {type(exc).__name__}: {exc}")
        return None
    finally:
        if mapdl is not None:
            try:
                mapdl.exit()
            except Exception:
                pass


# -------------------------------------------------------------------- batch
def try_batch():
    print("\n[2] MAPDL batch ...")
    if not MAPDL_EXE.exists():
        print(f"    solver not found: {MAPDL_EXE}")
        return None

    wd = Path.cwd() / "batch_run"
    wd.mkdir(exist_ok=True)
    inp = wd / "cu_bar.dat"
    inp.write_text(DECK.read_text() + "\n/EXIT,NOSAVE\n")
    for stale in ("results.txt", "cubar.out"):
        (wd / stale).unlink(missing_ok=True)

    cmd = [str(MAPDL_EXE), "-b", "-smp", "-np", "1",
           "-dir", str(wd), "-j", "cubar",
           "-i", str(inp), "-o", str(wd / "cubar.out")]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, errors="replace")
    solve_s = time.perf_counter() - t0

    res = wd / "results.txt"
    if not res.exists():
        out = (wd / "cubar.out")
        tail = out.read_text(errors="replace")[-1500:] if out.exists() else proc.stderr
        print(f"    FAILED (rc={proc.returncode}); last output:\n{tail}")
        return None

    m = re.search(r"RESULT\s+(\S+)\s+(\S+)", res.read_text())
    if not m:
        print(f"    could not parse: {res.read_text()[:300]}")
        return None
    return report(float(m.group(1)), float(m.group(2)), "MAPDL batch", solve_s)


if __name__ == "__main__":
    print(f"solver : {MAPDL_EXE}  (exists={MAPDL_EXE.exists()})")
    print(f"deck   : {DECK}")
    results = [r for r in (try_pymapdl(), try_batch()) if r is not None]
    if not results:
        print("\nNo connection method succeeded.")
        sys.exit(2)
    sys.exit(0 if all(results) else 1)
