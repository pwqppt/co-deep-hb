#!/usr/bin/env python
"""
Compare finished Stage 1 runs against a reference, to confirm that a restructure, a
different PC or a different checkout reproduces the same numbers.

    python src/tools/check_run.py R01 R26                 # vs docs/history/regression_reference_pre_restructure.json
    python src/tools/check_run.py R01 --ref git           # vs the committed results/<RUN>/record.json (HEAD)

Compared at the end of heating (400 C): protrusion u_center (nm), eps'_zz, Cu volume-averaged
sigma_z, pad lengthening, node/element/CE counts and the deck physics hash.
MATCH means every quantity agrees to within --rtol (default 1e-9 relative, 1e-12 absolute floor)
and the counts and hash are identical.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BASELINE = REPO / "docs" / "history" / "regression_reference_pre_restructure.json"
FLOATS = ("u_center_nm", "eps_dev_zz", "s_z", "padlen")
EXACT = ("n_nodes", "n_elems", "n_ce")


def summarize(rec):
    h = next(r for r in rec["rows"] if r["tag"] == "h400")
    return dict(u_center_nm=h["u_center_nm"], eps_dev_zz=h["eps_dev_zz"], s_z=h["s_z"],
                padlen=h["u_center_rel_padbottom_nm"], n_nodes=rec["n_nodes"], n_elems=rec["n_elems"],
                n_ce=rec["deck"]["n_ce"], sha=rec["deck"]["physics_sha256"])


def reference(rid, mode):
    if mode == "baseline":
        return json.loads(BASELINE.read_text())[rid]
    blob = subprocess.run(["git", "show", f"HEAD:results/{rid}/record.json"], cwd=REPO,
                          capture_output=True, text=True, check=True).stdout
    return summarize(json.loads(blob))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--ref", choices=["baseline", "git"], default="baseline")
    ap.add_argument("--rtol", type=float, default=1e-9)
    a = ap.parse_args(argv)
    all_ok = True
    for rid in a.runs:
        cur = summarize(json.loads((REPO / "results" / rid / "record.json").read_text()))
        ref = reference(rid, a.ref)
        ok = True
        parts = []
        for k in FLOATS:
            d = abs(cur[k] - ref[k])
            good = d <= max(a.rtol * abs(ref[k]), 1e-12)
            ok &= good
            parts.append(f"{k}={cur[k]:.9g} (ref {ref[k]:.9g}, |d|={d:.2e}{'' if good else ' !'})")
        for k in EXACT:
            good = cur[k] == ref[k]; ok &= good
            if not good:
                parts.append(f"{k}={cur[k]} (ref {ref[k]}) !")
        good = cur["sha"] == ref["sha"]; ok &= good
        parts.append(f"deck hash {'same' if good else 'DIFFERENT'}")
        print(f"{rid}: {'MATCH' if ok else 'MISMATCH'}\n    " + "\n    ".join(parts))
        all_ok &= ok
    print("ALL MATCH" if all_ok else "SOME RUNS DO NOT MATCH")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
