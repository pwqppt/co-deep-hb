"""Step 11 orientation test: model delta eps'_zz(30 -> 400 C after 2 h) vs Ayoub Fig. 2(c), per elastic pad.

Pre-registered rule (docs/stage1_step11_orientation_test.md):
    r_p = delta_model_p / delta_meas_p for p in {4, 5, 7, 9, 10}; r_mean = mean(r_p)
    tol_p = max(0.10, err_p / delta_meas_p)
    PASS  iff  |r_p / r_mean - 1| <= tol_p  for all five pads.

Usage (repo root):  .venv\\Scripts\\python.exe src\\tools\\step11_compare.py
Writes results/diagnostics/step11_orientation_test.csv and prints the verdict.
"""
import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
RUNS = {10: "R26", 4: "R31", 5: "R32", 7: "R33", 9: "R34"}
LS_START, LS_HOLD_END, LS_HEAT_END = 0, 8, 7      # 30 C (after load step 0), 400 C after 2 h hold, 400 C end of heating


def main():
    fig2 = json.loads((REPO / "inputs/ayoub2022_fig2ac_digitized.json").read_text(encoding="utf-8"))
    meas = {d["pad"]: d for d in fig2["fig2c_delta_30_400h"]}
    fig5 = json.loads((REPO / "inputs/ayoub2022_fig5_digitized.json").read_text(encoding="utf-8"))
    f5 = {p["pad"]: p for p in fig5["points"]}
    a, b = fig5["paper_fit"]["a"], fig5["paper_fit"]["b"]
    rows = []
    for pad, rid in RUNS.items():
        f = REPO / "results" / rid / "record.json"
        if not f.exists():
            print(f"{rid} (pad {pad}): record.json missing - run it first"); return 1
        rec = json.loads(f.read_text(encoding="utf-8"))
        by = {r["ls"]: r for r in rec["rows"]}
        d_model = (by[LS_HOLD_END]["eps_dev_zz"] - by[LS_START]["eps_dev_zz"]) * 1e3
        d_heat = (by[LS_HEAT_END]["eps_dev_zz"] - by[LS_START]["eps_dev_zz"]) * 1e3
        m = meas[pad]
        Ez = rec["deck"]["mat"]["Ez_GPa"]
        rows.append(dict(pad=pad, run=rid, Ez=Ez, E_table=m["E_table"], d_model=d_model, d_model_heating=d_heat,
                         d_meas=m["value"], err=m["err"], ratio=d_model / m["value"], tol=max(0.10, m["err"] / m["value"]),
                         u_model=by[LS_HEAT_END]["u_center_nm"], u_fig5_point=f5.get(pad, {}).get("u_nm", np.nan),
                         u_fig5_fit_Etable=a + b * m["E_table"]))
    r_mean = np.mean([r["ratio"] for r in rows])
    for r in rows:
        r["dev_from_mean"] = r["ratio"] / r_mean - 1
        r["pass"] = abs(r["dev_from_mean"]) <= r["tol"]
    verdict = all(r["pass"] for r in rows)
    Ez = np.array([r["Ez"] for r in rows]); u = np.array([r["u_model"] for r in rows])
    slope = np.polyfit(Ez, u, 1)[0]
    out = REPO / "results/diagnostics/step11_orientation_test.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    keys = list(rows[0].keys())
    with out.open("w", encoding="utf-8") as fh:
        fh.write(",".join(keys) + "\n")
        for r in sorted(rows, key=lambda r: r["Ez"]):
            fh.write(",".join(f"{r[k]:.5g}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
    print(f"{'pad':>3} {'run':>4} {'Ez':>7} {'d_model':>8} {'d_meas':>7} {'ratio':>6} {'dev':>7} {'tol':>6}  pass   u_model  u_fig5")
    for r in sorted(rows, key=lambda r: r["Ez"]):
        print(f"{r['pad']:>3} {r['run']:>4} {r['Ez']:7.2f} {r['d_model']:8.3f} {r['d_meas']:7.3f} {r['ratio']:6.3f} "
              f"{r['dev_from_mean']:+7.3f} {r['tol']:6.3f}  {str(r['pass']):5}  {r['u_model']:7.3f}  {r['u_fig5_point']:6.3f}")
    print(f"mean ratio = {r_mean:.3f};  model u-vs-E_z slope (elastic pads) = {slope:.4f} nm/GPa (Fig. 5 fit: {b})")
    print("STEP 2 VERDICT:", "PASS (single orientation-independent factor) -> Option B justified" if verdict
          else "FAIL -> Option B NOT justified, stop")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
