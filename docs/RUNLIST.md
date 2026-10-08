# Stage 1 run list (Ayoub 2022 reproduction)

Run status is in `docs/HANDOFF.md`; this file lists the runs and their commands.

Run these in **your** PowerShell on a PC that has the repo cloned at `$HOME\co-deep-hb`, a `.venv`
(see `README.md`) and ANSYS. Every command starts with the `cd`, so it can be pasted as is. If the Student
licence refuses 4 cores, append `--np 2` (or set `nproc` in `configs\machine.json`).

## Status at hand-off

| Step | Runs | Ready now | Blocked on |
|---|---|---|---|
| 1 | R01, R02 — pad 10 elastic, in-plane 0° / 45° | yes | — |
| 2 | R03, R04 — pad 3 elastic, 0° / 45° | yes | — |
| 3 | R05 — pad 3 plastic ② (default) | **yes** (G(T) from Chang & Himmel, Table I) | — |
| 4 | R06 — pad 3 plastic ① · R07 — plastic ③ | yes | — |
| 5 | R08 — pad 10 plastic (= elastic by rule) | yes; **no solve**: deck identical to R01 | — |
| 6 | R09 CTE(T) · R10 sliding · R11 no cap · R12 circular · R13 plastic 45° · R14 Si-slab check | R10–R14 yes | R09: NIST α(T) table (`inputs/cu_cte_vs_T.csv`) |
| 7 | R15 (L0), R16 (L1), R05 (L2), R17 (L3) — G2 mesh study | yes | — |
| 8 | R18–R23 — optional pads 1, 2, 6 (elastic + plastic ②) | yes | — |
| 9 | R01 re-run, R24 no cap, R25 sliding, R26 no cap + sliding — diagnostics for the R01 finding | yes | — |

`--step N` runs every run of that step and skips blocked ones with a message. Runs already done
are skipped unless `--force` is given.

## Commands

**Step 1 — pad 10 elastic (R01, R02).** About 5 min total (2–10).
```
cd "$HOME\co-deep-hb"; .\.venv\Scripts\python.exe src\stage1\stage1_pipeline.py --step 1
```

**Step 2 — pad 3 elastic (R03, R04).** About 5 min (2–10).
```
cd "$HOME\co-deep-hb"; .\.venv\Scripts\python.exe src\stage1\stage1_pipeline.py --step 2
```

**Step 3 — pad 3 plastic ② (R05).** About 20 min (10–45).
```
cd "$HOME\co-deep-hb"; .\.venv\Scripts\python.exe src\stage1\stage1_pipeline.py --step 3
```

**Step 4 — pad 3 plastic ① and ③ (R06, R07).** About 20 min per run. The first plastic run to
finish is the **first test of KEYOPT(2)=3 with ANEL + MISO** (open issue O8).
```
cd "$HOME\co-deep-hb"; .\.venv\Scripts\python.exe src\stage1\stage1_pipeline.py --step 4
```

**Step 5 — pad 10 plastic confirmation (R08).** Seconds, no solve. Writes `results\R08\R08_check.json`.
```
cd "$HOME\co-deep-hb"; .\.venv\Scripts\python.exe src\stage1\stage1_pipeline.py --step 5
```

**Step 6 — variants (R09–R14).** About 20–27 min per run; R09 is skipped until the NIST table exists.
```
cd "$HOME\co-deep-hb"; .\.venv\Scripts\python.exe src\stage1\stage1_pipeline.py --step 6
```

**Step 7 — G2 mesh study (R15, R16, R17; L2 = R05).** About 7 + 13 + 37 min.
```
cd "$HOME\co-deep-hb"; .\.venv\Scripts\python.exe src\stage1\stage1_pipeline.py --step 7
```

**Step 8 — optional pads 1/2/6 (R18–R23).** Elastic about 2–3 min each, plastic about 20 min each.
```
cd "$HOME\co-deep-hb"; .\.venv\Scripts\python.exe src\stage1\stage1_pipeline.py --step 8
```

**Step 9 — R01 diagnostics (R01 re-run + R24, R25, R26).** Elastic, pad 10. About 80 s each (R01 measured
78.6 s), so ≈ 5–6 min. The R01 re-run only adds the vertical force-split output; its physics is unchanged.
```
cd "$HOME\co-deep-hb"; .\.venv\Scripts\python.exe src\stage1\stage1_pipeline.py --case R01,R24,R25,R26 --force
```

**Any single run** (example R03):
```
cd "$HOME\co-deep-hb"; .\.venv\Scripts\python.exe src\stage1\stage1_pipeline.py --case R03
```

**Status table at any time:** the same command with `--list`.

## Expected output files (per run, in `results\<RUN>\`)

| File | Content |
|---|---|
| `<RUN>.dat` | the exact deck solved (≈ 4.4 MB, ≈ 101k lines at L2) |
| `<RUN>.out` | MAPDL log: check this first if a run fails |
| `hist_strain.txt` | per load step: Cu volume-averaged EPEL (6), ε<sub>p</sub> mean/max |
| `hist_stress.txt` | per load step: Cu volume-averaged stress (6), max von Mises, sidewall σ<sub>n</sub> max/min, τ max |
| `hist_disp.txt` | per load step: UZ at the probe nodes, CPU time, convergence flag |
| `meta.txt` | element/node counts and frame check (Cu element ESYS, element 1 material) |
| `deck_meta.json` | inputs: E<sub>z</sub>, T<sub>sf</sub> and source, plasticity rule, probe node IDs |
| `record.json` | everything above, parsed, with exit code and wall time |

Aggregated after every invocation: `results\results.csv` (one row per run × load step) and
`results\results.json`.

**Runtime caveat.** The minutes above are estimates (±2×) from model size and substep count, not
measurements. R01's actual wall time recalibrates them.

## Unblocking the remaining runs (data only you can fetch)

1. **G(T) — done.** `inputs/cu_elastic_constants_vs_T.csv` is transcribed from Chang & Himmel
   (UCRL-16697), Table I, report p. 14.
2. **NIST α(T) for R09.** Save the NIST SRM 736 certificate, or Hahn, *J. Appl. Phys.* 41, 5096
   (1970), into `inputs\`. I will build `cu_cte_vs_T.csv` (instantaneous CTE, 30–400 °C)
   from it.
