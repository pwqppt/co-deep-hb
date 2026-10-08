# co-deep-hb — Cu–Cu hybrid bonding: orientation and recess effects on pad contact

Computational study: when 300 nm single-crystal Cu/SiO2 pads vary in crystal orientation and CMP recess,
does an anneal condition chosen with average Cu properties still guarantee Cu–Cu contact?
Thermo-mechanical FEA in ANSYS Mechanical APDL (Student 2026 R1), then GP surrogates, Monte Carlo and
Bayesian optimisation.

**Design decisions live in the main planning document** (link in `docs/HANDOFF.md`), not in this repo.
`docs/stage1_model_spec.md` is the implementation spec derived from it.

## Layout

| Path | Content |
|---|---|
| `src/stage1/` | Stage 1 pipeline (`stage1_pipeline.py`: decks, runner, collector) and mesh (`s1_mesh.py`) |
| `src/validation/` | G1 Cu bar (`g1_cu_bar_test.py` + deck), G3 anisotropy (`g3_anisotropic_cu_check.py`), ANEL+MISO check (`cp_check.py`) |
| `src/tools/` | `check_run.py` — compare runs against a reference or the committed results |
| `inputs/` | digitized Ayoub 2022 data (ε′zz, Fig. 5), Cu constants vs T (Chang & Himmel), CTE template |
| `configs/` | `machine.example.json` → copy to `machine.json` (per PC, not in git) |
| `results/<run>/` | per-run text outputs (`.out`, `.mntr`, `hist_*.txt`, `meta.txt`, `deck_meta.json`, `record.json`); `results/results.csv` |
| `docs/` | `HANDOFF.md` (run status), `CHANGELOG.md`, `RUNLIST.md`, `stage1_model_spec.md`, `REFERENCES.md`, `GIT_SETUP.md`, `history/` |

Not in git: MAPDL `.rst` and scratch files, generated decks (`results/**/*.dat`), paper PDFs, `.venv/`, `configs/machine.json`.

## Setup on a PC (once)

```
cd $HOME; git clone https://github.com/pwqppt/co-deep-hb.git; cd co-deep-hb
py -3 -m venv .venv; .\.venv\Scripts\python.exe -m pip install -r requirements.txt
Copy-Item configs\machine.example.json configs\machine.json
```
Edit `configs\machine.json` if ANSYS is installed elsewhere. ANSYS must be installed and licensed on that PC.

## Running

From the repo root: `.\.venv\Scripts\python.exe src\stage1\stage1_pipeline.py --case R01` (see `docs/RUNLIST.md`).
ANSYS cannot run inside the Claude Science app; every run is executed in a PowerShell on a PC with ANSYS.
