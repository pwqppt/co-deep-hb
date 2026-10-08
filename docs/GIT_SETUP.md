# Git setup and multi-PC routine

Remote: `https://github.com/pwqppt/co-deep-hb.git` (private). Working copy on every PC: `$HOME\co-deep-hb`.

## Step 1 (PC 1) — check git, then the restructure regression check, **before** any commit

Runs R01 and R26 from the new layout in the staging folder and compares against the pre-restructure
numbers. About 3.5 min.
```
git --version; git config --global user.name; git config --global user.email; cd "C:\Users\USER\.claude-science\orgs\f50c262b-2754-4e2c-b4cd-c0deca628d34\workspaces\0b7dce47-624a-4a52-8bdc-8d464908c370\cu-hybrid-bonding"; & "C:\Users\USER\.claude-science\conda\envs\ansys\python.exe" src\stage1\stage1_pipeline.py --case R01,R26 --force; & "C:\Users\USER\.claude-science\conda\envs\ansys\python.exe" src\tools\check_run.py R01 R26
```
- **Expected first:** `git version 2.x.x.windows.x`, then your name and e-mail. If either is empty, set
  them once: `git config --global user.name "Your Name"; git config --global user.email "you@example.com"`.
- **Expected last line:** `ALL MATCH`. Commit only if you see it.

## Step 2 (PC 1) — init, first commit, push (only after ALL MATCH)

```
cd "C:\Users\USER\.claude-science\orgs\f50c262b-2754-4e2c-b4cd-c0deca628d34\workspaces\0b7dce47-624a-4a52-8bdc-8d464908c370\cu-hybrid-bonding"; git init -b main; git add -A; git status --short | Measure-Object -Line; git commit -m "Stage 1 repo: pipeline + mesh (incl. sliding-sidewall fix), inputs, docs, results R01 R24-R26, validation G1 G3"; git remote add origin https://github.com/pwqppt/co-deep-hb.git; git push -u origin main
```
- **Expected:** a line count of about 100 files.
  - The first push opens a browser window ("Sign in to GitHub", from Git Credential Manager).
  - It ends with `branch 'main' set up to track 'origin/main'`.
- **Check:** no `.rst` or `.pdf` file was added:
  `git ls-files | Select-String "\.rst$|\.pdf$"` should print nothing.

## Step 3 (PC 1) — permanent working copy (stop using the chat workspace)

```
cd $HOME; git clone https://github.com/pwqppt/co-deep-hb.git co-deep-hb; cd co-deep-hb; & "C:\Users\USER\.claude-science\conda\envs\ansys\python.exe" -m venv .venv; .\.venv\Scripts\python.exe -m pip install -r requirements.txt; Copy-Item configs\machine.example.json configs\machine.json; .\.venv\Scripts\python.exe src\stage1\stage1_pipeline.py --list
```
- **Expected:** the run table, with R01, R24, R25 and R26 marked `DONE`.
- From now on, run everything from `$HOME\co-deep-hb`. The staging folder in the chat workspace can be
  deleted later. Its `R01.rst` and `R25.rst` are the only binary result copies, so copy them first if
  you want them.

## Step 4 (each other PC) — clone, set up, run one case and compare with the committed result

ANSYS Student 2026 R1 must be installed and licensed there, and the teammate needs collaborator access
to the repository (GitHub → Settings → Collaborators). If `py` is missing, install Python first:
`winget install -e --id Python.Python.3.12`.
```
cd $HOME; git clone https://github.com/pwqppt/co-deep-hb.git co-deep-hb; cd co-deep-hb; py -3 -m venv .venv; .\.venv\Scripts\python.exe -m pip install -r requirements.txt; Copy-Item configs\machine.example.json configs\machine.json; .\.venv\Scripts\python.exe src\stage1\stage1_pipeline.py --case R01 --force; .\.venv\Scripts\python.exe src\tools\check_run.py R01 --ref git; git restore results
```
- **Expected:** about 2 min of solve, then `R01: MATCH` and `ALL MATCH`.
  - Tiny last-digit differences across PCs would show as `MISMATCH` with |d| ≈ 1e-12; report them.
  - `git restore results` discards the verification output, so nothing is committed from this check.
- If ANSYS is installed elsewhere, edit `mapdl_exe` in `configs\machine.json` before the run.

## Routine

- **Start of a session:** `cd "$HOME\co-deep-hb"; git pull`, then follow `docs/HANDOFF.md`.
- **End of a session** (gate passed or "wrap up"): update `docs/HANDOFF.md` and `docs/CHANGELOG.md`, then
  ```
  cd "$HOME\co-deep-hb"; git add -A; git status --short; git commit -m "<summary>"; git push
  ```
- **Model freeze** (e.g. end of Stage 1):
  ```
  cd "$HOME\co-deep-hb"; git tag -a v1.0 -m "Stage 1 model freeze"; git push origin v1.0
  ```
- Only one PC should write a given `results/<run>/` at a time. Pull before running and push right
  after, so two PCs never commit the same run.
