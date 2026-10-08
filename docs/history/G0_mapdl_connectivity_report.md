# ANSYS MAPDL connectivity check — Student 2026 R1 (v261), Windows

**Date:** 2026-10-08 · **Verdict: no connection method works from inside the Claude Science sandbox.**
The solver binary is reachable and runs, but it cannot obtain a license there. The test case
was therefore **not solved**; no FEA-vs-hand-calculation comparison is reported below.

## Installation found

| Item | Value |
|---|---|
| Install root | `C:\Program Files\ANSYS Inc\ANSYS Student\v261` |
| MAPDL launcher | `...\v261\ansys\bin\winx64\ANSYS261.exe` |
| Version banner | `Ansys MAPDL 2026 R1  Build 26.1  UP20260202  WINDOWS x64` |
| Licensing client | `...\v261\licensingclient\winx64\{ansyscl,ansysli_util,lmutil}.exe` |

The sandbox sees the **real** host filesystem for `C:\Program Files`, so discovery and file
staging work normally. No host-folder grant was required.

## What was attempted

| # | Method | Outcome | Wall clock |
|---|---|---|---|
| 1 | MAPDL batch, default launcher | Failed before solve — Intel MPI bootstrap | 0.3 s |
| 2 | MAPDL batch, `-smp` | Solver ran; **license checkout failed** | 257 s |
| 3 | PyMAPDL `launch_mapdl` (`ansys-mapdl-core`) | Process spawned; gRPC connect refused | 120 s |
| 4 | Ansys MCP connector | **Not attached to this session** | — |

## Root causes (three independent sandbox limits)

**1. Intel MPI bootstrap cannot spawn children.** The default launcher routes even `-np 1`
through `mpiexec`, which fails with `HYD_spawn ... unable to create stdout pipe`.
*Workaround found:* add `-smp`. This one is solved — keep `-smp` in all MAPDL calls.

**2. The licensing client cannot load (hard blocker).** Running `ansyscl.exe` or
`ansysli_util.exe` directly fails immediately with NTSTATUS **623
`STATUS_ILLEGAL_DLL_RELOCATION`** ("a system DLL was relocated in memory"). MAPDL therefore
reports, 15 times over ~4 minutes:

```
ANSYS LICENSE MANAGER ERROR:
ANSYSLI exited or could not read server port ansyscl.<HOST>.16488.96047
Unable to access ...\Packages\operon.kernel.<id>\AC\Temp\.ansys\ansyscl.<HOST>...log
*** ERROR - ANSYS license not available.
```

The `...\Packages\<pkg>\AC\Temp\...` path shows the sandbox is a Windows **AppContainer**.
The licensing helper cannot initialise inside it. This is not a license-file or firewall
problem — `obtaining a license: 0.0 seconds` and every solver phase reads `0.0 seconds`,
i.e. the model never entered `/PREP7`. Not fixable from inside the sandbox.

**3. gRPC to loopback is blocked (blocks PyMAPDL specifically).** PyMAPDL talks to MAPDL over
`127.0.0.1:50052`. The sandbox exports `HTTP_PROXY=http://127.0.0.1:52194`, which answers the
gRPC CONNECT with **403**. Setting `GRPC_ENABLE_HTTP_PROXY=0` does *not* help — a direct
gRPC channel to an in-sandbox loopback listener also times out, because connections to
private/reserved IPs are refused below the proxy layer. (A plain Python TCP socket to
in-sandbox loopback *does* succeed, so the block is specific to the routed path, not to
loopback in general.)

Consequence: **PyMAPDL cannot be used from the sandbox at all** — not by launching MAPDL
itself, and not by attaching to an MAPDL gRPC server you start natively.

## Permissions needed

Nothing grantable would fix this. No host-folder grant is needed (the install is already
visible), and `request_network_access` explicitly **cannot** authorise a private/reserved IP
target, so the loopback block is not a setting anyone can toggle.

## Recommended path: the Ansys MCP connector

An MCP server runs as a **native host process, outside the AppContainer** — so it has a
working licensing client and no proxy in front of loopback. That makes it the only viable
route, and it sidesteps all three limits above. It is currently **not attached**: the
connectors in this session are all biology/chemistry data sources.

To enable it: **Customize → Connectors**, add the Ansys MCP server, then ask for this test
again. If the connector exposes a generic "run APDL input file" tool, `cu_bar_thermal_test.dat`
runs as-is.

Interim fallback: run `run_cu_bar_test.py` yourself in a normal (non-sandbox) terminal. It
tries PyMAPDL, falls back to batch, and prints PASS/FAIL per quantity.

## Analytic reference for the test case

µMKS consistent units — length µm, mass kg, time s, temperature °C ⇒ force µN, stress **MPa**.
Cu bar: `L0 = 1000 µm`, `A = 1 µm²`, `E = 110 GPa = 110000 MPa`, `α = 16.5e-6 /°C` (secant,
ref 25 °C), `T: 25 → 400 °C`, so `ΔT = 375 °C`.

| Quantity | Closed form | Value |
|---|---|---|
| Thermal strain | `α·ΔT` | `6.187500e-03` |
| Free expansion `UX` | `α·ΔT·L0` | **`6.187500 µm`** |
| Constrained stress `Sxx` | `−E·α·ΔT` | **`−680.625000 MPa`** |
| Constrained reaction | `σ·A` | `−680.625 µN` |

The deck solves both load steps: step 1 free (checks CTE and length units), step 2 fully
constrained (checks that `E` in MPa is consistent with µm/µN — the unit pairing that matters
for your pad-recess and contact-margin models). The runner asserts both to 1 %.

## Reusable gotchas for this project

- Always pass `-smp` to MAPDL on this machine; the MPI launcher is unusable from the sandbox.
- Student limit is 128k nodes/elements — not approached here (2 nodes, 1 element).
- `MP,REFT` **and** `TREF` are both set to 25 °C. If they disagree, MAPDL uses `MP,REFT` for
  the material and silently gives a strain offset — a classic source of wrong recess values.
