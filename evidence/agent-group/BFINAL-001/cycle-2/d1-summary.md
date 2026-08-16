# BFINAL-001 / cycle-2 / D1 — re-anchor run results

Date: 2026-08-16 (session), run 11:35:36–11:43:36 (ExecutionTime 480.59 s, ClockTime 484 s, exit 0)
Mode: DIAGNOSTIC_ONLY (read-only wrt src/ and case inputs). Sanctioned B4 diagnostic path run once.

## Command

```bash
bash /home/ys/dsH/run_mtohf.sh /home/ys/dsH/b2_case_smoke > evidence/agent-group/BFINAL-001/cycle-2/rerun.log 2>&1
```

Fresh log: `evidence/agent-group/BFINAL-001/cycle-2/rerun.log` (1635 lines).
Reference: `/home/ys/b2_case_smoke/log.b4export` (1634 lines) — only benign diffs, see below.

## Extracted numbers (rerun.log vs reference log.b4export)

| Probe | rerun.log value | reference log.b4export | match |
|---|---|---|---|
| g_w thermal | 2.7025486001e-11 | 2.7025486001e-11 | exact |
| g_w pressure | 7.07793481308e-12 | 7.07793481308e-12 | exact |
| Rx-A (fixed phi, dR_U/dalpha) relL2 | 1.61351152954e-14 | 1.61351152954e-14 | exact |
| Rx-B \|dRp/dalpha\|L2 | 8.14730839032e-10 | 8.14730839032e-10 | exact |
| Stage B4 T1 (momentum U-U, dp=0) relL2 / cos | 0.00852323456072 / 0.999963676577 | same | exact |
| Stage B4 T-cont (continuity row) | 0.000648577427676 | same | exact |
| Stage B4 T2 (dp!=0 continuity) | 7.51176551272e-21 | same | exact |
| Stage B4 T-fvm (frozen-phi fvm matrix part) | 0.999999999866 | same | exact |
| Stage B4 T-dev (deviatoric part) | 5.36114041232e-05 | same | exact |
| Stage B4 T-conv (convection phi-variation) | 0.418119174867 | same | exact |
| GateH1 same-basis dH relL2 / cos | 0.117010638336 / 0.999433960379 | same | exact |

D1 acceptance values reproduced: g_w 2.7e-11/7.08e-12, Rx-A 1.6e-14, Rx-B 8.1e-10, T1 0.0085, T-fvm≈1.0, T-conv 0.418. ✓

## Full-log diff vs reference (only benign / expected-import-state lines)

`diff log.b4export rerun.log` shows ONLY:

1. Header: Date/Time/PID/Case (new run identity; Case = /home/ys/dsH/b2_case_smoke writable copy vs reference /home/ys/b2_case_smoke).
2. `Reduced discrete flow adjoint pressureDrop` block:
   - reference (export run): `GMRES iterations=0, relative residual=1`, `residual=1 > tolerance` , `UcFingerprint[afterSolve pressureDrop]` all zeros, `gsensPressureDrop` solve 0 iterations, `dDP=[-0, -0]`.
   - rerun (import state): `Loaded explicit solution "/home/ys/b2_case_smoke/explicitSol_pressureDrop.mtx" relRes=3.4765473443e-09`, nonzero UcFingerprint (L2=46984.37), `gsensPressureDrop` solved to 4.45992125251e-10 in 11 iterations, `dDP=[-0.0200726738411, 0.0194203165319]`.
   This is the expected B4-import behavior (optProperties `discreteUseExplicitSolution=true`, explicit files read from the read-only reference copy); the reference log.b4export was the export run before the explicit solution existed. All Stage-B4 probe numbers are identical between the two.
3. CSV log path (writable copy vs reference).
4. ExecutionTime (480.59 vs 491.5 s).

## No source/case-input modification

- `git status --short` after the run: identical tracked modifications to pre-state (src/solveDiscreteFlowAdjoint.H +222/-2, src/validateStageBRepeatability.H mode change) — no new tracked changes.
- Case inputs (0/, constant/, system/, constant/optProperties) NOT hand-edited. Solver wrote its own outputs into the writable copy during the sanctioned run (stageB2_*.mtx, explicitRhs_*.mtx, optimization_history.csv/log, time dir `1/`) — this is the sanctioned run's own export behavior, not a hand-edit.
- Reference case /home/ys/b2_case_smoke untouched (its optimization_history.csv/log.dat still dated Aug 15 20:27).

## Verdict

D1 ANCHOR achieved: current HEAD (df695c7 + B4 edits) reproduces log.b4export exactly on all B-final probe numbers. This is an anchor only, not a discriminator; (a) J!=R_w, (b) R_x incomplete, (c) RHS/sign/assembly remain OPEN per plan (steps D2–D4 to follow in subsequent executor steps).
