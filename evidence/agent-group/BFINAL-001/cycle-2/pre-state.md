# BFINAL-001 / cycle-2 / D1 — Pre-state baseline (current-HEAD re-anchor)

Date: 2025-08-16 (session)
Mode: DIAGNOSTIC_ONLY — read-only wrt src/ and case inputs.

## Workspace identity

- `pwd` = `/home/ys/dsH/TO-ANISOTROPIC`
- `git rev-parse --show-toplevel` = `/home/ys/dsH/TO-ANISOTROPIC`  (matches shared workspace)
- `git rev-parse HEAD` = `df695c7373297425e2d76767d2c415a9f7f55b18` (short `df695c7`)

## HEAD state

HEAD = `df695c7` — matches the approved plan's stated HEAD.

## Modified tracked files (uncommitted B4 edits)

| File | Change |
|---|---|
| `src/solveDiscreteFlowAdjoint.H` | content modified: +222 / -2 (B4 probe instrumentation: explicit frozen deviatoric-stress transpose export, Stage B4 Jacobian probes, etc.) |
| `src/validateStageBRepeatability.H` | mode change `100755 => 100644`, 0 content lines changed |

These are exactly the two modified src files expected by the approved plan ("HEAD df695c7 + uncommitted B4 edits (solveDiscreteFlowAdjoint.H, validateStageBRepeatability.H)").

## Full `git status --short` snapshot (pre-existing user changes preserved; nothing modified by this cycle)

- `D` CODEX_LOCAL_VALIDATION.md, LOCAL_FROZEN_GRADIENT_VALIDATION_REPORT.md, README.md (pre-existing deletions — preserved)
- `M` src/solveDiscreteFlowAdjoint.H, src/validateStageBRepeatability.H (above)
- Untracked: build_minifd.sh, build_oracle.sh, check_*.sh, crosscheck_pd.py, evidence/, peek_log.sh, probe_*.sh, project_pd.py, ps_check.sh, readme.md, run_b2diag.sh, run_b4_export.sh, run_b4_import.sh, show_cfg.sh, solve_oracle_pd.py, src/validation/, src/validation_stage_*.log, write_C.sh

No `git reset --hard` / `git clean -fd` used. User's dirty worktree fully preserved.

## File hashes (pre-run, current worktree)

- `src/solveDiscreteFlowAdjoint.H`  md5 = `881bc664a8acf5d97814eb46dfb08dce`
- `src/validateStageBRepeatability.H` md5 = `a348117873f48c85555eae8af96b98a8`

## Writable case copy state (read-only wrt inputs this cycle)

Case: `/home/ys/dsH/b2_case_smoke` (contains 0/ constant/ system/; optProperties in B4-import state):

- `stageB4JacobianProbe true`
- `discreteExportOnly true`
- `discreteUseExplicitSolution true`
- `adjointMode discrete`
- `mmaUpdateEnabled false`
- `frozenGradientValidated false`
- Explicit matrix/rhs/sol files point at the read-only reference copy `/home/ys/b2_case_smoke/` (explicitJT.mtx, explicitRhs_*, explicitSol_*) — not modified.

## Reference logs (read-only, for comparison)

- `/home/ys/b2_case_smoke/log.b4export` (the B4-export run this cycle re-anchors)
- `/home/ys/b2_case_smoke/log.b4import`
- `/home/ys/b2_case_smoke/log.minifd`

## Purpose of this re-anchor (D1 acceptance)

Run the sanctioned B4 diagnostic path once (bash /home/ys/dsH/run_mtohf.sh /home/ys/dsH/b2_case_smoke), capture fresh log to
`evidence/agent-group/BFINAL-001/cycle-2/rerun.log`, and extract g_w thermal/pressure, Rx-A, Rx-B, Stage B4 T1/T-cont/T2/T-fvm/T-dev/T-conv and GateH1 same-basis dH.
This is an ANCHOR ONLY, not a discriminator: it must reproduce `log.b4export` (g_w 2.7e-11/7.08e-12, Rx-A 1.6e-14, Rx-B 8.1e-10, T1 0.0085, T-fvm≈1.0, T-conv 0.418).
