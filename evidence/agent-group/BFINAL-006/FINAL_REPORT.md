# BFINAL-006 FINAL_REPORT — Pressure-Drop Tangent Oracle (cycle-2, DIAGNOSTIC_ONLY)

- Stage: B-final · Mode: **DIAGNOSTIC_ONLY** · Hypothesis: `BFINAL-006-TANGENT-CLOSURE`
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace) · HEAD `ca8a772b8339a36e7906ea32a7125061c47aa280` · branch `agent/dsH-stage-b-validation`
- **Decision: FAIL_PLAN_OR_HYPOTHESIS** (research-level: tangent closure FAILED at T1 — the
  tangent linear system cannot be closed; see §7)

> The corrected pressure pin (pRefRow=100800) is RIGHT, but it is NOT sufficient: the
> closed BFINAL-003 operator (transpose of `explicitJT.mtx`, re-pinned) is **numerically
> singular** with a non-constant, outlet-concentrated pressure null mode. The approved
> direct-solve T1 stalls at true relative residual ~3.6e-4 (target ≤1e-9). This refutes the
> hypothesis that the cycle-1 failure was *purely* the pin identification with no defect in
> J; it promotes handoff hypothesis #2 ("implemented J is self-consistent with its transpose
> but not the exact nonsingular derivative of the frozen primal residual").

---

## Provenance

| item | value |
|---|---|
| pwd | `/home/ys/dsH/TO-ANISOTROPIC` |
| git root | `/home/ys/dsH/TO-ANISOTROPIC` (== workspace) |
| HEAD | `ca8a772b8339a36e7906ea32a7125061c47aa280` (unchanged) |
| branch | `agent/dsH-stage-b-validation` |
| case | `/home/ys/dsH/b2_case_smoke` |
| source changes this round | **NONE** — `git status --porcelain` (38 lines) byte-identical to `cycle-2/pre-state/git_status_porcelain.txt` |
| locked-head sha256 | unchanged (`solveDiscreteFlowAdjoint.H` f0c81497…, `solveDiscreteFlowAdjointProduction.H` 8c4901bf…, `sensitivity.H` f0ed0726…, `rxPressureRowTranspose.H` fd648ae0…, `stageB6RxDesignOracle.H` 09832dd6…) |
| artifacts (sha256) | `explicitJT.mtx` 3dfcc9ab…, `stageB6_rxc_analytic.mtx` d30e1317…, `stageB6_dirs.mtx` 75722615…, `explicitRhs_pressureDrop.mtx` 2d7cf846… (all unchanged from cycle-1) |

Only evidence files were written this round (`evidence/agent-group/BFINAL-006/cycle-2/…`
and this report). No `stageB7TangentOracle.H` was created, `MTO_HF.C` untouched.

## Locked derivative status

- **R_w/J**: BFINAL-003 CLOSED for *linearization correctness* (FD oracle momentum relL2
  1.653e-4 / P-total 1.528e-4; J/J^T dot 6.42e-14; explicit/matrix-free 3.5e-16). **NOT**
  validated for *invertibility* — this round shows the exported operator is numerically
  singular (see below). Status: **closure-to-oracle still holds, but the operator is not
  directly invertible**.
- **R_x**: BFINAL-005 COMPLETE (momentum row + pressure/continuity row, machine-precision
  contraction). Used as-is in T1 (`stageB6_rxc_analytic.mtx`). No R_x defect implicated by
  this round's evidence.
- **g_w** (`pressureConstraintDerivative`): direct P-block of `explicitRhs_pressureDrop.mtx`;
  inlet-only support; `sum(g_w,P) = 2.357e-5` (not gauge-neutral because the outlet is
  fixedValue p=0).

## Pressure-Reference Gate (S1) — ALL 5 PASS (independently re-verified)

| # | check | result |
|---|---|---|
| C1 | 0/p outlet `fixedValue uniform 0` → `p.needReference()==false` → `setRefCell` no-op → pRefCell stays 0 (fvSolution pRefCell=5600 ineffective) | PASS |
| C2 | unique identity row in pressure block `[100800,134400)` via "single \|v\|>1e-12 AND diagonal==1" | **100800** (14 stored = 13 explicit zeros + diag 1.0) |
| C3 | identity row == `discretePIndex(0)` == 3·33600+0 | PASS |
| C4 | row 106400 (=discretePIndex(5600)) NOT identity (18 stored / 18 nonzero, physical continuity row, diag −1.899e-09) | PASS |
| C5 | pRefRow=100800 recorded | PASS |

## Tangent solve table (T1) — ACCEPTANCE NOT MET

`J = transpose(explicitJT.mtx)` re-pinned row 100800→identity, `rhs[100800]=0`, `rhs=-R_x·d`
(physical, no extra U/P scaling), scipy `splu` (COLAMD) + iterative refinement (30 iters).

| direction | \|\|R_x·d\|\| | true tangent residual (target ≤1e-9) | U-row residual | P-row residual | \|\|wPrime\|\| | wPrime_p[100800] |
|---|---|---|---|---|---|---|---|
| D1 | 5.663185e-01 | **3.647081e-04** | 3.647055e-04 | 4.021298e-03 | 2.084470e+18 | 0 |
| D2 | 8.615839e-01 | **4.044110e-05** | 4.044098e-05 | 1.926119e-04 | 5.539179e+17 | 0 |
| D3 | 8.516528e-01 | **3.092642e-03** | 3.092637e-03 | 1.843539e-02 | 3.796498e+19 | 0 |

Floor is 3–5 orders of magnitude above the 1e-9 target. `D_TAN` computed from these
solutions (+5.85 / +1.19 / −3.31) is **INVALID** (null-component-dominated, must not be used
for T2/T3).

## State tangent table (T2) / Functional tangent table (T3)

**NOT EXECUTED** — the round stopped at the T1 floor per the plan's escape clause
("If the exact solver cannot reach this level for a numerical reason, report the floor and
STOP rather than weakening the operator-validation logic"). No frozen-primal ±FD points were
generated, no eps sweep, no `D_FD`. (The re-converged-frozen-primal FD platform from
`validateStageB2GradientAmplitude.H` exists but was not exercised.)

## Root cause (verified)

The re-pinned J is **numerically singular**:
- `||J·n||/||n|| = 2.716846e-19` for the null vector n (independent matvec check);
- `sigma_min ≈ 7.4e-27 … 8.6e-26` (fresh random solves, independent LU);
- the null mode has `n[100800]=0` (NOT the constant-pressure mode) and is concentrated at the
  OUTLET plane (top-1% = 336 cells carry 89.5% of |p| mass; big-P cells at x∈[0.037,0.040],
  U-concentration at x=0.040 = outlet);
- outlet cells' P rows are physical continuity rows (16–18 nonzeros), NOT a Dirichlet p=0
  constraint row → the outlet `fixedValue p=0` BC is not enforced as a constraint row in the
  exported Jacobian, so an outlet pressure null mode survives the P(0) pin;
- RHS is nearly orthogonal to the null mode (`(-R_x d)^T n/||-R_x d||` = 1.3e-9 / 4.4e-11 /
  1.5e-9) and `g_w^T n = 2.8e-18`, so the system is consistent in exact arithmetic, but LU
  roundoff (cond ~1e26) sets the residual floor.

## Regression anchors

Zero production-math diff this round (git status byte-identical to pre-state; locked-head
sha256 unchanged). BFINAL-003/005 anchors therefore remain as locked (no re-run of `MTO_HF`
was required because no source changed; the artifacts used are byte-identical to cycle-1).

## Reviewer independent reproduction (this invocation)

1. **Identity-row re-scan** (`cycle-2/post-review/rev_nullcheck.py`, fresh scipy path):
   unique identity row = 100800; row 106400 = 18 nonzero physical row. Reproduces S1.
2. **Null-vector matvec** (same script): `||J·n||/||n|| = 2.716846e-19`; `n[PREF]=0`; P-block
   top-1% mass fraction 89.5%; RHS orthogonality 1.297e-9 / 4.421e-11 / 1.452e-9;
   `g_w^T n = 2.805e-18`. Reproduces the Executor's singularity claim exactly.
3. **Fresh LU tangent-solve rerun** (`cycle-2/post-review/rev_tangent_rerun.py`, background,
   fresh `splu` + fresh random sigma_min + fresh D1/D2/D3 solves + refinement):
   - sigma_min min = 7.382e-27;
   - D1 trueRelRes = 3.647081e-04 (Urel 3.647055e-04, Prel 4.021298e-03, ||w|| 2.084470e+18);
   - D2 = 4.044110e-05 (||w|| 5.539179e+17); D3 = 3.092642e-03 (||w|| 3.796498e+19);
   - w[PREF]=0 for all three.
   These reproduce the Executor's T1 numbers to the last digit — independent, deterministic.
4. **Null-mode geometry** (`cycle-2/post-review/rev_nullmode_geometry.log`): big-P cells at
   x∈[0.037,0.040] (outlet), U at x=0.040. Confirms outlet localization.

## Final decision

**FAIL_PLAN_OR_HYPOTHESIS** — the root-cause assumption ("cycle-1 failure was purely the
pressure-reference/pin identification, not a defect in J") is WRONG. The corrected pin is
necessary but insufficient: the exported closed J is numerically singular (outlet pressure
null mode), so the tangent system `J·w' = -R_x·d` cannot be solved to 1e-9 by the approved
direct solver. The research-level equivalent is **FAIL_TANGENT_CLOSURE** (task §12 / §20): the
tangent chain J + R_x + forward solve is NOT closed.

Do NOT let Flash patch J or R_x directly. Replanner must form a new single primary hypothesis
(e.g., enforce the outlet fixedValue p=0 as a hard Dirichlet constraint row in the exported J,
or solve the singular system null-space-consistently mirroring the C++ Schur-reduced adjoint).

STOP — no BFINAL-007.
