# BFINAL-006 cycle-2 — Post-Reviewer independent verification notes

## 1. Workspace / scope (verified personally)

- pwd `/home/ys/dsH/TO-ANISOTROPIC`; git root == workspace; HEAD `ca8a772b8339a36e7906ea32a7125061c47aa280`; branch `agent/dsH-stage-b-validation`.
- `git status --porcelain` (38 lines) **byte-identical** to `cycle-2/pre-state/git_status_porcelain.txt` (diff empty) → **zero source/case changes this round**.
- Locked-head sha256 (recomputed now) match the BFINAL-003/005 locked heads:
  - `solveDiscreteFlowAdjoint.H` = f0c81497… (BFINAL-003 locked)
  - `solveDiscreteFlowAdjointProduction.H` = 8c4901bf… (BFINAL-003 locked)
  - `sensitivity.H` = f0ed0726… / `rxPressureRowTranspose.H` = fd648ae0… (BFINAL-005)
  - `stageB6RxDesignOracle.H` = 09832dd6… / `stageB5BoundaryRelaxOracle.H` = 0b274ee2…
- `stageB7TangentOracle.H` was **never created**, `MTO_HF.C` untouched → S3/S4 were correctly NOT executed (Executor stopped at the T1 floor per the plan's escape clause). No scope violation.

## 2. Independent reruns performed (this invocation)

### (a) Identity-row / pressure-reference re-verification (rev_nullcheck.py)
Fresh scipy path (independent of the Executor's streaming s1_gate.py):
- Unique identity row in the pressure block `[100800,134400)` = **100800** (criterion "single |v|>1e-12 AND diagonal==1").
- Row 106400 = 18 stored / 18 nonzero, physical continuity row, NOT identity.
- `discretePIndex(celli)=3N+celli` confirmed from `solveDiscreteFlowAdjoint.H` L113-120
  (`discreteUIndex=3*celli+cmpt`, `discretePIndex=discreteNVelocity+celli`), so
  `discretePIndex(0)=100800`, `discretePIndex(5600)=106400`. S1 gate is CORRECT.

### (b) Null-vector / singularity confirmation (rev_nullcheck.py, no LU needed)
- Built `J = transpose(explicitJT.mtx)` re-pinned at row 100800, loaded the Executor's
  `wPrime_TAN_D1.mtx`, computed `||J w||/||w||`:
  - `||w|| = 2.084470e+18`, `||J w|| = 5.663185e-01`, **`||J n||/||n|| = 2.716846e-19`**.
  - This independently reproduces the Executor's "2.7e-19" null-vector ratio → **J is
    numerically singular**.
- `n[PREF] = 0` → the null mode is **NOT** the constant-pressure mode (a constant mode
  would have n[100800]≠0 after the P(0) pin; the pin does NOT remove this mode).
- P-block localization: top-1% (336 cells) carry **89.5%** of |p| mass; nullmode geometry
  rerun (rev_nullmode_geometry.log) locates them at x∈[0.037,0.040] (the OUTLET plane,
  xmax=0.0398), U-concentration at x=0.040. → outlet-concentrated pressure null mode.
- RHS orthogonality: `(-R_x d)^T n/||-R_x d||` = 1.297e-09 / 4.421e-11 / 1.452e-09 (D1/D2/D3)
  → the tangent system is consistent in exact arithmetic but the null mode is not removed.
- `g_w^T n = 2.805e-18` (inlet-only support, orthogonal to the null mode).

### (c) Fresh LU tangent-solve rerun (rev_tangent_rerun.py, background)
Gold-standard rerun: fresh `splu(J)` + fresh random sigma_min + fresh D1/D2/D3 solves with
iterative refinement. See `rev_tangent_rerun.log` (results filled in when it completes).

## 3. Ordering / R_x alignment (verified correct)

- `stageB6_rxc_analytic.mtx` export (stageB6RxDesignOracle.H L1093-1099) writes
  `[anRUd(3N cell-major); anRPd(N)]` per direction — matches the Executor's
  `rxc[off:off+3N]` (U) + `rxc[off+3N:off+4N]` (P) split.
- `explicitRhs_pressureDrop.mtx` is cell-major `[Ux,Uy,Uz,p]` (solveDiscreteFlowAdjoint.H
  L3031-3037); the Executor's `gw_cell[3::4]` correctly extracts the P block.
- The transpose `J = (J^T)^T` is correct (exported CSR is the physical J^T, pinned row
  pRefRow; transposing yields column-pinned J; the driver's row re-pin restores the gauge).
- Conclusion: the Executor's T1 implementation matches the approved plan exactly; the
  failure is NOT an implementation defect.

## 4. Decision rationale

The primary hypothesis ("the locked chain J + R_x + forward tangent + direct g_w is
physically closed, and the cycle-1 failure was purely the pressure-reference/pin
identification, not a defect in J or R_x") is **REFUTED**:

1. The cycle-1 pin mis-identification (100800 vs 106400) was real and is now corrected (S1 PASS).
2. BUT re-pinning at the correct pRefRow=100800 does NOT close T1: the exported J (transpose
   of explicitJT.mtx, re-pinned) is **numerically singular** with a non-constant,
   outlet-concentrated pressure null mode (||J n||/||n|| = 2.7e-19, sigma_min ~1e-26), so the
   direct splu solve stalls at true relative residual ~3.6e-4 (target ≤1e-9).
3. The BFINAL-003 "J closed" evidence (FD oracle relL2 ~1.5-1.7e-4, J/J^T dot 6.4e-14,
   explicit/matrix-free ~3.5e-16) never tested INVERTIBILITY; it is fully consistent with a
   singular operator whose null mode is also a flat direction of the frozen primal residual.
   This finding promotes handoff hypothesis #2 ("implemented J is self-consistent with its
   transpose but not the exact nonsingular derivative of the frozen primal residual").

→ The planned remedy (fix only the pin identification, then direct-solve) is insufficient;
the root-cause assumption ("no defect in J") is wrong. **FAIL_PLAN_OR_HYPOTHESIS** (do NOT
let Flash patch J/R_x directly).

## 5. Files (this review)

- `rev_nullcheck.py` / `rev_nullcheck.log`
- `rev_nullmode_geometry.log`
- `rev_tangent_rerun.py` / `rev_tangent_rerun.log` (background LU rerun)
- `POST_REVIEW_NOTES.md` (this file)
