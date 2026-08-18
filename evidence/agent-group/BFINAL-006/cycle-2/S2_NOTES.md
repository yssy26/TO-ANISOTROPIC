# BFINAL-006 cycle-2 — S2-T1 step notes (Executor)

## Status
**BLOCKED_BY_NEW_EVIDENCE** — T1 tangent linear-system closure could not meet acceptance.

## What was done (exactly per approved plan S2-T1)
1. Pressure-Reference Gate re-run in-driver from the ACTUAL case/source/artifact: all 5 checks PASS
   (unique identity row 100800 = discretePIndex(0); row 106400 physical continuity, NOT identity;
   pRefRow=100800 recorded).  See `T1_tangent_solve.log` lines 3-17.
2. T1: J = transpose(explicitJT.mtx) (physical J^T -> J; scales M/A recorded, no extra scaling);
   re-pin J row 100800 -> identity; rhs[100800]=0; rhs = -R_x*d (stageB6_rxc_analytic.mtx +
   stageB6_dirs.mtx, BFINAL-005 closed R_x); scipy splu (COLAMD) + iterative refinement (30 iters).
3. Reported ||rhs||, true relative residual, U-row and P-row residuals per direction; wrote
   wPrime_TAN_{D1,D2,D3}.mtx; recorded D_TAN = g_w^T wPrime (g_w = P block of
   explicitRhs_pressureDrop.mtx = production pressureConstraintDerivative; no adjoint quantity).

## Result vs acceptance
| dir | trueRelRes | U-rowRel | P-rowRel | wPrime_p[100800] | target |
|---|---|---|---|---|---|
| D1 | 3.647e-04 | 3.647e-04 | 4.021e-03 | 0 | <=1e-9 (pref 1e-10) |
| D2 | 4.044e-05 | 4.044e-05 | 1.926e-04 | 0 | <=1e-9 |
| D3 | 3.093e-03 | 3.093e-03 | 1.844e-02 | 0 | <=1e-9 |

Acceptance NOT met (3-5 orders above).  wPrime_p[100800]=0 satisfied.

## Root cause (evidence-backed; see S2_T1_FLOOR_EVIDENCE.md)
- The closed operator (transpose of explicitJT.mtx, column-pin; equivalently the C++ J^T with its
  row-P(0) pin) is NUMERICALLY SINGULAR: machine-precision null vector (||J n||/||n|| = 2.7e-19),
  sigma_min ~ 1e-26 physical / ~4e-20 with the closed adjoint's U/P scaling; the scaled tangent
  solve also stalls (3.9e-3).
- The null mode is a pressure mode concentrated at the OUTLET plane (84 outlet cells = 35% of |p|
  mass; outlet P rows are physical continuity rows, NOT p=0 Dirichlet rows).  The exported
  Jacobian does not enforce the outlet fixedValue p=0 as a constraint row, so a pressure null
  mode at the outlet survives the P(0) pin.
- The C++ adjoint converges because it iterates with a Schur-reduced block preconditioner
  (solveDiscreteFlowAdjoint.H L3126-3133) and its rhs g_w is exactly orthogonal to the null mode
  (inlet-only support; g_w^T n = 2.8e-18).  Direct splu on the full operator cannot exploit this.
- D_TAN values (+5.85 / +1.19 / -3.31) are INVALID (non-converged solve; solution dominated by
  the arbitrary null component) and must not be used for PASS/FAIL.

## Deviations
None from the approved procedure.  Additional read-only diagnostic work (conditioning, scaled
solve, null-mode structure) was performed to characterize the floor, per the plan's T1 escape
clause; no production file was touched.

## Files written (evidence only)
- cycle-2/tangent_solve.py, cycle-2/T1_tangent_solve.log
- cycle-2/t1_floor_diagnosis.py, cycle-2/t1_floor_diagnosis.log
- cycle-2/nullmode_analysis.py
- cycle-2/S2_T1_FLOOR_EVIDENCE.md (this analysis)
- cycle-2/artifacts/wPrime_TAN_D{1,2,3}.mtx (garbage solutions, NOT usable)
- cycle-2/artifacts/T1_summary.json
