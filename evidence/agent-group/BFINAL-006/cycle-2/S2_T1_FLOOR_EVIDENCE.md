# BFINAL-006 cycle-2 — S2/T1 result: TANGENT SOLVE FAILS ACCEPTANCE (numerical floor, structural)

Status: **BLOCKED_BY_NEW_EVIDENCE** — the closed BFINAL-003 operator (transpose of the
exported `explicitJT.mtx`, re-pinned per the approved plan) is numerically SINGULAR, so
the approved T1 procedure (scipy splu, no extra U/P scaling) cannot reach the required
true relative residual <= 1e-9 (floor 3.6e-4 / 4.0e-5 / 3.1e-3).  This is a structural
property of the closed operator (outlet fixedValue p=0 is NOT enforced as a Dirichlet
constraint row in the exported Jacobian), not a solver-tuning issue.

---

## 1. Provenance

- workspace : `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace, verified)
- branch    : `agent/dsH-stage-b-validation`
- HEAD      : `ca8a772b8339a36e7906ea32a7125061c47aa280` (unchanged)
- git status: dirty set identical to the locked BFINAL-003/005 baseline
  (M src/sensitivity.H, M src/solveDiscreteFlowAdjoint.H,
   M src/solveDiscreteFlowAdjointProduction.H + pre-existing untracked).  ZERO new
  production changes by this step (only evidence files under
  `evidence/agent-group/BFINAL-006/cycle-2/`).
- case      : `/home/ys/dsH/b2_case_smoke`
- artifacts (sha256, unchanged from cycle-1 pre-state):
  explicitJT.mtx `3dfcc9ab41f9344649d31638eac161c0b64e9603b2907b36527f42af88924940`
  stageB6_rxc_analytic.mtx `d30e13179e5dcd5065a4b5f509ad24d000efda7499c1c0ce5fc7dd4897bc984c`
  stageB6_dirs.mtx `7572261514f8a8eee2c352772c1b19961c7245ca9539340813f916e0cfd75bb5`
  explicitRhs_pressureDrop.mtx `2d7cf846829762427298f9a2bb732e5dcf83ce6d8d93b7a77ec13650ab6bdb5f`

## 2. Pressure-Reference Gate (re-verified in-driver from the ACTUAL case/source)

| # | check | result |
|---|---|---|
| C1 | 0/p outlet `fixedValue uniform 0` (7 other patches zeroGradient) -> `p.needReference()==false` -> `setRefCell` (readTransportProperties.H L66-68) no-op -> pRefCell stays 0 -> fvSolution pRefCell=5600 INEFFECTIVE | PASS |
| C2 | identity-row scan of explicitJT.mtx pressure rows [100800,134400) with CORRECTED criterion "single \|v\|>1e-12 AND diagonal==1" -> unique identity row **[100800]** (14 stored entries = 13 explicit zeros + diag 1.0; nrow==1 criterion would fail as in cycle-1) | PASS |
| C3 | identity row 100800 == discretePIndex(0) = 3*33600+0 | PASS |
| C4 | row 106400 (=discretePIndex(5600)): 18 stored / 18 nonzero, diag -1.899e-09, physical continuity row, NOT identity | PASS |
| C5 | pRefRow = 100800 recorded as T1 pinning row | PASS |

## 3. T1 solve as executed (per approved plan)

- J = transpose(explicitJT.mtx) (physical J^T -> J; exported csrVal is physical —
  scales M=0.00786409044975, A=2.47771624976e-07, M/A=31739.2697832 applied only in the
  C++ scaled wrapper; recorded, no extra scaling applied).
- Re-pin J row pRefRow=100800 -> identity; rhs[100800]=0 (wPrime_p[100800]=0 gauge;
  pPrime[5600]=106400 NOT zeroed).
- rhs = -R_x*d per direction (stageB6_rxc_analytic.mtx [anRUd(3N); anRPd(N)] per dir,
  stageB6_dirs.mtx 3xN); scipy splu (COLAMD) + iterative refinement (30 iters).

| dir | \|\|R_x*d\|\| | \|\|rhs(post-pin)\|\| (U / P) | \|\|J w'-rhs\|\| | trueRelRes | U-rowRel | P-rowRel | \|\|wPrime\|\| | wPrime_p[100800] | D_TAN (INVALID) |
|---|---|---|---|---|---|---|---|---|---|
| D1 | 5.663185e-01 | 5.663185e-01 (5.663184e-01 / 1.930627e-04) | 2.065409e-04 | **3.647e-04** | 3.647e-04 | 4.021e-03 | 2.084e+18 | 0 | +5.84756577e+00 |
| D2 | 8.615839e-01 | 8.615839e-01 (8.615837e-01 / 4.479779e-04) | 3.484340e-05 | **4.044e-05** | 4.044e-05 | 1.926e-04 | 5.539e+17 | 0 | +1.19303916e+00 |
| D3 | 8.516528e-01 | 8.516528e-01 (8.516527e-01 / 2.595428e-04) | 2.633857e-03 | **3.093e-03** | 3.093e-03 | 1.844e-02 | 3.796e+19 | 0 | -3.31183111e+00 |

Acceptance target true relative residual <= 1e-9 (prefer <= 1e-10): **NOT MET** (3-5 orders
of magnitude above).  wPrime_p[100800]=0 holds.  No adjoint quantity used in D_TAN
(load-bearing inputs: R_x + J + forward solve + direct g_w = P block of
explicitRhs_pressureDrop.mtx = production pressureConstraintDerivative, +rho*A/(pMax*A_inlet)
on the 84 inlet cells, U block exactly zero).

## 4. Root-cause diagnosis (numerical floor is STRUCTURAL)

1. **The closed operator is singular.**  `||J·n||/||n|| = 2.7e-19` for the solution vector n
   (machine-precision null mode).  Random-solve sigma_min estimate of the physical re-pinned
   J: **7.1e-27** (cond ~ 1e26).  The null mode is intrinsic to the exported operator: it is a
   null vector of the column-pin-only J (as-exported transpose, WITHOUT the row re-pin) with the
   same ratio 2.7e-19, so it is NOT created by the plan's row re-pin.

2. **U/P scaling (the closed adjoint's momentum/area scaling) does not fix it.**
   Scaled tangent Mtilde = Dout^-1 J Din (row: U/M, P/A; col: P*M/A — the exact scaling of the
   closed BFINAL-003 adjoint operator applyScaledDiscreteFlowJT): sigma_min ~ 3.8e-20
   (cond ~ 2e22); scaled solve + refinement stalls at relres **3.9e-3** (29 iters); physical-space
   residual of the scaled solution 1.5e-1/4.3e-4.  Both representations are numerically singular.

3. **The null mode is a pressure mode concentrated at the OUTLET plane.**
   - 84 outlet cells (0.25% of cells) carry **35.2%** of the total |p| null-mode mass;
     mean|p|null outlet = 1.0e17 vs 7.2e14 global mean (140x).
   - The outlet cells' P rows in explicitJT.mtx are PHYSICAL continuity rows (16-18 nonzeros,
     U couplings ~1.25e-7, P couplings ~1e-10), NOT the Dirichlet p=0 row.
   - The only identity row in the whole matrix is 100800 = P(0) (an inlet-plane cell, the C++
     pRef pin).  => the outlet fixedValue p=0 BC is NOT enforced as a constraint row in the
     exported Jacobian; a pressure null mode at the outlet survives the P(0) pin.
   - The null mode is zero at the inlet plane (g_w support), so `g_w^T n = 2.8e-18` (machine
     zero): D_TAN = g_w^T wPrime is gauge-invariant w.r.t. the null component.  But the
     non-converged particular part makes the computed D_TAN untrustworthy -> INVALID evidence.

4. **Why the C++ adjoint works anyway.**  solveDiscreteFlowAdjoint.H L3126-3133 solves J^T
   lambda = g with a Schur-reduced BLOCK PRECONDITIONER (S = SIMPLE pressure-Poisson
   laplacian(rAtU), GAMG; zU = (A^T)^-1 rU, PBiCGStab) — it never factorizes/solves the full
   singular operator directly, and its RHS g_w is exactly orthogonal to the null mode
   (inlet-only support).  The approved T1 (direct splu on the full physical J, no scaling)
   cannot exploit this.

5. **Consistency check.**  (-R_x d)^T n / ||-R_x d|| = 1.3e-9 (D1), 4.4e-11 (D2), 1.5e-9 (D3):
   the tangent RHS is *nearly* orthogonal to the null mode (consistent system in exact
   arithmetic), but double-precision LU amplifies the null component (cond ~1e22-1e26), so the
   residual floor is set by roundoff, not by inconsistency.

## 5. Not modified / out of scope

- No production source, case file, or forbidden module touched (git diff identical to the
  locked baseline).  No attempt to patch J or R_x (plan forbids it; failure is reported, not
  fixed).

## 6. Evidence files (this step)

- `T1_tangent_solve.log` — full driver log (gate + T1 + summary).
- `t1_floor_diagnosis.py` / `t1_floor_diagnosis.log` — sigma_min estimates (physical 7.1e-27,
  scaled 3.8e-20), physical & scaled solve tests.
- `nullmode_analysis.py` — null-mode spatial analysis (outlet concentration, inlet zero,
  alpha correlation -0.04).
- `artifacts/wPrime_TAN_D1.mtx`, `wPrime_TAN_D2.mtx`, `wPrime_TAN_D3.mtx` — T1 solutions
  (garbage; dominated by the arbitrary null component — do NOT use for T2/T3).
- `artifacts/T1_summary.json` — machine-readable summary.

## 7. Recommendation for Replanner (NOT executed)

- The tangent linear system J w' = -R_x d cannot be closed with the current operator by a
  direct LU solve.  Two options to discriminate (next hypothesis, planner decision):
  (a) Solve the SAME operator with an iterative method mirroring the C++ adjoint
      (GMRES + Schur-reduced block preconditioner; solveDiscreteFlowAdjoint.H L3126-3133) —
      this is not a new operator, but it deviates from the approved splu procedure;
  (b) Treat the outlet Dirichlet p=0 as a hard constraint row in J (or drop/replace the
      outlet-concentrated null mode) — this modifies the operator and needs the planner's
      authorization; it also speaks directly to handoff hypothesis #2 ("implemented J is
      self-consistent with its transpose but not the exact derivative of the actual frozen
      primal residual").
- Either way the T1 acceptance (residual <= 1e-9) requires a plan revision; the step stopped
  per the plan's T1 escape clause ("If the exact solver cannot reach this level for a
  numerical reason, report the floor and STOP rather than weakening the operator-validation
  logic").
