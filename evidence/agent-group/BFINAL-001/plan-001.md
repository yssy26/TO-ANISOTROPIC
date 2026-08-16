# BFINAL-001 Planner round — DIAGNOSTIC_ONLY plan notes

Workspace check (Planner, read-only):
- pwd        = /home/ys/dsH/TO-ANISOTROPIC
- toplevel   = /home/ys/dsH/TO-ANISOTROPIC  (matches required shared workspace)
- HEAD       = df695c7373297425e2d76767d2c415a9f7f55b18  (B2/B3 checkpoint proxy)
- uncommitted: M src/solveDiscreteFlowAdjoint.H, M src/validateStageBRepeatability.H,
  plus many untracked scripts/logs/evidence (B4 worktree ≈ B4 commit 292c711).

## Objective (this round)
Stage B-final pressure-gradient physical-closure audit; DIAGNOSTIC_ONLY; no PATCH;
no MMA. Per AGENT_GROUP_HANDOFF.md §8 (B-F1..B-F5) and §19.

## Key evidence gathered (read-only)
1. g_w (pressure objective state derivative) — src/validateDiscreteObjectiveDerivatives.H
   already passes on current case: thermal=2.7e-11, pressure=7.08e-12 (threshold 1e-9).
   This is the B-F1 anchor and confirms the adjoint RHS (pressureConstraintDerivative)
   equals dg/dp of the area-averaged pressure-drop functional.
2. RxProbe (src/NS.H, stageB4JacobianProbe=true), from /home/ys/b2_case_smoke/log.b4export:
   - Rx-A (fixed phi, Brinkman dR_U/dalpha): relL2=1.6e-14, cos=1  -> momentum R_x EXACT.
   - Rx-B (fixed U,p, phi(alpha) rebuilt, dRp/dalpha): |dRp/dalpha|L2=8.1e-10,
     max=2.4e-11  -> continuity/pressure-row design derivative ~= 0 (negligible).
   => hypothesis "R_x incomplete (missing pressure-row term)" is essentially EXCLUDED.
3. Stage B4 flow-Jacobian self-tests (src/solveDiscreteFlowAdjoint.H), log.b4export:
   - T1 (momentum U-U, dphi=flux(dU)): relL2=0.0085 (0.85%).
   - T-cont (continuity row): relL2=0.00065.
   - T2 (dp!=0 continuity): relL2=7.5e-21 (exact).
   - T-dev (deviatoric part): relL2=5.4e-5.
   - T-fvm (frozen-phi fvm matrix part): relL2=0.999999999866 (≈100%).
   - T-conv (convection phi-variation / velocityJump): relL2=0.418 (42%).
   - frozen-phi dH analytic vs FD: relL2=4.7e7; same-basis dH: relL2=0.117 (cos 0.9994).
   => the momentum U-U convection/flux-sensitivity (dphi/dU through SIMPLE phi=phiHbyA-pEqn.flux())
   is the failing block; Brinkman/deviatoric/pressure-coupling blocks are clean.
4. Transpose self-consistency (explicit CSR J^T == matrix-free J^T) = 3.3e-16 (handoff §5.3);
   direct SuperLU solve true residual = 3.48e-9 (handoff §5.4). Linear algebra ruled out.
5. Final assembly (src/sensitivity.H L65-68): gsenshPressureDrop = -dAlphaDxh*(U&Uc)*V
   (momentum Brinkman term only). Sign convention -lambda^T R_x is internally consistent
   given J^T lambda = g_w. No independent evidence implicates RHS/sign while T-conv fails directly.

## Primary hypothesis (hypothesis_id = BFINAL-J-NEQ-RW-CONV-DPHI)
The assembled B4 flow Jacobian J is NOT the exact state derivative R_w of the converged
frozen primal residual; the error is localized to the momentum U-U block convection/flux
sensitivity (dphi/dU through the SIMPLE flux reconstruction), so the adjoint velocity Uc
is wrong and the pressure gradient fails FD. g_w, R_x (both rows), transpose self-consistency,
and direct-solve accuracy are all verified, leaving J=R_w as the only unclosed link.

## Smallest missing diagnostic (next numerical Gate = B-F2)
A clean full-Jacobian column FD test:
  J*v  vs  [ R(w + eps*v) - R(w - eps*v) ] / (2*eps),
with R = the ACTUAL NS.H frozen-primal residual semantics (momentum fvMatrix::residual()
at the converged state for the U block; div(phi) for the continuity block), an epsilon
sweep, and block-wise relL2 report. Gate criterion: momentum U-U block must be < ~1e-3
for J=R_w; a reproducible relL2 ~0.4 in the convection sub-block confirms this hypothesis
and localizes the fix (one-pass/reduced flux tangent vs full fixed-point flux Jacobian).

## What is OUT OF SCOPE this round
Implementing SIMPLE-transpose or any J fix; enabling MMA; touching src/ or the case;
changing objective/constraint/gradient scales/thresholds; thermal-objective debugging.
