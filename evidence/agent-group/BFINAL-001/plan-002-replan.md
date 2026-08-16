# BFINAL-001 Replanner round — DIAGNOSTIC_ONLY (replan after REJECTED_HYPOTHESIS)

## Workspace check (re-verified this round)
- pwd = /home/ys/dsH/TO-ANISOTROPIC
- toplevel = /home/ys/dsH/TO-ANISOTROPIC  (== required shared workspace)
- HEAD = df695c7373297425e2d76767d2c415a9f7f55b18
- uncommitted: M src/solveDiscreteFlowAdjoint.H, M src/validateStageBRepeatability.H (+ untracked scripts/logs)
- case optProperties: stageB4JacobianProbe=true, discreteExportOnly=true,
  discreteUseExplicitSolution=true, adjointMode=discrete, mmaUpdateEnabled=false.

## Why the previous hypothesis was rejected (accepted by replanner)
1. It declared (a) "J != R_w, localized to momentum U-U convection block" the single
   surviving cause, but the same log also shows T-fvm relL2 = 0.999999999866 (~100%)
   and frozen-phi dH relL2 = 4.67e7 — i.e. the fvm part fails too, so "only convection
   fails" is internally inconsistent.
2. It excluded (b) R_x-incomplete using Rx-B |dRp/dalpha|L2 = 8.1e-10, but Rx-B rebuilds
   phi(alpha) with the SAME analytic SIMPLE formula it is testing, uses a sine direction
   (not the design direction), and mixes a relaxed-diag convention (A0 = UEqn.A() - alpha).

## New evidence gathered (read-only source audit, this round)
### Residual conventions in the B4 probe suite are mutually inconsistent
- T1 FD reference = momentumResidualB4 (solveDiscreteFlowAdjoint.H L1783-1818):
  R_U = mesh.V()*(M & UU) + 2.0*H_ldu(UU), which equals A*UU - source (diag + offdiag,
  no /V on the final result, deviatoric folded into source via M -= fvc::div(...dev2...)).
  T1 also perturbs phi: phiper = phibase + fvc::flux(dU) (ad-hoc, SIMPLE-free dphi).
  => T1 relL2 = 0.0085 (the only COMPLETE full-residual comparison; it is clean).
- T-fvm reference (L2074-2076): Afd = mesh.V()*(discretePrimalMomentum & dUField) + 2*H_ldu(dU)
  = A_raw*dU with NO boundary-diag and NO SIMPLE dphi, while Jfvm (L2042-2063) uses
  discreteDiagX = diag + boundary-diag adjustment (L30-71) plus raw upper/lower.
  => T-fvm relL2 = 0.999999999866 is a CONVENTION MISMATCH of the reference, not a proof
  the fvm block is wrong.
- T-conv (L2166-2204) uses a separate resWithHldu formula; frozen-phi dH analytic uses
  "-offdiag*dU - alpha*V*dU" (NS.H L729-753) vs the /V-scaled H() FD. Both use different
  hand-rolled conventions => T-conv=0.418, frozen-phi dH=4.67e7, same-basis dH=0.117 are
  decomposition artifacts, mutually inconsistent with T1.
=> Conclusion: the block-wise numbers DO NOT localize any single broken block. They show
   the decomposition REFERENCES are convention-inconsistent; only T1 is a valid full test.

### Rx-B is self-referential and contradicts first principles
- Rx-B (NS.H L267-354): phiMod = phiHbyA(alphaMod) - pEqnMod.flux(), with
  phiHbyAmod = fvc::flux(rAUmod*Hbase), rAUmod = 1/(A0 + alphaMod), A0 = UEqn.A() - alpha
  (RELAXED diag mixed with raw alpha). Same analytic SIMPLE path as the thing it tests.
- First principles at fixed (U,p): d(rAU)/dalpha = -1/(A0+alpha)^2 = -rAU^2 != 0, so
  dphiHbyA/dalpha and d[laplacian(rAU,p).flux()]/dalpha are NONZERO in general.
  => |dRp/dalpha|L2 = 8.1e-10 is a probe artifact (cancellation/roundoff/mixed convention),
     NOT evidence that R_P,x = 0. Hypothesis (b) stays OPEN (leading suspect per handoff §4.6, §8 B-F3).

### What is still genuinely unvalidated (all three hypotheses remain OPEN)
- (a) J != R_w: T1 clean at 0.85% argues the momentum row J ~ R_w for THAT direction, but
  J's reduced SIMPLE phi tangent dphi = Sf&(rAU*dH) + mobF*dp*dcf*maf is only compared against
  ad-hoc flux(dU), never against the true SIMPLE phi(alpha,U) reconstruction. OPEN.
- (b) R_x incomplete: Rx-B self-referential; no TRUE fixed-state design-direction FD of
  R_U,x and R_P,x. OPEN (leading suspect).
- (c) RHS/sign/final assembly: g_w=7.08e-12 and -lambda^T R_x convention are self-consistent,
  but the end-to-end 52x/2.3x/sign-flip is NOT explained and no independent FD validates the
  final gsenshPressureDrop = -dAlphaDxh*(U&Uc)*V (sensitivity.H L65-68) chain. OPEN.

### Ruled out (unchanged)
- linear solver (direct residual 3.48e-9), transpose self-consistency (3.3e-16),
  projection chain (0.1555186511), volume chain. g_w anchor = 2.7e-11 / 7.08e-12.

## Revised primary hypothesis (CANDIDATE, no survivor named)
hypothesis_id = BFINAL-INDETERMINATE-PROBE-CONVENTION
The B4 probe suite is internally convention-inconsistent (block-decomposition references mix
raw-coefficient / V-scaled / factor-2-H / boundary-diag / source conventions), and Rx-B is
self-referential with a non-design direction, so the suite cannot discriminate (a)/(b)/(c).
No hypothesis is excluded. The minimal missing diagnostic is a single-convention block-wise
B-F2 (J*v vs central FD of the actual residual, momentum + continuity rows, epsilon sweep)
plus design-direction B-F3 (fixed-state R_x*d vs design FD, split R_U,x and R_P,x).

## Next numerical Gate (for the following round, not patched now)
- B-F2: block-wise relL2(J*v, central FD of actual R_w) < 1e-3 for BOTH momentum and
  continuity rows over an epsilon sweep with a plateau.
- B-F3: relL2(R_x*d, design-direction fixed-state central FD) < 1e-3 for R_U,x AND R_P,x,
  using the REAL design direction d (not a sine).
