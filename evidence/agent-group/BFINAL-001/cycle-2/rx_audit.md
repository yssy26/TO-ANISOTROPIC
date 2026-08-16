# D3 audit — Rx-B / R_x completeness (fixed-state design residual derivative)

Stage: B-final (DIAGNOSTIC_ONLY, cycle 2)
Step: D3 of plan-002-replan (BFINAL-INDETERMINATE-PROBE-CONVENTION)
Author: Executor (flash), evidence-only; no src/ or case modification.
Date/HEAD: workspace /home/ys/dsH/TO-ANISOTROPIC @ df695c7 (uncommitted B4 edits
preserved); reference repo /home/ys/TO-ANISOTROPIC @ 292c711.
Evidence base: cycle-2 rerun.log (reproduces log.b4export), src/NS.H,
src/solveDiscreteFlowAdjoint.H, src/sensitivity.H, src/AdjNS_PD.H,
src/createFrozenHotRegionFields.H, src/validateDiscreteObjectiveDerivatives.H,
src/writeCSVLog.H, handoff §4.6-4.7/§5.2/§8 (B-F1..B-F3)/§19, and the
installed OpenFOAM-7 sources under /opt/openfoam7 (relax()/H()/Sp()
semantics verified from source, not from memory).

---

## 0. Verdict (one paragraph)

The code's pressure sensitivity assembles only the momentum-Brinkman row of
`R_x` (`gsenshPressureDrop = -dAlphaDxh*(U&Uc)*V`, sensitivity.H L67-68), i.e.
it *assumes* `R_P,x = d(div phi)/dalpha = 0` at fixed (U,p) and also assumes
`dphi/dalpha = 0` inside the momentum row. First principles say both
assumptions fail generically: `d rAU/dalpha = -rAU^2 != 0`, and the SIMPLE
flux satisfies `dphi_f/dalpha = -(rAU^2*H)_f*Sf_f - rAU_f^2*snGrad(p)_f*|Sf|_f
!= 0` wherever `H != 0` (i.e. `U != 0`) or `grad(p) != 0` (pressure-driven
flow: everywhere), with no mechanism forcing the two terms to cancel — the
converged SIMPLE state only enforces `div(phi(alpha0)) = 0` at the base
alpha. The existing Rx-B probe measures exactly this object and gets
`|dRp/dalpha|L2 = 8.1e-10` (cell level) and `|dPhi/dalpha|L2 = 2.1e-9`
(face level) — *nonzero*, triggering the probe's own "continuity dRp/dalpha
is NONZERO" message (rerun.log L146-148) and thereby contradicting the
"current analytic Rx_P = 0" assumption it prints. But Rx-B cannot settle
hypothesis (b) because it is self-referential (rebuilds phi with the same
analytic SIMPLE formula family it tests, boundary flux pinned to base by
construction), uses a sine direction (`sin(0.371(celli+1))*designMask`,
NS.H L156-161) rather than the design direction `dAlphaDxh`, and uses a
mixed relaxed-diagonal convention (`A0 = UEqn.A() - alpha` with UEqn
RELAXED; with U-relax = 0.4 this is `A0 = A_base/0.4 + 1.5*alpha`, verified
from fvMatrix::relax() `D /= alpha` and the run log max
`A0 = 1.50004681441e8 = 1e8/0.4 - 1e8`), so its absolute 8.1e-10 has no
normalization against the momentum term and cannot be read as "R_P,x = 0 is
fine". (b) and (c) remain OPEN: no pass criterion in the B-final sequence
asserts `R_P,x = 0` — handoff §8 B-F3 explicitly requires the FD to answer
"Is R_P,x*d actually zero?" — and the minimal discriminator is the B-F3
design-direction, fixed-state, split `R_U,x / R_P,x` FD specified in
gate_spec.md (D4).

---

## 1. Scope

Read-only audit of:

1. the Rx-B probe in NS.H (L136-410) — direction, step, A0 convention,
   phi rebuild, boundary handling;
2. the exact `R_x` the final sensitivity actually uses (sensitivity.H
   L63-68) and the adjoint RHS that produces `lambda = (Uc, pc)`
   (AdjNS_PD.H L1-29; solveDiscreteFlowAdjoint.H adjoint solve);
3. the first-principles content of `dR_P/dalpha` at fixed (U,p);
4. what hypothesis (c) (RHS/sign/final assembly) still needs.

No src/ or case file was modified; only this evidence note is written.

---

## 2. What the sensitivity actually uses as R_x (the (c) map)

### 2.1 Final assembly — src/sensitivity.H L63-68

```
if (adjointMode == "discrete")
{
    ...
    fsenshMeanT = -dAlphaDxh*(U & Ub);   fsenshMeanT *= mesh.V();
    gsenshPressureDrop = -dAlphaDxh*(U & Uc);  gsenshPressureDrop *= mesh.V();
```

- `dAlphaDxh = -alphaMax*(1+qu)*qu/sqr(qu+xh+SMALL)` (sensitivity.H L13-15)
  is the design-direction chain derivative `d(alpha)/d(xh)` (per cell,
  zeroed outside the design mask and at the alpha-clip, L25-35).
- `Uc` is the discrete adjoint velocity; `pc` (adjoint pressure) appears
  NOWHERE in the pressure sensitivity. Per-cell formula, momentum-Brinkman
  only.
- Sign convention: `gsens = -lambda^T * dR/dxh` with `dR_U/dxh =
  dAlphaDxh*U*V` — exactly the object Rx-A verifies
  (`dR_U/dalpha = dA*U*V`, relL2=1.6e-14, rerun.log L144). The assembly is
  therefore self-consistent *with the momentum-Brinkman-only R_x*; it is
  not self-consistent with the full reduced-system R_x (see §6).

### 2.2 Adjoint RHS — src/AdjNS_PD.H L1-29

```
vectorField discreteFlowVelocityRhs(mesh.nCells(), vector::zero);
scalarField discreteFlowPressureRhs(pressureConstraintDerivative);
```

- Momentum RHS = 0 (`dJ/dU = 0`), pressure RHS =
  `pressureConstraintDerivative = dJ/dp` of the implemented pressure-drop
  functional
  `J = rhoFluid*(mean p_in - mean p_out)/pressureDropMaxPa` (area-weighted;
  createFrozenHotRegionFields.H L624-662; runtime value writeCSVLog.H
  L85-92; pressureDropMaxPa = 25000 in optProperties).
- The g_w unit test (validateDiscreteObjectiveDerivatives.H L57-127)
  validates `dJ/dp` against a central FD of `J(p)` along
  `sin(0.371(celli+1))` with eps = 100: pressure = 7.08e-12 (threshold
  1e-9) — PASS (rerun.log L119). The adjoint RHS is therefore the correct
  `dJ/dx` of the implemented J; the RHS is anchored.
- Adjoint solve: `J^T lambda = dJ/dx` via applyScaledDiscreteFlowJT
  (momentum/area scaling cancels at the fixed point; pRef row zeroed,
  solveDiscreteFlowAdjoint.H L2229-2268); direct imported solution relRes
  = 3.48e-9 (rerun.log; handoff §5.4); transpose self-consistency 3.3e-16
  (handoff §5.3). Uc fingerprint nonzero after the pressureDrop solve and
  zero after the thermal solve (rerun.log L1281/L1603/L1617) — plumbing is
  wired, but no independent check that Uc solves J^T lambda = dJ/dx beyond
  the import residual.

### 2.3 What (c) still needs

All still OPEN:

1. B-F4 tangent oracle: solve `J*w' = -R_x*d` accurately and form
   `g_w^T*w'`; compare with the frozen-primal design FD (D1/D2/D3 ground
   truth -2.537 / +0.523 / +6.477, handoff §5.2). This bypasses the
   adjoint: if tangent != FD then the state/design model (J or R_x) is
   wrong ((a)/(b)); if tangent == FD but adjoint != tangent then it is a
   (c) issue (RHS sign / Lagrangian convention / final assembly sign).
2. B-F5 tangent-adjoint identity on the *same* J, R_x, g_w.
3. An audit that `gsenshPressureDrop` contracts the exact solved lambda
   (Uc) with the same scaling used in the adjoint solve (currently only the
   import relRes 3.48e-9 + fingerprint attest to this).

(c) cannot be closed by the existing evidence: g_w anchors the RHS, the
transpose is self-consistent, the solve is accurate, yet the end-to-end
projection still disagrees 52x / 2.3x / sign-flip (handoff §5.2) — the
unvalidated links are J vs R_w ((a)), R_x completeness ((b)), and the
B-F4/B-F5 assembly/sign checks ((c)).

---

## 3. Rx-B exact trace (NS.H L136-410)

### 3.1 Direction and step (L143-167)

```
dAlphaDir[celli] = sin(0.371*(celli+1)) * designMask[celli]     // SINE, masked
alphaMaxVal = gMax(mag(alpha));  hAlpha = 1e-6*max(alphaMaxVal,1e-12)  // = 100
alphaP = alpha + hAlpha*dAlphaDir;  alphaM = alpha - hAlpha*dAlphaDir
```

- Direction is a fixed-index sine wave (same seed as the g_w test), NOT the
  design direction `dAlphaDxh` (sensitivity.H L13-15). A JVP along an
  arbitrary sine direction does not test the direction that enters the
  sensitivity, and alternating signs can partially cancel in the L2 norm.
- `hAlpha = 100` on `alpha ~ 1e8` (max) — relative perturbation ~1e-6.

### 3.2 The A0 relaxed-diagonal mix (L45, L193-244)

```
UEqn.relax();                                  // NS.H L45, BEFORE the probe
A0 = UEqn.A() - alpha;                         // NS.H L243-244
rAUmod = 1/(A0 + alphaMod);                    // NS.H L287-291
```

- `UEqn.A()` is the RELAXED per-volume diagonal. Verified from
  /opt/openfoam7 fvMatrix.C relax(const scalar alpha):
  `D /= alpha;  S += (D - D0)*psi` — the diagonal is DIVIDED by the
  relaxation factor (U-relax = 0.4, fvSolution), so
  `A_relaxed = (A_base + alpha)/0.4` and
  `A0 = A_relaxed - alpha = A_base/0.4 + 1.5*alpha`.
- Run-log confirmation (rerun.log L145): max `UEqn.A() = 250004681.441 ~
  1e8/0.4`, max `A-alpha = 150004681.441 = 1e8/0.4 - 1e8` — all `A-alpha`
  values positive (avg +9.55e7, min +3.9e5).
- Hence `rAUmod = 1/(A0 + alphaMod) = 1/(A_base/0.4 + 2.5*alphaMod -
  alpha + ...)` is the rAU of NO physical matrix: in solid cells
  (alpha = 1e8) the ±100 perturbation is diluted against
  `A0 + alpha ~ 2.5e8`; in fluid cells the alpha-dependence is distorted by
  the `1/0.4` factor.

> **Correction to D2 audit.md §1.10 / §4.9.** D2 claimed
> `diag' = alpha*diag` (multiplication) and `A0 ~ (alpha_relax-1)*alpha ~
> -3e7` in solid cells. The OpenFOAM-7 source says `D /= alpha`, and the
> run log's all-positive `A-alpha` (min +3.9e5, max +1.50004681441e8)
> rules out the -3e7 formula. The relaxed-diag artifact is real but its
> sign/magnitude is `+1.5*alpha` (with U-relax = 0.4), not `-0.3*alpha`.
> The D2 conclusion that the artifact swamps the ±100 perturbation (and
> hence explains Rx-B ~ 0) survives; the formula and numbers are corrected
> here.

### 3.3 The phi rebuild (L287-355) — self-referential

```
Hbase  = UEqn.H()                      // relaxed-matrix H (incl. relax source)
phiHbyAmod = fvc::flux(rAUmod*Hbase)   // constrainHbyA + adjustPhi
if simple.consistent(): ... (NOT active: no "consistent" in fvSolution)
pEqnMod = fvm::laplacian(rAUmod, pbase)
phiMod  = phiHbyAmod - pEqnMod.flux()
rpMod   = div(phiMod) on internal faces + BASE phi on boundary faces
```

- This is the *same* analytic SIMPLE reconstruction family
  (`phi = phiHbyA - pEqn.flux()` with `rAU = 1/A`) whose alpha-dependence
  the "current analytic Rx_P = 0" assumption (L404) relies on. The probe
  therefore tests the analytic phi formula against itself; it cannot
  validate whether the analytic sensitivity NEEDS a `dR_P/dalpha` term.
- Boundary faces use the base phi ("U,p fixed semantics"), i.e. the
  boundary contribution to dRp/dalpha is zero BY CONSTRUCTION (admitted in
  the code comment NS.H L388-391). The inlet/outlet cells — exactly where
  the pressure-drop objective lives — are excluded from the measurement.
- Result (rerun.log L146-148):
  `|dRp/dalpha|L2 = 8.14730839032e-10`, max 2.39e-11, hAlpha = 100;
  face-level `|dPhi/dalpha|L2 = 2.13348655414e-09`; the probe's own
  `nRP > SMALL` branch fires: "continuity dRp/dalpha is NONZERO".

---

## 4. First principles: dRp/dalpha != 0 at fixed (U,p)

### 4.1 Setup

SIMPLE-condensed reduced system, state `w = (U, p)`, phi a dependent
variable `phi(w; alpha) = phiHbyA(w;alpha) - pEqn(w;alpha).flux()` with

```
rAU(alpha) = 1/(A_base + alpha)          // A_base: non-Brinkman per-volume diag
phiHbyA    = flux(rAU*H(U))              // H(U) frozen w.r.t. alpha (H() has no diag term)
pEqn.flux()_f = - rAU_f * snGrad(p)_f * |Sf|_f
phi_f      = (rAU*H)_f * Sf_f + rAU_f * snGrad(p)_f * |Sf|_f     // no consistent term in this case
```

(no `consistent` flag in fvSolution, so the rAtU branch of NS.H L308-318 is
not executed; `primalPressureMobility = rAU`.)

Residuals of the reduced system:

```
R_U(U,p;alpha) = A(alpha, phi)*U - source(U,p)      // momentum (frozen-coefficient J = dR_U/d(U,p))
R_P(U,p;alpha) = div( phi_SIMPLE(U,p;alpha) )       // continuity
```

### 4.2 The chain

```
d rAU/dalpha = -1/(A_base+alpha)^2 = -rAU^2  != 0        (rAU > 0 everywhere)
d phi_f/dalpha = - (rAU^2*H)_f * Sf_f  -  rAU_f^2 * snGrad(p)_f * |Sf|_f
```

- Term 1 (`-(rAU^2*H)*Sf`): nonzero wherever `H(U) != 0`. At the SIMPLE
  fixed point `HbyA = rAU*H ~ U`, so `H ~ A*U != 0` wherever `U != 0`
  (fluid cells, channels). (The relaxation inflation of H is cancelled by
  the relaxed rAU in the product `rAU*H ~ U`; the derivative term
  `rAU^2*H = rAU*(rAU*H) ~ rAU*U` is physical.)
- Term 2 (`-rAU^2*snGrad(p)*|Sf|`): nonzero wherever `grad(p) != 0` — a
  pressure-drop-constrained flow has `grad(p) != 0` essentially everywhere,
  including the Brinkman-dominated solid cells where `U ~ 0` (Darcy-like
  pressure gradient through the porous medium).
- There is no mechanism forcing Term 1 + Term 2 to vanish identically: the
  converged SIMPLE state only enforces `div(phi(alpha0)) = 0` at the base
  alpha; nothing constrains the alpha-derivative of the continuity residual
  at the fixed point.

### 4.3 Consequence

```
dR_P/dalpha = div( dphi/dalpha )  != 0   generically
dR_U/dalpha|fixed = dAlphaDxh*U*V  +  div( dphi/dalpha (*) U )   // phi-mediated part != 0
```

(the second momentum term is the variation of the convection operator
`fvm::div(phi, U)` through `dphi/dalpha`; nonzero wherever `U != 0` and
`dphi/dalpha != 0`).

### 4.4 Measured confirmation (this run)

- `|dPhi/dalpha|L2 = 2.13e-9` (internal faces, rerun.log L147) — nonzero
  face-level flux derivative;
- `|dRp/dalpha|L2 = 8.15e-10` (cell level, L146) — nonzero cell-level
  continuity derivative;
- probe prints "RxProbe Rx-B ... continuity dRp/dalpha is NONZERO"
  (L146-148 fires the NONZERO branch) — the probe's own output contradicts
  the "current analytic Rx_P = 0" label it prints.

The absolute magnitudes are small (rAU^2 is small in this case), but the
probe gives no normalization against the momentum term's contribution to
the sensitivity (`-Uc * R_U,x * d_hat`), so "small in absolute terms" says
nothing about "negligible for the gradient". That is exactly B-F3's job.

---

## 5. Why Rx-B's ~0 cannot close hypothesis (b)

1. **Self-referential.** `phiMod = phiHbyAmod - pEqnMod.flux()` with
   `rAUmod = 1/(A0+alphaMod)` is the same analytic SIMPLE formula family
   whose alpha-dependence the "analytic Rx_P = 0" assumption relies on; the
   boundary flux is pinned to the base phi by construction (NS.H L388-391),
   so inlet/outlet cells are excluded. A probe cannot validate an analytic
   assumption by re-running the same analytic formula.
2. **Non-design direction.** `dAlphaDir = sin(0.371*(celli+1))*designMask`
   (NS.H L156-161), not `dAlphaDxh = -alphaMax*(1+qu)*qu/sqr(qu+xh)` (the
   direction that enters the sensitivity). A JVP along a sine direction
   tests `dR_P/dalpha * d_hat_sine`, not `dR_P/dalpha * d_hat_design`;
   alternating-sign directions can partially cancel in the L2 norm.
3. **Mixed relaxed-diag convention.** `A0 = UEqn.A() - alpha` with UEqn
   relaxed (U-relax = 0.4) gives `A0 = A_base/0.4 + 1.5*alpha` (verified:
   `D /= alpha` in fvMatrix::relax; log max A0 = 1e8/0.4 - 1e8). The
   ±100*sin perturbation is diluted against `A0+alpha ~ 2.5e8` in solid
   cells, and `rAUmod` is the reciprocal of a non-physical matrix; in fluid
   cells the alpha-dependence is distorted by the 1/0.4 factor.
4. **No normalization.** `|dRp/dalpha|L2 = 8.1e-10` is an absolute L2 over
   cells; it is never compared with `|Uc * R_U,x * d_hat|` or
   `|pc * R_P,x * d_hat|`, so no significance statement follows.
5. **Incomplete object.** The probe measures only the continuity-row cell
   divergence; it does not measure the phi-mediated momentum row
   `div(dphi/dalpha * U)` nor the boundary-flux derivative at inlet/outlet
   — both of which are part of the full reduced-system `R_x` (see §6).

Therefore Rx-B cannot settle (b): the "≈0" is an artifact of a
self-referential, non-design-direction, relaxed-diag-mixed, boundary-pinned
construction — and even so it measures NONZERO.

---

## 6. Missing pieces of the code's R_x (the (b) statement)

Relative to the reduced-system `R_x = dR/dalpha` at fixed (U,p), the
sensitivity assembles only:

```
R_x(code) = ( dAlphaDxh*U*V ,  0 )      // momentum Brinkman only, continuity = 0
```

Missing by construction:

1. **Continuity row:** `R_P,x = div(dphi/dalpha)` — the `-pc * R_P,x *
   dAlphaDxh` term is absent (no pc anywhere in gsenshPressureDrop).
2. **Phi-mediated momentum row:** `div(dphi/dalpha * U)` (convection of U by
   the alpha-induced flux change) — absent; Rx-A's "fixed phi" assumption
   (NS.H L169-190) is contradicted by Rx-B's own face-level
   `|dPhi/dalpha|L2 = 2.1e-9 != 0`.
3. **Boundary-flux derivative** at the inlet/outlet cells — absent (the
   probe pins boundary phi by construction).

Whether the missing terms are numerically significant against the momentum
term is an open quantitative question — that is the B-F3 discriminator, not
an assumption to carry forward.

---

## 7. Status of (a)/(b)/(c) after D3

- **(a) J != R_w — OPEN.** Constrained by T1 = 0.85%, T-dev = 5.4e-5
  (dominant momentum blocks consistent for that direction). Unvalidated J
  links: the SIMPLE dphi/dU tangent mixes unrelaxed `rAUAdj` (dH part) with
  relaxed `primalPressureMobility = rAtU` (dp part) — factor ~1/0.4
  mismatch vs the primal phi built entirely with relaxed rAU/rAtU
  (NS.H L113-129, L1016-1033); the missing `-bc*dU_b` boundarySource term
  (D2 §3); upwind velocityJump lumping. Discriminator: B-F2 (single-
  convention J*v vs central FD of the actual frozen residual, momentum +
  continuity rows separately; spec in gate_spec.md / D4).
- **(b) R_x incomplete — OPEN and LEADING.** First-principles
  `dRp/dalpha != 0` (§4); the probe itself measures NONZERO (§4.4); the
  code's R_x misses the continuity row, the phi-mediated momentum row, and
  the boundary-flux derivative (§6). Discriminator: B-F3 — fixed-state,
  design-direction, split `R_U,x / R_P,x` central FD with the REAL design
  direction `dAlphaDxh`, epsilon sweep, and the missing terms
  `-pc*R_P,x*d_hat`, `-Uc*div(dphi/dalpha*U)*d_hat` compared against the
  momentum term `-Uc*(dAlphaDxh*U*V)*d_hat`.
- **(c) RHS/sign/final assembly — OPEN.** RHS anchored (g_w = 7.08e-12),
  transpose self-consistent (3.3e-16), direct solve accurate (3.48e-9);
  the end-to-end chain is still 52x / 2.3x / sign-flip vs frozen FD
  (handoff §5.2). Discriminators: B-F4 (tangent oracle) and B-F5
  (tangent-adjoint identity) after (a)/(b) — see §2.3.

**No pass criterion asserts `R_P,x = 0`.** The "current analytic Rx_P = 0"
(NS.H L404) is an assumption printed by a probe whose own measurement
contradicts it; handoff §8 B-F3 explicitly requires the fixed-state FD to
answer "Is R_P,x*d actually zero?" before the momentum-only assembly can be
trusted, and B-F3 is a mandatory step of the B-final diagnostic sequence
(handoff §8, §19).

---

## 8. Feed-forward to D4 (gate_spec.md)

B-F3 must use: fixed (U,p) (frozen turbulence), design direction
`d_hat = dAlphaDxh` (real chain derivative), central FD with epsilon sweep,
residual convention R_U = A*U - source of the unrelaxed matrix
`div - laplacian + Sp(alpha)` (with dev2 when rans) and
`R_P = div(phi_SIMPLE)` rebuilt through the SAME primal path (rAU from the
UNRELAXED A(), phiHbyA = flux(rAU*H), pEqn = laplacian(rAU,p),
phi = phiHbyA - pEqn.flux(), boundary flux REBUILT not pinned), rows split
(U index 3*cell+cmpt, P index 3N+cell), and the three candidate missing
terms of §6 quantified against the momentum term. Threshold proposal (per
plan D4): momentum and continuity block-wise J*v vs FD < 1e-3 for B-F2;
R_P,x*d vs FD < 1e-3 before declaring R_x complete for B-F3.

---

## 9. Files/lines consulted (read-only)

- src/NS.H: L45 (UEqn.relax), L113-129 (rAU/HbyA/phiHbyA/rAtU/
  primalPressureMobility), L136-190 (Rx-A), L156-167 (sine direction,
  hAlpha=100), L192-410 (Rx-B: A0 = UEqn.A()-alpha at L243-244; rAUmod at
  L287-291; phiMod = phiHbyAmod - pEqnMod.flux() at L319-332; boundary
  pinned at L341-349 and L388-391; NONZERO branch at L405-409),
  L1016-1033 (primal pEqn with rAtU, phi = phiHbyA - pEqn.flux()).
- src/sensitivity.H: L3-15 (dAlphaDxh), L25-35 (mask/clip), L40-94
  (discrete branch, fingerprint L43-62, final assembly L65-68).
- src/AdjNS_PD.H: L1-29 (discrete RHS: 0 / pressureConstraintDerivative),
  L33-131 (continuous branch, not used in discrete mode).
- src/solveDiscreteFlowAdjoint.H: L22-99 (discretePrimalMomentum UNRELAXED,
  discreteDiagX=diag+internalCoeffs, rAUAdj = 1/A()), L625-641 (SIMPLE dphi
  tangent: rAUAdj dH part + mobF=primalPressureMobility dp part),
  L646-657 (velocityJump, continuity row), L2229-2268 (adjoint RHS scaling,
  pRef row), plus D2 §3 for the full J row mapping.
- src/createFrozenHotRegionFields.H: L582, L624-662 (dJ/dp of the
  area-weighted pressure-drop functional).
- src/validateDiscreteObjectiveDerivatives.H: L57-127 (g_w pressure FD).
- src/writeCSVLog.H: L85-92 (runtime J = rhoFluid*(pbar_in - pbar_out)).
- case: system/fvSolution (U-relax 0.4, p-relax 0.2, no SIMPLE consistent,
  pressureDropMaxPa 15000 in fvSolution but 25000 in constant/optProperties
  — the code path uses the optProperties value), constant/optProperties.
- /opt/openfoam7: fvMatrix.C relax() (`D /= alpha; S += (D-D0)*psi`),
  H() (`boundaryDiag*psi + lduMatrix::H(psi) + source + boundarySource`,
  then `/V`), fvmSup.C (fvm::Sp: `diag += V*sp`).
- Evidence: cycle-2 rerun.log L144-148 (Rx-A/SpDiag/Rx-B raw+faces+NONZERO),
  L119 (g_w), L1281/L1603/L1617/L1630 (Uc fingerprints, dDP range),
  L1250-1262 (T1/T-cont/T2/T-fvm/T-dev/T-conv, GateH1 — see D1/D2);
  AGENT_GROUP_HANDOFF.md §4.6-4.7, §5.2-5.4, §6.3, §8 (B-F1..B-F5), §19.
- Companion notes: evidence/agent-group/BFINAL-001/cycle-2/pre-state.md,
  d1-summary.md, audit.md (D2), rerun.log.
