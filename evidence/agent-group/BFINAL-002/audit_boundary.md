# BFINAL-002 S1 — Boundary-consistency audit (Q-A map) at HEAD ca8a772b

Stage: B-final (DIAGNOSTIC_ONLY, BFINAL-002, step S1-audit-boundary)
Author: Executor (flash). READ-ONLY: no src/, case, or physics modification.
HEAD: ca8a772b8339a36e7906ea32a7125061c47aa280, branch agent/dsH-stage-b-validation
(git root /home/ys/dsH/TO-ANISOTROPIC == session workspace; dirty tree = untracked
src/validation_stage_*_build.log only).
Purpose: answer Q-A's audit question — map internal phi / physical-boundary phi /
fixedValue vs non-fixed U and p patches / adjustPhi / constrainHbyA /
constrainPressure / pRefCell / pEqn.flux / final reconstructed phi in NS.H and
the J boundary terms in solveDiscreteFlowAdjoint.H, cross-checked against the
case files — and state, per existing FD reference, whether it holds boundary
face fluxes at baseline.

---

## 0. Executive summary

1. **The production phi reconstruction at this HEAD is:**
   `UEqn.relax()` (NS.H L45) → `rAU = 1/UEqn.A()` (L113, RELAXED,
   avg 2.8234e-07) → `HbyA = constrainHbyA(rAU*UEqn.H(), U, p)` (L114) →
   `phiHbyA = fvc::flux(HbyA)` (L115) → `adjustPhi` (L116, **NO-OP** in this
   case) → `rAtU = rAU` (consistent=false, L118-127) → `constrainPressure`
   (L1016, **NO-OP** in this case) → `pEqn = fvm::laplacian(rAtU,p)` +
   `pEqn.setReference(pRefCell, pRefValue)` (L1020-1025, **NO-OP** in this
   case) → `phi = phiHbyA - pEqn.flux()` (L1032).
2. **Boundary semantics reduce dramatically in this case.** Because `p` has a
   fixedValue outlet, `p.needReference() == false`, which makes **adjustPhi
   and pEqn.setReference(pRefCell,pRefValue) no-ops**; because no patch is
   `fixedFluxPressure`, **constrainPressure is a no-op**. The only active
   boundary mechanisms are (i) **constrainHbyA**, which overwrites `HbyA_b =
   U_b` at the 7 non-assignable (fixedValue/noSlip) U patches, and (ii) the
   fixedValue-p outlet's `pEqn.flux()` boundary.
3. **Production boundary flux tangent (dphi_b/dU, dp=0): zero at 7 of 8
   patches** (inlet, hotInlet, hotOutlet, solidEndWalls, bottomWall, topWall,
   sideWalls — fixedValue/noSlip U, zeroGradient p: `phiHbyA_b = Sf_b&U_b`
   pinned by constrainHbyA, `pEqn.flux()_b = 0`), and
   `Sf_b & d(rAU_relaxed*H_relaxed)_b/dU` at the outlet (zeroGradient U,
   fixedValue p) only.
4. **J (applyDiscreteFlowJ L414-711) adds a boundary dphi at ALL 8 patches**
   (L689-697: `P(celli) += Sf_b & (rAUAdj*dH_b)`, extrapolated-flux
   convention, unrelaxed rAUAdj), i.e. it models a boundary flux tangent that
   production does NOT have at the 7 fixedValue-U patches. **This single term
   accounts for the entire P-row boundary mismatch** (|diff|L2 boundary =
   2.69e-06 ≈ |Jv_P| total = 3.59e-06, localize_prow.out).
5. **Every existing phi-rebuild FD reference is convention-blind to the
   production boundary semantics in a different way** (see the mechanism
   matrix, §6):
   - NS.H rebuilt-phi P-row FD (`stageB2_rpDudU_rebuilt.mtx`): fresh
     UNRELAXED matrix + **boundary flux PINNED to the stored base phi**
     (NS.H L676-681) → boundary dphi ≡ 0 by construction.
   - GatePR/GateOracle FD: fresh UNRELAXED matrix + boundary flux REBUILT
     but **without constrainHbyA/adjustPhi/pEqn.flux** (extrapolated at all
     patches) → self-consistent with J's boundary term by construction
     (relL2 ~ 1e-10), matching NEITHER production.
   - Frozen-phi FD (`stageB2_ruDudU.mtx`): phi not perturbed at all.
6. **Explicit answer to the acceptance question:** YES — the NS.H rebuilt-phi
   P-row FD (the reference behind the 0.824 result) **holds boundary face
   fluxes at baseline** (pins `phi.boundaryField()` into the manual
   divergence, L676-681); GatePR does not (rebuilds them, L1497-1502) but
   also omits constrainHbyA; production pins them at the 7 fixedValue-U
   patches via constrainHbyA and rebuilds only at the outlet.
7. **Q-A prediction (to be settled by the Gate-A oracle, S2-S4):** the 0.824
   is mechanically a reference-side artifact (J's boundary dphi term vs
   base-pinned FD), but a fully production-consistent FD (constrainHbyA
   pinning + relaxed rAU + relax-source H) is expected to keep the P-boundary
   mismatch O(1) at the 7 fixedValue-U patch cells because J's boundary dphi
   term is not the production boundary flux tangent there. Which outcome
   (collapse vs O(1) localization) is exactly what Gate A must measure.
8. **Q-B prediction:** the 2.5 = 1/U-relax ratio is exact for the mobility
   (verified numerically, §5.3), but the full dphi/dU tangent also involves
   the relax-source term `(D_relaxed-D0)*psi` inside `fvMatrix::H()`
   (1.5*D0*psi/V with U-relax 0.4) which contributes a `0.6*dU`-like diagonal
   to `d(rAU*H)/dU` — partially compensating the 0.4 factor. No existing
   probe includes it (all rebuild with fresh unrelaxed matrices). The oracle
   must rebuild phi through the exact NS.H sequence (relaxed H + relaxed rAU)
   to measure the net effect; the 2.5 ratio alone is NOT proof.

---

## 1. Repository / case identity

```
pwd                       /home/ys/dsH/TO-ANISOTROPIC
git rev-parse --show-toplevel  /home/ys/dsH/TO-ANISOTROPIC
git branch --show-current agent/dsH-stage-b-validation
git rev-parse HEAD        ca8a772b8339a36e7906ea32a7125061c47aa280
git status --short        only ?? src/validation_stage_*_build.log (untracked)
```
Matches the BFINAL-002 baseline. Case: `/home/ys/dsH/b2_case_smoke`
(writable; contains the B4-rerun stageB2_*.mtx exports). Read-only reference
exports: `/home/ys/b2_case_smoke/*.mtx` (explicitJT.mtx, stageB2_*.mtx).

Key run facts from `evidence/agent-group/BFINAL-001/cycle-2/rerun.log`
(produced by the same B4 code path):
```
g_w pressure FD error        = 7.07793481308e-12   (L119, PASS)
NS UEqn A()  avg             = 158965329.116        (RELAXED, ~1e8/0.4)
NS rAU = 1/A avg             = 2.82337769549e-07    (RELAXED)
Adj rAUAdj  avg              = 7.05844426338e-07    (UNRELAXED, discretePrimalMomentum.A() avg 63586131.6461)
rAUAdj/rAU                   = 2.5000000087  (= 1/0.4 = 2.5, U-relax 0.4)
|phiRb0 - phi|L2             = 3.49247878283e-11    (RELAXED rebuild == stored phi)
deltaA |dA_FD|L2             = 0                    (A is U-independent at frozen phi)
GatePR h=1e-3  relL2         = 1.12295958888e-10    (J P-row vs rebuilt-boundary UNRELAXED FD)
GatePRint h=1e-3 internal    = 1.46436377217e-10
crosscheck_bf2: momentum     = 1.653123e-04 cos=0.999999986
crosscheck_bf2: P-row total  = 8.240837e-01 cos=0.695444677
localize_prow: internal      = 1.090701e-07 cos=1.000000000
localize_prow: boundary      = 1.407930e+00 cos=0.246014468
localize_prow: |diff|L2      ALL=2.693799e-06 internal=2.890788e-13 boundary=2.693799e-06
face-level dphi: J vs GateOracle-FD relL2=9.003456e-10; J vs NS.H-rebuilt 4.717604e-08
```

---

## 2. Case configuration (cross-checked against the case files)

### 2.1 system/fvSolution
```
SIMPLE { nNonOrthogonalCorrectors 1; pRefCell 5600; pRefValue 0; }
relaxationFactors { fields p 0.2; equations U 0.4; (k|omega) 0.4; }
solvers: p -> PCG/DIC tol 1e-9 relTol 0; U -> smoothSolver GaussSeidel 2 sweeps.
NO "consistent" entry in SIMPLE -> simple.consistent() == false.
constraintFunctionDict_flow { constraintPatchesNames (inlet outlet);
                              pressureDropMaxPa 15000.0; rhoFluid 0.5892; }
```
`pRefCell = 5600` is read by `src/readTransportProperties.H` L66-68
(`setRefCell(p, simple.dict(), pRefCell, pRefValue)`), pRefValue = 0.

### 2.2 0/U (velocity boundary conditions)
| patch | type | assignable() | fixesValue() | faces |
|---|---|---|---|---|
| inlet | **fixedValue (193.97 0 0)** | false | true | 84 |
| outlet | **zeroGradient** | true | false | 84 |
| hotInlet | fixedValue (0 0 0) | false | true | 196 |
| hotOutlet | fixedValue (0 0 0) | false | true | 196 |
| solidEndWalls | noSlip (fixedValue subclass) | false | true | 280 |
| bottomWall | noSlip | false | true | 1120 |
| topWall | noSlip | false | true | 1120 |
| sideWalls | noSlip | false | true | 4800 |

Non-fixed U patches: **outlet only** (84 faces). Fixed U patches: 7 of 8.

### 2.3 0/p (pressure boundary conditions)
| patch | type | fixesValue() |
|---|---|---|
| outlet | **fixedValue 0** | true |
| inlet, hotInlet, hotOutlet, solidEndWalls, bottomWall, topWall, sideWalls | zeroGradient | false |

Non-fixed p patches: 7 of 8; fixed p patch: **outlet only**.

> **CORRECTION to the approved-plan case-config description.** The plan's
> "evidence" bullet says "0/U inlet zeroGradient/outlet fixedValue 0". The
> actual case (both `/home/ys/dsH/b2_case_smoke` and `/home/ys/b2_case_smoke`)
> is the opposite for U: **inlet = fixedValue (193.97 0 0), outlet =
> zeroGradient**. 0/p matches the plan (inlet zeroGradient / outlet fixedValue
> 0). This does not change the audit conclusions but is recorded because the
> fixedValue-inlet/noSlip-wall set is exactly the 7-patch set where
> constrainHbyA pins HbyA_b and where J's boundary dphi term is
> production-inconsistent (§4-§5).

### 2.4 system/fvSchemes
```
div(phi,U)      bounded Gauss upwind;      (upwind -> J's velocityJump downwind lumping)
laplacianSchemes default Gauss linear corrected;   (nNonOrthogonalCorrectors 1)
snGradSchemes   default corrected;
fluxRequired { default no; p; pa; pb; pc; }        (pEqn.flux() is legal for p)
```

### 2.5 constant/optProperties (relevant switches)
```
adjointMode discrete; flowModel incompressibleRANSFrozen;
freezeTurbulenceForValidation true; mmaUpdateEnabled false;
stageB4JacobianProbe true; discreteExportOnly true;
discreteExplicitMatrixFile "/home/ys/b2_case_smoke/explicitJT.mtx";
pressureDropMaxPa 25000.0 (used by the code; fvSolution's 15000 is not used);
pressureGradientScale 1.0; objectiveGradientScale 1.0.
```
No `constant/fvOptions` / `system/fvOptions` in the case → `fvOptions(U)` is
a zero source and `fvOptions.constrain`/`fvOptions.correct` are no-ops.

---

## 3. Production phi reconstruction — NS.H (exact lines at ca8a772b)

```
L26-34  UEqn = fvm::div(phi, U) - fvm::laplacian(nuEffFrozen, U) + fvm::Sp(alpha, U)
        == -fvc::grad(p) + fvOptions(U)
L37-43  if (ransFlowModel) UEqn -= fvc::div(nuEffFrozen*dev2(T(grad(U))))   // source += V*div(dev2)
L45     UEqn.relax();            // D /= 0.4 ; S += (D - D0)*U   (fvMatrix::relax)
L47     fvOptions.constrain(UEqn);       // no-op (no fvOptions)
L52     solve(UEqn == -fvc::grad(p));
L113    rAU = 1.0/UEqn.A();              // RELAXED mobility (avg 2.8234e-07)
L114    HbyA = constrainHbyA(rAU*UEqn.H(), U, p);
L115    phiHbyA = fvc::flux(HbyA);
L116    adjustPhi(phiHbyA, U, p);        // NO-OP: p.needReference()==false
L118-127 rAtU = rAU; if (simple.consistent()) {...}   // consistent==false -> rAtU == rAU
L129    primalPressureMobility = rAtU().primitiveField();   // RELAXED mobility
L1016   constrainPressure(p, U, phiHbyA, rAtU());   // NO-OP: no fixedFluxPressure patch
L1020-1023 pEqn = fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA);
L1025   pEqn.setReference(pRefCell, pRefValue);    // NO-OP: needReference()==false
L1032   phi = phiHbyA - pEqn.flux();       // final nonOrthogonal iteration only
L1038-1039 p.relax(); U = HbyA - rAtU()*fvc::grad(p);
L1040-1041 U.correctBoundaryConditions(); fvOptions.correct(U);
```

### 3.1 Why adjustPhi / setReference / constrainPressure are no-ops here
Verified from /opt/openfoam7 sources:
- `GeometricField::needReference()` (GeometricField.C ~L1100): returns false as
  soon as ANY patch `fixesValue()`. `p` has a fixedValue outlet → **false**.
- `adjustPhi` (cfdTools/general/adjustPhi/adjustPhi.C): entire body inside
  `if (p.needReference())` → skipped.
- `fvMatrix::setReference(celli, value, forceReference=false)`
  (fvMatrix.H L348-353, fvMatrix.C L505): applies only if
  `forceReference || psi_.needReference()` → skipped (forceReference defaults
  false).
- `constrainPressure` (cfdTools/general/constrainPressure/constrainPressure.C):
  only touches `fixedFluxPressureFvPatchScalarField` patches; none exist →
  skipped.
- `pEqn.flux()` (fvMatrix.C L863): legal (fluxRequired p). Internal faces:
  `lduMatrix::faceH(p)`. Boundary: `flux_b = internalCoeffs_b*p_internal_b -
  boundaryCoeffs_b*p_neighbour_b` (non-coupled; neighbour term only for
  coupled patches). For the laplacian(rAtU,p): zeroGradient-p patches → flux_b
  = 0; fixedValue-p outlet → `flux_b = -rAtU_b*dc_b*(p_b - p_int_b)`-type
  term, which depends on **internal p** and **rAtU**, not on U (dp=0 test).

### 3.2 Boundary semantics that ARE active in production
- **constrainHbyA** (cfdTools/general/constrainHbyA/constrainHbyA.C):
  `if (!U.boundaryField()[patchi].assignable() && p not fixedFluxExtrapolatedPressure)
   HbyA_b = U_b;` → at inlet, hotInlet, hotOutlet, solidEndWalls, bottomWall,
  topWall, sideWalls (7 patches): **HbyA_b pinned to U_b** (fixed base value).
  At the outlet (zeroGradient U, assignable): HbyA_b = extrapolated
  `rAU*H` (fvMatrix::H() returns an `extrapolatedCalculated`-type boundary).
- **phiHbyA_b = fvc::flux(HbyA)_b**: at the 7 pinned patches =
  `Sf_b & U_b` (base, U-perturbation-independent for the U-column); at the
  outlet = `Sf_b & (rAU*H)_b` (U-dependent).
- **pEqn.flux()_b**: nonzero only at the fixedValue-p outlet; U-independent at
  dA/dU = 0 with dp = 0 (depends on frozen p and rAtU).
- **phi_b = phiHbyA_b - pEqn.flux()_b**:
  - 7 fixedValue-U patches: phi_b = `Sf_b & U_b` → **dphi_b/dU = 0**.
  - outlet: phi_b = `Sf_b & (rAU_relaxed*H_relaxed)_b - pEqn.flux()_b` →
    dphi_b/dU = `Sf_b & d(rAU_relaxed*H_relaxed)_b/dU` (pEqn.flux()_b term
    U-independent).
- Walls (noSlip U_b = 0, zeroGradient p): phi_b = 0 (impermeable) ✓.

### 3.3 Relaxation inside H(): the relax-source term
`fvMatrix::relax()` (fvMatrix.C ~L672) does `D /= alpha` and
`S += (D - D0)*psi`; `fvMatrix::H()` (fvMatrix.C) = `[BdiagH*psi +
lduMatrix::H(psi) + source + boundarySource]/V` with `BdiagH =
-internalCoeffs_cmpt + cmptAv(internalCoeffs)` (per component, all patches).
Therefore, with U-relax = 0.4:
```
D_relaxed = D0/0.4 ;  S_relax = S0 + 1.5*D0*U
H_relaxed(U) = H_unrelaxed(U) + 1.5*(D0/V)*U        (per volume)
rAU_relaxed = 0.4*rAU_unrelaxed
rAU_relaxed*H_relaxed = 0.4*rAU_unrelaxed*H_unrelaxed + 0.6*U
d(rAU*H)/dU = 0.4*rAUu*dHu/dU + 0.6*I               (dA/dU = 0, verified)
```
So the true production dphi/dU (dH part) is NOT simply 0.4× the unrelaxed
one: the relax-source adds a `0.6*dU` diagonal contribution that partially
compensates the 0.4 mobility factor. Net difference vs J's unrelaxed tangent
is `flux(0.6*(I - rAUu*dHu/dU))` — nonzero only through the off-diagonal /
boundary structure of `dHu/dU`. **None of the existing probes include the
relax-source term** (all rebuild H from fresh unrelaxed matrices, §6). This is
the exact "relaxation factor could in principle cancel through H() details"
gap identified in gate_spec.md (D4) §2.5 — it must be measured, not assumed.

---

## 4. J side — solveDiscreteFlowAdjoint.H (exact lines at ca8a772b)

### 4.1 Assembly
```
L22-28  discretePrimalMomentum = fvm::div(phi,U) - fvm::laplacian(nuEffFrozen,U)
                                 + fvm::Sp(alpha,U)        // NO relax, NO dev2, NO grad(p)
L30-32  discreteDiagX/Y/Z = diag + internalCoeffs (all patches)   // = boundaryDiag
L71-72  discreteLower/Upper = raw ldu coefficients
L76-81  rAUAdj = 1/discretePrimalMomentum.A()             // UNRELAXED (avg 7.0584e-07)
```

### 4.2 applyDiscreteFlowJ (L414-711) — boundary-relevant terms
```
L445-451  deviatoric deltaU: zeroed on fixedValue U patches, then
          correctBoundaryConditions (extrapolates at zeroGradient outlet)
L476-484  U-row diagonal: discreteDiagX/Y/Z * dU_cell
L502-582  dH tangent: internal faces dH -= upper*dU_nei (own), -= lower*dU_own (nei);
          boundary dH += (-ic_cmpt + cmptAv(ic))*dU_cell  (ALL patches, L562-574);
          /= V (L576-582). (Matches fvMatrix::H() boundaryDiag semantics, frozen coeffs.)
L597-623  U-row off-diagonal action (upper/lower)
L625-641  SIMPLE phi tangent (internal faces):
          dphi_f = Sf_f & ( wf*rAUAdj_o*dH_o + (1-wf)*rAUAdj_n*dH_n )   // UNRELAXED dH part
                 + mobF*(dp_nei - dp_own)*dcf*maf                        // mobF = primalPressureMobility
                                                                         //        = rAtU RELAXED (dp part)
          -> J's dphi is INTERNALLY MIXED: unrelaxed rAUAdj (dH) + relaxed mobF (dp)
L646-653  velocityJump: U-row(downwind) += dphi*(U_nei - U_own)   (upwind lumping)
L656-657  continuity row: P(own) += dphi ; P(nei) -= dphi
L675-697  boundary loop (ALL 8 patches):
L689-697  boundary dphi: P(celli) += Sf_b & (rAUAdj[celli]*dH_b)   // extrapolated-flux
                                                                   // convention, UNRELAXED,
                                                                   // NO constrainHbyA pinning
L699-708  boundary dp -> U-row: += Sf_b*dp_cell  only if !pressureFixed (7 patches)
L418-419  pressureFluxCorrection computed (fixedValue-p boundary dp flux, outlet) but
          NEVER added to output -> DEAD CODE; J's P-row has NO boundary dp term at the
          fixedValue-p outlet (T2 = 7.5e-21 passes because its reference's outlet
          boundary term is small relative to internal faces; not a tested link)
```
Absent from J's U-row (unchanged from D2): the `-boundaryCoeffs*dU_b`
boundarySource term of the true residual Jacobian; nonzero only at the
zeroGradient outlet where `dU_b = dU_int` (bc_outlet ~ 1.57e-4, rerun.log).
Small but unvalidated.

### 4.3 Explicit J^T COO assembly (L2271-2485) replicates applyDiscreteFlowJT
Includes the same mixed conventions: U-U diag/offdiag (L2292-2325), direct
P<-U (L2328-2334), kf = mobF*dcf*maf dp path (L2338-2346), hAdjoint
Stage1/Stage2 with rAUAdj (L2349-2449), boundary-flux transpose with rAUAdj
and Bdiag (L2451-2485), deviatoric transpose (L2487+). pRef P-row pinned to
identity in CSR (L2758-2780); ExplicitJToracle maxRelL2 = 3.9e-16.

---

## 5. The three boundary-flux conventions in play (Q-A core)

| # | Convention | Boundary flux tangent dphi_b/dU (dp=0) | Who |
|---|---|---|---|
| P | production (NS.H state map) | 0 at the 7 fixedValue/noSlip U patches (constrainHbyA pins HbyA_b=U_b, pEqn.flux()_b=0 there); `Sf_b&d(rAU_rel*H_rel)_b/dU` at the outlet | NS.H L45/L113-116/L1020-1032 |
| J | J's reduced phi tangent | `Sf_b&(rAUAdj*dH_b)` at ALL 8 patches (extrapolated, unrelaxed); no constrainHbyA pinning | applyDiscreteFlowJ L689-697 |
| F1 | NS.H rebuilt-phi FD | **0 everywhere** (boundary flux pinned to stored base phi) | NS.H L676-681 |
| F2 | GatePR/GateOracle FD | `Sf_b&d(rAUf*H)_b/dU` at ALL patches (extrapolated, unrelaxed; no constrainHbyA/adjustPhi/pEqn.flux) | solveDiscreteFlowAdjoint.H L1486-1502 |

J (row "J") == F2 by construction (same extrapolated-unrelaxed convention) →
GatePR passes at 1e-10. J != P at the 7 fixedValue-U patches. F1 != J (and
F1 != P at the outlet). The 0.824 total / 1.41 boundary result is J vs F1.

---

## 6. Mechanism matrix — does each reference include the mechanism?

| Mechanism | Production NS.H | J (applyDiscreteFlowJ) | GatePR/GateOracle FD | NS.H rebuilt-phi FD | frozen-phi FD (ruDudU) |
|---|---|---|---|---|---|
| internal phi (upwind div) | yes (L28, phi stored) | yes (coefficients at frozen phi; dphi tangent L625-641) | yes (fvm::div at frozen phi) | yes (fvm::div at frozen phi) | yes (frozen phi, no dphi) |
| physical-boundary phi, fixedValue-U patches | `Sf_b&U_b`, pinned (constrainHbyA L114) | dphi_b = `Sf_b&(rAUAdj*dH_b)` ≠ 0 (L689-697) — **inconsistent with production** | rebuilt extrapolated dphi_b ≠ 0 (no constrainHbyA) | **pinned to stored base phi (L676-681) → 0** | n/a (frozen) |
| physical-boundary phi, non-fixed U (outlet) | `Sf_b&(rAU*H)_b` extrapolated − pEqn.flux()_b | dphi_b = `Sf_b&(rAUAdj*dH_b)` (scale 2.5× if compared to relaxed production) | same extrapolated form (unrelaxed) | pinned → 0 | n/a |
| fixedValue U patches (perturbation handling) | U_b fixed (BC) | dU zeroed on fixedValue patches for dev2 (L445-451); dH boundaryDiag acts on cell value at ALL patches (L562-574) | UmodR.correctBoundaryConditions() (L1476) → U_b fixed | Umod.correctBoundaryConditions() (L631) → U_b fixed | Umod.correctBoundaryConditions() (L499) |
| non-fixed U patches (outlet extrapolation) | zeroGradient → U_b = U_int | same via correctBoundaryConditions | same | same | same |
| fixedValue p (outlet) | pEqn.flux()_b = fvm flux at outlet; phi_b = phiHbyA_b − flux_b | P-row: boundary dp NOT applied at pressureFixed (pressureFluxCorrection dead, L418-419); dp part only internal faces (L640) | not modeled (dp=0 test; phiF has no pEqn.flux()) | pEqnRb.flux() present (L657) but boundary flux contribution overridden by pinned base (L676-681) | n/a |
| non-fixed p patches | zeroGradient → pEqn.flux()_b = 0 | n/a (no boundary dp term there) | n/a | n/a | n/a |
| adjustPhi | **no-op** (needReference false) — but if it ran, a global massCorr on outflow | NOT modeled (no massCorr derivative anywhere) | NOT modeled | called (L656) but no-op for same reason | NOT modeled |
| constrainHbyA | yes (L114) | NOT modeled in dphi (boundary term uses extrapolated flux) | NOT modeled | yes (L654) but boundary flux still overridden at L676-681 | NOT modeled |
| constrainPressure | **no-op** (no fixedFluxPressure patches) | NOT modeled (no-op in this case anyway) | NOT modeled | NOT modeled | NOT modeled |
| pRefCell / setReference | **no-op** (needReference false); pRefCell 5600 read (readTransportProperties.H L66-68) | P(pRef) pinned to identity in the exported/adjoint matrix (L2267, L2758-2780); NOT part of J's physics rows | pRef row excluded from P-row J-vs-FD by the export pin (crosscheck_bf2.out) | pRef row = 0 in FD (no pin needed; comparison excludes it) | n/a |
| pEqn.flux() (dp part) | yes, production | internal dphi dp part (L640) with relaxed mobF; boundary dp part missing (dead code) | dp = 0 test only | pEqnRb.flux() internal part (L657); boundary overridden | n/a |
| final reconstructed phi | phi = phiHbyA − pEqn.flux() (L1032) | dphi = dphiHbyA − d(pEqn.flux()) with d(pEqn.flux())/dU = 0 at dp=0 | phiF = fvc::flux(rAU*H) only (no pEqn.flux() term) | phiRb = phiHbyArb − pEqnRb.flux() (L668) | n/a |
| U relaxation (rAU convention) | RELAXED (L45 before L113) | UNRELAXED rAUAdj (L76-81, L628-629, L689) | UNRELAXED rAUc/rAUf (L1467, L1486) | UNRELAXED rAUrb (L634-640) | UNRELAXED fresh matrices |
| relax-source in H() | yes (relax modifies source; H() includes it) | NO (discretePrimalMomentum never relaxed) | NO (tEqR fresh) | NO (tUEqnRb fresh) | NO (fresh) |
| boundarySource −bc*dU_b | in residual (BC) | absent (D2; only matters at zeroGradient outlet) | absent | absent | absent |

**Answer to the acceptance question:** the only FD reference that holds
boundary face fluxes at baseline is the NS.H rebuilt-phi P-row FD — it
explicitly adds the STORED base `phi.boundaryField()` into the manual cell
divergence (NS.H L676-681), so its boundary-cell rows have zero boundary dphi
by construction (the same code pattern is admitted in NS.H L388-391 for the
Rx-B probe). GatePR rebuilds the boundary flux (L1497-1502) but without
constrainHbyA; J includes a boundary dphi term; production pins via
constrainHbyA. The three conventions (J, GatePR, production) coincide only at
the outlet; they differ at the 7 fixedValue/noSlip U patches.

---

## 7. Where the 0.824 lives (mechanism, not yet verdict)

From localize_prow.out: internal cells agree to 1.09e-7 (cos = 1.0); the
entire discrepancy |diff|L2 = 2.69e-06 sits in boundary cells (relL2 1.41,
cos 0.246). The boundary-cell rows differ ONLY in the boundary-face flux
tangent (internal-face dphi agrees because both sides share the unrelaxed
fresh-matrix convention). Per-cell estimate of J's boundary dphi
(Sf_b ~ 1e-4…1e-7 m², rAU*dH ~ 1e-1, aggregate over 7392 cells) matches
|diff|L2 ≈ 2.7e-06. Hence:

- **Mechanism:** the 0.824 = J's P-row boundary dphi term
  (solveDiscreteFlowAdjoint.H L689-697) present on the J side, absent on the
  NS.H rebuilt-phi FD side (boundary pinned, L676-681). It is a
  reference-convention artifact of F1. (This confirms the plan's Q-A framing:
  the FD reference pins boundary fluxes at baseline.)
- **But** J's boundary dphi is also NOT the production boundary tangent at the
  7 fixedValue-U patches (production: 0 there). So the "correct" oracle must
  decide which residual is the production state map's: R_P = div(phi_SIMPLE)
  with the production phi (constrainHbyA-pinned boundary). Under that
  semantics, J's P-boundary rows are expected to STILL be wrong at the 7
  fixedValue-U patch cells (J models a boundary dphi production does not
  have), i.e. P-boundary O(1) would survive a production-consistent oracle —
  localizing the genuine J boundary-term gap (candidate
  PATCH_P_BOUNDARY_JACOBIAN). Under the extrapolated-unrelaxed convention
  (F2), J matches by construction (GatePR 1e-10) and the 0.824 is fully a
  reference artifact (GENERAL_P_ROW_JACOBIAN_DEFECT = NOT_SUPPORTED). **Which
  of these two outcomes Gate A produces is the discriminator; this audit only
  maps the alternatives.**

---

## 8. Q-B note: the 2.5 ratio and the relax-source term

- `rAUAdj/rAU = 2.5000000087` (verified, §5.3) — exact 1/U-relax for the
  mobility. Confirms D2/D3/D4 facts at this HEAD.
- `discretePrimalMomentum` is built WITHOUT relax (L22-28) and its H is
  without the relax source; the production UEqn is relaxed (L45) and its H()
  carries the relax source `1.5*D0*psi/V`. The full production dphi/dU is
  `flux(d(0.4*rAUu*H_rel)/dU)` = `flux(0.4*rAUu*dHu/dU + 0.6*dU)` — the
  0.6*dU relax-source term partially compensates the 0.4 factor. Candidate B
  (relaxed mobility) must therefore be built with the FULL production
  sequence (relaxed H including the relax source), exactly as the Gate-B spec
  in the task file §6 requires ("equation relaxation ... rAU/rAtU ... HbyA
  ... phiHbyA ... pEqn.flux()"). A naive "0.4×rAUAdj" candidate that reuses
  the unrelaxed dH would measure the wrong object.
- The momentum residual R_U = A*U − source is relaxation-INVARIANT
  (S_relax = S0 + (D_relax−D0)*U exactly cancels D_relax*U − D0*U), so the
  momentum row's unrelaxed J*v-vs-FD closure (1.65e-4) is valid for the true
  residual regardless of relaxation; only the P-row's phi tangent is
  relaxation-sensitive. (This is why Q-B is a P-row/phi question, not a
  momentum question.)

---

## 9. Confirmation status of D2/D3/D4 facts at ca8a772b

| Claim (D2/D3/D4) | Status at ca8a772b |
|---|---|
| NS.H L45 relax before L113 rAU=1/A() (relaxed) | CONFIRMED (lines identical) |
| rAUAdj unrelaxed = 1/discretePrimalMomentum.A() (L76-81); matrix unrelaxed (L22-28) | CONFIRMED |
| J dphi mixes unrelaxed rAUAdj (dH) + relaxed primalPressureMobility (dp) (L625-641) | CONFIRMED |
| rAUAdj/rAU = 2.5 = 1/0.4 | CONFIRMED numerically (2.5000000087) |
| |phiRb0−phi| = 3.5e-11 ⇒ stored phi is the RELAXED rebuild | CONFIRMED |
| deltaA |dA_FD| = 0 | CONFIRMED |
| NS.H rebuilt FD: fresh unrelaxed matrix (L634-639), boundary PINNED (L676-681) | CONFIRMED |
| GatePR/GateOracle: unrelaxed rAUc (L1467), fresh tEqR (L1477-1482), rebuilt boundary flux (L1497-1502), no constrainHbyA/adjustPhi | CONFIRMED |
| P-row error 100% boundary-concentrated (7392 cells) | CONFIRMED (localize_prow.out) |
| Rx-B self-referential, A0 = UEqn.A()−alpha with RELAXED UEqn, boundary pinned (NS.H L243-244, L341-349, L388-391) | CONFIRMED |
| Momentum row 1.65e-4 / cos 0.99999999 | CONFIRMED (crosscheck_bf2.out) |
| constrainHbyA semantics (HbyA_b = U_b at non-assignable U patches) | NEW, verified from /opt/openfoam7 source |
| adjustPhi / setReference / constrainPressure are NO-OPs in this case | NEW, verified from /opt/openfoam7 source + case patch types |
| relax-source term in H() (1.5*D0*psi/V) and its 0.6*dU tangent contribution | NEW analysis (§3.3, §8) — untested by every existing probe |
| 0/U inlet fixedValue (193.97 0 0) / outlet zeroGradient (opposite of the plan text) | NEW correction (§2.2) |

---

## 10. Feed-forward to S2-S4 (oracle requirements this audit implies)

1. Gate-A FD must rebuild phi through the production sequence with the
   boundary flux REBUILT (not pinned), and must REPORT the split
   internal / fixedValue-U-patch cells / outlet-patch cells so the two
   candidate outcomes (§7) are separable. If the FD uses constrainHbyA
   (production), J's boundary dphi at the 7 fixedValue-U patches is the
   expected O(1) residual; if it uses the extrapolated convention (GatePR),
   the mismatch should collapse.
2. Gate-B candidates must use the FULL production semantics (relaxed H with
   the relax source + relaxed rAU + constrainHbyA + pEqn.flux with frozen p);
   a "relaxed mobility × unrelaxed dH" shortcut is NOT the production phi
   tangent (§8).
3. Existing anchors to re-verify in the oracle run: momentum relL2 ~ 1.65e-4,
   GatePR ~ 1e-10, |phiRb0 − phi| ~ 3.5e-11, deltaA = 0.

---

## 11. Files/lines consulted (read-only)

- src/NS.H: L26-45 (UEqn + dev2 + relax), L72-82 (IC stats), L113-129
  (rAU/HbyA/phiHbyA/adjustPhi/rAtU/primalPressureMobility), L136-410 (Rx-A/B),
  L465-577 (frozen-phi FD), L579-716 (rebuilt-phi FD: L591 relaxed baseline,
  L634-640 fresh unrelaxed rAUrb, L654-669 constrainHbyA/adjustPhi/pEqn/
  phiRb, L676-681 boundary PINNED), L1014-1041 (pEqn/phi/U update).
- src/solveDiscreteFlowAdjoint.H: L22-28 (discretePrimalMomentum unrelaxed),
  L30-32, L71-81 (rAUAdj), L250-304 (applyPressureFluxCorrection), L414-711
  (applyDiscreteFlowJ: dH L502-582, SIMPLE dphi L625-641, velocityJump
  L646-653, continuity L656-657, boundary L675-710), L1444-1744 (GatePR/
  GateOracle: rAUc L1467, tEqR L1477-1482, phiF L1486-1488, boundary rebuilt
  L1497-1502, face dphi export L1508-1513, GateDHsplit L1561-1699, GatePRint/
  GatePR L1700-1742), L1971-2039 (T2), L2271-2485 (explicit J^T COO),
  L2706-2780 (CSR + pRef pin), L2889-2907 (explicitJT export).
- src/readTransportProperties.H L66-68 (pRefCell via setRefCell).
- /opt/openfoam7: fvMatrix.C (relax ~L672, A() L738, H() L760+, flux() L863,
  setReference L505), fvMatrix.H L348-353; GeometricField.C L1100
  (needReference); cfdTools/general/{constrainHbyA,adjustPhi,constrainPressure}
  .C; fvcFlux.C (scheme dotInterpolate); fixedValueFvPatchField.H L148-165
  (fixesValue/assignable), fvPatchField.H L295-320 (defaults).
- Case: system/fvSolution, system/fvSchemes, system/controlDict, 0/U, 0/p,
  constant/optProperties, constant/polyMesh/boundary; no constant/fvOptions.
- Evidence: evidence/agent-group/BFINAL-001/cycle-2/{gate_spec.md, audit.md,
  rx_audit.md, crosscheck_bf2.out, localize_prow.out, crosscheck_bf2.py,
  localize_prow.py, rerun.log}.
