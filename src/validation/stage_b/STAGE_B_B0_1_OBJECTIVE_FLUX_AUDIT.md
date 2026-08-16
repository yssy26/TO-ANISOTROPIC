# Stage B0.1 — Mass-Flow-Weighted Outlet Objective Flux-Derivative Audit and Minimal Correction

- Repository: `/home/ys/TO-ANISOTROPIC`
- Branch: `agent/stage-a-anisotropic-validation`
- Baseline commit (audit start): `f3c4aa25e4ae46b4bed80143cc1ef81ab6c7c68f`
- OpenFOAM: OpenFOAM-7, gate4 case (`/home/ys/Twostream/validation/gate4`), frozen state at time `1`
- Audit date: 2026-08-06
- Scope: Stage B0.1 only. No B1, no MMA, no `gradientValidated=true`, no full SST Gate 6, no change to the objective/constraint/adjoint definitions, no FD-threshold relaxation.

## 1. Objective definition and complete differential

```text
M    = sum_f(phi_f)
N    = sum_f(phi_f*T_f)
Tmix = N/M
J    = -Tmix/Tref
```

Complete first variation:

```text
delta J = sum_f (dJ/dT_f)*delta T_f  +  sum_f (dJ/dphi_f)*delta phi_f
dJ/dT_f   = -phi_f/(M*Tref)
dJ/dphi_f = -(T_f - Tmix)/(M*Tref)
```

Sign conventions: phi is the cold-stream volumetric flux, positive out of the
domain at the outlet; M > 0 for outflow. If `sum_f delta phi_f = 0` (global mass
fixed), the term `sum_f -(T_f - Tmix)/(M*Tref) delta phi_f` is generally NOT
zero because the weights `-(T_f - Tmix)` are non-uniform when the outlet
temperature is not uniform (`Tmix` is the flux-weighted mean, so the weights
sum to zero only when weighted by phi, not by the raw perturbation). Special
conditions where `dJ/dphi` vanishes: all outlet T identical (then
`T_f - Tmix = 0`); a single independent outlet flux degree of freedom with the
total flux fixed (then the only allowed perturbation is zero); or every outlet
face flux strictly fixed by the BCs. **None holds for the gate4 case**: outlet T
spans 606-688 K and the outlet flux has 84 independent calculated-patch faces
(pinned only by U/p through the flow solve), so `dJ/dphi` is a genuine,
O(0.1) contribution that must flow into the discrete flow adjoint.

## 2. Current boundary conditions (gate4)

| Field | inlet | outlet |
|---|---|---|
| U | fixedValue (193.97,0,0) | zeroGradient |
| p | zeroGradient | fixedValue 0 |
| T | fixedValue 600 | zeroGradient |
| phi | calculated | calculated |

`phiThermal = coldFaceMask*phi + hotRhoCpRatio*hotFaceMask*phiHotFrozen`
(updateThermalFlux.H); on the cold outlet `coldFaceMask = 1`,
`hotFaceMask = 0`, so `phiThermal[outlet] = phi[outlet]` and
`d phiThermal/d phi = 1` there. The outlet face flux DOES depend on the design
through the flow solve (U outlet is zeroGradient, p outlet fixedValue, so the
outlet velocity/flux is a free response), and even with the global outlet mass
flow fixed the local flux distribution varies with the design.

## 3. dJ/dT: implemented and verified

`computeObjective.H:26-56` sets, for `maximizeColdOutletTemperature`,

```cpp
thermalObjectiveDerivative[outletCells[facei]] -= phiOutlet[facei]
    /max(coldOutletMassFlow*thermalReferenceTemperature, SMALL);
```

i.e. `dJ/dT_f = -phi_f/(M*Tref)`. Algebraic FD (per outlet face, 3 steps):
max relative error = 9.54e-9 <= 1e-8, sign consistent. The initialization-time
area-weighted variant (createFrozenHotRegionFields.H:584-601) is used only by
validateDiscreteObjectiveDerivatives.H at start-up; the runtime adjoint always
uses the phi-weighted version refreshed each outer iteration (MTO_HF.C:54).

## 4. dJ/dphi: MISSING before the B0.1 fix (confirmed)

No code path constructed or wrote `-(T_f - Tmix)/(M*Tref)`. The adjoint-source
chain for phi (AdjNS_HT.H discrete branch) filled `discreteExternalFaceFluxAdjoint`
only for INTERNAL cold faces (`-coldThermalFaceMask*adjointDownwind*
temperatureJump`, AdjNS_HT.H:11-26); the cold outlet boundary kept its initial
zero. Algebraic FD of the exact formula on the frozen gate4 outlet state (3
directions x {1e-3,1e-4,1e-5,1e-6}, Kahan summation, L2-normalized directions):

| direction | D_an | median rel err | max rel err | sign mismatches |
|---|---:|---:|---:|---:|
| d1 non-conserving | 6.40e-2 | - | 2.88e-9 | 0 |
| d2 conserving (zero sum) | -1.05e-1 | - | - | 0 |
| d3 T-correlated zero-sum | -4.75e-1 | - | - | 0 |
| overall | - | 1.65e-11 | 2.88e-9 | 0 |

Gate: median <= 1e-9, max <= 1e-8 -> PASS. Magnitude O(0.1), not negligible.

## 5. Thermal-residual flux derivative <Tb, dR_T/dphi>

The discrete residual is the SOLVED operator

```cpp
R_T = fvm::div(phiThermal,T,"div(phi,T)") - fvm::laplacian(DTEffective,T) - Q
```

Verified with a standalone OpenFOAM tool (`heatFluxDerivTest`, temporary
directory `/home/ys/heatflux_deriv_test`, NOT committed) that re-assembles
`M.residual()` (exact solver semantics, `fvMatrix::residual()` = source - A*T)
for perturbed cold flux `phi -> coldFaceMask*phiPert + hotPart`, and compares
`<Tb, r(phi+h d)>` centered differences against candidate analytic forms on
the frozen state.

Single-face results (h=1e-7):

| face | D_fd = <Tb, dr>/dphi | code internal expr | Tf*(Tb_o-Tb_n) |
|---|---:|---:|---:|
| internal 16330 (own 5600, nei 5601) | 1.83e-14 | 1.47e-4 | 0.109 |
| outlet face 0 (cell 5679) | 0 (dr vector identically zero, count=0) | 0 (not filled) | -248.3 |
| outlet faces 0-5 | 0 | - | -248 |
| inlet faces 0-5 | -2.0e-14 | - | -211 |

Key findings:

1. **Outlet boundary residual-flux term is EXACTLY zero** under the fvm
   discretization: perturbing the outflow zeroGradient-T outlet flux leaves
   the residual vector unchanged (`dr == 0`). The implicit upwind handling of
   the outflow face makes `dR_T/dphi_out = 0`. Therefore the code's zero
   outlet-boundary value is CORRECT; there is no missing boundary residual term.
2. **Internal-face residual-flux term is numerically negligible**: `dr` has a
   zero-sum "temperature-difference" structure
   (`dr[own] ~ T_own - T_nei`, `dr[nei] ~ T_nei - T_own`), so
   `<Tb, dr> ~ (T_nei - T_own)*(Tb_own - Tb_nei)*V ~ 1e-14` for this
   converged state, versus the objective term O(0.1).
3. The code's internal-face expression `-Tb_downwind*(T_nei - T_own)`
   (magnitude 1.5e-4) does NOT equal the exact fvm residual derivative
   (1.8e-14); it matches the fvc-explicit divergence derivative instead. This
   is a pre-existing inaccuracy in the discrete flow-adjoint thermal-coupling
   term, **never activated by any Stage A/B test** (`solveFlowAdjoints=false`
   in every validation case). It is listed as a finding and a B1 precondition;
   it is NOT part of the dJ/dphi gap fixed here.
4. Inlet faces: residual-flux term ~1e-14 (negligible); the inlet flux is
   design-fixed (U fixedValue) and the flow-adjoint boundary loop skips
   fixedValue-U patches, consistent.

Consequence: the combined flow-adjoint boundary flux RHS on the cold outlet is
exactly the objective direct term `dJ/dphi`; the thermal-residual boundary term
contributes nothing. **Case B of the task's decision tree applies** (dJ/dphi
missing, boundary residual term complete/zero), plus the internal-face
inaccuracy finding.

## 6. Flow-adjoint interface path (discreteExternalFaceFluxAdjoint)

`solveDiscreteFlowAdjoint.H`:
- internal faces: `discreteExternalFluxTranspose[U(own)] += weights*Sf*faceAdj`,
  `[U(nei)] += (1-weights)*Sf*faceAdj` (lines 398-412),
- boundary faces with U not fixedValue (outlet is zeroGradient):
  `discreteExternalFluxTranspose[U(celli)] += Sf*faceAdj[patchStart+facei]`
  (lines 413-434),
- pressure path: `applyPressureFluxCorrectionTranspose(faceAdjoint, ...)`
  adds the same face adjoint to the pressure RHS through the pressure-flux
  correction chain (lines 435-449, 301-397).

The discrete flow adjoint solves `A^T lambda = RHS` with
`RHS[U] = faceAdjoint^T (dphi/dU)` and `RHS[p]` from the pressure-flux
correction; the gradient then uses `-dAlphaDxh*(U & Uadj)*V` (sensitivity.H).
The Lagrangian derivation gives `faceAdjoint = dJ/dphi + <Tb, dR_T/dphi>` with
the POSITIVE sign as the `dJ/du` RHS (since `dJ/dx = dJ/dx_direct -
lambda^T dR/dx` and `A^T lambda = dJ/du`). With the boundary residual term
zero, the cold outlet face adjoint must be exactly `dJ/dphi_f =
-(T_f - Tmix)/(M*Tref)`.

## 7. Minimal correction applied (Case B)

Three files changed, 80 insertions, no deletions, no changes to the objective,
constraints, adjoint equations, gradient formulas, thresholds or MMA:

1. `src/createFrozenHotRegionFields.H`: declare `thermalObjectiveFluxDerivative`
   (surfaceScalarField, dimless, calculated patches, NO_READ/NO_WRITE) next to
   `thermalObjectiveDerivative`.
2. `src/computeObjective.H`: each outer iteration reset the field to zero and,
   for `maximizeColdOutletTemperature`, set the cold outlet boundary to
   `-(T_face - Tmix)/(M*Tref)` using the same M/N reduction as costfunction.H.
3. `src/AdjNS_HT.H`: in the discrete branch, after the internal-face loop, add
   the cold outlet boundary values of `thermalObjectiveFluxDerivative` into
   `discreteExternalFaceFluxAdjoint` (global face index
   `patchStart + facei`).

The internal-face thermal-coupling term was intentionally left unchanged (it is
a separate, pre-existing inaccuracy; fixing it is a B1 precondition, see
section 11). The pressure-drop adjoint was not touched. F1 (thermal-only)
gradients are unaffected because F1 fixes phi and does not solve Ub.

Clean rebuild on the working tree (wmake, 0 errors): log
`src/validation_stage_b01_build.log`.

## 8. Verification of the fix

1. dJ/dT algebraic FD: max 9.54e-9 (PASS).
2. dJ/dphi algebraic FD: median 1.65e-11, max 2.88e-9, sign 0 (PASS).
3. Internal-face residual-flux FD: <Tb,dr> = 1.8e-14 (numerically negligible).
4. Outlet-boundary residual-flux FD: dr identically zero (term is exactly 0).
5. Combined T-phi algebraic FD (3 directions, 4 steps):
   median 3.12e-11, max over h<=1e-4 = 6.28e-10 (PASS); the h=1e-3 point shows
   the expected second-order truncation (h^2 convergence, 6.3e-8 at 1e-3 ->
   6.3e-10 at 1e-4).
6. The complete objective flux term now enters the discrete flow-adjoint RHS:
   `discreteExternalFaceFluxAdjoint[outlet] = dJ/dphi` and the boundary loop
   feeds it to both U and p RHS (section 6). (End-to-end magnitude validation
   of the Ub chain is Stage B1 scope.)
7-8. F1 regression serial and 4-rank parallel: see below.
9. No sign mismatches in any significant direction.
10. `objectiveGradientScale` / `pressureGradientScale` unchanged (1.0 / 0.3,
    and the 1.0/1.0 recommendation of B0 for amplitude FD is untouched).
11. Objective definition unchanged.
12. All tests ran on copies (hf_case, f1_diag_serial/parallel); the original
    gate4 case was not modified; the frozen state was restored after every
    perturbation.

## 9. F1 regression (after the fix)

| metric | serial | 4-rank parallel | requirement |
|---|---:|---:|---:|
| medianRelError | 2.25e-4 | 9.95e-5 | <= 1% |
| maxRelError | 2.49e-3 | 1.51e-3 | <= 5% |
| failures | 0/9 | 0/9 | 0/9 |
| sign mismatch | 0 | 0 | 0 |
| D_ADJ dir 0 | -4.056556e-3 | -4.056556e-3 | - |
| D_ADJ dir 1 | 1.731376e-4 | 1.731376e-4 | - |
| D_ADJ dir 2 | 5.671516e-5 | 5.671517e-5 | - |
| dir 1,2 max serial/parallel rel diff | 3.65e-8 | - | <= 1e-6 |

D_ADJ values are identical to the pre-fix run (bit-for-bit), confirming that the
added dJ/dphi term does not enter the F1 thermal-only gradient (F1 fixes phi).

## 10. Combined derivative chain after the fix

Objective gradient now contains both paths:

```text
Path T:   x -> xp -> xh -> alpha, DTMolecular -> U,p,phi,T -> J -> dJ/dT
          -> Tb (adjoint) -> dJ/dxh (thermal diffusion + Brinkman) -> dfdx
Path phi: x -> xp -> xh -> alpha -> U,p,phi -> J -> dJ/dphi (new, cold outlet)
          -> Ub (flow adjoint RHS) -> dJ/dxh (Brinkman) -> dfdx
Coupling: x -> U,p,phi -> R_T -> T -> J (through Tb) and the thermal residual
          flux terms (numerically zero at outlet, negligible internally)
```

Ub now corresponds to the full `dJ/d(U,p)` (objective direct phi path + thermal
equation coupling), not only the thermal-equation indirect coupling. The
pressure-drop chain is untouched.

## 11. Findings and recommended follow-up (not applied)

1. [FIXED] `dJ/dphi` was missing; implemented and verified (sections 4, 7).
2. [PRE-EXISTING, NOT FIXED] The internal-face thermal-coupling term
   `-Tb_downwind*(T_nei - T_own)` in AdjNS_HT.H does not equal the exact fvm
   residual derivative `<Tb, dR_T/dphi_f>` (1.8e-14 measured vs 1.5e-4 coded).
   It matches the fvc-explicit divergence derivative. This path is inactive in
   every Stage A/B validation case (`solveFlowAdjoints=false`), and its
   magnitude (1.5e-4) is small relative to dJ/dphi (O(0.1)), but it should be
   corrected or re-derived before Stage B1 activates the discrete flow adjoint.
3. [INFO] `fvMatrix::residual()` must be used for residual-flux derivative
   checks; `fvc::div` gives a different (fvc-explicit) phi-derivative and is
   NOT the operator the discrete adjoint transposes.
4. [INFO] The `1/Tb` field read at time `1` does not satisfy
   `A0^T Tb = g` (relative residual ~2.4 by our manual assembly); the reason is
   that the gate4 `1/` state was produced with `solveFlowAdjoints=false`
   (Tb written as zero-field fill by the frozen-cold-flow branch). This does
   not affect the B0.1 conclusions: the residual-flux derivative measurements
   used the actual Tb field and the conclusions (dr structure, dr=0 at outlet)
   are independent of Tb.

## 12. B0.1 acceptance

All 12 acceptance criteria of the task are satisfied: (1) dJ/dT algebraic FD,
(2) dJ/dphi algebraic FD, (3) internal-face residual-flux discrete FD,
(4) outlet-boundary residual-flux discrete FD, (5) combined T-phi algebraic FD,
(6) objective flux term enters the discrete flow-adjoint RHS,
(7) F1 serial regression, (8) F1 parallel regression, (9) no sign mismatches,
(10) no gradient-scale adjustment, (11) objective definition unchanged,
(12) state fully restored.

**Final status: PASS** (Case B fix + two documented pre-existing findings).

## 13. Deliverables

- `/home/ys/TO-ANISOTROPIC/src/validation/stage_b/STAGE_B_B0_1_OBJECTIVE_FLUX_AUDIT.md`
- `/home/ys/TO-ANISOTROPIC/src/validation/stage_b/stage_b_b0_1_objective_flux_audit.json`
- `/home/ys/TO-ANISOTROPIC/src/validation/stage_b/stage_b_b0_1_direct_derivative.tsv`

Updated (B0 files): `STAGE_B_B0_AUDIT_REPORT.md`, `stage_b_b0_audit.json`,
`stage_b_b0_variable_trace.tsv` (B0-01/B0-17/B0-20/readyForB1/overallStatus).

## 14. B0 closure and B1 readiness

B0 can be closed with the B0.1 corrections included. Stage B1 formal
finite-difference validation must additionally:
- keep `objectiveGradientScale 1.0; pressureGradientScale 1.0;`,
- run serial (discrete flow adjoint is serial-only),
- use the flow-coupled configuration (`solveFlowAdjoints=true`,
  `freezeColdFlowForValidation=false`) and compare `dgdx[1]` against the
  normalized `gDP`,
- re-derive or remove the internal-face `-Tb_downwind*(T_nei-T_own)` term so
  the Ub chain matches the exact fvm residual derivatives.

`gradientValidated` remains `false`; no optimization was started.
