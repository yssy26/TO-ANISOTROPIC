# Stage B0 — Frozen-Turbulence Gradient Definition and Scaling Audit

## 1. Repository, branch, commit

| Item | Value |
|---|---|
| Repository | `/home/ys/TO-ANISOTROPIC` |
| Branch | `agent/stage-a-anisotropic-validation` |
| HEAD | `f3c4aa25e4ae46b4bed80143cc1ef81ab6c7c68f` |
| OpenFOAM | OpenFOAM-7 |
| Audit date | 2026-08-06 |
| Audit scope | Stage B0 only: definitions, normalization, scaling, array mapping of J / gV / gDP. No B1, no formal FD, no MMA, no SST Gate 6, no physics changes. |

## 2. Working-tree state

`git status --short` at audit start (unchanged after audit):

```
?? src/validation_stage_a_diagnostic_cleanup_build.log
?? src/validation_stage_a_processor_unit.log
?? src/validation_stage_a_rebuild.log
?? src/validation_stage_a_unit.log
```

Only the four Stage-A validation log files are untracked. They pre-date this audit, were not created by it, and were not modified, staged, committed, or pushed. This audit created three new files under `src/validation/stage_b/` (see section 15); they remain untracked by design (no commit per task instructions).

## 3. Objective J — actual code definition

```cpp
// costfunction.H (runtime value)
MeanT = sum(phi_out*T_out)/sum(phi_out);                       // Tmix, mass-flow weighted
normalizedThermalObjective = -MeanT/thermalReferenceTemperature; // minimization form
Objective = normalizedThermalObjective;                          // objectiveMode == "outletTemperature"
// or Objective = (1-w1)*normalizedThermalObjective + w1*PowerDiss (legacyComposite)
```

- Value: `J = -Tmix / Tref`, `Tref = thermalReferenceTemperature` (`createFrozenHotRegionFields.H:119-124`).
- Direct state derivative: `computeObjective.H:47-56` gives `dJ/dT[celli] = -phi_out[face]/(sum(phi_out)*Tref)` for cold-outlet-adjacent cells (zero elsewhere). `validateCommon.H` `evaluateObjective()` uses the identical formula `-meanT/Tref` (used by F1/F2/F3/Gate6 FD).
- Adjoint source sign: `AdjHeatTransfer.H:30` `thermalTransposeEqn.source() = thermalObjectiveDerivative` — the discrete transpose is solved with the dJ/dT source carrying the negative sign of the minimization form. Consistent.
- Normalized: YES, divided by `Tref` exactly once.
- `dfdx` maps to the normalized objective: `dfdx = objectiveGradientScale * fsensMeanT` (`sensitivity.H:122-127`), where `fsensMeanT` is the filtered per-cell dJ/dx of `J = -Tmix/Tref`.
- `objectiveGradientScale`: a single user-supplied positive multiplier, default 1.0 (`createFields.H:207-210`). With scale=1, `dfdx = raw dJ/dx`.

## 4. Volume constraint gV — actual code definition

```cpp
// costfunction.H:49-53
designFluidFraction = domainIntegrate(xh*designMask)/designVolume;   // = sum(xh*mask*V)/designVolume
designSolidFraction = 1.0 - designFluidFraction;
V = solidFractionMin - designSolidFraction;                          // gV
```

- `gx[0] = V` (`sensitivity.H:115`).
- `dgdx[0] = gsensVol[celli]` with NO scale (`sensitivity.H:128`).
- `gsenshVol = designMask*mesh.V()/designVolume` (`filter_chainrule.H:6-7`) — the direct derivative of gV w.r.t. xh, containing `mesh.V()/designVolume` exactly once (domainIntegrate carries the cell volume).
- Projection chain rule and eta5 implicit derivative are included (`filter_chainrule.H:33-187`): `gsenshVol *= drho` then `+= designMask*V_cell*(1-drho)*volumeEtaCorrection` with `dProjectionDeta` and the volume-preserving `eta5` correction. E5 (`validateGradientChain.H`) re-computes eta5 for every FD perturbation and validates this.
- Helmholtz filter adjoint output `gsensVol` is per-cell; NO additional mesh.V() (E4 dot-product test in `validateGradientChain.H:190-266` uses consistent V-weighting; filter_chainrule comments 221-227 state this explicitly).
- Fixed regions: zeroed by `designMask` multiplication (`filter_chainrule.H:10-13, 216-232`), by `dDTDxh` zeroing (`updateMaterialProperties.H:68-71`) and `dAlphaDxh` zeroing (`sensitivity.H:29-34`); E6 checks.

## 5. Pressure-drop constraint gDP — actual code definition

```cpp
// costfunction.H:88-118
pInletAvg  = sum(p_in*area)/sum(area);      // area-weighted
pOutletAvg = sum(p_out*area)/sum(area);
PressureDropPa = rhoFluid*(pInletAvg - pOutletAvg);
pressureDropConstraint = PressureDropPa/pressureDropMaxPa - 1.0;    // gDP
```

- `gx[1] = PressureDropPa/pressureDropMaxPa - 1.0` (`sensitivity.H:116`).
- `dgdx[1] = pressureGradientScale*gsensPressureDrop[celli]` (`sensitivity.H:129-130`).
- The direct derivative is `pressureConstraintDerivative` (`createFrozenHotRegionFields.H:603-641`):
  `+rhoFluid*A_face/(pressureDropMaxPa*A_inletTotal)` at non-fixed inlet cells, `-rhoFluid*A_face/(pressureDropMaxPa*A_outletTotal)` at non-fixed outlet cells. It is d(gDP)/dp — the `rhoFluid/pressureDropMaxPa` normalization is already inside.
- This feeds the pressure adjoint RHS: `AdjNS_PD.H:4` `discreteFlowPressureRhs(pressureConstraintDerivative)` → `solveDiscreteFlowAdjoint.H` (J^T solve) → `Uc`. Then `gsenshPressureDrop = -dAlphaDxh*(U & Uc)*mesh.V()` (`sensitivity.H:46-47`) → projection/filter chain → `gsensPressureDrop`.
- Validation FD compares against the NORMALIZED constraint: F2 (`validateFrozenGradient.H:902-911`) and Gate6-strict (`validateGate6SSTStrict.H:425-449`) both compute `(DPplus-DPminus)/(2h)/pressureDropMaxPa`; Gate6 comments (lines 6-7, 399) explicitly state `dgdx[1]` is not divided by pressureDropMaxPa again. So Stage B FD must compare `dgdx[1]` against `d(gDP)/dx` (not `d(PressureDropPa)/dx`).
- `pressureGradientScale`: single user multiplier, default 1.0 (`createFields.H:211-214`). With scale=1, `dgdx[1] = raw dgDP/dx`.

## 6. Three gradient chains (full trace)

### 6.1 Objective chain

| # | Step | File | Variable | Cell-quantity? | mesh.V()? | Normalized? | Mask? | eta5? | User scale? |
|---|---|---|---|---|---|---|---|---|---|
| 1 | design variable | opt_initialization.H / createFields.H | x | per-cell | no | no | (fixed enforced) | - | - |
| 2 | Helmholtz filter | filter_x.H:7-13; validateCommon.H:80-87 | xp | per-cell | matrix form (no explicit V) | no | yes (designFilterFaceMask + enforceFixed) | - | - |
| 3 | Heaviside projection | filter_x.H:62-86 | xh | per-cell | no | no | yes | eta5 adaptive | - |
| 4 | Brinkman penalty | updateMaterialProperties.H:26-30 | alpha | per-cell | no | no | yes | - | - |
| 5 | thermal diffusivity | updateMaterialProperties.H:17-23 | DTMolecular | per-cell tensor | no | no | yes | - | - |
| 6 | effective diffusivity | updateFrozenTurbulenceFields.H:71 | DTEffective | per-cell tensor | no | no | - | - | - |
| 7 | primal fields | NS.H, HeatTransfer.H | U,p,phi,T | - | - | - | - | - | - |
| 8 | objective value | costfunction.H:38-40,71-79 | J = -Tmix/Tref | global | no | YES (÷Tref) | - | - | - |
| 9 | dJ/dT | computeObjective.H:47-56 | thermalObjectiveDerivative | per-cell | no | YES (÷Tref) | - | - | - |
| 10 | thermal adjoint | AdjHeatTransfer.H:26-126 | Tb | per-cell | matrix | (source normalized) | - | - | - |
| 11 | diffusion derivative | AdjHeatTransfer.H:189-328 | thermalDiffusionDerivativeDTCell | per-cell face-sum | NO (already per-cell) | - | yes (dDTDxh) | - | - |
| 12 | Brinkman derivative | sensitivity.H:44-47,82-85 | -dAlphaDxh*(U&Ub)*V | per-cell | YES (once) | - | yes (dAlphaDxh) | - | - |
| 13 | projection chain rule | filter_chainrule.H:33-187 | drho, dProjectionDeta, etaCorr | per-cell | V only inside etaCorr factor | - | yes | YES (implicit) | - |
| 14 | Helmholtz filter adjoint | filter_chainrule.H:190-213 | fsensMeanT | per-cell | NO extra | - | yes (operator) | - | - |
| 15 | post-filter mask | filter_chainrule.H:216-232 | fsensMeanT | per-cell | - | - | yes | - | - |
| 16 | MMA gradient | sensitivity.H:122-127 | dfdx | per design cell | - | - | - | - | YES (objectiveGradientScale) |

### 6.2 Volume chain

| # | Step | File | Variable | V? | Normalized? | Mask? | eta5? | Scale? |
|---|---|---|---|---|---|---|---|---|
| 1-3 | x→xp→xh | filter_x.H | (as above) | - | - | yes | eta5 | - |
| 4 | solid fraction | costfunction.H:49-52 | designSolidFraction = 1 - ∫xh·mask/V_design | V in integral | yes (÷designVolume) | yes | - | - |
| 5 | constraint value | costfunction.H:53; sensitivity.H:115 | gV; gx[0] | - | yes | - | - | - |
| 6 | direct derivative | filter_chainrule.H:6-7 | gsenshVol = mask·V_cell/designVolume | YES (once) | yes (÷designVolume) | yes | - | - |
| 7 | projection chain rule | filter_chainrule.H:148,157 | ×drho, +eta correction | V only in correction | - | yes | YES | - |
| 8 | filter adjoint | filter_chainrule.H:210-213 | gsensVol | NO extra | - | operator | - | - |
| 9 | post-filter mask | filter_chainrule.H:218,231 | gsensVol | - | - | yes | - | - |
| 10 | MMA gradient | sensitivity.H:128 | dgdx[0] = gsensVol | - | - | - | - | NONE |

### 6.3 Pressure chain

| # | Step | File | Variable | V? | Normalized? | Mask? | eta5? | Scale? |
|---|---|---|---|---|---|---|---|---|
| 1-3 | x→xp→xh | filter_x.H | (as above) | - | - | yes | eta5 | - |
| 4 | pressure drop | costfunction.H:88-116 | PressureDropPa = rho·(pInAvg-pOutAvg) | - | no (raw Pa) | - | - | - |
| 5 | constraint value | costfunction.H:117-118; sensitivity.H:116 | gDP; gx[1] | - | YES (÷pressureDropMaxPa) | - | - | - |
| 6 | direct derivative | createFrozenHotRegionFields.H:603-641 | pressureConstraintDerivative (±rho·A/(pMax·A_tot)) | - | YES (rho/pMax inside) | - | - | - |
| 7 | pressure adjoint | AdjNS_PD.H:4; solveDiscreteFlowAdjoint.H | Uc (J^T solve) | - | - | - | - | - |
| 8 | Darcy derivative | sensitivity.H:46-47 | gsenshPressureDrop = -dAlphaDxh·(U&Uc)·V | YES (once) | (normalized source) | yes (dAlphaDxh) | - | - |
| 9 | projection + filter | filter_chainrule.H | gsensPressureDrop | NO extra | - | yes | YES | - |
| 10 | MMA gradient | sensitivity.H:129-130 | dgdx[1] | - | - | - | - | YES (pressureGradientScale) |

## 7. Normalization audit

- J: divided by `Tref` exactly once (value and derivative). No duplicate division found.
- gDP: divided by `pressureDropMaxPa` exactly once. The derivative contains `rhoFluid/pressureDropMaxPa` once (inside `pressureConstraintDerivative`). F2 and Gate6 FD normalize the raw-pressure FD once. No duplicate normalization anywhere (`validateGate6SSTStrict.H:6-7,399` explicitly documents this).
- gV: divided by `designVolume` once, in both value and direct derivative.
- `PowerDiss` is normalized by `PowerDiss0` once (legacyComposite objective only; not part of the Stage B outletTemperature objective).
- No dynamic normalization found: searched for `gradient/=gMax`, `/=Foam::max(gMax`, `gradient.normalize`, `/=maxGrad`, `/=currentObjective`, `/=currentPressureDrop`, `GradientScale`-like dynamic patterns across `src/`. Only the two static positive multipliers `objectiveGradientScale` and `pressureGradientScale` exist.

## 8. Gradient-scaling audit

- Mathematical raw gradients: `fsensMeanT` (objective), `gsensVol` (volume), `gsensPressureDrop` (pressure) — these are what F1/F2/F3 FD should ideally compare against; F1 uses `fsensThermalOnly`/`dfdxThermal`.
- Gradients passed to MMA: `dfdx = objectiveGradientScale*fsensMeanT`, `dgdx[0] = gsensVol` (no scale), `dgdx[1] = pressureGradientScale*gsensPressureDrop`.
- Gradients written to disk: raw fields `fsensMeanT`, `gsensVol`, `gsensPressureDrop`, `gsenshMeanT`, etc., plus `dfdxField` (= the SCALED dfdx) (`sensitivity.H:187-218`).
- The identities `dfdx = objectiveGradientScale*rawObjectiveGradient` and `dgdx[1] = pressureGradientScale*rawPressureGradient` hold exactly (single multiplication in `sensitivity.H:122,130`; verified by E7 `validateGradientChain.H:539-553`).
- **B1 formal amplitude validation must use `objectiveGradientScale 1.0; pressureGradientScale 1.0`** — this is the only safe configuration: F1/F2/F3/Gate6 compare `scale*gradient` against the FD of the normalized J/gDP. With scale≠1 the FD comparison would differ by the factor scale and fail (or be silently tuned). No suggestion is made to adjust the scales to make FD "pass".
- No `gradient /= max(...)` style dynamic normalization exists (see section 7).

## 9. mesh.V() audit

1. Darcy (Brinkman) term: `fvm::Sp(alpha,U)` assembles `alpha*V_cell` on the discrete diagonal (OpenFOAM Sp semantics), so `dJ/dxh` from the Darcy path must carry one cell volume. Code: `fsenshMeanT = -dAlphaDxh*(U&Ub); *= mesh.V()` (`sensitivity.H:44-47`) — exactly one `mesh.V()`. Same for `gsenshPressureDrop` and `gsenshPowerDiss`.
2. Thermal-diffusion term: `thermalDiffusionDerivativeDTCell` is assembled face-by-face as `(Tb_o-Tb_n)(T_o-T_n)*deltaCoeffs*magSf*e_f^T(dD/dxh)e_f` (`AdjHeatTransfer.H:189-328`). This is the derivative of the discrete Laplacian residual already integrated over the cell volume (coefficient `gamma*Sf/d`, no extra V factor needed); it is a per-cell quantity and is added WITHOUT mesh.V() (`sensitivity.H:69-72`). Correct.
3. Volume term: `gsenshVol = designMask*mesh.V()/designVolume` — exactly one `mesh.V()`, one `designVolume`. Correct.
4. Helmholtz filter adjoint output: used directly, no additional mesh.V() (`filter_chainrule.H:221-227` comments). Correct — the filter equation `(L - bV I) xp = -bV x` is symmetric and the E4 dot-product test (`validateGradientChain.H:246-265`) confirms consistency with V-weighted inner products.
5. Serial/parallel: the A7 fixes make the thermal adjoint (`AdjHeatTransfer.H` processor branch), the design-filter mask (`createFrozenHotRegionFields.H`) and the F1 directions partition-consistent. `solveDiscreteFlowAdjoint.H` explicitly aborts in parallel (`lines 1-8`) — the reduced discrete flow adjoint is serial-only by design. Both paths apply mesh.V() identically within their supported modes.

## 10. Frozen-turbulence scope audit

Fields frozen under `freezeTurbulenceForValidation=true`:
- `k`, `omega`, `nut` (RANS): not re-solved during design perturbations; `turbulence->correct()` is guarded by `!freezeTurbulenceForValidation` (`NS.H:106-109`, `validateCommon.H:246` only inside `solveFullSST` which sets the flag false).
- `nutFrozen`: retained (the `if (!freezeTurbulenceForValidation)` guards at `updateFrozenTurbulenceFields.H:13-21,51-55` skip recomputation in frozen mode). Fixed-region cells are set to 0 unconditionally (`lines 34-40`).
- `nuEffFrozen = nu + nutFrozen` (`line 60`): frozen.
- `alphaTurbulent = nutFrozen/turbulentPrandtl` (`line 61`): frozen.

Fields recomputed after a design perturbation:
- `x → xp → xh` (filter_x.H), `alpha`, `DTMolecular`, `dDTDxh` (updateMaterialProperties.H), `DTEffective = DTMolecular + alphaTurbulent*I` (`updateFrozenTurbulenceFields.H:71`) — DTEffective DOES update through the molecular part while the turbulent part stays frozen. This is the intended Stage B semantics: "fixed turbulent transport coefficients, design-dependent molecular part".

Perturbation paths (F1 `solveThermalOnly`, F2/F3 `solveFrozenPrimalFixed`, Gate6 `solveFullSST` with flag forced false) never call `turbulence->correct()` in the frozen validation loops. `evaluateCandidate.H` only enables turbulence correction for MMA acceptance evaluation (not used in B0/B1).

## 11. B0-01..B0-20 checklist

| # | Check | Result | Evidence |
|---|---|---|---|
| B0-01 | Objective value and direct derivative agree | PASS | costfunction.H:38-40 vs computeObjective.H:47-56 (both -Tmix/Tref form); validateDiscreteObjectiveDerivatives.H:21-55 |
| B0-02 | Pressure value and direct derivative agree | PASS | costfunction.H:116-118 vs createFrozenHotRegionFields.H:603-641; validateDiscreteObjectiveDerivatives.H:66-123 |
| B0-03 | Volume value and gradient agree | PASS | costfunction.H:49-53 vs filter_chainrule.H:6-7; E5 FD incl. eta5 recompute |
| B0-04 | gx[0], gx[1] ordering | PASS | sensitivity.H:115-116; MMA m=2 (opt_initialization.H:103-117) |
| B0-05 | dgdx[0], dgdx[1] ordering | PASS | sensitivity.H:128-130 |
| B0-06 | Pressure normalized exactly once | PASS | rho/pMax in pressureConstraintDerivative only; F2/Gate6 normalize FD once; Gate6 comments 6-7,399 |
| B0-07 | Thermal objective divided by Tref once | PASS | costfunction.H:73; computeObjective.H:53; no second division |
| B0-08 | objectiveGradientScale settable to 1 | PASS | createFields.H:207-210 default 1.0; required for B1 |
| B0-09 | pressureGradientScale settable to 1 | PASS | createFields.H:211-214 default 1.0; required for B1 |
| B0-10 | No hidden dynamic gradient normalization | PASS | full-tree search; only the two static scales |
| B0-11 | Darcy gradient has mesh.V() once | PASS | sensitivity.H:44-47,82-85 (Sp diagonal carries V) |
| B0-12 | Thermal-diffusion gradient no extra mesh.V() | PASS | sensitivity.H:69-72; AdjHeatTransfer.H per-cell face-sum |
| B0-13 | Volume gradient mesh.V()/designVolume once | PASS | filter_chainrule.H:6-7 |
| B0-14 | Filter adjoint output no extra mesh.V() | PASS | filter_chainrule.H:221-227; E4 dot test |
| B0-15 | Fixed-region gradients strictly zero | PASS | filter_chainrule.H:10-13,216-232; dDTDxh/dAlphaDxh zeroing; E6 |
| B0-16 | eta5 implicit derivative included | PASS | filter_chainrule.H:74-161 (dProjectionDeta + corrections); E5 recomputes eta5 |
| B0-17 | dfdx same sign convention as J | PASS | J=-Tmix/Tref; dfdx from adjoint with dJ/dT<0 source; F1 sign check |
| B0-18 | dgdx[1] same sign convention as gDP | PASS | pressureConstraintDerivative (+inlet/-outlet); F2/Gate6 FD |
| B0-19 | Frozen turbulence fields fixed under perturbation | PASS | updateFrozenTurbulenceFields.H guards; NS.H:106-109; no turbulence->correct() in F1/F2/F3 loops |
| B0-20 | Code ready for Stage B1 | PASS_WITH_WARNINGS | Chains consistent; see section 12 for the three B1 preconditions |

## 12. Findings

1. [WARNING - structural, not a defect] In `freezeColdFlowForValidation=true` mode (the Stage-A F1 configuration), `gsenshPressureDrop` is set to zero (`sensitivity.H:57-62`) and `Uc` is never solved (AdjNS_PD.H skipped, `MTO_HF.C:79-84`), so `dgdx[1] ≡ 0`. The F2 pressure-drop FD check (`validateFrozenGradient.H:902-924`) would then compare 0 against a non-zero FD of the actual pressure response and fail. The pressure chain is instead validated through the discrete flow adjoint dot-tests (`solveDiscreteFlowAdjoint.H:755-789`) and the objective-derivative check (`validateDiscreteObjectiveDerivatives.H`) in the non-frozen path. For Stage B1, run the pressure-gradient FD with `solveFlowAdjoints=true`, `freezeColdFlowForValidation=false`.
2. [WARNING] `solveDiscreteFlowAdjoint.H:1-8` aborts in parallel — the reduced discrete flow adjoint (Ub/Uc) is serial-only. Stage B formal FD must run serial.
3. [WARNING] The disk field `dfdx` (`sensitivity.H:200-218`) is the SCALED gradient, while `fsensMeanT`/`gsensVol`/`gsensPressureDrop` are raw. Report consumers must not mix the two.
4. [INFO] `thermalObjectiveDerivative` has two definitions: area-weighted at initialization (`createFrozenHotRegionFields.H:584-601`, used by `validateDiscreteObjectiveDerivatives.H` where phi may be zero) and phi-weighted at runtime (`computeObjective.H:26-56`, refreshed every outer iteration, `MTO_HF.C:54`). The runtime adjoint always uses the phi-weighted version. Not a conflict.
5. [INFO] `gsenshVol` multiplies `designMask` twice (definition + post-multiplication); harmless since mask ∈ {0,1}.
6. [INFO] `validateDiscreteObjectiveDerivatives.H` constructs its test direction from local cell indices (`sin(0.371*(celli+1))`) — the same serial-centric pattern as the removed Gate 6 local-index directions; harmless for a serial init-time check but not partition-invariant (serial-only usage documented).

## 13. Recommended minimal fixes (proposed, NOT applied)

1. In frozen-cold-flow mode, skip or mark N/A the F2/F3 pressure-drop FD comparison (or gate it on `solveFlowAdjoints && !freezeColdFlowForValidation`) to avoid a guaranteed false failure.
2. Write the B1 case configuration template documenting `objectiveGradientScale 1.0; pressureGradientScale 1.0;` as mandatory, and serial execution for the discrete flow adjoint.
3. Optionally rename the disk scalar `dfdx` to `dfdxScaled` (or add a `dfdxRaw` field) to disambiguate scaled vs raw output. Does not change any solve.
4. Document the dual definition of `thermalObjectiveDerivative` (init vs runtime) in the code header.

None of these change the objective, constraints, adjoint equations, gradient formulas, acceptance thresholds, or MMA. They were not applied in this audit.

## 14. Stage B1 readiness

**YES, with the following mandatory preconditions:**
- Run with `objectiveGradientScale = 1.0` and `pressureGradientScale = 1.0` (the unique safe configuration for amplitude FD).
- Run serial (discrete flow adjoint is serial-only).
- For the pressure-gradient FD, use the flow-coupled configuration (`freezeColdFlowForValidation=false`, `solveFlowAdjoints=true`); compare `dgdx[1]` against FD of the NORMALIZED `gDP = PressureDropPa/pressureDropMaxPa - 1`.
- For the thermal-only FD, the F1 gate already compares the scaled thermal gradient against FD of `-Tmix/Tref`.

`gradientValidated` remains `false`; no optimization was started.

## 15. Deliverables

- `/home/ys/TO-ANISOTROPIC/src/validation/stage_b/STAGE_B_B0_AUDIT_REPORT.md`
- `/home/ys/TO-ANISOTROPIC/src/validation/stage_b/stage_b_b0_audit.json`
- `/home/ys/TO-ANISOTROPIC/src/validation/stage_b/stage_b_b0_variable_trace.tsv`

## 16. Final status

**PASS (B0 + B0.1 + B0.2)** — the three function/gradient definitions are
mutually consistent, normalized exactly once, free of mesh.V() duplication or
omission, and free of hidden dynamic normalization. B0.1 added the objective
direct flux derivative `dJ/dphi` on the cold outlet; B0.2 proved the
internal-face thermal coupling term to be the exact discrete transpose for the
current temperature convection scheme and superseded the earlier
"inconsistent" conclusion (below). The only B1 preconditions are the three
documented ones: `objectiveGradientScale 1.0; pressureGradientScale 1.0;`,
serial execution (discrete flow adjoint is serial-only), and the flow-coupled
pressure FD configuration.

**Scheme-scope note:** every B0.2 transpose result applies to
`div(phi,T) bounded Gauss upwind` (the gate4/solver fvSchemes entry). If the
temperature convection scheme in `system/fvSchemes` is ever changed (e.g. to
plain `Gauss upwind`, `linearUpwind`, `limitedLinear`, or `vanLeer`), the
Stage B0.2 internal-face and boundary-face transpose verification MUST be
re-executed before any further gradient validation relies on it.


## Stage B0.1 — Mass-flow-weighted outlet objective flux derivative (2026-08-06)

See `STAGE_B_B0_1_OBJECTIVE_FLUX_AUDIT.md` for the full audit and the
`stage_b_b0_1_*` deliverables. Summary:

- J = -Tmix/Tref has the full direct derivative
  `dJ/dT_f = -phi_f/(M*Tref)` (was implemented) and
  `dJ/dphi_f = -(T_f - Tmix)/(M*Tref)` (was MISSING, now implemented in
  `computeObjective.H` + `AdjNS_HT.H`, declared in `createFrozenHotRegionFields.H`).
- Algebraic FD: dJ/dT max 9.5e-9; dJ/dphi median 1.65e-11 / max 2.88e-9;
  combined T-phi median 3.1e-11 (h<=1e-4) - all PASS.
- Thermal-residual flux derivative (fvMatrix::residual() semantics): outlet
  boundary term exactly zero (`dr == 0`), internal faces ~1.8e-14 (negligible),
  inlet ~-2e-14 (negligible).
- F1 regression after the fix: serial 0/9 failures (median 2.25e-4, max 2.49e-3),
  4-rank parallel 0/9 (median 9.95e-5, max 1.51e-3), dir1+2 serial/parallel max
  rel diff 3.65e-8 (<=1e-6). D_ADJ bit-identical to pre-fix (F1 fixes phi).
- B0-01, B0-17, B0-20 updated to PASS; readyForB1 = true; overallStatus =
  PASS (B0 + B0.1).
- ~~Pre-existing finding: the internal-face thermal-coupling term
  `-Tb_downwind*(T_nei-T_own)` does not equal the exact fvm residual
  derivative~~ **SUPERSEDED by Stage B0.2** (2026-08-06): that B0.1-era
  conclusion was an artifact of a cell-volume-weighted test; with the correct
  no-volume inner product the expression is the EXACT fvm transpose of
  `div(phi,T) bounded Gauss upwind` (median per-face FD 2.05e-10). It is NOT
  to be re-derived; no change is required.


## Stage B0.2 — Exact discrete thermal-to-flow coupling transpose (2026-08-06)

See `STAGE_B_B0_2_THERMAL_FLOW_TRANSPOSE_REPORT.md` and the
`stage_b_b0_2_*` deliverables. Summary:

- Actual discretization: `div(phi,T) bounded Gauss upwind`
  (= Gauss upwind fvmDiv - fvm::Sp(div(phi),T)).
- The B0.1 report's finding that the internal-face expression
  `-Tb_downwind*(T_nei-T_own)` is "inconsistent with the fvm residual
  derivative" is SUPERSEDED: the B0.1 test weighted the inner product by the
  cell volume although the fvMatrix residual is already volume-integrated.
  With the correct (no-volume) dot product the current expression is the
  EXACT fvm transpose (internal faces: median FD rel err 2.05e-10, max
  1.56e-9; sign 0).
- cold outlet: `dR_T/dphi_out == 0` exactly (upwind boundary internalCoeffs
  and the bounded `-V*div(phi)` correction cancel; measured `d(diagFull)/dh=0`,
  `lambda_fd ~ 7e-15`). faceAdjoint on the outlet = dJ/dphi only (B0.1 fix
  confirmed).
- cold inlet: `lambda_phi = (T_cell - T_bnd)*Tb_cell` (max rel 2.95e-8), but
  the inlet flux is fixed by fixedValue U and never enters the flow unknowns.
- Full-objective smoke (smoke_b02: solveFlowAdjoints=true, frozen cold flow
  false): F1 0/9 (median 2.25e-4), F2 max 6.7%, F3 objective 2.4% (same
  sign), discrete flow-adjoint dot-test 7.44e-14.
- F1 regression serial+parallel PASS (bit-identical D_ADJ to B0.1; dir1+2
  serial/parallel rel diff 3.65e-8 <= 1e-6).
- Source change: `src/AdjNS_HT.H` +60 lines of diagnostic output only (no
  physics change). B0.1 fix (computeObjective.H / createFrozenHotRegionFields.H
  / AdjNS_HT.H outlet term) remains correct.
- B0 checks B0-01/B0-17/B0-20 remain PASS; readyForB1 = true; overallStatus =
  PASS (B0 + B0.1 + B0.2).


## Stage B1 — Frozen-turbulence repeatability & noise floor (2026-08-07)

See `STAGE_B_B1_REPEATABILITY_NOISE_REPORT.md` and the `stage_b_b1_*`
deliverables. Status: **PASS** (16/16 criteria). Key results and finding:

- Two-phase baseline: full-SST converged in 246 iterations; frozen baseline
  saved with field statistics (all nBad=0).
- `solveFrozenPrimalConverged` added (replaces fixed-300-iteration practice):
  dJrel<=1e-9 && dgDPrel<=1e-8 for 5 consecutive iterations + a 20-iteration
  platform window; the 1e-10 window is never reached (J limit-cycle
  oscillation ~1e-7) and the run honestly reports the 1000-iteration platform
  value (Iter=-1), which is nonetheless bit-reproducible.
- Test A (5 identical starts) and Test B (3 perturbed starts): all fields,
  gradients and projections bit-identical (J rel range 0, gDP rel range 0,
  gV range 0, drift 0, recovery 3.4e-15).
- SNR probe h=1e-3: SNR_J = 8.3e25 (gate 100), SNR_DP = 1.3e32 (gate 20).
- FINDING: the stock GAMG pressure solver is non-deterministic across
  identical restarts (J platform noise ~1e-5, SNR_J=9.1); the B1 case uses
  deterministic PCG+DIC, after which everything is bit-identical. This
  solver choice is a B2/B3 precondition.
