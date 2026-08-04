# Frozen Hot-Side / Cold-Side Topology Optimization

## 1. Scope

This branch implements the corrected **scheme 2** for the `TO-HF` OpenFOAM 7 solver:

- the upper hot-side channel is fixed and is not designed;
- the hot-side velocity field is precomputed and frozen;
- the separator wall is a fixed solid region;
- only the lower cold-side region enters the MMA design vector;
- the existing cold-side primal and adjoint Navier–Stokes systems remain active;
- one conformal mesh and one global temperature field represent the hot fluid, separator wall and cold-side topology domain.

This is not a VOF or moving-interface model and it is not a `chtMultiRegionFoam` multi-mesh solver.

The existing `TO-HF` convention is retained:

```text
x = 1  -> fluid
x = 0  -> solid
```

## 2. Model structure

The mesh is divided by cell zones:

```text
hotFluidZone        fixed hot fluid
separatorWallZone   fixed solid separator
remaining cold cells
    cold design cells
    optional fixed-fluid cells
    optional fixed-solid cells
```

The following fields are used:

```text
U, p                 solved cold-side flow
UHotFrozen           prescribed hot-side velocity
T                    global temperature
Tb                   global thermal adjoint
Ua, pa                power-dissipation flow adjoint
Ub, pb                thermal-objective flow adjoint
Uc, pc                pressure-drop flow adjoint
x, xp, xh            raw, filtered and projected design variables
```

No hot-side pressure or velocity adjoint is introduced because `UHotFrozen` is prescribed and independent of the design.

## 3. Forward equations

The cold-side flow continues to use the existing Brinkman formulation. In `hotFluidZone` and `separatorWallZone`, the cold flow is blocked using `coldBlockAlpha`.

The normalized thermal equation is

```text
div(phiThermal,T) - laplacian(DT,T) = 0
```

in frozen-hot mode, with

```text
phiThermal
    = coldFaceMask*phi
    + (rhoCpHot/rhoCpCold)*hotFaceMask*phiHotFrozen
```

and

```text
DT = k/rhoCpCold
```

The fixed-region values are

```text
hot fluid:       DT = hotK/rhoCpCold
separator wall:  DT = wallK/rhoCpCold
cold topology:   existing interpolated DT(xh)
```

`phiThermal` is assembled **after** the current cold-side `NS.H` solve. This prevents the temperature equation from using the previous optimization iteration's `phi`.

The cell-to-face masks are binary. Internal faces require both owner and neighbour cells to belong to the same flowing region. Coupled/processor faces also use the neighbour-side mask. Therefore the convective term cannot cross the separator wall.

## 4. Thermal objective

Frozen-hot mode uses

```text
thermalObjectiveType maximizeColdOutletTemperature;
```

The thermal objective is

```text
J_T = -T_c,out,area / thermalReferenceTemperature
```

where `T_c,out,area` is the area-averaged cold-outlet temperature. The complete objective remains

```text
J = (1-w1)*J_T + w1*PowerDiss
```

with the existing fixed `PowerDiss0` normalization.

The area-averaged outlet-temperature objective was chosen instead of an enthalpy-flow objective because it does not add a separate explicit outlet-velocity derivative. The cold velocity still affects the objective through the temperature equation and the existing thermal flow-adjoint coupling.

Legacy mode retains the original mean-temperature objective:

```text
thermalObjectiveType legacyMeanTemperature;
```

The two objective types are intentionally mode-locked to avoid silently running the wrong optimization problem.

## 5. Design-domain enforcement

The hot region, separator wall and existing prescribed solid/fluid zones are not merely assigned zero sensitivities. They are removed from the MMA design vector.

`activeDesignCells` contains only cells for which

```text
designMask > 0.5
```

MMA is constructed with this active list only. Fixed cells are therefore unable to move due to MMA regularization terms or asymptote updates.

Prescribed values are also restored:

1. before material interpolation;
2. before the density PDE filter;
3. after the density PDE filter;
4. after Heaviside projection;
5. before writing solved states.

The Heaviside volume-preserving threshold is calculated only over active design cells.

## 6. Sensitivity treatment

The branch removes the previous per-iteration maximum-gradient normalization. Objective and constraint values now use fixed reference quantities, and the gradients retain the matching physical scaling.
Sensitivity densities are multiplied by the actual cell volume before entering MMA; no `maxCellVolume` normalization is applied.

Because the volume-preserving projection threshold depends on all active filtered variables, the chain rule includes the implicit `d eta/d xp` term. The bisection tolerance is controlled by `projectionEtaTolerance`.

The interpolation derivatives are evaluated as

```text
dAlphaDxh = d alpha/d xh
dDTDxh    = d DT/d xh
```

with

- zero derivative in fixed regions;
- zero `dAlphaDxh` where `alpha` is clipped by `alphamin`;
- design masking before and after the adjoint PDE filter.

The pressure-drop adjoint boundary conditions already include `rhoFluid/pressureDropMaxPa`; this factor is not applied a second time in `sensitivity.H`.

The old one-cell finite-difference block was removed because it perturbed `xh` while comparing against a gradient with respect to raw `x`, altered the primal state without fully restoring it, and was not parallel safe.

## 7. Discrete adjoint status

The thermal adjoint transposes the assembled primal temperature matrix. The
outlet-temperature objective derivative is distributed to the owner cells of
the cold outlet faces. The frozen-flow adjoint uses the same `nuEffFrozen`
Laplacian as the primal equation and explicitly transposes the deviatoric-stress
operator in reverse order.

Matrix/operator point-product tests are run for both paths. Full design
direction tests are still required because the frozen-RANS approximation omits
the derivatives of `k`, `omega`, and `nutFrozen`. Until those tests pass:

```text
gradientValidated false;
ransDirectionValidated false;
mmaUpdateEnabled false;
```

must be used initially.

The solver terminates if MMA is enabled while the required validation gate is
false. These flags may only be changed after case-level directional-derivative
tests pass.

## 8. Required dictionaries

### 8.1 `constant/thermalProperties`

```text
frozenHotRegion
{
    enabled                     true;
    strictValidation            true;

    hotZone                     hotFluidZone;
    wallZone                    separatorWallZone;
    hotVelocityField            UHotFrozen;

    hotInletPatch               hotInlet;
    hotOutletPatch              hotOutlet;
    coldInletPatch              inlet;
    coldOutletPatch             outlet;

    hotRhoCp
        hotRhoCp [1 -1 -2 -1 0 0 0] 1.20e3;

    hotK
        hotK [1 1 -3 -1 0 0 0] 4.0e-2;

    wallK
        wallK [1 1 -3 -1 0 0 0] 15.0;

    coldBlockAlpha
        coldBlockAlpha [0 0 -1 0 0 0 0] 1e10;

    hotMassBalanceTolerance      1e-6;
    hotDivergenceTolerance       1e-8;
}
```

The dimensions must match the existing `rhoc`, `kf`, `ks` and `alpha` definitions.

### 8.2 `constant/optProperties`

For initial forward and adjoint validation:

```text
flowModel                         incompressibleRANSFrozen;
frozenTurbulenceAdjoint           true;
turbulentPrandtl                  0.85;
turbulenceDampingExponent         3.0;
turbulenceMinCorrectors           40;
turbulenceMaxCorrectors           200;
turbulenceConvergenceTolerance    1e-6;

thermalObjectiveType             maximizeColdOutletTemperature;
thermalReferenceTemperature      300;

mmaUpdateEnabled                  false;
gradientValidated                 false;
ransDirectionValidated            false;
recomputeDesignMappingWhenMmaDisabled true;
validationProjectionBeta          8.0;
projectionEtaTolerance            1e-10;

thermalBalanceTolerance           1e-2;
coldLeakageTolerance              1e-8;
coldHotPatchVelocityTolerance     1e-12;
```

Only after a directional-derivative test passes:

```text
gradientValidated                 true;
ransDirectionValidated            true;
mmaUpdateEnabled                  true;
```

### 8.3 `constant/turbulenceProperties`

For `flowModel incompressibleRANSFrozen`:

```text
simulationType RAS;

RAS
{
    RASModel    kOmegaSST;
    turbulence  on;
    printCoeffs on;
}
```

The primal solver converges `U`, `p`, `k`, `omega`, and the effective frozen
eddy-viscosity field. The adjoint treats the converged turbulent transport
coefficients as constants; it does not include `dnut/dx`. Use
`simulationType laminar` with `flowModel laminar` to reproduce the laminar path.
LES and non-frozen RANS adjoints are not supported.

## 9. Boundary conditions

### Temperature `T`

```text
hotInlet       fixedValue
coldInlet      fixedValue
hotOutlet      zeroGradient
coldOutlet     zeroGradient
external walls zeroGradient unless another physical condition is intended
```

The hot-fluid/wall and wall/cold interfaces must be conformal internal faces, not external patches.

### Cold velocity `U`

`U` represents only the cold stream. The hot inlet and hot outlet patches must use fixed zero cold velocity. The solver checks these values at startup.

### Hot frozen velocity `UHotFrozen`

Provide `UHotFrozen` in the starting time directory. It should contain the precomputed hot velocity in `hotFluidZone`; internal values outside that zone are reset to zero.

The precomputed field must be conservative. The solver checks:

- hot inlet/outlet volume-flow balance;
- maximum `div(phiHotFrozen)`;
- values outside the hot region.

## 10. Mesh checks

At startup, the solver verifies:

- hot and wall zones exist and are non-empty;
- hot and wall zones do not overlap;
- no hot-fluid cell directly neighbours a cold-fluid cell;
- hot/cold inlet and outlet patch owner cells belong to the correct regions;
- coupled/processor face masks agree with neighbour-side region masks;
- serial pressure reference cells lie in the cold region.

The separator must form a continuous cell layer between the hot and cold flowing regions.

## 11. Runtime health gates

Before an MMA update, the solver checks:

```text
|Qhot - Qcold|/max(|Qhot|,|Qcold|) <= thermalBalanceTolerance
max |Ucold| in hot/wall             <= coldLeakageTolerance
```

`Qhot` and `Qcold` are reported in watts by multiplying the normalized flux integrals by `rhoc`.

Non-finite objective, constraint or gradient values terminate the run.

## 12. Output and restart behaviour

Fields are written after the primal and adjoint equations are solved for the current design and before MMA creates the next design. This keeps `x/xp/xh`, `U`, `p`, `T` and the adjoints synchronized in each output directory.

The final solved iteration is written even when it does not coincide with `writeInterval`.

The following diagnostics are explicitly written:

```text
designMask
hotRegionMask
wallRegionMask
coldRegionMask
hotFaceMask
coldFaceMask
UHotFrozen
phiHotFrozen
phiThermal
dAlphaDxh
dDTDxh
filtered and unfiltered sensitivity fields
```

`x`, `xp` and `xh` use `READ_IF_PRESENT`, allowing a warm restart from a saved design and field state. The MMA asymptote history and continuation counters are not checkpointed, so this is not an exact optimizer-state restart.

## 13. Validation sequence

Use the following sequence without skipping stages.

### Stage A: compile

```bash
source /opt/openfoam7/etc/bashrc
wclean
wmake
```

### Stage B: legacy regression

Disable `frozenHotRegion` and confirm that the legacy single-stream case still runs with `thermalObjectiveType legacyMeanTemperature`.

### Stage C: masks and frozen velocity

Run with

```text
mmaUpdateEnabled false;
gradientValidated false;
```

Inspect all cell and face masks, `UHotFrozen`, `phiHotFrozen`, and the startup conservation report.

### Stage D: three-layer conduction

Set both velocities to zero and compare the numerical heat flux with

```text
q = DeltaT/(Lh/kh + Lw/kw + Lc/kc)
```

This test is mandatory because large conductivity jumps are sensitive to the face interpolation selected by the laplacian scheme.

### Stage E: fixed-topology conjugate heat transfer

Verify:

- hot enthalpy decrease and cold enthalpy increase agree;
- cold velocity is negligible in hot/wall regions;
- temperature and heat flux are continuous through the separator;
- results converge with mesh refinement.

### Stage F: directional derivative

For a direction `delta x` containing active design cells only, compare

```text
gFD = [J(x + eps*delta x) - J(x - eps*delta x)]/(2*eps)
```

with

```text
gAdj = sum_i dJdx_i*delta x_i
```

Use several `eps` values and several directions. Set `recomputeDesignMappingWhenMmaDisabled true` and use the same fixed `validationProjectionBeta` in the baseline and both perturbed cases. The solver recomputes `xp/xh` from raw `x` at startup while keeping `alphaMax` and `qu` fixed. The entire forward problem must be converged for both perturbed designs. Do not perturb `xh` directly when comparing against `dJ/dx`.

A practical initial acceptance target is a relative error below 5%, followed by tighter checking as the discrete adjoint is refined.

### Stage G: optimization

Only after all previous stages pass, set

```text
gradientValidated true;
ransDirectionValidated true;
mmaUpdateEnabled true;
```

## 14. Known limitations

The following are intentionally not hidden by the implementation:

1. The thermal and frozen-flow transpose tests pass algebraically, but the RANS pressure-drop direction fails when SST is reconverged after design perturbations. MMA remains gated.
2. The conductivity face interpolation is still controlled by the case laplacian scheme. A three-layer analytic test is required; a dedicated harmonic face-conductance implementation may be needed for very large conductivity ratios.
3. The cold velocity is blocked from the hot/wall cells with a finite Brinkman coefficient rather than a separate cold-flow mesh. Leakage is monitored and gated.
4. Laminar and incompressible frozen-`kOmegaSST` modes are supported. The RANS adjoint omits turbulence-state and turbulence-design derivatives.
5. The M1/M2 outlet-temperature comparison is not mesh independent; current RANS thermal results are for algorithm development only.
6. Contact resistance, thermal radiation, compressibility, buoyancy and temperature-dependent hot-side momentum are not represented.
7. Warm restart does not restore MMA asymptote history or continuation state exactly.

## 15. Modified modules

```text
MTO_HF.C
createFields.H
readTransportProperties.H
readThermalProperties.H
createFrozenHotRegionFields.H
createFrozenTurbulenceFields.H
validateFrozenHotCase.H
enforceFixedDesignRegions.H
applyFrozenHotRegionProperties.H
updateFrozenTurbulenceFields.H
updateThermalFlux.H
update.H
HeatTransfer.H
AdjHeatTransfer.H
solveDiscreteThermalAdjoint.H
solveDiscreteFlowAdjoint.H
AdjNS_HT.H
costfunction.H
sensitivity.H
filter_chainrule.H
filter_x.H
diff.c
writeOptimizationState.H
opt_initialization.H
```

## 16. Build status

The GitHub connector environment does not contain OpenFOAM 7 or `wmake`. The branch has undergone a static source review against OpenFOAM 7 interfaces, but a real OpenFOAM compilation and case-level validation remain mandatory before merging or using the solver for publication results.
