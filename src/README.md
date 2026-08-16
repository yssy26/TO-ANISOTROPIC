# TO-HF OpenFOAM Solver

OpenFOAM 7 thermal-fluid topology optimization solver with primal and adjoint paths for:

- thermal objectives;
- power dissipation;
- pressure-drop constraints.

The `agent/frozen-hot-region` development line adds a single-mesh conjugate
heat-transfer configuration with:

- a precomputed, frozen upper hot-side velocity field;
- a fixed solid separator wall;
- a lower cold-side topology design region;
- one global temperature field;
- active cold-side primal and adjoint Navier-Stokes equations.

The branch `agent/stage-a-anisotropic-validation` adds a fixed solid
conductivity tensor, strict validation gates and the Stage A acceptance suite.

## Repository scope

This repository tracks solver source only. It is not intended to store large
meshes, processor outputs, runtime logs or production result directories.

## Main entry point

- `MTO_HF.C`

## Main modules

- `NS.H`: cold-side primal Navier-Stokes solve;
- `HeatTransfer.H`: global primal temperature solve;
- `AdjHeatTransfer.H`: thermal adjoint and tensor thermal sensitivity;
- `AdjNS_HT.H`: cold-flow adjoint for the thermal objective;
- `AdjNS_FF.H`: power-dissipation flow adjoint;
- `AdjNS_PD.H`: pressure-drop flow adjoint;
- `createFrozenHotRegionFields.H`: region masks, frozen velocity and startup checks;
- `createFrozenTurbulenceFields.H`: density-damped frozen-RANS transport fields;
- `updateFrozenTurbulenceFields.H`: frozen eddy-viscosity and thermal-diffusivity update;
- `updateMaterialProperties.H`: scalar/tensor material interpolation and derivatives;
- `validateFrozenHotCase.H`: boundary, flow-model and validation gates;
- `validateGate6SSTStrict.H`: strict full-SST directional-derivative gate;
- `updateThermalFlux.H`: current-iteration combined hot/cold thermal flux;
- `costfunction.H`: objective, constraints and forward-health diagnostics;
- `sensitivity.H`: active-design sensitivity assembly and MMA update;
- `enforceFixedDesignRegions.H`: fixed hot/wall/prescribed-zone enforcement.

## Build

Use an OpenFOAM Foundation v7 environment:

```bash
source /home/ys/OpenFOAM-7/etc/bashrc
wclean
wmake
```

## Anisotropic-conductivity safety gate

The production Stage A scope is an axis-aligned diagonal solid conductivity
tensor. Off-diagonal components are disabled by default.

Initial validation settings:

```text
useAnisotropicConductivity       true;
anisotropicConductivityValidated false;
allowOffDiagonalAnisotropy       false;
mmaUpdateEnabled                 false;
```

When anisotropic conductivity is enabled, MMA is blocked until
`anisotropicConductivityValidated` is explicitly set to `true`. That flag may
only be enabled after the Stage A automatic report is `PASS`.

Run the algebraic tests with:

```bash
python3 tests/test_anisotropic_material_model.py
```

Read and execute:

- `ANISOTROPIC_STAGE_A_VALIDATION.md`;
- `CODEX_STAGE_A_EXECUTION_PROMPT.md`;
- `tests/stage_a_acceptance.py`.

## Frozen-RANS safety and validation gates

Frozen-hot mode supports `laminar` and `incompressibleRANSFrozen`. The RANS
path uses `kOmegaSST` for the primal flow and freezes the density-damped eddy
viscosity for the temperature and flow-adjoint solves. The solver does not
allow ordinary MMA updates until the relevant validation gates permit them.

The strict Gate 6 implementation:

- removes duplicate pressure-gradient normalization;
- requires both sign and magnitude agreement with full-SST central differences;
- requires at least four valid volume-neutral directions;
- checks repeated-baseline objective and pressure-constraint noise;
- restores the complete design, thermal and turbulence state after validation.

These checks do not create a full SST adjoint. The adjoint still omits the
`k`, `omega` and turbulent-viscosity design derivatives.

The first runs should normally use:

```text
gradientValidated false;
ransDirectionValidated false;
mmaUpdateEnabled false;
```

The thermal and frozen-flow adjoints include matrix/operator transpose tests.
These algebraic tests do not replace full finite-difference checks in which the
relevant primal equations are reconverged after each design perturbation.

## Detailed documentation

See:

- `FROZEN_HOT_REGION_IMPLEMENTATION.md` for the physical model and known limitations;
- `ANISOTROPIC_CONDUCTIVITY_PLAN.md` for the tensor implementation plan;
- `ANISOTROPIC_STAGE_A_VALIDATION.md` for mandatory anisotropic validation.
