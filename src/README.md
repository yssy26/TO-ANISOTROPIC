# TO-HF OpenFOAM Solver

OpenFOAM 7 thermal-fluid topology optimization solver with primal and adjoint paths for:

- thermal objectives;
- power dissipation;
- pressure-drop constraints.

The `agent/frozen-hot-region` branch adds a single-mesh conjugate heat-transfer configuration with:

- a precomputed, frozen upper hot-side velocity field;
- a fixed solid separator wall;
- a lower cold-side topology design region;
- one global temperature field;
- active cold-side primal and adjoint Navier–Stokes equations.

## Repository scope

This repository tracks solver source only. It is not intended to store large meshes, processor outputs, runtime logs or production result directories.

## Main entry point

- `MTO_HF.C`

## Main modules

- `NS.H`: cold-side primal Navier–Stokes solve;
- `HeatTransfer.H`: global primal temperature solve;
- `AdjHeatTransfer.H`: thermal adjoint;
- `AdjNS_HT.H`: cold-flow adjoint for the thermal objective;
- `AdjNS_FF.H`: power-dissipation flow adjoint;
- `AdjNS_PD.H`: pressure-drop flow adjoint;
- `createFrozenHotRegionFields.H`: region masks, frozen velocity and startup checks;
- `createFrozenTurbulenceFields.H`: density-damped frozen-RANS transport fields;
- `updateFrozenTurbulenceFields.H`: frozen eddy-viscosity and thermal-diffusivity update;
- `validateFrozenHotCase.H`: boundary, flow-model and validation gates;
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

## Safety and validation gates

Frozen-hot mode supports `laminar` and `incompressibleRANSFrozen`. The RANS path
uses `kOmegaSST` for the primal flow and freezes the density-damped eddy viscosity
for the temperature and discrete adjoint solves. The solver does not allow MMA
updates until the case-level direction gates explicitly permit them:

```text
gradientValidated true;
ransDirectionValidated true;
mmaUpdateEnabled true;
```

The first runs should use both values as `false`, followed by mask, energy-balance, leakage, mesh-convergence and directional-derivative checks.
With MMA disabled, the solver recomputes `xp/xh` from raw `x` once at startup, keeps material-continuation parameters fixed, and uses `validationProjectionBeta` as the fixed projection strength. This makes perturbed raw-design cases evaluate the full design mapping without advancing MMA.

The thermal and frozen-flow adjoints include matrix/operator transpose tests.
These algebraic tests do not replace full finite-difference checks in which the
RANS equations are reconverged after each design perturbation.

## Detailed documentation

See [`FROZEN_HOT_REGION_IMPLEMENTATION.md`](FROZEN_HOT_REGION_IMPLEMENTATION.md) for:

- physical equations and assumptions;
- required cell zones and boundary conditions;
- dictionary entries;
- design-domain enforcement;
- objective and sensitivity definitions;
- validation sequence;
- known limitations.
