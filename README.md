# TO-ANISOTROPIC

OpenFOAM-7 implementation for two-stream heat-exchanger topology optimization with axis-aligned anisotropic solid conduction, frozen-RANS (k-omega SST baseline), reduced discrete cold-flow adjoints, Helmholtz filtering / Heaviside projection, and guarded MMA updates.

The intended production route is deliberately conservative:

1. converge a full SST baseline;
2. freeze `k`, `omega`, `nut`, `nuEffFrozen` and turbulent thermal transport;
3. solve the frozen cold-side `U-p-T` problem;
4. compute discrete thermal and cold-flow adjoints;
5. verify objective / pressure-drop / volume gradients against multi-step central finite differences;
6. unlock MMA only after the frozen-model amplitude gate passes;
7. use full-SST re-convergence as a higher-level candidate acceptance / periodic correction layer.

## Source layout

```text
src/
├── MTO_HF.C
├── NS.H
├── HeatTransfer.H
├── AdjHeatTransfer.H
├── AdjNS_HT.H / AdjNS_PD.H / AdjNS_FF.H
├── solveDiscreteFlowAdjointProduction.H   # production J^T + restarted FGMRES
├── solveDiscreteFlowAdjoint.H             # Stage-B4 diagnostic/oracle solver
├── validateMmaUnlockGate.H                # hard optimization safety gate
├── validateStageBRepeatability.H
├── validateStageB2GradientAmplitude.H
├── validateGate6SSTStrict.H
├── sensitivity.H
├── filter_x.H / filter_chainrule.H
└── MMA/
```

## Discrete cold-flow adjoint

Two solver paths are intentionally separated.

### Production path

When `stageB4JacobianProbe false`, `AdjNS_HT.H` and `AdjNS_PD.H` use `solveDiscreteFlowAdjointProduction.H`.

The production solver:

- uses the same reduced `J^T` terms already closed by the Stage-B4 dot tests;
- includes the frozen explicit deviatoric-stress transpose;
- preserves the corrected `fvMatrix::H()` `/V` semantics and boundary diagonal;
- solves the scaled coupled system with restarted right-preconditioned FGMRES;
- uses a signed Jacobi preconditioner (no experimental ILU0);
- checks the true matrix-free residual;
- aborts if the adjoint does not meet `discreteFlowAdjointTolerance` instead of passing a partial solution to `sensitivity.H`.

Recommended initial controls:

```text
discreteFlowAdjointRestart    80;
discreteFlowAdjointMaxIter    4000;
discreteFlowAdjointTolerance  1e-9;
```

### Stage-B4 diagnostic path

When `stageB4JacobianProbe true`, the legacy large `solveDiscreteFlowAdjoint.H` is retained because it contains the H1-H4, P-row, block-dot, Rx and explicit-operator diagnostics.

Do **not** run Stage B2/B3 and B4 in the same execution. `validateMmaUnlockGate.H` enforces this separation.

Important: in frozen RANS the legacy `explicitJT.mtx` is not, by itself, the complete operator; the deviatoric-stress transpose is applied matrix-free by the diagnostic code. An external solution is therefore rejected unless the case explicitly declares that it was generated from the full operator.

## Frozen-gradient validation

`validateStageB2GradientAmplitude.H` performs D1/D2/D3 multi-step central finite differences.

Formal acceptance is per direction, not a global best-step metric:

```text
J:    sign correct + best relative error <= 5%
gDP:  sign correct + best relative error <= 10%
gV:   sign correct + best relative error <= 1e-6
FD:    at least one adjacent-step plateau <= 30% for each metric/direction
+/-h: every frozen primal must converge
```

A PASS writes:

```text
stageB2/FROZEN_GRADIENT_UNLOCK.txt
```

After reviewing the report, copy the generated line into `constant/optProperties`:

```text
frozenGradientValidated true;
```

Do not set this flag by hand before a PASS.

## MMA safety gate

Production MMA is blocked unless either:

```text
frozenGradientValidated true;
```

or the explicit direction-only feasibility safeguard is enabled together with full-SST candidate acceptance.

Validation/probe modes (`stageBEnabled`, `stageB2Enabled`, `stageB4JacobianProbe`) are mutually incompatible with active MMA updates.

Recommended validation configuration:

```text
mmaUpdateEnabled              false;
stageB2Enabled                true;
stageB4JacobianProbe          false;
frozenGradientValidated       false;
objectiveGradientScale        1.0;
pressureGradientScale         1.0;
freezeTurbulenceForValidation true;
solveFlowAdjoints             true;
```

Recommended first production smoke test after a PASS:

```text
mmaUpdateEnabled          true;
stageB2Enabled            false;
stageB4JacobianProbe      false;
frozenGradientValidated   true;
enableSSTAcceptanceCheck  true;
```

Use a small move limit for the first optimization run.

## Build

```bash
source /opt/openfoam7/etc/bashrc
cd src
wclean
wmake
```

## Validation status and scope

- Stage A: axis-aligned diagonal anisotropic conduction path validated; off-diagonal anisotropy remains a separate validation task.
- Stage B1: frozen-primal repeatability/noise framework implemented.
- Stage B2/B3: objective, normalized pressure-drop and volume gradient amplitude gate implemented with per-direction acceptance.
- Stage B4: reduced `J/J^T` diagnostic/oracle framework retained for operator debugging.
- Reduced discrete flow adjoint remains serial-only until processor-patch transpose terms are implemented.

Full-SST Gate 6 is a physical consistency/correction check; it is not a substitute for the frozen-model B2/B3 amplitude gate.
