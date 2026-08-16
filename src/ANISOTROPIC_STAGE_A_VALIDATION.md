# Stage A: Anisotropic-Conduction Validation

## 1. Purpose

Stage A validates the anisotropic molecular thermal-diffusivity path before it
is used in topology updates. It covers the material interpolation, tensor
Laplacian, physical boundary contributions, thermal adjoint sensitivity and
serial/parallel consistency.

The solver now blocks MMA when all of the following are true:

```text
useAnisotropicConductivity      true;
anisotropicConductivityValidated false;
mmaUpdateEnabled                true;
```

Set `anisotropicConductivityValidated true` only after the final acceptance
script reports `PASS`.

The production scope of Stage A is an axis-aligned diagonal tensor:

```text
solidK solidK [1 1 -3 -1 0 0 0]
    (Kxx 0 0 Kyy 0 Kzz);
```

Off-diagonal components remain disabled by default. They require a separate
rotated-tensor flux and adjoint validation.

## 2. Required environment

- OpenFOAM Foundation v7
- Repository branch: `agent/stage-a-anisotropic-validation`
- A known working TO-HF validation case containing all required fields and
  dictionaries
- Python 3.8 or newer

Recommended working directory:

```bash
/home/ys/Twostream/validation/stage_a
```

## 3. Acceptance sequence

Do not skip or reorder the checks.

### A0. Build and algebraic unit test

```bash
cd /home/ys/TO-ANISOTROPIC/src
source /home/ys/OpenFOAM-7/etc/bashrc
wclean
wmake
python3 tests/test_anisotropic_material_model.py
```

Required result:

```text
ANISOTROPIC MATERIAL MODEL UNIT TESTS: PASS
```

The test verifies:

- scalar/tensor isotropic regression;
- analytical `dD/dxh` against central finite differences;
- x/y/z directional conductivity;
- pure-solid and pure-fluid material endpoints.

### A1. Isotropic solver regression

Run the same fixed design twice:

1. scalar fallback:

```text
useAnisotropicConductivity false;
ks 19.1;
```

2. tensor path:

```text
useAnisotropicConductivity true;
solidK (19.1 0 0 19.1 0 19.1);
anisotropicConductivityValidated false;
mmaUpdateEnabled false;
```

Use identical initial fields, mesh, decomposition and solver settings. Compare
all common scalar, vector and tensor output fields after reconstruction.

Acceptance:

```text
maxRelativeFieldError <= 1e-12
filesCompared >= 1
```

A bit-for-bit result is preferred when the two cases use the same compiler and
execution mode.

### A2. One-dimensional x/y/z conduction

Create three conduction-only cases using the same uniform orthogonal hexahedral
mesh. In each case:

- set all velocities to zero;
- use `freezeColdFlowForValidation true`;
- use `solveFlowAdjoints false` for the forward-only check;
- set `mmaUpdateEnabled false`;
- remove volumetric heat generation;
- impose fixed temperatures on the two faces normal to the tested axis;
- set the remaining faces to `zeroGradient`;
- use a uniform pure-solid design;
- use `solidK (21.6 0 0 21.6 0 16.6)`.

For axis `i`, the analytical heat rate is

```text
Q_i = K_ii A DeltaT / L_i
```

Check the domain-integrated heat flux and the linearity of the centerline
temperature profile.

Acceptance for every axis:

```text
relativeFluxError <= 0.01
temperatureLinearityR2 >= 0.999
```

### A3. Three-layer thermal-resistance case

Use a one-dimensional stack with three conformal cell layers representing hot
fluid, separator solid and cold material. Set both velocities to zero.

The analytical heat rate is

```text
Q = A DeltaT /
    (L_hot/k_hot + L_wall/k_wall,normal + L_cold/k_cold,normal)
```

Use the conductivity component normal to the layers. Compare the numerical
heat rate with the analytical value and sample temperatures on both sides of
each internal interface.

Acceptance:

```text
relativeFluxError <= 0.02
maxInterfaceTemperatureJump <= 1e-6 K
```

If the flux error exceeds the limit, test a harmonic face-conductance scheme
before changing the tensor sensitivity implementation.

### A4. Physical directionality

For the diagonal tensor

```text
solidK (21.6 0 0 21.6 0 16.6)
```

run otherwise identical x- and z-direction conduction cases. Confirm:

```text
Q_x > Q_z
Q_x/Q_z approximately equals 21.6/16.6
```

Acceptance:

```text
orderingPass = true
relativeConductivityRatioError <= 0.05
```

This check must use the numerical equation flux, not a post-processed gradient
that bypasses the assembled tensor Laplacian.

### A5. Boundary sensitivity finite difference

Use a conduction-only case with designable cells adjacent to a fixed-temperature
boundary. This specifically exercises the physical-boundary contribution in
`AdjHeatTransfer.H`.

For at least three deterministic raw-design directions:

1. perturb raw `x`, not `xh`;
2. recompute filter, projection and material properties;
3. solve the complete temperature equation for `x+h d` and `x-h d`;
4. compare the central finite difference with the adjoint derivative;
5. use an epsilon ladder such as `1e-3`, `3e-4`, `1e-4`;
6. verify a stable finite-difference plateau.

Acceptance:

```text
directions >= 3
maximumRelativeError <= 0.05
```

### A6. Full thermal-adjoint finite difference

Freeze cold velocity and turbulent transport coefficients so that this test
isolates the anisotropic thermal operator and the complete raw-design chain.
Use at least four volume-neutral deterministic directions distributed through
the active design domain.

For each direction compare

```text
gFD = [J(x+h d)-J(x-h d)]/(2h)
gAD = sum_i (dJ/dx_i) d_i
```

Use the same epsilon ladder for all directions. The `+h` and `-h` cases must
start from the same converged baseline.

Acceptance:

```text
directions >= 4
maximumRelativeError <= 0.05
fdPlateauPass = true
```

### A7. Serial/parallel consistency

Run the same fixed anisotropic case in serial and with four MPI ranks. Reconstruct
parallel results and compare at least:

```text
T
DTMolecular
DTEffective
dDTDxh
fsensMeanT
objective
```

Acceptance:

```text
maxRelativeFieldError <= 1e-10
objectiveRelativeError <= 1e-10
gradientRelativeError <= 1e-8
```

The flow discrete adjoint remains serial-only. This Stage A check concerns the
thermal tensor path and thermal sensitivity assembly.

## 4. Result file and automatic judgement

Copy the schema:

```bash
cp src/tests/stage_a_results.example.json \
   /home/ys/Twostream/validation/stage_a/stage_a_results.json
```

Replace every `null` and placeholder with measured values. Then run:

```bash
cd /home/ys/TO-ANISOTROPIC/src
python3 tests/stage_a_acceptance.py \
    /home/ys/Twostream/validation/stage_a/stage_a_results.json \
    --report /home/ys/Twostream/validation/stage_a/STAGE_A_VALIDATION_REPORT.md
```

The script exits with status zero only when all mandatory checks pass.

## 5. Unlock rule

Only after the report contains

```text
Overall result: PASS
```

may the case dictionary be changed to:

```text
anisotropicConductivityValidated true;
```

Keep the following production restriction unless a separate general-tensor
validation is completed:

```text
allowOffDiagonalAnisotropy false;
```

## 6. Evidence to preserve

Keep the following under the validation directory:

- build log;
- unit-test log;
- each case dictionary and mesh summary;
- solver logs;
- analytical calculations;
- raw field-comparison tables;
- finite-difference tables for every direction and epsilon;
- `stage_a_results.json`;
- `STAGE_A_VALIDATION_REPORT.md`.

Do not mark Stage A complete from screenshots or qualitative temperature plots
alone.
