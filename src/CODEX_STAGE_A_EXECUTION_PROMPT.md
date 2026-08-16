# Codex Execution Prompt: Stage A Anisotropic-Conduction Validation

Use the following prompt without weakening the acceptance criteria.

---

You are working on the OpenFOAM 7 repository:

```text
GitHub: yssy26/TO-ANISOTROPIC
branch: agent/stage-a-anisotropic-validation
solver source: /home/ys/TO-ANISOTROPIC/src
working cases: /home/ys/Twostream
```

Your task is to execute and document Stage A anisotropic-conduction validation.
Read these files first:

```text
src/ANISOTROPIC_STAGE_A_VALIDATION.md
src/readThermalProperties.H
src/updateMaterialProperties.H
src/HeatTransfer.H
src/AdjHeatTransfer.H
src/validateGate6SSTStrict.H
src/tests/test_anisotropic_material_model.py
src/tests/stage_a_acceptance.py
src/tests/stage_a_results.example.json
```

Do not start topology optimization. Do not set
`anisotropicConductivityValidated true` unless the automatic acceptance script
reports PASS.

## Required work

1. Confirm the checked-out branch and record the exact commit SHA.
2. Build in the real OpenFOAM Foundation v7 environment:

```bash
cd /home/ys/TO-ANISOTROPIC/src
source /home/ys/OpenFOAM-7/etc/bashrc
wclean
wmake 2>&1 | tee validation_stage_a_build.log
```

3. Run:

```bash
python3 tests/test_anisotropic_material_model.py \
    2>&1 | tee validation_stage_a_unit.log
```

4. If compilation fails, make only the minimum OpenFOAM-7 compatibility fix.
Do not alter the mathematical model or relax validation thresholds. Record the
exact error and patch.

5. Create a self-contained directory:

```text
/home/ys/Twostream/validation/stage_a/
```

Reuse a known working validation case as the template. Do not modify the
production case in place.

6. Execute all checks A1-A7 exactly as specified in
`ANISOTROPIC_STAGE_A_VALIDATION.md`:

- scalar/tensor isotropic regression;
- x-, y- and z-direction one-dimensional conduction;
- three-layer thermal resistance;
- physical directionality and conductivity-ratio check;
- fixed-temperature-boundary sensitivity finite difference;
- full thermal-adjoint raw-design finite difference in at least four directions;
- serial/four-rank thermal-field and gradient consistency.

7. Every finite-difference comparison must perturb raw `x`, then recompute the
filter, projection and material interpolation. Do not perturb `xh` directly
when comparing with `dJ/dx`.

8. For every `+h` and `-h` pair, restore the same converged baseline before the
solve. Use at least:

```text
h = 1e-3, 3e-4, 1e-4
```

Report whether a finite-difference plateau exists. Do not report only the best
step size.

9. Use the assembled temperature-equation flux for heat-rate comparisons.
Do not substitute a separate post-processing gradient formula.

10. Fill a copy of:

```text
src/tests/stage_a_results.example.json
```

and save it as:

```text
/home/ys/Twostream/validation/stage_a/stage_a_results.json
```

Every `null` must be replaced by a measured finite number.

11. Run the mandatory automatic judgement:

```bash
cd /home/ys/TO-ANISOTROPIC/src
python3 tests/stage_a_acceptance.py \
    /home/ys/Twostream/validation/stage_a/stage_a_results.json \
    --report /home/ys/Twostream/validation/stage_a/STAGE_A_VALIDATION_REPORT.md \
    2>&1 | tee /home/ys/Twostream/validation/stage_a/stage_a_acceptance.log
```

Do not manually override a FAIL exit status.

## Required return package

Return all of the following:

1. exact tested commit SHA;
2. `git status --short`;
3. complete build result and any compile patch;
4. unit-test result;
5. one table covering A1-A7 with measured values and thresholds;
6. all finite-difference tables, including every direction and every `h`;
7. `stage_a_results.json`;
8. `STAGE_A_VALIDATION_REPORT.md`;
9. paths to all solver logs and comparison files;
10. a clear final status: `PASS`, `FAIL`, or `BLOCKED`.

Use `BLOCKED` only for an environmental problem that prevents execution, such
as a missing OpenFOAM installation or missing template case. A numerical or
validation failure is `FAIL`, not `BLOCKED`.

If any check fails, identify whether the likely cause is:

- tensor Laplacian discretization;
- face interpolation across conductivity jumps;
- physical-boundary sensitivity assembly;
- projection/filter chain;
- serial/processor-face sensitivity assembly;
- insufficient nonlinear convergence;
- or a different source supported by evidence.

Do not proceed to SST topology optimization or publication runs.

---
