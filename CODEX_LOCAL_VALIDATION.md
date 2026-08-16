# Codex Local Validation Runbook

This runbook validates the post-`ae3a252` frozen-RANS discrete-adjoint fixes on a local OpenFOAM-7 installation.  The validation order is intentionally strict: compile -> production adjoint smoke test -> Stage B4 regression -> Stage B2/B3 amplitude gate -> guarded MMA smoke test.

## Non-negotiable rules for Codex

Codex may fix compile/runtime defects introduced by the current branch, but it must **not**:

- change the objective or constraints;
- change `objectiveGradientScale` or `pressureGradientScale` away from `1.0`;
- weaken Stage B2/B3 error or plateau thresholds;
- set `frozenGradientValidated true` before a formal PASS;
- enable MMA during Stage B validation;
- use `explicitJT.mtx` alone as a full frozen-RANS `J^T` oracle;
- add a speculative `alpha -> phi` design derivative unless a same-basis finite-difference oracle demonstrates a non-negligible missing term;
- hide a non-converged adjoint by loosening the production tolerance beyond the stated validation limits.

Preserve the already-validated Stage-B4 `J/J^T` operator while diagnosing the production linear solver.

## 0. Synchronize and record provenance

```bash
cd /home/ys/TO-ANISOTROPIC
git fetch origin
git checkout agent/stage-a-anisotropic-validation
git pull --ff-only
git status --short
git log -8 --oneline
```

Record the exact HEAD SHA in every report.

## 1. Clean build

Locate the actual OpenFOAM-7 environment first.  Then perform a clean build:

```bash
source /opt/openfoam7/etc/bashrc   # replace only if the local OF7 path differs
cd /home/ys/TO-ANISOTROPIC/src
wclean
wmake
```

Acceptance:

```text
BUILD = PASS
```

If compilation fails, Codex should make the smallest source-level compatibility fix, rebuild from clean state, and document the exact compiler error and patch.  Do not alter equations or validation thresholds merely to compile.

## 2. Production adjoint smoke test

Use a copy of the frozen-validation case.  Before the full finite-difference sweep, verify that both production coupled adjoints converge with the clean solver path.

Required `constant/optProperties` settings:

```text
mmaUpdateEnabled              false;
stageBEnabled                 false;
stageB2Enabled                false;
stageB4JacobianProbe          false;
frozenGradientValidated       false;
adjointMode                   discrete;
flowModel                     incompressibleRANSFrozen;
frozenTurbulenceAdjoint       true;
solveFlowAdjoints             true;
objectiveGradientScale        1.0;
pressureGradientScale         1.0;
discreteFlowAdjointRestart    80;
discreteFlowAdjointMaxIter    4000;
discreteFlowAdjointTolerance  1e-9;
discreteUseExplicitSolution   false;
```

Run one solved frozen state and collect these lines:

```text
Production reduced discrete flow adjoint thermalCoupling
Production reduced discrete flow adjoint pressureDrop
```

Acceptance for both systems:

```text
true relative residual <= 1e-9
finite Ub/pb/Uc/pc
no FatalError
```

If FGMRES stalls, first try solver-only changes such as restart 120 or 160 and `discreteFlowAdjointMaxIter 8000`; keep the operator and validation tolerance unchanged.  Save the true residual history.  If simple Jacobi remains insufficient, implement a better right preconditioner while keeping `applyProdFlowJT` unchanged and prove the same converged solution by true residual.

Do **not** proceed to B2/B3 with a non-converged production adjoint.

## 3. Stage B4 operator regression (separate run)

B4 remains the diagnostic/oracle path and must be run separately from B2/B3.

Use:

```text
mmaUpdateEnabled             false;
stageBEnabled                false;
stageB2Enabled               false;
stageB4JacobianProbe         true;
frozenGradientValidated      false;
discreteUseExplicitSolution  false;
```

Check that the earlier operator closure has not regressed.  Collect at least:

- full `J/J^T` dot-test error;
- PU / PP / UU / UP block-dot errors;
- H1/H4 reduced-`H()` diagnostics;
- P-row/reduced-flux finite-difference diagnostics;
- explicit-core vs matrix-free oracle diagnostics.

Expected qualitative result: the previously closed dot/oracle checks remain at near-machine precision / their established tolerances.  B4's legacy experimental linear solver is diagnostic-only and is not used to approve gradients.

Important: `explicitJT.mtx` omits the frozen deviatoric-stress transpose that the diagnostic code adds matrix-free.  Do not solve that file alone and label the result a full frozen-RANS adjoint.

## 4. Stage B2/B3 formal frozen-gradient amplitude validation

Use a fresh copy of the baseline validation case and set:

```text
mmaUpdateEnabled              false;
stageBEnabled                 false;
stageB2Enabled                true;
stageB4JacobianProbe          false;
frozenGradientValidated       false;
adjointMode                   discrete;
flowModel                     incompressibleRANSFrozen;
frozenTurbulenceAdjoint       true;
solveFlowAdjoints             true;
objectiveGradientScale        1.0;
pressureGradientScale         1.0;
discreteFlowAdjointRestart    80;
discreteFlowAdjointMaxIter    4000;
discreteFlowAdjointTolerance  1e-9;
discreteUseExplicitSolution   false;
```

Run the solver in serial.  Inspect:

```text
stageB2/stage_b2_fd_scan.tsv
stageB2/stage_b2_summary.tsv
```

Formal PASS requires **every direction D1/D2/D3** to satisfy:

```text
all +/-h frozen primals converge
J:    sign PASS and bestRel <= 0.05
 gDP: sign PASS and bestRel <= 0.10
 gV:  sign PASS and bestRel <= 1e-6
J/gDP/gV: adjacent-step plateau PASS (spread <= 0.30)
```

The final summary must contain:

```text
status  PASS
```

Only then should the run create:

```text
stageB2/FROZEN_GRADIENT_UNLOCK.txt
```

If B2/B3 fails, Codex must decompose the failure instead of changing thresholds.  Report separately:

1. production adjoint true residual;
2. D1/D2/D3 FD plateau quality;
3. objective thermal/direct-flux contribution;
4. pressure-drop adjoint projection;
5. Brinkman design derivative;
6. filter/projection chain;
7. any same-basis `R_x` oracle evidence.

## 5. Unlock only after reviewed PASS

After the formal report is reviewed, copy the generated setting into the production case:

```text
frozenGradientValidated true;
```

Keep all Stage-B modes off during optimization.

## 6. Guarded MMA smoke test

First optimization run:

```text
mmaUpdateEnabled          true;
stageBEnabled             false;
stageB2Enabled            false;
stageB4JacobianProbe      false;
frozenGradientValidated   true;
enableSSTAcceptanceCheck  true;
```

Use a conservative move limit, initially about `0.01-0.02`, and only 5-10 design updates.

For each accepted/rejected candidate record:

- frozen predicted objective change;
- frozen realized objective change;
- predicted and realized normalized pressure-drop change;
- volume constraint;
- full-SST candidate objective and pressure drop;
- acceptance/rejection reason;
- rollback integrity if rejected.

Acceptance of the smoke phase is qualitative plus safety-based: gradients should consistently move the design in the intended direction, constraints remain controlled, full-SST acceptance does not systematically reject the frozen-model direction, and rollback is exact when rejection occurs.

## 7. Final report requested from Codex

Create `LOCAL_FROZEN_GRADIENT_VALIDATION_REPORT.md` containing:

1. repository HEAD and case provenance;
2. build result;
3. production thermal-coupling and pressure-drop adjoint convergence histories;
4. Stage B4 regression table;
5. Stage B2/B3 D1/D2/D3 table for every formal step;
6. explicit PASS/FAIL for J, gDP, gV and plateau per direction;
7. whether `FROZEN_GRADIENT_UNLOCK.txt` was generated;
8. any code changes Codex made locally, with diffs and reasons;
9. recommendation: `BLOCK_MMA`, `READY_FOR_MMA_SMOKE`, or `READY_FOR_GUARDED_OPTIMIZATION`.

Do not mark the code ready for MMA if the production adjoint or any formal frozen-gradient gate fails.
