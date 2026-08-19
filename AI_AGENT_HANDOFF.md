# TO-ANISOTROPIC — AI Agent Project Handoff

> **Purpose:** This is the authoritative onboarding document for a new AI coding agent taking over development of this repository.
>
> **Repository:** `yssy26/TO-ANISOTROPIC`  
> **Active development branch:** `agent/dsH-stage-b-validation`  
> **Source snapshot used for this handoff:** `4c499038b9dffe2a90e0e92dbdb734ab36f82242` (`4c49903`, "Add production pressure-GAMG preconditioner")  
> **Handoff date:** 2026-08-19  
> **Current macro-stage:** **B-final — frozen-turbulence gradient / production-adjoint closure**  
> **MMA status:** **LOCKED** (`mmaUpdateEnabled=false`, `frozenGradientValidated=false`)

> **Status note (2026-08-19, post-handoff):** **BFINAL-010 implemented and run.** The pressure-GAMG preconditioner matrix is now hard-gated to the scaled J^T P–P block by a four-group equivalence decomposition (`PRODPRECGAMGCHECK`: interior/boundary/effective/offDiag relL2 ≈ 1e-16 on both labels — the BFINAL-009 `diagRelL2=0.041` stop is retired as a bare-vs-effective-diagonal comparison artifact), plus a `laplacianSchemes` coefficient-interpolation guard (`PRODPRECGAMGSCHEME`; the interpolation word comes from the laplacianSchemes entry ITstream, not interpolationSchemes). Production-path result: **thermalCoupling converges for the first time** (FGMRES 751 iters, trueRelRes 9.49e-10 ≤ 1e-9; diagonal baseline stalled at 5.08e-5/4000 iters); **pressureDrop NaNs inside its first GAMG inner solve** (iter=2 exit, MMA safety gate aborts, MTO_RC=134) — this is the BFINAL-011 entry problem. See `evidence/agent-group/BFINAL-010/FINAL_REPORT.md`.

> **Status note (2026-08-19 late, BFINAL-011 cycle-1 PASS):** pressureDrop NaN root-caused and fixed. **The trigger was the BFINAL-010 gate's own non-const `lduMatrix::lower()` call, which materialises the symmetric preconditioner matrix into asymmetric storage** (PCG then rejected; GAMG switched to its asymmetric branch with scaleCorrection=false, diverging on the pressureDrop rhs). Fix: the gate now reads through a `const lduMatrix&` alias (values bit-identical; regression-tested). With storage kept symmetric, **pure default GAMG converges BOTH adjoint labels in the production iterative path — the first such run ever**: thermalCoupling 750 iters 9.4957e-10, pressureDrop 902 iters 9.9282e-10 (tol 1e-9), GRADPROXY plateaus +2217.480 / -205194.141, four-group gate values unchanged at 1e-16, full program MTO_RC=0 (~120 s). Diagnostic switches added (per-label adjoint isolation; `discreteProdPressurePrecSolver` GAMG|PCG; GAMG coefficient overrides; inner-psi non-finite fast-fail). Next: BFINAL-012 (frozen-gradient FD amplitude gate). See `evidence/agent-group/BFINAL-011/cycle-1/FINAL_REPORT.md`.

---

## 0. Read this first

This project has accumulated many diagnostic files and historical Stage-B reports. Do **not** assume that an older report is current merely because it is detailed.

When documents disagree, use this precedence:

1. **current branch source code**;
2. this handoff document;
3. newest `BFINAL-00N` final report/evidence corresponding to the source being discussed;
4. older `src/validation/stage_b/STAGE_B_*` reports only as historical context.

The current development philosophy is simple:

> **Gradient validation is a gate, not the research objective.**
>
> The goal is not to build a universal exact RANS/SST adjoint framework. The goal is to obtain a numerically stable frozen-turbulence topology-optimization method whose gradient direction and useful amplitude are verified by finite differences and which can reliably evolve a physically meaningful heat-exchanger topology.

A new agent should therefore **preserve validated mathematics, avoid widening scope, close B-final, then move to MMA as soon as the final gates pass.**

---

# 1. What this project is

`MTO_HF` / `TO-ANISOTROPIC` is an **OpenFOAM Foundation v7** solver and optimization framework for a **two-stream turbulent conjugate-heat-transfer topology-optimization problem**.

The physical structure contains:

- a **hot stream / hot-side region**;
- a separating / conducting solid structure;
- a **cold stream** whose design region can become fluid or solid through density-based topology optimization;
- anisotropic solid heat conduction with a fixed material principal direction in the current production scope.

The main engineering goal is to improve heat-exchange performance — preferably total heat transfer `Q`, or alternatively the cold-side mixed outlet temperature — while respecting:

- a **cold-side pressure-drop constraint**;
- a **solid-volume-fraction constraint**;
- manufacturable / numerically regularized topology through filtering and projection.

The final method is intended for a scientific study in which the optimized structure will be checked by full SST primal simulations and can later be compared with experiment. The research contribution is **not** "we implemented a general exact SST discrete adjoint".

---

# 2. Scope lock: what the final algorithm MUST and MUST NOT do

## 2.1 Required final physics / optimization scope

### Hot side

The intended production algorithm uses:

- frozen hot-side velocity `U_h`;
- frozen hot-side pressure `p_h`;
- frozen hot-side turbulence state (`k_h`, `omega_h`, `nu_t,h` and derived frozen transport quantities);
- **live hot-side temperature**, re-solved as the topology changes.

This is appropriate for the intended forced-convection regime in which the topology change is on the cold side and temperature feedback does not materially alter the hot-side flow field.

Relevant implementation documentation/files include:

- `src/FROZEN_HOT_REGION_IMPLEMENTATION.md`
- `src/createFrozenHotRegionFields.H`
- `src/applyFrozenHotRegionProperties.H`
- `src/validateFrozenHotCase.H`

### Cold side

At each outer optimization state, the cold side should use a **full SST-RANS primal state**. During sensitivity evaluation, turbulence is frozen:

- solve/update cold-side `U`, `p`, `k`, `omega`, `nu_t` with SST-RANS;
- freeze `k`, `omega`, turbulent viscosity and derived turbulent transport coefficients for the sensitivity stage;
- re-solve the frozen-turbulence `U/p/T` state as required by the validated optimization path;
- solve thermal and flow adjoints without introducing turbulence-variable adjoints.

Relevant files include:

- `src/NS.H`
- `src/createFrozenTurbulenceFields.H`
- `src/updateFrozenTurbulenceFields.H`
- `src/validateGate6SSTStrict.H`
- `src/validateSSTDirection.H`

### Material / topology model

The raw design variable `x` (legacy convention in this code: approximately `0 = solid`, `1 = fluid`) goes through:

```text
raw design x
    -> Helmholtz/density filtering
    -> Heaviside projection
    -> projected material variable
    -> Brinkman resistance alpha(xh)
    -> thermal conductivity / diffusivity tensor K(xh)
```

The design therefore controls both:

1. flow blockage / fluid-solid topology via **Brinkman penalization**;
2. solid/fluid thermal transport via the material interpolation, including the validated **anisotropic conductivity** path.

Relevant files:

- `src/filter_x.H`
- `src/filter_chainrule.H`
- `src/updateMaterialProperties.H`
- `src/enforceFixedDesignRegions.H`
- `src/harmonicSymmTensor.H/.C`
- `src/ANISOTROPIC_CONDUCTIVITY_PLAN.md`
- `src/ANISOTROPIC_STAGE_A_VALIDATION.md`

### Objective and constraints

The production optimization problem should remain conceptually simple:

```text
maximize heat-exchange performance (prefer Q; mixed cold-outlet temperature is an alternative)
subject to:
    cold-side pressure drop <= prescribed limit
    solid volume fraction <= prescribed limit
    0 <= design variable <= 1
```

The exact sign conventions and normalization used by the code must be read from the current `costfunction.H`, `computeObjective.H`, `sensitivity.H`, and case dictionaries rather than guessed from this conceptual statement.

### Gradient path

The intended final gradient chain is:

```text
thermal objective / constraint sources
    -> discrete thermal adjoint

pressure-drop source
    -> reduced discrete frozen-flow adjoint

solid volume
    -> analytic derivative

all derivatives w.r.t. projected design
    -> Heaviside chain rule
    -> filter transpose / chain rule
    -> raw design gradient
    -> validation gates
    -> guarded MMA update
```

Main files:

- `src/AdjHeatTransfer.H`
- `src/solveDiscreteFlowAdjointProduction.H`
- `src/sensitivity.H`
- `src/filter_chainrule.H`
- `src/validateFrozenGradient.H`
- `src/validateMmaUnlockGate.H`
- `src/MMA/`

## 2.2 Explicit non-goals

Do **not** expand this project into any of the following unless the user explicitly changes the scope:

- no full `k-omega SST` turbulence-variable adjoint;
- no `lambda_k` / `lambda_omega` production system;
- no derivatives of turbulence state with respect to design as a required final feature;
- no universal OpenFOAM adjoint framework;
- no support requirement for arbitrary turbulence models;
- no support requirement for arbitrary boundary-condition families;
- no general SIMPLE/PISO/PIMPLE adjoint abstraction;
- no machine-precision closure of every possible OpenFOAM residual as the research objective;
- no new physical model before the current frozen-turbulence path is closed and produces an optimizing topology.

Local operator-equivalence tests may legitimately use `1e-8` to `1e-12` tolerances, but this must not be confused with the final scientific requirement on a turbulent topology gradient. For the latter, stable sign, direction, amplitude and finite-difference plateau are what matter.

---

# 3. Intended production optimization loop

The end-state algorithm should look like this.

## Step 1 — Baseline / periodic full SST primal

Obtain a physically converged SST-RANS state. Hot-side flow is initialized/frozen according to the project model. Cold-side SST state is refreshed at the beginning of an optimization state and periodically as needed.

## Step 2 — Filter and project the design

```text
x_n -> filtered x -> projected xh
```

Apply fixed-design masks/regions before evaluating material properties.

## Step 3 — Update material properties

Build:

- Brinkman flow resistance `alpha(xh)`;
- molecular / effective thermal transport;
- anisotropic conductivity/diffusivity tensor in the validated fixed-direction form.

## Step 4 — Frozen-turbulence primal solve

Freeze turbulence quantities and solve the state used by the derivative definition. The finite-difference validation cases must perturb the **raw design** and re-run the same filter/projection/material/primal chain.

## Step 5 — Evaluate objective / constraints

At minimum track:

- heat-exchange objective;
- cold-side pressure drop;
- solid volume fraction.

## Step 6 — Solve adjoints

- thermal adjoint for the heat-transfer objective/coupling;
- production reduced discrete frozen-flow adjoint for flow-dependent gradients, especially pressure drop;
- analytic solid-volume derivative.

## Step 7 — Assemble raw-design gradients

Apply all direct terms plus state-adjoint terms, then propagate through projection and filter.

## Step 8 — Validation gate

Before MMA is allowed, require repeatability, noise-floor and finite-difference amplitude evidence on the current source version.

## Step 9 — MMA

Only after the gate is explicitly unlocked:

```text
x_n -> MMA -> candidate x_(n+1)
```

Use the existing accepted-state / candidate-evaluation / rollback safeguards rather than introducing an unguarded design update.

## Step 10 — SST refresh and final validation

Periodically refresh with a full SST primal solution. The final topology must be validated with full SST-RANS, not only the frozen sensitivity state.

---

# 4. Code architecture: where to look

The solver uses the OpenFOAM style of including many `.H` implementation fragments into `MTO_HF.C`. A new agent should understand the include/control flow before editing isolated files.

| File / area | Role |
|---|---|
| `src/MTO_HF.C` | Main program and high-level orchestration |
| `src/createFields.H` | Main fields and field initialization |
| `src/opt_initialization.H` | Optimization state and controls |
| `src/readTransportProperties.H` | Flow/material transport parameters |
| `src/readThermalProperties.H` | Thermal/material parameters |
| `src/NS.H` | Primal SIMPLE/RANS flow path and frozen-turbulence logic |
| `src/HeatTransfer.H` | Forward heat-transfer solve |
| `src/updateMaterialProperties.H` | Brinkman / thermal material interpolation |
| `src/filter_x.H` | Design filtering / projection-side path |
| `src/filter_chainrule.H` | Gradient back-propagation through design regularization |
| `src/AdjHeatTransfer.H` | Thermal adjoint / thermal sensitivity path |
| `src/sensitivity.H` | Main gradient assembly; do not change casually after validation |
| `src/solveDiscreteFlowAdjoint.H` | Large Stage-B diagnostic/oracle implementation; explicit CSR, matrix-free tests, many probes |
| `src/solveDiscreteFlowAdjointProduction.H` | Compact **production** reduced discrete frozen-flow adjoint; current BFINAL-009/010 focus |
| `src/rxPressureRowTranspose.H` | Pressure/continuity-row design derivative contribution established in BFINAL-005 |
| `src/stageB5BoundaryRelaxOracle.H` | Historical/diagnostic Stage-B oracle |
| `src/stageB6RxDesignOracle.H` | R_x / tangent diagnostic support |
| `src/stageB7PressureBCDiagnostic.H` | Pressure-boundary diagnostic |
| `src/stageB8JPPActualResidualFD.H` | Actual frozen-primal pressure residual/J_PP FD probe |
| `src/validateStageB2GradientAmplitude.H` | Stage-B directional FD amplitude gate |
| `src/validateStageBRepeatability.H` | Repeatability/noise validation |
| `src/validateFrozenGradient.H` | Frozen-gradient validation orchestration |
| `src/validateMmaUnlockGate.H` | Hard MMA safety gate |
| `src/evaluateCandidate.H` | Candidate design evaluation |
| `src/saveAcceptedState.H`, `restoreAcceptedState.H` | Rollback / accepted-state transaction mechanism |
| `src/saveOptimizerState.H`, `loadOptimizerState.H` | Restart/persistence |
| `src/MMA/` | MMA implementation |
| `src/tests/test_stage_b_safety_gates.py` | Static Stage-B safety regressions |
| `src/validation/stage_b/` | Human-readable Stage-B reports and current task specs |
| `evidence/agent-group/BFINAL-*` | Detailed per-round evidence, logs, independent recomputations |

### Important distinction: diagnostic vs production adjoint

`solveDiscreteFlowAdjoint.H` is intentionally large. It is the Stage-B diagnostic/oracle environment used to inspect discrete Jacobian blocks, export matrices, perform dot tests and compare against finite differences.

`solveDiscreteFlowAdjointProduction.H` is the path that must eventually run every optimization iteration. It should remain substantially smaller and scalable. A new agent must **not** copy all diagnostic machinery into production merely to make a test convenient.

---

# 5. Stage A — anisotropic thermal-conduction path

The Stage-A validation design is documented in `src/ANISOTROPIC_STAGE_A_VALIDATION.md` and covers:

- algebraic material-model checks;
- isotropic regression;
- 1D x/y/z conduction;
- layered thermal resistance;
- physical directionality;
- boundary sensitivity FD;
- full thermal-adjoint FD;
- serial/parallel consistency.

Current project policy treats the **axis-aligned diagonal anisotropic conductivity path** as the validated production scope. The validation document explicitly keeps off-diagonal tensor components disabled unless a separate rotated-tensor validation is performed.

Do not broaden anisotropy to a general rotated tensor while Stage B is still open.

A previously completed Stage-A retest established strong serial/parallel consistency for the tested gradient directions. Preserve this layer unless a later regression specifically points to it.

---

# 6. Stage B — why it exists

Stage B is not an attempt to prove a universal adjoint theorem. It exists to answer the practical question:

> **For the exact frozen-turbulence state definition used by the optimizer, do the objective, pressure-drop and volume gradients have repeatable, low-noise, correct-direction and useful-amplitude agreement with re-converged finite differences?**

The pressure-drop path became the difficult part because a hand-built reduced SIMPLE/RANS discrete operator has to match OpenFOAM's actual pressure/velocity residual semantics closely enough for the adjoint gradient.

Key historical reports live in:

- `src/validation/stage_b/STAGE_B_B0_AUDIT_REPORT.md`
- `STAGE_B_B0_1_OBJECTIVE_FLUX_AUDIT.md`
- `STAGE_B_B0_2_THERMAL_FLOW_TRANSPOSE_REPORT.md`
- `STAGE_B_B1_REPEATABILITY_NOISE_REPORT.md`
- `STAGE_B_B2B3_GRADIENT_AMPLITUDE_REPORT.md`

However, those older reports predate several BFINAL corrections and are **historical evidence**, not current numerical truth.

The authoritative detailed chronology is under `evidence/agent-group/BFINAL-001` through `BFINAL-008`, plus the current `src/validation/stage_b/BFINAL_009_PRODUCTION_PRECONDITIONER.md`.

---

# 7. BFINAL chronology: what has already been learned/fixed

The following summary reflects the development history captured in the repo.

## BFINAL-001 — diagnostic inconsistency exposed

The existing B4 probe suite was using inconsistent residual conventions. `J`, residual derivatives, RHS/sign and assembly could not yet be treated as one closed system.

## BFINAL-002 — P<-U relaxation semantics identified

The pressure-row velocity tangent was using the wrong unrelaxed `rAUAdj` semantics. The correct path had to follow production relaxed-SIMPLE behavior.

## BFINAL-003 — `J_PU` closed

The P<-U tangent was patched to the actual production relaxed-SIMPLE mapping. This reduced the major P-row discrepancy from order unity to approximately the `1e-4` numerical derivative floor.

**Do not casually rewrite this block.**

## BFINAL-004 — missing pressure-row design derivative identified

The fixed-state pressure residual derivative with respect to design, `R_P,x`, is nonzero. A momentum-only pressure-drop sensitivity was incomplete.

## BFINAL-005 — `R_x` closed

The pressure/continuity-row contribution `-lambda_P^T R_P,x` was added via `rxPressureRowTranspose.H`.

The complete pressure-drop design derivative is therefore not simply the Brinkman/momentum contribution.

**Do not remove this term when simplifying production code.**

## BFINAL-006 — structural pressure null mode discovered

The coupled diagnostic Jacobian contained a non-constant, outlet-localized pressure null mode. This was not merely the normal constant-pressure gauge.

## BFINAL-007 — root cause localized to `J_PP`

The pressure-pressure state derivative omitted the fixedValue outlet boundary `internalCoeffs` contribution and used the wrong interior residual sign. A candidate correction matched actual-residual FD and removed the pathological mode.

## BFINAL-008 — `J_PP` patched and validated

The corrected production/diagnostic operator uses the actual frozen-primal pressure residual convention:

- corrected interior P-P sign;
- fixedValue outlet boundary contribution included;
- consistent forward `J`, matrix-free `J^T`, explicit CSR and production transpose.

The archived BFINAL-008 evidence reports:

- actual pressure-residual FD closure around `1e-4` to `1e-8`, depending on the probe/noise scale;
- transpose consistency at approximately machine precision;
- removal of the old outlet-localized structural singularity;
- tangent solves reaching approximately `1e-14` true relative residual for tested directions;
- no regression in the already-closed `J_PU` and `R_x` layers.

Treat **BFINAL-003 (`J_PU`) + BFINAL-005 (`R_x`) + BFINAL-008 (`J_PP`) as locked mathematical milestones** unless a new independent test actually disproves one of them.

---

# 8. Why BFINAL-009 exists: production scalability

After the discrete operator was substantially closed, a separate production issue remained:

> The production preconditioner setup could not be allowed to scale as `O(N * nnz)` by applying `J^T` once per unknown during topology-optimization iterations.

BFINAL-009 therefore changes **only the production solve layer**. It introduces a sparse pressure-Poisson / GAMG preconditioner in `solveDiscreteFlowAdjointProduction.H` while leaving the validated physical `J^T`, pressure-reference semantics, `R_x` and gradient assembly unchanged.

Current production modes are conceptually:

- `pressureGAMG` — intended production default;
- `diagonal` — rollback/simple mode;
- `bruteForceL1` — diagnostic-only, not scalable production setup.

The current implementation constructs a scaled pressure preconditioner using OpenFOAM `GAMG`, then applies the existing block correction for velocity.

This is intended to close the **P0 scalability problem**: preconditioner construction must be `O(N)` / `O(nnz)`-class rather than `O(N*nnz)`.

---

# 9. Latest runtime status: BFINAL-009 STOPPED at first checkpoint

The latest user-run validation was performed after updating to HEAD `4c49903`.

Static checks/build passed:

```text
python3 src/tests/test_stage_b_safety_gates.py
-> Ran 5 tests ... OK

wclean; wmake
-> WMAKE_EXIT=0
```

Latest local validation case:

```text
/home/ys/dsH/b8_verify_gamg
```

It was copied from the previous production verification case and only the preconditioner controls were changed to:

```text
discreteProdPreconditionerSetup       pressureGAMG;
discreteProdPressurePrecTolerance     1e-3;
discreteProdPressurePrecMaxIter       50;
discreteProdSolverType                fgmres;
discreteFlowAdjointRestart            80;
discreteFlowAdjointMaxIter            4000;
discreteFlowAdjointTolerance          1e-9;
discreteProdEnableRitzPilot           false;
discreteProdConvergeFatal             false;
mmaUpdateEnabled                      false;
frozenGradientValidated               false;
```

The run printed:

```text
PRODPRECSETUP (thermalCoupling): mode=pressureGAMG avg=27.7275430014 min=1.00309714066 max=725.954987504
PRODPRECGAMGSETUP (thermalCoupling): diagRelL2=0.0408489769827 innerTol=0.001 innerMaxIter=50
```

The BFINAL-009 task specification required stopping if this setup mismatch exceeded `1e-8`. The process was therefore stopped before accepting FGMRES behavior.

**This stop was correct.** It prevented an unverified preconditioner from being declared successful merely because an inner GAMG solve or an imported direct solution converged.

Local log supplied by the user:

```text
/home/ys/dsH/b8_verify_gamg/Log.verify_gamg.txt
```

This local runtime result may not yet be archived under `evidence/` in the repo. A new agent should preserve/commit the raw log when continuing BFINAL-010/011 evidence.

---

# 10. Critical correction to the first BFINAL-009 root-cause hypothesis

The first diagnosis suggested that OpenFOAM `interpolationSchemes default linear` uses a fixed arithmetic `0.5/0.5` average whereas the production block uses `mesh.weights()`.

**Do not proceed from that assumption. It is not generally correct in OpenFOAM 7.**

OpenFOAM 7's `linear` surface interpolation obtains its weights from `mesh().surfaceInterpolation::weights()`. Therefore, the internal-face interpolation may already be consistent with the production `mesh.weights()` path.

There is a more immediate issue in the BFINAL-009 diagnostic itself:

- `prodJacobiDiag` / `prodPDiagPhysical` includes the fixedValue-pressure boundary diagonal contribution explicitly;
- `PRODPRECGAMGSETUP` currently compares that against raw `prodPressurePrecMatrix.diag()`;
- OpenFOAM's scalar matrix solver temporarily calls `addBoundaryDiag(diag(), 0)` before creating/solving the LDU system, which adds patch `internalCoeffs` into the effective solver diagonal.

Therefore the reported `diagRelL2 ~= 0.04085` may be a **false negative caused by comparing a raw fvMatrix diagonal with an effective solver diagonal that includes boundary `internalCoeffs`**.

This is currently a strong diagnosis, **not yet a validated result**. BFINAL-010 must test it directly.

Do not "fix" GAMG interpolation or alter the validated `J^T` before this diagnostic is resolved.

---

# 11. Immediate next task for the incoming Agent: BFINAL-010

## Goal

Strengthen the pressure-GAMG matrix-equivalence diagnostic **without changing the physical production operator or the mathematical pressure-preconditioner definition**.

## Hard scope restrictions

For BFINAL-010, do **not**:

- modify validated production `J^T`;
- modify `J_PU`;
- modify `J_PP`;
- modify `R_x` / `rxPressureRowTranspose.H`;
- modify objective/constraint definitions;
- modify filter/projection;
- modify MMA;
- tune FGMRES to hide the setup discrepancy;
- add a custom interpolation scheme merely because `diagRelL2` was 0.041;
- declare success from GAMG's own inner residual.

## Required diagnostic decomposition

Before the first outer Krylov iteration, compare the production pressure block and preconditioner by separating:

```text
interiorDiagRelL2
boundaryDiagRelL2
effectiveDiagRelL2
offDiagRelL2
```

The effective preconditioner diagonal must reproduce the same boundary treatment that OpenFOAM's scalar solver actually uses (`raw diag + boundary internalCoeffs` in solver semantics).

At minimum verify:

1. internal face coefficient correspondence;
2. fixedValue pressure boundary coefficient correspondence;
3. effective full diagonal correspondence;
4. internal off-diagonal correspondence (`upper/lower` vs expected pressure-face coefficient/sign/scaling).

## First acceptance checkpoint

Prefer:

```text
effectiveDiagRelL2 < 1e-10
offDiagRelL2       < 1e-10
```

Hard acceptance ceiling for continuing the outer solve:

```text
effectiveDiagRelL2 < 1e-8
offDiagRelL2       < 1e-8
```

If the gate fails, **stop the solve and localize the coefficient mismatch**. Do not proceed to a long FGMRES run.

A useful log block would be:

```text
PRODPRECGAMGCHECK (...):
    interiorDiagRelL2=...
    boundaryDiagRelL2=...
    effectiveDiagRelL2=...
    offDiagRelL2=...
```

Add a static regression test if practical so future refactors cannot silently revert to comparing incompatible diagonal definitions.

---

# 12. After BFINAL-010: BFINAL-011

Only after the pressure-GAMG matrix-equivalence checkpoint passes should the agent run the production outer solve.

Use the initial controls already specified by BFINAL-009 unless evidence requires otherwise:

```text
pressureGAMG inner tolerance = 1e-3
pressureGAMG inner maxIter   = 50
FGMRES restart               = 80
outer maxIter                = 4000
outer tolerance              = 1e-9
Ritz pilot                   = false
```

Acceptance must use the **true outer residual** of the physical matrix-free production operator:

```text
||b - J^T lambda|| / ||b||
```

Do not use these as substitutes for the outer true residual:

- GAMG's inner `Final residual`;
- an imported SuperLU solution;
- a favorable gradient sign without convergence.

Also require:

- `gradProxy` / checkpoint behavior to stabilize as the solve converges;
- thermal-coupling and pressure-drop adjoints both to converge under the configured production criteria;
- no regression in BFINAL-008 operator dot tests / explicit-vs-matrix-free checks;
- final derivative outputs consistent with the previously validated operator, within the relevant numerical floor.

If this succeeds, the **production preconditioner P0** can be considered substantially closed.

---

# 13. Then BFINAL-012: final current-source gradient gate

A scalable linear solve is not enough to unlock topology updates.

On the exact current source used for optimization, perform the final directional finite-difference validation for:

1. heat-transfer objective;
2. cold-side pressure-drop constraint;
3. solid-volume constraint.

The perturbation must be applied to the **raw design variable**, and each `x +/- h d` case must re-run the same filter/projection/material/frozen-primal chain.

Use multiple deterministic directions and an epsilon ladder to establish a stable FD plateau.

For the difficult turbulent flow-dependent gradients, prioritize:

- matching sign;
- stable direction;
- repeatability;
- useful amplitude.

Target amplitude error should preferably be below ~5%; a stable 5–10% band may be acceptable for the intended frozen-turbulence method if supported by a clear step-size plateau and consistent directions. Do not accept factor-of-two/factor-of-five errors or sign reversal.

The analytic volume derivative should be much tighter.

Only after this final gate is explicitly passed may:

```text
frozenGradientValidated = true
```

be set for an MMA-enabled case.

`validateMmaUnlockGate.H` is a hard safety feature. Do not add a direction-only or "temporary" bypass.

---

# 14. Stage C onward: stop debugging the Jacobian and optimize

Once BFINAL-012 passes, the project should change character from "adjoint debugging" to "topology optimization".

## C0 — serial MMA smoke

Run a small/representative case for approximately 5–10 MMA iterations to prove the whole loop:

```text
design
 -> filter/projection
 -> primal
 -> frozen turbulence
 -> objectives/constraints
 -> adjoints
 -> gradients
 -> MMA
 -> candidate evaluation / safeguards
 -> accepted design
```

Success means the topology begins to evolve in a physically sensible direction while constraints remain controlled. This is more important to the research goal than continuing to improve already-adequate Jacobian micro-closure.

## C1 — production parallelization

The current `solveDiscreteFlowAdjointProduction.H` intentionally aborts under `Pstream::parRun()` because processor-patch transpose terms are not implemented.

This is the current **P2 engineering limitation**.

Parallel production adjoint is desirable for final large meshes, but it should **not block the first serial MMA smoke once B-final is validated**. First prove that the optimizer can evolve the structure; then implement/validate processor-patch transpose terms without changing the serial mathematics.

## C2 — full optimization runs

After serial workflow validity and adequate production scalability are established:

- run representative topology optimization;
- monitor objective, pressure-drop constraint, solid-volume constraint and design change;
- preserve accepted-state rollback and optimizer state;
- perform periodic full SST refresh/correction.

## D — SST refresh study

Establish a robust policy for how frequently the complete SST primal state is refreshed relative to frozen-turbulence sensitivity/MMA updates.

## E — final validation

For final optimized structures:

- run full SST-RANS primal verification;
- compare baseline vs optimized heat-transfer performance;
- verify pressure drop and solid fraction;
- inspect topology and flow/thermal fields;
- if possible, compare against later experimental results.

---

# 15. Current open problems, ranked

## P0 — Production adjoint linear-solver scalability / robustness

**Status:** active; BFINAL-009/010/011.

The old brute-force preconditioner setup is unacceptable in topology iterations. `pressureGAMG` is the intended scalable route. The immediate issue is determining whether the current `diagRelL2=0.04085` is a genuine coefficient mismatch or a raw-vs-effective-diagonal diagnostic mistake.

**Do this now.**

## P1 — Final current-source FD amplitude gate

**Status:** blocked behind a trustworthy production solve.

Even though individual discrete blocks were closed in BFINAL-003/005/008, the source used for MMA still needs final end-to-end directional FD evidence.

**Do this after P0.**

## P2 — Production flow adjoint is serial-only

**Status:** known limitation.

Processor-patch transpose terms are not implemented in `solveDiscreteFlowAdjointProduction.H`.

**Do this after the first validated serial MMA smoke unless mesh size makes it unavoidable earlier.**

## P3 — Optimization-loop behavior not yet proven on the final frozen-turbulence production path

**Status:** MMA deliberately locked.

After B-final, run a short guarded serial optimization before any large production campaign.

## P4 — Documentation drift

Several older README/report statements describe obsolete commits and solver behavior. This handoff document is intended to correct that. When a new BFINAL stage closes, update this file or add a short dated status note rather than letting the next agent infer state from old reports.

---

# 16. Safety invariants and development rules

A new AI agent must preserve the following unless the user explicitly authorizes a change.

1. **MMA remains disabled while `frozenGradientValidated=false`.**
2. Do not add any direction-only bypass to `validateMmaUnlockGate.H`.
3. A solve-layer change must not silently alter validated `J^T` mathematics.
4. A preconditioner is allowed to approximate the inverse, but it must not redefine the physical operator or residual used for the outer true-residual test.
5. Imported SuperLU/direct solutions are diagnostic oracles only; they are not proof that the production iterative solver works.
6. Checkpoint failures are stop conditions. Do not continue expensive runs and rationalize the failed gate afterward.
7. Preserve BFINAL-003 `J_PU`, BFINAL-005 `R_x`, and BFINAL-008 `J_PP` unless new independent evidence demonstrates an actual regression.
8. Keep diagnostic/oracle code separate from production code whenever possible.
9. Do not widen the physics scope to full SST adjoints.
10. Do not tune tolerances to make a failed mathematical-equivalence check disappear.
11. Every major BFINAL change should preserve build logs, run logs, exact controls, commit SHA and independent recomputation evidence.
12. Before editing any sensitive layer, identify which prior validation gate it can invalidate and plan the necessary regression tests.

---

# 17. Validation philosophy

There are two different standards in this repo. Do not mix them.

## Local discrete-operator / assembly tests

Examples:

- transpose dot tests;
- matrix-free vs explicit CSR;
- pressure preconditioner coefficient equivalence;
- analytic material interpolation derivative.

These can and should use very tight numerical thresholds when the relationship is algebraic, often `1e-8` to machine precision.

## End-to-end turbulent topology-gradient tests

These compare an adjoint derivative against re-converged finite differences of the frozen/full workflow. The meaningful requirements are:

- repeatability above the numerical noise floor;
- correct sign;
- direction cosine close to 1 where a vector comparison is meaningful;
- stable amplitude across an epsilon plateau;
- practically acceptable amplitude error (prefer <5%; stable <10% may be defensible for the deliberately approximate frozen-turbulence method).

The project must **not** spend unlimited time trying to force end-to-end turbulent gradients to `1e-10` if they are already stable and sufficiently accurate to drive the topology in the correct direction.

---

# 18. Recommended reading order for a new Agent

Read in this order before proposing a large patch:

1. **`AI_AGENT_HANDOFF.md`** — this document.
2. `src/MTO_HF.C` — understand include/control flow.
3. `src/README.md` and `src/SRC_CODE_GUIDE.md` — architecture context, noting that some descriptions may predate BFINAL work.
4. `src/ANISOTROPIC_STAGE_A_VALIDATION.md` — understand the validated anisotropic production boundary.
5. `src/FROZEN_HOT_REGION_IMPLEMENTATION.md` — understand hot-side freeze/live-temperature semantics.
6. `src/validation/stage_b/STAGE_B_BFINAL_SUMMARY.md` — BFINAL chronology, but check its recorded commit against current HEAD.
7. `evidence/agent-group/BFINAL-003/FINAL_REPORT.md` — `J_PU` closure.
8. `evidence/agent-group/BFINAL-005/FINAL_REPORT.md` — `R_x` closure.
9. `evidence/agent-group/BFINAL-008/FINAL_REPORT.md` and `POST_REVIEWER_FINAL_REPORT.md` — `J_PP` closure and independent validation.
10. `src/validation/stage_b/BFINAL_009_PRODUCTION_PRECONDITIONER.md` — current production preconditioner intent.
11. `src/solveDiscreteFlowAdjointProduction.H` — current production implementation.
12. Only then inspect `src/solveDiscreteFlowAdjoint.H` and Stage-B diagnostic oracles for details required by the current task.
13. `src/tests/test_stage_b_safety_gates.py` and `src/validateMmaUnlockGate.H` before changing any gates.

For BFINAL-010 specifically, also inspect OpenFOAM 7's implementation of:

- `linear` surface interpolation weights;
- `fvScalarMatrix` solver construction/solve;
- `fvMatrix::addBoundaryDiag`;
- laplacian scheme coefficient assembly;
- the actual pressure boundary types in the validation case.

Do not reason about these from memory if the source can be inspected directly.

---

# 19. Local development / validation locations seen in the current work

Recent work has used:

```text
repository/workspace:
    /home/ys/dsH/TO-ANISOTROPIC

latest BFINAL-009 local case:
    /home/ys/dsH/b8_verify_gamg

latest reported log:
    /home/ys/dsH/b8_verify_gamg/Log.verify_gamg.txt
```

Historical evidence also references scratch cases around `b2_case_smoke` / `b8_scratch_*`.

These paths are environment-specific and may not exist on another machine. Never hard-code them into production mathematics.

Build environment is OpenFOAM Foundation v7. Rebuild after source changes; preserve the full build log and verify the binary is newer than the changed sources.

---

# 20. Definition of done for the whole project

The project is not done when a Jacobian probe reaches machine precision. It is done when the following chain is demonstrated:

1. anisotropic material/thermal path is validated;
2. frozen-turbulence objective, pressure-drop and volume gradients are repeatable and FD-validated strongly enough to drive optimization;
3. production adjoint solve is computationally usable at topology-optimization scales;
4. guarded MMA can run without bypassing validation gates;
5. optimization produces a physically interpretable cold-side fluid-solid topology and anisotropic heat-conduction structure;
6. heat-transfer performance improves relative to baseline while pressure-drop and solid-volume constraints are respected;
7. periodic/full SST primal refresh does not invalidate the improvement trend;
8. final topology is re-verified with full SST-RANS and, where available, later experimental comparison.

A scientifically successful result can therefore use a **frozen-turbulence adjoint approximation** as long as its limitations are stated and the optimization direction/performance are validated. A full turbulence-model adjoint is explicitly outside the required endpoint.

---

# 21. Immediate assignment template for the incoming Agent

If no newer handoff/status commit exists, the new agent's first implementation task should be:

> **BFINAL-010 — pressure-GAMG effective-matrix equivalence diagnostic**
>
> Starting from branch `agent/dsH-stage-b-validation` at or after `4c49903`, inspect the current production pressure-GAMG assembly and OpenFOAM 7 scalar-matrix solver semantics. Do not modify the validated physical `J^T`, `J_PU`, `J_PP`, `R_x`, objective, filter/projection or MMA. Replace/extend the current raw-diagonal-only setup check with a decomposition that independently reports internal diagonal, boundary `internalCoeffs`, effective solver diagonal and off-diagonal closure against the production P-P block. Confirm from OpenFOAM 7 source that `linear` interpolation uses mesh surface-interpolation weights rather than assuming fixed 0.5/0.5. Rebuild and run the first checkpoint only. Continue to the outer FGMRES solve only if `effectiveDiagRelL2 < 1e-8` and `offDiagRelL2 < 1e-8` (prefer <1e-10). If the checkpoint fails, stop and report localized coefficient evidence instead of tuning the solver.

After that:

```text
BFINAL-010 matrix-equivalence gate
    -> BFINAL-011 GAMG+FGMRES true-residual convergence/regression
    -> BFINAL-012 current-source FD amplitude gate
    -> C0 guarded serial MMA smoke
    -> C1 production parallel adjoint
    -> C2 full optimization
    -> D periodic SST refresh
    -> E final full-SST validation
```

---

# 22. Final note to future Agents

This codebase contains many sophisticated diagnostic experiments because Stage B had to localize several subtle SIMPLE/pressure-discretization errors. Those experiments were useful, but they are **not the final product**.

When deciding what to do next, ask:

> Does this change help us obtain a validated, scalable frozen-turbulence gradient and then optimize the topology?

If the answer is no — especially if the change instead expands toward a universal exact SST adjoint — it is probably outside the intended project path.
