# BFINAL-011 cycle-1 FINAL REPORT — pressureDrop GAMG inner NaN: root cause found and fixed

## Verdict

**PASS — root cause identified and fixed; both adjoint labels now converge in
the production iterative path with pure default GAMG settings.**

```
thermalCoupling: FGMRES 750 iters, trueRelRes = 9.4956818492e-10  <= 1e-9  CONVERGED
pressureDrop:    FGMRES 902 iters, trueRelRes = 9.92818139401e-10 <= 1e-9  CONVERGED
Full program completed (sensitivity.H finished), MTO_RC = 0.
```

This is the first production-path run with BOTH adjoints converged and the
full solver program exiting cleanly.

## Provenance

- Repo: `/home/ys/dsH/TO-ANISOTROPIC`, branch `agent/dsH-stage-b-validation`, base `b01999c` + this task's changes (committed together with this report)
- Binary: rebuilt twice; final `build/bin/MTO_HF` 2026-08-19 20:48 (wmake exit 0, `wmake_bfinal011.log` at repo root)
- Case: `/home/ys/dsH/b11_diag` (`cp -a b10_verify_gamg`, production path, `discreteProdPreconditionerSetup pressureGAMG`, `discreteProdSolverType fgmres`, tol 1e-9, restart 80, maxIter 4000, `discreteProdConvergeFatal false`)
- Winning run: `Log.verify_b011_E0p.txt` (2026-08-19 20:55–20:57, ~120 s wall, sha256 `33f4b855e181168fae52dfb13830e0d8da55d7ee63caab31e24e77f29e4113c8`)

## Root cause

The pressureDrop inner-GAMG NaN observed in BFINAL-010 was **not** an
inherent GAMG weakness on the pressureDrop rhs class, and **not** cross-label
state leakage. It was a matrix-storage side effect introduced by the
BFINAL-010 equivalence gate itself:

1. `-fvm::laplacian(mobility, psi)` assembles the preconditioner matrix in
   **symmetric lduMatrix storage** (only `upper` allocated; `lower()` is the
   `upper` alias). `matrix.symmetric() == true`.
2. The BFINAL-010 gate read `prodPressurePrecMatrix.upper()` / `.lower()` on
   the **non-const** matrix. OF7's non-const `lduMatrix::lower()`
   (lduMatrix.C:168-182) **materialises the lower array as a separate copy**,
   converting the matrix to asymmetric storage — silently changing the solver
   path: `lduMatrix::solver::New` then treats PCG as invalid ("Unknown
   asymmetric matrix solver PCG") and `GAMGSolver` takes its asymmetric
   branch with `scaleCorrection_ = matrix.symmetric() = false`
   (GAMGSolver.C:76).
3. On that corrupted path GAMG diverged on the pressureDrop rhs (all 33600
   cells non-finite in the first inner solve) while still converging
   thermalCoupling — which is why the b010 run looked like a
   label-specific solver robustness problem.

**Fix:** read the matrix through a `const lduMatrix&` alias in the gate
(const `lower()` returns the upper alias without materialising —
lduMatrix.C:262-281). Gate values are unchanged (bit-identical four-group
relL2); the solver path is unchanged (symmetric GAMG, scaleCorrection=true).
A static regression assertion forbids the materialising access pattern.

## Experiment matrix and decision table

| run | binary / matrix state | config | result | conclusion |
|---|---|---|---|---|
| E0 (pre-fix) | b011 build #1, lower() materialised (asymmetric GAMG path) | defaults | gate values bit-identical to b010; TC converged 750 iters 9.58e-10; PD fast-FatalError at preconditioner apply 1 (33600/33600 non-finite cells) | baseline reproduced; determinism; fast-fail works |
| E1 (pre-fix) | same | `solveThermalCouplingFlowAdjoint false` (PD alone, no HT anywhere in the log) | identical PD NaN at apply 1; PRODRHS PD sumRhsP=95.1198507991 unchanged | cross-label leakage (H2) excluded under this binary |
| E2 (pre-fix) | same | `discreteProdPressurePrecSolver PCG` | `FOAM FATAL IO ERROR: Unknown asymmetric matrix solver PCG` | **discovery**: the matrix was asymmetric — gate side effect found |
| **E0' (post-fix)** | b011 build #2 (const alias), symmetric storage restored | pure defaults (no new keys) | **TC 750/9.496e-10, PD 902/9.928e-10, both CONVERGED, MTO_RC=0**; TC trajectory bit-for-bit equal to the BFINAL-009 b8 partial run (e.g. @80 3.04473267442e-4, @240 4.90659573665e-5) | **root cause confirmed = storage materialisation; fix = const access; no PCG / no GAMG parameter sweep needed** |

Evidence-loss note: E0/E1 logs were deleted by the run script's original
`rm -f Log*` (fixed to delete only the current run's log); their load-bearing
lines are quoted above and in EXECUTOR_SUMMARY.md. E2's log survives verbatim
(`Log.verify_b011_E2.txt`). E1'/E2' (post-fix TC-skip / PCG variants) were not
run — the decision tree resolved at E0' and the acceptance criteria do not
require them; they remain cheap future checks.

## Acceptance criteria (plan fix-round)

| criterion | result |
|---|---|
| both labels true residual <= 1e-9 | ✓ 9.4957e-10 / 9.9282e-10 |
| GRADPROXY plateau | ✓ TC 2217.4777→2217.48004 (8 digits); PD -205190.7→-205194.141 (converged plateau, see checkpoint tsv) |
| PRODPRECGAMGCHECK gate values unchanged (1e-16 level) | ✓ bit-identical on both labels: interior 1.13142348768e-16, boundary 8.25782114104e-17, effective 1.42213008543e-16, offDiag 1.34302013134e-16 |
| static tests pass | ✓ 7 tests OK |
| no FATAL / clean exit | ✓ MTO_RC=0, sensitivity.H completed |

## Production-path comparison table (tolerance 1e-9)

| label | precond | matrix state | iterations | final trueRelRes | status |
|---|---|---|---|---|---|
| thermalCoupling | diagonal (BFINAL-008 fresh) | n/a | 4000 | 5.083e-5 | stalled |
| thermalCoupling | GAMG (b010/b011-E0, materialised) | asymmetric | 750-751 | 9.5-9.6e-10 | converged (corrupted-path variant) |
| **thermalCoupling** | **GAMG (E0', symmetric)** | symmetric | **750** | **9.4957e-10** | **CONVERGED** |
| pressureDrop | diagonal (BFINAL-008 fresh) | n/a | 4000 | 1.161e-2 | stalled (finite) |
| pressureDrop | GAMG (b010/E0/E1, materialised) | asymmetric | 2 | -nan | inner GAMG divergence (side effect) |
| **pressureDrop** | **GAMG (E0', symmetric)** | symmetric | **902** | **9.9282e-10** | **CONVERGED** |

## What this task implemented (diagnostic + fix, all defaults = baseline)

1. `src/createFields.H` + `src/MTO_HF.C`: per-label isolation switches
   `solveThermalCouplingFlowAdjoint` / `solvePressureDropFlowAdjoint`
   (default true) guarding the `AdjNS_HT.H` / `AdjNS_PD.H` includes.
2. `src/solveDiscreteFlowAdjointProduction.H`:
   - `discreteProdPressurePrecSolver` (word, default `"GAMG"`; `"PCG"` →
     PCG+DIC; else FatalError);
   - GAMG robustness overrides from optProperties (defaults = OF7/current:
     `discreteProdPressurePrecScaleCorrection` true,
     `discreteProdPressurePrecNPreSweeps` 0,
     `discreteProdPressurePrecDirectSolveCoarsest` false);
   - inner-solve psi finiteness fast-fail (`non-finite psi` FatalError with
     label / apply index / hint) replacing the silent nan→MMA-gate path;
   - **the fix**: gate reads the matrix via `const lduMatrix&
     prodPrecLduMatrix` (no storage materialisation);
   - `PRODPRECGAMGSETUP` line now also prints `solver=`.
3. `src/tests/test_stage_b_safety_gates.py`: new
   `test_bfinal011_label_isolation_switches_and_inner_solver_controls`
   (7 tests OK), including an assertion that forbids the materialising
   non-const `.lower()` access pattern.

## Operational notes

- BFINAL-010's FINAL_REPORT statement "pressureDrop NaNs inside its first
  GAMG inner solve — solver-layer problem handed to BFINAL-011" is hereby
  corrected: the trigger was the BFINAL-010 gate's own storage side effect,
  not the GAMG/rhs pair. The BFINAL-010 equivalence RESULT (four-group gate
  at 1e-16) stands unchanged — the materialisation does not alter values,
  only the storage/solver path.
- PCG fallback and the GAMG coefficient overrides remain available as
  robustness insurance; they are NOT needed for the current case.
- Known issue for future parallel work (unchanged): serial matrix-free J^T
  only (see OPENFOAM_USAGE §8.5).

## Files

| file | note |
|---|---|
| `Log.verify_b011_E0p.txt` | winning run (both labels converged, MTO_RC=0) |
| `adjointCheckpoint_E0p_thermalCoupling.tsv` | TC 80→750 checkpoint series |
| `adjointCheckpoint_E0p_pressureDrop.tsv` | PD 80→902 checkpoint series |
| `Log.verify_b011_E2.txt` | pre-fix PCG attempt → asymmetric-solver rejection (side-effect discovery artifact) |
| `optProperties.b11_diag.txt` | case controls snapshot (defaults only) |
| `../../../src/MTO_HF.C`, `../../../src/createFields.H` | isolation switches |
| `../../../src/solveDiscreteFlowAdjointProduction.H` | solver options + fast-fail + const-access fix |
| `../../../src/tests/test_stage_b_safety_gates.py` | static regression (7 OK) |
