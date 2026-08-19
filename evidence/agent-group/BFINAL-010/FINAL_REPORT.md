# BFINAL-010 FINAL REPORT — pressure-GAMG four-group matrix equivalence decomposition + hard gate

## Provenance

- Repo: `/home/ys/dsH/TO-ANISOTROPIC`, branch `agent/dsH-stage-b-validation`
- Base HEAD: `6a31ed1` (BFINAL-008/009 evidence commit) + this task's working changes (committed together with this report; see git log for the exact SHA of the commit that adds this file)
- Binary: `build/bin/MTO_HF` built 2026-08-19 18:34 (wmake exit 0; first attempt failed on two misspelled gate identifiers `effectiveDiagRel2`/`offDiagRel2`, fixed to `effectiveDiagRelL2`/`offDiagRelL2`, static tests updated accordingly; build log `wmake_bfinal010.log` at repo root)
- Case: `/home/ys/dsH/b10_verify_gamg` (`cp -a` of `b8_verify_gamg`, old logs/checkpoints removed) — production path, `discreteProdPreconditionerSetup pressureGAMG`, `discreteProdSolverType fgmres`, `discreteProdEnableRitzPilot false`, `discreteProdConvergeFatal false`, `restart 80`, `maxIter 4000`, tolerance 1e-9
- Run: 2026-08-19 18:34:36 (PID 15478, host DESKTOP-DSG34PK), log `Log.verify_b010.txt` archived here verbatim (sha256 `ff197ded03ae1fdb5ad8500d32880d0e2803902312cb86407c21e574f018e131`)

## What was implemented (src/solveDiscreteFlowAdjointProduction.H)

1. **laplacianSchemes guard** (before `tProdPressurePrecMatrix` construction, pressureGAMG mode only). Fact correction vs. the earlier handoff wording: the coefficient face interpolation of `fvm::laplacian(volScalarField, …)` does NOT read `interpolationSchemes`; the interpolation word comes from the `laplacianSchemes` entry ITstream (`Gauss <interp> <snGrad>`: `Gauss` is consumed by `laplacianScheme::New`, the next word selects the mobility face interpolation — verified against OF7 `gaussLaplacianScheme.H`, `laplacianScheme.H:120-135`, `laplacianSchemes.C`). The guard walks `default` plus any `laplacian(prodPressurePrecMobility…` named entry, extracts that word, and FatalErrors unless it is `linear` (the J^T kf assumption `linear == mesh.weights()`). Pass marker: `PRODPRECGAMGSCHEME`.
2. **Four-group equivalence decomposition + hard gate** (replaces the old bare `diagRelL2` check). Operator side is independently recomputed from `primalPressureMobility`/`prodWeights`/`deltaCoeffs`/`magSf` times `prodPressurePrecScale`; preconditioner side reads `diag()` (bare, negSumDiag), patch-accumulated `internalCoeffs()` (addBoundaryDiag semantics), `upper()`/`lower()`. Four global relL2 values printed as `PRODPRECGAMGCHECK`; FatalError on any non-finite value or `effectiveDiagRelL2 > 1e-8` or `offDiagRelL2 > 1e-8` (hardcoded safety invariant, not configurable). `PRODPRECGAMGSETUP` line retained without the misleading `diagRelL2=` term.
3. **Static regression**: `src/tests/test_stage_b_safety_gates.py::test_pressure_gamg_equivalence_decomposition_gate` asserts the four metric names, `PRODPRECGAMGCHECK`, the FatalError fragment, the `laplacianSchemes` guard marker, the internalCoeffs accumulation marker, and that the old bare comparison string `prodPressurePrecMatrix.diag()[celli] - refDiag` cannot return. 6 tests OK.

## Gate results (both adjoint labels, identical to all printed digits)

```
PRODPRECGAMGSCHEME (thermalCoupling): laplacianSchemes coefficient interpolation 'linear' == mesh.weights()
PRODPRECGAMGCHECK (thermalCoupling): interiorDiagRelL2=1.13142348768e-16 boundaryDiagRelL2=8.25782114104e-17 effectiveDiagRelL2=1.42213008543e-16 offDiagRelL2=1.34302013134e-16
PRODPRECGAMGSCHEME (pressureDrop):    laplacianSchemes coefficient interpolation 'linear' == mesh.weights()
PRODPRECGAMGCHECK (pressureDrop):      interiorDiagRelL2=1.13142348768e-16 boundaryDiagRelL2=8.25782114104e-17 effectiveDiagRelL2=1.42213008543e-16 offDiagRelL2=1.34302013134e-16
```

**GATE PASS at machine precision.** The pressure-GAMG preconditioner matrix (`-fvm::laplacian` of the scaled mobility) reproduces the scaled J^T P–P block term by term: interior diag +kf, fixedValue-p boundary diag +mob_c·δ_b·|Sf|_b, upper/lower −kf, and interior+boundary == `prodJacobiDiag[P]`.

**Consequence:** BFINAL-009's `PRODPRECGAMGSETUP diagRelL2=0.0408 > 1e-8` stop is confirmed as a comparison artifact — it compared the boundary-less fvMatrix diagonal against the boundary-inclusive effective diagonal (addBoundaryDiag semantics). The operator and the preconditioner were never inconsistent.

## Solver-layer run results

### thermalCoupling — FIRST-EVER production-path convergence

```
iter=80  trueRelRes=3.0447e-4      GRADPROXY=1348.55   PRODDEFL cos=0.2199
iter=160 trueRelRes=1.2072e-4      GRADPROXY=1934.84   cos=0.1538
iter=240 trueRelRes=4.9182e-5      GRADPROXY=2143.99   cos=0.1416
iter=320 trueRelRes=1.9048e-5      GRADPROXY=2207.37   cos=0.1027
iter=400 trueRelRes=5.4579e-6      GRADPROXY=2216.84   cos=0.1456
iter=480 trueRelRes=1.4277e-6      GRADPROXY=2217.38   cos=0.1191
iter=560 trueRelRes=2.8050e-7      GRADPROXY=2217.45   cos=0.1566
iter=640 trueRelRes=2.0681e-8      GRADPROXY=2217.48
iter=720 trueRelRes=2.5603e-9      GRADPROXY=2217.48
iter=751 trueRelRes=9.48919598855e-10  <= 1e-9  CONVERGED
PRODRESID (thermalCoupling): |rU|L2=6.43e-09 |rP|L2=1.27e-09 (finite, no NaN)
GRADPROXY-FINAL = 2217.48004323 (plateaued)
```

- 751 FGMRES iterations, one GAMG apply each (5–8 V-cycle iterations, inner final residual ~4e-4–1e-3 at innerTol=1e-3), whole run wall time ≈ 1–2 min.
- The ~5e-5 stall floor observed with the `diagonal` preconditioner (BFINAL-008 fresh run: 5.083e-5 after **4000** iterations) is **broken** — it was a preconditioning artifact, not an operator defect. This is consistent with (and now explained by) the gate result: the GAMG matrix is exactly the scaled P–P block.
- BFINAL-009's partial 240-iteration trend reproduces (4.9066e-5 → 4.9182e-5 @240, ~0.2% difference from run-to-run GAMG ordering effects), confirming that log was healthy partial data.

### pressureDrop — GAMG inner solve NaN (new blocker, BFINAL-011 scope)

```
PRODDIAG (pressureDrop): jacobiUavg=1 jacobiPavg=103.910 jacobiMin=0.003097 jacobiMax=724.955 jacobiNearZero=0   (sane)
GAMG: Solving for prodPressurePrecPsi_pressureDrop, Initial residual = 1, Final residual = nan, No Iterations 50   <- first apply diverges
GAMG: Solving for prodPressurePrecPsi_pressureDrop, Initial residual = 0, Final residual = 0, No Iterations 0
[FGMRES-PROD pressureDrop] iter=2 trueRelRes=-nan
Production reduced discrete flow adjoint pressureDrop: FGMRES iterations=2, true relative residual=-nan
--> Warning: did not converge (DIAGNOSTIC MODE) ... Continuing for cross-objective discrimination.
```

- The divergence happens **inside the first GAMG inner solve** (50 V-cycle iterations → nan), before any FGMRES restart logic matters. The preconditioner matrix is identical to thermalCoupling's (gate values identical to all digits; same mobility/scale), and the same matrix solved fine 751 times for thermalCoupling — so the trigger is the pressureDrop rhs/label interaction with the GAMG solver, not the assembled operator.
- Downstream: `gsensPressureDrop` solve residuals nan → the existing MMA safety gate fires (`sensitivity.H:263`, "Non-finite MMA data at active design index 0, cell 5610") → clean abort, MTO_RC=134. This is the pre-existing safety architecture working as designed; `discreteProdConvergeFatal=false` kept the solver layer in diagnostic mode so the data above is complete and valid.
- Per the task's no-masking rule, no parameter was tuned and no gate was relaxed to push past this; the NaN is reported as the BFINAL-011 entry problem.

## Baseline comparison table (production path, tolerance 1e-9)

| label | precond | iterations | final trueRelRes | status |
|---|---|---|---|---|
| thermalCoupling | diagonal (BFINAL-008 fresh) | 4000 | 5.083e-5 | stalled |
| thermalCoupling | pressureGAMG (BFINAL-009 partial) | 240 (killed) | 4.907e-5 | partial |
| **thermalCoupling** | **pressureGAMG (this run)** | **751** | **9.489e-10** | **CONVERGED** |
| pressureDrop | diagonal (BFINAL-008 fresh) | 4000 | 1.161e-2 | stalled (finite) |
| pressureDrop | pressureGAMG (this run) | 2 | −nan | GAMG inner diverges (BFINAL-011) |

## Files

| file | note |
|---|---|
| `Log.verify_b010.txt` | full run log (2431 lines) |
| `adjointCheckpoint_thermalCoupling.tsv` | 80→720 checkpoint series (final 751 in log) |
| `adjointCheckpoint_pressureDrop.tsv` | 2 nan rows (documentational) |
| `../../../src/solveDiscreteFlowAdjointProduction.H` | guard + decomposition + hard gate |
| `../../../src/tests/test_stage_b_safety_gates.py` | static regression (6 tests OK) |
| `../../../wmake_bfinal010.log` (repo root) | successful rebuild log |

## Verdict

- **PASS (equivalence gate)**: the pressure-GAMG preconditioner matrix is the scaled J^T P–P block to machine precision on both labels; the BFINAL-009 stop is retired as a comparison artifact.
- **Operational consequence**: thermalCoupling now converges in the production iterative path (751 iters, 9.49e-10); the exact/oracle route is no longer the only trustworthy path for this label.
- **Open problem handed to BFINAL-011**: pressureDrop first GAMG inner solve NaN — solver-layer (inner GAMG robustness for this label's rhs), not matrix equivalence.
