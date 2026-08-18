# BFINAL-009 — Production pressure-GAMG preconditioner

## Scope

This change affects only the production solve layer in
`solveDiscreteFlowAdjointProduction.H`. It does not change the validated
matrix-free `J^T`, the pressure-reference semantics, `R_x`, or gradient
assembly.

The former default divided the pressure block by its diagonal. The new
`pressureGAMG` mode approximately inverts the scaled pressure Poisson block
with OpenFOAM GAMG, then applies the existing block upper-triangular U
correction. Both setup and each multigrid cycle are sparse/near-linear in mesh
size. `diagonal` remains available as a rollback mode and
`bruteForceL1` remains diagnostic-only.

## Required case controls

```
discreteProdPreconditionerSetup       pressureGAMG;
discreteProdPressurePrecTolerance     1e-3;
discreteProdPressurePrecMaxIter       50;
discreteProdEnableRitzPilot           false;
discreteProdSolverType                fgmres;
discreteProdConvergeFatal             false;
mmaUpdateEnabled                      false;
frozenGradientValidated               false;
```

Keep the existing outer controls for the first run:
`discreteFlowAdjointRestart 80`,
`discreteFlowAdjointMaxIter 4000`, and
`discreteFlowAdjointTolerance 1e-9`.

## Acceptance sequence

1. Rebuild with OpenFOAM 7 and run
   `python3 src/tests/test_stage_b_safety_gates.py`.
2. Run the same serial `b8_verify_prod` case used at source HEAD
   `6320943`, changing only
   `discreteProdPreconditionerSetup pressureGAMG`.
3. Require `PRODPRECGAMGSETUP ... diagRelL2` to be finite and preferably
   below `1e-10`. If it is larger than `1e-8`, stop and report the log.
4. Compare thermalCoupling and pressureDrop convergence against the recorded
   diagonal baselines `5.08309124088e-05` and `1.16060277262e-02`.
5. The production gate passes only if both true relative residuals reach the
   configured outer tolerance without using an imported solution.
6. After convergence, rerun the P4 dot test and P5 explicit/matrix-free test
   to confirm that the solve-layer change did not alter `J^T`.

Do not enable MMA from this preconditioner test. T2/T3 and the current-HEAD
B2/B3 finite-difference amplitude gates remain required.
