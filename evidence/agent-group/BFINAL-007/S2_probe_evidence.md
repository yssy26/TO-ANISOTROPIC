# BFINAL-007 S2 — Actual-pressure-residual central FD probe (Q1 + raw Q2/Q3/Q4b/Q5 data)

- Stage: B-final · Mode: **DIAGNOSTIC_ONLY** · Step: S2-probe-actual-residual-FD
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace, verified) · HEAD `ca8a772b8339a36e7906ea32a7125061c47aa280` · branch `agent/dsH-stage-b-validation`
- Case: `/home/ys/dsH/b2_case_smoke` — runs executed in scratch copies `b7_scratch_{pre,off,on}` (the new switch cannot be added to the real case's `constant/optProperties`, which is files_forbidden; scratch copies are byte-identical case snapshots with only `discreteExplicitMatrixFile` re-pointed into the scratch and `stageB7PressureBCDiagnostic true;` appended in the "on" copy)
- Date: 2026-08-17 (Executor, BFINAL-007)

## 0. Scope (audited)

- New files (files_allowed): `src/stageB7PressureBCDiagnostic.H` (switch-guarded, default off, read-only probe); `src/MTO_HF.C` (+5 lines: guarded include after the AdjNS_PD.H flow-adjoint block).
- `git status`: newly modified tracked file = ONLY `src/MTO_HF.C`. The other 3 modified tracked files (`sensitivity.H`, `solveDiscreteFlowAdjoint.H`, `solveDiscreteFlowAdjointProduction.H`) are the pre-existing locked BFINAL-003/005 changes (untouched). New untracked: `src/stageB7PressureBCDiagnostic.H` + evidence. No forbidden file touched (NS.H, J/J^T, R_x, MMA, filter, case 0/constant/system all clean; real case `constant/optProperties` mtime unchanged).
- No production math/field modified; probe is a complete no-op when the switch is absent/false (verified byte-identical below).

## 1. What the probe does (per approved plan)

1. Read-only reduced-SIMPLE rebuild (rxPressureRowTranspose.H pattern): UEqn -> rAU=1/A() -> HbyA -> phiHbyA=fvc::flux(HbyA)+adjustPhi. Rebuild fidelity vs persisted production rAtU: **relL2=1.88e-09, maxRel=1.52e-08**.
2. Mobility field: internal = persisted `primalPressureMobility` (production rAtU), boundary = rebuilt rAU boundary values.
3. `pEqnDbg = fvm::laplacian(mobilityField, p)` — production-equivalent matrix (same schemes, incl. "corrected" laplacian + fluxRequired face-flux correction).
4. (a) Export outlet internalCoeffs/boundaryCoeffs/deltaCoeffs/magSf/faceCells + full diag/upper/source + addBoundaryDiag diag (+ all-patch boundary file).
5. (b) Read BFINAL-006 null vector `wPrime_TAN_D1.mtx` (134400, plain text), normalize (||w||=2.08447046292e+18, matches BFINAL-006), P block `n_P` (||n_P||=0.99998; n_P[PREF=cell0]=0). Evaluate the ACTUAL frozen-primal pressure residual `R_P(p)=div(phiHbyA)-div(pEqnDbg.flux())` (production fvMatrix::flux() semantics, raw face-flux sum, NO 1/V, boundary owner +) at p±eps·n_P, eps∈{1e-2,1e-4,1e-6}; export central FD **and** the exact `FD_lin = -div(flux(psi=n_P))` (exact linearity; no central-difference roundoff).
6. (c) Export candidate boundary diagonal `+rAtU_c·deltaCoeffs_b·|Sf|_b` per outlet cell (residual-Jacobian convention) + raw internalCoeffs_b (laplacian convention).
7. In-probe self-checks: FD_lin vs -L_prod·n_P (L_prod = addBoundaryDiag diag + upper), outlet mass fraction, FD vs candidate·n_P closure, eps robustness.

## 2. Runs (deterministic serial; PCG+DIC p solver / GMRES adjoint, per case)

| run | binary | switch | scratch | RC | wall |
|---|---|---|---|---|---|
| pre | pre-probe `744d8e74…` | n/a | b7_scratch_pre | 0 | 480 s |
| off | probe `…final` | absent (default false) | b7_scratch_off | 0 | 506.6 s |
| on  | probe `…final` | true | b7_scratch_on | 0 | 505.9 s |

- Environment: `source /opt/openfoam7/etc/bashrc` THEN `unset FOAM_SIGFPE`. (OpenFOAM-7 `Foam::env()` tests variable PRESENCE; bashrc `export FOAM_SIGFPE=` (present-but-empty) still enables SIGFPE trapping → the RxProbe `1.0/A()` division crashed with SIGFPE when the unset was done before sourcing. Verified empirically.)
- All runs wrote byte-identical physics outputs (see §5).

## 3. Q1 — actual frozen-primal pressure/SIMPLE residual central FD along n_P: **NON-ZERO** (outlet-concentrated)

```
FD_lin = -div(flux(psi=n_P))   (production pEqn.flux semantics)
  max|FD| = 3.8984e-10
  L2      = 1.7683e-09
  nnz(|FD|>1e-30) = 33600 (nonzero on all cells; interior content = -L_int·n_P = J_PU·n_U
                           part of the null mode, concentrated at the outlet plane)
  |FD| mass in the 84 outlet cells = 68.6%
  |FD| mass in the top-1% (336 cells) = 94.6%
  mean|FD| outlet / global = 274.5x
```

- The production pressure residual derivative is ~9 orders above the exported-J null floor (||J_exp·n||/||n|| = 2.7e-19): the production frozen-primal pressure mapping does **NOT** share the exported-J null mode.
- FD is eps-independent (exact linearity): |FD(1e-2)−FD(1e-4)|/|FD(1e-4)| = 6.1e-7, |FD(1e-6)−FD(1e-4)|/|FD(1e-4)| = 2.3e-5 (absolute roundoff level; FD values O(1e-10)).
- Central-difference FD matches the exact FD_lin to 6.1e-7 relative (roundoff of the non-orthogonal correction evaluation inside flux()).

## 4. Q2/Q3 structural evidence (raw data for S3)

- **Outlet (patch 1, fixedValue p=0) laplacian boundary coefficients** (from `pEqnDbg.internalCoeffs/boundaryCoeffs`):
  - `internalCoeffs_b = -rAtU_b·|Sf|_b·deltaCoeffs_b`: min **−1.7063e-09**, max **−5.4800e-10**, mean −1.1574e-09 (negative, per gaussLaplacianScheme internalCoeffs = pGamma·gradientInternalCoeffs, fixedValue gradientInternalCoeffs = −δ).
  - `boundaryCoeffs_b = 0` exactly (p_b = 0 → gradientBoundaryCoeffs = δ·p_b = 0).
  - candidate `+rAtU_c·δ·|Sf|` per outlet cell: min +5.4800e-10, max +1.7063e-09; |internalCoeffs_b − (−candidate)|/|candidate| ≤ 1.4e-08 (cell vs boundary-interpolated mobility).
- **All-patch audit**: patch 1 (outlet) is the ONLY patch with nonzero internalCoeffs; the other 7 patches (inlet, hotInlet, hotOutlet, solidEndWalls, bottomWall, topWall, sideWalls — all zeroGradient p) have internalCoeffs = 0 EXACTLY (zeroGradientFvPatchField::gradientInternalCoeffs = 0) and boundaryCoeffs = 0. `|diag_bnd_full − diag_bnd| = 0`.
- **Q3 (preview, full numerical diff is S3)**: exported `explicitJT.mtx` P-P diagonal == probe interior-only laplacian diagonal at **1.70e-16** (excluding the pRef identity row P(0)); production matrix diagonal (interior + addBoundaryDiag) minus exported J P-P diagonal: **0 on all off-outlet cells, exactly internalCoeffs_b (−1.706e-9..−5.480e-10) on the 84 outlet cells**. I.e. the ONLY P-P difference between the production pressure matrix and the exported J is the outlet internalCoeffs_b diagonal (residual-convention: +rAtU_c·δ_b·|Sf|_b per outlet cell).

## 5. Byte-identical checks (probe is a no-op when off)

- `diff -rq pre off` (excluding Logs): ONLY `constant/optProperties` (the sed'd explicitMatrixFile path — harness), `optimization_history.csv`, `optimization_log.dat` differ; the two CSVs differ ONLY in the trailing run-time field (478 vs 506 s). **All physics outputs byte-identical** (time dir `1` fields, explicitJT.mtx, explicitRhs_*, stageB5_*, stageB6_*, rxpr_*, …).
- `diff -rq off on`: same three files (optProperties additionally has the appended switch); probe's own artifacts go to `evidence/agent-group/BFINAL-007/artifacts/` only (none in the case). Log diff off-vs-on = exactly the 8 probe Info lines.
- Log pre-vs-off differs only in Time/PID/Case/ExecutionTime metadata.

## 6. Artifacts (evidence/agent-group/BFINAL-007/artifacts/, sha256 recorded)

| file | content |
|---|---|
| stageB7_outlet_patch.mtx | per outlet face: cell, deltaCoeffs_b, magSf_b, internalCoeffs_b, boundaryCoeffs_b, mobility_b, mobility_cell, candidate |
| stageB7_boundary_all_patches.mtx | all 7840 boundary faces: patch, cell, deltaCoeffs, magSf, ic, bc |
| stageB7_laplacian_diag.mtx / _upper.mtx / _source.mtx | interior-only laplacian diag / upper / source |
| stageB7_laplacian_diag_bnd.mtx (=_bnd_full.mtx) | diag + addBoundaryDiag (all patches; identical to outlet-only) |
| stageB7_laplacian_Lnp.mtx | L_prod·n_P (diag_bnd + upper action) |
| stageB7_nP.mtx | normalized null-vector P block |
| stageB7_RP_full.mtx | R_P(p−eps n), R_P(p), R_P(p+eps n) (full residual incl. div(phiHbyA)) |
| stageB7_RP_divFluxN.mtx | div(flux(psi=n_P)) |
| stageB7_RP_FD.mtx / _FD_lin.mtx / _FD_all_eps.mtx | central FD (1e-4) / exact FD_lin / eps sweep |
| stageB7_candidate_diag.mtx | candidate +rAtU_c·δ·|Sf| per outlet cell (residual convention) |
| stageB7_rebuilt_rAU.mtx | rebuilt 1/UEqn.A() internal (fidelity check) |
| stageB7_divPhiHbyA.mtx | div(rebuilt phiHbyA) |

## 7. Notes for S3 / Post-Reviewer

- The exported FD in this (final, corrected) run is the **true** residual FD `dR_P/dp·n_P = -(L_prod·n_P)` (residual convention). An earlier probe version had a sign bug (FD = +L_prod·n_P) — caught by the in-probe self-check, fixed, and the ON run re-executed; the artifacts above are from the corrected run.
- `|FD_lin − (−L_prod·n_P)|/|L_prod·n_P| = 2.4e-09` — residual of the numerically-evaluated non-orthogonal correction inside flux() (mesh is perfectly orthogonal, correction ≈ 1e-19 absolute; diag/upper reconstruction cannot include it). Does not affect any conclusion (FD dominated by the outlet boundary term + interior null-mode part at O(1e-10)).
- For Q4b (S3): candidate J_PP·n_P (residual convention) = −(L_int + bnd)·n_P; FD = −L_prod·n_P; the two differ only by the exported-J vs probe-L_int agreement (1.7e-16) — closure expected at ~1e-9 (correction residue) or better.
- For Q5 (S3): the boundary term here is the STATE Jacobian term (dR_P/dp). The R_P,x design-row operator (BFINAL-005 rxPressureRowTranspose.H) uses the design field drAU·w whose boundary value is 0 (calculated BC) → dflux_b = 0 there; whether a boundary term enters R_P,x is quantified in S3 using deltaCoeffs_b/magSf_b/drAU exported here and from the case.
- **No production file modified; no tangent entered; stopped after the three verification runs; Post-Reviewer to independently re-verify (see §5 files to re-run).**
