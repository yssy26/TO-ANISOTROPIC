# BFINAL-008 S3 — Probe (run-1) + production (run-2) MTO_HF evidence runs

- Stage: B-final · Mode: PATCH · Step: S3 (run + archive; P1-P8 gate *evaluation* is S4+)
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace, verified) · HEAD
  `ca8a772b8339a36e7906ea32a7125061c47aa280` · branch `agent/dsH-stage-b-validation`
- Binary: `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF` (S1 patch + S2 probe build,
  mtime 2026-08-18 11:27, newer than every src file: solveDiscreteFlowAdjoint.H 09:05,
  solveDiscreteFlowAdjointProduction.H 08:43, stageB8JPPActualResidualFD.H 10:56)
- Env: `source /opt/openfoam7/etc/bashrc` THEN `unset FOAM_SIGFPE`;
  `FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin` (BFINAL-007 documented order)
- Date: 2026-08-18 (Executor, BFINAL-008 S3)

## 0. Scratch cases (NEVER the real case)

| run | scratch | optProperties changes vs b2_case_smoke |
|---|---|---|
| run-1 probe | `/home/ys/dsH/b8_scratch_r1` | `discreteExplicitMatrixFile` -> b8_scratch_r1; `stageB8JPPActualResidualFD true;` appended (stageB4/stageB5/stageB6 already true in the base case) |
| run-2 production | `/home/ys/dsH/b8_scratch_prod` | `stageB4JacobianProbe false`, `stageB5BoundaryRelaxOracle false`, `stageB6RxDesignOracle false`, `discreteExplicitMatrixFile` -> b8_scratch_prod (stageB8 absent -> false) |

Real case `/home/ys/dsH/b2_case_smoke` untouched: `constant/optProperties` mtime
`2026-08-17 00:22:38` unchanged, `stageB8JPPActualResidualFD` count = 0 (grep), sha256
`6d7e4de1ee10c9f4968acc2eecf546dc207d45bc56a526a23b7670a43d680ed7`.

## 1. Run-1 (probe path): MTO_RC=0, 0 crash markers, wall 496.7 s

All load-bearing lines (deterministic; identical to the S2 smoke run, which doubles
as an internal repeatability check):

### P4 — J/J^T transpose closure (BlockDot + reduced cold-flow dot test), BOTH adjoint labels
```
Reduced cold-flow operator transpose dot-test (thermalCoupling) max relative error=6.86412558997e-14
BlockDot (thermalCoupling): PU=2.07803300324e-13 PP=2.29255815574e-14 UU=7.3095183677e-14 UP=1.79609743792e-12 boundaryU=1.80391921425e-13
Reduced cold-flow operator transpose dot-test (pressureDrop) max relative error=6.86412558997e-14
BlockDot (pressureDrop):      PU=2.07803300324e-13 PP=2.29255815574e-14 UU=7.3095183677e-14 UP=1.79609743792e-12 boundaryU=1.80391921425e-13
```
- The **PP block dot = 2.29e-14** exercises the corrected P-P sign + boundary diagonal
  (J and J^T both carry the corrected operator; PP symmetric).

### P5 — explicit / matrix-free consistency (corrected explicitJT.mtx vs corrected matrix-free J^T)
```
ExplicitJToracle: maxRelL2=3.75731177462e-16 maxCosErr=2.22044604925e-16 maxRelL_U=3.53382916868e-16 maxRelL_P=3.75831856052e-16
```
- 5 random test vectors; ~3.8e-16 both blocks. The corrected explicit CSR == corrected
  matrix-free J^T to machine precision.

### P6 — StageB5 G1/G2 anchors (no regression; identical to BFINAL-003)
```
GateA eps=0.001 momentum:   relL2=1.65311363075e-4 cos=0.999999986345   (BFINAL-003 1.65311e-4)
GateB eps=0.001 candA:      relL2=2.98506496345e-05 cos=0.999999999556 (BFINAL-003 2.98506496e-05)
GateB eps=0.001 candB:      relL2=1.05921688465e-11 cos=1              (BFINAL-003 1.05921670e-11)
```
- eps sweep stable (candA ~2.985e-5 all eps; candB 1.06e-11..1.06e-9; momentum 1.6531e-4).

### P7 — StageB6 R_x dot tests + RxPressureRow transpose dot test (no regression vs BFINAL-005)
```
StageB6 weighted D1: |R_U,x d|L2=0.566318445188 |R_P,x d|L2=0.000193062746532
  D_momentum=-0.0486867554277 D_pressure=-0.00754887795106 D_total=-0.0562356333787
  prodGsenDPressDrop*d=-0.0486867553288 |D_mom-prod|/|prod|=2.03023945995e-09
StageB6 weighted D2: ... |D_mom-prod|/|prod|=5.39814148733e-10
StageB6 weighted D3: ... |D_mom-prod|/|prod|=7.65196471109e-10
RxPressureRowTranspose: sum(T*w)="-2.696485451749351e-08" sum(pc*(J_P*w))="-2.696485451749358e-08"
  relErr=2.57680498181e-15 |T|L2=0.000256656999292 |J_P*w|L2=6.22344977639e-10 w-support=5040
```
- BFINAL-005 references: |D_mom-prod|/|prod| = 2.4e-09 / 5.4e-10 / 1.0e-09 — identical.

### P1/P3 — stageB8 actual-residual FD (production residual semantics; probe ran in BOTH
adjoint blocks with byte-identical metrics — internal repeatability)
```
stageB8: fixesValue(p) patch = outlet nFaces=84; cells N=33600 internal=26208 otherBoundary=7308 outletAdjacent=84
stageB8: null vector ||w||=2.08447046292e+18 ||n_P||=0.999978325577 n_P[PREF=cell0]=0
```
Independent Python recomputation from artifacts (s3_gate_metrics.py), eps=1e-2:

| metric | relL2 | cos | normRatio |
|---|---|---|---|
| P1 full-P | 1.4429e-08 | 1.0 | 0.9999999992 |
| P1 outlet-adjacent | 5.4549e-09 | 1.0 | 0.9999999993 |
| P1 off-outlet | 5.2956e-08 | 1.0 | 0.9999999974 |
| P3 dir1 P-total | 1.4429e-08 | 1.0 | 1.0 |
| P3 dir2 U rows / P-total | 1.6537e-04 / 1.5158e-04 | 0.99999999 / 0.99999999 | — |
| P3 dir3 U rows / P-total | 1.6532e-04 / 1.5280e-04 | 0.99999999 / 0.99999999 | — |

- P1 eps plateau (full-P): 1e-2: 1.44e-08, 1e-3: 1.17e-07, 1e-4: 1.04e-06, 1e-5: 8.23e-06
  (exact linearity => best closure at largest eps; consistent with the verified FD).
- P1 sign/structure (eps=1e-2): sameSignFraction=0.7845, |FD| mass outlet=0.686248479375,
  |Jv| mass outlet=0.686248477495 (identical outlet concentration).
- P3 dir2/dir3 P-total ~1.5e-4 == known O(1e-4) derivative floor (BFINAL-003 momentum 1.653e-4).

### Corrected explicitJT.mtx (P5 + P2 raw material)
- Written to `/home/ys/dsH/b8_scratch_r1/explicitJT.mtx` (134400x134400); archived
  `artifacts/s3/explicitJT.mtx`.
- **Independent P-P diff vs OLD export** (`b2_case_smoke/explicitJT.mtx`, pre-patch):
  - P-P structure identical (227320 nnz both; 0 only-in-OLD / 0 only-in-NEW);
  - all 193,717 off-diagonal P-P entries: NEW == -OLD exactly (rel dev = 0.0);
  - off-outlet diagonals: NEW == -OLD to 2.4e-18;
  - 84 outlet-cell diagonals: NEW == -OLD + [+5.479951e-10, +1.706289e-09], mean +1.157420e-09
    == BFINAL-007 Q2 internalCoeffs_b magnitude (mean -1.1574e-09) exactly;
  - pRef row 100800: byte-identical identity pin (unchanged).
  → the corrected explicit operator carries exactly the patch: interior P-P sign flip
  (-L_int) + fixedValue-outlet boundary internalCoeffs diagonal, reference row untouched.

## 2. Run-2 (production path): MTO_RC=0, 0 crash markers, wall 3136 s

Scratch `/home/ys/dsH/b8_scratch_prod` (stageB4/stageB5/stageB6 false, stageB8 absent).
Production site = `solveDiscreteFlowAdjointProduction.H` (FGMRES + L1 row-sum +
block-lower-triangular preconditioners) — the FIRST production-path exercise on this
case/HEAD (no prior "Production reduced discrete flow adjoint" log exists anywhere).

- **MTO_RC=0**, 0 crash markers (`FOAM FATAL|sigFpe|Segmentation|abort|Floating point`
  grep count = 0), `ExecutionTime = 3136.4 s`.
- **Production-site dot test / state-derivative validation (always-on in the
  production run path):**
  - `Thermal matrix transpose dot-test max relative error=1.09689091013e-13` (same
    value as run-1 — deterministic);
  - `Discrete objective derivative errors: thermal=2.7025486001e-11, pressure=7.07793481308e-12`
    (g_w unit test `validateDiscreteObjectiveDerivatives.H`, abort threshold 1e-9 — PASS);
  - `E4` Helmholtz dot test does NOT run in this case: `gradientValidated true` in
    optProperties gates Stage-E off (pre-existing case configuration, unchanged).
- **Production adjoint solves with the corrected operator (both labels):**
  - thermalCoupling: `FGMRES iterations=4000, true relative residual=4.55432043137e-05`
    → `did not converge ... (DIAGNOSTIC MODE). trueRelRes=4.55432043137e-05, tolerance=1e-09. Continuing...`
  - pressureDrop: `FGMRES iterations=4000, true relative residual=0.0110569777489`
    → same DIAGNOSTIC MODE warning (discreteProdConvergeFatal=false in the case).
  - NOTE (context, NOT a gate and NOT a regression claim): production FGMRES
    convergence is a preconditioner-strength limitation, not operator validity —
    the project's documented position is that production Krylov convergence is NOT
    the primary blocker (exact offline SuperLU solves reach 3.48e-9, handoff §5.4)
    and solver/preconditioner tuning is explicitly NOT authorized by this task.
    P8 (exact direct tangent solve) is the structural-solvability gate.
- **Corrected-operator structural diagnostics at the production site:**
  - `PRODSOLV (thermalCoupling): sum(J^T*1)_P=12454.1904658 |J^T*1_P|L2=1429.48391438`
    (and identical for pressureDrop) — with the corrected P-P boundary diagonal the
    all-ones P vector is NO LONGER in the J^T null space (row sums no longer ~0),
    i.e. the outlet-localized non-constant null mode structure is gone;
  - `PRODDIAG ... jacobiNearZero=0` (no near-zero P Jacobi diagonal);
  - `PRODRESID (thermalCoupling): |rU|L2=2.48e-04 |rP|L2=1.94e-04` (partial adjoint).

## 3. Files (this step)

- `src/` — unchanged this step (S1 patch + S2 probe verified in place; incremental
  hashes in `cycle-1/pre-state-s2/` and `post-state/`).
- New: `cycle-1/s3_gate_metrics.py`, `cycle-1/s3_gate_metrics.json`,
  `cycle-1/s3_explicit_PP_diff.py`, `cycle-1/s3_archive.sh`, `cycle-1/S3_run_evidence.md`,
  `cycle-1/VALIDATION_SUMMARY_S3.json`
- Artifacts: `cycle-1/artifacts/` (stageB8_* probe exports, regenerated by run-1,
  sha256-identical to S2) + `cycle-1/artifacts/s3/` (both run logs, corrected
  explicitJT.mtx, explicitRhs_*, stageB5_*, stageB6_* exports); sha256 manifest
  `cycle-1/artifacts/sha256_artifacts.txt` (57 entries).

## 4. Gate status (evaluation is S4+; this step delivers the raw evidence)

| gate | evidence in this step | value |
|---|---|---|
| P1 | stageB8 P1 tables + JvP vs BFINAL-007 FD_lin | relL2 1.44e-8 (eps=1e-2); **2.946e-9 vs verified FD_lin** |
| P3 | stageB8 P3 dir1/2/3 tables | P-total 1.44e-8 / 1.52e-4 / 1.53e-4 @eps=1e-2 |
| P4 | BlockDot + reduced cold-flow dot | PP=2.29e-14, PU=2.08e-13, UU=7.31e-14, UP=1.80e-12; dot-test 6.86e-14 |
| P5 | ExplicitJToracle + explicit P-P diff | 3.76e-16; diff = sign flip + outlet boundary term, pRef row identical |
| P6 | StageB5 G1/G2 | candA 2.985e-5, candB 1.059e-11, momentum 1.653e-4 (== BFINAL-003) |
| P7 | StageB6 D1/D2/D3 + RxPressureRow dot | 2.03e-9/5.40e-10/7.65e-10; relErr 2.58e-15 (== BFINAL-005) |
| P8 | (later step: exact direct tangent solve) | — |
