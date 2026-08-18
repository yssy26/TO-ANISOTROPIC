# BFINAL-008 — Post-Reviewer independent verification & final decision

- Stage: B-final · Mode: **PATCH** · Hypothesis: `BFINAL-008-JPP-ACTUAL-PRIMAL-CLOSURE`
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace, verified) · HEAD `ca8a772b8339a36e7906ea32a7125061c47aa280` · branch `agent/dsH-stage-b-validation`
- Reviewer: independent Post-Reviewer invocation (this file written by the Post-Reviewer, not the Executor)
- Date: 2026-08-18
- **Final decision: PASS_JPP_PRESSURE_BC_PATCH**

---

## 0. What I personally verified (not from Executor summaries)

1. **Workspace gate** — `pwd`, `git rev-parse --show-toplevel`, `git rev-parse HEAD`,
   `git branch --show-current`, `git status --short` run by me: root == workspace,
   HEAD == ca8a772b, branch == agent/dsH-stage-b-validation. Dirty set = the 4 expected
   tracked files (MTO_HF.C / sensitivity.H pre-existing, solveDiscreteFlowAdjoint*.H =
   BFINAL-003 heads + this round's patch) + untracked locked probes/evidence; `build/`
   untracked as expected.
2. **Scope / diff audit** — I read the S1 incremental diffs
   (`cycle-1/incremental_diff_solveDiscreteFlowAdjoint.H.diff`,
   `...Production.H.diff`) which isolate the round's changes vs a reconstructed
   pre-state worktree: ONLY the two authorized files touched. Re-verified against the
   live worktree: `sensitivity.H` = `f0ed0726…`, `rxPressureRowTranspose.H` = `fd648ae0…`,
   `MTO_HF.C` = `1f190c0c…`, `stageB7PressureBCDiagnostic.H` = `d7f7ecf0…` (all locked
   hashes recomputed by me, matching BFINAL-003/005/007 records). Pre-state
   `solveDiscreteFlowAdjoint.H` = `f0c81497…` / `solveDiscreteFlowAdjointProduction.H` =
   `8c4901bf…` == BFINAL-007 locked heads. Post-state: `8390cfa4…` / `fdd5ea2f…` (current
   live hashes == FINAL_REPORT records). `git diff --check` clean. Real case
   `b2_case_smoke/constant/optProperties` mtime `2026-08-17 00:22:38`, sha256
   `6d7e4de1…`, `stageB8` count = 0 (untouched).
3. **No outlet-cell identity hack** — I read the actual patched loops (forward
   L~760-775, transpose L~920-945, CSR export L~2660-2685, production L~475-490): the
   fixedValue-p boundary term is a boundary-FACE flux linearization
   `+primalPressureMobility[celli]*deltaCoeffs_b*magSf_b*input[P(celli)]` gated ONLY on
   `p.boundaryField()[patchi].fixesValue()`; no identity row is inserted anywhere.
   Independent offline scan of the corrected `explicitJT.mtx` P block: unique identity
   row = [100800] (= discretePIndex(0), the existing reference pin), row 106400 =
   physical (34 nnz); my own scan confirmed `row100800.nnz == 1` with single 1.0.
4. **Corrected PP interior sign** — I read the forward interior loop: the P-P term in
   `deltaPhiFacePU` is now `- mobF*(input[P(nei)]-input[P(own)])*dcfF*mafF` (exported
   +L_int → −L_int); the relaxation U-part `Sf&(alphaRel*rAUdH+(1-alphaRel)*interpDU)`
   and the OLD `deltaPhiFace` (velocityJump J_UU/J_UP) are byte-identical. Transpose kf
   path: `phiAdjointPP = lambdaPdiff - convTerm` with `P(own) += kf*phiAdjointPP;
   P(nei) -= kf*phiAdjointPP` (only the lambdaPdiff coefficient flips; convTerm/
   phiAdjointRel/hA/directRel unchanged). CSR export: four P-P addCoo flipped
   (diag +kf, off-diag −kf) with convTerm velocityJump addCoo (both own and nei rows)
   byte-identical. Production `applyProdFlowJT` same flip; `prodPDiagPhysical` −kf→+kf
   plus boundary term; production pins byte-identical.
5. **Forbidden hardcodes** — grep of added lines: no 84 / 100800 / 106400 / 5600 /
   7.4e-27 / 7.9e-12 in production logic; the only 0.4/0.6 hits are comments and
   pre-existing diagnostic code (sin-seeded test vectors, Info strings), not the patch.
6. **Probe FD semantics** — I read `stageB8JPPActualResidualFD.H` in full: it rebuilds
   the ACTUAL frozen-primal residual `R(w)=[R_U;R_P]` from production semantics
   (`resU = -UEqn.residual()` NS.H sign; `phiHbyA` rebuilt from PRODUCTION p —
   p-independent; `pEqn = fvm::laplacian(rAtU(),pmod)` on the perturbed psi so
   `pEqn.flux()` carries the P-P block incl. fixedValue outlet internalCoeffs;
   `resP = div(phiHbyA - pEqn.flux())` raw face sums, no 1/V). This is the BFINAL-007
   Q1-verified actual-residual semantics, NOT the old diagnostic residual definition.

---

## 1. Independent reproduction performed (this round)

I executed my own fresh runs/computations (no reuse of Executor scripts):

1. **Fresh probe binary run** (production binary, corrected J + probe) on a new scratch
   `b8_review_r1` (byte-identical copy of the validated scratch with only the explicit
   matrix path re-pointed): MTO_RC=0, ExecutionTime 508.7 s. **Determinism**: the fresh
   run log differs from the archived run-1 log in EXACTLY 3 lines (Exec/Case/CSV path
   strings); every numerical line is byte-identical. All 22 stageB8 artifacts
   byte-identical (sha256) to the archive; fresh `explicitJT.mtx` sha256
   `5ca1f592…` == archived corrected export.
2. **Fresh offline P1/P3** from raw artifacts with my own code
   (`reviewer_independent.py`, fresh implementation): P1 full-P relL2 = **1.4429e-08**
   (cos 1.0, normRatio 0.99999999915), outlet-adjacent 5.4549e-09, off-outlet 5.2956e-08;
   P3 dir1/2/3 P-total = 1.443e-08 / 1.516e-04 / 1.528e-04 — all matching the
   Executor's fresh validate JSON bit-for-bit (e.g. Executor full-P 1.4428509283701247e-08).
3. **Fresh P2 (splu sigma_min)** — my own splu of the corrected re-pinned J (row
   100800 → identity): sigma_min trials 8.074e-12 / 7.222e-12 / 6.802e-12 → **min
   6.80e-12, mean 7.37e-12** (Executor 8.62e-12–1.22e-11; BFINAL-007 candidate 4.6e-12–
   1.2e-11; random-vector spread expected). `||J_corrected·n||/||n||` = **3.019933e-09**
   (Executor 3.01993291393e-09 — identical). Null-mode geometry: old n_P top-1% mass
   0.895 / outlet 0.352; corrected v_min top-1% mass 0.029 / outlet mass ~6.5e-07,
   max |v_P| at cell 32495 (NOT an outlet cell) — no outlet-concentrated near-null.
4. **Fresh P5** — my own explicit P-P block vs −L_prod rebuilt from BFINAL-007
   production laplacian exports: PP diag relL2 = **2.424e-10** (known cell-vs-boundary
   mobility interpolation residue); outlet boundary entry vs −internalCoeffs_b relL2 =
   **5.933e-09**; P-P block exactly symmetric (max|PP−PPᵀ| excl cell0 = 0.0); pRef row
   nnz=1 single 1.0.
5. **Fresh P8** — my own splu solve `J wPrime = -R_x·d` D1: true relative residual =
   **5.330576e-14** (≤1e-9), Urel 4.63e-14, Prel 7.75e-11, ||wPrime|| 2.89e+07
   (physical scale; pre-patch floor D1 3.647e-4 with ~1e18–19 null-dominated wPrime).
6. **P7 direct case-field checks** — `1/designMask` on the 84 outlet cells: min=max=
   sum=0 (verified from the case field myself); `stageB6_deltaAlpha.mtx` on outlet
   cells: nnz=0. Executor's independent D1/D2/D3 contraction vs locked BFINAL-005:
   relErr ≤ 1.4e-12; `RxPressureRowTranspose` dot relErr 2.58e-15 (in my fresh run log).

---

## 2. Gate verdicts (my independent evaluation)

| gate | acceptance | my independent value | verdict |
|---|---|---|---|
| P1 | relL2 ≤ 1e-6 (pref 1e-8); O(1)=FAIL | full-P 1.44e-8, outlet 5.45e-9, off-outlet 5.30e-8 (fresh run + fresh offline) | **PASS** |
| P2 | sigma_min ~1e-11 (not 7.4e-27); no outlet near-null; gauge intact | 6.80e-12–8.07e-12 (my splu); ||Jn||/||n|| 3.019933e-09; v_min outlet mass ~0; identity row [100800] only | **PASS** |
| P3 | P-total ≤ 1e-3 (pref 1e-4 / O(1e-4) floor) | 1.44e-8 / 1.52e-4 / 1.53e-4 | **PASS** |
| P4 | full dot ≤ 1e-10 (pref machine) | PP 2.29e-14, all blocks ≤ 1.8e-12, reduced cold-flow 6.86e-14 (fresh run byte-identical) | **PASS** |
| P5 | explicit == matrix-free ~1e-16 | ExplicitJToracle 3.76e-16; explicitJ vs matrix-free Jv 4.0e-16/3.9e-16 (clean dirs); PP vs −L_prod 2.42e-10; bnd 5.93e-9; sym 0.0 | **PASS** |
| P6 | no J_PU/momentum regression | GateA momentum 1.65311363075e-4, candA 2.98506496345e-5, candB 1.05921688465e-11, T1 dp=0 0.00852323456072 — all == BFINAL-003 (fresh run) | **PASS** |
| P7 | R_P,x negligible; no R_X_REOPEN | designMask outlet 0, deltaAlpha outlet 0, boundary L2 5.36e-13, D1/D2/D3 relErr ≤ 1.4e-12 | **PASS** |
| P8 | true relRes ≤ 1e-9 (≥1 strong dir) | D1 5.33e-14 (mine); D2/D3 4.29e-14/9.17e-14 (Executor, same splu class) | **PASS** |

---

## 3. Scope / stop-rule compliance

- No change to NS.H / sensitivity.H / rxPressureRowTranspose.H / J_PU relaxed-SIMPLE /
  MMA / filter / projection / objective / constraint / outlet BC / relaxation / scales /
  tolerances — all forbidden-file hashes unchanged (recomputed by me).
- No empirical diagonal regularization, no epsilon·I, no pinning of outlet cells, no
  pressure mean removal, no pRefCell change.
- The corrected transpose is the exact adjoint of the corrected forward (derived, not
  tuned); explicit CSR export mirrors the same corrected operator.
- **No tangent oracle entered**: the run contains only the pre-existing Stage-B4
  momentum/continuity diagnostics (T1/T2 labels are BFINAL-003-era anchors required by
  P6, not the BFINAL-006 tangent oracle); P8 is a T1-style structural smoke only.
  T2/T3 NOT executed.
- The three mandatory statements are present verbatim in the Executor's FINAL_REPORT:
  "The flow Jacobian now represents the actual frozen-primal pressure residual
  including the fixedValue outlet boundary contribution and the correct P-P residual
  sign convention." / "J_PU remains closed from BFINAL-003." / "R_x remains complete
  from BFINAL-005 for the current design region."

---

## 4. Minor notes (non-blocking, for the record)

1. P3 dir3 P-outlet-adjacent relL2 = 1.03e-3 slightly exceeds the 1e-3 preferred bound,
   but it is the tiny-FD-noise block (|FD| ~7.3e-8) at the same intrinsic O(1e-4)
   derivative floor; P-total for dir3 = 1.53e-4 passes. Documented and reproduced.
2. P1 sameSignFraction = 0.7845 (not ~1): the direction has significant near-zero
   entries where sign is roundoff-sensitive; the decisive metrics are relL2 (1.44e-8),
   cos=1.0, and the identical outlet mass concentration (0.68625 both Jv and FD).
3. P3 dir1 U-rows relL2 = 0.123: the P-only null direction has a tiny U FD signal
   (|FD| 3.28e-7) dominated by the roundoff-level UEqn response — a U-block statement,
   not a P-row statement; P-block values are the gate metrics.
4. sigma_min spread across trials (6.8e-12–8.1e-12 mine; 8.6e-12–1.2e-11 Executor) is
   the expected random-vector variability of the splu-based estimate; all trials are
   ≥1e-13 and consistent with the BFINAL-007 candidate ~7.9e-12 (×1e15 above the
   7.4e-27 pre-patch floor).

---

## 5. Decision

**PASS_JPP_PRESSURE_BC_PATCH.** Every gate P1–P8 was independently reproduced by the
Post-Reviewer from the raw artifacts, the fresh binary run, and fresh offline
computations (my own splu / explicit-matrix / FD arithmetic). Scope is exactly the two
authorized files; no outlet identity hack; corrected P-P sign and fixedValue boundary
internalCoeffs present in forward / transpose / explicit CSR / production; no tangent
oracle entered. Per task §22: STOP after Post-Review; the next user-authorized
operation is a corrected rerun of the tangent oracle.

Raw commands (my independent path):
```bash
# fresh probe binary run (deterministic; log differs from archive only in paths)
source /opt/openfoam7/etc/bashrc && unset FOAM_SIGFPE && \
FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin \
/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF -case /home/ys/dsH/b8_review_r1 \
  > /home/ys/dsH/b8_review_r1/Log.stageB8.review.txt 2>&1   # MTO_RC=0
# fresh offline P1/P2/P3/P5/P7/P8 (my own implementation)
PY=/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-002/.venv_s4/bin/python3
$PY evidence/agent-group/BFINAL-008/reviewer_independent.py   # -> reviewer_independent.json/.log
```
Evidence files: `evidence/agent-group/BFINAL-008/reviewer_independent.py`,
`reviewer_independent.log`, `reviewer_independent.json`; fresh run logs under
`/home/ys/dsH/b8_review_r1/`.
