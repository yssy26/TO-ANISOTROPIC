# BFINAL-008 FINAL_REPORT — Production J_PP := dR_P/dp actual-primal closure (−L_int interior sign + fixedValue outlet boundary internalCoeffs) across forward J / matrix-free J^T / explicit CSR / production transpose

- Stage: B-final · Mode: **PATCH** · Hypothesis: `BFINAL-008-JPP-ACTUAL-PRIMAL-CLOSURE`
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace, verified) · HEAD `ca8a772b8339a36e7906ea32a7125061c47aa280` · branch `agent/dsH-stage-b-validation`
- Case: `/home/ys/dsH/b2_case_smoke` (real case untouched; runs on scratch `b8_scratch_r1` / `b8_scratch_prod`)
- Evidence: `evidence/agent-group/BFINAL-008/cycle-1/` (pre-state, post-state, incremental diffs, S3 run evidence, S4 offline validation, validate-step independent re-run)
- Date: 2026-08-18 (Executor, BFINAL-008 S4 + validate step)
- **Final decision: PASS_JPP_PRESSURE_BC_PATCH**

> Validate step (Executor, independent re-run of every approved validation entry on
> fresh scratch cases `b8_scratch_v1`/`b8_scratch_v2` and fresh offline
> computations): **all entries reproduce the archived S2/S3/S4 numbers
> bit-for-bit** — run-1/run-2 logs differ from the archived logs ONLY in Time/PID
> lines; all 22 stageB8 artifacts byte-identical; fresh `explicitJT.mtx` sha256
> `5ca1f592…` == archived corrected export; fresh P2/P8 json identical
> (sigma_min 8.62e-12–1.22e-11, ||Jn||/||n|| 3.019933e-09, P8 D1/D2/D3
> trueRelRes 5.33e-14/4.29e-14/9.17e-14). Consolidated per-entry PASS table in
> `cycle-1/VALIDATION_SUMMARY.json` (validate step).

> The flow Jacobian now represents the actual frozen-primal pressure residual
> including the fixedValue outlet boundary contribution and the correct P-P
> residual sign convention.
>
> J_PU remains closed from BFINAL-003.
>
> R_x remains complete from BFINAL-005 for the current design region.

---

## Provenance

| item | value |
|---|---|
| pwd / git root | `/home/ys/dsH/TO-ANISOTROPIC` (== workspace) |
| HEAD | `ca8a772b8339a36e7906ea32a7125061c47aa280` (unchanged through the round) |
| branch | `agent/dsH-stage-b-validation` |
| tracked dirty files | `src/MTO_HF.C` (pre-existing), `src/sensitivity.H` (BFINAL-005 locked), `src/solveDiscreteFlowAdjoint.H`, `src/solveDiscreteFlowAdjointProduction.H` (BFINAL-003 locked heads **+ this round's authorized J_PP patch**) |
| untracked (preserved) | `src/rxPressureRowTranspose.H`, `src/stageB5BoundaryRelaxOracle.H`, `src/stageB6RxDesignOracle.H`, `src/stageB7PressureBCDiagnostic.H`, `src/stageB8JPPActualResidualFD.H`, `build/`, `evidence/…` |
| files changed this round | S1 (PATCH): ONLY `src/solveDiscreteFlowAdjoint.H` + `src/solveDiscreteFlowAdjointProduction.H` (see `cycle-1/incremental_diff_*.diff`). S2 (probe wiring): `src/solveDiscreteFlowAdjoint.H` +25-line switch-guarded include of `stageB8JPPActualResidualFD.H` (default false, static once-guard, never modifies production state — same proven pattern as stageB5/B6/B7) + new untracked `src/stageB8JPPActualResidualFD.H`. `sensitivity.H`/`MTO_HF.C` byte-identical to pre-state (S1-verified). Evidence under `cycle-1/`. |
| post-state sha256 | S1 post: `solveDiscreteFlowAdjoint.H` `f41d1753…`, `solveDiscreteFlowAdjointProduction.H` `fdd5ea2f…`. S2 final (current): `solveDiscreteFlowAdjoint.H` `8390cfa4…` (= S1 post + guarded probe include), `solveDiscreteFlowAdjointProduction.H` `fdd5ea2f…` (unchanged — no stageB8 reference). Records in `cycle-1/post-state/` and `cycle-1/pre-state-s2/`. |
| binary | `build/bin/MTO_HF` (S1 patch + S2 probe build, mtime 2026-08-18 11:27, newer than all src files) |
| real case | `/home/ys/dsH/b2_case_smoke` untouched: `constant/optProperties` mtime `2026-08-17 00:22:38`, sha256 `6d7e4de1…`, no `stageB8JPPActualResidualFD` switch |
| runs | run-1 probe `b8_scratch_r1` (stageB8 true, MTO_RC=0, 496.7 s) · run-2 production `b8_scratch_prod` (stageB4/5/6 false, stageB8 absent, MTO_RC=0, 3136.4 s) |

Forbidden modules (`NS.H`, `sensitivity.H`, `rxPressureRowTranspose.H`, MMA, filter/projection,
objective/constraint, case `{0,constant,system}`) have zero diff this round (S1 verified;
`git diff --check` clean). No tangent oracle entered; T2/T3 not run.

---

## Defect repaired

Old J_PP convention (BFINAL-007, independently reproduced by Executor + Post-Reviewer):
the coupled J (forward, transpose, exported `explicitJT.mtx` alike) assembled the P-P block
ONLY from interior faces with the sign `+L_int` and omitted the fixedValue outlet boundary
contribution of the production `fvm::laplacian(rAtU,p)`; the exported interior P-P block was
the NEGATIVE of the interior residual derivative (`dR_P/dp|_int = −L_int`).

Actual-primal convention (this patch): `J_PP_actual = dR_P/dp = −L_prod = −(L_int + bnd)`.

1. **Interior sign** — forward interior loop (solveDiscreteFlowAdjoint.H): the P-P term in
   `deltaPhiFacePU` flipped `+ mobF*(p_nei−p_own)*dcfF*mafF` → `− mobF*(…)` (exported
   `+L_int` → `−L_int`). The relaxation U-part `Sf&(alphaRel*rAUdH + (1−alphaRel)*interpDU)`
   and the OLD `deltaPhiFace` (velocityJump J_UU/J_UP) are byte-identical (verified in the
   incremental diff).
2. **fixedValue-p boundary contribution** — for every patch with `p.boundaryField()[patchi].fixesValue()`
   (this case: the outlet only), the cell P-P diagonal gains
   `+ primalPressureMobility[celli]*deltaCoeffs_b*magSf_b` (= `−internalCoeffs_b` of the
   production laplacian, residual convention; BFINAL-007 Q2). Derived from actual mesh boundary
   fields, gated only on `fixesValue()` (zeroGradient p patches contribute 0); **no source term**
   (boundaryCoeffs = 0 since p_b = 0) and **no outlet identity rows**; the only reference/pin row
   remains the existing system reference treatment (row 100800 = P(0), unchanged).
3. **Transpose (exact adjoint of the corrected forward, not tuned)** — the kf path becomes
   `phiAdjointPP = lambdaPdiff − convTerm; P(own) += kf*phiAdjointPP; P(nei) −= kf*phiAdjointPP`
   (only the `lambdaPdiff` coefficient flips; `convTerm`/`phiAdjointRel`/`hA`/`directRel`
   (J_PU) byte-identical); the boundary P-P diagonal is symmetric and carried identically.
   Same in the production transpose `applyProdFlowJT`.
4. **Explicit CSR export** — the four P-P `addCoo` flip (`diag +kf`, `off-diag −kf`), the
   convTerm velocityJump `addCoo` unchanged, plus boundary `addCoo(P(celli),P(celli),+term)` on
   `fixesValue()` patches; pRef pin row byte-identical. Production `prodPDiagPhysical`
   `−kf → +kf` plus the same boundary term.

No hard-coded 84 / 100800 / 106400 / 5600 / 0.4 / 0.6 / 7.4e-27 / 7.9e-12 in production logic
(S1 grep-verified); all values derived from the actual case/operator.

Boundary term magnitudes (this case; from actual fields): `δ_b = 4000`, `|Sf|_b = 2.5e-7`,
`rAtU_b ∈ [5.48e-7, 1.71e-6]` → `internalCoeffs_b ∈ [−1.706e-9, −5.480e-10]`, residual-convention
`+rAtU_c·δ_b·|Sf|_b ∈ [5.48e-10, 1.71e-9]` on the 84 outlet cells.

---

## Nullspace table (Gate P2)

Re-pinned corrected J = transpose of corrected `explicitJT.mtx`, row 100800 → identity
(BFINAL-006 S1 gate convention); n = normalized `wPrime_TAN_D1` (BFINAL-006 null vector,
`||w|| = 2.084470e18`, `n[100800] = 0`). splu (COLAMD) + fresh random solves (3 seeds).

| metric | exported J (pre-patch) | BFINAL-007 candidate | **corrected production J (this round)** |
|---|---|---|---|
| sigma_min | 7.382e-27 (BFINAL-006) | 7.9e-12 (Post-Reviewer mean) | **8.62e-12 min / 1.088e-11 mean** (3 fresh random solves: 1.187e-11, 8.620e-12, 1.216e-11) |
| \|\|J·n\|\|/\|\|n\|\| | 2.717e-19 | 3.020e-09 | **3.019933e-09** |
| unique identity row in P block | [100800] | [100800] | **[100800]** (row 106400 physical, 34 nnz) |
| null-vector localization (old n) | outlet-concentrated (top-1% 0.895, outlet-cell 0.352) | — | old n under corrected J: **no longer a null vector** |
| corrected-J near-null direction geometry | — | — | v_min P-block top-1% mass = **0.031** (old 0.895), outlet-cell mass = **6.5e-07** (~0), max \|v_P\| at cell 32480 (NOT an outlet cell) |

**P2 verdict: PASS** — sigma_min raised 7.4e-27 → 1.09e-11 (×1.5e15, consistent with the
BFINAL-007 candidate ~7.9e-12), no outlet-localized near-null vector comparable to the old mode
(v_min outlet mass ≈ 0 vs old 0.352; top-1% 0.031 vs 0.895), reference/gauge semantics intact
(unique identity row = 100800; row 106400 physical). The all-ones P vector is also no longer in
the J^T null space at the production site (PRODSOLV run-2: `sum(J^T·1)_P = 12454.19`,
`|J^T·1_P|L2 = 1429.48`; PRODDIAG `jacobiNearZero = 0`).

---

## Actual residual FD table (Gates P1 / P3)

FD semantics = ACTUAL frozen-primal pressure residual (production semantics:
`R_P = div(phiHbyA) − div(pEqn.flux())`, `pEqn = fvm::laplacian(rAtU,p)` incl. fixedValue outlet
boundary coefficients; probe `stageB8JPPActualResidualFD.H`; NOT the old diagnostic residual
definition). Jv = corrected matrix-free forward `applyDiscreteFlowJ`. Central FD `(R(w+eps v) −
R(w−eps v))/(2eps)`, `eps = 1e-2` authoritative (exact linearity: best closure at largest eps).

### P1 — targeted null-vector actual-residual closure (dir1 = n_P, dp-only)

| block | relL2 | cos | normRatio | \|Jv\| / \|FD\| |
|---|---|---|---|---|
| full P | **1.4429e-08** | 1.0 | 0.9999999992 | 1.76831225e-9 / 1.76831225e-9 |
| outlet-adjacent (84) | 5.4549e-09 | 1.0 | 0.9999999993 | 1.71051053e-9 / 1.71051053e-9 |
| off-outlet | 5.2956e-08 | 1.0 | 0.9999999974 | 4.48421620e-10 / 4.48421621e-10 |

sameSignFraction = 0.7845; |FD| mass in the 84 outlet cells = 0.68625, |Jv| mass = 0.68625
(identical outlet concentration — same sign/structure). eps plateau (full-P relL2): 1e-2: 1.44e-8,
1e-3: 1.17e-7, 1e-4: 1.04e-6, 1e-5: 8.23e-6 (linear, roundoff-dominated at small eps).
**Acceptance relL2 ≤ 1e-6 (pref. 1e-8): PASS.**

### P3 — general actual-residual B-F2 (deterministic directions, eps=1e-2)

| dir | U rows | P internal | P outlet-adj | P other-bnd | **P total** |
|---|---|---|---|---|---|
| dir1 (n_P) | 1.231e-1* | 6.343e-8 | 5.455e-9 | 3.276e-8 | **1.443e-8** |
| dir2 (dp-heavy outlet-supported) | 1.6537e-4 | 1.2816e-4 | 1.3574e-4 | 1.6670e-4 | **1.5158e-4** |
| dir3 (seeded dU+dp) | 1.6532e-4 | 1.2815e-4 | 1.0310e-3 | 1.6671e-4 | **1.5280e-4** |

all cos ≥ 0.9999994; normRatio ≈ 1. *dir1 U rows: the P-only null direction has a tiny U FD
signal (|FD| 3.28e-7) dominated by the roundoff-level UEqn response — matches the S3/C++ value
exactly and is not a P-row statement. P-total ≤ 1e-3 hard, ~1e-4 preferred or consistent with the
known O(1e-4) derivative floor (momentum 1.653e-4): **PASS** for dir1/dir2; dir3 P-total 1.53e-4
(PASS); dir3 P-outlet-adjacent 1.03e-3 is the tiny-FD-noise block (|FD| 7.3e-8) at the same
intrinsic floor, identical to the S3 recompute — not a defect.

---

## Operator consistency (Gates P4 / P5)

### P4 — J/J^T transpose closure (run log, both adjoint labels; deterministic)

| block | value (thermalCoupling == pressureDrop) | threshold |
|---|---|---|
| PP | **2.29255815574e-14** | ≤1e-10 |
| PU | 2.07803300324e-13 | ≤1e-10 |
| UU | 7.3095183677e-14 | ≤1e-10 |
| UP | 1.79609743792e-12 | ≤1e-10 |
| boundaryU | 1.80391921425e-13 | ≤1e-10 |
| reduced cold-flow dot (both labels) | 6.86412558997e-14 | ≤1e-10 |
| Thermal matrix transpose dot-test | 1.09689091013e-13 | ≤1e-10 |

The PP block (2.29e-14) exercises the corrected P-P sign + boundary diagonal in both J and J^T.
**PASS.** Independent structural cross-check (explicit matrices): corrected P-P block is exactly
symmetric (max|PP − PPᵀ| excl. cell 0 = **0.0**).

### P5 — explicit / matrix-free consistency

- In-log ExplicitJToracle (explicit J^T vs matrix-free J^T, 5 random vectors): maxRelL2 =
  **3.757e-16**, maxRelL_U = 3.534e-16, maxRelL_P = 3.758e-16 (BFINAL-003 reference ~3.3e-16). **PASS.**
- Independent forward check (corrected exportᵀ vs matrix-free Jv exports, clean directions with
  v[pRef]=0): dir1 relL2 = **4.00e-16** (max|diff| 3.6e-23), dir2 relL2 = **3.93e-16** (max|diff| 7.6e-19).
  dir3 has v[PREF] = 0.172 ≠ 0: differs on exactly the 7 pRef-coupled rows (0,1,2,3,241,3362,100800)
  by `J_phys[i,PREF]·v[PREF]` (ratios ≈ 1.2e-7 physical couplings; row PREF ratio = 1.0) — the
  documented BFINAL-003 G4 pin-column convention of the export, NOT a patch defect.
- Corrected explicit P-P block vs `−L_prod` rebuilt from the BFINAL-007 production laplacian
  exports (same design point; primal untouched): diag relL2 = **2.42e-10** (known
  cell-vs-boundary mobility interpolation residue, same as the BFINAL-007 candidate);
  outlet-cell boundary entry `B[c,c] − (−d[c])` vs `−internalCoeffs_b` relL2 = **5.93e-9**
  (same residue class); pRef row = identity pin, byte-identical (14 stored, single 1.0).
- S3 explicit P-P diff vs OLD export: all 193,717 off-diagonal P-P entries NEW == −OLD exactly
  (rel dev 0.0); off-outlet diagonals NEW == −OLD to 2.4e-18; 84 outlet diagonals
  NEW == −OLD + [5.48e-10, 1.71e-9] (the boundary term); pRef row byte-identical.

**PASS** (machine precision; corrected P-P sign + outlet boundary internalCoeffs + reference row
all present; no old explicit matrix reused).

---

## No-regression (Gates P6 / P7)

### P6 — J_PU / momentum anchors (== BFINAL-003 byte-for-byte)

| anchor | this round | BFINAL-003 |
|---|---|---|
| GateA momentum relL2 (eps=1e-3) | **1.65311363075e-4** (cos 0.999999986345) | 1.65311e-4 |
| GateB candA relL2 (eps=1e-3) | **2.98506496345e-5** (cos 0.999999999556) | 2.98506496e-5 |
| GateB candB relL2 (eps=1e-3) | **1.05921688465e-11** (cos 1) | 1.05921670e-11 |
| eps sweep | momentum 1.6531e-4 stable; candA 2.985e-5 stable; candB 1.06e-11 → 1.06e-9 | identical |
| P<-U tangent anchor | P-U / P-P unaffected by design (dp=0 directions byte-identical; candidate J == exported J on dp=0) | BFINAL-003 |

### P7 — R_x re-audit (unchanged artifacts; R_x NOT modified by this patch)

- **designMask on outlet cells: max = 0.0** (case `1/designMask`).
- **deltaAlpha support on outlet cells: 0 cells** (`stageB6_deltaAlpha.mtx` nnz on the 84 outlet cells = 0).
- **pressure-row boundary R_x contribution**: `d(phi_b)/d(alpha) = +drAtU_b·δ_b·|Sf|_b·p_c` per
  outlet cell ∈ [−1.4388e-13, +5.3862e-14], **L2 = 5.358e-13** (BFINAL-007: 5.36e-13) =
  8.6e-4 of |J_P·w| (6.223e-10); design-weighted (deltaAlpha/designMask zero on outlet) it is
  **identically 0** — the BFINAL-005 `dflux_b = 0` assumption is exact for this case.
- **BFINAL-005 D1/D2/D3 contraction** (independent recompute, `D = −λᵀ R_x d`):
  D_mom −0.0486867554 / +0.2290796624 / −0.0301314016 (relErr vs locked ≤ 9e-13);
  D_pre −0.00754887795 / +0.06859276939 / −0.00621002822 (≤ 4e-13);
  D_tot −0.0562356334 / +0.2976724318 / −0.03634142984 (≤ 1.4e-12) — **locked anchors reproduced exactly**.
- Production contractions: momentum field·d = −0.0486867553 / +0.2290796625 / −0.0301314016
  (relErr 2.03e-9 / 5.40e-10 / 7.65e-10 == BFINAL-005 G2/G3 scale); volume field·d =
  0.1555186511 / −0.0844966774 / 0.0639623768 (relErr ≤ 3.2e-15 == BFINAL-005 G5 anchors).
- RxPressureRowTranspose dot test (run log): relErr = **2.57680498181e-15** (== BFINAL-005);
  StageB6 weighted D1/D2/D3 |D_mom−prod|/|prod| = 2.03e-9 / 5.40e-10 / 7.65e-10 (identical).

**PASS — no R_X_REOPEN required.** The corrected pressure-BC representation creates NO
non-negligible R_P,x term at this design point.

---

## Direct tangent solvability smoke (Gate P8)

`J_corrected · wPrime = −R_x·d` (re-pinned row 100800 → identity; rhs[100800] = 0),
splu (COLAMD) + iterative refinement (30 iters). This is a T1-style structural smoke ONLY
(not the full BFINAL-006 tangent oracle; T2/T3 not run).

| dir | \|\|R_x·d\|\| | true relative residual (target ≤ 1e-9) | U rel | P rel | \|\|wPrime\|\| | w[PREF] |
|---|---|---|---|---|---|---|
| D1 (strongest) | 5.663185e-01 | **5.330576e-14** | 4.631e-14 | 7.745e-11 | 2.891e+07 | 0 |
| D2 | 8.615839e-01 | **4.291321e-14** | 3.681e-14 | 4.242e-11 | 3.707e+07 | 0 |
| D3 | 8.516528e-01 | **9.166669e-14** | 7.964e-14 | 1.489e-10 | 9.004e+07 | 0 |

(BFINAL-006 pre-patch floor for comparison: D1 3.647e-4, D2 4.044e-5, D3 3.093e-3.)

**P8 verdict: PASS** — the structural singularity that blocked T1 has disappeared: all three
directions solve to true relative residual ~1e-13–1e-14 (5+ orders below the 1e-9 acceptance),
and \|\|wPrime\|\| dropped from ~1e18–1e19 (null-component-dominated) to ~1e7–1e8 (physical). No
new structural nullspace appeared. T2/T3 NOT executed (out of scope for this task).

---

## Reviewer independent reproduction

Post-Reviewer must independently re-verify (task §18):
1. production diff scope / no outlet identity hack — `cycle-1/incremental_diff_*.diff`, `post-state/`;
2. corrected PP interior sign + fixedValue boundary internalCoeffs — `s4_p1_p3_p4_p5_p7_fast.py` P5 + `s3_explicit_PP_diff.py`;
3. null-vector actual-residual FD — `s4_p1_p3_p4_p5_p7_fast.py` P1 (from `stageB8_Jv1/FD1_eps_*`);
4. smallest-singular/nullspace — `s4_p2_p8_nullspace_tangent.py` + `s4_p2_p8.json`;
5. PP/full dot test — run logs `Log.stageB8.r1.txt` BlockDot / reduced cold-flow dot;
6. explicit vs matrix-free — `ExplicitJToracle` log lines + `s4_*_fast` P5;
7. J_PU no-regression anchor — `StageB5 GateA/B` log lines (momentum 1.653e-4);
8. R_x no-regression contraction — `s4_p1_p3_p4_p5_p7_fast.py` P7 (D1/D2/D3 vs locked);
9. direct tangent solvability smoke — `s4_p2_p8.json`.

Raw commands:
```bash
cd /home/ys/dsH/TO-ANISOTROPIC
PY=/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-002/.venv_s4/bin/python3
$PY evidence/agent-group/BFINAL-008/cycle-1/s4_p1_p3_p4_p5_p7_fast.py   # P1/P3/P4/P5/P7 (fast)
$PY evidence/agent-group/BFINAL-008/cycle-1/s4_p2_p8_nullspace_tangent.py  # P2/P8 (heavy, ~35 min)
```

All load-bearing numbers were recomputed by fresh S4 code paths from the raw artifacts and run
logs (see `s4_*.json`, `S3_run_evidence.md`, `s3_gate_metrics.json`).

---

## Gate verdict summary

| gate | acceptance | result | verdict |
|---|---|---|---|
| P1 | relL2 ≤ 1e-6 (pref 1e-8), O(1)=FAIL | 1.44e-8 (full-P), 5.45e-9 (outlet), 5.30e-8 (off-outlet) | **PASS** |
| P2 | sigma_min ≥ 1e-13 ≈ 7.9e-12; no outlet near-null; gauge intact | sigma_min 8.62e-12–1.22e-11; v_min outlet mass 6.5e-7; identity row [100800] | **PASS** |
| P3 | P-total ≤ 1e-3 (pref 1e-4 / O(1e-4) floor) | 1.44e-8 / 1.52e-4 / 1.53e-4 | **PASS** |
| P4 | full dot ≤ 1e-10 (pref machine) | PP 2.29e-14, all blocks ≤ 1.8e-12, dot 6.86e-14 | **PASS** |
| P5 | explicit == matrix-free ~1e-16 | 3.76e-16 (J^T); 4.0e-16/3.9e-16 (J, clean dirs) | **PASS** |
| P6 | no J_PU/momentum regression | momentum 1.653e-4, candA 2.985e-5, candB 1.06e-11 (== BFINAL-003) | **PASS** |
| P7 | R_P,x negligible; designMask/deltaAlpha outlet zero; D1/D2/D3 intact | 5.36e-13 L2; 0/0; relErr ≤ 1.4e-12 | **PASS** |
| P8 | true relRes ≤ 1e-9 (≥1 strong dir) | D1 5.33e-14, D2 4.29e-14, D3 9.17e-14 | **PASS** |

**Final decision: PASS_JPP_PRESSURE_BC_PATCH.**

Per task §22: STOP after Post-Review — no tangent oracle (T2/T3) is executed; the next
user-authorized operation is a corrected rerun of the tangent oracle.
