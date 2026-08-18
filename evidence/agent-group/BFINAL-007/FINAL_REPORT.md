# BFINAL-007 FINAL_REPORT — Outlet fixedValue pressure-BC boundary term: derived, proven missing from J_PP, closed by a diagnostic-only candidate J

- Stage: B-final · Mode: **DIAGNOSTIC_ONLY** · Step: S4-final-report (Executor; builds on S1/S2/S3, every load-bearing number re-derived independently in this step)
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace, verified) · HEAD `ca8a772b8339a36e7906ea32a7125061c47aa280` · branch `agent/dsH-stage-b-validation`
- Case: `/home/ys/dsH/b2_case_smoke` · Evidence: `evidence/agent-group/BFINAL-007/`
- Date: 2026-08-18 (Executor, BFINAL-007 S4)

---

## 0. Provenance and scope audit (S4 re-verified)

| item | value |
|---|---|
| pwd / git root | `/home/ys/dsH/TO-ANISOTROPIC` (== workspace) |
| HEAD | `ca8a772b8339a36e7906ea32a7125061c47aa280` (unchanged through the round) |
| branch | `agent/dsH-stage-b-validation` |
| tracked dirty files | ONLY the pre-existing LOCKED BFINAL-003/005 heads + this round's guarded probe include |
| locked-head sha256 (this round, re-verified) | `solveDiscreteFlowAdjoint.H` = `f0c814975cf8bc7f98b9f69efa2f0bc01bd251579920bbedda7d266b72927507`; `solveDiscreteFlowAdjointProduction.H` = `8c4901bf8da4edbc7cc1efe73a7583963e4688bba93c8c58539c4c39cae1ad91`; `sensitivity.H` = `f0ed0726a89950b2c0ff777f9708feabb13d57f4814754edadaf82a30873af09`; `stageB6RxDesignOracle.H` = `09832dd6ae5a3c7f37201b200bfb19dc0c48f42a3989f9593979c73ca098baec` — byte-identical to the BFINAL-005/006 records |
| this round's production diff | `src/MTO_HF.C` +5 lines only: switch-guarded `#include "stageB7PressureBCDiagnostic.H"` after the `AdjNS_PD.H` block (no math/field change) |
| new untracked | `src/stageB7PressureBCDiagnostic.H` (switch `stageB7PressureBCDiagnostic`, default false; verified no-op when off by byte-identical runs) |
| forbidden files | `NS.H`, `sensitivity.H` (unchanged hash), `solveDiscreteFlowAdjoint*.H` (unchanged hash), `rxPressureRowTranspose.H`, MMA, filter/projection, objective/constraint, gradient-scale files, case `{0,constant,system}` — all untouched (git-diff/`git diff HEAD` verified; real case `constant/optProperties` mtime unchanged) |

- DIAGNOSTIC_ONLY: no production J/J^T / R_x / sensitivity / NS.H / MMA / filter / objective change; no tangent entered; no identity-row substitution on outlet cells; candidate J is an offline matrix only, never written to production.
- Runs (S2, byte-identical verification): pre-probe binary vs probe-binary switch OFF vs switch ON on scratch copies `b7_scratch_{pre,off,on}` — all `MTO_RC=0`, all physics outputs byte-identical (only `optProperties` harness path + CSV wall-time field + the 8 probe Info lines differ; `diff -rq` verified). Probe is a read-only no-op when the switch is false/absent.

---

## 1. Q1 — actual frozen-primal pressure/SIMPLE residual central FD along n_P: **NON-ZERO** (outlet-concentrated)

Production semantics rebuilt independently in the probe (NOT the diagnostic J residual definition): `R_P(p) = div(phiHbyA) − div(pEqn.flux())` with `pEqn = fvm::laplacian(rAtU(), p)` (production `fvMatrix::flux()` semantics incl. the outlet `fixedValue p=0` boundary coefficients and the `fluxRequired` face-flux correction; raw face-flux sum, owner `+`/neighbour `−`, NO `/V`; `phiHbyA` independent of `p`). Central FD along the BFINAL-006 null vector `n_P` (P block of `wPrime_TAN_D1.mtx`, normalized; `||w||=2.08447046292e+18`, `||n_P||=0.99998`, `n_P[PREF=cell0]=0`):

| quantity | value (S2 run) | S4 independent re-derivation (from artifacts) |
|---|---|---|
| `FD_lin = −div(flux(psi=n_P))` max\|.\| | **3.8984e-10** | 3.8984e-10 |
| L2 | **1.7683e-09** | 1.7683e-09 |
| \|FD\| mass in the 84 outlet cells | **68.6%** | 0.6862 |
| \|FD\| mass in top-1% (336 cells) | 94.6% | 0.9462 |
| mean\|FD\| outlet / global | 274.5× | 274.5× |
| FD_lin vs `−(L_prod·n_P)` (independent rebuild from diag_bnd+upper) | 2.86e-10 (C++ in-probe) / 2.41e-9 (python) | **2.414e-09** |
| eps-spread, central FD (1e-2/1e-6 vs 1e-4) | 6.1e-7 / 2.3e-5 (exact linearity, roundoff of the non-orthogonal correction) | — |

**Q1 answer: the production frozen-primal pressure residual derivative along n_P is NON-ZERO and outlet-concentrated**, ~9 orders of magnitude above the exported-J null floor (`||J·n||/||n|| = 2.7e-19`): **the production pressure mapping does NOT share the exported-J null mode.** (The 2.86e-10 vs 2.41e-9 self-check spread is the roundoff-level residue of the numerically-evaluated non-orthogonal correction inside `flux()`, which the diag/upper reconstruction cannot include; mesh is perfectly orthogonal, correction ≈ 1e-19 absolute — immaterial to every conclusion.)

---

## 2. Q2 — exact outlet fixedValue boundary coefficient/source contribution (S1 derivation; OpenFOAM-7 sources re-verified)

Production pEqn (`src/NS.H` L1016-1034): `fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)`; `pEqn.setReference(pRefCell,pRefValue)`; `phi = phiHbyA − pEqn.flux()`. Case `0/p`: **outlet = fixedValue uniform 0** (the ONLY `fixesValue()` patch of p; other 7 patches zeroGradient). `p.needReference()==false` (fixedValue outlet present) ⇒ `setRefCell` no-op ⇒ effective `pRefCell=0`, `pRefValue=0` (fvSolution `pRefCell 5600` ineffective — BFINAL-006 S2, re-verified).

Per-face exact values (gaussLaplacianScheme.C L63-64,L82-83; fixedValueFvPatchField.C L130,L138; fvMatrixSolve.C L134,L148 addBoundarySource/addBoundaryDiag; fvMatrix.C L110-125,L144-158; fvMatrix.C L903-936 `flux()`):

```
pGamma        = rAtU_b · |Sf|_b                      (> 0)
internalCoeffs_b = pGamma · gradientInternalCoeffs    = (rAtU_b·|Sf|_b)·(−δ_b) = −rAtU_b·|Sf|_b·δ_b   (< 0)
boundaryCoeffs_b = −pGamma · gradientBoundaryCoeffs   = −(rAtU_b·|Sf|_b)·(δ_b·p_b) = 0                 (p_b = 0)
laplacian diag  += internalCoeffs_b   (addBoundaryDiag, on the 84 outlet cells)
laplacian source+= boundaryCoeffs_b   = 0
pEqn.flux()_b   = internalCoeffs_b·p_c − boundaryCoeffs_b = −rAtU_b·|Sf|_b·δ_b·p_c
d(phi_b)/dp_c   = +rAtU_b·|Sf|_b·δ_b
```

So in the pressure-residual-Jacobian convention `J_PP = dR_P/dp` (`R_P = div(phiHbyA) − div(pEqn.flux())` ⇒ `dR_P/dp = −L_prod`):

> **the outlet-cell diagonal of J_PP must gain `+rAtU_c · δ_b · |Sf|_b`, with NO source term (boundaryCoeffs vanishes, p_b = 0) and NO cell-center identity row** (the term is a boundary-FACE flux linearization, not a Dirichlet row).

Magnitudes this case (runtime probe export; authoritative): `δ_b = 4000.0` (uniform), `|Sf|_b = 2.5e-07` (uniform), `rAtU_b = rAtU_c ∈ [5.48e-7, 1.71e-6]` ⇒ `internalCoeffs_b ∈ [−1.7063e-09, −5.4800e-10]`, candidate `+rAtU_c·δ_b·|Sf|_b ∈ [5.48e-10, 1.71e-9]`. (Note: S1's offline mesh estimate δ≈5.95e3/|Sf|=1.25e-7 differs from the runtime probe values δ=4000/|Sf|=2.5e-7 by a constant factor; runtime values are authoritative and used throughout. Magnitude conclusions unchanged: boundary term ≈ 1–2e-10 per outlet cell vs exported-J outlet diag −1.36e-9..−4.04e-9.) The flux linearization `d(phi_b)/dp_c = +rAtU_b·|Sf|_b·δ_b` matches the dead-code `applyPressureFluxCorrection` (`solveDiscreteFlowAdjoint.H` L307-312: `corr_f = +mobility*delta*Area*p_cell`).

**Q2 answer: production pEqn adds `internalCoeffs_b = −rAtU_b·|Sf|_b·δ_b` to the laplacian diagonal of each of the 84 outlet cells and `boundaryCoeffs_b = 0` to the source** — i.e. `+rAtU_c·δ_b·|Sf|_b` on the outlet-cell J_PP diagonal, no source, no identity row.

---

## 3. Q3 — what the current J_PP / coupled J is missing (numerically proven)

Exported `explicitJT.mtx` (transposed → J, re-pinned row `100800 = P(0)`) P-P block vs the S2-exported production laplacian (`L_int` = interior-only, `L_prod = L_int + addBoundaryDiag`):

| comparison | result (S3) | S4 independent re-derivation |
|---|---|---|
| exported J P-P diag vs L_int diag (excl. pRef row) | **1.70e-16** | 1.700e-16 |
| exported J P-P upper(own,nei) vs L upper(+kf), 96857 faces excl. cell-0 | **1.38e-16**, max\|d\|=1.03e-24 | — (same artifacts) |
| L_prod − J_PP, off-outlet cells | max\|d\| = **2.48e-24** (zero) | 2.482e-24 |
| L_prod − J_PP, the 84 outlet cells | **[−1.7063e-09, −5.4800e-10] == internalCoeffs_b exactly** | [−1.7063e-09, −5.4800e-10] |
| probe laplacian source max\|.\| | 8.9e-20 (boundaryCoeffs = 0) | — |
| identity rows in P block | exactly `[100800]`; row `106400` physical (18 stored/18 nz) | — (rev_nullcheck.log) |

**Q3 answer: relative to the PRODUCTION pressure matrix, the exported J P-P is missing exactly ONE thing — the outlet `internalCoeffs_b` diagonal on the 84 outlet cells** (laplacian convention `−rAtU_b·|Sf|_b·δ_b`; residual convention `+rAtU_c·δ_b·|Sf|_b`); boundaryCoeffs source is 0 (p_b=0); nothing else differs (off-diagonal and off-outlet diagonals agree to 1e-16/1e-24). Source audit (this round, re-read): all three representations (forward `applyDiscreteFlowJ`, transpose `applyDiscreteFlowJT`, exported `explicitJT.mtx`) assemble the P-P block ONLY from interior faces (`kf = mobF·dcfF·mafF`, `solveDiscreteFlowAdjoint.H` L642-667); the boundary loops (L701-753 forward, L2550-2604 export) add only P-U (`uAssignable` dphi) and U-P (`!pressureFixed`) terms, never a P-P term on the `pressureFixed` outlet. The one function encoding the correct term (`applyPressureFluxCorrection` L262-316) is computed-and-discarded dead code (called L431, result never consumed).

**Sign-convention result (Q4b empirical arbitration, S4 re-verified):** the exported interior P-P block equals **+L_int**, i.e. the NEGATIVE of the interior residual derivative `dR_P/dp|_int = −L_int`, AND the boundary diagonal is missing:

```
|exported J_PP·n_P − FD|/|FD|            = 1.708  (as-is fails)
|−exported J_PP·n_P − FD|/|FD|           = 0.307  (bnd still missing)
|(exported+diag(cand))·n_P − FD|/|FD|    = 1.422  (literal "J+term", sign NOT flipped: FAILS)
|(−exported+diag(cand))·n_P − FD|/|FD|   = 2.95e-9  (the closing candidate = −L_prod)
|(−L_prod)·n_P − FD|/|FD|                = 2.41e-9  (independent probe-exports build)
```

The correct residual-convention candidate is **`J_cand_PP = −L_prod = −exported_J_PP + diag(cand)`**, not `+exported_J_PP + diag(cand)`.

---

## 4. Q4 — diagnostic-only candidate J = re-pinned exported J with P-P block := −L_prod (offline only, never written to production)

Candidate built offline (S3/S4): same U rows/cols as the exported J; P-P block REPLACED by `−L_prod` (built two ways, identical to 2.4e-10: from `−expPP + diag(cand)` and directly from probe exports); row `100800` re-pinned to identity. Candidate differs from the exported J ONLY in the P-P block (`J_cand·v == J_exp·v` for any `dp=0` direction: **relL2 = 0.000e+00 exactly**).

### Q4a — null mode gone: **YES**

| metric | exported J | candidate J |
|---|---|---|
| \|\|J·n\|\|/\|\|n\|\| (n = normalized `wPrime_TAN_D1`) | **2.717e-19** | **3.020e-09** (×1.1e10) — S4 re-derived: 2.717e-19 → 3.020e-09 |
| sigma_min (splu random solves, 3 trials) | 7.382e-27 (BFINAL-006 Post-Reviewer) | **8.02e-12 / 1.06e-11 / 9.32e-12** (mean **9.30e-12**, ×1.3e15) |

### Q4b — candidate J_PP·n_P vs S2 actual-residual FD: **CLOSES**

| check | relL2 |
|---|---|
| candidate J_PP·n_P vs exact FD_lin (`−L_prod·n_P`) | **2.95e-9** (candPP=−expPP+diag(cand)) / **2.41e-9** (candPP=−L_prod direct) — both S4-confirmed |
| candidate J_PP·n_P vs central-difference FD (eps=1e-4) | 6.14e-7 (= the central-diff FD's own roundoff vs FD_lin) |
| sign identity: exported J_PP·n_P + FD_lin == diag(cand)·n_P | 9.59e-9 |
| candidate P-P diag vs −(production L diag) (excl pRef) | 2.42e-10 (cell-vs-boundary mobility interpolation residue) |

The candidate J (P-P = −L_prod) reproduces the independently-built production-residual FD along n_P at ~1e-9. The literal "exported J + diag(cand)" (without the interior sign flip) FAILS at relL2=1.42 and is NOT the residual Jacobian; the plan's `J_PP = −laplacian` convention (S1 §3) is the correct one.

### Q4c — BFINAL-003 G1/G2 anchors: **NO REGRESSION** (reproduced to the digit from raw artifacts)

- `J_cand·v == J_exp·v` for the Gate-A direction `v=[dUdir; 0]` (dp=0): **relL2 = 0.000e+00 exactly** — the candidate changes ONLY the P-P block, which no dp=0 direction exercises.
- `J_exp·v` vs the matrix-free oracle `stageB5_gateA_Jv.mtx` (excl. pRef row): **3.966e-16** (BFINAL-003 G4 ref 3.98e-16). (Including the pinned pRef row the ratio is 2.9e-8 — that row is a harness pin, not physics; the ex-PREF comparison is the BFINAL-003 anchor.)
- **G1 (Gate B, per-face dphi)**: candA relL2 = **2.985065e-5** cos=1 (BFINAL-003: 2.98506496e-5); candB = **1.059217e-11** cos=1 (1.05921670e-11) at eps=1e-3.
- **G2 (Gate A) blocks** @ eps=1e-3: momentum **1.653114e-4**/cos 0.9999999863; P-internal **1.281877e-4**/0.9999999924; P-boundary **1.682302e-4**/0.9999999863; P-total **1.528289e-4**/0.9999999883 — all identical to BFINAL-003's accepted anchors (1.65311e-4 / 1.28188e-4 / 1.68230e-4 / 1.52829e-4) and stable across eps ∈ {1e-3, 3e-4, 1e-4, 3e-5, 1e-5}.

---

## 5. Q5 — does the SAME boundary term affect R_P,x? **Structurally YES; numerically ZERO here (design support excludes the outlet); no R_P,x change needed at this design point**

- R_P,x (BFINAL-004/005 closed, `rxPressureRowTranspose.H` L8-13): `J_P·w = div(dphiHbyA_w − dflux_w)`, `dflux_w = flux(fvm::laplacian(drAU·w, p))`, `drAU = −rAU²/alphaRel` (alphaRel = 0.4, from `mesh.equationRelaxationFactor("U")` — verified).
- The design-perturbed mobility field `drAU·w` is created with default **calculated** BC value **0** (`rxPressureRowTranspose.H` L36-39; `calculatedFvPatchField::updateCoeffs` empty ⇒ boundary value stays 0), so the interpolated boundary `gamma_b = 0` and **`dflux_b == 0` exactly** in the operator as implemented — runtime dot test `relErr = 2.57680498181e-15` (Log.stageB7.on.txt L1900), `|T|L2=0.000256656999292`, `|J_P·w|L2=6.22344977639e-10`, `w-support=5040`.
- The SAME outlet fixedValue internalCoeffs mechanism has a design-derivative counterpart (frozen p_c): `d(phi_b)/d(alpha) = +drAtU_b·δ_b·|Sf|_b·p_c` per outlet cell. Magnitude (S3, re-derived in S4): `drAtU_b ∈ [−7.51e-12, −7.28e-13]`; `p_c ∈ [−65.06, +36.78]` at the outlet cells; term per cell ∈ **[−1.44e-13, +5.39e-14]**, L2 = **5.36e-13**; ratio to `|J_P·w|L2` = **8.6e-4** (~0.1%).
- **Decisive:** the weighted design directions actually validated (`stageB6_deltaAlpha` for RX-A, and `1/designMask` from the case) are **identically ZERO on all 84 outlet cells** (S3/S4: `deltaAlpha` outlet nnz = 0; `designMask` outlet max = 0.0), so the design-weighted boundary term is **identically 0** in the RX-A/B validated directions; BFINAL-005 closure (dot test 2.6e-15; R_P,alpha relL2 ~1e-6 vs FD) is **unaffected**, and the `dflux_b=0` assumption is **exact for this case**. If the design domain ever extended to the outlet plane, the operator would need the small (~0.1%) `+drAtU_b·δ_b·|Sf|_b·p_c` term — beyond this design point and NOT authorized this round.

**Q5 answer: the same boundary-face term does have a small design-derivative counterpart in R_P,x, but the current R_P,x is exact at this design point because the design direction has zero support on the outlet cells; no R_P,x change is authorized or needed here.**

---

## 6. Final conclusion — the exact missing term in J

The coupled J (forward, transpose, and exported `explicitJT.mtx` alike) is missing, relative to the PRODUCTION pressure matrix, exactly **the outlet fixedValue `fvm::laplacian(rAtU,p)` boundary-face contribution: the `internalCoeffs_b` diagonal on the 84 outlet cells** (= `−rAtU_b·|Sf|_b·δ_b` laplacian-matrix convention; `+rAtU_c·δ_b·|Sf|_b` pressure-residual-Jacobian convention), with **NO source term** (boundaryCoeffs = 0 because p_b = 0) and **NO cell-center identity row**. This missing term is the source of the outlet-concentrated pressure null mode (BFINAL-006: sigma_min ≈ 7.4e-27, `||J·n||/||n|| ≈ 2.7e-19`, 89.5% of P-mode mass at the outlet plane).

The diagnostic-only candidate `J_cand_PP = −L_prod = −(L_int + bnd)` (= `−exported_J_PP + diag(cand)`, the residual Jacobian):
1. **removes the null mode** — `||J·n||/||n||`: 2.7e-19 → 3.0e-9 (×1.1e10); sigma_min: 7.4e-27 → 9.3e-12 (×1.3e15);
2. **closes against the independently-built actual-residual FD** at ~1e-9 (Q4b), where the FD was built from the production `NS.H`/`pEqn` boundary/SIMPLE semantics — not from the diagnostic J's residual definition;
3. **leaves every BFINAL-003 anchor byte-for-byte unchanged** (all validated directions have dp=0; the candidate touches only the P-P block).

The exported interior P-P block is the NEGATIVE of the interior residual derivative (empirically arbitrated in Q4b): the correct residual-convention fix is `−L_int + diag(cand)`, NOT `+L_int + diag(cand)` (which fails closure at relL2=1.42). R_P,x (Q5) does not need the boundary term at this design point (zero design support on the outlet cells).

This is **DIAGNOSTIC ONLY**: no production J/J^T / R_x / sensitivity / NS.H / MMA / filter / objective / gradient-scale change was made; the candidate matrix is offline evidence only. Per the round's STOP rule: **Post-Reviewer must independently re-verify (re-run at least one actual-residual FD check and one candidate-J null-mode/Jv-FD check, review the boundary derivation and the cited OpenFOAM source lines, and inspect the raw logs/matrices) before this round is closed; no subsequent task is auto-started.**

---

## 7. Post-Reviewer rerun recipe (everything re-derivable from artifacts)

```bash
cd /home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-007
PY=/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-002/.venv_s4/bin/python3
# Q1+Q3+Q4a+Q4b+Q4c+Q5 from raw artifacts (fresh code path, ~8 s):
$PY s4_executor_verify.py          # -> s4_executor_verify.log/.json
# Q3/Q4 detail + G1/G2 anchors from raw artifacts:
$PY s3_q3_q4.py                    # Q3 diff, Q4a nullvec, Q4b closure, Q4c G1/G2
$PY s3_verify_independent.py       # independent rebuild path (Q3/Q4a/Q4b/Q4c G1)
$PY s3_q5_rpx.py                   # Q5 boundary-term magnitude
# sigma_min of candidate J (splu, ~30 min; 3 random trials):
$PY s3_q4a_sigma.py                # -> s3_q4a_sigma.json (mean 9.30e-12)
# actual-residual FD probe run (rebuild + MTO_HF on a scratch case with
# stageB7PressureBCDiagnostic true; ~8 min): see S2_probe_evidence.md §2/§5.
```

Raw artifacts: `evidence/agent-group/BFINAL-007/artifacts/` (`stageB7_outlet_patch.mtx`, `stageB7_laplacian_{diag,diag_bnd,upper,source}.mtx`, `stageB7_RP_FD{,_lin,_all_eps}.mtx`, `stageB7_candidate_diag.mtx`, `stageB7_nP.mtx`, `stageB7_rebuilt_rAU.mtx`, `stageB7_boundary_all_patches.mtx`, `logs/Log.stageB7.{pre,off,on}.txt`); case anchors `stageB5_*.mtx`, `stageB6_deltaAlpha.mtx`; null vector `BFINAL-006/cycle-2/artifacts/wPrime_TAN_D1.mtx`; BFINAL-006 Post-Reviewer null/sigma evidence in `BFINAL-006/cycle-2/post-review/`.

## 8. Files (this round)

- Evidence: `S1_boundary_derivation.md`, `s1_verify_boundary.py|.log`, `S2_probe_evidence.md`, `s2_compare_runs.sh|.log`, `s2_verify_artifacts.py|.log`, `artifacts/*` (17 .mtx + 3 logs), `s3_q3_q4.py|.log|.json`, `s3_q4a_sigma.py|.log|.json`, `s3_q4c_anchors.py|.log`, `s3_q5_rpx.py|.log`, `s3_verify_independent.py|.log`, `s3_summary.json`, and this report's S4 verification: `s4_executor_verify.py|.log|.json`.
- Source: `src/stageB7PressureBCDiagnostic.H` (new, switch-guarded default-off read-only probe); `src/MTO_HF.C` (+5 guarded include lines). No forbidden file touched.
