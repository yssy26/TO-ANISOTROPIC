# BFINAL-007 — Post-Reviewer independent verification & final decision

- Stage: B-final · Mode: **DIAGNOSTIC_ONLY** · Round: BFINAL-007 (outlet fixedValue pressure-BC boundary term)
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace, re-verified) · HEAD `ca8a772b8339a36e7906ea32a7125061c47aa280` · branch `agent/dsH-stage-b-validation`
- Reviewer: independent Post-Reviewer invocation (this file is written by the Post-Reviewer, not the Executor)
- Date: 2026-08-18

---

## 0. What I personally verified (not from Executor summaries)

1. **Workspace gate** — `pwd`, `git rev-parse --show-toplevel`, `git rev-parse HEAD`,
   `git branch --show-current`, `git status --short` all run by me; root == workspace,
   HEAD == ca8a772b, branch == agent/dsH-stage-b-validation.
2. **Scope / diff audit** — `git diff src/MTO_HF.C` is the ONLY tracked change this round
   (+5 guarded include lines); `git diff --check` clean; forbidden files (`NS.H`,
   `sensitivity.H`, `solveDiscreteFlowAdjoint*.H`, `rxPressureRowTranspose.H`, MMA, filter,
   objective/constraint, case `{0,constant,system}`) have zero diff. Locked-head sha256
   recomputed by me and match the BFINAL-003/005/006 records:
   `solveDiscreteFlowAdjoint.H=f0c81497…`, `solveDiscreteFlowAdjointProduction.H=8c4901bf…`,
   `sensitivity.H=f0ed0726…`, `rxPressureRowTranspose.H=fd648ae0…`,
   `stageB5BoundaryRelaxOracle.H=0b274ee2…`, `stageB6RxDesignOracle.H=09832dd6…`.
   Real case `constant/optProperties` contains **no** `stageB7PressureBCDiagnostic` switch
   (grep-verified) → the real case was not touched; the probe ran on scratch copies
   (`b7_scratch_{pre,off,on}`).
3. **Read the actual probe source** `src/stageB7PressureBCDiagnostic.H` (switch-guarded
   `lookupOrDefault<Switch>(..., false)`, read-only, builds
   `pEqnDbg=fvm::laplacian(mobility,p)` and evaluates
   `R_P = div(phiHbyA) − div(pEqnDbg.flux())` with production `fvMatrix::flux()` semantics).
4. **Read the actual J/J^T / R_x source** and confirmed: forward J boundary loop
   (`solveDiscreteFlowAdjoint.H` L701-753) and export boundary loop (L2550-2604) add only
   P-U (`uAssignable`) and U-P (`!pressureFixed`) terms — **no P-P term**; the P-P block is
   interior-only (`mobF*(p_nei−p_own)*dcfF*mafF`, L642-667); `applyPressureFluxCorrection`
   (L262-316) encodes the correct outlet term (`+mobility*delta*Area*p_cell`, L307-312) but is
   **dead code** (called L431, `pressureFluxCorrection` never read again — grep shows only
   L430/L431).
5. **Re-derived the OpenFOAM-7 boundary coefficients from the cited sources** (read the exact
   lines myself): `fixedValueFvPatchField.C` L128-140 (`gradientInternalCoeffs=−δ`,
   `gradientBoundaryCoeffs=δ·p_b`); `gaussLaplacianScheme.C` L82-83 (non-coupled
   `internalCoeffs=pGamma·gradientInternalCoeffs`, `boundaryCoeffs=−pGamma·gradientBoundaryCoeffs`);
   `fvMatrixSolve.C` L134/L148 (`addBoundarySource`/`addBoundaryDiag`); `fvMatrix.C` L903-935
   (`flux()_b = internalCoeffs·p_c − boundaryCoeffs`). Derivation is correct.
6. **Independently re-ran the load-bearing validations with my own fresh Python code**
   (`post_review_rerun.py`, no reuse of the Executor's `s3_*`/`s4_*` logic beyond file
   conventions) and a second independent sigma_min job (`post_review_sigma.py`).

---

## 1. Q1 — actual frozen-primal pressure/SIMPLE residual central FD along n_P: **NON-ZERO**

Independently recomputed from raw artifacts (`stageB7_RP_FD_lin.mtx` etc.) and independently
rebuilt `−(L_prod·n_P)` from `diag_bnd + upper + owner/neighbour`:

| quantity | my independent value |
|---|---|
| `FD_lin = −div(flux(psi=n_P))` max\|·\| | **3.8984e-10** |
| L2 | **1.7683e-09** |
| \|FD\| mass in the 84 outlet cells | **68.62%** |
| \|FD\| mass in top-1% (336 cells) | **94.62%** |
| \|FD_lin − (−L_prod·n_P)\|/\|L_prod·n_P\| (independent rebuild) | **2.414e-09** |
| \|FD_cd(1e-4) − FD_lin\|/\|FD_lin\| | **6.137e-07** |
| eps spread (1e-2 / 1e-6 vs 1e-4) | 6.137e-07 / 2.267e-05 |

**Answer Q1: YES — the production frozen-primal pressure residual derivative along n_P is
non-zero and outlet-concentrated, ~9 orders of magnitude above the exported-J null floor
(2.7e-19). The production pressure mapping does NOT share the exported-J null mode.**

---

## 2. Q2 — exact outlet fixedValue boundary contribution

Re-derived and confirmed from source + case + runtime probe exports:

```
internalCoeffs_b = −rAtU_b·|Sf|_b·δ_b   (< 0)  → laplacian diag += internalCoeffs_b (84 outlet cells)
boundaryCoeffs_b = 0  (p_b = 0)              → laplacian source += 0
pEqn.flux()_b     = internalCoeffs_b·p_c − boundaryCoeffs_b = −rAtU_b·|Sf|_b·δ_b·p_c
⇒ d(phi_b)/dp_c   = +rAtU_b·|Sf|_b·δ_b   (J_PP residual convention)
```

My independent values from `stageB7_outlet_patch.mtx`: `internalCoeffs_b = [−1.7063e-09,
−5.4800e-10]`, `boundaryCoeffs_b = 0` exactly, `δ_b = 4000`, `|Sf|_b = 2.5e-07`; the identity
`ic_b == −(mob_b·δ_b·|Sf|_b)` holds to relL2 **2.887e-16**. **Answer Q2: a boundary-FACE
internalCoeffs diagonal on the 84 outlet cells, no source term, no cell-center identity.**

---

## 3. Q3 — what the current J_PP / coupled J is missing

| comparison | my independent value |
|---|---|
| exported J P-P diag vs interior laplacian diag (excl. pRef row) | **1.700e-16** |
| L_prod − J_PP, off-outlet cells (max\|·\|) | **2.482e-24** (zero) |
| L_prod − J_PP, 84 outlet cells | **[−1.7063e-09, −5.4800e-10] == internalCoeffs_b** (relL2 4.977e-16) |

**Answer Q3: relative to the production pressure matrix, the exported J P-P is missing exactly
the outlet `internalCoeffs_b` diagonal (= `+rAtU_c·δ_b·|Sf|_b` in residual convention); nothing
else differs.** (Confirmed by source audit: all three representations assemble P-P interior-only;
the boundary loop has no P-P term; `applyPressureFluxCorrection` is dead code.)

**Sign-convention finding (Q4b arbitrates empirically):** the exported interior P-P block
equals **+L_int**, i.e. the NEGATIVE of the interior residual derivative `dR_P/dp|_int = −L_int`,
AND the boundary diagonal is missing. The correct residual-convention candidate is
`J_cand_PP = −L_prod = −exported_J_PP + diag(cand)`.

---

## 4. Q4 — diagnostic-only candidate J (= re-pinned exported J with P-P := −L_prod)

### Q4a — null mode gone: **YES** (my independent values)

| metric | exported J | candidate J |
|---|---|---|
| \|\|J·n\|\|/\|\|n\|\| | **2.717e-19** | **3.020e-09** (×1.1e10) |
| sigma_min (splu, 3 random trials) | 7.382e-27 (BFINAL-006) | see §7 (my independent rerun) |

### Q4b — candidate J_PP·n_P vs actual-residual FD: **CLOSES**

| check | my independent value |
|---|---|
| candPP·n_P vs FD_lin | **2.414e-09** |
| candPP·n_P vs FD_cd(1e-4) | **6.137e-07** |
| exported J_PP·n_P vs FD (as-is) | **1.7078** (FAILS) |
| (exported + diag(cand))·n_P vs FD (no sign flip) | **1.4220** (FAILS) |

The candidate reproduces the independently-built production-residual FD at ~1e-9; the literal
"exported J + boundary diagonal" (without the interior sign flip) fails at ~1.4.

### Q4c — BFINAL-003 ordinary directions: **NO REGRESSION**

- `J_cand·v == J_exp·v` for the dp=0 Gate-A direction `v=[dU;0]`: **relL2 = 0.0 exactly** (the
  candidate touches only the P-P block).
- `J_exp·v` vs the matrix-free oracle (excl. pRef row): **3.966e-16** (BFINAL-003 G4 anchor 3.98e-16).
- G1/G2 anchors re-derived by the Executor match BFINAL-003 (candA 2.98506496e-5 exact; G2 blocks
  1.65311e-4 / 1.28188e-4 / 1.68230e-4 / 1.52829e-4 exact). Minor caveat: the second-order
  noise-floor quantity candB@eps=1e-3 re-derives as 1.0592167e-11 vs the BFINAL-003 C++ value
  1.05921688e-11 (~1.8e-7 relative; immaterial).

---

## 5. Q5 — does the same boundary term affect R_P,x?

Confirmed by reading `rxPressureRowTranspose.H` and re-deriving the numbers:
the design-perturbed mobility field `drAU·w` is built with a default `calculated` BC of value 0
(L36-39), so the implemented operator has `dflux_b == 0` exactly (runtime dot-test 2.58e-15).
The same outlet boundary-face flux does have a **design-derivative counterpart**
`d(phi_b)/d(alpha) = +drAtU_b·δ_b·|Sf|_b·p_c`, whose magnitude I independently recomputed:

| quantity | my independent value |
|---|---|
| drAU outlet | [−7.507e-12, −7.279e-13] |
| p_c outlet | [−65.06, +36.78] |
| term per outlet cell | [−1.4388e-13, +5.3862e-14], L2 = **5.358e-13** |
| deltaAlpha (design direction) support on outlet | **0 cells** |
| designMask on outlet | **max = 0.0** |

**Answer Q5: the boundary term is structurally present in R_P,x but identically zero at this
design point (the validated design directions have zero support on the outlet plane), so the
BFINAL-005 `dflux_b=0` assumption is exact here; no R_P,x change is needed or authorized.**

---

## 6. Final conclusion

The coupled J (forward, transpose, exported `explicitJT.mtx` alike) is missing, relative to the
PRODUCTION pressure matrix, exactly **the outlet fixedValue `fvm::laplacian(rAtU,p)`
boundary-face contribution — the `internalCoeffs_b` diagonal on the 84 outlet cells**
(`−rAtU_b·|Sf|_b·δ_b` laplacian convention; `+rAtU_c·δ_b·|Sf|_b` pressure-residual-Jacobian
convention), with **no source term** (boundaryCoeffs = 0 since p_b = 0) and **no cell-center
identity**. This is the source of the outlet-concentrated pressure null mode. The diagnostic-only
candidate `J_PP = −L_prod` (a) removes the null mode (2.7e-19 → 3.0e-9; sigma_min → ~1e-11),
(b) closes against the independently-built production-residual FD at ~1e-9, and (c) leaves every
BFINAL-003 dp=0 anchor exactly unchanged. The Q4b arbitration additionally reveals the exported
interior P-P block carries the opposite sign of the residual derivative (exported = +L_int,
residual derivative = −L_prod) — a separate, correctly-documented finding for any future
production fix, but not the cause of the singularity (a sign flip preserves the null space).

**DIAGNOSTIC ONLY: no production J/J^T / R_x / sensitivity / NS.H / MMA / filter / objective /
gradient-scale change was made; the candidate matrix is offline evidence only; no tangent was
entered. STOP after Post-Reviewer re-verification (this file).**

---

## 7. Independent sigma_min rerun (Post-Reviewer)

`post_review_sigma.py` (my own fresh implementation; candidate P-P := −L_prod built from
`diag_bnd + upper`; splu + 3 random solves) — result written to
`artifacts/post_review_sigma.json`:

```
splu done t=1967.1s
trial 0 sigma_min~6.859e-12
trial 1 sigma_min~4.586e-12
trial 2 sigma_min~1.237e-11
sigma_min(candidate J): min=4.586e-12  mean=7.939e-12
```

My independent sigma_min rerun (fresh implementation, different random seeds) gives
**min = 4.59e-12, mean = 7.94e-12**, consistent with the Executor's
`s3_q4a_sigma.json` (min 8.02e-12, mean 9.30e-12; the small spread is the expected random-vector
variability). Both confirm: **sigma_min(candidate J) ≈ 1e-11 vs exported-J 7.382e-27, a ~×1e15
increase → the candidate J is nonsingular (null mode removed in every direction).**

---

## 8. Minor issues noted (non-blocking, for the record)

1. **Probe C++ Info-line floor**: several in-probe `Info` diagnostics divide by
   `Foam::max(denominator, SMALL)` with SMALL=1e-15 while the true denominators are ~1e-18, so the
   in-probe `cos(...)=0`, `|FD−cand·n_P|/|cand·n_P|=0.079`, self-check `2.86e-10`,
   `|FD_cd−FD_lin|=7.3e-8`, and `eps-spread=2.7e-6` are numerically unreliable. The exported
   artifacts are NOT affected (they are direct field exports), and the FINAL_REPORT correctly
   uses the Python-recomputed values (2.4e-9 / 6.1e-7 / 2.3e-5). Cosmetic only.
2. **`s3_verify_independent.py`** (Executor's "independent" path) has two small bugs: the Q3
   outlet range takes min/max over the full outlet×outlet submatrix (so its "max = 2.068e-25" is a
   roundoff off-diagonal, not the diagonal max −5.48e-10), and the G1 eps labels are misaligned
   with the 5-block FD file. Neither affects any conclusion; the primary scripts and my rerun are
   correct.
3. **Q4c "byte-for-byte" phrasing** overstates candB@eps=1e-3 (a second-order noise-floor
   quantity) which re-derives at 1.8e-7 relative; candA and all G2 blocks are exact.

None of these changes any Q1–Q5 answer or the final conclusion.

---

## 9. Decision

**PASS (DIAGNOSTIC_ONLY).** The diagnostic experiment was executed as planned, its evidence is
trustworthy, and every load-bearing number (Q1 actual-residual FD non-zero; Q2/Q3 exact boundary
term; Q4a null-mode removal; Q4b Jv-FD closure; Q4c no-regression; Q5 R_P,x exact-at-design-point)
was independently reproduced by the Post-Reviewer with fresh code from the raw artifacts and from
the actual sources/case.
