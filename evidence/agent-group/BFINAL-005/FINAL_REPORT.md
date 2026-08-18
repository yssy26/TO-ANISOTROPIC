# BFINAL-005 FINAL_REPORT — Post-Reviewer (independent, post-Fixer)

- Stage: B-final · Mode: **PATCH** · Hypothesis: `BFINAL-005-MISSING-PRESSURE-RX`
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace, verified)
- HEAD: `ca8a772b8339a36e7906ea32a7125061c47aa280` · branch `agent/dsH-stage-b-validation`
- Reviewer: independent Post-Reviewer invocation (post bounded-Fixer round)
- **Decision: PASS_RX_PRESSURE_ROW_PATCH**

> Production R_x now contains both momentum-row and pressure/continuity-row
> design derivatives. No further R_x patch is authorized unless an independent
> fixed-state R_x regression appears.
>
> R_w/J remains frozen from BFINAL-003.

---

## Provenance

| item | value |
|---|---|
| pwd | `/home/ys/dsH/TO-ANISOTROPIC` |
| git root | `/home/ys/dsH/TO-ANISOTROPIC` (== workspace) |
| HEAD | `ca8a772b8339a36e7906ea32a7125061c47aa280` (unchanged) |
| branch | `agent/dsH-stage-b-validation` |
| tracked dirty files | `src/sensitivity.H` (this patch) + `src/solveDiscreteFlowAdjoint.H` + `src/solveDiscreteFlowAdjointProduction.H` (the two LOCKED BFINAL-003 heads) |
| new untracked | `src/rxPressureRowTranspose.H` (authorized helper) |

Locked-head hashes (verified by Post-Reviewer, unchanged from BFINAL-003/004):

- `src/solveDiscreteFlowAdjoint.H` = `f0c814975cf8bc7f98b9f69efa2f0bc01bd251579920bbedda7d266b72927507`
- `src/solveDiscreteFlowAdjointProduction.H` = `8c4901bf8da4edbc7cc1efe73a7583963e4688bba93c8c58539c4c39cae1ad91`
- `src/stageB6RxDesignOracle.H` = `09832dd6ae5a3c7f37201b200bfb19dc0c48f42a3989f9593979c73ca098baec` (BFINAL-004 locked, zero diff)

`cycle-1/` is marked `ABORTED_PREAUTHORIZATION` (preserved, not deleted) and is
not mixed with formal acceptance evidence. Binary `build/bin/MTO_HF` sha256
`744d8e7433fd1ae0f62d65a0d6667072378946adc42116be92e71d2d85c705ac` (built by the
Fixer round, WMAKE_EXIT=0).

---

## Production diff summary

`git diff --stat HEAD` = 3 tracked files: `src/sensitivity.H` (+106/−2, this
patch) and the two LOCKED BFINAL-003 heads (+260/+69, pre-existing). New untracked
`src/rxPressureRowTranspose.H`. `git diff --check` clean. No forbidden module
changed (checked every forbidden entry: only the two locked heads appear).

### `src/sensitivity.H` (DISCRETE branch only; continuous branch byte-identical)

1. Momentum term **UNCHANGED** (renamed for transparency):
   `gsenshPressureDropMomentum = -dAlphaDxh*(U & Uc); *= mesh.V();`
2. New pressure/continuity-row term:
   `gsenshPressureDropPressureRow = -rxPressureRowT*dAlphaDxh;` (NO extra V),
   where `rxPressureRowT = J_P^T pc` is computed by the helper.
3. Total: `gsenshPressureDrop = momentum + pressureRow` (single field feeds the
   EXISTING `filter_chainrule.H` path).
4. `freezeColdFlowForValidation` zeroes BOTH components + total.
5. Guarded (`stageB6RxDesignOracle`) exports of the three xh-level fields
   (`rxpr_prod_gsensh_*.mtx`) and the post-chain raw-design field
   (`rxpr_prod_gsens.mtx`) for G1–G3.

### `src/rxPressureRowTranspose.H` (new, authorized helper)

Computes `rxPressureRowT = J_P^T pc`, the exact transpose of the BFINAL-004
validated forward operator `J_P*w = div(dphiHbyA_w - dflux_w)`, using the locked
reduced-SIMPLE semantics: `drAU = -rAU²/alphaRel`,
`dHbyA = (rAU/alphaRel)((1-alphaRel)U - HbyA)` (NO extra V); `g0 =
flux(fvm::laplacian(one,p))` (scheme-exact face gradient); T1 (HbyA) own/nei +
assignable-boundary pinning; T2 (flux) `-wf*g0*(pc_o-pc_n)*drAU`; boundary
`dflux` derivative = 0 (this case: p fixedValue 0 on outlet, zeroGradient
elsewhere — BFINAL-004-verified).

### Fixer change (bounded, same files)

The prior Post-Reviewer flagged one BLOCKING implementation defect: the entire
`rxPressureRowT` computation was gated behind `stageB6RxDesignOracle`, so
production (switch OFF) had `rxPressureRowT==0` and R_x stayed momentum-only. The
bounded Fixer moved the rebuild + drAU/dHbyA + g0 + T1/T2 loops + the
`rxPressureRowT = T` assignment OUTSIDE the switch guard, keeping ONLY the runtime
dot-test self-check and the diagnostic Info/exports gated. The math is byte-for-byte
unchanged (only guard placement + a comment in sensitivity.H). This is confirmed by
direct code inspection: the assignment loop is unconditional, and the `if
(stageB6RxDesignOracle)` block wraps only the dot-test re-assembly + `Info` +
`FatalError`.

---

## Sensitivity decomposition table (Post-Reviewer independent recompute)

Independent pure-Python (no numpy, `math.fsum`) recomputation from the raw
`.mtx` artifacts — a code path separate from the C++ probe and from the
executor/fixer scripts. `nC = 33600`.

| dir | D_momentum oracle | D_pressure oracle | D_total oracle | G2 D_momentum prod | G1 D_pressure prod | G2 D_total prod | G3 raw-design prod | relErr G1 | relErr G2 | relErr G3 | \|D_pressure/D_total\| |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| D1 | −4.8686755428e-02 | −7.5488779511e-03 | −5.6235633379e-02 | −4.8686755428e-02 | −7.5488779511e-03 | −5.6235633379e-02 | −5.6235633243e-02 | 3.3e-15 | 3.7e-16 | 2.4e-09 | 0.134237 |
| D2 | +2.2907966236e-01 | +6.8592769391e-02 | +2.9767243175e-01 | +2.2907966236e-01 | +6.8592769391e-02 | +2.9767243175e-01 | +2.9767243191e-01 | 8.1e-16 | 1.9e-16 | 5.4e-10 | 0.230430 |
| D3 | −3.0131401621e-02 | −6.2100282187e-03 | −3.6341429840e-02 | −3.0131401621e-02 | −6.2100282187e-03 | −3.6341429840e-02 | −3.6341429802e-02 | 1.4e-14 | 4.0e-15 | 1.0e-09 | 0.170880 |

All signs identical to the oracle. Decomposition self-check
`max|gT − (gM + gP)| = 5.204e-18` (cell-wise). BFINAL-004 locked anchors
(`D_total` −0.0562356333787 / +0.297672431753 / −0.0363414298395; `D_pressure`
−0.00754887795106 / +0.0685927693909 / −0.00621002821874; `D_momentum`
−0.0486867554277 / +0.229079662362 / −0.0301314016208) reproduced to the last
digit (relErr ≤ 1.4e-12 = the rounding of the 13-digit quoted anchors).

Sign convention: `gsenshPressureDropPressureRow = -rxPressureRowT*dAlphaDxh`
contracts to `-pcᵀ J_P(dAlphaDxh·z) = -pcᵀ R_P,x d = -λ_Pᵀ R_P,x d`, exactly the
task §9 requirement `-λ_Pᵀ R_P,x`. The negative sign follows the `J^T λ = g_w`
Lagrangian convention; no sign was flipped to fit D1/D2/D3.

---

## Gate table

| Gate | verdict | evidence |
|---|---|---|
| G1 local pressure-row contraction | **PASS** (relErr 3.3e-15 / 8.1e-16 / 1.4e-14) | production `gsensh_pressurerow·z` == `−pcᵀR_P,x d` |
| G2 full R_x contraction | **PASS** (relErr 3.7e-16 / 1.9e-16 / 4.0e-15) | `gsensh_total·z` == `D_momentum + D_pressure` |
| G3 raw-design projection/filter closure | **PASS** (relErr 2.4e-9 / 5.4e-10 / 1.0e-9) | `rxpr_prod_gsens·d` == `D_total` through EXISTING chain |
| G4 R_w/J no-regression | **PASS** (all BFINAL-003 anchors byte-identical) | Post-Reviewer full rerun (below) |
| G5 trusted non-flow anchors | **PASS** (unchanged) | Post-Reviewer full rerun (below) |
| **production activation (switch OFF)** | **PASS** | Post-Reviewer switch-OFF rerun (below) |

G4 anchors (Post-Reviewer rerun, byte-identical): momentum relL2 =
`1.65311363075e-4`; P-total relL2 = `1.52828911949e-4`; J/J^T dot =
`6.42244770078e-14`; ExplicitJToracle maxRelL_U = `3.52693970502e-16`,
maxRelL_P = `3.98062510668e-16`; GatePR = `relL2=1.03546813403
cos=0.306822860167` (pre-existing R_w-level anchor, identical in both diagnostic
occurrences — NOT a regression, NOT touched by this R_x patch).

G5 anchors (Post-Reviewer rerun): volume proj D1 `0.1555186511144686` /
D2 `-0.08449667741658891` / D3 `0.06396237684432686` (byte-identical);
Anisotropic conductivity `solidK=(21.6 0 0 21.6 0 16.6)`; discrete objective
`thermal=2.7025486001e-11, pressure=7.07793481308e-12`; projection
`eta5=0.748005161603 del=8`.

---

## Reviewer independent reproduction

Two independent load-bearing reruns were performed in THIS invocation, plus an
independent recomputation.

1. **Independent recompute** (`post_reviewer_fixer/independent_recompute.py`):
   pure-Python (no numpy, `math.fsum`), fresh code path reading only the raw
   `cycle-2/fixer/artifacts/*.mtx`. Reproduced `D_momentum` / `D_pressure` /
   `D_total` from `stageB6_lambda.mtx` × `stageB6_rxc_analytic.mtx`, and the
   G1/G2/G3 production contractions. All close to machine precision; BFINAL-004
   locked anchors reproduced (table above). `-λ_Pᵀ R_P,x d` (G1) and
   `-λᵀ R_x d` (G2) independently recomputed for all three directions.

2. **Independent full rerun — switch ON** (`post_reviewer_fixer/run_switchON_postreviewer.log`):
   re-executed `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF` (sha256
   `744d8e74…`) on `/home/ys/dsH/b2_case_smoke` with the exact clean env
   (`clean PATH; source /opt/openfoam7/etc/bashrc; FOAM_USER_APPBIN=…/build/bin;
   unset FOAM_SIGFPE`). `MTO_RC=0`, 0 crash markers. Dot-test
   `relErr=2.57680498181e-15`; xh decomposition `|momentum|L2=0.0373267050442
   |pressureRow|L2=0.0145680567705 |total|L2=0.0516444470747`; all G4/G5 anchors
   byte-identical to the executor's run. Regenerated case-dir `.mtx` exports
   sha256-identical (9/9) to the executor's archived `artifacts/` — determinism
   confirmed. This is the mandated R_w/J no-regression rerun.

3. **Independent full rerun — switch OFF (production)**
   (`post_reviewer_fixer/run_switchOFF_postreviewer.log`): same binary on
   `/home/ys/dsH/b2_case_smoke_prod_fix` (`stageB6RxDesignOracle false`).
   `MTO_RC=0`, 0 crash markers. Verified the fix:
   - dot-test `Info` and xh-fields `Info` lines are ABSENT (gated as intended);
   - the gated `rxpr_prod_*.mtx` exports are NOT regenerated (mtimes remain
     11:15 = pre-run copies);
   - `dDP` range = `[-0.0277320682043, 0.0259232358427]` (the PATCHED total),
     NOT the momentum-only `[-0.0200726738411, 0.0194203165319]`;
   - written `1/gsenshPressureDrop` and `1/gsensPressureDrop` are BYTE-IDENTICAL
     between switch-ON and switch-OFF runs (production computation path
     identical, pressure-row term active);
   - **G3-in-production**: projection of the PRODUCTION-written raw-design
     `1/gsensPressureDrop` onto D1/D2/D3 = −0.056235633243 / +0.297672431912 /
     −0.036341429802, matching the oracle TOTAL at relErr 2.4e-9 / 5.4e-10 /
     1.0e-9 (NOT the momentum-only −0.0486867554 / +0.2290796624 /
     −0.0301314016). This proves production R_x now contains BOTH rows with the
     switch OFF.

4. **Source inspection**: full diff read; `git diff --check` clean; continuous
   branch byte-identical; forbidden modules zero diff; no empirical coefficient
   and no hard-coded D1/D2/D3 in production code; the transpose is the exact
   adjoint of the BFINAL-004 forward operator (T1/T2 face loops + boundary
   pinning match `stageB6RxDesignOracle.H` `assembleWeightedRPa`), pinned by the
   runtime dot-test at 2.58e-15.

---

## Optional exact-adjoint smoke

Not performed. Task §18 makes this optional and explicitly defers the full
frozen-primal FD comparison to the next independent diagnostic (BFINAL-006
tangent oracle). BFINAL-005 only validates that production implements the
BFINAL-004-validated R_x model; it does not claim full total-derivative closure.

---

## Final decision

**PASS_RX_PRESSURE_ROW_PATCH**

G1–G5 all PASS (G1/G2 at machine precision, G3 ~1e-9 through the existing
filter/projection chain), the production activation of the new pressure-row term
is independently confirmed with the diagnostic switch OFF, and the Post-Reviewer
independently reran both load-bearing validations (switch-ON R_w/J no-regression
+ switch-OFF production activation). Only `src/sensitivity.H` (discrete branch)
and the new authorized helper `src/rxPressureRowTranspose.H` changed; forbidden
modules have zero diff; no empirical coefficient or D1/D2/D3 fitting was
introduced; the sign follows `J^T λ = g_w` with no convention flip.

Production R_x now contains both momentum-row and pressure/continuity-row design
derivatives. No further R_x patch is authorized unless an independent fixed-state
R_x regression appears.

R_w/J remains frozen from BFINAL-003.

STOP — per task §24, no BFINAL-006 / tangent oracle / FD campaign / MMA is
automatically started.
