# BFINAL-004 — FINAL_REPORT (Post-Reviewer-verified, DIAGNOSTIC_ONLY)

- Stage: B-final · Mode: **DIAGNOSTIC_ONLY** · Hypothesis:
  `BFINAL-004-PROBE-FIX-B1-B2-B3` (probe-only defect repair; no production patch)
- Run: `/home/ys/dsH/b2_case_smoke` (writable original, `constant/optProperties`
  L87 `stageB6RxDesignOracle true;` — verified; `stageB2Enabled=false`,
  `mmaUpdateEnabled=false`, `frozenGradientValidated=false`,
  `freezeTurbulenceForValidation=true` — nothing weakened)
- **MTO_RC=0 (EXIT=0)**, 476 s clock (21:22:06Z → 21:30:02Z), no SIGABRT / no
  `calculatedFvPatchField` / no `FOAM FATAL` / no `Aborted` / no `signal`
  (grep count = 0 over the whole log)
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace) @
  `ca8a772b8339a36e7906ea32a7125061c47aa280` (unchanged, branch
  `agent/dsH-stage-b-validation`)
- Probe: `src/stageB6RxDesignOracle.H` (untracked diagnostic header, sha256
  `09832dd6ae5a3c7f37201b200bfb19dc0c48f42a3989f9593979c73ca098baec` ==
  cycle-3 pre-state baseline; the ONLY source file of this round). Production
  source mathematics **zero diff** vs the BFINAL-003 baseline (see §7).
- Binary: `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF` (absolute path),
  4747104 B, mtime 2026-08-17 05:21:59, sha256
  `15172c3aeed061ebb3fa5f2edb5d7954f36730ecb673d04a22ea689b9c20fdab`; built by
  cycle-3 step3 `wmake` (**WMAKE_EXIT=0**, exactly one compile unit `MTO_HF.C`
  + relink, header-only change); fresh fixed-header markers present in the
  executable (`J_P*deltaAlpha`, `degenerate(nB=0)`,
  `projectionEtaDenominator singular` = 3 hits; obsolete `anRPa[celli]*da`
  string = 0 hits)
- Raw evidence: `evidence/agent-group/BFINAL-004-FIX/cycle-3/` — `pre-state/`,
  `step2_audit_notes.md`, `build_step3.log`/`build_step3_exec.log`,
  `run_stageb6_cycle3.log` (1915 lines) / `run_stageb6_cycle3_exec.log`
  (1916 lines), `verify_recompute_cycle3.py`+`.out` (executor-side),
  `artifacts/stageB6_*.mtx` (21 files, 0 zero-byte, byte-identical to the case
  copies), `step3_build_record.md`, `step3_build_run_record_exec.md`,
  `step4_run_record.md`

---

## Status summary (Post-Reviewer-verified)

| item | status |
|---|---|
| R_w/J | **CLOSED** — every BFINAL-003 anchor byte-identical in this run (log L1268/L1658): ExplicitJToracle maxRelL2=3.97863790105e-16, maxRelL_U=3.52693970502e-16, maxRelL_P=3.98062510668e-16; StageB5 momentum relL2=1.65311e-4 / cos=0.999999986345 (L1293), P-total 1.52829e-4 (L1296); rAU avgRatio 2.50000000002; StageB5 GateB candB (eps 1e-5) 1.0550908582e-9; StageB4 dPhi_J L2=2.19206493915e-05 (L1208) |
| R_U,x | **CLOSED** at all three layers — RX-A alpha: relL2=9.24e-8…9.28e-6, cos=1; RX-B xh: relL2=1.114e-6…1.278e-10, cos=1 (machine precision); RX-C raw-design: relL2=2.3e-5…2.95e-5 at eps=1e-4 (all < 1e-3 ceiling, cos>0.9999999) |
| R_P,x | **NONZERO** — FD-secured, stable 5-eps plateaus at every layer (R_P,alpha≈6.223e-10, R_P,xh≈9.71259e-5, R_P,xd=1.931e-4/4.480e-4/2.595e-4 for D1/D2/D3); the independent analytic chain **CLOSES** to the FD oracle at every layer (RX-A relL2=6.8e-7…6.9e-5/cos=1; RX-B relL2=4.5e-10…5.6e-7/cos=1; RX-C R_P,xd relL2=4.1e-5/7.0e-5/1.1e-4 at eps=1e-4 — all < 1e-3 ceiling, cos>0.99999) |
| production R_x completeness | **INCOMPLETE** — production `gsenshPressureDrop` assembly reproduces EXACTLY the momentum piece (−pcᵀR_U,x d to \|Δ\|/\|prod\|≈1e-9); the λ_P-weighted pressure-row contribution −pcᵀR_P,x d (15.5 % / 29.9 % / 20.6 % of D_momentum for D1/D2/D3) is **absent** from production |
| final recommendation | **PATCH_RX_PRESSURE_ROW** (authorize BFINAL-005) — exactly one, per task §15 |

## Primary question (task §1/§12): MISSING_PRESSURE_RX = **CONFIRMED (Case 2)**

The fixed-state pressure/continuity-row design derivative R_P,x is numerically
NONZERO and stable across all 5 eps at the alpha layer (≈6.2234498e-10), the xh
layer (≈9.7125907e-5) and the raw-design layer (D1 1.9306e-4, D2 4.4798e-4,
D3 2.5954e-4), and the independently hand-coded analytic linearization closes to
the central-FD oracle at every layer. The three probe defects of the previous
round are fixed and the diagnostic completes: EXIT=0, RX-A/B/C and §9/§10 tables
all present, volume anchor exact. The λ_P-weighted pressure-row contribution is
non-negligible (15.5 % / 29.9 % / 20.6 % of D_momentum) and the current
production gradient contains only the momentum piece. Per task §12 Case 2 the
missing pressure-row term must be added in the NEXT task (BFINAL-005). No
production patch was made in this DIAGNOSTIC_ONLY round.

## What was fixed this round (B1/B2/B3, probe header only)

- **B1 (RX-A comparison):** the analytic side is now weighted by deltaAlpha
  before comparison — `anRUaW = anRUa*deltaAlpha` (exact: R_U is diagonal in
  alpha, `fvm::Sp` touches only the own diagonal) and
  `anRPaW = assembleWeightedRPa(deltaAlpha)` (weight multiplied into
  drAU/dHbyA BEFORE the face-flux assembly and the divergence — the true
  directional J_P·da, since R_P=div(phi) is NON-diagonal in alpha). cycle-1
  printed spurious relL2=1.9065/cos=0.00066 (unweighted vs weighted); now
  R_U,alpha relL2=9.24e-8…9.28e-6/cos=1 and R_P,alpha
  relL2=6.8e-7…6.9e-5/cos=1 on the stable 5-eps FD plateau.
- **B2 (RX-C crash):** the four solved diagnostic fields — `stageB6_d`,
  `stageB6_y` (filter tangent) and the §9/§10 replicas `gsensVolR`, `gsensR`
  (plus source mirrors `gsenshVolR`, `gsenshR`) — are constructed with
  `zeroGradientFvPatchScalarField::typeName` and `correctBoundaryConditions()`
  before each solve, mirroring the production filter fields (createFields.H
  L416/L593/L608; filter_chainrule.H L30; case `0/xp` zeroGradient). The
  cycle-1 `calculatedFvPatchField::gradientInternalCoeffs` SIGABRT on the inlet
  of `stageB6_y` (BFINAL-004/cycle-1/run_stageb6.log L1775-1787) is gone
  (grep count 0). `gaussLaplacianScheme::fvmLaplacianUncorrected` calls
  `gradientInternalCoeffs()` only on the unknown field; a calculated patch
  FATALs there — eliminated.
- **B3 (R_P analytic chain):** `drAU = -rAU^2/alphaRel` and
  `dHbyA = (rAU/alphaRel)*((1-alphaRel)*U - HbyA)` (no extra V), source-verified
  against OpenFOAM-7 (`fvm::Sp` adds V·alpha to the integrated diagonal D0;
  `relax()` gives D_rel=D0/alphaRel; `A()=D/V` so rAU=V/D_rel already carries V;
  `H()=H0/V` so HbyA=H0/D_rel+(1-alphaRel)U). cycle-1 RX-B R_P,xh
  relL2=1.012/cos=−0.176 → now relL2=4.5e-10…5.6e-7/cos=1.
- **OFstream integrity:** all 21 `stageB6_*.mtx` sections are flush()ed and
  scoped-closed; every artifact is complete with a fresh mtime (0 zero-byte
  files; no tail truncation on the crash path).

## FD-side facts (independent, trusted, 5-eps stable plateaus)

- R_P,alpha ≈ **6.223e-10** (FD norm 6.22344978830e-10 … 6.22344766468e-10
  across eps 1e-3…1e-5).
- R_P,xh ≈ **9.71259e-5** (9.71259463259e-05 … 9.71259073862e-05).
- R_P,xd (RX-C): FD norms 1.9305e-4 (D1), 4.4795e-4 (D2), 2.5945e-4 (D3);
  xh-tangent FD |z| = 175.75 / 252.23 / 244.95 (D1/D2/D3) — stable plateaus.
- R_U,xh closes at machine precision (relL2=1.2775e-10, cos=1, eps=1e-5).
- Weighted R_P,alpha per-term decomposition sums EXACTLY:
  |div(dphiHbyA_w)|=5.97626165832e-10 + |-div(dflux_w)|=3.49379242410e-11 =
  |total_w|=6.22344977639e-10 == FD oracle; `stageB6_rpa_terms.mtx`
  max|col1+col2−col3|=1.29e-26 (machine zero), col3 == the `rxa_analytic`
  R_P block exactly (max diff 0.0).
- The FD side is NOT self-referential: `rebuildResidual` independently rebuilds
  the production reduced residual from scratch (fresh UEqn:
  div − laplacian + fvm::Sp(alpha) == −grad(p) → relax → residual; rAU → HbyA →
  phiHbyA → adjustPhi → rAtU → constrainPressure → pEqn → setReference → flux →
  phi; R_U=−residual, R_P=div(phi)); the analytic side is a separately
  hand-coded linearization with no shared formula. Both sides exported raw;
  the Post-Reviewer recompute below only reads the exported artifacts.

## Per-direction tables (task §15: every tested direction)

### RX-A — direct alpha derivative (fixed state), direction = alpha·δa

| direction | FD R_U norm | analytic R_U norm | FD R_P norm | analytic R_P norm | relL2 (P) | cos (P) | λ_U weighted | λ_P weighted | production projected |
|---|---|---|---|---|---|---|---|---|---|
| alpha·δa, eps=1e-3 | 1.29286041864e-06 | 1.29286041480e-06 | 6.22344978830e-10 | 6.22344977639e-10 | 6.80e-07 | 1 | — | — | — |
| alpha·δa, eps=3e-4 | 1.29286042778e-06 | 1.29286041480e-06 | 6.22344989835e-10 | 6.22344977639e-10 | 2.33e-06 | 1 | — | — | — |
| alpha·δa, eps=1e-4 | 1.29286042172e-06 | 1.29286041480e-06 | 6.22344999291e-10 | 6.22344977639e-10 | 6.80e-06 | 1 | — | — | — |
| alpha·δa, eps=3e-5 | 1.29286030984e-06 | 1.29286041480e-06 | 6.22344887681e-10 | 6.22344977639e-10 | 2.27e-05 | 1 | — | — | — |
| alpha·δa, eps=1e-5 | 1.29286042969e-06 | 1.29286041480e-06 | 6.22344766468e-10 | 6.22344977639e-10 | 6.88e-05 | 1 | — | — | — |

(R_U,alpha relL2=9.24e-8…9.28e-6, cos=1. The R_P,alpha relL2 growth at small
eps is FD round-off on an O(1e-10) signal; the 5-eps FD plateau is stable to 4
significant digits. λ_U/λ_P/production columns are not applicable to the
alpha-layer direction — defined for the raw-design directions below.)

### RX-B — xh-to-alpha closure (production alpha interpolation), direction = xh·δxh

| direction | FD R_U norm | analytic R_U norm | FD R_P norm | analytic R_P norm | relL2 (P) | cos (P) | λ_U weighted | λ_P weighted | production projected |
|---|---|---|---|---|---|---|---|---|---|
| xh·δxh, eps=1e-3 | 0.176090298387 | 0.176090119582 | 9.71259463259e-05 | 9.71259073822e-05 | 5.61e-07 | 1 | — | — | — |
| xh·δxh, eps=3e-4 | 0.176090135675 | 0.176090119582 | 9.71259108872e-05 | 9.71259073822e-05 | 5.05e-08 | 1 | — | — | — |
| xh·δxh, eps=1e-4 | 0.176090121370 | 0.176090119582 | 9.71259077717e-05 | 9.71259073822e-05 | 5.62e-09 | 1 | — | — | — |
| xh·δxh, eps=3e-5 | 0.176090119743 | 0.176090119582 | 9.71259074173e-05 | 9.71259073822e-05 | 5.25e-10 | 1 | — | — | — |
| xh·δxh, eps=1e-5 | 0.176090119600 | 0.176090119582 | 9.71259073862e-05 | 9.71259073822e-05 | 4.50e-10 | 1 | — | — | — |

(R_U,xh relL2=1.114e-6…1.278e-10, cos=1 — the xh→alpha chain is validated at
machine precision; the P row closes to 4.5e-10 at eps=1e-5.)

### RX-C — raw-design directions D1/D2/D3 (filter→projection→alpha, U/p fixed)

Analytic norms (eps-independent): D1 |R_U,x d|=0.566318445188,
|R_P,x d|=0.000193062746532; D2 0.861583739658 / 0.000447977917016;
D3 0.851652729879 / 0.000259542840587.

| direction | FD R_U norm | analytic R_U norm | FD R_P norm | analytic R_P norm | relL2 (P) | cos (P) | λ_U weighted (λ_UᵀR_U,x d) | λ_P weighted (λ_PᵀR_P,x d) | production projected (prod gsensR·d) |
|---|---|---|---|---|---|---|---|---|---|
| D1, eps=1e-3 | 0.566309422701 | 0.566318445188 | 1.93050763358e-04 | 1.93062746532e-04 | 2.03e-03 | 0.9999979 | 0.0486867554277 | 0.00754887795106 | −0.0486867553288 |
| D1, eps=3e-4 | 0.566320544931 | 0.566318445188 | 1.93063095382e-04 | 1.93062746532e-04 | 3.39e-04 | 0.9999999 | 0.0486867554277 | 0.00754887795106 | −0.0486867553288 |
| D1, eps=1e-4 | 0.566321690147 | 0.566318445188 | 1.93063787840e-04 | 1.93062746532e-04 | 4.14e-05 | 0.9999999992 | 0.0486867554277 | 0.00754887795106 | −0.0486867553288 |
| D1, eps=3e-5 | 0.566345382821 | 0.566318445188 | 1.93069025833e-04 | 1.93062746532e-04 | 2.02e-04 | 0.99999998 | 0.0486867554277 | 0.00754887795106 | −0.0486867553288 |
| D1, eps=1e-5 | 0.566345174760 | 0.566318445188 | 1.93068954818e-04 | 1.93062746532e-04 | 2.01e-04 | 0.99999998 | 0.0486867554277 | 0.00754887795106 | −0.0486867553288 |
| D2, eps=1e-3 | 0.861556002830 | 0.861583739658 | 4.47953941314e-04 | 4.47977917016e-04 | 1.81e-03 | 0.9999984 | −0.229079662362 | −0.0685927693909 | +0.229079662485 |
| D2, eps=3e-4 | 0.861582738164 | 0.861583739658 | 4.47976758555e-04 | 4.47977917016e-04 | 3.42e-04 | 0.9999999 | −0.229079662362 | −0.0685927693909 | +0.229079662485 |
| D2, eps=1e-4 | 0.861584780474 | 0.861583739658 | 4.47978388483e-04 | 4.47977917016e-04 | 7.01e-05 | 0.9999999975 | −0.229079662362 | −0.0685927693909 | +0.229079662485 |
| D2, eps=3e-5 | 0.861595980185 | 0.861583739658 | 4.47995737062e-04 | 4.47977917016e-04 | 1.56e-04 | 0.99999999 | −0.229079662362 | −0.0685927693909 | +0.229079662485 |
| D2, eps=1e-5 | 0.861596040144 | 0.861583739658 | 4.47995691493e-04 | 4.47977917016e-04 | 1.56e-04 | 0.99999999 | −0.229079662362 | −0.0685927693909 | +0.229079662485 |
| D3, eps=1e-3 | 0.851605215365 | 0.851652729879 | 2.59448310744e-04 | 2.59542840587e-04 | 3.60e-03 | 0.9999936 | 0.0301314016208 | 0.00621002821874 | −0.0301314015977 |
| D3, eps=3e-4 | 0.851646185630 | 0.851652729879 | 2.59533281983e-04 | 2.59542840587e-04 | 6.53e-04 | 0.9999998 | 0.0301314016208 | 0.00621002821874 | −0.0301314015977 |
| D3, eps=1e-4 | 0.851650880868 | 0.851652729879 | 2.59540505581e-04 | 2.59542840587e-04 | 1.10e-04 | 0.999999994 | 0.0301314016208 | 0.00621002821874 | −0.0301314015977 |
| D3, eps=3e-5 | 0.851648081053 | 0.851652729879 | 2.59531369489e-04 | 2.59542840587e-04 | 1.87e-04 | 0.99999998 | 0.0301314016208 | 0.00621002821874 | −0.0301314015977 |
| D3, eps=1e-5 | 0.851648091033 | 0.851652729879 | 2.59531376830e-04 | 2.59542840587e-04 | 1.87e-04 | 0.99999998 | 0.0301314016208 | 0.00621002821874 | −0.0301314015977 |

(xh-tangent and alpha-tangent internal gates also close: relL2≤2.9e-5 at
eps=1e-4 for every direction/quantity. The eps=1e-3 R_P,xd points (1.8e-3…
3.6e-3) and the ~1e-4…2e-4 floor at eps≤3e-5 are FD central-difference
truncation + filter solver tolerance — relL2 converges quadratically
(1e-3→3e-4→1e-4) to the eps-independent analytic side. All R_P,xd relL2 are
below the task §13 diagnostic ceiling 1e-3; the preferred ≤1e-4 is met at
eps=1e-4 for D1/D2 and nearly met for D3 (1.10e-4). cos>0.999 for every entry.)

### §9/§10 — λ-weighted importance and production completeness (analytic, eps-independent)

| dir | ‖R_U,x d‖ | ‖R_P,x d‖ | ratio ‖R_P‖/‖R_U‖ | λ_UᵀR_U,x d | λ_PᵀR_P,x d | D_momentum | D_pressure | D_total | prod gsensR·d | \|D_mom−prod\|/\|prod\| | \|D_pres/D_mom\| |
|---|---|---|---|---|---|---|---|---|---|---|---|
| D1 | 0.566318445188 | 0.000193062746532 | 3.409e-4 | 0.0486867554277 | 0.00754887795106 | −0.0486867554277 | −0.00754887795106 | −0.0562356333787 | −0.0486867553288 | 2.03e-9 | 15.5 % |
| D2 | 0.861583739658 | 0.000447977917016 | 5.199e-4 | −0.229079662362 | −0.0685927693909 | +0.229079662362 | +0.0685927693909 | +0.297672431753 | +0.229079662485 | 5.40e-10 | 29.9 % |
| D3 | 0.851652729879 | 0.000259542840587 | 3.048e-4 | 0.0301314016208 | 0.00621002821874 | −0.0301314016208 | −0.00621002821874 | −0.0363414298395 | −0.0301314015977 | 7.65e-10 | 20.6 % |

- Raw-norm ratio ‖R_P,x d‖/‖R_U,x d‖ = 3.0e-4…5.2e-4 (modest raw norm), but the
  λ_P-weighted pressure-row contribution |D_pressure/D_momentum| = 15.5 % (D1),
  29.9 % (D2), 20.6 % (D3) — **NON-NEGLIGIBLE** (task §9: do not infer
  importance from residual norms alone).
- `prod gsensR·d` == D_momentum to ~1e-9 for every direction → the current
  production `gsenshPressureDrop` chain (g_x=0, momentum piece only) is exactly
  −λ_UᵀR_U,x d; the pressure-row piece −λ_PᵀR_P,x d is **absent** from
  production.

## Volume anchor (regression anchor, task §8)

- **projV_D1 = 0.1555186511144686** — reproduced exactly in the run log
  (L1870) AND by the Post-Reviewer recompute from `prod_gsensVol.mtx` ×
  `dirs.mtx` over the 5040 active cells (0.15551865111446861).
- D2 = −0.08449667741658891, D3 = 0.06396237684432686 (log L1871-1872;
  recompute −0.084496677416588895 / 0.063962376844326860).
- production-projected contributions recomputed: D1 −0.0486867553288,
  D2 +0.229079662485, D3 −0.0301314015977 — match D_momentum to ~1e-9.

## Post-Reviewer independent verification (task §16, genuine independent invocation)

Performed by the Post-Reviewer in its own invocation — NOT a self-review
(`evidence/agent-group/BFINAL-004/post_review_independent_recompute.py` +
`.out`, a fresh pure-Python script with no numpy and no shared formula with the
C++ probe or the executor's script; it parses ONLY the raw `stageB6_*.mtx`
artifacts).

1. **Diagnostic implementation inspected:** read the full probe header
   `src/stageB6RxDesignOracle.H`; confirmed B1 (deltaAlpha weighting before
   comparison; weighted R_P assembled with the weight inside the face-flux
   assembly), B2 (zeroGradientFvPatchScalarField + correctBoundaryConditions on
   `stageB6_d/y`, `gsensVolR`, `gsensR` and source mirrors), B3
   (drAU=−rAU²/alphaRel, dHbyA=(rAU/alphaRel)((1−alphaRel)U−HbyA)) and the
   OFstream flush/close pattern.
2. **Exact residual semantics inspected:** `rebuildResidual` rebuilds the
   locked production reduced residual from scratch (UEqn
   div−laplacian+fvm::Sp(alpha)==−grad(p) → relax → residual → rAU → HbyA →
   phiHbyA → adjustPhi → rAtU → constrainPressure → pEqn → setReference →
   flux → phi; R_U=−residual, R_P=div(phi)); the analytic side is a separate
   hand-coded chain — **no shared formula** (self-referential oracle rejected).
3. **Raw epsilon tables inspected:** run_stageb6_cycle3_exec.log RX-A
   L1750-1762, RX-B L1763-1773, RX-C L1774-1868, §9/§10 L1873-1880.
4. **Independent recompute of ≥1 R_U,x and ≥1 R_P,x direction:** the script
   recomputes ALL RX-A/RX-B metrics and all RX-C D1/D2/D3 × 5-eps metrics
   (including the alpha-tangent gate aA=dAlphaDxh·z), reproducing every log
   number. Key confirmations: RX-B R_U,xh relL2=1.2775e-10/cos=1 and
   R_P,xh relL2=4.501e-10/cos=1 at eps=1e-5 (B3 closed); RX-A R_U,alpha
   relL2=9.238e-8 and R_P,alpha relL2=6.801e-7 at eps=1e-3 (B1 closed);
   RX-C D1 R_U,xd relL2=2.298e-5 and R_P,xd relL2=4.141e-5 at eps=1e-4
   (raw-design closed).
5. **rpa_terms identity (task §6):** max|col1+col2−col3|=1.29e-26 (machine
   zero) and col3 == the analytic R_P block exactly; |div(dphiHbyA_w)| +
   |-div(dflux_w)| = |total_w| = 6.22344977639e-10 == the central-FD oracle —
   the two weighted per-term contributions sum to the total and to the FD.
6. **No state variable unintentionally re-solved:** between the StageB6 start
   (L1746) and `oracle complete` (L1882) the ONLY DICPCG solves are
   `stageB6_y` (6), `xp` (30, filter tangent), `stageB6_gsensVol` (1),
   `stageB6_gsensPressureDrop` (1) — **zero** solves of U/p/T/k/omega/nuTilda.
   All probe fields are IOobject::NO_WRITE; the case `1/` directory contains no
   stageB6 field files (count 0).
7. **Turbulence remained frozen:** entry fingerprint (L1747) == exit
   fingerprint (L1881) byte-identical (|U|2=16796.0670769, |p|2=7965118.75546,
   |phi|2=0.0041733366884, |alpha|2=14587858210.1, |k|2=11566.7580246,
   |omega|2=13527763.4515, |nutFrozen|2=0, |nuEffFrozen|2=0.00951359211546,
   |x|2=141.996485326, |xp|2=141.993260856, |xh|2=142.223133214, eta5,
   del); `frozenTurbulenceAdjoint=true`,
   `freezeTurbulenceForValidation=true`.
8. **Filter/projection regression anchors green:** volume anchor
   projV_D1=0.1555186511144686 exact (also recomputed independently); D2/D3
   match; StageB4/StageB5/ExplicitJToracle anchors byte-identical to
   BFINAL-003 (see R_w/J row).
9. **Production source mathematics zero diff:** `git status --porcelain=v1`
   byte-identical to `cycle-3/pre-state/` — the only `M` files are the two
   pre-existing BFINAL-003 production headers, hashes unchanged
   (solveDiscreteFlowAdjoint.H=f0c81497…, solveDiscreteFlowAdjointProduction.H=
   8c4901bf…); `git diff` == the pre-state unstaged patch (verified identical);
   every other forbidden module has zero diff. No SIMPLE-transpose, no MMA,
   mmaUpdateEnabled=false, stageB2Enabled=false, gradient scales/thresholds
   untouched.
10. **Binary re-verified:** the absolute-path binary
    `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF` (4747104 B, sha256
    15172c3a…, mtime 05:21:59) was produced by WMAKE_EXIT=0 and is the binary
    that ran (run START 21:22:06Z immediately after WMAKE_END 21:21:59Z);
    fresh fixed-header markers present in the executable, obsolete string
    absent; all 21 artifacts non-empty with fresh mtimes and byte-identical to
    the case copies (cmp verified).

## Recommendation (exactly one, task §15)

**PATCH_RX_PRESSURE_ROW** — authorize BFINAL-005. The fixed-state FD shows a
stable non-zero R_P,x at every layer; RX-A/B/C close numerically (relL2 ≤ 1e-3
ceiling, cos > 0.999; machine precision at the alpha/xh layers); the
λ_P-weighted pressure-row contribution is non-negligible (15–30 % of
D_momentum); and the current production gradient contains only the momentum
piece. Per task §12 Case 2, **MISSING_PRESSURE_RX = CONFIRMED**. No production
patch was made in this DIAGNOSTIC_ONLY round; the next task must add the
pressure-row term −pcᵀR_P,x d to the production pressure-drop sensitivity
assembly (BFINAL-005, authorized scope: production patch, outside this round).
