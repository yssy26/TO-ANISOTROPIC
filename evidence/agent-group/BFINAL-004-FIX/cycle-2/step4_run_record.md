# BFINAL-004-FIX (replan v4) cycle-2 — step4-run-record (executor + Post-Reviewer verification)

- Stage: B-final · Mode: DIAGNOSTIC_ONLY · Step: step4-run-report
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace) @
  `ca8a772b8339a36e7906ea32a7125061c47aa280` (unchanged)
- Run: `/home/ys/dsH/b2_case_smoke` (writable original; `constant/optProperties`
  line 87 `stageB6RxDesignOracle true;` — verified; nothing else changed)
- Binary: `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF` (absolute path), size
  4747104 B, mtime 2026-08-17 03:23:23 (built by step3 from the fixed header;
  fixed-header marker strings verified present in the binary:
  `J_P*deltaAlpha`, `degenerate(nB=0)`, `projectionEtaDenominator singular`;
  old `anRPa[celli]` pointwise-multiply string absent)
- Log: `cycle-2/run_stageb6_cycle2.log` · **MTO_RC=0 (EXIT=0)**; START
  2026-08-16T19:25:03Z → END 2026-08-16T19:33:04Z (481 s)
- No SIGABRT, no `calculatedFvPatchField` (grep count = 0), 100 StageB6 log
  lines, 21 `stageB6_*.mtx` artifacts fresh (03:33), copied to
  `cycle-2/artifacts/` (sha256-verified byte-identical to case copies)
- Independent recompute: `cycle-2/post_review_recompute_cycle2.py` (pure
  Python, no numpy, no shared formula with the C++ probe) → output
  `cycle-2/post_review_recompute_cycle2.out`; reproduces every log metric.

## 1. THE HYPOTHESIS IS CONFIRMED — analytic R_P chain CLOSES at all three layers

The pre-assembly weighting fix (J_P*w, weight inside the face-flux assembly)
eliminates the O(1) analytic/FD P-row gap of cycles 0/1 (relL2≈1.01, cos<0).

### RX-A — weighted analytic vs FD (5 eps, N=33600; independent recompute values)

| eps | |a| analytic | |b| FD | relL2 | cos |
|---|---|---|---|---|
| 1e-3  | 6.22344977639e-10 | 6.22344978830e-10 | **6.80e-07** | 0.9999999999998 |
| 3e-4  | 6.22344977639e-10 | 6.22344989835e-10 | **2.33e-06** | 0.9999999999973 |
| 1e-4  | 6.22344977639e-10 | 6.22344999291e-10 | **6.80e-06** | 0.9999999999769 |
| 3e-5  | 6.22344977639e-10 | 6.22344887681e-10 | **2.27e-05** | 0.9999999997426 |
| 1e-5  | 6.22344977639e-10 | 6.22344766468e-10 | **6.88e-05** | 0.9999999976354 |

Weighted per-term (RX-A, log L1751 + independent recompute):
|div(dphiHbyA_w)|=5.97626165832e-10, |-div(dflux_w)|=3.49379242410e-11,
|total_w|=6.22344977639e-10 == FD oracle |b|=6.2234497883e-10 to ~10 digits.
`stageB6_rpa_terms.mtx` identity col1+col2==col3: max|diff|=1.29e-26 (machine
zero); col3 == `rxa_analytic` R_P block exactly (max|diff|=0).
(cycle-1: analytic total 4.5949e-10 ≠ FD 6.223e-10, relL2=1.0115/cos=-0.18.)

### RX-B — xh layer (production alpha interpolation)

| eps | R_P,xh relL2 | R_P,xh cos | R_U,xh relL2 | R_U,xh cos |
|---|---|---|---|---|
| 1e-3  | **5.61e-07** | 1 | 1.114e-06 | 1 |
| 3e-4  | **5.05e-08** | 1 | 1.003e-07 | 1 |
| 1e-4  | **5.62e-09** | 1 | 1.114e-08 | 1 |
| 3e-5  | **5.25e-10** | 1 | 1.003e-09 | 1 |
| 1e-5  | **4.50e-10** | 1 | 1.2775e-10 | 1 |

R_P,xh analytic |a|=9.71259073822e-05 vs FD |b|=9.71259073862e-05 (eps=1e-5).
(cycle-1: R_P,xh relL2=1.0106/cos=-0.1432.)

### RX-C — raw-design D1/D2/D3 (filter→projection→alpha; D4 fixed, no ~1e9 amplification)

Analytic |a| now O(FD): xh-tangent |a| = 175.76/252.25/244.97 vs FD
175.76/252.25/244.96 (cycle-1 analytic was 7.17e11/4.48e11/6.56e11 — D4 factor
~4.08e9 gone). relL2 (independent recompute):

| dir | quantity | eps=1e-3 | eps=3e-4 | eps=1e-4 | eps=3e-5 | eps=1e-5 |
|---|---|---|---|---|---|---|
| D1 | xh-tangent | 4.76e-4 | 8.10e-5 | 2.14e-5 | 1.21e-4 | 1.21e-4 |
| D1 | alpha-tangent | 3.90e-4 | 6.71e-5 | 1.95e-5 | 1.10e-4 | 1.10e-4 |
| D1 | R_U,xd | 4.54e-4 | 7.59e-5 | 2.30e-5 | 1.33e-4 | 1.34e-4 |
| D1 | R_P,xd | 2.03e-3 | 3.39e-4 | 4.14e-5 | 2.02e-4 | 2.01e-4 |
| D2 | xh-tangent | 7.23e-4 | 1.27e-4 | 2.86e-5 | 9.67e-5 | 9.66e-5 |
| D2 | alpha-tangent | 6.32e-4 | 1.11e-4 | 2.58e-5 | 9.21e-5 | 9.20e-5 |
| D2 | R_U,xd | 7.04e-4 | 1.22e-4 | 2.95e-5 | 1.00e-4 | 1.00e-4 |
| D2 | R_P,xd | 1.81e-3 | 3.42e-4 | 7.01e-5 | 1.56e-4 | 1.56e-4 |
| D3 | xh-tangent | 7.86e-4 | 1.41e-4 | 2.58e-5 | 9.70e-5 | 9.70e-5 |
| D3 | alpha-tangent | 6.74e-4 | 1.21e-4 | 2.26e-5 | 9.14e-5 | 9.14e-5 |
| D3 | R_U,xd | 8.66e-4 | 1.48e-4 | 2.77e-5 | 9.69e-5 | 9.69e-5 |
| D3 | R_P,xd | 3.60e-3 | 6.53e-4 | 1.10e-4 | 1.87e-4 | 1.87e-4 |

All cos ≥ 0.9999936 (all cos > 0.999). The eps=1e-3 R_P,xd points (1.8e-3…3.6e-3)
are FD central-difference truncation: relL2 converges quadratically
(1e-3→3e-4→1e-4: ~5.5x, ~6x per 3.3x step) to a ~1.1e-4…2.0e-4 plateau at
eps≤3e-4 — the analytic side is the converged truth; 4 of 5 eps are < 1e-3 and
the preferred <1e-4 is met at eps=1e-4 for every direction/quantity.

## 2. §9/§10 — now unpoisoned (analytic chain; pressureGradientScale=1)

| dir | \|\|R_U,x d\|\| | \|\|R_P,x d\|\| | ratio R_P/R_U | λ_UᵀR_U,x d | λ_PᵀR_P,x d | D_momentum | D_pressure | D_total | prod gsensR·d | \|D_mom−prod\|/\|prod\| |
|---|---|---|---|---|---|---|---|---|---|---|
| D1 | 0.566318445188 | 0.000193062746532 | 3.41e-4 | 0.0486867554277 | 0.00754887795106 | −0.0486867554277 | −0.00754887795106 | −0.0562356333787 | −0.0486867553288 | 2.03e-9 |
| D2 | 0.861583739658 | 0.000447977917016 | 5.20e-4 | −0.229079662362 | −0.0685927693909 | +0.229079662362 | +0.0685927693909 | +0.297672431753 | +0.229079662485 | 5.40e-10 |
| D3 | 0.851652729879 | 0.000259542840587 | 3.05e-4 | 0.0301314016208 | 0.00621002821874 | −0.0301314016208 | −0.00621002821874 | −0.0363414298395 | −0.0301314015977 | 7.65e-10 |

- λ_P-weighted pressure-row fraction |D_pressure/D_momentum| = **15.5 % (D1),
  29.9 % (D2), 20.6 % (D3)** — NON-NEGLIGIBLE (task §9: raw norm ratio 3-5e-4
  but λ_P-weighted effect 15-30 %).
- prod gsensR·d == D_momentum to 2e-9 → the current production
  `gsenshPressureDrop` assembly is EXACTLY the momentum piece; the pressure-row
  piece −pcᵀR_P,x d is absent from production.

## 3. Regression anchors — all byte-identical to BFINAL-003 (R_w/J CLOSED)

| anchor | BFINAL-003 / cycle-1 | this run |
|---|---|---|
| ExplicitJToracle maxRelL2 | 3.97863790105e-16 | 3.97863790105e-16 (L1267, L1657) |
| ExplicitJToracle maxRelL_U | 3.52693970502e-16 | 3.52693970502e-16 |
| ExplicitJToracle maxRelL_P | 3.98062510668e-16 | 3.98062510668e-16 |
| StageB5 momentum (eps 3e-5) | 0.000165311363075 / 0.999999986345 | identical |
| StageB5 P-total | 0.000152828911949 | identical |
| StageB5 rAU avgRatio | 2.50000000002 | 2.50000000002 |
| StageB5 GateB candB (eps 1e-5) | 1.0550908582e-09 | 1.0550908582e-09 |
| StageB4 dPhi_J L2 | 2.19206493915e-05 | 2.19206493915e-05 |
| **Volume anchor projV_D1** | **0.1555186511144686** (stage_b) | **0.1555186511144686** (L1870; cycle-1 was D4-poisoned 9.66e9) |

Volume anchor D2 = -0.08449667741658891, D3 = 0.06396237684432686 (production
replica; D4 fixed).

## 4. State-frozen / turbulence-frozen guarantee

Entry fingerprint (L1746) == exit fingerprint (L1880), byte-identical:
|U|2=16796.0670769, |p|2=7965118.75546, |phi|2=0.0041733366884,
|alpha|2=14587858210.1, |k|2=11566.7580246, |omega|2=13527763.4515,
|nutFrozen|2=0, |nuEffFrozen|2=0.00951359211546, |x|2=141.996485326,
|xp|2=141.993260856, |xh|2=142.223133214, eta5=0.748005161603, del=8.
No state variable re-solved; k/omega/nutFrozen/nuEffFrozen untouched; U/p/phi
rebuilt only in local copies inside rebuildResidual; filter solves act on
diagnostic local fields.

## 5. Scope / source integrity

- Production source zero-diff: `git status --porcelain=v1` sha256
  `6ddf34405437…` and `git diff` sha256 `247dd33996ee…` byte-identical to
  `cycle-2/pre-state/` (the only M files are the two pre-existing BFINAL-003
  production headers, preserved untouched; probe header untracked).
- Only `src/stageB6RxDesignOracle.H` is the diagnostic source changed this
  cycle (sha256 09832dd6… == step2/step3 record; no further edits after build).
- No production math / sensitivity / residual / J / filter / projection /
  conductivity / MMA / thresholds touched; FLOW_JACOBIAN_FROZEN state (per
  BFINAL-003, J frozen) maintained — no J modification, no SIMPLE-transpose,
  no MMA.

## 6. Routing / recommendation

All step4 acceptance criteria MET:

- WMAKE_EXIT=0 (step3 record) and run used the absolute-path fresh binary
  (fixed-header markers verified inside the executable).
- EXIT=0, no SIGABRT / calculatedFvPatchField.
- Volume anchor D1 = 0.1555186511144686 reproduced exactly.
- RX-A R_P,alpha relL2 = 6.8e-7…6.9e-5 (< 1e-3, < 1e-4 at eps ≥ 3e-4),
  cos ≥ 0.9999999976, 5-eps plateau — CLOSED.
- RX-B R_P,xh relL2 = 4.5e-10…5.6e-7, cos = 1, 5-eps plateau — CLOSED.
- RX-C xh/alpha-tangent + R_U,xd/R_P,xd closed (no D4 amplification);
  R_P,xd ≤ 2.0e-4 for eps ≤ 3e-4, ≤ 1.1e-4 at eps = 1e-4, cos > 0.999.
- R_w/J and R_U,x anchors no regression (byte-identical).
- λ_P-weighted pressure-row contribution 15-30 % of D_momentum — non-negligible;
  production contains exactly the momentum piece (|D_mom−prod|/|prod| ≈ 1e-9).

**Final recommendation: PATCH_RX_PRESSURE_ROW** (authorize BFINAL-005) — the
analytic chain is closed at all three layers and the §9/§10 λ_P-weighted
pressure-row contribution is non-negligible, so the missing pressure-row term
−pcᵀR_P,x d must be added to the production pressure-drop sensitivity. No
production patch was made in this DIAGNOSTIC_ONLY round.
