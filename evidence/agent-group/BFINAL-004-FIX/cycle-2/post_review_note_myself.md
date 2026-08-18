# Post-Reviewer independent note (BFINAL-004-FIX cycle-2, DIAGNOSTIC_ONLY)

Reviewer re-executed load-bearing validations in its own invocation (not from
executor reports).

## State check
- pwd == git root == /home/ys/dsH/TO-ANISOTROPIC ; HEAD == ca8a772b… (unchanged).
- `git status --porcelain=v1` sha256 == cycle-2/pre-state == 6ddf3440… (no new
  source changes this round).
- `git diff` sha256 == 247dd339… (production heads zero new diff; only the two
  pre-existing BFINAL-003 M files remain).
- Production heads sha256 == cycle-2/pre-state: solveDiscreteFlowAdjoint.H
  f0c81497…, solveDiscreteFlowAdjointProduction.H 8c4901bf….
- Only source changed = src/stageB6RxDesignOracle.H (untracked diagnostic probe),
  sha 09832dd6… (== step2/step3 records).
- Production `gsenshPressureDrop = -dAlphaDxh*(U&lambda_U)*V` (momentum-only)
  confirmed unchanged; no R_P/pressure-row/assembleWeighted code in production.

## Independent recomputation (own script, not executor's)
`cycle-2/post_review_independent_myself.py` reads raw `stageB6_*.mtx` (N=33600)
and recomputes every metric. Results (reproduce executor/log to the last digit):
- RX-A R_P,alpha relL2 = 6.80e-7 / 2.33e-6 / 6.80e-6 / 2.27e-5 / 6.88e-5, cos=1;
  R_U,alpha relL2=9.24e-8, cos=1.
- rpa_terms identity col1+col2==col3 to 1.29e-26; col3 == rxa_analytic R_P block
  exactly (0.0); weighted terms |div(dphiHbyA_w)|=5.976e-10, |-div(dflux_w)|=
  3.494e-11, total=6.22344977639e-10 == FD oracle.
- RX-B R_P,xh relL2 = 5.61e-7…4.50e-10, cos=1; R_U,xh relL2=1.28e-10.
- RX-C D1/D2/D3 (eps=1e-4): xh-tangent/alpha-tangent/R_U,xd/R_P,xd relL2
  1.95e-5…1.10e-4, cos=1; analytic |a| (175.76/252.25/244.97) matches FD (no D4
  ~1e9 amplification).
- Volume anchor projV_D1 = 0.1555186511144686 (exact), D2=-0.0844966774165889,
  D3=0.0639623768443269.
- §9/§10: λ_P-weighted D_pressure = 15.50%/29.94%/20.61% of D_momentum;
  prod gsensR·d == D_momentum to 2e-9…7e-10.
- FD non-triviality: FD R_P,alpha |b| varies with eps (6.22344979e-10 →
  6.22344766e-10) while analytic is fixed → genuine central FD, not a copy of the
  analytic linearization (no self-reference).

## Independent full binary rerun (load-bearing)
Reran `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF` (absolute path, source
bashrc → FOAM_USER_APPBIN → unset FOAM_SIGFPE) in /home/ys/dsH/b2_case_smoke.
Result: REVIEWER_MTO_RC=0 (EXIT=0), no SIGABRT / calculatedFvPatchField, 477s.
Fresh run reproduces executor log byte-identically:
- StageB5 anchors (R_w/J CLOSED): ExplicitJToracle maxRelL2=3.97863790105e-16;
  momentum 0.000165311363079; P-total 0.000152828914432; rAU avgRatio
  2.50000000002; GateB candB 1.0550908582e-09.
- StageB6 entry fingerprint == exit fingerprint (|U|2=16796.0670769,
  |p|2=7965118.75546, |k|2=11566.7580246, |omega|2=13527763.4515, |nutFrozen|2=0,
  |nuEffFrozen|2=0.00951359211546) → no state re-solve, turbulence frozen.
- RX-A weighted terms + R_P,alpha relL2=6.8011824783e-07/cos=1 at eps=1e-3;
  RX-B/RX-C/§9/§10 all identical to the recorded cycle-2 run.

## Conclusion
Primary hypothesis (pre-assembly weighting order defect, J_P*w vs (J_P·1)*w) is
CONFIRMED by independent code inspection, independent raw-mtx recomputation, and
an independent full binary rerun. Diagnostic evidence is trustworthy. Decision:
PASS (DIAGNOSTIC_ONLY). Recommendation stands: PATCH_RX_PRESSURE_ROW (BFINAL-005).
