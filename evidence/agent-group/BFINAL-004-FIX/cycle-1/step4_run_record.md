# BFINAL-004-FIX step4 — run record (MTO_HF on writable /home/ys/dsH/b2_case_smoke)

- Stage: B-final · Mode: DIAGNOSTIC_ONLY · Step: step4-run (executor)
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace) @ `ca8a772b8339a36e7906ea32a7125061c47aa280`
- Run: `/home/ys/dsH/b2_case_smoke` (writable original; `constant/optProperties` diff =
  exactly `+ stageB6RxDesignOracle true;`, verified vs `cycle-1/optProperties.pre-append`)
- Binary: `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF` (fresh 4696344 bytes, mtime
  2026-08-17 00:44:47, built from the B1/B2/B3-fixed header; invoked by absolute path)
- Log: `cycle-1/run_stageb6_fix.log` · **MTO_RC=0 (EXIT=0)** — no SIGABRT, no
  `calculatedFvPatchField` errors (grep count = 0).  Step acceptance (EXIT=0 + tables)
  MET; numerical pass criteria for the analytic P-row chain NOT met (see §4).

## 1. BFINAL-003 regression anchors — BYTE-IDENTICAL (R_w/J still CLOSED)

| anchor | BFINAL-003 / cycle-1 | this run |
|---|---|---|
| ExplicitJToracle maxRelL2 | 3.97863790105e-16 | 3.97863790105e-16 (L1266, L1656) |
| ExplicitJToracle maxRelL_U | 3.52693970502e-16 | 3.52693970502e-16 |
| ExplicitJToracle maxRelL_P | 3.98062510668e-16 | 3.98062510668e-16 |
| StageB5 momentum (eps 3e-5) | 1.65311e-4 / cos 0.999999986345 | 0.000165311363075 / 0.999999986345 |
| StageB5 P-total (eps 3e-5) | 1.52829e-4 | 0.000152828911949 |
| StageB5 rAU avgRatio | 2.50000000002 | 2.50000000002 |
| StageB5 GateB candB (eps 1e-5) | 1.0550908582e-9 | 1.0550908582e-9 |
| StageB4 dPhi_J L2 | 2.19206493915e-05 | 2.19206493915e-05 |

## 2. B1 (RX-A deltaAlpha weighting) — CONFIRMED FIXED

- RX-A R_U,alpha (weighted `anRUaW` vs FD): |a|=1.2928604148e-06, |b|=1.29286041864e-06,
  **relL2=9.23785693523e-08, cos=1**, stable across all 5 eps.  (cycle-1 unweighted artifact
  relL2=1.906/cos=0.00066 is gone.)
- RX-A R_P,alpha (weighted `anRPaW` vs FD): |a|=3.47017865559e-11 vs |b|=6.2234497883e-10.
  Independent mtx recomputation (N=33600): **relL2=1.0115, cos=-0.1799** — NOT closed.
  (The printed relL2=0.0422/cos=-1.75e-5 on the log line is geometrically impossible
  given its own |a|,|b|,|err|: |a-b|/|b| from the printed norms is 1.0115 and cos>1 would
  be required; the mtx-derived numbers are authoritative.  Either way relL2>1e-3.)

## 3. B3 (analytic drAU/dHbyA) — did NOT close the analytic P-row chain

- RX-A R_P,alpha analytic terms: |div(dphiHbyA)|=4.60721246188e-10, |-div(dflux)|=3.29773587666e-11,
  total=4.59494018896e-10; weighted comparison to FD 6.223e-10: relL2=1.0115, cos=-0.18.
- RX-B R_P,xh: analytic |a|=5.97326244089e-06 vs FD |b|=9.71259073862e-05 (5-eps stable plateau),
  **relL2=1.0106, cos=-0.143** — essentially unchanged from cycle-1 (1.012 / -0.176).
- RX-B R_U,xh: |a|=0.176090119582 vs FD 0.1760901196, relL2=1.2775e-10, cos=1 — CLOSED (no regression).
- Conclusion: the B3 formula correction (source-verified no-V drAU=-rAU^2/alphaRel,
  dHbyA=(rAU/alphaRel)((1-alphaRel)U-HbyA)) is in place and correct, but the analytic R_P
  chain does NOT reproduce the FD P-row.  The remaining gap is a probe-oracle analytic-chain
  incompleteness (see plan stop condition #1), not a B3 implementation error.

## 4. NEW defect discovered at runtime (D4) — max(projEtaDenom, SMALL) sign bug

- `projEtaDenom = -3.1385996818e-06` (production E5 DIAG prints the SAME negative value).
- Probe guards `num/Foam::max(projEtaDenom, SMALL)` (L782 projectionTangent etaResp,
  L1035 §9 volEtaCorr, L1121 §10 pressEtaCorr).  With projEtaDenom<0,
  `max(projEtaDenom, SMALL) = SMALL = 1e-15` (since 1e-15 > -3.14e-6), i.e. division by
  +1e-15 instead of -3.14e-6 → ~3.14e9 amplification AND sign flip.
- Production divides by the raw (possibly negative) denominator with a `mag()<=SMALL`
  singularity guard (filter_chainrule.H L81/L108-111) — no such bug there.
- Effect: RX-C analytic chain (xh-tangent |a|=7.17e11 vs FD 175.7, alpha-tangent 1.26e17 vs
  FD 3.26e7, R_U,xd 2.42e9 vs FD 0.566, R_P,xd 9.24e4 vs FD 1.93e-4 — uniform ~4e9 ratio);
  §9 volume anchor projV=9.66e9 (D1) instead of the production-replica value (production
  E5 DIAG gsenshVol sum=1 max=0.0001984; stage_b anchor 0.1555186511144686);
  §10 gsensR/gsenshR magnitudes inflated (max 8.3e5 / 1.9e7).
- The FD sides of RX-C / §9 / §10 are independent of D4 and remain valid:
  xh-tangent FD |z|=175.75 (D1) / 252.25 (D2) / 244.96 (D3); R_U,xd FD 0.566/0.862/0.852;
  R_P,xd FD 1.93e-4 / 4.48e-4 / 2.60e-4 — all stable 5-eps plateaus → **R_P,x numerically
  NONZERO at the raw-design layer too** (ratio to R_U,xd ≈ 3.4e-4, consistent with the
  alpha/xh-layer ratios).
- D4 is a NEW diagnostic-probe defect not in the approved B1/B2/B3 fix list.  Per executor
  rules it is recorded here, NOT patched in this round (no improvisation beyond the plan).

## 5. RX-A/B/C + §9/§10 tables

- RX-A/B/C all ran to completion (EXIT=0); §9/§10 printed (volume-anchor, weighted
  D1/D2/D3 lines).  All `stageB6_*.mtx` non-empty (21 files, incl. rxc_* 1.0-10.7 MB).
- Fingerprint entry == exit (|U|2=16796.0670769, |p|2=7965118.75546, |x|2=141.996485326,
  |xh|2=142.223133214 ...) — no state variable was re-solved; turbulence frozen
  (|nutFrozen|2=0, |nuEffFrozen|2=0.00951359211546 unchanged).

## 6. Artifacts

- `cycle-1/run_stageb6_fix.log` (1909 lines)
- `cycle-1/run_stageb6_fix.sh` (absolute-path env script)
- `cycle-1/optProperties.pre-append` + appended case optProperties (diff = 1 line)
- `cycle-1/artifacts/stageB6_*.mtx` (21 files, copied from the run case)
- `cycle-1/step2_step3_record.md` (fix + build record)
- `evidence/agent-group/BFINAL-004/FINAL_REPORT.md` (this round's report)

## 7. Routing

Step acceptance EXIT=0 + tables + no-crash MET; the numeric pass criteria for the analytic
P-row chain are NOT met (RX-A/B R_P relL2~1.01, cos<0) and a NEW probe defect (D4) was
found.  Status = BLOCKED_BY_NEW_EVIDENCE; final report recommends REPLAN_DIAGNOSTICS per
plan stop condition #1.  FD-side facts (R_P,x NONZERO at all three layers; R_U,xh CLOSED;
R_w/J CLOSED) are secured and survive this round.
