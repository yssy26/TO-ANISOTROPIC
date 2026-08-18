# BFINAL-004 step-3-run record — MTO_HF run with stageB6RxDesignOracle

- Stage: B-final · Mode: DIAGNOSTIC_ONLY · Step: step-3-run (executor)
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace) @ `ca8a772b8339a36e7906ea32a7125061c47aa280`
- Run: `/home/ys/dsH/b2_case_smoke_b6` (writable rsync copy of `/home/ys/dsH/b2_case_smoke`,
  old `*.mtx`/`optimization_history.*` excluded; `constant/optProperties` diff = exactly
  `+ stageB6RxDesignOracle true;` — verified with `diff`).
- Binary: `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF` (fresh 4691928 bytes, 22:42,
  48 `stageB6` symbols; invoked by absolute path — PATH resolution would pick the STALE
  `/home/ys/OpenFOAM/ys-7/.../MTO_HF` Aug-15 binary with 0 stageB6 symbols).
- **Outcome: EXIT=134 (SIGABRT) — the run CRASHED in RX-C** (diagnostic-probe defect, see §5).
  Step acceptance (EXIT=0 + RX-A/B/C + §9/§10 tables) **NOT MET**.  Status = TEST_FAILED.

## 1. Environment notes (two launch bugs fixed)

1. `set -u` + `source /opt/openfoam7/etc/bashrc` aborts the whole shell (bashrc references
   unset vars) → first two background attempts produced empty logs.  Removed `set -u`.
2. `which MTO_HF` after sourcing bashrc resolves to the STALE
   `/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin/MTO_HF` (Aug 15, 0 stageB6
   symbols).  The fresh build lives in `FOAM_USER_APPBIN` which is NOT on PATH.  Fixed by
   invoking the absolute path.  Run script: `cycle-1/run_stageb6.sh`.

## 2. BFINAL-003 anchors — BYTE-IDENTICAL (regression green, R_w/J still CLOSED)

| anchor | BFINAL-003 | this run (run_stageb6.log) |
|---|---|---|
| ExplicitJToracle maxRelL2 | 3.97863790105e-16 | 3.97863790105e-16 (L1266, L1656) |
| ExplicitJToracle maxRelL_U | 3.52693970502e-16 | 3.52693970502e-16 |
| ExplicitJToracle maxRelL_P | 3.98062510668e-16 | 3.98062510668e-16 |
| StageB5 momentum (eps 3e-5) | 1.65311e-4 / 0.999999986345 | 0.000165311363077 / 0.999999986345 |
| StageB5 P-total (eps 3e-5) | 1.52829e-4 | 0.000152828916722 |
| StageB5 rAU avgRatio | 2.50000000002 | 2.50000000002 |
| StageB5 GateB candB (eps 1e-5) | 1.05509085820e-09 | 1.0550908582e-09 |
| StageB5 GateB candB (eps 1e-4) | 1.07093812803e-10 | 1.07093812803e-10 |
| StageB4 dPhi_J stats | 2.19206493915e-05 | 2.19206493915e-05 |

No production-math regression.  `R_w/J = CLOSED` stands.

## 3. StageB6 probe entry state (fixed-state guarantee)

```
StageB6 fingerprint: |U|2=16796.0670769 |p|2=7965118.75546 |phi|2=0.0041733366884
 |alpha|2=14587858210.1 |k|2=11566.7580246 |omega|2=13527763.4515 |nutFrozen|2=0
 |nuEffFrozen|2=0.00951359211546 |x|2=141.996485326 |xp|2=141.993260856
 |xh|2=142.223133214 eta5=0.748005161603 del=8
StageB6 support: active=5040 unclipped=5040 deltaAlpha-support=5040 deltaXh-support=5040
```
Exit fingerprint never printed (crash in RX-C).  No `stageB6*` field files were written to
the case (all probe fields are `IOobject::NO_WRITE`; verified `1/` contains 0 stageB6 files).

## 4. Gate results (raw, from run_stageb6.log)

### RX-A (alpha layer, fixed state) — FD values valid; printed analytic-vs-FD metrics INVALID (see §5 bug B1)

- Analytic R_U,alpha = `V*U` (norm 2.099e-6); FD directional norm 1.293e-6.
  relL2=1.906, cos≈0.00066 — **comparison artifact, not physics** (analytic side not
  weighted by deltaAlpha; deltaAlpha = sin(0.271*(celli+1)) signed pattern).
- **R_P,alpha FD = 6.223e-10, stable plateau across all 5 eps**
  (6.2234497883e-10, 6.22344989835e-10, 6.22344999291e-10, 6.22344887681e-10,
  6.22344766468e-10).  → **alpha-layer P-row directional derivative is numerically NON-ZERO
  (tiny: 4.8e-4 × R_U,alpha FD).**
- Analytic R_P,alpha term decomposition: |div(dphiHbyA)|=4.598e-10,
  |-div(dflux)|=5.276e-12, total 4.596e-10; alphaRel=0.4.

### RX-B (xh layer, production alpha interpolation) — THE DECISIVE ROW-SPLIT

- **R_U,xh: CLOSED** — relL2 1.114e-6 → 1.003e-9 → **1.278e-10** (eps 1e-3→1e-5), cos=1,
  |a|=0.176090119582 vs FD 0.1760901196.  xh→alpha chain (production interpolation +
  dAlphaDxh) validated for the U row at machine precision.
- **R_P,xh: FD = 9.71259e-5, stable plateau across all 5 eps**
  (9.71259463259e-05, 9.71259108872e-05, 9.71259077717e-05, 9.71259074173e-05,
  9.71259073862e-05) vs analytic 5.7046e-6: **relL2=1.012, cos=-0.176 — analytic chain
  FAILS to close.**  FD R_P,xh / FD R_U,xh = 9.71e-5 / 0.1761 ≈ **5.5e-4**.

→ **Primary qualitative answer: the fixed-state pressure/continuity-row design derivative
R_P,x is numerically NON-ZERO** (stable multi-eps plateau at both the alpha layer and the
xh layer).  The probe's analytic R_P chain (alpha→relaxed diag→rAU/rAtU→HbyA→phiHbyA→
pEqn.flux→div) does NOT reproduce it — the analytic side is either incomplete (missing
adjustPhi/constrainPressure/setReference/boundary-flux-rebuild derivatives) or has a sign/
formula defect.  Both interpretations are probe-oracle issues; the FD side (full production
rebuild) is the trusted measurement.

### RX-C (raw design D1/D2/D3) — NOT RUN (crash, §5 bug B2)

### §9/§10 (lambda_U/lambda_P weighting, production replica, volume anchor) — NOT RUN

## 5. Probe defects found at runtime (files_allowed diagnostic header only)

- **B1 (RX-A comparison bug):** analytic side compares unweighted `V*U` / `div(dphiHbyA-dflux)`
  against the deltaAlpha-weighted FD.  The plan's own RX-A acceptance is
  `analytic -V*U*deltaAlpha` (i.e. weighted).  Affects only the printed RX-A closure metrics
  (relL2/cos); the FD norms remain valid measurements.  RX-B/RX-C already apply the
  `dAlphaDxh*deltaXh` / `dAlphaDxh*z` weighting correctly.
- **B2 (RX-C crash):** `stageB6_y` (and `stageB6_d`) are created with default `calculated`
  boundary conditions; the Helmholtz tangent solve
  `fvm::laplacian(designFilterFaceMask, yfield) - fvm::Sp(b, yfield) + dfield*b` then aborts:
  "cannot be called for a calculatedFvPatchField on patch inlet of field stageB6_y"
  (FOAM FATAL, run_stageb6.log L1775-1777).  Production `filter_x.H` solves the same
  operator on `xp` which carries proper BCs.  Fix direction: give `stageB6_y`/`stageB6_d`
  the production filter's BCs (or zeroGradient/fixedValue mirror of the `xp` solve) before
  solving; then RX-C and §9/§10 can complete.
- **B3 (RX-B P-row analytic, candidate):** the analytic R_P chain fails (relL2=1.01,
  cos=-0.176).  Likely missing derivatives of production operations that the FD rebuild
  includes: `adjustPhi`, `constrainPressure`, `pEqn.setReference(pRefCell,...)`, and the
  boundary-flux rebuild in R_P.  Needs the Fixer to audit against the FD term decomposition
  (task §6) before any conclusion is drawn from the analytic side.

## 6. Artifacts preserved

- `cycle-1/run_stageb6.log` (1798 lines, full log incl. crash stack)
- `cycle-1/run_stageb6.sh` (fixed env script)
- `cycle-1/optProperties.stageB6.enabled` (case optProperties with the switch)
- `cycle-1/artifacts/stageB6_*.mtx` — rxa_FD_all_eps, rxa_analytic, rpa_terms, rxb_FD_all_eps,
  rxb_analytic populated; rxc_* empty (crash before export)
- `cycle-1/artifacts/stageB5_*.mtx` — full set (anchors)
- Case copy `/home/ys/dsH/b2_case_smoke_b6` left in place for the Fixer/rerun.

## 7. Recommendation routing

`TEST_FAILED` → **Fixer** (diagnostic probe only: fix B2 stageB6_y/d BCs → RX-C + §9/§10
completable; fix B1 RX-A analytic deltaAlpha weighting; audit B3 analytic R_P chain against
the FD decomposition) → rebuild (incremental, header-only) → rerun step-3 → step-4 report.
FD-side facts already secured: R_w/J CLOSED (anchors byte-identical); R_P,x numerically
NON-ZERO (alpha-layer 6.2e-10, xh-layer 9.7e-5, stable plateaus); R_U,x chain closed at
1.28e-10 through the xh→alpha layer; production completeness (§9/§10) UNRESOLVED until the
probe is fixed and re-run.

## 8. Independent mtx sanity check (pure Python, no numpy, 23:0x)

From `artifacts/stageB6_rxa_analytic.mtx` (4N) and `stageB6_rxa_FD_all_eps.mtx` (5×4N), N=33585:
- |anU|=2.099479e-06, |fdU(eps0..4)|=1.292860e-06 (stable) — cos(anU,fdU)=6.6e-4
- |anP|=4.596511e-10, |fdP(eps0..4)|=6.223448..6223450e-10 (stable plateau) — cos≈-2.6e-8
Confirms: (a) FD R_P,alpha plateau is real and non-zero; (b) the printed RX-A cos≈0 is the
unweighted-vs-weighted comparison artifact (bug B1), not physics.
