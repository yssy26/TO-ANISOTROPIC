# BFINAL-002 Post-Reviewer evidence (independent invocation)

Post-Reviewer: deepseek-v4-pro (independent of Executor/S4).

## Workspace check
- pwd = /home/ys/dsH/TO-ANISOTROPIC; git root == same; HEAD ca8a772b8339a36e7906ea32a7125061c47aa280 (== task baseline); branch agent/dsH-stage-b-validation. No mismatch.

## Scope audit (git)
- Tracked modifications: exactly 1 file — `src/solveDiscreteFlowAdjoint.H` (+20 insertions, the guarded `stageB5BoundaryRelaxOracle` invocation at EOF).
- New untracked: `src/stageB5BoundaryRelaxOracle.H` (781 lines, switch-guarded).
- Forbidden files (NS.H / sensitivity.H / AdjNS_PD.H / AdjNS_HT.H / costfunction.H / computeObjective.H / filter_x.H / filter_chainrule.H / updateMaterialProperties.H) — `git diff --name-only` on all of them = EMPTY. No production math touched. No SIMPLE-transpose / MMA / R_x patch.
- `git diff --cached` empty (nothing staged).

## Source verification (direct read)
- Probe `rebuildProduction` lambda reproduces NS.H SIMPLE sequence exactly (UEqn.relax → fvOptions.constrain → rAU=1/A → constrainHbyA → fvc::flux → adjustPhi → rAtU(consistent=false) → constrainPressure → pEqn=laplacian(rAtU,p)+setReference → phi=phiHbyA−pEqn.flux()).
- J dphi tangent (solveDiscreteFlowAdjoint.H L628-644 internal, L689-697 boundary) uses UNRELAXED rAUAdj, no relax-source, no constrainHbyA pinning — Candidate A. dp term uses relaxed primalPressureMobility (mixed convention confirmed).
- `discretePrimalMomentum` assembled without `.relax()` (unrelaxed) — rAUAdj = 1/A_unrelaxed.
- `collectDPhiJ`/`collectedDPhiJ` pre-existing (L411-412); probe only reuses them.

## Raw-log verification
- stageB5_run_s4_independent.log StageB5 lines match FINAL_REPORT.md to the printed digit.
- stageB5_run_off.log: 0 "StageB5" lines; anchors (NS UEqn A() avg 158965329.116; |phiRb0−phi| 3.49247878283e-11) reproduced → no-op when disabled.

## Independent recomputation (MY OWN code path)
- Wrote `postreview_recompute.py` (fresh numpy, independent of recompute_stageB5.py / recompute_s4_independent.py) and ran it on artifacts_s3/.
- Reproduced exactly: Gate A momentum relL2 1.653114e-4 / cos 0.999999986; P-int 0.858360 / cos 0.612374; P-bnd 1.145370 / cos 0.152835; P-tot 1.035515 / cos 0.370431 (eps-stable over 5 eps). Gate B Candidate A relL2 0.794617 / normRatio 0.590838 (int 0.790385, bnd 1.387680); Candidate B relL2 1.059e-11→1.055e-9 / cos 1.0 / normRatio 1.0.
- Verified FD eps-independence (affine phi(U) → central FD exact; dphiFD relL2 5e-11→1.5e-9 across eps, cos 1.0).
- Verified relaxation fingerprint: |A|/|FD|=0.590838 (not 0.4); |FD−0.4A|/|FD|=0.8768 (a naive 0.4 rescale would NOT fix it).

## Determinism
- All 8 stageB5_*.mtx md5-identical across artifacts/ (S2), artifacts_s3/ (S3), and /home/ys/dsH/b2_case_smoke (fresh rerun).

## Minor non-blocking notes
- Probe runs twice (solveDiscreteFlowAdjoint.H included from both AdjNS_PD.H L16 and AdjNS_HT.H L101 → two function-scope `static bool`s). Benign: identical output, no state change.
- `evidence/agent-group/BFINAL-002/.venv_s4/` is a full Python venv inside the evidence dir (hygiene only; inside allowed dir).
- X3 relaxation-algebra cross-check is coarse (relL2 4.6e-2, cos 0.9989) but secondary/confirmatory.

## Verdict
PASS (DIAGNOSTIC_ONLY). Diagnostic executed correctly and its evidence is trustworthy. Q-A: boundary-oracle-mismatch REJECTED (O(1) P-row error is a real J defect — unrelaxed rAUAdj + missing constrainHbyA pinning). Q-B: CONFIRMED (unrelaxed rAU causes wrong dphi/dU; Candidate B matches at roundoff). R_x DEFERRED. Recommended next action: A PATCH_J_RELAXATION.
