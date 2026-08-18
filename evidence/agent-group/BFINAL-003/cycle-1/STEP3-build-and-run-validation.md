# BFINAL-003 step: build-and-run-validation — execution record

Status: EXECUTED (EXIT=0, all Gate markers present, artifacts exported)
Date: 2026-08-16 (session clock)
Workspace: /home/ys/dsH/TO-ANISOTROPIC, git root == workspace, HEAD ca8a772b8339a36e7906ea32a7125061c47aa280
branch agent/dsH-stage-b-validation (unchanged by this step).

## Build

- Re-verified wmake clean/current: sources mtime <= 19:50, binary
  build/bin/MTO_HF 20:17:50 (4237328 B, newer than all sources). WMAKE_EXIT=0,
  nothing to relink. Log: `cycle-1/build3-validation.log` (previous steps'
  build.log/build2.log preserved; build2.log already recorded the successful
  relink with the correct env order: source bashrc FIRST, then export
  FOAM_USER_APPBIN).

## Run

- Command (background, cwd = /home/ys/dsH/b2_case_smoke):
  `export PATH=...; source /opt/openfoam7/etc/bashrc; unset FOAM_SIGFPE;
  /home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF | tee cycle-1/run.log`
- Case /home/ys/dsH/b2_case_smoke with `stageB4JacobianProbe true;
  stageB5BoundaryRelaxOracle true;` (constant/optProperties lines 79-80,
  unchanged from BFINAL-002 s3run).
- fvSolution `equations { U 0.4; ... }` -> alphaRel = 0.4 read from the
  actual equation relaxation factor (`mesh.equationRelaxationFactor("U")`
  guarded by `mesh.relaxEquation("U")`); no hard-coded 0.4/0.6/2.5 in
  production math (only comments).
- Result: RUN_EXIT=0, ExecutionTime = 465.16 s (single B4 path).
- Determinism: the StageB5 oracle block runs once per inclusion site
  (AdjNS_HT.H lines 1286-1348 and AdjNS_PD.H lines 1676-1738 of run.log);
  the two blocks are byte-identical (`diff` clean) — no state drift.
- Raw log: `cycle-1/run.log` (1770 lines). No FOAM FATAL / errors (only
  routine primal continuity diagnostics and the pre-existing
  "residual > tolerance (continuing, diagnostic)" reduced-adjoint notices).
- Reduced-adjoint notices: identical structure to BFINAL-002 (4 occurrences in
  each log). thermalCoupling relRes=1 unchanged; pressureDrop relRes changed
  3.4765473443e-09 (BFINAL-002) -> 36.8066580805 because the diagnostic loads
  the stale case-dir file explicitSol_pressureDrop.mtx (computed for the OLD
  operator) and measures its residual under the CORRECTED (patched) operator —
  expected consequence of the J^T change, and corroborates that the operator
  changed. These are load-and-report diagnostics, not gates.

## Exported artifacts (cycle-1/artifacts/, copied from the case dir after the run)

stageB5_cell_type.mtx, stageB5_dUdir.mtx, stageB5_face_type.mtx,
stageB5_gateA_FD_all_eps.mtx, stageB5_gateA_Jv.mtx, stageB5_gateB_dphi_A.mtx,
stageB5_gateB_dphi_B_all_eps.mtx, stageB5_gateB_dphi_FD_all_eps.mtx,
explicitJT.mtx (467508907 B). md5 (same names in the case dir):
cell_type 00a728a2dd1b8d543927eda9afe0eacb, dUdir c4630c928f4b913ce142f56ef5e34f86,
face_type 7f0f0ce89c4eb193d11a08c4f385ee20, gateA_FD_all_eps
db67751fd5571ba4f81c60dc028b97eb, gateA_Jv 0343172f49744e99fd5968da50de6dd9,
gateB_dphi_A 6e57baee2be077b85aa6379d4f9b17b0, gateB_dphi_B_all_eps
6b33d9580080a79ea87f40363e13f4dd, gateB_dphi_FD_all_eps
c0c3c59605876d9f2a16a230f2105876.
stageB5_gateA_Jv.mtx and stageB5_gateB_dphi_A.mtx sizes differ from BFINAL-002
(2943624 vs 2943858; 2186030 vs 2344860) — expected: they export the
CORRECTED (Candidate-B) operator dphi.

## Gate results (all from run.log; both inclusion sites identical)

### G1 — Candidate-B regression (StageB5 Gate B; "Candidate A" = patched J dphi vs unchanged production FD; "Candidate B" = unchanged BFINAL-002 reconstruction)

| eps | candA relL2 | candA cos | candA normRatio | candB relL2 | candB cos |
|---|---|---|---|---|---|
| 1e-3  | 2.98506496345e-05 | 0.999999999556 | 1.00000150089 | 1.05921688465e-11 | 1 |
| 3e-4  | 2.98506498901e-05 | 0.999999999556 | 1.00000150089 | 3.47047269226e-11 | 1 |
| 1e-4  | 2.98506506028e-05 | 0.999999999556 | 1.00000150089 | 1.07093812803e-10 | 1 |
| 3e-5  | 2.98506477103e-05 | 0.999999999556 | 1.00000150090 | 3.53850848932e-10 | 1 |
| 1e-5  | 2.98506483110e-05 | 0.999999999556 | 1.00000150088 | 1.05509085820e-09 | 1 |

Gate verdict: PASS (hard gate relL2 < 1e-3 met with 33x margin; cos >
0.999999; normRatio ~1.000). "Prefer <1e-6" not met — explanation below.

G1 degradation note vs the BFINAL-002 candB reference (1.06e-11): the patched
J (candA) matches the production FD at relL2 = 2.985e-5, stable to 10
significant digits across all 5 eps (converged plateau, not FD truncation).
candB (the unchanged reconstruction) still matches the same FD at roundoff
(1.06e-11..1.06e-9) — the reconstruction is built from the probe's own
relaxed-H FD basis, so it agrees with the FD by construction. The J's
matrix-free dphi differs from that reconstruction at the ~3e-5 level, which is
the same order as the J's momentum-row closure (1.653e-4) and the Gate-A
P-total closure (1.528e-4) — i.e., the intrinsic frozen-turbulence
J-vs-production-rebuild assembly accuracy of this code, not the O(1) relaxation
defect (which is gone: 0.795 -> 2.985e-5). If sub-1e-6 is later required, the
next step would be a micro-audit of the J's rAU*dH interpolation vs the
probe reconstruction basis; it is not needed for any gate here.

### G2 — full B-F2 R_w closure (StageB5 Gate A, J*v vs boundary-consistent central FD of the production residual; 5 eps)

| block | BFINAL-002 (defective) | this run (eps=1e-3) | all 5 eps |
|---|---|---|---|
| momentum | relL2 1.653e-4, cos 0.999999986 | relL2 1.65311e-4, cos 0.999999986345 | identical at all eps (no regression) |
| P-internal | relL2 0.85836, cos 0.612 | relL2 1.28188e-4, cos 0.999999992443 | 1.2819e-4 @ all eps |
| P-boundary | relL2 1.14537, cos 0.153 | relL2 1.68230e-4, cos 0.999999986268 | 1.6823e-4 @ all eps |
| P-total | relL2 1.03551, cos 0.370 | relL2 1.52829e-4, cos 0.999999988323 | 1.5283e-4 @ all eps |

Gate verdict: PASS. P-total relL2 = 1.528e-4 < 1e-3 (hard gate; preferred
<1e-4 not met by 1.5x — at the same scale as the momentum closure, i.e. the
probe's FD conditioning floor), cos = 0.99999999 > 0.999. Momentum not
regressed (1.653e-4, byte-identical to BFINAL-002). The O(1) P-row mismatch is
GONE (1.036 -> 1.528e-4); per the plan an O(1) P-row mismatch would be an
automatic FAIL — it is absent.

### G3 — J/J^T transpose consistency

- BlockDot (thermalCoupling and pressureDrop, identical): PU=2.07803300324e-13
  PP=2.6655422681e-14 UU=7.3095183677e-14 UP=1.79609743792e-12
  boundaryU=1.80391921425e-13 — all <= 1.8e-12 <= 1e-10.
- Reduced cold-flow operator transpose dot-test: (thermalCoupling)
  max relative error = 6.42244770078e-14; (pressureDrop) = 6.42244770078e-14.
- Thermal (T-adjoint, untouched) matrix transpose dot-test:
  1.09689091013e-13 (byte-identical to BFINAL-002).
Gate verdict: PASS. The transpose (applyDiscreteFlowJT and the production
applyProdFlowJT) corresponds to the corrected Candidate-B forward operator —
BlockDot PU now validates the new P-U transpose block (was 7.99e-14 with the
old unrelaxed block, now 2.08e-13 with the corrected block; still roundoff).

### G4 — explicit/matrix-free oracle

ExplicitJToracle: maxRelL2=3.97863790105e-16, maxCosErr=2.22044604925e-16,
maxRelL_U=3.52693970502e-16, maxRelL_P=3.98062510668e-16 (near machine
precision; ~3.3e-16 reference). Gate verdict: PASS. The explicit CSR oracle
carries the identical alphaRel + direct (1-alphaRel) corrections, so the CSR
J^T and the matrix-free J^T agree at roundoff.

### G5 — trust anchors

- phiRb0: |phiRb0-phi|L2 = 3.49247878283e-11 (byte-identical to BFINAL-002).
- StageB5 momentum anchor (frozen-phi FD): relL2=1.65312e-4, cos=0.999999986345.
- StageB5 rAU: relaxFactorU(fvSolution)=0.4; rAUAdj avg=7.05844426338e-07;
  rAUrel avg=2.82337770533e-07; avgRatio(rAUAdj/rAUrel)=2.50000000002
  (byte-identical; the 0.4 mobility factor and 2.5 ratio are DERIVED from the
  equation relaxation factor, not hard-coded).
- StageB5 selfcheck (recomputed dH boundary + internal div vs outJ[P]):
  relL2=0 (the probe's corrected-semantics self-consistency check).
- Stage B4 T2 (P-P pressure-perturbation continuity, J_P vs
  div(mob*snGrad(dp)*magSf)): relL2 = 7.51176551272e-21 (machine precision,
  byte-identical to BFINAL-002) — J_PP mathematically unchanged.
- Stage B4 T1 (momentum U-U, dp=0): relL2=0.00852323456072,
  cos=0.999963676577 (byte-identical to BFINAL-002) — J_UU unchanged.
- Candidate-B reconstruction byte-identical to BFINAL-002 across all 5 eps.
- Gate-A FD plateau stable to >=9 digits across all 5 eps.

### Legacy probe notes (expected, NOT production regressions) — for Post-Reviewer

Two pre-existing BFINAL-001-era probes compare the J against an UNRELAXED
fvm reconstruction (no fvMatrix::relax, no relax-source, no constrainHbyA):
- GatePR (run.log 1230-1246 etc.): |J_P|=4.33741818206e-06 vs
  |FD_P|=6.2437647062e-06 -> relL2=1.035, cos=0.307 (was relL2=1.12e-10 in
  BFINAL-002). The FD here is `fvm::div(phi,UmodR)-laplacian+Sp(alpha)` without
  relax (solveDiscreteFlowAdjoint.H ~L1544-1567); it encoded the OLD unrelaxed
  semantics, which is exactly the defect BFINAL-002 proved wrong. After the
  patch the J implements the relaxed (production) tangent, so this probe
  intentionally no longer closes: its FD reference is the outdated one.
- GateOracle full(dPhiJ vs dPhiFD): relL2=0.790, cos=0.356 (was 2.59e-11).
  Same cause (unrelaxed FD reference, dPhiJ = corrected collectedDPhiJ).
- Stage B4 T-cont (continuity row vs div(flux(dU))): 0.0006486 -> 0.0002530
  (the J's P-row changed; this partial-basis probe is not a gate).
- By contrast, the PRODUCTION-semantics closures (the actual gates) all
  improved: StageB5 Gate-A P-total 1.036 -> 1.528e-4; Gate-B candA 0.795 ->
  2.985e-5. The GatePR/GateOracle mismatch is the expected signature of the
  J having switched semantics from unrelaxed to relaxed.

## Acceptance checklist (step criteria)

- EXIT=0: YES (RUN_EXIT=0).
- Log contains StageB5 Gate-A 5-eps sweep: YES (eps 1e-3..1e-5, both sites).
- Gate-B Candidate A/B: YES (candA/candB, 5 eps).
- Transpose dot test full + PU/UP/UU/PP: YES (BlockDot PU/PP/UU/UP/boundaryU
  + Reduced cold-flow operator transpose dot-test + thermal dot-test).
- ExplicitJToracle: YES (2 sites).
- Anchors: YES (phiRb0, momentum anchor, GatePR legacy, GatePRint, GateOracle,
  Stage B4 T1/T2, StageB5 rAU, selfcheck).
- Artifacts exported: YES (9 files in cycle-1/artifacts/, incl. explicitJT.mtx).

## Files

- cycle-1/run.log (raw run log, 1770 lines)
- cycle-1/build3-validation.log (wmake re-verification)
- cycle-1/artifacts/ (9 mtx files)

Pre-state baseline (cycle-1/pre-state/) untouched; BFINAL-002 dirty include and
untracked probe preserved; HEAD unchanged. No forbidden file modified
(git status: only M src/solveDiscreteFlowAdjoint.H,
M src/solveDiscreteFlowAdjointProduction.H).
