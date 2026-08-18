# BFINAL-003 cycle-1 — Step: patch-diagnostic-operator (EXECUTED)

- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace), HEAD
  `ca8a772b8339a36e7906ea32a7125061c47aa280` (unchanged), branch
  `agent/dsH-stage-b-validation`.  Pre-existing dirty state preserved:
  `M src/solveDiscreteFlowAdjoint.H` (BFINAL-002 guarded include) + untracked
  `src/stageB5BoundaryRelaxOracle.H` remain; nothing staged/reset/cleaned.

## Files changed (this step)

- `src/solveDiscreteFlowAdjoint.H` (tracked; the only tracked modification):
  full diff saved as `patch-diagnostic-operator.diff` (354 lines).
- `src/stageB5BoundaryRelaxOracle.H` (untracked probe; copy saved as
  `stageB5BoundaryRelaxOracle.H.patched`).
- `evidence/agent-group/BFINAL-003/cycle-1/{build.log, patch-diagnostic-operator.diff, stageB5BoundaryRelaxOracle.H.patched}`.
- `src/solveDiscreteFlowAdjointProduction.H` intentionally UNTOUCHED (next
  step `patch-production-transpose`), sha256 still
  `135b299f...` == pre-state baseline.

## Changes implemented

1. **alphaRel** (new, after the rAUAdj block):
   `const scalar alphaRel = mesh.relaxEquation("U") ? mesh.equationRelaxationFactor("U") : scalar(1.0);`
   — read from fvSolution `relaxationFactors/equations/U` (0.4 in
   `b2_case_smoke/system/fvSolution`), exactly mirroring `fvMatrix::relax()`
   (OF-7 `solution::relaxEquation`/`equationRelaxationFactor` semantics
   verified in `/opt/openfoam7/src/OpenFOAM/matrices/solution/solution.C` and
   the `fvMatrix::relax()` call site in `fvMatrix.C` L678-685).  Default 1.0
   when U relaxation is absent.  No hard-coded 0.4/0.6/2.5 anywhere in the
   production math (only comments and pre-existing diagnostic vectors).

2. **Forward `applyDiscreteFlowJ` internal faces**: kept the existing
   `deltaPhiFace` (unrelaxed `Sf&rAUdH + kf*(p_nei-p_own)`) **only** for the
   momentum velocityJump (conv) term (J_UU, byte-for-byte unchanged).  New
   `deltaPhiFacePU = Sf&(alphaRel*rAUdH + (1-alphaRel)*interpDU)
   + mobF*(p_nei-p_own)*dcfF*mafF` (P-P part identical to the old expression)
   routed **only** to the P-row `+=`/`-=` and to `collectedDPhiJ` (so every
   J-dphi diagnostic reflects the corrected operator).

3. **Forward boundary loop**: `dphi_b = 0` when
   `!U.boundaryField()[patchi].assignable()` (constrainHbyA pinning on the 7
   fixedValue/noSlip U patches), else
   `Sf_b&(alphaRel*rAUb*dHb + (1-alphaRel)*dU_cell)` (extrapolated outlet).
   Direct U-P (pressure-gradient) boundary term unchanged.

4. **Transpose `applyDiscreteFlowJT` internal faces**: kf (P-P) path keeps
   `phiAdjoint = lambdaPdiff + convTerm` unchanged; hA path now uses
   `phiAdjointRel = convTerm + alphaRel*lambdaPdiff`; new direct
   `(1-alphaRel)*lambdaPdiff*(wf/1-wf)*Sf` contributions to
   `output[U(own/nei)]`.  U-U offdiag transpose, direct pressure-gradient
   transpose, dev2 transpose unchanged.

5. **Transpose boundary loop**: hB now `lambdaPc*(alphaRel*rAUAdj)*patchSf`
   and **assignable-only**; new direct `(1-alphaRel)*lambdaPc*patchSf` term on
   assignable faces only; direct U-P (P row <- U) term unchanged.

6. **Explicit CSR oracle (J^T, physical)**: internal direct
   `(1-alphaRel)` COO entries (row U, col P); hB accumulation assignable-only
   and alphaRel-scaled; internal hA-path P-part coeffs `coP_o2/coP_n2`
   alphaRel-scaled (`convCo` untouched); boundary hB alphaRel-scaled +
   assignable-only + direct `(1-alphaRel)` boundary term.  U-U/P-P/direct-U-P
   entries unchanged.  (ExplicitJToracle check at the next run will verify
   explicit == matrix-free at ~1e-16.)

7. **`stageB5BoundaryRelaxOracle.H`**: Candidate-A boundary extraction and its
   self-consistency check updated to the corrected semantics (pinning +
   alphaRel + direct dU term) so Gate-A P rows and Gate-B "Candidate A"
   measure the patched operator.  Candidate-B reconstruction, Gate-A FD,
   anchors, and exports unchanged.

## Transpose consistency (derivation, recorded for the reviewer)

Forward P-row = `dphi_PU` (U part relaxed, P-P part unchanged) while the conv
(velocityJump) term keeps the old unrelaxed `deltaPhiFace`; the adjoint scalar
therefore splits into `convTerm + alphaRel*lambdaPdiff` (hA) + direct
`(1-alphaRel)*lambdaPdiff*Sf` (U rows), with kf keeping the full
`lambdaPdiff + convTerm` — verified term-by-term; UU/UP/PP block dot tests are
structurally unchanged (alphaRel terms vanish when the respective block's
delta/lambda are zero), so the runtime BlockDot and full dot tests are the
checks that the patch is an exact transpose.

## Build (acceptance)

- Command: `cd src && (PATH=clean; source /opt/openfoam7/etc/bashrc; export
  FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin; unset FOAM_SIGFPE;
  wmake)` — note: `FOAM_USER_APPBIN` must be exported AFTER sourcing bashrc
  (bashrc L163 unconditionally resets it to the default user path, which the
  sandbox cannot write).
- Result: `WMK_EXIT=0`; binary rebuilt
  `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF` (Aug 16 19:48, 4233232 B).
- Compile warnings are all pre-existing (verified against
  `evidence/agent-group/BFINAL-002/stageB5_build.log`: unused-variable
  `nei`/`rAUo`/`Vinvo`/`nBS`/`numT1Conv`/`stageB2FDStepsMulti` etc.).

## Acceptance checklist (patch-diagnostic-operator)

- [x] Builds clean (WMK_EXIT=0).
- [x] No literal 0.4/0.6/2.5 in production math (grep: only comments /
      pre-existing diagnostic constants).
- [x] Momentum velocityJump term byte-for-byte unchanged
      (`output[U_downwind] += deltaPhiFace*velocityJump.component(cmpt)`
      is a context line in `git diff`).
- [x] P-P kf/laplacian path byte-for-byte unchanged (forward P-P expression
      identical in `deltaPhiFacePU`; transpose kf path untouched).
- [x] `collectedDPhiJ[facei] = deltaPhiFacePU` == the corrected P-row dphi
      (same expression routed to `output[P(own)] +=` / `output[P(nei)] -=`).
- [x] Only files_allowed touched; production file byte-identical to baseline;
      pre-state baseline dir not modified; BFINAL-002 dirty include preserved.

## Handoff to next step (build-and-run-validation)

Run `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF` on `/home/ys/dsH/b2_case_smoke`
(writable copy) with `stageB4JacobianProbe true;` and
`stageB5BoundaryRelaxOracle true;` (~8 min).  Expected: Gate-A P-total
relL2 < 1e-3/cos > 0.999 (was 1.036), Gate-B Candidate A == Candidate B at
roundoff, momentum stays 1.653e-4, BlockDot/full dot tests green, and the
explicit/matrix-free oracle check at ~1e-16.  Note: the pre-existing GatePR
diagnostic in the Stage-B4 block compares J's P-row against an UNRELAXED
phi-FD basis; after the patch it will print an O(1) mismatch BY DESIGN (the
unrelaxed basis is not the production semantics — BFINAL-002 §1.1/2.1); the
StageB5 Gate-A P-row closure uses the production relaxed FD and is the
authoritative G2 evidence.
