# BFINAL-003 step: patch-production-transpose — execution record

Status: EXECUTED (build clean)
Date: 2026-08-16 (session clock)
Workspace: /home/ys/dsH/TO-ANISOTROPIC, git root == workspace, HEAD ca8a772b8339a36e7906ea32a7125061c47aa280
branch agent/dsH-stage-b-validation (unchanged by this step).

## What was done

Patched the production matrix-free transpose `applyProdFlowJT` in
`src/solveDiscreteFlowAdjointProduction.H` with the identical J^T correction
already applied to the diagnostic `applyDiscreteFlowJT` (BFINAL-003
patch-diagnostic-operator step), so the production pressure adjoint uses the
Candidate-B (relaxed-SIMPLE) P<-U operator.

Three targeted edits, nothing else:

1. **alphaRel read** (after `prodRAU`):
   `const scalar alphaRel = mesh.relaxEquation("U") ? mesh.equationRelaxationFactor("U") : scalar(1.0);`
   Same semantics as the diagnostic file (`fvMatrix::relax()` convention; for the
   b2_case_smoke case fvSolution `equations/U = 0.4`, so alphaRel = 0.4).
   No hard-coded 0.4/0.6/2.5 anywhere in production math (only explanatory
   comments). The production file is included in a mutually exclusive branch
   from the diagnostic file (AdjNS_HT.H / AdjNS_PD.H), so its own alphaRel
   definition is required and cannot collide.

2. **Internal faces**: `phiAdjoint` split into `convTerm = lambdaDown & velocityJump`
   and `lambdaPdiff = lambdaPown - lambdaPnei`.
   - kf (P-P) path: unchanged, still uses `phiAdjoint = lambdaPdiff + convTerm`.
   - hA path (P-U): `phiAdjointRel = convTerm + alphaRel*lambdaPdiff` replaces
     `phiAdjoint` in `hOwn/hNei = phiAdjointRel*(weight|1-weight)*prodRAU*prodSf`.
   - Direct relax-source transpose: `output[prodUIndex] += (1-alphaRel)*lambdaPdiff*(weight*prodSf)`
     (own weighted by `weight`, nei by `1-weight`).

3. **Boundary faces**: `uAssignable = U.boundaryField()[patchi].assignable()`;
   hB is assignable-only and alphaRel-scaled (`alphaRel*prodRAU`), plus a direct
   `(1-alphaRel)*lambdaPc*patchSf` term on assignable faces. Non-assignable
   (fixedValue/noSlip) patches -> dphi_b/dU = 0 (constrainHbyA pinning; the 7
   fixedValue-U faces). The direct U-P transpose (`!fixesValue()` lambdaCell
   term) is untouched.

## Unchanged blocks (verified by diff)

- prodDiag (U-U diagonal), prodLower/prodUpper (U-U offdiag A^T)
- direct momentum pressure-gradient transpose (U-P)
- kf (P-P) path
- deltaH^T Stage-2 (`-upper/V` hA routing) and boundary H() diagonal transpose
- applyProdFrozenDeviatoricJT (dev2), prodMomentumScale/prodAreaScale,
  preconditioners, FGMRES solve layer

## Mirror check vs diagnostic applyDiscreteFlowJT

Identical formula, only naming differs: weight == discreteWeights,
prodRAU == rAUAdj, prodSf == discreteSf, prodUIndex/prodPIndex ==
discreteUIndex/discretePIndex, convTerm/lambdaPdiff identical definitions,
boundary guard identical (`assignable()`), alphaRel identical read.

## Build

- Attempt 1 (env facts order: export FOAM_USER_APPBIN BEFORE sourcing
  /opt/openfoam7/etc/bashrc): MTO_HF.C compiled cleanly (new MTO_HF.o, no
  compile errors), but the link failed with "Permission denied" because
  etc/bashrc resets FOAM_USER_APPBIN to $WM_PROJECT_USER_DIR/... . Recorded in
  cycle-1/build.log (the single 'error:' line is the linker output-file
  permission error, not a compile error).
- Attempt 2 (source bashrc FIRST, then export FOAM_USER_APPBIN=build/bin):
  WMAKE_EXIT=0, binary relinked at
  /home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF (20:17, 4237328 bytes).
  cycle-1/build2.log: 0 errors.

## Artifacts

- cycle-1/patch-production-transpose.diff  (git diff of the production file)
- cycle-1/build.log / cycle-1/build2.log
- cycle-1/file_sha256_production_patched.txt:
  sha256(src/solveDiscreteFlowAdjointProduction.H) =
  8c4901bf8da4edbc7cc1efe73a7583963e4688bba93c8c58539c4c39cae1ad91

## Acceptance (step criteria)

- Builds clean: YES (WMAKE_EXIT=0, 0 errors).
- Mirrors the validated diagnostic transpose formula exactly (only prodX naming
  differs): YES (verified line-by-line).
- No change to any U-U/U-P/P-P block in the production file: YES (diff shows
  only the J_PU transpose hA/direct/boundary terms plus the alphaRel read).

Pre-existing BFINAL-002 dirty include and untracked probe preserved; pre-state
baseline untouched; HEAD unchanged. Next step: build-and-run-validation (single
B4 path on /home/ys/dsH/b2_case_smoke).
