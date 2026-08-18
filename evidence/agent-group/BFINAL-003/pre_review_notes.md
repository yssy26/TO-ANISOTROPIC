# BFINAL-003 Pre-Reviewer notes (independent invocation, deepseek-v4-pro)

## Workspace verification
- pwd = /home/ys/dsH/TO-ANISOTROPIC; git root == workspace; HEAD ca8a772b8339a36e7906ea32a7125061c47aa280 (== task baseline); branch agent/dsH-stage-b-validation. No mismatch.
- Pre-existing dirty state confirmed: `M src/solveDiscreteFlowAdjoint.H` (20-line guarded stageB5 include, EOF) + untracked `src/stageB5BoundaryRelaxOracle.H`. `git diff --stat` = 1 file / 20 insertions. BFINAL-002 tooling preserved.

## Source verification (independent read)
- `applyDiscreteFlowJ` (solveDiscreteFlowAdjoint.H L414-711): internal P-U dphi = `Sf&(wf*rAUAdj_o*dH_o + (1-wf)*rAUAdj_n*dH_n)` (UNRELAXED), while the P-P term in the SAME `deltaPhiFace` already uses `mobF` = interpolated `primalPressureMobility` (RELAXED). Mixed convention = the defect. Boundary dphi (L689-697) extrapolates `Sf_b&(rAUAdj*dH)` at ALL patches (no constrainHbyA pinning).
- `applyDiscreteFlowJT` (L713-867): hA path `phiAdjoint*wf*rAUAdj*Sf` (unrelaxed), boundary hB `lambdaPc*rAUAdj*patchSf` (unrelaxed, all patches). Mirrors the same defect.
- CSR oracle (L2270-2485): `coP_o2/coP_n2 = coefP*rAUo2` (unrelaxed mobility), boundary hB `rAUAdj*pSf` all patches. Same defect.
- `applyProdFlowJT` (solveDiscreteFlowAdjointProduction.H L312-458): identical structure with prodX naming; `prodMomentum` built WITHOUT `.relax()` so `prodRAU` is unrelaxed — same defect in the production adjoint transpose.
- No other J_PU/dphi representation outside these two files (grep of rAUAdj/prodRAU/deltaPhiFace/collectedDPhiJ/dphi across src/*.H).

## Independent algebra check of the proposed tangent
- OpenFOAM-7 `fvMatrix::relax`: D/=alpha, S+=(D-D0)*psi (invariant residual). Derivation: HbyA_rel = alpha*rAU_u*H_u + (1-alpha)*psi, hence dphi/dU = flux(alpha*rAU_u*dH_u + (1-alpha)*dU). Matches BFINAL-002 Candidate B exactly (probe computes rAUrel*dHrel which equals the same expression). Naive rAUAdj*=0.4 rejected (misses relax-source; BFINAL-002 X3 already disproved).
- Relax factor: case fvSolution `equations.U = 0.4` (confirmed), rAUAdj/rAUrel = 2.5 = 1/0.4 (confirmed in s3 log). `solution::relaxEquation`/`equationRelaxationFactor` are methods of `fvSolution` (accessed via `mesh.solutionDict()`), NOT fvMesh — plan's "mesh.equationRelaxationFactor" shorthand needs `mesh.solutionDict()`.

## Gate assessment
- SCOPE: pass. J_PU-only within B-final R_w closure; forbidden modules protected; no SIMPLE-transpose.
- PHYSICS: pass. Tangent = validated production relaxed-SIMPLE mapping + constrainHbyA pinning.
- ROOT_CAUSE: pass. BFINAL-002 discriminated decisively; hypothesis_id BFINAL-003-PU-RELAXATION-PATCH is new (patch counts empty).
- CHANGE_SCOPE: pass. VelocityJump (momentum), kf (P-P), direct grad-p (U-P) byte-for-byte unchanged → J_UU/J_UP/J_PP preserved; P-row rerouted to corrected flux.
- VALIDATION: pass with clarifications (G5 anchor coverage + exact direct-(1-alpha) COO pattern + solutionDict API).

## Decision
APPROVED (with 3 non-blocking required clarifications).
