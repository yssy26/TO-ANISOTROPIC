# BFINAL-002 S2 — stageB5BoundaryRelaxOracle probe design (record)

Stage: B-final (DIAGNOSTIC_ONLY, BFINAL-002, step S2-add-oracle-probe)
Author: Executor (flash).
HEAD: ca8a772b (branch agent/dsH-stage-b-validation), workspace /home/ys/dsH/TO-ANISOTROPIC.
Purpose: record the design and implementation of the Gate-A / Gate-B guarded
diagnostic probe added in S2, plus build/run environment notes, so the
Post-Reviewer can re-verify without re-deriving the code.

---

## 1. Files changed (only files_allowed)

1. `src/stageB5BoundaryRelaxOracle.H` (NEW, 774 lines) — the whole probe.
2. `src/solveDiscreteFlowAdjoint.H` (+20 lines at the end) — guarded invocation:

```cpp
if (optProperties.lookupOrDefault<Switch>("stageB5BoundaryRelaxOracle", false))
{
    static bool stageB5BoundaryRelaxOracleDone = false;
    if (!stageB5BoundaryRelaxOracleDone)
    {
        stageB5BoundaryRelaxOracleDone = true;
        #include "stageB5BoundaryRelaxOracle.H"
    }
}
```

`git diff --stat` at HEAD shows only `solveDiscreteFlowAdjoint.H | 20 +`; the
new header is untracked. No production math line was touched (diff above).

No-op when disabled: with the switch absent/false the guard is skipped
(one dictionary lookup only), so the production path is byte-identical
(no field created/modified, no file written, no residual/value change).
Enabled: runs once (static guard) on the FIRST inclusion of
solveDiscreteFlowAdjoint.H (the thermalCoupling adjoint block in
AdjNS_HT.H), i.e. after the converged NS.H state exists and after
applyDiscreteFlowJ / rAUAdj / discretePrimalMomentum are in scope.
The probe never writes U/p/phi/alpha/k/omega/nutFrozen/nuEffFrozen/
primalPressureMobility/rAUAdj; it only reads them and creates local
temporaries + writes stageB5_*.mtx exports to the case CWD.

## 2. Gate A (Q-A) construction

- Deterministic direction `dUdir` = sin(0.271*(celli+1)) pattern, constructed
  EXACTLY like NS.H L466-485 (so stageB5_dUdir.mtx should match the BFINAL-001
  stageB2_dUdirPhys.mtx).
- v = [dU; dp=0] (4N state vector, P part zero).
- `Jv = applyDiscreteFlowJ(v)` (existing matrix-free reduced Jacobian),
  `collectDPhiJ=true` during the call so `collectedDPhiJ` holds J's
  internal-face dphi tangent (Candidate A internal).
- FD reference = central difference of the FULL production reduced residual
  R = (R_U, R_P) with phi rebuilt through the exact NS.H SIMPLE sequence:
  fresh UEqn (div(phi,U)-laplacian(nuEffFrozen,U)+Sp(alpha,U) == -grad(p)+
  fvOptions) -> dev2 source (ransFlowModel) -> UEqn.relax() ->
  fvOptions.constrain -> rAU=1/A -> constrainHbyA(rAU*H,Umod,p) ->
  fvc::flux(HbyA) -> adjustPhi -> rAtU(=rAU, consistent=false) ->
  constrainPressure -> pEqn=laplacian(rAtU,p) == div(phiHbyA) ->
  setReference(pRefCell) -> phi = phiHbyA - pEqn.flux(); boundary flux
  REBUILT from the rebuilt phi (NOT pinned).  R_U = -(UEqn.residual())
  (NS.H sign convention); R_P = div(rebuilt phi) with rebuilt boundary.
- Per eps in {1e-3,3e-4,1e-4,3e-5,1e-5}: report momentum rows (3N),
  P-internal cells, P-boundary cells (cells with >=1 boundary face),
  P-total: |Jv|L2, |FD|L2, |err|L2, relL2 (=|err|/|FD|, same convention as
  crosscheck_bf2.py), cosine, max-abs-diff + location (cell/patch).
- Anchor reproduction: frozen-phi momentum FD (exact NS.H L501-526 recipe,
  hU=1e-3*max|U|, no dev2, no RHS) vs Jv momentum -> expect
  relL2~1.65e-4 cos~0.999999986 (the BFINAL-001 anchor).  The existing
  GatePR block (1e-10) and RxProbe |phiRb0-phi|~3.5e-11 anchors are emitted
  by the pre-existing stageB4JacobianProbe code in the same run.

## 3. Gate B (Q-B) construction

- p / design (alpha) / k / omega / nutFrozen / nuEffFrozen fixed; dU = dUdir.
- True FD dphi = [phi_rebuilt(U+eps dU) - phi_rebuilt(U-eps dU)]/(2eps) per
  face (internal then boundary faces), production rebuild.
- Candidate A (current production/J semantics, unrelaxed rAUAdj):
  internal faces = collectedDPhiJ; boundary faces = Sf_b & (rAUAdj*dH_J)
  with dH_J recomputed as the frozen unrelaxed H-tangent (copy of
  applyDiscreteFlowJ L502-582); self-consistency of the recomputed dH is
  verified against outJ[P] (printed as StageB5 selfcheck).
- Candidate B (diagnostic-only relaxed primal-equivalent mobility):
  dphiB = Sf & (wf*rAU_rel_o*dH_rel_o + (1-wf)*rAU_rel_n*dH_rel_n) internal,
  boundary: 0 at !assignable U patches (constrainHbyA pinning), extrapolated
  Sf_b & (rAU_rel*dH_rel) at the assignable U patch (outlet).  rAU_rel =
  production relaxed rAU (1/A of the relaxed UEqn, == primalPressureMobility);
  dH_rel = FD of the production RELAXED H (reassembled production UEqn at
  U±eps dU; affine in U at frozen p/coeffs since dA/dU=0, so the H-FD is the
  exact H-tangent including the relax source and dev2 source).  This is the
  audit's §8 requirement (relaxed H incl. relax source + relaxed rAU +
  constrainHbyA), not a naive "0.4*rAUAdj x unrelaxed dH" shortcut.
- Per candidate per eps in {1e-3,1e-4,1e-5} (subset of Gate A eps): |cand|L2,
  |FD|L2, norm ratio, relL2, cos, internal-face rel err, boundary-face rel err
  (per-face vectors over internal faces then boundary faces in patch order).
- rAU stats: relaxFactorU from fvSolution relaxationFactors/equations/U,
  rAUAdj (unrelaxed) avg/min/max, relaxed rAU (primalPressureMobility)
  avg/min/max, avg ratio, norm ratio.

## 4. Exports written to the case CWD (stageB5_*.mtx)

- stageB5_gateA_Jv.mtx           (4N: J*v)
- stageB5_dUdir.mtx              (3N: direction)
- stageB5_gateA_FD_all_eps.mtx   (5 blocks of 4N; NOTE: each 4N block is
  written PER CELL interleaved as Ux Uy Uz P, i.e. index 4*celli+k for the k-th
  U component and 4*celli+3 for P — NOT [U(3N), P(N)]; recompute_stageB5.py
  unpacks this layout)
- stageB5_gateB_dphi_A.mtx       (nInt+nBnd: Candidate A dphi, eps-independent)
- stageB5_gateB_dphi_B_all_eps.mtx (5 blocks of nInt+nBnd: Candidate B dphi)
- stageB5_gateB_dphi_FD_all_eps.mtx (5 blocks of nInt+nBnd: FD dphi)
- stageB5_face_type.mtx          (nInt+nBnd: 0 internal, 1..8 patch index+1)
- stageB5_cell_type.mtx          (N: 0 internal cell, 1 boundary cell)

Copies are archived under evidence/agent-group/BFINAL-002/artifacts/ (the run
writes them into the case CWD; the archive is the stable record).

## 4a. Verified behavior notes (S2 runs)

- The probe runs once per #include site of solveDiscreteFlowAdjoint.H: the
  block-scope `static bool` is per-inclusion-site, so the probe executes twice
  per run (AdjNS_HT thermalCoupling block and AdjNS_PD pressureDrop block).
  Both executions produce IDENTICAL output (verified by diff of the two
  StageB5 sections in stageB5_run_on.log), because the flow state
  (U/p/phi/alpha) is unchanged between the two adjoint blocks.  Harmless,
  deterministic, read-only; kept as-is rather than adding a marker-file hack.
- S2 acceptance results (stageB5_run_on2.log, switch ON):
  - Anchor (frozen-phi momentum FD): relL2=1.653123e-4 cos=0.999999986
    (reproduces the BFINAL-001 crosscheck_bf2 anchor).
  - Gate A eps=1e-3 (stable across all 5 eps): momentum relL2=1.6531e-4 /
    cos=0.999999986; P-internal relL2=0.8584 / cos=0.6124; P-boundary
    relL2=1.1454 / cos=0.1528 (max @ P cell 7761 patch sideWalls);
    P-total relL2=1.0355 / cos=0.3704.
  - Gate B eps=1e-3 (>=3 eps plateau): Candidate A |A|/|FD|=0.5908,
    relL2=0.7946 (int 0.7904, bnd 1.3877), cos=0.6073; Candidate B
    normRatio=1.0000, relL2=1.06e-11 (int 1.06e-11, bnd 1.33e-11), cos=1.
  - rAU: relaxFactorU(fvSolution)=0.4; rAUAdj avg 7.0584e-7 [min 1.0e-8,
    max 5.13e-6]; relaxed rAU avg 2.8234e-7 [min 4.0e-9, max 2.05e-6];
    avgRatio=2.50000000002, normRatio=2.49999999999.
  - Independent recomputation with evidence/agent-group/BFINAL-002/
    recompute_stageB5.py from the raw artifacts reproduces every log number.
  - All BFINAL-001 anchors reproduced in the same run: |phiRb0-phi|
    =3.49247878283e-11, GatePR=1.12295958888e-10, NS UEqn A() avg
    =158965329.116, Frozen-RANS converged after 44 correctors.
- Switch OFF (stageB5_run_off.log): zero StageB5 lines; every anchor
  IDENTICAL to rerun.log (NS A() avg, rAtU, |phiRb0-phi|, GatePR, GatePRint,
  convergence) => production path byte-identical.

## 5. Build / run environment notes (sandbox adaptations, all documented)

- The DSH file sandbox (workspace-write) only permits writes under
  /home/ys/dsH.  wmake's link target /home/ys/OpenFOAM/ys-7/.../bin/MTO_HF is
  outside the workspace, so the link step was executed manually with the same
  g++ command (from the wmake log) and `-o` redirected to
  /home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF, run with the OpenFOAM bashrc
  environment sourced (needed for the linker to resolve dependent libs).
  Object files compiled fine inside src/Make/... (workspace).
- The run writes explicitJT.mtx (455 MB) to the configured
  discreteExplicitMatrixFile; the BFINAL-001 reference location
  /home/ys/b2_case_smoke/ is outside the workspace (read-only here), so the
  writable case copy's discreteExplicitMatrixFile was redirected to
  /home/ys/dsH/b2_case_smoke/explicitJT.mtx (same content, new location).
  discreteExplicitSolutionFile (READ) stays at /home/ys/b2_case_smoke/explicitSol.
- Runs are executed with the same environment as run_mtohf.sh (clean PATH,
  source /opt/openfoam7/etc/bashrc, unset FOAM_SIGFPE, cd case) but invoking
  /home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF directly instead of via
  `exec MTO_HF` (PATH would resolve the old binary in FOAM_USER_APPBIN).
- Case state: b2_case_smoke optProperties baseline backed up to
  case_optProperties.baseline; run-1 copy at case_optProperties.run1.
  stageB4JacobianProbe remains true (as in the baseline).

## 6. Verification plan for this step

- Run 1 (switch OFF): anchors reproduced (GatePR ~1e-10, |phiRb0-phi|
  ~3.5e-11, NS UEqn A() avg ~1.59e8, NS rAtU avg ~2.82e-7, momentum anchor
  ~1.65e-4) and NO "StageB5" lines in the log => production path unchanged.
- Run 2 (switch ON): probe emits Gate-A per-block metrics for 5 eps, Gate-B
  per-candidate metrics for >=3 eps, rAU stats; exports appear; anchors still
  reproduced (probe does not disturb them).
