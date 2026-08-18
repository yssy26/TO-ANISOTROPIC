# BFINAL-002 S2 — executor completion summary

Step: S2-add-oracle-probe (DIAGNOSTIC_ONLY). Status: EXECUTED.
HEAD ca8a772b, branch agent/dsH-stage-b-validation, git root
/home/ys/dsH/TO-ANISOTROPIC == session workspace.

## Deliverables

1. `src/stageB5BoundaryRelaxOracle.H` (NEW, ~780 lines) — the guarded
   Gate-A / Gate-B diagnostic oracle (design: stageB5_probe_design.md).
2. `src/solveDiscreteFlowAdjoint.H` (+20 lines, guarded invocation only;
   the ONLY tracked modification; `git diff --stat` = 20 insertions).
3. Evidence under `evidence/agent-group/BFINAL-002/`:
   - stageB5_probe_design.md (design + verified-behavior notes)
   - recompute_stageB5.py (independent recomputation tool for the reviewer)
   - crosscheck_bf2_runoff.py (BFINAL-001 crosscheck pointed at the run-off
     exports)
   - stageB5_build.log / stageB5_build2.log (wmake; link step sandbox-adapted)
   - stageB5_run_off.log (switch OFF), stageB5_run_on.log and
     stageB5_run_on2.log (switch ON; on2 = fixed probe)
   - artifacts/ (all 8 stageB5_*.mtx raw exports)
   - case_optProperties.{baseline,run1,run2} (case config records)

## Acceptance verification

- Probe compiles (wmake compile OK; link done manually into
  /home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF because the DSH sandbox
  forbids writes to /home/ys/OpenFOAM/.../bin — documented adaptation).
- Switch OFF: production path byte-identical — zero StageB5 lines; every
  BFINAL-001 anchor reproduced EXACTLY (NS UEqn A() avg 158965329.116,
  |phiRb0-phi| 3.49247878283e-11, GatePR 1.12295958888e-10, GatePRint,
  44 correctors); re-running crosscheck_bf2.py on the run-off exports gives
  momentum 1.653123e-4 / cos 0.999999986, P-row 0.824, J-matrix-free
  1.239622e-15 — identical to crosscheck_bf2.out.
- Switch ON: emits per-block Gate-A metrics (momentum / P-internal /
  P-boundary / P-total: |Jv|, |FD|, |err|, relL2, cos, max-abs-diff +
  location) for 5 eps {1e-3,3e-4,1e-4,3e-5,1e-5} with a stable FD plateau;
  Gate-B per-candidate dphi metrics (|analytic|, |FD|, norm ratio, relL2,
  cos, internal/boundary-face rel err) for 5 eps; relax factor (0.4) and
  rAUAdj vs relaxed-rAU avg/min/max/norm-ratio (2.50000000002); exports 8
  stageB5_*.mtx files.  Independent recomputation (recompute_stageB5.py)
  from the raw artifacts reproduces every log number.
- No production math touched; no patch; no SIMPLE-transpose; no MMA; no R_x.

## Key numbers (for the S4 report; interpretation is the S4 step's job)

Gate A (production-consistent FD, stable across eps):
  momentum relL2=1.653e-4 cos=0.999999986
  P-internal relL2=0.858 cos=0.612 | P-boundary relL2=1.145 cos=0.153
  (max @ P cell 7761 patch sideWalls) | P-total relL2=1.036 cos=0.370
Gate B:
  Candidate A (J, unrelaxed rAUAdj): normRatio=0.591 relL2=0.795
  (int 0.790, bnd 1.388) cos=0.607
  Candidate B (production relaxed mobility x relaxed H tangent): normRatio
  =1.000 relL2=1.06e-11 (int 1.06e-11, bnd 1.33e-11) cos=1.000
  -> Candidate B reproduces the true production dphi/dU to roundoff; the
  O(1) mismatch is on Candidate A (the current J formulation).  The 2.5x
  rAU ratio does NOT appear 1:1 in dphi: the relax-source term in the
  relaxed H() partially compensates (normRatio 0.59, not 0.4 or 2.5).

## Hand-off note for S3 (run step)

- The case `/home/ys/dsH/b2_case_smoke/constant/optProperties` currently has
  `stageB4JacobianProbe true;` and `stageB5BoundaryRelaxOracle true;` (see
  case_optProperties.run2), and `discreteExplicitMatrixFile` redirected to
  `/home/ys/dsH/b2_case_smoke/explicitJT.mtx` (the sandbox forbids writing to
  /home/ys/b2_case_smoke).
- `bash /home/ys/dsH/run_mtohf.sh <case>` resolves `MTO_HF` from PATH
  (FOAM_USER_APPBIN = /home/ys/OpenFOAM/ys-7/.../bin), which holds the OLD
  pre-BFINAL-002 binary — the sandbox prevents installing the new binary
  there.  S3 must invoke the workspace binary directly with the same
  environment as run_mtohf.sh:
  `export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin";
   source /opt/openfoam7/etc/bashrc;
   export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH";
   unset FOAM_SIGFPE; cd /home/ys/dsH/b2_case_smoke;
   /home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF`
  (This exact invocation produced stageB5_run_on2.log.)
- `build/` (untracked) holds the workspace-linked binary + link scripts; it is
  a build byproduct, not production code.
