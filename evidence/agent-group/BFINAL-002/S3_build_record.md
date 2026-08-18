# BFINAL-002 S3 — rebuild record (Executor's own run)

Date: 2026-08-16 (~17:12 CST). Workspace /home/ys/dsH/TO-ANISOTROPIC @ ca8a772b.

## Rebuild

- Fresh rebuild of MTO_HF from current source (HEAD ca8a772b + S2 guarded
  plumbing in src/solveDiscreteFlowAdjoint.H + new src/stageB5BoundaryRelaxOracle.H).
- Sandbox adaptation (same as S2, documented in stageB5_s2_summary.md): the
  DSH sandbox forbids writing to the normal FOAM_USER_APPBIN
  (/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin), so the link
  target was redirected by exporting
  `FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin` before `wmake`.
- wmake exit code 0; binary rebuilt:
  `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF`
  size 4213040 bytes, mtime 2026-08-16 17:12:13 (+0800).
- Only MTO_HF.o / MTO_HF.C.dep were regenerated (other objects unchanged);
  no production source line touched during this step.

## Run (this S3 execution)

- Case: /home/ys/dsH/b2_case_smoke (writable copy; constant/optProperties has
  `stageB4JacobianProbe true;` and `stageB5BoundaryRelaxOracle true;`,
  snapshot saved as case_optProperties.s3run).
- Invocation (same environment as /home/ys/dsH/run_mtohf.sh, but resolving the
  workspace binary directly because run_mtohf.sh PATH would pick the OLD
  pre-BFINAL-002 binary in FOAM_USER_APPBIN):
  ```
  cd /home/ys/dsH/b2_case_smoke
  export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
  source /opt/openfoam7/etc/bashrc
  export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
  unset FOAM_SIGFPE
  /home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF
  ```
- Raw log: evidence/agent-group/BFINAL-002/stageB5_run_s3.log
- Expected anchors (BFINAL-001 / BFINAL-002 S2): momentum relL2 ~1.65e-4
  cos ~0.999999986; GatePR relL2 ~1e-10; |phiRb0-phi| ~3.49e-11; Gate-A eps
  sweep {1e-3, 3e-4, 1e-4, 3e-5, 1e-5}; Gate-B >=3 eps per candidate.
