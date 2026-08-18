# BFINAL-004 step-2-build record

- Stage: B-final · Mode: DIAGNOSTIC_ONLY · Step: step-2-build (executor)
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace) @ `ca8a772b8339a36e7906ea32a7125061c47aa280` · branch `agent/dsH-stage-b-validation`
- Date: 2026-08-16

## Acceptance criteria

1. **wmake exits 0** — PASS (WMAKE_EXIT=0, two passes, see below).
2. **Build log recorded** — PASS: `evidence/agent-group/BFINAL-004/cycle-1/build_step2.log`
   (full forced recompile + link, 239 lines, no errors, no stageB6-related warnings).
3. **Pre-state baseline captured before any source edit** — PASS: step-1 captured
   `cycle-1/pre-state/` (git_head.txt, git_status_porcelain.txt, git_status_full.txt,
   git_diff_unstaged.patch, file_sha256.txt) at 21:33 before the stageB6 probe edit;
   verified intact and NOT overwritten by this step.

## Build environment (clean PATH)

OpenFOAM-7 `etc/bashrc` → `config.sh/settings:163` **unconditionally overrides**
`FOAM_USER_APPBIN` to `$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin`, so the
export MUST come AFTER `source /opt/openfoam7/etc/bashrc` (documented in
step-1 design note; the plan's literal `export && source` ordering links into the
wrong dir → Permission denied). Used:

```bash
env -i HOME=$HOME PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin bash -c '
  source /opt/openfoam7/etc/bashrc
  export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin
  unset FOAM_SIGFPE
  cd /home/ys/dsH/TO-ANISOTROPIC/src
  wmake'
```

## Results

- Pass 1 (incremental, after step-1 build): WMAKE_EXIT=0, empty log (all objects
  up-to-date w.r.t. current sources; binary already contained stageB6 symbols).
- Pass 2 (forced full recompile: `touch src/MTO_HF.C` — mtime only, zero content
  change, git diff unchanged — then wmake): WMAKE_EXIT=0, `MTO_HF.o` recompiled,
  full link line executed:
  `... -o /home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF`.
- Fresh binary: `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF`
  4691928 bytes, mtime 2026-08-16 22:42:24 +0800; contains 7 `stageB6` string
  symbols (`stageB6RxDesignOracle`, `stageB6_d`, `stageB6_y`, `stageB6_phiHbyA`,
  ...) → probe compiled in.
- Warnings: 63 total, all pre-existing benign (`unused variable` in production
  code, `-Wsign-compare` in pre-existing includes); ZERO warnings reference
  stageB6/stageB5. The `could not open file ompi/mpi/cxx/*.h` lines at the top of
  the log are the wmake makedepends scanner probing missing optional MPI headers
  (pre-existing, non-fatal; link succeeded).

## Source-state verification (no scope violation)

- `git status --short` = `M src/solveDiscreteFlowAdjoint.H`,
  `M src/solveDiscreteFlowAdjointProduction.H` (pre-existing BFINAL-003 patch,
  untouched by this step) + untracked probes/evidence. No new tracked content.
- `git diff --stat` = the same 2 files (276+/53-), identical to the pre-state
  patch record — the forced recompile added no source delta.
- Forbidden files byte-identical to the step-1 pre-state baseline (sha256 match
  for NS.H, sensitivity.H, filter_x.H, filter_chainrule.H, costfunction.H,
  createFrozenHotRegionFields.H, updateMaterialProperties.H, AdjNS_PD.H,
  solveDiscreteFlowAdjointProduction.H, stageB5BoundaryRelaxOracle.H).
- The only `solveDiscreteFlowAdjoint.H` delta vs the reconstructed BFINAL-003
  pre-state (HEAD + pre-state patch, sha verified = baseline) is exactly the
  31-line guarded stageB6 include block at L4055-4085 (switch
  `stageB6RxDesignOracle`, default false, pressureDrop block only).

## Next (step-3-run)

Case `/home/ys/dsH/b2_case_smoke` currently has `stageB4JacobianProbe true;`
`stageB5BoundaryRelaxOracle true;` in `constant/optProperties`; the step-3
executor must add `stageB6RxDesignOracle true;` (switch already compiled into the
binary) and run `MTO_HF` on a writable copy (~8 min single B4 path), then export
the *.mtx artifacts for step-4 recompute/report.
