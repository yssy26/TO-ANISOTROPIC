# BFINAL-006 cycle-2 pre-state notes (Executor S1)

Date: recorded during BFINAL-006 cycle-2 execution.
Task: BFINAL-006 Pressure-Drop Tangent Oracle (DIAGNOSTIC_ONLY) — cycle-2 is the
plan correction after cycle-1 FAIL_PLAN_OR_HYPOTHESIS (wrong pRefRow=106400 pin +
too-strict identity test nrow==1). hypothesis_id stays BFINAL-006-TANGENT-CLOSURE.
This file records the immutable pre-state baseline for cycle-2. No source/case
files were modified by this step (S1). See git_status_porcelain.txt,
git_diff_stat.txt, git_head.txt, git_branch.txt, pwd.txt, git_toplevel.txt
alongside.

## 1. Provenance (run from /home/ys/dsH/TO-ANISOTROPIC)

- pwd                        : /home/ys/dsH/TO-ANISOTROPIC
- git rev-parse --show-toplevel : /home/ys/dsH/TO-ANISOTROPIC
- git rev-parse HEAD         : ca8a772b8339a36e7906ea32a7125061c47aa280
- git branch --show-current  : agent/dsH-stage-b-validation
- Workspace == git root      : YES (matches task expectation)

## 2. git status (pre-existing dirty set, preserved, NOT modified by S1)

Tracked modified (pre-existing, locked, must be preserved):
- M src/sensitivity.H                        (BFINAL-005 locked head)
- M src/solveDiscreteFlowAdjoint.H           (BFINAL-003 locked head)
- M src/solveDiscreteFlowAdjointProduction.H (BFINAL-003 locked head)

Untracked (pre-existing; preserved):
- build/
- evidence/agent-group/BFINAL-002/ ... BFINAL-006/
- src/rxPressureRowTranspose.H  (BFINAL-005 helper)
- src/stageB5BoundaryRelaxOracle.H (BFINAL-002 probe)
- src/stageB6RxDesignOracle.H   (BFINAL-004 probe)
- src/validation_stage_*.log (pre-existing build logs in src/)

Full listing saved to git_status_porcelain.txt; git diff --stat HEAD saved to
git_diff_stat.txt (only the 3 locked tracked files).

## 3. Key artifacts in /home/ys/dsH/b2_case_smoke (existence + sha256, cycle-2)

| file                          | sha256                                                             | size (bytes) |
|-------------------------------|--------------------------------------------------------------------|--------------|
| explicitJT.mtx                | 3dfcc9ab41f9344649d31638eac161c0b64e9603b2907b36527f42af88924940  | 467508907    |
| stageB6_rxc_analytic.mtx      | d30e13179e5dcd5065a4b5f509ad24d000efda7499c1c0ce5fc7dd4897bc984c  | 2221113      |
| stageB6_dirs.mtx              | 7572261514f8a8eee2c352772c1b19961c7245ca9539340813f916e0cfd75bb5  | 463561       |
| explicitRhs_pressureDrop.mtx  | 2d7cf846829762427298f9a2bb732e5dcf83ce6d8d93b7a77ec13650ab6bdb5f  | 270505       |
| stageB6_rxc_FD_all_eps.mtx    | ab4efc3a04fdf8311bd3959901c3da42d0d237b53f3c76c9543a0a24888baf84  | 10671064     |

All five hashes are byte-identical to the cycle-1 recorded values (unchanged
artifacts; no re-export happened between cycles).

## 4. Matrix layout facts (verified this round)

- explicitJT.mtx : MatrixMarket "matrix coordinate real general", 134400 x 134400,
  nnz-header 13558441 coordinate lines (COO with duplicates; unique ~6209576).
- N = discreteNCells = 134400/4 = 33600.
- discretePIndex(celli) = discreteNVelocity + celli = 3*33600 + celli.
  - discretePIndex(0)    = 100800
  - discretePIndex(5600) = 106400
- fvSolution (b2_case_smoke/system/fvSolution L76-83): pRefCell 5600; pRefValue 0;
  paRefCell/pbRefCell/pcRefCell 5600.
- 0/p boundaryField: outlet { type fixedValue; value uniform 0; }; all other 8
  patches zeroGradient -> p.needReference()==false (OpenFOAM-7
  GeometricField::needReference: any patch fixesValue() -> needRef=false).

## 5. Python environment (verified working)

- venv: evidence/agent-group/BFINAL-002/.venv_s4
- numpy 2.4.6, scipy 1.17.1, scipy.sparse.linalg.splu import OK.

## 6. S1 gate outcome (filled by the S1 script run)

See cycle-2/S1_pressure_reference_gate/ directory:
- s1_gate_scan.log        (full programmatic scan output)
- s1_raw_scan_evidence.md (extracted raw evidence)
- S1_GATE_PASSED marker file when all 5 checks pass.
