# BFINAL-006 cycle-1 pre-state notes (Executor S1)

Date: recorded during BFINAL-006 cycle-1 execution.
Task: BFINAL-006 Pressure-Drop Tangent Oracle (DIAGNOSTIC_ONLY).
This file records the immutable pre-state baseline for cycle-1. No source/case files were
modified by this step (S1). See git_status_porcelain.txt and git_diff_stat.txt alongside.

## 1. Provenance (run from /home/ys/dsH/TO-ANISOTROPIC)

- pwd                       : /home/ys/dsH/TO-ANISOTROPIC
- git rev-parse --show-toplevel : /home/ys/dsH/TO-ANISOTROPIC
- git rev-parse HEAD        : ca8a772b8339a36e7906ea32a7125061c47aa280
- git branch --show-current : agent/dsH-stage-b-validation
- Workspace == git root     : YES (matches task expectation)

## 2. git status (pre-existing dirty set, preserved, NOT modified by S1)

Tracked modified (pre-existing, must be preserved):
- M src/sensitivity.H
- M src/solveDiscreteFlowAdjoint.H
- M src/solveDiscreteFlowAdjointProduction.H

Untracked (pre-existing; must be preserved):
- build/
- evidence/agent-group/BFINAL-002/
- evidence/agent-group/BFINAL-003/
- evidence/agent-group/BFINAL-004-FIX/
- evidence/agent-group/BFINAL-004/
- evidence/agent-group/BFINAL-005/
- src/rxPressureRowTranspose.H
- src/stageB5BoundaryRelaxOracle.H
- src/stageB6RxDesignOracle.H
- src/validation_stage_*.log (a set of pre-existing untracked build logs in src/)

Full listing saved to git_status_porcelain.txt (38 lines).
git diff --stat saved to git_diff_stat.txt (only the 3 tracked modified files; no new tracked changes).

## 3. Key artifacts in /home/ys/dsH/b2_case_smoke (existence + sha256)

All five artifacts referenced by the approved plan are PRESENT:

| file                          | sha256                                                             | size (bytes) |
|-------------------------------|--------------------------------------------------------------------|--------------|
| explicitJT.mtx                | 3dfcc9ab41f9344649d31638eac161c0b64e9603b2907b36527f42af88924940  | 467508907    |
| stageB6_rxc_analytic.mtx      | d30e13179e5dcd5065a4b5f509ad24d000efda7499c1c0ce5fc7dd4897bc984c  | 2221113      |
| stageB6_dirs.mtx              | 7572261514f8a8eee2c352772c1b19961c7245ca9539340813f916e0cfd75bb5  | 463561       |
| explicitRhs_pressureDrop.mtx  | 2d7cf846829762427298f9a2bb732e5dcf83ce6d8d93b7a77ec13650ab6bdb5f  | 270505       |
| stageB6_rxc_FD_all_eps.mtx    | ab4efc3a04fdf8311bd3959901c3da42d0d237b53f3c76c9543a0a24888baf84  | 10671064     |

NOTE: sizes are exact bytes from `ls -la` (verified 2026-08 run). stageB6_rxc_analytic.mtx has
no MatrixMarket header — it is a bare list (see §5).

## 4. pRef gauge (system/fvSolution in b2_case_smoke)

- Line 78: pRefCell  5600;
- Line 82: pRefValue 0;
- SIMPLE dict also sets paRefCell/pbRefCell/pcRefCell = 5600.

pRefRow (BFINAL-003 convention): discretePIndex(pRefCell) = 3*N + pRefCell with N = 33600
(134400/4) => pRefRow = 100800 + 5600 = 106400. Cross-check against explicitJT.mtx identity-row
detection is part of S2.

## 5. Artifact matrix headers / layout (recorded for S2 driver)

- explicitJT.mtx : MatrixMarket "matrix coordinate real general", 134400 x 134400, 13558441 nnz.
- explicitRhs_pressureDrop.mtx : MatrixMarket "matrix array real general", 134400 x 1 (cell-major).
- stageB6_rxc_analytic.mtx : bare list (no MM header), first values 0/0/-0 -> 403200 entries
  expected = 3 directions x [anRUd(3N); anRPd(N)] block-major bare list.
- stageB6_dirs.mtx : bare list, first values 0/0/0 -> 100800 entries expected = 3 x N bare list.
- stageB6_rxc_FD_all_eps.mtx : bare list (optional FD-Rx oracle; layout per BFINAL-005 report).

## 6. Python environment (scipy absence re-confirmed)

- venv: evidence/agent-group/BFINAL-002/.venv_s4
- python: 3.11.15 (symlink to /home/ys/.local/bin/python3.11)
- numpy: 2.4.6 PRESENT
- scipy: ABSENT — `import scipy` -> ModuleNotFoundError: No module named 'scipy'

This matches the approved plan's S2 assumption: scipy must be pip-installed into this venv
before the offline tangent solve (planned, NOT executed in S1).

## 7. S1 acceptance checklist

- [x] pre-state files written (this dir)
- [x] HEAD == ca8a772b8339a36e7906ea32a7125061c47aa280
- [x] artifact hashes recorded (5/5)
- [x] scipy-absence re-confirmed and recorded
- [x] zero new tracked/untracked changes beyond the pre-existing dirty set (git status before
      and after S1 identical; only evidence/agent-group/BFINAL-006/ was created)
