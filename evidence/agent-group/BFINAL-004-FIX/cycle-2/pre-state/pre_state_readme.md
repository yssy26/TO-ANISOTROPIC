# BFINAL-004-FIX cycle-2 pre-state (step1-pre-state)

Captured 2025-08-17 02:18 (workspace /home/ys/dsH/TO-ANISOTROPIC).

## Files in this directory

- `git_head.txt`              — `git rev-parse HEAD`
- `git_status_porcelain.txt`  — `git status --porcelain=v1`
- `git_diff_unstaged.patch`   — `git diff -- src/solveDiscreteFlowAdjoint.H src/solveDiscreteFlowAdjointProduction.H`
- `git_diff_staged.patch`     — `git diff --cached` (empty, as in cycle-1)
- `file_sha256.txt`           — sha256 of the probe + 2 production headers + stageB5BoundaryRelaxOracle.H
- `pre_state_readme.md`       — this file

## Verified facts

- git root: `/home/ys/dsH/TO-ANISOTROPIC` (exact workspace, checked via `git rev-parse --show-toplevel`).
- HEAD: `ca8a772b8339a36e7906ea32a7125061c47aa280` (matches cycle-1 pre-state).
- Dirty worktree: `M src/solveDiscreteFlowAdjoint.H`, `M src/solveDiscreteFlowAdjointProduction.H`
  (BFINAL-003 pre-existing user/production changes — FORBIDDEN to touch this cycle).
- Untracked: `build/`, `evidence/agent-group/BFINAL-00{2,3,4,4-FIX}/`, `src/stageB5BoundaryRelaxOracle.H`,
  `src/stageB6RxDesignOracle.H`, `src/validation_stage_*.log`.

## Acceptance checks (all PASS)

| check | cycle-2 | cycle-1 | match |
|---|---|---|---|
| HEAD | ca8a772b8339a36e7906ea32a7125061c47aa280 | ca8a772b8339a36e7906ea32a7125061c47aa280 | IDENTICAL |
| src/solveDiscreteFlowAdjoint.H sha256 | f0c814975cf8bc7f98b9f69efa2f0bc01bd251579920bbedda7d266b72927507 | f0c814975cf8bc7f98b9f69efa2f0bc01bd251579920bbedda7d266b72927507 | IDENTICAL |
| src/solveDiscreteFlowAdjointProduction.H sha256 | 8c4901bf8da4edbc7cc1efe73a7583963e4688bba93c8c58539c4c39cae1ad91 | 8c4901bf8da4edbc7cc1efe73a7583963e4688bba93c8c58539c4c39cae1ad91 | IDENTICAL |
| src/stageB5BoundaryRelaxOracle.H sha256 | 0b274ee26dc93c73d614f8c2af623085f25f9d9f6b008261f4efa3cda4de2d6b | 0b274ee26dc93c73d614f8c2af623085f25f9d9f6b008261f4efa3cda4de2d6b | IDENTICAL |
| git_diff_unstaged.patch | sha256 247dd33996eeef37bc194928c49653d428db732757df2f011b8ae806bf5cce5b | same | byte-IDENTICAL (production zero new diff) |
| git_status_porcelain.txt | sha256 6ddf344054373ea5a640d0aee32bc156219923208ee2733beafa9084b5619204 | same | byte-IDENTICAL |

## Probe state (expected, differs from cycle-1 pre-state by design)

`src/stageB6RxDesignOracle.H` sha256 = `51a02fc6288da612fc7e9b775cfd05ac757c137b492d5e863185f53c74c4b31d`
(cycle-1 pre-state had `f0b88ba9...`). This is EXPECTED: cycle-1 landed B1/B2/B3 fixes inside the
probe (only allowed source file). Current probe contains the cycle-2 defect surface to be fixed:

- L374: B3 no-V formulas present (`drAU = -rAU^2/alphaRel; dHbyA = (rAU/alphaRel)*((1-alphaRel)*U - HbyA)`).
- L521: `anRPaW[celli] = anRPa[celli]*da;`          — RX-A post-multiplication (weighting-order defect)
- L599: `anRPxh[celli] = anRPa[celli]*f;`           — RX-B post-multiplication (weighting-order defect)
- L857: `anRPd[celli] = anRPa[celli]*f;`            — RX-C post-multiplication (weighting-order defect)
- L782 / L1035 / L1121: `Foam::max(projEtaDenom, SMALL)` — D4 sign-clamp defect

No source file was modified by this capture step (evidence files only).
