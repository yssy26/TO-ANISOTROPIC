# BFINAL-009 evidence — GAMG preconditioner run（240 迭代部分收敛数据 / partial, 240-iteration data）

> **注意：本目录日志为 240 迭代部分收敛数据 —— 运行中途被终止，日志截断，求解未达到容差。
> NOT a converged run. The solver was killed mid-run; the log is truncated at
> cycleIter=240 k=8 with no completion marker.**

## Provenance

- Case: `/home/ys/dsH/b8_verify_gamg` (production verify case, `discreteProdPreconditionerSetup` = pressure-GAMG, commit 4c49903 "Add production pressure-GAMG preconditioner")
- Exec: `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF`, OpenFOAM 7 (Build 7-3bcbaf946ae9)
- Repo HEAD at run time: `6a0004b` (`agent/dsH-stage-b-validation`), tracked tree clean
- Run started: 2026-08-19 00:35:03 (PID 388624, host DESKTOP-DSG34PK); log last write 00:35, process no longer alive → run terminated externally
- Archived: 2026-08-19, copied verbatim from `b8_verify_gamg/Log.verify_gamg.txt`

## Files

| file | sha256 | note |
|---|---|---|
| `Log.verify_gamg.txt` | `c944bb213fcedd34814793a7ecc724d2d9168d78e986e8d7dfea1caa11b7475b` | 128 KB, truncated partial log (ends mid restart-cycle at cycleIter=240 k=8) |
| `adjointCheckpoint_thermalCoupling.tsv` | — | 0 bytes at kill time (checkpoint never written) |

## Key numbers (thermalCoupling production FGMRES + pressure-GAMG)

FGMRES-PROD trueRelRes checkpoints (tolerance 1e-9):

| iter | trueRelRes | GRADPROXY | PRODDEFL cos(Au,u) |
|---|---|---|---|
| 80  | 3.0447e-4 | 1348.76 | 0.2199 |
| 160 | 1.2072e-4 | 1935.18 | 0.1532 |
| 240 | 4.9066e-5 | 2144.25 | 0.1410 |

- GAMG inner solves (`prodPressurePrecPsi_thermalCoupling`): 249 applications visible, each 5–8 V-cycle iterations, inner final residual ~4e-4–1e-3 (inner tolerance appears not tightened).
- Comparison, fresh diagonal-preconditioner run (BFINAL-008 cycle-1, `VERIFY_FRESH_6320943.md` §7): thermalCoupling stalled at 5.08e-5 after **4000** iterations. GAMG reaches ~4.9e-5 by **240** iterations → same stall floor, reached ~17x sooner, but still NOT converged to 1e-9.
- Deflation cosine degrades slowly (0.22 → 0.14); gradient proxy grows — consistent with the cond~1e19 stall already documented.
- `adjointCheckpoint_thermalCoupling.tsv` empty → no gradient checkpoint was produced by this run.

## Interpretation (for BFINAL-009 analysis)

Pressure-GAMG accelerates the early phase dramatically but does not break the
~5e-5 trueRelRes stall floor of the ill-conditioned physical operator on the
33600-cell case. The 240-iteration data above is partial evidence only; any
GAMG conclusion must be re-anchored on a run that either converges or is
stopped at a pre-declared iteration budget.
