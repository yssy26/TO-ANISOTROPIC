# BFINAL-004-FIX cycle-3 — Pre-edit baseline (step1-pre-state)

- Role: Executor (PATCH cycle first step, read-only)
- Created: 2026-08-17 (before any source edit in cycle-3)
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` == git root (verified by
  `git rev-parse --show-toplevel`; pwd matches).
- HEAD: `ca8a772b8339a36e7906ea32a7125061c47aa280` (branch
  `agent/dsH-stage-b-validation`) — recorded in `git_head.txt`.

## Files in this directory

| file | content |
|---|---|
| `git_head.txt` | `git rev-parse HEAD` + `--abbrev-ref HEAD` |
| `git_status_porcelain_v1.txt` | `git status --porcelain=v1` (34 lines) |
| `git_status_full.txt` | `git status` verbose snapshot |
| `git_diff_unstaged.patch` | `git diff` (497 lines: the 2 BFINAL-003 production files) |
| `git_diff_staged.patch` | `git diff --cached` (0 bytes — nothing staged) |
| `file_sha256.txt` | sha256 of the probe header + every files_forbidden module (22 entries) |
| `pre_state_readme.md` | this file |

## Acceptance verification (step1)

1. **pre-state dir exists** — yes, populated above.
2. **HEAD = ca8a772b** — `ca8a772b8339a36e7906ea32a7125061c47aa280` ✓
   (matches approved-plan current_state and BFINAL-003 FINAL_REPORT HEAD).
3. **Probe header pre-hash recorded** — `src/stageB6RxDesignOracle.H` =
   `09832dd6ae5a3c7f37201b200bfb19dc0c48f42a3989f9593979c73ca098baec`
   (matches the plan's stated on-disk hash `09832dd6ae5a...`; 61807 bytes).
   Note: this differs from cycle-1 (`f0b88ba9...`) and cycle-2
   (`51a02fc6...`) pre-state hashes — the on-disk candidate already carries
   the self-reported cycle-2 B1/B2/B3 edits + OFstream flush calls. This
   hash is the authoritative pre-edit baseline for cycle-3.
4. **Production heads match BFINAL-003 locked hashes** —
   - `src/solveDiscreteFlowAdjoint.H` =
     `f0c814975cf8bc7f98b9f69efa2f0bc01bd251579920bbedda7d266b72927507`
     == BFINAL-004-FIX cycle-1/cycle-2 pre-state record (post-BFINAL-003 state) ✓
   - `src/solveDiscreteFlowAdjointProduction.H` =
     `8c4901bf8da4edbc7cc1efe73a7583963e4688bba93c8c58539c4c39cae1ad91`
     == BFINAL-003 `cycle-1/file_sha256_production_patched.txt` (locked
     post-patch hash) and cycle-1/cycle-2 pre-state records ✓
   - Both are the BFINAL-003-locked dirty-worktree versions (git `M`),
     preserved untouched; the unstaged diff patch is archived in
     `git_diff_unstaged.patch`.
5. **All other forbidden files show no diff** — `git diff --quiet` against
   every other files_forbidden module (AdjNS_PD.H, AdjNS_HT.H, AdjNS_FF.H,
   NS.H, sensitivity.H, costfunction.H, computeObjective.H, filter_x.H,
   filter_chainrule.H, updateMaterialProperties.H, update.H,
   updateFrozenTurbulenceFields.H, HeatTransfer.H,
   createFrozenHotRegionFields.H, validateStageB2GradientAmplitude.H,
   validateFrozenGradient.H, MMA.h, MTO_HF.C, stageB5BoundaryRelaxOracle.H)
   exits 0 — zero diff vs HEAD for all tracked forbidden files.
   `src/stageB5BoundaryRelaxOracle.H` =
   `0b274ee26dc93c73d614f8c2af623085f25f9d9f6b008261f4efa3cda4de2d6b`
   == cycle-1/cycle-2 pre-state record (BFINAL-002 probe, preserved) ✓

## Scope note

Nothing was modified in this step: all commands were read-only except
creating this evidence directory under `evidence/agent-group/BFINAL-004-FIX/`
(permitted by files_allowed). No source file, no optProperties, no build
artifact was touched.
