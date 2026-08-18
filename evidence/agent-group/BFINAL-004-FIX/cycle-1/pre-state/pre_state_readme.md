# BFINAL-004-FIX cycle-1 pre-state baseline (step1, read-only)

Captured by Executor (deepseek-v4-flash) before ANY source modification.
No file outside `evidence/agent-group/BFINAL-004-FIX/` was created or changed.

## Files

- `git_head.txt` — HEAD = `ca8a772b8339a36e7906ea32a7125061c47aa280`
- `git_status_porcelain.txt` — `git status --porcelain=v1` (34 lines, below)
- `git_diff_unstaged.patch` — full unstaged diff (25229 bytes)
- `git_diff_staged.patch` — empty (0 bytes: no staged changes)
- `file_sha256.txt` — sha256 of probe + adjacent sources:
  - `src/stageB6RxDesignOracle.H` = `f0b88ba9c0c885d01e140ff3b958cdcc990a4792947f08eaa36e234481f35fe1`
  - `src/stageB5BoundaryRelaxOracle.H` = `0b274ee26dc93c73d614f8c2af623085f25f9d9f6b008261f4efa3cda4de2d6b`
  - `src/solveDiscreteFlowAdjoint.H` = `f0c814975cf8bc7f98b9f69efa2f0bc01bd251579920bbedda7d266b72927507`
  - `src/solveDiscreteFlowAdjointProduction.H` = `8c4901bf8da4edbc7cc1efe73a7583963e4688bba93c8c58539c4c39cae1ad91`

## Verification performed (read-only, all OK)

1. `pwd` = `/home/ys/dsH/TO-ANISOTROPIC`; `git rev-parse --show-toplevel` = same → repo root == workspace.
2. HEAD matches the approved plan (`ca8a772b`).
3. Saved `git status --porcelain=v1` diffed byte-identical against a fresh `git status --porcelain=v1` run after capture (STATUS_MATCH=OK).
4. Saved `git_diff_unstaged.patch` diffed byte-identical against a fresh `git diff` (UNSTAGED_DIFF_MATCH=OK); staged diff empty on both sides.
5. `sha256sum src/stageB6RxDesignOracle.H` after capture == value in `file_sha256.txt` (SHA_MATCH=OK, no drift).

## Worktree state at baseline (matches plan `current_state`)

- Modified (BFINAL-003, FORBIDDEN to touch): `src/solveDiscreteFlowAdjoint.H`, `src/solveDiscreteFlowAdjointProduction.H`
- Untracked: `build/`, `evidence/agent-group/BFINAL-002|003|004|004-FIX/`, `src/stageB5BoundaryRelaxOracle.H`, `src/stageB6RxDesignOracle.H`, `src/validation_stage_*.log`

## Confirmed defect sites in `src/stageB6RxDesignOracle.H` (for step2)

- B3 (RX-B P-row): L383 `drAU[celli] = -alphaRel*rAUc*rAUc;` and L384-389
  `dHbyA[celli] = rAUc*((alphaRel-1.0)*U[celli] - alphaRel*HbyABase[celli]);`
  — matches plan's wrong-code description (alphaRel vs 1/alphaRel mixing; U-term sign).
  Header comment L34-36 carries the same wrong formulas.
- B1 (RX-A): L552-553 `stageB6Metrics(anRUa, fdRUa, ...)` / `stageB6Metrics(anRPa, fdRPa, ...)`
  compare UNWEIGHTED analytic vs deltaAlpha-weighted FD; `anRUa` built at L364-370 (V*U),
  `anRPa` built at L453-483.
- B2 (RX-C): `stageB6_d` (L707) and `stageB6_y` (L714) constructed with default
  (calculated) BCs before Helmholtz tangent solve at L717-726; `gsensVolR` (L1009) and
  `gsensR` (L1092) likewise for §9/§10 filter solves (L1023-1024, L1107-1108).

## No changes made

`git status` after capture is byte-identical to the saved baseline (nothing added/modified
outside evidence/). Production sources untouched.
