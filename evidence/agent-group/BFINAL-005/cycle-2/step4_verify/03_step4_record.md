# BFINAL-005 cycle-2 — step4_verify record (executor)

- Task: BFINAL-005 PATCH_RX_PRESSURE_ROW — step4_verify: independent G1-G5
  recomputation from raw artifacts + decomposition/gate tables.
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` == git root (verified); HEAD
  `ca8a772b8339a36e7906ea32a7125061c47aa280`; branch `agent/dsH-stage-b-validation`.
- Method: `verify_gates.py` — a PURE PYTHON script that reads ONLY the raw
  `.mtx` artifacts (45 files in `cycle-2/step3_build_run/artifacts/`); NO code
  shared with the C++ probe. The oracle anchors are RE-COMPUTED here from
  `stageB6_lambda.mtx` × `stageB6_rxc_analytic.mtx` (D_pressure = -Σ pc·R_P,xd,
  D_momentum = -Σ Uc·R_U,xd), then compared against the production contractions.

## 1. What was verified

- **G1** D_pressure_production = Σ(gsensh_pressurerow·z) vs D_pressure_oracle:
  relErr 4.4e-15 / 2.0e-15 / 4.6e-14 (D1/D2/D3), sign same → PASS (≪1e-4).
- **G2** D_total_production = Σ(gsensh_total·z) vs D_total_oracle (recomputed
  D_momentum+D_pressure): relErr 1.2e-15 / 2.4e-15 / 1.3e-13, sign same,
  |D_pressure/D_total| = 0.134237 / 0.230430 / 0.170880 → PASS.
  Momentum component closes too (2.7e-15 / 4.8e-16 / 1.2e-13).
- **G3** raw-design projection Σ(rxpr_prod_gsens·d) vs D_total_oracle:
  relErr 2.4e-9 / 5.4e-10 / 1.0e-9, sign same → PASS. Volume projection
  unchanged (G5).
- **G4** R_w/J no regression — re-grepped raw log: momentum relL2
  1.65311363075e-4, P-total 1.52828911949e-4, J/J^T dot 6.42244770078e-14,
  explicit maxRelL_U 3.52693970502e-16 / maxRelL_P 3.98062510668e-16, GatePR
  identical across both occurrences → PASS.
- **G5** volume proj D1 0.1555186511144686 (independent dot identical),
  filter/projection RX-C tangent closures (cos>0.9999 at eps=1e-4),
  Stage-A anisotropy solidK=(21.6 0 0 21.6 0 16.6), discrete objective
  derivative errors thermal=2.70e-11 pressure=7.08e-12 → PASS.

## 2. Cross-check vs run-log §9/§10 anchors

Log-printed D_pressure (D1 -0.00754887795106, D2 +0.0685927693909,
D3 -0.00621002821874) and D_total (D1 -0.0562356333787, D2 +0.297672431753,
D3 -0.0363414298395) are reproduced to the last printed digit by the
independent recomputation (see `01_decomposition_table.md`). The production
contraction Σ(gsensh_pressurerow·z) equals the recomputed oracle to 4.6e-14
worst case — i.e. the production pressure-row term IS -λ_P^T R_P,x d.

## 3. Evidence written

- `verify_gates.py` (independent script)
- `verify_gates_output.txt` (raw stdout)
- `01_decomposition_table.md` (sensitivity decomposition D1/D2/D3)
- `02_gate_table.md` (G1-G5 gate table with raw log quotes)

## 4. Acceptance (step4)

| criterion | result |
|---|---|
| G1-G5 all PASS per pass criteria | PASS (all; relErr ≤ 4.6e-14 for G1/G2, ≤ 2.4e-9 for G3) |
| decomposition table for D1/D2/D3 with sign + rel err + \|D_pressure/D_total\| | PASS (01_decomposition_table.md) |
| gate table written to cycle-2 | PASS (02_gate_table.md) |
| independent script, no shared formula | PASS (pure python from .mtx; oracle recomputed, not copied) |
| G4/G5 anchors from raw log re-grep | PASS (quoted raw lines) |

## 5. Scope integrity at step4 exit

- No source file modified in this step (evidence writes only under
  `evidence/agent-group/BFINAL-005/cycle-2/step4_verify/`).
- `git status` unchanged from step3: only `src/sensitivity.H` (authorized
  patch) + the two pre-existing LOCKED BFINAL-003 heads `M`;
  `src/rxPressureRowTranspose.H` untracked (authorized new helper).
- No forbidden module touched; no J/J^T change; no empirical coefficients;
  no D1/D2/D3 hard-coding (anchors only used as comparison targets in the
  gate table, never in production code).
