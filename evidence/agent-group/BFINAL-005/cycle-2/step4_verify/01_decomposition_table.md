# BFINAL-005 cycle-2 — step4_verify: sensitivity decomposition table (independent)

Independent recomputation from raw `.mtx` artifacts in
`cycle-2/step3_build_run/artifacts/` (script: `verify_gates.py`, no shared
formula with the C++ probe — oracle anchors are RE-COMPUTED here from
`stageB6_lambda.mtx` × `stageB6_rxc_analytic.mtx`, not copied from log lines).

- nC = 33600 cells; lambda file = [3nC Uc-lambda ; nC pc-lambda];
  rxc_analytic = 3 dirs × [3nC R_U,xd ; nC R_P,xd]; z = dxh/dx·d per dir.
- D_pressure = -Σ pc·R_P,xd ; D_momentum = -Σ Uc·R_U,xd ; D_total = sum.
- G1 production = Σ gsensh_pressurerow·z ; G2 total = Σ gsensh_total·z ;
  G3 raw = Σ rxpr_prod_gsens·d (post-chain, patched production field).

## 1. Decomposition table (D1/D2/D3)

| dir | D_momentum oracle | D_pressure oracle | D_total oracle | G2 D_momentum prod | G1 D_pressure prod | G2 D_total prod | G3 raw-design prod | relErr G1 | relErr G2 | relErr G3 | \|D_pressure/D_total\| |
|---|---|---|---|---|---|---|---|---|---|---|---|
| D1 | -0.0486867554277 | -0.00754887795106 | -0.0562356333787 | -0.0486867554277 | -0.00754887795106 | -0.0562356333787 | -0.0562356332429 | 4.37e-15 | 1.23e-15 | 2.42e-09 | 0.134237 |
| D2 | +0.229079662362 | +0.0685927693909 | +0.297672431753 | +0.229079662362 | +0.0685927693909 | +0.297672431753 | +0.297672431913 | 2.02e-15 | 2.42e-15 | 5.38e-10 | 0.230430 |
| D3 | -0.0301314016208 | -0.00621002821874 | -0.0363414298395 | -0.0301314016208 | -0.00621002821874 | -0.0363414298395 | -0.0363414298017 | 4.55e-14 | 1.26e-13 | 1.04e-09 | 0.170880 |

Full precision (from `verify_gates_output.txt`):

- D1: mom -4.8686755427656533e-02, press -7.5488779510572059e-03, total -5.6235633378713741e-02
  - G1 prod -7.5488779510571730e-03 (relErr 4.366e-15), G2 mom prod -4.8686755427656665e-02 (relErr 2.708e-15),
    G2 total prod -5.6235633378713810e-02 (relErr 1.234e-15), G3 raw prod -5.6235633242873749e-02 (relErr 2.416e-09)
- D2: mom +2.2907966236169880e-01, press +6.8592769390886446e-02, total +2.9767243175258523e-01
  - G1 prod +6.8592769390886585e-02 (relErr 2.023e-15), G2 mom prod +2.2907966236169891e-01 (relErr 4.846e-16),
    G2 total prod +2.9767243175258451e-01 (relErr 2.424e-15), G3 raw prod +2.9767243191274928e-01 (relErr 5.381e-10)
- D3: mom -3.0131401620773057e-02, press -6.2100282187409978e-03, total -3.6341429839514057e-02
  - G1 prod -6.2100282187407150e-03 (relErr 4.553e-14), G2 mom prod -3.0131401620776741e-02 (relErr 1.223e-13),
    G2 total prod -3.6341429839518644e-02 (relErr 1.262e-13), G3 raw prod -3.6341429801673833e-02 (relErr 1.041e-09)

Decomposition self-check: max|gsensh_total - (gsensh_momentum + gsensh_pressurerow)| = 5.204e-18 (cell-wise).

## 2. Sign convention

All production projections have IDENTICAL sign to the oracle for D1/D2/D3 on
G1/G2/G3 (J^T λ = g_w convention preserved; no sign flip to fit anchors).

## 3. Notes

- G1 (pressure-row local contraction) closes to **4.4e-15 / 2.0e-15 / 4.6e-14**
  — the transpose operator T = J_P^T pc applied via `-T*dAlphaDxh` contracts
  exactly to -pc^T R_P,x d. This is the DIRECT proof that the production
  pressure-row term is the true -λ_P^T R_P,x (and not the cycle-1 forward
  +pointwise form, which would fail here).
- G2 (total) closes to **1.2e-15 / 2.4e-15 / 1.3e-13**; |D_pressure/D_total| =
  13.4% / 23.0% / 17.1% — the pressure row is a substantial fraction of the
  total design derivative, confirming BFINAL-004's INCOMPLETE finding.
- G3 (raw-design projection through filter+projection chain) closes to
  **2.4e-9 / 5.4e-10 / 1.0e-9** — the patched production raw-design field
  `rxpr_prod_gsens.mtx` (which goes through the EXISTING filter_chainrule.H
  path) reproduces the total design derivative to ~1e-9.
- G5 volume projection (independent dot from `stageB6_prod_gsensVol.mtx` × dirs):
  D1 +0.1555186511144686, D2 -0.08449667741658891, D3 +0.06396237684432686 —
  byte-identical to the log anchors and to BFINAL-004-FIX cycle-3 (volume chain
  unchanged by the patch).
