# BFINAL-030 PREREGISTRATION — T3/(1-alphaRel) 语义审判轮

Date: 2026-08-21.  Written BEFORE any fingerprint data is computed.
Decisions made in Stage 0 derivation (see b30_NOTEBOOK.md); numbers below are
the pre-registered expectations against which the measured values will be
judged.

## Adjudicated correction coefficient

- **c* = 1.0** (primary, FIXED-POINT semantics): at SIMPLE outer-loop
  convergence the map design -> (U,p) solves the steady residual; the relax
  factor is a per-iteration preconditioner and its single-iteration derivative
  `(1-alphaRel)*dU` disappears at the fixed point.  Correct steady-state flux
  tangent = `Sf&dU`, i.e. T3 coefficient should be 1.0, not 0.6.
- Alternates (kept for verdict sensitivity):
  - c* = 1.09  (empirical U-row fit: 0.6*1.82)
  - c* = 1.23  (gDP-only full-explanation fit: 0.6*2.05)

Patch recipe (fingerprint 1, coefficient-level, ZERO production changes):
  M_corrected = M + (c* - 0.6) * T3_base
  T3_base = sparse pattern with coefficient 1 in the U-row/P-column block:
    internal: [U(own,c),P(own)] += w*sf_c ; [U(own,c),P(nei)] -= w*sf_c
              [U(nei,c),P(own)] += (1-w)*sf_c ; [U(nei,c),P(nei)] -= (1-w)*sf_c
    assignable boundary: [U(cell,c),P(cell)] += patchSf_c

## Fingerprint 2 — corrected-T3 route vs exported b_TC (thermalCoupling)

- Route rebuild: b29_slot_instrument.py T3 slot coefficient `(1.0-ALPHA_REL)`
  = 0.6 -> c* = 1.0 (candidate change; T6 boundary slot left at 0.6 first,
  re-tested with 1.0 as a secondary check).
- Baseline (existing): route("Tmix", FULL-{"T8"}) relL2 = 0.39701093640754337.
- Pre-registered: relL2 drops from 0.397 to < 5% (0.05), U row amplitude gap
  1.82x -> ~1.0.  This is the direction-of-the-hypothesis statement; a
  residual mismatch that stays > 30% falsifies the T3-only attribution.

## Fingerprint 3 — corrected b_TC x existing w_true -> M1 three-direction ratios

- b_TC^T w_true / RHS with RHS = {D1:0.07443678641907196,
  D2:0.0006387577668902512, D3:0.01094029857921596}.
- Pre-registered: ratios (b_TC^T w_true / RHS) must satisfy
  |ratio - 1| <= 15% on ALL THREE directions simultaneously, with the primary
  c* = 1.0.  (D2 is a known large-cancellation small-residue row; noted as
  possible instability, not a verdict blamer by itself.)

## Fingerprint 1 — coefficient-patched explicitJT_H7 re-solve

- Patch U-row/P-column block of `/home/ys/dsH/b24_diag/explicitJT_H7.mtx`
  (shape 134400x134400, nnz 6395416), rescale T3 entries by c*/0.6.
- Re-solve w_true with the b26_wstar_h7.py splu recipe in background
  (~17 min, poll).
- Main criterion: corrected b_TC^T (corrected w_true) closure on all three
  directions (|ratio-1| <= 15%).
- Secondary: solve lambda_PD against `/home/ys/dsH/b25_qgate/b18rhs_pressureDrop.mtx`
  (B13 decomposition), predict gDP factor.
- Pre-registered gDP prediction: ADJ/FD direction-uniform factor goes from
  2.05x -> 1.0 +/- 0.15 (i.e. |factor-1| <= 0.15).  FD_gDP values read from
  the existing `b25 stage_b2_fd_scan.tsv`; NO FD re-run.
  - c*=1.0 => gDP factor predicted ~1.0 (primary).
  - If only partial (T3-direct-only with code-path 0.6), predicted factor
    ~2.05/1.667 ~ 1.23 (fails the |factor-1|<=0.15 criterion).

## Verdict mapping (pre-registered)

- 定罪 (conviction): fingerprint 2 < 5%, fingerprint 3 all-three within 15%,
  fingerprint 1 main criterion within 15% AND gDP within 1.0 +/- 0.15, all
  with the same c* = 1.0.
- 部分 (partial): at least one fingerprint moves in the preregistered
  direction but not all three meet their thresholds; report per-fingerprint
  and quantify residual.
- 无罪 (innocent): no fingerprint moves toward the preregistered numbers.

## Zero-production-change guard

`src/` is NOT touched this round.  All corrections are offline, coefficient-
level (route instrument script and matrix/RHS text files only).
