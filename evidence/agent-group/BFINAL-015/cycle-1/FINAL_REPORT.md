# BFINAL-015 cycle-1 FINAL REPORT — T2/T3 tangent-vs-state comparison

## Verdict

**P-A confirmed in essence; preregistered sub-prediction falsified and
reported; one defect beyond the P-A/P-B framework discovered.** This round
delivered the first external validation of λ/w′ and split the BFINAL-012
failure into TWO independent defects:

1. **gDP misfit (≈2.13x) = J-operator defect.** With the source and the
   state truth machine-precision closed, the coded tangent w′ = J⁻¹(−R_x d)
   differs from the true state response by direction-dependent factors
   1.705 / 0.809 / 2.200 (D1/D2/D3) — J (both the production matrix-free
   AND the exported explicit variant, which differ from each other) is NOT
   the Jacobian of the phi-eliminated fixed-point map. Mismatch
   concentrates in the **P rows** (relP 0.72–1.13 vs relU 0.18–0.53).
2. **J-objective misfit = source/elimination defect (b_TC), independent of
   the operator.** Even with the ground-truth state response w_true, the
   thermalCoupling source contraction fails: b_TC^T w_true = 0.1125 /
   −0.0256 / −0.0448 vs FD_J = 0.0429 / −0.00347 / −0.000283 (factors
   2.6x / 7.4x / 158x). The thermal-adjoint elimination chain that folds
   dJ/dT and dJ/dphi into (U,p) rows is itself defective.

## Measurement chain and integrity proofs

- **Pre-registration** (`PREREGISTRATION.md`) committed before any w′/w_true
  data existed (timestamped by this commit chain).
- **Export run** (`b15_export`, export-only diagnostic path, MTO_RC=0):
  ExplicitJT-vs-matrix-free oracle **3.757e-16** (bit-identical to the
  BFINAL-008 fresh value), transpose dot tests 6.86e-14, RX-A design-FD
  ~1e-7 — export integrity proven despite the (irrelevant here) degenerate
  diagnostic SOLVE path (export-only skips it).
- **State-truth run** (`b13_probe3`, B2 module + new switch-gated
  `stageB15StateExport`): full FD scan reproduced BFINAL-012 bit-for-bit
  (FROZEN_GRADIENT_STATUS=FAIL, MTO_RC=0); 30 state files (baseline + 14
  probes at h ∈ {3e-4, 1e-3, 3e-3}).
- **Offline** (`b15_t2t3.py`, venv_s4): splu tangent solves trueRelRes
  5.3e-14 / 4.3e-14 / 9.2e-14 (≤1e-10 gate ✓); rhs norms bit-identical to
  the BFINAL-008 P8 historical values (5.663185e-01 / 8.615839e-01 /
  8.516528e-01).
- **Layout pitfall (documented)**: the rhs export is per-cell interleaved
  (solveDiscreteFlowAdjoint.H:3113-3119) while the matrix export and
  rxc_analytic are block-major — first-pass contractions were scrambled and
  are preserved in the console log history only as the debugging trail; all
  numbers below use the corrected de-interleaved loader.

## Decisive numbers (h = 1e-3 formal step; ladder stable, see table)

### Truth-chain closure (machine precision)

| dir | b_PD^T w_true | BFINAL-012 FD_gDP | verdict |
|---|---|---|---|
| D1 | −2.537460 | −2.53745977185 | ✓ closes |
| D2 | +0.523409 | +0.523409452158 | ✓ closes |
| D3 | +6.477216 | +6.47721554707 | ✓ closes |

Source (P1) + state truth (this round) + objective FD (BFINAL-012) are one
consistent triple. Everything below is measured against this closed truth.

### Defect 1 — J operator (gDP side)

| dir | b_PD^T w′ (exported J) | b_PD^T w_true | ratio | production ADJ | ADJ/FD |
|---|---|---|---|---|---|
| D1 | −4.326003 | −2.537460 | **1.705** | −5.405059 | 2.130 |
| D2 | +0.426151 | +0.523409 | **0.814** | +1.105927 | 2.113 |
| D3 | +14.252071 | +6.477216 | **2.200** | +14.044714 | 2.168 |

- w′ ≠ w_true with direction-dependent ratios → the exported J is not the
  true fixed-point Jacobian (P-A core confirmed; NOT a scalar error).
- The production λ differs from the exported-J adjoint as well (−4.326 vs
  −5.405 on D1; 1.5% agreement on D3) → **two different, differently-wrong
  operators**; the production matrix-free and the diagnostic/export path
  disagree with each other by up to 61% (D2).
- State-level mismatch profile (D1): cos 0.9981, relL2 0.725;
  **relP 0.7247 vs relU 0.1802** — mismatch concentrated in the P rows
  (preregistered "U-row/velocityJump concentration" is FALSIFIED; reported
  per the pre-registration discipline). Top mismatch cells come in
  face-adjacent pairs (7598/7599, 8718/8719) — face-flux/pressure-coupling
  structure.

### Defect 2 — thermalCoupling source/elimination (J side)

| dir | b_TC^T w_true | FD_J | factor | b_TC^T w′ | w′/w_true |
|---|---|---|---|---|---|
| D1 | +0.112506 | +0.042856 | 2.63 | +0.152532 | 1.356 |
| D2 | −0.025593 | −0.003469 | 7.38 | −0.035935 | 1.404 |
| D3 | −0.044814 | −0.000283 | 158 | −0.065202 | 1.455 |

Even the ground-truth state response fails to close the J contraction →
**b_TC (the eliminated dJ/d(U,p) assembled in AdjNS_HT + the flux-transpose
folding) is defective independently of the operator**. This also explains
the end-to-end J sign flips: a wrong source contracted through a wrong
operator.

### eps ladder stability (w_true contractions, D1)

gDP: −2.537102 / −2.537460 / −2.535369 (3e-4 / 1e-3 / 3e-3, spread <0.1%);
J: 0.112511 / 0.112506 / 0.112245. All mismatches are systematic, not
noise.

## Decision table

| preregistered prediction | outcome |
|---|---|
| P-A core: J not the true fixed-point Jacobian → w′ ≠ w_true | **CONFIRMED** (ratios 0.81–2.20) |
| P-A sub: mismatch concentrated in U rows, velocityJump pattern | **FALSIFIED** — P rows dominate (relP 0.72–1.13 vs relU 0.18–0.53); face-adjacent-pair structure |
| P-B: w′ ≈ w_true → contradiction stop | NOT triggered |
| Beyond both: b_TC source defect | **DISCOVERED** — J-side misfit persists with ground-truth w |

## Localization status for the next (user-authorized) rounds

- **Operator round (gDP)**: the P-row blocks of J — and the production-vs-
  exported operator divergence (they disagree with each other by up to 61%,
  so at most one can be right; both are wrong vs truth) — specifically the
  pressure-row / flux-coupling semantics (rAtU/relaxed-SIMPLE composition).
  The 97.6%-carrier momentum contraction from BFINAL-013 inherits whatever
  λ this operator produces.
- **Source round (J)**: the thermalCoupling rhs assembly / thermal-adjoint
  elimination (AdjHeatTransfer + discreteExternalFaceFluxAdjoint folding in
  AdjNS_HT + the diagnostic rhs export path) — externally measurable now
  via b_TC^T w_true vs FD_J (this round's new instrument, permanent).
- No fixes were attempted in this round (pure measurement, per discipline).

## Files

| file | note |
|---|---|
| `PREREGISTRATION.md` | pre-registered predictions (before data) |
| `b15_t2t3.py` | offline tangent + comparison instrument |
| `t2t3_console.log`, `t2t3_results.json`, `wprime_cache.npz` | analysis output |
| `Log.verify_b015_export.txt` | export run (integrity markers) |
| `Log.verify_b015_states.txt`, `stage_b2_fd_scan.states_run.tsv` | state-truth run (B12 bit-reproduction) |
| `stageB6_rxc_analytic.mtx`, `explicitRhs_*.mtx` | tangent rhs / sources (checksummed) |
| `sha256_artifacts.txt` | checksums incl. explicitJT.mtx (467 MB, not committed) |
| `optProperties.*.txt` | case snapshots |
