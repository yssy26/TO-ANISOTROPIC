# Stage B2/B3 — Frozen-Turbulence Gradient-Amplitude Validation (FAIL)

- Repository: `/home/ys/TO-ANISOTROPIC`
- Branch: `agent/stage-a-anisotropic-validation`
- Baseline commit: `f3c4aa25e4ae46b4bed80143cc1ef81ab6c7c68f`
- OpenFOAM-7; case `/home/ys/b2_case` (gate4-derived, serial, deterministic PCG+DIC, 33600 cells); date 2026-08-07
- Config: `objectiveGradientScale 1.0; pressureGradientScale 1.0;` `freezeTurbulenceForValidation true`; `solveFlowAdjoints true`; `div(phi,T) bounded Gauss upwind`; MMA off; `gradientValidated false`.
- Scope: Stage B1.5 closure (direction-consistency self-check, limit-cycle phase determinism), B2 single-direction multi-step FD diagnostic, chain decomposition, B3 multi-direction formal acceptance. No MMA, no design update, no full-SST Gate 6.

## 1. Method

Two-phase baseline (full-SST 246 iterations -> frozen turbulence), identical to B1. Baseline: `J=-1.1673569`, `PressureDropPa=48336.3`, `gDP=0.933452`, `gV=-7.933e-9`. Directions D1/D2/D3: physical-coordinate low-frequency modes, designMask-only, bounds-safe, L-inf normalized to 1 (identical construction to B1). Adjoint projections `<dfdx,Dk>`, `<dgdx[1],Dk>`, `<dgdx[0],Dk>` computed on the baseline. Central FD along each direction at h in {1e-4,3e-4,1e-3,3e-3,1e-2} (D1) and {1e-4,3e-4,1e-3} (D2/D3), each side a full restore->design-chain recompute->frozen solve->objective evaluation.

## 2. B1.5 closure

- Direction-consistency self-check: the same Dk field object is used for the FD probes and the projections (no direction/ sign mismatch in the harness).
- Limit-cycle phase determinism: the deterministic PCG solver gives bit-reproducible plus/minus states for identical restarts (repeated 3 times in B1, all bit-identical); the frozen iterate limit cycle does not pollute the FD comparison (FD is stable across h, see Section 3).
- Pressure-normalization identity: `dgDP/dh == (1/PDmax) dPD/dh` verified to relDiff=2.5e-14 (PASS).

## 3. B2 single-direction multi-step FD (D1)

| h | FD_J | ADJ_J | relJ | signJ | FD_gDP | ADJ_gDP | relDP | FD_gV | relV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1e-4 | +0.0417192 | -0.0363597 | 2.147 | **0** | -2.53758 | -5.45652 | 0.535 | 0.155512 | 3.96e-5 |
| 3e-4 | +0.0417190 | -0.0363597 | 2.147 | **0** | -2.53784 | -5.45652 | 0.535 | 0.155513 | 3.81e-5 |
| 1e-3 | +0.0417182 | -0.0363597 | 2.147 | **0** | -2.53792 | -5.45652 | 0.535 | 0.155518 | 1.80e-6 |
| 3e-3 | +0.0416817 | -0.0363597 | 2.146 | **0** | -2.53606 | -5.45652 | 0.535 | 0.155518 | 1.84e-6 |
| 1e-2 | +0.0410765 | -0.0363597 | 2.130 | **0** | -2.51909 | -5.45652 | 0.538 | 0.155519 | 1.34e-7 |

FD_J is flat at +0.04172 across h=1e-4..3e-3 (asymptotic), so the D1 sign flip is not an FD artifact. ADJ_J = -0.03636 has the opposite sign; magnitudes differ by 14.7%.

## 4. B3 multi-direction formal acceptance

| dir | best-step FD_J | ADJ_J | relJ | signJ | best-step FD_gDP | ADJ_gDP | relDP | FD_gV relV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| D1 | +0.0417 | -0.0364 | 2.147 | **0** | -2.538 | -5.457 | 0.535 | 1.3e-7 (h=1e-2) |
| D2 | -0.0038 | -0.0196 | 0.806 | 1 | +0.507 | +1.179 | 0.570 | 3.4e-6 |
| D3 | -0.0003..-0.0004 | -0.0191 | 0.979-0.986 | 1 | +6.48 | +14.05 | 0.539 | 7.6e-6 |

Verdicts: J bestRel=0.806 (gate <=5%) FAIL; signOK=0 FAIL (D1 flip); gDP bestRel=0.535 (gate <=10%) FAIL; gV bestRel=1.3e-7 (gate <=1e-6) PASS; plateauOK=1.

## 5. Findings

1. **[FAIL] Objective-function adjoint gradient is wrong in amplitude and in D1 sign.**
   - D1: sign flip (FD +0.0417 vs ADJ -0.0364), all five steps.
   - D2: same sign, ADJ/FD ratio 5.2.
   - D3: same sign, ADJ/FD ratio 50-70 (D3 FD is small and h-dependent; not asymptotic).
   - The error is direction-dependent, so it is NOT a single global sign or a single global scale in the sensitivity assembly; it points to the thermal objective derivative chain (dJ/dT/dJ/dphi, the Tb adjoint, or the fsensMeanT assembly) rather than the filter/projection chain.
2. **[FAIL] Pressure-drop adjoint is systematically ~2.15-2.33x the FD value** on all three directions (D1 2.15, D2 2.33, D3 2.17), same sign. A uniform multiplicative error in gsensPressureDrop (or its normalization) is indicated. (Known since B0.2 F2: PD maxRelError 0.58, 9/9 fail - now quantified as a consistent factor.)
3. **[PASS] Volume-constraint chain** (x->xp->xh->gV) is correct: FD_gV vs <gsensVol,Dk> rel error 1.3e-7..7.6e-6 depending on h/direction, sign consistent. The filter/projection chain and the design derivative are trustworthy.
4. **[PASS] The frozen-flow FD machinery is reliable**: deterministic (bit-reproducible), pressure identity holds to 2.5e-14, FD stable across h (D1), so the FAIL is in the adjoint assembly, not in the primal/FD side.
5. **[INFO] D3 objective FD is not asymptotic** (FD_J changes by 1.5x from h=1e-4 to h=1e-3); D3 objective sensitivity is weak. Does not affect the FAIL conclusion.

## 6. Recommendation (for review)

The objective and pressure adjoints must be debugged before any gradient-based step:
- Objective: audit the full chain J -> (dJ/dT, dJ/dphi) -> Tb adjoint source -> Tb solve -> fsensMeanT assembly, focusing on the thermalObjectiveMask/domain-integrate terms and any missing factor; a single-direction "slice" test (one design cell) comparing ADJ vs FD cell-by-cell would localize the error.
- Pressure: audit gsensPressureDrop assembly for a uniform factor (~2.15x) - candidates: 1/rho, area factor, PDmax normalization, or a missing -1.
- Volume chain and the primal FD harness need no changes.

Status: **FAIL** (J FAIL amplitude+sign, gDP FAIL amplitude, gV PASS). Next action pending external review.

## 7. Deliverables

- `STAGE_B_B2B3_GRADIENT_AMPLITUDE_REPORT.md`
- `stage_b_b2b3_fd_scan.tsv` (all 11 rows)
- `stage_b_b2b3_summary.json`
