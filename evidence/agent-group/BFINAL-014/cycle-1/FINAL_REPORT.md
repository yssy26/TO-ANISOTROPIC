# BFINAL-014 cycle-1 FINAL REPORT — derivation round: STOP (fix requires operator layer)

## Verdict

**STOP per the task's rule.** The phase-0 derivation (from NS.H's actual
forward path, `DERIVATION.md` in this directory) does **not** support any
change to the `sensitivity.H` contraction terms: the current
`−dAlphaDxh·(U & U_adj)·V` formula is exactly `−λ_U^T(∂R_U/∂α)` of the
true fixed-point system. The remaining defect converges on **λ itself —
i.e., the J operator's semantics as the Jacobian of the phi-eliminated
fixed-point system — which is a verified operator layer outside this
round's authorized scope.**

## Phase-0 derivation outcome (hypothesis adjudication)

| hypothesis | verdict | basis |
|---|---|---|
| H-F1 relaxation-semantics scalar mismatch | **REJECTED** | OF7 `fvMatrix::relax` (fvMatrix.C:521-670) updates the source with `(D_rel−D0)·ψ`; at the fixed point (ψ_prevIter=U) the relaxation cancels exactly and the fixed-point residual is the UNRELAXED momentum equation — independent of α_rel. No 1/α_rel factor can enter. Additionally the measured ADJ/FD ratio is not a pure scalar (2.113–2.169, 1.3% spread across directions). |
| H-F2 alpha-dependence beyond Sp(alpha,U) | **REJECTED** | `nuEffFrozen = nu + nutFrozen` (updateFrozenTurbulenceFields.H:60, alpha-independent); fvOptions empty (log-confirmed); the pressure-row mobility derivative `drAU=−rAU²/α_rel` already matches the derived form (rxPressureRowTranspose.H:15). The dominance clamp (`D=max(|D0|,sumOff)`) is quantitatively irrelevant: alpha dominates the diagonal by 2-3 orders (|alpha|₂≈1.46e10, RMS≈659). |
| H-F3 fixed-point / linearization system mismatch | **SURVIVES, narrowed** | The one element never externally validated is λ (J⁻ᵀsource) as the Jacobian of the phi-eliminated fixed-point map (convection-coefficient state sensitivity + relaxed-SIMPLE P←U tangent composition). This lives in the verified operator layer (J_PU/J_PP/J_U,phi). |

## New adjudication experiment (zero-compile: B4 diagnostic solver + B6 oracle)

`Log.verify_b014_b4b6.txt` (case snapshot `optProperties.b14_b4b6.txt`):

1. **RX-A fixed-state design FD**: `∂R_U/∂α` analytic vs FD relL2
   9.2e-8 → 3.0e-6 (eps 1e-3→3e-5), cos=1; `∂R_P/∂α` ~1e-6.
   **The R_x operator is externally correct.**
2. **StageB6 weighted D1/D2/D3 rerun**: the historical
   "|D_mom-prod|/|prod| ≈ 2e-9" validation is exposed as **λ-relative
   internal consistency** (the same λ on both sides). In today's rerun the
   diagnostic-path λ_U is degenerate (λ_U^T R_U,x d ≈ 7.4e-47) and the run
   aborts at the MMA non-finite gate — the diagnostic solver path is broken
   on the current HEAD (production path unaffected; BFINAL-010/011 fixed
   only the production solver). Historical anchors built on that path never
   validated λ externally.
3. **External-coverage table (final)**: source ✓ (BFINAL-013 P1, machine
   precision), R_x ✓ (this round, ~1e-7), design chain ✓ (gV end-to-end
   1e-7), **λ ✗ never externally validated** — BFINAL-006's T2/T3
   (tangent vs state FD) were explicitly deferred and never run; P8 proved
   only that J·w′=−R_x·d is exactly solvable (6e-14), not that w′ is the
   true state response.

## Why this cannot be fixed under this round's authorization

- The contraction formula is derived-correct given a correct λ (so no
  `sensitivity.H` edit is derivable — and an empirical rescale would be the
  forbidden "fitting 2.13");
- every candidate mechanism that survives adjudication lives in the
  verified operator layer (the J system's P←U / phi-elimination blocks);
- the task's stop rule for exactly this case applies.

## Recommended next round (cheapest first)

1. **Run T2/T3 at last** (no new mathematics): solve `J·w′ = −R_x·d`
   (D1/D2/D3) via the existing BFINAL-008 P8 machinery and compare w′
   against the state FD of the fully re-converged primal under x±h·d —
   one run yields the direction-dependent λ-coded/λ-true spectrum and
   localizes the defective J block.
2. Independent external anchor for the J_PU tangent (phiHbyA U/p-FD vs
   analytic tangent).
3. Separate small task: repair the B4 diagnostic solver path (currently
   degenerate λ) so λ-relative oracles become usable again.

## Provenance

- Repo `/home/ys/dsH/TO-ANISOTROPIC`, branch `agent/dsH-stage-b-validation`, base `5aaae31`
- No code changes this round (stop before implementation); static tests 8/8 unchanged
- Runs: `b14_b4b6` (B4+B6 oracle adjudication, MTO_RC=134 after the oracle
  block completed — abort at the MMA gate due to the degenerate
  diagnostic-path adjoint; all load-bearing oracle data printed before it)

## Files

| file | note |
|---|---|
| `DERIVATION.md` | phase-0 derivation with line numbers + hypothesis adjudication |
| `Log.verify_b014_b4b6.txt` | RX-A external FD + StageB6 exposure run |
| `optProperties.b14_b4b6.txt` | case snapshot |
