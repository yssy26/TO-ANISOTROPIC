# BFINAL-012 cycle-1 FINAL REPORT — current-source end-to-end FD amplitude gate

## Verdict

**FAIL (checkpoint rule: amplitude ratio >= 2x on gDP; sign flips on J).**
Executed on the current source `2aaf649` (+ the D1-ladder extension committed
with this report), the current validation case, with BOTH adjoint labels
converged to <= 1e-9 at every adjoint solve in the run. Per the task's stop
rule the gate stopped here; no parameter was tuned, no threshold relaxed, no
re-tuning attempts were made.

```
FROZEN_GRADIENT_STATUS = FAIL   (in-run Stage-B3 verdict)
formalConvergenceOK    = 1      (every +/- FD primal point converged, windowed)
plateauOK              = 1      (all FD ladders show stable plateaus)
```

## Provenance

- Repo `/home/ys/dsH/TO-ANISOTROPIC`, branch `agent/dsH-stage-b-validation`, base `2aaf649`
- Binary rebuilt 2026-08-19 21:54 (wmake exit 0, `wmake_bfinal012.log` at repo root); only code change = D1 diagnostic ladder extended to {1e-5 … 1e-2} (four decades, task spec); formal acceptance steps unchanged {1e-4, 3e-4, 1e-3}
- Case `/home/ys/dsH/b12_fdgate` (= `cp -a b11_diag`, same 33600-cell / 104740-face mesh as all B-final runs and as the historical stage-B2/B3 case), `stageB2Enabled true`, all probes off, `mmaUpdateEnabled false`, production adjoint `pressureGAMG`, GAMG default
- Run A: `Log.verify_b012_A.txt` (sha256 `f2e237003cba3c008d3c9998f46ee6eb2ce52d19bf4a4d074f26335532fea506`), MTO_RC=0
- Method: existing Stage-B2/B3 machinery (`validateStageB2GradientAmplitude.H`): raw-x perturbation x±h·Dk → full `recomputeDesignChain()` (filter → projection → material) → windowed-converged frozen-primal re-solve per point → objective re-evaluation; adjoint projections from the PRODUCTION adjoints; deterministic physical directions D1/D2/D3

## Adjoint convergence gate (task requirement) — PASS

Every adjoint solve in the run, true relative residual <= 1e-9, zero
"did not converge" warnings:

```
main loop   thermalCoupling: 750 iters  9.4956818492e-10   (bit-identical to b11 E0')
main loop   pressureDrop:    902 iters  9.92818139401e-10  (bit-identical to b11 E0')
B2 baseline thermalCoupling: 711 iters  9.17308458009e-10
B2 baseline pressureDrop:    679 iters  9.80164486971e-10
```

## Baseline (post full-SST + refreeze)

```
J = -1.16745040317   PD = 48353.7815132 Pa   gDP = 0.934151260527   gV = -7.887e-09
(PDmax = 25000 Pa exactly; historical stage-B2/B3 baseline on the same mesh family:
 J = -1.16736, PD = 48336.3, gDP = 0.93345 — same state family)
```

## Amplitude table (formal steps; full ladder in stage_b2_fd_scan.tsv)

| dir | metric | FD (h=1e-3) | ADJ | ratio ADJ/FD | sign | in-run rel |
|---|---|---|---|---|---|---|
| D1 | J    | +0.0428559 | −0.0159322 | **−0.372** | **FLIP** | 3.69 |
| D2 | J    | −0.0034689 | −0.0241647 | **+6.97** | same | 0.856 |
| D3 | J    | −0.0002827 | +0.0031606 | **−11.18** | **FLIP** | 1.089 |
| D1 | gDP  | −2.53746 | −5.40506 | **2.130** | same | 0.531 |
| D2 | gDP  | +0.523409 | +1.10593 | **2.113** | same | 0.527 |
| D3 | gDP  | +6.47722 | +14.0447 | **2.169** | same | 0.539 |
| D1 | gV   | 0.155518 | 0.155519 | 1.0000 | same | 1.58e-6 |
| D2 | gV   | −0.0844964 | −0.0844967 | 1.0000 | same | 3.60e-6 |
| D3 | gV   | 0.0639628 | 0.0639624 | 1.0000 | same | 6.83e-6 |

FD quality: D1 ladder spans 1e-5 → 1e-2; FD_gDP = −2.5363/−2.5363/−2.5371/−2.5371/−2.5375/−2.5354/−2.5181 (spread < 0.15% over the interior; 0.75% including h=1e-2); FD_J plateau 0.0425–0.0437 (< 3%). Plateaus are unambiguous.

Pressure identity: dgDP/dh vs (1/PDmax)·dPD/dh agree to 2.2e-14 (FD-side
normalization consistent). Volume chain FD vs <gsensVol,Dk> matches the
probeFD rows (flow solve does not affect gV).

## Failure localization (facts from this run)

1. **Shared design chain exonerated.** gV (volume fraction) projections match FD with correct signs on all three directions; rel error 1.4e-7 at h=1e-2, noise-floor-limited (∝1/h: 3.9e-5 at 1e-4 → 1.4e-7 at 1e-2). The filter/projection/material chain-rule used by ALL gradients is correct. (The in-code formal gV criterion 1e-6 is marginally missed at formal steps 1.6e-6–6.8e-6 — FD noise floor, not an adjoint error; far tighter than the flow terms as the task expects.)
2. **gDP: a single systematic scalar.** Sign correct on all directions, FD plateaus rock-solid, but ADJ = (2.113 … 2.169)×FD (mean ≈ 2.137, ±1.6%). This localizes to the pressureDrop adjoint→gradient assembly (`sensitivity.H` PD column / its R_x contraction into dgdx[1]) — NOT the solver (adjoints converged ≤1e-9), NOT the matrix/operator (BFINAL-010 gate 1e-16), NOT the rhs constructions (BFINAL-005 anchors).
3. **J: structural, not scalar.** Sign flips on D1 and D3 plus 6.7–11.2× amplification on D2/D3 and −0.37× on D1. A thermal-objective adjoint-chain defect (AdjHeatTransfer / thermalCoupling chain or its contraction into dfdx), not explainable by one scale factor.
4. **Historical continuity.** The pre-correction stage-B2/B3 (commit f3c4aa2, same mesh family) failed with gDP factors 2.15/2.33/2.17 and J D1 sign-flip / D2 5.17 / D3 −20 — essentially the SAME pattern and magnitudes. The BFINAL-002…008 solve-layer corrections and the BFINAL-010/011 convergence fixes were each verified against solve-layer references (J·v vs residual FD, oracle contractions, dot tests) and remain valid at that layer; the end-to-end objective-FD discrepancy they were expected to remove is still present. The end-to-end FD gate is the binding authority, exactly as the task book states.
5. GRADPROXY plateaus (TC +2217 / PD −205194) demonstrate adjoint solve stability only and do not substitute for this gate (task wording).

## Candidate hypotheses for the next task (labeled, unverified)

- **H-PD-scale**: a ≈2.13 scalar bookkeeping error in the dgdx[1] assembly (double-counted term or wrong scale; note PDmax = 25000 exactly and 1+gDPbase = 1.934 — neither matches 2.13, so it is not a simple gDP normalization slip).
- **H-J-struct**: thermal adjoint chain defect — objective sign/scale convention in the dfdx thermal term, or a missing/incorrect contribution path (e.g., frozen-hot Q coupling) that differs between the adjoint linearization and the FD-re-evaluated objective.
- Suggested first probes for the follow-up task: per-term decomposition of dfdx·D / dgdx[1]·D (gsens components) vs FD; rerun the stageB6 R_x design oracle on the current binary; a small-mesh reproduction to iterate quickly.

## Repeatability (task requirement) — resolved without a second run

Skipped per the stop rule, with justification recorded: the failure is
deterministic and systematic — (a) main-loop adjoints are bit-identical to
the b11 E0' run (750/902 iterations, identical residuals), (b) FD plateaus
spread < 0.15% across three decades on D1, (c) the ADJ/FD ratio is a
constant ~2.13 across three independent directions. Run-to-run noise cannot
produce any of these; a repeat would reproduce the failure bit-for-bit.

## Files

| file | note |
|---|---|
| `Log.verify_b012_A.txt` | full run (baseline, projections, all FD rows, verdict) |
| `stage_b2_fd_scan.tsv` | machine-readable FD scan (all directions/steps) |
| `stage_b2_summary.tsv` | machine-readable verdict/projections |
| `optProperties.b12_fdgate.txt` | case controls snapshot |
| `../../../src/validateStageB2GradientAmplitude.H` | D1 ladder extension (only code change) |
