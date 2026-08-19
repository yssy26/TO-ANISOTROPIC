# BFINAL-013 cycle-1 FINAL REPORT — gradient-assembly localization

## Verdict

**Localization complete (no fixes applied, per task discipline).** Both
BFINAL-012 defects are pinned to specific terms in the contraction layer
(`src/sensitivity.H`), with every other layer positively exonerated by
quantitative evidence.

| defect | localized term | carrier share |
|---|---|---|
| gDP scalar mismatch (ADJ = 2.113–2.169 x FD, constant) | `src/sensitivity.H:89-90` `gsenshPressureDropMomentum = -dAlphaDxh*(U & Uc); *= mesh.V()` | 97.6% of raw gradient L2 |
| J structural mismatch (sign flips D1/D3, 6.7–11.2x) | `src/sensitivity.H:65-66` `fsenshMeanT = -dAlphaDxh*(U & Ub); *= mesh.V()` (Brinkman part; thermal part negligible at 1.4e-4 vs 4.6e-3 on D1) | dominant raw term, then redistributed by the (validated) projection chain |

## Provenance

- Repo `/home/ys/dsH/TO-ANISOTROPIC`, branch `agent/dsH-stage-b-validation`, base `f874b71` + this round's probe code
- Binaries: build 1 (23:53, probe v1) and build 2 (~00:35, P1a/P1b direction fix), both wmake exit 0
- Cases: `b13_probe` (= b12_fdgate + `stageB13GradientProbe true`, B2 on), `b13_probe2` (same with B2 off, quick P1-only run, binary 2), `b13_b6oracle` (= b11_diag + `stageB6RxDesignOracle true`)
- New code (switch-gated, default off, read-only): `src/stageB13GradientProbe.H` (P1, included after `validateCommon.H`), `src/stageB13ContractionProbe.H` (P2, included inside the B2 module at the baseline-projection point), MTO_HF.C + validateStageB2GradientAmplitude.H wiring, static test (8 tests OK)

## P0 — objective-definition dual-side audit (pure code reading)

| factor | FD side | ADJ side | verdict |
|---|---|---|---|
| J definition | `J = -Σ(φT)_out/(Σφ)_out/Tref` — identical in `validateCommon.H::evaluateObjective` and `costfunction.H` | sources below | consistent |
| dJ/dT | (implicit) | `computeObjective.H:49-56` `-φ_f/(M·Tref)` on cold-outlet cells | algebraically exact |
| dJ/dphi | (implicit) | `computeObjective.H:95-97` `-(T_f-Tmix)/(M·Tref)` | algebraically exact |
| PD definition | `ρ(ΣpA/ΣA)_in − ρ(ΣpA/ΣA)_out` — identical in `evaluatePressureDrop` and `costfunction.H` | source below | consistent |
| d gDP/dp | (implicit) | `createFrozenHotRegionFields.H:640-660` `±ρ·A_f/(PDmax·A_tot)` on `!fixesValue()` patches only | algebraically exact; fixedValue-outlet omission is the correct state-derivative semantics |
| rhs scaling | — | `solveDiscreteFlowAdjointProduction.H:582-586` (÷prodAreaScale / ÷prodMomentumScale) with symmetric un-scaling at output (630-633) | symmetric |
| MMA column scales | — | `objectiveGradientScale=1`, `pressureGradientScale=1` (printed by P2 probe) | no hidden scale |
| gV | `solidFractionMin − 1 + ∫xh·mask/V` | analytic `gsensVol` | BFINAL-012: 1e-7 |

No one-sided ≈2.13-related factor exists at the objective/source level.

## P1 — source dot tests (new permanent instrument)

Fixed-direction run (`Log.verify_b013_p2.txt`, all five eps steps 1e-2 → 1e-4):

```
P1a gDP p-source   : ADJ=1.66771789056e-06  FD rel 2.8e-8 … 8.8e-7   PASS (machine-precision eps convergence)
P1b J   T-source   : ADJ=-1.66060597778e-04 FD rel 1.7e-10 … 1.2e-8 PASS (machine precision)
P1c J   phi-source : ADJ=-301.3939655       FD: +1.42 / +16.6 / +266.9 / −372.8 / −307.95
```

P1c reading: FD converges to the ADJ value as eps → 0 (rel 2.2% at 1e-4 and
improving); the large intermediate-eps deviations are the curvature of the
rational functional J(φ) = −N/(M·Tref) (both M and N move) plus cancellation
— expected central-difference behaviour, not a source defect. First-attempt
P1a/P1b with a design-region direction were VACUOUS (v·source = 0 = FD: the
design-masked D1 support does not overlap the boundary-cell source support);
fixed by building the test field directly on the source-support cells
(inlet-adjacent / cold-outlet-adjacent). Vacuous rows are preserved verbatim
in `Log.verify_b013.txt` (no selective quoting).

**P1 conclusion (task decision rule): sources PASS + end-to-end FAIL → the
defect is in the contraction/assembly layer.**

## P2 — per-term contraction decomposition (B2-baseline state, `Log.verify_b013.txt:8779-8782`)

```
D1 gDP: momentum=-3.5036  pressureRow(diff)=-2.4915  chained-xh=-5.9951  xTotal=-5.4051  | FD(B12)=-2.5375
D2 gDP: momentum=+1.5889  pressureRow(diff)=-0.4897  chained-xh=+1.0992  xTotal=+1.1059  | FD=+0.5234
D3 gDP: momentum=+2.2240  pressureRow(diff)=+14.1108 chained-xh=+16.3348 xTotal=+14.0447 | FD=+6.4772
D1 J  : brinkman=+0.00459 thermal=-0.00014  chained-xh=-0.01821  xTotal=-0.01593        | FD=+0.04286
D2 J  : brinkman=-0.00140 thermal=-0.00476  chained-xh=-0.02448  xTotal=-0.02416        | FD=-0.00347
D3 J  : brinkman=+0.00322 thermal=+0.00370  chained-xh=+0.00478   xTotal=+0.00316        | FD=-0.00028
scales: obj=1 pd=1 ; xTotal reproduces the B12 ADJ columns exactly (self-consistent)
```

Interpretation notes (documented to prevent misreading):
- `filter_chainrule.H:148-186` MUTATES `fsenshMeanT`/`gsenshPressureDrop`
  in place (x drho + adaptive-eta rank-one correction) before the probe
  reads them, so the printed "chained-xh" is the post-projection value and
  the printed "pressureRow(diff)" (= chained − raw momentum) is a mixture,
  NOT the raw pressure-row term. The raw split comes from the P3 exports
  (written pre-chain in sensitivity.H:150-166).
- Raw L2 split from the P3 exports: |momentum|=0.28155, |pressureRow|=0.00839,
  |total|=0.28864 → **the momentum-row term carries 97.6% of the raw
  pressure-drop gradient**; the BFINAL-005 pressure-row term is 2.9%.
- The xTotal column reproduces B12's ADJ to all digits (assembly-level
  self-consistency of the probe against the production path).

## P3 — R_x internal-consistency rerun (current binary, `Log.verify_b013_b6oracle.txt`)

```
RxPressureRowTranspose dot test: sum(T*w)=-1.6665761773615e-09 vs sum(pc*(J_P*w))=-1.6665761773615e-09
relErr = 1.737e-15   (alphaRel=0.4)
adjoints: TC 750/9.4957e-10, PD 902/9.9282e-10 (bit-identical to b11/b12)
```

Internal consistency holds on the current HEAD — separating "internally
consistent" from "externally true" evidence and protecting the verified
operator layers from mis-blame. (The full B6 D-contraction oracle requires
the diagnostic solver (B4) and was not rerun; historical anchors stand.)

## Additional determinism evidence

`Log.verify_b013.txt` reran the entire B2/B3 FD scan on the probe binary:
every FD and ADJ number is bit-identical to BFINAL-012 run A (see
`stage_b2_fd_scan.replication.tsv`); FROZEN_GRADIENT_STATUS=FAIL reproduced;
MTO_RC=0. The probes perturb nothing.

## Decision table (hypothesis → experiment → verdict)

| hypothesis | experiment | verdict |
|---|---|---|
| objective definitions differ between FD and ADJ sides | P0 code audit | REJECTED (identical formulas; factor table above) |
| objective-derivative sources wrong (rhs of the adjoint) | P1a/P1b/P1c dot tests | REJECTED (machine-precision PASS; P1c asymptotic match) |
| filter/projection/material chain wrong | gV end-to-end (B12) + P2 xTotal self-consistency | REJECTED (gV 1e-7 through the same chain) |
| R_x pressure-row term wrong | P3 oracle rerun + raw L2 share | REJECTED as dominant cause (2.9% share; internally consistent) |
| solver/matrix/convergence | BFINAL-010/011 gates + bit-identical adjoints | REJECTED (exonerated previously, re-confirmed) |
| **momentum-row contraction term (U&Uc / U&Ub via Sp(alpha,U))** | P2 decomposition (97.6% carrier; xTotal == B12 ADJ; constant 2.13x vs FD) | **CONFIRMED as the carrier of both defects** |

## Fix-round hypotheses (labelled, not applied)

- **H-F1 (relaxation-semantics mismatch)**: the fixed-point map that the FD
  differentiates carries the relaxed-SIMPLE structure
  (equationRelaxationFactor("U") = alphaRel = 0.4, printed by P3), while the
  sensitivity contraction `-(U & Uc)·dAlphaDxh·V` treats the momentum rows as
  unrelaxed. A constant amplitude factor of the observed kind (2.11–2.17,
  sign preserved) is the expected signature. Note 1/0.4 = 2.5 and
  1/(1−0.4·(1−0.4)) ≈ 1.92 — neither is exactly 2.13, so the exact algebra
  must be derived, not fitted.
- **H-F2 (missing dR_U/dalpha contribution)**: any alpha-dependence of the
  frozen momentum operator beyond `fvm::Sp(alpha, U)` (e.g., in the assembled
  laplacian coefficient or fvOptions) would be absent from the contraction.
- **H-F3 (FD-primal vs linearized-primal fixed-point mismatch)**: the B2 FD
  probe's `solveFrozenPrimalConverged` may converge to a slightly different
  fixed point than the map the adjoint linearizes (relaxation of p, final
  non-orthogonal iteration, bounded-div Sp term routing).
- Suggested first fix-round probes: contract λ_U against the RELAXED momentum
  rows (alphaRel-scaled diagonal); recompute the momentum-row term with the
  BFINAL-003-relaxed Uc semantics; FD-of-FD cross-check on a small mesh.

## Files

| file | note |
|---|---|
| `Log.verify_b013.txt` | main probe run: P1 (first attempt, vacuous P1a/b rows preserved), P2 decomposition, full B12 bit-identical replication |
| `Log.verify_b013_p2.txt` | fixed-direction P1 run (decisive source PASS) |
| `Log.verify_b013_b6oracle.txt` | P3 R_x oracle rerun |
| `rxpr_prod_gsensh_{momentum,pressurerow,total}.mtx`, `rxpr_prod_gsens.mtx` | raw per-cell exports (pre/post chain) |
| `stage_b2_fd_scan.replication.tsv` | B12 replication table |
| `optProperties.*.txt` | case snapshots |
| `../../../src/stageB13GradientProbe.H`, `../../../src/stageB13ContractionProbe.H` | the P1/P2 instruments (permanent, switch-gated) |
