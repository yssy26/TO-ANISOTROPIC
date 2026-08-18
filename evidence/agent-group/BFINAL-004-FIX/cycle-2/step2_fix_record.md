# BFINAL-004-FIX (replan v4) cycle-2 — step2-fix-probe record

- Stage: B-final · Mode: DIAGNOSTIC_ONLY · Step: step2-fix-probe (executor)
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace) @ `ca8a772b8339a36e7906ea32a7125061c47aa280`
- File changed: `src/stageB6RxDesignOracle.H` ONLY (untracked diagnostic header; 1280 -> 1499 lines)
- Pre-edit sha256 (cycle-2/pre-state capture): `51a02fc6288da612fc7e9b775cfd05ac757c137b492d5e863185f53c74c4b31d`
- Post-edit sha256: `09832dd6ae5a3c7f37201b200bfb19dc0c48f42a3989f9593979ca73ca098baec`
- Pre-fix byte copy: `cycle-2/pre-state/stageB6RxDesignOracle.H.pre-fix-copy`

## What was changed (exactly the approved plan, hypothesis BFINAL-004-RP-PREASSEMBLY-WEIGHTING)

### (a) Weighting order — J_P*w, NOT (J_P*1)*w (the primary fix)

New helper lambda `assembleWeightedRPa(const scalarField& w, scalarField& anRPwOut,
scalarField& anRPaHbyAwOut, scalarField& anRPaFluxwOut)` (inserted after the unweighted
R_P,alpha chain, before the RX-A block):

- `drAUw[celli] = drAU[celli]*w[celli]`, `dHbyAw[celli] = dHbyA[celli]*w[celli]`
  (drAU/dHbyA = the already-verified unweighted per-cell derivatives, B3 formulas).
- `dphiHbyAw[face] = Sf&(wf*dHbyAw[own] + (1-wf)*dHbyAw[nei])`; boundary faces:
  `uAssignable ? dHbyAw[celli] : 0` (identical boundary semantics to the unweighted chain).
- `dflux_w = flux(fvm::laplacian(drAUwField, p))` with `drAUwField` constructed exactly
  like the existing `drAUfield` (dimensionedScalar + correctBoundaryConditions), materialized
  over ALL faces (internal + boundary) as in the unweighted chain.
- Divergence with per-term decomposition: `anRPwOut = div(dphiHbyAw - dfluxWFaces)`,
  `anRPaHbyAwOut = div(dphiHbyAw)`, `anRPaFluxwOut = -div(dfluxWFaces)`,
  so `anRPaHbyAwOut + anRPaFluxwOut == anRPwOut` per cell (col1+col2 == col3).

Call sites (replacing the three pointwise post-multiplies `anRPa*da` / `anRPa*f`):

- RX-A: `assembleWeightedRPa(deltaAlpha, anRPaW, anRPaHbyAW, anRPaFluxW)` (was L521
  `anRPaW[celli] = anRPa[celli]*da`).
- RX-B: `assembleWeightedRPa(wXh, anRPxh, ...)` with `wXh = dAlphaDxh*deltaXh` (was L599
  `anRPxh[celli] = anRPa[celli]*f`).
- RX-C `analyticRxDirection`: `assembleWeightedRPa(w, anRPd, ...)` with
  `w = dAlphaDxh*z` (was L857 `anRPd[celli] = anRPa[celli]*f`).

R_U side UNCHANGED (per-cell product `anRUa*w`; R_U is diagonal in alpha — Sp touches only
the own diagonal — already closed to relL2=1.28e-10).

Export change (RX-A): `stageB6_rxa_analytic.mtx` now writes the WEIGHTED fields
`anRUaW` / `anRPaW` (directly comparable to `rxa_FD_all_eps`, no further deltaAlpha
multiplication); `stageB6_rpa_terms.mtx` now writes the WEIGHTED per-term columns
`[N div(dphiHbyA_w) ; N -div(dflux_w) ; N total_w]` (col1+col2 == col3 == FD oracle).

### (b) D4 — sign-preserving eta-denominator division (mirror of production)

`Foam::max(projEtaDenom, SMALL)` at the three sites (projectionTangent etaResp, §9
volEtaCorr, §10 pressEtaCorr) replaced by `safeEtaDivide(num)`:

```cpp
auto safeEtaDivide = [&](const scalar num) -> scalar {
    if (mag(projEtaDenom) <= SMALL) { Info<< "...singular..." ; return 0.0; }
    return num/projEtaDenom;
};
```

Mirrors production filter_chainrule.H L81/L108-111 (mag()<=SMALL singularity guard + raw
division); a negative denominator (here projEtaDenom=-3.1386e-6) is now divided directly
instead of being clamped to +SMALL=1e-15 (the ~3.14e9 sign-flipping amplification).

### (c) stageB6Metrics display clamp

`rel = sqrt(nD/nB)` / `cos = dot/sqrt(nA*nB)` computed with EXACT-zero degenerate guard
only (`nB2<=0 || nA2*nB2<=0` -> print "degenerate(nB=0)" with raw |a|,|b|,|err|); no SMALL
clamping of the denominator, so the printed relL2/cos is geometrically consistent with the
printed |a|,|b|,|err| (cycle-1 RX-A R_P,alpha printed 0.0422/-1.75e-5 while |a-b|/|b|=1.0115).

### (d) Header comment (L33-41 area) + inline comments

Header RX-A bullet now documents the directional J_P*w chain (weight inside the face-flux
assembly), the NON-diagonal R_P vs DIAGONAL R_U asymmetry, and the weighted per-term export.
RX-B/RX-C section comments and the B1 inline comment updated accordingly. OFstream
per-section flush/close retained unchanged.

## Acceptance-criteria verification (step2)

| criterion | result |
|---|---|
| `grep -n 'Foam::max(projEtaDenom' src/stageB6RxDesignOracle.H` empty | PASS (EMPTY) |
| `grep -nE 'anRPa\[celli\]\*da\|anRPa\[celli\]\*f\|anRPd\[celli\] *= *anRPa'` empty | PASS (EMPTY) |
| `assembleWeightedRPa` defined and called at RX-A/B/C | PASS (L557 def; L686 RX-A, L795 RX-B, L1077 RX-C) |
| only this header has a (new) diff | PASS — production headers sha256 identical to cycle-2/pre-state (f0c81497... / 8c4901bf...); git diff --stat still only the two pre-existing BFINAL-003 M files |
| rpa_terms weighted columns col1+col2==col3==total | PASS by construction (per-cell += of the same dphiHbyAw/dfluxWFaces split); verified at runtime by mtx in step4 |
| stageB6Metrics no SMALL-clamped relL2/cos | PASS (exact-zero guard only) |
| brace balance | 169 open / 169 close |

## Scope

- No production file touched (solveDiscreteFlowAdjoint.H, solveDiscreteFlowAdjointProduction.H
  untouched — sha256 unchanged; BFINAL-003 diffs preserved).
- No filter/projection/conductivity/MMA/threshold changes; no SIMPLE-transpose; no MMA.
- Compile sanity check: header-only incremental wmake (FOAM_USER_APPBIN=build/bin, unset
  FOAM_SIGFPE) run in this step; WMAKE_EXIT recorded in `step2_build_check.log` if run here,
  otherwise completed by step3-build.
