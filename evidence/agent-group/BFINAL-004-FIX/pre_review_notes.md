# BFINAL-004-FIX (replan v2) — Pre-Reviewer independent verification

- Role: Pre-Reviewer (independent). Decision: **APPROVED**.
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace, verified).
- HEAD `ca8a772b8339a36e7906ea32a7125061c47aa280`; `git status --short`: `M solveDiscreteFlowAdjoint.H`,
  `M solveDiscreteFlowAdjointProduction.H` (BFINAL-003 production diffs, untouched this round),
  untracked `src/stageB6RxDesignOracle.H`, `build/`, `evidence/*`. `git diff --stat` = 2 files (the two
  forbidden production headers) → plan's files_forbidden correctly excludes them.

## What I independently verified (not relying on the plan/summaries)

1. **B1 confirmed** — `src/stageB6RxDesignOracle.H` L552-553 compare unweighted
   `anRUa`/`anRPa` (per-cell `dR/dalpha`) against the deltaAlpha-weighted central-FD
   `fdRUa`/`fdRPa`. Raw log L1750-1759: `R_U,alpha relL2=1.906 cos=0.00066`, `R_P,alpha`
   `|a|=4.596e-10` vs `|b|=6.223e-10 cos=-1.4e-7` — both are the unweighted-vs-weighted artifact.
   Fix direction (new weighted fields `anRUaW`/`anRPaW`, keep `anRUa`/`anRPa` for RX-B/C reuse) is correct.

2. **B2 confirmed** — `filterTangent` (L705-737) constructs `stageB6_d`/`stageB6_y` with the
   `dimensionedScalar(...)` constructor (default `calculated` BC) and solves
   `fvm::laplacian(designFilterFaceMask, yfield)-fvm::Sp(b,yfield)+dfield*b`. Raw log L1775-1781:
   `calculatedFvPatchField::gradientInternalCoeffs` FATAL on `stageB6_y` inlet → SIGABRT (EXIT=134).
   The same calculated-BC defect also affects the solved replicas `gsensVolR` (L1009-1027) and
   `gsensR` (L1092-1111) (§9/§10) which crash once RX-C is reached. Production mirrors use
   `zeroGradientFvPatchScalarField::typeName` (createFields.H L416 xp, L593 gsensVol, L608 gsenshVol;
   filter_chainrule.H L30 dProjectionDeta; adjoint-filter solves filter_chainrule.H L189-213).
   Fix direction (zeroGradient on the solved fields + source-field mirrors, `correctBoundaryConditions`)
   is correct.

3. **B3 confirmed and independently re-derived** —
   - Wrong code at L380-389: `drAU = -alphaRel*rAUc*rAUc` and
     `dHbyA = rAUc*((alphaRel-1.0)*U - alphaRel*HbyA)`; wrong comments L34-36 / L380-382.
   - From OpenFOAM-7 semantics: `fvm::Sp(alpha,U)` adds `alpha*V` to the integrated diagonal D0
     (`fvmSup.C:118`); `relax()` → `D_rel = D0/alphaRel`, `S_rel = S0 + (D_rel-D0)*U`;
     `A() = D_rel/V` (`fvMatrix.C:751`); `H() = H0/V` with `H0 = S0 - offDiag*U + boundary`
     (alpha-independent). Therefore:
     `rAU = 1/A = V/D_rel`; `HbyA = H0/D_rel + (1-alphaRel)U`.
     → `drAU/dalpha = -rAU^2/alphaRel` (V cancels),
       `dHbyA/dalpha = (rAU/alphaRel)*((1-alphaRel)U - HbyA)` (V cancels).
   - Plan's B3 formula is exactly this; the "no extra V" conclusion is correct (the momentum-row
     anchor `dR_U/dalpha = V*U` carries V because `R_U = -residual()` is integrated, while the P-row
     `R_P = div(phi)` V-dependence cancels through `HbyA = rAU*H()`). v1's explicit `mesh.V()` was a
     category error and is correctly abandoned.
   - Wrong-code coefficient errors are consistent with the observed RX-B `R_P,xh` gap
     (raw log L1762-1770: `|a|=5.705e-6` vs `|b(FD)|=9.713e-5`, `relL2=1.012 cos=-0.176`, stable 5-eps plateau).

4. **FD side trusted (not self-referential)** — `rebuildResidual` (L240-310) is a full independent
   production rebuild (fresh UEqn → relax → residual, rAU → HbyA → phiHbyA → adjustPhi → rAtU →
   constrainPressure → pEqn → setReference → flux → div). The analytic side is a separate hand-coded
   linearization; no shared formula. FD facts already secured by the prior Post-Reviewer (independent
   recompute): `R_P,alpha ≈ 6.223e-10` and `R_P,xh ≈ 9.713e-5` stable 5-eps plateaus; `R_U,xh`
   `relL2=1.28e-10`. `R_w/J = CLOSED` (BFINAL-003 anchors byte-identical in cycle-1 log).

5. **Switch guard confirmed** — `solveDiscreteFlowAdjoint.H` L4069-4085: single
   `optProperties.lookupOrDefault<Switch>` + static once-guard + `discreteAdjointLabel=="pressureDrop"`
   gate; switch-off is a no-op. `build/bin/MTO_HF` exists (fresh 4691928 bytes, absolute-path target).

## Gate assessment

- SCOPE_GATE: PASS (DIAGNOSTIC_ONLY, only the diagnostic header + evidence + case optProperties; no production math).
- PHYSICS_GATE: PASS (B1/B2/B3 all physically/numerically correct; B3 derivation verified from source).
- ROOT_CAUSE_GATE: PASS (3 defects map 1:1 to the observed failures; `{}` patch count → ceiling not reached).
- CHANGE_SCOPE_GATE: PASS (1 source file, switch-guarded, frozen modules protected).
- VALIDATION_GATE: PASS (build+run+anchors+relL2/cos+Post-Reviewer independent recompute; stop conditions
  cover the residual-failure path → REPLAN_DIAGNOSTICS).

## Residual notes (non-blocking)

- `rebuildResidual` calls `pEqn.setReference(pRefCell,pRefValue)` unconditionally (L289), while production
  gates on `p.needReference()`; BFINAL-002 established these are no-ops and the rebuild is byte-identical
  to the validated stageB5 oracle, so this is a shared-oracle property, not an analytic-chain completeness
  defect for the alpha derivative. Whether the analytic R_P chain fully closes is empirically arbitrated by
  the run; the plan's stop condition (relL2>1e-3 → REPLAN_DIAGNOSTICS) covers the failure path.
- The FD-only facts (`R_P,x` NONZERO at both alpha and xh layers) survive even if the analytic chain needs
  further work, so the round's primary question (§1) is answerable regardless.
