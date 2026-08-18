# BFINAL-004-FIX cycle-3 — step2 audit-fix-probe (executor, independent)

- Stage: B-final · Mode: DIAGNOSTIC_ONLY · Step: step2-audit-fix-probe
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace, verified)
- HEAD: `ca8a772b8339a36e7906ea32a7125061c47aa280`; dirty tree: `M
  src/solveDiscreteFlowAdjoint.H` + `M src/solveDiscreteFlowAdjointProduction.H`
  (BFINAL-003 production diffs, LOCKED/forbidden), untracked probe header +
  build/ + evidence/ (unchanged from cycle-3/pre-state).
- Probe header pre-hash (cycle-3 baseline): `09832dd6ae5a3c7f37201b200bfb19dc0c48f42a3989f9593979ca73ca098baec`
  (identical to cycle-2 post-fix hash; the on-disk candidate already carries the
  B1/B2/B3 fixes — nothing is accepted from the self-reported cycle; every fix is
  re-audited below against OpenFOAM-7 source and the approved plan).

## Audit method

Read the full 1499-line probe header and re-derived every claimed formula from the
OpenFOAM-7 sources and the case files. No production file was modified; no source
edit was needed this cycle because the audit confirms all three fixes are present
and correct (see below). The empirical closure is arbitrated by the fresh build+run
(step3/step4); this note documents the source-level verification.

## B1 — RX-A analytic side deltaAlpha-weighting: VERIFIED CORRECT

- `anRUaW[3c+i] = anRUa[3c+i]*deltaAlpha[c]` (L675-682); RX-A compares
  `stageB6Metrics(anRUaW, fdRUa, ...)` and `stageB6Metrics(anRPaW, fdRPa, ...)`
  (L762-763), i.e. weighted-vs-weighted, matching the FD directional derivative
  `[R(alpha+eps*da)-R(alpha-eps*da)]/(2eps)`.
- `anRPaW` is NOT the pointwise product `anRPa*da`: it is
  `assembleWeightedRPa(deltaAlpha, ...)` (L686) which multiplies the weight into
  `drAU`/`dHbyA` BEFORE the face-flux assembly and the divergence (L563-569, L637-662):
  `dphiHbyA_w[f] = Sf&(wf*(dHbyA*da)_own+(1-wf)*(dHbyA*da)_nei)` and
  `dflux_w = flux(fvm::laplacian(drAU*da, p))`. This is the true directional
  derivative `J_P*da`; the pointwise `(J_P*1)*da` equals it only for uniform `da`
  and is the documented wrong form (header comment L665-674; approved plan states
  the weight must enter BEFORE the divergence). R_U is diagonal in alpha (`fvm::Sp`
  touches only the own diagonal), so `anRUaW = anRUa*da` per cell is exact.
- Unweighted `anRUa`/`anRPa`/`drAU`/`dHbyA` are retained and reused by RX-B/RX-C via
  `assembleWeightedRPa` with weights `dAlphaDxh*deltaXh` / `dAlphaDxh*z` (L784-795,
  L1064-1077). ✓ acceptance item "weighted comparison fields" satisfied.

## B2 — zeroGradient BC on the four solved fields: VERIFIED CORRECT

- `filterTangent` constructs `stageB6_d` (dfield) and `stageB6_y` (yfield) with
  `zeroGradientFvPatchScalarField::typeName` (L947-962) and calls
  `correctBoundaryConditions()` on both before the Helmholtz tangent solve
  (L967-968). The solve `fvm::laplacian(designFilterFaceMask, yfield)
  - fvm::Sp(b, yfield) + dfield*b` (L974-979) can no longer hit
  `calculatedFvPatchField::gradientInternalCoeffs` (SIGABRT of cycle-1):
  `gaussLaplacianScheme::fvmLaplacianUncorrected` calls
  `pvf.gradientInternalCoeffs()` on the UNKNOWN field's boundary patch
  (gaussLaplacianScheme.C L82-83); a `calculated` patch FATALs there.
- §9 replica `gsensVolR` (L1264-1271) and §10 replica `gsensR` (L1351-1358) are
  constructed with `zeroGradientFvPatchScalarField::typeName` +
  `correctBoundaryConditions()` before their adjoint-filter solves (L1272, L1359);
  their sources `gsenshVolR` (L1233-1240) and `gsenshR` (L1316-1323) carry the same
  BC. Production mirrors: createFields.H L416 (xp), L593 (gsensVol), L608
  (gsenshVol); filter_chainrule.H L30 (dProjectionDeta); the case field `0/xp` has
  `zeroGradient` on every patch. ✓ acceptance item "zeroGradient on all four solved
  fields" satisfied (stageB6_y, stageB6_d, gsensVolR, gsensR).
- Coefficient fields `drAUfield`/`drAUwField` (RX-A/RX-B/RX-C `dflux` assembly)
  keep the default (calculated, boundary=0) BC. This is intentional and exact for
  this case: the laplacian coefficient is only interpolated (never differentiated),
  and the boundary of `dflux = flux(fvm::laplacian(drAU,p))` is
  `p_bnd * d(boundaryCoeffs)/dalpha = -p_bnd*drAU_bnd*gradientBoundaryCoeffs(p)`.
  The case `0/p` is `fixedValue uniform 0` ONLY on `outlet` (p_bnd=0 => term=0) and
  `zeroGradient` elsewhere (gradientBoundaryCoeffs=0 => term=0), so `dflux` boundary
  = 0 reproduces the FD oracle's boundary-flux derivative exactly.

## B3 — drAU/dHbyA formulas and R_P per-term decomposition: VERIFIED CORRECT

Re-derived from OpenFOAM-7 source (independent of the plan's replan note):

1. `fvm::Sp(alpha,U)` adds `mesh.V()*alpha` to the INTEGRATED diagonal D0
   (`fvmSup.C`: `fvm.diag() += mesh.V()*sp.field()`), so `dD0/dalpha = V`.
2. `relax(alphaRel)`: `D /= alphaRel` and `S += (D-D0)*U` (fvMatrix.C L521 region:
   `D /= alpha;` `S += (D - D0)*psi_.primitiveField();`), so
   `D_rel = D0/alphaRel` and the residual is relaxation-invariant
   (`residual = S_rel - D_rel*U - offdiag*U - bnd = S0 - D0*U - offdiag*U - bnd`).
3. `A() = D()/mesh.V()` (fvMatrix.C L738: `primitiveFieldRef() = D()/psi_.mesh().V()`),
   hence `rAU = 1/A = V/D_rel` ALREADY carries V (units s).
4. `H() = (bndDiag*U + lduMatrix::H(U) + source_ + bndSource)/V` (fvMatrix.C L760),
   so with the relaxed source `rAU*H() = H0/D_rel + (1-alphaRel)*U`
   (`H0` = alpha-independent integrated part); the probe's `HbyABase =
   constrainHbyA(rAU*UEqn.H(), U, p)` (rebuildResidual L282) is exactly this in the
   interior.

Therefore, with U,p fixed:
- `drAU/dalpha = d(V/D_rel)/dalpha = -(V/D_rel)^2/alphaRel = -rAU^2/alphaRel`
  (V cancels) — code L418: `drAU[celli] = -rAUc*rAUc/alphaRel;` ✓
- `dHbyA/dalpha = d(H0/D_rel)/dalpha = -(H0/D_rel)*(V/D_rel)/alphaRel
  = (rAU/alphaRel)*((1-alphaRel)*U - HbyA)` (no extra V) — code L419-424 ✓

R_P per-term decomposition: `assembleWeightedRPa` accumulates
`anRPaHbyAwOut = div(dphiHbyA_w)` and `anRPaFluxwOut = -div(dflux_w)` through the
SAME face loops as the total `anRPwOut = div(dphiHbyA_w - dflux_w)` (L637-662), so
`col1 + col2 == col3` per cell EXACTLY (machine zero); `stageB6_rpa_terms.mtx`
exports `[N div(dphiHbyA_w) ; N -div(dflux_w) ; N total_w]` (L720-723), i.e. the
two weighted terms sum to the total, which the run compares against the
independent central-FD oracle. ✓ acceptance item "correct drAU/dHbyA" satisfied.

The remaining chain members that the FD rebuild includes but the analytic side
does not model explicitly — `adjustPhi`, `constrainPressure` (no fixedFluxPressure
patches in this case => no-op), `pEqn.setReference` (affects only the solve, not
`pEqn.flux()`), boundary-flux rebuild — are empirically arbitrated by the run
(relL2<1e-3/cos>0.999 => CLOSED; else REPLAN_DIAGNOSTICS per plan stop condition).
For this case the outlet `phiHbyA` derivative is captured by the assignable-patch
extrapolation (uAssignable => `dphiHbyA_b = pSf & dHbyA_cell`) and `dflux_bnd=0`
is exact (p=0 outlet), so closure is expected; the run decides.

## OFstream flush/close audit: VERIFIED (21 streams, all flushed + scoped)

| section | streams | flush | scope-close |
|---|---|---|---|
| RX-A | rxaFD, rxaAn, rpaTerms (L710-712) | L724, L767 | block L709-769 |
| RX-B | rxbFD, rxbAn (L798-799) | L803, L874 | block L797-876 |
| RX-C | rxcFD, rxcAn, rxcZ, rxcXhFD, rxcXhAn, rxcAlphFD (L1084-1089) | L1102, L1222 | block L1083-1225 |
| §9 vol | vOut (L1305) | L1308 | block L1232-1310 |
| §10 | ghOut, gOut (L1376-1377) | L1381 | block L1375-1382 |
| §9/§10 | lamOut (L1390) | L1400 | block L1388-1495 |
| §9/§10 | daOut, dxOut, dadxOut, dirOut, ctOut, prOut (L1402-1407) | L1430-1431 | same block |

All OFstreams are function-local objects in scoped blocks (destructor closes the
file) AND flush()ed after each write batch, so no `.mtx` tail can be truncated if a
later gate aborts. ✓

## Grep acceptance checks

- `zeroGradientFvPatchScalarField::typeName` present on: stageB6_d (L953),
  stageB6_y (L961), gsenshVolR (L1239), gsensVolR (L1270), gsenshR (L1322),
  gsensR (L1357) — 6 sites, all 4 solved fields + 2 source mirrors.
- `assembleWeightedRPa` call sites: L686 (RX-A, w=deltaAlpha), L795 (RX-B,
  w=dAlphaDxh*deltaXh), L1077 (RX-C, w=dAlphaDxh*z) — the weighted comparison
  fields.
- `drAU[celli] = -rAUc*rAUc/alphaRel` (L418) and
  `dHbyA[celli] = (rAUc/alphaRel)*((1.0 - alphaRel)*U[celli] - HbyABase[celli])`
  (L419-424) — the correct no-V formulas.
- No `Foam::max(projEtaDenom,` (D4 raw-division guard `safeEtaDivide` L931-941).

## Scope integrity

`git status --porcelain=v1` unchanged vs cycle-3/pre-state (only the two
BFINAL-003 production heads modified, hashes identical to pre-state record:
solveDiscreteFlowAdjoint.H=f0c81497..., solveDiscreteFlowAdjointProduction.H=
8c4901bf...). No production file touched; only this evidence file was written.

## Outcome

step2 acceptance MET: B1/B2/B3 edits present and source-verified against
OpenFOAM-7 semantics; grep confirms zeroGradient on all four solved fields,
weighted comparison fields, and correct drAU/dHbyA; every OFstream has a matching
flush() and scope-close; no production file touched. Proceeding to step3
(header-only incremental wmake) + step4 (run on /home/ys/dsH/b2_case_smoke) +
step5 (FINAL_REPORT).
