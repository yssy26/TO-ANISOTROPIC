# BFINAL-008 S2 — stageB8JPPActualResidualFD probe: written, built, smoke-verified

- Stage: B-final · Mode: PATCH · Step: S2 (probe write + build + probe smoke)
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace, verified) · HEAD
  `ca8a772b8339a36e7906ea32a7125061c47aa280` · branch `agent/dsH-stage-b-validation`
- Case (scratch, switch ON): `/home/ys/dsH/b8_scratch_on` (byte-identical copy of
  `/home/ys/dsH/b2_case_smoke` with only `discreteExplicitMatrixFile` re-pointed and
  `stageB8JPPActualResidualFD true;` appended — same pattern as BFINAL-007's `b7_scratch_on`)

## 0. Files (this step)

- `src/stageB8JPPActualResidualFD.H` (NEW, 539 lines): switch-guarded
  (`stageB8JPPActualResidualFD`, default false), read-only actual-residual FD probe.
- `src/solveDiscreteFlowAdjoint.H`: +25-line guarded include block at END of header
  (stageB5/stageB6 pattern; static guard; NO MTO_HF.C change).
- No other file touched (verified hashes: `solveDiscreteFlowAdjointProduction.H`,
  `MTO_HF.C`, `NS.H`, `sensitivity.H`, `rxPressureRowTranspose.H` all unchanged).

## 1. Acceptance

- `WMAKE_EXIT=0` (3 builds: S1 binary, probe v1, probe v2 after metric fix).
- Probe switch default-false no-op: structural (entire probe body inside the
  `lookupOrDefault<Switch>("stageB8JPPActualResidualFD", false)` guard; absent/false
  => one dictionary lookup only, same proven pattern as stageB5/stageB6/stageB7).
  Byte-identical run optional per acceptance.
- Probe reads `wPrime_TAN_D1.mtx` (134400 values): `||w||=2.08447046292e+18`
  (matches BFINAL-006/007 exactly), `||n_P||=0.999978325577`, `n_P[PREF=cell0]=0`.
- Probe writes ONLY under `evidence/agent-group/BFINAL-008/cycle-1/artifacts/`
  (23 files; none in the case dir).
- MTO_RC=0 both runs; probe executes in both adjoint blocks (thermalCoupling +
  pressureDrop) with byte-identical metrics — double-run matches stageB5's
  established pattern and doubles as an internal repeatability check.

## 2. Probe semantics (per approved plan S2)

- Rebuilds the ACTUAL frozen-primal residual pair R(w) = [R_U; R_P] from production
  semantics:
  - U rows: full NS.H UEqn (`fvm::div(phi,U) - fvm::laplacian(nuEffFrozen,U)
    + fvm::Sp(alpha,U) == -fvc::grad(p) + fvOptions(U)`, +dev2 when ransFlowModel,
    relax, constrain), `resU = -UEqn.residual()` (NS.H sign, stageB5-validated).
  - P rows: `rAU=1/A() -> HbyA=constrainHbyA -> phiHbyA=fvc::flux -> adjustPhi ->
    rAtU(consistent branch) -> constrainPressure -> pEqn=fvm::laplacian(rAtU,p) ->
    setReference -> phi = phiHbyA - pEqn.flux()`, `resP = div(phi)` raw face sums
    (owner +, neighbour -, boundary +, NO 1/V) — production `flux()` semantics incl.
    the fixedValue outlet internalCoeffs boundary contribution.
  - **phiHbyA is p-independent (rebuilt from PRODUCTION p); only pEqn.flux() responds
    to dp** — this is the BFINAL-007 Q1-verified semantics the corrected J implements.
    (Initial probe version rebuilt phiHbyA from pmod, injecting a spurious
    d(phiHbyA)/dp term: its FD mismatched the verified FD at relL2=0.37. Fixed.)
  - U-row residual responds to pmod through `-fvc::grad(pmod)` (J's U-row dp
    contribution); pEqn built on psi=pmod so `flux()` carries the P-P block.
- P1: v = [0; n_P] (dp-only); P3: dir1=n_P, dir2=dp-heavy outlet-supported,
  dir3=seeded smooth dU+dp; central FD [R(w+eps v)-R(w-eps v)]/(2eps), eps in
  {1e-2,1e-3,1e-4,1e-5}; block metrics relL2/cos/normRatio/maxDiff per U rows /
  P-internal / P-outlet-adjacent / P-other-boundary / P-total.
- Metric denominators use VSMALL (not SMALL=1e-15): P-block squared norms are
  O(1e-18) << SMALL, so `max(nB,SMALL)` corrupted in-probe relL2/cos prints
  (raw artifacts were always correct). Fixed.

## 3. Load-bearing smoke numbers (authoritative Python recomputation from artifacts)

- P1 (dir1, eps=1e-2): **Jv-P vs FD relL2 = 1.44e-08, cos = 1.0, normRatio = 1.0**;
  outlet-adjacent 5.46e-09, off-outlet 5.30e-08; maxDiff 3.99e-18.
- **Jv-P vs BFINAL-007 verified FD_lin (stageB7_RP_FD_lin.mtx): relL2 = 2.95e-09**
  — the corrected forward J closes the independently-verified actual-residual FD.
- Probe FD vs BFINAL-007 FD_lin: 1.42e-08; vs BFINAL-007 central FD: 6.14e-07
  (the known central-diff roundoff of the non-orthogonal correction).
- P1 eps-plateau (full-P relL2): 1e-2: 1.44e-08, 1e-3: 1.17e-07, 1e-4: 1.04e-06,
  1e-5: 8.23e-06 (exact linearity => best closure at largest eps).
- P3 (eps=1e-2): dir2 P-total relL2=1.52e-04, dir3 P-total relL2=1.53e-04,
  U-rows relL2 ~1.65e-04 (consistent with the known O(1e-4) derivative floor and
  BFINAL-003 momentum anchor 1.653e-4).

## 4. Gates

P1-P8 evaluation is S3+ (not part of S2). S2 delivers the probe + build + raw
evidence artifacts.

Raw artifacts: `evidence/agent-group/BFINAL-008/cycle-1/artifacts/` (23 files:
`stageB8_dir{1,2,3}.mtx`, `stageB8_Jv{1,2,3}.mtx`, `stageB8_FD{1,2,3}_eps_*.mtx`,
`stageB8_nP.mtx`, `stageB8_celltype.mtx`, `stageB8_outlet_cells.mtx`,
`stageB8_metrics.txt`, `sha256_artifacts.txt`); run log
`/home/ys/dsH/b8_scratch_on/Log.stageB8.on.txt` (2 identical probe executions).
