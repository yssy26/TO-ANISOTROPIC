# stageB6RxDesignOracle.H — design note (BFINAL-004, step-1)

- Stage: B-final · Mode: DIAGNOSTIC_ONLY · Step: step-1-rx-probe (executor)
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace) @ `ca8a772b`
  branch `agent/dsH-stage-b-validation`

## What was added (files_allowed only)

1. `src/stageB6RxDesignOracle.H` (NEW, switch-guarded diagnostic probe,
   optProperties switch `stageB6RxDesignOracle`, default false).
2. `src/solveDiscreteFlowAdjoint.H` — ONLY a guarded include block appended
   AFTER the existing stageB5 block (end of file), mirroring the stageB5
   pattern: `if (lookupOrDefault<Switch>("stageB6RxDesignOracle", false)) {
   static done; if (!done) { done=true; if (discreteAdjointLabel ==
   "pressureDrop") { #include "stageB6RxDesignOracle.H" } } }`.

Verified: the stageB6-only delta vs the BFINAL-003 pre-state is exactly this
guarded include block (32 diff lines, 6 mentions of stageB6RxDesignOracle).
No production-mathematics line changed; forbidden files untouched
(`git diff --stat` = the two pre-existing BFINAL-003 modified files only;
`src/stageB6RxDesignOracle.H` is new/untracked).

## Why the probe runs in the pressureDrop block only

`solveDiscreteFlowAdjoint.H` is included twice per outer iteration:
- AdjNS_HT.H (label `thermalCoupling`): `discreteAdjointVelocity=Ub`, `pb`.
- AdjNS_PD.H (label `pressureDrop`): `discreteAdjointVelocity=Uc`, `pc`.

The §9/§10 adjoint weighting needs the EXACT PRESSURE-DROP adjoint (Uc/pc),
so the probe runs only when `discreteAdjointLabel == "pressureDrop"` (static
guard: once per program run, first outer iteration, after the adjoint solve).

## Locked residual semantics (FD side, identical to BFINAL-003)

`rebuildResidual(alphaMod, phiOut, resU, rp, rAUrel, Hrel, HbyA)` replicates
stageB5BoundaryRelaxOracle.H `rebuildProduction` verbatim with alpha as a
parameter:
- UEqn: `fvm::div(phi,U) - fvm::laplacian(nuEffFrozen,U) + fvm::Sp(alphaMod,U)
  == -fvc::grad(p) + fvOptions(U)`; `-= fvc::div(nuEffFrozen*dev2(T(grad U)))`
  when ransFlowModel; `relax()`; `fvOptions.constrain()`.
- `R_U = -UEqn.residual()` (relaxation-invariant).
- `rAU = 1/A`; `HbyA = constrainHbyA(rAU*H, U, p)`; `phiHbyA = fvc::flux(HbyA)`;
  `adjustPhi`; `rAtU` (consistent=false => rAU in this case); `constrainPressure`;
  `pEqn = laplacian(rAtU,p) == div(phiHbyA)`; `setReference(pRefCell)` (no-op,
  p.needReference()==false); `phi = phiHbyA - pEqn.flux()`; boundary flux
  REBUILT.
- `R_P = div(phi)` per cell = raw face-flux sum (owner +, neighbour -, boundary +),
  NO 1/V.

U/p/phi/k/omega/nutFrozen/nuEffFrozen are NEVER modified; UEqn/pEqn are NEVER
solved (only residual()/A()/H()/flux()).

## Analytic side (separately implemented linearization, no shared formula)

- RX-A:
  - `dR_U/dalpha = +V*U` per cell (residual relaxation-invariant => exact
    affine in alpha; FD == analytic to machine precision).
  - `dR_P/dalpha = div(dphiHbyA - dflux)` with:
    - `drAU = -alphaRel*rAU^2` (A()=D/V so the fvm::Sp V factor cancels; NOTE
      the earlier V-factor bug was fixed),
    - `dHbyA = rAU*((alphaRel-1)*U - alphaRel*HbyA)` (relax-source + mobility),
    - `dphiHbyA_f = Sf&(w*dHbyA_o+(1-w)*dHbyA_n)`; boundary 0 on non-assignable
      U patches (constrainHbyA), extrapolated (cell value) on assignable,
    - `dflux = flux(fvm::laplacian(drAU, p))` materialized over ALL faces
      (surfaceScalarField operator[] covers internal faces only — fixed).
    - per-term decomposition (div(dphiHbyA), -div(dflux)) logged + exported.
- RX-B: `alpha(xh) = max(alphamin, alphaMax*qu*(1-xh)/(qu+xh+SMALL))` exact
  production interpolation; analytic = `R_alpha*(dAlphaDxh*deltaXh)` with the
  production `dAlphaDxh` (mask/alphamin-clip zeroed).
- RX-C: FD = exact chain x±eps*d -> Helmholtz filter (production equation +
  PCG/DIC tol 1e-9) -> adaptive-Heaviside (bisected eta) -> alpha -> residual;
  analytic = `R_alpha*(dAlphaDxh*z)` with `z = dxh/dx*d` via the exact filter
  tangent solve `(M-B)y = -b d` and the exact projection Jacobian
  `z = drho*y + dProjEta*[sum(designMask*(1-drho)*V*y)/sum(designMask*dProjEta*V)]`
  (true tangent incl. implicit adaptive-eta derivative — the adjoint of the
  production filter_chainrule formula).
- Internal gates: xh-tangent FD vs analytic z; alpha-tangent FD vs
  dAlphaDxh*z.

## §9/§10

- `||R_P,x d||/||R_U,x d||` (analytic, eps-independent).
- `D_momentum = -Uc^T R_U,x d`, `D_pressure = -pc^T R_P,x d`,
  `D_total = D_momentum + D_pressure`.
- Production replica (momentum-only): `gsenshPressureDrop =
  -dAlphaDxh*(U&Uc)*V` -> designMask -> drho + eta correction -> adjoint filter
  (PCG/DIC 1e-9) -> designMask -> `gsensR`; `prodContrib = gsensR*d`;
  `prodDgdx1*d = pressureGradientScale*prodContrib`.
- Volume regression anchor: `gsenshVol = designMask*V/designVolume` -> same
  projection chain -> `projV_D1` must equal 0.1555186511144686.

## Exports (written to the case working directory)

stageB6_rxa_FD_all_eps.mtx (5x4N: U3N+P N per eps), stageB6_rxa_analytic.mtx
(4N), stageB6_rpa_terms.mtx (3N: div(dphiHbyA), -div(dflux), total),
stageB6_rxb_FD_all_eps.mtx / stageB6_rxb_analytic.mtx (5x4N / 4N),
stageB6_rxc_FD_all_eps.mtx (15x4N: D1..D3 x 5 eps), stageB6_rxc_analytic.mtx
(3x4N), stageB6_rxc_z_analytic.mtx (3N), stageB6_rxc_xh_FD_all_eps.mtx (15N),
stageB6_rxc_xh_analytic.mtx (3N), stageB6_rxc_alpha_FD_all_eps.mtx (15N),
stageB6_lambda.mtx (Uc 3N + pc N), stageB6_deltaAlpha.mtx / stageB6_deltaXh.mtx
(N), stageB6_dalphadxh.mtx (N), stageB6_dirs.mtx (3N), stageB6_celltype.mtx
(N), stageB6_proj.mtx (drho N), stageB6_prod_gsensh.mtx / stageB6_prod_gsens.mtx
(N), stageB6_prod_gsensVol.mtx (N).

## Acceptance for step-1 (this step)

- Header compiles: wmake exit 0 (build_stageb6.log). Precondition only.
- The only solveDiscreteFlowAdjoint.H delta = guarded include (+ new header).
- Switch-off path: single dictionary lookup, no-op (verified by construction).
- FD side = full production rebuild (verbatim locked semantics); analytic side
  = independent linearization (no shared formula).

## Notes for the run step (step-3, next executor)

- Run on a writable copy of /home/ys/dsH/b2_case_smoke; enable in
  constant/optProperties: `stageB4JacobianProbe true;`
  `stageB5BoundaryRelaxOracle true;` `stageB6RxDesignOracle true;`
  (everything else unchanged).
- Expect BFINAL-003 anchors unchanged (StageB5 momentum 1.653e-4,
  P-total 1.528e-4, rAU avgRatio 2.5, volume projection 0.1555186511144686).
- State fingerprints at probe entry/exit must be identical.

## CRITICAL build-environment finding (cost two ~20-min builds)

`source /opt/openfoam7/etc/bashrc` UNCONDITIONALLY overrides
`FOAM_USER_APPBIN` (etc/config.sh/settings L163:
`export FOAM_USER_APPBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin`).
Exporting FOAM_USER_APPBIN BEFORE sourcing is silently lost, and the link then
fails with `cannot open output file
/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin/MTO_HF: Permission
denied` (exit 2).  The correct order is:

```
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin   # AFTER source
unset FOAM_SIGFPE
cd /home/ys/dsH/TO-ANISOTROPIC/src && wmake
```

This produced WMAKE_EXIT=0 and a fresh
`/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF` (22:15, 4691928 bytes).  If only
the link step needs redoing, MTO_HF.o is already fresh (22:10) so wmake
re-links in seconds.
