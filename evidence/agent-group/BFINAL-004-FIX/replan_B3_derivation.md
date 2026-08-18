# BFINAL-004-FIX — Replanner note v2: B3 corrected analytic dHbyA/dalpha & drAU/dalpha

v2 correction of the Pre-Reviewer REJECTED_PHYSICS verdict. v1 (superseded) wrongly
added an explicit `mesh.V()` factor. The OpenFOAM-7 source proves `rAU = 1/UEqn.A()`
ALREADY carries V, so the correct alpha-derivatives have NO extra V.

## Facts locked by source inspection

- `fvm::Sp(alpha,U)` adds `mesh.V()*alpha` to the fvMatrix **integrated** diagonal
  (`/opt/openfoam7/src/finiteVolume/finiteVolume/fvm/fvmSup.C:118`):
  => unrelaxed integrated diagonal `D0 = D_rest + alpha*V`, so `dD0/dalpha = V`.
- `fvMatrix::A()` returns `D()/mesh.V()` (per-volume diagonal),
  `/opt/openfoam7/src/finiteVolume/fvMatrices/fvMatrix/fvMatrix.C:751`:
  => `A() = D_rel/V`, and `rAU = 1/UEqn.A() = V/D_rel` (rAU ALREADY carries V, units s).
- `fvMatrix::H()` returns the per-volume H: `H() = (S - offdiag*U + boundary)/V`
  (fvMatrix.C `Hphi.primitiveFieldRef() /= psi_.mesh().V()`).
- `UEqn.relax()` (fvMatrix.C:521): `D_rel = D0/alphaRel`, `S_rel = S0 + (D_rel-D0)*U`
  (assuming the diagonal-dominance `max(D,sumOff)` is a no-op, standard here).
- In the probe `rebuildResidual` (src/stageB6RxDesignOracle.H):
  `rAUrelBase = rAU = 1/UEqn.A()` (relaxed, L268/L293) and
  `HbyABase = constrainHbyA(rAU*UEqn.H(), U, p)` (relaxed, L269/L296).
  `alphaRel` = the actual fvSolution U relaxation factor (0.4 for this case).

## Correct alpha-derivatives (U, p fixed)

With `H0 = S0 - offdiag*U + boundary` (integrated, alpha-independent; alpha enters ONLY
via `fvm::Sp` which touches the diagonal, not source/offdiag/boundary):

```
HbyA = H0/D_rel + (1 - alphaRel)*U          (V cancels: rAU*H() = (V/D_rel)*(H0/V) )
drAU/dalpha = -rAU^2 / alphaRel
dHbyA/dalpha = (rAU/alphaRel) * ((1-alphaRel)*U - HbyA)
```

Check:
- `drAU/dalpha = d(V/D_rel)/dalpha = -V*(dD_rel/dalpha)/D_rel^2
   = -V*(V/alphaRel)/D_rel^2 = -(V/D_rel)^2/alphaRel = -rAU^2/alphaRel`.  (V cancels.)
- `dHbyA/dalpha = d(H0/D_rel)/dalpha = H0*(-1/D_rel^2)*(V/alphaRel)
   = -(H0/D_rel)*(V/D_rel)/alphaRel
   = -(HbyA-(1-alphaRel)U)*rAU/alphaRel
   = (rAU/alphaRel)*((1-alphaRel)U - HbyA)`.  (No extra V.)

## Why v1's explicit V was a regression (category error)

The momentum anchor `dR_U/dalpha = V*U` carries V because `R_U = -residual()` is an
INTEGRATED quantity (`fvMatrix::residual()` does NOT divide by V). The P-row
`R_P = div(phi)` is already integrated (raw face-flux sum, no 1/V), and its V-dependence
flows through `HbyA = rAU*H()` where `rAU=V/D_rel` and `H()=H0/V`, so V cancels. The V in
the momentum row does NOT imply an explicit V belongs in the P-row rAU derivative.

## What was in the probe (wrong, pre-fix)

```
drAU[celli]  = -alphaRel*rAUc*rAUc;                                  // alphaRel instead of 1/alphaRel
dHbyA[celli] = rAUc*((alphaRel-1.0)*U[celli] - alphaRel*HbyABase[celli]); // coeffs wrong, U-term sign wrong
```
(6.25x from alphaRel^2=0.16 plus a U-term sign error -> RX-B R_P,xh relL2=1.012/cos=-0.176.)

## Required edit (B3, no V)

```
drAU[celli]  = -rAUc*rAUc/alphaRel;
dHbyA[celli] = (rAUc/alphaRel)*((1.0 - alphaRel)*U[celli] - HbyABase[celli]);
```

with `rAUc = rAUrelBase[celli]`, `HbyABase` = relaxed HbyA, `alphaRel` = 0.4 here.
Also fix the header comment (L34-36) and body comment (L380-382) to the no-V formulas.
`dflux = flux(fvm::laplacian(drAUfield,p))` assembly is unchanged; only the `drAU` values
change. RX-A `R_P,alpha` (weighted by deltaAlpha) and RX-B `R_P,xh` relL2<1e-3 / cos>0.999
against the independent central-FD oracle are the arbitrating gates.
