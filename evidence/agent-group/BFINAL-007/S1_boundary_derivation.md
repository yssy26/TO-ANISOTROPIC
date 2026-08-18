# BFINAL-007 S1 — Outlet fixedValue pressure-BC boundary term: exact derivation (Q2)

- Stage: B-final · Mode: **DIAGNOSTIC_ONLY** · Step: S1-derive-boundary-term (Q2)
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace, verified) · HEAD `ca8a772b8339a36e7906ea32a7125061c47aa280` · branch `agent/dsH-stage-b-validation`
- Case: `/home/ys/dsH/b2_case_smoke`
- Date: 2026-08-17 (Executor, BFINAL-007)
- Status: derivation complete; no production file modified (evidence-only)

---

## 0. Question answered here (Q2)

> Production pEqn 中 outlet fixedValue 对 pressure residual Jacobian 的精确 boundary
> coefficient/source contribution 是什么？(fvm::laplacian(rAtU,p) 的
> internalCoeffs/boundaryCoeffs 在 outlet fixedValue patch 上的对角/源项贡献；
> 从 OpenFOAM-7 源码与 case 0/p 边界逐项推导。)

**Exact per-face answer (per the approved plan's S1 formula):**

```
laplacian fvMatrix diagonal   +=  -rAtU_b · |Sf|_b · deltaCoeffs_b      (per outlet face, owner cell)
laplacian fvMatrix source     +=  +rAtU_b · |Sf|_b · deltaCoeffs_b · p_b  = 0   (p_b = 0)
```
i.e. in the pressure residual Jacobian convention `J_PP = -laplacian` (the P-P block of the
pressure residual derivative `dR_P/dp` with `R_P = div(phi) = div(phiHbyA) - div(pEqn.flux())`),
**the outlet-cell diagonal of J_PP must gain `+rAtU_c · deltaCoeffs_b · |Sf|_b`**
with **NO source term** (boundaryCoeffs source vanishes because `p_b = 0`) and **NO
cell-center identity** (the term is a boundary-FACE flux linearization, not a Dirichlet row).

84 outlet cells are affected (mesh boundary: `outlet { nFaces 84; startFace 96944; }`; all 84
boundary faces are owned by 84 distinct cells — verified from `owner`).

---

## 1. Production pressure equation (exact semantics, NS.H)

`src/NS.H` L1016-1034 (production pEqn):

```cpp
constrainPressure(p, U, phiHbyA, rAtU());          // L1016
while (simple.correctNonOrthogonal())               // L1018
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)   // L1020-1023
    );
    pEqn.setReference(pRefCell, pRefValue);        // L1025
    SolverPerformance<scalar> perfP = pEqn.solve();// L1026
    ...
    phi = phiHbyA - pEqn.flux();                   // L1032
}
```

Key facts:
- `rAtU() = rAU` in this case (SIMPLE `consistent=false`; `simple.consistent()` is false —
  verified in `stageB5BoundaryRelaxOracle.H` L286-291 and `rxPressureRowTranspose.H` L17-23).
  `primalPressureMobility = rAtU().primitiveField()` is recorded at `NS.H` L129.
- `phiHbyA` is independent of `p` (frozen-U SIMPLE reconstruction); all p-dependence of the
  residual is through `pEqn.flux()`.
- The continuity residual is `R_P = div(phi) = div(phiHbyA) - div(pEqn.flux())` (raw face-flux
  sum, owner `+` / neighbour `-`, NO `/V` — the convention used by the FD oracle
  `stageB5BoundaryRelaxOracle.H` L476-501 and by `rxPressureRowTranspose.H` L6).
- `p.needReference() == false` because `p` has a `fixedValue` outlet (`GeometricField::needReference()`
  returns false as soon as any patch `fixesValue()`; BFINAL-002 audit_boundary.md L209-210), so
  `setRefCell` (`readTransportProperties.H` L66-68 → OpenFOAM-7 `findRefCell.C` L40)
  is a no-op: `pRefCell` stays `0`, `pRefValue` stays `0`. The `fvSolution` `pRefCell 5600`
  is INEFFECTIVE. (BFINAL-006 S2 gate C1, re-verified this round.)

Case `0/p` boundaryField (L33626-33661): `inlet` zeroGradient, **`outlet` fixedValue uniform 0**,
hotInlet/hotOutlet/solidEndWalls/bottomWall/topWall/sideWalls zeroGradient. So the ONLY patch
with `fixesValue()` on `p` is the outlet.

---

## 2. fvm::laplacian(rAtU,p) matrix assembly — boundary coefficients

### 2.1 gaussLaplacianScheme.C — `fvmLaplacianUncorrected` (L44-88)

`/opt/openfoam7/src/finiteVolume/finiteVolume/laplacianSchemes/gaussLaplacianScheme/gaussLaplacianScheme.C`:

```cpp
fvm.upper() = deltaCoeffs.primitiveField()*gammaMagSf.primitiveField();  // L63  (internal)
fvm.negSumDiag();                                                         // L64  diag = -sum(upper) per row
...
forAll(vf.boundaryField(), patchi)                                         // L66
{
    ...
    else   // non-coupled patch (outlet is non-coupled)
    {
        fvm.internalCoeffs()[patchi] = pGamma*pvf.gradientInternalCoeffs();        // L82
        fvm.boundaryCoeffs()[patchi] = -pGamma*pvf.gradientBoundaryCoeffs();       // L83
    }
}
```

where `pGamma = gammaMagSf.boundaryField()[patchi]` and for scalar `gamma = rAtU`:
`SfGamma = mesh.Sf() & gamma = gamma*Sf`, `SfGammaSn = SfGamma & Sn = gamma*|Sf|`, i.e.
**`pGamma = rAtU_b · |Sf|_b > 0`**.

### 2.2 fixedValueFvPatchField.C — gradient coefficients (L128-139)

`/opt/openfoam7/src/finiteVolume/fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchField.C`:

```cpp
gradientInternalCoeffs()  { return -pTraits<Type>::one*this->patch().deltaCoeffs(); }  // L130  = -deltaCoeffs_b
gradientBoundaryCoeffs()  { return this->patch().deltaCoeffs()*(*this); }              // L138  = deltaCoeffs_b*p_b
```

### 2.3 Outlet fixedValue p_b = 0 — exact boundary coefficients

Combining 2.1 + 2.2 for the outlet patch:

```text
internalCoeffs[outlet] = pGamma · gradientInternalCoeffs
                       = (rAtU_b·|Sf|_b) · (−deltaCoeffs_b)
                       = − rAtU_b · |Sf|_b · deltaCoeffs_b          < 0
boundaryCoeffs[outlet] = − pGamma · gradientBoundaryCoeffs
                       = − (rAtU_b·|Sf|_b) · (deltaCoeffs_b · p_b)
                       = − rAtU_b · |Sf|_b · deltaCoeffs_b · p_b
                       = 0                                          (p_b = 0)
```

### 2.4 How they enter the solved linear system — fvMatrixSolve.C `solveSegregated` (L127-148)

`/opt/openfoam7/src/finiteVolume/fvMatrices/fvMatrix/fvMatrixSolve.C`:

```cpp
Field<Type> source(source_);
addBoundarySource(source);                      // L134  source[cell] += boundaryCoeffs  (= 0 for outlet)
...
addBoundaryDiag(diag(), cmpt);                  // L148  diag[cell] += internalCoeffs   (= -rAtU_b|Sf|_b δ_b)
```

`fvMatrix.C` L110-125 (`addBoundaryDiag`): `diag[cell] += internalCoeffs_[patchi][face]`.
`fvMatrix.C` L144-158 (`addBoundarySource`, non-coupled): `source[cell] += boundaryCoeffs_[patchi][face]`.

**Therefore the production pressure matrix has, on each of the 84 outlet cells:**
`diag += −rAtU_b · |Sf|_b · deltaCoeffs_b` and `source += 0`.

### 2.5 pEqn.flux() boundary — the flux linearization (fvMatrix.C L903-936)

`fvMatrix::flux()` (L894-901) internal faces: `faceH(p) = upper·p_nei − lower·p_own = kf·(p_nei−p_own)`.
Boundary faces (L903-936):
```cpp
InternalContrib = internalCoeffs * patchInternalField;   // internalCoeffs_b * p_c
NeighbourContrib = boundaryCoeffs;                        // (non-coupled: not multiplied by p_b)
ffbf[patchi] = InternalContrib − NeighbourContrib;        // L935
```
So `pEqn.flux()_b = internalCoeffs_b · p_c − boundaryCoeffs_b = − rAtU_b·|Sf|_b·δ_b · p_c`.
Then `d(phi_b)/dp_c = −d(pEqn.flux()_b)/dp_c = +rAtU_b·|Sf|_b·δ_b` (POSITIVE), and the P-row
residual `R_P(c) = V·div(phi)(c)` gains `+phi_b` (owner `+`), i.e. **the boundary contribution
to `dR_P/dp` at the outlet-cell diagonal is `+rAtU_b·|Sf|_b·δ_b`** — exactly the plan's
`+rAtU_c · deltaCoeffs_b · |Sf|_b` and the dead-code `applyPressureFluxCorrection`
(`solveDiscreteFlowAdjoint.H` L307-312: `corr_f = +mobility*delta*Area*p_cell`).

### 2.6 Exact source line map (Q2 citations)

| fact | source |
|---|---|
| production pEqn (laplacian == div(phiHbyA); phi = phiHbyA − pEqn.flux()) | `src/NS.H` L1016-1034 (L1020-1023, L1032) |
| upper = deltaCoeffs·gammaMagSf; negSumDiag | gaussLaplacianScheme.C L63-64 |
| internalCoeffs/boundaryCoeffs non-coupled | gaussLaplacianScheme.C L82-83 |
| pGamma = rAtU_b·|Sf|_b | gaussLaplacianScheme.C L66-69 (SfGammaSn for scalar gamma) |
| fixedValue gradientInternalCoeffs = −δ_b | fixedValueFvPatchField.C L128-131 |
| fixedValue gradientBoundaryCoeffs = δ_b·p_b | fixedValueFvPatchField.C L136-139 |
| solveSegregated addBoundarySource/addBoundaryDiag | fvMatrixSolve.C L127-148 (L134, L148) |
| addBoundaryDiag / addBoundarySource semantics | fvMatrix.C L110-125, L144-158 |
| flux() boundary = internalCoeffs·p_c − boundaryCoeffs | fvMatrix.C L903-936 (L935) |
| outlet fixedValue uniform 0; 7 others zeroGradient | `case 0/p` L33626-33661 |
| outlet nFaces=84 startFace=96944 | `case constant/polyMesh/boundary` L26-31 |
| p.needReference()==false → setRefCell no-op; pRefCell stays 0 | GeometricField::needReference (BFINAL-002 audit L209-210); findRefCell.C L40-113; readTransportProperties.H L66-68 |

---

## 3. What the CURRENT J's P-P block contains, and what it is missing (Q3 preview)

Current J P-P (all three representations — forward `applyDiscreteFlowJ`, transpose
`applyDiscreteFlowJT`, exported `explicitJT.mtx`):
- **Interior faces only**: `kf = mobF·dcfF·mafF` with `mobF = wf·primalPressureMobility[own] +
  (1−wf)·primalPressureMobility[nei]` (`solveDiscreteFlowAdjoint.H` L642-644, L666-667).
  `J_PP(own,own) += −kf`, `J_PP(own,nei) += +kf`, `J_PP(nei,own) += +kf`, `J_PP(nei,nei) += −kf`
  — this equals the laplacian fvMatrix interior (diag −Σkf, upper/lower +kf). Verified against
  the exported matrix: row P(10) diag = −7.775e-10 < 0, off-diag +1.9e-10 > 0; row P(5600)
  diag = −1.899e-09; T2 `relT2 = 7.5e-21` (`solveDiscreteFlowAdjoint.H` L2048-2116).
- **Boundary loop** (L701-753 forward; L2550-2604 export): adds only P-U (`uAssignable`
  boundary dphi) and U-P (`!pressureFixed` direct `Sf_b·p_c`) terms. **It never adds a P-P
  diagonal term on the pressureFixed outlet** — confirmed by direct source read and by the
  exported outlet-cell P-P diagonals (all equal to the interior-only sum, no boundary term).

**Missing term (Q3 answer, preview):** the production `addBoundaryDiag` contribution
`internalCoeffs_b = −rAtU_b·|Sf|_b·deltaCoeffs_b` on the 84 outlet-cell diagonals. In the
`J_PP = -laplacian` (pressure residual Jacobian) convention used by the plan/Q4b, the outlet-cell
diagonal of the candidate J must gain **`+rAtU_c · deltaCoeffs_b · |Sf|_b`** (which is
`−internalCoeffs_b`), with no source term. The full Q3 numerical diff (explicitJT.mtx P-P block
vs the S2-exported production laplacian diag/upper/lower) is executed in S3.

**Sign-convention note (for the Post-Reviewer):** the T2-validated code P-P block equals the
laplacian fvMatrix `L` (diag −Σkf, off-diag +kf), i.e. `J_PP = +L` in matrix terms, while the
plan's candidate convention `J_PP = -laplacian` is the pressure residual derivative
`dR_P/dp = -L` in the `R_P = div(phi)` sense (since `R_P = div(phiHbyA) − div(pEqn.flux())` and
`div(pEqn.flux()) = L·p/V`). Both statements describe the SAME boundary-face flux; the relative
sign between the exported interior block and the production residual derivative is exactly what
Q4b (candidate `J_PP·n_P` vs the independently-built actual-residual FD, S2) arbitrates
empirically. The boundary term's magnitude and sign derived here
(`+rAtU_c·δ_b·|Sf|_b` in the J_PP-diagonal convention) is convention-independent in magnitude
and matches both the production `d(phi_b)/dp_c` and the dead-code correction.

---

## 4. NOT a cell-center identity — and no source term

- The outlet `fixedValue p = 0` enters the pressure equation as a **boundary-face flux
  linearization** (internalCoeffs diagonal of the adjacent cell equation), NOT as a Dirichlet
  row: OpenFOAM never replaces the outlet cell row with `p = 0`; it adds
  `internalCoeffs_b = −rAtU_b·|Sf|_b·δ_b` to the cell's diagonal (via `addBoundaryDiag` during
  `solveSegregated`) and `boundaryCoeffs_b·p_b = 0` to the source. The identity-row scan of
  `explicitJT.mtx` confirms the ONLY identity row is `100800 = discretePIndex(0)` (the C++ pRef
  pin); row `106400 = discretePIndex(5600)` is a physical continuity row (18 stored / 18
  nonzero) — BFINAL-006 S2 gate C2-C4, re-verified this round.
- `boundaryCoeffs` source vanishes because `p_b = 0`: `boundaryCoeffs_b = −rAtU_b·|Sf|_b·δ_b·p_b = 0`.
- Therefore the plan's requirement holds: add the boundary-face-flux diagonal term, do NOT
  replace outlet-adjacent pressure cells with cell-center identity rows.

---

## 5. Magnitude (this case, mesh + runtime data)

- outlet: 84 faces = 84 unique cells; `deltaCoeffs_b`: min 5.333e+03, max 6.000e+03, avg 5.952e+03
  (computed from `constant/polyMesh` via `patch.delta() = nHat·(nHat·(Cf−Cn))`, fvPatch.C L142-147);
  `magSf_b = 1.250e-07` (all 84); sum |Sf|_b = 1.0500e-05.
- `rAtU` (primalPressureMobility): avg 2.823e-07, min 4.0e-09, max 2.05e-06 (run.log `NS rAtU`).
- typical boundary term `rAtU_c·δ_b·|Sf|_b ≈ 2.82e-07 · 5.95e+03 · 1.25e-07 ≈ 2.1e-10` per outlet cell.
- current exported J P-P diagonal on the 84 outlet cells: −4.04e-09 .. −1.36e-09 (interior-only
  sum; all 84 negative). The missing boundary term (≈ 2e-10) is the same order as the interior
  diagonal — consistent with it being the term that removes the outlet-concentrated pressure
  null mode (BFINAL-006: sigma_min ≈ 7.4e-27, ||J·n||/||n|| ≈ 2.7e-19, 89.5% of P-mode mass at
  the outlet plane).

---

## 6. Q5 preview — does the SAME boundary term affect R_P,x?

The boundary internalCoeffs term derived here is the **state**-Jacobian contribution
(`dR_P/dp`, the P-P block of J). The design/row derivative `R_P,x` (BFINAL-005) is a different
operator: `R_P,x = J_P·(dAlphaDxh·z)` with `J_P*w = div(dphiHbyA_w − dflux_w)` and
`dflux_w = flux(fvm::laplacian(drAU·w, p))` (`rxPressureRowTranspose.H` L8-13). There, the
design coefficient field `drAU·w` carries a default `calculated` BC with boundary value **0**
(`rxPressureRowTranspose.H` L36-39), so the interpolated boundary `gamma_b = 0` and
`dflux_b == 0` exactly — the boundary-face flux derivative of `R_P,x` vanishes in this case
(runtime dot-test arbitration at L200+). Whether the state-boundary term derived here has a
non-zero counterpart in `R_P,x` through the drAU field is quantified in S3 (Q5), but the
`dflux_b = 0` assumption of BFINAL-005 is about the design field `drAU` (BC value 0), NOT about
the state Jacobian boundary term derived in this document.

---

## 7. Files / provenance

- This document: `evidence/agent-group/BFINAL-007/S1_boundary_derivation.md`
- No production file modified; `git status` unchanged from the locked BFINAL-003/005 baseline
  (M `src/sensitivity.H`, M `src/solveDiscreteFlowAdjoint.H`,
  M `src/solveDiscreteFlowAdjointProduction.H` + pre-existing untracked files).
- Helper numeric checks (read-only python on mesh + exported `explicitJT.mtx`):
  outlet face/cell count, deltaCoeffs/magSf, rAtU stats, exported P-P diagonal signs,
  identity-row scan.
