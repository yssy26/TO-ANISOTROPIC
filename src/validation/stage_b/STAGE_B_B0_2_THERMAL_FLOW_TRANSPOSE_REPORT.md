# Stage B0.2 — Exact Discrete Thermal-to-Flow Coupling Transpose Audit and Correction

- Repository: `/home/ys/TO-ANISOTROPIC`
- Branch: `agent/stage-a-anisotropic-validation`
- Baseline commit: `f3c4aa25e4ae46b4bed80143cc1ef81ab6c7c68f` (with the B0.1 uncommitted fix)
- OpenFOAM-7, gate4 case (`/home/ys/Twostream/validation/gate4`), frozen state at time `1`
- Audit date: 2026-08-06
- Scope: Stage B0.2 only. No B1, no MMA, no `gradientValidated=true`, no full SST Gate 6 acceptance, no change to objective/constraints/turbulence/thresholds. commit/push/PR: none.

## 1. Actual convection discretization of `div(phi,T)`

`system/fvSchemes`:
```text
div(phi,T)      bounded Gauss upwind;
```
OpenFOAM-7 implements `bounded Gauss upwind` as
```cpp
boundedConvectionScheme::fvmDiv(faceFlux, vf)
    = scheme_().fvmDiv(faceFlux, vf)          // Gauss upwind matrix
      - fvm::Sp(fvc::surfaceIntegrate(faceFlux), vf);   // diag -= V*div(phi)
```
i.e. the upwind operator MINUS a `Sp(div(phi))` diagonal correction. This
correction is the reason the boundary behavior differs from plain upwind and
must be included in any exact transpose derivation (Section 4).

## 2. Actual discrete thermal residual

The solver operator (HeatTransfer.H) is
```cpp
fvScalarMatrix TEqn(fvm::div(activeThermalFluxF, T, "div(phi,T)")
                  - fvm::laplacian(DTEffective, T) == Q);
```
`fvMatrix::residual()` computes `r = source_ + bC - (diag + iC)*T - upper/lower*T`
(physical patches; `addBoundarySource` adds `bC`, `addBoundaryDiag` adds `iC`).
This is the exact residual used for the discrete adjoint transpose. All FD
tests below re-assemble `M.residual()`; a long-double hand assembly of the
same expression matches `residual()` to `2.86e-17` (machine precision),
eliminating summation round-off so that per-face FD reaches `1e-10` accuracy.

`phiThermal = coldFaceMask*phi + hotRhoCpRatio*hotFaceMask*phiHotFrozen`
(updateThermalFlux.H); on cold faces `d phiThermal/d phi = 1`, on hot faces 0.

## 3. Why the B0.1 "inconsistency" conclusion was wrong

B0.1 measured `<Tb, dr>` with an inner product weighted by the CELL VOLUME:
`I = sum_cells V_cell*Tb_cell*r_cell`. The fvMatrix residual `r` is already the
volume-INTEGRATED residual, so the volume weight is spurious and shrinks the
measured value by `V ~ 1.25e-10`, giving `1.8e-14` instead of the true `1.5e-4`.
The exact transpose inner product for a volume-integrated operator is the plain
vector dot product `<Tb, dr> = sum_cells Tb_cell*dr_cell` (no V). With the
correct inner product:

| internal face 16330 (phi>0) | value |
|---|---:|
| `<Tb, dr>` (correct, no V) | **+1.4668e-4** |
| current code `-Tb_downwind*(T_nei-T_own)` | **+1.4668e-4** |
| B0.1 reported (spurious V-weighted) | 1.8e-14 |

**The current internal-face expression is the exact fvm transpose**; the B0.1
"inconsistency" finding was an artifact of the volume-weighted test.

## 4. Exact internal-face formula and matrix extraction path

From the actual matrix (bounded Gauss upwind), for internal face `f`
(owner `o`, neighbour `n`, upwind sign held):

```text
lower[f] = -w*phi_f            (A[o][n])
upper[f] = (1-w)*phi_f         (A[n][o])
negSumDiag:  diag[o] -= upper[f]; diag[n] -= lower[f]
bounded:     diag[o] -= V_o*div(phi)_o ;  diag[n] -= V_n*div(phi)_n
```
Measured coefficient sensitivities (per-face FD of the matrix, face 16330):
`d(diag)[n]/dphi = +1`, `d(diag)[o]/dphi = 0`, `d(lower)/dphi = -1`,
`d(upper)/dphi = 0` (upwind w=1). With `R = s - A*T`:

```text
dr[o] = -d(diag[o]*T_o)/dphi - d(upper*T_n)/dphi = 0
dr[n] = -d(diag[n]*T_n)/dphi - d(lower*T_o)/dphi = -(T_n - T_o)
lambda_phi_f = sum_cells dr_cell*Tb_cell = -Tb_downwind*(T_n - T_o)
```
This is EXACTLY the code's `-coldThermalFaceMask*adjointDownwind*
temperatureJump` (AdjNS_HT.H). For `phi<0` the same result holds with
`downwind = owner`.

**Extraction path (priority order used)**: (1) fvm::div coefficient derivation,
(2) the actual `lower/upper/internalCoeffs/boundaryCoeffs` of `fvMatrix`,
(3) per-face FD of `residual()` — all three agree.

## 5. Boundary-face results (all explained, not just "≈0")

### 5.1 cold outlet (zeroGradient T, outflow, `coldFaceMask=1`)
The upwind boundary `internalCoeffs = phi` enters `diag` via `addBoundaryDiag`;
the bounded correction `-V*div(phi)` also enters `diag`; their sensitivities
w.r.t. `phi_out` are `+1` and `-1` (measured: `d(diag)/dh = -1`,
`d(iC)/dh = +1`, `d(diagFull)/dh = 0`). Therefore the residual is independent
of `phi_out`: **`dR_T/dphi_out == 0` exactly** (diag+iC cancellation, not a
physical independence). Measured `lambda_fd ~ 7e-15` (`relToTbT ~ 3e-17`).
The code's zero value is correct. B0.1's numeric conclusion stands; the
reason is the bounded-cancel, now proven from the actual coefficients.

### 5.2 cold inlet (fixedValue T=600, inflow, `coldFaceMask=1`)
`bC = -phi*T_bnd`, `iC = 0`; bounded `diag -= V*div(phi)` with
`V*div(phi)_cell` containing `+phi_in`. Measured `dr[cell] = -(T_bnd - T_cell)`
(= `T_cell - T_bnd`, here ~`4.4e-4`), so
`lambda_phi(inlet) = (T_cell - T_bnd)*Tb_cell ~ -1.5e-4` (max rel err over
faces/steps `2.95e-8`). This face flux is NOT a design-free flow variable:
inlet U is `fixedValue`, so `solveDiscreteFlowAdjoint.H` skips fixedValue-U
boundaries and the inlet flux never enters the flow unknowns. No inlet
faceAdjoint fill is required or performed.

### 5.3 hot inlet / hot outlet (`coldFaceMask = 0`)
`phiThermal` does not depend on the cold flux there
(`hotPart` design-fixed); measured `lambda_fd = 0`, mask = 0.

### 5.4 walls / solidEndWalls
`coldFaceMask = 0` (wall cells not in the design cold flow) or zero flux;
no contribution.

## 6. Internal-face per-face FD table (summary; full table in TSV)

20 significant cold internal faces (9 phi>0, 9 phi<0, 2 near-zero skipped for
sign-flip protection), steps `{1e-3,3e-4,1e-4,3e-5} * |phi_f|`:

| metric | value | gate |
|---|---:|---:|
| median relative error | **2.05e-10** | <= 1e-8 |
| maximum relative error | **1.56e-9** | <= 1e-6 |
| sign mismatches | **0** | 0 |

The reference is the current-code expression `-Tb_downwind*(T_nei-T_own)`.
Near-zero faces (`|phi| < 1e-6`, where a perturbation flips the upwind sign and
the operator is non-differentiable) are reported as SKIP with the reason.

## 7. Boundary-face per-face FD table (summary; full in TSV)

| patch | lambda_fd | explanation | gate |
|---|---:|---|---|
| outlet face 0-5 | 6.9e-15 .. 0 | bounded iC/diag cancel, exact zero | relToTbT <= 3e-17 |
| inlet face 0-5 | matches `(T_cell-T_bnd)*Tb_cell` | fixedValue bC source | max rel 2.95e-8 |
| hotInlet/hotOutlet | 0 | coldFaceMask = 0 | exact |

## 8. Matrix-level transpose dot-product tests

Per-face FD (Section 6) is the strongest form: it verifies every column of
`(dR_T/dphi)^T Tb` individually at `1e-10`. Additionally:

1. Combined face directions (`lambda_phi^T dphi` vs FD of `<Tb,R>`), amplitude
   `1e-3`, steps `{1e-3..3e-5}*1e-3`: all significant directions agree in
   magnitude and sign; best-step relative errors: internal_sin `8.2e-7`,
   internal_outlet `5.8e-6`, internal_inlet `1.7e-5` (FD round-off dominated);
   outlet-only directions give `D_an = 0 = D_fd ~ 1e-13` (verifies the outlet
   residual term is zero); median over all combined samples `1.5e-6`.
2. End-to-end `delta U -> delta phi = fvc::flux(deltaU)` (5 smooth
   physical-coordinate directions; SIMPLE pressure-flux coupling not
   transposed here — that is the `solveDiscreteFlowAdjoint.H` dot-test, which
   passes at `7.44e-14` in the smoke run): dU_x rel `3.5e-4`, dU_xy rel
   `1.2e-6`, dU_y a zero direction (`D_an ~ -1e-13`, `D_fd ~ 1e-10` noise,
   verifies y-flux independence). No sign mismatches in any direction.

## 9. Combination of dJ/dphi and the residual flux term

Current Lagrangian sign convention (verified by construction): the discrete
flow adjoint solves `A^T lambda = RHS` with
`RHS = (dJ/dq) + (dR_T/dq)^T Tb`, and `dR_T/dq = dR_T/dphi * dphi/dq`.
Therefore the face adjoint is

```text
faceAdjoint = dJ/dphi + (dR_T/dphi)^T Tb
internal faces:  dJ/dphi = 0            -> faceAdjoint = -Tb_downwind*(T_n-T_o)
cold outlet:     (dR_T/dphi)^T Tb = 0   -> faceAdjoint = dJ/dphi = -(T_f-Tmix)/(M*Tref)
```
Both terms carry the positive sign into `discreteExternalFaceFluxAdjoint`
(`+=` on the outlet; direct assignment on internal faces). The sign convention
is confirmed end-to-end by the full-objective smoke test (Section 11): a sign
error would flip `D_ADJ` relative to `D_FD`, which does not happen
(`F3 obj_err = 2.4%`, same sign).

## 10. Modified files

`src/AdjNS_HT.H` (+60 lines, diagnostics and comments ONLY; no physics change):
- updated the outlet-boundary comment to document the proven bounded-cancel
  (dR_T/dphi_out = 0);
- added a diagnostic block printing
  `maxInternalThermalResidualFluxAdjoint`,
  `maxBoundaryThermalResidualFluxAdjoint` (0),
  `maxObjectiveFluxDerivative`, `maxCombinedFluxAdjoint`.

The internal-face expression itself was NOT changed: it is the exact
transpose (Section 4). B0.1's `computeObjective.H` /
`createFrozenHotRegionFields.H` changes remain as-is (correct). Clean rebuild:
`src/validation_stage_b02_build.log` (0 errors).

## 11. Full-objective single-direction smoke test

Case `/home/ys/smoke_b02` (gate4 copy, `flowModel incompressibleRANSFrozen`,
`adjointMode discrete`, `freezeTurbulenceForValidation true`,
`freezeColdFlowForValidation false`, `solveFlowAdjoints true`,
`mmaUpdateEnabled false`, `objectiveGradientScale 1.0`,
`pressureGradientScale 1.0`, 1 outer iteration; F1/F2/F3 run automatically):

| item | value |
|---|---:|
| discrete flow-adjoint dot-test (thermalCoupling) | 7.44e-14 |
| F1 [GATE4] (fixed-phi thermal only) | median 2.25e-4, max 2.49e-3, 0/9, sign 0 |
| F2 full frozen U-p-T | maxRelError 6.66e-2 (3/9 fails, informational) |
| **F3 combined full-objective direction** | **obj_err 2.44e-2 (2.4%)** |
| F3 sign | same sign (no mismatch) |

F3 re-solves the FULL frozen U-p-T primal for `x +/- h*d` (i.e. it includes
the dJ/dphi path through the perturbed outlet flux) and compares against the
adjoint `D_ADJ(J) = -0.013899`; the `2.4%` agreement at every step is far
inside the `<=10%` smoke gate and confirms the complete chain
`dJ/dphi -> discrete flow adjoint (Ub) -> dfdx` with the correct sign.
(The subsequent GATE6-STRICT baseline repeatability FAIL is expected in this
configuration: full-SST Gate 6 is out of scope and requires converged SST
primals; it is not part of the B0.2 acceptance.)

## 12. F1 regression (after B0.2, serial + 4 ranks)

| metric | serial | 4-rank parallel | requirement |
|---|---:|---:|---:|
| medianRelError | 2.25e-4 | 9.95e-5 | <= 1% |
| maxRelError | 2.49e-3 | 1.51e-3 | <= 5% |
| failures | 0/9 | 0/9 | 0/9 |
| sign mismatch | 0 | 0 | 0 |
| D_ADJ dir 0/1/2 | -4.056556e-3 / 1.731376e-4 / 5.671516e-5 | identical | - |
| dir 1,2 max serial/parallel rel diff | 3.65e-8 | - | <= 1e-6 |

D_ADJ values are bit-identical to B0.1 (the B0.2 diagnostics do not touch the
thermal-only gradient; F1 fixes phi).

## 13. B0 closure

B0 is closable. The three B0 checks that B0.1 re-classified (B0-01, B0-17,
B0-20) remain PASS; the B0.1 report's incorrect "internal-face expression
inconsistent" finding is hereby superseded by the exact-transpose proof of
Section 4 (the B0.1 test multiplied the inner product by cell volume).
B0.1's dJ/dphi fix on the cold outlet is confirmed correct (Section 9, 11).

## 14. B1 readiness

`readyForB1 = true`. Preconditions (unchanged from B0/B0.1):
- `objectiveGradientScale 1.0; pressureGradientScale 1.0;`
- serial run (reduced discrete flow adjoint is serial-only);
- pressure-gradient FD in the flow-coupled configuration
  (`solveFlowAdjoints=true`, `freezeColdFlowForValidation=false`), comparing
  `dgdx[1]` against the normalized `gDP`;
- B1 must also re-validate the internal-face coupling term with the
  no-volume inner product used here, and the pressure-drop chain separately.

`gradientValidated` remains `false`; no optimization started.

## 15. Deliverables

- `src/validation/stage_b/STAGE_B_B0_2_THERMAL_FLOW_TRANSPOSE_REPORT.md`
- `src/validation/stage_b/stage_b_b0_2_transpose_results.json`
- `src/validation/stage_b/stage_b_b0_2_face_fd.tsv`

## 16. Final status

**PASS** — all 10 acceptance criteria satisfied; the only source change is the
non-physical diagnostic block in `AdjNS_HT.H`.
