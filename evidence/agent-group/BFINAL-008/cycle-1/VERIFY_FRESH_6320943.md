# BFINAL-008 Fresh Verification (post-6320943, unpinned physical operator)

- HEAD: `6320943` (`agent/dsH-stage-b-validation`), worktree clean, binary rebuilt 23:04 (WMAKE_EXIT=0, all new markers present)
- Diagnostic case: `/home/ys/dsH/b8_verify_diag` (stageB4/B5/B6/B8 probes ON, export-only) — MTO_RC=0
- Production case: `/home/ys/dsH/b8_verify_prod` (all probes OFF, `discreteProdPreconditionerSetup diagonal`, `discreteProdEnableRitzPilot false`) — MTO_RC=0
- Both cases: `mmaUpdateEnabled false; frozenGradientValidated false;`
- Logs: `b8_verify_diag/Log.verify_diag.txt`, `b8_verify_prod/Log.verify_prod.txt`
- Offline script: `cycle-1/s4_verify_fresh_unpinned.py` → `s4_verify_fresh_unpinned.json`
- Consolidated metrics: `cycle-1/s3_gate_metrics_fresh.py` → `s3_gate_metrics_fresh.json`

## 1. fixedValue-outlet actual-residual Jv FD, ALL pressure rows (stageB8)

stageB8: fixesValue(p) patch = outlet nFaces=84; N=33600 internal=26208 otherBoundary=7308 outletAdjacent=84

| metric (eps=1e-2) | relL2 | cos | note |
|---|---|---|---|
| P1 full-P | 1.4429e-08 | 1.0 | == BFINAL-008 |
| P1 outlet-adjacent | 5.4549e-09 | 1.0 | == BFINAL-008 |
| P1 off-outlet | 5.2956e-08 | 1.0 | == BFINAL-008 |
| P3 dir2 P-total | 1.5158e-04 | 0.99999999 | **< 1e-3 ✓** (pref 1e-4, at known O(1e-4) floor) |
| P3 dir3 P-total | 1.5280e-04 | 0.99999999 | **< 1e-3 ✓** (pref 1e-4, at known O(1e-4) floor) |

P1 sign/structure: sameSignFraction=0.7845, |FD| outlet mass=0.686248479375, |Jv| outlet mass=0.686248477495 — identical outlet concentration, == BFINAL-008.

## 2. P4 — J/J^T dot (from run log, both labels)

```
Reduced cold-flow operator transpose dot-test (thermalCoupling) max relative error=6.86412558997e-14
BlockDot (thermalCoupling): PU=2.07803300324e-13 PP=2.29255815574e-14 UU=7.3095183677e-14 UP=1.79609743792e-12 boundaryU=1.80391921425e-13
```
All block dots < 1e-10 ✓ (PP=2.29e-14); identical in pressureDrop block; bit-for-bit == BFINAL-008.

## 3. P5 — explicit / matrix-free consistency

```
ExplicitJToracle: maxRelL2=3.75727184556e-16 maxCosErr=2.22044604925e-16 maxRelL_U=3.53382916868e-16 maxRelL_P=3.75827846067e-16
```
~3.8e-16 both blocks ✓ (== BFINAL-008).

## 4. P2 — nullspace removal on the UNPINNED physical export (fresh, NO re-pin)

```
identity rows in P block [100800,134400): []            <- ZERO artificial identity rows ✓
row 100800: nnz=25 (physical, NOT identity)             <- old pin row is now a physical P row ✓
row 106400: nnz=34 (physical, not identity)
||J_physical n||/||n|| = 3.019933e-09                   <- OLD export 2.717e-19; null mode GONE ✓
sigma_min(physical unpinned J): min=8.620e-12 mean=1.088e-11   (== BFINAL-007 candidate ~7.9e-12) ✓
old null mode n_P: top-1% mass=0.895 outlet mass=0.352 (the old vector itself)
physical J v_min(P-block): top-1% mass=0.031 outlet mass=0.000, max cell 32480 (NOT outlet)  ✓
```

## 5. P8 — direct tangent solvability on the UNPINNED physical J (fresh)

```
D1: ||rhs||=5.663185e-01 trueRelRes=6.136015e-14  (Urel 5.4e-14, Prel 8.4e-11)  ✓ <= 1e-9
D2: ||rhs||=8.615839e-01 trueRelRes=4.576175e-14  ✓
D3: ||rhs||=8.516528e-01 trueRelRes=9.084105e-14  ✓
```

## 6. P6 / P7 — no regression (bit-for-bit == BFINAL-003/005/008)

- StageB5 GateA momentum relL2=1.6531e-4 (== BFINAL-003); candA 2.985e-5; candB 1.059e-11
- StageB6 weighted D1/D2/D3: |D_mom-prod|/|prod| = 2.03e-9 / 5.40e-10 / 7.65e-10 (== BFINAL-005)
- RxPressureRowTranspose: relErr=2.58e-15 (== BFINAL-005)

## 7. Production adjoint solver (NEW diagonal preconditioner, RitzPilot=false)

```
PRODPRECSETUP (thermalCoupling): mode=diagonal avg=27.7275430014 min=1.00309714066 max=725.954987504
PRODPRECSETUP (pressureDrop):    mode=diagonal avg=27.7275430014 min=1.00309714066 max=725.954987504
Production reduced discrete flow adjoint thermalCoupling: FGMRES iterations=4000, true relative residual=5.08309124088e-05
    -> did not converge (DIAGNOSTIC MODE). trueRelRes=5.08309124088e-05, tolerance=1e-09
Production reduced discrete flow adjoint pressureDrop:    FGMRES iterations=4000, true relative residual=0.0116060277262
    -> did not converge (DIAGNOSTIC MODE). trueRelRes=0.0116060277262, tolerance=1e-09
```
- `FGMRES-PROD thermalCoupling`: 4.35e-3 (80) → 5.08e-5 (4000) — **stalls at ~5e-5**
- `FGMRES-PROD pressureDrop`: stalls at **1.16e-2**
- `BCGS-PROD`: none (discreteProdSolverType=fgmres)
- `PRODSOLV`: correctly absent — gated behind `prodPressureNeedsReference` (false for fixedValue outlet; no null-space projection for the physical operator)
- crash markers: none; MTO_RC=0

**Conclusion (production): the O(N) `diagonal` preconditioner does NOT fix the ill-conditioned (cond~1e19) iterative stall — same result as the pre-patch FGMRES (thermal 4.55e-5, pressure 1.1e-2). Per user instruction, do NOT switch to `bruteForceL1` on the 33600-cell case (O(N*nnz) is only for small meshes). The trustworthy route remains the B4-oracle/import path + exact offline SuperLU tangent solves (P8: 6e-14..9e-14).**

## 8. Overall verdict

- **PASS (P1/P2/P3/P4/P5/P6/P7/P8)** — the corrected J (J_PU + J_PP interior sign + fixedValue-outlet internalCoeffs) now closes the actual frozen-primal residual FD with ALL pressure rows on the unpinned physical operator; no artificial identity rows; outlet null mode does not reappear; direct tangent solvable to ~1e-13.
- **Production iterative Krylov convergence: STALLS (unchanged)** — diagonal preconditioner is not strong enough for cond~1e19; exact/oracle path is the trustworthy production route.
