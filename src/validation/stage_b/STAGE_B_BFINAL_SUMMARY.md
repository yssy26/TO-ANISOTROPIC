# Stage B — BFINAL Summary (pressure-drop gradient closure)

- **Branch:** `agent/dsH-stage-b-validation`
- **HEAD:** `5116f88` (after BFINAL-003/005/007/008 production patches)
- **Stage:** B-final (in progress — pressure-drop gradient gate not yet closed)
- **Date:** 2026-08-18

> The older reports in this directory (`STAGE_B_B0_*.md`, `STAGE_B_B1_*`,
> `STAGE_B_B2B3_*`) are **historical evidence** produced on an older branch /
> commit (`agent/stage-a-anisotropic-validation`, August 6-7). They predate the
> J_PU / R_x / J_PP corrections below and must not be treated as current
> numerical ground truth. The full per-round evidence lives under
> `evidence/agent-group/BFINAL-0NN/` (final reports + raw logs).

---

## Roadmap position

```
B-final -> C-MMA -> D-SST-refresh -> E-final-validation
```

Current status: **B-final**, pressure-drop gradient closure in progress. MMA and
SIMPLE-transpose remain **blocked** until the B-final pressure gate passes.

## BFINAL rounds

| Round | Mode | What was established / fixed |
|---|---|---|
| BFINAL-001 | DIAGNOSTIC | Existing B4 probe suite uses mutually inconsistent residual conventions; (a) J≠R_w, (b) R_x incomplete, (c) RHS/sign/assembly all still open. |
| BFINAL-002 | DIAGNOSTIC | P-row B-F2 error is **not** a boundary-oracle artifact; unrelaxed `rAUAdj` in J's `dphi/dU` tangent is wrong at O(1) (Candidate B relaxed semantics closes to relL2≈1e-11). |
| BFINAL-003 | **PATCH** | **J_PU closed**: P<-U tangent now uses production relaxed-SIMPLE semantics (`flux(alphaRel·rAU_u·dH_u + (1−alphaRel)·dU)` + constrainHbyA pinning). P-total relL2 1.036 → 1.528e-4; momentum 1.653e-4. |
| BFINAL-004 | DIAGNOSTIC | Fixed-state R_P,x is **NONZERO** at all layers; production pressure sensitivity was momentum-only (λ_P-weighted pressure-row contribution 15–30% of D_momentum). |
| BFINAL-005 | **PATCH** | **R_x closed**: added `−pcᵀ R_P,x` (pressure/continuity-row design derivative) to `gsenshPressureDrop` via `rxPressureRowTranspose.H`. G1–G5 pass; D_total = D_momentum + D_pressure to ~1e-9. |
| BFINAL-006 | DIAGNOSTIC | Tangent oracle exposed a **structural defect**: coupled J has an outlet-localized non-constant pressure null mode (sigma_min≈7.4e-27, ||J·n||/||n||≈2.7e-19). |
| BFINAL-007 | DIAGNOSTIC | Root cause located: J_PP omits the fixedValue outlet boundary internalCoeffs contribution and the interior P-P block has the opposite residual sign. Candidate J closes actual-residual FD to ~2.4e-9 and removes the null mode (sigma_min→~7.9e-12). |
| BFINAL-008 | **PATCH** | **J_PP closed**: interior sign `+L_int → −L_prod` + fixedValue outlet boundary internalCoeffs across forward J / J^T / explicit CSR / production. P1–P8 pass: null mode gone, direct tangent solve converges to ~5e-14, no J_PU/momentum/R_x regression. |

## Current derivative status

| Derivative layer | Status |
|---|---|
| J_PU (momentum P<-U) | **CLOSED** (BFINAL-003) |
| R_x (design derivative, momentum + pressure row) | **COMPLETE** (BFINAL-005) |
| J_PP (pressure state derivative, interior sign + outlet boundary) | **CLOSED** (BFINAL-008) |
| J/J^T transpose + explicit/matrix-free | machine-precision consistent |
| Structural singularity (outlet null mode) | **removed** (tangent solvable to ~1e-14) |
| B-final pressure-drop gradient gate | **NOT YET** — requires corrected tangent oracle (full J·w′=−R_x·d vs re-converged frozen-primal FD) |

## Required declarations (BFINAL-008 §19)

- The flow Jacobian now represents the actual frozen-primal pressure residual,
  including the fixedValue outlet boundary contribution and the correct P-P
  residual sign convention.
- J_PU remains closed from BFINAL-003.
- R_x remains complete from BFINAL-005 for the current design region.

## Next step

A corrected rerun of the pressure-drop **tangent oracle** (BFINAL-006 route with
the corrected J): solve `J·w′ = −R_x·d`, compare `g_wᵀ·w′` with re-converged
frozen-primal central finite differences (D1/D2/D3), and only then apply the
B-final pressure-drop gradient Gate (sign + amplitude <5–10% + step-size
plateau). This is a separate, user-authorized task.
