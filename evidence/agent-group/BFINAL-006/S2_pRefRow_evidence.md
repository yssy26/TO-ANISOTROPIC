# BFINAL-006 S2 — NEW EVIDENCE: pRefRow is 100800 (P cell 0), NOT 106400 (P cell 5600)

Status: BLOCKED_BY_NEW_EVIDENCE (executor S2 stopped; see structured result).

The approved plan's S2 assumed `pRefRow = discretePIndex(pRefCell) = 3*33600 + 5600 = 106400`
(pRefCell=5600 from `system/fvSolution`).  Direct inspection of the actual exported
artifact contradicts this: the ONLY identity row of `explicitJT.mtx` is
**row 100800 = `discretePIndex(0)`**, i.e. the P row of cell 0.

## 1. Raw matrix evidence (explicitJT.mtx, sha256 3dfcc9ab…, b2_case_smoke)

```
identity rows (single nonzero 1.0 on diagonal, after scipy duplicate-sum): [100800]
expected P(0)  = 3*33600 +    0 = 100800
expected P(5600) = 3*33600 + 5600 = 106400
row 100800 entries: [(0,0.0),(1,0.0),(2,0.0),(3,0.0),(4,0.0),(5,0.0),(240,0.0),(241,0.0),
                     (242,0.0),(3362,0.0),(100800,1.0),(100801,0.0),(100880,0.0),(101920,0.0)]
row 106400 entries: 18 columns, 18 nonzero, physical continuity row
                    [(13440,6.97e-08),(13441,3.47e-11),(13442,1.224e-07),(16800,-1.25e-07),...]
```

- Row 100800 = the exported pin signature (all original columns zeroed, diagonal = 1.0).
- Row 106400 = the physical continuity row of cell 5600 (NOT pinned).
- The raw MatrixMarket file has 13,558,441 coordinate lines but only 6,209,576 unique
  entries (the C++ COO->CSR merge at solveDiscreteFlowAdjoint.H L2861-2870 does not merge
  duplicates because `csrRowPtr[r+1]==csrCol.size()` is stale during row processing;
  scipy sums duplicates -> identical operator content, different nnz count).

## 2. Root cause (verified in source + case)

OpenFOAM 7 `setRefCell` (findRefCell.C) assigns `pRefCell` only when
`field.needReference() || forceReference` is true:

```cpp
if (fieldRef.needReference() || forceReference) { ... refCelli = readLabel(...); ... }
```

`b2_case_smoke/0/p` has `outlet { type fixedValue; value uniform 0; }` (all other patches
zeroGradient) -> `p.needReference() == false` -> `setRefCell(p, simple.dict(), pRefCell,
pRefValue)` in `readTransportProperties.H:68` is a **no-op**: `label pRefCell = 0;` stays 0
and `scalar pRefValue = 0.0;` stays 0.

Consequences (all verified):
- The closed BFINAL-003 J^T export (solveDiscreteFlowAdjoint.H L2877-2899 pin +
  L3015-3026 export) pins `discretePIndex(pRefCell) = discretePIndex(0) = 100800`.
- The in-code operator pins the same row (`applyExplicitScaled` L2938
  `output[discretePIndex(pRefCell)] = input[...]`); run log confirms:
  `P-cell0: vP=0.999579529469 JTvP[P]=0.999579529469 MvP[P]=0.999579529469` (identity pin on P(0)).
- The primal NS.H:1025 `pEqn.setReference(pRefCell=0, pRefValue=0)` weakly biases p[0]
  (doubles the diagonal); the STRONG pressure reference is the outlet fixedValue p=0 BC.
- Converged p fields (0/p and 1/p) have NO cell at exactly 0: p[0]=74550.087,
  p[5599]=3.234, p[5600]=74809.058 (consistent with weak setReference(0,0) + outlet BC).

## 3. Impact on the approved S2 procedure

- The plan's S2 cross-check `pRefRow == 3*N+pRefCell == 106400` FAILS (actual 100800).
- The plan's gauge `wPrime_p[5600] = 0` does not match the closed Jacobian's reference
  convention (P(0)) nor the primal's effective reference (outlet BC + weak p[0] bias).
- After transpose, the exported matrix has column 100800 = e_100800 (cell 0's pressure is
  decoupled from momentum/continuity rows) and row 100800 is the corrupted
  (physical-column-with-diagonal-1) row; the physical continuity equation of cell 0 is
  destroyed in the file.  The tangent built from this file is the tangent of the closed
  PINNED operator with its P(0) reference — which is what BFINAL-003 closed, but the T2
  (state tangent) gauge must be handled deliberately (FD gauge = outlet BC + weak p[0]
  bias), and T3's D_TAN is gauge-invariant only w.r.t. the constant-pressure mode
  (g_w P-block sums to zero: +rho*area/(pMax*inletArea) inlet - ... outlet).

## 4. Not modified

- No source file, no case file, no forbidden module touched.
- `evidence/agent-group/BFINAL-006/tangent_solve.py` was written (read-only driver) and
  executed once; it stopped at the pRefRow cross-check (candidates empty because the
  pinned row carries explicit zeros; the identity row is row 100800).  No wPrime files
  were written.  scipy 1.17.1 was pip-installed into
  `evidence/agent-group/BFINAL-002/.venv_s4` (numpy 2.4.6 unchanged; `splu` import OK);
  install log: `scipy_install.log`.

## 5. Recommended re-plan direction (for Planner/Pre-Reviewer; NOT executed)

- Use the ACTUALLY detected pRefRow = 100800 (the closed Jacobian's reference), OR
  re-export the matrix with the physical row preserved; decide the T2 gauge protocol
  given the primal's effective reference is the outlet fixedValue BC (p_outlet=0), with
  p[0] weakly biased by setReference(0,0).
