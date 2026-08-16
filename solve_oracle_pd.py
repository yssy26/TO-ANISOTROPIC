#!/usr/bin/env python3
"""
Exact pressure-drop J^T oracle: direct sparse solve (SuperLU) + iterative
refinement to push the linear residual toward machine precision.

The frozen-RANS J^T is a saddle-point system (ill-conditioned), so a single
SuperLU solve leaves a residual ~1e-5.  Iterative refinement reuses the LU
factors and drives ||A x - b||/||b|| down to ~1e-13 (or lower).

Ordering: matrix is BLOCK-MAJOR (U-block rows 0..3N-1 = 3*celli+c, then
P-block rows 3N..4N-1 = 3N+celli); RHS export is CELL-MAJOR (per celli:
Ux,Uy,Uz,P).  We reorder the RHS to block-major.  Solution is emitted in
block-major order (what the OpenFOAM import path expects).
"""
import sys, os
import numpy as np
import scipy.io
import scipy.sparse as sp
import scipy.sparse.linalg as spla

BASE = r"\\wsl.localhost\Ubuntu-20.04\home\ys\b2_case_smoke"
MAT = os.path.join(BASE, "explicitJT.mtx")
RHS = os.path.join(BASE, "explicitRhs_pressureDrop.mtx")
SOL = os.path.join(BASE, "explicitSol_pressureDrop.mtx")

print("reading matrix ...", flush=True)
A = scipy.io.mmread(MAT).tocsr()
n = A.shape[0]
print("matrix:", A.shape, "nnz=", A.nnz, flush=True)

print("reading rhs (cell-major) ...", flush=True)
b_cell = np.asarray(scipy.io.mmread(RHS), dtype=np.float64).ravel()
N = n // 4

# reorder rhs cell-major -> block-major
b = np.zeros(n, dtype=np.float64)
for celli in range(N):
    for c in range(3):
        b[3 * celli + c] = b_cell[4 * celli + c]
    b[3 * N + celli] = b_cell[4 * celli + 3]
bnorm = np.linalg.norm(b)
print("||b|| =", bnorm, flush=True)

print("factorizing with SuperLU (splu) ...", flush=True)
lu = spla.splu(A)
print("initial solve ...", flush=True)
x = lu.solve(b)
print("||x|| =", np.linalg.norm(x), flush=True)

def relres(x):
    r = A @ x - b
    return np.linalg.norm(r) / max(bnorm, 1e-300)

print("iterative refinement ...", flush=True)
rr = relres(x)
print("  iter 0 relres=%.6e" % rr, flush=True)
for it in range(1, 25):
    if rr < 1e-13:
        break
    r = b - A @ x
    dx = lu.solve(r)
    x = x + dx
    rr = relres(x)
    print("  iter %d relres=%.6e" % (it, rr), flush=True)

print("final linear residual ||A x - b|| / ||b|| =", rr, flush=True)

print("writing solution (block-major) ...", flush=True)
with open(SOL, "w") as f:
    for v in x:
        f.write("%.17g\n" % v)
print("wrote", SOL, flush=True)
print("DONE relres=%.6e" % rr, flush=True)
