#!/usr/bin/env python3
"""
BFINAL-006 S2 -- Offline T1 tangent solve (read-only, no C++ rebuild).

Gate T1 (task BFINAL-006 §8): assemble rhs = -R_x*d with the BFINAL-005 CLOSED
production R_x model (stageB6_rxc_analytic.mtx), solve the closed BFINAL-003
Jacobian J*wPrime = rhs to high accuracy, and report true relative residuals
split by U rows / P rows.

Conventions (LOCKED by BFINAL-003 / BFINAL-005, verified against source):
  * Matrix: explicitJT.mtx = MatrixMarket "matrix coordinate real general",
    134400 x 134400, physical J^T CSR.  The pRef ROW (discretePIndex(pRefCell)
    = 3*N + pRefCell = 106400, N=33600, pRefCell=5600, pRefValue=0) is pinned
    to identity in the exported J^T (solveDiscreteFlowAdjoint.H L2877-2899,
    export L3015-3026).  Ordering is BLOCK-MAJOR: U rows 0..3N-1 =
    3*celli+cmpt, then P rows 3N..4N-1 = 3N+celli.
  * Tangent: J*wPrime = -R_x*d  with J = dR/dw (R_U = -UEqn.residual(),
    R_P = div(phi)), R_x*d = [anRUd(3N); anRPd(N)] block-major physical
    (stageB6_rxc_analytic.mtx layout: per direction di, block
    [anRUd(3N); anRPd(N)]).
  * Gauge: transposing the pinned J^T gives J with pRef COLUMN pinned to
    identity.  J's pRef ROW is the natural continuity row; it is RE-PINNED to
    identity and rhs[pRefRow]=0 so that wPrime_p[pRefCell]=0, matching the
    primal pEqn.setReference(pRefCell, pRefValue=0) gauge.  No arbitrary
    pressure mean is subtracted anywhere.
  * g_w: explicitRhs_pressureDrop.mtx (cell-major [Ux,Uy,Uz,P] per cell); its
    P block IS the physical pressureConstraintDerivative (direct g_w).
    D_TAN = g_w^T*wPrime is computed as an informational T3 output (allowed:
    R_x + J + forward solve + direct g_w; NO adjoint quantity).

Read-only: never modifies the case, the matrix, or any source file.
"""
import os
import sys
import time
import numpy as np
import scipy.io
import scipy.sparse as sp
import scipy.sparse.linalg as spla

CASE = "/home/ys/dsH/b2_case_smoke"
EVID = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-006"
MAT = os.path.join(CASE, "explicitJT.mtx")
RXC_AN = os.path.join(CASE, "stageB6_rxc_analytic.mtx")
DIRS = os.path.join(CASE, "stageB6_dirs.mtx")
RHS_GW = os.path.join(CASE, "explicitRhs_pressureDrop.mtx")

N = 33600
NV = 3 * N          # discreteNVelocity
NUNK = 4 * N        # 134400
PREF_CELL = 5600
PREF_ROW = NV + PREF_CELL   # 100800 + 5600 = 106400
DIR_NAMES = ["D1", "D2", "D3"]

t0 = time.time()

def log(msg):
    print(msg, flush=True)

# ---------------------------------------------------------------------------
# 1. Read J^T (physical CSR, pRef row pinned).
# ---------------------------------------------------------------------------
log("[1] reading explicitJT.mtx (%.1f MB) ..." % (os.path.getsize(MAT) / 1e6))
Jt = scipy.io.mmread(MAT).tocsr()
n = Jt.shape[0]
assert n == NUNK, "unexpected size %d" % n
log("    matrix %dx%d nnz=%d  (t=%.1fs)" % (n, n, Jt.nnz, time.time() - t0))

# ---------------------------------------------------------------------------
# 2. Identify the pRef row: the unique row of J^T equal to the identity row
#    (single 1.0 on the diagonal).  Cross-check == 3*N + pRefCell == 106400.
# ---------------------------------------------------------------------------
log("[2] locating pRef row (unique identity row of J^T) ...")
nrow = np.diff(Jt.indptr)
cand = []
for r in range(n):
    if nrow[r] == 1:
        c0 = Jt.indices[Jt.indptr[r]]
        v0 = Jt.data[Jt.indptr[r]]
        if c0 == r and abs(v0 - 1.0) < 1e-14:
            cand.append(r)
log("    identity-row candidates: %s" % cand)
if len(cand) == 1:
    pRefRow = cand[0]
elif PREF_ROW in cand:
    pRefRow = PREF_ROW
    log("    WARNING: %d identity-row candidates; using expected pRefRow %d"
        % (len(cand), PREF_ROW))
else:
    raise SystemExit("BLOCKED_BY_NEW_EVIDENCE: pRef identity row not found "
                     "(candidates=%s, expected=%d)" % (cand, PREF_ROW))
assert pRefRow == PREF_ROW, "pRefRow %d != expected %d" % (pRefRow, PREF_ROW)
# explicit check of the pinned row content
pr_s = Jt.indptr[pRefRow]; pr_e = Jt.indptr[pRefRow + 1]
log("    pRefRow=%d (=3*N+pRefCell=%d) row nnz=%d entries=%s"
    % (pRefRow, PREF_ROW, pr_e - pr_s,
       list(zip(Jt.indices[pr_s:pr_e].tolist(), Jt.data[pr_s:pr_e].tolist()))))

# ---------------------------------------------------------------------------
# 3. Transpose J^T -> J via scatter: J[col, row] += val  for every (row,col,val).
#    After this, J has the pRef COLUMN pinned to identity.
# ---------------------------------------------------------------------------
log("[3] transposing J^T -> J (vectorized scatter) ...")
rows = np.repeat(np.arange(n, dtype=np.int64), nrow)
J = sp.coo_matrix((Jt.data, (Jt.indices.astype(np.int64), rows)),
                  shape=(n, n)).tocsr()
del Jt
log("    J nnz=%d  (t=%.1fs)" % (J.nnz, time.time() - t0))
# verify the pRef column is identity (inherited from the J^T row pin)
col = J[:, pRefRow].toarray().ravel()
col_ok = (np.count_nonzero(col) == 1) and (abs(col[pRefRow] - 1.0) < 1e-14)
log("    J[:,pRefRow] identity-column check: %s" % col_ok)
assert col_ok, "J pRef column is not identity"

# ---------------------------------------------------------------------------
# 4. Re-pin J's pRef ROW to identity (gauge wPrime_p[pRefCell]=0, matching
#    primal pRefValue=0); the solve then sets wPrime[pRefRow]=rhs[pRefRow]=0.
# ---------------------------------------------------------------------------
log("[4] re-pinning J pRef row to identity ...")
J = J.tolil()
J[pRefRow, :] = 0.0
J[pRefRow, pRefRow] = 1.0
J = J.tocsr()

# ---------------------------------------------------------------------------
# 5. Factorize once (SuperLU), iterative refinement per direction.
# ---------------------------------------------------------------------------
log("[5] SuperLU factorization of pinned J (n=%d) ..." % n)
lu = spla.splu(J)
log("    factorized (t=%.1fs)" % (time.time() - t0))

def solve_refine(rhs, label):
    bnorm = np.linalg.norm(rhs)
    x = lu.solve(rhs)
    rr = np.linalg.norm(J @ x - rhs) / max(bnorm, 1e-300)
    log("    %s: iter 0 relres=%.6e" % (label, rr))
    for it in range(1, 30):
        if rr < 1e-13:
            break
        dx = lu.solve(rhs - J @ x)
        x = x + dx
        rr = np.linalg.norm(J @ x - rhs) / max(bnorm, 1e-300)
        if it % 5 == 0 or rr < 1e-12:
            log("    %s: iter %d relres=%.6e" % (label, it, rr))
    return x, rr

# ---------------------------------------------------------------------------
# 6. Load directions (stageB6_dirs.mtx: 3 x N bare list) and analytic R_x*d
#    (stageB6_rxc_analytic.mtx: per direction [anRUd(3N); anRPd(N)]).
# ---------------------------------------------------------------------------
log("[6] loading directions and R_x*d ...")
dirs_raw = np.loadtxt(DIRS, dtype=np.float64)
assert dirs_raw.size == 3 * N, "dirs size %d" % dirs_raw.size
dirs = [dirs_raw[di * N:(di + 1) * N] for di in range(3)]

rxc = np.loadtxt(RXC_AN, dtype=np.float64)
assert rxc.size == 3 * 4 * N, "rxc size %d" % rxc.size

# g_w (physical direct pressure functional derivative): explicitRhs P block.
gw_cell = np.asarray(scipy.io.mmread(RHS_GW), dtype=np.float64).ravel()
assert gw_cell.size == 4 * N
gw = gw_cell[3::4].copy()          # P block of the cell-major export
gwU = np.concatenate([gw_cell[0::4], gw_cell[1::4], gw_cell[2::4]])
log("    ||g_w(P)||=% .6e  ||g_w(U)||=% .6e (U block should be ~0)"
    % (np.linalg.norm(gw), np.linalg.norm(gwU)))

# ---------------------------------------------------------------------------
# 7. Per direction: rhs = -R_x*d, solve J*wPrime=rhs, report residuals.
# ---------------------------------------------------------------------------
results = {}
for di in range(3):
    off = di * 4 * N
    anRUd = rxc[off:off + 3 * N]
    anRPd = rxc[off + 3 * N:off + 4 * N]
    d = dirs[di]
    rxd = np.concatenate([anRUd, anRPd])
    rhs = -rxd.copy()
    rhs[PREF_ROW] = 0.0
    rnorm = np.linalg.norm(rhs)
    rnormU = np.linalg.norm(rhs[:3 * N])
    rnormP = np.linalg.norm(rhs[3 * N:])
    log("[7] direction %s: ||R_x*d||=% .6e  ||rhs(post-pin)||=% .6e"
        % (DIR_NAMES[di], np.linalg.norm(rxd), rnorm))
    wP, rr = solve_refine(rhs, DIR_NAMES[di])
    res = J @ wP - rhs
    resU = np.linalg.norm(res[:3 * N])
    resP = np.linalg.norm(res[3 * N:])
    true_rel = np.linalg.norm(res) / max(rnorm, 1e-300)
    relU = resU / max(rnormU, 1e-300)
    relP = resP / max(rnormP, 1e-300)
    log("    %s: ||J*wPrime-rhs||=% .6e  trueRelres=% .6e  "
        "U-rowRel=% .6e  P-rowRel=% .6e  ||wPrime||=% .6e"
        % (DIR_NAMES[di], np.linalg.norm(res), true_rel, relU, relP,
           np.linalg.norm(wP)))
    # wPrime gauge check
    assert abs(wP[PREF_ROW]) < 1e-12, "wPrime[pRef] not zero"
    # D_TAN = g_w^T * wPrime (informational T3 output; direct g_w, no adjoint)
    pPrime = wP[3 * N:]
    dtan = float(np.dot(gw, pPrime))
    log("    %s: D_TAN = g_w^T*wPrime = % .8e" % (DIR_NAMES[di], dtan))
    # write wPrime_TAN (block-major, full 4N)
    out = os.path.join(EVID, "wPrime_TAN_%s.mtx" % DIR_NAMES[di])
    with open(out, "w") as f:
        for v in wP:
            f.write("%.17g\n" % v)
    log("    wrote %s (%.0f bytes)" % (out, os.path.getsize(out)))
    results[DIR_NAMES[di]] = dict(
        rxd_norm=float(np.linalg.norm(rxd)),
        rhs_norm=rnorm, rhs_normU=rnormU, rhs_normP=rnormP,
        abs_res=float(np.linalg.norm(res)),
        true_relres=true_rel, relU=relU, relP=relP,
        wnorm=float(np.linalg.norm(wP)), dtan=dtan)

log("[8] done in %.1fs" % (time.time() - t0))
print("\n=== SUMMARY (T1 tangent solve) ===")
hdr = "dir     ||R_x*d||       ||Jw-rhs||      trueRel   Urel      Prel      D_TAN"
print(hdr)
for di in range(3):
    r = results[DIR_NAMES[di]]
    print("%-6s % .6e  % .6e  % .3e  % .3e  % .3e  % .8e"
          % (DIR_NAMES[di], r["rxd_norm"], r["abs_res"], r["true_relres"],
             r["relU"], r["relP"], r["dtan"]))
