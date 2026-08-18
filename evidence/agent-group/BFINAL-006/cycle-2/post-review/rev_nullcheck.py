#!/usr/bin/env python3
"""
Post-Reviewer INDEPENDENT check (cheap, no LU): re-verify identity row,
build J=(JT)^T re-pinned at pRefRow, and confirm the T1 solution is a
numerical null vector (||J n||/||n|| ~ 1e-19) plus its spatial structure
and RHS orthogonality.  Read-only.
"""
import numpy as np, scipy.io, scipy.sparse as sp

CASE = "/home/ys/dsH/b2_case_smoke"
EVID = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-006/cycle-2"
N = 33600; NV = 3 * N; NUNK = 4 * N; PREF = NV  # 100800

M = scipy.io.mmread(CASE + "/explicitJT.mtx").tocsr()
print("loaded JT %dx%d nnz=%d" % (M.shape[0], M.shape[1], M.nnz))

# --- re-verify identity row in the pressure block (independent of s1_gate.py) ---
ident = []
for r in range(NV, NUNK):
    s, e = M.indptr[r], M.indptr[r + 1]
    nz = [(M.indices[j], M.data[j]) for j in range(s, e) if abs(M.data[j]) > 1e-12]
    if len(nz) == 1 and nz[0][0] == r and abs(nz[0][1] - 1.0) <= 1e-12:
        ident.append(r)
print("identity rows in pressure block:", ident)
row5600 = NV + 5600
s, e = M.indptr[row5600], M.indptr[row5600 + 1]
nz5600 = sum(1 for j in range(s, e) if abs(M.data[j]) > 1e-12)
print("row106400 stored=%d nz>1e-12=%d (must be physical, NOT identity)" % (e - s, nz5600))

# --- build J = M^T, re-pin row PREF ---
J = M.T.tocsr()
J = J.tolil(); J[PREF, :] = 0.0; J[PREF, PREF] = 1.0; J = J.tocsr()
print("J built nnz=%d" % J.nnz)

# --- null-vector check on the T1 solution (D1) ---
w = np.loadtxt(EVID + "/artifacts/wPrime_TAN_D1.mtx", dtype=np.float64)
wn = np.linalg.norm(w)
rw = J @ w
print("D1 wPrime: ||w||=%.6e  ||J w||=%.6e  ||J w||/||w||=%.6e" %
      (wn, np.linalg.norm(rw), np.linalg.norm(rw) / wn))
n = w / wn
Jn = J @ n
print("normalized n: ||J n||=%.6e  (== ||J w||/||w||, null-vector ratio)" % np.linalg.norm(Jn))
print("n[PREF]=%.3e (gauge pin, must be ~0)" % n[PREF])

# --- P-block spatial structure of the null vector ---
p = np.abs(w[NV:])
pmass = np.sum(p)
thr = np.percentile(p, 99)
big = p > thr
print("P-block |p|: total mass=%.3e  top-1%% (%d cells) mass fraction=%.4f  max=%.3e mean=%.3e"
      % (pmass, int(big.sum()), float(np.sum(p[big]) / pmass), p.max(), p.mean()))
print("n[PREF]=0 => null mode is NOT the constant-pressure mode (constant would have n[PREF]=1)")

# --- RHS orthogonality (consistent system check) ---
rxc = np.loadtxt(CASE + "/stageB6_rxc_analytic.mtx", dtype=np.float64)
for di, name in enumerate(["D1", "D2", "D3"]):
    off = di * 4 * N
    rxd = rxc[off:off + 4 * N]
    rhs = -rxd.copy(); rhs[PREF] = 0.0
    orth = float(np.dot(rhs, n) / np.linalg.norm(rhs))
    print("%s: (-R_x d)^T n / ||-R_x d|| = %.3e" % (name, orth))

# --- g_w orthogonality to null mode ---
gw_cell = np.asarray(scipy.io.mmread(CASE + "/explicitRhs_pressureDrop.mtx")).ravel()
gw = gw_cell[3::4].copy()  # P block (cell-major export)
print("g_w^T n = %.3e   sum(g_w,P)=%.3e" % (float(np.dot(gw, n[NV:])), float(np.sum(gw))))
print("REV_NULLCHECK_DONE")
