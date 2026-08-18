#!/usr/bin/env python3
"""
BFINAL-006 cycle-2 -- T1 numerical-floor diagnosis.

Purpose: explain why the PHYSICAL re-pinned tangent system J*wPrime=-R_x*d
cannot be solved by SuperLU to the T1 acceptance (true relative residual
<= 1e-9; observed floor 3.6e-4 / 4.0e-5 / 3.1e-3 with ||wPrime||~1e17-1e19).

Checks performed:
  A. sigma_min estimate of the physical re-pinned J via random solves
     (sigma_min ~= ||z||/||x|| with x = J^{-1} z).
  B. Solve the physical system for D1 rhs; report ||x||, residual.
  C. Solve the SCALED tangent system (Mtilde = Dout^-1 J Din, the exact
     momentum/area scaling used by the closed BFINAL-003 adjoint operator
     applyScaledDiscreteFlowJT), re-pinned the same way; report residual in
     PHYSICAL space and ||wPrime||.  This isolates whether the floor is
     purely the unbalanced physical scaling or a genuine matrix defect.
  D. Condition-number estimate from the scaled solve (for the record).
"""
import os, time
import numpy as np
import scipy.io
import scipy.sparse as sp
import scipy.sparse.linalg as spla

CASE = "/home/ys/dsH/b2_case_smoke"
N = 33600
NV = 3 * N
NUNK = 4 * N
PREF = NV           # 100800
M = 0.00786409044975
A = 2.47771624976e-07
MA = M / A

t0 = time.time()
def log(m): print(m, flush=True)

log("=== T1 numerical-floor diagnosis (physical vs scaled tangent system) ===")
Jt = scipy.io.mmread(os.path.join(CASE, "explicitJT.mtx")).tocsr()
J = Jt.T.tocsr()
del Jt
J = J.tolil(); J[PREF, :] = 0.0; J[PREF, PREF] = 1.0; J = J.tocsr()
log("physical re-pinned J built (nnz=%d, t=%.1fs)" % (J.nnz, time.time() - t0))

# ---- A. sigma_min estimate via random solves (x = J^{-1} z) ----
log("[A] estimating sigma_min of physical re-pinned J ...")
lu = spla.splu(J, permc_spec="COLAMD")
log("  splu done (t=%.1fs)" % (time.time() - t0))
sigmas = []
for k in range(3):
    rng = np.random.default_rng(1000 + k)
    z = rng.standard_normal(NUNK)
    x = lu.solve(z)
    sig = np.linalg.norm(z) / np.linalg.norm(x)
    sigmas.append(sig)
    log("  trial %d: ||z||=%.4e ||x||=%.4e  sigma_min~%.3e" % (k, np.linalg.norm(z), np.linalg.norm(x), sig))
log("  sigma_min estimate: min=%.3e mean=%.3e" % (min(sigmas), float(np.mean(sigmas))))
log("  ||J||_2 estimate ~ max row norm = %.3e -> cond ~ %.3e"
    % (np.sqrt(np.asarray(J.multiply(J).sum(axis=1)).max()),
       np.sqrt(np.asarray(J.multiply(J).sum(axis=1)).max()) / min(sigmas)))

# ---- B. physical-system solve for D1 ----
rxc = np.loadtxt(os.path.join(CASE, "stageB6_rxc_analytic.mtx"), dtype=np.float64)
d1 = rxc[0:4 * N]
rhs = -d1.copy()
rhs[PREF] = 0.0
wB = lu.solve(rhs)
resB = np.linalg.norm(J @ wB - rhs)
log("[B] physical solve D1: ||wPrime||=%.4e ||Jw-rhs||=%.4e relres=%.3e"
    % (np.linalg.norm(wB), resB, resB / np.linalg.norm(rhs)))

# ---- C. scaled-system solve ----
# row scale: U rows /M, P rows /A ; col scale: P cols *M/A
log("[C] building scaled tangent system Mtilde = Dout^-1 J Din ...")
Mtilde = J.copy()
Mtilde = Mtilde.tolil()
for r in range(NUNK):
    if r < NV:
        Mtilde[r, :] = Mtilde[r, :] / M
    else:
        Mtilde[r, :] = Mtilde[r, :] / A
Mtilde = Mtilde.tocsr()
# col scale P cols by M/A
Dcol = np.ones(NUNK)
Dcol[NV:] = MA
Mtilde = Mtilde @ sp.diags(Dcol)
Mtilde = Mtilde.tolil()
# re-pin row PREF (scaled space gauge  w'_P[PREF]=0 <=> z[PREF]=0)
Mtilde[PREF, :] = 0.0
Mtilde[PREF, PREF] = 1.0
Mtilde = Mtilde.tocsr()
log("  Mtilde built (nnz=%d, t=%.1fs)" % (Mtilde.nnz, time.time() - t0))
lu2 = spla.splu(Mtilde, permc_spec="COLAMD")
log("  splu(scaled) done (t=%.1fs)" % (time.time() - t0))

bscale = rhs.copy()
bscale[:NV] /= M
bscale[NV:] /= A
bscale[PREF] = 0.0
z = lu2.solve(bscale)
wC = z.copy()
wC[NV:] *= MA          # wPrime = Din z
resC = np.linalg.norm(J @ wC - rhs)
log("[C] scaled solve D1: ||z||=%.4e ||wPrime||=%.4e  PHYSICAL relres=%.3e  w'[PREF]=%.3e"
    % (np.linalg.norm(z), np.linalg.norm(wC), resC / np.linalg.norm(rhs), wC[PREF]))
# iterative refinement on the scaled system
rr = np.linalg.norm(Mtilde @ z - bscale) / np.linalg.norm(bscale)
log("  scaled-system relres (iter0)=%.3e" % rr)
for it in range(1, 30):
    if rr < 1e-13: break
    dz = lu2.solve(bscale - Mtilde @ z)
    z = z + dz
    rr = np.linalg.norm(Mtilde @ z - bscale) / np.linalg.norm(bscale)
log("  scaled-system relres after refinement=%.3e (iter %d)" % (rr, it))
wC2 = z.copy(); wC2[NV:] *= MA
resC2 = np.linalg.norm(J @ wC2 - rhs)
log("  physical-space residual with refined scaled solution: relres=%.3e  ||wPrime||=%.4e"
    % (resC2 / np.linalg.norm(rhs), np.linalg.norm(wC2)))

# ---- D. cond estimate of scaled system (sigma_min via random solve) ----
sig2 = []
for k in range(3):
    rng = np.random.default_rng(2000 + k)
    zz = rng.standard_normal(NUNK)
    xx = lu2.solve(zz)
    sig2.append(np.linalg.norm(zz) / np.linalg.norm(xx))
log("[D] scaled-system sigma_min estimate: min=%.3e -> cond~%.1e"
    % (min(sig2), np.sqrt(np.asarray(Mtilde.multiply(Mtilde).sum(axis=1)).max()) / min(sig2)))
log("done in %.1fs" % (time.time() - t0))
