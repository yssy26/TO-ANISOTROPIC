#!/usr/bin/env python3
"""
Post-Reviewer INDEPENDENT tangent-solve rerun (gold standard, ~35 min):
fresh splu of the re-pinned J, fresh random sigma_min estimates, and a
fresh solve of J wPrime = -R_x d for D1/D2/D3 with iterative refinement.
Reports the true relative residual floor to test the T1 acceptance
(<= 1e-9).  Read-only.
"""
import os, time
import numpy as np, scipy.io, scipy.sparse as sp, scipy.sparse.linalg as spla

CASE = "/home/ys/dsH/b2_case_smoke"
N = 33600; NV = 3 * N; NUNK = 4 * N; PREF = NV

t0 = time.time()
def log(m): print(m, flush=True)

log("=== Post-Reviewer independent tangent-solve rerun ===")
M = scipy.io.mmread(os.path.join(CASE, "explicitJT.mtx")).tocsr()
J = M.T.tocsr(); del M
J = J.tolil(); J[PREF, :] = 0.0; J[PREF, PREF] = 1.0; J = J.tocsr()
log("re-pinned J built nnz=%d (t=%.1fs)" % (J.nnz, time.time() - t0))

log("[LU] splu(J, COLAMD) ...")
lu = spla.splu(J, permc_spec="COLAMD")
log("splu done (t=%.1fs)" % (time.time() - t0))

log("[sigma_min] 2 fresh random solves")
sig = []
for k in range(2):
    rng = np.random.default_rng(424242 + k)
    z = rng.standard_normal(NUNK)
    x = lu.solve(z)
    s = np.linalg.norm(z) / np.linalg.norm(x)
    sig.append(s)
    log("  trial %d: sigma_min~%.3e (||x||=%.3e)" % (k, s, np.linalg.norm(x)))
log("  sigma_min: min=%.3e" % min(sig))

rxc = np.loadtxt(os.path.join(CASE, "stageB6_rxc_analytic.mtx"), dtype=np.float64)
log("[solves] D1/D2/D3 with iterative refinement (30 iters)")
for di, name in enumerate(["D1", "D2", "D3"]):
    off = di * 4 * N
    rxd = rxc[off:off + 4 * N]
    rhs = -rxd.copy(); rhs[PREF] = 0.0
    rnorm = np.linalg.norm(rhs)
    x = lu.solve(rhs)
    rr = np.linalg.norm(J @ x - rhs) / rnorm
    for it in range(1, 30):
        if rr < 1e-13: break
        x = x + lu.solve(rhs - J @ x)
        rr = np.linalg.norm(J @ x - rhs) / rnorm
    res = J @ x - rhs
    relU = np.linalg.norm(res[:3 * N]) / np.linalg.norm(rhs[:3 * N])
    relP = np.linalg.norm(res[3 * N:]) / np.linalg.norm(rhs[3 * N:])
    log("  %s: ||rhs||=%.6e trueRelRes=%.6e Urel=%.6e Prel=%.6e ||wPrime||=%.6e w[PREF]=%.3e"
        % (name, rnorm, rr, relU, relP, np.linalg.norm(x), x[PREF]))
log("REV_TANGENT_RERUN_DONE (t=%.1fs)" % (time.time() - t0))
