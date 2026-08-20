#!/usr/bin/env python3
"""BFINAL-021 E2c: clean unpinned operator solve w* = J^{-1}(-rxd).

The b15 wprime_cache was produced with a pinned row (B16 documented the
pinning artifact); the control r_P(w') != 0 confirms contamination.  Here:
one unpinned splu of the exported M (= J^T), three solves for -rxd, plus
validation r_P(w*) ~ 0 and clean |w_state - w*| block metrics at h=1e-3.
~25 min runtime.
"""
import time, json
import numpy as np
import scipy.io
import scipy.sparse as sp
import scipy.sparse.linalg as spla

EXPORT_JT = "/home/ys/dsH/b8_verify_diag/explicitJT.mtx"
E = "/home/ys/dsH/b15_export"
STATE_CASE = "/home/ys/dsH/b13_probe3"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-021/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N; P0 = NV
t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-021 E2c: clean unpinned w* ===")

M = scipy.io.mmread(EXPORT_JT).tocsc()
M.sort_indices()
log("M loaded nnz=%d (t=%.1fs)" % (M.nnz, time.time()-t0))
lu = spla.splu(M)
log("splu done (t=%.1fs)" % (time.time()-t0))

rxc = np.loadtxt(E + "/stageB6_rxc_analytic.mtx", dtype=np.float64)

def load_state(dname, tag, sign):
    Uu = np.loadtxt(STATE_CASE + "/stageB2/wstate_%s_h%s_%s_U.mtx" % (dname, tag, sign))
    pp = np.loadtxt(STATE_CASE + "/stageB2/wstate_%s_h%s_%s_p.mtx" % (dname, tag, sign))
    wst = np.zeros(NUNK)
    wst[0:NV:3] = Uu[:, 0]; wst[1:NV:3] = Uu[:, 1]; wst[2:NV:3] = Uu[:, 2]
    wst[P0:] = pp
    return wst

wstar = {}
res = {}
for di, nm in enumerate(["D1", "D2", "D3"]):
    rxd = rxc[di*NUNK:(di+1)*NUNK]
    w = lu.solve(-rxd, trans="T")     # M^T w = -rxd  <=>  J w = -rxd
    wstar[nm] = w
    # validation: r_P(w*) = (J w*)[P] + rxd_P
    Jw = (M.T @ w)
    r_P = Jw[P0:] + rxd[P0:]
    # clean state mismatch at h=1e-3
    wt = (load_state(nm, "0.001", "p") - load_state(nm, "0.001", "m"))/2e-3
    d = wt - w
    res[nm] = {
        "|r_P(w*)|/|rxdP|": float(np.linalg.norm(r_P)/np.linalg.norm(rxd[P0:])),
        "cos_w_wt": float(w @ wt/(np.linalg.norm(w)*np.linalg.norm(wt))),
        "relL2": float(np.linalg.norm(d)/np.linalg.norm(wt)),
        "relU": float(np.linalg.norm(d[:NV])/np.linalg.norm(wt[:NV])),
        "relP": float(np.linalg.norm(d[P0:])/np.linalg.norm(wt[P0:])),
        "|J_P dw|/|rxdP|": float(np.linalg.norm((M.T @ d)[P0:])
                                 / np.linalg.norm(rxd[P0:])),
    }
    log("  %s: |r_P(w*)|/|rxdP|=%.3e  cos(w*,w_state)=%.6f relL2=%.4f "
        "(relU=%.4f relP=%.4f) |J_P dw|/|rxdP|=%.4f"
        % (nm, res[nm]["|r_P(w*)|/|rxdP|"], res[nm]["cos_w_wt"],
           res[nm]["relL2"], res[nm]["relU"], res[nm]["relP"],
           res[nm]["|J_P dw|/|rxdP|"]))

np.savez(OUT + "/b21_wstar_cache.npz", **wstar)
json.dump(res, open(OUT + "/b21_wstar_clean.json", "w"), indent=1)
log("DONE t=%.1fs" % (time.time()-t0))
