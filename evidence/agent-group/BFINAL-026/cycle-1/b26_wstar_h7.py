#!/usr/bin/env python3
"""BFINAL-026 instrument: w_true on the H7 operator (coordinator-executed salvage).

Replicates b21_wstar_clean.py exactly, with three changes:
  1. M = explicitJT_H7.mtx (b24_diag export, B24 ExplicitJToracle 5.48e-16) --
     the b8_verify_diag explicitJT.mtx is PRE-H7 and must not be used for w_true.
  2. splu ordering MMD_AT_PLUS_A + SymmetricMode (structurally symmetric
     saddle-point pattern) with fallback to the b21 default (COLAMD).
  3. State-response validation against b25_qgate wstate exports (Q-functional
     state, bit-identical fields to b8_verify_diag per NOTEBOOK).

Outputs (BFINAL-026/cycle-1/): b26_wstar_h7.npz + b26_wstar_h7.json
Closure acceptance: |r_P(w*)|/|rxdP| ~ 1e-11 (b21 level); BLOCK the round if > 1e-9.
"""
import os, time, json
import numpy as np
import scipy.io
import scipy.sparse as sp
import scipy.sparse.linalg as spla

EXPORT_JT = "/home/ys/dsH/b24_diag/explicitJT_H7.mtx"
RXC = "/home/ys/dsH/b8_verify_diag/stageB6_rxc_analytic.mtx"
STATE_CASE = "/home/ys/dsH/b25_qgate"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-026/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N; P0 = NV
t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-026: clean unpinned w* on the H7 operator ===")

M = scipy.io.mmread(EXPORT_JT).tocsc()
M.sort_indices()
log("M loaded nnz=%d shape=%s (t=%.1fs)" % (M.nnz, M.shape, time.time()-t0))

# structural symmetry sanity (informational; old scipy's sparse != returns int)
with open(EXPORT_JT) as fh:
    for _ in range(2):
        hdr = fh.readline()
        if not hdr.startswith("%%"):
            declared = int(hdr.split()[2])
            break
log("declared coordinate entries: %d -> unique CSC nnz: %d (duplicates summed)"
    % (declared, M.nnz))
asy = (abs(M) != abs(M.T))
patt_asym = asy.nnz if hasattr(asy, "nnz") else int(asy)
log("pattern asymmetry entries: %d / %d (%.2f%%)"
    % (patt_asym, M.nnz, 100.0*patt_asym/M.nnz))

lu = None
try:
    lu = spla.splu(M, permc_spec="MMD_AT_PLUS_A",
                   options=dict(SymmetricMode=True, diag_pivot_thresh=0.0))
    log("splu MMD_AT_PLUS_A+SymmetricMode done (t=%.1fs)" % (time.time()-t0))
except Exception as e:
    log("MMD splu failed (%s); falling back to default COLAMD" % e)
    lu = spla.splu(M)
    log("splu default done (t=%.1fs)" % (time.time()-t0))

rxc = np.loadtxt(RXC, dtype=np.float64)
log("rxc loaded %s (t=%.1fs)" % (rxc.shape, time.time()-t0))

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
    Jw = (M.T @ w)
    r_P = Jw[P0:] + rxd[P0:]
    entry = {"|r_P(w*)|/|rxdP|": float(np.linalg.norm(r_P)/np.linalg.norm(rxd[P0:]))}
    try:
        wt = (load_state(nm, "0.001", "p") - load_state(nm, "0.001", "m"))/2e-3
        d = wt - w
        entry.update({
            "cos_w_wt": float(w @ wt/(np.linalg.norm(w)*np.linalg.norm(wt))),
            "relL2": float(np.linalg.norm(d)/np.linalg.norm(wt)),
            "relU": float(np.linalg.norm(d[:NV])/np.linalg.norm(wt[:NV])),
            "relP": float(np.linalg.norm(d[P0:])/np.linalg.norm(wt[P0:])),
        })
    except IOError as e:
        log("  %s: wstate files missing (%s) -- closure only" % (nm, e))
    res[nm] = entry
    log("  %s: |r_P(w*)|/|rxdP|=%.3e  cos=%.6f relL2=%.4f (relU=%.4f relP=%.4f)"
        % (nm, entry["|r_P(w*)|/|rxdP|"], entry.get("cos_w_wt", float("nan")),
           entry.get("relL2", float("nan")), entry.get("relU", float("nan")),
           entry.get("relP", float("nan"))))
    log("  (t=%.1fs)" % (time.time()-t0))

np.savez(OUT + "/b26_wstar_h7.npz", **wstar)
json.dump(res, open(OUT + "/b26_wstar_h7.json", "w"), indent=1)
log("DONE t=%.1fs" % (time.time()-t0))
