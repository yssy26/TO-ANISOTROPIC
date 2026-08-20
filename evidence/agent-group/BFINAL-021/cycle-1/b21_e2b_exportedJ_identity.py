#!/usr/bin/env python3
"""BFINAL-021 E2b: P-row identity with the EXACT exported J rows (zero-code).

Improvement over b21_e2: instead of the python-rebuilt PsiU/PsiP (which is
only a ~0.15% approximation of the operator), use the exported explicitJT
matrix directly: J = M^T, so (J w)[P-block] = (M.T @ w)[P-block].  This
removes all rebuild-representation error from the identity test.

Per direction and eps:
  r_P(w_state) = (J w_state)[P] + rxd_P      (exported J, exported rxd_P)
  cos(div part, -rxd_P), |r_P|/|rxd_P|, eps-dependence
Controls:
  r_P(w') for the cached b15 operator solve (quantifies the b15 pinning
  artifact), and the decomposition r_P(w_state) = J_P (w_state - w').
"""
import time, json
import numpy as np
import scipy.io
import scipy.sparse as sp

EXPORT_JT = "/home/ys/dsH/b8_verify_diag/explicitJT.mtx"
E = "/home/ys/dsH/b15_export"
STATE_CASE = "/home/ys/dsH/b13_probe3"
WPC = ("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/"
       "BFINAL-015/cycle-1/wprime_cache.npz")
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-021/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N; P0 = NV
t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-021 E2b: identity via exported J ===")

M = scipy.io.mmread(EXPORT_JT).tocsr()
M.sort_indices()
MT = M.T.tocsr()
log("M loaded nnz=%d (t=%.1fs)" % (M.nnz, time.time()-t0))

rxc = np.loadtxt(E + "/stageB6_rxc_analytic.mtx", dtype=np.float64)

def load_state(dname, tag, sign):
    Uu = np.loadtxt(STATE_CASE + "/stageB2/wstate_%s_h%s_%s_U.mtx" % (dname, tag, sign))
    pp = np.loadtxt(STATE_CASE + "/stageB2/wstate_%s_h%s_%s_p.mtx" % (dname, tag, sign))
    wst = np.zeros(NUNK)
    wst[0:NV:3] = Uu[:, 0]; wst[1:NV:3] = Uu[:, 1]; wst[2:NV:3] = Uu[:, 2]
    wst[P0:] = pp
    return wst

cz = np.load(WPC)
wprime = {k: cz[k] for k in cz.files}
ladder = {"D1": ["0.0003", "0.001", "0.003"],
          "D2": ["0.0003", "0.001"],
          "D3": ["0.0003", "0.001"]}

def metrics(w, rxdP):
    Jw_P = (MT @ w)[P0:]
    dv = Jw_P                        # = J_PU dU + J_PP dp (unclosed tangent div)
    r_P = dv + rxdP
    cosv = float(dv @ (-rxdP)/(np.linalg.norm(dv)*np.linalg.norm(rxdP)))
    return {"cos": cosv,
            "|rP|/|rxdP|": float(np.linalg.norm(r_P)/np.linalg.norm(rxdP)),
            "|J_P w|/|rxdP|": float(np.linalg.norm(dv)/np.linalg.norm(rxdP))}

res = {"ladder": {}, "wprime_control": {}, "decomp": {}}
log("\n[E2b] identity ladder (exported J P-rows + exported rxd_P):")
for nm in ["D1", "D2", "D3"]:
    di = ["D1", "D2", "D3"].index(nm)
    rxdP = rxc[di*NUNK + P0: (di+1)*NUNK]
    res["ladder"][nm] = {}
    for tag in ladder[nm]:
        h = float(tag)
        wt = (load_state(nm, tag, "p") - load_state(nm, tag, "m"))/(2*h)
        m = metrics(wt, rxdP)
        res["ladder"][nm][tag] = m
        log("  %s h=%s: cos=%.4f |r_P|/|rxdP|=%.4f |J_P w|/|rxdP|=%.4f"
            % (nm, tag, m["cos"], m["|rP|/|rxdP|"], m["|J_P w|/|rxdP|"]))

log("\n[control] r_P(w') (b15 cached solve):")
for nm in ["D1", "D2", "D3"]:
    di = ["D1", "D2", "D3"].index(nm)
    rxdP = rxc[di*NUNK + P0: (di+1)*NUNK]
    m = metrics(wprime[nm], rxdP)
    res["wprime_control"][nm] = m
    log("    %s: |r_P(w')|/|rxdP|=%.4e cos=%.6f"
        % (nm, m["|rP|/|rxdP|"], m["cos"]))

log("\n[decomp] r_P(w_state) = J_P (w_state - w') (h=1e-3):")
for nm in ["D1", "D2", "D3"]:
    di = ["D1", "D2", "D3"].index(nm)
    rxdP = rxc[di*NUNK + P0: (di+1)*NUNK]
    wt = (load_state(nm, "0.001", "p") - load_state(nm, "0.001", "m"))/2e-3
    d = wt - wprime[nm]
    Jd_P = (MT @ d)[P0:]
    r_P = (MT @ wt)[P0:] + rxdP
    res["decomp"][nm] = {
        "|J_P dw|/|rxdP|": float(np.linalg.norm(Jd_P)/np.linalg.norm(rxdP)),
        "|r_P|/|rxdP|": float(np.linalg.norm(r_P)/np.linalg.norm(rxdP)),
        "|dw_U|/|w_U|": float(np.linalg.norm(d[:NV])/np.linalg.norm(wt[:NV])),
        "|dw_P|/|w_P|": float(np.linalg.norm(d[P0:])/np.linalg.norm(wt[P0:])),
        "resid_vs_identity": float(np.linalg.norm(Jd_P - r_P)
                                   /np.linalg.norm(r_P))}
    log("    %s: |J_P dw|/|rxdP|=%.4f  |r_P|/|rxdP|=%.4f  identity-resid=%.3e"
        "  (|dw_U|/|w_U|=%.3f |dw_P|/|w_P|=%.3f)"
        % (nm, res["decomp"][nm]["|J_P dw|/|rxdP|"],
           res["decomp"][nm]["|r_P|/|rxdP|"],
           res["decomp"][nm]["resid_vs_identity"],
           res["decomp"][nm]["|dw_U|/|w_U|"],
           res["decomp"][nm]["|dw_P|/|w_P|"]))

json.dump(res, open(OUT + "/b21_e2b_exportedJ_identity.json", "w"), indent=1)
log("\nDONE t=%.1fs" % (time.time()-t0))
