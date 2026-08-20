#!/usr/bin/env python3
"""
BFINAL-015 T2/T3 offline — tangent (w') vs state-FD (w_true) comparison.

Pure measurement (no source/case modification beyond reading exports).

  w'      : solves J_phys w' = -R_x d   (splu + iterative refinement)
            J_phys = transpose(explicitJT.mtx)   [physical, COMPLETE J^T]
            R_x d  = stageB6_rxc_analytic.mtx    [raw-x D1/D2/D3, full chain]
  w_true  : (U+ - U-)/(2h), (p+ - p-)/(2h) from the B2 state-export run
            (stageB15StateExport, x +/- h*D full-chain re-converged states)
  b       : explicitRhs_pressureDrop.mtx / explicitRhs_thermalCoupling.mtx
            [physical objective sources, P1-verified]

Checks:
  A. export integrity cross-check: b^T w' vs the B2-run ADJ projections
     (gDP D1 = -5.4050587339 etc. from BFINAL-012) -> linearization-point
     equivalence measure;
  B. w' vs w_true: cos, relL2, U/P block split, top-mismatch cells;
  C. double-objective attribution: b_gDP^T w' / b_gDP^T w_true (expect the
     end-to-end 2.11-2.17 band) and the same for b_J (expect the J misfit);
  D. eps ladder (3e-4, 1e-3, 3e-3) stability of w_true.
"""
import os, sys, time, json, glob
import numpy as np
import scipy.io, scipy.sparse as sp, scipy.sparse.linalg as spla

EXPORT_CASE = "/home/ys/dsH/b15_export"
# The case's discreteExplicitMatrixFile is an absolute path inherited from
# the b8_verify_diag lineage, so today's explicitJT.mtx landed there:
EXPORT_JT   = "/home/ys/dsH/b8_verify_diag/explicitJT.mtx"
STATE_CASE  = "/home/ys/dsH/b13_probe3"       # B2 + stageB15StateExport run
CYCLE = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-015/cycle-1"

N = 33600
NV = 3 * N
NUNK = 4 * N
PREF = NV
P0 = NV

# BFINAL-012 run A reference values (same case family, formal step h=1e-3)
REF = {
    "gDP": {"ADJ": {"D1": -5.4050587339, "D2": 1.10592698054, "D3": 14.0447142088},
             "FD":  {"D1": -2.53745977185, "D2": 0.523409452158, "D3": 6.47721554707}},
    "J":   {"ADJ": {"D1": -0.0159321971595, "D2": -0.0241647083211, "D3": 0.00316057765059},
             "FD":  {"D1": 0.0428559144263, "D2": -0.00346892642156, "D3": -0.000282658160322}},
}

t0 = time.time()
def log(m): print(m, flush=True)

log("=== BFINAL-015 T2/T3 offline ===")

# ---------------- 1. J from explicit export ----------------
M = scipy.io.mmread(EXPORT_JT).tocsr()
log("explicitJT loaded nnz=%d (t=%.1fs)" % (M.nnz, time.time() - t0))
J = M.T.tocsr()
del M
J = J.tolil()
J[PREF, :] = 0.0
J[PREF, PREF] = 1.0
J = J.tocsr()
log("J built (pinned P0 row for splu) nnz=%d" % J.nnz)

# ---------------- 2. tangent solves ----------------
rxc = np.loadtxt(EXPORT_CASE + "/stageB6_rxc_analytic.mtx", dtype=np.float64)
log("rxc_analytic loaded shape=%s" % (rxc.shape,))
lu = spla.splu(J, permc_spec="COLAMD")
log("splu done (t=%.1fs)" % (time.time() - t0))

wprime = {}
cache_f = CYCLE + "/wprime_cache.npz"
if os.path.exists(cache_f):
    cz = np.load(cache_f)
    wprime = {k: cz[k] for k in cz.files}
    log("w' cache loaded: %s" % list(wprime))
else:
    for di, name in enumerate(["D1", "D2", "D3"]):
        off = di * NUNK
        rxd = rxc[off:off + NUNK]
        rhs = -rxd.copy()
        rhs[PREF] = 0.0
        rnorm = np.linalg.norm(rhs)
        x = lu.solve(rhs)
        rr = np.linalg.norm(J @ x - rhs) / rnorm
        for it in range(60):
            if rr < 1e-12:
                break
            x = x + lu.solve(rhs - J @ x)
            rr = np.linalg.norm(J @ x - rhs) / rnorm
        log("  w' %s: ||rhs||=%.6e trueRelRes=%.3e ||w'||=%.6e" %
            (name, rnorm, rr, np.linalg.norm(x)))
        assert rr <= 1e-10, "tangent true residual above 1e-10 gate"
        wprime[name] = x
    np.savez(cache_f, **wprime)

# ---------------- 3. objective sources ----------------
# The rhs EXPORT layout is per-cell interleaved (Ux,Uy,Uz,p per cell,
# solveDiscreteFlowAdjoint.H:3113-3119) while the matrix export and the
# rxc_analytic files are block-major [all-U, all-P].  De-interleave here.
def load_rhs_export(path):
    v = np.asarray(scipy.io.mmread(path)).ravel()
    assert v.size == NUNK
    m = v.reshape(N, 4)
    b = np.zeros(NUNK)
    b[0:NV:3] = m[:, 0]
    b[1:NV:3] = m[:, 1]
    b[2:NV:3] = m[:, 2]
    b[P0:P0 + N] = m[:, 3]
    return b

bPD = load_rhs_export(EXPORT_CASE + "/explicitRhs_pressureDrop.mtx")
bTC = load_rhs_export(EXPORT_CASE + "/explicitRhs_thermalCoupling.mtx")
log("b loaded: |bPD|=%.6e |bTC|=%.6e" % (np.linalg.norm(bPD), np.linalg.norm(bTC)))

# integrity cross-check A: b^T w' vs the B2-run ADJ projections
log("[A] linearization-point cross-check (b^T w' vs BFINAL-012 ADJ):")
for name in ["D1", "D2", "D3"]:
    cPD = float(bPD @ wprime[name])
    cTC = float(bTC @ wprime[name])
    log("  %s: bPD^T w'=%.10f (ref ADJ gDP %.10f, rel %.2e) | "
        "bTC^T w'=%.10f (ref ADJ J %.10f, rel %.2e)" %
        (name, cPD, REF["gDP"]["ADJ"][name], abs(cPD - REF["gDP"]["ADJ"][name]) / abs(REF["gDP"]["ADJ"][name]),
         cTC, REF["J"]["ADJ"][name], abs(cTC - REF["J"]["ADJ"][name]) / max(abs(REF["J"]["ADJ"][name]), 1e-12)))

# ---------------- 4. state-FD truth ----------------
def load_state(dname, h, sign):
    upat = STATE_CASE + "/stageB2/wstate_%s_h%s_%s_U.mtx" % (dname, h, sign)
    ppat = STATE_CASE + "/stageB2/wstate_%s_h%s_%s_p.mtx" % (dname, h, sign)
    U = np.loadtxt(upat, dtype=np.float64)
    p = np.loadtxt(ppat, dtype=np.float64)
    return U, p

def state_to_w(U, p):
    w = np.zeros(NUNK)
    w[0:NV:3] = U[:, 0]
    w[1:NV:3] = U[:, 1]
    w[2:NV:3] = U[:, 2]
    w[P0:P0 + N] = p
    return w

# discover available h tags per direction from filenames
hmap = {"0.0003": 3e-4, "0.001": 1e-3, "0.003": 3e-3}

results = {"tangent_residuals": {}, "compare": {}, "contractions": {}}
for name in ["D1", "D2", "D3"]:
    tags = sorted(set(
        os.path.basename(f).split("_h")[1].split("_")[0]
        for f in glob.glob(STATE_CASE + "/stageB2/wstate_%s_h*_p.mtx" % name)))
    log("%s state-export h tags: %s" % (name, tags))
    for tag in tags:
        h = hmap[tag]
        Up, pp = load_state(name, tag, "p")
        Um, pm = load_state(name, tag, "m")
        w_true = (state_to_w(Up, pp) - state_to_w(Um, pm)) / (2.0 * h)
        wp = wprime[name]
        mask = np.ones(NUNK, dtype=bool)
        mask[PREF] = False           # pinned-row artifact cell
        d = wp - w_true
        cosv = float(wp @ w_true / (np.linalg.norm(wp) * np.linalg.norm(w_true)))
        rel = float(np.linalg.norm(d) / np.linalg.norm(w_true))
        relU = float(np.linalg.norm(d[0:NV]) / np.linalg.norm(w_true[0:NV]))
        relP = float(np.linalg.norm(d[P0:]) / np.linalg.norm(w_true[P0:]))
        # spatial pattern: top mismatch cells in U block
        du = np.sqrt(d[0:NV:3]**2 + d[1:NV:3]**2 + d[2:NV:3]**2)
        top = np.argsort(du)[::-1][:10]
        cPDp = float(bPD @ wp); cPDt = float(bPD @ w_true)
        cTCp = float(bTC @ wp); cTCt = float(bTC @ w_true)
        log("[B] %s h=%s: cos=%.6f relL2=%.4e relU=%.4e relP=%.4e" %
            (name, tag, cosv, rel, relU, relP))
        log("     top mismatch cells: %s" % list(map(int, top[:5])))
        log("[C] contractions gDP: b^T w'=%.6f b^T w_true=%.6f ratio=%.4f "
            "(end-to-end ref %.4f)" %
            (cPDp, cPDt, cPDp / cPDt if cPDt != 0 else float("nan"),
             REF["gDP"]["ADJ"][name] / REF["gDP"]["FD"][name]))
        log("            J  : b^T w'=%.6ef b^T w_true=%.6ef ratio=%.4f" %
            (cTCp, cTCt, cTCp / cTCt if cTCt != 0 else float("nan")))
        results["compare"]["%s_h%s" % (name, tag)] = {
            "cos": cosv, "relL2": rel, "relU": relU, "relP": relP,
            "top_cells": [int(x) for x in top]}
        results["contractions"]["%s_h%s" % (name, tag)] = {
            "gDP_wp": cPDp, "gDP_wt": cPDt, "J_wp": cTCp, "J_wt": cTCt}

with open(CYCLE + "/t2t3_results.json", "w") as f:
    json.dump(results, f, indent=1)
log("T2T3_DONE t=%.1fs -> t2t3_results.json" % (time.time() - t0))
