#!/usr/bin/env python3
"""BFINAL-030 fingerprint 1 (operator side): patch explicitJT_H7 T3 entries
and re-solve w_corr.

Patches the exported J^T matrix M (b24_diag/explicitJT_H7.mtx) at the
OFF-DIAGONAL U-row/P-column (1-alphaRel) direct-relax entries ONLY:
    M_corrected = M + (c* - 0.6) * T3_offdiag
T3_offdiag (coefficient 1, INTERNAL faces only, off-diagonal slots):
    [U(own,c),P(nei)] -= w*sf_c        (tag 'on')
    [U(nei,c),P(own)] += (1-w)*sf_c    (tag 'no')
DIAGONAL U-row/P-col positions [U(cell,c),P(cell)] are NOT pure T3 slots:
the hA/hB Stage-2 path (solveDiscreteFlowAdjoint.H 2719-2807) writes to the
same positions, so production values there are composite.  Boundary
U-row/P-col is likewise composite (hB alphaRel*rAU path + direct
(1-alphaRel)*patchSf + internal accumulation).  Both are skipped.

Position-level purity mask: within the off-diagonal slots, production M must
equal BASE_COEF*T3_offdiag = 0.6*T3_offdiag exactly (|rel| < 1e-6).  The
remaining ~5% of off-diagonal positions (boundary-adjacent cells where the
hA[own]/hA[nei] closed-face sum over internal faces does not cancel to zero,
so hA contamination is present) are left at their production composite
values and are documented as excluded.  Verified empirically: 'on'/'no'
median rel = 1e-14 at pure positions; max rel = 0.398 at contaminated ones.

Solves M_corrected^T w = -rxd for D1/D2/D3 with the b26-proven default
COLAMD splu and reports closure |r_P|/|rxdP| per direction.

Usage: b30_fp1_patch_resolve.py [cstar]   (default 1.0)
Outputs: b30_fp1_wcorr_c<cstar>.npz  +  b30_fp1_resolve_c<cstar>.json
"""
import os, sys, time, json
import numpy as np
import scipy.io
import scipy.sparse as sp
import scipy.sparse.linalg as spla

CSTAR = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0

EXPORT_JT = "/home/ys/dsH/b24_diag/explicitJT_H7.mtx"
RXC = "/home/ys/dsH/b8_verify_diag/stageB6_rxc_analytic.mtx"
STATE_CASE = "/home/ys/dsH/b25_qgate"
MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-030/cycle-1"

N = 33600; NV = 3 * N; NUNK = 4 * N; P0 = NV
ALPHA_REL = 0.4
BASE_COEF = 1.0 - ALPHA_REL   # 0.6

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-030 fp1 patch-resolve (off-diag only) c*=%.3f ===" % CSTAR)

# ---------------- mesh ----------------
pts = []
intxt = False
for line in open(POLY + "/points"):
    s = line.strip()
    if not intxt:
        if s == "(": intxt = True
        continue
    if s == ")": break
    pts.append([float(t) for t in s.strip("()").split()])
pts = np.array(pts)
faces = []
intxt = False
for line in open(POLY + "/faces"):
    s = line.strip()
    if not intxt:
        if s == "(": intxt = True; continue
    if s == ")": break
    if "(" in s and s.endswith(")"):
        faces.append([int(v) for v in s[s.index("(") + 1:-1].split()])
def rl(p):
    t = open(p).read(); i = t.index("\n("); j = t.index("\n)", i)
    return np.array([int(x) for x in t[i + 2:j].split()], dtype=np.int64)
oa = rl(POLY + "/owner"); na = rl(POLY + "/neighbour")
nIF = len(na); nF = len(faces)
Sf = np.zeros((nF, 3))
for fi in range(nF):
    fv = pts[faces[fi]]
    Sf[fi] = 0.5 * np.cross(fv, np.roll(fv, -1, axis=0)).sum(axis=0)
owner = oa[:nIF]; nei = na[:nIF]
w = np.loadtxt(MESH + "/b16mesh_weights.mtx")
bnd = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
bc_cell = bnd[:, 0].astype(np.int64); bSf = bnd[:, 1:4]
Ufix = bnd[:, 7] > 0.5
bidx = np.where(~Ufix)[0]   # assignable boundary faces (outlet, zeroGradient)
log("mesh: nIF=%d nB=%d assignable=%d (t=%.1fs)"
    % (nIF, len(bnd), len(bidx), time.time() - t0))

def U(c, cm): return 3 * c + cm
def P(c): return NV + c

# ---------------- T3_offdiag (internal faces, off-diagonal slots) ----------------
ri = []; ci = []; vi = []; ti = []
SfI = Sf[:nIF]
for f in range(nIF):
    o = owner[f]; n = nei[f]
    wf = w[f]
    sfc = SfI[f]
    for cm in range(3):
        if sfc[cm] == 0:
            continue
        # off-diagonal U-row/P-col slots only (pure (1-alphaRel) T3)
        ri += [U(o, cm), U(n, cm)]
        ci += [P(n), P(o)]
        vi += [-wf * sfc[cm], (1.0 - wf) * sfc[cm]]
        ti += ['on', 'no']
T3b = sp.coo_matrix((np.array(vi), (np.array(ri), np.array(ci))),
                    shape=(NUNK, NUNK))
T3b.sum_duplicates()
log("T3_offdiag nnz=%d (internal off-diagonal only) (t=%.1fs)"
    % (T3b.nnz, time.time() - t0))

# ---------------- load M + purity gate (position-level mask) ----------------
M = scipy.io.mmread(EXPORT_JT)
log("M raw coo nnz=%d (t=%.1fs)" % (M.nnz, time.time() - t0))
M_csc = M.tocsc(); M_csc.sort_indices()
ii, jj, vv = T3b.row, T3b.col, T3b.data
orig = np.asarray(M_csc[ii, jj]).reshape(-1)
# NOTE: uniform mesh (w=0.5) gives exact cancellations at diagonal slots only;
# off-diagonal slots are unique per face so vv != 0 whenever sf_c != 0.
vmax = float(np.abs(vv).max())
sig = np.abs(vv) > 1e-9 * vmax
rel = np.abs(orig[sig] - BASE_COEF * vv[sig]) / np.maximum(np.abs(orig[sig]), np.abs(vv[sig]))
gate_max = float(rel.max())
gate_med = float(np.median(rel))
frac_exact = float((rel < 1e-6).mean())
pure = np.zeros(len(vv), dtype=bool)
pure[sig] = rel < 1e-6
n_pure = int(pure.sum())
log("OFFDIAG purity gate: sig=%d/%d medianrel=%.3e maxrel=%.3e frac_exact=%.6f pure=%d"
    % (int(sig.sum()), len(vv), gate_med, gate_max, frac_exact, n_pure))
if not (gate_med < 1e-9 and frac_exact > 0.9):
    log("FATAL: offdiag purity gate FAILED; aborting before solve")
    sys.exit(2)
excluded = int(sig.sum()) - n_pure
log("positions excluded (hA/hB-contaminated off-diagonal): %d (%.4f%%)"
    % (excluded, 100.0 * excluded / max(1, int(sig.sum()))))

# ---------------- build M_corrected (patch pure positions only) ----------------
scale = CSTAR - BASE_COEF
Mt = sp.coo_matrix((scale * vv[pure], (ii[pure], jj[pure])),
                   shape=(NUNK, NUNK))
Mt.sum_duplicates()
Mc = (M_csc + Mt).tocsc()
Mc.sort_indices()
del M, Mt
log("M_corrected csc nnz=%d (t=%.1fs)" % (Mc.nnz, time.time() - t0))

# ---------------- splu ----------------
# b26-proven path: default COLAMD (~17 min).  MMD_AT_PLUS_A is pathological
# on this near-singular matrix (never finishes); skipped deliberately.
lu = spla.splu(Mc)
log("splu default COLAMD done (t=%.1fs)" % (time.time() - t0))

rxc = np.loadtxt(RXC, dtype=np.float64)
log("rxc loaded %s (t=%.1fs)" % (rxc.shape, time.time() - t0))

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
    rxd = rxc[di * NUNK:(di + 1) * NUNK]
    w = lu.solve(-rxd, trans="T")     # M_corr^T w = -rxd  <=>  J_corr w = -rxd
    wstar[nm] = w
    Jw = (Mc.T @ w)
    r_P = Jw[P0:] + rxd[P0:]
    entry = {"|r_P(w*)|/|rxdP|": float(np.linalg.norm(r_P) / np.linalg.norm(rxd[P0:]))}
    try:
        wt = (load_state(nm, "0.001", "p") - load_state(nm, "0.001", "m")) / 2e-3
        d = wt - w
        entry.update({
            "cos_w_wt": float(w @ wt / (np.linalg.norm(w) * np.linalg.norm(wt))),
            "relL2": float(np.linalg.norm(d) / np.linalg.norm(wt)),
            "relU": float(np.linalg.norm(d[:NV]) / np.linalg.norm(wt[:NV])),
            "relP": float(np.linalg.norm(d[P0:]) / np.linalg.norm(wt[P0:])),
        })
    except IOError as e:
        log("  %s: wstate files missing (%s) -- closure only" % (nm, e))
    res[nm] = entry
    log("  %s: |r_P(w*)|/|rxdP|=%.3e  cos=%.6f relL2=%.4f (relU=%.4f relP=%.4f)"
        % (nm, entry["|r_P(w*)|/|rxdP|"], entry.get("cos_w_wt", float("nan")),
           entry.get("relL2", float("nan")), entry.get("relU", float("nan")),
           entry.get("relP", float("nan"))))
    log("  (t=%.1fs)" % (time.time() - t0))

meta = {
    "cstar": CSTAR, "base_coef": BASE_COEF, "scale": scale,
    "nnz_M": Mc.nnz, "nnz_T3_offdiag": T3b.nnz,
    "gate_offdiag_medianrel": gate_med, "gate_offdiag_maxrel": gate_max,
    "gate_offdiag_frac_exact": frac_exact,
    "n_pure_patched": n_pure, "n_excluded_composite": excluded,
    "patch": "internal off-diagonal U-row/P-col slots, |rel|<1e-6 position mask",
    "solver": "scipy splu default COLAMD",
}
res["meta"] = meta
tag = "c%.2f" % CSTAR
np.savez(OUT + "/b30_fp1_wcorr_%s.npz" % tag, **wstar)
with open(OUT + "/b30_fp1_resolve_%s.json" % tag, "w") as f:
    json.dump(res, f, indent=1)
log("wrote %s/b30_fp1_wcorr_%s.npz and b30_fp1_resolve_%s.json (t=%.1fs)"
    % (OUT, tag, tag, time.time() - t0))
