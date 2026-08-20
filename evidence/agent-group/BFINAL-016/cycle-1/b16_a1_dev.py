#!/usr/bin/env python3
"""
BFINAL-016 A1 - deviatoric-stress alignment gauge.

Replicates applyProdFrozenDeviatoricJT (solveDiscreteFlowAdjointProduction.H:105-231)
in numpy from the one-shot mesh export (b16mesh_*), adds it to the
exported-matrix tangent via the fixed-point iteration

    w_{k+1} = splu(J_core) * (rhs - J_dev * w_k)

and tests whether the resulting b_PD^T w' moves onto the PRODUCTION ADJ
values (-5.4050587339 / +1.10592698054 / +14.0447142088) to >= 3
significant digits.  Also re-measures the P-row concentration of the
mismatch on the COMPLETE operator.
"""
import os, time, json
import numpy as np
import scipy.io, scipy.sparse as sp, scipy.sparse.linalg as spla

MESH = "/home/ys/dsH/b16_mesh"
EXPORT_JT = "/home/ys/dsH/b8_verify_diag/explicitJT.mtx"
EXPORT_CASE = "/home/ys/dsH/b15_export"
CYCLE = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-016/cycle-1"

N = 33600
NV = 3 * N
NUNK = 4 * N
PREF = NV
P0 = NV

REF_ADJ_PD = {"D1": -5.4050587339, "D2": 1.10592698054, "D3": 14.0447142088}
REF_FD_PD = {"D1": -2.537460, "D2": 0.523409, "D3": 6.477216}

t0 = time.time()
def log(m): print(m, flush=True)

log("=== BFINAL-016 A1: deviatoric alignment ===")

# ---------------- mesh data ----------------
owner = np.loadtxt(MESH + "/b16mesh_owner.mtx", dtype=np.int64)
nei = np.loadtxt(MESH + "/b16mesh_neighbour.mtx", dtype=np.int64)
Sf = np.loadtxt(MESH + "/b16mesh_sf.mtx")
w = np.loadtxt(MESH + "/b16mesh_weights.mtx")
Vn = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")
V = Vn[:, 0]; nue = Vn[:, 1]
bnd = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
log("mesh loaded: nIF=%d cells=%d nBnd=%d (t=%.1fs)" %
    (len(owner), N, len(bnd), time.time() - t0))
bc = bnd[:, 0].astype(np.int64)
bS = bnd[:, 1:4]; bmag = bnd[:, 4]; bdel = bnd[:, 5]
bnue = bnd[:, 6]; bfix = bnd[:, 7] > 0.5
bn = bS / bmag[:, None]

def dev2T(A):
    # dev2(T(A)) = A^T - (2/3) tr(A) I
    out = np.swapaxes(A, 1, 2).copy()
    tr = np.trace(A, axis1=1, axis2=2)
    out[:, 0, 0] -= 2.0/3.0*tr
    out[:, 1, 1] -= 2.0/3.0*tr
    out[:, 2, 2] -= 2.0/3.0*tr
    return out

Vown = V[owner]; Vnei = V[nei]
wcol = w[:, None, None]; owcol = (1.0 - w)[:, None, None]

def apply_Jdev(x):
    U = x[:NV].reshape(N, 3)
    dl = U[nei] - U[owner]
    outer_f = Sf[:, :, None]*dl[:, None, :]          # (nF,3,3)
    st = np.zeros((N, 3, 3))
    np.add.at(st, owner, wcol*outer_f)
    np.add.at(st, nei, owcol*outer_f)
    # boundary stress + tangential-projected gradient contribution
    bg = bnue[:, None, None]*dev2T(-bS[:, :, None]*U[bc][:, None, :])
    projterm = bg - bn[:, :, None]*(np.einsum('fi,fij->fj', bn, bg))[:, None, :]
    gr = nue[:, None, None]*dev2T(st)
    # note: boundary unprojected bg feeds BOTH the projection (gr) and the
    # fixed-U normal output term below (per production.H:158-172)
    gr_pre = np.zeros((N, 3, 3))
    np.add.at(gr_pre, bc, projterm)
    gr += gr_pre
    # interior face output
    gO = gr[owner]/np.maximum(Vown, 1e-300)[:, None, None]
    gN = gr[nei]/np.maximum(Vnei, 1e-300)[:, None, None]
    faceVel = np.einsum('fi,fij->fj', Sf, gO - gN)   # Sf & (gO - gN)
    out = np.zeros(NUNK)
    dU = out[:NV].reshape(N, 3)
    np.add.at(dU, owner, w[:, None]*faceVel)
    np.add.at(dU, nei, (1.0 - w)[:, None]*faceVel)
    # non-fixed-U boundary: Sf & (gr_c/V_c)
    nf = ~bfix
    bvel = np.einsum('fi,fij->fj', bS[nf], gr[bc[nf]]/np.maximum(V[bc[nf]], 1e-300)[:, None, None])
    np.add.at(dU, bc[nf], bvel)
    # fixed-U boundary: -delta*(normal & bg) componentwise
    fx = bfix
    ngrad = -bdel[fx, None]*np.einsum('fi,fij->fj', bn[fx], bg[fx])
    np.add.at(dU, bc[fx], ngrad)
    return out

Jdev = spla.LinearOperator((NUNK, NUNK), matvec=apply_Jdev)

# ---------------- core matrix + splu ----------------
M = scipy.io.mmread(EXPORT_JT).tocsr()
J = M.T.tocsr()
del M
J = J.tolil()
J[PREF, :] = 0.0
J[PREF, PREF] = 1.0
J = J.tocsr()
lu = spla.splu(J, permc_spec="COLAMD")
log("splu done (t=%.1fs)" % (time.time() - t0))

# ---------------- sources ----------------
rxc = np.loadtxt(EXPORT_CASE + "/stageB6_rxc_analytic.mtx", dtype=np.float64)

def load_rhs_export(path):
    v = np.asarray(scipy.io.mmread(path)).ravel()
    m = v.reshape(N, 4)
    b = np.zeros(NUNK)
    b[0:NV:3] = m[:, 0]; b[1:NV:3] = m[:, 1]; b[2:NV:3] = m[:, 2]
    b[P0:P0 + N] = m[:, 3]
    return b

bPD = load_rhs_export(EXPORT_CASE + "/explicitRhs_pressureDrop.mtx")

# ---------------- tangent solves with J_dev ----------------
results = {}
for di, name in enumerate(["D1", "D2", "D3"]):
    off = di * NUNK
    rhs = -rxc[off:off + NUNK].copy()
    rhs[PREF] = 0.0
    x = lu.solve(rhs)
    for it in range(200):
        xn = lu.solve(rhs - Jdev @ x)
        dw = np.linalg.norm(xn - x)/max(np.linalg.norm(xn), 1e-300)
        x = xn
        if dw < 1e-13:
            break
    resid = np.linalg.norm((J @ x + Jdev @ x) - rhs)/np.linalg.norm(rhs)
    c = float(bPD @ x)
    ref = REF_ADJ_PD[name]
    # significant-digit agreement
    def sigdig(a, b):
        if a == b: return 16
        d = abs(a - b)/max(abs(b), 1e-300)
        return max(0, int(-np.floor(np.log10(d))) ) if d < 1 else 0
    log("A1 %s: iters=%d dw=%.2e combinedResid=%.2e bPD^T w'=%.10f "
        "(prod ADJ %.10f, sigdig=%d, rel=%.3e)" %
        (name, it, dw, resid, c, ref, sigdig(c, ref), abs(c-ref)/abs(ref)))
    results[name] = {"contraction": c, "prod_ADJ": ref,
                     "rel": abs(c-ref)/abs(ref),
                     "FD": REF_FD_PD[name],
                     "ratio_to_FD": c/REF_FD_PD[name]}

with open(CYCLE + "/a1_dev_alignment.json", "w") as f:
    json.dump(results, f, indent=1)
log("A1_DONE t=%.1fs" % (time.time() - t0))
