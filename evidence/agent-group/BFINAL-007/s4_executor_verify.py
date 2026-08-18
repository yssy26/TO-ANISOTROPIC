#!/usr/bin/env python3
"""
BFINAL-007 S4 (Executor, final-report step) INDEPENDENT re-derivation.
Fresh code path (written for S4, not reusing s3_* logic beyond matrix
conventions documented in S1) that re-derives every load-bearing number of
the FINAL_REPORT from the RAW artifacts, so the Post-Reviewer can compare.

Checks:
  Q1 : FD_lin (artifacts) max/L2/outlet-mass-fraction/top-1%  [+ independent
       rebuild of -(L_prod . n_P) from diag_bnd + upper + owner/neighbour]
  Q2 : outlet internalCoeffs_b / boundaryCoeffs_b ranges; candidate diag
  Q3 : exported J P-P diag/upper vs probe laplacian L_int; L_prod - J_PP
  Q4a: ||J_exp . n||/||n|| and ||J_cand . n||/||n|| (candidate = -L_prod)
  Q4b: candidate J_PP . n_P vs FD_lin (residual-convention closure)
  Q4c: J_cand . v == J_exp . v for the dp=0 Gate-A direction (regression)
  Q5 : design boundary term d(phi_b)/d(alpha) magnitude + design support
"""
import numpy as np, scipy.io, scipy.sparse as sp, json, time
CASE = "/home/ys/dsH/b2_case_smoke"
ART  = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-007/artifacts"
OUT  = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-007"
N = 33600; NV = 3*N; NUNK = 4*N; PREF = NV
t0 = time.time()
def log(m): print(m, flush=True)
res = {}

def col(p): return np.loadtxt(p, dtype=np.float64)
def read_of_ints(path, n):
    txt = open(path).read(); i = txt.index('(')
    return np.array([int(x) for x in txt[i+1:txt.rindex(')')].split()[:n]], dtype=np.int64)

diagL  = col(ART + "/stageB7_laplacian_diag.mtx")       # interior-only diag
diagB  = col(ART + "/stageB7_laplacian_diag_bnd.mtx")   # + addBoundaryDiag
upperL = col(ART + "/stageB7_laplacian_upper.mtx")
srcL   = col(ART + "/stageB7_laplacian_source.mtx")
nP     = col(ART + "/stageB7_nP.mtx")
fdLin  = col(ART + "/stageB7_RP_FD_lin.mtx")
fdCD   = col(ART + "/stageB7_RP_FD.mtx")
cand   = col(ART + "/stageB7_candidate_diag.mtx")
op     = np.loadtxt(ART + "/stageB7_outlet_patch.mtx")
ocells = op[:, 0].astype(int)
icB    = op[:, 3]
own = read_of_ints(CASE + "/constant/polyMesh/owner", len(upperL))
nei = read_of_ints(CASE + "/constant/polyMesh/neighbour", len(upperL))

# ============ Q1 ============
log("== Q1: actual-residual FD_lin from artifacts ==")
mx = np.abs(fdLin).max()
l2 = np.linalg.norm(fdLin)
om = np.zeros(N, bool); om[ocells] = True
outFrac = np.abs(fdLin[om]).sum() / np.abs(fdLin).sum()
top = np.argsort(np.abs(fdLin))[-336:]
topFrac = np.abs(fdLin[top]).sum() / np.abs(fdLin).sum()
res["Q1_FDlin_maxabs"] = float(mx)
res["Q1_FDlin_L2"] = float(l2)
res["Q1_outlet_mass_fraction"] = float(outFrac)
res["Q1_top1pct_mass_fraction"] = float(topFrac)
log("  max|FD|=%.4e L2=%.4e outlet-frac=%.4f top1pct-frac=%.4f"
    % (mx, l2, outFrac, topFrac))
# independent rebuild of -(L_prod . n_P) from diag_bnd + upper
Lnp = -(diagB * nP)
for f in range(len(upperL)):
    Lnp[own[f]] += -upperL[f]*nP[nei[f]]
    Lnp[nei[f]] += -upperL[f]*nP[own[f]]
rel = np.linalg.norm(fdLin - Lnp)/np.linalg.norm(Lnp)
res["Q1_FDlin_vs_minusLprod_nP"] = float(rel)
log("  |FD_lin - (-L_prod.n_P)|/|L_prod.n_P| = %.3e (independent rebuild)" % rel)

# ============ Q2 ============
log("== Q2: outlet boundary coefficients ==")
res["Q2_icB_range"] = [float(icB.min()), float(icB.max())]
res["Q2_bcB_maxabs"] = float(np.abs(op[:, 4]).max())
res["Q2_candidate_range"] = [float(cand[ocells].min()), float(cand[ocells].max())]
log("  internalCoeffs_b=[%.4e,%.4e] boundaryCoeffs_b max|.|=%.2e cand=[%.4e,%.4e]"
    % (icB.min(), icB.max(), np.abs(op[:,4]).max(), cand[ocells].min(), cand[ocells].max()))

# ============ Q3 ============
log("== Q3: exported J P-P vs production laplacian ==")
Jt = scipy.io.mmread(CASE + "/explicitJT.mtx").tocsr()
J = Jt.T.tocsr(); del Jt
J = J.tolil(); J[PREF, :] = 0.0; J[PREF, PREF] = 1.0; J = J.tocsr()
rows, cols, vals = [], [], []
for c in range(N):
    r = NV + c; s, e = J.indptr[r], J.indptr[r+1]
    for j in range(s, e):
        cc = J.indices[j]
        if NV <= cc < NUNK:
            rows.append(c); cols.append(cc - NV); vals.append(J.data[j])
expPP = sp.coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsr()
expdiag = np.asarray(expPP.diagonal()).copy()
mask = np.ones(N, bool); mask[0] = False
relD = np.linalg.norm((diagL - expdiag)[mask]) / np.linalg.norm(expdiag[mask])
res["Q3_diag_JPP_vs_Lint"] = float(relD)
dProd = diagB - expdiag
outMask = np.zeros(N, bool); outMask[ocells] = True
res["Q3_prod_minus_export_offoutlet_max"] = float(np.abs(dProd[mask & ~outMask]).max())
res["Q3_prod_minus_export_outlet_range"] = [float(dProd[outMask].min()), float(dProd[outMask].max())]
log("  J_PP diag vs L_int relL2=%.3e ; L_prod-J_PP off-outlet max=%.3e ; outlet=[%.4e,%.4e]"
    % (relD, np.abs(dProd[mask & ~outMask]).max(), dProd[outMask].min(), dProd[outMask].max()))

# ============ Q4a ============
log("== Q4a: null-vector ratio (exported vs candidate) ==")
w = np.loadtxt("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-006/cycle-2/artifacts/wPrime_TAN_D1.mtx", dtype=np.float64)
nvec = w / np.linalg.norm(w)
res["Q4a_w_norm"] = float(np.linalg.norm(w))
res["Q4a_n_PREF"] = float(nvec[PREF])
rn_exp = np.linalg.norm(J @ nvec)
# candidate: P-P := -L_prod, built from probe exports directly (not from expPP)
candPP = sp.dok_matrix((N, N), dtype=np.float64)
for c in range(N):
    candPP[c, c] = -diagB[c]
for f in range(len(upperL)):
    candPP[own[f], nei[f]] += -upperL[f]
    candPP[nei[f], own[f]] += -upperL[f]
candPP = candPP.tocsr()
Jc = J.tolil()
for c in range(N):
    r = NV + c; rr, dd = [], []
    for k, cc in enumerate(Jc.rows[r]):
        if not (NV <= cc < NUNK):
            rr.append(cc); dd.append(Jc.data[r][k])
    Jc.rows[r] = rr; Jc.data[r] = dd
s_e = candPP.indptr
for c in range(N):
    r = NV + c
    for j in range(s_e[c], s_e[c+1]):
        Jc[r, NV + candPP.indices[j]] = candPP.data[j]
Jc = Jc.tocsr(); Jc = Jc.tolil(); Jc[PREF, :] = 0.0; Jc[PREF, PREF] = 1.0; Jc = Jc.tocsr()
rn_cand = np.linalg.norm(Jc @ nvec)
res["Q4a_nullvec_Jexp"] = float(rn_exp)
res["Q4a_nullvec_Jcand"] = float(rn_cand)
log("  ||J_exp.n||/||n||=%.3e  ||J_cand.n||/||n||=%.3e  (ratio %.1e)"
    % (rn_exp, rn_cand, rn_cand/max(rn_exp,1e-300)))

# ============ Q4b ============
log("== Q4b: candidate J_PP . n_P vs actual-residual FD ==")
candPPn = candPP @ nP
relB = np.linalg.norm(candPPn - fdLin)/np.linalg.norm(fdLin)
relBcd = np.linalg.norm(candPPn - fdCD)/np.linalg.norm(fdCD)
res["Q4b_candPPn_vs_FDlin"] = float(relB)
res["Q4b_candPPn_vs_FDcd"] = float(relBcd)
log("  candPP.n_P vs FD_lin relL2=%.3e ; vs FD_cd(1e-4) relL2=%.3e" % (relB, relBcd))

# ============ Q4c ============
log("== Q4c: dp=0 Gate-A direction unchanged (regression) ==")
dU = col(CASE + "/stageB5_dUdir.mtx")
v = np.zeros(NUNK); v[:NV] = dU
Jv_exp = J @ v
Jv_cand = Jc @ v
relV = np.linalg.norm(Jv_cand - Jv_exp)/np.linalg.norm(Jv_exp)
res["Q4c_Jv_cand_vs_exp"] = float(relV)
Jv_orc = col(CASE + "/stageB5_gateA_Jv.mtx")
idx = np.arange(NUNK) != PREF
relO = np.linalg.norm((Jv_exp - Jv_orc)[idx])/np.linalg.norm(Jv_orc[idx])
res["Q4c_Jexp_vs_oracle_exclPREF"] = float(relO)
log("  J_cand.v == J_exp.v relL2=%.3e ; J_exp.v vs oracle (excl PREF) relL2=%.3e" % (relV, relO))

# ============ Q5 ============
log("== Q5: R_P,x boundary design term ==")
rAU = np.loadtxt(ART + "/stageB7_rebuilt_rAU.mtx")
def read_of(path, n):
    t = open(path).read(); i = t.index('(')
    return np.array([float(x) for x in t[i+1:t.rindex(')')].split()[:n]])
p = read_of(CASE + "/1/p", N)
alphaRel = 0.4
drAU = -rAU**2/alphaRel
term = drAU[ocells]*op[:, 1]*op[:, 2]*p[ocells]
res["Q5_term_range"] = [float(term.min()), float(term.max())]
res["Q5_term_L2"] = float(np.linalg.norm(term))
da = np.loadtxt(CASE + "/stageB6_deltaAlpha.mtx")
res["Q5_deltaAlpha_outlet_nnz"] = int((np.abs(da[ocells]) > 1e-300).sum())
dm = read_of(CASE + "/1/designMask", N)
res["Q5_designMask_outlet_max"] = float(dm[ocells].max())
log("  term range=[%.3e,%.3e] L2=%.3e ; deltaAlpha outlet nnz=%d ; designMask outlet max=%.1f"
    % (term.min(), term.max(), np.linalg.norm(term),
       res["Q5_deltaAlpha_outlet_nnz"], res["Q5_designMask_outlet_max"]))

res["_runtime_s"] = time.time() - t0
with open(OUT + "/s4_executor_verify.json", "w") as f:
    json.dump(res, f, indent=1, default=float)
log("\nS4_EXECUTOR_VERIFY_DONE (t=%.1fs)" % res["_runtime_s"])
