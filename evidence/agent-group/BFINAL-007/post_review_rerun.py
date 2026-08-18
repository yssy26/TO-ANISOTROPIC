#!/usr/bin/env python3
"""
BFINAL-007 POST-REVIEWER independent re-verification.
Written fresh (independent code path, not derived from the executor's s3_*/s4_*
scripts) to re-derive the load-bearing numbers of Q1..Q5 from the RAW artifacts.

Load-bearing checks re-run here (per the STOP rule):
  * Q1  actual frozen-primal pressure/SIMPLE residual central FD along n_P
        (non-zero? outlet-concentrated?) + independent rebuild of -(L_prod.n_P)
  * Q2  outlet fixedValue internalCoeffs_b / boundaryCoeffs_b
  * Q3  exported J P-P block vs production laplacian (what is missing)
  * Q4a candidate J null-mode (||J.n||/||n|| exported vs candidate)
  * Q4b candidate J_PP.n_P vs actual-residual FD (Jv-FD closure)
  * Q4c BFINAL-003 dp=0 Gate-A direction unchanged
  * Q5  R_P,x boundary design term magnitude + design support on the outlet
"""
import json
import time
import numpy as np
import scipy.io
import scipy.sparse as sp

ART = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-007/artifacts"
CASE = "/home/ys/dsH/b2_case_smoke"
WVEC = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-006/cycle-2/artifacts/wPrime_TAN_D1.mtx"

N = 33600
NV = 3 * N
NUNK = 4 * N
PREF = NV  # discretePIndex(0) = 100800

t0 = time.time()
R = {}


def col(path):
    return np.loadtxt(path, dtype=np.float64)


def of_list(path, n, dtype):
    txt = open(path).read()
    i = txt.index("(")
    toks = txt[i + 1:txt.rindex(")")].split()
    return np.array([dtype(x) for x in toks[:n]])


print("== loading raw artifacts ==", flush=True)
diagL = col(ART + "/stageB7_laplacian_diag.mtx")       # interior-only laplacian diag
diagB = col(ART + "/stageB7_laplacian_diag_bnd.mtx")   # + addBoundaryDiag (production)
upper = col(ART + "/stageB7_laplacian_upper.mtx")
nP = col(ART + "/stageB7_nP.mtx")                      # P block of normalized null vec
fdLin = col(ART + "/stageB7_RP_FD_lin.mtx")            # -div(flux(psi=n_P)) from probe
fdCD = col(ART + "/stageB7_RP_FD.mtx")                 # central FD @ eps=1e-4
cand = col(ART + "/stageB7_candidate_diag.mtx")        # +rAtU_c*delta_b*|Sf|_b per cell
op = np.loadtxt(ART + "/stageB7_outlet_patch.mtx")     # cell delta magSf ic bc mobB mobC cand
ocells = op[:, 0].astype(int)
icB = op[:, 3]
bcB = op[:, 4]
deltaB = op[:, 1]
magSfB = op[:, 2]

nIntF = upper.shape[0]
own = of_list(CASE + "/constant/polyMesh/owner", nIntF, np.int64)
nei = of_list(CASE + "/constant/polyMesh/neighbour", nIntF, np.int64)
assert own.shape[0] == nIntF and nei.shape[0] == nIntF

# ---------- helper: apply a (diag, upper) laplacian-like matrix to a field ----------
def matvec(diag_, upper_, x):
    y = diag_ * x
    for f in range(upper_.shape[0]):
        y[own[f]] += upper_[f] * x[nei[f]]
        y[nei[f]] += upper_[f] * x[own[f]]
    return y

# ---------- Q1: actual-residual FD along n_P ----------
print("\n== Q1: actual-residual central FD along n_P ==", flush=True)
fdMax = float(np.abs(fdLin).max())
fdL2 = float(np.linalg.norm(fdLin))
outMask = np.zeros(N, bool); outMask[ocells] = True
outFrac = float(np.abs(fdLin[outMask]).sum() / np.abs(fdLin).sum())
topIdx = np.argsort(np.abs(fdLin))[-336:]
topFrac = float(np.abs(fdLin[topIdx]).sum() / np.abs(fdLin).sum())
# independent rebuild of the production residual derivative -(L_prod . n_P):
# L_prod = diag (with addBoundaryDiag) + upper (symmetric), so
# dR_P/dp . n_P = -L_prod . n_P
Lprod_nP = matvec(diagB, upper, nP)
minusLprod_nP = -Lprod_nP
relQ1 = float(np.linalg.norm(fdLin - minusLprod_nP) / np.linalg.norm(minusLprod_nP))
# central-diff vs exact-linearity FD
relCD = float(np.linalg.norm(fdCD - fdLin) / np.linalg.norm(fdLin))
R["Q1_FDlin_maxabs"] = fdMax
R["Q1_FDlin_L2"] = fdL2
R["Q1_outlet_mass_fraction"] = outFrac
R["Q1_top1pct_mass_fraction"] = topFrac
R["Q1_FDlin_vs_minusLprod_nP_relL2"] = relQ1
R["Q1_FDcd_vs_FDlin_relL2"] = relCD
print("  max|FD_lin|=%.4e  L2=%.4e" % (fdMax, fdL2), flush=True)
print("  outlet mass frac=%.4f  top-1%% mass frac=%.4f" % (outFrac, topFrac), flush=True)
print("  |FD_lin - (-L_prod.n_P)|/|L_prod.n_P| = %.3e  (indep rebuild)" % relQ1, flush=True)
print("  |FD_cd(1e-4) - FD_lin|/|FD_lin| = %.3e" % relCD, flush=True)

# ---------- Q2: boundary coefficients ----------
print("\n== Q2: outlet fixedValue boundary coefficients ==", flush=True)
R["Q2_icB_range"] = [float(icB.min()), float(icB.max())]
R["Q2_bcB_maxabs"] = float(np.abs(bcB).max())
R["Q2_deltaB_range"] = [float(deltaB.min()), float(deltaB.max())]
R["Q2_magSfB"] = float(magSfB[0])
# consistency: internalCoeffs_b should == -(mobB * deltaB * magSfB)
mobB = op[:, 5]
recon_icB = -mobB * deltaB * magSfB
relRecon = float(np.linalg.norm(recon_icB - icB) / np.linalg.norm(icB))
R["Q2_icB_reconstruction_relL2"] = relRecon
print("  internalCoeffs_b=[%.4e, %.4e]  boundaryCoeffs_b max|.|=%.2e" %
      (icB.min(), icB.max(), np.abs(bcB).max()), flush=True)
print("  delta_b=[%.4e,%.4e] magSf_b=%.3e" % (deltaB.min(), deltaB.max(), magSfB[0]), flush=True)
print("  ic_b vs -(mob_b*delta_b*magSf_b) relL2=%.3e" % relRecon, flush=True)

# ---------- load exported J and re-pin ----------
print("\n== loading explicitJT.mtx (exported J^T) ==", flush=True)
Jt = scipy.io.mmread(CASE + "/explicitJT.mtx").tocsr()
print("  J^T shape=%s nnz=%d" % (Jt.shape, Jt.nnz), flush=True)
J = Jt.T.tocsr().tolil()
J[PREF, :] = 0.0
J[PREF, PREF] = 1.0
J = J.tocsr()
del Jt

# extract P-P block (rows NV..NUNK, cols NV..NUNK) of exported J
rows, cols, vals = [], [], []
for c in range(N):
    r = NV + c
    for j in range(J.indptr[r], J.indptr[r + 1]):
        cc = J.indices[j]
        if NV <= cc < NUNK:
            rows.append(c); cols.append(cc - NV); vals.append(J.data[j])
expPP = sp.coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsr()
expdiag = np.asarray(expPP.diagonal()).copy()

# ---------- Q3: what is missing ----------
print("\n== Q3: exported J P-P vs production laplacian ==", flush=True)
mask = np.ones(N, bool); mask[0] = False  # exclude pRef pin row
relDiag = float(np.linalg.norm((diagL - expdiag)[mask]) / np.linalg.norm(expdiag[mask]))
dProd = diagB - expdiag
offOut = float(np.abs(dProd[mask & ~outMask]).max())
outRange = [float(dProd[outMask].min()), float(dProd[outMask].max())]
R["Q3_diag_JPP_vs_Lint_relL2"] = relDiag
R["Q3_prod_minus_export_offoutlet_maxabs"] = offOut
R["Q3_prod_minus_export_outlet_range"] = outRange
print("  J_PP diag vs L_int diag relL2 = %.3e" % relDiag, flush=True)
print("  L_prod - J_PP diag: off-outlet max|.|=%.3e ; outlet=[%.4e, %.4e]" %
      (offOut, outRange[0], outRange[1]), flush=True)
# confirm outlet diff == internalCoeffs_b (laplacian convention)
dOut = dProd[outMask]
relIcmatch = float(np.linalg.norm(dOut - icB) / np.linalg.norm(icB))
R["Q3_outlet_diff_vs_icB_relL2"] = relIcmatch
print("  (L_prod - J_PP) outlet == internalCoeffs_b  relL2=%.3e" % relIcmatch, flush=True)

# ---------- build candidate J: P-P := -L_prod ----------
# candPP = -(diagB on diag) + (-upper on off-diag), symmetric
candPP = sp.dok_matrix((N, N), dtype=np.float64)
for c in range(N):
    candPP[c, c] = -diagB[c]
for f in range(nIntF):
    candPP[own[f], nei[f]] += -upper[f]
    candPP[nei[f], own[f]] += -upper[f]
candPP = candPP.tocsr()

Jc = J.tolil()
for c in range(N):
    r = NV + c
    keep = [(cc, v) for cc, v in zip(Jc.rows[r], Jc.data[r]) if not (NV <= cc < NUNK)]
    Jc.rows[r] = [cc for cc, _ in keep]
    Jc.data[r] = [v for _, v in keep]
for c in range(N):
    r = NV + c
    for j in range(candPP.indptr[c], candPP.indptr[c + 1]):
        Jc[r, NV + candPP.indices[j]] = candPP.data[j]
Jc = Jc.tocsr().tolil()
Jc[PREF, :] = 0.0; Jc[PREF, PREF] = 1.0
Jc = Jc.tocsr()

# ---------- Q4a: null mode ----------
print("\n== Q4a: null mode (exported vs candidate) ==", flush=True)
w = np.loadtxt(WVEC, dtype=np.float64)
wnorm = float(np.linalg.norm(w))
nvec = w / wnorm
nPREF = float(nvec[PREF])
rnExp = float(np.linalg.norm(J @ nvec))
rnCand = float(np.linalg.norm(Jc @ nvec))
R["Q4a_w_norm"] = wnorm
R["Q4a_n_PREF"] = nPREF
R["Q4a_nullvec_Jexp"] = rnExp
R["Q4a_nullvec_Jcand"] = rnCand
R["Q4a_ratio"] = rnCand / max(rnExp, 1e-300)
print("  ||w||=%.6e  n[PREF]=%.3e" % (wnorm, nPREF), flush=True)
print("  ||J_exp.n||/||n||=%.3e  ||J_cand.n||/||n||=%.3e  ratio=%.1e" %
      (rnExp, rnCand, rnCand / max(rnExp, 1e-300)), flush=True)

# ---------- Q4b: candidate J_PP . n_P vs FD ----------
print("\n== Q4b: candidate J_PP.n_P vs actual-residual FD ==", flush=True)
candPPn = candPP @ nP
relBlin = float(np.linalg.norm(candPPn - fdLin) / np.linalg.norm(fdLin))
relBcd = float(np.linalg.norm(candPPn - fdCD) / np.linalg.norm(fdCD))
R["Q4b_candPPn_vs_FDlin_relL2"] = relBlin
R["Q4b_candPPn_vs_FDcd_relL2"] = relBcd
print("  candPP.n_P vs FD_lin relL2=%.3e ; vs FD_cd(1e-4) relL2=%.3e" % (relBlin, relBcd), flush=True)
# also the "literal exported + candidate" (without interior sign flip) to show it FAILS
literalPPn = expPP @ nP + cand * nP
relLit = float(np.linalg.norm(literalPPn - fdLin) / np.linalg.norm(fdLin))
relExp = float(np.linalg.norm(expPP @ nP - fdLin) / np.linalg.norm(fdLin))
R["Q4b_expPPn_vs_FDlin_relL2"] = relExp
R["Q4b_expPPn_plus_cand_vs_FDlin_relL2"] = relLit
print("  (exported J_PP).n_P vs FD relL2=%.4f ; (exported+diag(cand)).n_P vs FD relL2=%.4f" %
      (relExp, relLit), flush=True)

# ---------- Q4c: dp=0 Gate-A direction unchanged ----------
print("\n== Q4c: dp=0 Gate-A direction (regression) ==", flush=True)
dU = col(CASE + "/stageB5_dUdir.mtx")
v = np.zeros(NUNK); v[:NV] = dU
JvExp = J @ v
JvCand = Jc @ v
relV = float(np.linalg.norm(JvCand - JvExp) / np.linalg.norm(JvExp))
R["Q4c_Jv_cand_vs_exp_relL2"] = relV
oracle = col(CASE + "/stageB5_gateA_Jv.mtx")
idx = np.arange(NUNK) != PREF
relO = float(np.linalg.norm((JvExp - oracle)[idx]) / np.linalg.norm(oracle[idx]))
R["Q4c_Jexp_vs_oracle_exclPREF_relL2"] = relO
print("  J_cand.v == J_exp.v relL2=%.3e ; J_exp.v vs oracle (excl PREF) relL2=%.3e" %
      (relV, relO), flush=True)

# ---------- Q5: R_P,x boundary design term ----------
print("\n== Q5: R_P,x boundary design term ==", flush=True)
rAU = col(ART + "/stageB7_rebuilt_rAU.mtx")
p = of_list(CASE + "/1/p", N, float)
alphaRel = 0.4
drAU = -rAU ** 2 / alphaRel
term = drAU[ocells] * deltaB * magSfB * p[ocells]
R["Q5_drAU_outlet_range"] = [float(drAU[ocells].min()), float(drAU[ocells].max())]
R["Q5_pc_outlet_range"] = [float(p[ocells].min()), float(p[ocells].max())]
R["Q5_term_range"] = [float(term.min()), float(term.max())]
R["Q5_term_L2"] = float(np.linalg.norm(term))
da = col(CASE + "/stageB6_deltaAlpha.mtx")
dm = of_list(CASE + "/1/designMask", N, float)
R["Q5_deltaAlpha_outlet_nnz"] = int((np.abs(da[ocells]) > 1e-300).sum())
R["Q5_designMask_outlet_max"] = float(dm[ocells].max())
print("  drAU outlet=[%.3e,%.3e]  p_c outlet=[%.4f,%.4f]" %
      (drAU[ocells].min(), drAU[ocells].max(), p[ocells].min(), p[ocells].max()), flush=True)
print("  d(phi_b)/d(alpha) per outlet cell: range=[%.4e, %.4e]  L2=%.4e" %
      (term.min(), term.max(), np.linalg.norm(term)), flush=True)
print("  deltaAlpha outlet nnz=%d ; designMask outlet max=%.3e" %
      (R["Q5_deltaAlpha_outlet_nnz"], R["Q5_designMask_outlet_max"]), flush=True)

R["_runtime_s"] = time.time() - t0
outp = ART + "/post_review_rerun.json"
with open(outp, "w") as f:
    json.dump(R, f, indent=1, default=float)
print("\nPOST_REVIEW_RERUN_DONE (t=%.1fs) -> %s" % (R["_runtime_s"], outp), flush=True)
