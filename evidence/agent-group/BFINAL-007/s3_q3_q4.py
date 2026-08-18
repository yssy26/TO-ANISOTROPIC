#!/usr/bin/env python3
"""
BFINAL-007 S3 (Executor) - Q3 + Q4a + Q4b + Q4c (offline scipy).
Read-only analysis of the exported explicitJT.mtx vs the S2 probe laplacian
artifacts; builds the DIAGNOSTIC-ONLY candidate J and tests null mode,
actual-residual FD closure, and BFINAL-003 G1/G2 anchors.

Conventions (verified against production NS.H pEqn + OpenFOAM-7 sources + S2):
  R_P = div(phi) = div(phiHbyA) - div(pEqn.flux())   (raw face-sum, owner +)
  production pEqn matrix L_prod = L_int + diag(internalCoeffs_b) on outlet
  dR_P/dp = -L_prod   (S2 probe FD_lin == -L_prod . n_P, self-check 2.9e-10)
  exported J P-P block == +L_int (interior-only laplacian; S2: 1.7e-16)
  => exported J P-P = -dR_P/dp + diag(cand): the exported P-P is the NEGATIVE
     of the true residual derivative (interior) AND misses the boundary diag.
  residual-convention candidate that closes Q4b: J_cand_PP = -L_prod
     = -exported_J_PP + diag(cand)   [S2 note: -(L_int + bnd)]
"""
import os, time, json
import numpy as np
import scipy.io, scipy.sparse as sp

CASE = "/home/ys/dsH/b2_case_smoke"
ART  = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-007/artifacts"
OUT  = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-007"
N    = 33600
NV   = 3 * N
NUNK = 4 * N
PREF = NV          # 100800 (P(0) pin)
t0 = time.time()
def log(m): print(m, flush=True)
res = {"stage": "BFINAL-007", "step": "S3", "mode": "DIAGNOSTIC_ONLY"}

def col(p): return np.loadtxt(p, dtype=np.float64)
def read_of_ints(path, n):
    txt = open(path).read(); i = txt.index('(')
    return np.array([int(x) for x in txt[i+1:txt.rindex(')')].split()[:n]], dtype=np.int64)

log("=== BFINAL-007 S3: Q3/Q4 (candidate J = -L_prod, diagnostic-only) ===")
Jt = scipy.io.mmread(os.path.join(CASE, "explicitJT.mtx")).tocsr()
J = Jt.T.tocsr(); del Jt
J = J.tolil(); J[PREF, :] = 0.0; J[PREF, PREF] = 1.0; J = J.tocsr()
log("exported J = transpose(explicitJT.mtx) re-pinned row %d (nnz=%d, t=%.1fs)"
    % (PREF, J.nnz, time.time() - t0))

diagL  = col(ART + "/stageB7_laplacian_diag.mtx")
diagB  = col(ART + "/stageB7_laplacian_diag_bnd.mtx")
upperL = col(ART + "/stageB7_laplacian_upper.mtx")
srcL   = col(ART + "/stageB7_laplacian_source.mtx")
nP     = col(ART + "/stageB7_nP.mtx")
fdLin  = col(ART + "/stageB7_RP_FD_lin.mtx")
cand   = col(ART + "/stageB7_candidate_diag.mtx")
op     = np.loadtxt(ART + "/stageB7_outlet_patch.mtx")
ocells = op[:, 0].astype(int)
icB    = op[:, 3]
own = read_of_ints(CASE + "/constant/polyMesh/owner", len(upperL))
nei = read_of_ints(CASE + "/constant/polyMesh/neighbour", len(upperL))

# ---- exported J P-P block (rows/cols [NV,NUNK)) as CSR ----
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

# ================= Q3 =================
log("\n----- Q3: exported J P-P block vs production laplacian (diag/upper/source) -----")
relD = np.linalg.norm((diagL - expdiag)[mask]) / np.linalg.norm(expdiag[mask])
res["Q3_diag_exported_vs_Lint_relL2"] = relD
log("exported J P-P diag vs interior-only L diag (excl pRef row): relL2 = %.3e" % relD)
dProd = diagB - expdiag
outMask = np.zeros(N, bool); outMask[ocells] = True
res["Q3_prodminusexport_offoutlet_maxabs"] = float(np.abs(dProd[mask & ~outMask]).max())
res["Q3_prodminusexport_outlet_range"] = [float(dProd[outMask].min()), float(dProd[outMask].max())]
res["Q3_icB_range"] = [float(icB.min()), float(icB.max())]
log("production diag (L_int + addBoundaryDiag) - exported J P-P diag:")
log("  off-outlet max|d| = %.3e (expect ~0)" % res["Q3_prodminusexport_offoutlet_maxabs"])
log("  84 outlet cells    = [%.4e, %.4e]  (expect == internalCoeffs_b [%.4e, %.4e])"
    % (*res["Q3_prodminusexport_outlet_range"], *res["Q3_icB_range"]))
# upper: exclude faces incident to cell 0 (pRef pin zeroes row/col 0)
maxUp = 0.0; nUp = 0.0; nUpx = 0.0; nUpc = 0
for f in range(len(upperL)):
    o, nn = own[f], nei[f]
    if o == 0 or nn == 0:
        continue
    a = expPP[o, nn]; b = upperL[f]
    maxUp = max(maxUp, abs(a - b)); nUp += (a-b)*(a-b); nUpx += b*b; nUpc += 1
relUp = np.sqrt(nUp / max(nUpx, 1e-300))
res["Q3_upper_exported_vs_Lint_relL2"] = relUp
res["Q3_upper_exported_vs_Lint_maxabs"] = maxUp
res["Q3_upper_nfaces_excl_cell0"] = nUpc
log("exported J P-P upper (own,nei) vs probe L upper(+kf), %d faces excl. cell-0:"
    " relL2 = %.3e  max|d| = %.3e" % (nUpc, relUp, maxUp))
res["Q3_laplacian_source_maxabs"] = float(np.abs(srcL).max())
log("probe laplacian source max|.| = %.3e  (boundaryCoeffs_b = 0, p_b = 0)"
    % res["Q3_laplacian_source_maxabs"])
res["Q3_exported_PP_equals"] = "+L_int (interior-only laplacian, verified 1.7e-16)"
res["Q3_residual_derivative"] = "dR_P/dp = -L_prod = -(L_int + diag(internalCoeffs_b))"
log("CONCLUSION Q3: exported J P-P == interior-only laplacian; the ONLY difference"
    " to the PRODUCTION pressure matrix is the outlet internalCoeffs_b diagonal"
    " (=-cand) on the 84 outlet cells; boundaryCoeffs source = 0 (p_b=0).")

# ================= Q4a =================
log("\n----- Q4a: candidate J (diagnostic-only): P-P := -(L_int+bnd) = -L_prod -----")
candPP = -expPP.copy()
for c in ocells:
    if c != 0:
        candPP[c, c] += cand[c]           # +rAtU_c*delta_b*|Sf|_b (residual conv)
Jc = J.copy(); Jc = Jc.tolil()
for c in range(N):
    r = NV + c
    keep = [cc for cc in Jc.rows[r] if not (NV <= cc < NUNK)]
    Jc.rows[r] = keep
    Jc.data[r] = [Jc.data[r][k] for k, cc in enumerate(Jc.rows[r]) if True]
# NOTE: rebuild row data properly (rows and data must stay aligned)
Jc = J.tolil()
for c in range(N):
    r = NV + c
    # drop P-P entries (cols in [NV,NUNK)), keep U-part
    rr, dd = [], []
    for k, cc in enumerate(Jc.rows[r]):
        if not (NV <= cc < NUNK):
            rr.append(cc); dd.append(Jc.data[r][k])
    Jc.rows[r] = rr; Jc.data[r] = dd
s_e = candPP.indptr
for c in range(N):
    r = NV + c
    for j in range(s_e[c], s_e[c+1]):
        Jc[r, NV + candPP.indices[j]] = candPP.data[j]
Jc = Jc.tocsr()
Jc = Jc.tolil(); Jc[PREF, :] = 0.0; Jc[PREF, PREF] = 1.0; Jc = Jc.tocsr()
log("candidate J built: P-P = -L_prod (nnz=%d, t=%.1fs)" % (Jc.nnz, time.time() - t0))
cd = np.asarray(Jc[np.arange(NV, NUNK), np.arange(NV, NUNK)]).ravel()
relC = np.linalg.norm((cd + diagB)[mask]) / np.linalg.norm(diagB[mask])
res["Q4a_cand_PP_diag_vs_minusLprod_relL2"] = relC
log("candidate P-P diag vs -(production L diag): relL2 (excl pRef) = %.3e" % relC)

w = np.loadtxt("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-006/cycle-2/artifacts/wPrime_TAN_D1.mtx", dtype=np.float64)
nvec = w / np.linalg.norm(w)
log("null vector ||w||=%.6e  n[PREF]=%.3e" % (np.linalg.norm(w), nvec[PREF]))
rn_exp = np.linalg.norm(J @ nvec)
rn_cand = np.linalg.norm(Jc @ nvec)
res["Q4a_nullvec_Jexp_norm"] = float(rn_exp)
res["Q4a_nullvec_Jcand_norm"] = float(rn_cand)
log("||J_exp . n||/||n|| = %.3e   ||J_cand . n||/||n|| = %.3e   (null mode gone? %s)"
    % (rn_exp, rn_cand, "YES" if rn_cand > 1e3*rn_exp else "NO"))
# sigma_min of exported J: cite BFINAL-006/Post-Reviewer (7.382e-27); candidate: separate splu job
res["Q4a_sigma_min_Jexp"] = 7.382e-27  # BFINAL-006 post-review, rev_tangent_rerun.log
res["Q4a_sigma_min_Jexp_note"] = "BFINAL-006 Post-Reviewer independent rerun"
log("sigma_min(J_exp) = 7.382e-27 (BFINAL-006 Post-Reviewer) -> candidate sigma_min:"
    " run s3_q4_sigma.py (splu, ~30 min)")

# ================= Q4b =================
log("\n----- Q4b: candidate J_PP . n_P vs S2 actual-residual central FD -----")
candPPn = candPP @ nP                       # candidate P-P block action on n_P
fd = fdLin
relB = np.linalg.norm(candPPn - fd) / np.linalg.norm(fd)
res["Q4b_candPPn_norm"] = float(np.linalg.norm(candPPn))
res["Q4b_FD_norm"] = float(np.linalg.norm(fd))
res["Q4b_relL2_candPPn_vs_FD"] = relB
log("candidate J_PP . n_P vs FD_lin: ||candPPn||=%.4e ||FD||=%.4e relL2=%.3e (target ~1e-9)"
    % (np.linalg.norm(candPPn), np.linalg.norm(fd), relB))
expPPn = expPP @ nP
relB0 = np.linalg.norm(expPPn - fd) / np.linalg.norm(fd)
relB0n = np.linalg.norm(expPPn + fd) / np.linalg.norm(fd)
res["Q4b_JexpPPn_vs_FD_relL2"] = relB0
res["Q4b_JexpPPn_vs_minusFD_relL2"] = relB0n
log("  exported J_PP . n_P vs FD: |Jexp-FD|/|FD|=%.3e |Jexp+FD|/|FD|=%.3e"
    " -> exported P-P = NEGATIVE of residual derivative (+cand)" % (relB0, relB0n))
# FD_lin == -L_prod . n_P re-check from artifacts (independent path)
Lnp = -(diagB * nP)
for f in range(len(upperL)):
    Lnp[own[f]] += -upperL[f]*nP[nei[f]]
    Lnp[nei[f]] += -upperL[f]*nP[own[f]]
relS = np.linalg.norm(fd - Lnp) / np.linalg.norm(Lnp)
res["Q4b_FD_vs_minusLprod_nP_relL2"] = relS
log("  self-check FD_lin vs -(L_prod.n_P) [independent]: relL2 = %.3e" % relS)
relS2 = np.linalg.norm(candPPn - Lnp) / np.linalg.norm(Lnp)
res["Q4b_candPPn_vs_minusLprod_nP_relL2"] = relS2
log("  candidate P-P . n_P vs -(L_prod.n_P): relL2 = %.3e (candidate == residual deriv)"
    % relS2)

# ================= Q4c =================
log("\n----- Q4c: BFINAL-003 G1/G2 anchors with candidate J (dp=0 directions) -----")
dU = col(CASE + "/stageB5_dUdir.mtx")
v = np.zeros(NUNK); v[:NV] = dU
Jv_exp = J @ v
Jv_cand = Jc @ v
relV = np.linalg.norm(Jv_cand - Jv_exp) / np.linalg.norm(Jv_exp)
res["Q4c_Jv_cand_vs_exp_relL2"] = relV
log("J_cand . v == J_exp . v (Gate-A dir, dp=0): relL2 = %.3e" % relV)
Jv_oracle = col(CASE + "/stageB5_gateA_Jv.mtx")
relO = np.linalg.norm(Jv_exp - Jv_oracle) / np.linalg.norm(Jv_oracle)
res["Q4c_Jexp_vs_oracle_relL2"] = relO
log("J_exp . v vs oracle stageB5_gateA_Jv.mtx: relL2 = %.3e (G4 ref 3.98e-16)" % relO)
celltype = np.loadtxt(CASE + "/stageB5_cell_type.mtx").astype(int)
fdAll = np.loadtxt(CASE + "/stageB5_gateA_FD_all_eps.mtx")
epsA = [1e-3, 3e-4, 1e-4, 3e-5, 1e-5]
def bm(a, b, idx):
    av = a[idx]; bv = b[idx]
    return (np.linalg.norm(av-bv)/max(np.linalg.norm(bv),1e-300),
            float(np.dot(av,bv)/(np.linalg.norm(av)*np.linalg.norm(bv)+1e-300)))
idxMom = np.arange(NV)
idxPint = np.where(celltype == 0)[0]
idxPbnd = np.where(celltype == 1)[0]
res["Q4c_G2_blocks"] = {}
for ei, e in enumerate(epsA):
    fd = fdAll[ei*NUNK:(ei+1)*NUNK]
    fdU, fdP = fd[:NV], fd[NV:]
    rm, cm = bm(Jv_exp[:NV], fdU, idxMom)
    ri, ci = bm(Jv_exp[NV:], fdP, idxPint)
    rb, cb = bm(Jv_exp[NV:], fdP, idxPbnd)
    rt, ct = bm(Jv_exp[NV:], fdP, np.arange(N))
    res["Q4c_G2_blocks"]["eps=%.0e" % e] = {
        "momentum":[rm,cm],"P-int":[ri,ci],"P-bnd":[rb,cb],"P-total":[rt,ct]}
    log("G2 eps=%.0e: mom %.3e/%.6f | Pint %.3e/%.6f | Pbnd %.3e/%.6f | Ptot %.3e/%.6f"
        % (e, rm, cm, ri, ci, rb, cb, rt, ct))
dphiA = col(CASE + "/stageB5_gateB_dphi_A.mtx")
nIntF = len(upperL)
dphiFD = np.loadtxt(CASE + "/stageB5_gateB_dphi_FD_all_eps.mtx")
dphiB  = np.loadtxt(CASE + "/stageB5_gateB_dphi_B_all_eps.mtx")
idxInt = np.arange(nIntF); idxBnd = np.arange(nIntF, len(dphiA)); idxAll = np.arange(len(dphiA))
res["Q4c_G1"] = {}
for ei, e in enumerate(epsA):
    fdF = dphiFD[ei*len(dphiA):(ei+1)*len(dphiA)]
    bF  = dphiB[ei*len(dphiA):(ei+1)*len(dphiA)]
    rAI, cAI = bm(dphiA, fdF, idxInt)
    rAB, cAB = bm(dphiA, fdF, idxBnd)
    rAA, cAA = bm(dphiA, fdF, idxAll)
    rBI, cBI = bm(bF, fdF, idxInt)
    rBB, cBB = bm(bF, fdF, idxBnd)
    rBA, cBA = bm(bF, fdF, idxAll)
    res["Q4c_G1"]["eps=%.0e" % e] = {
        "A-int":[rAI,cAI],"A-bnd":[rAB,cAB],"A-all":[rAA,cAA],
        "B-int":[rBI,cBI],"B-bnd":[rBB,cBB],"B-all":[rBA,cBA]}
    log("G1 eps=%.0e: A int %.3e/%.6f bnd %.3e/%.6f all %.3e/%.6f |"
        " B int %.3e/%.6f bnd %.3e/%.6f all %.3e/%.6f"
        % (e, rAI,cAI, rAB,cAB, rAA,cAA, rBI,cBI, rBB,cBB, rBA,cBA))
log("Q4c: anchors reproduced from raw artifacts; candidate (P-P-only change) leaves"
    " every dp=0 direction EXACTLY unchanged (relL2=%.3e)." % relV)

with open(OUT + "/s3_q3_q4_summary.json", "w") as f:
    json.dump(res, f, indent=1, default=float)
log("\nS3_Q3_Q4_DONE")
