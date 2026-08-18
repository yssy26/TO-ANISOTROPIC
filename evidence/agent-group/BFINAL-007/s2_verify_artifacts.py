#!/usr/bin/env python3
"""
BFINAL-007 S2 analysis (read-only, offline): verify the (corrected) probe
artifacts.

Checks (S2 acceptance):
  1. outlet internalCoeffs_b ~ -rAtU*delta*|Sf| (<0), boundaryCoeffs_b ~ 0
  2. candidate diagonal +rAtU_c*deltaCoeffs_b*|Sf|_b per outlet cell
  3. Q1: actual-residual FD along n_P (production pEqn.flux() semantics,
     raw face-flux sum) is NON-ZERO and outlet-concentrated
  4. self-check: FD_lin == -div(flux(psi=n_P)) == -(L_prod . n_P) at machine
     precision (L_prod = addBoundaryDiag diag + upper; zeroGradient ic = 0,
     orthogonal mesh -> no non-orthogonal correction)
  5. all-patch boundary diag == outlet-only boundary diag (zeroGradient
     gradientInternalCoeffs = 0 verified from the export)
  6. eps robustness of the central-difference FD
  7. Q3 preview: production laplacian P-P block (diag/upper + outlet bnd diag)
     vs the exported explicitJT.mtx P-P block -> difference localized to the
     84 outlet cells
"""
import numpy as np, scipy.io, re

ART = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-007/artifacts"
CASE = "/home/ys/dsH/b2_case_smoke"
N = 33600
NV = 3 * N

def col(p):
    return np.loadtxt(p, dtype=np.float64)

diag = col(ART + "/stageB7_laplacian_diag.mtx")          # interior-only diag
diagB = col(ART + "/stageB7_laplacian_diag_bnd.mtx")     # + outlet addBoundaryDiag
diagBF = col(ART + "/stageB7_laplacian_diag_bnd_full.mtx")  # + all patches
upper = col(ART + "/stageB7_laplacian_upper.mtx")        # symmetric: lower == upper
nP = col(ART + "/stageB7_nP.mtx")
fd = col(ART + "/stageB7_RP_FD.mtx")                     # central FD (eps=1e-4)
fdLin = col(ART + "/stageB7_RP_FD_lin.mtx")              # exact -div(flux(n_P))
divFluxN = col(ART + "/stageB7_RP_divFluxN.mtx")         # div(flux(n_P))
fdEps = np.loadtxt(ART + "/stageB7_RP_FD_all_eps.mtx")   # [eps1e-2, eps1e-4, eps1e-6]
cand = col(ART + "/stageB7_candidate_diag.mtx")          # residual-conv bnd diag
op = np.loadtxt(ART + "/stageB7_outlet_patch.mtx")       # cell delta magSf ic bc mobB mobC candB
ocells = op[:, 0].astype(int)
icB, bcB = op[:, 3], op[:, 4]
bp = np.loadtxt(ART + "/stageB7_boundary_all_patches.mtx")  # patch cell delta magSf ic bc

print("== 1. outlet boundary coefficients ==")
print("internalCoeffs_b: min=%.4e max=%.4e mean=%.4e" %
      (icB.min(), icB.max(), icB.mean()))
print("boundaryCoeffs_b: min=%.4e max=%.4e mean=%.4e  (p_b=0 -> ~0)" %
      (bcB.min(), bcB.max(), bcB.mean()))
print("candidate rAtU_c*delta*|Sf|: min=%.4e max=%.4e mean=%.4e" %
      (cand[ocells].min(), cand[ocells].max(), cand[ocells].mean()))
rel = np.abs(icB + cand[ocells]) / np.maximum(np.abs(cand[ocells]), 1e-300)
print("|icB - (-cand)|/|cand|: max=%.3e  (expect ~0: icB = -rAtU*delta*|Sf|)"
      % rel.max())

print("\n== 1b. all-patch boundary audit ==")
for p in np.unique(bp[:, 0]):
    m = bp[:, 0] == p
    print("patch %d: nFaces=%d  ic: min=%.3e max=%.3e  bc: min=%.3e max=%.3e"
          % (p, int(m.sum()), bp[m, 4].min(), bp[m, 4].max(),
             bp[m, 5].min(), bp[m, 5].max()))
print("|diag_bnd_full - diag_bnd| (must be 0: only outlet contributes) = %.3e"
      % np.linalg.norm(diagBF - diagB))

print("\n== 2. Q1: actual-residual FD along n_P ==")
print("FD_lin: max|.|=%.4e  L2=%.4e  nnz(|FD|>1e-30)=%d" %
      (np.abs(fdLin).max(), np.linalg.norm(fdLin), int((np.abs(fdLin) > 1e-30).sum())))
fm = np.abs(fdLin)
frac_out = fm[ocells].sum() / fm.sum()
thr = np.percentile(fm, 99)
big = fm > thr
print("outlet (84 cells) |FD| mass fraction = %.4f" % frac_out)
print("top-1%% |FD| cells = %d, mass fraction = %.4f" % (int(big.sum()), fm[big].sum() / fm.sum()))
print("mean|FD| outlet=%.4e  global=%.4e  ratio=%.1f" %
      (fm[ocells].mean(), fm.mean(), fm[ocells].mean() / fm.mean()))

print("\n== 3. self-check: FD_lin == -div(flux(n_P)) == -(L_prod . n_P) ==")
def data_block(path):
    txt = open(path).read()
    m = re.search(r'\n(\d+)\s*\n\(', txt)
    return int(m.group(1)), txt[m.end():]

def read_ints(path):
    n, body = data_block(path)
    return np.array(
        [int(t) for t in body.replace('(', ' ').replace(')', ' ').split()[:n]]
    )

Lnp = diagB * nP
own = read_ints(CASE + "/constant/polyMesh/owner")[:len(upper)]
nei = read_ints(CASE + "/constant/polyMesh/neighbour")[:len(upper)]
np.add.at(Lnp, own, upper * nP[nei])
np.add.at(Lnp, nei, upper * nP[own])
print("|FD_lin - (-div(fluxN))|/|div(fluxN)| = %.3e" %
      (np.linalg.norm(fdLin + divFluxN) / np.linalg.norm(divFluxN)))
print("|FD_lin - (-L_prod.n_P)|/|L_prod.n_P| = %.3e  (expect ~1e-12..1e-14)"
      % (np.linalg.norm(fdLin + Lnp) / np.linalg.norm(Lnp)))
print("|FD_cd(1e-4) - FD_lin|/|FD_lin| = %.3e" %
      (np.linalg.norm(fd - fdLin) / np.linalg.norm(fdLin)))

print("\n== 4. eps robustness (central-difference FD) ==")
for i, e in enumerate([1e-2, 1e-4, 1e-6]):
    r = np.linalg.norm(fdEps[:, i] - fd) / np.linalg.norm(fd)
    print("eps=%.0e: |FD(eps)-FD(1e-4)|/|FD(1e-4)| = %.3e" % (e, r))

print("\n== 5. Q3 preview: production laplacian vs exported J P-P block ==")
M = scipy.io.mmread(CASE + "/explicitJT.mtx").tocsr()
expPPdiag = np.zeros(N)
for c in range(N):
    r = NV + c
    s, e = M.indptr[r], M.indptr[r + 1]
    for j in range(s, e):
        if M.indices[j] == r:
            expPPdiag[c] = M.data[j]
            break
# exclude the pRef identity row P(0) (=cell 0, diag 1.0) from the comparison
mask = np.ones(N, bool); mask[0] = False
relDiag = np.linalg.norm((diag - expPPdiag)[mask]) / np.linalg.norm(expPPdiag[mask])
print("|diag_laplacian - diag_exportedJ_PP|/|diag_exported| (excl. pRef row) = %.3e"
      % relDiag)
diff = np.abs(diag - expPPdiag)
outMask = np.zeros(N, bool); outMask[ocells] = True
interiorDiff = diff[mask & ~outMask]
outletDiff = diff[mask & outMask]
print("off-outlet cells: max|diag diff| = %.3e  (expect ~0: P-P blocks match)"
      % interiorDiff.max())
print("84 outlet cells: max|diag diff| = %.3e  (expect = |internalCoeffs_b|)"
      % outletDiff.max())
print("  |icB| on outlet: min=%.3e max=%.3e" % (np.abs(icB).min(), np.abs(icB).max()))
print("production bnd diag on outlet cells (internalCoeffs, laplacian conv): "
      "min=%.4e max=%.4e  -> exported J P-P missing this" %
      ((diagB[ocells] - diag[ocells]).min(), (diagB[ocells] - diag[ocells]).max()))
print("residual-convention missing term (candidate) on outlet cells: "
      "min=%.4e max=%.4e" % (cand[ocells].min(), cand[ocells].max()))
# PRODUCTION matrix (interior + addBoundaryDiag) vs exported J P-P:
dProd = diagB - expPPdiag
print("production diag (with addBoundaryDiag) minus exported J P-P diag:")
print("  off-outlet cells: max|diff| = %.3e  (expect ~0)" %
      np.abs(dProd[mask & ~outMask]).max())
print("  84 outlet cells : min=%.4e max=%.4e  (expect == internalCoeffs_b "
      "-1.7e-9..-5.5e-10)" % (dProd[outMask].min(), dProd[outMask].max()))
print("  => Q3: the ONLY P-P difference between the production pressure "
      "matrix and the exported J is the outlet internalCoeffs_b diagonal")

print("\nS2_ANALYSIS_DONE")
