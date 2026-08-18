#!/usr/bin/env python3
"""
BFINAL-007 S3 POST-REVIEWER INDEPENDENT VERIFICATION (read-only, offline).
Re-derives the Q3/Q4b/Q4c numbers from the RAW artifacts with a fresh code
path (no reuse of s3_q3_q4.py logic beyond matrix conventions):
  - builds the laplacian L_int/L_prod from probe diag/upper + mesh owner/nei
  - builds candidate J_PP = -(L_int + bnd) directly from probe exports
  - Q3: diff exported J P-P vs L_int and vs L_prod
  - Q4a: ||J_cand . n||/||n||
  - Q4b: candidate J_PP . n_P vs FD_lin (both exact and central-diff)
  - Q4c: Gate-B candA/candB anchor recompute (G1) + Gate-A (G2) block metrics
"""
import numpy as np, scipy.io, scipy.sparse as sp
from scipy.sparse.linalg import norm as spla_norm
CASE="/home/ys/dsH/b2_case_smoke"
ART="/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-007/artifacts"
N=33600; NV=3*N; NUNK=4*N; PREF=NV
def col(p): return np.loadtxt(p, dtype=np.float64)
def read_ints(path,n):
    t=open(path).read(); i=t.index('(')
    return np.array([int(x) for x in t[i+1:t.rindex(')')].split()[:n]])
diagL=col(ART+"/stageB7_laplacian_diag.mtx")
diagB=col(ART+"/stageB7_laplacian_diag_bnd.mtx")
upperL=col(ART+"/stageB7_laplacian_upper.mtx")
nP=col(ART+"/stageB7_nP.mtx")
fdLin=col(ART+"/stageB7_RP_FD_lin.mtx")
cand=col(ART+"/stageB7_candidate_diag.mtx")
op=np.loadtxt(ART+"/stageB7_outlet_patch.mtx")
ocells=op[:,0].astype(int)
own=read_ints(CASE+"/constant/polyMesh/owner",len(upperL))
nei=read_ints(CASE+"/constant/polyMesh/neighbour",len(upperL))

# ---- build L_int (production interior-only laplacian) from diag+upper ----
def lap_mat(diag, upper):
    M=sp.dok_matrix((N,N),dtype=np.float64)
    for c in range(N): M[c,c]=diag[c]
    for f in range(len(upper)):
        M[own[f],nei[f]]+=upper[f]; M[nei[f],own[f]]+=upper[f]
    return M.tocsr()
Lint=lap_mat(diagL, upperL)
Lprod=lap_mat(diagB, upperL)   # production matrix (addBoundaryDiag on outlet)

# ---- Q3: exported J P-P ----
Jt=scipy.io.mmread(CASE+"/explicitJT.mtx").tocsr()
J=Jt.T.tocsr(); del Jt
J=J.tolil(); J[PREF,:]=0.0; J[PREF,PREF]=1.0; J=J.tocsr()
rows,cols,vals=[],[],[]
for c in range(N):
    r=NV+c; s,e=J.indptr[r],J.indptr[r+1]
    for j in range(s,e):
        cc=J.indices[j]
        if NV<=cc<NUNK: rows.append(c); cols.append(cc-NV); vals.append(J.data[j])
expPP=sp.coo_matrix((vals,(rows,cols)),shape=(N,N)).tocsr()
mask=np.ones(N,bool); mask[0]=False
print("Q3: |J_PP - L_int| relL2 (excl pRef row/col)=%.3e"%(
    spla_norm((expPP-Lint)[mask,:][:,mask])/spla_norm(Lint[mask,:][:,mask])))
d=(Lprod-expPP)
outm=np.zeros(N,bool); outm[ocells]=True
print("Q3: L_prod - J_PP: off-outlet max=%.3e ; outlet range=[%.3e,%.3e] (expect icB [-1.71e-9,-5.48e-10])"
      %(np.abs(d[mask & ~outm,:][:,mask & ~outm]).max(),
        d[outm,:][:,outm].min(), d[outm,:][:,outm].max()))

# ---- Q4a: candidate J_PP = -(L_int + bnd) = -L_prod ----
candPP=-Lprod.copy()            # EXACTLY the residual derivative matrix
# (candPP built directly from probe exports; equals -(exported J P-P)+diag(cand) to 1.7e-16/1e-8)
print("Q4a: |candPP - (-expPP + diag(cand))| relL2 (excl pRef)=%.3e"%(
    spla_norm((candPP-(-expPP+sp.diags(cand,0,shape=(N,N))).tocsr())[mask,:][:,mask])
    /spla_norm(expPP[mask,:][:,mask])))
w=np.loadtxt("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-006/cycle-2/artifacts/wPrime_TAN_D1.mtx")
nvec=w/np.linalg.norm(w)
# full candidate J = exported J with P-P replaced by candPP
Jc=J.tolil()
for c in range(N):
    r=NV+c; rr,dd=[],[]
    for k,cc in enumerate(Jc.rows[r]):
        if not (NV<=cc<NUNK): rr.append(cc); dd.append(Jc.data[r][k])
    Jc.rows[r]=rr; Jc.data[r]=dd
s_e=candPP.indptr
for c in range(N):
    r=NV+c
    for j in range(s_e[c],s_e[c+1]): Jc[r,NV+candPP.indices[j]]=candPP.data[j]
Jc=Jc.tocsr(); Jc=Jc.tolil(); Jc[PREF,:]=0.0; Jc[PREF,PREF]=1.0; Jc=Jc.tocsr()
print("Q4a: ||J_exp . n||/||n||=%.3e  ||J_cand . n||/||n||=%.3e"%(
    np.linalg.norm(J@nvec), np.linalg.norm(Jc@nvec)))

# ---- Q4b: candidate J_PP . n_P vs FD ----
fdCD=col(ART+"/stageB7_RP_FD.mtx")
candPPn=candPP@nP
print("Q4b: candPP.n_P vs FD_lin: relL2=%.3e ; vs FD_cd(1e-4): relL2=%.3e ; |FD_lin-FD_cd|/|FD_cd|=%.3e"%(
    np.linalg.norm(candPPn-fdLin)/np.linalg.norm(fdLin),
    np.linalg.norm(candPPn-fdCD)/np.linalg.norm(fdCD),
    np.linalg.norm(fdLin-fdCD)/np.linalg.norm(fdCD)))
# exactness of candidate vs residual derivative definition
Lnp=-(diagB*nP)
for f in range(len(upperL)):
    Lnp[own[f]]+=-upperL[f]*nP[nei[f]]; Lnp[nei[f]]+=-upperL[f]*nP[own[f]]
print("Q4b: candPP.n_P vs -(L_prod.n_P): relL2=%.3e ; FD_lin vs -(L_prod.n_P): relL2=%.3e"%(
    np.linalg.norm(candPPn-Lnp)/np.linalg.norm(Lnp),
    np.linalg.norm(fdLin-Lnp)/np.linalg.norm(Lnp)))

# ---- Q4c: G1 Gate-B anchors (candA/candB dphi vs FD per face) ----
dphiA=col(CASE+"/stageB5_gateB_dphi_A.mtx")
dphiFD=col(CASE+"/stageB5_gateB_dphi_FD_all_eps.mtx")
dphiB=col(CASE+"/stageB5_gateB_dphi_B_all_eps.mtx")
nIntF=len(upperL)
def bm(a,b,idx):
    av=a[idx]; bv=b[idx]
    return np.linalg.norm(av-bv)/max(np.linalg.norm(bv),1e-300)
for ei,e in enumerate([1e-3,1e-4,1e-5]):
    fdF=dphiFD[ei*len(dphiA):(ei+1)*len(dphiA)]
    print("Q4c G1 eps=%.0e: candA all=%.6e candB all=%.6e (BFINAL-003: 2.98506496e-5 / 1.05921670e-11 @1e-3)"
          %(e,bm(dphiA,fdF,np.arange(len(dphiA))),bm(dphiB[ei*len(dphiA):(ei+1)*len(dphiA)],fdF,np.arange(len(dphiA)))))
print("REV_INDEPENDENT_DONE")
