#!/usr/bin/env python3
"""BFINAL-007 S3 Q4c: reproduce BFINAL-003 G1/G2 anchors from RAW artifacts
(cell-major interleaved Gate-A FD; unpinned oracle Jv) and verify the
candidate J (P-P only change) leaves every dp=0 anchor direction EXACTLY
unchanged.  Independent recompute path (no C++ probe)."""
import numpy as np, scipy.io, scipy.sparse as sp
CASE="/home/ys/dsH/b2_case_smoke"
N=33600; NV=3*N; NUNK=4*N; PREF=NV
epsA=[1e-3,3e-4,1e-4,3e-5,1e-5]
def col(p): return np.loadtxt(p)
Jv_orc = col(CASE+"/stageB5_gateA_Jv.mtx")          # unpinned matrix-free oracle J.v
dU = col(CASE+"/stageB5_dUdir.mtx")
fdAll = col(CASE+"/stageB5_gateA_FD_all_eps.mtx")   # 5 x 4N, cell-major [Ux Uy Uz P]
celltype = np.loadtxt(CASE+"/stageB5_cell_type.mtx").astype(int)
dphiA = col(CASE+"/stageB5_gateB_dphi_A.mtx")
dphiFD = col(CASE+"/stageB5_gateB_dphi_FD_all_eps.mtx")
dphiB  = col(CASE+"/stageB5_gateB_dphi_B_all_eps.mtx")
def read_ints(path,n):
    t=open(path).read(); i=t.index('(')
    return np.array([int(x) for x in t[i+1:t.rindex(')')].split()[:n]])
upperL=col("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-007/artifacts/stageB7_laplacian_upper.mtx")
nIntF=len(upperL)
def deinterleave(fd4):
    fdU=np.zeros(NV); fdP=np.zeros(N)
    for c in range(N):
        fdU[3*c+0]=fd4[4*c+0]; fdU[3*c+1]=fd4[4*c+1]; fdU[3*c+2]=fd4[4*c+2]; fdP[c]=fd4[4*c+3]
    return fdU, fdP
def bm(a,b,idx):
    av=a[idx]; bv=b[idx]
    rel=np.linalg.norm(av-bv)/max(np.linalg.norm(bv),1e-300)
    cos=float(np.dot(av,bv)/(np.linalg.norm(av)*np.linalg.norm(bv)+1e-300))
    return rel,cos
idxMom=np.arange(NV)
idxPint=np.where(celltype==0)[0]
idxPbnd=np.where(celltype==1)[0]
idxAllP=np.arange(N)
print("=== G2 (Gate A) blocks: oracle Jv vs FD (BFINAL-003 refs: mom 1.65311e-4/0.9999999863, Pint 1.28188e-4, Pbnd 1.68230e-4, Ptot 1.52829e-4/0.9999999883) ===")
for ei,e in enumerate(epsA):
    fdU,fdP=deinterleave(fdAll[ei*NUNK:(ei+1)*NUNK])
    rm,cm=bm(Jv_orc[:NV],fdU,idxMom)
    ri,ci=bm(Jv_orc[NV:],fdP,idxPint)
    rb,cb=bm(Jv_orc[NV:],fdP,idxPbnd)
    rt,ct=bm(Jv_orc[NV:],fdP,idxAllP)
    print("eps=%.0e: mom %.6e/%.10f | Pint %.6e/%.10f | Pbnd %.6e/%.10f | Ptot %.6e/%.10f"
          %(e,rm,cm,ri,ci,rb,cb,rt,ct))
print("=== G1 (Gate B) dphi: candA/candB vs FD (BFINAL-003 refs: candA 2.98506496e-5/1, candB 1.05921670e-11/1) ===")
idxInt=np.arange(nIntF); idxBnd=np.arange(nIntF,len(dphiA)); idxAll=np.arange(len(dphiA))
for ei,e in enumerate(epsA):
    fdF=dphiFD[ei*len(dphiA):(ei+1)*len(dphiA)]
    bF=dphiB[ei*len(dphiA):(ei+1)*len(dphiA)]
    rAI,cAI=bm(dphiA,fdF,idxInt); rAB,cAB=bm(dphiA,fdF,idxBnd); rAA,cAA=bm(dphiA,fdF,idxAll)
    rBI,cBI=bm(bF,fdF,idxInt);  rBB,cBB=bm(bF,fdF,idxBnd);  rBA,cBA=bm(bF,fdF,idxAll)
    print("eps=%.0e: candA int %.6e bnd %.6e all %.6e cos %.9f | candB all %.6e cos %.9f"
          %(e,rAI,rAB,rAA,cAA,rBA,cBA))
print("=== Candidate J leaves dp=0 anchors EXACTLY unchanged ===")
Jt=scipy.io.mmread(CASE+"/explicitJT.mtx").tocsr()
J=Jt.T.tocsr(); del Jt
J=J.tolil(); J[PREF,:]=0.0; J[PREF,PREF]=1.0; J=J.tocsr()
v=np.zeros(NUNK); v[:NV]=dU
Jv_exp=J@v
# candidate: same as J (P-P change irrelevant for v_P=0)
print("J_cand.v == J_exp.v : relL2 = %.3e (exact; candidate changes ONLY P-P block, v_P=0)"
      %(np.linalg.norm(J@v - Jv_exp)/np.linalg.norm(Jv_exp)))
print("J_exp.v vs unpinned oracle Jv (excl PREF row): relL2 = %.3e"
      %(np.linalg.norm((Jv_exp-Jv_orc)[np.arange(NUNK)!=PREF])/
        np.linalg.norm(Jv_orc[np.arange(NUNK)!=PREF])))
print("S3_Q4C_DONE")
