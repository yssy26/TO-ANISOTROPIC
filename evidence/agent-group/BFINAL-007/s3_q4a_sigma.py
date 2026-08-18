#!/usr/bin/env python3
"""BFINAL-007 S3 Q4a: sigma_min of the DIAGNOSTIC-ONLY candidate J (-L_prod P-P)
via random solves (splu), same method as BFINAL-006 t1_floor_diagnosis.py."""
import os, time, json
import numpy as np, scipy.io, scipy.sparse as sp, scipy.sparse.linalg as spla
CASE="/home/ys/dsH/b2_case_smoke"
ART="/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-007/artifacts"
OUT="/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-007"
N=33600; NV=3*N; NUNK=4*N; PREF=NV
t0=time.time()
def log(m): print(m,flush=True)
def col(p): return np.loadtxt(p)
log("=== S3 Q4a: sigma_min(candidate J) via random solves (diagnostic-only) ===")
Jt=scipy.io.mmread(CASE+"/explicitJT.mtx").tocsr()
J=Jt.T.tocsr(); del Jt
J=J.tolil(); J[PREF,:]=0.0; J[PREF,PREF]=1.0; J=J.tocsr()
cand=col(ART+"/stageB7_candidate_diag.mtx")
op=np.loadtxt(ART+"/stageB7_outlet_patch.mtx")
ocells=op[:,0].astype(int)
# candidate: P-P := -L_prod = -(L_int + bnd) ; build from exported J P-P
rows,cols,vals=[],[],[]
for c in range(N):
    r=NV+c; s,e=J.indptr[r],J.indptr[r+1]
    for j in range(s,e):
        cc=J.indices[j]
        if NV<=cc<NUNK: rows.append(c); cols.append(cc-NV); vals.append(J.data[j])
expPP=sp.coo_matrix((vals,(rows,cols)),shape=(N,N)).tocsr()
candPP=-expPP.copy()
for c in ocells:
    if c!=0: candPP[c,c]+=cand[c]
Jc=J.tolil()
for c in range(N):
    r=NV+c; rr,dd=[],[]
    for k,cc in enumerate(Jc.rows[r]):
        if not (NV<=cc<NUNK): rr.append(cc); dd.append(Jc.data[r][k])
    Jc.rows[r]=rr; Jc.data[r]=dd
s_e=candPP.indptr
for c in range(N):
    r=NV+c
    for j in range(s_e[c],s_e[c+1]):
        Jc[r,NV+candPP.indices[j]]=candPP.data[j]
Jc=Jc.tocsr()
Jc=Jc.tolil(); Jc[PREF,:]=0.0; Jc[PREF,PREF]=1.0; Jc=Jc.tocsr()
log("candidate J built (nnz=%d t=%.1fs)"%(Jc.nnz,time.time()-t0))
lu=spla.splu(Jc,permc_spec="COLAMD")
log("splu done (t=%.1fs)"%(time.time()-t0))
sig=[]
for k in range(3):
    rng=np.random.default_rng(1000+k)
    z=rng.standard_normal(NUNK)
    x=lu.solve(z)
    s=np.linalg.norm(z)/np.linalg.norm(x)
    sig.append(s)
    log("trial %d: ||z||=%.4e ||x||=%.4e sigma_min~%.3e"%(k,np.linalg.norm(z),np.linalg.norm(x),s))
log("sigma_min(candidate J): min=%.3e mean=%.3e"%(min(sig),float(np.mean(sig))))
res={"Q4a_sigma_min_Jcand_min":min(sig),"Q4a_sigma_min_Jcand_mean":float(np.mean(sig)),
     "Q4a_sigma_min_Jcand_trials":sig,"note":"diagnostic-only candidate, no production change"}
with open(OUT+"/s3_q4a_sigma.json","w") as f: json.dump(res,f,indent=1)
log("S3_Q4A_SIGMA_DONE")
