#!/usr/bin/env python3
# Post-Reviewer's OWN independent recomputation (not the executor's script).
# Reads raw stageB6_*.mtx from cycle-2/artifacts and recomputes:
#   (1) RX-A R_U,alpha / R_P,alpha relL2+cos across 5 eps
#   (2) rpa_terms identity col1+col2==col3 and col3==rxa_analytic R_P block
#   (3) RX-B R_U,xh / R_P,xh closure
#   (4) RX-C D1/D2/D3 xh/alpha tangent + R_U,xd/R_P,xd closure
#   (5) FD non-triviality: FD R_P,alpha must VARY with eps (true FD, not a copy)
#   (6) volume anchor projV_D1/D2/D3 from gsensVolR . dirs
#   (7) §9/§10 lambda-weighted D_momentum/D_pressure + prod gsensR.d
import sys, math

def load(p):
    v=[]
    with open(p) as f:
        for line in f:
            s=line.strip()
            if s: v.append(float(s))
    return v

def m(a,b):
    nA=sum(x*x for x in a); nB=sum(x*x for x in b)
    nD=sum((x-y)**2 for x,y in zip(a,b)); dot=sum(x*y for x,y in zip(a,b))
    nA=math.sqrt(nA); nB=math.sqrt(nB); nD=math.sqrt(nD)
    if nB==0 or nA*nB==0:
        return (None,None,nA,nB,nD)
    return (nD/nB, dot/(nA*nB), nA, nB, nD)

d = sys.argv[1] if len(sys.argv)>1 else "evidence/agent-group/BFINAL-004-FIX/cycle-2/artifacts"

da = load(f"{d}/stageB6_deltaAlpha.mtx")
N = len(da)
an = load(f"{d}/stageB6_rxa_analytic.mtx")
fd = load(f"{d}/stageB6_rxa_FD_all_eps.mtx")
terms = load(f"{d}/stageB6_rpa_terms.mtx")
assert len(an)==4*N and len(fd)==20*N and len(terms)==3*N, (len(an),len(fd),len(terms),N)
print(f"N={N}")
anRU = an[0:3*N]; anRP = an[3*N:4*N]

# (2) rpa_terms identity
c1=terms[0:N]; c2=terms[N:2*N]; c3=terms[2*N:3*N]
print("rpa_terms max|col1+col2-col3| =", max(abs(c1[i]+c2[i]-c3[i]) for i in range(N)))
print("rpa_terms col3 == rxa_analytic R_P max|diff| =", max(abs(c3[i]-anRP[i]) for i in range(N)))
n1=math.sqrt(sum(x*x for x in c1)); n2=math.sqrt(sum(x*x for x in c2)); n3=math.sqrt(sum(x*x for x in c3))
print(f"weighted terms |div(dphiHbyA_w)|={n1:.12e} |-div(dflux_w)|={n2:.12e} |total_w|={n3:.12e}")

print("\n== RX-A ==")
fdp_prev=None
for ei in range(5):
    eps=(1e-3,3e-4,1e-4,3e-5,1e-5)[ei]
    fRU=fd[ei*4*N:ei*4*N+3*N]; fRP=fd[ei*4*N+3*N:(ei+1)*4*N]
    ru=m(anRU,fRU); rp=m(anRP,fRP)
    print(f"  eps={eps:g} R_U,alpha relL2={ru[0]:.6e} cos={ru[1]:.12f} |a|={ru[2]:.8e} |b|={ru[3]:.8e}")
    print(f"            R_P,alpha relL2={rp[0]:.6e} cos={rp[1]:.12f} |a|={rp[2]:.8e} |b|={rp[3]:.8e}")
    # FD non-triviality: |b| of R_P must change across eps
    if fdp_prev is not None:
        print(f"            FD R_P |b| delta vs prev eps = {abs(rp[3]-fdp_prev):.3e}")
    fdp_prev=rp[3]

# (3) RX-B
anb=load(f"{d}/stageB6_rxb_analytic.mtx"); fdb=load(f"{d}/stageB6_rxb_FD_all_eps.mtx")
anRUb=anb[0:3*N]; anRPb=anb[3*N:4*N]
print("\n== RX-B ==")
for ei in range(5):
    eps=(1e-3,3e-4,1e-4,3e-5,1e-5)[ei]
    fRU=fdb[ei*4*N:ei*4*N+3*N]; fRP=fdb[ei*4*N+3*N:(ei+1)*4*N]
    ru=m(anRUb,fRU); rp=m(anRPb,fRP)
    print(f"  eps={eps:g} R_U,xh relL2={ru[0]:.6e} cos={ru[1]:.12f} ; R_P,xh relL2={rp[0]:.6e} cos={rp[1]:.12f} |a|={rp[2]:.8e} |b|={rp[3]:.8e}")

# (4) RX-C
anr=load(f"{d}/stageB6_rxc_analytic.mtx"); fdr=load(f"{d}/stageB6_rxc_FD_all_eps.mtx")
zAn=load(f"{d}/stageB6_rxc_z_analytic.mtx"); zFD=load(f"{d}/stageB6_rxc_xh_FD_all_eps.mtx")
aFD=load(f"{d}/stageB6_rxc_alpha_FD_all_eps.mtx"); dadx=load(f"{d}/stageB6_dalphadxh.mtx")
print("\n== RX-C (eps=1e-4 representative) ==")
for di,name in enumerate(("D1","D2","D3")):
    anRU=anr[di*4*N:di*4*N+3*N]; anRP=anr[di*4*N+3*N:(di+1)*4*N]
    zA=zAn[di*N:(di+1)*N]; aA=[dadx[i]*zA[i] for i in range(N)]
    for ei in range(5):
        eps=(1e-3,3e-4,1e-4,3e-5,1e-5)[ei]
        if eps!=1e-4: continue
        fRU=fdr[(di*5+ei)*4*N:(di*5+ei)*4*N+3*N]; fRP=fdr[(di*5+ei)*4*N+3*N:(di*5+ei+1)*4*N]
        zF=zFD[(di*5+ei)*N:(di*5+ei+1)*N]; aF=aFD[(di*5+ei)*N:(di*5+ei+1)*N]
        zt=m(zA,zF); at=m(aA,aF); ru=m(anRU,fRU); rp=m(anRP,fRP)
        print(f"  {name}: xh-tan relL2={zt[0]:.3e} cos={zt[1]:.6f} |a|={zt[2]:.6e} |b|={zt[3]:.6e}")
        print(f"         al-tan relL2={at[0]:.3e} cos={at[1]:.6f}")
        print(f"         R_U,xd relL2={ru[0]:.3e} cos={ru[1]:.6f} ; R_P,xd relL2={rp[0]:.3e} cos={rp[1]:.6f} |a|={rp[2]:.6e} |b|={rp[3]:.6e}")

# (6) volume anchor
gvol=load(f"{d}/stageB6_prod_gsensVol.mtx"); dirs=load(f"{d}/stageB6_dirs.mtx"); ct=load(f"{d}/stageB6_celltype.mtx")
print("\n== volume anchor (sum gsensVolR * Dk over active) ==")
for di,name in enumerate(("D1","D2","D3")):
    Dk=dirs[di*N:(di+1)*N]
    s=sum(gvol[i]*Dk[i] for i in range(N) if ct[i]>=1)
    print(f"  projV_{name} = {s:.16f}")

# (7) §9/§10 lambda weighted
lam=load(f"{d}/stageB6_lambda.mtx"); gsensR=load(f"{d}/stageB6_prod_gsens.mtx")
Uc=lam[0:3*N]; pc=lam[3*N:4*N]
print("\n== §9/§10 lambda-weighted (analytic chain) ==")
for di,name in enumerate(("D1","D2","D3")):
    anRU=anr[di*4*N:di*4*N+3*N]; anRP=anr[di*4*N+3*N:(di+1)*4*N]
    Dk=dirs[di*N:(di+1)*N]
    lU=sum(Uc[i]*anRU[i] for i in range(3*N))
    lP=sum(pc[i]*anRP[i] for i in range(N))
    prod=sum(gsensR[i]*Dk[i] for i in range(N) if ct[i]>=1)
    Dmom=-lU; Dpres=-lP
    print(f"  {name}: lU={lU:.10e} lP={lP:.10e} D_mom={Dmom:.10e} D_pres={Dpres:.10e}")
    print(f"         prod gsensR.d={prod:.10e}  |D_mom-prod|/|prod|={abs(Dmom-prod)/abs(prod):.3e}")
    print(f"         |D_pres/D_mom|={abs(Dpres/Dmom)*100:.2f}%")
