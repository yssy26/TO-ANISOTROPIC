#!/usr/bin/env python3
# BFINAL-004-FIX cycle-3 — executor verification: recompute every StageB6 metric
# from the raw exported stageB6_*.mtx artifacts with a pure-Python path
# (no numpy, no shared formula with the C++ probe).
#
# Usage: python3 verify_recompute_cycle3.py <artifact_dir>
#
# Verifies:
#   RX-A : weighted analytic vs FD (5 eps) for R_U,alpha and R_P,alpha;
#          rpa_terms identity col1+col2==col3 and col3==analytic R_P block;
#          weighted per-term norms.
#   RX-B : R_U,xh and R_P,xh (5 eps).
#   RX-C : per direction (D1/D2/D3) xh-tangent, alpha-tangent, R_U,xd, R_P,xd
#          (5 eps), recomputing the analytic alpha-tangent as dAlphaDxh*z.
#   §9/§10: per direction ||R_U,x d||, ||R_P,x d||, lambda-weighted
#          contributions (from stageB6_lambda.mtx), ratio.
import sys, math, os

def load(p):
    vals = []
    with open(p) as f:
        for line in f:
            s = line.strip()
            if s:
                vals.append(float(s))
    return vals

def metrics(a, b, name):
    nA = sum(x*x for x in a); nB = sum(x*x for x in b)
    nD = sum((x-y)*(x-y) for x, y in zip(a, b))
    dot = sum(x*y for x, y in zip(a, b))
    nA2 = math.sqrt(nA); nB2 = math.sqrt(nB)
    if nB2 <= 0.0 or nA2*nB2 <= 0.0:
        print(f"{name}: |a|={nA2:.16e} |b|={nB2:.16e} |err|={math.sqrt(nD):.16e} relL2=degenerate cos=degenerate")
        return None
    rel = math.sqrt(nD)/nB2
    cosv = dot/(nA2*nB2)
    print(f"{name}: |a|={nA2:.16e} |b|={nB2:.16e} |err|={math.sqrt(nD):.16e} relL2={rel:.12e} cos={cosv:.12e}")
    return rel

def main(d):
    da = load(os.path.join(d, "stageB6_deltaAlpha.mtx"))
    N = len(da)
    print(f"# N={N} (deltaAlpha.mtx)")

    # ---------- RX-A ----------
    an = load(os.path.join(d, "stageB6_rxa_analytic.mtx"))
    fd = load(os.path.join(d, "stageB6_rxa_FD_all_eps.mtx"))
    terms = load(os.path.join(d, "stageB6_rpa_terms.mtx"))
    assert len(an) == 4*N and len(fd) == 20*N and len(terms) == 3*N, \
        (len(an), len(fd), len(terms), N)
    print("== RX-A ==")
    c1 = terms[0:N]; c2 = terms[N:2*N]; c3 = terms[2*N:3*N]
    maxid = max(abs(c1[i]+c2[i]-c3[i]) for i in range(N))
    anRP = an[3*N:4*N]
    maxid2 = max(abs(c3[i]-anRP[i]) for i in range(N))
    n1 = math.sqrt(sum(x*x for x in c1)); n2 = math.sqrt(sum(x*x for x in c2))
    n3 = math.sqrt(sum(x*x for x in c3))
    print(f"rpa_terms: max|col1+col2-col3|={maxid:.6e} (machine zero expected)")
    print(f"rpa_terms: max|col3 - rxa_analytic R_P block|={maxid2:.6e} (0 expected)")
    print(f"weighted per-term: |div(dphiHbyA_w)|={n1:.16e} |-div(dflux_w)|={n2:.16e} |total_w|={n3:.16e}")
    anRU = an[0:3*N]
    for ei, eps in enumerate((1e-3, 3e-4, 1e-4, 3e-5, 1e-5)):
        fRU = fd[ei*4*N: ei*4*N+3*N]
        fRP = fd[ei*4*N+3*N: (ei+1)*4*N]
        metrics(anRU, fRU, f"  RX-A R_U,alpha eps={eps:g}")
        metrics(anRP, fRP, f"  RX-A R_P,alpha eps={eps:g}")

    # ---------- RX-B ----------
    print("== RX-B ==")
    anb = load(os.path.join(d, "stageB6_rxb_analytic.mtx"))
    fdb = load(os.path.join(d, "stageB6_rxb_FD_all_eps.mtx"))
    assert len(anb) == 4*N and len(fdb) == 20*N
    anRUb = anb[0:3*N]; anRPb = anb[3*N:4*N]
    for ei, eps in enumerate((1e-3, 3e-4, 1e-4, 3e-5, 1e-5)):
        fRU = fdb[ei*4*N: ei*4*N+3*N]
        fRP = fdb[ei*4*N+3*N: (ei+1)*4*N]
        metrics(anRUb, fRU, f"  RX-B R_U,xh eps={eps:g}")
        metrics(anRPb, fRP, f"  RX-B R_P,xh eps={eps:g}")

    # ---------- RX-C ----------
    print("== RX-C ==")
    anr = load(os.path.join(d, "stageB6_rxc_analytic.mtx"))
    fdr = load(os.path.join(d, "stageB6_rxc_FD_all_eps.mtx"))
    zAn = load(os.path.join(d, "stageB6_rxc_z_analytic.mtx"))
    zFD = load(os.path.join(d, "stageB6_rxc_xh_FD_all_eps.mtx"))
    aFD = load(os.path.join(d, "stageB6_rxc_alpha_FD_all_eps.mtx"))
    dadx = load(os.path.join(d, "stageB6_dalphadxh.mtx"))
    assert len(anr) == 12*N and len(fdr) == 60*N and len(zAn) == 3*N \
        and len(zFD) == 15*N and len(aFD) == 15*N, (len(anr), len(fdr), len(zAn), len(zFD), len(aFD))
    for di, name in enumerate(("D1", "D2", "D3")):
        anRU = anr[di*4*N: di*4*N+3*N]
        anRP = anr[di*4*N+3*N: (di+1)*4*N]
        zA = zAn[di*N: (di+1)*N]
        aA = [dadx[i]*zA[i] for i in range(N)]
        for ei, eps in enumerate((1e-3, 3e-4, 1e-4, 3e-5, 1e-5)):
            fRU = fdr[(di*5+ei)*4*N: (di*5+ei)*4*N+3*N]
            fRP = fdr[(di*5+ei)*4*N+3*N: (di*5+ei+1)*4*N]
            zF = zFD[(di*5+ei)*N: (di*5+ei+1)*N]
            aF = aFD[(di*5+ei)*N: (di*5+ei+1)*N]
            metrics(zA, zF, f"  RX-C {name} xh-tangent eps={eps:g}")
            metrics(aA, aF, f"  RX-C {name} alpha-tangent eps={eps:g}")
            metrics(anRU, fRU, f"  RX-C {name} R_U,xd eps={eps:g}")
            metrics(anRP, fRP, f"  RX-C {name} R_P,xd eps={eps:g}")

    # ---------- §9/§10 weighted contributions (analytic chain, eps-independent) ----------
    print("== §9/§10 ==")
    lam = load(os.path.join(d, "stageB6_lambda.mtx"))
    assert len(lam) == 4*N, len(lam)
    lamU = lam[0:3*N]; lamP = lam[3*N:4*N]
    for di, name in enumerate(("D1", "D2", "D3")):
        anRU = anr[di*4*N: di*4*N+3*N]
        anRP = anr[di*4*N+3*N: (di+1)*4*N]
        nRU = math.sqrt(sum(x*x for x in anRU))
        nRP = math.sqrt(sum(x*x for x in anRP))
        lU = sum(lamU[i]*anRU[i] for i in range(3*N))
        lP = sum(lamP[i]*anRP[i] for i in range(N))
        ratio = nRP/max(nRU, 1e-300)
        print(f"  {name}: ||R_U,x d||={nRU:.16e} ||R_P,x d||={nRP:.16e} "
              f"ratio={ratio:.6e} lambda_U^T R_U,x d={lU:.16e} "
              f"lambda_P^T R_P,x d={lP:.16e} D_momentum={-lU:.16e} "
              f"D_pressure={-lP:.16e} D_total={-(lU+lP):.16e} "
              f"|D_pressure/D_momentum|={abs(lP)/max(abs(lU),1e-300):.6e}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
