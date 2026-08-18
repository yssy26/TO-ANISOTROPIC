#!/usr/bin/env python3
# BFINAL-004-FIX cycle-2 — Post-Reviewer independent recomputation (pure Python,
# no numpy, no shared formula with the C++ probe).  Reads ONLY the raw exported
# stageB6_*.mtx artifacts and recomputes relL2/cos/per-term decomposition for
# RX-A, RX-B, RX-C (xh-tangent + R_U,xd + R_P,xd) and the rpa_terms identity
# col1+col2==col3.
#
# Usage: python3 post_review_recompute_cycle2.py <artifact_dir>
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
    """relL2 = |a-b|_2/|b|_2 ; cos = <a,b>/(|a||b|). Independent of the probe's
    stageB6Metrics (which reduces over ranks and prints via Foam::Info)."""
    nA = sum(x*x for x in a)
    nB = sum(x*x for x in b)
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
    # ---------- RX-A ----------
    da = load(os.path.join(d, "stageB6_deltaAlpha.mtx"))
    N = len(da)
    an = load(os.path.join(d, "stageB6_rxa_analytic.mtx"))
    fd = load(os.path.join(d, "stageB6_rxa_FD_all_eps.mtx"))
    terms = load(os.path.join(d, "stageB6_rpa_terms.mtx"))
    assert len(an) == 4*N and len(fd) == 20*N and len(terms) == 3*N, \
        (len(an), len(fd), len(terms), N)
    print(f"# N={N}")
    print("== RX-A ==")
    # rpa_terms identity: [N divHbyA_w ; N -divflux_w ; N total_w] ; col1+col2==col3
    c1 = terms[0:N]; c2 = terms[N:2*N]; c3 = terms[2*N:3*N]
    maxid = max(abs(c1[i]+c2[i]-c3[i]) for i in range(N))
    print(f"rpa_terms identity max|col1+col2-col3|={maxid:.6e}")
    # col3 should equal the analytic R_P,alpha block of rxa_analytic
    anRP = an[3*N:4*N]
    maxid2 = max(abs(c3[i]-anRP[i]) for i in range(N))
    print(f"rpa_terms col3 == rxa_analytic R_P block max|diff|={maxid2:.6e}")
    # per-term norms
    n1 = math.sqrt(sum(x*x for x in c1)); n2 = math.sqrt(sum(x*x for x in c2))
    n3 = math.sqrt(sum(x*x for x in c3))
    print(f"weighted per-term: |div(dphiHbyA_w)|={n1:.16e} |-div(dflux_w)|={n2:.16e} |total_w|={n3:.16e}")
    anRU = an[0:3*N]
    for ei in range(5):
        eps = (1e-3, 3e-4, 1e-4, 3e-5, 1e-5)[ei]
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
    for ei in range(5):
        eps = (1e-3, 3e-4, 1e-4, 3e-5, 1e-5)[ei]
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
        # analytic alpha-tangent: dAlphaDxh*z (independent reconstruction)
        aA = [dadx[i]*zA[i] for i in range(N)]
        for ei in range(5):
            eps = (1e-3, 3e-4, 1e-4, 3e-5, 1e-5)[ei]
            fRU = fdr[(di*5+ei)*4*N: (di*5+ei)*4*N+3*N]
            fRP = fdr[(di*5+ei)*4*N+3*N: (di*5+ei+1)*4*N]
            zF = zFD[(di*5+ei)*N: (di*5+ei+1)*N]
            aF = aFD[(di*5+ei)*N: (di*5+ei+1)*N]
            metrics(zA, zF, f"  RX-C {name} xh-tangent eps={eps:g}")
            metrics(aA, aF, f"  RX-C {name} alpha-tangent eps={eps:g}")
            metrics(anRU, fRU, f"  RX-C {name} R_U,xd eps={eps:g}")
            metrics(anRP, fRP, f"  RX-C {name} R_P,xd eps={eps:g}")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "."))
