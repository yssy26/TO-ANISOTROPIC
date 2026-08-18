#!/usr/bin/env python3
# BFINAL-004 Post-Reviewer INDEPENDENT recompute (not *_myself.*, not the
# executor's verify_recompute_cycle3.py).  Written fresh for this invocation;
# parses ONLY the raw stageB6_*.mtx artifacts and recomputes the key closure
# metrics from first principles (pure Python, no numpy):
#
#   (a) RX-A  : weighted analytic vs central-FD for R_U,alpha and R_P,alpha
#               (B1 check: the analytic side exported to rxa_analytic.mtx IS
#                the deltaAlpha-weighted directional derivative).
#   (b) rpa_terms identity: col1+col2 == col3 per cell, and col3 == the
#               analytic R_P block of rxa_analytic.mtx.
#   (c) RX-B  : R_U,xh and R_P,xh analytic vs FD (B3 check).
#   (d) RX-C  : per direction D1/D2/D3, recompute the alpha-tangent gate as
#               aA[i]=dAlphaDxh[i]*zA[i] and compare to the FD alpha tangent,
#               then R_U,xd and R_P,xd analytic vs FD (>=1 R_U,x and >=1 R_P,x
#               direction per the task).
#   (e) §9/§10 : ||R_U,x d||, ||R_P,x d||, lambda_U^T R_U,x d,
#               lambda_P^T R_P,x d, |D_pressure/D_momentum| per direction.
#   (f) volume anchor: projV_Dk = sum_active gsensVolR*Dk recomputed from
#               prod_gsensVol.mtx x dirs.mtx over the celltype active mask;
#               and production projected contribution prod_gsens*d.
#
# Usage: python3 post_review_independent_recompute.py <artifact_dir>
import sys, os, math

def load(p):
    with open(p) as f:
        return [float(s) for s in f.read().split()]

def l2(v):
    return math.sqrt(sum(x*x for x in v))

def metrics(a, b, label):
    # relL2 = |a-b|_2/|b|_2 (b = FD reference), cosine = a.b/(|a||b|)
    nA, nB, nD, dot = 0.0, 0.0, 0.0, 0.0
    for x, y in zip(a, b):
        nA += x*x; nB += y*y; d = x-y; nD += d*d; dot += x*y
    nA2, nB2, nD2 = math.sqrt(nA), math.sqrt(nB), math.sqrt(nD)
    if nA2*nB2 == 0.0:
        cos = float('nan')
    else:
        cos = dot/(nA2*nB2)
    rel = (nD2/nB2) if nB2 > 0 else float('nan')
    print("  %-34s |a|L2=%.16e |b(FD)|L2=%.16e |err|L2=%.16e relL2=%.6e cos=%.12f"
          % (label, nA2, nB2, nD2, rel, cos))
    return rel, cos

def main(d):
    N = len(load(os.path.join(d, "stageB6_deltaAlpha.mtx")))
    print("# BFINAL-004 Post-Reviewer independent recompute")
    print("# artifact dir: %s   N=%d" % (d, N))
    EPS = (1e-3, 3e-4, 1e-4, 3e-5, 1e-5)

    # ---------------- (a) RX-A ----------------
    print("== RX-A (B1: analytic side must be the deltaAlpha-weighted directional derivative) ==")
    an = load(os.path.join(d, "stageB6_rxa_analytic.mtx"))
    fd = load(os.path.join(d, "stageB6_rxa_FD_all_eps.mtx"))
    terms = load(os.path.join(d, "stageB6_rpa_terms.mtx"))
    assert len(an) == 4*N and len(fd) == 20*N and len(terms) == 3*N, \
        (len(an), len(fd), len(terms), N)
    anRU = an[0:3*N]; anRP = an[3*N:4*N]
    # sanity: B1 weighting is NOT uniform -> weighted analytic must differ from
    # the unweighted R_U row of a unit-weight comparison; here we simply assert
    # the exported analytic matches the FD in direction (cos) and relL2.
    for ei, eps in enumerate(EPS):
        fRU = fd[ei*4*N: ei*4*N+3*N]
        fRP = fd[ei*4*N+3*N: (ei+1)*4*N]
        metrics(anRU, fRU, "RX-A R_U,alpha eps=%g" % eps)
        metrics(anRP, fRP, "RX-A R_P,alpha eps=%g" % eps)

    # ---------------- (b) rpa_terms identity ----------------
    print("== rpa_terms identity (col1+col2==col3 ; col3==analytic R_P block) ==")
    c1 = terms[0:N]; c2 = terms[N:2*N]; c3 = terms[2*N:3*N]
    maxid = max(abs(c1[i]+c2[i]-c3[i]) for i in range(N))
    maxid2 = max(abs(c3[i]-anRP[i]) for i in range(N))
    print("  max|col1+col2-col3|            = %.6e" % maxid)
    print("  max|col3 - analytic R_P block| = %.6e" % maxid2)
    print("  |div(dphiHbyA_w)|=%.16e  |-div(dflux_w)|=%.16e  |total_w|=%.16e"
          % (l2(c1), l2(c2), l2(c3)))
    print("  sum-check: |c1|+|c2| vs |c3| (triangle only, not the identity)")

    # ---------------- (c) RX-B ----------------
    print("== RX-B (B3: R_U,xh / R_P,xh closure) ==")
    anb = load(os.path.join(d, "stageB6_rxb_analytic.mtx"))
    fdb = load(os.path.join(d, "stageB6_rxb_FD_all_eps.mtx"))
    assert len(anb) == 4*N and len(fdb) == 20*N
    anRUb = anb[0:3*N]; anRPb = anb[3*N:4*N]
    for ei, eps in enumerate(EPS):
        fRU = fdb[ei*4*N: ei*4*N+3*N]
        fRP = fdb[ei*4*N+3*N: (ei+1)*4*N]
        metrics(anRUb, fRU, "RX-B R_U,xh eps=%g" % eps)
        metrics(anRPb, fRP, "RX-B R_P,xh eps=%g" % eps)

    # ---------------- (d) RX-C ----------------
    print("== RX-C (raw-design D1/D2/D3; alpha-tangent gate + R_U,xd/R_P,xd) ==")
    anr = load(os.path.join(d, "stageB6_rxc_analytic.mtx"))
    fdr = load(os.path.join(d, "stageB6_rxc_FD_all_eps.mtx"))
    zAn = load(os.path.join(d, "stageB6_rxc_z_analytic.mtx"))
    zFD = load(os.path.join(d, "stageB6_rxc_xh_FD_all_eps.mtx"))
    aFD = load(os.path.join(d, "stageB6_rxc_alpha_FD_all_eps.mtx"))
    dadx = load(os.path.join(d, "stageB6_dalphadxh.mtx"))
    assert len(anr) == 12*N and len(fdr) == 60*N and len(zAn) == 3*N \
        and len(zFD) == 15*N and len(aFD) == 15*N, \
        (len(anr), len(fdr), len(zAn), len(zFD), len(aFD))
    for di, name in enumerate(("D1", "D2", "D3")):
        anRU = anr[di*4*N: di*4*N+3*N]
        anRP = anr[di*4*N+3*N: (di+1)*4*N]
        zA = zAn[di*N: (di+1)*N]
        # alpha-tangent analytic: aA = dAlphaDxh * z (production interpolation
        # linearization), recomputed here independently.
        aA = [dadx[i]*zA[i] for i in range(N)]
        for ei, eps in enumerate(EPS):
            fRU = fdr[(di*5+ei)*4*N: (di*5+ei)*4*N+3*N]
            fRP = fdr[(di*5+ei)*4*N+3*N: (di*5+ei+1)*4*N]
            zF = zFD[(di*5+ei)*N: (di*5+ei+1)*N]
            aF = aFD[(di*5+ei)*N: (di*5+ei+1)*N]
            metrics(zA, zF, "RX-C %s xh-tangent eps=%g" % (name, eps))
            metrics(aA, aF, "RX-C %s alpha-tangent eps=%g" % (name, eps))
            metrics(anRU, fRU, "RX-C %s R_U,xd eps=%g" % (name, eps))
            metrics(anRP, fRP, "RX-C %s R_P,xd eps=%g" % (name, eps))

    # ---------------- (e) §9/§10 lambda-weighted ----------------
    print("== §9/§10 lambda-weighted (analytic chain, eps-independent) ==")
    lam = load(os.path.join(d, "stageB6_lambda.mtx"))
    assert len(lam) == 4*N, len(lam)
    lamU = lam[0:3*N]; lamP = lam[3*N:4*N]
    for di, name in enumerate(("D1", "D2", "D3")):
        anRU = anr[di*4*N: di*4*N+3*N]
        anRP = anr[di*4*N+3*N: (di+1)*4*N]
        lU = sum(lamU[i]*anRU[i] for i in range(3*N))
        lP = sum(lamP[i]*anRP[i] for i in range(N))
        ratio = l2(anRP)/l2(anRU)
        print("  %s: ||R_U,x d||=%.16e ||R_P,x d||=%.16e ratio=%.6e"
              % (name, l2(anRU), l2(anRP), ratio))
        print("      lambda_U^T R_U,x d=%.16e  lambda_P^T R_P,x d=%.16e"
              % (lU, lP))
        print("      D_momentum=%.16e D_pressure=%.16e D_total=%.16e"
              " |D_pres/D_mom|=%.6f"
              % (-lU, -lP, -(lU+lP), abs(lP)/max(abs(lU), 1e-300)))

    # ---------------- (f) volume anchor + production projected ----------------
    print("== volume anchor (prod_gsensVol x dirs over active cells) ==")
    gv = load(os.path.join(d, "stageB6_prod_gsensVol.mtx"))
    gr = load(os.path.join(d, "stageB6_prod_gsens.mtx"))
    ct = load(os.path.join(d, "stageB6_celltype.mtx"))
    dirs = load(os.path.join(d, "stageB6_dirs.mtx"))
    assert len(gv) == N and len(gr) == N and len(ct) == N and len(dirs) == 3*N
    active = [i for i in range(N) if ct[i] > 0.5]
    print("  active cells (celltype>0.5): %d" % len(active))
    for di, name in enumerate(("D1", "D2", "D3")):
        dvec = dirs[di*N: (di+1)*N]
        pV = sum(gv[i]*dvec[i] for i in active)
        pG = sum(gr[i]*dvec[i] for i in active)
        print("  %s: projV=%.16e   prod gsensR*d=%.16e" % (name, pV, pG))

    print("# END independent recompute")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "."))
