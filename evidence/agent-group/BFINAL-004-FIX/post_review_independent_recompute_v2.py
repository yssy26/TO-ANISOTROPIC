#!/usr/bin/env python3
# Post-Reviewer INDEPENDENT recompute (fresh invocation, no shared code).
# Reads ONLY the raw stageB6_*.mtx artifacts exported by the C++ probe and
# recomputes the load-bearing closure metrics from scratch:
#   - RX-A R_U,alpha / R_P,alpha  (B1: deltaAlpha-weighted closure)
#   - RX-B R_U,xh    / R_P,xh     (B3: xh->alpha closure, machine precision)
#   - RX-C D1/D2/D3  R_U,xd/R_P,xd, xh-tangent, alpha-tangent
#   - rpa_terms identity col1+col2==col3 (B3 per-term decomposition)
#   - lambda-weighted D_momentum/D_pressure and prod gsensR.d (sec9/10)
#   - volume anchor projV (sec8 regression anchor)
# Artifact layout (nC = total cells = 33600, nActive = 5040):
#   per-cell scalar : nC lines
#   U block (3N)    : 3*nC lines
#   lambda.mtx      : 3*nC (Uc) then nC (pc)
#   dirs.mtx        : 3*nC (D1,D2,D3)
#   rxa/rxb analytic: 3*nC (R_U) + nC (R_P)
#   rxa/rxb FD      : 5 eps * (3*nC + nC)
#   rxc analytic    : 3 dirs * (3*nC + nC)
#   rxc FD          : 3 dirs * 5 eps * (3*nC + nC)
#   rxc_xh/alpha FD : 3 dirs * 5 eps * nC
#   rxc_z/xh an     : 3 dirs * nC
#   rpa_terms       : nC (col1) + nC (col2) + nC (col3)
import sys, math, os

ART = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cycle-3", "artifacts")

def read_floats(path):
    with open(path) as f:
        return [float(x) for x in f.read().split()]

def l2(v):
    return math.sqrt(sum(x*x for x in v))

def dot(a, b):
    return sum(x*y for x, y in zip(a, b))

def rel_l2(a, b):
    return l2([x-y for x, y in zip(a, b)]) / l2(b)

def cosine(a, b):
    na, nb = l2(a), l2(b)
    return dot(a, b) / (na*nb) if na*nb > 0 else float("nan")

nC = 33600
N3 = 3*nC
EPS = [1e-3, 3e-4, 1e-4, 3e-5, 1e-5]

# ---- load artifacts ----
lam = read_floats(f"{ART}/stageB6_lambda.mtx")          # 4*nC
Uc = [lam[i] for i in range(N3)]
pc = [lam[N3+i] for i in range(nC)]

dirs = read_floats(f"{ART}/stageB6_dirs.mtx")           # 3*nC
D = [dirs[i*nC:(i+1)*nC] for i in range(3)]

rxa_an = read_floats(f"{ART}/stageB6_rxa_analytic.mtx") # 4*nC
rxa_fd = read_floats(f"{ART}/stageB6_rxa_FD_all_eps.mtx") # 5*4*nC
rxb_an = read_floats(f"{ART}/stageB6_rxb_analytic.mtx")
rxb_fd = read_floats(f"{ART}/stageB6_rxb_FD_all_eps.mtx")
rxc_an = read_floats(f"{ART}/stageB6_rxc_analytic.mtx")   # 3*4*nC
rxc_fd = read_floats(f"{ART}/stageB6_rxc_FD_all_eps.mtx") # 3*5*4*nC
rxc_xh_fd = read_floats(f"{ART}/stageB6_rxc_xh_FD_all_eps.mtx")    # 3*5*nC
rxc_alpha_fd = read_floats(f"{ART}/stageB6_rxc_alpha_FD_all_eps.mtx")
rxc_z_an = read_floats(f"{ART}/stageB6_rxc_z_analytic.mtx")        # 3*nC
rpa = read_floats(f"{ART}/stageB6_rpa_terms.mtx")          # 3*nC
gsens = read_floats(f"{ART}/stageB6_prod_gsens.mtx")       # nC
gsensVol = read_floats(f"{ART}/stageB6_prod_gsensVol.mtx") # nC

def split4(arr):
    """split a [3N U, N P] block into (Ublock, Pblock)."""
    return arr[:N3], arr[N3:]

def block_of(arr, idx, per=4*nC):
    return arr[idx*per:(idx+1)*per]

print("=== RX-A (alpha layer, deltaAlpha-weighted) ===")
for ei, eps in enumerate(EPS):
    au, ap = split4(rxa_an)
    fu, fp = split4(block_of(rxa_fd, ei))
    print(f" eps={eps:g}  R_U relL2={rel_l2(au,fu):.3e} cos={cosine(au,fu):.12f}"
          f" | R_P relL2={rel_l2(ap,fp):.3e} cos={cosine(ap,fp):.12f}"
          f" FD|R_P|={l2(fp):.6e}")

print("\n=== RX-B (xh layer) ===")
for ei, eps in enumerate(EPS):
    bu, bp = split4(rxb_an)
    gu, gp = split4(block_of(rxb_fd, ei))
    print(f" eps={eps:g}  R_U relL2={rel_l2(bu,gu):.3e} cos={cosine(bu,gu):.12f}"
          f" | R_P relL2={rel_l2(bp,gp):.3e} cos={cosine(bp,gp):.12f}"
          f" FD|R_P|={l2(gp):.6e}")

print("\n=== rpa_terms identity (B3 per-term decomposition) ===")
c1, c2, c3 = rpa[:nC], rpa[nC:2*nC], rpa[2*nC:3*nC]
maxdiff = max(abs(c1[i]+c2[i]-c3[i]) for i in range(nC))
print(f" max|col1+col2-col3| = {maxdiff:.3e}   |col1|={l2(c1):.6e}"
      f" |col2|={l2(c2):.6e} |col3|={l2(c3):.6e} (col1+col2 in L2? "
      f" {l2([c1[i]+c2[i] for i in range(nC)]):.6e})")

print("\n=== RX-C (raw-design D1/D2/D3) ===")
print(f"{'dir':<3}{'eps':>8} | {'xh relL2':>10} {'alpha relL2':>12} "
      f"{'R_U relL2':>11} {'R_U cos':>11} {'R_P relL2':>11} {'R_P cos':>11}")
for di in range(3):
    anU, anP = split4(block_of(rxc_an, di))
    for ei, eps in enumerate(EPS):
        fU, fP = split4(block_of(rxc_fd, di*5+ei))
        xh_fd = rxc_xh_fd[di*5*nC+ei*nC : di*5*nC+(ei+1)*nC]
        al_fd = rxc_alpha_fd[di*5*nC+ei*nC : di*5*nC+(ei+1)*nC]
        z_an = rxc_z_an[di*nC:(di+1)*nC]
        print(f"D{di+1} {eps:8g} | {rel_l2(z_an,xh_fd):10.3e}"
              f" {rel_l2(z_an,al_fd):12.3e}"
              f" {rel_l2(anU,fU):11.3e} {cosine(anU,fU):11.6f}"
              f" {rel_l2(anP,fP):11.3e} {cosine(anP,fP):11.6f}")

print("\n=== sec9/10: lambda-weighted + production completeness (eps-independent) ===")
for di in range(3):
    anU, anP = split4(block_of(rxc_an, di))
    d = D[di]
    nRU, nRP = l2(anU), l2(anP)
    lU = dot(Uc, anU)
    lP = dot(pc, anP)
    Dmom, Dpres = -lU, -lP
    prod = dot(gsens, d)
    print(f" D{di+1}: ||R_U||={nRU:.9f} ||R_P||={nRP:.9f} ratio={nRP/nRU:.3e}")
    print(f"   lambda_U^T R_U = {lU:+.10f}  lambda_P^T R_P = {lP:+.10f}")
    print(f"   D_momentum={Dmom:+.10f} D_pressure={Dpres:+.10f}"
          f" |D_pres/D_mom|={abs(Dpres/Dmom)*100:.2f}%")
    print(f"   prod gsensR.d = {prod:+.10f}  |Dmom-prod|/|prod|="
          f"{abs(Dmom-prod)/abs(prod):.3e}")

print("\n=== sec8 volume anchor (prod gsensVol . D) ===")
for di in range(3):
    projV = dot(gsensVol, D[di])
    print(f" D{di+1}: projV = {projV:.16f}")
