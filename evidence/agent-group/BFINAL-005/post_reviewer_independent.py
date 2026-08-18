#!/usr/bin/env python3
"""
BFINAL-005 Post-Reviewer INDEPENDENT recomputation (fresh invocation).

Reads ONLY the raw .mtx artifacts exported by the run and recomputes the
load-bearing G1/G2/G3 contractions from scratch. This is a separate code path
from both the C++ probe and the executor's verify_gates.py (different variable
naming, block indexing via explicit slices, and an independent dot-product
implementation using math.fsum for the vector inner products).

Semantics (from the BFINAL-004-locked reduced-SIMPLE oracle):
  lambda file : [3*nC  Uc-lambda ; nC  pc-lambda]
  rxc_analytic: 3 dirs x [3*nC  R_U,xd ; nC  R_P,xd]   (R_P,xd = J_P*(dAlphaDxh*z))
  z           : per dir, nC z = dxh/dx*d  (xh-tangent)
  dirs        : per dir, nC raw-design direction d
  oracle      : D_momentum = -sum(Uc . R_U,xd) ; D_pressure = -sum(pc * R_P,xd)
  production  : gsensh_pressurerow = -(J_P^T pc)*dAlphaDxh  (exported .mtx)
  G1          : sum(gsensh_pressurerow * z)   vs  D_pressure
  G2          : sum(gsensh_total * z)         vs  D_total
  G3          : sum(rxpr_prod_gsens * d)      vs  D_total
"""
import math

ART = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-005/cycle-2/step3_build_run/artifacts"
nC = 33600

def load(path, expect):
    vals = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.strip()
            if ln:
                vals.append(float(ln))
    assert len(vals) == expect, f"{path}: {len(vals)} != {expect}"
    return vals

lam  = load(f"{ART}/stageB6_lambda.mtx", 4*nC)
rxc  = load(f"{ART}/stageB6_rxc_analytic.mtx", 3*4*nC)
z    = load(f"{ART}/stageB6_rxc_z_analytic.mtx", 3*nC)
dirs = load(f"{ART}/stageB6_dirs.mtx", 3*nC)
gM   = load(f"{ART}/rxpr_prod_gsensh_momentum.mtx", nC)
gP   = load(f"{ART}/rxpr_prod_gsensh_pressurerow.mtx", nC)
gT   = load(f"{ART}/rxpr_prod_gsensh_total.mtx", nC)
gRaw = load(f"{ART}/rxpr_prod_gsens.mtx", nC)
gVol = load(f"{ART}/stageB6_prod_gsensVol.mtx", nC)

UcL = lam[0:3*nC]
pcL = lam[3*nC:4*nC]

def fdot(a, b):
    return math.fsum(x*y for x, y in zip(a, b))

def vecdot3(a, b):
    # a is 3*nC (interleaved x,y,z), b is 3*nC (interleaved x,y,z)
    return math.fsum(a[i]*b[i] for i in range(len(a)))

print("="*100)
print("BFINAL-005 Post-Reviewer independent recomputation (G1/G2/G3 from raw .mtx)")
print("="*100)

results = []
for di in range(3):
    base = di*4*nC
    RU = rxc[base:base+3*nC]
    RP = rxc[base+3*nC:base+4*nC]
    zDi = z[di*nC:(di+1)*nC]
    dDi = dirs[di*nC:(di+1)*nC]

    lU = vecdot3(UcL, RU)
    lP = fdot(pcL, RP)
    D_mom = -lU
    D_pre = -lP
    D_tot = D_mom + D_pre

    g1 = fdot(gP, zDi)          # production pressure-row (xh level)
    g2 = fdot(gT, zDi)          # production total (xh level)
    g3 = fdot(gRaw, dDi)        # production total (raw-design, post chain)
    gv = fdot(gVol, dDi)        # volume anchor

    def rel(a, b):
        return abs(a-b)/max(abs(b), 1e-300)

    results.append((di, D_mom, D_pre, D_tot, g1, g2, g3, gv, rel(g1, D_pre),
                    rel(g2, D_tot), rel(g3, D_tot)))
    print(f"--- D{di+1} ---")
    print(f"  D_momentum oracle = {D_mom:.16e}")
    print(f"  D_pressure oracle = {D_pre:.16e}")
    print(f"  D_total   oracle = {D_tot:.16e}")
    print(f"  G1 pressure-row prod = {g1:.16e}   relErr={rel(g1,D_pre):.3e}  signOK={math.copysign(1,g1)==math.copysign(1,D_pre)}")
    print(f"  G2 total prod       = {g2:.16e}   relErr={rel(g2,D_tot):.3e}  signOK={math.copysign(1,g2)==math.copysign(1,D_tot)}")
    print(f"  G3 raw-design prod  = {g3:.16e}   relErr={rel(g3,D_tot):.3e}  signOK={math.copysign(1,g3)==math.copysign(1,D_tot)}")
    print(f"  |D_pressure/D_total|= {abs(D_pre/D_tot):.6f}")
    print(f"  G5 volume proj     = {gv:.16e}")
    print()

print("="*100)
print("Decomposition self-check: max|gT-(gM+gP)| =",
      max(abs(gT[i]-(gM[i]+gP[i])) for i in range(nC)))
print("="*100)

# Locked BFINAL-004 anchors for comparison
anchors = [
    ("D1", -0.0562356333787, -0.00754887795106),
    ("D2", +0.297672431753, +0.0685927693909),
    ("D3", -0.0363414298395, -0.00621002821874),
]
print("\nAnchor reproduction check (D_total / D_pressure vs BFINAL-004 locked):")
for (di, Dm, Dp, Dt, g1, g2, g3, gv, re1, re2, re3), (nm, at, ap) in zip(results, anchors):
    print(f"  {nm}: total={Dt:.13e} (anchor {at})  pressure={Dp:.13e} (anchor {ap})")

# Gate verdicts
print("\nGATE VERDICTS (Post-Reviewer independent):")
allok = True
for (di, Dm, Dp, Dt, g1, g2, g3, gv, re1, re2, re3) in results:
    nm = f"D{di+1}"
    for gname, re_, sg in (("G1", re1, math.copysign(1,g1)==math.copysign(1,Dp)),
                           ("G2", re2, math.copysign(1,g2)==math.copysign(1,Dt)),
                           ("G3", re3, math.copysign(1,g3)==math.copysign(1,Dt))):
        ok = re_ <= 1e-3 and sg
        allok = allok and ok
        print(f"  {gname} {nm}: relErr={re_:.3e} signOK={sg} -> {'PASS' if ok else 'FAIL'}")
print(f"\nALL G1/G2/G3 PASS (independent): {allok}")
