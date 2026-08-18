#!/usr/bin/env python3
"""
BFINAL-005 Post-Reviewer (fixer round) INDEPENDENT recomputation.
Fresh pure-python code path (math.fsum), no shared code with C++ probe or
executor script. Reads ONLY the raw archived .mtx artifacts.
"""
import os, math

ART = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-005/cycle-2/fixer/artifacts"
nC = 33600

def load(name):
    vals = []
    with open(os.path.join(ART, name)) as f:
        for line in f:
            line = line.strip()
            if line:
                vals.append(float(line))
    return vals

lam  = load("stageB6_lambda.mtx")          # [3*nC Uc ; nC pc]
assert len(lam) == 4*nC, len(lam)
Uc = lam[0:3*nC]                           # adjoint velocity lambda_U (x,y,z interleaved)
pc = lam[3*nC:4*nC]                        # adjoint pressure lambda_P

rxc = load("stageB6_rxc_analytic.mtx")     # per dir: [3*nC R_U,xd ; nC R_P,xd]
assert len(rxc) == 12*nC, len(rxc)
z   = load("stageB6_rxc_z_analytic.mtx")   # nC z per dir (xh tangent)
assert len(z) == 3*nC, len(z)
dirs= load("stageB6_dirs.mtx")             # nC raw-design d per dir
assert len(dirs) == 3*nC, len(dirs)

gM = load("rxpr_prod_gsensh_momentum.mtx")
gP = load("rxpr_prod_gsensh_pressurerow.mtx")
gT = load("rxpr_prod_gsensh_total.mtx")
gRaw= load("rxpr_prod_gsens.mtx")
gVol= load("stageB6_prod_gsensVol.mtx")

# BFINAL-004 locked anchors (eps-independent analytic)
LOCK_D_mom = [-0.0486867554277, +0.229079662362, -0.0301314016208]
LOCK_D_pre = [-0.00754887795106, +0.0685927693909, -0.00621002821874]
LOCK_D_tot = [-0.0562356333787, +0.297672431753, -0.0363414298395]
LOCK_vol   = [0.1555186511144686, -0.08449667741658891, 0.06396237684432686]

def dot3(Uc, RU, nC):
    # Uc: list of 3*nC (x,y,z interleaved); RU: list of 3*nC
    return math.fsum(Uc[i]*RU[i] for i in range(3*nC))

print("=== (E) decomposition cellwise self-check ===")
dd = max(abs(gT[i] - (gM[i]+gP[i])) for i in range(nC))
print(f"max|total-(mom+pr)| = {dd:.3e}")

print("\n=== per-direction independent recomputation ===")
all_ok = True
for di in range(3):
    base = di*4*nC
    RU = rxc[base : base+3*nC]
    RP = rxc[base+3*nC : base+4*nC]
    lU = dot3(Uc, RU, nC)          # lambda_U^T R_U,x d
    lP = math.fsum(pc[i]*RP[i] for i in range(nC))   # lambda_P^T R_P,x d
    D_mom = -lU
    D_pre = -lP
    D_tot = -(lU + lP)
    zD = z[di*nC:(di+1)*nC]
    dD = dirs[di*nC:(di+1)*nC]
    g1 = math.fsum(gP[i]*zD[i] for i in range(nC))
    g2 = math.fsum(gT[i]*zD[i] for i in range(nC))
    g2m= math.fsum(gM[i]*zD[i] for i in range(nC))
    g3 = math.fsum(gRaw[i]*dD[i] for i in range(nC))
    gv = math.fsum(gVol[i]*dD[i] for i in range(nC))
    def rel(a,b): return abs(a-b)/max(abs(b),1e-300)
    checks = [
      ("D_momentum anchor", D_mom, LOCK_D_mom[di]),
      ("D_pressure anchor", D_pre, LOCK_D_pre[di]),
      ("D_total    anchor", D_tot, LOCK_D_tot[di]),
      ("G1 pressurerow vs oracle", g1, D_pre),
      ("G2 momentum   vs oracle", g2m, D_mom),
      ("G2 total      vs oracle", g2, D_tot),
      ("G3 rawdesign  vs oracle", g3, D_tot),
      ("G5 volume     vs anchor", gv, LOCK_vol[di]),
    ]
    nm = ["D1","D2","D3"][di]
    print(f"\n--- {nm} ---")
    for k,a,b in checks:
        r = rel(a,b)
        ok = r <= 1e-3
        all_ok = all_ok and ok
        print(f"  {k:28s} val={a:+.16e} ref={b:+.16e} relErr={r:.3e} {'OK' if ok else 'FAIL'}")
    print(f"  |D_pressure/D_total| = {abs(D_pre/D_tot):.6f}")

print("\n=== BFINAL-004 locked anchor reproduction (D_pressure oracle) ===")
for di in range(3):
    base = di*4*nC
    RP = rxc[base+3*nC : base+4*nC]
    lP = math.fsum(pc[i]*RP[i] for i in range(nC))
    print(f"  D{di+1}: D_pressure oracle = {-lP:+.16e}  (locked {-LOCK_D_pre[di]:+.16e})")

print(f"\nALL CHECKS PASS (<=1e-3): {all_ok}")
