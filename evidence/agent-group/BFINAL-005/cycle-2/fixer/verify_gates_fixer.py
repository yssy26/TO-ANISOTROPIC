#!/usr/bin/env python3
"""
BFINAL-005 cycle-2 step4_verify: independent G1-G5 recomputation.

INDEPENDENT script: no shared formula with the C++ probe. It reads ONLY the
raw .mtx artifacts dumped by the run and recomputes every quantity from
scratch:

  oracle  D_pressure = -sum_cells( pc * R_P,xd )   [lambda_P^T R_P,x d]
  oracle  D_momentum = -sum_cells( Uc . R_U,xd )   [lambda_U^T R_U,x d]
  oracle  D_total    = D_momentum + D_pressure
  G1: D_pressure_production = sum_cells( gsensh_pressurerow * z )
  G2: D_total_production    = sum_cells( gsensh_total * z )
  G3: D_total_raw           = sum_cells( rxpr_prod_gsens * d )   (raw-design)
      volume projection     = sum_cells( gsensVolR * d )
  G4/G5: re-grep raw run log for the locked anchors.

File layouts (per cell count nC = 33600):
  stageB6_lambda.mtx        : [3*nC Uc-lambda ; nC pc-lambda]  (4*nC rows)
  stageB6_rxc_analytic.mtx  : for each direction di in 0..2:
                               [3*nC R_U,xd ; nC R_P,xd]       (4*nC rows)
  stageB6_rxc_z_analytic.mtx: nC z = dxh/dx*d per direction     (3*nC rows)
  stageB6_dirs.mtx          : nC D1, nC D2, nC D3               (3*nC rows)
  rxpr_prod_gsensh_momentum/pressurerow/total.mtx : nC each (xh-level, pre-chain)
  rxpr_prod_gsens.mtx       : nC raw-design post-chain gsensPressureDrop
  stageB6_prod_gsensVol.mtx : nC raw-design volume gradient
"""
import sys, math, os

ART = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-005/cycle-2/fixer/artifacts"
RUNLOG = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-005/cycle-2/fixer/run_switchON_fixer.log"
OUTDIR = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-005/cycle-2/fixer"

def load(name, n):
    p = os.path.join(ART, name)
    vals = []
    with open(p) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            vals.append(float(line))
    if len(vals) != n:
        raise SystemExit(f"{name}: expected {n} values, got {len(vals)}")
    return vals

def block(vals, start, n):
    return vals[start:start+n]

nC = 33600
lam = load("stageB6_lambda.mtx", 4*nC)
UcL = block(lam, 0, 3*nC)          # 3*nC
pcL = block(lam, 3*nC, nC)         # nC
rxc = load("stageB6_rxc_analytic.mtx", 3*4*nC)
z   = load("stageB6_rxc_z_analytic.mtx", 3*nC)
dirs= load("stageB6_dirs.mtx", 3*nC)
gM  = load("rxpr_prod_gsensh_momentum.mtx", nC)
gP  = load("rxpr_prod_gsensh_pressurerow.mtx", nC)
gT  = load("rxpr_prod_gsensh_total.mtx", nC)
gRaw= load("rxpr_prod_gsens.mtx", nC)
gVol= load("stageB6_prod_gsensVol.mtx", nC)

# sanity: production decomposition total == momentum + pressurerow (cell-wise)
maxDecompDiff = max(abs(gT[i] - (gM[i] + gP[i])) for i in range(nC))

def oracle_dir(di):
    """Independent recomputation of D_momentum/D_pressure/D_total for direction di."""
    base = di*4*nC
    RU = block(rxc, base, 3*nC)
    RP = block(rxc, base + 3*nC, nC)
    lU = sum(UcL[3*i+c]*RU[3*i+c] for i in range(nC) for c in range(3))
    lP = sum(pcL[i]*RP[i] for i in range(nC))
    return -lU, -lP, -(lU+lP)

def dot(a, b):
    return sum(a[i]*b[i] for i in range(len(a)))

def relerr(a, b):
    return abs(a-b)/max(abs(b), 1e-300)

rows = []
print("="*100)
print("BFINAL-005 cycle-2 step4_verify — independent G1-G5 (python, from raw .mtx)")
print("="*100)
print(f"decomposition self-check: max|total-(momentum+pressurerow)| = {maxDecompDiff:.3e}\n")

names = ["D1", "D2", "D3"]
table = []
for di in range(3):
    dD, dP, dT = oracle_dir(di)                 # independent oracle (recomputed)
    zDi = block(z, di*nC, nC)
    dDi = block(dirs, di*nC, nC)
    # G1: production pressure-row contraction with the analytic z (xh-level)
    g1_prod = dot(gP, zDi)
    # G2: production total contraction with z; also momentum component
    g2_prod_total = dot(gT, zDi)
    g2_prod_mom   = dot(gM, zDi)
    # G3: raw-design projection of the patched production field
    g3_raw = dot(gRaw, dDi)
    # G5: volume projection
    gVolProj = dot(gVol, dDi)
    table.append((names[di], dD, dP, dT, g2_prod_mom, g1_prod, g2_prod_total, g3_raw, gVolProj))
    print(f"--- {names[di]} (independent oracle recomputation) ---")
    print(f"  D_momentum oracle   = {dD: .16e}")
    print(f"  D_pressure oracle   = {dP: .16e}")
    print(f"  D_total oracle      = {dT: .16e}")
    print(f"  G1 D_pressure prod  = {g1_prod: .16e}  relErr={relerr(g1_prod,dP):.3e}  signOK={math.copysign(1,g1_prod)==math.copysign(1,dP)}")
    print(f"  G2 D_momentum prod  = {g2_prod_mom: .16e}  relErr={relerr(g2_prod_mom,dD):.3e}")
    print(f"  G2 D_total prod     = {g2_prod_total: .16e}  relErr={relerr(g2_prod_total,dT):.3e}  signOK={math.copysign(1,g2_prod_total)==math.copysign(1,dT)}")
    print(f"  G3 raw-design prod  = {g3_raw: .16e}  relErr={relerr(g3_raw,dT):.3e}  signOK={math.copysign(1,g3_raw)==math.copysign(1,dT)}")
    print(f"  |D_pressure/D_total|= {abs(dP/dT):.6f}")
    print(f"  G5 volume proj      = {gVolProj: .16e}")
    print()

# ---------- G4 / G5: re-grep raw run log for locked anchors ----------
def grep1(pattern, default=None):
    with open(RUNLOG) as f:
        for line in f:
            if pattern in line:
                return line.strip()
    return default

g4 = {}
g4["momentum_relL2"] = grep1("relL2=0.000165311363075")
g4["ptotal_relL2"]   = grep1("P-total: |a|L2=4.33741818206e-06")
g4["jdot"]           = grep1("Reduced cold-flow operator transpose dot-test (pressureDrop)")
g4["explicit"]       = grep1("ExplicitJToracle:")
g4["gatepr"]         = grep1("GatePR h=0.001 ")
g4["gatea_mom"]      = grep1("GateA eps=0.001")
g4["discrete"]       = grep1("Discrete objective derivative errors")
g4["anisok"]         = grep1("Anisotropic conductivity:")
g4["projeta"]        = grep1("StageB6 RX-C projection:")
g4["rxA_rp"]         = grep1("StageB6 RX-A eps=0.0001:StageB6  R_P,alpha:")
g4["rxB_rp"]         = grep1("StageB6 RX-B eps=0.0001:StageB6  R_P,xh:")
g4["rxc_d1_rp"]      = grep1("StageB6 RX-C D1 eps=0.0001:StageB6  R_P,xd:")
g4["dot_test"]       = grep1("RxPressureRowTranspose:")

print("="*100)
print("G4/G5 — locked anchors re-grepped from raw run log")
print("="*100)
for k, v in g4.items():
    print(f"  {k}: {v}")

# volume anchors (G5) — independent dot-product vs log prints
print()
print("G5 volume anchors (independent dot vs log):")
for di in range(3):
    logline = grep1(f"StageB6 volume-anchor {names[di]}:")
    print(f"  {names[di]}: independent={table[di][8]: .16e}  log={logline}")

# ---------- pass/fail summary ----------
print()
print("="*100)
print("GATE TABLE")
print("="*100)
ok = True
print(f"{'Gate':5s} {'dir':3s} {'production':>22s} {'oracle':>22s} {'relErr':>12s} {'sign':5s} {'PASS':5s}")
for di in range(3):
    n_, dD, dP, dT, gm, gp, gt, gr, gv = table[di]
    for gname, prod, orac in (("G1", gp, dP), ("G2", gt, dT), ("G3", gr, dT)):
        re_ = relerr(prod, orac)
        sg = math.copysign(1, prod) == math.copysign(1, orac)
        pass_ = (re_ <= 1e-3) and sg
        ok = ok and pass_
        print(f"{gname:5s} {n_:3s} {prod:22.16e} {orac:22.16e} {re_:12.3e} {str(sg):5s} {str(pass_):5s}")
    print(f"{'G5':5s} {n_:3s} volume proj = {gv:.16e} (anchor D1=0.1555186511144686)")

# G1 preferred 1e-4 / G2 preferred 1e-4 / G3 preferred 1e-4
pref_ok = True
for di in range(3):
    n_, dD, dP, dT, gm, gp, gt, gr, gv = table[di]
    if max(relerr(gp,dP), relerr(gt,dT), relerr(gr,dT)) > 1e-4:
        pref_ok = False
print(f"\npreferred (<=1e-4) achieved for all G1/G2/G3: {pref_ok}")
print(f"ALL GATES PASS (relErr<=1e-3 & sign): {ok}")
print(f"max|decomp-(mom+pr)| = {maxDecompDiff:.3e}")
