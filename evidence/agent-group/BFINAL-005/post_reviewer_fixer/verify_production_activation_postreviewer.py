#!/usr/bin/env python3
"""
BFINAL-005 Post-Reviewer (fixer round): INDEPENDENT production activation check.
Reads the SWITCH-OFF production-written 1/gsensPressureDrop (raw-design field)
and 1/gsenshPressureDrop (xh-level total), projects onto the deterministic
D1/D2/D3 (from the archived stageB6_dirs.mtx), and compares against the
BFINAL-004-locked D_total oracle. Proves production R_x now contains BOTH rows
with stageB6RxDesignOracle=false.
"""
import os, math

CASE = "/home/ys/dsH/b2_case_smoke_prod_fix"
ART = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-005/cycle-2/fixer/artifacts"
nC = 33600

def load_mtx(name):
    vals=[]
    with open(os.path.join(ART,name)) as f:
        for line in f:
            line=line.strip()
            if line: vals.append(float(line))
    return vals

def load_field(path):
    # parse OpenFOAM ascii volScalarField internalField
    vals=[]
    with open(path) as f:
        lines=f.read().splitlines()
    # find "internalField" and the count line, then read next nC numeric lines
    i=0
    while i < len(lines) and "internalField" not in lines[i]:
        i+=1
    i+=1  # internalField line
    # next non-empty lines: possibly "nonuniform List<scalar>" then count then "("
    while i < len(lines):
        ln=lines[i].strip()
        i+=1
        if ln.startswith("("):
            break
    cnt=0
    while i < len(lines) and cnt < nC:
        ln=lines[i].strip(); i+=1
        if ln.startswith(")"): break
        if ln=="": continue
        # a line may contain multiple numbers
        for tok in ln.split():
            if tok.startswith("("): tok=tok[1:]
            if tok.endswith(")"): tok=tok[:-1]
            if tok=="": continue
            vals.append(float(tok)); cnt+=1
    assert len(vals)==nC, f"expected {nC}, got {len(vals)}"
    return vals

dirs = load_mtx("stageB6_dirs.mtx")            # 3*nC
gRaw = load_field(CASE+"/1/gsensPressureDrop") # raw-design post-chain
gTot = load_field(CASE+"/1/gsenshPressureDrop")# xh-level total

# note: z (xh tangent) is needed for the xh-level contraction; load archived z
z = load_mtx("stageB6_rxc_z_analytic.mtx")

LOCK_D_tot = [-0.0562356333787, +0.297672431753, -0.0363414298395]
LOCK_D_mom = [-0.0486867554277, +0.229079662362, -0.0301314016208]

print("=== production (switch OFF) activation: projection onto D1/D2/D3 ===")
all_ok=True
for di in range(3):
    dD = dirs[di*nC:(di+1)*nC]
    zD = z[di*nC:(di+1)*nC]
    proj_raw = math.fsum(gRaw[i]*dD[i] for i in range(nC))     # raw-design
    proj_xh  = math.fsum(gTot[i]*zD[i] for i in range(nC))     # xh total . z
    def rel(a,b): return abs(a-b)/max(abs(b),1e-300)
    r_raw = rel(proj_raw, LOCK_D_tot[di])
    r_xh  = rel(proj_xh, LOCK_D_tot[di])
    ok = r_raw<=1e-3 and r_xh<=1e-3
    all_ok = all_ok and ok
    print(f"D{di+1}: raw-design proj={proj_raw:+.16e} (oracle total {LOCK_D_tot[di]:+.16e}) relErr={r_raw:.3e}")
    print(f"       xh-level   proj={proj_xh:+.16e} (oracle total {LOCK_D_tot[di]:+.16e}) relErr={r_xh:.3e}")
    print(f"       [momentum-only would be {LOCK_D_mom[di]:+.16e}]")

print(f"\nPRODUCTION CONTAINS BOTH ROWS (switch OFF): {all_ok}")
