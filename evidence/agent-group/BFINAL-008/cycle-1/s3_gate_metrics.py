#!/usr/bin/env python3
# BFINAL-008 S3 — independent recomputation of P1/P3/P4/P5/P6/P7 gate metrics
# from run-1 (probe path) raw artifacts + log, and run-2 (production path) log.
# Authoritative Python path (matches the S2 pattern): reads the .mtx exports and
# the raw run logs; writes a compact JSON summary.
import json, os, re, math, sys

BASE = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-008/cycle-1/artifacts"
LOG1 = "/home/ys/dsH/b8_scratch_r1/Log.stageB8.r1.txt"
LOG2 = "/home/ys/dsH/b8_scratch_prod/Log.stageB8.prod.txt"

def read_vec(path, n=None):
    with open(path) as f:
        vals = [float(x) for x in f.read().split()]
    return vals if n is None else vals[:n]

def read_celltype(path):
    with open(path) as f:
        return [int(x) for x in f.read().split()]

def block_metrics(a, b, idx):
    nA = nB = nD = dot = 0.0
    for i in idx:
        av, bv = a[i], b[i]
        d = av - bv
        nA += av*av; nB += bv*bv; nD += d*d; dot += av*bv
    if nB == 0.0:
        return {"relL2": None, "cos": None, "normRatio": None, "|Jv|": math.sqrt(nA), "|FD|": 0.0}
    return {
        "relL2": math.sqrt(nD/nB),
        "cos": dot/math.sqrt(nA*nB) if nA > 0 else None,
        "normRatio": math.sqrt(nA/nB),
        "|Jv|": math.sqrt(nA),
        "|FD|": math.sqrt(nB),
    }

out = {"HEAD": None, "P1": {}, "P3": {}, "P4": {}, "P5": {}, "P6": {}, "P7": {}, "run2": {}}

# ---------------- cell classification ----------------
ct = read_celltype(os.path.join(BASE, "stageB8_celltype.mtx"))
N = len(ct)
outlet = [i for i,c in enumerate(ct) if c == 2]
other = [i for i,c in enumerate(ct) if c == 1]
interior = [i for i,c in enumerate(ct) if c == 0]
offout = [i for i,c in enumerate(ct) if c != 2]
tot = list(range(N))
NV = 3*N
print(f"cells: N={N} internal={len(interior)} otherBoundary={len(other)} outletAdjacent={len(outlet)}")

# ---------------- P1 / P3 (dir1, eps=1e-2 authoritative) ----------------
for name, jvf, fdf, dirf in [("dir1","stageB8_Jv1.mtx","stageB8_FD1_eps_1e-2.mtx","stageB8_dir1.mtx"),
                             ("dir2","stageB8_Jv2.mtx","stageB8_FD2_eps_1e-2.mtx","stageB8_dir2.mtx"),
                             ("dir3","stageB8_Jv3.mtx","stageB8_FD3_eps_1e-2.mtx","stageB8_dir3.mtx")]:
    Jv = read_vec(os.path.join(BASE, jvf))
    FD = read_vec(os.path.join(BASE, fdf))
    Pidx = [NV + i for i in tot]
    Uidx = list(range(NV))
    Pint = [NV + i for i in interior]
    Pout = [NV + i for i in outlet]
    Poth = [NV + i for i in other]
    Poff = [NV + i for i in offout]
    m = {
        "U_rows": block_metrics(Jv, FD, Uidx),
        "P_internal": block_metrics(Jv, FD, Pint),
        "P_outlet_adjacent": block_metrics(Jv, FD, Pout),
        "P_other_boundary": block_metrics(Jv, FD, Poth),
        "P_total": block_metrics(Jv, FD, Pidx),
    }
    if name == "dir1":
        out["P1"]["full_P"] = m["P_total"]
        out["P1"]["outlet_adjacent"] = m["P_outlet_adjacent"]
        out["P1"]["off_outlet"] = block_metrics(Jv, FD, Poff)
        # sign/structure
        sPos = sum(1 for i in Pidx if Jv[i]*FD[i] > 0)
        fdOut = sum(abs(FD[i]) for i in Pout); fdTot = sum(abs(FD[i]) for i in Pidx)
        jvOut = sum(abs(Jv[i]) for i in Pout); jvTot = sum(abs(Jv[i]) for i in Pidx)
        out["P1"]["sameSignFraction"] = sPos/len(Pidx)
        out["P1"]["FD_mass_outlet_fraction"] = fdOut/fdTot
        out["P1"]["Jv_mass_outlet_fraction"] = jvOut/jvTot
        nP = read_vec(os.path.join(BASE, "stageB8_nP.mtx"))
        out["P1"]["||n_P||"] = math.sqrt(sum(v*v for v in nP))
        w = read_vec("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-006/cycle-2/artifacts/wPrime_TAN_D1.mtx", 4*N)
        out["P1"]["||w_TAN_D1||"] = math.sqrt(sum(v*v for v in w))
    out["P3"][name + "_eps1e-2"] = m

# ---------------- P1 eps plateau (full-P relL2 from probe metrics file) ----------------
metrics_txt = os.path.join(BASE, "stageB8_metrics.txt")
plateau = {}
with open(metrics_txt) as f:
    for line in f:
        p = line.split()
        if p and p[0] == "P1":
            plateau[p[1]] = float(p[2])
out["P1"]["full_P_relL2_eps_plateau"] = plateau

# ---------------- P4 BlockDot (from run-1 log) ----------------
log1 = open(LOG1).read() if os.path.exists(LOG1) else ""
block = re.findall(r"BlockDot \((\w+)\): PU=(\S+) PP=(\S+) UU=(\S+) UP=(\S+) boundaryU=(\S+)", log1)
out["P4"]["BlockDot"] = [{"label": b[0], "PU": b[1], "PP": b[2], "UU": b[3], "UP": b[4], "boundaryU": b[5]} for b in block]
red = re.findall(r"Reduced cold-flow operator transpose dot-test \((\w+)\) max relative error=(\S+)", log1)
out["P4"]["reduced_coldflow_dot"] = [{"label": r[0], "maxRelErr": r[1]} for r in red]

# ---------------- P5 ExplicitJToracle ----------------
tor = re.findall(r"ExplicitJToracle: maxRelL2=(\S+) maxCosErr=(\S+) maxRelL_U=(\S+) maxRelL_P=(\S+)", log1)
out["P5"]["ExplicitJToracle"] = [{"maxRelL2": t[0], "maxCosErr": t[1], "maxRelL_U": t[2], "maxRelL_P": t[3]} for t in tor]

# ---------------- P6 StageB5 G1/G2 ----------------
g1 = re.findall(r"StageB5 GateB eps=(\S+) candA: \|A\|L2=\S+ \|FD\|L2=\S+ normRatio=\S+ relL2=(\S+) cos=\S+ intRel=\S+ bndRel=\S+", log1)
g2 = re.findall(r"StageB5 GateB eps=(\S+) candB: \|B\|L2=\S+ \|FD\|L2=\S+ normRatio=\S+ relL2=(\S+) cos=\S+ intRel=\S+ bndRel=\S+", log1)
mom = re.findall(r"StageB5 GateA eps=(\S+):StageB5  momentum: \|a\|L2=\S+ \|b\(FD\)\|L2=\S+ \|err\|L2=\S+ relL2=(\S+) cos=\S+ maxDiff=\S+", log1)
out["P6"]["GateB_candA"] = [{"eps": g[0], "relL2": g[1]} for g in g1]
out["P6"]["GateB_candB"] = [{"eps": g[0], "relL2": g[1]} for g in g2]
out["P6"]["GateA_momentum"] = [{"eps": m[0], "relL2": m[1]} for m in mom]

# ---------------- P7 stageB6 R_x dot tests ----------------
rx = re.findall(r"RxPressureRowTranspose: sum\(T\*w\)=\"(\S+)\" sum\(pc\*\(J_P\*w\)\)=\"(\S+)\" relErr=(\S+)", log1)
out["P7"]["rxPressureRow_dot"] = [{"sumTw": r[0], "sumPcJPw": r[1], "relErr": r[2]} for r in rx]
d = re.findall(r"StageB6 weighted (D\d): \|\|R_U,x d\|\|L2=(\S+) \|\|R_P,x d\|\|L2=(\S+) ratio\|\|R_P\|\|/\|\|R_U\|\|=(\S+) lambda_U\^T R_U,x d=(\S+) lambda_P\^T R_P,x d=(\S+) D_momentum=(\S+) D_pressure=(\S+) D_total=(\S+) prodGsenDPressDrop\*d=(\S+) prodDgdx1\*d\(scaled\)=(\S+) \|D_mom-prod\|/\|prod\|=(\S+)", log1)
out["P7"]["stageB6_weighted"] = [{"dir": x[0], "|R_U,x d|L2": x[1], "|R_P,x d|L2": x[2], "ratio": x[3],
                                  "lambdaU_RUx": x[4], "lambdaP_RPx": x[5], "D_momentum": x[6],
                                  "D_pressure": x[7], "D_total": x[8], "prod": x[9], "prodDgdx1": x[10],
                                  "|D_mom-prod|/|prod|": x[11]} for x in d]

# ---------------- run-2 production path ----------------
if os.path.exists(LOG2):
    log2 = open(LOG2).read()
    out["run2"]["log_exists"] = True
    out["run2"]["thermal_dot"] = re.findall(r"Thermal matrix transpose dot-test max relative error=(\S+)", log2)
    out["run2"]["E4_max_dot"] = re.findall(r"E4: Helmholtz filter transpose dot-product test", log2)
    out["run2"]["prod_fgmres"] = re.findall(r"Production reduced discrete flow adjoint (\w+): FGMRES iterations=(\S+), true relative residual=(\S+)", log2)
    out["run2"]["prodsolv"] = re.findall(r"PRODSOLV \((\w+)\): sum\(J\^T\*1\)_P=(\S+) \|J\^T\*1_P\|L2=(\S+) rhsPmean=(\S+)", log2)
    out["run2"]["proddiag"] = re.findall(r"PRODDIAG \((\w+)\): jacobiUavg=(\S+) jacobiPavg=(\S+) jacobiMin=(\S+) jacobiMax=(\S+) jacobiNearZero=(\S+)", log2)
    out["run2"]["crash_markers"] = [m for m in ["FOAM FATAL", "FATAL ERROR", "sigFpe", "Segmentation", "FOAM aborting", "Floating point exception"] if m.lower() in log2.lower()]
else:
    out["run2"]["log_exists"] = False

with open(os.path.join(BASE, "..", "s3_gate_metrics.json"), "w") as f:
    json.dump(out, f, indent=1)
print(json.dumps(out, indent=1)[:6000])
