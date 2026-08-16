#!/usr/bin/env python3
"""Read-only localization of the J P-row vs NS.H rebuilt-phi FD discrepancy.

DIAGNOSTIC_ONLY (cycle-2, D4). Reads ONLY existing artifacts:
  explicitJT.mtx (J^T), stageB2_dUdirPhys.mtx, stageB2_rpDudU_rebuilt.mtx,
  stageB2_dphiDudU_face.mtx (NS.H rebuilt face dphi),
  stageB2_dphiDudU_FD.mtx  (solveDiscreteFlowAdjoint GateOracle face dphi FD,
                            same convention as J),
  stageB2_dphiDudU_J.mtx   (J's collectedDPhiJ face dphi),
  stageB2_Jdq_P_forward.mtx (J matrix-free P-row),
  constant/polyMesh/{owner,neighbour,boundary}  (boundary-cell identification)
"""
import sys

CASE = "/home/ys/b2_case_smoke/"
N = 33600
NU = 3 * N
TOT = 4 * N

def load_vec(name, count):
    with open(CASE + name, "r") as f:
        vals = [float(x) for x in f.read().split()]
    assert len(vals) == count, f"{name}: {len(vals)} != {count}"
    return vals

def l2(x):
    return sum(a*a for a in x) ** 0.5

def rel_cos(a, b, idxs, label):
    n = na = nb = dot = 0.0
    for i in idxs:
        d = a[i]-b[i]; n += d*d; na += a[i]*a[i]; nb += b[i]*b[i]; dot += a[i]*b[i]
    r = (n/nb)**0.5 if nb else float('nan')
    c = dot/((na*nb)**0.5) if na and nb else float('nan')
    print(f"  {label}: relL2={r:.6e} cos={c:.9f} |a|L2={na**0.5:.6e} |b|L2={nb**0.5:.6e}")

# ---- boundary cells from polyMesh ----
# owner file: header ends with a blank line, then "count\n(\nval\n..."
def read_label_list(fname):
    with open(fname) as f:
        raw = f.read()
    # first standalone integer token anywhere (the count), then values follow
    toks = raw.split()
    cnt = None
    for i, t in enumerate(toks):
        try:
            v = int(t)
        except ValueError:
            continue
        # the count is the first token that is a pure integer AND has at least
        # v following tokens (values are also ints); the header contains no
        # standalone integers except "7" (version) inside braces - guard by
        # requiring the token immediately after to be "("
        if i + 1 < len(toks) and toks[i + 1] == "(":
            cnt = v
            vals = [int(x) for x in toks[i + 2:i + 2 + v]]
            assert len(vals) == v, f"{fname}: parsed {len(vals)} != {v}"
            return vals
    raise SystemExit("no count found in " + fname)

owner = read_label_list(CASE + "constant/polyMesh/owner")
neigh = read_label_list(CASE + "constant/polyMesh/neighbour")
nFaces = len(owner)
nInt = len(neigh)          # 96860
print(f"nFaces={nFaces} nInternal={nInt} nCells={N}")

# boundary face ranges from the boundary dict (parsed manually below)
bnd = [
    ("inlet", 84, 96860), ("outlet", 84, 96944),
    ("hotInlet", 196, 97028), ("hotOutlet", 196, 97224),
    ("solidEndWalls", 280, 97420), ("bottomWall", 1120, 97700),
    ("topWall", 1120, 98820), ("sideWalls", 4800, 99940),
]
bndCells = set()
for name, nf, start in bnd:
    for fi in range(start, start+nf):
        bndCells.add(owner[fi])
print(f"boundary cells: {len(bndCells)}")

# ---- J*v P-row from explicitJT (scatter) ----
dU = load_vec("stageB2_dUdirPhys.mtx", NU)
v = [0.0]*TOT
for i in range(NU): v[i] = dU[i]
out = [0.0]*TOT
with open(CASE + "explicitJT.mtx") as f:
    f.readline(); f.readline()
    for line in f:
        p = line.split()
        if len(p) < 3: continue
        r = int(p[0])-1; c = int(p[1])-1; val = float(p[2])
        if val != 0.0: out[c] += val*v[r]
outP = out[NU:]
rp = load_vec("stageB2_rpDudU_rebuilt.mtx", N)
jq = load_vec("stageB2_Jdq_P_forward.mtx", N)

internal = [i for i in range(N) if i not in bndCells]
print("\nP-row J*v vs NS.H rebuilt-phi FD (stageB2_rpDudU_rebuilt):")
rel_cos(outP, rp, list(range(N)), "ALL cells")
rel_cos(outP, rp, internal, "internal cells only")
rel_cos(outP, rp, sorted(bndCells), "boundary cells only")
print("\nP-row J*v vs J matrix-free P-row (Jdq_P_forward):")
rel_cos(outP, jq, list(range(N)), "ALL cells")

# ---- face-level dphi comparisons ----
nF = nInt  # internal faces 96860
dF = load_vec("stageB2_dphiDudU_FD.mtx", nF)     # GateOracle FD (J convention)
dJ = load_vec("stageB2_dphiDudU_J.mtx", nF)      # J's collectedDPhiJ
dR = load_vec("stageB2_dphiDudU_face.mtx", nF)   # NS.H rebuilt face dphi
print("\nface-level dphi/dU (internal faces):")
rel_cos(dJ, dF, list(range(nF)), "J-collected vs GateOracle-FD (same conv)")
rel_cos(dJ, dR, list(range(nF)), "J-collected vs NS.H-rebuilt face dphi")
rel_cos(dF, dR, list(range(nF)), "GateOracle-FD vs NS.H-rebuilt face dphi")

# ---- difference concentration: per-region contribution to |Jv_P - rp| ----
diff = [outP[i]-rp[i] for i in range(N)]
print(f"\n|diff|L2 ALL={l2(diff):.6e}  internal={l2([diff[i] for i in internal]):.6e}  "
      f"boundary={l2([diff[i] for i in bndCells]):.6e}")
