#!/usr/bin/env python3
"""BFINAL-002 S4 (Post-Reviewer) — INDEPENDENT recomputation of Gate A / Gate B
metrics from the raw stageB5_*.mtx artifacts, plus two cross-checks the probe
does NOT print:

  (X1) per-patch boundary-face analysis: which patches carry |dphi_FD| (must be
       ~0 at the 7 constrainHbyA-pinned fixedValue-U patches, nonzero only at
       the zeroGradient-U outlet), and per-patch Candidate A vs FD relL2.
  (X2) divergence consistency: reconstruct div(dphi_FD) per cell from the mesh
       (owner/neighbour + Sf from faces/points) and compare to the Gate-A
       per-cell P-row FD from the artifacts — i.e., the two oracles measure the
       SAME derivative object.

Implementation is deliberately independent of recompute_stageB5.py (numpy,
different metric code path).  Reads ONLY exported artifacts + the case mesh.

Usage: python3 recompute_s4_independent.py <caseDir>
"""
import sys
import numpy as np

CASE = sys.argv[1] if len(sys.argv) > 1 else "."
EPS = [1e-3, 3e-4, 1e-4, 3e-5, 1e-5]
PATCH_NAMES = ["inlet", "outlet", "hotInlet", "hotOutlet",
               "solidEndWalls", "bottomWall", "topWall", "sideWalls"]


def load(name, count=None):
    v = np.fromfile(CASE + "/" + name, sep=" ", dtype=np.float64)
    if count is not None:
        assert v.size == count, f"{name}: expected {count} got {v.size}"
    return v


def metrics(a, b, idx, label):
    a = a[idx]
    b = b[idx]
    d = a - b
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    nd = np.linalg.norm(d)
    rel = nd / nb if nb > 0 else float("nan")
    cos = float(a @ b / (na * nb)) if (na > 0 and nb > 0) else float("nan")
    k = int(np.argmax(np.abs(d)))
    print(f"  {label}: |a|L2={na:.6e} |FD|L2={nb:.6e} |err|L2={nd:.6e} "
          f"relL2={rel:.6e} cos={cos:.9f} maxDiff={np.abs(d[k]):.6e} @idx={idx[k]}")
    return rel, cos


N = 33600
NU = 3 * N
TOT = 4 * N
cell_type = load("stageB5_cell_type.mtx", N)
face_type = load("stageB5_face_type.mtx", None)
nF = face_type.size
nInt = int((face_type == 0).sum())
nBnd = nF - nInt
assert nF == 104740 and nInt == 96860 and nBnd == 7880, (nF, nInt, nBnd)
print(f"N={N} nFaces={nF} nInternal={nInt} nBoundary={nBnd}")

jv = load("stageB5_gateA_Jv.mtx", TOT)
jvU = jv[:NU]
jvP = jv[NU:]
fdAll = load("stageB5_gateA_FD_all_eps.mtx", 5 * TOT)
dphiA = load("stageB5_gateB_dphi_A.mtx", nF)
dphiBAll = load("stageB5_gateB_dphi_B_all_eps.mtx", 5 * nF)
dphiFDAll = load("stageB5_gateB_dphi_FD_all_eps.mtx", 5 * nF)

int_cells = np.where(cell_type == 0)[0]
bnd_cells = np.where(cell_type == 1)[0]
int_faces = np.where(face_type == 0)[0]
bnd_faces = np.where(face_type != 0)[0]
all_faces = np.arange(nF)

# ---------------- Gate A ----------------
print("\n=== [independent] Gate A: J*v vs production-consistent FD, per eps ===")
print("eps        | momentum relL2/cos | P-int relL2/cos | P-bnd relL2/cos | P-tot relL2/cos")
gateA = {}
for e, eps in enumerate(EPS):
    blk = fdAll[e * TOT:(e + 1) * TOT]
    # artifact layout: per cell interleaved Ux Uy Uz P
    fdU = blk.reshape(N, 4)[:, :3].ravel()
    fdP = blk.reshape(N, 4)[:, 3]
    rM = metrics(jvU, fdU, np.arange(NU), f"eps={eps} momentum")
    rI = metrics(jvP, fdP, int_cells, f"eps={eps} P-internal")
    rB = metrics(jvP, fdP, bnd_cells, f"eps={eps} P-boundary")
    rT = metrics(jvP, fdP, np.arange(N), f"eps={eps} P-total")
    print(f"  {eps:.1e}   {rM[0]:.3e}/{rM[1]:.6f}   {rI[0]:.3e}/{rI[1]:.6f}"
          f"   {rB[0]:.3e}/{rB[1]:.6f}   {rT[0]:.3e}/{rT[1]:.6f}")
    gateA[eps] = dict(mom=rM, pint=rI, pbnd=rB, ptot=rT, fdP=fdP)

# ---------------- Gate B ----------------
print("\n=== [independent] Gate B: dphi candidates vs FD, per eps ===")
for e, eps in enumerate(EPS):
    dphiFD = dphiFDAll[e * nF:(e + 1) * nF]
    dphiB = dphiBAll[e * nF:(e + 1) * nF]
    print(f"eps={eps}: Candidate A (J, unrelaxed):")
    metrics(dphiA, dphiFD, int_faces, "  A internal faces")
    metrics(dphiA, dphiFD, bnd_faces, "  A boundary faces")
    rA = metrics(dphiA, dphiFD, all_faces, "  A all faces")
    print(f"   A normRatio={np.linalg.norm(dphiA)/np.linalg.norm(dphiFD):.9f}")
    print(f"eps={eps}: Candidate B (production relaxed):")
    metrics(dphiB, dphiFD, int_faces, "  B internal faces")
    metrics(dphiB, dphiFD, bnd_faces, "  B boundary faces")
    rB = metrics(dphiB, dphiFD, all_faces, "  B all faces")
    print(f"   B normRatio={np.linalg.norm(dphiB)/np.linalg.norm(dphiFD):.9f}")

# ---------------- X1: per-patch boundary analysis ----------------
print("\n=== [X1] per-patch boundary dphi (eps=1e-3 block) ===")
dphiFD0 = dphiFDAll[0:nF]
print(f"{'patch':12s} {'nFaces':>6s} {'|dphi_FD|L2':>12s} {'|dphi_A|L2':>12s} "
      f"{'A-vs-FD relL2':>13s} {'A-vs-FD cos':>12s}")
for p in range(1, 9):
    f = np.where(face_type == p)[0]
    a = dphiA[f]
    b = dphiFD0[f]
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    rel = np.linalg.norm(a - b) / nb if nb > 0 else float("nan")
    cos = float(a @ b / (na * nb)) if (na > 0 and nb > 0) else float("nan")
    print(f"{PATCH_NAMES[p-1]:12s} {f.size:6d} {nb:12.4e} {na:12.4e} "
          f"{rel:13.4e} {cos:12.6f}")
# check the production-pinning expectation: |dphi_FD| ~ 0 on the 7 pinned patches
pinned = np.isin(face_type, [1, 3, 4, 5, 6, 7, 8])
outlet = face_type == 2
print(f"|dphi_FD| on 7 pinned patches: {np.linalg.norm(dphiFD0[pinned]):.4e} "
      f"(expect ~1e-11, constrainHbyA pinning)")
print(f"|dphi_FD| on outlet patch   : {np.linalg.norm(dphiFD0[outlet]):.4e} "
      f"(the only assignable-U patch)")

# ---------------- X2: divergence consistency ----------------
print("\n=== [X2] div(dphi_FD) per cell vs Gate-A P-row FD per cell ===")
def parse_mesh():
    # points
    with open(f"{CASE}/constant/polyMesh/points") as f:
        txt = f.read()
    head = txt.split("(\n", 1)[1]
    body = head.split("\n)\n", 1)[0]
    pts = np.array([[float(x) for x in line.strip(" ()").split()]
                    for line in body.splitlines() if line.strip()])
    # faces
    with open(f"{CASE}/constant/polyMesh/faces") as f:
        txt = f.read()
    head = txt.split("(\n", 1)[1]
    body = head.split("\n)\n", 1)[0]
    faces = []
    for line in body.splitlines():
        line = line.strip()
        if not line:
            continue
        # OpenFOAM compact face format: 4(1 2 3 4)
        if "(" in line:
            rest = line.split("(", 1)[1]
            idx = [int(x) for x in rest.strip(") ").split()]
        else:
            idx = [int(x) for x in line.split()]
        faces.append(idx)
    # owner / neighbour
    def label_list(path):
        with open(path) as f:
            txt = f.read()
        body = txt.split("(\n", 1)[1].split("\n)\n", 1)[0]
        return np.array([int(x) for x in body.split()], dtype=np.int64)
    owner = label_list(f"{CASE}/constant/polyMesh/owner")
    nb = label_list(f"{CASE}/constant/polyMesh/neighbour")
    return pts, faces, owner, nb

pts, faces, owner, nb = parse_mesh()
nPts, nFacesM, nCells = pts.shape[0], len(faces), 33600
assert nFacesM == 104740 and nPts == 37665
# face area vectors (planar polygons): Sf = 0.5 * sum_i (p_i x p_{i+1})
Sf = np.zeros((nFacesM, 3))
Cf = np.zeros((nFacesM, 3))
for fi, fv in enumerate(faces):
    P = pts[fv]
    Cf[fi] = P.mean(axis=0)
    cross = np.cross(P, np.roll(P, -1, axis=0)).sum(axis=0)
    Sf[fi] = 0.5 * cross
# orientation sanity check: internal faces owner->neighbour, boundary outward
# cell volumes via divergence theorem: V = (1/3) sum_f Sf_f . Cf_f (outward)
V = np.zeros(nCells)
cell_faces = [[] for _ in range(nCells)]
for fi in range(96860):
    o, n = int(owner[fi]), int(nb[fi])
    cell_faces[o].append((fi, +1))
    cell_faces[n].append((fi, -1))
for fi in range(96860, nFacesM):
    o = int(owner[fi])
    cell_faces[o].append((fi, +1))
for c in range(nCells):
    V[c] = (1.0 / 3.0) * sum(s * float(Sf[fi] @ Cf[fi]) for fi, s in cell_faces[c])
    assert V[c] > 0, c
# cell centres via the face-based formula (exact for linear tets/quads):
#   Cc = (1/(4V)) sum_f (Sf_f . Cf_f)(Cf_f + Cf_f)  -- for hexahedra use
#   the standard (1/(2V)) sum (Sf.Cf) Cf; adequate for orientation checks
Cc = np.zeros((nCells, 3))
for c in range(nCells):
    num = np.zeros(3)
    for fi, s in cell_faces[c]:
        num += s * float(Sf[fi] @ Cf[fi]) * Cf[fi]
    Cc[c] = num / (2.0 * V[c])
# orientation: internal face Sf should point owner->neighbour
bad = 0
for fi in range(96860):
    o, n = int(owner[fi]), int(nb[fi])
    if np.dot(Sf[fi], Cc[n] - Cc[o]) < 0:
        bad += 1
print(f"Sf orientation check: {bad}/96860 internal faces point neighbour->owner")
# boundary outward
badb = 0
for fi in range(96860, nFacesM):
    o = int(owner[fi])
    if np.dot(Sf[fi], Cf[fi] - Cc[o]) < 0:
        badb += 1
print(f"Sf orientation check: {badb}/{nFacesM-96860} boundary faces point inward")

# reconstruct div(dphi_FD) per cell (eps=1e-3 block) -- raw face sum (the
# probe's P-row residual convention is the raw sum, NOT (1/V) sum; the 1/V
# factor observed = 8.09e9 = 1/V_cell avg)
dph = dphiFD0
divphi = np.zeros(nCells)   # raw sum per cell
for fi in range(96860):
    o, n = int(owner[fi]), int(nb[fi])
    divphi[o] += dph[fi]
    divphi[n] -= dph[fi]
for fi in range(96860, nFacesM):
    divphi[int(owner[fi])] += dph[fi]

fdP1 = gateA[1e-3]["fdP"]  # Gate-A P-row FD (production-consistent) per cell
print(f"|sum_f dphi_FD|L2(all)      = {np.linalg.norm(divphi):.6e}")
print(f"|GateA P-row FD|L2(all)    = {np.linalg.norm(fdP1):.6e}")
rel_all = np.linalg.norm(divphi - fdP1) / np.linalg.norm(fdP1)
cos_all = float(divphi @ fdP1 / (np.linalg.norm(divphi) * np.linalg.norm(fdP1)))
print(f"all cells:      relL2={rel_all:.3e} cos={cos_all:.9f}")
rel_i = np.linalg.norm(divphi[int_cells] - fdP1[int_cells]) / np.linalg.norm(fdP1[int_cells])
rel_b = np.linalg.norm(divphi[bnd_cells] - fdP1[bnd_cells]) / np.linalg.norm(fdP1[bnd_cells])
print(f"internal cells: relL2={rel_i:.3e}")
print(f"boundary cells: relL2={rel_b:.3e}")

# also reconstruct div of Candidate B and Candidate A for reference
for name, vec in (("A", dphiA), ("B", dphiBAll[0:nF])):
    dv = np.zeros(nCells)
    for fi in range(96860):
        o, n = int(owner[fi]), int(nb[fi])
        dv[o] += vec[fi]
        dv[n] -= vec[fi]
    for fi in range(96860, nFacesM):
        dv[int(owner[fi])] += vec[fi]
    dv /= V
    print(f"div(Candidate {name}): |.|L2={np.linalg.norm(dv):.6e}")

# ---------------- X3: relaxation algebra on face vectors ----------------
print("\n=== [X3] production-tangent algebra check ===")
# theory (audit §8, verified from fvMatrix::relax: D/=alpha, S+=(D-D0)*psi):
#   rAU_rel*H_rel = 0.4*rAU_u*H_u + 0.6*U/V
#   =>  dphi_prod = flux(0.4*rAU_u*dH_u) + flux(0.6*dU/V)
#   =>  FD = 0.4*A + X, X = flux(0.6*dU/V)
a, b = dphiA, dphiFD0
res = b - 0.4 * a
print(f"|FD|={np.linalg.norm(b):.6e} |A|={np.linalg.norm(a):.6e} "
      f"|FD-0.4A|={np.linalg.norm(res):.6e}")
print(f"|A|/|FD|={np.linalg.norm(a)/np.linalg.norm(b):.6f} "
      f"(0.591, not 0.4: the 0.6*dU/V relax-source term is NOT small)")

# numerical test of X = flux(0.6*dU/V): build the face flux of the per-cell
# field 0.6*dUdir/V with distance-based linear weights on internal faces
dU = load("stageB5_dUdir.mtx", 3 * N)
dU3 = dU.reshape(N, 3)
w = np.zeros(96860)
for fi in range(96860):
    o, n = int(owner[fi]), int(nb[fi])
    do = np.linalg.norm(Cf[fi] - Cc[o])
    dn = np.linalg.norm(Cc[n] - Cf[fi])
    w[fi] = dn / (do + dn)          # weight of owner cell value
X = np.zeros(nF)
for fi in range(96860):
    o, n = int(owner[fi]), int(nb[fi])
    vo = 0.6 * dU3[o]
    vn = 0.6 * dU3[n]
    X[fi] = float(Sf[fi] @ (w[fi] * 0.6 * dU3[o] + (1.0 - w[fi]) * 0.6 * dU3[n]))
# boundary: at the 7 pinned patches dphi_FD=0 (constrainHbyA replaces HbyA_b
# with U_b, so the 0.6*dU/V term is absent); at the outlet (extrapolated) it
# enters with the cell value (zeroGradient => U_b = U_cell)
for fi in range(96860, nFacesM):
    o = int(owner[fi])
    if face_type[fi] == 2.0:        # outlet, the assignable-U patch
        X[fi] = float(Sf[fi] @ (0.6 * dU3[o]))
    else:
        X[fi] = 0.0
relX = np.linalg.norm(res - X) / np.linalg.norm(res)
cosX = float(res @ X / (np.linalg.norm(res) * np.linalg.norm(X)))
print(f"|X_theory=flux(0.6 dU/V)|L2 = {np.linalg.norm(X):.6e}")
print(f"X vs (FD-0.4A): relL2={relX:.3e} cos={cosX:.6f}")
ri = np.linalg.norm((res - X)[int_faces]) / np.linalg.norm(res[int_faces])
print(f"  internal faces only: relL2={ri:.3e}")

print("\n[done]")
