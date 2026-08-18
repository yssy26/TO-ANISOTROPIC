#!/usr/bin/env python3
"""Independent recomputation of the BFINAL-002 Gate-A / Gate-B metrics from the
raw stageB5_*.mtx artifacts (Post-Reviewer / evidence tool).

DIAGNOSTIC_ONLY (BFINAL-002, S2/S3).  Reads ONLY the exported artifacts:
  stageB5_gateA_Jv.mtx                  (4N: J*v, U part 3N then P part N)
  stageB5_gateA_FD_all_eps.mtx          (5 blocks of 4N: FD_U(3N)+FD_P(N))
  stageB5_dUdir.mtx                     (3N direction)
  stageB5_gateB_dphi_A.mtx              (nInt+nBnd: Candidate A dphi)
  stageB5_gateB_dphi_B_all_eps.mtx      (5 blocks of nInt+nBnd)
  stageB5_gateB_dphi_FD_all_eps.mtx     (5 blocks of nInt+nBnd)
  stageB5_cell_type.mtx                 (N: 0 internal cell, 1 boundary cell)
  stageB5_face_type.mtx                 (nInt+nBnd: 0 internal, 1..8 patch+1)
Layout conventions are documented in evidence/agent-group/BFINAL-002/
stageB5_probe_design.md.  relL2 = |a-b|_2/|b|_2 (b = FD reference), cosine =
(a.b)/(|a||b|), exactly the conventions of BFINAL-001 crosscheck_bf2.py.

Usage: python3 recompute_stageB5.py <caseDir>
"""
import sys

CASE = sys.argv[1] if len(sys.argv) > 1 else "."
EPS = [1e-3, 3e-4, 1e-4, 3e-5, 1e-5]

def load_vec(name, count=None):
    with open(CASE + "/" + name) as f:
        vals = [float(x) for x in f.read().split()]
    if count is not None and len(vals) != count:
        raise SystemExit(f"{name}: expected {count} values, got {len(vals)}")
    return vals

def l2(x, idx=None):
    s = 0.0
    if idx is None:
        for v in x: s += v * v
    else:
        for i in idx: s += x[i] * x[i]
    return s ** 0.5

def metrics(a, b, idx, label):
    na = nb = nd = dot = 0.0
    mx = 0.0; arg = -1
    for i in idx:
        av, bv = a[i], b[i]
        d = av - bv
        na += av*av; nb += bv*bv; nd += d*d; dot += av*bv
        if abs(d) > mx: mx = abs(d); arg = i
    rel = (nd/nb)**0.5 if nb > 0 else float("nan")
    cos = dot/((na*nb)**0.5) if (na > 0 and nb > 0) else float("nan")
    print(f"  {label}: |a|L2={na**0.5:.6e} |FD|L2={nb**0.5:.6e} "
          f"|err|L2={nd**0.5:.6e} relL2={rel:.6e} cos={cos:.9f} "
          f"maxDiff={mx:.6e} @idx={arg}")
    return rel, cos

# ---- load ----
N = len(load_vec("stageB5_cell_type.mtx"))
NU = 3 * N
TOT = 4 * N
cell_type = load_vec("stageB5_cell_type.mtx", N)
face_type = load_vec("stageB5_face_type.mtx")
nF = len(face_type)
nInt = sum(1 for t in face_type if t == 0)
print(f"N={N} nFaces(internal+boundary)={nF} nInternal={nInt} "
      f"nBoundary={nF-nInt}")

jv = load_vec("stageB5_gateA_Jv.mtx", TOT)
jvU = jv[:NU]
jvP = jv[NU:]
fdAll = load_vec("stageB5_gateA_FD_all_eps.mtx", 5*TOT)
dphiA = load_vec("stageB5_gateB_dphi_A.mtx", nF)
dphiBAll = load_vec("stageB5_gateB_dphi_B_all_eps.mtx", 5*nF)
dphiFDAll = load_vec("stageB5_gateB_dphi_FD_all_eps.mtx", 5*nF)

int_cells = [c for c in range(N) if cell_type[c] == 0]
bnd_cells = [c for c in range(N) if cell_type[c] == 1]
int_faces = [f for f in range(nF) if face_type[f] == 0]
bnd_faces = [f for f in range(nF) if face_type[f] != 0]
all_faces = list(range(nF))

print("\n=== Gate A (J*v vs production-consistent FD), per eps ===")
print("eps        | momentum relL2/cos | P-int relL2/cos | P-bnd relL2/cos | P-tot relL2/cos")
for e, eps in enumerate(EPS):
    blk = fdAll[e*TOT:(e+1)*TOT]
    # NOTE: the probe writes the Gate-A FD per cell INTERLEAVED as
    #   celli: Ux Uy Uz P   (4N values per eps block)
    fdU = [0.0]*NU
    fdP = [0.0]*N
    for c in range(N):
        for k in range(3):
            fdU[3*c + k] = blk[4*c + k]
        fdP[c] = blk[4*c + 3]
    rM = metrics(jvU, fdU, list(range(NU)), f"eps={eps} momentum")
    rI = metrics(jvP, fdP, int_cells, f"eps={eps} P-internal")
    rB = metrics(jvP, fdP, bnd_cells, f"eps={eps} P-boundary")
    rT = metrics(jvP, fdP, list(range(N)), f"eps={eps} P-total")
    print(f"  {eps:.1e}   {rM[0]:.3e}/{rM[1]:.6f}   {rI[0]:.3e}/{rI[1]:.6f}"
          f"   {rB[0]:.3e}/{rB[1]:.6f}   {rT[0]:.3e}/{rT[1]:.6f}")

print("\n=== Gate B (dphi candidates vs FD), per eps ===")
for e, eps in enumerate(EPS):
    dphiFD = dphiFDAll[e*nF:(e+1)*nF]
    dphiB = dphiBAll[e*nF:(e+1)*nF]
    print(f"eps={eps}: Candidate A (J, unrelaxed):")
    metrics(dphiA, dphiFD, int_faces, "  A internal faces")
    metrics(dphiA, dphiFD, bnd_faces, "  A boundary faces")
    metrics(dphiA, dphiFD, all_faces, "  A all faces")
    print(f"eps={eps}: Candidate B (production relaxed):")
    metrics(dphiB, dphiFD, int_faces, "  B internal faces")
    metrics(dphiB, dphiFD, bnd_faces, "  B boundary faces")
    metrics(dphiB, dphiFD, all_faces, "  B all faces")

# face-level dphi norms for context
print("\nface-dphi context (eps=1e-3):")
dphiFD0 = dphiFDAll[0:nF]
dphiB0 = dphiBAll[0:nF]
print(f"  |dphi_FD|L2={l2(dphiFD0):.6e} int={l2(dphiFD0, int_faces):.6e} "
      f"bnd={l2(dphiFD0, bnd_faces):.6e}")
print(f"  |dphi_A|L2={l2(dphiA):.6e} int={l2(dphiA, int_faces):.6e} "
      f"bnd={l2(dphiA, bnd_faces):.6e}")
print(f"  |dphi_B|L2={l2(dphiB0):.6e} int={l2(dphiB0, int_faces):.6e} "
      f"bnd={l2(dphiB0, bnd_faces):.6e}")
