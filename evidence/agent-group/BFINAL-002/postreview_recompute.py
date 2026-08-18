#!/usr/bin/env python3
"""BFINAL-002 Post-Reviewer INDEPENDENT recomputation (this invocation).

Deliberately separate code path from recompute_stageB5.py and
recompute_s4_independent.py.  Loads ONLY the exported stageB5_*.mtx artifacts
and re-derives, with fresh code, the load-bearing Gate A and Gate B metrics.

Checks performed:
  (1) Gate A per-block: J*v vs production-consistent FD, momentum/P-int/P-bnd/P-tot.
  (2) Gate B per eps: Candidate A (J, unrelaxed rAUAdj) and Candidate B
      (production relaxed mobility) vs true FD dphi, internal/boundary/all faces.
  (3) FD eps-independence (linearity): the reconstructed phi(U) is affine in U
      (momentum matrix is linear in U with fixed frozen phi/coefficients), so the
      central-difference dphi_FD must be eps-independent to roundoff.
  (4) The "0.591 not 0.4" relaxation fingerprint at face level: |FD - 0.4*A|/|FD|
      is NOT ~0 => the relax-source term (0.6*dU/V) is present in production.
"""
import numpy as np

ART = "artifacts_s3"
N = 33600
NU = 3 * N
TOT = 4 * N
EPS = [1e-3, 3e-4, 1e-4, 3e-5, 1e-5]


def load(name, expect=None):
    v = np.fromfile(f"{ART}/{name}", sep=" ", dtype=np.float64)
    if expect is not None:
        assert v.size == expect, f"{name}: size {v.size} != {expect}"
    return v


cell_type = load("stageB5_cell_type.mtx", N).astype(int)
face_type = load("stageB5_face_type.mtx").astype(int)
nF = face_type.size
nInt = int((face_type == 0).sum())
nBnd = nF - nInt
assert nF == 104740 and nInt == 96860 and nBnd == 7880, (nF, nInt, nBnd)

jv = load("stageB5_gateA_Jv.mtx", TOT)
jvU = jv[:NU]
jvP = jv[NU:]
fdAll = load("stageB5_gateA_FD_all_eps.mtx", 5 * TOT)
dphiA = load("stageB5_gateB_dphi_A.mtx", nF)
dphiB = load("stageB5_gateB_dphi_B_all_eps.mtx", 5 * nF)
dphiFD = load("stageB5_gateB_dphi_FD_all_eps.mtx", 5 * nF)

int_cells = np.where(cell_type == 0)[0]
bnd_cells = np.where(cell_type == 1)[0]
int_faces = np.where(face_type == 0)[0]
bnd_faces = np.where(face_type != 0)[0]
all_faces = np.arange(nF)


def m(a, b):
    """return (relL2, cos, |a|, |b|, |err|)"""
    d = a - b
    na = float(np.sqrt(a @ a))
    nb = float(np.sqrt(b @ b))
    nd = float(np.sqrt(d @ d))
    rel = nd / nb if nb > 0 else float("nan")
    cos = float(a @ b) / (na * nb) if (na > 0 and nb > 0) else float("nan")
    return rel, cos, na, nb, nd


print(f"N={N} nFaces={nF} nInternal={nInt} nBoundary={nBnd}")
print(f"internal cells={int_cells.size} boundary cells={bnd_cells.size}")

print("\n=== [post-reviewer] Gate A: J*v vs production FD ===")
for e, eps in enumerate(EPS):
    blk = fdAll[e * TOT:(e + 1) * TOT].reshape(N, 4)
    fdU = blk[:, :3].ravel()
    fdP = blk[:, 3]
    rM = m(jvU, fdU)
    rI = m(jvP[None, :], fdP[None, :]) if False else m(jvP[int_cells], fdP[int_cells])
    rB = m(jvP[bnd_cells], fdP[bnd_cells])
    rT = m(jvP, fdP)
    print(f"eps={eps:.0e}: mom relL2={rM[0]:.6e} cos={rM[1]:.9f} | "
          f"P-int relL2={rI[0]:.6e} cos={rI[1]:.9f} | "
          f"P-bnd relL2={rB[0]:.6e} cos={rB[1]:.9f} | "
          f"P-tot relL2={rT[0]:.6e} cos={rT[1]:.9f}")

print("\n=== [post-reviewer] Gate B: candidates vs FD ===")
for e, eps in enumerate(EPS):
    f = dphiFD[e * nF:(e + 1) * nF]
    b = dphiB[e * nF:(e + 1) * nF]
    a_int = m(dphiA[int_faces], f[int_faces])
    a_bnd = m(dphiA[bnd_faces], f[bnd_faces])
    a_all = m(dphiA, f)
    b_int = m(b[int_faces], f[int_faces])
    b_bnd = m(b[bnd_faces], f[bnd_faces])
    b_all = m(b, f)
    print(f"eps={eps:.0e}: A relL2={a_all[0]:.6e} cos={a_all[1]:.9f} "
          f"normRatio={a_all[2]/a_all[3]:.9f} (int {a_int[0]:.6e}, bnd {a_bnd[0]:.6e})")
    print(f"          B relL2={b_all[0]:.6e} cos={b_all[1]:.9f} "
          f"normRatio={b_all[2]/b_all[3]:.9f} (int {b_int[0]:.6e}, bnd {b_bnd[0]:.6e})")

print("\n=== [post-reviewer] FD eps-independence (affine phi(U) => central FD exact) ===")
ref = dphiFD[0:nF]
for e in range(1, 5):
    cur = dphiFD[e * nF:(e + 1) * nF]
    rel, cos, _, _, _ = m(ref, cur)
    print(f"  dphiFD[eps={EPS[e]:.0e}] vs dphiFD[eps=1e-3]: relL2={rel:.3e} cos={cos:.9f}")

print("\n=== [post-reviewer] relaxation fingerprint: |FD - 0.4*A| / |FD| ===")
f = dphiFD[0:nF]
res = f - 0.4 * dphiA
print(f"  |A|/|FD| = {np.linalg.norm(dphiA)/np.linalg.norm(f):.9f}  (0.591 => relax-source term present)")
print(f"  |FD - 0.4*A| / |FD| = {np.linalg.norm(res)/np.linalg.norm(f):.9f}  "
      f"(would be ~0 if J were a pure 0.4-rescale of production)")
print(f"  |FD - 1.0*A| / |FD| = {np.linalg.norm(f - dphiA)/np.linalg.norm(f):.9f}  "
      f"(this is Candidate A relL2, the actual J mismatch)")

print("\n[done]")
