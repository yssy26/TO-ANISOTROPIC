#!/usr/bin/env python3
# Independent Post-Reviewer recomputation from the exported BFINAL-003 artifacts.
# Pure Python (no numpy), deliberately a different code path from the C++ probe.
import math, os

ART = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-003/cycle-1/artifacts"

def load(path):
    with open(path) as f:
        return [float(x) for x in f.read().split()]

def l2(a):
    return math.sqrt(sum(x*x for x in a))

def rel_l2(a, b):
    # relL2 = |a-b| / |b|  (b = reference/FD)
    nb = l2(b)
    d = [a[i]-b[i] for i in range(len(a))]
    return l2(d)/nb if nb > 0 else float('nan')

def cosim(a, b):
    na, nb = l2(a), l2(b)
    if na == 0 or nb == 0:
        return float('nan')
    return sum(a[i]*b[i] for i in range(len(a)))/(na*nb)

def norm_ratio(a, b):
    nb = l2(b)
    return l2(a)/nb if nb > 0 else float('nan')

def block_stats(a, b, idx):
    aa = [a[i] for i in idx]; bb = [b[i] for i in idx]
    return rel_l2(aa, bb), cosim(aa, bb), norm_ratio(aa, bb), l2(aa), l2(bb)

dphiA  = load(os.path.join(ART, "stageB5_gateB_dphi_A.mtx"))
dphiB  = load(os.path.join(ART, "stageB5_gateB_dphi_B_all_eps.mtx"))
dphiFD = load(os.path.join(ART, "stageB5_gateB_dphi_FD_all_eps.mtx"))
jv     = load(os.path.join(ART, "stageB5_gateA_Jv.mtx"))
fdAll  = load(os.path.join(ART, "stageB5_gateA_FD_all_eps.mtx"))
cellt  = load(os.path.join(ART, "stageB5_cell_type.mtx"))
facet  = load(os.path.join(ART, "stageB5_face_type.mtx"))

nFace = len(dphiA)
N     = len(cellt)
print(f"nFace={nFace}  Ncells={N}")
assert len(dphiB) == 5*nFace, f"dphiB size {len(dphiB)} != 5*{nFace}"
assert len(dphiFD) == 5*nFace, f"dphiFD size {len(dphiFD)} != 5*{nFace}"
assert len(jv) == 4*N, f"jv size {len(jv)} != 4*{N}"
assert len(fdAll) == 5*4*N, f"fdAll size {len(fdAll)} != 5*4*{N}"

nIntF = sum(1 for t in facet if t == 0)
nBndF = nFace - nIntF
print(f"internal faces={nIntF}  boundary faces={nBndF}")

# internal vs boundary face indices
int_idx = [i for i in range(nFace) if facet[i] == 0]
bnd_idx = [i for i in range(nFace) if facet[i] != 0]

print("\n=== G1: dphi oracle (per face) ===")
for ei, eps in enumerate([1e-3, 3e-4, 1e-4, 3e-5, 1e-5]):
    fd = dphiFD[ei*nFace:(ei+1)*nFace]
    b  = dphiB[ei*nFace:(ei+1)*nFace]
    # Candidate A is eps-independent
    ra, ca, na, laa, lfa = block_stats(dphiA, fd, list(range(nFace)))
    rb, cb, nb, lab, lfb = block_stats(b, fd, list(range(nFace)))
    ri_a, ci_a, _, _, _ = block_stats(dphiA, fd, int_idx)
    rbnd_a, cbnd_a, _, _, _ = block_stats(dphiA, fd, bnd_idx)
    print(f"eps={eps:g}: candA relL2={ra:.8e} cos={ca:.12f} normRatio={na:.12f} "
          f"(int={ri_a:.8e}, bnd={rbnd_a:.8e})")
    print(f"         candB relL2={rb:.8e} cos={cb:.12f} normRatio={nb:.12f}")

print("\n=== G2: full residual J*v vs FD (cell blocks) ===")
pint = [c for c in range(N) if cellt[c] == 0]
pbnd = [c for c in range(N) if cellt[c] == 1]
ptot = list(range(N))
# Layout: Jv is block-ordered U(3N) then P(N).
# FD export is interleaved per cell [Ux,Uy,Uz,P] -> de-interleave.
def deinterleave_fd(fd):
    fu = [0.0]*(3*N); fp = [0.0]*N
    for c in range(N):
        fu[3*c+0] = fd[4*c+0]
        fu[3*c+1] = fd[4*c+1]
        fu[3*c+2] = fd[4*c+2]
        fp[c] = fd[4*c+3]
    return fu, fp

ju = jv[:3*N]; jp = jv[3*N:]
for ei, eps in enumerate([1e-3, 3e-4, 1e-4, 3e-5, 1e-5]):
    fd = fdAll[ei*4*N:(ei+1)*4*N]
    fu, fp = deinterleave_fd(fd)
    ru, cu, _, _, _ = block_stats(ju, fu, list(range(3*N)))
    rpi, cpi, _, _, _ = block_stats(jp, fp, pint)
    rpb, cpb, _, _, _ = block_stats(jp, fp, pbnd)
    rpt, cpt, _, _, _ = block_stats(jp, fp, ptot)
    print(f"eps={eps:g}: U relL2={ru:.8e} cos={cu:.12f} | P-int={rpi:.8e}/{cpi:.12f} "
          f"P-bnd={rpb:.8e}/{cpb:.12f} P-tot={rpt:.8e}/{cpt:.12f}")

# === Localize the candA residual: A - B structure (eps=1e-3) ===
print("\n=== Localize candA residual: A - B (eps=1e-3) ===")
fd = dphiFD[0:nFace]
b  = dphiB[0:nFace]
diff_AB = [dphiA[i]-b[i] for i in range(nFace)]
diff_AFD = [dphiA[i]-fd[i] for i in range(nFace)]
diff_BFD = [b[i]-fd[i] for i in range(nFace)]
print(f"|A-B|L2 / |A|L2  = {l2(diff_AB)/l2(dphiA):.8e}")
print(f"|A-FD|L2/|FD|L2  = {rel_l2(dphiA, fd):.8e}")
print(f"|B-FD|L2/|FD|L2  = {rel_l2(b, fd):.8e}")
# per-face-type A-B magnitude
iab = l2([diff_AB[i] for i in int_idx])
bab = l2([diff_AB[i] for i in bnd_idx])
print(f"|A-B|L2 internal faces={iab:.8e}  boundary faces={bab:.8e}")
print(f"(so the A-vs-B gap is {'boundary-' if bab>iab else 'internal-'}dominated)")

# dominant faces in A-B gap
top = sorted(range(nFace), key=lambda i: abs(diff_AB[i]), reverse=True)[:5]
for i in top:
    typ = "internal" if facet[i]==0 else f"boundary patch {int(facet[i])-1}"
    print(f"  face {i} ({typ}): A={dphiA[i]:.8e} B={b[i]:.8e} diff={diff_AB[i]:.8e}")
