#!/usr/bin/env python3
"""BFINAL-023 cycle-1 T8: relax-source boundary post-piece (preregistered
in PREREGISTRATION_T7_ADDENDUM.md T8 section).

H8a value gate, H8b anchors (du+dp+press configuration primary; full
3-channel also reported).
"""
import time, json
import numpy as np

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
SB = "/home/ys/dsH/b22_gauge/stageB2"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-023/cycle-1"

N = 33600
ALPHA_REL = 0.4
NU = 5.19009e-05

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-023 cycle-1 T9: A/B-form with corrected map ===")

pts = []
intxt = False
for line in open(POLY + "/points"):
    s = line.strip()
    if not intxt:
        if s == "(": intxt = True
        continue
    if s == ")": break
    pts.append([float(t) for t in s.strip("()").split()])
pts = np.array(pts)
faces = []
intxt = False
for line in open(POLY + "/faces"):
    s = line.strip()
    if not intxt:
        if s == "(": intxt = True; continue
    if s == ")": break
    if "(" in s and s.endswith(")"):
        faces.append([int(v) for v in s[s.index("(")+1:-1].split()])
def rl(p):
    t = open(p).read(); i = t.index("\n("); j = t.index("\n)", i)
    return np.array([int(x) for x in t[i+2:j].split()], dtype=np.int64)
oa = rl(POLY + "/owner"); na = rl(POLY + "/neighbour")
nIF = len(na); nF = len(faces)
Sf = np.zeros((nF, 3))
for fi in range(nF):
    fv = pts[faces[fi]]
    Sf[fi] = 0.5*np.cross(fv, np.roll(fv, -1, axis=0)).sum(axis=0)
magSf = np.linalg.norm(Sf, axis=1)
cv = {}
for fi in range(nF):
    cs = (oa[fi], na[fi]) if fi < nIF else (oa[fi],)
    for c in cs: cv.setdefault(c, set()).update(faces[fi])
C = np.zeros((N, 3))
for c, vs in cv.items(): C[c] = pts[sorted(vs)].mean(axis=0)
owner = oa[:nIF]; nei = na[:nIF]
dvec = C[nei] - C[owner]
magSfI = magSf[:nIF]
dc = magSfI/np.einsum("ij,ij->i", Sf[:nIF], dvec)
V = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")[:, 0]
w = np.loadtxt(MESH + "/b16mesh_weights.mtx")
bnd = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
bc = bnd[:, 0].astype(np.int64); bSf = bnd[:, 1:4]; bmag = bnd[:, 4]; bdel = bnd[:, 5]
Ufix = bnd[:, 7] > 0.5
nB = len(bc); bidx = np.where(~Ufix)[0]
sfI = Sf[:nIF]
import re
bfile = open(POLY + "/boundary").read()
patch_names = re.findall(r"^\s{4}(\w+)\s*$", bfile, re.M)
patch_sizes = [int(m) for m in re.findall(r"nFaces\s+(\d+);", bfile)]
patch_of_bface = np.concatenate(
    [np.full(s, i, dtype=np.int64) for i, s in enumerate(patch_sizes)])
OUTLET_IDX = patch_names.index("outlet")

def build_basis(phi_, phiB_, alpha_, nuEff_):
    qp = (phi_ >= 0).astype(float); qn = 1 - qp
    divp = np.zeros(N)
    np.add.at(divp, owner, phi_); np.add.at(divp, nei, -phi_)
    np.add.at(divp, bc, phiB_)
    gam = w*nuEff_[owner] + (1 - w)*nuEff_[nei]
    kfv = gam*magSfI*dc
    lower = -qp*phi_ - kfv; upper = qn*phi_ - kfv
    diag = np.zeros(N)
    np.add.at(diag, owner, -lower); np.add.at(diag, nei, -upper)
    diag -= divp; diag += alpha_*V
    sumOff = np.zeros(N)
    np.add.at(sumOff, owner, np.abs(upper)); np.add.at(sumOff, nei, np.abs(lower))
    ab = nuEff_[bc]*bmag*bdel
    ic = np.where(Ufix, -ab, phiB_)
    Db = diag.copy(); np.add.at(Db, bc, np.abs(ic))
    Drel = np.maximum(np.abs(Db), sumOff)/ALPHA_REL
    # relax-internal diagonal = Drel - cmptMin(ic): Ufixed: +ab ; outlet: -phi_b
    bnd_post = np.zeros(N)
    np.add.at(bnd_post, bc, np.where(Ufix, ab, -phiB_))
    return dict(lower=lower, upper=upper, diag=diag, D_rel=Drel, Db=Db,
                mob=V/Drel, bnd_post=bnd_post)

def grad_p(p_):
    pf = w*p_[owner] + (1 - w)*p_[nei]
    g = np.zeros((N, 3))
    np.add.at(g, owner, sfI*pf[:, None])
    np.add.at(g, nei, -(sfI*pf[:, None]))
    pb = p_[bc].copy()
    pb[patch_of_bface == OUTLET_IDX] = 0.0
    np.add.at(g, bc, bSf*pb[:, None])
    return g/V[:, None]

def interp_vec(F):
    return w[:, None]*F[owner] + (1-w)[:, None]*F[nei]

def phi_map(U_, p_, bs, press=True, bndpost=True):
    mob = bs["mob"]
    relaxD = bs["D_rel"] + (bs["bnd_post"] if bndpost else 0.0)
    src = (relaxD - bs["diag"])[:, None]*U_
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, bs["upper"][:, None]*U_[nei])
    np.add.at(offU, nei, bs["lower"][:, None]*U_[owner])
    HbyA = mob[:, None]*(src - offU)/V[:, None]
    if press:
        HbyA = HbyA - mob[:, None]*grad_p(p_)
    mobF = w*mob[owner] + (1 - w)*mob[nei]
    kf = mobF*dc*magSfI
    phi_int = np.einsum("fi,fi->f", sfI, interp_vec(HbyA)) - kf*(p_[nei] - p_[owner])
    phi_b = np.zeros(nB)
    for bf in bidx:
        c = bc[bf]
        phi_b[bf] = bSf[bf] @ HbyA[c] + mob[c]*bdel[bf]*bmag[bf]*p_[c]
    return np.concatenate([phi_int, phi_b])


U_B = np.loadtxt(SB + "/wstate_baseline_U.mtx")
p_B = np.loadtxt(SB + "/wstate_baseline_p.mtx")
phi_B = np.loadtxt(SB + "/wstate_baseline_phi.mtx")
phiB_B = np.loadtxt(SB + "/wstate_baseline_phiB.mtx")
nut = np.loadtxt(SB + "/wstate_baseline_nutFrozen.mtx")
al = np.loadtxt(SB + "/wstate_baseline_alpha.mtx")
nuE = NU + nut
mob_ex = np.loadtxt(SB + "/wstate_baseline_primalPressureMobility.mtx")
bsB = build_basis(phi_B, phiB_B, al, nuE)
assert np.linalg.norm(bsB["mob"] - mob_ex)/np.linalg.norm(mob_ex) < 1e-12

def grad_p(p_):
    pf = w*p_[owner] + (1 - w)*p_[nei]
    g = np.zeros((N, 3))
    np.add.at(g, owner, sfI*pf[:, None])
    np.add.at(g, nei, -(sfI*pf[:, None]))
    pb = p_[bc].copy()
    pb[patch_of_bface == OUTLET_IDX] = 0.0
    np.add.at(g, bc, bSf*pb[:, None])
    return g/V[:, None]

def interp_vec(F):
    return w[:, None]*F[owner] + (1-w)[:, None]*F[nei]

def phi_map(U_, p_, bs):
    mob = bs["mob"]
    relaxD = bs["D_rel"] + bs["bnd_post"]
    src = (relaxD - bs["diag"])[:, None]*U_
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, bs["upper"][:, None]*U_[nei])
    np.add.at(offU, nei, bs["lower"][:, None]*U_[owner])
    HbyA = mob[:, None]*(src - offU)/V[:, None] - mob[:, None]*grad_p(p_)
    mobF = w*mob[owner] + (1 - w)*mob[nei]
    kf = mobF*dc*magSfI
    phi = np.einsum("fi,fi->f", sfI, interp_vec(HbyA)) - kf*(p_[nei] - p_[owner])
    pb = np.zeros(nB)
    for bf in bidx:
        c = bc[bf]
        pb[bf] = bSf[bf] @ HbyA[c] + mob[c]*bdel[bf]*bmag[bf]*p_[c]
    return np.concatenate([phi, pb])

res = {}
ph = phi_map(U_B, p_B, bsB)
res["value_gate"] = float(np.linalg.norm(ph[:nIF] - phi_B)/np.linalg.norm(phi_B))
log("[T9] corrected-map value gate: %.4f" % res["value_gate"])
res["AB"] = {}
for nm in ("D1", "D2", "D3"):
    pre = SB + "/wstate_%s_h0.001_" % nm
    U_p = np.loadtxt(pre + "p_U.mtx"); p_p = np.loadtxt(pre + "p_p.mtx")
    U_m = np.loadtxt(pre + "m_U.mtx"); p_m = np.loadtxt(pre + "m_p.mtx")
    phi_p = np.loadtxt(pre + "p_phi.mtx"); phi_m = np.loadtxt(pre + "m_phi.mtx")
    phiBp = np.loadtxt(pre + "p_phiB.mtx"); phiBm = np.loadtxt(pre + "m_phiB.mtx")
    dph = np.concatenate([(phi_p - phi_m)/2e-3, (phiBp - phiBm)/2e-3])
    phA = (phi_map(U_p, p_p, bsB) - phi_map(U_m, p_m, bsB))/2e-3
    alp = np.loadtxt(pre + "p_alpha.mtx"); alm = np.loadtxt(pre + "m_alpha.mtx")
    bsp = build_basis(phi_p, phiBp, alp, nuE)
    bsm = build_basis(phi_m, phiBm, alm, nuE)
    phB_ = (phi_map(U_p, p_p, bsp) - phi_map(U_m, p_m, bsm))/2e-3
    bs_pa = build_basis(phi_B, phiB_B, alp, nuE)
    bs_ma = build_basis(phi_B, phiB_B, alm, nuE)
    ph_alpha = (phi_map(U_p, p_p, bs_pa) - phi_map(U_m, p_m, bs_ma))/2e-3
    def m(a):
        return (float(a @ dph/(np.linalg.norm(a)*np.linalg.norm(dph))),
                float(np.linalg.norm(a - dph)/np.linalg.norm(dph)))
    ca, ra = m(phA); cb, rb = m(phB_); _, ral = m(ph_alpha)
    res["AB"][nm] = {"A_cos": ca, "A_relL2": ra, "B_cos": cb, "B_relL2": rb,
                     "alphaMatrix_relL2": ral}
    log("[T9] %s: A-form cos=%.4f relL2=%.4f ; B-form(+motion) cos=%.4f "
        "relL2=%.4f ; alpha-matrix-only relL2=%.4f" % (nm, ca, ra, cb, rb, ral))

json.dump(res, open(OUT + "/b23_t9_abform_corrected.json", "w"), indent=1)
log("\nDONE t=%.1fs" % (time.time()-t0))
