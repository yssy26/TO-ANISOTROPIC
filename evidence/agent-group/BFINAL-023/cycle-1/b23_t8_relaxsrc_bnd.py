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
log("=== BFINAL-023 cycle-1 T8: relax-source boundary post-piece ===")

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
res = {}

ph = phi_map(U_B, p_B, bsB, press=True, bndpost=True)
relI = float(np.linalg.norm(ph[:nIF] - phi_B)/np.linalg.norm(phi_B))
ph0 = phi_map(U_B, p_B, bsB, press=True, bndpost=False)
res["H8a"] = {"with_bndpost": relI, "without": float(
    np.linalg.norm(ph0[:nIF] - phi_B)/np.linalg.norm(phi_B))}
log("[H8a] value gate: without bnd-post %.4f -> with %.4f"
    % (res["H8a"]["without"], res["H8a"]["with_bndpost"]))
err = ph[:nIF] - phi_B
sfh = sfI/magSfI[:, None]; axx = np.abs(sfh[:, 0])
bndc = np.zeros(N, bool); bndc[bc] = True
adj = bndc[owner] | bndc[nei]
res["H8a"]["err_structure"] = {
    "xfrac": float((err[axx > 0.9]**2).sum()/(err**2).sum()),
    "bndadj_frac": float((err[adj]**2).sum()/(err**2).sum())}
log("    err structure: x %.3f bndadj %.3f"
    % (res["H8a"]["err_structure"]["xfrac"],
       res["H8a"]["err_structure"]["bndadj_frac"]))

def channels(bs, U_, p_, dU, dp, da, press=True, bndpost=True, da_on=True):
    A_u = np.abs(bs["Db"])/V
    mob = bs["mob"]
    mobF = w*mob[owner] + (1 - w)*mob[nei]
    kf = mobF*dc*magSfI
    dphi_p = -kf*(dp[nei] - dp[owner])
    if press:
        dphi_p -= np.einsum("fi,fi->f", sfI, interp_vec(mob[:, None]*grad_p(dp)))
    dphi_p_b = np.zeros(nB)
    for bf in bidx:
        c = bc[bf]
        hb = mob[c]*bdel[bf]*bmag[bf]*dp[c]
        if press:
            hb -= bSf[bf] @ (mob[c]*grad_p(dp)[c])
        dphi_p_b[bf] = hb
    offdU = np.zeros((N, 3))
    np.add.at(offdU, owner, bs["upper"][:, None]*dU[nei])
    np.add.at(offdU, nei, bs["lower"][:, None]*dU[owner])
    dH = -offdU/V[:, None]
    fac = ALPHA_REL/A_u
    ccoef = np.full(N, 1.0 - ALPHA_REL) \
        + (mob*bs["bnd_post"]/V if bndpost else 0.0)
    dHbyA = fac[:, None]*dH + ccoef[:, None]*dU
    dphi_u = np.einsum("fi,fi->f", sfI, interp_vec(dHbyA))
    dphi_u_b = np.zeros(nB)
    for bf in bidx:
        c = bc[bf]
        dphi_u_b[bf] = bSf[bf] @ dHbyA[c]
    if not da_on:
        return (np.concatenate([dphi_u, dphi_u_b]),
                np.concatenate([dphi_p, dphi_p_b]),
                np.zeros_like(np.concatenate([dphi_p, dphi_p_b])))
    gp = grad_p(p_)
    drAU = -(mob**2/ALPHA_REL)*da
    relaxD = bs["D_rel"] + (bs["bnd_post"] if bndpost else 0.0)
    src = (relaxD - bs["diag"])[:, None]*U_
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, bs["upper"][:, None]*U_[nei])
    np.add.at(offU, nei, bs["lower"][:, None]*U_[owner])
    Hcur = (src - offU)/V[:, None] - (gp if press else 0.0)
    dHsrc = ((1.0/ALPHA_REL - 1.0)*da)[:, None]*U_
    dHbyA_a = drAU[:, None]*Hcur + mob[:, None]*dHsrc
    dphi_a = np.einsum("fi,fi->f", sfI, interp_vec(dHbyA_a))
    dphi_a += (drAU[owner]*w + drAU[nei]*(1-w))*dc*magSfI*(p_[nei] - p_[owner])
    if press:
        dphi_a -= np.einsum("fi,fi->f", sfI,
            (drAU[owner]*w + drAU[nei]*(1-w))[:, None]*interp_vec(gp))
    dphi_a_b = np.zeros(nB)
    for bf in bidx:
        c = bc[bf]
        dphi_a_b[bf] = bSf[bf] @ dHbyA_a[c] + drAU[c]*bdel[bf]*bmag[bf]*p_[c]
    return (np.concatenate([dphi_u, dphi_u_b]),
            np.concatenate([dphi_p, dphi_p_b]),
            np.concatenate([dphi_a, dphi_a_b]))

res["H8b"] = {}
for nm in ("D1", "D2", "D3"):
    pre = SB + "/wstate_%s_h0.001_" % nm
    U_p = np.loadtxt(pre + "p_U.mtx"); p_p = np.loadtxt(pre + "p_p.mtx")
    U_m = np.loadtxt(pre + "m_U.mtx"); p_m = np.loadtxt(pre + "m_p.mtx")
    phi_p = np.loadtxt(pre + "p_phi.mtx"); phi_m = np.loadtxt(pre + "m_phi.mtx")
    phiBp = np.loadtxt(pre + "p_phiB.mtx"); phiBm = np.loadtxt(pre + "m_phiB.mtx")
    dphi_true = np.concatenate([(phi_p - phi_m)/2e-3, (phiBp - phiBm)/2e-3])
    dU = (U_p - U_m)/2e-3; dp_ = (p_p - p_m)/2e-3
    ap = np.loadtxt(pre + "p_alpha.mtx"); am = np.loadtxt(pre + "m_alpha.mtx")
    da = (ap - am)/2e-3
    out = {}
    for cfg, kw in (("du_dp_press", dict(da_on=False, bndpost=True)),
                    ("full3", dict(da_on=True, bndpost=True)),
                    ("full3_nobndpost", dict(da_on=True, bndpost=False))):
        duE, dpE, daE = channels(bsB, U_B, p_B, dU, dp_, da, **kw)
        tot = duE + dpE + daE
        out[cfg] = {"cos": float(tot @ dphi_true/(np.linalg.norm(tot)*np.linalg.norm(dphi_true))),
                    "relL2": float(np.linalg.norm(tot - dphi_true)/np.linalg.norm(dphi_true))}
    # direction self-check for the bndpost correction (du_dp config)
    duA, dpA, _ = channels(bsB, U_B, p_B, dU, dp_, da, da_on=False, bndpost=False)
    duB, dpB, _ = channels(bsB, U_B, p_B, dU, dp_, da, da_on=False, bndpost=True)
    corr = float(np.corrcoef((duB - duA)[:nIF], (duA + dpA - dphi_true)[:nIF])[0, 1])
    out["corr_bndpost_vs_residual"] = corr
    res["H8b"][nm] = out
    log("[H8b] %s: du_dp_press %.4f (cos %.4f) ; full3 %.4f ; full3_nobndpost "
        "%.4f ; corr(bndpost-corr, resid)=%.3f"
        % (nm, out["du_dp_press"]["relL2"], out["du_dp_press"]["cos"],
           out["full3"]["relL2"], out["full3_nobndpost"]["relL2"], corr))

json.dump(res, open(OUT + "/b23_t8_relaxsrc_bnd.json", "w"), indent=1)
log("\nDONE t=%.1fs" % (time.time()-t0))
