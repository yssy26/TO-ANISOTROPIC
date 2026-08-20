#!/usr/bin/env python3
"""BFINAL-023 cycle-1 T7: phiHbyA pressure channel -interp(rAtU*grad p)&Sf.
Preregistered in PREREGISTRATION_T7_ADDENDUM.md.

H7a  value gate (+/- sign discrimination).
H7b  anchors with the dp-channel extension -interp(rAtU*grad dp)&Sf.
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
log("=== BFINAL-023 cycle-1 T7: phiHbyA pressure channel ===")

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
log("mesh ready (t=%.1fs)" % (time.time()-t0))

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
    return dict(lower=lower, upper=upper, diag=diag, D_rel=Drel, Db=Db, mob=V/Drel)

def grad_p(p_):
    """Green-Gauss gradient of p (Gauss linear). Boundary: outlet
    fixedValue 0, everything else zeroGradient (cell value)."""
    pf = w*p_[owner] + (1 - w)*p_[nei]
    g = np.zeros((N, 3))
    np.add.at(g, owner, sfI*pf[:, None])
    np.add.at(g, nei, -(sfI*pf[:, None]))
    pb = p_[bc].copy()
    pb[patch_of_bface == OUTLET_IDX] = 0.0
    np.add.at(g, bc, bSf*pb[:, None])
    return g/V[:, None]

def phi_map(U_, p_, bs, extra_p_channel=0.0):
    mob = bs["mob"]
    src = (bs["D_rel"] - bs["diag"])[:, None]*U_
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, bs["upper"][:, None]*U_[nei])
    np.add.at(offU, nei, bs["lower"][:, None]*U_[owner])
    HbyA = mob[:, None]*(src - offU)/V[:, None]
    if extra_p_channel != 0.0:
        HbyA = HbyA - extra_p_channel*mob[:, None]*grad_p(p_)
    mobF = w*mob[owner] + (1 - w)*mob[nei]
    kf = mobF*dc*magSfI
    phi_int = np.einsum("fi,fi->f", sfI, w[:, None]*HbyA[owner]
                        + (1-w)[:, None]*HbyA[nei]) - kf*(p_[nei] - p_[owner])
    phi_b = np.zeros(nB)
    for bf in bidx:
        c = bc[bf]
        Hb = HbyA[c] - (extra_p_channel*mob[c]*grad_p(p_)[c]
                        if extra_p_channel != 0.0 else 0.0)
        phi_b[bf] = bSf[bf] @ Hb + mob[c]*bdel[bf]*bmag[bf]*p_[c]
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

# ---------------- H7a/H7c: value gate with sign discrimination ----------------
def relI(ph): return float(np.linalg.norm(ph[:nIF] - phi_B)/np.linalg.norm(phi_B))
ph0 = phi_map(U_B, p_B, bsB, 0.0)
phm = phi_map(U_B, p_B, bsB, -1.0)   # HbyA -= mob*grad p  (minus channel)
php = phi_map(U_B, p_B, bsB, +1.0)   # wrong-sign control
res["H7a"] = {"no_channel": relI(ph0), "minus_channel": relI(phm),
              "plus_channel": relI(php)}
log("[H7a] value: no-channel %.4f ; with -interp(mob*grad p).Sf %.4f ; "
    "wrong-sign + %.4f" % (relI(ph0), relI(phm), relI(php)))
# magnitude of the channel
ch = (phm - ph0)[:nIF]
res["H7a"]["channel_over_phi"] = float(np.linalg.norm(ch)/np.linalg.norm(phi_B))
res["H7a"]["channel_over_residual"] = float(
    np.linalg.norm(ch)/np.linalg.norm(ph0[:nIF] - phi_B))
log("    |channel|/|phi|=%.3f ; |channel|/|residual|=%.3f"
    % (res["H7a"]["channel_over_phi"], res["H7a"]["channel_over_residual"]))

# ---------------- H7b: anchors with dp-channel extension ----------------
def channels(bs, U_, p_, dU, dp, da, dp_ext=False):
    A_u = np.abs(bs["Db"])/V
    mob = bs["mob"]
    mobF = w*mob[owner] + (1 - w)*mob[nei]
    kf = mobF*dc*magSfI
    dphi_p = -kf*(dp[nei] - dp[owner])
    if dp_ext:
        dH_p = -mob[:, None]*grad_p(dp)
        dphi_p += np.einsum("fi,fi->f", sfI, w[:, None]*dH_p[owner]
                            + (1-w)[:, None]*dH_p[nei])
    dphi_p_b = np.zeros(nB)
    for bf in bidx:
        c = bc[bf]
        hb = mob[c]*bdel[bf]*bmag[bf]*dp[c]
        if dp_ext:
            hb -= bSf[bf] @ (mob[c]*grad_p(dp)[c])
        dphi_p_b[bf] = hb
    offdU = np.zeros((N, 3))
    np.add.at(offdU, owner, bs["upper"][:, None]*dU[nei])
    np.add.at(offdU, nei, bs["lower"][:, None]*dU[owner])
    dH = -offdU/V[:, None]
    fac = ALPHA_REL/A_u
    dHbyA = fac[:, None]*dH + (1.0 - ALPHA_REL)*dU
    dphi_u = np.einsum("fi,fi->f", sfI, w[:, None]*dHbyA[owner]
                       + (1-w)[:, None]*dHbyA[nei])
    dphi_u_b = np.zeros(nB)
    for bf in bidx:
        c = bc[bf]
        dphi_u_b[bf] = bSf[bf] @ dHbyA[c]
    drAU = -(mob**2/ALPHA_REL)*da
    src = (bs["D_rel"] - bs["diag"])[:, None]*U_
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, bs["upper"][:, None]*U_[nei])
    np.add.at(offU, nei, bs["lower"][:, None]*U_[owner])
    Hcur = (src - offU)/V[:, None]
    dHsrc = ((1.0/ALPHA_REL - 1.0)*da)[:, None]*U_
    dHbyA_a = drAU[:, None]*Hcur + mob[:, None]*dHsrc
    dphi_a = np.einsum("fi,fi->f", sfI, w[:, None]*dHbyA_a[owner]
                       + (1-w)[:, None]*dHbyA_a[nei])
    dphi_a += (drAU[owner]*w + drAU[nei]*(1-w))*dc*magSfI*(p_[nei] - p_[owner])
    dphi_a_b = np.zeros(nB)
    for bf in bidx:
        c = bc[bf]
        dphi_a_b[bf] = bSf[bf] @ dHbyA_a[c] + drAU[c]*bdel[bf]*bmag[bf]*p_[c]
    return (np.concatenate([dphi_u, dphi_u_b]),
            np.concatenate([dphi_p, dphi_p_b]),
            np.concatenate([dphi_a, dphi_a_b]))

res["H7b"] = {}
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
    du0, dp0, daE = channels(bsB, U_B, p_B, dU, dp_, da, dp_ext=False)
    du1, dp1, _ = channels(bsB, U_B, p_B, dU, dp_, da, dp_ext=True)
    tot0 = du0 + dp0 + daE; tot1 = du1 + dp1 + daE
    def mets(a):
        return {"cos": float(a @ dphi_true/(np.linalg.norm(a)*np.linalg.norm(dphi_true))),
                "relL2": float(np.linalg.norm(a - dphi_true)/np.linalg.norm(dphi_true))}
    res["H7b"][nm] = {"no_ext": mets(tot0), "with_ext": mets(tot1)}
    corr = float(np.corrcoef((tot1 - tot0)[:nIF], (tot0 - dphi_true)[:nIF])[0, 1])
    res["H7b"][nm]["corr_correction_vs_residual"] = corr
    res["H7b"][nm]["ext_over_residual"] = float(
        np.linalg.norm((tot1 - tot0)[:nIF])/np.linalg.norm((tot0 - dphi_true)[:nIF]))
    log("[H7b] %s: no-ext relL2=%.4f cos=%.4f -> with-ext relL2=%.4f cos=%.4f ; "
        "corr(corr, resid)=%.3f ; |ext|/|resid|=%.3f"
        % (nm, mets(tot0)["relL2"], mets(tot0)["cos"], mets(tot1)["relL2"],
           mets(tot1)["cos"], corr, res["H7b"][nm]["ext_over_residual"]))

json.dump(res, open(OUT + "/b23_t7_phihbya_pressure.json", "w"), indent=1)
log("\nDONE t=%.1fs" % (time.time()-t0))
