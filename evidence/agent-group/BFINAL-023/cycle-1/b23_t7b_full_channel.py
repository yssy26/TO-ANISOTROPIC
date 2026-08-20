#!/usr/bin/env python3
"""BFINAL-023 cycle-1 T7b: combined exact instrument (corrected basis +
minus phiHbyA pressure channel, full linearization incl. the da-channel
piece and the interp-product term) -> final value gate + anchors +
residual structure.  Within preregistered T7 scope (complete linearization
of the confirmed channel; no new hypothesis).
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
log("=== BFINAL-023 cycle-1 T7b: combined exact instrument ===")

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
    return dict(lower=lower, upper=upper, diag=diag, D_rel=Drel, Db=Db, mob=V/Drel)

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

def phi_map(U_, p_, bs, press=True):
    mob = bs["mob"]
    src = (bs["D_rel"] - bs["diag"])[:, None]*U_
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, bs["upper"][:, None]*U_[nei])
    np.add.at(offU, nei, bs["lower"][:, None]*U_[owner])
    HbyA = mob[:, None]*(src - offU)/V[:, None]
    gp = grad_p(p_) if press else None
    if press:
        HbyA = HbyA - mob[:, None]*gp
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

# value gate
ph1 = phi_map(U_B, p_B, bsB, press=True)
relI = float(np.linalg.norm(ph1[:nIF] - phi_B)/np.linalg.norm(phi_B))
res["value_gate_with_channel"] = relI
log("[T7b] value gate with full channel: %.4f" % relI)
err = ph1[:nIF] - phi_B
sfh = sfI/magSfI[:, None]
axx = np.abs(sfh[:, 0])
bndc = np.zeros(N, bool); bndc[bc] = True
adj = bndc[owner] | bndc[nei]
res["value_err_structure"] = {
    "xfrac": float((err[axx > 0.9]**2).sum()/(err**2).sum()),
    "bndadj_frac": float((err[adj]**2).sum()/(err**2).sum()),
    "fluidowner_frac": float((err[al[owner] <= 1e3]**2).sum()/(err**2).sum()),
}
log("    err structure: x-normal %.3f ; bndadj %.3f ; fluid-owner %.3f"
    % (res["value_err_structure"]["xfrac"],
       res["value_err_structure"]["bndadj_frac"],
       res["value_err_structure"]["fluidowner_frac"]))

# combined tangent
def channels_full(bs, U_, p_, dU, dp, da, press=True):
    A_u = np.abs(bs["Db"])/V
    mob = bs["mob"]
    mobF = w*mob[owner] + (1 - w)*mob[nei]
    kf = mobF*dc*magSfI
    dphi_p = -kf*(dp[nei] - dp[owner])
    dphi_a_press = np.zeros(nIF)
    if press:
        gp_dp = mob[:, None]*grad_p(dp)
        dphi_p -= np.einsum("fi,fi->f", sfI, interp_vec(gp_dp))
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
    dHbyA = fac[:, None]*dH + (1.0 - ALPHA_REL)*dU
    dphi_u = np.einsum("fi,fi->f", sfI, interp_vec(dHbyA))
    dphi_u_b = np.zeros(nB)
    for bf in bidx:
        c = bc[bf]
        dphi_u_b[bf] = bSf[bf] @ dHbyA[c]
    # da channel: H now includes the -grad p source; kf interp as before;
    # plus the design derivative of the NEW channel: -interp(drAU*grad p).Sf
    gp = grad_p(p_)
    drAU = -(mob**2/ALPHA_REL)*da
    src = (bs["D_rel"] - bs["diag"])[:, None]*U_
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, bs["upper"][:, None]*U_[nei])
    np.add.at(offU, nei, bs["lower"][:, None]*U_[owner])
    Hcur = (src - offU)/V[:, None] - (gp if press else 0.0)
    dHsrc = ((1.0/ALPHA_REL - 1.0)*da)[:, None]*U_
    dHbyA_a = drAU[:, None]*Hcur + mob[:, None]*dHsrc
    dphi_a = np.einsum("fi,fi->f", sfI, interp_vec(dHbyA_a))
    dphi_a += (drAU[owner]*w + drAU[nei]*(1-w))*dc*magSfI*(p_[nei] - p_[owner])
    if press:
        dphi_a_press = -np.einsum(
            "fi,fi->f",
            sfI, (drAU[owner]*w + drAU[nei]*(1-w))[:, None]*interp_vec(gp))
        dphi_a += dphi_a_press
    dphi_a_b = np.zeros(nB)
    for bf in bidx:
        c = bc[bf]
        dphi_a_b[bf] = bSf[bf] @ dHbyA_a[c] + drAU[c]*bdel[bf]*bmag[bf]*p_[c]
    return (np.concatenate([dphi_u, dphi_u_b]),
            np.concatenate([dphi_p, dphi_p_b]),
            np.concatenate([dphi_a, dphi_a_b]))

res["anchors"] = {}
res["residual_structure"] = {}
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
    duE, dpE, daE = channels_full(bsB, U_B, p_B, dU, dp_, da, press=True)
    tot = duE + dpE + daE
    m = {"cos": float(tot @ dphi_true/(np.linalg.norm(tot)*np.linalg.norm(dphi_true))),
         "relL2": float(np.linalg.norm(tot - dphi_true)/np.linalg.norm(dphi_true)),
         "chan": {"du": float(np.linalg.norm(duE)/np.linalg.norm(tot)),
                  "dp": float(np.linalg.norm(dpE)/np.linalg.norm(tot)),
                  "da": float(np.linalg.norm(daE)/np.linalg.norm(tot))}}
    res["anchors"][nm] = m
    log("[T7b] %s: relL2=%.4f cos=%.4f (chan du/dp/da %.2f/%.2f/%.2f)"
        % (nm, m["relL2"], m["cos"], m["chan"]["du"], m["chan"]["dp"], m["chan"]["da"]))
    if nm == "D2":
        r = (tot - dphi_true)[:nIF]
        proj = {k: float((v[:nIF] @ r)/(np.linalg.norm(r)**2))
                for k, v in (("du", duE), ("dp", dpE), ("da", daE))}
        loo = {k: float(np.linalg.norm((tot - v - dphi_true)[:nIF])
                        /np.linalg.norm(dphi_true[:nIF]))
               for k, v in (("none", 0*duE), ("du", duE), ("dp", dpE), ("da", daE))}
        energy = r**2
        thr = np.sort(energy)[-max(nIF//10, 1)]
        carrier = (energy >= thr) & (axx > 0.9)
        res["residual_structure"]["D2"] = {
            "proj": proj, "leave_one_out": loo,
            "carrier_relL2": float(np.linalg.norm(r[carrier])
                                   /np.linalg.norm(dphi_true[:nIF][carrier])),
            "carrier_frac_energy": float(energy[carrier].sum()/energy.sum()),
            "carrier_xfrac_energy": float((energy[axx > 0.9]).sum()/energy.sum()),
        }
        log("    D2 residual: proj du/dp/da = %.3f/%.3f/%.3f ; LOO none/-du/-dp/-da "
            "= %.4f/%.4f/%.4f/%.4f ; carrier relL2=%.4f (energy %.2f, x-faces %.2f)"
            % (proj["du"], proj["dp"], proj["da"], loo["none"], loo["du"],
               loo["dp"], loo["da"],
               res["residual_structure"]["D2"]["carrier_relL2"],
               res["residual_structure"]["D2"]["carrier_frac_energy"],
               res["residual_structure"]["D2"]["carrier_xfrac_energy"]))

json.dump(res, open(OUT + "/b23_t7b_full_channel.json", "w"), indent=1)
log("\nDONE t=%.1fs" % (time.time()-t0))
