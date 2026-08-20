#!/usr/bin/env python3
"""BFINAL-023 cycle-1 T5: matrix-motion discrimination (preregistered in
PREREGISTRATION_T5_ADDENDUM.md).

  H5a  value-level gate: corrected-basis phi_map(U_B,p_B,a_B;matrix(phi_B))
       reproduces phi_B internal relL2 < 1%.
  H5b  A-form (states at +/-, matrix frozen at B) vs B-form (matrix
       re-assembled at +/- states) against dphi_true.
  H5c  channel attribution at fixed matrix on the D2 carrier faces
       (residual projections + leave-one-out).
  H5d  upwind donor-flip census on the carrier.

All on the corrected (production-sign) basis validated to 2.2e-15 against
the exact exported mobility in T1.
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
log("=== BFINAL-023 cycle-1 T5: matrix-motion discrimination ===")

# ---------------- mesh ----------------
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
def read_labels(path):
    txt = open(path).read(); i = txt.index("\n("); j = txt.index("\n)", i)
    return np.array([int(t) for t in txt[i+2:j].split()], dtype=np.int64)
owner_all = read_labels(POLY + "/owner"); nei_all = read_labels(POLY + "/neighbour")
nIF = len(nei_all); nF = len(faces)
Sf_all = np.zeros((nF, 3))
for fi in range(nF):
    fv = pts[faces[fi]]
    Sf_all[fi] = 0.5*np.cross(fv, np.roll(fv, -1, axis=0)).sum(axis=0)
magSf_all = np.linalg.norm(Sf_all, axis=1)
cellverts = {}
for fi in range(nF):
    cells = (owner_all[fi], nei_all[fi]) if fi < nIF else (owner_all[fi],)
    for c in cells:
        cellverts.setdefault(c, set()).update(faces[fi])
Ccen = np.zeros((N, 3))
for c, vs in cellverts.items(): Ccen[c] = pts[sorted(vs)].mean(axis=0)
owner = owner_all[:nIF]; nei = nei_all[:nIF]
d_vec = Ccen[nei] - Ccen[owner]
deltaCoeffs = magSf_all[:nIF]/np.einsum("ij,ij->i", Sf_all[:nIF], d_vec)
V = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")[:, 0]
w = np.loadtxt(MESH + "/b16mesh_weights.mtx")
sf = Sf_all[:nIF]
bnd = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
bc_cell = bnd[:, 0].astype(np.int64)
bSf = bnd[:, 1:4]; bmag = bnd[:, 4]; bdel = bnd[:, 5]
Ufixed_b = bnd[:, 7] > 0.5
nBnd = len(bc_cell)
bidx = np.where(~Ufixed_b)[0]
log("mesh ready (t=%.1fs)" % (time.time()-t0))

def read_of_internal(path):
    txt = open(path).read()
    i = txt.index("internalField")
    j = txt.index("\n(", i); k = txt.index("\n)", j)
    vals = []
    for line in txt[j+1:k].splitlines():
        for tok in line.replace(";", " ").split():
            for tt in tok.replace("(", " ").replace(")", " ").split():
                vals.append(float(tt))
    return np.array(vals)

# ---------------- corrected-basis machinery (T1-validated) ----------------
def build_basis(phi_, phiB_, alpha_, nuEff_, nuEff_is_full=True):
    upw_ = phi_ >= 0.0
    qp_ = upw_.astype(float); qn_ = 1.0 - qp_
    divphi_ = np.zeros(N)
    np.add.at(divphi_, owner, phi_); np.add.at(divphi_, nei, -phi_)
    np.add.at(divphi_, bc_cell, phiB_)
    gam = w*nuEff_[owner] + (1.0 - w)*nuEff_[nei]
    kfv = gam*magSf_all[:nIF]*deltaCoeffs
    lower_ = -qp_*phi_ - kfv
    upper_ = qn_*phi_ - kfv
    diag_ = np.zeros(N)
    np.add.at(diag_, owner, -lower_)
    np.add.at(diag_, nei, -upper_)
    diag_ -= divphi_
    diag_ += alpha_*V
    sumOff_ = np.zeros(N)
    np.add.at(sumOff_, owner, np.abs(upper_)); np.add.at(sumOff_, nei, np.abs(lower_))
    ab_ = nuEff_[bc_cell]*bmag*bdel
    ic_ = np.where(Ufixed_b, -ab_, phiB_)
    Db_ = diag_.copy()
    np.add.at(Db_, bc_cell, np.abs(ic_))
    D_rel_ = np.maximum(np.abs(Db_), sumOff_)/ALPHA_REL
    return dict(lower=lower_, upper=upper_, diag=diag_, D_rel=D_rel_,
                Db=Db_, sumOff=sumOff_, ic=ic_, qp=qp_, qn=qn_,
                mob=V/D_rel_)

def phi_map(U_, p_, bs):
    """SIMPLE flux reconstruction: interp(rAtU*H)&Sf - kf*grad(p)."""
    mob = bs["mob"]
    src = (bs["D_rel"] - bs["diag"])[:, None]*U_
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, bs["upper"][:, None]*U_[nei])
    np.add.at(offU, nei, bs["lower"][:, None]*U_[owner])
    HbyA = mob[:, None]*(src - offU)/V[:, None]
    mobF = w*mob[owner] + (1.0 - w)*mob[nei]
    kf = mobF*deltaCoeffs*magSf_all[:nIF]
    phi_int = np.einsum("fi,fi->f", sf, w[:, None]*HbyA[owner]
                        + (1-w)[:, None]*HbyA[nei])
    phi_int -= kf*(p_[nei] - p_[owner])
    phi_bnd = np.zeros(nBnd)
    for bf in bidx:
        c = bc_cell[bf]
        phi_bnd[bf] = bSf[bf] @ HbyA[c] + mob[c]*bdel[bf]*bmag[bf]*p_[c]
    return np.concatenate([phi_int, phi_bnd])

# ---------------- state B ----------------
U_B = np.loadtxt(SB + "/wstate_baseline_U.mtx")
p_B = np.loadtxt(SB + "/wstate_baseline_p.mtx")
phi_B = np.loadtxt(SB + "/wstate_baseline_phi.mtx")
phiB_B = np.loadtxt(SB + "/wstate_baseline_phiB.mtx")
nut_B = np.loadtxt(SB + "/wstate_baseline_nutFrozen.mtx")
alpha_B = np.loadtxt(SB + "/wstate_baseline_alpha.mtx")
nuEff_B = NU + nut_B
mob_exact = np.loadtxt(SB + "/wstate_baseline_primalPressureMobility.mtx")
res = {}

bsB = build_basis(phi_B, phiB_B, alpha_B, nuEff_B)
assert np.linalg.norm(bsB["mob"] - mob_exact)/np.linalg.norm(mob_exact) < 1e-12

# ---------------- H5a: value-level gate ----------------
phi_re = phi_map(U_B, p_B, bsB)
errI = phi_re[:nIF] - phi_B
errO = phi_re[nIF:][bidx] - phiB_B[bidx]
res["H5a"] = {
    "internal_relL2": float(np.linalg.norm(errI)/np.linalg.norm(phi_B)),
    "outlet_relL2": float(np.linalg.norm(errO)/np.linalg.norm(phiB_B[bidx])),
    "internal_maxRelTop": float(np.percentile(
        np.abs(errI)/np.maximum(np.abs(phi_B), 1e-300), 99)),
}
log("[H5a] value rebuild @B: internal relL2=%.4e outlet relL2=%.3e p99rel=%.3e"
    % (res["H5a"]["internal_relL2"], res["H5a"]["outlet_relL2"],
       res["H5a"]["internal_maxRelTop"]))

# ---------------- H5b: A-form / B-form ----------------
def load_state(pre, sgn):
    U = np.loadtxt(pre + sgn + "_U.mtx")
    p = np.loadtxt(pre + sgn + "_p.mtx")
    phi = np.loadtxt(pre + sgn + "_phi.mtx")
    phiB = np.loadtxt(pre + sgn + "_phiB.mtx")
    alpha = np.loadtxt(pre + sgn + "_alpha.mtx")
    nut = np.loadtxt(pre + sgn + "_nutFrozen.mtx")
    return U, p, phi, phiB, alpha, NU + nut

res["H5b"] = {}
sfhat = sf/np.maximum(np.linalg.norm(sf, axis=1, keepdims=True), 1e-300)
absx = np.abs(sfhat[:, 0])
carrier = None
tang_cache = {}
for nm, tag in (("D1", "0.001"), ("D2", "0.001"), ("D3", "0.001")):
    h = float(tag)
    pre = SB + "/wstate_%s_h%s_" % (nm, tag)
    U_p, p_p, phi_p, phiB_p, alpha_p, nuEff_p = load_state(pre, "p")
    U_m, p_m, phi_m, phiB_m, alpha_m, nuEff_m = load_state(pre, "m")
    dphi_true = np.concatenate(
        [(phi_p - phi_m)/(2*h), (phiB_p - phiB_m)/(2*h)])
    # A-form: states +/-, matrix frozen at B
    phiA = (phi_map(U_p, p_p, bsB) - phi_map(U_m, p_m, bsB))/(2*h)
    # B-form: matrix re-assembled at +/- states
    bs_p = build_basis(phi_p, phiB_p, alpha_p, nuEff_p)
    bs_m = build_basis(phi_m, phiB_m, alpha_m, nuEff_m)
    phiBform = (phi_map(U_p, p_p, bs_p) - phi_map(U_m, p_m, bs_m))/(2*h)
    def mets(a):
        return {"cos": float(a @ dphi_true
                /(np.linalg.norm(a)*np.linalg.norm(dphi_true))),
                "relL2": float(np.linalg.norm(a - dphi_true)
                               /np.linalg.norm(dphi_true))}
    # matrix-motion contribution in isolation (B - A)
    mot = phiBform - phiA
    res["H5b"][nm] = {
        "Aform": mets(phiA), "Bform": mets(phiBform),
        "motion_over_truth": float(np.linalg.norm(mot[:nIF])
                                   /np.linalg.norm(dphi_true[:nIF])),
        "residA_over_truth": float(np.linalg.norm(
            (phiA - dphi_true)[:nIF])/np.linalg.norm(dphi_true[:nIF])),
        "corr_motion_residA": float(np.corrcoef(
            mot[:nIF], (phiA - dphi_true)[:nIF])[0, 1]),
    }
    log("[H5b] %s: A-form cos=%.4f relL2=%.4f ; B-form cos=%.4f relL2=%.4f ; "
        "|motion|/|truth|=%.3f corr(motion, residA)=%.3f"
        % (nm, res["H5b"][nm]["Aform"]["cos"], res["H5b"][nm]["Aform"]["relL2"],
           res["H5b"][nm]["Bform"]["cos"], res["H5b"][nm]["Bform"]["relL2"],
           res["H5b"][nm]["motion_over_truth"],
           res["H5b"][nm]["corr_motion_residA"]))
    if nm == "D2":
        energy = (phiA - dphi_true)[:nIF]**2
        thr = np.sort(energy)[-max(nIF//10, 1)]
        carrier = (energy >= thr) & (absx > 0.9)
        res["H5b"]["D2_motion_map"] = {
            "carrier_frac_energy": float(energy[carrier].sum()/energy.sum()),
        }
        # channel-level: which piece of matrix motion (upwind weights vs
        # mobility vs kf) carries it? partial B-forms:
        #   B1: matrix upwind/diag motion ONLY in H (mob from B)
        #   B2: full B-form (already)
        # cheap attribution: rebuild with alpha_+- and phi_B (design-driven
        # matrix motion only) vs phi_+- and alpha_B (flux-driven only)
        bs_pa = build_basis(phi_B, phiB_B, alpha_p, nuEff_B)
        bs_ma = build_basis(phi_B, phiB_B, alpha_m, nuEff_B)
        phiB_alphaOnly = (phi_map(U_p, p_p, bs_pa)
                          - phi_map(U_m, p_m, bs_ma))/(2*h)
        bs_pf = build_basis(phi_p, phiB_p, alpha_B, nuEff_B)
        bs_mf = build_basis(phi_m, phiB_m, alpha_B, nuEff_B)
        phiB_fluxOnly = (phi_map(U_p, p_p, bs_pf)
                         - phi_map(U_m, p_m, bs_mf))/(2*h)
        res["H5b"]["D2_partial"] = {
            "alphaDrivenMatrix": mets(phiB_alphaOnly),
            "fluxDrivenMatrix": mets(phiB_fluxOnly),
            "alphaDriven_motion_over_truth": float(np.linalg.norm(
                (phiB_alphaOnly - phiA)[:nIF])/np.linalg.norm(dphi_true[:nIF])),
            "fluxDriven_motion_over_truth": float(np.linalg.norm(
                (phiB_fluxOnly - phiA)[:nIF])/np.linalg.norm(dphi_true[:nIF])),
        }
        log("    D2 partial-forms: alpha-driven matrix cos=%.4f relL2=%.4f "
            "(motion|x %.3f) ; flux-driven cos=%.4f relL2=%.4f (motion|x %.3f)"
            % (res["H5b"]["D2_partial"]["alphaDrivenMatrix"]["cos"],
               res["H5b"]["D2_partial"]["alphaDrivenMatrix"]["relL2"],
               res["H5b"]["D2_partial"]["alphaDriven_motion_over_truth"],
               res["H5b"]["D2_partial"]["fluxDrivenMatrix"]["cos"],
               res["H5b"]["D2_partial"]["fluxDrivenMatrix"]["relL2"],
               res["H5b"]["D2_partial"]["fluxDriven_motion_over_truth"]))

# ---------------- H5c: channel attribution at fixed matrix (D2) ----------------
# channels from the T3 run (VC == production semantics here): rebuild them
def channels_fixed(bs, mob, U_, p_, dU, dp, da):
    A_u = np.abs(bs["Db"])/V
    mobF = w*mob[owner] + (1.0 - w)*mob[nei]
    kf = mobF*deltaCoeffs*magSf_all[:nIF]
    dphi_p = -kf*(dp[nei] - dp[owner])
    dphi_p_bnd = np.zeros(nBnd)
    for bf in bidx:
        c = bc_cell[bf]
        dphi_p_bnd[bf] = mob[c]*bdel[bf]*bmag[bf]*dp[c]
    offdU = np.zeros((N, 3))
    np.add.at(offdU, owner, bs["upper"][:, None]*dU[nei])
    np.add.at(offdU, nei, bs["lower"][:, None]*dU[owner])
    dH = -offdU/V[:, None]
    fac = ALPHA_REL/A_u
    dHbyA = fac[:, None]*dH + (1.0 - ALPHA_REL)*dU
    dphi_u = np.einsum("fi,fi->f", sf, w[:, None]*dHbyA[owner]
                       + (1-w)[:, None]*dHbyA[nei])
    dphi_u_bnd = np.zeros(nBnd)
    for bf in bidx:
        c = bc_cell[bf]
        dphi_u_bnd[bf] = bSf[bf] @ dHbyA[c]
    drAU = -(mob**2/ALPHA_REL)*da
    src = (bs["D_rel"] - bs["diag"])[:, None]*U_
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, bs["upper"][:, None]*U_[nei])
    np.add.at(offU, nei, bs["lower"][:, None]*U_[owner])
    Hcur = (src - offU)/V[:, None]
    dHsrc = ((1.0/ALPHA_REL - 1.0)*da)[:, None]*U_
    dHbyA_a = drAU[:, None]*Hcur + mob[:, None]*dHsrc
    dphi_a = np.einsum("fi,fi->f", sf, w[:, None]*dHbyA_a[owner]
                       + (1-w)[:, None]*dHbyA_a[nei])
    dphi_a += (drAU[owner]*w + drAU[nei]*(1-w))*deltaCoeffs*magSf_all[:nIF]*(
        p_[nei] - p_[owner])
    dphi_a_bnd = np.zeros(nBnd)
    for bf in bidx:
        c = bc_cell[bf]
        dphi_a_bnd[bf] = bSf[bf] @ dHbyA_a[c] + drAU[c]*bdel[bf]*bmag[bf]*p_[c]
    return (np.concatenate([dphi_u, dphi_u_bnd]),
            np.concatenate([dphi_p, dphi_p_bnd]),
            np.concatenate([dphi_a, dphi_a_bnd]))

pre = SB + "/wstate_D2_h0.001_"
U_p = np.loadtxt(pre + "p_U.mtx"); p_p = np.loadtxt(pre + "p_p.mtx")
U_m = np.loadtxt(pre + "m_U.mtx"); p_m = np.loadtxt(pre + "m_p.mtx")
phi_p = np.loadtxt(pre + "p_phi.mtx"); phi_m = np.loadtxt(pre + "m_phi.mtx")
phiBp = np.loadtxt(pre + "p_phiB.mtx"); phiBm = np.loadtxt(pre + "m_phiB.mtx")
dphi_true = np.concatenate([(phi_p - phi_m)/2e-3, (phiBp - phiBm)/2e-3])
dU = (U_p - U_m)/2e-3; dp_ = (p_p - p_m)/2e-3
alpha_p = np.loadtxt(pre + "p_alpha.mtx"); alpha_m = np.loadtxt(pre + "m_alpha.mtx")
da = (alpha_p - alpha_m)/2e-3
duE, dpE, daE = channels_fixed(bsB, mob_exact, U_B, p_B, dU, dp_, da)
tot = duE + dpE + daE
r = tot - dphi_true
mF = carrier
res["H5c"] = {
    "chan_norms": {k: float(np.linalg.norm(v[:nIF]))
                   for k, v in (("du", duE), ("dp", dpE), ("da", daE))},
    "proj_resid": {k: float((v[:nIF] @ r[:nIF])/(np.linalg.norm(r[:nIF])**2))
                   for k, v in (("du", duE), ("dp", dpE), ("da", daE))},
    "proj_resid_carrier": {
        k: float((v[:nIF][mF] @ r[:nIF][mF])/(np.linalg.norm(r[:nIF][mF])**2))
        for k, v in (("du", duE), ("dp", dpE), ("da", daE))},
    "leave_one_out_relL2": {
        k: float(np.linalg.norm((tot - v - dphi_true)[:nIF])
                 /np.linalg.norm(dphi_true[:nIF]))
        for k, v in (("none", 0*duE), ("du", duE), ("dp", dpE), ("da", daE))},
}
log("[H5c] D2 fixed-matrix channel attribution:")
log("    residual projections (all faces): du=%.3f dp=%.3f da=%.3f"
    % (res["H5c"]["proj_resid"]["du"], res["H5c"]["proj_resid"]["dp"],
       res["H5c"]["proj_resid"]["da"]))
log("    residual projections (carrier):   du=%.3f dp=%.3f da=%.3f"
    % (res["H5c"]["proj_resid_carrier"]["du"],
       res["H5c"]["proj_resid_carrier"]["dp"],
       res["H5c"]["proj_resid_carrier"]["da"]))
log("    leave-one-out relL2: none=%.4f -du=%.4f -dp=%.4f -da=%.4f"
    % tuple(res["H5c"]["leave_one_out_relL2"][k]
            for k in ("none", "du", "dp", "da")))

# ---------------- H5d: upwind donor-flip census ----------------
flip_p = (phi_B >= 0) != (phi_p >= 0)
flip_m = (phi_B >= 0) != (phi_m >= 0)
flip_any = flip_p | flip_m
# sensitivity: |phi| small relative to |kfv| -> upper structure dominated by
# convection sign locally
small_phi = np.abs(phi_B) < 0.1*(w*nuEff_B[owner] + (1-w)*nuEff_B[nei]) \
    *magSf_all[:nIF]*deltaCoeffs
res["H5d"] = {
    "flip_faces_frac": float(flip_any.mean()),
    "carrier_flip_frac": float(flip_any[carrier].mean()),
    "residA_energy_on_flipfaces": float(
        ((phiA - dphi_true)[:nIF][flip_any]**2).sum()
        /((phiA - dphi_true)[:nIF]**2).sum()),
    "carrier_smallphi_frac": float(small_phi[carrier].mean()),
    "truth_energy_on_flipfaces": float(
        (dphi_true[:nIF][flip_any]**2).sum()/(dphi_true[:nIF]**2).sum()),
}
log("[H5d] upwind flips (B vs +/-): all %.4f ; carrier %.4f ; "
    "residA energy on flip faces %.3f ; truth energy on flip faces %.3f ; "
    "carrier small-|phi| frac %.3f"
    % (res["H5d"]["flip_faces_frac"], res["H5d"]["carrier_flip_frac"],
       res["H5d"]["residA_energy_on_flipfaces"],
       res["H5d"]["truth_energy_on_flipfaces"],
       res["H5d"]["carrier_smallphi_frac"]))

json.dump(res, open(OUT + "/b23_t5_matrix_motion.json", "w"), indent=1)
log("\nDONE t=%.1fs" % (time.time()-t0))
