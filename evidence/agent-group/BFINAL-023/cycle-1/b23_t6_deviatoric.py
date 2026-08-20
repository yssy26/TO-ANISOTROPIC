#!/usr/bin/env python3
"""BFINAL-023 cycle-1 T6: explicit deviatoric-stress slot (preregistered in
PREREGISTRATION_T6_ADDENDUM.md).

H6a  value gate: phi_map + S_phys = -div(nuEff*dev2(T(grad U)))  vs phi_B.
H6b  anchors D1/D2/D3 with the deviatoric term added to the du channel
     (production semantics VP + dev).
H6c  correlation of the correction with the T5 residual.
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
U_INLET = np.array([193.97, 0.0, 0.0])

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-023 cycle-1 T6: deviatoric slot ===")

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
# boundary patch id per boundary face (from mesh boundary file ordering):
# read patch names+sizes from constant/polyMesh/boundary
bfile = open(POLY + "/boundary").read()
import re
patch_names = re.findall(r"^\s{4}(\w+)\s*$", bfile, re.M)
patch_sizes = [int(m) for m in re.findall(r"nFaces\s+(\d+);", bfile)]
patch_of_bface = np.concatenate(
    [np.full(s, i, dtype=np.int64) for i, s in enumerate(patch_sizes)])
log("patches: %s sizes %s" % (patch_names, patch_sizes))
log("mesh ready (t=%.1fs)" % (time.time()-t0))

# patch boundary U values (fixed) / zeroGradient flag
INLET_IDX = patch_names.index("inlet")
OUTLET_IDX = patch_names.index("outlet")

# ---------------- corrected basis machinery ----------------
def build_basis(phi_, phiB_, alpha_, nuEff_):
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
                Db=Db_, sumOff=sumOff_, ic=ic_, mob=V/D_rel_)

def dev2T(G):
    """dev2(T(G)) for cell gradient array G (N,3,3): T = transpose;
    dev2(A) = A - (2/3) tr(A) I."""
    A = np.transpose(G, (0, 2, 1))
    tr = np.trace(A, axis1=1, axis2=2)
    return A - (2.0/3.0)*tr[:, None, None]*np.eye(3)[None]

def green_gauss_grad(Uf_int, Uf_bnd):
    """G_c = (1/V) [sum_int Sf U_f + sum_bnd Sf_b U_b]; U_f = w U_o+(1-w)U_n."""
    Uf = w[:, None]*Uf_int[owner] + (1-w)[:, None]*Uf_int[nei]
    acc = np.zeros((N, 3, 3))
    np.add.at(acc, owner, sf[:, None, :]*Uf[:, :, None])  # Sf outer Uf
    np.add.at(acc, nei, -(sf[:, None, :]*Uf[:, :, None]))
    np.add.at(acc, bc_cell, bSf[:, None, :]*Uf_bnd[:, :, None])
    return acc/V[:, None, None]

def divDev(Uf_int, Uf_bnd, nuEff_, wall_onesided=False):
    """(1/V)[sum_int Sf.X_f + sum_bnd Sf_b.X_b], X = nuEff*dev2(T(grad U)).
    Gauss linear: X_f = w X_o + (1-w) X_n; X_b = X_c (calculated) or
    one-sided wall gradient variant."""
    G = green_gauss_grad(Uf_int, Uf_bnd)
    X = nuEff_[:, None, None]*dev2T(G)
    Xf = w[:, None, None]*X[owner] + (1-w)[:, None, None]*X[nei]
    acc = np.einsum("fi,fij->fj", sf, Xf) if False else np.einsum("fi,fij->fj", sf, Xf)
    div = np.zeros((N, 3))
    np.add.at(div, owner, np.einsum("fi,fij->fj", sf, Xf))
    np.add.at(div, nei, -np.einsum("fi,fij->fj", sf, Xf))
    Xb = X[bc_cell]
    if wall_onesided:
        n_hat = bSf/np.maximum(bmag[:, None], 1e-300)
        # one-sided wall gradient: G_b = delta*(n (x) (U_b - U_c))
        dU_wall = Uf_bnd - Uf_int[bc_cell]
        Gb = bdel[:, None, None]*(n_hat[:, :, None]*dU_wall[:, None, :])
        Xb = nuEff_[bc_cell][:, None, None]*dev2T(Gb)
    np.add.at(div, bc_cell, np.einsum("fi,fij->fj", bSf, Xb))
    return div/V[:, None]

def phi_map(U_, p_, bs, dev_div=None):
    mob = bs["mob"]
    src = (bs["D_rel"] - bs["diag"])[:, None]*U_
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, bs["upper"][:, None]*U_[nei])
    np.add.at(offU, nei, bs["lower"][:, None]*U_[owner])
    HbyA = mob[:, None]*(src - offU)/V[:, None]
    if dev_div is not None:
        HbyA -= mob[:, None]*dev_div
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

def U_bnd_values(U_):
    Uf_bnd = U_[bc_cell].copy()
    Uf_bnd[Ufixed_b] = 0.0                      # noSlip walls + hot in/outlets
    inlet_m = patch_of_bface == INLET_IDX
    Uf_bnd[inlet_m] = U_INLET                   # fixed inlet
    # outlet (zeroGradient): keep cell value
    return Uf_bnd

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

UfB_bnd = U_bnd_values(U_B)
devB = divDev(U_B, UfB_bnd, nuEff_B)
devB_wall = divDev(U_B, UfB_bnd, nuEff_B, wall_onesided=True)

# ---------------- H6a: value gate ----------------
phi_re0 = phi_map(U_B, p_B, bsB)
phi_re1 = phi_map(U_B, p_B, bsB, dev_div=devB)
phi_re2 = phi_map(U_B, p_B, bsB, dev_div=devB_wall)
def relI(ph):
    return float(np.linalg.norm(ph[:nIF] - phi_B)/np.linalg.norm(phi_B))
res["H6a"] = {
    "no_dev": relI(phi_re0), "dev_calc": relI(phi_re1),
    "dev_wall_onesided": relI(phi_re2),
    "dev_norm_over_phiHbyA": float(np.linalg.norm(
        bsB["mob"][:, None]*devB)/np.linalg.norm(phi_B)),
}
log("[H6a] value gate: no-dev %.4f -> dev(calc) %.4f ; dev(wall-1s) %.4f ; "
    "|mob*dev|/|phi| %.3f"
    % (res["H6a"]["no_dev"], res["H6a"]["dev_calc"],
       res["H6a"]["dev_wall_onesided"], res["H6a"]["dev_norm_over_phiHbyA"]))

# ---------------- H6b: anchors with dev in du channel ----------------
def channels(bs, U_, p_, dU, dp, da, dev_dU=None):
    A_u = np.abs(bs["Db"])/V
    mob = bs["mob"]
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
    if dev_dU is not None:
        dH = dH - dev_dU
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
    Hcur = (src - offU)/V[:, None] - devB
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

res["H6b"] = {}
res["H6c"] = {}
for nm in ("D1", "D2", "D3"):
    pre = SB + "/wstate_%s_h0.001_" % nm
    U_p = np.loadtxt(pre + "p_U.mtx"); p_p = np.loadtxt(pre + "p_p.mtx")
    U_m = np.loadtxt(pre + "m_U.mtx"); p_m = np.loadtxt(pre + "m_p.mtx")
    phi_p = np.loadtxt(pre + "p_phi.mtx"); phi_m = np.loadtxt(pre + "m_phi.mtx")
    phiBp = np.loadtxt(pre + "p_phiB.mtx"); phiBm = np.loadtxt(pre + "m_phiB.mtx")
    dphi_true = np.concatenate([(phi_p - phi_m)/2e-3, (phiBp - phiBm)/2e-3])
    dU = (U_p - U_m)/2e-3; dp_ = (p_p - p_m)/2e-3
    alpha_p = np.loadtxt(pre + "p_alpha.mtx"); alpha_m = np.loadtxt(pre + "m_alpha.mtx")
    da = (alpha_p - alpha_m)/2e-3
    dU_bnd = U_bnd_values(dU)   # walls/inlet zero, outlet cell
    devdU = divDev(dU, dU_bnd, nuEff_B)
    du0, dpE, daE = channels(bsB, U_B, p_B, dU, dp_, da)
    du1, _, _ = channels(bsB, U_B, p_B, dU, dp_, da, dev_dU=devdU)
    tot0 = du0 + dpE + daE; tot1 = du1 + dpE + daE
    def mets(a):
        return {"cos": float(a @ dphi_true/(np.linalg.norm(a)*np.linalg.norm(dphi_true))),
                "relL2": float(np.linalg.norm(a - dphi_true)/np.linalg.norm(dphi_true))}
    res["H6b"][nm] = {"no_dev": mets(tot0), "with_dev": mets(tot1)}
    corr = float(np.corrcoef((tot1 - tot0)[:nIF], (tot0 - dphi_true)[:nIF])[0, 1])
    res["H6c"][nm] = {
        "corr_correction_vs_residual": corr,
        "correction_over_residual": float(np.linalg.norm((tot1 - tot0)[:nIF])
                                          /np.linalg.norm((tot0 - dphi_true)[:nIF])),
        "devchan_norm_over_tot": float(np.linalg.norm((du1 - du0)[:nIF])
                                       /np.linalg.norm(tot1[:nIF])),
    }
    log("[H6b] %s: no-dev relL2=%.4f cos=%.4f -> with-dev relL2=%.4f cos=%.4f"
        % (nm, res["H6b"][nm]["no_dev"]["relL2"], res["H6b"][nm]["no_dev"]["cos"],
           res["H6b"][nm]["with_dev"]["relL2"], res["H6b"][nm]["with_dev"]["cos"]))
    log("[H6c]    corr(correction, residual)=%.3f ; |corr|/|resid|=%.3f ; "
        "devchan/|tot|=%.3f"
        % (corr, res["H6c"][nm]["correction_over_residual"],
           res["H6c"][nm]["devchan_norm_over_tot"]))

# also: value-form A/B-form with dev (D2 quick check, value FD of full map)
nm = "D2"
pre = SB + "/wstate_%s_h0.001_" % nm
U_p = np.loadtxt(pre + "p_U.mtx"); p_p = np.loadtxt(pre + "p_p.mtx")
U_m = np.loadtxt(pre + "m_U.mtx"); p_m = np.loadtxt(pre + "m_p.mtx")
phi_p = np.loadtxt(pre + "p_phi.mtx"); phi_m = np.loadtxt(pre + "m_phi.mtx")
phiBp = np.loadtxt(pre + "p_phiB.mtx"); phiBm = np.loadtxt(pre + "m_phiB.mtx")
dphi_true = np.concatenate([(phi_p - phi_m)/2e-3, (phiBp - phiBm)/2e-3])
dev_p = divDev(U_p, U_bnd_values(U_p), nuEff_B)
dev_m = divDev(U_m, U_bnd_values(U_m), nuEff_B)
phiA_dev = (phi_map(U_p, p_p, bsB, dev_div=dev_p)
            - phi_map(U_m, p_m, bsB, dev_div=dev_m))/2e-3
res["H6b"]["D2_Aform_with_dev"] = {
    "cos": float(phiA_dev @ dphi_true/(np.linalg.norm(phiA_dev)*np.linalg.norm(dphi_true))),
    "relL2": float(np.linalg.norm(phiA_dev - dphi_true)/np.linalg.norm(dphi_true))}
log("[H6b] D2 A-form (du+dp, no da) with dev: relL2=%.4f cos=%.4f"
    % (res["H6b"]["D2_Aform_with_dev"]["relL2"],
       res["H6b"]["D2_Aform_with_dev"]["cos"]))

json.dump(res, open(OUT + "/b23_t6_deviatoric.json", "w"), indent=1)
log("\nDONE t=%.1fs" % (time.time()-t0))
