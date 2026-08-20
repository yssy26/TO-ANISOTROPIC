#!/usr/bin/env python3
"""BFINAL-023 cycle-1: H-CLAMP discriminating gauge (preregistered).

Tests (see PREREGISTRATION.md, same directory):
  T0  reproduction control: b22-lineage instrument (old basis) reproduces
      mobility 12.9% and D1/D2/D3 face anchors 8.68/44.9/8.99% (rebuild
      basis, h=1e-3).
  T1  H-WB: corrected basis (production diffusion sign + symmetric interp
      gamma) rebuilds exact mobility to <1% (verdict rule) + non-clamp
      identity + true clamp set.
  T2  coincidence: D2 top-10% x-normal mismatch faces vs true clamp cells.
  T3  H-CLAMP variants on the face anchor: V0old / V0c / VP / VC / VC1 / VC2
      + W-C channel attribution + A_u sensitivity.
  T4  verdict table numbers.

Zero production changes; offline only, reads only existing exports.
"""
import time, json, os
import numpy as np

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
SA = "/home/ys/dsH/b16_states/1"
SB = "/home/ys/dsH/b22_gauge/stageB2"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-023/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N; P0 = NV
ALPHA_REL = 0.4
NU = 5.19009e-05

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-023 cycle-1: H-CLAMP gauge ===")

# ---------------- mesh (verbatim b22 lineage loaders) ----------------
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
assign_b = ~Ufixed_b
bidx = np.where(assign_b)[0]
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

# ---------------- state B fields ----------------
U_B = np.loadtxt(SB + "/wstate_baseline_U.mtx")
p_B = np.loadtxt(SB + "/wstate_baseline_p.mtx")
phi_B = np.loadtxt(SB + "/wstate_baseline_phi.mtx")
phiB_B = np.loadtxt(SB + "/wstate_baseline_phiB.mtx")
nut_B = np.loadtxt(SB + "/wstate_baseline_nutFrozen.mtx")
alpha_B = np.loadtxt(SB + "/wstate_baseline_alpha.mtx")
mob_B_exact = np.loadtxt(SB + "/wstate_baseline_primalPressureMobility.mtx")
nuEff_B = NU + nut_B
res = {}

# ---------------- bases ----------------
def build_common(nuEff_, alpha_, phi_, phiB_):
    upw_ = phi_ >= 0.0
    qp_ = upw_.astype(float); qn_ = 1.0 - qp_
    divphi_ = np.zeros(N)
    np.add.at(divphi_, owner, phi_); np.add.at(divphi_, nei, -phi_)
    np.add.at(divphi_, bc_cell, phiB_)
    ab_ = nuEff_[bc_cell]*bmag*bdel
    ic_ = np.where(Ufixed_b, -ab_, phiB_)
    return upw_, qp_, qn_, divphi_, ab_, ic_

def build_basis_old(U_, p_, phi_, phiB_, alpha_, nuEff_):
    """b22-lineage instrument: diffusion diag -= kfv, offdiag += kfv,
    one-sided gamma (owner at lower, nei at upper).  VERBATIM semantics."""
    upw_, qp_, qn_, divphi_, ab_, ic_ = build_common(nuEff_, alpha_, phi_, phiB_)
    kfv_o = nuEff_[owner]*magSf_all[:nIF]*deltaCoeffs
    kfv_n = nuEff_[nei]*magSf_all[:nIF]*deltaCoeffs
    lower_ = -qp_*phi_ + kfv_o
    upper_ = qn_*phi_ + kfv_n
    diag_ = np.zeros(N)
    np.add.at(diag_, owner, qp_*phi_)
    np.add.at(diag_, nei, -qn_*phi_)
    diag_ -= divphi_
    np.add.at(diag_, owner, -kfv_o); np.add.at(diag_, nei, -kfv_n)
    diag_ += alpha_*V
    sumOff_ = np.zeros(N)
    np.add.at(sumOff_, owner, np.abs(upper_)); np.add.at(sumOff_, nei, np.abs(lower_))
    Db_ = diag_.copy()
    np.add.at(Db_, bc_cell, np.abs(ic_))
    D_rel_ = np.maximum(np.abs(Db_), sumOff_)/ALPHA_REL
    return dict(lower=lower_, upper=upper_, diag=diag_, D_rel=D_rel_,
                Db=Db_, sumOff=sumOff_, ic=ic_, divphi=divphi_,
                mob_re=V/D_rel_)

def build_basis_corr(U_, p_, phi_, phiB_, alpha_, nuEff_):
    """Corrected: production equation signs (fvm::div - fvm::laplacian):
    upper = qn*phi - kfv_f, lower = -qp*phi - kfv_f (negSumDiag gives
    diag += kfv_f), gamma_f = linear-interpolated (symmetric)."""
    upw_, qp_, qn_, divphi_, ab_, ic_ = build_common(nuEff_, alpha_, phi_, phiB_)
    gam = w*nuEff_[owner] + (1.0 - w)*nuEff_[nei]
    kfv = gam*magSf_all[:nIF]*deltaCoeffs
    lower_ = -qp_*phi_ - kfv
    upper_ = qn_*phi_ - kfv
    diag_ = np.zeros(N)
    np.add.at(diag_, owner, -lower_)     # negSumDiag: diag[own] -= lower
    np.add.at(diag_, nei, -upper_)       # diag[nei] -= upper
    diag_ -= divphi_
    diag_ += alpha_*V
    sumOff_ = np.zeros(N)
    np.add.at(sumOff_, owner, np.abs(upper_)); np.add.at(sumOff_, nei, np.abs(lower_))
    Db_ = diag_.copy()
    np.add.at(Db_, bc_cell, np.abs(ic_))
    D_rel_ = np.maximum(np.abs(Db_), sumOff_)/ALPHA_REL
    return dict(lower=lower_, upper=upper_, diag=diag_, D_rel=D_rel_,
                Db=Db_, sumOff=sumOff_, ic=ic_, divphi=divphi_,
                mob_re=V/D_rel_)

bsB_old = build_basis_old(U_B, p_B, phi_B, phiB_B, alpha_B, nuEff_B)
bsB_cor = build_basis_corr(U_B, p_B, phi_B, phiB_B, alpha_B, nuEff_B)

# ---------------- T0: reproduction control ----------------
r_old = bsB_old["mob_re"]
num = np.linalg.norm(mob_B_exact - r_old); den = np.linalg.norm(mob_B_exact)
res["T0_mobility_old_basis"] = float(num/den)
log("[T0] old-basis mobility relL2=%.4e (b22 read 1.2932e-01)"
    % res["T0_mobility_old_basis"])

# ---------------- T1: H-WB corrected basis ----------------
r_cor = bsB_cor["mob_re"]
clamp_old = bsB_old["sumOff"] > np.abs(bsB_old["Db"])
clamp_cor = bsB_cor["sumOff"] > np.abs(bsB_cor["Db"])
err_cor = np.abs(mob_B_exact - r_cor)/np.maximum(np.abs(mob_B_exact), 1e-300)
res["T1"] = {
    "mobility_corr_relL2": float(np.linalg.norm(mob_B_exact - r_cor)/den),
    "mobility_corr_frac_gt1pct": float(np.mean(err_cor > 0.01)),
    "clamp_frac_old": float(clamp_old.mean()),
    "clamp_frac_corr": float(clamp_cor.mean()),
    # non-clamp identity: mob_exact == alphaRel*V/|Db| exactly on non-clamp
    "nonclamp_identity_relL2": float(
        np.linalg.norm(mob_B_exact[~clamp_cor]
                       - ALPHA_REL*V[~clamp_cor]/np.abs(bsB_cor["Db"][~clamp_cor]))
        /np.linalg.norm(mob_B_exact[~clamp_cor])),
    # clamp self-check: mob_exact*sumOff/(alphaRel*V) ~ 1 on clamp cells
    "clamp_selfcheck_median": float(np.median(
        mob_B_exact[clamp_cor]*bsB_cor["sumOff"][clamp_cor]
        /(ALPHA_REL*V[clamp_cor]))),
    # production mobility overshoot factor on clamp cells: (alphaRel/A_u)/rAtU
    "A_rel_over_A_u_clamp_median": float(np.median(
        bsB_cor["sumOff"][clamp_cor]/np.abs(bsB_cor["Db"][clamp_cor]))
        if clamp_cor.any() else float("nan")),
    "A_rel_over_A_u_clamp_p90": float(np.percentile(
        bsB_cor["sumOff"][clamp_cor]/np.abs(bsB_cor["Db"][clamp_cor]), 90)
        if clamp_cor.any() else float("nan")),
}
nut_active = nut_B > 1e-6
res["T1"]["clamp_nut_active_frac"] = float(
    np.mean(clamp_cor & nut_active)/max(np.mean(clamp_cor), 1e-300))
# margin structure of the corrected basis: how close is |D| to sumOff?
dom_margin = np.abs(bsB_cor["Db"])/np.maximum(bsB_cor["sumOff"], 1e-300)
res["T1"]["dominance_margin_pctiles"] = [
    float(np.percentile(dom_margin, q)) for q in (0.1, 1, 5, 50)]
res["T1"]["dominance_margin_min_nutactive"] = float(
    dom_margin[nut_active].min() if nut_active.any() else float("nan"))
log("[T1] corrected-basis mobility relL2=%.4e frac>1%%=%.4f"
    % (res["T1"]["mobility_corr_relL2"], res["T1"]["mobility_corr_frac_gt1pct"]))
log("    clamp frac: old %.3f -> corr %.3f ; non-clamp identity relL2=%.3e"
    % (res["T1"]["clamp_frac_old"], res["T1"]["clamp_frac_corr"],
       res["T1"]["nonclamp_identity_relL2"]))
log("    clamp self-check median=%.4f ; A_rel/A_u on clamp: med=%.2f p90=%.2f"
    % (res["T1"]["clamp_selfcheck_median"],
       res["T1"]["A_rel_over_A_u_clamp_median"],
       res["T1"]["A_rel_over_A_u_clamp_p90"]))

# where do the residual errors sit (structure check, corrected basis)
for nm, msk in (("solid(a>1e6)", alpha_B > 1e6),
                ("transition", (alpha_B > 1e2) & (alpha_B <= 1e6)),
                ("fluid(a<=1e2)", alpha_B <= 1e2),
                ("clamp_corr", clamp_cor), ("nonclamp", ~clamp_cor),
                ("nut_active", nut_active)):
    if msk.sum() == 0: continue
    log("    corr-basis err %-14s n=%5d mean=%.3e frac>1%%=%.3f"
        % (nm, msk.sum(), err_cor[msk].mean(), np.mean(err_cor[msk] > 0.01)))

# ---------------- tangent machinery ----------------
def H_field(bs, mob, U_):
    """H = [(D_rel - D0)*U - off*U]/V  (HbyA = mob*H)."""
    src = (bs["D_rel"] - bs["diag"])[:, None]*U_
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, bs["upper"][:, None]*U_[nei])
    np.add.at(offU, nei, bs["lower"][:, None]*U_[owner])
    return (src - offU)/V[:, None]

def channels(bs, mob, U_, p_, dU, dp, da, clamp,
             fac_mode="clamped", dU_coef_mode="plain", da_mode="old"):
    """dphi channels on basis bs with mobility field `mob`.
    fac_mode: 'clamped' -> dH factor = mob (true rAtU);
              'unclamped' -> alphaRel/A_u (production prodRAU*alphaRel).
    dU_coef_mode: 'plain' -> (1-alphaRel); 'clamped' -> 1 - mob*A_u.
    da_mode: 'old' -> drAU=-mob^2/alphaRel*da, dHsrc=(1/alphaRel-1)*da*U;
             'clamped' -> drAU=0 & dHsrc=-da*U on clamp cells."""
    A_u = np.abs(bs["Db"])/V
    mobF = w*mob[owner] + (1.0 - w)*mob[nei]
    kf = mobF*deltaCoeffs*magSf_all[:nIF]
    # --- dp channel ---
    dphi_p = -kf*(dp[nei] - dp[owner])
    dphi_p_bnd = np.zeros(nBnd)
    for bf in bidx:
        c = bc_cell[bf]
        dphi_p_bnd[bf] = mob[c]*bdel[bf]*bmag[bf]*dp[c]
    # --- dU channel ---
    offdU = np.zeros((N, 3))
    np.add.at(offdU, owner, bs["upper"][:, None]*dU[nei])
    np.add.at(offdU, nei, bs["lower"][:, None]*dU[owner])
    dH = -offdU/V[:, None]
    if fac_mode == "clamped":   fac = mob
    else:                       fac = ALPHA_REL/A_u
    if dU_coef_mode == "plain": ccoef = np.full(N, 1.0 - ALPHA_REL)
    else:                       ccoef = 1.0 - mob*A_u
    dHbyA = fac[:, None]*dH + ccoef[:, None]*dU
    dphi_u = np.einsum("fi,fi->f", sf, w[:, None]*dHbyA[owner]
                       + (1-w)[:, None]*dHbyA[nei])
    dphi_u_bnd = np.zeros(nBnd)
    for bf in bidx:
        c = bc_cell[bf]
        dphi_u_bnd[bf] = bSf[bf] @ dHbyA[c]
    # --- da channel ---
    if da_mode == "old":
        drAU = -(mob**2/ALPHA_REL)*da
        dHsrc = ((1.0/ALPHA_REL - 1.0)*da)[:, None]*U_
    else:
        drAU = -(mob**2/ALPHA_REL)*da*(~clamp)
        dHsrc = (((1.0/ALPHA_REL - 1.0)*(~clamp) + (-1.0)*clamp)*da)[:, None]*U_
    Hcur = H_field(bs, mob, U_)
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

def mets(a, b):
    return {"cos": float(a @ b/(np.linalg.norm(a)*np.linalg.norm(b))),
            "relL2": float(np.linalg.norm(a - b)/np.linalg.norm(b))}

dA = read_of_internal(SA + "/dAlphaDxh")
zdirs = np.loadtxt("/home/ys/dsH/b15_export/stageB6_rxc_z_analytic.mtx").reshape(3, N)

VARIANTS = {
    "V0old": dict(bs=bsB_old, mob=None,        fac="clamped",   dU="plain",    da="old"),
    "V0c":   dict(bs=bsB_cor, mob=mob_B_exact, fac="clamped",   dU="plain",    da="old"),
    "VP":    dict(bs=bsB_cor, mob=mob_B_exact, fac="unclamped", dU="plain",    da="old"),
    "VC1":   dict(bs=bsB_cor, mob=mob_B_exact, fac="clamped",   dU="clamped",  da="old"),
    "VC2":   dict(bs=bsB_cor, mob=mob_B_exact, fac="clamped",   dU="plain",    da="clamped"),
    "VC":    dict(bs=bsB_cor, mob=mob_B_exact, fac="clamped",   dU="clamped",  da="clamped"),
}

# ---------------- T3: anchors ----------------
ladder = {"D1": ["0.001", "0.0003"], "D2": ["0.001", "0.0003"],
          "D3": ["0.001", "0.0003"]}
res["T3"] = {}
d2_sets = {}
log("\n[T3] face anchors (internal+assignable-bnd concat):")
for di, nm in enumerate(["D1", "D2", "D3"]):
    res["T3"][nm] = {}
    for tag in ladder[nm]:
        h = float(tag)
        pre = SB + "/wstate_%s_h%s_" % (nm, tag)
        U_p = np.loadtxt(pre + "p_U.mtx"); p_p = np.loadtxt(pre + "p_p.mtx")
        U_m = np.loadtxt(pre + "m_U.mtx"); p_m = np.loadtxt(pre + "m_p.mtx")
        phi_p = np.loadtxt(pre + "p_phi.mtx"); phi_m = np.loadtxt(pre + "m_phi.mtx")
        phiBp = np.loadtxt(pre + "p_phiB.mtx"); phiBm = np.loadtxt(pre + "m_phiB.mtx")
        dphi_true = np.concatenate([(phi_p - phi_m)/(2*h),
                                    (phiBp - phiBm)/(2*h)])
        dU = (U_p - U_m)/(2*h); dp_ = (p_p - p_m)/(2*h)
        da_analytic = dA*zdirs[di]
        alpha_p = np.loadtxt(pre + "p_alpha.mtx")
        alpha_m = np.loadtxt(pre + "m_alpha.mtx")
        da_export = (alpha_p - alpha_m)/(2*h)
        da_xcheck = float(np.linalg.norm(da_export - da_analytic)
                          /np.linalg.norm(da_analytic))
        row = {"da_export_vs_analytic": da_xcheck}
        for vn, cfg in VARIANTS.items():
            mob = cfg["mob"]
            if mob is None: mob = cfg["bs"]["mob_re"]   # V0old self-consistent
            duE, dpE, daE = channels(cfg["bs"], mob, U_B, p_B, dU, dp_,
                                     da_export, clamp_cor,
                                     fac_mode=cfg["fac"], dU_coef_mode=cfg["dU"],
                                     da_mode=cfg["da"])
            tot = duE + dpE + daE
            row[vn] = mets(tot, dphi_true)
            row[vn+"_chan"] = {
                "du": float(np.linalg.norm(duE)/np.linalg.norm(tot)),
                "dp": float(np.linalg.norm(dpE)/np.linalg.norm(tot)),
                "da": float(np.linalg.norm(daE)/np.linalg.norm(tot))}
            if nm == "D2" and tag == "0.001" and vn == "V0old":
                d2_sets["V0old"] = (tot - dphi_true, dphi_true, tot)
            if nm == "D2" and tag == "0.001":
                d2_sets[vn] = (tot - dphi_true, dphi_true, tot)
                d2_sets[vn+"_chan"] = (duE, dpE, daE)
        res["T3"][nm][tag] = row
        log("  %s h=%s: da xcheck=%.2e" % (nm, tag, da_xcheck))
        for vn in VARIANTS:
            m = row[vn]
            log("    %-6s cos=%.4f relL2=%.4f  (chan du/dp/da %.2f/%.2f/%.2f)"
                % (vn, m["cos"], m["relL2"], row[vn+"_chan"]["du"],
                   row[vn+"_chan"]["dp"], row[vn+"_chan"]["da"]))

# ---------------- T2: D2 carrier set vs clamp ----------------
mismF, truthF, totF = d2_sets["V0old"]
mismF = mismF[:nIF]; truthF = truthF[:nIF]
sfhat = sf/np.maximum(np.linalg.norm(sf, axis=1, keepdims=True), 1e-300)
absx = np.abs(sfhat[:, 0])
energy = mismF**2
top_thr = np.sort(energy)[-max(nIF//10, 1)]
carrier = (energy >= top_thr) & (absx > 0.9)
clamp_adj = clamp_cor[owner] | clamp_cor[nei]
res["T2"] = {
    "carrier_faces": int(carrier.sum()),
    "carrier_frac_of_mismatch_energy": float(
        energy[carrier].sum()/energy.sum()),
    "carrier_frac_clamp_adjacent": float(clamp_adj[carrier].mean()),
    "allfaces_frac_clamp_adjacent": float(clamp_adj.mean()),
    # relL2 split by clamp adjacency (V0old and VC)
}
for vn in ["V0old", "V0c", "VP", "VC1", "VC2", "VC"]:
    mm, tt, _ = d2_sets[vn]
    mmF = mm[:nIF]; ttF = tt[:nIF]
    res["T2"][vn] = {
        "relL2_carrier": float(np.linalg.norm(mmF[carrier])
                               /np.linalg.norm(ttF[carrier])),
        "relL2_clampadj": float(np.linalg.norm(mmF[clamp_adj])
                                /np.linalg.norm(ttF[clamp_adj])),
        "relL2_nonclampadj": float(np.linalg.norm(mmF[~clamp_adj])
                                   /np.linalg.norm(ttF[~clamp_adj])),
        "cos_carrier": float(mmF[carrier] @ -ttF[carrier]
            /(np.linalg.norm(mmF[carrier])*np.linalg.norm(ttF[carrier]))*-1.0),
    }
log("\n[T2] D2 carrier (top-10%% energy, x-normal): %d faces, %.1f%% of "
    "mismatch energy; %.3f clamp-adjacent (all faces: %.3f)"
    % (res["T2"]["carrier_faces"],
       100*res["T2"]["carrier_frac_of_mismatch_energy"],
       res["T2"]["carrier_frac_clamp_adjacent"],
       res["T2"]["allfaces_frac_clamp_adjacent"]))
for vn in ["V0old", "V0c", "VP", "VC1", "VC2", "VC"]:
    r = res["T2"][vn]
    log("    %-6s carrier relL2=%.4f ; clampadj=%.4f ; nonclamp=%.4f"
        % (vn, r["relL2_carrier"], r["relL2_clampadj"], r["relL2_nonclampadj"]))

# ---------------- W-C: channel attribution on the carrier ----------------
wc = {}
duO, dpO, daO = d2_sets["V0old_chan"]
duC, dpC, daC = d2_sets["VC_chan"]
fix_tot = np.linalg.norm((d2_sets["V0old"][2] - d2_sets["VC"][2])[:nIF][carrier])
wc["carrier_fix_L2"] = float(fix_tot)
wc["fix_from_dUchan"] = float(np.linalg.norm((duO - duC)[:nIF][carrier])/max(fix_tot, 1e-300))
wc["fix_from_dpchan"] = float(np.linalg.norm((dpO - dpC)[:nIF][carrier])/max(fix_tot, 1e-300))
wc["fix_from_dachan"] = float(np.linalg.norm((daO - daC)[:nIF][carrier])/max(fix_tot, 1e-300))
# channel magnitudes on the carrier (VC basis)
wc["VC_chan_shares_carrier"] = {
    "du": float(np.linalg.norm(duC[:nIF][carrier])/np.linalg.norm(duC[:nIF][carrier] + dpC[:nIF][carrier] + daC[:nIF][carrier])),
    "dp": float(np.linalg.norm(dpC[:nIF][carrier])/np.linalg.norm(duC[:nIF][carrier] + dpC[:nIF][carrier] + daC[:nIF][carrier])),
    "da": float(np.linalg.norm(daC[:nIF][carrier])/np.linalg.norm(duC[:nIF][carrier] + dpC[:nIF][carrier] + daC[:nIF][carrier])),
}
res["WC"] = wc
log("\n[WC] carrier fix decomposition |V0old-VC|: dU-chan %.2f, dp-chan %.2f, "
    "da-chan %.2f (of fix L2)" % (wc["fix_from_dUchan"], wc["fix_from_dpchan"],
                                  wc["fix_from_dachan"]))
log("    VC channel shares on carrier: du/dp/da = %.2f/%.2f/%.2f"
    % (wc["VC_chan_shares_carrier"]["du"], wc["VC_chan_shares_carrier"]["dp"],
       wc["VC_chan_shares_carrier"]["da"]))

# ---------------- A_u sensitivity for VC ----------------
bsS = dict(bsB_cor)
for sgn, nm in ((1.15, "plus15"), (0.85, "minus15")):
    bsS = dict(bsB_cor)
    bsS["Db"] = bsB_cor["Db"]*sgn
    bsS["D_rel"] = bsB_cor["D_rel"]*sgn
    pre = SB + "/wstate_D2_h0.001_"
    U_p = np.loadtxt(pre + "p_U.mtx"); p_p = np.loadtxt(pre + "p_p.mtx")
    U_m = np.loadtxt(pre + "m_U.mtx"); p_m = np.loadtxt(pre + "m_p.mtx")
    phi_p = np.loadtxt(pre + "p_phi.mtx"); phi_m = np.loadtxt(pre + "m_phi.mtx")
    phiBp = np.loadtxt(pre + "p_phiB.mtx"); phiBm = np.loadtxt(pre + "m_phiB.mtx")
    dphi_true = np.concatenate([(phi_p - phi_m)/2e-3, (phiBp - phiBm)/2e-3])
    dU = (U_p - U_m)/2e-3; dp_ = (p_p - p_m)/2e-3
    alpha_p = np.loadtxt(pre + "p_alpha.mtx"); alpha_m = np.loadtxt(pre + "m_alpha.mtx")
    da_ex = (alpha_p - alpha_m)/2e-3
    duE, dpE, daE = channels(bsS, mob_B_exact, U_B, p_B, dU, dp_, da_ex,
                             clamp_cor, fac_mode="clamped",
                             dU_coef_mode="clamped", da_mode="clamped")
    m = mets(duE + dpE + daE, dphi_true)
    res["T3"]["D2"]["0.001"].setdefault("VC_sens", {})[nm] = m
    log("[sens] VC with A_u %s: D2 cos=%.4f relL2=%.4f"
        % (nm, m["cos"], m["relL2"]))

json.dump(res, open(OUT + "/b23_hclamp_gauge.json", "w"), indent=1)
log("\nDONE t=%.1fs" % (time.time()-t0))
