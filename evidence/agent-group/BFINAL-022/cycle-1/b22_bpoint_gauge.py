#!/usr/bin/env python3
"""BFINAL-022 W4.3/W5: B-point consistent offline gauge.

All four ingredients at the SAME state B (the stageB2 full-SST frozen
baseline):
  operator : offline flux-map rebuild at B, with the kf/rAU channels taken
             from the EXACT exported primalPressureMobility(B) (the W1
             refresh) -- removes the rAU-rebuild uncertainty that limited
             the B21 instrument (3.35% internal fidelity, rebuild basis);
  rxd      : alpha-channel reconstruction on the same exact basis;
  rhs/FD   : the exported wstate_* probe states (FD truth) and b18rhs_*.

Tests:
  0. exact-vs-rebuild mobility: per-cell relL2 of the exported
     primalPressureMobility(B) vs the offline D_rel rebuild at B
     (quantifies the rAU-rebuild instrument error directly, and the
     SLOT-2 magnitude with the true A-capture vs true B-capture).
  1. flux-rebuild fidelity at B: exact-kf basis vs rebuild basis
     (does the exact kf push the internal fidelity below the 3.35%
     floor?).
  2. face-level flux-response anchor at B (exact basis): D1/D2/D3
     total-tangent vs dphi_state-FD  (B21 read 8.7/45/9.0% with the
     rebuild basis).
  3. P-row identity at B (exact basis, E2b form):
     |div(UV channels) + rxd_P| / |rxd_P|.
  4. (W5) D2 localization: spatial structure of the D2 face-level
     mismatch (by cell coordinate, alpha value, design mask).
"""
import time, json, os
import numpy as np

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
SA = "/home/ys/dsH/b16_states/1"          # state A fields (mobility compare)
SB = "/home/ys/dsH/b22_gauge/stageB2"     # state B exports (this round)
OLD = "/home/ys/dsH/b21_anchor/stageB2"   # pre-fix B exports (control)
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-022/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N; P0 = NV
ALPHA_REL = 0.4
NU = 5.19009e-05

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-022 W4.3/W5: B-point consistent gauge ===")

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
nFaceTot = nIF + nBnd
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
import re
def read_boundary_vals(path):
    txt = open(path).read()
    resd = {}
    for m in re.finditer(r"\b(inlet|outlet|hotInlet|hotOutlet|solidEndWalls|"
                         r"bottomWall|topWall|sideWalls)\s*\n?\s*\{", txt):
        name = m.group(1); i = m.end(); depth = 1; j = i
        while depth > 0:
            if txt[j] == "{": depth += 1
            elif txt[j] == "}": depth -= 1
            j += 1
        blk = txt[i:j]
        vm = re.search(r"value\s+nonuniform[^0-9\-]*\d+\s*\(", blk)
        if vm:
            k = blk.index("(", vm.end()-1); kk = blk.index(")", k)
            vals = []
            for tok in blk[k+1:kk].replace("(", " ").replace(")", " ").split():
                try: vals.append(float(tok))
                except ValueError: pass
            resd[name] = np.array(vals); continue
        vm = re.search(r"value\s+uniform\s+([^\n;]+)", blk)
        if vm:
            u = vm.group(1).strip().rstrip(";").strip()
            if u.startswith("("):
                resd[name] = np.array([float(x) for x in u.strip("()").split()])
            else:
                resd[name] = np.array([float(u)])
        else:
            resd[name] = None
    return resd

# ---------------- state B fields (this round's run) ----------------
U_B = np.loadtxt(SB + "/wstate_baseline_U.mtx")
p_B = np.loadtxt(SB + "/wstate_baseline_p.mtx")
phi_B = np.loadtxt(SB + "/wstate_baseline_phi.mtx")
phiB_B = np.loadtxt(SB + "/wstate_baseline_phiB.mtx")
nut_B = np.loadtxt(SB + "/wstate_baseline_nutFrozen.mtx")
nuEff_B_raw = np.loadtxt(SB + "/wstate_baseline_nuEffFrozen.mtx")
alpha_B = np.loadtxt(SB + "/wstate_baseline_alpha.mtx")
mob_B_exact = np.loadtxt(SB + "/wstate_baseline_primalPressureMobility.mtx")
nuEff_B = NU + nut_B
res = {}

# ---------------- bases ----------------
def build_basis(U_, p_, phi_, phiB_, alpha_, nuEff_, mob_override=None):
    upw_ = phi_ >= 0.0
    qp_ = upw_.astype(float); qn_ = 1.0 - qp_
    kfv_o = nuEff_[owner]*magSf_all[:nIF]*deltaCoeffs
    kfv_n = nuEff_[nei]*magSf_all[:nIF]*deltaCoeffs
    lower_ = -qp_*phi_ + kfv_o
    upper_ = qn_*phi_ + kfv_n
    diag_ = np.zeros(N)
    np.add.at(diag_, owner, qp_*phi_)
    np.add.at(diag_, nei, -qn_*phi_)
    divphi_ = np.zeros(N)
    np.add.at(divphi_, owner, phi_); np.add.at(divphi_, nei, -phi_)
    np.add.at(divphi_, bc_cell, phiB_)
    diag_ -= divphi_
    np.add.at(diag_, owner, -kfv_o); np.add.at(diag_, nei, -kfv_n)
    diag_ += alpha_*V
    ab_ = nuEff_[bc_cell]*bmag*bdel
    ic_ = np.where(Ufixed_b, -ab_, phiB_)
    sumOff_ = np.zeros(N)
    np.add.at(sumOff_, owner, np.abs(upper_))
    np.add.at(sumOff_, nei, np.abs(lower_))
    Db_ = diag_.copy()
    np.add.at(Db_, bc_cell, np.abs(ic_))
    D_rel_ = np.maximum(np.abs(Db_), sumOff_)/ALPHA_REL
    rAU_rel_re = 1.0/(D_rel_/V)
    if mob_override is None:
        rAU_rel_ = rAU_rel_re
        mob_src = "rebuild"
    else:
        rAU_rel_ = mob_override
        mob_src = "exact"
    rAU_u_ = rAU_rel_/ALPHA_REL
    mobF_ = w*rAU_rel_[owner] + (1.0-w)*rAU_rel_[nei]
    kf_ = mobF_*deltaCoeffs*magSf_all[:nIF]
    return dict(upw=upw_, qp=qp_, qn=qn_, lower=lower_, upper=upper_,
                diag=diag_, D_rel=D_rel_, rAU_rel=rAU_rel_, rAU_u=rAU_u_,
                mobF=mobF_, kf=kf_, divphi=divphi_,
                rAU_rel_rebuild=rAU_rel_re, mob_src=mob_src)

def HbyA_field(U_, bs):
    src = (bs["D_rel"] - bs["diag"])[:, None]*U_
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, bs["upper"][:, None]*U_[nei])
    np.add.at(offU, nei, bs["lower"][:, None]*U_[owner])
    H = (src - offU)/V[:, None]
    return bs["rAU_rel"][:, None]*H

def phi_rebuild(U_, p_, bs):
    HbyA = HbyA_field(U_, bs)
    phi_int = np.einsum("fi,fi->f", sf, w[:, None]*HbyA[owner]
                        + (1-w)[:, None]*HbyA[nei])
    phi_int -= bs["kf"]*(p_[nei] - p_[owner])
    phi_bnd = np.zeros(nBnd)
    for bf in bidx:
        c = bc_cell[bf]
        phi_bnd[bf] = bSf[bf] @ HbyA[c] + bs["rAU_rel"][c]*bdel[bf]*bmag[bf]*p_[c]
    return np.concatenate([phi_int, phi_bnd])

# ---------------- 0. exact vs rebuild mobility ----------------
bsB_re = build_basis(U_B, p_B, phi_B, phiB_B, alpha_B, nuEff_B)
rAU_re = bsB_re["rAU_rel_rebuild"]
num = np.linalg.norm(mob_B_exact - rAU_re)
den = np.linalg.norm(mob_B_exact)
res["mobility_exact_vs_rebuild_B"] = {
    "relL2": float(num/den),
    "maxRel": float(np.max(np.abs(mob_B_exact - rAU_re)
                           /np.maximum(np.abs(mob_B_exact), 1e-300))),
    "frac_out_1pct": float(np.mean(
        np.abs(mob_B_exact - rAU_re)/np.maximum(np.abs(mob_B_exact), 1e-300)
        > 0.01)),
}
log("[0] exact vs rebuild mobility @B: relL2=%.4e maxRel=%.3e frac>1%%=%.3f"
    % (res["mobility_exact_vs_rebuild_B"]["relL2"],
       res["mobility_exact_vs_rebuild_B"]["maxRel"],
       res["mobility_exact_vs_rebuild_B"]["frac_out_1pct"]))

# SLOT-2 magnitude with the TRUE captures: A-capture (main-loop, molecular
# system, = what the module adjoint used before W1) vs B-capture (refresh).
try:
    # A-state true rAU: rebuild from A fields (molecular nuEff) -- the b21
    # instrument; report as reference only.
    def read_A():
        U_A = read_of_internal(SA + "/U").reshape(N, 3)
        p_A = read_of_internal(SA + "/p")
        phi_A = read_of_internal(SA + "/phi")
        alpha_A = read_of_internal(SA + "/alpha")
        phi_bv = read_boundary_vals(SA + "/phi")
        phiB_A = np.zeros(nBnd); o = 0
        for pn, nf in [("inlet", 84), ("outlet", 84), ("hotInlet", 196),
                       ("hotOutlet", 196), ("solidEndWalls", 280),
                       ("bottomWall", 1120), ("topWall", 1120),
                       ("sideWalls", 4800)]:
            v = phi_bv[pn]
            phiB_A[o:o+nf] = (v[0] if v.size == 1 else v)
            o += nf
        return U_A, p_A, phi_A, phiB_A, alpha_A
    U_A, p_A, phi_A, phiB_A, alpha_A = read_A()
    bsA = build_basis(U_A, p_A, phi_A, phiB_A, alpha_A, np.full(N, NU))
    ratio = np.where(bsA["rAU_rel"] > 0, mob_B_exact/bsA["rAU_rel"], 1.0)
    res["mobility_trueB_vs_Arebuild"] = {
        "median": float(np.median(ratio)), "p10": float(np.percentile(ratio, 10)),
        "p90": float(np.percentile(ratio, 90)), "min": float(ratio.min()),
        "max": float(ratio.max()),
        "frac_out_0.9_1.1": float(np.mean((ratio < 0.9) | (ratio > 1.1)))}
    log("    true B-capture vs A-rebuild: median=%.4f frac_out_10%%=%.3f "
        "min=%.3f max=%.3f"
        % (res["mobility_trueB_vs_Arebuild"]["median"],
           res["mobility_trueB_vs_Arebuild"]["frac_out_0.9_1.1"],
           res["mobility_trueB_vs_Arebuild"]["min"],
           res["mobility_trueB_vs_Arebuild"]["max"]))
except Exception as e:
    log("    A-compare failed: %s" % e)

# ---------------- 1. flux-rebuild fidelity at B ----------------
bsB_ex = build_basis(U_B, p_B, phi_B, phiB_B, alpha_B, nuEff_B,
                     mob_override=mob_B_exact)
def fidelity(bs):
    re = phi_rebuild(U_B, p_B, bs)
    errI = re[:nIF] - phi_B
    errB = re[nIF:][bidx] - phiB_B[bidx]
    return (float(np.linalg.norm(errI)/np.linalg.norm(phi_B)),
            float(np.linalg.norm(errB)/max(np.linalg.norm(phiB_B[bidx]), 1e-300)))
fi_re, fo_re = fidelity(bsB_re)
fi_ex, fo_ex = fidelity(bsB_ex)
res["rebuild_fidelity_B"] = {"rebuild_basis": {"internal": fi_re, "outlet": fo_re},
                             "exact_kf_basis": {"internal": fi_ex, "outlet": fo_ex}}
log("[1] flux-rebuild fidelity @B: rebuild %.4e internal / %.3e outlet"
    " -> exact-kf %.4e / %.3e" % (fi_re, fo_re, fi_ex, fo_ex))

# ---------------- channels ----------------
dA = read_of_internal(SA + "/dAlphaDxh")
zdirs = np.loadtxt("/home/ys/dsH/b15_export/stageB6_rxc_z_analytic.mtx").reshape(3, N)

def psi_channels(bs, U_, p_, dU, dp, da):
    dphi_p = -bs["kf"]*(dp[nei] - dp[owner])
    dphi_p_bnd = np.zeros(nBnd)
    for bf in bidx:
        c = bc_cell[bf]
        dphi_p_bnd[bf] = bs["rAU_rel"][c]*bdel[bf]*bmag[bf]*dp[c]
    offdU = np.zeros((N, 3))
    np.add.at(offdU, owner, bs["upper"][:, None]*dU[nei])
    np.add.at(offdU, nei, bs["lower"][:, None]*dU[owner])
    dH = -offdU/V[:, None]
    dHbyA = ALPHA_REL*bs["rAU_u"][:, None]*dH + (1.0-ALPHA_REL)*dU
    dphi_u = np.einsum("fi,fi->f", sf, w[:, None]*dHbyA[owner]
                       + (1-w)[:, None]*dHbyA[nei])
    dphi_u_bnd = np.zeros(nBnd)
    for bf in bidx:
        c = bc_cell[bf]
        dphi_u_bnd[bf] = bSf[bf] @ dHbyA[c]
    drAU = -bs["rAU_rel"]**2/ALPHA_REL*da
    Hcur = HbyA_field(U_, bs)/bs["rAU_rel"][:, None]
    dHsrc = (1.0/ALPHA_REL - 1.0)*da[:, None]*U_
    dHbyA_a = drAU[:, None]*Hcur + bs["rAU_rel"][:, None]*dHsrc
    dphi_a = np.einsum("fi,fi->f", sf, w[:, None]*dHbyA_a[owner]
                       + (1-w)[:, None]*dHbyA_a[nei])
    dphi_a += (drAU[owner]*w + drAU[nei]*(1-w))*deltaCoeffs*magSf_all[:nIF]*(
        p_[nei] - p_[owner])
    dphi_a_bnd = np.zeros(nBnd)
    for bf in bidx:
        c = bc_cell[bf]
        dphi_a_bnd[bf] = bSf[bf] @ dHbyA_a[c] \
            + drAU[c]*bdel[bf]*bmag[bf]*p_[c]
    return (np.concatenate([dphi_u, dphi_u_bnd]),
            np.concatenate([dphi_p, dphi_p_bnd]),
            np.concatenate([dphi_a, dphi_a_bnd]))

def divF(x):
    d = np.zeros(N)
    np.add.at(d, owner, x[:nIF]); np.add.at(d, nei, -x[:nIF])
    np.add.at(d, bc_cell, x[nIF:])
    return d

ladder = {"D1": ["0.0003", "0.001", "0.003"], "D2": ["0.0003", "0.001"],
          "D3": ["0.0003", "0.001"]}
res["anchor"] = {}
log("\n[2] face-level anchor at B (exact-kf basis) + [3] P-row identity:")
d2_mismatch = None
for di, nm in enumerate(["D1", "D2", "D3"]):
    res["anchor"][nm] = {}
    for tag in ladder[nm]:
        h = float(tag)
        pre = SB + "/wstate_%s_h%s_" % (nm, tag)
        if not os.path.exists(pre + "p_phi.mtx"):
            log("  %s h=%s: phi exports missing, skip" % (nm, tag))
            continue
        U_p = np.loadtxt(pre + "p_U.mtx"); p_p = np.loadtxt(pre + "p_p.mtx")
        U_m = np.loadtxt(pre + "m_U.mtx"); p_m = np.loadtxt(pre + "m_p.mtx")
        phi_p = np.loadtxt(pre + "p_phi.mtx")
        phi_m = np.loadtxt(pre + "m_phi.mtx")
        phiBp = np.loadtxt(pre + "p_phiB.mtx")
        phiBm = np.loadtxt(pre + "m_phiB.mtx")
        dphi_true = np.concatenate([(phi_p - phi_m)/(2*h),
                                    (phiBp - phiBm)/(2*h)])
        dU = (U_p - U_m)/(2*h); dp_ = (p_p - p_m)/(2*h)
        da = dA*zdirs[di]
        duE, dpE, daE = psi_channels(bsB_ex, U_B, p_B, dU, dp_, da)
        totE = duE + dpE + daE
        duR, dpR, daR = psi_channels(bsB_re, U_B, p_B, dU, dp_, da)
        totR = duR + dpR + daR
        def mets(a, b):
            return {"cos": float(a @ b/(np.linalg.norm(a)*np.linalg.norm(b))),
                    "relL2": float(np.linalg.norm(a - b)/np.linalg.norm(b))}
        # P-row identity (E2b form) on the exact basis:
        rxdP_re = divF(daE)
        rP = divF(duE + dpE) + rxdP_re
        row = {
            "exact_total": mets(totE, dphi_true),
            "rebuild_total": mets(totR, dphi_true),
            "identP_exact": {
                "|rP|/|rxdP|": float(np.linalg.norm(rP)/np.linalg.norm(rxdP_re)),
                "cos": float(divF(duE + dpE) @ (-rxdP_re)
                    /(np.linalg.norm(divF(duE + dpE))*np.linalg.norm(rxdP_re)))},
        }
        res["anchor"][nm][tag] = row
        log("  %s h=%s: exact cos=%.4f rel=%.4f (rebuild rel=%.4f) ; "
            "identP |rP|=%.3f cos=%.4f"
            % (nm, tag, row["exact_total"]["cos"], row["exact_total"]["relL2"],
               row["rebuild_total"]["relL2"],
               row["identP_exact"]["|rP|/|rxdP|"],
               row["identP_exact"]["cos"]))
        if nm == "D2" and tag == "0.001":
            # W5 localization uses the SELF-CONSISTENT rebuild basis (the
            # exact-kf hybrid mixes two rAU conventions and is only a
            # bracket; the rebuild basis is the calibrated instrument).
            d2_mismatch = (totR - dphi_true, dphi_true, totR)

# ---------------- 4. W5: D2 mismatch spatial structure ----------------
if d2_mismatch is not None:
    mism, truth, tot = d2_mismatch
    mismF = mism[:nIF]; truthF = truth[:nIF]
    res["D2_localization"] = {}
    # face -> cell coordinate bins (owner cell)
    xo = Ccen[owner]
    frac = np.abs(mismF)/max(np.linalg.norm(mismF)/np.sqrt(nIF), 1e-300)
    # by x-decile
    xdec = np.clip((xo[:, 0]/np.ptp(pts[:, 0])*10).astype(int), 0, 9)
    byx = []
    for d in range(10):
        msk = xdec == d
        byx.append(float(np.linalg.norm(mismF[msk])/max(np.linalg.norm(truthF[msk]), 1e-300)))
    res["D2_localization"]["relL2_by_x_decile"] = byx
    # by solid(α) vs fluid region of owner cell
    a_o = alpha_B[owner]
    solid = a_o > 1e3
    fluid = ~solid
    res["D2_localization"]["relL2_solid_faces"] = float(
        np.linalg.norm(mismF[solid])/max(np.linalg.norm(truthF[solid]), 1e-300))
    res["D2_localization"]["relL2_fluid_faces"] = float(
        np.linalg.norm(mismF[fluid])/max(np.linalg.norm(truthF[fluid]), 1e-300))
    res["D2_localization"]["frac_energy_top10pct_faces"] = float(
        np.sum(np.sort(mismF**2)[-max(nIF//10, 1):])/np.sum(mismF**2))
    # correlation of the mismatch with the D2 direction support
    dz = zdirs[1][owner]
    res["D2_localization"]["corr_mismatch_vs_zdir_owner"] = float(
        np.corrcoef(mismF, dz)[0, 1])
    # directional decomposition of the mismatch: face-normal orientation
    sfhat = sf/np.maximum(np.linalg.norm(sf, axis=1, keepdims=True), 1e-300)
    absx = np.abs(sfhat[:, 0]); absy = np.abs(sfhat[:, 1]); absz = np.abs(sfhat[:, 2])
    emsk = mismF**2
    res["D2_localization"]["mismatch_energy_by_normal"] = {
        "x-norm_faces_frac": float(emsk[absx > 0.9].sum()/emsk.sum()),
        "y-norm_faces_frac": float(emsk[absy > 0.9].sum()/emsk.sum()),
        "z-norm_faces_frac": float(emsk[absz > 0.9].sum()/emsk.sum()),
        "mixed_frac": float(emsk[(absx <= 0.9)*(absy <= 0.9)*(absz <= 0.9)].sum()/emsk.sum()),
    }
    # where is the truth response itself large?
    tmsk = truthF**2
    res["D2_localization"]["truth_energy_by_normal"] = {
        "x": float(tmsk[absx > 0.9].sum()/tmsk.sum()),
        "y": float(tmsk[absy > 0.9].sum()/tmsk.sum()),
        "z": float(tmsk[absz > 0.9].sum()/tmsk.sum()),
    }
    log("\n[4] D2 localization: relL2 by x-decile=%s" % ["%.2f" % v for v in byx])
    log("    solid(owner α>1e3) faces: %.3f ; fluid faces: %.3f ; "
        "top-10%% faces carry %.1f%% of mismatch energy"
        % (res["D2_localization"]["relL2_solid_faces"],
           res["D2_localization"]["relL2_fluid_faces"],
           100*res["D2_localization"]["frac_energy_top10pct_faces"]))
    log("    mismatch energy by face normal (x/y/z/mixed): %s"
        % res["D2_localization"]["mismatch_energy_by_normal"])

json.dump(res, open(OUT + "/b22_bpoint_gauge.json", "w"), indent=1)
log("\nDONE t=%.1fs" % (time.time()-t0))
