#!/usr/bin/env python3
"""BFINAL-024 cycle-1: corrected-instrument re-attribution (H7 live).

Layers (all on the SAME state B, b24_gauge exports, FD side bit-identical
to B22 - verified by cmp):
  L1 face-level production P-row tangent vs true wstate FD, per direction:
     (a) pre-H7  instrument (du + kf)             [B23 corrected basis]
     (b) post-H7 instrument (du + kf + H7)        [the implemented channel]
     (c) + da (true alpha-FD channel, reference only - production has no
         design-row H7/da slot in J; quantifies the T9 A/B gap)
  L2 row-level amplification: |div(tangent) - div(FD)| / |div(FD)| per
     direction (the B21 mechanism: face-level mismatch amplified through
     div/pressure inversion).
  L3 system-level: FD-gate ADJ columns B22 (pre-H7) vs B24 (post-H7),
     gDP factor and J sign table (from the solver logs; recorded here for
     the attribution table).
No fitting factors; corrected-sign basis only (b23 lineage).
"""
import time, json
import numpy as np

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
SB = "/home/ys/dsH/b24_gauge/stageB2"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-024/cycle-1"

N = 33600
ALPHA_REL = 0.4
NU = 5.19009e-05

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-024: corrected-instrument re-attribution (H7 live) ===")

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
nB = len(bc); bidx = np.where(~Ufix)[0]     # assignable-U (outlet) boundary
import re
bfile = open(POLY + "/boundary").read()
patch_names = re.findall(r"^\s{4}(\w+)\s*$", bfile, re.M)
patch_sizes = [int(m) for m in re.findall(r"nFaces\s+(\d+);", bfile)]
patch_of_bface = np.concatenate(
    [np.full(s, i, dtype=np.int64) for i, s in enumerate(patch_sizes)])
OUTLET_IDX = patch_names.index("outlet")
sfI = Sf[:nIF]

# ---------------- corrected basis (b23 lineage, production sign) ---------
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
    return dict(lower=lower, upper=upper, diag=diag, D_rel=Drel, Db=Db,
                mob=V/Drel, kfV=kfv)

def grad_p(p_):
    pf = w*p_[owner] + (1 - w)*p_[nei]
    g = np.zeros((N, 3))
    np.add.at(g, owner, sfI*pf[:, None])
    np.add.at(g, nei, -(sfI*pf[:, None]))
    pb = p_[bc].copy()
    pb[patch_of_bface == OUTLET_IDX] = 0.0     # fixedValue p (outlet) -> 0
    np.add.at(g, bc, bSf*pb[:, None])
    return g/V[:, None]

def interp_vec(F):
    return w[:, None]*F[owner] + (1-w)[:, None]*F[nei]

U_B = np.loadtxt(SB + "/wstate_baseline_U.mtx")
p_B = np.loadtxt(SB + "/wstate_baseline_p.mtx")
phi_B = np.loadtxt(SB + "/wstate_baseline_phi.mtx")
phiB_B = np.loadtxt(SB + "/wstate_baseline_phiB.mtx")
nuE = NU + np.loadtxt(SB + "/wstate_baseline_nutFrozen.mtx")
al = np.loadtxt(SB + "/wstate_baseline_alpha.mtx")
mob_ex = np.loadtxt(SB + "/wstate_baseline_primalPressureMobility.mtx")
bs = build_basis(phi_B, phiB_B, al, nuE)
mob_err = float(np.linalg.norm(bs["mob"] - mob_ex)/np.linalg.norm(mob_ex))
log("[0] corrected-basis mobility reproduction vs exact export: %.3e "
    "(b23: 2.19e-15)" % mob_err)
assert mob_err < 1e-10
mob = bs["mob"]
mobF = w*mob[owner] + (1-w)*mob[nei]
kf = mobF*dc*magSfI

res = {"mobility_reproduction": mob_err, "anchor": {}, "rowamp": {}}

def divF(x):
    d = np.zeros(N)
    np.add.at(d, owner, x[:nIF]); np.add.at(d, nei, -x[:nIF])
    np.add.at(d, bc, x[nIF:])
    return d

ladder = {"D1": ["0.0003", "0.001"], "D2": ["0.0003", "0.001"],
          "D3": ["0.0003", "0.001"]}
for nm in ("D1", "D2", "D3"):
    res["anchor"][nm] = {}
    for tag in ladder[nm]:
        h = float(tag)
        pre = SB + "/wstate_%s_h%s_" % (nm, tag)
        U_p = np.loadtxt(pre + "p_U.mtx"); p_p = np.loadtxt(pre + "p_p.mtx")
        U_m = np.loadtxt(pre + "m_U.mtx"); p_m = np.loadtxt(pre + "m_p.mtx")
        phi_p = np.loadtxt(pre + "p_phi.mtx"); phi_m = np.loadtxt(pre + "m_phi.mtx")
        phiBp = np.loadtxt(pre + "p_phiB.mtx"); phiBm = np.loadtxt(pre + "m_phiB.mtx")
        ap = np.loadtxt(pre + "p_alpha.mtx"); am = np.loadtxt(pre + "m_alpha.mtx")
        dphi_true = np.concatenate([(phi_p - phi_m)/(2*h),
                                    (phiBp - phiBm)/(2*h)])
        dU = (U_p - U_m)/(2*h); dp_ = (p_p - p_m)/(2*h)
        da = (ap - am)/(2*h)

        # du channel: production relaxed mapping on the corrected basis
        offdU = np.zeros((N, 3))
        np.add.at(offdU, owner, bs["upper"][:, None]*dU[nei])
        np.add.at(offdU, nei, bs["lower"][:, None]*dU[owner])
        dH = -offdU/V[:, None]
        A_u = np.abs(bs["Db"])/V
        dHbyA = (ALPHA_REL/A_u)[:, None]*dH + (1.0 - ALPHA_REL)*dU
        dphi_u = np.einsum("fi,fi->f", sfI, interp_vec(dHbyA))
        dphi_u_b = np.zeros(nB)
        for bf in bidx:
            dphi_u_b[bf] = bSf[bf] @ dHbyA[bc[bf]]
        dphi_u_all = np.concatenate([dphi_u, dphi_u_b])

        # dp channel: kf (pEqn.flux) + H7 (phiHbyA pressure channel)
        dphi_kf = -kf*(dp_[nei] - dp_[owner])
        dphi_kf_b = np.zeros(nB)
        for bf in bidx:
            c = bc[bf]
            dphi_kf_b[bf] = mob[c]*bdel[bf]*bmag[bf]*dp_[c]
        gdp = grad_p(dp_)
        h7 = -np.einsum("fi,fi->f", sfI, interp_vec(mob[:, None]*gdp))
        h7_b = np.zeros(nB)
        for bf in bidx:
            c = bc[bf]
            h7_b[bf] = -(bSf[bf] @ (mob[c]*gdp[c]))
        dphi_kf_all = np.concatenate([dphi_kf, dphi_kf_b])
        h7_all = np.concatenate([h7, h7_b])

        # da channel (true alpha-FD, reference only)
        src = (bs["D_rel"] - bs["diag"])[:, None]*U_B
        offU = np.zeros((N, 3))
        np.add.at(offU, owner, bs["upper"][:, None]*U_B[nei])
        np.add.at(offU, nei, bs["lower"][:, None]*U_B[owner])
        Hcur = (src - offU)/V[:, None] - grad_p(p_B)
        drAU = -(mob**2/ALPHA_REL)*da
        dHsrc = ((1.0/ALPHA_REL - 1.0)*da)[:, None]*U_B
        dHbyA_a = drAU[:, None]*Hcur + mob[:, None]*dHsrc
        dphi_a = np.einsum("fi,fi->f", sfI, interp_vec(dHbyA_a))
        dphi_a += (drAU[owner]*w + drAU[nei]*(1-w))*dc*magSfI*(p_B[nei] - p_B[owner])
        dphi_a -= np.einsum(
            "fi,fi->f", sfI,
            (drAU[owner]*w + drAU[nei]*(1-w))[:, None]*interp_vec(grad_p(p_B)))
        dphi_a_b = np.zeros(nB)
        for bf in bidx:
            c = bc[bf]
            dphi_a_b[bf] = (bSf[bf] @ dHbyA_a[c]) \
                + drAU[c]*bdel[bf]*bmag[bf]*p_B[c] \
                - (bSf[bf] @ (drAU[c]*grad_p(p_B)[c]))
        dphi_a_all = np.concatenate([dphi_a, dphi_a_b])

        def mets(a, b):
            return {"cos": float(a @ b/(np.linalg.norm(a)*np.linalg.norm(b))),
                    "relL2": float(np.linalg.norm(a - b)/np.linalg.norm(b))}
        pre_h7 = dphi_u_all + dphi_kf_all
        post_h7 = pre_h7 + h7_all
        post_da = post_h7 + dphi_a_all
        row = {"pre_h7": mets(pre_h7, dphi_true),
               "post_h7": mets(post_h7, dphi_true),
               "post_h7_da": mets(post_da, dphi_true),
               "h7_share": float(np.linalg.norm(h7_all)
                                 /np.linalg.norm(dphi_true))}
        res["anchor"][nm][tag] = row
        # row-level: the P-row continuity residual the pressure equation
        # must absorb (the truth dphi is divergence-free to ~1e-13, so the
        # meaningful metric is |div(tangent)| itself, pre vs post H7)
        rP_pre = divF(pre_h7)
        rP_post = divF(post_h7)
        rP_truth = divF(dphi_true)
        res["rowamp"][nm] = {
            "face_pre": row["pre_h7"]["relL2"], "face_post": row["post_h7"]["relL2"],
            "row_pre": float(np.linalg.norm(rP_pre)),
            "row_post": float(np.linalg.norm(rP_post)),
            "row_truth": float(np.linalg.norm(rP_truth))}
        log("  %s h=%s: pre-H7 relL2=%.4f -> post-H7 %.4f (+da ref %.4f) "
            "|H7|/|FD|=%.3f ; |div(tangent)| %.3e -> %.3e (truth %.3e)"
            % (nm, tag, row["pre_h7"]["relL2"], row["post_h7"]["relL2"],
               row["post_h7_da"]["relL2"], row["h7_share"],
               res["rowamp"][nm]["row_pre"], res["rowamp"][nm]["row_post"],
               res["rowamp"][nm]["row_truth"]))

# L3: system-level deltas recorded from the solver logs (no recomputation)
res["fdgate"] = {
    "ADJ_gDP": {"D1": {"b22": -5.44409679856, "b24": -5.23029384627,
                        "FD_h1e-3": -2.53745977185},
                 "D2": {"b22": 1.13007807765, "b24": 1.0505249811,
                        "FD_h1e-3": 0.523409452158},
                 "D3": {"b22": 14.1276593948, "b24": 13.6034225815,
                        "FD_h1e-3": 6.47721554707}},
    "ADJ_J": {"D1": {"b22": -0.0349244139673, "b24": -0.0358350400528,
                      "FD_h1e-3": 0.0428559144263},
               "D2": {"b22": -0.0233184637414, "b24": -0.0249259903913,
                      "FD_h1e-3": -0.00346892642156},
               "D3": {"b22": -0.0139740506764, "b24": -0.0120030122811,
                      "FD_h1e-3": -0.000282658160322}},
}
f22 = {}; f24 = {}
for d in ("D1", "D2", "D3"):
    f22[d] = res["fdgate"]["ADJ_gDP"][d]["b22"]/res["fdgate"]["ADJ_gDP"][d]["FD_h1e-3"]
    f24[d] = res["fdgate"]["ADJ_gDP"][d]["b24"]/res["fdgate"]["ADJ_gDP"][d]["FD_h1e-3"]
res["fdgate"]["gDP_factor_ADJoverFD"] = {"b22": f22, "b24": f24}
log("[L3] gDP factor ADJ/FD: B22 %.3f/%.3f/%.3f -> B24 %.3f/%.3f/%.3f"
    % (f22["D1"], f22["D2"], f22["D3"], f24["D1"], f24["D2"], f24["D3"]))

json.dump(res, open(OUT + "/b24_reattribution.json", "w"), indent=1)
log("\nDONE t=%.1fs" % (time.time() - t0))
