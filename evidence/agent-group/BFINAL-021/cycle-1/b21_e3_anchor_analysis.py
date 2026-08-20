#!/usr/bin/env python3
"""BFINAL-021 E3: face-level flux-response anchor + B-basis identity.

Inputs: the b21_anchor re-run (b13_probe3 config + the new gated exports
phi/phiB/nutFrozen/nuEffFrozen/alpha per probe state and at the stageB2
baseline).

Tests:
  0. bit-reproducibility of the FD scan vs the archived b13_probe3 run.
  1. E1(b) direct: |div(phi_pm)| of the exported final fluxes (windowed
     convergence leaves phi div-free to machine level?).
  2. rebuild fidelity at B: phi*(U_B,p_B) rebuilt from the exported
     ingredients vs the exported phi_B (validates the offline flux-map
     representation = the solver's actual flux construction).
  3. E3 anchor: dphi_state-FD = (phi+ - phi-)/2h face field vs the
     operator reconstruction at B (Psi_U dU + Psi_p dp + PsiA da),
     and at A (same machinery at the old linearization point) --
     quantifies how much of the mismatch is the point/system error.
  4. P-row identity at B vs at A (exported-J A numbers from E2b).
  5. mobility staleness: mobF(A, molecular) vs mobF(B, turbulent) --
     the kf-basis error of the module round-2 adjoint.
"""
import time, json, os
import numpy as np
import scipy.sparse as sp

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
SA = "/home/ys/dsH/b16_states/1"          # state A (old linearization point)
SB = "/home/ys/dsH/b21_anchor/stageB2"    # state B exports (new run)
OLD = "/home/ys/dsH/b13_probe3/stageB2"   # archived states (bit-repro check)
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-021/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N; P0 = NV
ALPHA_REL = 0.4
NU = 5.19009e-05

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-021 E3: anchor analysis ===")

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

# ---------------- 0. bit-repro of the FD scan ----------------
res = {}
try:
    a = open(SB + "/stage_b2_fd_scan.tsv").read()
    b = open(OLD + "/stage_b2_fd_scan.tsv").read()
    res["fd_scan_bit_identical"] = bool(a == b)
    log("[0] FD scan bit-identical to archived b13_probe3: %s" % (a == b))
except Exception as e:
    res["fd_scan_bit_identical"] = "ERR %s" % e
    log("[0] FD scan compare failed: %s" % e)
# state exports identical to archived?
try:
    same = True
    for nm, tag, sgn, f in [("D1", "0.001", "p", "U"), ("D1", "0.001", "m", "p"),
                            ("D2", "0.0003", "p", "U"), ("D3", "0.001", "m", "p")]:
        x1 = np.loadtxt(SB + "/wstate_%s_h%s_%s_%s.mtx" % (nm, tag, sgn, f))
        x2 = np.loadtxt(OLD + "/wstate_%s_h%s_%s_%s.mtx" % (nm, tag, sgn, f))
        if x1.shape != x2.shape or not np.array_equal(x1, x2):
            same = False
            log("    state %s %s differs: maxabs=%.3e rel=%.3e"
                % (nm + f, tag, np.abs(x1 - x2).max() if x1.shape == x2.shape else -1,
                   np.linalg.norm(x1 - x2)/np.linalg.norm(x2)
                   if x1.shape == x2.shape else -1))
    res["wstate_bit_identical"] = same
    log("[0] wstate exports bit-identical: %s" % same)
except Exception as e:
    res["wstate_bit_identical"] = "ERR %s" % e

# ---------------- B-state fields ----------------
U_B = np.loadtxt(SB + "/wstate_baseline_U.mtx")
p_B = np.loadtxt(SB + "/wstate_baseline_p.mtx")
phi_B = np.loadtxt(SB + "/wstate_baseline_phi.mtx")
phiB_B = np.loadtxt(SB + "/wstate_baseline_phiB.mtx")
nut_B = np.loadtxt(SB + "/wstate_baseline_nutFrozen.mtx")
alpha_B = np.loadtxt(SB + "/wstate_baseline_alpha.mtx")
# A-state (old point)
def read_of_boundary_vals(path):
    import re
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
U_A = read_of_internal(SA + "/U").reshape(N, 3)
p_A = read_of_internal(SA + "/p")
phi_A = read_of_internal(SA + "/phi")
alpha_A = read_of_internal(SA + "/alpha")
phi_bv = read_of_boundary_vals(SA + "/phi")
phiB_A = np.zeros(nBnd); o = 0
for pn, nf in [("inlet", 84), ("outlet", 84), ("hotInlet", 196), ("hotOutlet", 196),
               ("solidEndWalls", 280), ("bottomWall", 1120), ("topWall", 1120),
               ("sideWalls", 4800)]:
    v = phi_bv[pn]
    phiB_A[o:o+nf] = (v[0] if v.size == 1 else v)
    o += nf
log("A-B recheck: U relL2=%.4e p relL2=%.4e phi relL2=%.4e"
    % (np.linalg.norm(U_A - U_B)/np.linalg.norm(U_B),
       np.linalg.norm(p_A - p_B)/np.linalg.norm(p_B),
       np.linalg.norm(phi_A - phi_B)/np.linalg.norm(phi_B)))

# ---------------- E1(b) direct: div(phi_pm) ----------------
log("\n[1] E1(b): continuity of exported probe fluxes")
divB = np.zeros(N)
np.add.at(divB, owner, phi_B); np.add.at(divB, nei, -phi_B)
np.add.at(divB, bc_cell, phiB_B)
log("    baseline B: |div(phi_B)|inf=%.3e L2=%.3e"
    % (np.abs(divB).max(), np.linalg.norm(divB)))

# ---------------- operator machinery at an arbitrary point ---------------
def build_basis(U_, p_, phi_, phiB_, alpha_, nuEff_):
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
    # boundary ic: fixedValue U -> -nuEff_c*magSf*deltaCoeffs (laplacian),
    # zeroGradient U -> phi_b (convective vic=1)
    ab_ = nuEff_[bc_cell]*bmag*bdel
    ic_ = np.where(Ufixed_b, -ab_, phiB_)
    sumOff_ = np.zeros(N)
    np.add.at(sumOff_, owner, np.abs(upper_))
    np.add.at(sumOff_, nei, np.abs(lower_))
    Db_ = diag_.copy()
    np.add.at(Db_, bc_cell, np.abs(ic_))
    D_rel_ = np.maximum(np.abs(Db_), sumOff_)/ALPHA_REL
    rAU_rel_ = 1.0/(D_rel_/V)
    rAU_u_ = rAU_rel_/ALPHA_REL
    mobF_ = w*rAU_rel_[owner] + (1.0-w)*rAU_rel_[nei]
    kf_ = mobF_*deltaCoeffs*magSf_all[:nIF]
    return dict(upw=upw_, qp=qp_, qn=qn_, lower=lower_, upper=upper_,
                diag=diag_, D_rel=D_rel_, rAU_rel=rAU_rel_, rAU_u=rAU_u_,
                mobF=mobF_, kf=kf_, divphi=divphi_)

def HbyA_field(U_, bs):
    src = (bs["D_rel"] - bs["diag"])[:, None]*U_
    # own row gets upper*U_nei ; nei row gets lower*U_own
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, bs["upper"][:, None]*U_[nei])
    np.add.at(offU, nei, bs["lower"][:, None]*U_[owner])
    H = (src - offU)/V[:, None]
    return bs["rAU_rel"][:, None]*H

def phi_rebuild(U_, p_, bs, phiB_fixed_):
    """absolute flux map: phiHbyA - pEqnFlux at (U_,p_) given coefficient
    basis bs (assembled with some phi).  Returns internal+boundary fluxes."""
    HbyA = HbyA_field(U_, bs)
    phi_int = np.einsum("fi,fi->f", sf, w[:, None]*HbyA[owner]
                        + (1-w)[:, None]*HbyA[nei])
    phi_int -= bs["kf"]*(p_[nei] - p_[owner])
    phi_bnd = np.zeros(nBnd)
    for bf in bidx:
        c = bc_cell[bf]
        phi_bnd[bf] = bSf[bf] @ HbyA[c] + bs["rAU_rel"][c]*bdel[bf]*bmag[bf]*p_[c]
    return np.concatenate([phi_int, phi_bnd])

# ---------------- 2. rebuild fidelity ----------------
log("\n[2] flux-map rebuild fidelity (internal faces + outlet boundary only;")
log("    non-assignable boundary fluxes are not rebuilt by the machinery):")
def fidelity(U_, p_, phi_, phiB_, alpha_, nuEff_):
    bs = build_basis(U_, p_, phi_, phiB_, alpha_, nuEff_)
    re = phi_rebuild(U_, p_, bs, phiB_)
    errI = re[:nIF] - phi_
    errB = re[nIF:][bidx] - phiB_[bidx]
    return (float(np.linalg.norm(errI)/np.linalg.norm(phi_)),
            float(np.linalg.norm(errB)/max(np.linalg.norm(phiB_[bidx]), 1e-300)))
nuEff_B = NU + nut_B
bsB = build_basis(U_B, p_B, phi_B, phiB_B, alpha_B, nuEff_B)
fiB, foB = fidelity(U_B, p_B, phi_B, phiB_B, alpha_B, nuEff_B)
res["rebuild_fidelity_B"] = {"internal_relL2": fiB, "outletBnd_relL2": foB}
log("    B(turbulent): internal relL2=%.4e outlet-boundary relL2=%.4e" % (fiB, foB))

# A-point rebuild fidelity (calibrates the offline machinery error itself)
nuEff_A_full = np.full(N, NU)
bsA_pre = build_basis(U_A, p_A, phi_A, phiB_A, alpha_A, nuEff_A_full)
fiA, foA = fidelity(U_A, p_A, phi_A, phiB_A, alpha_A, nuEff_A_full)
res["rebuild_fidelity_A"] = {"internal_relL2": fiA, "outletBnd_relL2": foA}
log("    A(molecular): internal relL2=%.4e outlet-boundary relL2=%.4e"
    " (instrument floor)" % (fiA, foA))

# ---------------- probe states: face-level anchor ----------------
dA = read_of_internal(SA + "/dAlphaDxh")
# design-chain tangent z = projection(filter(d)) (state-independent, valid at
# B too): rxc = R_alpha*(dAlphaDxh*z), NOT dAlphaDxh*rawD (B21 fix; b20's
# augmented check missed this chain).
zdirs = np.loadtxt("/home/ys/dsH/b15_export/stageB6_rxc_z_analytic.mtx").reshape(3, N)
dirs = zdirs

# A-basis (old point, molecular nuEff)
nuEff_A = np.full(N, NU)
bsA = build_basis(U_A, p_A, phi_A, phiB_A, alpha_A, nuEff_A)

def psi_channels(bs, U_, p_, dU, dp, da, dAlphaDxh_dirs=None):
    """dphi = Psi_U dU + Psi_p dp + Psi_alpha da at basis bs."""
    # Psi_p: -kf (dp_n - dp_o) ; boundary +kf_b dp_c
    dphi_p = -bs["kf"]*(dp[nei] - dp[owner])
    dphi_p_bnd = np.zeros(nBnd)
    for bf in bidx:
        c = bc_cell[bf]
        dphi_p_bnd[bf] = bs["rAU_rel"][c]*bdel[bf]*bmag[bf]*dp[c]
    # dH from offdiag + boundary diag (approx: same structure as b20)
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
    # alpha channel (design)
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
log("\n[3] face-level flux-response anchor (dphi_state-FD vs operator):")
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
        # E1(b): continuity of probe fluxes
        dp_ = np.zeros(N)
        np.add.at(dp_, owner, phi_p - phi_m)
        np.add.at(dp_, nei, -(phi_p - phi_m))
        np.add.at(dp_, bc_cell, phiBp - phiBm)
        cont_res = np.linalg.norm(dp_)
        # FD response
        dphi_true = np.concatenate([(phi_p - phi_m)/(2*h), (phiBp - phiBm)/(2*h)])
        dU = (U_p - U_m)/(2*h); dp_ = (p_p - p_m)/(2*h)
        da = dA*dirs[di]
        # B-basis channels
        duB, dpB, daB = psi_channels(bsB, U_B, p_B, dU, dp_, da)
        totB = duB + dpB + daB
        # A-basis channels
        duA, dpA, daA = psi_channels(bsA, U_A, p_A, dU, dp_, da)
        totA = duA + dpA + daA
        def mets(a, b):
            return {"cos": float(a @ b/(np.linalg.norm(a)*np.linalg.norm(b))),
                    "relL2": float(np.linalg.norm(a - b)/np.linalg.norm(b))}
        row = {
            "probe_cont_res_L2": float(cont_res),
            "B": {"total": mets(totB, dphi_true),
                  "UV": mets(duB + dpB, dphi_true),
                  "alpha_vs_resid": mets(daB, dphi_true - duB - dpB)},
            "A": {"total": mets(totA, dphi_true)},
            "identB": {"|rP|/|rxdP_re|":
                       float(np.linalg.norm(divF(duB + dpB) + divF(daB))
                             /np.linalg.norm(divF(daB))),
                       "cos": float(divF(duB + dpB) @ (-divF(daB))
                           /(np.linalg.norm(divF(duB + dpB))
                             *np.linalg.norm(divF(daB))))},
            "identA_re": {"|rP|/|rxdP_re|":
                       float(np.linalg.norm(divF(duA + dpA) + divF(daA))
                             /np.linalg.norm(divF(daA))),
                       "cos": float(divF(duA + dpA) @ (-divF(daA))
                           /(np.linalg.norm(divF(duA + dpA))
                             *np.linalg.norm(divF(daA))))},
            "alpha_share_B": float(np.linalg.norm(daB)
                                   /np.linalg.norm(totB)),
        }
        res["anchor"][nm][tag] = row
        log("  %s h=%s: cont=%.1e | B-tot cos=%.4f rel=%.4f | A-tot cos=%.4f"
            " rel=%.4f | identB |rP|=%.3f cos=%.4f | identA_re |rP|=%.3f cos=%.4f"
            " | alpha_share=%.3f"
            % (nm, tag, cont_res, row["B"]["total"]["cos"],
               row["B"]["total"]["relL2"], row["A"]["total"]["cos"],
               row["A"]["total"]["relL2"],
               row["identB"]["|rP|/|rxdP_re|"], row["identB"]["cos"],
               row["identA_re"]["|rP|/|rxdP_re|"], row["identA_re"]["cos"],
               row["alpha_share_B"]))

# ---------------- 5. mobility staleness ----------------
mobA = bsA["mobF"]; mobB = bsB["mobF"]
ratio = np.where(mobA > 0, mobB/mobA, 1.0)
res["mobility_staleness"] = {
    "median_ratio": float(np.median(ratio)),
    "p10": float(np.percentile(ratio, 10)),
    "p90": float(np.percentile(ratio, 90)),
    "min": float(ratio.min()), "max": float(ratio.max()),
    "frac_out_0.9_1.1": float(np.mean((ratio < 0.9) | (ratio > 1.1)))}
log("\n[5] mobility staleness mobB/mobA: median=%.4f [p10=%.3f p90=%.3f]"
    " min=%.3f max=%.3f frac_out_1pm10%%=%.3f"
    % (res["mobility_staleness"]["median_ratio"],
       res["mobility_staleness"]["p10"], res["mobility_staleness"]["p90"],
       res["mobility_staleness"]["min"], res["mobility_staleness"]["max"],
       res["mobility_staleness"]["frac_out_0.9_1.1"]))

json.dump(res, open(OUT + "/b21_e3_anchor_analysis.json", "w"), indent=1)
log("\nDONE t=%.1fs" % (time.time()-t0))
