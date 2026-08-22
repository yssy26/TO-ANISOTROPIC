#!/usr/bin/env python3
"""BFINAL-026 attempt-2 M1: source-layer identity.

    lhs = b_TC^T w_true(H7)        [b_TC: PHYSICAL/unscaled Q-type B-state export]
    rhs = FD_J(Q) - C - Gx         [C/Gx per PREREGISTRATION-pinned semantics]
    check lhs/rhs, sign, rel-dev per D1/D2/D3 + pre-H7 control (H7 sensitivity).

b18_folding_lab.py is the template; DEFECTS corrected per b26a_NOTEBOOK §4:
  D1  exact pyramid centroids (V-gate vs exported V) instead of vertex-mean
  D2  FD pins from B25 stage_b2_fd_scan.tsv (h=1e-3 criterion rows)
  D3  z from b8_verify_diag/stageB6_rxc_z_analytic.mtx (chain-rule Dxh dirs)
  D4  SB repointed to /home/ys/dsH/b25_qgate/stageB2
  D5  patch sizes read from polyMesh/boundary (regex)
  D6  b_TC read np.loadtxt + route_V0 layout
  +NEW mixed-BC snGrad correction for C_field boundary (bottomWall mixed
   vf=0.00175607, refGradient=0 => snGrad = vf*(refValue-T_cell)*bdel, not
   (T_bvals-T_cell)*bdel which overcounts ~500x).

Gates: mobility gate <1e-12 (b23 recipe, SB=b25 stageB2); V-gate ~1e-12;
route_V0 self-gate vs exported b_TC (relL2 reported, no hard threshold).
"""
import time, json, re
import numpy as np

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
SB = "/home/ys/dsH/b25_qgate/stageB2"
Q1 = "/home/ys/dsH/b25_qgate/1"
E26 = "/home/ys/dsH/b8_verify_diag"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-026/cycle-1"

N = 33600; NV = 3 * N; NUNK = 4 * N
ALPHA_REL = 0.4
NU = 5.19009e-05
SMALL = 1e-15

# FD_J(Q) criterion rows (h=1e-3) from stage_b2_fd_scan.tsv (platform stability
# reported separately): D1 +0.04252469623637622, D2 -0.003528972517971574,
# D3 -0.0004449511895182612.  Main refs: D1(h=1e-5) +0.04244654087520727,
# D2(h=1e-4) -0.003699503574594587, D3(h=1e-4) -0.0004211702528400529.
FD_J = {"D1": 0.04252469623637622,
        "D2": -0.003528972517971574,
        "D3": -0.0004449511895182612}
FD_MAINREF = {"D1": 0.04244654087520727,
              "D2": -0.003699503574594587,
              "D3": -0.0004211702528400529}
ADJ_J = {"D1": -0.0358350400485765,
         "D2": -0.02492599038981423,
         "D3": -0.01200301226991504}
ADJ_GDP = {"D1": -5.23029384626718,
           "D2": 1.0505249811003,
           "D3": 13.60342258150915}
PRODH7SHARE = 0.701947825646

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-026 attempt-2 M1 instrument ===")

# ---------------- field parsing (b18 read_field, validated) ----------------
def read_field(path):
    txt = open(path).read()
    i = txt.index("internalField")
    j = txt.index("\n(", i)
    k = txt.index("\n)", j)
    if txt[i:j].split()[1] == "uniform":
        seg = txt[i:j]
        u = seg[seg.index("uniform") + len("uniform"):].strip()
        toks = []
        for tok in u.replace(";", " ").replace("(", " ").replace(")", " ").split():
            try: toks.append(float(tok))
            except ValueError: break
        val = np.array(toks) if len(toks) > 1 else float(toks[0])
        internal = None
    else:
        vals = []
        for line in txt[j + 1:k].splitlines():
            for tok in line.replace(";", " ").split():
                for tt in tok.replace("(", " ").replace(")", " ").split():
                    vals.append(float(tt))
        internal = np.array(vals)
        val = None
    bnd = {}
    b = txt.index("boundaryField")
    lines = txt[b:].splitlines()
    p = 0
    while p < len(lines):
        s = lines[p].strip()
        if s in ("{", "}", ""):
            p += 1; continue
        if s.endswith("{"):
            name = s[:-1].strip(); hdr = p
        elif (p + 1 < len(lines) and lines[p + 1].strip() == "{"
              and " " not in s and not s.startswith(("type", "value", "internalField"))):
            name = s; hdr = p + 1
        else:
            p += 1; continue
        if name == "boundaryField":
            p += 1; continue
        depth = 1; q = hdr + 1; bval = None; btype = None
        while q < len(lines) and depth > 0:
            t = lines[q].strip()
            if t.endswith("{") and t != "{": depth += 1
            if t == "{": depth += 1
            if t == "}": depth -= 1
            if t.startswith("type"): btype = t.split()[1].rstrip(";")
            if t.startswith("uniform") and depth == 1:
                u = t[len("uniform"):].rstrip(";").strip()
                if u.startswith("("):
                    bval = ("u", np.array([float(x) for x in u.strip("()").split()]))
                else:
                    bval = ("u", float(u))
            if t.startswith("value") and depth == 1:
                if " nonuniform" in t or t.rstrip(";").endswith("nonuniform"):
                    bval = ("list", q)
                elif "uniform" in t:
                    u = t[t.index("uniform") + len("uniform"):].rstrip(";").strip()
                    if u.startswith("("):
                        bval = ("u", np.array([float(x) for x in u.strip("()").split()]))
                    else:
                        bval = ("u", float(u))
                else:
                    bval = ("list", q)
            q += 1
        if isinstance(bval, tuple) and bval[0] == "list":
            q0 = bval[1]; qq = q0
            while "(" not in lines[qq]: qq += 1
            kk = qq
            while ")" not in lines[kk]: kk += 1
            vals = []
            for line in lines[qq:kk + 1]:
                for tok in line.strip("()").replace("(", " ").replace(")", " ").split():
                    try: vals.append(float(tok))
                    except ValueError: pass
            bnd[name] = (btype, np.array(vals))
        elif bval is not None:
            bnd[name] = (btype, bval[1])
        else:
            bnd[name] = (btype, None)
        p = q
    return internal, val, bnd

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
        faces.append([int(v) for v in s[s.index("(") + 1:-1].split()])
def rl(p):
    t = open(p).read(); i = t.index("\n("); j = t.index("\n)", i)
    return np.array([int(x) for x in t[i + 2:j].split()], dtype=np.int64)
oa = rl(POLY + "/owner"); na = rl(POLY + "/neighbour")
nIF = len(na); nF = len(faces)
Sf = np.zeros((nF, 3))
for fi in range(nF):
    fv = pts[faces[fi]]
    Sf[fi] = 0.5 * np.cross(fv, np.roll(fv, -1, axis=0)).sum(axis=0)
magSf = np.linalg.norm(Sf, axis=1)
fCtrs = np.array([pts[fv].mean(axis=0) for fv in faces])
owner = oa[:nIF]; nei = na[:nIF]
bc = oa[nIF:]                       # boundary face -> cell
nB = nF - nIF
bfile = open(POLY + "/boundary").read()
patch_names = re.findall(r"^\s{4}(\w+)\s*$", bfile, re.M)
patch_sizes = [int(m) for m in re.findall(r"nFaces\s+(\d+);", bfile)]
patch_of_bface = np.concatenate(
    [np.full(s, i, dtype=np.int64) for i, s in enumerate(patch_sizes)])
bface_start = {}
o = nIF
for pn, sz in zip(patch_names, patch_sizes):
    bface_start[pn] = o; o += sz

# --- exact pyramid centroids (primitiveMeshCellCentresAndVols.C) ---
cEst = np.zeros((N, 3)); nFc = np.zeros(N, dtype=np.int64)
np.add.at(cEst, oa, fCtrs); np.add.at(nFc, oa, 1)
np.add.at(cEst, na, fCtrs[:nIF]); np.add.at(nFc, na, 1)
cEst /= nFc[:, None]
num = np.zeros((N, 3)); vol3 = np.zeros(N)
SfI = Sf[:nIF]; fCtrsI = fCtrs[:nIF]
p3o = np.einsum("ij,ij->i", SfI, fCtrsI - cEst[owner])
pc_o = 0.75 * fCtrsI + 0.25 * cEst[owner]
np.add.at(num, owner, pc_o * p3o[:, None]); np.add.at(vol3, owner, p3o)
p3n = np.einsum("ij,ij->i", -SfI, fCtrsI - cEst[nei])
pc_n = 0.75 * fCtrsI + 0.25 * cEst[nei]
np.add.at(num, nei, pc_n * p3n[:, None]); np.add.at(vol3, nei, p3n)
SfB = Sf[nIF:]; fCtrsB = fCtrs[nIF:]
p3b = np.einsum("ij,ij->i", SfB, fCtrsB - cEst[bc])
pc_b = 0.75 * fCtrsB + 0.25 * cEst[bc]
np.add.at(num, bc, pc_b * p3b[:, None]); np.add.at(vol3, bc, p3b)
Cp = num / vol3[:, None]
V_pyr = vol3 / 3.0
V = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")[:, 0]
V_gate = float(np.linalg.norm(V_pyr - V) / np.linalg.norm(V))
log("V gate (pyramid) relL2=%.3e" % V_gate)
# vertex-mean centroid (b23 recipe, for comparison)
cv = {}
for fi in range(nF):
    cs = (oa[fi], na[fi]) if fi < nIF else (oa[fi],)
    for c in cs: cv.setdefault(c, set()).update(faces[fi])
Cv = np.zeros((N, 3))
for c, vs in cv.items(): Cv[c] = pts[sorted(vs)].mean(axis=0)
C = Cp  # exact pyramid (corrected per NOTEBOOK DEFECT 1)
dvec = C[nei] - C[owner]
magSfI = magSf[:nIF]
dc = magSfI / np.einsum("ij,ij->i", SfI, dvec)
# vertex-mean dc for comparison
dv2 = Cv[nei] - Cv[owner]
dc2 = magSfI / np.einsum("ij,ij->i", SfI, dv2)
log("dc pyramid vs vertex-mean: relL2=%.3e  (deltaCoeffs centroid sensitivity)"
    % (np.linalg.norm(dc - dc2) / np.linalg.norm(dc)))
w = np.loadtxt(MESH + "/b16mesh_weights.mtx")
bnd = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
bc_cell = bnd[:, 0].astype(np.int64); bSf = bnd[:, 1:4]
bmag = bnd[:, 4]; bdel = bnd[:, 5]
Ufix = bnd[:, 7] > 0.5
bidx = np.where(~Ufix)[0]

# ---------------- fields (b25_qgate/1 exclusively) ----------------
U_i, _, U_b = read_field(Q1 + "/U"); U = U_i.reshape(N, 3)
p_i, _, p_b = read_field(Q1 + "/p"); p = p_i
phi_i, _, phi_b = read_field(Q1 + "/phi"); phi_int = phi_i
phTh_i, _, phTh_b = read_field(Q1 + "/phiThermal"); phiTh_int = phTh_i
cfm_i, _, _ = read_field(Q1 + "/coldFaceMask"); coldmask_int = cfm_i
Tb_i, _, Tb_b = read_field(Q1 + "/Tb"); Tb = Tb_i
T_i, _, T_b = read_field(Q1 + "/T"); T = T_i
alpha_i, _, _ = read_field(Q1 + "/alpha"); alpha = alpha_i
dAlphaDxh = read_field(Q1 + "/dAlphaDxh")[0]
dtdxh_i, _, _ = read_field(Q1 + "/dDTDxh")
dDTDxh = dtdxh_i.reshape(N, 6)
assert len(phi_int) == nIF, "phi internal len %d != nIF %d" % (len(phi_int), nIF)
assert len(phiTh_int) == nIF, "phiThermal internal len %d" % len(phiTh_int)
assert len(coldmask_int) == nIF, "coldFaceMask internal len %d" % len(coldmask_int)
p_btype = {pn: p_b[pn][0] for pn in patch_names}   # for H7 stage-2 fixedValue-p skip
log("fields loaded (t=%.1fs)" % (time.time() - t0))

# phi/p boundary values
phi_bvals = np.zeros(nB); p_bvals = np.zeros(nB)
for pn, sz in zip(patch_names, patch_sizes):
    o0 = bface_start[pn] - nIF
    btp, bvp = phi_b[pn]
    if isinstance(bvp, float): phi_bvals[o0:o0 + sz] = bvp
    elif bvp is not None: phi_bvals[o0:o0 + sz] = bvp
    btp2, bvp2 = p_b[pn]
    if isinstance(bvp2, float): p_bvals[o0:o0 + sz] = bvp2
    elif bvp2 is not None: p_bvals[o0:o0 + sz] = bvp2
    else: p_bvals[o0:o0 + sz] = p[bc_cell[o0:o0 + sz]]
# T patch types + values + mixed params
T_btype = {}; T_bval = {}
for pn in patch_names:
    bt, bv = T_b[pn]; T_btype[pn] = bt; T_bval[pn] = bv
Ttxt = open(Q1 + "/T").read()
mixed_params = {}
for pn in patch_names:
    m = re.search(r"^\s{4}%s\s*\{(.*?)\}" % pn,
                  Ttxt[Ttxt.index("boundaryField"):], re.S | re.M)
    seg = m.group(1) if m else ""
    vf = re.search(r"valueFraction\s+uniform\s+([\d.eE+-]+)", seg)
    rv = re.search(r"refValue\s+uniform\s+([\d.eE+-]+)", seg)
    rg = re.search(r"refGradient\s+uniform\s+([\d.eE+-]+)", seg)
    mixed_params[pn] = (float(vf.group(1)) if vf else None,
                        float(rv.group(1)) if rv else None,
                        float(rg.group(1)) if rg else None)
log("T patches: %s" % {pn: T_btype[pn] for pn in patch_names})
log("mixed params: %s" % mixed_params)

# ---------------- g face functional (B0.2-verified formulas) ----------------
g = np.zeros(nIF + nB)
jumpT = T[nei] - T[owner]
down_nei = phiTh_int >= 0.0
adjDown = np.where(down_nei, Tb[nei], Tb[owner])
g[:nIF] = -coldmask_int * adjDown * jumpT
Mflow, Tref = np.loadtxt(SB + "/wstate_outletMeta.mtx")
o0 = bface_start["outlet"] - nIF
ophi = phi_bvals[o0:o0 + patch_sizes[patch_names.index("outlet")]]
Toc = T[bc_cell[o0:o0 + patch_sizes[patch_names.index("outlet")]]]
Tmix = float(np.sum(ophi * Toc) / Mflow)
dJdphi_out = -(Toc - Tmix) / (Mflow * Tref)
g[nIF + o0:nIF + o0 + len(Toc)] += dJdphi_out
log("g: |internal|max=%.6e nz=%d ; |outlet dJ/dphi|max=%.6e Tmix=%.6f M=%.6e"
    % (np.abs(g[:nIF]).max(), int((g[:nIF] != 0).sum()),
       np.abs(dJdphi_out).max(), Tmix, Mflow))

# ---------------- mobility (b23 build_basis, SB=b25 stageB2) ----------------
def build_basis(phi_, phiB_, alpha_, nuEff_):
    qp = (phi_ >= 0).astype(float); qn = 1 - qp
    divp = np.zeros(N)
    np.add.at(divp, owner, phi_); np.add.at(divp, nei, -phi_)
    np.add.at(divp, bc_cell, phiB_)
    gam = w * nuEff_[owner] + (1 - w) * nuEff_[nei]
    kfv = gam * magSfI * dc
    lower = -qp * phi_ - kfv; upper = qn * phi_ - kfv
    diag = np.zeros(N)
    np.add.at(diag, owner, -lower); np.add.at(diag, nei, -upper)
    diag -= divp; diag += alpha_ * V
    sumOff = np.zeros(N)
    np.add.at(sumOff, owner, np.abs(upper)); np.add.at(sumOff, nei, np.abs(lower))
    ab = nuEff_[bc_cell] * bmag * bdel
    ic = np.where(Ufix, -ab, phiB_)
    Db = diag.copy(); np.add.at(Db, bc_cell, np.abs(ic))
    Drel = np.maximum(np.abs(Db), sumOff) / ALPHA_REL
    return dict(lower=lower, upper=upper, diag=diag, D_rel=Drel, Db=Db,
                mob=V / Drel)

U_B = np.loadtxt(SB + "/wstate_baseline_U.mtx")
p_B = np.loadtxt(SB + "/wstate_baseline_p.mtx")
phi_B = np.loadtxt(SB + "/wstate_baseline_phi.mtx")
phiB_B = np.loadtxt(SB + "/wstate_baseline_phiB.mtx")
nut = np.loadtxt(SB + "/wstate_baseline_nutFrozen.mtx")
al_B = np.loadtxt(SB + "/wstate_baseline_alpha.mtx")
nuE = NU + nut
mob_ex = np.loadtxt(SB + "/wstate_baseline_primalPressureMobility.mtx")
bs = build_basis(phi_B, phiB_B, al_B, nuE)
mob_gate = float(np.linalg.norm(bs["mob"] - mob_ex) / np.linalg.norm(mob_ex))
log("MOBILITY GATE (pyramid dc): relL2=%.3e  (need <1e-12)" % mob_gate)
if mob_gate >= 1e-12:
    bs2 = dict(bs)
    dc2v = magSfI / np.einsum("ij,ij->i", SfI, dv2)
    bs2["dc"] = dc2v
    # rebuild with vertex-mean dc to isolate
    qp = (phi_B >= 0).astype(float); qn = 1 - qp
    divp = np.zeros(N)
    np.add.at(divp, owner, phi_B); np.add.at(divp, nei, -phi_B)
    np.add.at(divp, bc_cell, phiB_B)
    gam = w * nuE[owner] + (1 - w) * nuE[nei]
    kfv = gam * magSfI * dc2v
    lower = -qp * phi_B - kfv; upper = qn * phi_B - kfv
    diag = np.zeros(N)
    np.add.at(diag, owner, -lower); np.add.at(diag, nei, -upper)
    diag -= divp; diag += al_B * V
    sumOff = np.zeros(N)
    np.add.at(sumOff, owner, np.abs(upper)); np.add.at(sumOff, nei, np.abs(lower))
    ab = nuE[bc_cell] * bmag * bdel
    ic = np.where(Ufix, -ab, phiB_B)
    Db = diag.copy(); np.add.at(Db, bc_cell, np.abs(ic))
    Drel = np.maximum(np.abs(Db), sumOff) / ALPHA_REL
    mob2 = V / Drel
    g2 = float(np.linalg.norm(mob2 - mob_ex) / np.linalg.norm(mob_ex))
    log("  vertex-mean dc gate: relL2=%.3e" % g2)
    # pick the better one for the run
    if g2 < mob_gate:
        C = Cv; dc = dc2v
        bs = build_basis(phi_B, phiB_B, al_B, nuE)
        mob_gate = g2
        log("  -> using vertex-mean centroids (gate improved)")
log("MOBILITY GATE final: relL2=%.3e" % mob_gate)
assert mob_gate < 1e-12, "BLOCKED: mobility gate failed (%.3e)" % mob_gate
mob = bs["mob"]                       # == primalPressureMobility (rAtU = rAU, SIMPLE no-consistent)

# ---------------- route_V1 self-gate vs exported b_TC ----------------
# Production V1 folding mirror (solveDiscreteFlowAdjointProduction.H L351-497,
# fix-3 route; exported b_TC is this production route, written 12:57 same state
# as the B25 baselines):
#   internal:  hA += alphaRel*prodRAU * w * (Sf*g)      [alphaRel*prodRAU==mob]
#              U(own/nei) += (1-alphaRel)*w/(1-w)*(Sf*g)          [direct]
#              P(own) += kf*g ; P(nei) -= kf*g          [operator J_PP sign]
#              h7G += mob*w/(1-w)*g*Sf                  [H7 stage-1]
#   boundary (assignable-U only = outlet): hA/direct/H7 stage-1 +
#              pressureFixed(outlet) -> P(c) += mob*dc_b*|Sf|_b*g
#   deltaH^T:  U(nei) -= upper/V_o*hA_o ; U(own) -= lower/V_n*hA_n
#   boundary H-diagonal: (-internalCoeffs+cav)*Vinv*hA(c) -- identically zero
#              here (outlet zeroGradient-U => internalCoeffs=0; fixedValue-U
#              patches have hA=0); kept as documented parity no-op
#   H7 stage-2: q = Sf&(h7G_o/V_o - h7G_n/V_n) ; P(own) -= w*q ; P(nei) -= (1-w)*q
#              + zeroGradient-p boundary self -(h7G_c&S_b)/V_c
# prodUpper/prodLower = unrelaxed momentum off-diagonals (bs upper/lower);
# prodRAU = mob/ALPHA_REL (alphaRel*prodRAU == rAtU == mob, machine-verified).
kf_int = (w * mob[owner] + (1 - w) * mob[nei]) * magSfI * dc
kf_b = mob[bc_cell] * bmag * bdel
mobc = mob[:, None]

def route_V0():
    """OLD defective route (pre fix-3): plain w*Sf*g interpolation on the U
    rows + kf with sign OPPOSITE to the operator's own J_PP block.  Kept only as
    a diagnostic proving the 0.32 relL2 mismatch is the old-vs-V1 routing
    difference, not a state mismatch (DERIVATION_FIX3 defect list)."""
    b = np.zeros(NUNK)
    gU = np.zeros((N, 3))
    ft = SfI * g[:nIF][:, None]
    np.add.at(gU, owner, w[:, None] * ft)
    np.add.at(gU, nei, (1.0 - w)[:, None] * ft)
    b[0:NV:3] = gU[:, 0]; b[1:NV:3] = gU[:, 1]; b[2:NV:3] = gU[:, 2]
    gP = np.zeros(N)
    np.add.at(gP, owner, -kf_int * g[:nIF])
    np.add.at(gP, nei, kf_int * g[:nIF])
    for f in range(nB):
        if Ufix[f]: continue
        c = bc_cell[f]
        gf = g[nIF + f]
        b[3 * c:3 * c + 3] += bSf[f] * gf
        gP[c] += kf_b[f] * gf
    b[NV:] = gP
    return b

def route_V1():
    b = np.zeros(NUNK)
    gf_i = g[:nIF]; gf_b = g[nIF:]
    hA = np.zeros((N, 3)); h7G = np.zeros((N, 3))
    gU = np.zeros((N, 3)); gP = np.zeros(N)
    ft = SfI * gf_i[:, None]
    # internal: hA (alphaRel*prodRAU == mob)
    np.add.at(hA, owner, mobc[owner] * w[:, None] * ft)
    np.add.at(hA, nei, mobc[nei] * (1.0 - w)[:, None] * ft)
    # internal: direct (1-alphaRel) U rows
    np.add.at(gU, owner, (1.0 - ALPHA_REL) * w[:, None] * ft)
    np.add.at(gU, nei, (1.0 - ALPHA_REL) * (1.0 - w)[:, None] * ft)
    # internal: kf P rows (+own/-nei, operator J_PP sign)
    np.add.at(gP, owner, kf_int * gf_i)
    np.add.at(gP, nei, -kf_int * gf_i)
    # internal: H7 stage-1 -> h7G
    np.add.at(h7G, owner, mobc[owner] * w[:, None] * ft)
    np.add.at(h7G, nei, mobc[nei] * (1.0 - w)[:, None] * ft)
    # boundary assignable-U only (outlet): hA/direct/H7 stage-1 + pressureFixed
    for f in bidx:
        c = bc_cell[f]
        ftb = bSf[f] * gf_b[f]
        hA[c] += mob[c] * ftb
        gU[c] += (1.0 - ALPHA_REL) * ftb
        h7G[c] += mob[c] * ftb
        gP[c] += mob[c] * bdel[f] * bmag[f] * gf_b[f]
    # deltaH^T
    Vinvo = 1.0 / V[owner]; Vinn = 1.0 / V[nei]
    gU[nei] -= (bs["upper"] * Vinvo)[:, None] * hA[owner]
    gU[owner] -= (bs["lower"] * Vinn)[:, None] * hA[nei]
    # boundary H-diagonal: identically zero here (parity note above)
    # H7 stage-2 internal
    q = np.einsum("fi,fi->f", SfI,
                  h7G[owner] * Vinvo[:, None] - h7G[nei] * Vinn[:, None])
    np.add.at(gP, owner, -w * q)
    np.add.at(gP, nei, -(1.0 - w) * q)
    # H7 stage-2 boundary: zeroGradient-p patches only (fixedValue-p dp_b = 0)
    for pn, sz in zip(patch_names, patch_sizes):
        if p_btype[pn] == "fixedValue": continue
        o0 = bface_start[pn] - nIF
        for f in range(o0, o0 + sz):
            c = bc_cell[f]
            gP[c] -= float(h7G[c] @ bSf[f]) / V[c]
    b[0:NV:3] = gU[:, 0]; b[1:NV:3] = gU[:, 1]; b[2:NV:3] = gU[:, 2]
    b[NV:] = gP
    return b

m_exp = np.loadtxt("/home/ys/dsH/b25_qgate/b18rhs_thermalCoupling.mtx")
bexp = np.zeros(NUNK)
bexp[0:NV:3] = m_exp[:, 0]; bexp[1:NV:3] = m_exp[:, 1]
bexp[2:NV:3] = m_exp[:, 2]; bexp[NV:] = m_exp[:, 3]
bV0 = route_V0()
bV1 = route_V1()
r0 = float(np.linalg.norm(bV0 - bexp) / np.linalg.norm(bexp))
c0 = float(bV0 @ bexp / (np.linalg.norm(bV0) * np.linalg.norm(bexp)))
route_gate = float(np.linalg.norm(bV1 - bexp) / np.linalg.norm(bexp))
route_cos = float(bV1 @ bexp / (np.linalg.norm(bV1) * np.linalg.norm(bexp)))
log("route_V0 self-gate vs exported b_TC: relL2=%.3e cos=%.6f (old defective route)"
    % (r0, c0))
log("route_V1 self-gate vs exported b_TC: relL2=%.3e cos=%.6f (production mirror)"
    % (route_gate, route_cos))

# ---------------- C_field (production formula + mixed-BC correction) ------
def dir_deriv(D6, ex, ey, ez):
    return (ex * ex * D6[:, 0] + ex * ey * D6[:, 1] + ex * ez * D6[:, 2]
            + ey * ex * D6[:, 1] + ey * ey * D6[:, 3] + ey * ez * D6[:, 4]
            + ez * ex * D6[:, 2] + ez * ey * D6[:, 4] + ez * ez * D6[:, 5])

dvecC = C[nei] - C[owner]
ef = dvecC / np.maximum(np.linalg.norm(dvecC, axis=1), 1e-30)[:, None]
D6o = dDTDxh[owner]; D6n = dDTDxh[nei]
dd_o = dir_deriv(D6o, ef[:, 0], ef[:, 1], ef[:, 2])
dd_n = dir_deriv(D6n, ef[:, 0], ef[:, 1], ef[:, 2])
base = (Tb[owner] - Tb[nei]) * (T[owner] - T[nei]) * dc * magSfI
Cfield = np.zeros(N)
np.add.at(Cfield, owner, -w * dd_o * base)
np.add.at(Cfield, nei, -(1.0 - w) * dd_n * base)
# boundary
bCtrs = fCtrs[nIF:]
bDvec = bCtrs - C[bc_cell]
bef = bDvec / np.maximum(np.linalg.norm(bDvec, axis=1), 1e-30)[:, None]
dd_b = dir_deriv(dDTDxh[bc_cell], bef[:, 0], bef[:, 1], bef[:, 2])
# patch-aware snGrad (corrected mixed-BC formula; zeroGradient -> 0 skip)
snb = np.zeros(nB)
T_bvals = np.zeros(nB)
for pn, sz in zip(patch_names, patch_sizes):
    o0 = bface_start[pn] - nIF
    bv = T_bval[pn]
    if isinstance(bv, float): T_bvals[o0:o0 + sz] = bv
    elif bv is not None: T_bvals[o0:o0 + sz] = bv
    else: T_bvals[o0:o0 + sz] = T[bc_cell[o0:o0 + sz]]
for pn, sz in zip(patch_names, patch_sizes):
    o0 = bface_start[pn] - nIF
    sl = slice(o0, o0 + sz)
    t_ = T_btype[pn]
    if t_ == "fixedValue":
        snb[sl] = (T_bvals[sl] - T[bc_cell[sl]]) * bdel[sl]
    elif t_ == "mixed":
        vf, rv, rg = mixed_params[pn]
        snb[sl] = (vf * (rv - T[bc_cell[sl]]) * bdel[sl]
                   + (1.0 - vf) * (rg if rg else 0.0))
    else:
        snb[sl] = 0.0
np.add.at(Cfield, bc_cell, dd_b * Tb[bc_cell] * snb * bmag)
log("C_field: maxabs=%.6e ; |internal|=%.6e ; |boundary|=%.6e (mixed-corrected)"
    % (np.abs(Cfield).max(), np.linalg.norm(Cfield * (np.arange(N) < 0)),
       np.linalg.norm(Cfield)))
nzc = int((np.abs(Cfield) > 0).sum())
log("C_field nnz=%d" % nzc)

# ---------------- Gx (production dHbyA/drAU + dphi formula) ----------------
def grad_p(p_):
    pf = w * p_[owner] + (1 - w) * p_[nei]
    g0 = np.zeros((N, 3))
    np.add.at(g0, owner, SfI * pf[:, None])
    np.add.at(g0, nei, -(SfI * pf[:, None]))
    pb_ = p_[bc_cell].copy()
    out_start = bface_start["outlet"] - nIF
    pb_[out_start:out_start + patch_sizes[patch_names.index("outlet")]] = 0.0
    np.add.at(g0, bc_cell, bSf * pb_[:, None])
    return g0 / V[:, None]

mob = bs["mob"]
src = (bs["D_rel"] - bs["diag"])[:, None] * U_B
offU = np.zeros((N, 3))
np.add.at(offU, owner, bs["upper"][:, None] * U_B[nei])
np.add.at(offU, nei, bs["lower"][:, None] * U_B[owner])
HbyA_base = mob[:, None] * (src - offU) / V[:, None] - mob[:, None] * grad_p(p_B)
drAU = -mob * mob / ALPHA_REL
dHbyA = (mob / ALPHA_REL)[:, None] * ((1.0 - ALPHA_REL) * U_B - HbyA_base)
g0_int = magSfI * dc * (p_B[nei] - p_B[owner])

def Gx_of(zdir):
    wc = dAlphaDxh * zdir
    dHw = dHbyA * wc[:, None]
    phx_int = np.einsum("fi,fi->f", SfI,
                        w[:, None] * dHw[owner] + (1 - w)[:, None] * dHw[nei])
    gf_int = w * drAU[owner] * wc[owner] + (1 - w) * drAU[nei] * wc[nei]
    phx_int -= gf_int * g0_int
    val = float(np.dot(g[:nIF], phx_int))
    am = bidx  # U-assignable boundary faces (non-fixed-U)
    phx_b = np.einsum("fi,fi->f", bSf[am], dHw[bc_cell[am]])
    val += float(np.dot(g[nIF:][am], phx_b))
    return val

# ---------------- M1 identity ----------------
z = np.loadtxt(E26 + "/stageB6_rxc_z_analytic.mtx").reshape(3, N)
wH7 = np.load("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/"
              "BFINAL-026/cycle-1/b26_wstar_h7.npz")
wPre = np.load("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/"
               "BFINAL-021/cycle-1/b21_wstar_cache.npz")
rows = []
for d, zi in zip(("D1", "D2", "D3"), z):
    w_ = wH7[d]
    lhs = float(bexp @ w_)
    lhs_pre = float(bexp @ wPre[d])
    Ck = float(np.dot(Cfield, zi))
    Gk = Gx_of(zi)
    fd = FD_J[d]
    rhs = fd - Ck - Gk
    ratio = lhs / rhs
    rel = abs(ratio - 1.0)
    rows.append(dict(dir=d, lhs=lhs, rhs=rhs, ratio=ratio, sign_lhs=np.sign(lhs),
                     sign_rhs=np.sign(rhs), reldev=rel,
                     C=Ck, Gx=Gk, FD=fd,
                     lhs_pre=lhs_pre,
                     lhs_pre_frac=lhs_pre / lhs if lhs != 0 else float("nan"),
                     lhs_diff_frac=abs(lhs - lhs_pre) / abs(lhs_pre),
                     share_anchor=PRODH7SHARE))
    log("[M1] %s: lhs=%.8e rhs=%.8e ratio=%.5f sign=%+d/%+d reldev=%.4f "
        "| C=%.6e Gx=%.6e FD=%.8e | lhs_pre=%.8e"
        % (d, lhs, rhs, ratio, np.sign(lhs), np.sign(rhs), rel, Ck, Gk, fd, lhs_pre))

# verdict vs discrimination matrix (<=10% source layer clean)
for r in rows:
    v = "SOURCE-CLEAN (lesion in lambda_T solve / thermalC contraction a/b)"
    if r["reldev"] > 0.30:
        v = "O(1)-MISMATCH (lesion in T-elimination folding c segment)"
    elif r["reldev"] > 0.10:
        v = "BORDERLINE (between a/b and c)"
    r["verdict"] = v
    log("[M1] %s verdict: %s" % (r["dir"], v))

json.dump(dict(rows=rows, gates=dict(mob=mob_gate, V=V_gate,
                                     route_V1=route_gate, route_V1_cos=route_cos,
                                     route_V0=r0, route_V0_cos=c0),
               FD_J=FD_J, FD_MAINREF=FD_MAINREF, ADJ_J=ADJ_J, PRODH7SHARE=PRODH7SHARE),
          open(OUT + "/b26a_m1_results.json", "w"), indent=1)

# tsv for reviewer recomputation
with open(OUT + "/b26a_m1_identity.tsv", "w") as f:
    f.write("dir\tlhs\tFD\tC\tGx\trhs=FD-C-Gx\tlhs/rhs\tsign_lhs\tsign_rhs"
            "\treldev|lhs/rhs-1|\tlhs_pre\tlhs_pre/lhs\t|lhs-lhs_pre|/|lhs_pre|"
            "\tH7share_anchor\tverdict\n")
    for r in rows:
        f.write("%s\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.8f\t%d\t%d\t%.6f"
                "\t%.10e\t%.6f\t%.6f\t%.6f\t%s\n"
                % (r["dir"], r["lhs"], r["FD"], r["C"], r["Gx"], r["rhs"],
                   r["ratio"], int(r["sign_lhs"]), int(r["sign_rhs"]),
                   r["reldev"], r["lhs_pre"], r["lhs_pre_frac"],
                   r["lhs_diff_frac"], r["share_anchor"], r["verdict"]))
f2 = open(OUT + "/b26a_m1_gates.tsv", "w")
f2.write("gate\tvalue\nmobility_relL2\t%.6e\nV_pyramid_relL2\t%.6e\n"
         "route_V1_vs_b_TC_relL2\t%.6e\nroute_V1_cos\t%.8f\n"
         "route_V0_vs_b_TC_relL2\t%.6e\nroute_V0_cos\t%.8f\n"
         % (mob_gate, V_gate, route_gate, route_cos, r0, c0))
f2.close()
log("\nDONE t=%.1fs" % (time.time() - t0))
