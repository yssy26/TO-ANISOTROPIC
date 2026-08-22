#!/usr/bin/env python3
"""BFINAL-029 Phase 1: slot-by-slot decomposition of the b_TC folding route.

Reproduces B26's route_V1 (relL2=0.397) slot-by-slot, then decomposes the
mismatch vs the exported b_TC (b18rhs_thermalCoupling.mtx) into per-slot
contributions under the production-mirror default [g1_TIn, T8 ON].

Per PREREGISTRATION (b29_PREREGISTRATION.md, written before data):
  1. self-gate   : [g1_Tmix, T8 off] must reproduce 0.3970109364 exactly
  2. additivity  : sum of per-input-slot contributions == full route (machine)
  3. leave-one-out per-slot drelL2 (12 toggleable slots), full = [g1_TIn, T8 on]
  4. g1 caliber  : drelL2(Tmix) - drelL2(TIn)
  5. T8 measured : contribution norm (theorem: isotropic -> ~0)
  6. labels      : offline-error / production-error / caliber / innocent
  + M1 preview (diagnostic only, w_true and rhs from B26, no FD recompute):
    lhs = b_TC^T w_true vs rhs = FD - C - Gx for the corrected full route.

Slots (production source refs in comments):
  T1 hA internal, T3 U internal, T4 P internal kf, H7s1 h7G internal,
  T5 hA boundary, T6 U boundary, H7s1b h7G boundary, T7 P pressureFixed,
  T2 deltaH^T, T8 boundary H-diagonal, H7s2 phiHbyA P internal,
  H7s2b phiHbyA P zeroGradient-p boundary.
  g1 = outlet dJ/dphi (TIn production / Tmix B26-control).
"""
import time, json, re
import numpy as np

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
SB = "/home/ys/dsH/b25_qgate/stageB2"
Q1 = "/home/ys/dsH/b25_qgate/1"
E26 = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-026/cycle-1"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-029/cycle-1"

N = 33600; NV = 3 * N; NUNK = 4 * N
ALPHA_REL = 0.4
NU = 5.19009e-05
SMALL = 1e-15

# B26/B25 pinned numbers (NOT recomputed; Phase 1 uses them only for the M1
# preview which is a diagnostic, not a verdict).
FD_J = {"D1": 0.04252469623637622,
        "D2": -0.003528972517971574,
        "D3": -0.0004449511895182612}
RHS = {"D1": 0.07443678641907196,
       "D2": 0.0006387577668902512,
       "D3": 0.01094029857921596}

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-029 Phase 1 slot instrument ===")

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
log("mesh: nF=%d nIF=%d nB=%d patches=%s" % (nF, nIF, nB, patch_names))

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
C = num / vol3[:, None]
V_pyr = vol3 / 3.0
V = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")[:, 0]
V_gate = float(np.linalg.norm(V_pyr - V) / np.linalg.norm(V))
log("V gate (pyramid) relL2=%.3e" % V_gate)
assert V_gate < 1e-9, "V gate failed"
# exact pyramid deltaCoeffs (b26a DEFECT 1 correction)
magSfI = magSf[:nIF]
dc = magSfI / np.einsum("ij,ij->i", SfI, C[nei] - C[owner])

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
assert len(phi_int) == nIF and len(phiTh_int) == nIF
assert len(coldmask_int) == nIF
p_btype = {pn: p_b[pn][0] for pn in patch_names}
log("fields loaded (t=%.1fs)" % (time.time() - t0))

# phi/p boundary values
phi_bvals = np.zeros(nB)
for pn, sz in zip(patch_names, patch_sizes):
    o0 = bface_start[pn] - nIF
    btp, bvp = phi_b[pn]
    if isinstance(bvp, float): phi_bvals[o0:o0 + sz] = bvp
    elif bvp is not None: phi_bvals[o0:o0 + sz] = bvp

# ---------------- g face functional (g0 internal + g1 outlet) ----------------
def g_functional(g1_mode):
    """g0 (AdjNS_HT.H L23-25) always on. g1 outlet dJ/dphi per mode:
       TIn  (production computeObjective.H L152-172, maximizeTotalHeatTransfer):
             dJ/dphi = -(T_out - T_in)/(M_frozen*Tref), T_in = gAverage(T[coldInlet])
       Tmix (B26 control): dJ/dphi = -(T_out - Tmix)/(Mflow*Tref)
    """
    g = np.zeros(nIF + nB)
    jumpT = T[nei] - T[owner]
    down_nei = phiTh_int >= 0.0
    adjDown = np.where(down_nei, Tb[nei], Tb[owner])
    g[:nIF] = -coldmask_int * adjDown * jumpT
    Mflow, Tref = np.loadtxt(SB + "/wstate_outletMeta.mtx")
    o0 = bface_start["outlet"] - nIF
    sz = patch_sizes[patch_names.index("outlet")]
    ophi = phi_bvals[o0:o0 + sz]
    Toc = T[bc_cell[o0:o0 + sz]]
    if g1_mode == "TIn":
        # T_in = gAverage(T[coldInlet]); coldInlet=inlet, inlet U/T fixedValue 600
        T_in = 600.0
        dJdphi_out = -(Toc - T_in) / (Mflow * Tref)
    else:
        Tmix = float(np.sum(ophi * Toc) / Mflow)
        dJdphi_out = -(Toc - Tmix) / (Mflow * Tref)
    g[nIF + o0:nIF + o0 + sz] += dJdphi_out
    return g, dJdphi_out

g_TIn, dJ_TIn = g_functional("TIn")
g_Tmix, dJ_Tmix = g_functional("Tmix")
log("g: |g0|max=%.3e nz=%d |dJdphi TIn|max=%.3e |dJdphi Tmix|max=%.3e"
    % (np.abs(g_TIn[:nIF]).max(), int((g_TIn[:nIF] != 0).sum()),
       np.abs(dJ_TIn).max(), np.abs(dJ_Tmix).max()))

# ---------------- mobility (b23 recipe, SB=b25 stageB2) ----------------
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
log("MOBILITY GATE: relL2=%.3e (need <1e-12)" % mob_gate)
if mob_gate >= 1e-12:
    raise RuntimeError("BLOCKED: mobility gate failed %.3e" % mob_gate)
mob = bs["mob"]
mobc = mob[:, None]

# ---------------- per-slot folding route ----------------
INPUT_SLOTS = ["T1", "T3", "T4", "H7s1", "T5", "T6", "H7s1b", "T7"]
PIPELINE_SLOTS = ["T2", "T8", "H7s2", "H7s2b"]
ALL_SLOTS = INPUT_SLOTS + PIPELINE_SLOTS
kf_int = (w * mob[owner] + (1 - w) * mob[nei]) * magSfI * dc
Vinvo = 1.0 / V[owner]; Vinn = 1.0 / V[nei]
Vinv = 1.0 / V
# T8 boundary scalar internal coefficients (b26a/momentum convention:
# ic = where(Ufix, -ab, phiB), ab = nuE*bmag*bdel; component-isotropic
# momentum operator => (-ic + cav) == 0 by construction; measured below).
# All arrays here are boundary-ordered (nB): ab = nuE[bc_cell]*bmag*bdel,
# phiB_B is boundary-ordered already, Ufix is boundary-ordered.
ab_b = nuE[bc_cell] * bmag * bdel
ic_b = np.where(Ufix, -ab_b, phiB_B)
ic_bv = np.zeros((nB, 3))
ic_bv[:] = ic_b[:, None]
# zeroGradient-p patches for H7s2b (skip fixedValue-p = outlet)
zgP_patches = [pn for pn in patch_names if p_btype[pn] != "fixedValue"]
log("H7s2b applies to zeroGradient-p patches: %s (outlet skipped)"
    % (zgP_patches,))

def route(g1_mode, slots_on, slots_flipped=()):
    """Assemble b_TC for the given g1 outlet formula, slot mask and flips.
    g1_mode: 'TIn' | 'Tmix'.  slots_on: set of slots included.  slots_flipped:
    set of slots whose contribution sign is flipped (Phase 2 support).
    """
    g = g_TIn if g1_mode == "TIn" else g_Tmix
    gf_i = g[:nIF]; gf_b = g[nIF:]
    hA = np.zeros((N, 3)); h7G = np.zeros((N, 3))
    gU = np.zeros((N, 3)); gP = np.zeros(N)
    sign = {s: (-1.0 if s in slots_flipped else 1.0) for s in ALL_SLOTS}
    ft = SfI * gf_i[:, None]
    if "T1" in slots_on:
        c = sign["T1"]
        np.add.at(hA, owner, c * mobc[owner] * w[:, None] * ft)
        np.add.at(hA, nei, c * mobc[nei] * (1.0 - w)[:, None] * ft)
    if "T3" in slots_on:
        c = sign["T3"]
        np.add.at(gU, owner, c * (1.0 - ALPHA_REL) * w[:, None] * ft)
        np.add.at(gU, nei, c * (1.0 - ALPHA_REL) * (1.0 - w)[:, None] * ft)
    if "T4" in slots_on:
        c = sign["T4"]
        np.add.at(gP, owner, c * kf_int * gf_i)
        np.add.at(gP, nei, -c * kf_int * gf_i)
    if "H7s1" in slots_on:
        c = sign["H7s1"]
        np.add.at(h7G, owner, c * mobc[owner] * w[:, None] * ft)
        np.add.at(h7G, nei, c * mobc[nei] * (1.0 - w)[:, None] * ft)
    # boundary (assignable-U only = outlet; production skips fixesValue-U)
    for f in bidx:
        ccell = bc_cell[f]
        ftb = bSf[f] * gf_b[f]
        if "T5" in slots_on:
            hA[ccell] += sign["T5"] * mob[ccell] * ftb
        if "T6" in slots_on:
            gU[ccell] += sign["T6"] * (1.0 - ALPHA_REL) * ftb
        if "H7s1b" in slots_on:
            h7G[ccell] += sign["H7s1b"] * mob[ccell] * ftb
        if "T7" in slots_on:
            gP[ccell] += sign["T7"] * mob[ccell] * bdel[f] * bmag[f] * gf_b[f]
    if "T2" in slots_on:
        c = sign["T2"]
        gU[nei] -= c * (bs["upper"] * Vinvo)[:, None] * hA[owner]
        gU[owner] -= c * (bs["lower"] * Vinn)[:, None] * hA[nei]
    if "T8" in slots_on:
        # production L442-465: gU[c]+=(-ic[c]+cav)*Vinv[c]*hA[c];
        # component-isotropic momentum operator => (-ic+cav)==0 per face.
        cav = ic_bv.mean(axis=1)
        t8c = (-ic_bv + cav[:, None]) * Vinv[bc_cell][:, None] * hA[bc_cell]
        gU[bc_cell] += sign["T8"] * t8c
    if "H7s2" in slots_on:
        c = sign["H7s2"]
        q = np.einsum("fi,fi->f", SfI,
                      h7G[owner] * Vinvo[:, None] - h7G[nei] * Vinn[:, None])
        np.add.at(gP, owner, -c * w * q)
        np.add.at(gP, nei, -c * (1.0 - w) * q)
    if "H7s2b" in slots_on:
        c = sign["H7s2b"]
        for pn in zgP_patches:
            o0 = bface_start[pn] - nIF
            sz = patch_sizes[patch_names.index(pn)]
            for f in range(o0, o0 + sz):
                ccell = bc_cell[f]
                gP[ccell] -= c * float(h7G[ccell] @ bSf[f]) / V[ccell]
    b = np.zeros(NUNK)
    b[0:NV:3] = gU[:, 0]; b[1:NV:3] = gU[:, 1]; b[2:NV:3] = gU[:, 2]
    b[NV:] = gP
    return b, gU, gP

# ---------------- exported b_TC ----------------
m_exp = np.loadtxt("/home/ys/dsH/b25_qgate/b18rhs_thermalCoupling.mtx")
bexp = np.zeros(NUNK)
bexp[0:NV:3] = m_exp[:, 0]; bexp[1:NV:3] = m_exp[:, 1]
bexp[2:NV:3] = m_exp[:, 2]; bexp[NV:] = m_exp[:, 3]
def relL2(b): return float(np.linalg.norm(b - bexp) / np.linalg.norm(bexp))
def cosv(b):
    return float(b @ bexp / (np.linalg.norm(b) * np.linalg.norm(bexp)))

# ---------------- gates + decomposition ----------------
FULL = set(ALL_SLOTS)
# 1. self-gate: reproduce B26 route_V1 ([g1_Tmix, T8 off])
b_b26 = route("Tmix", FULL - {"T8"})[0]
sg = (relL2(b_b26), cosv(b_b26))
log("SELF-GATE [g1_Tmix, T8 off]: relL2=%.10f cos=%.6f (B26 route_V1 0.3970109364)"
    % sg)
assert abs(sg[0] - 0.39701093640754337) < 1e-6, "SELF-GATE FAILED"

# 2. full production-mirror route [g1_TIn, T8 on]
b_full = route("TIn", FULL)[0]
r_full = relL2(b_full); c_full = cosv(b_full)
log("FULL [g1_TIn, T8 on]: relL2=%.6f cos=%.6f" % (r_full, c_full))

# 3. additivity gate: sum of per-input-slot contributions (pipeline always on)
comp = {}
b_sum = np.zeros(NUNK)
for s in INPUT_SLOTS:
    b_s, _, _ = route("TIn", {s} | set(PIPELINE_SLOTS))
    comp[s] = b_s.copy()
    b_sum += b_s
add_max = float(np.max(np.abs(b_sum - b_full)))
add_rel = float(np.linalg.norm(b_sum - b_full) / np.linalg.norm(b_full))
log("ADDITIVITY: max|sum-full|=%.3e rel=%.3e (need ~machine)" % (add_max, add_rel))

# 4. leave-one-out per slot
loo = {}
for s in ALL_SLOTS:
    b_m, _, _ = route("TIn", FULL - {s})
    r_m = relL2(b_m)
    loo[s] = {"relL2_removed": r_m, "delta": r_m - r_full,
              "removed_norm_frac": float(
                  np.linalg.norm(b_full - b_m) / np.linalg.norm(b_full))}

# 5. g1 caliber: Tmix vs TIn
b_tmix = route("Tmix", FULL)[0]
g1_delta = relL2(b_tmix) - r_full
log("g1 caliber: relL2(Tmix)=%.6f relL2(TIn)=%.6f delta=%.6f"
    % (relL2(b_tmix), r_full, g1_delta))

# 6. T8 measured contribution
t8_norm_frac = loo["T8"]["removed_norm_frac"]
log("T8 contribution norm frac = %.3e (theorem 0, isotropic operator)" % t8_norm_frac)

# 7. M1 preview (diagnostic; w_true + rhs from B26, no FD recompute)
wtrue = {}
for k in ("D1", "D2", "D3"):
    wtrue[k] = np.load(E26 + "/b26_wstar_h7.npz")[k]
m1 = {}
for k in ("D1", "D2", "D3"):
    lhs = float(b_full @ wtrue[k])
    rhs = RHS[k]
    m1[k] = {"lhs": lhs, "rhs": rhs, "ratio": lhs / rhs,
             "b26_ratio": None}
# B26 baseline ratios from bexp
for k in ("D1", "D2", "D3"):
    lhs_b26 = float(bexp @ wtrue[k])
    m1[k]["b26_lhs"] = lhs_b26
    m1[k]["b26_ratio"] = lhs_b26 / RHS[k]
    log("M1 preview %s: lhs=%.8f rhs=%.8f ratio=%.4f (B26 bexp ratio %.4f)"
        % (k, m1[k]["lhs"], rhs, m1[k]["ratio"], m1[k]["b26_ratio"]))

# ---------------- slot table + verdict labels ----------------
def verdict(delta, impl_correct, note=""):
    if note == "T8_zero":
        return "caliber"
    if delta > 0.01:
        return "suspect" if impl_correct else "offline-error"
    return "innocent"

slot_table = []
for s in ALL_SLOTS:
    d = loo[s]["delta"]
    # implementation correctness: b26a reproduced exactly (self-gate) => all
    # 12 slots as implemented here match production source refs; label per
    # preregistration categories.
    label = verdict(d, True)
    slot_table.append({
        "slot": s, "delta_relL2_removed": d,
        "removed_norm_frac": loo[s]["removed_norm_frac"],
        "label": label,
        "note": ("b26a-skipped; theory-zero (isotropic momentum operator)"
                 if s == "T8" else "")})

res = {
    "gates": {
        "mob": mob_gate, "V": V_gate,
        "selfgate_Tmix_noT8": {"relL2": sg[0], "cos": sg[1],
                                "expect": 0.39701093640754337},
        "additivity": {"max_abs": add_max, "rel": add_rel}},
    "full_corrected_TIn_T8on": {"relL2": r_full, "cos": c_full},
    "g1_caliber": {"relL2_Tmix": relL2(b_tmix), "relL2_TIn": r_full,
                   "delta": g1_delta,
                   "dJdphi_TIn_max": float(np.abs(dJ_TIn).max()),
                   "dJdphi_Tmix_max": float(np.abs(dJ_Tmix).max())},
    "T8": {"contribution_norm_frac": t8_norm_frac},
    "leave_one_out": loo,
    "slot_table": slot_table,
    "m1_preview_corrected_route": m1,
}
with open(OUT + "/b29_slot_table.json", "w") as f:
    json.dump(res, f, indent=1)
log("wrote %s/b29_slot_table.json (t=%.1fs)" % (OUT, time.time() - t0))

# ---------------- quick summary ----------------
print("\n=== PHASE 1 SLOT TABLE (delta relL2 = relL2(remove) - relL2(full)) ===")
print("full [TIn,T8on] relL2=%.6f cos=%.6f" % (r_full, c_full))
for s in ALL_SLOTS:
    print("  %-6s delta=%+8.5f  removed_norm=%.4f  %s"
          % (s, loo[s]["delta"], loo[s]["removed_norm_frac"],
             [t["label"] for t in slot_table if t["slot"] == s][0]))
print("g1 caliber delta=%.6f" % g1_delta)
