#!/usr/bin/env python3
"""BFINAL-030 fingerprint 1 (RHS side): corrected b_TC x corrected w_corr.

Main criterion (preregistered):
    lhs = b_TC_corr^T w_corr  vs  rhs (FD-derived, pinned from B26a)
    rhs = {D1:0.07443678641907196, D2:0.0006387577668902512,
           D3:0.01094029857921596}
    PASS iff |ratio - 1| <= 0.15 on ALL THREE directions simultaneously,
    where b_TC_corr uses the SAME c* as the operator patch on BOTH the T3
    (internal) and T6 (boundary) slots (the operator boundary U<-P direct
    term is part of T3_base, so the RHS must be corrected on both).

Secondary outputs:
    - gDP factor prediction (b20 convention):
        gdp_factor = (bPD @ w_corr) / (bPD @ w_true)
      with the defcheck caveat (bPD^T w_true does NOT reproduce ADJ_gDP:
      0.8178 / 0.3660 / 1.0362; bPD has zero U-block).
    - cross-reference: fp3-style ratios with w_true (uncorrected adjoint).
"""
import time, json
import numpy as np

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
SB = "/home/ys/dsH/b25_qgate/stageB2"
Q1 = "/home/ys/dsH/b25_qgate/1"
E26 = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-026/cycle-1"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-030/cycle-1"
import sys
CSTAR = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0

N = 33600; NV = 3 * N; NUNK = 4 * N
ALPHA_REL = 0.4
NU = 5.19009e-05
RHS = {"D1": 0.07443678641907196,
       "D2": 0.0006387577668902512,
       "D3": 0.01094029857921596}

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-030 fp1 analyze  c*=%.3f ===" % CSTAR)

# ---------------- field parsing (b29) ----------------
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
import re
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
bc = oa[nIF:]
nB = nF - nIF
bfile = open(POLY + "/boundary").read()
patch_names = re.findall(r"^\s{4}(\w+)\s*$", bfile, re.M)
patch_sizes = [int(m) for m in re.findall(r"nFaces\s+(\d+);", bfile)]
bface_start = {}
o = nIF
for pn, sz in zip(patch_names, patch_sizes):
    bface_start[pn] = o; o += sz

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
magSfI = magSf[:nIF]
dc = magSfI / np.einsum("ij,ij->i", SfI, C[nei] - C[owner])

w = np.loadtxt(MESH + "/b16mesh_weights.mtx")
bnd = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
bc_cell = bnd[:, 0].astype(np.int64); bSf = bnd[:, 1:4]
bmag = bnd[:, 4]; bdel = bnd[:, 5]
Ufix = bnd[:, 7] > 0.5
bidx = np.where(~Ufix)[0]

# ---------------- fields ----------------
U_i, _, U_b = read_field(Q1 + "/U"); U = U_i.reshape(N, 3)
p_i, _, p_b = read_field(Q1 + "/p"); p = p_i
phi_i, _, phi_b = read_field(Q1 + "/phi"); phi_int = phi_i
phTh_i, _, phTh_b = read_field(Q1 + "/phiThermal"); phiTh_int = phTh_i
cfm_i, _, _ = read_field(Q1 + "/coldFaceMask"); coldmask_int = cfm_i
Tb_i, _, Tb_b = read_field(Q1 + "/Tb"); Tb = Tb_i
T_i, _, T_b = read_field(Q1 + "/T"); T = T_i
alpha_i, _, _ = read_field(Q1 + "/alpha"); alpha = alpha_i
p_btype = {pn: p_b[pn][0] for pn in patch_names}
log("fields loaded (t=%.1fs)" % (time.time() - t0))

phi_bvals = np.zeros(nB)
for pn, sz in zip(patch_names, patch_sizes):
    o0 = bface_start[pn] - nIF
    btp, bvp = phi_b[pn]
    if isinstance(bvp, float): phi_bvals[o0:o0 + sz] = bvp
    elif bvp is not None: phi_bvals[o0:o0 + sz] = bvp

def g_functional(g1_mode):
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
        T_in = 600.0
        dJdphi_out = -(Toc - T_in) / (Mflow * Tref)
    else:
        Tmix = float(np.sum(ophi * Toc) / Mflow)
        dJdphi_out = -(Toc - Tmix) / (Mflow * Tref)
    g[nIF + o0:nIF + o0 + sz] += dJdphi_out
    return g

g_TIn = g_functional("TIn")

# ---------------- mobility (gated) ----------------
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
mob = bs["mob"]; mobc = mob[:, None]

# ---------------- slots + route ----------------
INPUT_SLOTS = ["T1", "T3", "T4", "H7s1", "T5", "T6", "H7s1b", "T7"]
PIPELINE_SLOTS = ["T2", "T8", "H7s2", "H7s2b"]
ALL_SLOTS = INPUT_SLOTS + PIPELINE_SLOTS
kf_int = (w * mob[owner] + (1 - w) * mob[nei]) * magSfI * dc
Vinvo = 1.0 / V[owner]; Vinn = 1.0 / V[nei]
Vinv = 1.0 / V
ab_b = nuE[bc_cell] * bmag * bdel
ic_b = np.where(Ufix, -ab_b, phiB_B)
ic_bv = np.zeros((nB, 3)); ic_bv[:] = ic_b[:, None]
zgP_patches = [pn for pn in patch_names if p_btype[pn] != "fixedValue"]

def route(g1_mode, slots_on, t3coef=0.6, t6coef=0.6):
    g = g_TIn if g1_mode == "TIn" else g_Tmix
    gf_i = g[:nIF]; gf_b = g[nIF:]
    hA = np.zeros((N, 3)); h7G = np.zeros((N, 3))
    gU = np.zeros((N, 3)); gP = np.zeros(N)
    ft = SfI * gf_i[:, None]
    if "T1" in slots_on:
        np.add.at(hA, owner, mobc[owner] * w[:, None] * ft)
        np.add.at(hA, nei, mobc[nei] * (1.0 - w)[:, None] * ft)
    if "T3" in slots_on:
        np.add.at(gU, owner, t3coef * w[:, None] * ft)
        np.add.at(gU, nei, t3coef * (1.0 - w)[:, None] * ft)
    if "T4" in slots_on:
        np.add.at(gP, owner, kf_int * gf_i)
        np.add.at(gP, nei, -kf_int * gf_i)
    if "H7s1" in slots_on:
        np.add.at(h7G, owner, mobc[owner] * w[:, None] * ft)
        np.add.at(h7G, nei, mobc[nei] * (1.0 - w)[:, None] * ft)
    for f in bidx:
        ccell = bc_cell[f]
        ftb = bSf[f] * gf_b[f]
        if "T5" in slots_on:
            hA[ccell] += mob[ccell] * ftb
        if "T6" in slots_on:
            gU[ccell] += t6coef * ftb
        if "H7s1b" in slots_on:
            h7G[ccell] += mob[ccell] * ftb
        if "T7" in slots_on:
            gP[ccell] += mob[ccell] * bdel[f] * bmag[f] * gf_b[f]
    if "T2" in slots_on:
        gU[nei] -= (bs["upper"] * Vinvo)[:, None] * hA[owner]
        gU[owner] -= (bs["lower"] * Vinn)[:, None] * hA[nei]
    if "T8" in slots_on:
        cav = ic_bv.mean(axis=1)
        t8c = (-ic_bv + cav[:, None]) * Vinv[bc_cell][:, None] * hA[bc_cell]
        gU[bc_cell] += t8c
    if "H7s2" in slots_on:
        q = np.einsum("fi,fi->f", SfI,
                      h7G[owner] * Vinvo[:, None] - h7G[nei] * Vinn[:, None])
        np.add.at(gP, owner, -w * q)
        np.add.at(gP, nei, -(1.0 - w) * q)
    if "H7s2b" in slots_on:
        for pn in zgP_patches:
            o0 = bface_start[pn] - nIF
            sz = patch_sizes[patch_names.index(pn)]
            for f in range(o0, o0 + sz):
                ccell = bc_cell[f]
                gP[ccell] -= float(h7G[ccell] @ bSf[f]) / V[ccell]
    b = np.zeros(NUNK)
    b[0:NV:3] = gU[:, 0]; b[1:NV:3] = gU[:, 1]; b[2:NV:3] = gU[:, 2]
    b[NV:] = gP
    return b

FULL = set(ALL_SLOTS)

# ---------------- adjoint states ----------------
wcorr = {}
for k in ("D1", "D2", "D3"):
    wcorr[k] = np.load(OUT + "/b30_fp1_wcorr_c%.2f.npz" % CSTAR)[k]
wtrue = {}
for k in ("D1", "D2", "D3"):
    wtrue[k] = np.load(E26 + "/b26_wstar_h7.npz")[k]
log("w_corr and w_true loaded (t=%.1fs)" % (time.time() - t0))

# ---------------- exported b_TC / b_PD ----------------
m_exp = np.loadtxt("/home/ys/dsH/b25_qgate/b18rhs_thermalCoupling.mtx")
bexp = np.zeros(NUNK)
bexp[0:NV:3] = m_exp[:, 0]; bexp[1:NV:3] = m_exp[:, 1]
bexp[2:NV:3] = m_exp[:, 2]; bexp[NV:] = m_exp[:, 3]

m_pd = np.loadtxt("/home/ys/dsH/b25_qgate/b18rhs_pressureDrop.mtx")
bPD = np.zeros(NUNK)
bPD[0:NV:3] = m_pd[:, 0]; bPD[1:NV:3] = m_pd[:, 1]
bPD[2:NV:3] = m_pd[:, 2]; bPD[NV:] = m_pd[:, 3]
log("bPD U-block norm=%.3e total=%.3e (defcheck: U-block must be ~0)"
    % (np.linalg.norm(bPD[:NV]), np.linalg.norm(bPD)))

res = {}

# ---------------- MAIN criterion: corrected b_TC x corrected w_corr --------
for label, t6c in (("T3corr_T6corr", CSTAR), ("T3corr_T6keep", 0.6)):
    b = route("TIn", FULL, CSTAR, t6c)
    ratios = {}
    for k in ("D1", "D2", "D3"):
        lhs = float(b @ wcorr[k])
        ratios[k] = {"lhs": lhs, "rhs": RHS[k], "ratio": lhs / RHS[k]}
    maxdev = max(abs(ratios[k]["ratio"] - 1.0) for k in ("D1", "D2", "D3"))
    log("MAIN %-14s ratios D1=%.4f D2=%.4f D3=%.4f  maxdev=%.4f"
        % (label, ratios["D1"]["ratio"], ratios["D2"]["ratio"],
           ratios["D3"]["ratio"], maxdev))
    res["main_%s" % label] = ratios

# ---------------- fp3 cross-reference: corrected b_TC x w_true ------------
b = route("TIn", FULL, CSTAR, CSTAR)
fp3x = {}
for k in ("D1", "D2", "D3"):
    lhs = float(b @ wtrue[k])
    fp3x[k] = {"lhs": lhs, "rhs": RHS[k], "ratio": lhs / RHS[k]}
res["fp3x_T3corr_T6corr_x_wtrue"] = fp3x
log("FP3X (corr b_TC x w_true) D1=%.4f D2=%.4f D3=%.4f"
    % (fp3x["D1"]["ratio"], fp3x["D2"]["ratio"], fp3x["D3"]["ratio"]))

# ---------------- gDP prediction (b20 convention) ----------------
gdp = {}
ADJ_gDP = {"D1": -5.23029384626718, "D2": 1.0505249811003, "D3": 13.60342258150915}
FD_gDP = {"D1": -2.53745977415909, "D2": 0.5234094521754384, "D3": 6.477215547197601}
for k in ("D1", "D2", "D3"):
    num = float(bPD @ wcorr[k])
    den = float(bPD @ wtrue[k])
    gdp[k] = {
        "bPD_wcorr": num, "bPD_wtrue": den,
        "gdp_factor": num / den if den != 0 else None,
        "defcheck_ratio_bPD_wtrue_ADJ": den / ADJ_gDP[k],
    }
res["gdp"] = gdp
log("GDP factor (bPD@wcorr)/(bPD@wtrue): D1=%.4f D2=%.4f D3=%.4f"
    % (gdp["D1"]["gdp_factor"], gdp["D2"]["gdp_factor"],
       gdp["D3"]["gdp_factor"]))
log("defcheck bPD@wtrue/ADJ_gDP: D1=%.4f D2=%.4f D3=%.4f  (preregistered"
    " target: production ADJ/FD 2.05; prediction: 2.05->1.0 => factor ~0.49)"
    % (gdp["D1"]["defcheck_ratio_bPD_wtrue_ADJ"],
       gdp["D2"]["defcheck_ratio_bPD_wtrue_ADJ"],
       gdp["D3"]["defcheck_ratio_bPD_wtrue_ADJ"]))

# ---------------- gDP via production-consistent b24 formula ----------------
# b20 convention gdp_factor=(bPD@x)/(bPD@wt); also report the raw ADJ-like
# value using the same bPD normalization as production ADJ_gDP denominator:
# bPD@wtrue reproduces ADJ_gDP only to 0.37-1.04 (defcheck), so factor is the
# meaningful relative prediction.
res["cstar"] = CSTAR
res["preregistered_criterion"] = (
    "MAIN: |ratio-1|<=0.15 all three (b_TC_corr^T w_corr vs RHS); "
    "gDP: ADJ/FD 2.05 -> 1.0 +- 0.15 (i.e. gdp_factor ~ 1/2.05 ~= 0.49)")

with open(OUT + "/b30_fp1_analyze_c%.2f.json" % CSTAR, "w") as f:
    json.dump(res, f, indent=1)
log("wrote %s/b30_fp1_analyze_c%.2f.json (t=%.1fs)" % (OUT, CSTAR, time.time() - t0))
