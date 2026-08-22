#!/usr/bin/env python3
"""BFINAL-029 Phase 1 state-mismatch test.

The exported b18rhs_thermalCoupling.mtx is the STAGE B2 baseline export
(mtime 12:57:06; scale-factor arithmetic matches Stage B2 PRODRHS).  The b29
slot instrument reconstructs b at the /1 (main-loop) state: g-layer fields
(T, Tb, phiThermal, phi, coldFaceMask) from /1, basis/mobility from Stage B2.
This script tests how much of the relL2=0.386174 residual is explained by the
state difference between /1 (main loop) and Stage B2 (full-SST converged).

Decisive cheap test first: g1 outlet dJ/dphi at the Stage B2 state must
reproduce the Stage B2 B0.2 diagnostic maxObjectiveFluxDerivative=108.815961611
(log line 6811), exactly as the /1-state g1 reproduces 150.300586467 (log 940).

Then: rebuild the full route with Stage B2 baseline fields (T, phi, phiB,
alpha, nuEffFrozen, nutFrozen, primalPressureMobility) plus the reconstructed
Stage B2 phiThermal = coldFaceMask*phi_B + hotRhoCpRatio*hotFaceMask*phiHotFrozen
(recipe validated to 2e-13 at /1).  Tb at the Stage B2 state is NOT exported,
so g0 uses the /1 Tb as the sole proxy (limitation reported).

Field-swap ladder: start from /1 state (relL2=0.386174), swap each field
group to Stage B2 one at a time, and observe relL2 vs the exported b_TC.
"""
import time, json, re
import numpy as np

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
SB = "/home/ys/dsH/b25_qgate/stageB2"
Q1 = "/home/ys/dsH/b25_qgate/1"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-029/cycle-1"
B25 = "/home/ys/dsH/b25_qgate"

N = 33600; NV = 3 * N; NUNK = 4 * N
ALPHA_REL = 0.4
NU = 5.19009e-05
HOT_RHO_CP = 0.425273056469

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-029 Phase 1 state-mismatch test ===")

def read_field(path):
    txt = open(path).read()
    i = txt.index("internalField")
    nl = txt.find("\n", i)
    is_uniform = txt[i:nl].strip().split()[1] == "uniform"
    if is_uniform:
        seg_end = txt.find(";", i)
        u = txt[txt.index("uniform", i) + len("uniform"):seg_end].strip()
        toks = []
        for tok in u.replace("(", " ").replace(")", " ").split():
            try: toks.append(float(tok))
            except ValueError: break
        val = np.array(toks) if len(toks) > 1 else float(toks[0])
        internal = None
    else:
        j = txt.index("\n(", i)
        k = txt.index("\n)", j)
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
bc = oa[nIF:]
nB = nF - nIF
bfile = open(POLY + "/boundary").read()
patch_names = re.findall(r"^\s{4}(\w+)\s*$", bfile, re.M)
patch_sizes = [int(m) for m in re.findall(r"nFaces\s+(\d+);", bfile)]
bface_start = {}
o = nIF
for pn, sz in zip(patch_names, patch_sizes):
    bface_start[pn] = o; o += sz
log("mesh: nF=%d nIF=%d nB=%d patches=%s" % (nF, nIF, nB, patch_names))

# --- exact pyramid centroids ---
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
V = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")[:, 0]
magSfI = magSf[:nIF]
dc = magSfI / np.einsum("ij,ij->i", SfI, C[nei] - C[owner])

w = np.loadtxt(MESH + "/b16mesh_weights.mtx")
bnd = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
bc_cell = bnd[:, 0].astype(np.int64); bSf = bnd[:, 1:4]
bmag = bnd[:, 4]; bdel = bnd[:, 5]
Ufix = bnd[:, 7] > 0.5
bidx = np.where(~Ufix)[0]

# ---------------- /1 fields ----------------
U1_i, _, U1_b = read_field(Q1 + "/U"); U1 = U1_i.reshape(N, 3)
p1_i, _, p1_b = read_field(Q1 + "/p"); p1 = p1_i
phi1_i, _, phi1_b = read_field(Q1 + "/phi"); phi1 = phi1_i
phTh1_i, _, phTh1_b = read_field(Q1 + "/phiThermal"); phTh1 = phTh1_i
cfm1_i, _, _ = read_field(Q1 + "/coldFaceMask"); cfm1 = cfm1_i
hfm1_i, _, _ = read_field(Q1 + "/hotFaceMask"); hfm1 = hfm1_i
phihf1_i, _, _ = read_field(Q1 + "/phiHotFrozen"); phihf1 = phihf1_i
Tb1_i, _, Tb1_b = read_field(Q1 + "/Tb"); Tb1 = Tb1_i
T1_i, _, T1_b = read_field(Q1 + "/T"); T1 = T1_i
alpha1_i, _, _ = read_field(Q1 + "/alpha"); alpha1 = alpha1_i
nuE1_raw, _, _ = read_field(Q1 + "/nuEffFrozen")
nuE1 = nuE1_raw if nuE1_raw is not None else np.full(N, NU)
if nuE1_raw is None:
    log("nuEffFrozen at /1 is uniform = %.6e (NU)" % NU)
# phi / p boundary values (/1)
phi_bvals1 = np.zeros(nB)
for pn, sz in zip(patch_names, patch_sizes):
    o0 = bface_start[pn] - nIF
    btp, bvp = phi1_b[pn]
    if isinstance(bvp, float): phi_bvals1[o0:o0 + sz] = bvp
    elif bvp is not None: phi_bvals1[o0:o0 + sz] = bvp
p_btype = {pn: p1_b[pn][0] for pn in patch_names}

# ---------------- Stage B2 fields ----------------
phi_B = np.loadtxt(SB + "/wstate_baseline_phi.mtx")
phiB_B = np.loadtxt(SB + "/wstate_baseline_phiB.mtx")
T_B = np.loadtxt(SB + "/wstate_baseline_T.mtx")
p_B = np.loadtxt(SB + "/wstate_baseline_p.mtx")
U_B = np.loadtxt(SB + "/wstate_baseline_U.mtx")
al_B = np.loadtxt(SB + "/wstate_baseline_alpha.mtx")
nut_B = np.loadtxt(SB + "/wstate_baseline_nutFrozen.mtx")
nuE_B = NU + nut_B
mob_ex = np.loadtxt(SB + "/wstate_baseline_primalPressureMobility.mtx")
phTh_B = cfm1 * phi_B + HOT_RHO_CP * hfm1 * phihf1
log("Stage B2 phiThermal reconstructed: relL2 vs /1 phiThermal = %.3e"
    % (np.linalg.norm(phTh_B - phTh1) / np.linalg.norm(phTh1)))
log("phi_B len=%d phiB_B len=%d T_B len=%d (nIF=%d nB=%d)"
    % (len(phi_B), len(phiB_B), len(T_B), nIF, nB))
assert len(phi_B) == nIF and len(phiB_B) == nB and len(T_B) == N

# ---------------- mobility (per state) ----------------
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

bsB = build_basis(phi_B, phiB_B, al_B, nuE_B)
mob_gate = float(np.linalg.norm(bsB["mob"] - mob_ex) / np.linalg.norm(mob_ex))
log("MOBILITY GATE (B2): relL2=%.3e (need <1e-12)" % mob_gate)
if mob_gate >= 1e-12:
    raise RuntimeError("BLOCKED: mobility gate failed %.3e" % mob_gate)
bs1 = build_basis(phi1, phi_bvals1, alpha1, nuE1)

Vinvo = 1.0 / V[owner]; Vinn = 1.0 / V[nei]; Vinv = 1.0 / V
zgP_patches = [pn for pn in patch_names if p_btype[pn] != "fixedValue"]

# ---------------- g face functional (state-parameterized) ----------------
def g_functional(T, phiTh, Tb, phi_bvals, g1_mode="TIn"):
    g = np.zeros(nIF + nB)
    jumpT = T[nei] - T[owner]
    down_nei = phiTh >= 0.0
    adjDown = np.where(down_nei, Tb[nei], Tb[owner])
    g[:nIF] = -cfm1 * adjDown * jumpT
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
    return g, dJdphi_out

# ---------------- full route, state-parameterized ----------------
def route(f):
    """f: dict with T, phiTh, Tb, phi_bvals, bs (basis dict: lower/upper/mob),
    nuEff, phiB_b (boundary phi for ic_b)."""
    g, _ = g_functional(f["T"], f["phiTh"], f["Tb"], f["phi_bvals"])
    gf_i = g[:nIF]; gf_b = g[nIF:]
    bs = f["bs"]; mob = bs["mob"]; mobc = mob[:, None]
    nuEff = f["nuEff"]
    phiB_b = f["phiB_b"]
    hA = np.zeros((N, 3)); h7G = np.zeros((N, 3))
    gU = np.zeros((N, 3)); gP = np.zeros(N)
    kf_int = (w * mob[owner] + (1 - w) * mob[nei]) * magSfI * dc
    ab_b = nuEff[bc_cell] * bmag * bdel
    ic_b = np.where(Ufix, -ab_b, phiB_b)
    ic_bv = np.zeros((nB, 3)); ic_bv[:] = ic_b[:, None]
    ft = SfI * gf_i[:, None]
    # internal input slots
    np.add.at(hA, owner, mobc[owner] * w[:, None] * ft)
    np.add.at(hA, nei, mobc[nei] * (1.0 - w)[:, None] * ft)
    np.add.at(gU, owner, (1.0 - ALPHA_REL) * w[:, None] * ft)
    np.add.at(gU, nei, (1.0 - ALPHA_REL) * (1.0 - w)[:, None] * ft)
    np.add.at(gP, owner, kf_int * gf_i)
    np.add.at(gP, nei, -kf_int * gf_i)
    np.add.at(h7G, owner, mobc[owner] * w[:, None] * ft)
    np.add.at(h7G, nei, mobc[nei] * (1.0 - w)[:, None] * ft)
    # boundary input slots (assignable-U = outlet only)
    for f in bidx:
        ccell = bc_cell[f]
        ftb = bSf[f] * gf_b[f]
        hA[ccell] += mob[ccell] * ftb
        gU[ccell] += (1.0 - ALPHA_REL) * ftb
        h7G[ccell] += mob[ccell] * ftb
        gP[ccell] += mob[ccell] * bdel[f] * bmag[f] * gf_b[f]
    # pipeline slots
    gU[nei] -= (bs["upper"] * Vinvo)[:, None] * hA[owner]
    gU[owner] -= (bs["lower"] * Vinn)[:, None] * hA[nei]
    cav = ic_bv.mean(axis=1)
    t8c = (-ic_bv + cav[:, None]) * Vinv[bc_cell][:, None] * hA[bc_cell]
    gU[bc_cell] += t8c
    q = np.einsum("fi,fi->f", SfI,
                  h7G[owner] * Vinvo[:, None] - h7G[nei] * Vinn[:, None])
    np.add.at(gP, owner, -w * q)
    np.add.at(gP, nei, -(1.0 - w) * q)
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

# ---------------- exported b_TC ----------------
m_exp = np.loadtxt(B25 + "/b18rhs_thermalCoupling.mtx")
bexp = np.zeros(NUNK)
bexp[0:NV:3] = m_exp[:, 0]; bexp[1:NV:3] = m_exp[:, 1]
bexp[2:NV:3] = m_exp[:, 2]; bexp[NV:] = m_exp[:, 3]
def relL2(b): return float(np.linalg.norm(b - bexp) / np.linalg.norm(bexp))
def cosv(b): return float(b @ bexp / (np.linalg.norm(b) * np.linalg.norm(bexp)))

# ================= TEST 1: g1 at Stage B2 state vs log 108.815961611 =================
dJ_1 = g_functional(T1, phTh1, Tb1, phi_bvals1)[1]
dJ_B = g_functional(T_B, phTh_B, Tb1, phiB_B)[1]
log("\nTEST 1: g1 outlet dJ/dphi, state comparison")
log("  /1 state  : max|dJdphi TIn| = %.12f  (log 940 = 150.300586467)"
    % np.abs(dJ_1).max())
log("  B2 state  : max|dJdphi TIn| = %.12f  (log 6811 = 108.815961611)"
    % np.abs(dJ_B).max())

# ================= TEST 2: full route state ladder =================
# state-mix anchor: B2 basis + /1 g-layer == instrument full route (0.386174)
fMix = dict(T=T1, phiTh=phTh1, Tb=Tb1, phi_bvals=phi_bvals1, bs=bsB,
            nuEff=nuE_B, phiB_b=phiB_B)
f1 = dict(T=T1, phiTh=phTh1, Tb=Tb1, phi_bvals=phi_bvals1, bs=bs1,
          nuEff=nuE1, phiB_b=phi_bvals1)
fB = dict(T=T_B, phiTh=phTh_B, Tb=Tb1, phi_bvals=phiB_B, bs=bsB,
          nuEff=nuE_B, phiB_b=phiB_B)
log("\nTEST 2: field-swap ladder (relL2 vs exported b_TC)")
r_mix = route(fMix)
log("  [state-mix (B2 basis, /1 g-layer)] relL2=%.10f cos=%.6f  (instrument full = 0.386174)"
    % (relL2(r_mix), cosv(r_mix)))
r_all1 = route(f1)
log("  [all /1          ] relL2=%.10f cos=%.6f" % (relL2(r_all1), cosv(r_all1)))
def swp(**kw):
    d = dict(fMix); d.update(kw); return route(d)
for name, kw in [("T -> B2", dict(T=T_B)),
                 ("phiTh -> B2", dict(phiTh=phTh_B)),
                 ("phi bvals -> B2", dict(phi_bvals=phiB_B)),
                 ("Tb -> (proxy)", dict(Tb=Tb1))]:
    b = swp(**kw)
    log("  [swap %-14s] relL2=%.10f cos=%.6f" % (name, relL2(b), cosv(b)))
r_allB = route(fB)
log("  [all B2 (Tb proxy /1)] relL2=%.10f cos=%.6f" % (relL2(r_allB), cosv(r_allB)))

# ================= TEST 3: component norms vs export =================
log("\nTEST 3: component norms")
log("  export |U|L2=5.210372e-04 |P|L2=1.932741e-06 (sumRhsU=3.9845 sumRhsP=-20.6076)")
for nm, b in (("/1", r_all1), ("B2", r_allB)):
    log("  %s route |U|L2=%.6e |P|L2=%.6e sumU=%.6e sumP=%.6e"
        % (nm, np.linalg.norm(b[:NV]), np.linalg.norm(b[NV:]),
           b[:NV].sum(), b[NV:].sum()))

# ================= TEST 4: g0 layer state sensitivity =================
g0_1 = np.linalg.norm(-cfm1 * np.where(phTh1 >= 0, Tb1[nei], Tb1[owner])
                      * (T1[nei] - T1[owner]))
g0_B = np.linalg.norm(-cfm1 * np.where(phTh_B >= 0, Tb1[nei], Tb1[owner])
                      * (T_B[nei] - T_B[owner]))
sel1 = (phTh1 >= 0).astype(float); selB = (phTh_B >= 0).astype(float)
log("\nTEST 4: g0 layer (T3-dominant) state sensitivity")
log("  |g0@/1|L2=%.6e  |g0@B2(T,Tb proxy)|L2=%.6e  ratio=%.4f"
    % (g0_1, g0_B, g0_B / g0_1))
log("  downwind-selection flips between /1 and B2: %d faces of %d cold"
    % (int((sel1 != selB).sum()), int((cfm1 != 0).sum())))

json.dump({"T1_dJdphi_max_1": float(np.abs(dJ_1).max()),
           "T1_dJdphi_max_B2": float(np.abs(dJ_B).max()),
           "log_main_940": 150.300586467,
           "log_B2_6811": 108.815961611,
           "route_relL2_all1": relL2(r_all1),
           "route_cos_all1": cosv(r_all1),
           "route_relL2_allB2": relL2(r_allB),
           "route_cos_allB2": cosv(r_allB),
           "export_U_L2": 5.210372e-04, "export_P_L2": 1.932741e-06,
           "route_U_L2_all1": float(np.linalg.norm(r_all1[:NV])),
           "route_U_L2_allB2": float(np.linalg.norm(r_allB[:NV])),
           "route_P_L2_all1": float(np.linalg.norm(r_all1[NV:])),
           "route_P_L2_allB2": float(np.linalg.norm(r_allB[NV:])),
           "route_sumU_all1": float(r_all1[:NV].sum()),
           "route_sumP_all1": float(r_all1[NV:].sum()),
           "route_sumU_allB2": float(r_allB[:NV].sum()),
           "route_sumP_allB2": float(r_allB[NV:].sum()),
           "g0_L2_1": float(g0_1), "g0_L2_B2": float(g0_B),
           "downwind_flips": int((sel1 != selB).sum())},
          open(OUT + "/b29_state_test.json", "w"), indent=1)
log("\nwrote %s/b29_state_test.json (t=%.1fs)" % (OUT, time.time() - t0))
