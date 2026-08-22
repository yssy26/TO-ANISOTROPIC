#!/usr/bin/env python3
"""BFINAL-032 b32_decomp: same-round decomposition of the gated exports.

Phase-3 deliverable (runs immediately after the b32_decomp run, per B26 lesson).
Loads the round-tagged fields written by the stageB31WriteDecomposition gated
block in sensitivity.H / rxPressureRowTranspose.H and recomputes the B31 §7
accounting with SAME-ROUND data (no more /1/ cross-round directory).

Gates (see PREREGISTRATION.md):
  G1  same-round closure (bit-exact target): gsenshPD_r2 == mom_r2 + prow_r2
      (elementwise relL2 <= 1e-12); also the J-label 4-segment sum pre-chain vs
      chained dfdx.
  G2  T-ring mirror: production rxPressureRowT_r2 (pc contraction) vs offline
      T_row(pc_r2); relL2 < 2% + corr.  rxPressureRowTb_r2 vs T_row(pb_r2).
  G3  cross-round artifact dissolution: same-round mom/prow projections vs B31
      §7 cross-round numbers; chained raw segment sum vs ADJ_gDP (post-chain);
      residual after dissolution = true lesion (or none).
  G4  HbyA calibration: production rxDHbyA_r2/rxrAU_r2/rxG0Field_r2 vs
      coordinator reconstruction dHbyA/mob/g0.
Cross-round attribution evidence: untagged /1/Uc pc pb Tb == _r1 copies
(bit-exact), untagged sensitivity family == r2 chain outputs.
"""
import time, json, re
import numpy as np

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
SB = "/home/ys/dsH/b25_qgate/stageB2"
Q1 = "/home/ys/dsH/b32_decomp/1"
B8 = "/home/ys/dsH/b8_verify_diag"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-032/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N
ALPHAREL = 0.4
NU = 5.19009e-05

# pinned truth values (b25 / B31 §7, NOT recomputed here)
ADJ_J  = {"D1": -0.0358350400485765, "D2": -0.02492599038981423, "D3": -0.01200301226991504}
ADJ_gDP= {"D1": -5.23029384626718,   "D2": 1.0505249811003,     "D3": 13.60342258150915}
FD_J   = {"D1": 0.04252469623637622, "D2": -0.003528972517971574, "D3": -0.0004449511895182612}
bPD_w  = {"D1": -4.277461237487605,  "D2": 0.3844628061126605,  "D3": 14.096353328798937}
bTC_w  = {"D1": -0.04884018700713266,"D2": -0.006440998953954891,"D3": 0.011566617518607991}
C_pinned = {"D1": -0.030341656277042862, "D2": -0.008210889966675959, "D3": -0.01979195745200000}
Gx_pinned= {"D1": -0.0015704339056528767,"D2": 0.004043159681714134,"D3": 0.00840670768310000}
# B31 §7 cross-round numbers (from b31b_m2_run.log) for direct comparison
B31_XR = {  # (mom_z, prow_prod_z, closure, prod_total, true_prow, prow_lesion, total_lesion)
    "D1": (-4.253786, -0.077828, -4.331615, -16.116120, -0.023675, -0.054153, -11.838659),
    "D2": (None, None, None, None, None, None, None),
    "D3": (None, None, None, None, None, None, None),
}
DIRS = ["D1", "D2", "D3"]

t0 = time.time()
def log(m): print("[%5.1fs] %s" % (time.time()-t0, m), flush=True)
log("=== b32_decomp: same-round gated-field decomposition ===")

# ---------------- read_field (b29/b27 validated, uniform-safe) -------------
def read_field(path):
    txt = open(path).read()
    i = txt.index("internalField")
    j = txt.index("\n(", i)
    if txt[i:j].split()[1] == "uniform":
        return None, None, {}
    k = txt.index("\n)", j)
    vals = []
    for line in txt[j+1:k].splitlines():
        for tok in line.replace(";", " ").split():
            for tt in tok.replace("(", " ").replace(")", " ").split():
                vals.append(float(tt))
    internal = np.array(vals)
    bnd = {}
    b = txt.index("boundaryField")
    lines = txt[b:].splitlines()
    p = 0
    while p < len(lines):
        s = lines[p].strip()
        if s.endswith("{") and s != "{":
            name = s[:-1].strip()
            depth, q, btype, bval = 1, p+1, None, None
            while q < len(lines) and depth > 0:
                t = lines[q].strip()
                if t.endswith("{") and t != "{": depth += 1
                if t == "{": depth += 1
                if t == "}": depth -= 1
                if t.startswith("type") and btype is None: btype = t.split()[1].rstrip(";")
                if (t.startswith("uniform") or t.startswith("value")) and depth == 1 and "uniform" in t:
                    u = t[t.index("uniform")+len("uniform"):].rstrip(";").strip()
                    bval = np.array([float(x) for x in u.strip("()").split()]) if u.startswith("(") else float(u)
                q += 1
            bnd[name] = (btype, bval)
        p += 1
    return internal, None, bnd

# ---------------- mesh (b29 validated) --------------------------------------
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
nIF = len(na); nF = len(faces); nB = nF - nIF
Sf = np.zeros((nF, 3))
for fi in range(nF):
    fv = pts[faces[fi]]
    Sf[fi] = 0.5*np.cross(fv, np.roll(fv, -1, axis=0)).sum(axis=0)
magSf = np.linalg.norm(Sf, axis=1)
fCtrs = np.array([pts[fv].mean(axis=0) for fv in faces])
owner = oa[:nIF]; nei = na[:nIF]; bc = oa[nIF:]
bfile = open(POLY + "/boundary").read()
patch_names = re.findall(r"^\s{4}(\w+)\s*$", bfile, re.M)
patch_sizes = [int(m) for m in re.findall(r"nFaces\s+(\d+);", bfile)]
bface_start = {}; o = nIF
for pn, sz in zip(patch_names, patch_sizes):
    bface_start[pn] = o; o += sz
cEst = np.zeros((N, 3)); nFc = np.zeros(N, dtype=np.int64)
np.add.at(cEst, oa, fCtrs); np.add.at(nFc, oa, 1)
np.add.at(cEst, na, fCtrs[:nIF]); np.add.at(nFc, na, 1)
cEst /= nFc[:, None]
num = np.zeros((N, 3)); vol3 = np.zeros(N)
SfI = Sf[:nIF]; fCtrsI = fCtrs[:nIF]; magSfI = magSf[:nIF]
p3o = np.einsum("ij,ij->i", SfI, fCtrsI - cEst[owner])
np.add.at(num, owner, (0.75*fCtrsI + 0.25*cEst[owner])*p3o[:, None]); np.add.at(vol3, owner, p3o)
p3n = np.einsum("ij,ij->i", -SfI, fCtrsI - cEst[nei])
np.add.at(num, nei, (0.75*fCtrsI + 0.25*cEst[nei])*p3n[:, None]); np.add.at(vol3, nei, p3n)
p3b = np.einsum("ij,ij->i", Sf[nIF:], fCtrs[nIF:] - cEst[bc])
np.add.at(num, bc, (0.75*fCtrs[nIF:] + 0.25*cEst[bc])*p3b[:, None]); np.add.at(vol3, bc, p3b)
C = num/vol3[:, None]
V = vol3/3.0
dc = magSfI/np.einsum("ij,ij->i", SfI, C[nei] - C[owner])
w = np.loadtxt(MESH + "/b16mesh_weights.mtx")
bnd = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
bc_cell = bnd[:, 0].astype(np.int64); bSf = bnd[:, 1:4]
bmag = bnd[:, 4]; bdel = bnd[:, 5]
Ufix = bnd[:, 7] > 0.5
am_assign = np.where(~Ufix)[0]
log("mesh ready: nIF=%d nB=%d" % (nIF, nB))

# ---------------- fields ----------------------------------------------------
def scalar(q, fb=None):
    v, _, _ = read_field(Q1 + "/" + q)
    return v
def vector(q):
    v, _, _ = read_field(Q1 + "/" + q)
    return v.reshape(N, 3)
# primal baseline (b25-gated; b32 is bit-identical in the FD columns)
U = np.loadtxt(SB + "/wstate_baseline_U.mtx").reshape(N, 3)
p = np.loadtxt(SB + "/wstate_baseline_p.mtx")
phi_int = np.loadtxt(SB + "/wstate_baseline_phi.mtx")
phi_bvals = np.loadtxt(SB + "/wstate_baseline_phiB.mtx")
alpha = np.loadtxt(SB + "/wstate_baseline_alpha.mtx")
dAlphaDxh = scalar("dAlphaDxh")
nutF = np.loadtxt(SB + "/wstate_baseline_nutFrozen.mtx"); nuEff = NU + nutF

# boundary values for reconstruction (mirrors b31b_m2.py bvals())
def bvals(q, scalarf=True, pdefault=None):
    _, _, bd = read_field(Q1 + "/" + q)
    out = np.zeros(nB if scalarf else (nB, 3))
    for pn, sz in zip(patch_names, patch_sizes):
        o0 = bface_start[pn] - nIF
        bt, bv = bd.get(pn, (None, None))
        if bv is None:
            if scalarf: out[o0:o0+sz] = pdefault[o0:o0+sz] if pdefault is not None else 0.0
            else: out[o0:o0+sz] = (pdefault[o0:o0+sz] if pdefault is not None else 0.0)
        else:
            out[o0:o0+sz] = bv if scalarf else bv
    return out
U_bvals = bvals("U", scalarf=False)
p_bvals = bvals("p", pdefault=p[bc_cell])

# round-tagged gated exports
R = {}
for tag in ("_r1", "_r2"):
    R[tag] = dict(
        mom_gDP = scalar("gsenshPressureDropMomentum" + tag),
        prow_gDP = scalar("gsenshPressureDropPressureRow" + tag),
        raw_gDP = scalar("gsenshPressureDropRaw" + tag),
        rxT = scalar("rxPressureRowT" + tag),
        rxTb = scalar("rxPressureRowTb" + tag),
        rxFx = scalar("rxFluxDirectT" + tag),
        gmtMom = scalar("gsenshMeanTMomentum" + tag),
        gmtPRow = scalar("gsenshMeanTPressureRow" + tag),
        gmtFx = scalar("gsenshMeanTFluxDirect" + tag),
        gmtC = scalar("gsenshMeanTThermalC" + tag),
        Uc = vector("Uc" + tag), Ub = vector("Ub" + tag),
        pc = scalar("pc" + tag), pb = scalar("pb" + tag),
        Tb = scalar("Tb" + tag),
        rxrAU = scalar("rxrAU" + tag),
        rxDHbyA = vector("rxDHbyA" + tag),
    )
    g0v, _, _ = read_field(Q1 + "/rxG0Field" + tag)
    R[tag]["rxG0"] = g0v
log("round-tagged fields loaded (r1+r2)")

# untagged /1/ fields for cross-round attribution evidence
untag = dict(
    Uc = vector("Uc"), Ub = vector("Ub"),
    pc = scalar("pc"), pb = scalar("pb"), Tb = scalar("Tb"),
    gsenshPD = scalar("gsenshPressureDrop"),    # post-chain r2 (mutated by filter)
    fsenshMT = scalar("fsenshMeanT"),           # post-chain r2
    fsensMT = scalar("fsensMeanT"),             # post-chain solved
    gsensPD = scalar("gsensPressureDrop"),      # post-chain solved
    dfdx = scalar("dfdx"),
    xp = scalar("xp"), designMask = scalar("designMask"),
)
log("untagged /1/ fields loaded")

# ---------------- G0: cross-round attribution evidence -----------------------
g0 = {}
for q, tag in [("Uc", "_r1"), ("Ub", "_r1"), ("pc", "_r1"), ("pb", "_r1"), ("Tb", "_r1")]:
    rel = float(np.linalg.norm(untag[q] - R[tag][q])/max(np.linalg.norm(R[tag][q]), 1e-300))
    g0[q + "_untag_eq_r1"] = rel
    log("G0: /1/%s == %s%s relL2=%.3e (0 => untagged is round-1)" % (q, q, tag, rel))
log("G0 attribution: %s" % g0)

# ---------------- G1: same-round closure -------------------------------------
g1 = {}
for tag in ("_r1", "_r2"):
    mom = R[tag]["mom_gDP"]; prow = R[tag]["prow_gDP"]; raw = R[tag]["raw_gDP"]
    closure = mom + prow
    denom = np.linalg.norm(raw)
    rel = float(np.linalg.norm(raw - closure)/max(denom, 1e-300))
    mx = float(np.max(np.abs(raw - closure)))
    g1[tag] = dict(relL2=rel, maxAbs=mx, norm_raw=denom)
    log("G1 %s: gsenshPD_raw == mom+prow relL2=%.3e maxAbs=%.3e (target <=1e-12)"
        % (tag, rel, mx))
    # same-round Uc/pc independent recompute (offline ingredients)
    mom_off = -dAlphaDxh*np.einsum("ci,ci->c", U, R[tag]["Uc"])*V
    prow_off = -R[tag]["rxT"]*dAlphaDxh
    rel_m = float(np.linalg.norm(mom_off - mom)/max(np.linalg.norm(mom), 1e-300))
    rel_p = float(np.linalg.norm(prow_off - prow)/max(np.linalg.norm(prow), 1e-300))
    g1[tag]["mom_offline_relL2"] = rel_m; g1[tag]["prow_offline_relL2"] = rel_p
    log("G1 %s: mom(offline Uc) relL2=%.3e | prow(offline rxT*dA) relL2=%.3e"
        % (tag, rel_m, rel_p))
    # J label: 4-segment sum
    js = R[tag]["gmtMom"] + R[tag]["gmtPRow"] + R[tag]["gmtFx"] + R[tag]["gmtC"]
    g1[tag]["J_seg_sum_L2"] = float(np.linalg.norm(js))
    log("G1 %s: J 4-seg sum |L2|=%.4e" % (tag, np.linalg.norm(js)))

# ---------------- dHbyA/drAU/HbyA (coordinator reconstruction) ---------------
qp = (phi_int >= 0).astype(float); qn = 1 - qp
divp = np.zeros(N)
np.add.at(divp, owner, phi_int); np.add.at(divp, nei, -phi_int)
np.add.at(divp, bc_cell, phi_bvals)
gam = w*nuEff[owner] + (1-w)*nuEff[nei]
kfv = gam*magSfI*dc
lower_m = -qp*phi_int - kfv
upper_m = qn*phi_int - kfv
diag_m = np.zeros(N)
np.add.at(diag_m, owner, -lower_m); np.add.at(diag_m, nei, -upper_m)
diag_m -= divp; diag_m += alpha*V
sumOff_m = np.zeros(N)
np.add.at(sumOff_m, owner, np.abs(upper_m)); np.add.at(sumOff_m, nei, np.abs(lower_m))
ab_m = nuEff[bc_cell]*bmag*bdel
ic_m = np.where(Ufix, -ab_m, phi_bvals)
Db_m = diag_m.copy(); np.add.at(Db_m, bc_cell, np.abs(ic_m))
D_rel = np.maximum(np.abs(Db_m), sumOff_m)/ALPHAREL
mob = V/D_rel
mob_ex = np.loadtxt(SB + "/wstate_baseline_primalPressureMobility.mtx")
g_mob = float(np.linalg.norm(mob - mob_ex)/np.linalg.norm(mob_ex))
log("G-mob: relL2=%.3e (b29 gate)" % g_mob)
assert g_mob < 1e-12, "mobility gate failed"

offU = np.zeros((N, 3))
np.add.at(offU, owner, upper_m[:, None]*U[nei])
np.add.at(offU, nei, lower_m[:, None]*U[owner])
pf_int = w*p[owner] + (1-w)*p[nei]
gradp = np.zeros((N, 3))
tmp = SfI*pf_int[:, None]
np.add.at(gradp, owner, tmp); np.add.at(gradp, nei, -tmp)
np.add.at(gradp, bc_cell, bSf*p_bvals[:, None])
gradp /= V[:, None]
src_exp = -gradp
Uf_int = w[:, None]*U[owner] + (1-w)[:, None]*U[nei]
gU = np.zeros((N, 3, 3))
tmpU = SfI[:, None]*Uf_int[:, :, None]
np.add.at(gU, owner, tmpU); np.add.at(gU, nei, -tmpU)
np.add.at(gU, bc_cell, bSf[:, None]*U_bvals[:, :, None])
gU /= V[:, None, None]
dev2T = np.swapaxes(gU, 1, 2).copy()
tr = np.trace(gU, axis1=1, axis2=2)
dev2T[:, 0, 0] -= 2.0/3.0*tr; dev2T[:, 1, 1] -= 2.0/3.0*tr; dev2T[:, 2, 2] -= 2.0/3.0*tr
nuint = 0.5*(nuEff[owner] + nuEff[nei])
Tf_int = w[:, None, None]*(nuint[:, None, None]*dev2T[owner]) \
       + (1-w)[:, None, None]*(nuint[:, None, None]*dev2T[nei])
Tf_bnd = nuEff[bc_cell][:, None, None]*dev2T[bc_cell]
divdev = np.zeros((N, 3))
tmpD = np.einsum("fi,fij->fj", SfI, Tf_int)
np.add.at(divdev, owner, tmpD); np.add.at(divdev, nei, -tmpD)
np.add.at(divdev, bc_cell, np.einsum("fi,fij->fj", bSf, Tf_bnd))
divdev /= V[:, None]
src_exp += -divdev
ob = phi_bvals >= 0.0
src_bcomp = np.zeros((N, 3))
ib = ~ob
np.add.at(src_bcomp, bc_cell[ib], -phi_bvals[ib, None]*U_bvals[ib])
np.add.at(src_bcomp, bc_cell[Ufix], ab_m[Ufix, None]*U_bvals[Ufix])
src_vec = src_exp*V[:, None] + src_bcomp
src_vec += (1.0-ALPHAREL)*D_rel[:, None]*U
H = (src_vec - offU)/V[:, None]
HbyA = mob[:, None]*H
drAU = -mob*mob/ALPHAREL
dHbyA = (mob/ALPHAREL)[:, None]*((1.0-ALPHAREL)*U - HbyA)
g0_int = magSfI*dc*(p[nei] - p[owner])
log("coordinator dHbyA/drAU built (|dHbyA|L2=%.4e |drAU|L2=%.4e)"
    % (np.linalg.norm(dHbyA), np.linalg.norm(drAU)))

SfDhA_o = np.einsum("ij,ij->i", SfI, dHbyA[owner])
SfDhA_n = np.einsum("ij,ij->i", SfI, dHbyA[nei])
def T_row(lam):
    """production T = T1+T2 for adjoint pressure field lam (per cell)."""
    pcDiff = lam[owner] - lam[nei]
    T1 = np.zeros(N); T2 = np.zeros(N)
    np.add.at(T1, owner,  w*SfDhA_o*pcDiff)
    np.add.at(T1, nei,  (1-w)*SfDhA_n*pcDiff)
    np.add.at(T2, owner, -w*g0_int*pcDiff*drAU[owner])
    np.add.at(T2, nei,  -(1-w)*g0_int*pcDiff*drAU[nei])
    for bf in am_assign:
        c = bc_cell[bf]
        T1[c] += float(np.dot(bSf[bf], dHbyA[c]))*lam[c]
    return T1 + T2

# ---------------- G2: T-ring mirror ------------------------------------------
g2 = {}
for tag in ("_r1", "_r2"):
    for lbl, prodT, lam in [("pc", "rxT", R[tag]["pc"]),
                            ("pb", "rxTb", R[tag]["pb"])]:
        T_off = T_row(lam)
        T_prod = R[tag][prodT]
        rel = float(np.linalg.norm(T_prod - T_off)/max(np.linalg.norm(T_prod), 1e-300))
        corr = float(np.corrcoef(T_prod, T_off)[0, 1])
        g2[tag + "_" + lbl] = dict(relL2=rel, corr=corr,
                                   norm_prod=float(np.linalg.norm(T_prod)),
                                   norm_off=float(np.linalg.norm(T_off)))
        log("G2 %s %s: production rxPressureRowT vs offline T_row relL2=%.4e "
            "corr=%.6f (target relL2<2%%)" % (tag, lbl, rel, corr))
    # prow rebuilt from production T vs gated prow export
    prow_re = -R[tag]["rxT"]*dAlphaDxh
    rel = float(np.linalg.norm(R[tag]["prow_gDP"] - prow_re)
                / max(np.linalg.norm(R[tag]["prow_gDP"]), 1e-300))
    g2[tag + "_prow_rebuilt"] = rel
    log("G2 %s: gated prow == -rxT*dAlphaDxh relL2=%.3e" % (tag, rel))

# ---------------- G4: HbyA calibration ---------------------------------------
g4 = {}
for tag in ("_r1", "_r2"):
    rel_dh = float(np.linalg.norm(R[tag]["rxDHbyA"] - dHbyA)
                   / max(np.linalg.norm(dHbyA), 1e-300))
    corr_dh = float(np.corrcoef(R[tag]["rxDHbyA"].ravel(), dHbyA.ravel())[0, 1])
    # rxrAU = 1/A(): compare against coordinator mobility mob (same diag scale)
    rel_r = float(np.linalg.norm(R[tag]["rxrAU"] - mob)
                  / max(np.linalg.norm(mob), 1e-300))
    corr_r = float(np.corrcoef(R[tag]["rxrAU"], mob)[0, 1])
    rel_g0 = float(np.linalg.norm(R[tag]["rxG0"] - g0_int)
                   / max(np.linalg.norm(g0_int), 1e-300))
    g4[tag] = dict(rxDHbyA_vs_dHbyA_relL2=rel_dh, rxDHbyA_vs_dHbyA_corr=corr_dh,
                   rxrAU_vs_mob_relL2=rel_r, rxrAU_vs_mob_corr=corr_r,
                   rxG0_vs_g0_relL2=rel_g0)
    log("G4 %s: rxDHbyA vs dHbyA relL2=%.4e corr=%.6f | rxrAU vs mob relL2=%.4e "
        "corr=%.6f | rxG0Field vs g0 relL2=%.4e"
        % (tag, rel_dh, corr_dh, rel_r, corr_r, rel_g0))

# ---------------- z (solver-exported chain^T directions) ---------------------
zx = np.loadtxt(B8 + "/stageB6_rxc_z_analytic.mtx").reshape(3, N)
Z = {d: zx[k] for k, d in enumerate(DIRS)}
Dall = np.loadtxt(B8 + "/stageB6_dirs.mtx").reshape(3, N)
DIRS_F = {d: Dall[k] for k, d in enumerate(DIRS)}

# ---------------- chain_field (filter_chainrule.H transpose) -----------------
xp = untag["xp"]; dm = untag["designMask"]
ETA5 = 0.7480051615837531
DEL = 8.0
drho = np.zeros(N); dPE = np.zeros(N)
lo = xp <= ETA5
pe_lo = np.exp(-DEL*(1.0 - xp[lo]/ETA5))
drho[lo] = DEL*pe_lo + np.exp(-DEL)
dPE[lo] = pe_lo*(1.0 - DEL*xp[lo]/ETA5) - np.exp(-DEL)
hi = ~lo
pe_hi = np.exp(-DEL*(xp[hi]-ETA5)/(1.0-ETA5))
drho[hi] = DEL*pe_hi + np.exp(-DEL)
dPE[hi] = pe_hi*(1.0 - DEL*(1.0-xp[hi])/(1.0-ETA5)) - np.exp(-DEL)
# production b coefficient is a constant 1.0 (createFields.H dimensionedScalar),
# so the discrete Helmholtz is (Lm - diag(V)) phi = -V*g  (b*V = V).
dfm_v, _, _ = read_field(Q1 + "/designFilterFaceMask")
dfm_int = dfm_v[:nIF]
kf_f = dfm_int*magSfI*dc
import scipy.sparse as sp
import scipy.sparse.linalg as spla
Lrows = np.concatenate([owner, nei, owner, nei])
Lcols = np.concatenate([owner, nei, nei, owner])
Lvals = np.concatenate([kf_f, kf_f, -kf_f, -kf_f])
Lm = sp.coo_matrix((Lvals, (Lrows, Lcols)), shape=(N, N)).tocsc()
Amat = Lm - sp.diags(V)
luF = spla.splu(Amat)
def chain_field(f):
    c = float((dm*dPE*f).sum())/float((dm*dPE*V).sum())
    g = drho*f + dm*V*(1.0-drho)*c
    return dm*luF.solve(-V*g)   # dm**2 == dm (production masks output twice)
log("chain_field ready")

# ---------------- G3: same-round projections vs B31 cross-round --------------
rows = []
g3 = {}
for d in DIRS:
    z = Z[d]; Dv = DIRS_F[d]
    for tag in ("_r1", "_r2"):
        mom = R[tag]["mom_gDP"]; prow = R[tag]["prow_gDP"]
        raw = R[tag]["raw_gDP"]
        mom_z = float(mom @ z); prow_z = float(prow @ z)
        closure = mom_z + prow_z
        raw_z = float(raw @ z)
        ch_total = float(chain_field(raw) @ Dv)
        ch_mom = float(chain_field(mom) @ Dv)
        ch_prow = float(chain_field(prow) @ Dv)
        # untagged post-chain production for the same-round ADJ anchor
        r = dict(dir=d, tag=tag,
                 mom_z=mom_z, prow_z=prow_z, closure_z=closure, raw_z=raw_z,
                 ch_mom=ch_mom, ch_prow=ch_prow, ch_closure=ch_mom+ch_prow,
                 ch_total=ch_total,
                 ch_valid=False,   # offline chain solve = near-singular artifact, see gc.isolation
                 ADJ_gDP=ADJ_gDP[d], ADJ_J=ADJ_J[d],
                 true_flow=bPD_w[d],
                 residual_closure_vs_raw=abs(closure-raw_z)/max(abs(raw_z),1e-300),
                 total_lesion=ch_total - bPD_w[d],
                 flow_lesion=(ch_mom+ch_prow) - bPD_w[d])
        rows.append(r)
        log("%s %s gDP: mom_z=%+.6f prow_z=%+.6f closure=%+.6f raw_z=%+.6f | "
            "chain mom=%+.6f prow=%+.6f closure=%+.6f total=%+.6f | "
            "ADJ=%+.6f true=%+.6f | flow_lesion=%+.6f total_lesion=%+.6f"
            % (d, tag, mom_z, prow_z, closure, raw_z,
               ch_mom, ch_prow, ch_mom+ch_prow, ch_total,
               ADJ_gDP[d], bPD_w[d], r["flow_lesion"], r["total_lesion"]))
    # B31 §7 cross-round comparison (D1 only, others None)
    if B31_XR[d][0] is not None:
        xr = B31_XR[d]
        log("B31-XR %s: mom=%+.6f prow_prod=%+.6f closure=%+.6f prod=%+.6f | "
            "true_prow=%+.6f prow_lesion=%+.6f total_lesion=%+.6f"
            % (d, *xr))

# ---------------- G-c: chain sanity (pre-chain -> post-chain disk) -----------
ch_disk_gDP = chain_field(R["_r2"]["raw_gDP"])
gc2 = float(np.linalg.norm(ch_disk_gDP - untag["gsensPD"])
            / max(np.linalg.norm(untag["gsensPD"]), 1e-300))
gc2_corr = float(np.corrcoef(ch_disk_gDP, untag["gsensPD"])[0, 1])
log("G-c2: chain(raw_gDP_r2) vs disk gsensPressureDrop relL2=%.3e corr=%.4f (INVALID offline solve)"
    % (gc2, gc2_corr))
js_r2 = R["_r2"]["gmtMom"] + R["_r2"]["gmtPRow"] + R["_r2"]["gmtFx"] + R["_r2"]["gmtC"]
ch_disk_J = chain_field(js_r2)
# compare against the unscaled post-chain solved field fsensMeanT (dfdx on disk
# is objectiveGradientScale-scaled and zero outside active design cells)
gc1 = float(np.linalg.norm(ch_disk_J - untag["fsensMT"])
            / max(np.linalg.norm(untag["fsensMT"]), 1e-300))
log("G-c1: chain(J-seg-sum_r2) vs disk fsensMeanT relL2=%.3e (INVALID offline solve)" % gc1)
gc1_dfdx = float(np.linalg.norm(ch_disk_J - untag["dfdx"])
                 / max(np.linalg.norm(untag["dfdx"]), 1e-300))
log("G-c1b: chain(J-seg-sum_r2) vs disk dfdx (scaled) relL2=%.3e (info)" % gc1_dfdx)
# ---- G-c isolation evidence (this session's /tmp investigations) ----
gc_iso = dict(
    chain_valid=False,
    raw_z_valid=True,               # same-round raw_z is a direct gate export, bit-exact
    source_bit_exact_relL2=2.911e-10,   # which_round.py: chain_src(raw_gDP_r2) vs disk gsenshPD
    forward_filter_corr=0.999870,       # scale_test.py: solve(Lm+VI, Vx) vs production xp
    forward_filter_norm_exact=1.4199e+02,  # |X1| == |xp| to full precision
    disk_vs_source_corr=0.9796,         # adj_rounds.py: disk gsensPD vs its source gsenshPD
    disk_null_mode_Lm_max_design=4.678e-06,  # source_id2.py: max |(Lm@disk)| on design cells
    null_mode_scale_Vg_over_V=6.3e-1,   # |V*g|/V amplification scale, matches |disk|=5.17e-1
    disk_norm=5.17e-1,
    prod_dicpcg_final=2.07e-10,         # Log.b32_decomp.txt L10420-10423 (r2, 10 iter)
    note=(
        "OFFLINE chain_field uses solve(Lm+VI) with Lm = production graph Laplacian "
        "(gaussLaplacianScheme verified; forward filter reproduced bit-exact corr 0.999870). "
        "The solve is near-singular: zeroGradient BCs everywhere -> Lm has a constant null "
        "mode, V~1.25e-10 on diagonal vs Lm scale ~1e-4 -> null-mode amplification 1/V ~ 8e9. "
        "Production DICPCG (tol 1e-9, final 2.07e-10) and offline splu diverge in the null "
        "direction, so ch_total/ch_mom/ch_prow projections are INVALID (disk is a near-null-"
        "mode field: corr 0.98 with its source, (Lm@disk)~0). G3 verdict must anchor on "
        "raw_z/closure_z (direct gated exports, bit-exact), NOT on chain projections."
    ),
)

# ---------------- dump -------------------------------------------------------
res = dict(g0=g0, g1=g1, g2=g2, g4=g4, g3=rows, gc=dict(c1=gc1, c1_dfdx=gc1_dfdx,
        c2=gc2, c2_corr=gc2_corr, mob=g_mob, isolation=gc_iso))
json.dump(res, open(OUT + "/b32_decomp.json", "w"), indent=1, default=str)
with open(OUT + "/b32_decomp.tsv", "w") as f:
    keys = list(rows[0].keys())
    f.write("\t".join(keys) + "\n")
    for r in rows:
        f.write("\t".join(str(r[k]) for k in keys) + "\n")
log("B32 DECOMP WRITTEN. DONE t=%.1fs" % (time.time()-t0))
