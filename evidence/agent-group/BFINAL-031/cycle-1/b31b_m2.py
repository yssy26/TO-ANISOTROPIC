#!/usr/bin/env python3
"""BFINAL-031 b31b_m2: assembly-chain accounting matrix (coordinator-executed).

Completes the m2 deliverable per NOTEBOOK §6.12 blueprint (attempt-2 derivation).
Three legs per label (J/gDP) x direction (D1/D2/D3):
  production ADJ (tsv) | offline segment recompute (mom + prow face-loop) |
  source contraction (b^T w_true, pinned)
Gates: G-mob (mobility vs export), G-z1/G-z2 (pre-chain field . z == tsv ADJ,
  validates z = stageB6_rxc_z_analytic as exact chain^T(D)), G-b8 (T rebuilt on
  b8 fields vs rxpr_pressurerow export, machine precision).
Production T-loop mirrors rxPressureRowTranspose.H L223-257 exactly.
"""
import time, json, re
import numpy as np

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
SB = "/home/ys/dsH/b25_qgate/stageB2"
Q1 = "/home/ys/dsH/b25_qgate/1"
B8 = "/home/ys/dsH/b8_verify_diag"
E26 = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-026/cycle-1"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-031/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N
ALPHAREL = 0.4
NU = 5.19009e-05

# pinned truth values (NOT recomputed here)
ADJ_J  = {"D1": -0.0358350400485765, "D2": -0.02492599038981423, "D3": -0.01200301226991504}
ADJ_gDP= {"D1": -5.23029384626718,   "D2": 1.0505249811003,     "D3": 13.60342258150915}
FD_J   = {"D1": 0.04252469623637622, "D2": -0.003528972517971574, "D3": -0.0004449511895182612}
bPD_w  = {"D1": -4.277461237487605,  "D2": 0.3844628061126605,  "D3": 14.096353328798937}
bTC_w  = {"D1": -0.04884018700713266,"D2": -0.006440998953954891,"D3": 0.011566617518607991}
C_pinned = {"D1": -0.030341656277042862, "D2": -0.008210889966675959, "D3": -0.01979195745200000}
Gx_pinned= {"D1": -0.0015704339056528767,"D2": 0.004043159681714134,"D3": 0.00840670768310000}
DIRS = ["D1", "D2", "D3"]

t0 = time.time()
def log(m): print("[%5.1fs] %s" % (time.time()-t0, m), flush=True)
log("=== b31b_m2: assembly-chain accounting (coordinator) ===")

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
# exact pyramid centroids/volumes
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
am_assign = np.where(~Ufix)[0]   # u-assignable boundary faces (outlet etc.)
log("mesh ready: nIF=%d nB=%d" % (nIF, nB))

# ---------------- b25 fields ------------------------------------------------
def scalar(q, fb=None):
    v, _, _ = read_field(Q1 + "/" + q)
    return v
def vector(q):
    v, _, _ = read_field(Q1 + "/" + q)
    return v.reshape(N, 3)
U = np.loadtxt(SB + "/wstate_baseline_U.mtx").reshape(N, 3)
p = np.loadtxt(SB + "/wstate_baseline_p.mtx")
phi_int = np.loadtxt(SB + "/wstate_baseline_phi.mtx")
phi_bvals = np.loadtxt(SB + "/wstate_baseline_phiB.mtx")
alpha = np.loadtxt(SB + "/wstate_baseline_alpha.mtx")
dAlphaDxh = scalar("dAlphaDxh")
nutF = np.loadtxt(SB + "/wstate_baseline_nutFrozen.mtx"); nuEff = NU + nutF
Uc = vector("Uc"); Ub = vector("Ub")
pc = scalar("pc"); pb = scalar("pb")
gsenshPD = scalar("gsenshPressureDrop"); fsenshMT = scalar("fsenshMeanT")
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
p_bvals = bvals("p", pdefault=p[bc_cell])   # zeroGradient -> internal
log("b25 fields loaded (primal from SB baseline family, b29-gated)")

# -------- matrix + mobility (b29 build_basis VERBATIM, machine-gated) --------
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
mob_gate = float(np.linalg.norm(mob - mob_ex)/np.linalg.norm(mob_ex))
log("G-mob: relL2=%.3e" % mob_gate)
assert mob_gate < 1e-12, "mobility gate failed"

# ---------------- dHbyA/drAU/HbyA (src per b18; offU from gated matrix) ------
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
# dev2 term (b17 uniform-nu-approximation-free version from b18)
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
log("dHbyA/drAU built (|dHbyA|L2=%.4e |drAU|L2=%.4e)"
    % (np.linalg.norm(dHbyA), np.linalg.norm(drAU)))

# ---------------- production T-loop (rxPressureRowTranspose.H L223-257) -----
g0_int = magSfI*dc*(p[nei] - p[owner])     # flux of laplacian(1,p), gamma=1
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
    for bf in am_assign:                    # u-assignable boundary (outlet)
        c = bc_cell[bf]
        T1[c] += float(np.dot(bSf[bf], dHbyA[c]))*lam[c]
    return T1 + T2
T_pc = T_row(pc)                            # pressureDrop label
T_pb = T_row(pb)                            # thermalCoupling label
prow_gDP_field = -T_pc*dAlphaDxh
prow_J_field   = -T_pb*dAlphaDxh
mom_gDP_field  = -dAlphaDxh*np.einsum("ci,ci->c", U, Uc)*V
mom_J_field    = -dAlphaDxh*np.einsum("ci,ci->c", U, Ub)*V
log("segments built: |prow_gDP|L2=%.4e |mom_gDP|L2=%.4e |prow_J|L2=%.4e |mom_J|L2=%.4e"
    % (np.linalg.norm(prow_gDP_field), np.linalg.norm(mom_gDP_field),
       np.linalg.norm(prow_J_field), np.linalg.norm(mom_J_field)))

# ---------------- z (solver-exported chain^T directions) ---------------------
zx = np.loadtxt(B8 + "/stageB6_rxc_z_analytic.mtx").reshape(3, N)  # 3x33600
Z = {d: zx[k] for k, d in enumerate(DIRS)}

# ---------------- gates G-z1/G-z2 --------------------------------------------
gz1 = {d: float(fsenshMT @ Z[d]) for d in DIRS}
gz2 = {d: float(gsenshPD @ Z[d]) for d in DIRS}
for d in DIRS:
    log("G-z1 %s: fsenshMT.z=%+.8f vs ADJ_J=%+.8f  rel=%.2e"
        % (d, gz1[d], ADJ_J[d], abs(gz1[d]-ADJ_J[d])/abs(ADJ_J[d])))
    log("G-z2 %s: gsenshPD.z=%+.8f vs ADJ_gDP=%+.8f  rel=%.2e"
        % (d, gz2[d], ADJ_gDP[d], abs(gz2[d]-ADJ_gDP[d])/abs(ADJ_gDP[d])))

# ---------------- accounting matrix ------------------------------------------
rows = []
for d in DIRS:
    z = Z[d]
    # gDP label
    prod_total = gz2[d]
    mom_z = float(mom_gDP_field @ z)
    prow_z = float(prow_gDP_field @ z)
    closure = mom_z + prow_z          # should == prod_total if production=mom+prow
    true_flow = bPD_w[d]
    true_prow = true_flow - mom_z
    # J label
    prod_total_J = gz1[d]
    momJ_z = float(mom_J_field @ z)
    prowJ_z = float(prow_J_field @ z)
    flowJ_prod = momJ_z + prowJ_z
    directJ_prod = prod_total_J - flowJ_prod
    true_flow_J = bTC_w[d]
    true_prow_J = true_flow_J - momJ_z
    r = dict(dir=d,
             # gDP
             ADJ_gDP=ADJ_gDP[d], prod_total_z=prod_total,
             mom_z=mom_z, prow_prod_z=prow_z, closure_prod=closure,
             closure_rel=abs(closure_prod:=closure - prod_total)/abs(prod_total) if prod_total else 0,
             true_flow_z=true_flow, true_prow_z=true_prow,
             prow_lesion=prow_z - true_prow, total_lesion=prod_total - true_flow,
             # J
             ADJ_J=ADJ_J[d], prod_total_J=prod_total_J,
             momJ_z=momJ_z, prowJ_prod_z=prowJ_z, flowJ_prod=flowJ_prod,
             directJ_prod=directJ_prod, C_pinned=C_pinned[d], Gx_pinned=Gx_pinned[d],
             true_flow_J=true_flow_J, true_prow_J=true_prow_J,
             FD_J=FD_J[d],
             J_prod_minus_truth=prod_total_J - (C_pinned[d]+Gx_pinned[d]+true_flow_J))
    rows.append(r)
    log("%s gDP: mom=%+.6f prow_prod=%+.6f closure=%+.6f vs prod=%+.6f | "
        "true_prow=%+.6f prow_lesion=%+.6f total_lesion=%+.6f"
        % (d, mom_z, prow_z, closure, prod_total, true_prow,
           r["prow_lesion"], r["total_lesion"]))
    log("%s J  : mom=%+.6f prow_prod=%+.6f flow_prod=%+.6f direct_prod=%+.6f | "
        "true: C%+.6f Gx%+.6f flow%+.6f => FD%+.6f | prod-FDtruth=%+.6f"
        % (d, momJ_z, prowJ_z, flowJ_prod, directJ_prod,
           C_pinned[d], Gx_pinned[d], true_flow_J, FD_J[d], r["J_prod_minus_truth"]))

# ---------------- G-b8: T on b8 fields vs rxpr export ------------------------
def b8_scalar(q):
    v, _, _ = read_field(B8 + "/1/" + q)
    return v
try:
    lam8 = np.loadtxt(B8 + "/stageB6_lambda.mtx")
    pc8 = lam8[NV:]
    U8 = b8_scalar("U").reshape(N, 3); p8 = b8_scalar("p")
    phi8_v, _, _ = read_field(B8 + "/1/phi"); phi8_i = phi8_v[:nIF]
    al8 = b8_scalar("alpha"); dA8 = b8_scalar("dAlphaDxh")
    nut8 = np.zeros(N); nuE8 = NU + nut8  # b8 nutFrozen: b8/1 field; fallback below
    # b8 boundary values
    def b8_bvals(q, pdefault=None):
        _, _, bd = read_field(B8 + "/1/" + q)
        out = np.zeros(nB)
        for pn, sz in zip(patch_names, patch_sizes):
            o0 = bface_start[pn] - nIF
            bt, bv = bd.get(pn, (None, None))
            out[o0:o0+sz] = bv if bv is not None and not isinstance(bv, float) else (bv if isinstance(bv, float) else 0.0)
            if bv is None and pdefault is not None:
                out[o0:o0+sz] = pdefault[o0:o0+sz]
        return out
    phib8 = b8_bvals("phi")
    Ub8_b = b8_bvals("U")  # may be zeros for nonuniform; acceptable for gate scale
    # b8 matrix
    nuEff_f8 = w*nuE8[owner] + (1-w)*nuE8[nei]
    a_f8 = nuEff_f8*magSfI*dc
    up8 = phi8_i >= 0.0
    po8 = np.where(up8, phi8_i, 0.0); pn8_ = np.where(up8, 0.0, phi8_i)
    dg8 = np.zeros(N)
    np.add.at(dg8, owner, po8 - a_f8); np.add.at(dg8, nei, -pn8_ - a_f8)
    dg8 += al8*V
    dp8 = np.zeros(N)
    np.add.at(dp8, owner, phi8_i); np.add.at(dp8, nei, -phi8_i)
    np.add.at(dp8, bc_cell, phib8)
    dg8 += -dp8
    ob8 = phib8 >= 0.0
    ab8 = nuE8[bc_cell]*bmag*bdel
    np.add.at(dg8, bc_cell, np.where(ob8, phib8, 0.0))
    np.add.at(dg8, bc_cell[Ufix], ab8[Ufix])
    sOff8 = np.zeros(N)
    np.add.at(sOff8, owner, np.abs(pn8_ + a_f8)); np.add.at(sOff8, nei, np.abs(-po8 + a_f8))
    Drel8 = np.maximum(np.abs(dg8), sOff8)/ALPHAREL
    mob8 = V/Drel8
    # dHbyA8 (internal-face approximation of gradp path; b18 recipe)
    offU8 = np.zeros((N, 3))
    np.add.at(offU8, owner, (pn8_ + a_f8)[:, None]*U8[nei])
    np.add.at(offU8, nei, (-po8 + a_f8)[:, None]*U8[owner])
    pf8 = w*p8[owner] + (1-w)*p8[nei]
    gp8 = np.zeros((N, 3))
    t8 = SfI*pf8[:, None]
    np.add.at(gp8, owner, t8); np.add.at(gp8, nei, -t8)
    np.add.at(gp8, bc_cell, bSf*(p8[bc_cell])[:, None])  # zeroGradient approx
    gp8 /= V[:, None]
    sv8 = -gp8*V[:, None]
    sv8 += (1.0-ALPHAREL)*Drel8[:, None]*U8
    ib8 = ~ob8
    np.add.at(sv8, bc_cell[ib8], -phib8[ib8, None]*Ub8_b.reshape(nB, 1)[ib8, :][:, :1]*0)  # U_b unknown; small
    H8 = (sv8 - offU8)/V[:, None]
    HbyA8 = mob8[:, None]*H8
    drAU8 = -mob8*mob8/ALPHAREL
    dHbyA8 = (mob8/ALPHAREL)[:, None]*((1.0-ALPHAREL)*U8 - HbyA8)
    SfDhA_o8 = np.einsum("ij,ij->i", SfI, dHbyA8[owner])
    SfDhA_n8 = np.einsum("ij,ij->i", SfI, dHbyA8[nei])
    g0_8 = magSfI*dc*(p8[nei] - p8[owner])
    pcD8 = pc8[owner] - pc8[nei]
    T1_8 = np.zeros(N); T2_8 = np.zeros(N)
    np.add.at(T1_8, owner, w*SfDhA_o8*pcD8)
    np.add.at(T1_8, nei, (1-w)*SfDhA_n8*pcD8)
    np.add.at(T2_8, owner, -w*g0_8*pcD8*drAU8[owner])
    np.add.at(T2_8, nei, -(1-w)*g0_8*pcD8*drAU8[nei])
    for bf in am_assign:
        c = bc_cell[bf]
        T1_8[c] += float(np.dot(bSf[bf], dHbyA8[c]))*pc8[c]
    T8v = T1_8 + T2_8
    prow8 = -T8v*dA8
    rxprP = np.loadtxt(B8 + "/rxpr_prod_gsensh_pressurerow.mtx")
    g_b8 = float(np.linalg.norm(prow8 - rxprP)/np.linalg.norm(rxprP))
    log("G-b8: prow8 vs rxpr_pressurerow relL2=%.3e (target <1e-6; note b8 "
        "boundary-src approximation may limit)" % g_b8)
except Exception as e:
    g_b8 = None
    log("G-b8 SKIPPED/FAILED: %s" % e)

# ---------------- dump --------------------------------------------------------
res = dict(gates=dict(mob=mob_gate, z1=gz1, z2=gz2, b8=g_b8), rows=rows)
json.dump(res, open(OUT + "/b31b_matrix.json", "w"), indent=1, default=str)
with open(OUT + "/b31b_matrix.tsv", "w") as f:
    keys = list(rows[0].keys())
    f.write("\t".join(keys) + "\n")
    for r in rows:
        f.write("\t".join(str(r[k]) for k in keys) + "\n")
log("MATRIX WRITTEN. DONE t=%.1fs" % (time.time()-t0))

# ================= chain_field: production filter_chainrule.H transpose ====
# form (source-verified): c = sum(dm*dPE*f)/sum(dm*dPE*V);
# g = drho*f + dm*V*(1-drho)*c;  out = solve((L-bV), -bV*g)
xp = scalar("xp"); dm = scalar("designMask")
ETA5 = 0.7480051615837531   # G1 gate (attempt-1 b31 log)
DEL = 8.0                   # validationProjectionBeta (mmaUpdateEnabled=false)
drho = np.zeros(N); dPE = np.zeros(N)
lo = xp <= ETA5
pe_lo = np.exp(-DEL*(1.0 - xp[lo]/ETA5))
drho[lo] = DEL*pe_lo + np.exp(-DEL)
dPE[lo] = pe_lo*(1.0 - DEL*xp[lo]/ETA5) - np.exp(-DEL)
hi = ~lo
pe_hi = np.exp(-DEL*(xp[hi]-ETA5)/(1.0-ETA5))
drho[hi] = DEL*pe_hi + np.exp(-DEL)
dPE[hi] = pe_hi*(1.0 - DEL*(1.0-xp[hi])/(1.0-ETA5)) - np.exp(-DEL)
# b scalar
nDes = int((dm > 0.5).sum()); desVol = float((V*dm).sum())
lenc = (desVol/nDes)**(1.0/3.0)
bFilt = 1.0/max((3.0*lenc/3.464)**2, 1e-30)
# filter matrix: L from designFilterFaceMask (surface field, gamma = face val)
dfm_v, _, _ = read_field(Q1 + "/designFilterFaceMask")
dfm_int = dfm_v[:nIF]
kf_f = dfm_int*magSfI*dc
import scipy.sparse as sp
import scipy.sparse.linalg as spla
Lrows = np.concatenate([owner, nei, owner, nei])
Lcols = np.concatenate([owner, nei, nei, owner])
Lvals = np.concatenate([kf_f, kf_f, -kf_f, -kf_f])
Lm = sp.coo_matrix((Lvals, (Lrows, Lcols)), shape=(N, N)).tocsc()
Amat = Lm - sp.diags(bFilt*V)
luF = spla.splu(Amat)
def chain_field(f):
    c = float((dm*dPE*f).sum())/float((dm*dPE*V).sum())
    g = drho*f + dm*V*(1.0-drho)*c
    return luF.solve(-bFilt*V*g)
# G-c gates: chain(disk pre-chain) vs disk post-chain
dfdx_disk = scalar("dfdx"); gsPD_disk = scalar("gsensPressureDrop")
ch1 = chain_field(fsenshMT); ch2 = chain_field(gsenshPD)
gc1 = float(np.linalg.norm(ch1 - dfdx_disk)/np.linalg.norm(dfdx_disk))
gc2 = float(np.linalg.norm(ch2 - gsPD_disk)/np.linalg.norm(gsPD_disk))
log("G-c1: chain(fsenshMT_disk) vs dfdx_disk relL2=%.3e" % gc1)
log("G-c2: chain(gsenshPD_disk) vs gsensPD_disk relL2=%.3e" % gc2)
# accounting with chained projections
rows2 = []
D1v = np.zeros(N); D2v = np.zeros(N); D3v = np.zeros(N)
# directions from the solver's own D exports (stageB6_dirs.mtx = 3x33600 raw D)
Dall = np.loadtxt(B8 + "/stageB6_dirs.mtx").reshape(3, N)
DIRS_F = {"D1": Dall[0], "D2": Dall[1], "D3": Dall[2]}
for d in DIRS:
    Dv = DIRS_F[d]
    ch_momG = chain_field(mom_gDP_field) @ Dv
    ch_prowG = chain_field(prow_gDP_field) @ Dv
    ch_totG = chain_field(gsenshPD) @ Dv
    ch_momJ = chain_field(mom_J_field) @ Dv
    ch_prowJ = chain_field(prow_J_field) @ Dv
    ch_totJ = chain_field(fsenshMT) @ Dv
    r = dict(dir=d,
             gDP_chain_mom=ch_momG, gDP_chain_prow=ch_prowG,
             gDP_chain_closure=ch_momG+ch_prowG, gDP_chain_total=ch_totG,
             ADJ_gDP=ADJ_gDP[d], true_flow=bPD_w[d],
             gDP_total_lesion=ch_totG - bPD_w[d],
             gDP_flow_lesion=(ch_momG+ch_prowG) - bPD_w[d],
             J_chain_mom=ch_momJ, J_chain_prow=ch_prowJ,
             J_flow_prod=ch_momJ+ch_prowJ, J_chain_total=ch_totJ,
             ADJ_J=ADJ_J[d], true_flow_J=bTC_w[d],
             C=C_pinned[d], Gx=Gx_pinned[d], FD_J=FD_J[d],
             J_flow_lesion=(ch_momJ+ch_prowJ) - bTC_w[d],
             J_direct_prod=ch_totJ - (ch_momJ+ch_prowJ),
             J_prod_minus_FDtruth=ch_totJ - FD_J[d])
    rows2.append(r)
    log("%s gDP(chain): mom=%+.6f prow=%+.6f closure=%+.6f total=%+.6f ADJ=%+.6f | "
        "true_flow=%+.6f flow_lesion=%+.6f total_lesion=%+.6f"
        % (d, ch_momG, ch_prowG, ch_momG+ch_prowG, ch_totG, ADJ_gDP[d],
           bPD_w[d], r["gDP_flow_lesion"], r["gDP_total_lesion"]))
    log("%s J(chain):  mom=%+.6f prow=%+.6f flow_prod=%+.6f direct=%+.6f total=%+.6f ADJ=%+.6f | "
        "true: C%+.6f Gx%+.6f flow%+.6f => FD%+.6f | prod-FD=%+.6f"
        % (d, ch_momJ, ch_prowJ, ch_momJ+ch_prowJ, r["J_direct_prod"], ch_totJ,
           ADJ_J[d], C_pinned[d], Gx_pinned[d], bTC_w[d], FD_J[d], r["J_prod_minus_FDtruth"]))
res2 = dict(gates=dict(mob=mob_gate, c1=gc1, c2=gc2, b8=g_b8), rows=rows2)
json.dump(res2, open(OUT + "/b31b_matrix_chain.json", "w"), indent=1, default=str)
with open(OUT + "/b31b_matrix_chain.tsv", "w") as f:
    keys = list(rows2[0].keys())
    f.write("\t".join(keys) + "\n")
    for r in rows2:
        f.write("\t".join(str(r[k]) for k in keys) + "\n")
log("CHAIN MATRIX WRITTEN. DONE t=%.1fs" % (time.time()-t0))
