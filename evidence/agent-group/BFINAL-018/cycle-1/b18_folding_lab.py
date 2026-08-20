#!/usr/bin/env python3
"""
BFINAL-018 fix-3 offline derivation lab: reconstruct the thermalCoupling
adjoint rhs b_TC from the frozen b16_states baseline fields, test routing
variants against the B1 gauge identity, and compute the direct-thermal part C.

Definitions (B16/B17 gauge data, frozen validation state):
  FD_J      : true FD of the heat objective (3 directions)
  thermal_g : dJ/dT . (T+ - T-)/(2h)   [TOTAL T response; B1 gauge]
  A         := FD_J - thermal_g         [outlet-flux-only "true flow-mediated"]
  C         : direct thermal xh-part (dDTDxh laplacian folding) contracted
              with the forward projected-design tangent z_k
  identity (production-consistent T-elimination):
      b_TC^T w_true == FD_J - C   (= A + B, B = T-response-to-flux part)
  identity (outlet-only decomposition, matches B1 gauge as implemented):
      b_TC^T w_true == A

Routing variants:
  V0 : current code (transpose of dphi = interp(dU)&Sf + kf(p_n-p_o) with
       INTERNAL kf sign flipped vs the operator)
  V1 : operator-consistent transpose (BFINAL-003 relaxed map):
       dphi_f = Sf&( alphaRel*rAU_u*dH + (1-alphaRel)*interp(dU) )
                - kf_f*(dp_n - dp_o)
       -> hA path (alphaRel*rAU_u) + direct ((1-alphaRel)) + deltaH^T + kf
          with P(own) += kf*g, P(nei) -= kf*g (matches operator J^T kf slot)
  V1k: V1 but with the OLD internal kf sign (isolate the kf-sign effect)
"""
import time, json
import numpy as np

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
E15 = "/home/ys/dsH/b15_export"
S = "/home/ys/dsH/b16_states/1"
ST = "/home/ys/dsH/b16_states/stageB2"
N = 33600; NV = 3*N; NUNK = 4*N
ALPHAREL = 0.4
FD_J = {"D1": 0.0428559144263, "D2": -0.00346892642156, "D3": -0.000282658160322}

t0 = time.time()
def log(m): print(m, flush=True)

# ---------------- field parsing (from b17_rebuildT.py, validated) ----------
def read_field(path):
    txt = open(path).read()
    i = txt.index("internalField")
    j = txt.index("\n(", i)
    k = txt.index("\n)", j)
    if txt[i:j].split()[1] == "uniform":
        seg = txt[i:j]
        u = seg[seg.index("uniform")+len("uniform"):].strip()
        toks = []
        for tok in u.replace(";", " ").replace("(", " ").replace(")", " ").split():
            try:
                toks.append(float(tok))
            except ValueError:
                break
        val = np.array(toks) if len(toks) > 1 else float(toks[0])
        internal = None
    else:
        vals = []
        for line in txt[j+1:k].splitlines():
            for tok in line.replace(";", " ").split():
                for tt in tok.replace("(", " ").replace(")", " ").split():
                    vals.append(float(tt))
        internal = np.array(vals)
        val = None
    bnd = {}
    b = txt.index("boundaryField")
    seg = txt[b:]
    lines = seg.splitlines()
    p = 0
    while p < len(lines):
        s = lines[p].strip()
        if s in ("{", "}", ""):
            p += 1; continue
        if s.endswith("{"):
            name = s[:-1].strip(); hdr = p
        elif (p+1 < len(lines) and lines[p+1].strip() == "{"
              and " " not in s and not s.startswith(("type", "value", "internalField"))):
            name = s; hdr = p+1
        else:
            p += 1; continue
        if name == "boundaryField":
            p += 1; continue
        depth = 1; q = hdr+1; bval = None; btype = None
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
                    u = t[t.index("uniform")+len("uniform"):].rstrip(";").strip()
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
            for line in lines[qq:kk+1]:
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
        faces.append([int(v) for v in s[s.index("(")+1:-1].split()])
def read_labels(path):
    txt = open(path).read(); i = txt.index("\n("); j = txt.index("\n)", i)
    return np.array([int(t) for t in txt[i+2:j].split()], dtype=np.int64)
owner_all = read_labels(POLY + "/owner"); nei = read_labels(POLY + "/neighbour")
owner = owner_all[:len(nei)]; nIF = len(nei)
Sf = np.zeros((nIF, 3))
for fi in range(nIF):
    fv = pts[faces[fi]]
    Sf[fi] = 0.5*np.cross(fv, np.roll(fv, -1, axis=0)).sum(axis=0)
magSf = np.linalg.norm(Sf, axis=1)
cellverts = {}
for fi in range(len(faces)):
    cells = (owner[fi], nei[fi]) if fi < nIF else (owner_all[fi],)
    for c in cells:
        cellverts.setdefault(c, set()).update(faces[fi])
Ccen = np.zeros((N, 3))
for c, vs in cellverts.items(): Ccen[c] = pts[sorted(vs)].mean(axis=0)
d_vec = Ccen[nei] - Ccen[owner]
deltaCoeffs = magSf/np.einsum("ij,ij->i", Sf, d_vec)
Vnueff = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")
V = Vnueff[:, 0]; nuEff = Nnueff = NUEFF = Vnueff[:, 1]
bnd_tbl = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
bc_cell = bnd_tbl[:, 0].astype(np.int64)
bSf = bnd_tbl[:, 1:4]; bmag = bnd_tbl[:, 4]; bdel = bnd_tbl[:, 5]
nBndF = len(bnd_tbl)
weights = np.loadtxt(MESH + "/b16mesh_weights.mtx")
# validate my geometry vs exported Sf/owner
sf_chk = np.loadtxt(MESH + "/b16mesh_sf.mtx")
own_chk = np.loadtxt(MESH + "/b16mesh_owner.mtx")
log("mesh geometry check: Sf maxdiff=%.3e  owner maxdiff=%d"
    % (np.abs(sf_chk - Sf).max(), int(np.abs(own_chk - owner).max())))

patch_names = ["inlet", "outlet", "hotInlet", "hotOutlet",
               "solidEndWalls", "bottomWall", "topWall", "sideWalls"]
nFacesPatch = {"inlet": 84, "outlet": 84, "hotInlet": 196, "hotOutlet": 196,
               "solidEndWalls": 280, "bottomWall": 1120, "topWall": 1120,
               "sideWalls": 4800}
patch_start = {}
o = 0
for pn in patch_names:
    patch_start[pn] = o; o += nFacesPatch[pn]

# ---------------- fields ----------------
U_i, _, U_b = read_field(S + "/U"); U = U_i.reshape(N, 3)
p_i, _, p_b = read_field(S + "/p"); p = p_i
phi_i, _, phi_b = read_field(S + "/phi"); phi_int = phi_i
alpha_i, _, _ = read_field(S + "/alpha"); alpha = alpha_i
Tb_i, _, Tb_b = read_field(S + "/Tb"); Tb = Tb_i
T_i, _, T_b = read_field(S + "/T"); T = T_i
phTh_i, _, phTh_b = read_field(S + "/phiThermal"); phiTh_int = phTh_i
cfm_i, _, _ = read_field(S + "/coldFaceMask"); coldmask_int = cfm_i
dAlphaDxh = read_field(S + "/dAlphaDxh")[0]
dtdxh_i, _, _ = read_field(S + "/dDTDxh")
dDTDxh = dtdxh_i.reshape(N, 6)  # symmTensor storage order xx xy xz yy yz zz

phi_bvals = np.zeros(nBndF); p_bvals = np.zeros(nBndF)
U_fixed = np.zeros(nBndF, dtype=bool); U_assign = np.zeros(nBndF, dtype=bool)
U_bvals = np.zeros((nBndF, 3)); T_bvals = np.zeros(nBndF)
o = 0
for pn in patch_names:
    nFp = nFacesPatch[pn]
    btypeU, bvu = U_b[pn]; btypeP, bvp = p_b[pn]
    btypePhi, bvphi = phi_b[pn]; btypeT, bvt = T_b[pn]
    for arr, src in ((phi_bvals, bvphi), (p_bvals, bvp)):
        if isinstance(src, float): arr[o:o+nFp] = src
        elif src is not None: arr[o:o+nFp] = src
    if bvp is None: p_bvals[o:o+nFp] = p[bc_cell[o:o+nFp]]
    if isinstance(bvu, np.ndarray):
        U_bvals[o:o+nFp] = bvu; U_fixed[o:o+nFp] = True
    else:
        U_bvals[o:o+nFp] = U[bc_cell[o:o+nFp]]
        U_assign[o:o+nFp] = btypeU in ("zeroGradient", "calculated",
                                       "inletOutlet", "outletInlet")
    if isinstance(bvt, float): T_bvals[o:o+nFp] = bvt
    elif bvt is not None: T_bvals[o:o+nFp] = bvt
    else: T_bvals[o:o+nFp] = T[bc_cell[o:o+nFp]]
    o += nFp
log("fields ready (t=%.1fs)" % (time.time()-t0))

# ---------------- g face functional (B0.2-verified formulas) --------------
g = np.zeros(nIF + nBndF)   # internal faces first, then boundary faces
jumpT = T[nei] - T[owner]
down_nei = phiTh_int >= 0.0
adjDown = np.where(down_nei, Tb[nei], Tb[owner])
g[:nIF] = -coldmask_int*adjDown*jumpT
# outlet dJ/dphi
Mflow, Tref = np.loadtxt(ST + "/wstate_outletMeta.mtx")
ophi = phi_bvals[patch_start["outlet"]:patch_start["outlet"]+84]
Toc = T[bc_cell[patch_start["outlet"]:patch_start["outlet"]+84]]
Tmix = float(np.sum(ophi*Toc)/Mflow)
dJdphi_out = -(Toc - Tmix)/(Mflow*Tref)
g[nIF + patch_start["outlet"]:nIF + patch_start["outlet"]+84] += dJdphi_out
log("g: |internal|max=%.6e nz=%d ; |outlet dJ/dphi|max=%.6e Tmix=%.6f M=%.6e"
    % (np.abs(g[:nIF]).max(), int((g[:nIF] != 0).sum()),
       np.abs(dJdphi_out).max(), Tmix, Mflow))

# ---------------- momentum matrix (per-cell nuEff, unrelaxed basis) -------
nuEff_f = weights*nuEff[owner] + (1.0-weights)*nuEff[nei]
a_f = nuEff_f*magSf*deltaCoeffs                     # laplacian face coeff
up_own = phi_int >= 0.0
po = np.where(up_own, phi_int, 0.0)   # diag(own) += po   ; row(nei)-=po*U(own)
pn_ = np.where(up_own, 0.0, phi_int)  # row(own) += pn*U(nei); diag(nei)-=pn
diag_u = np.zeros(N)
np.add.at(diag_u, owner, po - a_f)
np.add.at(diag_u, nei, -pn_ - a_f)
diag_u += alpha*V
# bounded Sp: diag += -divphi*1 * V? (Sp(-fvc::surfaceIntegrate(phi),U):
# the matrix entry is -divphi (integrated: -divphi*V/V = per fvMatrix the Sp
# term adds -surfaceIntegrate(phi) to the DIAGONAL integrated -> -divphi_int
# where divphi_int = sum(phi_f) over faces (no V). see b17 rebuild: diag+=-divphi
divphi = np.zeros(N)
np.add.at(divphi, owner, phi_int)
np.add.at(divphi, nei, -phi_int)
np.add.at(divphi, bc_cell, phi_bvals)
diag_u += -divphi
# boundary convection/laplacian internalCoeffs (per component, equal):
iC_b = np.zeros(nBndF)
ob = phi_bvals >= 0.0
iC_b += np.where(ob, phi_bvals, 0.0)          # outflow conv -> diag
a_b = (nuEff[bc_cell])*bmag*bdel              # laplacian boundary (all faces)
iC_b += np.where(U_fixed, a_b, 0.0)           # only fixed-U contributes to iC
np.add.at(diag_u, bc_cell, np.where(ob, phi_bvals, 0.0))
np.add.at(diag_u, bc_cell[U_fixed], a_b[U_fixed])
# UNRELAXED momentum A (production operator basis: no relax on prodMomentum)
D_u = diag_u.copy()                            # diag + boundary iC included
A_u = D_u/V
rAU_u = 1.0/A_u
# RELAXED diagonal (primalPressureMobility basis, with dominance clamp)
sumOff = np.zeros(N)
np.add.at(sumOff, owner, np.abs(po - a_f) + np.abs(pn_ + a_f)*0)
sumOff = np.zeros(N)
np.add.at(sumOff, owner, np.abs(pn_) + a_f + np.where(up_own, 0.0, 0.0))
# simpler & exact: sumMagOffDiag = sum over faces |offdiag| shared by both rows
offMag = np.abs(pn_) + a_f    # |A[own,nei]| = |pn_| + a_f ... (pn_ and a_f both in the
#  same off-diagonal entry A[own,nei] = pn_ + a_f ; A[nei,own] = -po + a_f)
np.add.at(sumOff, owner, np.abs(pn_ + a_f))
np.add.at(sumOff, nei, np.abs(-po + a_f))
for pn_2 in ():
    pass
D0 = diag_u.copy()
D_rel = np.maximum(np.abs(D0), sumOff)/ALPHAREL
mob = 1.0/(D_rel/V)      # primalPressureMobility = rAtU = 1/A_rel
log("rAU_u avg=%.6e ; mobility(relaxed) avg=%.6e (b17 log rAU avg 2.84e-7)"
    % (rAU_u.mean(), mob.mean()))
# fvMatrix upper/lower of the momentum matrix (unrelaxed; convection+laplacian)
upper = pn_ + a_f    # A[own,nei]
lower = -po + a_f    # A[nei,own]

# ---------------- routing variants ----------------------------------------
kf_int = (weights*mob[owner] + (1.0-weights)*mob[nei])*magSf*deltaCoeffs
kf_b = mob[bc_cell]*bmag*bdel

def route_V0():
    b = np.zeros(NUNK)
    gU = np.zeros((N, 3))
    ft = Sf*g[:nIF][:, None]
    np.add.at(gU, owner, weights[:, None]*ft)
    np.add.at(gU, nei, (1.0-weights)[:, None]*ft)
    b[0:NV:3] = gU[:, 0]; b[1:NV:3] = gU[:, 1]; b[2:NV:3] = gU[:, 2]
    gP = np.zeros(N)
    np.add.at(gP, owner, -kf_int*g[:nIF])
    np.add.at(gP, nei, kf_int*g[:nIF])
    # boundary
    for f in range(nBndF):
        if U_fixed[f]: continue
        c = bc_cell[f]
        gf = g[nIF+f]
        b[3*c:3*c+3] += bSf[f]*gf
        if True:  # kf_b branch for non-fixed-U
            gP[c] += kf_b[f]*gf
    b[NV:] = gP
    return b

def route_V1(kf_sign_new=True):
    b = np.zeros(NUNK)
    hA = np.zeros((N, 3))
    gU = np.zeros((N, 3))
    gP = np.zeros(N)
    ft = Sf*g[:nIF][:, None]
    np.add.at(hA, owner, (ALPHAREL*rAU_u[owner]*weights)[:, None]*ft)
    np.add.at(hA, nei, (ALPHAREL*rAU_u[nei]*(1.0-weights))[:, None]*ft)
    np.add.at(gU, owner, ((1.0-ALPHAREL)*weights)[:, None]*ft)
    np.add.at(gU, nei, ((1.0-ALPHAREL)*(1.0-weights))[:, None]*ft)
    s = -1.0 if kf_sign_new else +1.0   # V1: b[P(o)] += kf*g ; old: b[P(o)] -= kf*g
    np.add.at(gP, owner, s*kf_int*g[:nIF])
    np.add.at(gP, nei, -s*kf_int*g[:nIF])
    # boundary: non-fixed-U patches
    am = ~U_fixed
    hA[bc_cell[am]] += (ALPHAREL*rAU_u[bc_cell[am]])[:, None]*(bSf[am]*g[nIF:][am][:, None])
    gU[bc_cell[am]] += (1.0-ALPHAREL)*(bSf[am]*g[nIF:][am][:, None])
    # kf_b on fixed-p (outlet) patches, non-fixed-U: sign + (matches operator)
    gP[bc_cell[am]] += np.where(np.ones(am.sum()), kf_b[am]*g[nIF:][am], 0.0)
    # deltaH^T (transpose of dH = -(upper/lower applied to dU)/V)
    outU = np.zeros((N, 3))
    Vinvo = 1.0/V[owner]; Vinn = 1.0/V[nei]
    for cix in range(3):
        hAc = hA[:, cix]
        add_o = -(upper*Vinvo)*hAc[owner]   # out[U(nei)] -= upper/V_o * hA[o]
        add_n = -(lower*Vinn)*hAc[nei]      # out[U(own)] -= lower/V_n * hA[n]
        np.add.at(outU[:, cix], nei, add_o)
        np.add.at(outU[:, cix], owner, add_n)
    b[0:NV:3] = gU[:, 0] + outU[:, 0]
    b[1:NV:3] = gU[:, 1] + outU[:, 1]
    b[2:NV:3] = gU[:, 2] + outU[:, 2]
    b[NV:] = gP
    return b

# ---------------- w_true + contractions ------------------------------------
hmap = {("D1", 3e-4): ("D1_h0.0003"), ("D1", 1e-3): ("D1_h0.001"),
        ("D2", 3e-4): ("D2_h0.0003"), ("D2", 1e-3): ("D2_h0.001"),
        ("D3", 3e-4): ("D3_h0.0003"), ("D3", 1e-3): ("D3_h0.001")}
def wtrue(tag, h):
    Up = np.loadtxt(ST + "/wstate_%s_p_U.mtx" % tag)
    Um = np.loadtxt(ST + "/wstate_%s_m_U.mtx" % tag)
    pp = np.loadtxt(ST + "/wstate_%s_p_p.mtx" % tag)
    pm_ = np.loadtxt(ST + "/wstate_%s_m_p.mtx" % tag)
    w = np.zeros(NUNK)
    w[0:NV:3] = ((Up-Um)/(2*h))[:, 0]
    w[1:NV:3] = ((Up-Um)/(2*h))[:, 1]
    w[2:NV:3] = ((Up-Um)/(2*h))[:, 2]
    w[NV:] = (pp-pm_)/(2*h)
    return w

# thermal_gauge and C
dJdT = np.loadtxt(ST + "/wstate_dJdT.mtx")
z3 = np.loadtxt(E15 + "/stageB6_rxc_z_analytic.mtx").reshape(3, N)

def C_field():
    """thermalDiffusionDerivativeDTCell per-cell field (code formula)."""
    out = np.zeros(N)
    # e_f direction derivative uses delta vector (cell-center diff)
    ef = d_vec/np.maximum(np.linalg.norm(d_vec, axis=1), 1e-30)[:, None]
    def dir_deriv(cells):
        D6 = dDTDxh[cells]  # (M,6) xx xy xz yy yz zz
        ex, ey, ez = ef[:, 0], ef[:, 1], ef[:, 2]
        return (ex*ex*D6[:, 0] + ex*ey*D6[:, 1] + ex*ez*D6[:, 2]
                + ey*ex*D6[:, 1] + ey*ey*D6[:, 3] + ey*ez*D6[:, 4]
                + ez*ex*D6[:, 2] + ez*ey*D6[:, 4] + ez*ez*D6[:, 5])
    dd_o = dir_deriv(owner); dd_n = dir_deriv(nei)
    base = (Tb[owner]-Tb[nei])*(T[owner]-T[nei])*deltaCoeffs*magSf
    np.add.at(out, owner, -weights*dd_o*base)
    np.add.at(out, nei, -(1.0-weights)*dd_n*base)
    # boundary faces: fixedValue T (inlet) has snGrad != 0; e_f from
    # boundary delta vector (cell -> face centre)
    bCf = np.zeros((nBndF, 3))
    for f in range(nBndF):
        fv = pts[faces[nIF + f]]
        bCf[f] = fv.mean(axis=0)
    bDvec = bCf - Ccen[bc_cell]
    bef = bDvec/np.maximum(np.linalg.norm(bDvec, axis=1), 1e-30)[:, None]
    def bdir_deriv():
        D6 = dDTDxh[bc_cell]
        ex, ey, ez = bef[:, 0], bef[:, 1], bef[:, 2]
        return (ex*ex*D6[:, 0] + ex*ey*D6[:, 1] + ex*ez*D6[:, 2]
                + ey*ex*D6[:, 1] + ey*ey*D6[:, 3] + ey*ez*D6[:, 4]
                + ez*ex*D6[:, 2] + ez*ey*D6[:, 4] + ez*ez*D6[:, 5])
    snb = (T_bvals - T[bc_cell])*bdel
    ddb = bdir_deriv()
    np.add.at(out, bc_cell, ddb*Tb[bc_cell]*snb*bmag)
    return out

def C_direct(z):
    return float(np.dot(C_field(), z))

# ---------------- validation V0 vs exported rhs ----------------------------
import scipy.io as sio
vexp = np.asarray(sio.mmread(E15 + "/explicitRhs_thermalCoupling.mtx")).ravel()
m = vexp.reshape(N, 4)
bexp = np.zeros(NUNK)
bexp[0:NV:3] = m[:, 0]; bexp[1:NV:3] = m[:, 1]; bexp[2:NV:3] = m[:, 2]
bexp[NV:] = m[:, 3]
bV0 = route_V0()
dn = np.linalg.norm(bV0-bexp); rn = np.linalg.norm(bexp)
log("V0 reconstruction vs exported rhs: relL2=%.3e  maxabs=%.3e"
    % (dn/rn, np.abs(bV0-bexp).max()))

bV1 = route_V1(True)
bV1k = route_V1(False)
log("variants ready (t=%.1fs)" % (time.time()-t0))

# ---------------- contractions vs gauge ------------------------------------
res = {}
for d, h in [("D1", 1e-3), ("D2", 1e-3), ("D3", 1e-3), ("D1", 3e-4),
             ("D2", 3e-4), ("D3", 3e-4)]:
    tag = hmap[(d, h)]
    w = wtrue(tag, h)
    Tp = np.loadtxt(ST + "/wstate_%s_p_T.mtx" % tag)
    Tm = np.loadtxt(ST + "/wstate_%s_m_T.mtx" % tag)
    tg = float(np.dot(dJdT, (Tp-Tm)/(2*h)))
    fd = FD_J[d]
    A = fd - tg
    Ck = C_direct(z3[["D1", "D2", "D3"].index(d)])
    row = {
        "V0": float(bV0 @ w), "V1": float(bV1 @ w), "V1k": float(bV1k @ w),
        "fd": fd, "thermal_g": tg, "A": A, "C": Ck,
        "errA_V0": abs(bV0 @ w - A)/abs(A),
        "errA_V1": abs(bV1 @ w - A)/abs(A),
        "errFDC_V1": abs(bV1 @ w - (fd - Ck))/abs(fd - Ck),
        "errFDC_V0": abs(bV0 @ w - (fd - Ck))/abs(fd - Ck),
        "errFDC_V1k": abs(bV1k @ w - (fd - Ck))/abs(fd - Ck),
    }
    res["%s_h%g" % (d, h)] = row
    log("%s h=%g: V0=%.6f V1=%.6f V1k=%.6f | tg=%.6f A=%.6f C=%.6f "
        "FD-C=%.6f | errA(V0/V1)=%.3f/%.3f errFD-C(V1/V1k)=%.3f/%.3f"
        % (d, h, row["V0"], row["V1"], row["V1k"], tg, A, Ck, fd-Ck,
           row["errA_V0"], row["errA_V1"], row["errFDC_V1"], row["errFDC_V1k"]))

json.dump(res, open("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/"
                    "BFINAL-018/cycle-1/b18_folding_lab.json", "w"), indent=1)
log("LAB_DONE t=%.1fs" % (time.time()-t0))

# ---------------- part 2: Gx (flux direct-alpha), pb pressure-row, crosschecks
HbyA_lab = None
def build_dHbyA_drAU():
    # rxPressureRowT basis: RELAXED matrix (mobility semantics), validated by B17
    rAU_rel = mob  # 1/A_rel
    # HbyA = rAU_rel*H_rel; H_rel = (S0 + (D_rel-D0)*U + off*U)/V with S0 incl -gradp, dev2
    # off*U part:
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, (pn_ + a_f)[:, None]*U[nei])
    np.add.at(offU, nei, (-po + a_f)[:, None]*U[owner])
    # source S0*V: -grad(p)*V + convective inflow BC + lap fixed-U + (-div dev2)*V
    w_ = weights
    pf_int = w_*p[owner] + (1-w_)*p[nei]
    gradp = np.zeros((N, 3))
    tmp = Sf*pf_int[:, None]
    np.add.at(gradp, owner, tmp); np.add.at(gradp, nei, -tmp)
    np.add.at(gradp, bc_cell, bSf*p_bvals[:, None])
    gradp /= V[:, None]
    src_exp = -gradp
    # dev2 term (uniform-nu approximation as b17; small)
    Uf_int = w_[:, None]*U[owner] + (1-w_)[:, None]*U[nei]
    gU = np.zeros((N, 3, 3))
    tmpU = Sf[:, None]*Uf_int[:, :, None]
    np.add.at(gU, owner, tmpU); np.add.at(gU, nei, -tmpU)
    np.add.at(gU, bc_cell, bSf[:, None]*U_bvals[:, :, None])
    gU /= V[:, None, None]
    dev2T = np.swapaxes(gU, 1, 2).copy()
    tr = np.trace(gU, axis1=1, axis2=2)
    dev2T[:, 0, 0] -= 2.0/3.0*tr; dev2T[:, 1, 1] -= 2.0/3.0*tr; dev2T[:, 2, 2] -= 2.0/3.0*tr
    nuf_b = nuEff[bc_cell]
    tens = np.einsum("c,cij->cij", np.ones(N), dev2T)*0  # uniform nu handled below
    nuint = 0.5*(nuEff[owner] + nuEff[nei])
    Tf_int = w_[:, None, None]*(nuint[:, None, None]*dev2T[owner]) \
           + (1-w_)[:, None, None]*(nuint[:, None, None]*dev2T[nei])
    Tf_bnd = nuEff[bc_cell][:, None, None]*dev2T[bc_cell]
    divdev = np.zeros((N, 3))
    tmpD = np.einsum("fi,fij->fj", Sf, Tf_int)
    np.add.at(divdev, owner, tmpD); np.add.at(divdev, nei, -tmpD)
    np.add.at(divdev, bc_cell, np.einsum("fi,fij->fj", bSf, Tf_bnd))
    divdev /= V[:, None]
    src_exp += -divdev
    src_bcomp = np.zeros((N, 3))
    ib = ~ob
    np.add.at(src_bcomp, bc_cell[ib], -phi_bvals[ib, None]*U_bvals[ib])
    fb = U_fixed
    a_b2 = nuEff[bc_cell]*bmag*bdel
    np.add.at(src_bcomp, bc_cell[fb], a_b2[fb, None]*U_bvals[fb])
    src_vec = src_exp*V[:, None] + src_bcomp
    src_vec += (1.0-ALPHAREL)*D_rel[:, None]*U
    H = (src_vec - offU)/V[:, None]
    HbyA = rAU_rel[:, None]*H
    drAU = -rAU_rel*rAU_rel/ALPHAREL
    dHbyA = (rAU_rel/ALPHAREL)[:, None]*((1.0-ALPHAREL)*U - HbyA)
    return dHbyA, drAU, HbyA

dHbyA, drAU, HbyA = build_dHbyA_drAU()
# sanity: HbyA ~ U + rAU*gradp  (fixed point)
err = np.abs(HbyA - (U + mob[:, None]*0)).max()
log("HbyA-U maxabs=%.4e (should be ~ rAU*|grad p| scale, not 0)" % err)

g0_int = magSf*deltaCoeffs*(p[nei] - p[owner])

def Gx_of(zdir):
    wc = dAlphaDxh*zdir
    dHw = dHbyA*wc[:, None]
    phx_int = np.einsum("fi,fi->f", Sf, weights[:, None]*dHw[owner]
                        + (1.0-weights)[:, None]*dHw[nei])
    gf_int = weights*drAU[owner]*wc[owner] + (1.0-weights)*drAU[nei]*wc[nei]
    phx_int -= gf_int*g0_int
    val = float(np.dot(g[:nIF], phx_int))
    am = U_assign
    phx_b = np.einsum("fi,fi->f", bSf[am], dHw[bc_cell[am]])
    val += float(np.dot(g[nIF:][am], phx_b))
    return val

# pb pressure-row term: -pb^T anRPd_k  (rxc P-part = J_P*w = R_P,x)
pb_i, _, _ = read_field(S + "/pb"); pb = pb_i
Ub_i, _, _ = read_field(S + "/Ub"); Ub = Ub_i.reshape(N, 3)
rxc = np.loadtxt(E15 + "/stageB6_rxc_analytic.mtx")
projJ_b16 = {"D1": -0.01593219715949202, "D2": -0.02416470832109031,
             "D3": 0.003160577650589706}

log("\n=== part 2: full decomposition per direction (h-independent pieces) ===")
tbl = {}
for k, d in enumerate(["D1", "D2", "D3"]):
    anRUa = rxc[k*NUNK:(k+1)*NUNK][:NV].reshape(N, 3)
    anRPd = rxc[k*NUNK:(k+1)*NUNK][NV:]
    zk = z3[k]
    Ck = C_direct(zk)
    Gxk = Gx_of(zk)
    Umom = float(np.einsum("ij,ij->", Ub, anRUa))     # -lambda_U^T R_U,x
    Prow = float(np.dot(pb, anRPd))                    # -lambda_P^T R_P,x
    fd = FD_J[d]
    row = dict(C=Ck, Gx=Gxk, Umom=-Umom, Prow=-Prow,
               current_prod=Ck + (-Umom), b16_projJ=projJ_b16[d], fd=fd)
    tbl[d] = row
    log("%s: C=%.6f Gx=%.6f | -Ub.ru=%.6f -pb.rp=%.6f | C-Ub.ru=%.6f "
        "(b16 projJ=%.6f) | full=C+Ub+pb+Gx=%.6f vs FD=%.6f"
        % (d, Ck, Gxk, -Umom, -Prow, Ck-Umom, projJ_b16[d],
           Ck - Umom - Prow + Gxk, fd))

json.dump(tbl, open("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/"
                    "BFINAL-018/cycle-1/b18_decomp.json", "w"), indent=1)
log("PART2_DONE t=%.1fs" % (time.time()-t0))

# ---------------- part 3: validate C against the written fsenshMeanT --------
fsensh_w = read_field(S + "/fsenshMeanT")[0]
thermDiff_written = fsensh_w + dAlphaDxh*np.einsum("ij,ij->i", U, Ub)*V
log("\n=== part 3: C validation vs written thermalDiffusionDerivative ===")
myC = C_field()
log("per-cell thermalDiff: written-vs-mine maxabs=%.4e relL2=%.4e (|written|=%.4e)"
    % (np.abs(thermDiff_written - myC).max(),
       np.linalg.norm(thermDiff_written - myC)/np.linalg.norm(thermDiff_written),
       np.linalg.norm(thermDiff_written)))
for k, d in enumerate(["D1", "D2", "D3"]):
    log("  C_written(%s)=%.8f  C_mine=%.8f"
        % (d, float(np.dot(thermDiff_written, z3[k])), C_direct(z3[k])))

# ---------------- part 4: C validation on the single-pass b15_export case ---
S15 = "/home/ys/dsH/b15_export/1"
fsensh15 = read_field(S15 + "/fsenshMeanT")[0]
Ub15_i, uval15, _ = read_field(S15 + "/Ub")
Ub15 = (Ub15_i.reshape(N, 3) if Ub15_i is not None
        else np.tile(np.atleast_1d(uval15), (N, 1)))
thermDiff15 = fsensh15 + dAlphaDxh*np.einsum("ij,ij->i", U, Ub15)*V
log("\n=== part 4: C validation vs b15 single-pass fsenshMeanT ===")
log("b15 thermalDiff vs mine: maxabs=%.4e relL2=%.4e (|b15|=%.4e)"
    % (np.abs(thermDiff15 - myC).max(),
       np.linalg.norm(thermDiff15 - myC)/np.linalg.norm(thermDiff15),
       np.linalg.norm(thermDiff15)))
for k, d in enumerate(["D1", "D2", "D3"]):
    log("  C_b15(%s)=%.8f  C_mine=%.8f  ratio=%.4f"
        % (d, float(np.dot(thermDiff15, z3[k])), C_direct(z3[k]),
           float(np.dot(thermDiff15, z3[k]))/C_direct(z3[k])))

# ---------------- part 5: per-cell diagnosis of the C mismatch -------------
supp = np.abs(thermDiff15) > 1e-12
log("\n=== part 5: per-cell C mismatch diagnosis ===")
log("support b15: %d cells ; mine: %d cells ; overlap: %d"
    % (int(supp.sum()), int((np.abs(myC) > 1e-12).sum()),
       int((supp & (np.abs(myC) > 1e-12)).sum())))
d_ = thermDiff15 - myC
log("b15-field stats: max=%.4e L2=%.4e ; myC: max=%.4e L2=%.4e ; diff L2=%.4e"
    % (np.abs(thermDiff15).max(), np.linalg.norm(thermDiff15),
       np.abs(myC).max(), np.linalg.norm(myC), np.linalg.norm(d_)))
both = supp & (np.abs(myC) > 1e-12)
ratio = thermDiff15[both]/myC[both]
log("ratio stats: median=%.4f p10=%.4f p90=%.4f" %
    (float(np.median(ratio)), float(np.percentile(ratio, 10)),
     float(np.percentile(ratio, 90))))
# internal-only and boundary-only reconstruction
out_int = np.zeros(N)
ef = d_vec/np.maximum(np.linalg.norm(d_vec, axis=1), 1e-30)[:, None]
def dir_deriv2(cells, evec):
    D6 = dDTDxh[cells]
    ex, ey, ez = evec[:, 0], evec[:, 1], evec[:, 2]
    return (ex*ex*D6[:, 0] + ex*ey*D6[:, 1] + ex*ez*D6[:, 2]
            + ey*ex*D6[:, 1] + ey*ey*D6[:, 3] + ey*ez*D6[:, 4]
            + ez*ex*D6[:, 2] + ez*ey*D6[:, 4] + ez*ez*D6[:, 5])
ddo = dir_deriv2(owner, ef); ddn = dir_deriv2(nei, ef)
base = (Tb[owner]-Tb[nei])*(T[owner]-T[nei])*deltaCoeffs*magSf
np.add.at(out_int, owner, -weights*ddo*base)
np.add.at(out_int, nei, -(1.0-weights)*ddn*base)
log("internal-only vs b15 on support: relL2=%.4e"
    % (np.linalg.norm((out_int - thermDiff15)[supp])
       /np.linalg.norm(thermDiff15[supp])))
log("PART5_DONE")

# ---------------- part 6: reverse-engineer the production C formula --------
log("\n=== part 6: top cells forensics ===")
o = np.argsort(np.abs(thermDiff15))[::-1][:8]
# faces of each cell
cellfaces = {}
for fi in range(nIF):
    cellfaces.setdefault(owner[fi], []).append(fi)
    cellfaces.setdefault(nei[fi], []).append(fi)
for c in o:
    log("cell %d Ccen=%s b15=%.6e mine=%.6e nFaces=%d"
        % (c, np.round(Ccen[c], 6), thermDiff15[c], myC[c],
           len(cellfaces.get(c, []))))
    for fi in cellfaces.get(c, [])[:8]:
        ow, nn = owner[fi], nei[fi]
        d6o = dDTDxh[ow]; d6n = dDTDxh[nn]
        log("   f%d o=%d n=%d w=%.4f Tb=(%.3e,%.3e) T=(%.1f,%.1f) "
            "dDo=(%.2e %.2e %.2e %.2e %.2e %.2e)"
            % (fi, ow, nn, weights[fi], Tb[ow], Tb[nn], T[ow], T[nn],
               d6o[0], d6o[1], d6o[2], d6o[3], d6o[4], d6o[5]))

# ---------------- part 7: was b15's in-memory Ub the stale lambda15? --------
lam15 = np.loadtxt(E15 + "/stageB6_lambda.mtx")
U15 = lam15[:NV].reshape(N, 3)
mom15 = -dAlphaDxh*np.einsum("ij,ij->i", U, U15)*V
log("\n=== part 7: b15 fsensh decomposition test ===")
for k, d in enumerate(["D1", "D2", "D3"]):
    a = float(np.dot(myC, z3[k]))
    b = float(np.dot(fsensh15, z3[k]))
    c = float(np.dot(mom15, z3[k]))
    log("%s: fsensh15=%.6f | myC=%.6f mom(lambda15)=%.6f | myC+mom=%.6f | "
        "fsensh-myC=%.6f (vs mom %.6f)"
        % (d, b, a, c, a+c, b-a, c))

# ---------------- part 8: rebuild thermal matrix A_T, adjudicate C ---------
import scipy.sparse as sp
import scipy.sparse.linalg as spl
DT6 = dDTDxh
DTe6 = read_field(S + "/DTEffective")[0].reshape(N, 6)
# face direction conductivity projection
efl = d_vec/np.maximum(np.linalg.norm(d_vec, axis=1), 1e-30)[:, None]
def proj(D6c):
    ex, ey, ez = efl[:, 0], efl[:, 1], efl[:, 2]
    return (ex*ex*D6c[:, 0] + ex*ey*D6c[:, 1] + ex*ez*D6c[:, 2]
            + ey*ex*D6c[:, 1] + ey*ey*D6c[:, 3] + ey*ez*D6c[:, 4]
            + ez*ex*D6c[:, 2] + ez*ey*D6c[:, 4] + ez*ez*D6c[:, 5])
DDo = proj(DTe6[owner]); DDn = proj(DTe6[nei])
kfT = (weights*DDo + (1.0-weights)*DDn)*magSf*deltaCoeffs
# convection (bounded upwind, phiThermal internal + boundary)
poT = np.where(phiTh_int >= 0.0, phiTh_int, 0.0)
pnT = np.where(phiTh_int >= 0.0, 0.0, phiTh_int)
diagT = np.zeros(N); offT_own = np.zeros(nIF)  # offT: A[own,nei]
np.add.at(diagT, owner, poT)
np.add.at(diagT, nei, -pnT)
offT_p = pn_.shape  # unused
A_off_own = pnT.copy()   # A[own,nei] = pnT (from conv)
A_off_nei = -poT.copy()  # A[nei,own] = -poT
# laplacian (subtract: TEqn = div - laplacian == Q)
# laplacian matrix: diag += kfT, off -= kfT -> with the minus in TEqn:
# A += -lap: diag -= kfT, off += kfT
np.add.at(diagT, owner, kfT)
np.add.at(diagT, nei, kfT)
A_off_own -= kfT
A_off_nei -= kfT
# bounded Sp: diag += -div(phiThermal)
divphiT = np.zeros(N)
np.add.at(divphiT, owner, phiTh_int)
np.add.at(divphiT, nei, -phiTh_int)
phTh_b = np.zeros(nBndF)
o = 0
for pn_3 in patch_names:
    nFp = nFacesPatch[pn_3]
    btypePhi, bvphi = phTh_b_ = read_field(S + "/phiThermal")[2][pn_3]
    if isinstance(bvphi, float): phTh_b[o:o+nFp] = bvphi
    elif bvphi is not None: phTh_b[o:o+nFp] = bvphi
    o += nFp
np.add.at(divphiT, bc_cell, phTh_b)
diagT += -divphiT
# boundary: outflow conv -> diag += phi_b ; inflow -> source -phi_b*T_b
obT = phTh_b >= 0.0
np.add.at(diagT, bc_cell[obT], phTh_b[obT])
srcT = np.zeros(N)
np.add.at(srcT, bc_cell[~obT], -phTh_b[~obT]*T_bvals[~obT])
# laplacian boundary: fixedValue T (inlet + hotInlet?): diag += kf_b, src += kf_b*T_b
DTb6 = DTe6[bc_cell]
efb_all = np.zeros((nBndF, 3))
for f in range(nBndF):
    fv = pts[faces[nIF + f]]
    efb_all[f] = (fv.mean(axis=0) - Ccen[bc_cell[f]])
efb_all /= np.maximum(np.linalg.norm(efb_all, axis=1), 1e-30)[:, None]
def projb(D6c):
    ex, ey, ez = efb_all[:, 0], efb_all[:, 1], efb_all[:, 2]
    return (ex*ex*D6c[:, 0] + ex*ey*D6c[:, 1] + ex*ez*D6c[:, 2]
            + ey*ex*D6c[:, 1] + ey*ey*D6c[:, 3] + ey*ez*D6c[:, 4]
            + ez*ex*D6c[:, 2] + ez*ey*D6c[:, 4] + ez*ez*D6c[:, 5])
kfTb = projb(DTb6)*bmag*bdel
T_fixed = np.zeros(nBndF, dtype=bool)
o = 0
for pn_4 in patch_names:
    nFp = nFacesPatch[pn_4]
    btypeT, _ = T_b[pn_4]
    T_fixed[o:o+nFp] = (btypeT == "fixedValue")
    o += nFp
np.add.at(diagT, bc_cell[T_fixed], kfTb[T_fixed])
np.add.at(srcT, bc_cell[T_fixed], kfTb[T_fixed]*T_bvals[T_fixed])
# assemble CSR
rows = np.concatenate([owner, nei, np.arange(N)])
cols = np.concatenate([nei, owner, np.arange(N)])
vals = np.concatenate([A_off_own, A_off_nei, diagT])
AT = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
log("\n=== part 8: thermal operator adjudication ===")
# (a) residual check A_T T - src vs Q
Qv = np.zeros(N)  # Q field is uniform 0 (verified in file)
resT = AT @ T - srcT - Qv
log("A_T*T - src - Q: maxabs=%.4e  |Q|max=%.4e  relL2=%.3e"
    % (np.abs(resT).max(), np.abs(Qv).max(),
       np.linalg.norm(resT)/np.linalg.norm(Qv)))
# (b) adjoint check: A_T^T Tb_written vs dJdT
resTb = (AT.T @ Tb) - dJdT
log("A_T^T*Tb_written - dJdT: maxabs=%.4e |dJdT|max=%.4e relL2=%.3e"
    % (np.abs(resTb).max(), np.abs(dJdT).max(),
       np.linalg.norm(resTb)/np.linalg.norm(dJdT)))
# (c) definitive C by direct solve: A_T T' = -(A_x dx) T
ddD_o = proj(DT6[owner]); ddD_n = proj(DT6[nei])
for k, d in enumerate(["D1", "D2", "D3"]):
    zk = z3[k]
    dkf = (weights*ddD_o*zk[owner] + (1.0-weights)*ddD_n*zk[nei]) \
          *magSf*deltaCoeffs
    rhs = np.zeros(N)
    np.add.at(rhs, owner, +dkf*(T[owner] - T[nei]))
    np.add.at(rhs, nei, -dkf*(T[owner] - T[nei]))
    # boundary fixed-T faces: d(kf_b)* (T_c - T_b) enters rhs; dD_b nonzero?
    dkb = projb(DT6[bc_cell])*zk[bc_cell]*bmag*bdel
    np.add.at(rhs, bc_cell[T_fixed], -dkb[T_fixed]*(T[bc_cell] - T_bvals)[T_fixed])
    Tp = spl.spsolve(AT.tocsc(), rhs)
    Cfd = float(np.dot(dJdT, Tp))
    log("%s: C_direct-solve=%.8f | myC=%.8f | fsensh15=%.8f"
        % (d, Cfd, C_direct(zk), float(np.dot(fsensh15, z3[k]))))

# ---------------- part 9: pairing sanity + per-cell ratio correlation -------
gPDw = read_field(S15 + "/gsenshPressureDrop")[0]
log("\n=== part 9: pairing sanity ===")
for k, d in enumerate(["D1", "D2", "D3"]):
    log("%s: gsenshPD(b15)·z=%.6f (B17 XID -0.0562/0.2977/-0.0363) | "
        "thermDiff15·z=%.6f myC·z=%.6f ratio=%.3f"
        % (d, float(np.dot(gPDw, z3[k])), float(np.dot(thermDiff15, z3[k])),
           C_direct(z3[k]),
           float(np.dot(thermDiff15, z3[k]))/C_direct(z3[k])))
# deterministic D directions (mirror validateStageB2GradientAmplitude.H)
C_ = mesh_cc = np.zeros((N, 3))
Cx = np.zeros(N); Cy = np.zeros(N); Cz = np.zeros(N)
# recompute cell centers already in Ccen
Lx = Ccen[:, 0].max() - Ccen[:, 0].min()
Ly = Ccen[:, 1].max() - Ccen[:, 1].min()
Lz = Ccen[:, 2].max() - Ccen[:, 2].min()
xx = Ccen[:, 0]/Lx; yy = Ccen[:, 1]/Ly; zz = Ccen[:, 2]/Lz
xr = read_field(S + "/x")[0]
Ddir = np.zeros((3, N))
Ddir[0] = np.sin(2*np.pi*xx)*np.cos(np.pi*yy) + 0.5*np.sin(3*np.pi*zz)
Ddir[1] = np.sin(4*np.pi*xx) + 0.3*np.cos(2*np.pi*xx)
Ddir[2] = np.cos(2*np.pi*yy)*np.sin(2*np.pi*zz)
for j in range(3):
    m_ = Ddir[j].max()
    Ddir[j] /= m_
log("D-pairing: myC·D = %.6e %.6e %.6e (B13 probe thermal: "
    "-1.442e-4 -4.759e-3 3.698e-3)"
    % tuple(float(np.dot(myC, Ddir[j])) for j in range(3)))
log("D-pairing: thermDiff15·D = %.6e %.6e %.6e"
    % tuple(float(np.dot(thermDiff15, Ddir[j])) for j in range(3)))
# per-cell ratio correlation
both = (np.abs(thermDiff15) > 1e-10) & (np.abs(myC) > 1e-10)
rat = thermDiff15[both]/myC[both]
log("ratio stats: median=%.3f mean=%.3f min=%.3f max=%.3f n=%d"
    % (float(np.median(rat)), float(rat.mean()), float(rat.min()),
       float(rat.max()), int(both.sum())))

# ---------------- part 10: continuity test of flux-map variants ------------
def dH_apply(wU):
    """momentum H-tangent: dH_c = -(sum_f coeff*wU_other)/V (deltaH^T pair)."""
    out = np.zeros((N, 3))
    Vinvo = 1.0/V[owner]; Vinn = 1.0/V[nei]
    for cix in range(3):
        wUc = wU[:, cix]
        add_o = -(upper*Vinvo)*wUc[nei]   # dH[own] -= upper/V_o * wU[nei]
        add_n = -(lower*Vinn)*wUc[owner]
        np.add.at(out[:, cix], owner, add_o)
        np.add.at(out[:, cix], nei, add_n)
    return out

def flux_fwd(wU, wp, variant):
    dH = dH_apply(wU)
    if variant == "relaxed":      # BFINAL-003 operator map
        cellf_o = ALPHAREL*rAU_u[owner][:, None]*dH[owner] \
            + (1.0-ALPHAREL)*wU[owner]
        cellf_n = ALPHAREL*rAU_u[nei][:, None]*dH[nei] \
            + (1.0-ALPHAREL)*wU[nei]
        dphi = np.einsum("fi,fi->f", Sf, weights[:, None]*cellf_o
                         + (1.0-weights)[:, None]*cellf_n)
        dphi += kf_int*(wp[owner] - wp[nei])
    elif variant == "unrelaxed":
        cellf_o = rAU_u[owner][:, None]*dH[owner]
        cellf_n = rAU_u[nei][:, None]*dH[nei]
        dphi = np.einsum("fi,fi->f", Sf, weights[:, None]*cellf_o
                         + (1.0-weights)[:, None]*cellf_n)
        dphi += kf_int*(wp[owner] - wp[nei])
    elif variant == "v0code":
        interp = weights[:, None]*wU[owner] + (1.0-weights)[:, None]*wU[nei]
        dphi = np.einsum("fi,fi->f", Sf, interp)
        dphi -= kf_int*(wp[nei] - wp[owner])
    # boundary (assignable = outlet): dphi_b (v0code uses its own sign)
    db = np.zeros(nBndF)
    am = U_assign
    if variant == "relaxed":
        cellf_c = (ALPHAREL*rAU_u[bc_cell][:, None]*dH[bc_cell]
                   + (1.0-ALPHAREL)*wU[bc_cell])
        db[am] = np.einsum("fi,fi->f", bSf[am], cellf_c[am]) \
                 + kf_b[am]*wp[bc_cell[am]]
    elif variant == "unrelaxed":
        cellf_c = rAU_u[bc_cell][:, None]*dH[bc_cell]
        db[am] = np.einsum("fi,fi->f", bSf[am], cellf_c[am]) \
                 + kf_b[am]*wp[bc_cell[am]]
    elif variant == "v0code":
        db[am] = np.einsum("fi,fi->f", bSf[am], wU[bc_cell[am]]) \
                 + kf_b[am]*wp[bc_cell[am]]
    return dphi, db

def div_of(dphi, db):
    d = np.zeros(N)
    np.add.at(d, owner, dphi)
    np.add.at(d, nei, -dphi)
    np.add.at(d, bc_cell, db)
    return d

def phix_faces(zdir):
    wc = dAlphaDxh*zdir
    dHw = dHbyA*wc[:, None]
    phx = np.einsum("fi,fi->f", Sf, weights[:, None]*dHw[owner]
                    + (1.0-weights)[:, None]*dHw[nei])
    gf = weights*drAU[owner]*wc[owner] + (1.0-weights)*drAU[nei]*wc[nei]
    phx -= gf*g0_int
    phxb = np.zeros(nBndF)
    am = U_assign
    phxb[am] = np.einsum("fi,fi->f", bSf[am], dHw[bc_cell[am]])
    return phx, phxb

log("\n=== part 10: continuity test div(M w_true + phi_x) ===")
tag = hmap[("D1", 1e-3)]
w1 = wtrue(tag, 1e-3)
wU1 = np.stack([w1[0:NV:3], w1[1:NV:3], w1[2:NV:3]], axis=1)
wp1 = w1[NV:]
phx, phxb = phix_faces(z3[0])
div_phix = div_of(phx, phxb)
for variant in ["relaxed", "unrelaxed", "v0code"]:
    dphi, db = flux_fwd(wU1, wp1, variant)
    dmw = div_of(dphi, db)
    tot = dmw + div_phix
    log("%s: |div(Mw)+div(phix)|L2=%.4e  |div(phix)|L2=%.4e "
        "|div(Mw)|L2=%.4e  max|tot|=%.3e"
        % (variant, np.linalg.norm(tot), np.linalg.norm(div_phix),
           np.linalg.norm(dmw), np.abs(tot).max()))
log("PART10_DONE")

# ---------------- part 11: routing block decomposition ---------------------
def route_blocks(kf_sign_new=True, use_hA=True, use_direct=True,
                 arel=ALPHAREL, rau=None):
    if rau is None: rau = rAU_u
    b = np.zeros(NUNK)
    hA = np.zeros((N, 3)); gU = np.zeros((N, 3)); gP = np.zeros(N)
    ft = Sf*g[:nIF][:, None]
    if use_hA:
        np.add.at(hA, owner, (arel*rau[owner]*weights)[:, None]*ft)
        np.add.at(hA, nei, (arel*rau[nei]*(1.0-weights))[:, None]*ft)
    if use_direct:
        np.add.at(gU, owner, ((1.0-arel)*weights)[:, None]*ft)
        np.add.at(gU, nei, ((1.0-arel)*(1.0-weights))[:, None]*ft)
    s = 1.0 if kf_sign_new else -1.0
    np.add.at(gP, owner, s*kf_int*g[:nIF])
    np.add.at(gP, nei, -s*kf_int*g[:nIF])
    am = ~U_fixed
    if use_hA:
        hA[bc_cell[am]] += (arel*rau[bc_cell[am]])[:, None]*(
            bSf[am]*g[nIF:][am][:, None])
    if use_direct:
        gU[bc_cell[am]] += (1.0-arel)*(bSf[am]*g[nIF:][am][:, None])
    gP[bc_cell[am]] += kf_b[am]*g[nIF:][am]
    outU = np.zeros((N, 3))
    if use_hA:
        Vinvo = 1.0/V[owner]; Vinn = 1.0/V[nei]
        for cix in range(3):
            hAc = hA[:, cix]
            np.add.at(outU[:, cix], nei, -(upper*Vinvo)*hAc[owner])
            np.add.at(outU[:, cix], owner, -(lower*Vinn)*hAc[nei])
    b[0:NV:3] = gU[:, 0] + outU[:, 0]
    b[1:NV:3] = gU[:, 1] + outU[:, 1]
    b[2:NV:3] = gU[:, 2] + outU[:, 2]
    b[NV:] = gP
    return b, np.linalg.norm(gU)*0+np.array([0.0]), np.linalg.norm(outU), np.linalg.norm(gP)

log("\n=== part 11: block decomposition (hA/direct/kf on/off) ===")
w1 = wtrue(hmap[("D1", 1e-3)], 1e-3)
w2 = wtrue(hmap[("D2", 1e-3)], 1e-3)
w3 = wtrue(hmap[("D3", 1e-3)], 1e-3)
tgt = np.array([FD_J["D1"] - C_direct(z3[0]) - Gx_of(z3[0]),
                FD_J["D2"] - C_direct(z3[1]) - Gx_of(z3[1]),
                FD_J["D3"] - C_direct(z3[2]) - Gx_of(z3[2])])
log("target FD-C-Gx = %s" % np.round(tgt, 5))
for name, kw in [
    ("V1 full       ", dict()),
    ("no kf         ", dict(use_hA=True, use_direct=True)),
    ("hA only       ", dict(use_direct=False)),
    ("direct+kf only", dict(use_hA=False)),
    ("kf only       ", dict(use_hA=False, use_direct=False)),
    ("arel=1 unrelax", dict(arel=1.0)),
    ("hA+kf, nodir  ", dict(use_direct=False)),
]:
    if name == "no kf":
        continue
    b, _, nH, nP = route_blocks(**kw)
    if name == "kf only":
        b, _, nH, nP = route_blocks(use_hA=False, use_direct=False)
    vals = np.array([b @ w1, b @ w2, b @ w3])
    log("%s: D1=%+.5f D2=%+.5f D3=%+.5f" % (name, *vals))
# separate: hA-transpose output only
b_hA = route_blocks(use_direct=False)[0]
b_dir = route_blocks(use_hA=False)[0]
b_kf = route_blocks(use_hA=False, use_direct=False)[0]
for nm, b in [("hA^T part", b_hA), ("direct^T part", b_dir), ("kf^T part", b_kf)]:
    log("   %s: D1=%+.5f D2=%+.5f D3=%+.5f"
        % (nm, b @ w1, b @ w2, b @ w3))

# ---------------- part 12: round-2 state b_TC reconstruction vs export ----
ST18 = "/home/ys/dsH/b18_fix/stageB2"
U2 = np.loadtxt(ST18 + "/wstate_baseline_U.mtx")
p2 = np.loadtxt(ST18 + "/wstate_baseline_p.mtx")
T2 = np.loadtxt(ST18 + "/wstate_baseline_T.mtx")
dJdT2 = np.loadtxt(ST18 + "/wstate_dJdT.mtx")

def rebuild_flux(Uin, pin):
    """converged SIMPLE flux: phi = flux(rAU_rel*H_rel) - pflux (validated)."""
    offU_ = np.zeros((N, 3))
    np.add.at(offU_, owner, (pn_ + a_f)[:, None]*Uin[nei])
    np.add.at(offU_, nei, (-po + a_f)[:, None]*Uin[owner])
    pf = weights*pin[owner] + (1-weights)*pin[nei]
    gp = np.zeros((N, 3))
    tmp = Sf*pf[:, None]
    np.add.at(gp, owner, tmp); np.add.at(gp, nei, -tmp)
    np.add.at(gp, bc_cell, bSf*p_bvals[:, None])
    gp /= V[:, None]
    se = -gp
    Uf = weights[:, None]*Uin[owner] + (1-weights)[:, None]*Uin[nei]
    gU_ = np.zeros((N, 3, 3))
    tU = Sf[:, None]*Uf[:, :, None]
    np.add.at(gU_, owner, tU); np.add.at(gU_, nei, -tU)
    np.add.at(gU_, bc_cell, bSf[:, None]*U_bvals[:, :, None])
    gU_ /= V[:, None, None]
    d2 = np.swapaxes(gU_, 1, 2).copy()
    tr_ = np.trace(gU_, axis1=1, axis2=2)
    d2[:, 0, 0] -= 2.0/3.0*tr_; d2[:, 1, 1] -= 2.0/3.0*tr_; d2[:, 2, 2] -= 2.0/3.0*tr_
    ni = 0.5*(nuEff[owner] + nuEff[nei])
    Ti = weights[:, None, None]*(ni[:, None, None]*d2[owner]) \
       + (1-weights)[:, None, None]*(ni[:, None, None]*d2[nei])
    Tbnd = nuEff[bc_cell][:, None, None]*d2[bc_cell]
    dd = np.zeros((N, 3))
    tD = np.einsum("fi,fij->fj", Sf, Ti)
    np.add.at(dd, owner, tD); np.add.at(dd, nei, -tD)
    np.add.at(dd, bc_cell, np.einsum("fi,fij->fj", bSf, Tbnd))
    dd /= V[:, None]
    se += -dd
    sb = np.zeros((N, 3))
    np.add.at(sb, bc_cell[~ob], -phi_bvals[~ob, None]*U_bvals[~ob])
    np.add.at(sb, bc_cell[U_fixed], a_b[U_fixed, None]*U_bvals[U_fixed])
    sv = se*V[:, None] + sb
    sv += (1.0-ALPHAREL)*D_rel[:, None]*Uin
    Hh = (sv - offU_)/V[:, None]
    HbyA_ = mob[:, None]*Hh
    ph_int = np.einsum("fi,fi->f", Sf, weights[:, None]*HbyA_[owner]
                       + (1-weights)[:, None]*HbyA_[nei])
    ph_int -= kf_int*(pin[nei] - pin[owner])
    phb = np.zeros(nBndF)
    am = ~U_fixed
    phb[am] = np.einsum("fi,fi->f", bSf[am], HbyA_[bc_cell[am]])
    # fixed-p boundary (outlet): phi = phiHbyA_b - kf_b*(p_b - p_c), p_b=0
    phb[am] += np.where(kf_b[am] > 0, kf_b[am]*pin[bc_cell[am]], 0.0)
    return ph_int, phb

phi_rec, phib_rec = rebuild_flux(U, p)
log("\n=== part 12: flux rebuild validation (main state) ===")
log("phi internal: |rec - written|max=%.4e relL2=%.4e ; boundary outlet: "
    "rec sum=%.6e written sum=%.6e"
    % (np.abs(phi_rec - phi_int).max(),
       np.linalg.norm(phi_rec - phi_int)/np.linalg.norm(phi_int),
       phib_rec[U_assign].sum(), phi_bvals[U_assign].sum()))

# round-2 flux and thermal adjoint
phi2, phib2 = rebuild_flux(U2, p2)
# phiThermal ~ phi on cold faces (hot additions live in the hot region only)
phTh2 = phi2.copy(); phTh2b = phib2.copy()
# rebuild A_T at round-2 and solve A_T^T Tb = dJdT2
import scipy.sparse as sp
import scipy.sparse.linalg as spl
poT2 = np.where(phTh2 >= 0.0, phTh2, 0.0)
pnT2 = np.where(phTh2 >= 0.0, 0.0, phTh2)
diagT2 = np.zeros(N)
np.add.at(diagT2, owner, poT2)
np.add.at(diagT2, nei, -pnT2)
divp2 = np.zeros(N)
np.add.at(divp2, owner, phTh2); np.add.at(divp2, nei, -phTh2)
np.add.at(divp2, bc_cell, phTh2b)
diagT2 += -divp2
Aoff_o2 = pnT2 + kfT; Aoff_n2 = -poT2 + kfT
np.add.at(diagT2, owner, kfT); np.add.at(diagT2, nei, kfT)
obT2 = phTh2b >= 0.0
np.add.at(diagT2, bc_cell[obT2], phTh2b[obT2])
np.add.at(diagT2, bc_cell[T_fixed], kfTb[T_fixed])
rows2 = np.concatenate([owner, nei, np.arange(N)])
cols2 = np.concatenate([nei, owner, np.arange(N)])
vals2 = np.concatenate([Aoff_o2, Aoff_n2, diagT2])
AT2 = sp.csr_matrix((vals2, (rows2, cols2)), shape=(N, N))
Tb2 = spl.spsolve(AT2.tocsc().T, dJdT2)
log("round-2 Tb: |L2|=%.4e (written main Tb L2=%.4e); A_T2^T Tb2-dJdT2 "
    "maxabs=%.3e"
    % (np.linalg.norm(Tb2), np.linalg.norm(Tb), np.abs(AT2.T@Tb2-dJdT2).max()))
# round-2 g
g2 = np.zeros(nIF + nBndF)
jump2 = T2[nei] - T2[owner]
adjD2 = np.where(phTh2 >= 0.0, Tb2[nei], Tb2[owner])
g2[:nIF] = -coldmask_int*adjD2*jump2
M2 = phib2[patch_start["outlet"]:patch_start["outlet"]+84].sum()
Toc2 = T2[bc_cell[patch_start["outlet"]:patch_start["outlet"]+84]]
Tmix2 = float(np.sum(phib2[patch_start["outlet"]:patch_start["outlet"]+84]
                     * Toc2)/M2)
dJdp2 = -(Toc2 - Tmix2)/(M2*Tref)
g2[nIF+patch_start["outlet"]:nIF+patch_start["outlet"]+84] += dJdp2
log("round-2 g: |int|max=%.4e outlet dJ/dphi max=%.4e Tmix=%.4f M=%.6e "
    "(main: 107.47 114.36 687.84 4.073e-3)"
    % (np.abs(g2[:nIF]).max(), np.abs(dJdp2).max(), Tmix2, M2))
# route with V1 (round-2 rAU_u: alpha same -> rAU_u same; mobility same)
def route_V1_g(gv):
    b = np.zeros(NUNK)
    hA = np.zeros((N, 3)); gU_ = np.zeros((N, 3)); gP = np.zeros(N)
    ft = Sf*gv[:nIF][:, None]
    np.add.at(hA, owner, (ALPHAREL*rAU_u[owner]*weights)[:, None]*ft)
    np.add.at(hA, nei, (ALPHAREL*rAU_u[nei]*(1-weights))[:, None]*ft)
    np.add.at(gU_, owner, ((1-ALPHAREL)*weights)[:, None]*ft)
    np.add.at(gU_, nei, ((1-ALPHAREL)*(1-weights))[:, None]*ft)
    np.add.at(gP, owner, kf_int*gv[:nIF])
    np.add.at(gP, nei, -kf_int*gv[:nIF])
    am = ~U_fixed
    hA[bc_cell[am]] += (ALPHAREL*rAU_u[bc_cell[am]])[:, None]*(
        bSf[am]*gv[nIF:][am][:, None])
    gU_[bc_cell[am]] += (1-ALPHAREL)*(bSf[am]*gv[nIF:][am][:, None])
    gP[bc_cell[am]] += kf_b[am]*gv[nIF:][am]
    outU = np.zeros((N, 3))
    Vinvo = 1.0/V[owner]; Vinn = 1.0/V[nei]
    for cix in range(3):
        hAc = hA[:, cix]
        np.add.at(outU[:, cix], nei, -(upper*Vinvo)*hAc[owner])
        np.add.at(outU[:, cix], owner, -(lower*Vinn)*hAc[nei])
    b[0:NV:3] = gU_[:, 0]+outU[:, 0]; b[1:NV:3] = gU_[:, 1]+outU[:, 1]
    b[2:NV:3] = gU_[:, 2]+outU[:, 2]; b[NV:] = gP
    return b

bV1_2 = route_V1_g(g2)
m18 = np.loadtxt("/home/ys/dsH/b18_fix/b18rhs_thermalCoupling.mtx")
bexp18 = np.zeros(NUNK)
bexp18[0:NV:3] = m18[:, 0]; bexp18[1:NV:3] = m18[:, 1]
bexp18[2:NV:3] = m18[:, 2]; bexp18[NV:] = m18[:, 3]
d18 = bV1_2 - bexp18
log("round-2 V1 reconstruction vs exported b18 rhs: relL2=%.4e maxabs=%.4e"
    % (np.linalg.norm(d18)/np.linalg.norm(bexp18), np.abs(d18).max()))
for k, d in enumerate(["D1", "D2", "D3"]):
    w = wtrue({"D1": "D1_h0.001", "D2": "D2_h0.001", "D3": "D3_h0.001"}[d], 1e-3)
    log("  %s: b18export^T w=%.8f  myV1(round2)^T w=%.8f  labV1(main)^T w="
        "%.8f" % (d, bexp18 @ w, bV1_2 @ w,
                  {"D1": 0.05305, "D2": -0.00113, "D3": 0.03889}[d]))
log("PART12_DONE")

# ---------------- part 13: implementation fidelity: exported main-state b_TC
m18m = np.loadtxt("/home/ys/dsH/b18_main/b18rhs_thermalCoupling.mtx")
bexp18m = np.zeros(NUNK)
bexp18m[0:NV:3] = m18m[:, 0]; bexp18m[1:NV:3] = m18m[:, 1]
bexp18m[2:NV:3] = m18m[:, 2]; bexp18m[NV:] = m18m[:, 3]
bV1_main = route_V1(True) if False else None
# rebuild V1 (correct kf sign: s=+1 on owner) at MAIN state
bV1m = route_V1_g(g)  # g is main-state, route_V1_g uses correct kf sign
dm = bV1m - bexp18m
log("\n=== part 13: main-state V1 vs exported single-pass b18 rhs ===")
log("relL2=%.4e maxabs=%.4e  (|export|L2=%.4e |mine|L2=%.4e)"
    % (np.linalg.norm(dm)/np.linalg.norm(bexp18m), np.abs(dm).max(),
       np.linalg.norm(bexp18m), np.linalg.norm(bV1m)))
# per-block comparison
for nm, sl in [("U", slice(0, NV)), ("P", slice(NV, NUNK))]:
    a = bV1m[sl]; b_ = bexp18m[sl]
    log("  %s-block: relL2=%.4e  |mine|=%.4e |export|=%.4e"
        % (nm, np.linalg.norm(a-b_)/max(np.linalg.norm(b_), 1e-300),
           np.linalg.norm(a), np.linalg.norm(b_)))
for k, d in enumerate(["D1", "D2", "D3"]):
    w = wtrue({"D1": "D1_h0.001", "D2": "D2_h0.001", "D3": "D3_h0.001"}[d], 1e-3)
    log("  %s: export^T w=%.8f  mine^T w=%.8f"
        % (d, bexp18m @ w, bV1m @ w))
log("PART13_DONE")

# ---------------- part 14: new export vs OLD export (same main state) ------
d_old_new = bexp18m - bexp
log("\n=== part 14: b18 main export vs OLD b15 export (same state) ===")
log("relL2=%.4e maxabs=%.4e (|old|=%.4e)"
    % (np.linalg.norm(d_old_new)/np.linalg.norm(bexp),
       np.abs(d_old_new).max(), np.linalg.norm(bexp)))
log("OLD^T w: D1=%.6f D2=%.6f D3=%.6f"
    % (bexp @ wtrue("D1_h0.001", 1e-3),
       bexp @ wtrue("D2_h0.001", 1e-3),
       bexp @ wtrue("D3_h0.001", 1e-3)))
log("PART14_DONE")

# ---------------- part 15: block regression of the exported rhs ------------
# my per-block fields at main state
def route_parts():
    hA = np.zeros((N, 3)); gU = np.zeros((N, 3)); gP = np.zeros(N)
    ft = Sf*g[:nIF][:, None]
    np.add.at(hA, owner, (ALPHAREL*rAU_u[owner]*weights)[:, None]*ft)
    np.add.at(hA, nei, (ALPHAREL*rAU_u[nei]*(1-weights))[:, None]*ft)
    np.add.at(gU, owner, ((1-ALPHAREL)*weights)[:, None]*ft)
    np.add.at(gU, nei, ((1-ALPHAREL)*(1-weights))[:, None]*ft)
    np.add.at(gP, owner, kf_int*g[:nIF])
    np.add.at(gP, nei, -kf_int*g[:nIF])
    am = ~U_fixed
    hA[bc_cell[am]] += (ALPHAREL*rAU_u[bc_cell[am]])[:, None]*(
        bSf[am]*g[nIF:][am][:, None])
    gU[bc_cell[am]] += (1-ALPHAREL)*(bSf[am]*g[nIF:][am][:, None])
    gP[bc_cell[am]] += kf_b[am]*g[nIF:][am]
    hU = np.zeros((N, 3))
    Vinvo = 1.0/V[owner]; Vinn = 1.0/V[nei]
    for cix in range(3):
        hAc = hA[:, cix]
        np.add.at(hU[:, cix], nei, -(upper*Vinvo)*hAc[owner])
        np.add.at(hU[:, cix], owner, -(lower*Vinn)*hAc[nei])
    return hU, gU, gP

hU_m, gU_m, gP_m = route_parts()
# regression on U block: export_U = a*hU + b*gU (+ resid)
A_U = np.stack([hU_m.reshape(-1), gU_m.reshape(-1)], axis=1)
y_U = bexp18m[:NV]
coef_U, res_U, _, _ = np.linalg.lstsq(A_U, y_U, rcond=None)
pred_U = A_U @ coef_U
log("\n=== part 15: block regression vs export ===")
log("U-block: coef hA=%.6f direct=%.6f  residual relL2=%.4e"
    % (coef_U[0], coef_U[1],
       np.linalg.norm(y_U-pred_U)/np.linalg.norm(y_U)))
A_P = gP_m.reshape(-1, 1)
y_P = bexp18m[NV:]
coef_P, _, _, _ = np.linalg.lstsq(A_P, y_P, rcond=None)
log("P-block: coef kf=%.6f residual relL2=%.4e"
    % (coef_P[0], np.linalg.norm(y_P-A_P@coef_P)/np.linalg.norm(y_P)))
# also regress hA-only with alphaRel absorbed: export_U vs
# alpha*prodRAU*rAU... try variant: hA with rAU= mob (relaxed) instead
hA2 = np.zeros((N, 3))
ft = Sf*g[:nIF][:, None]
np.add.at(hA2, owner, (mob[owner]*weights)[:, None]*ft)
np.add.at(hA2, nei, (mob[nei]*(1-weights))[:, None]*ft)
am = ~U_fixed
hA2[bc_cell[am]] += mob[bc_cell[am]][:, None]*(bSf[am]*g[nIF:][am][:, None])
Vinvo = 1.0/V[owner]; Vinn = 1.0/V[nei]
hU2 = np.zeros((N, 3))
for cix in range(3):
    hAc = hA2[:, cix]
    np.add.at(hU2[:, cix], nei, -(upper*Vinvo)*hAc[owner])
    np.add.at(hU2[:, cix], owner, -(lower*Vinn)*hAc[nei])
A2 = np.stack([hU2.reshape(-1), gU_m.reshape(-1)], axis=1)
coef2, _, _, _ = np.linalg.lstsq(A2, y_U, rcond=None)
log("U-block (rAU=mob variant): coef hA=%.6f direct=%.6f residual=%.4e"
    % (coef2[0], coef2[1],
       np.linalg.norm(y_U-A2@coef2)/np.linalg.norm(y_U)))
log("PART15_DONE")
