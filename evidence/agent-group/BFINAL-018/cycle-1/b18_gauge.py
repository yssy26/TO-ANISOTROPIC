#!/usr/bin/env python3
"""
BFINAL-018 fix-3 gauge: b_TC^T w_true vs (FD_J - direct-thermal part).

Uses ONLY same-pass b18_fix artifacts:
  - b18rhs_thermalCoupling.mtx   (production-path physical rhs, fixed folding)
  - stageB2/wstate_*             (fresh FD states, same run)
  - 1/{Tb,T,dDTDxh,alpha,...}    (same-run fields)

Identity decomposition (BFINAL-018 derivation; see DERIVATION_FIX3.md):
  FD_J = A + B + C  with
    A = J_phi^T phi'_out        (outlet-flux explicit part) = FD - thermal_g
    B = -Tb^T (dR_T/dphi) phi'  (T response to flux change)  = thermal_g - C
    C = -Tb^T R_T,x             (direct thermal, dDTDxh)     [production term]
  exact T-elimination source:  b_TC = M^T g  with g = J_phi - (dR_T/dphi)^T Tb
  gauge identity:  b_TC^T w_true == FD_J - C - Gx
    Gx = g^T phi_x  (flux direct-alpha part; NOT representable as (U,p) source)
The B16-preregistered comparator FD-thermal_g (=A) matches only an
outlet-only folding and is reported for reference.
"""
import sys, json
import numpy as np

CASE = sys.argv[1] if len(sys.argv) > 1 else "/home/ys/dsH/b18_fix"
S = CASE + "/1"
ST = CASE + "/stageB2"
MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
N = 33600; NV = 3*N; NUNK = 4*N
ALPHAREL = 0.4
FD_J = {"D1": 0.0428559144263, "D2": -0.00346892642156, "D3": -0.000282658160322}

_BND_CACHE = {}

def read_field(path):
    """robust parser returning (internal, uniformVal); boundary blocks
    parsed once into _BND_CACHE[path] = {patch: (type, values|None)}"""
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
        internal = None
        val = np.array(toks) if len(toks) > 1 else float(toks[0])
    else:
        vals = []
        for line in txt[j+1:k].splitlines():
            for tok in line.replace(";", " ").split():
                for tt in tok.replace("(", " ").replace(")", " ").split():
                    vals.append(float(tt))
        internal = np.array(vals)
        val = None
    # boundary blocks
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
                    bval = np.array([float(x) for x in u.strip("()").split()])
                else:
                    bval = float(u)
            if t.startswith("value") and depth == 1:
                if " nonuniform" in t or t.rstrip(";").endswith("nonuniform"):
                    bval = ("list", q)
                elif "uniform" in t:
                    u = t[t.index("uniform")+len("uniform"):].rstrip(";").strip()
                    if u.startswith("("):
                        bval = np.array([float(x) for x in u.strip("()").split()])
                    else:
                        bval = float(u)
                else:
                    bval = ("list", q)
            q += 1
        if isinstance(bval, tuple) and bval[0] == "list":
            q0 = bval[1]; qq = q0
            while "(" not in lines[qq]: qq += 1
            kk = qq
            while ")" not in lines[kk]: kk += 1
            vals2 = []
            for line in lines[qq:kk+1]:
                for tok in line.strip("()").replace("(", " ").replace(")", " ").split():
                    try: vals2.append(float(tok))
                    except ValueError: pass
            bnd[name] = (btype, np.array(vals2))
        elif bval is not None:
            bnd[name] = (btype, bval)
        else:
            bnd[name] = (btype, None)
        p = q
    _BND_CACHE[path] = bnd
    return internal, val

# mesh (validated in the lab: Sf 1.4e-19, owner exact)
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
V = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")[:, 0]
nuEff = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")[:, 1]
bnd_tbl = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
bc_cell = bnd_tbl[:, 0].astype(np.int64)
bSf = bnd_tbl[:, 1:4]; bmag = bnd_tbl[:, 4]; bdel = bnd_tbl[:, 5]
nBndF = len(bnd_tbl)
weights = np.loadtxt(MESH + "/b16mesh_weights.mtx")
patch_names = ["inlet", "outlet", "hotInlet", "hotOutlet",
               "solidEndWalls", "bottomWall", "topWall", "sideWalls"]
nFacesPatch = {"inlet": 84, "outlet": 84, "hotInlet": 196, "hotOutlet": 196,
               "solidEndWalls": 280, "bottomWall": 1120, "topWall": 1120,
               "sideWalls": 4800}
patch_start = {}
o = 0
for pn in patch_names:
    patch_start[pn] = o; o += nFacesPatch[pn]

# fields
U = read_field(S + "/U")[0].reshape(N, 3)
p = read_field(S + "/p")[0]
alpha = read_field(S + "/alpha")[0]
Tb = read_field(S + "/Tb")[0]
T = read_field(S + "/T")[0]
phi_int = read_field(S + "/phi")[0]
phTh_int = read_field(S + "/phiThermal")[0]
dAlphaDxh = read_field(S + "/dAlphaDxh")[0]
dDTDxh = read_field(S + "/dDTDxh")[0].reshape(N, 6)
p_b = read_field(S + "/p")[0]  # placeholder
# boundary values needed: p_b (zeroGradient->cell), phi_b for conv diag
def read_bnd(path, name):
    if path not in _BND_CACHE:
        read_field(path)
    ent = _BND_CACHE[path].get(name)
    if ent is None:
        return None
    return ent[1]
phi_bv = np.zeros(nBndF)
p_bv = np.zeros(nBndF)
U_fixed = np.zeros(nBndF, dtype=bool)
for pn in patch_names:
    nFp = nFacesPatch[pn]
    vb = read_bnd(S + "/phi", pn)
    if isinstance(vb, float): phi_bv[patch_start[pn]:patch_start[pn]+nFp] = vb
    elif vb is not None: phi_bv[patch_start[pn]:patch_start[pn]+nFp] = vb
    vb = read_bnd(S + "/p", pn)
    if isinstance(vb, float): p_bv[patch_start[pn]:patch_start[pn]+nFp] = vb
    elif vb is not None: p_bv[patch_start[pn]:patch_start[pn]+nFp] = vb
    else: p_bv[patch_start[pn]:patch_start[pn]+nFp] = p[bc_cell[patch_start[pn]:patch_start[pn]+nFp]]
    U_fixed[patch_start[pn]:patch_start[pn]+nFp] = (pn in ("inlet",))
phTh_b = np.zeros(nBndF)
for pn in patch_names:
    nFp = nFacesPatch[pn]
    vb = read_bnd(S + "/phiThermal", pn)
    if isinstance(vb, float): phTh_b[patch_start[pn]:patch_start[pn]+nFp] = vb
    elif vb is not None: phTh_b[patch_start[pn]:patch_start[pn]+nFp] = vb
T_b = read_field(S + "/T")[0]

# ---- g face functional (B0.2 formulas) ----
jumpT = T[nei] - T[owner]
cfm = np.zeros(nIF)
vb = read_bnd(S + "/coldFaceMask", "inlet")  # not needed for internal
cfm_i = read_field(S + "/coldFaceMask")[0]
adjDown = np.where(phTh_int >= 0.0, Tb[nei], Tb[owner])
g = np.zeros(nIF + nBndF)
g[:nIF] = -cfm_i*adjDown*jumpT
Mflow, Tref = np.loadtxt(ST + "/wstate_outletMeta.mtx")
ophi = phi_bv[patch_start["outlet"]:patch_start["outlet"]+84]
Toc = T[bc_cell[patch_start["outlet"]:patch_start["outlet"]+84]]
Tmix = float(np.sum(ophi*Toc)/Mflow)
dJdphi_out = -(Toc - Tmix)/(Mflow*Tref)
g[nIF + patch_start["outlet"]:nIF + patch_start["outlet"]+84] += dJdphi_out

# ---- mobility / rAU (lab-validated rebuilds) ----
nuEff_f = weights*nuEff[owner] + (1.0-weights)*nuEff[nei]
a_f = nuEff_f*magSf*deltaCoeffs
up_own = phi_int >= 0.0
po = np.where(up_own, phi_int, 0.0)
pn_ = np.where(up_own, 0.0, phi_int)
diag_u = np.zeros(N)
np.add.at(diag_u, owner, po - a_f)
np.add.at(diag_u, nei, -pn_ - a_f)
diag_u += alpha*V
divphi = np.zeros(N)
np.add.at(divphi, owner, phi_int)
np.add.at(divphi, nei, -phi_int)
np.add.at(divphi, bc_cell, phi_bv)
diag_u += -divphi
ob = phi_bv >= 0.0
np.add.at(diag_u, bc_cell, np.where(ob, phi_bv, 0.0))
a_b = nuEff[bc_cell]*bmag*bdel
np.add.at(diag_u, bc_cell[U_fixed], a_b[U_fixed])
sumOff = np.zeros(N)
np.add.at(sumOff, owner, np.abs(pn_ + a_f))
np.add.at(sumOff, nei, np.abs(-po + a_f))
np.add.at(sumOff, bc_cell, np.abs(np.where(ob, phi_bv, 0.0)))
D_rel = np.maximum(np.abs(diag_u), sumOff)/ALPHAREL
mob = 1.0/(D_rel/V)
rAU_u = 1.0/(diag_u/V)
upper = pn_ + a_f
lower = -po + a_f
kf_int = (weights*mob[owner] + (1.0-weights)*mob[nei])*magSf*deltaCoeffs
kf_b = mob[bc_cell]*bmag*bdel

# ---- HbyA rebuild for Gx (rxPressureRowT basis) ----
w_ = weights
offU = np.zeros((N, 3))
np.add.at(offU, owner, (pn_ + a_f)[:, None]*U[nei])
np.add.at(offU, nei, (-po + a_f)[:, None]*U[owner])
pf_int = w_*p[owner] + (1-w_)*p[nei]
gradp = np.zeros((N, 3))
tmp = Sf*pf_int[:, None]
np.add.at(gradp, owner, tmp); np.add.at(gradp, nei, -tmp)
np.add.at(gradp, bc_cell, bSf*p_bv[:, None])
gradp /= V[:, None]
src_exp = -gradp
Uf_int = w_[:, None]*U[owner] + (1-w_)[:, None]*U[nei]
U_bvals = np.zeros((nBndF, 3))
for pn in patch_names:
    nFp = nFacesPatch[pn]
    sl = slice(patch_start[pn], patch_start[pn]+nFp)
    vb = read_bnd(S + "/U", pn)
    if isinstance(vb, np.ndarray): U_bvals[sl] = vb
    else: U_bvals[sl] = U[bc_cell[sl]]
gU = np.zeros((N, 3, 3))
tmpU = Sf[:, None]*Uf_int[:, :, None]
np.add.at(gU, owner, tmpU); np.add.at(gU, nei, -tmpU)
np.add.at(gU, bc_cell, bSf[:, None]*U_bvals[:, :, None])
gU /= V[:, None, None]
dev2T = np.swapaxes(gU, 1, 2).copy()
tr = np.trace(gU, axis1=1, axis2=2)
dev2T[:, 0, 0] -= 2.0/3.0*tr; dev2T[:, 1, 1] -= 2.0/3.0*tr; dev2T[:, 2, 2] -= 2.0/3.0*tr
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
np.add.at(src_bcomp, bc_cell[ib], -phi_bv[ib, None]*U_bvals[ib])
fb = U_fixed
np.add.at(src_bcomp, bc_cell[fb], a_b[fb, None]*U_bvals[fb])
src_vec = src_exp*V[:, None] + src_bcomp
src_vec += (1.0-ALPHAREL)*D_rel[:, None]*U
H = (src_vec - offU)/V[:, None]
HbyA = mob[:, None]*H
drAU = -mob*mob/ALPHAREL
dHbyA = (mob/ALPHAREL)[:, None]*((1.0-ALPHAREL)*U - HbyA)
g0_int = magSf*deltaCoeffs*(p[nei] - p[owner])
U_assign = ~U_fixed

z3 = np.loadtxt("/home/ys/dsH/b15_export/stageB6_rxc_z_analytic.mtx").reshape(3, N)

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

ef = d_vec/np.maximum(np.linalg.norm(d_vec, axis=1), 1e-30)[:, None]
def dir_deriv(cells):
    D6 = dDTDxh[cells]
    ex, ey, ez = ef[:, 0], ef[:, 1], ef[:, 2]
    return (ex*ex*D6[:, 0] + ex*ey*D6[:, 1] + ex*ez*D6[:, 2]
            + ey*ex*D6[:, 1] + ey*ey*D6[:, 3] + ey*ez*D6[:, 4]
            + ez*ex*D6[:, 2] + ez*ey*D6[:, 4] + ez*ez*D6[:, 5])
ddo = dir_deriv(owner); ddn = dir_deriv(nei)
base = (Tb[owner]-Tb[nei])*(T[owner]-T[nei])*deltaCoeffs*magSf
Cfield = np.zeros(N)
np.add.at(Cfield, owner, -weights*ddo*base)
np.add.at(Cfield, nei, -(1.0-weights)*ddn*base)
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
T_bv = np.zeros(nBndF)
for pn in patch_names:
    nFp = nFacesPatch[pn]
    sl = slice(patch_start[pn], patch_start[pn]+nFp)
    vb = read_bnd(S + "/T", pn)
    if isinstance(vb, float): T_bv[sl] = vb
    elif vb is not None: T_bv[sl] = vb
    else: T_bv[sl] = T[bc_cell[sl]]
snb = (T_bv - T[bc_cell])*bdel
np.add.at(Cfield, bc_cell, bdir_deriv()*Tb[bc_cell]*snb*bmag)

# ---- load the exported b_TC (fixed folding) ----
bTC = np.zeros(NUNK)
m = np.loadtxt(CASE + "/b18rhs_thermalCoupling.mtx")
bTC[0:NV:3] = m[:, 0]; bTC[1:NV:3] = m[:, 1]; bTC[2:NV:3] = m[:, 2]
bTC[NV:] = m[:, 3]

dJdT = np.loadtxt(ST + "/wstate_dJdT.mtx")
res = {}
print("dir  h       b_TC^T w_true   FD_J        thermal_g   C_direct    Gx     "
      "| FD-C-Gx   err    | A=FD-tg   errA")
for d, h in [("D1", 1e-3), ("D2", 1e-3), ("D3", 1e-3),
             ("D1", 3e-4), ("D2", 3e-4), ("D3", 3e-4)]:
    tag = "%s_h%s" % (d, ("0.001" if h == 1e-3 else "0.0003"))
    Up = np.loadtxt(ST + "/wstate_%s_p_U.mtx" % tag)
    Um = np.loadtxt(ST + "/wstate_%s_m_U.mtx" % tag)
    pp = np.loadtxt(ST + "/wstate_%s_p_p.mtx" % tag)
    pm_ = np.loadtxt(ST + "/wstate_%s_m_p.mtx" % tag)
    w = np.zeros(NUNK)
    w[0:NV:3] = ((Up-Um)/(2*h))[:, 0]
    w[1:NV:3] = ((Up-Um)/(2*h))[:, 1]
    w[2:NV:3] = ((Up-Um)/(2*h))[:, 2]
    w[NV:] = (pp-pm_)/(2*h)
    Tp = np.loadtxt(ST + "/wstate_%s_p_T.mtx" % tag)
    Tm = np.loadtxt(ST + "/wstate_%s_m_T.mtx" % tag)
    tg = float(np.dot(dJdT, (Tp-Tm)/(2*h)))
    k = ["D1", "D2", "D3"].index(d)
    Ck = float(np.dot(Cfield, z3[k]))
    Gxk = Gx_of(z3[k])
    fd = FD_J[d]
    A = fd - tg
    tgt = fd - Ck - Gxk
    val = float(bTC @ w)
    err = abs(val - tgt)/abs(tgt)
    errA = abs(val - A)/abs(A)
    print("%s h=%g  %+.8f  %+.8f  %+.8f  %+.8f  %+.6f | %+.8f  %.3f | "
          "%+.8f  %.3f"
          % (d, h, val, fd, tg, Ck, Gxk, tgt, err, A, errA))
    res["%s_h%g" % (d, h)] = dict(val=val, fd=fd, thermal=tg, C=Ck, Gx=Gxk,
                                  target=tgt, A=A, err=err, errA=errA)
json.dump(res, open("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/"
                    "BFINAL-018/cycle-1/b18_gauge.json", "w"), indent=1)
print("GAUGE_DONE")
