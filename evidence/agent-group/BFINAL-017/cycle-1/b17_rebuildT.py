#!/usr/bin/env python3
"""
BFINAL-017 phase-0 AD2-deep: independent offline rebuild of the production
rxPressureRowT (J_P^T pc) pipeline:

  UEqn = fvm::div(phi,U) [bounded Gauss upwind]
       - fvm::laplacian(nuEffFrozen,U) [Gauss linear corrected]
       + fvm::Sp(alpha,U)
      == -fvc::grad(p) + fvOptions(U)  [- fvc::div(nuEff*dev2(T(grad U))) if RANS]
  relax(0.4) -> A()=D/V, H()=(S - offdiag*U)/V, rAU=1/A(), HbyA=rAU*H()
  drAU  = -rAU^2/alphaRel
  dHbyA = (rAU/alphaRel)*((1-alphaRel)*U - HbyA)
  g0    = flux(fvm::laplacian(one,p))
  T1[own]+=wf*(Sf&dHbyA[own])*(pc[own]-pc[nei]); T1[nei]+=(1-wf)*...
  T1[b]  += (Sfb&dHbyA[c])*pc[c]   (assignable-U patches only)
  T2[own]-=wf*g0*(pc_o-pc_n)*drAU[own]; ...

Validation ladder (no production code on the oracle side):
  V1: T(pc15)  == -g15_pr/dAlphaDxh        (b15 raw pressure-row export)
  V2: J_P*w_k  == anRPd_k (rxc export) and J_P*sin == anRPaW (rxa export)
Adjudication:
  T(pc16) vs the b16-implied production transpose field (-gpr16/dAlphaDxh)
"""
import time
import numpy as np
import scipy.sparse as sp

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
E = "/home/ys/dsH/b15_export"
S = "/home/ys/dsH/b16_states/1"
N = 33600; NV = 3*N; NUNK = 4*N
ALPHAREL = 0.4

t0 = time.time()
def log(m): print(m, flush=True)

# ---------------- field parsing (internal + boundary values) ----------------
def read_field(path):
    """returns (internalValues, {patchName: np.array of boundary values})"""
    txt = open(path).read()
    i = txt.index("internalField")
    j = txt.index("\n(", i)
    k = txt.index("\n)", j)
    if txt[i:j].split()[1] == "uniform":
        # uniform <val>;
        toks = txt[i:j].split()
        # e.g. 'internalField uniform 5.19e-05 ;' possibly with (x y z)
        seg = txt[i:j]
        u = seg[seg.index("uniform")+len("uniform"):].strip().rstrip(";").strip()
        if u.startswith("("):
            val = np.array([float(t) for t in u.strip("()").split()])
        else:
            val = float(u)
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
            p += 1
            continue
        # block header may be "name {" on one line or "name" + "{" on next
        if s.endswith("{"):
            name = s[:-1].strip()
            hdr = p
        elif (p+1 < len(lines) and lines[p+1].strip() == "{"
              and " " not in s and not s.startswith(("type", "value", "internalField"))):
            name = s
            hdr = p+1
        else:
            p += 1
            continue
        if name == "boundaryField":
            p += 1
            continue
        # scan block for 'value' or 'uniform' (block body starts after hdr)
        depth = 1; q = hdr+1; bval = None; btype = None
        while q < len(lines) and depth > 0:
            t = lines[q].strip()
            if t.endswith("{") and t != "{": depth += 1
            if t == "{": depth += 1
            if t == "}": depth -= 1
            if t.startswith("type"):
                btype = t.split()[1].rstrip(";")
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
            q0 = bval[1]
            qq = q0
            while "(" not in lines[qq]: qq += 1
            kk = qq
            while ")" not in lines[kk]: kk += 1
            vals = []
            for line in lines[qq:kk+1]:
                for tok in line.strip().strip("()").replace("(", " ").replace(")", " ").split():
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
V = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")[:, 0]
# boundary table (patch order from polyMesh/boundary)
bnd_tbl = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
bc_cell = bnd_tbl[:, 0].astype(np.int64)
bSf = bnd_tbl[:, 1:4]; bmag = bnd_tbl[:, 4]; bdel = bnd_tbl[:, 5]
nBndF = len(bnd_tbl)
weights = np.loadtxt(MESH + "/b16mesh_weights.mtx")
log("mesh ready (t=%.1fs)" % (time.time()-t0))

# ---------------- fields ----------------
U_i, _, U_b = read_field(S + "/U")
U = U_i.reshape(N, 3)
p_i, _, p_b = read_field(S + "/p"); p = p_i
phi_i, _, phi_b = read_field(S + "/phi"); phi_int = phi_i
alpha_i, _, _ = read_field(S + "/alpha"); alpha = alpha_i
pc16_i, _, _ = read_field(S + "/pc"); pc16 = pc16_i
Uc16_i, _, _ = read_field(S + "/Uc"); Uc16 = Uc16_i.reshape(N, 3)
daD = read_field(S + "/dAlphaDxh")[0]
NU = 5.19009e-05

patch_names = ["inlet", "outlet", "hotInlet", "hotOutlet",
               "solidEndWalls", "bottomWall", "topWall", "sideWalls"]
nFacesPatch = {"inlet": 84, "outlet": 84, "hotInlet": 196, "hotOutlet": 196,
               "solidEndWalls": 280, "bottomWall": 1120, "topWall": 1120,
               "sideWalls": 4800}
phi_bvals = np.zeros(nBndF); p_bvals = np.zeros(nBndF); U_bvals = np.zeros((nBndF, 3))
U_fixed = np.zeros(nBndF, dtype=bool)   # fixesValue (fixedValue/noSlip)
U_assign = np.zeros(nBndF, dtype=bool)  # assignable (zeroGradient etc.)
o = 0
for pn in patch_names:
    nFp = nFacesPatch[pn]
    btypeU, bvu = U_b[pn]; btypeP, bvp = p_b[pn]
    btypePhi, bvphi = phi_b[pn]
    if isinstance(bvphi, float):
        phi_bvals[o:o+nFp] = bvphi
    elif bvphi is not None:
        phi_bvals[o:o+nFp] = bvphi
    if isinstance(bvp, float):
        p_bvals[o:o+nFp] = bvp
    elif bvp is not None:
        p_bvals[o:o+nFp] = bvp
    else:
        p_bvals[o:o+nFp] = p[bc_cell[o:o+nFp]]  # zeroGradient
    if isinstance(bvu, np.ndarray):
        U_bvals[o:o+nFp] = bvu
        U_fixed[o:o+nFp] = True
    elif bvu is not None:
        U_bvals[o:o+nFp] = bvu
        U_fixed[o:o+nFp] = True
    else:
        U_bvals[o:o+nFp] = U[bc_cell[o:o+nFp]]
        U_assign[o:o+nFp] = btypeU in ("zeroGradient", "calculated")
    o += nFp
log("fields ready (t=%.1fs)" % (time.time()-t0))

# ---------------- UEqn rebuild ----------------
def build_ueqn(bounded_corr=True):
    # integrated diagonal / source (scalar, same for each component)
    diag = np.zeros(N); offp = np.zeros(N); offn = np.zeros(N); src = np.zeros(N)
    # convection: upwind implicit
    #   internal: coefficient on U_up: +phi into owner row, -phi into neighbour row
    up_own = phi_int >= 0.0
    # contributions: A[own,own]+=phi; A[nei,own]-=phi (phi>=0)
    #                A[own,nei]+=phi; A[nei,nei]-=phi (phi<0)
    po = np.where(up_own, phi_int, 0.0)
    pn_ = np.where(up_own, 0.0, phi_int)
    np.add.at(diag, owner, po)
    np.add.at(offp, owner, pn_)          # A[own,nei]
    np.add.at(offp, nei, -po)            # A[nei,own]
    np.add.at(diag, nei, -pn_)
    # boundary convection: outflow -> diag; inflow -> source (-phi*U_b)
    ob = phi_bvals >= 0.0
    np.add.at(diag, bc_cell[ob], phi_bvals[ob])
    np.add.at(src, bc_cell[~ob], -phi_bvals[~ob]*1.0)  # per component: -phi*Ub_c
    src_bcomp = np.zeros((N, 3))
    np.add.at(src_bcomp, bc_cell[~ob],
              -phi_bvals[~ob, None]*U_bvals[~ob])
    if bounded_corr:
        # fvm::Sp(-fvc::div(phi), U): diag += -divphi*V
        divphi = np.zeros(N)
        np.add.at(divphi, owner, phi_int)
        np.add.at(divphi, nei, -phi_int)
        np.add.at(divphi, bc_cell, phi_bvals)
        diag += -divphi
    # laplacian(nuEffFrozen, U): nu uniform
    a_f = NU*magSf*deltaCoeffs
    np.add.at(diag, owner, -a_f)
    np.add.at(diag, nei, -a_f)
    np.add.at(offp, owner, a_f)
    np.add.at(offp, nei, a_f)
    # boundary laplacian: fixedValue U -> diag += a_b, src += a_b*U_b
    a_b = NU*bmag*bdel
    fb = U_fixed
    np.add.at(diag, bc_cell[fb], a_b[fb])
    np.add.at(src_bcomp, bc_cell[fb], a_b[fb, None]*U_bvals[fb])
    # Sp(alpha,U)
    diag += alpha*V
    # relax: D/alphaRel ; src += (1-alphaRel)*D_rel*U
    diag_rel = diag/ALPHAREL
    # RHS explicit: -grad(p)  (Gauss linear), and -div(nu*dev2(T(gradU))) if RANS
    # face p (linear interp): internal
    w = weights
    pf_int = w*p[owner] + (1.0-w)*p[nei]
    pf_bnd = p_bvals
    gradp = np.zeros((N, 3))
    tmp = Sf*pf_int[:, None]
    np.add.at(gradp, owner, tmp)
    np.add.at(gradp, nei, -tmp)
    np.add.at(gradp, bc_cell, bSf*pf_bnd[:, None])
    gradp /= V[:, None]
    # src += V*(-gradp)   (A U = b with b = -grad(p) ... sign fixed by validation)
    src_exp = -gradp
    # dev2 term: -fvc::div(nu*dev2(T(grad U)))
    # grad U (Gauss linear): internal face interp of U
    Uf_int = w[:, None]*U[owner] + (1.0-w)[:, None]*U[nei]
    Uf_bnd = U_bvals
    gU = np.zeros((N, 3, 3))
    tmpU = Sf[:, None]*Uf_int[:, :, None]     # Sf_j * U_k  (outer: (Sf x U))
    np.add.at(gU, owner, tmpU)
    np.add.at(gU, nei, -tmpU)
    np.add.at(gU, bc_cell, bSf[:, None]*Uf_bnd[:, :, None])
    gU /= V[:, None, None]
    dev2T = np.swapaxes(gU, 1, 2).copy()
    tr = np.trace(gU, axis1=1, axis2=2)
    dev2T[:, 0, 0] -= 2.0/3.0*tr; dev2T[:, 1, 1] -= 2.0/3.0*tr; dev2T[:, 2, 2] -= 2.0/3.0*tr
    tens = NU*dev2T
    Tf_int = w[:, None, None]*tens[owner] + (1.0-w)[:, None, None]*tens[nei]
    Tf_bnd = tens[bc_cell]  # linear boundary extrapolation for calculated BCs
    divdev = np.zeros((N, 3))
    tmpD = np.einsum("fi,fij->fj", Sf, Tf_int)
    np.add.at(divdev, owner, tmpD)
    np.add.at(divdev, nei, -tmpD)
    np.add.at(divdev, bc_cell, np.einsum("fi,fij->fj", bSf, Tf_bnd))
    divdev /= V[:, None]
    src_exp += -divdev
    src_vec = src_exp*V[:, None] + src_bcomp
    # HbyA per component: H = (src_vec + relax-adjust - offdiag*U)/V
    # relax source adjustment: src += (1-alphaRel)*diag_rel*U_current
    src_vec += (1.0-ALPHAREL)*diag_rel[:, None]*U
    A = diag_rel/V
    rAU = 1.0/A
    # offdiag*U:  for each cell: sum over faces: offp contribution from neighbours
    # (owner row gets A[own,nei]*U_nei; neighbour row gets A[nei,own]*U_own)
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, pn_[:, None]*U[nei] - po[:, None]*0.0)
    # careful: A[own,nei]=pn_ multiplies U[nei]; A[nei,own]=-po multiplies U[own]
    offU2 = np.zeros((N, 3))
    np.add.at(offU2, owner, pn_[:, None]*U[nei])
    np.add.at(offU2, nei, -po[:, None]*U[owner])
    # laplacian offdiag: A[own,nei]=a_f * U[nei]; A[nei,own]=a_f*U[own]
    np.add.at(offU2, owner, a_f[:, None]*U[nei])
    np.add.at(offU2, nei, a_f[:, None]*U[owner])
    H = (src_vec - offU2)/V[:, None]
    HbyA = rAU[:, None]*H
    return rAU, HbyA

rAU, HbyA = build_ueqn(bounded_corr=True)
log("UEqn rebuilt: rAU avg=%.6e min=%.6e max=%.6e (t=%.1fs)"
    % (rAU.mean(), rAU.min(), rAU.max(), time.time()-t0))

drAU = -rAU*rAU/ALPHAREL
dHbyA = (rAU/ALPHAREL)[:, None]*((1.0-ALPHAREL)*U - HbyA)

# g0 = flux(laplacian(one,p)); boundary gamma=1
g0_int = magSf*deltaCoeffs*(p[nei] - p[owner])
g0_bnd = bmag*bdel*(p_bvals - p[bc_cell])
log("g0: |bnd|max=%.6e (b15 log: 0.0650620425279)" % np.max(np.abs(g0_bnd)))

def build_T(pc):
    T1 = np.zeros(N); T2 = np.zeros(N)
    pcdiff = pc[owner] - pc[nei]
    SfdH_own = np.einsum("fi,fi->f", Sf, dHbyA[owner])
    SfdH_nei = np.einsum("fi,fi->f", Sf, dHbyA[nei])
    np.add.at(T1, owner, weights*SfdH_own*pcdiff)
    np.add.at(T1, nei, (1.0-weights)*SfdH_nei*pcdiff)
    np.add.at(T2, owner, -weights*g0_int*pcdiff*drAU[owner])
    np.add.at(T2, nei, -(1.0-weights)*g0_int*pcdiff*drAU[nei])
    # boundary T1 on assignable-U patches (outlet)
    am = U_assign
    np.add.at(T1, bc_cell[am],
              np.einsum("fi,fi->f", bSf[am], dHbyA[bc_cell[am]])*pc[bc_cell[am]])
    return T1 + T2, T1, T2

def JP_w(w):
    """assembleWeightedRPa replica: J_P*w = div(dphiHbyA_w - dflux_w)"""
    dHbyAw = dHbyA*w[:, None]
    drAUw = drAU*w
    dphi = np.einsum("fi,fi->f", Sf, weights[:, None]*dHbyAw[owner]
                     + (1.0-weights)[:, None]*dHbyAw[nei])
    # boundary: assignable patches use cell value; else 0
    dbnd = np.zeros(nBndF)
    am = U_assign
    dbnd[am] = np.einsum("fi,fi->f", bSf[am], dHbyAw[bc_cell[am]])
    # dflux_w: laplacian(drAU*w, p) flux; boundary gamma = 0 (w field zero at bnd)
    dfl = drAUw[owner]*0.0
    gf = (weights*drAUw[owner] + (1.0-weights)*drAUw[nei])
    dfl = gf*magSf*deltaCoeffs*(p[nei] - p[owner])
    out = np.zeros(N)
    np.add.at(out, owner, dphi - dfl)
    np.add.at(out, nei, -(dphi - dfl))
    np.add.at(out, bc_cell, dbnd)
    return out

# ---------------- validation ladder ----------------
lam15 = np.loadtxt(E + "/stageB6_lambda.mtx")
pc15 = lam15[NV:]; Uc15 = lam15[:NV].reshape(N, 3)
g15_pr = np.loadtxt(E + "/rxpr_prod_gsensh_pressurerow.mtx")
T15_export = np.zeros(N)
supp = daD != 0
T15_export[supp] = -g15_pr[supp]/daD[supp]

T15_mine, T1_15, T2_15 = build_T(pc15)
d = T15_mine - T15_export
nrm = np.linalg.norm(T15_export[supp])
log("V1 T(pc15) vs b15 raw export: |diff|=%.3e |export|=%.3e rel=%.3e  maxabs=%.3e"
    % (np.linalg.norm(d[supp]), nrm, np.linalg.norm(d[supp])/nrm,
       np.max(np.abs(d[supp]))))

z15 = np.loadtxt(E + "/stageB6_rxc_z_analytic.mtx").reshape(3, N)
rxc = np.loadtxt(E + "/stageB6_rxc_analytic.mtx")
for k, nm in enumerate(["D1", "D2", "D3"]):
    wk = daD*z15[k]
    a = JP_w(wk)
    b = rxc[k*NUNK:(k+1)*NUNK][NV:]
    log("V2 J_P*w(%s) vs anRPd: relL2=%.3e cos=%.15f"
        % (nm, np.linalg.norm(a-b)/np.linalg.norm(b),
           float(np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)))))
da_pat = np.loadtxt(E + "/stageB6_deltaAlpha.mtx")
rxa = np.loadtxt(E + "/stageB6_rxa_analytic.mtx")
a = JP_w(da_pat); b = rxa[3*N:]
log("V2 J_P*sin vs anRPaW: relL2=%.3e"
    % (np.linalg.norm(a-b)/np.linalg.norm(b)))

# ---------------- adjudication ----------------
T16_mine, T1_16, T2_16 = build_T(pc16)
# implied production T16 from the b16 written fields (eta inversion)
maskb = np.zeros(N, dtype=bool)
# designMask
dm_i, _, _ = read_field(S + "/designMask"); mask = dm_i; maskb = mask > 0.5
xp = read_field(S + "/xp")[0]
DEL = 8.0
def heaviside(g, eta):
    e = np.exp(-DEL); lo = g <= eta; out = np.empty_like(g)
    gl = g[lo]/eta; out[lo] = eta*(np.exp(-DEL*(1.0-gl))-(1.0-gl)*e)
    gh = (g[~lo]-eta)/(1.0-eta); out[~lo] = eta+(1.0-eta)*(1.0-np.exp(-DEL*gh)+gh*e)
    return out
def diffvol(g, eta): return np.sum((g[maskb]-heaviside(g[maskb], eta))*V[maskb])
e0, e1 = 1e-4, 0.9999; y0, y1 = diffvol(xp, e0), diffvol(xp, e1)
while (e1-e0) > 1e-10:
    e5 = 0.5*(e0+e1); y5 = diffvol(xp, e5)
    if y0*y5 < 0: e1, y1 = e5, y5
    else: e0, y0 = e5, y5
eta5 = 0.5*(e0+e1)
idx = np.where(maskb)[0]; xpm = xp[idx]; ex = np.exp(-DEL); lo = xpm <= eta5
drho = np.zeros(N); dPde = np.zeros(N)
pe = np.exp(-DEL*(1.0-xpm[lo]/eta5)); drho[idx[lo]] = DEL*pe+ex; dPde[idx[lo]] = pe*(1.0-DEL*xpm[lo]/eta5)-ex
pe2 = np.exp(-DEL*(xpm[~lo]-eta5)/(1.0-eta5)); drho[idx[~lo]] = DEL*pe2+ex; dPde[idx[~lo]] = pe2*(1.0-DEL*(1.0-xpm[~lo])/(1.0-eta5))-ex
Dden = np.sum(mask*dPde*V)
w16 = read_field(S + "/gsenshPressureDrop")[0]
rw = w16[idx]/drho[idx]
Cinv = np.sum((mask*dPde)[idx]*rw)/(Dden+np.sum((mask*dPde*V*(1.0-drho))[idx]/drho[idx]))
gtot = np.zeros(N); gtot[idx] = rw-(V*(1.0-drho)*Cinv)[idx]/drho[idx]
gmom = -daD*np.einsum("ij,ij->i", U, Uc16)*V
gpr = gtot-gmom
T16_implied = np.zeros(N); T16_implied[supp] = -gpr[supp]/daD[supp]
d16 = (T16_mine - T16_implied)[supp]
log("ADJ T(pc16) vs implied production T16: |diff|=%.3e |implied|=%.3e rel=%.3e maxabs=%.3e"
    % (np.linalg.norm(d16), np.linalg.norm(T16_implied[supp]),
       np.linalg.norm(d16)/np.linalg.norm(T16_implied[supp]), np.max(np.abs(d16))))
# per-cell top offenders
o = np.argsort(np.abs(d16))[::-1][:10]
for c in o:
    log("   cell %d: mine=%.6e implied=%.6e dA=%.3e alpha=%.3e"
        % (c, T16_mine[c], T16_implied[c], daD[c], alpha[c]))
# contraction with z for the mine side
for k, nm in enumerate(["D1", "D2", "D3"]):
    lhs = float(np.dot(-T16_mine*daD, z15[k]))
    rhs = -float(np.dot(pc16, rxc[k*NUNK:(k+1)*NUNK][NV:]))
    log("ADJ %s: <-T16_mine*dA, z>=%.10f  -pc16*anRPd=%.10f  diff=%.3e"
        % (nm, lhs, rhs, lhs-rhs))
log("REBUILD_DONE t=%.1fs" % (time.time()-t0))
