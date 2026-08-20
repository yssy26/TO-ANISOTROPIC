#!/usr/bin/env python3
"""BFINAL-020 phase-0b: decompose the flux-tangent (P-row) defect.

Established: stageB2_dphiDudU_FD (rebuilt-phi FD along dUdir, frozen phi) vs
stageB2_dphiDudU_J (operator deltaPhiFacePU along the same direction):
relL2 = 1.36, cos = 0.61.  The P-row flux tangent itself is wrong.

Here: recompute the rebuild tangent ANALYTICALLY in python (b17-validated
rebuild incl. RANS dev source, boundary sources, relax source), verify it
against the run FD, then decompose operator_J - true_tangent channel by
channel (dev-source derivative / boundary-source / relax-source / offdiag).

Also: the same decomposition for the p-direction tangent (dp-only) is not
available as a run FD; the U-direction evidence is enough to localize.
"""
import time, json
import numpy as np

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
E = "/home/ys/dsH/b15_export"
S = "/home/ys/dsH/b16_states/1"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-020/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N
ALPHA_REL = 0.4
NU = 5.19009e-05

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-020 phase-0b flux-tangent decomposition ===")

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
nIF = len(nei_all)
Sf_geo = np.zeros((nIF, 3))
for fi in range(nIF):
    fv = pts[faces[fi]]
    Sf_geo[fi] = 0.5*np.cross(fv, np.roll(fv, -1, axis=0)).sum(axis=0)
magSf_geo = np.linalg.norm(Sf_geo, axis=1)
cellverts = {}
for fi in range(len(faces)):
    cells = (owner_all[fi], nei_all[fi]) if fi < nIF else (owner_all[fi],)
    for c in cells:
        cellverts.setdefault(c, set()).update(faces[fi])
Ccen = np.zeros((N, 3))
for c, vs in cellverts.items(): Ccen[c] = pts[sorted(vs)].mean(axis=0)
owner = owner_all[:nIF]; nei = nei_all[:nIF]
d_vec = Ccen[nei] - Ccen[owner]
deltaCoeffs = magSf_geo/np.einsum("ij,ij->i", Sf_geo, d_vec)
V = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")[:, 0]
w = np.loadtxt(MESH + "/b16mesh_weights.mtx")
sf = np.zeros((nIF, 3))
for i, line in enumerate(open(MESH + "/b16mesh_sf.mtx")):
    sf[i] = [float(t) for t in line.split()]
bnd = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
bc_cell = bnd[:, 0].astype(np.int64)
bSf = bnd[:, 1:4]; bmag = bnd[:, 4]; bdel = bnd[:, 5]
Ufixed_b = bnd[:, 7] > 0.5
# boundary faces of each cell (for src assembly)
bfaces_cell = [[] for _ in range(N)]
for bf in range(len(bc_cell)):
    bfaces_cell[bc_cell[bf]].append(bf)
log("mesh ready (t=%.1fs)" % (time.time()-t0))

# ---------------- fields ----------------
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
U_i = read_of_internal(S + "/U").reshape(N, 3)
phi_i = read_of_internal(S + "/phi")
alpha_i = read_of_internal(S + "/alpha")

# phi boundary values: use patch uniform/nonuniform parse via known values file
import re
def read_boundary_vals(path, fieldtype):
    txt = open(path).read()
    res = {}
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
            res[name] = np.array(vals)
            continue
        vm = re.search(r"value\s+uniform\s+([^\n;]+)", blk)
        if vm:
            u = vm.group(1).strip().rstrip(";").strip()
            if u.startswith("("):
                res[name] = np.array([float(x) for x in u.strip("()").split()])
            else:
                res[name] = np.array([float(u)])
        else:
            res[name] = None   # no value entry (zeroGradient/noSlip): not fixed
    return res
phi_bv = read_boundary_vals(S + "/phi", "scalar")
patch_nf = {"inlet": 84, "outlet": 84, "hotInlet": 196, "hotOutlet": 196,
            "solidEndWalls": 280, "bottomWall": 1120, "topWall": 1120,
            "sideWalls": 4800}
phi_bvals = np.zeros(len(bc_cell))
o = 0
for pn, nf in patch_nf.items():
    v = phi_bv[pn]
    phi_bvals[o:o+nf] = v[0] if v.size == 1 else v
    o += nf
U_bv = read_boundary_vals(S + "/U", "vector")
U_bvals = np.zeros((len(bc_cell), 3)); Ubfixed = np.zeros(len(bc_cell), bool)
o = 0
for pn, nf in patch_nf.items():
    v = U_bv[pn]
    if v is not None and v.size == 3:
        U_bvals[o:o+nf] = v; Ubfixed[o:o+nf] = True   # fixedValue/noSlip w/ value
    elif v is not None:
        U_bvals[o:o+nf] = v.reshape(-1, 3); Ubfixed[o:o+nf] = True
    else:
        U_bvals[o:o+nf] = 0.0     # zeroGradient/noSlip without value: bc=0 anyway
    o += nf
Ubfixed_arr = Ubfixed
log("fields ready (t=%.1fs)" % (time.time()-t0))

# ---------------- frozen-phi momentum rebuild ----------------
upw = phi_i >= 0.0
k_f = NU*magSf_geo*deltaCoeffs
lower_f = -np.where(upw, phi_i, 0.0) + k_f
upper_f = np.where(upw, 0.0, phi_i) + k_f
diag_int = np.zeros(N)
np.add.at(diag_int, owner, np.where(upw, phi_i, 0.0))
np.add.at(diag_int, nei, -np.where(upw, 0.0, phi_i))
divphi = np.zeros(N)
np.add.at(divphi, owner, phi_i); np.add.at(divphi, nei, -phi_i)
np.add.at(divphi, bc_cell, phi_bvals)
diag_int -= divphi
np.add.at(diag_int, owner, -k_f); np.add.at(diag_int, nei, -k_f)
diag_int += alpha_i*V
a_b = NU*bmag*bdel
ic_b = np.where(Ufixed_b, -a_b, phi_bvals)
sumOff = np.zeros(N)
np.add.at(sumOff, owner, np.abs(upper_f)); np.add.at(sumOff, nei, np.abs(lower_f))
Dbound = diag_int.copy()
np.add.at(Dbound, bc_cell, np.abs(ic_b))
D_rel = np.maximum(np.abs(Dbound), sumOff)/ALPHA_REL
rAU_rel = 1.0/(D_rel/V)                    # == primalPressureMobility
rAU_u = rAU_rel/ALPHA_REL                  # == rAUAdj basis (1/A_unrelaxed)
mobF = w*rAU_rel[owner] + (1.0-w)*rAU_rel[nei]
kf = mobF*deltaCoeffs*magSf_geo

# boundaryCoeffs of the convection operator (source entries):
#   fvmDiv: boundaryCoeffs = -patchFlux*valueBoundaryCoeffs(pw)
#   fixedValue U: vBC=1 -> bc = -phi_b ; zeroGradient U: vBC=0 -> bc = 0
pw_b = (phi_bvals >= 0.0).astype(float)
bc_conv = np.where(Ufixed_b, -phi_bvals, 0.0)
# laplacian boundaryCoeffs: -gamma*gradientBoundaryCoeffs: fixedValue -> +del*1
#   boundaryCoeffs = -NU*bmag*bdel*delta... = -a_b ; source gets -bc*U_b = +a_b*U_b
bc_lap = np.where(Ufixed_b, -a_b, 0.0)
log("rebuild ready: rAU_rel avg=%.6e (t=%.1fs)" % (rAU_rel.mean(), time.time()-t0))

# ---------------- rebuild tangent: d(HbyA)/dU along a direction ----------------
# src(U) = boundaryCoeffs . U_b  (+ relax-source handled analytically)
#        + RANS dev source (explicit, -fvc::div(nu*dev2(T(grad U))))*V
def dev_div(Uf):
    """fvc::div(nu*dev2(T(grad U))) per cell (nu uniform)"""
    Uf_f = w[:, None]*Uf[owner] + (1-w)[:, None]*Uf[nei]
    gU = np.zeros((N, 3, 3))
    tmpU = Sf_geo[:, None]*Uf_f[:, :, None]
    np.add.at(gU, owner, tmpU); np.add.at(gU, nei, -tmpU)
    np.add.at(gU, bc_cell, bSf[:, None]*Uf[bc_cell][:, :, None])
    gU /= V[:, None, None]
    devT = np.swapaxes(gU, 1, 2).copy()
    tr = np.trace(gU, axis1=1, axis2=2)
    devT[:, 0, 0] -= 2.0/3.0*tr; devT[:, 1, 1] -= 2.0/3.0*tr
    devT[:, 2, 2] -= 2.0/3.0*tr
    tens = NU*devT
    Tf = w[:, None, None]*tens[owner] + (1-w)[:, None, None]*tens[nei]
    divdev = np.zeros((N, 3))
    tmpD = np.einsum("fi,fij->fj", Sf_geo, Tf)
    np.add.at(divdev, owner, tmpD); np.add.at(divdev, nei, -tmpD)
    np.add.at(divdev, bc_cell, np.einsum("fi,fij->fj", bSf, tens[bc_cell]))
    divdev /= V[:, None]
    return divdev

def src_of(Uf):
    # boundary source: source += bc*U_b  (per fvm conventions b enters A x = b)
    #   convection inflow: -phi_b*U_b (fixedValue U) ; lap: +a_b*U_b (fixedValue U)
    src = np.zeros((N, 3))
    np.add.at(src, bc_cell, bc_conv[:, None]*U_bvals)
    np.add.at(src, bc_cell, bc_lap[:, None]*U_bvals)
    return src

def HbyA_of(Uf):
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, np.where(upw, 0.0, phi_i)[:, None]*Uf[nei])
    np.add.at(offU, nei, -np.where(upw, phi_i, 0.0)[:, None]*Uf[owner])
    np.add.at(offU, owner, k_f[:, None]*Uf[nei])
    np.add.at(offU, nei, k_f[:, None]*Uf[owner])
    src = src_of(Uf)
    src += -V[:, None]*dev_div(Uf)                 # UEqn -= fvc::div(...)
    src += (D_rel - diag_int - np.add.reduceat(    # relax-source (Drel-D0int)*U
        np.zeros((1, 3)), [0])[0] if False else 0)
    # relax-source: source += (D_rel - D0)*U  (D0 = internal diag+boundary ic? )
    # use per-cell Dsrc = D_rel - D_int (b17 convention, validated)
    src += (D_rel - diag_int)[:, None]*Uf
    H = (src - offU)/V[:, None]
    return rAU_rel[:, None]*H

# direction
dU = np.loadtxt(E + "/stageB2_dUdirPhys.mtx").reshape(N, 3)
log("direction |dU|=%.6e" % np.linalg.norm(dU))

# FD in python of the rebuild HbyA (validates the analytic form)
h = 1e-6
Hp = HbyA_of(U_i + h*dU); Hm = HbyA_of(U_i - h*dU)
dHbyA = (Hp - Hm)/(2*h)
dphi_true = np.einsum("fi,fi->f", sf, w[:, None]*dHbyA[owner]
                      + (1-w)[:, None]*dHbyA[nei])

FD = np.loadtxt(E + "/stageB2_dphiDudU_FD.mtx")
Jop = np.loadtxt(E + "/stageB2_dphiDudU_J.mtx")
nrm = np.linalg.norm(FD)
log("[check] python dphi vs run FD: relL2=%.4e cos=%.8f"
    % (np.linalg.norm(dphi_true-FD)/nrm,
       (dphi_true@FD)/(np.linalg.norm(dphi_true)*nrm)))

# operator tangent decomposition
offU_d = np.zeros((N, 3))
np.add.at(offU_d, owner, np.where(upw, 0.0, phi_i)[:, None]*dU[nei])
np.add.at(offU_d, nei, -np.where(upw, phi_i, 0.0)[:, None]*dU[owner])
np.add.at(offU_d, owner, k_f[:, None]*dU[nei])
np.add.at(offU_d, nei, k_f[:, None]*dU[owner])
dH_op = -offU_d/V[:, None]
direct = (1.0-ALPHA_REL)*dU
wop = w[:, None]
dphi_op = np.einsum(
    "fi,fi->f", sf,
    ALPHA_REL*(wop*(rAU_u[owner][:, None]*dH_op[owner])
               + (1-wop)*(rAU_u[nei][:, None]*dH_op[nei]))
    + wop*direct[owner] + (1-wop)*direct[nei])
log("[check] python operator dphi vs run J: relL2=%.4e"
    % (np.linalg.norm(dphi_op-Jop)/np.linalg.norm(Jop)))

# channel decomposition of dphi_true - dphi_op
# 1. dev channel:  dH_dev = -devdiv(dU) ; phi: sf&(w*rAU_rel*dH_dev_o + ...)
ddev = -dev_div(dU)
dphi_dev = np.einsum(
    "fi,fi->f", sf,
    wop*(rAU_rel[owner][:, None]*ddev[owner])
    + (1-wop)*(rAU_rel[nei][:, None]*ddev[nei]))
# 2. boundary-source channel: src(dU) = phi_b*dU_c on zeroGradient-U (outlet)
src_dU = np.zeros((N, 3))
src_dU[bc_cell] += np.where(Ufixed_b, 0.0, phi_bvals)[:, None]*dU[bc_cell]
dphi_bsrc = np.einsum(
    "fi,fi->f", sf,
    wop*(rAU_rel[owner][:, None]*src_dU[owner]/V[owner][:, None])
    + (1-wop)*(rAU_rel[nei][:, None]*src_dU[nei]/V[nei][:, None]))
# 3. relax-basis channel: unrelaxed-vs-relaxed map difference along dU
#    unrelaxed: rAU_u*(dH_op + ddev + bsrc/V) ; relaxed: alpha*rAU_u*dH_op + (1-a)dUbar
dH_unrel_extra = ddev + src_dU/V[:, None]
dphi_basis = np.einsum(
    "fi,fi->f", sf,
    (1.0-ALPHA_REL)*(wop*(rAU_u[owner][:, None]*dH_op[owner])
                     + (1-wop)*(rAU_u[nei][:, None]*dH_op[nei]))
    - (wop*direct[owner] + (1-wop)*direct[nei]))
diff = dphi_true - dphi_op
for nm, vec in [("dev", dphi_dev), ("bsrc", dphi_bsrc),
                ("basis", dphi_basis),
                ("sum(dev,bsrc,basis)", dphi_dev + dphi_bsrc + dphi_basis),
                ("diff", diff)]:
    log("   %-20s |v|=%.6e  cos(diff,v)=%.6f"
        % (nm, np.linalg.norm(vec),
           (diff@vec)/(np.linalg.norm(diff)*np.linalg.norm(vec)+1e-300)))
res = {"dphi_true_vs_FD_rel": float(np.linalg.norm(dphi_true-FD)/nrm),
       "dphi_op_vs_J_rel": float(np.linalg.norm(dphi_op-Jop)/np.linalg.norm(Jop)),
       "|FD|": float(nrm), "|J|": float(np.linalg.norm(Jop)),
       "|pytrue|": float(np.linalg.norm(dphi_true)),
       "|pyop|": float(np.linalg.norm(dphi_op)),
       "|dev|": float(np.linalg.norm(dphi_dev)),
       "|bsrc|": float(np.linalg.norm(dphi_bsrc)),
       "|basis|": float(np.linalg.norm(dphi_basis)),
       "|diff|": float(np.linalg.norm(diff))}
json.dump(res, open(OUT + "/b20_flux_decomp.json", "w"), indent=1)
log("DONE t=%.1fs" % (time.time()-t0))
