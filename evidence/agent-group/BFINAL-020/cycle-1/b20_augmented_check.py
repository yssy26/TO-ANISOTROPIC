#!/usr/bin/env python3
"""BFINAL-020 phase-0d: augmented-system consistency test on w_true.

Hypothesis chain so far (all offline):
  D1-D4 (U-row convective feedback closure)  : insufficient (relL2 unchanged)
  R3 (unrelaxed basis + extras)              : WORSE (gDPf 2.03/0.55/1.97)
The remaining structural candidate: the phi-closure (I - Psi_phi)^-1 of the
eliminated flux map.  Test WITHOUT any big LU: given w_true = (dU, dp) from
the B15 state FD, the augmented system rows must close:

  phi-row : (I - Psi_phi) dphi = Psi_U dU + Psi_p dp + Psi_alpha da   -> solve dphi
  P-row   : div(dphi) = 0                                            -> check
  U-row   : J_UU dU + J_UP dp + J_Uphi dphi + rxd_U = 0               -> check

Psi = production flux map  phiHbyA_relaxed(U; phi) - pEqnFlux(phi; p)
with FULL coefficient feedback Psi_phi (rAU/H/kf channels + relax-source).
"""
import time, json
import numpy as np
import scipy.io
import scipy.sparse as sp
import scipy.sparse.linalg as spla

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
E = "/home/ys/dsH/b15_export"
S = "/home/ys/dsH/b16_states/1"
STATE_CASE = "/home/ys/dsH/b13_probe3"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-020/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N; P0 = NV
ALPHA_REL = 0.4
NU = 5.19009e-05

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-020 phase-0d augmented consistency test ===")

# ---------------- mesh & fields (same loaders as sweep) ----------------
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
nIF = len(nei_all); nF = len(faces)
Sf_all = np.zeros((nF, 3))
for fi in range(nF):
    fv = pts[faces[fi]]
    Sf_all[fi] = 0.5*np.cross(fv, np.roll(fv, -1, axis=0)).sum(axis=0)
magSf_all = np.linalg.norm(Sf_all, axis=1)
cellverts = {}
for fi in range(nF):
    cells = (owner_all[fi], nei_all[fi]) if fi < nIF else (owner_all[fi],)
    for c in cells:
        cellverts.setdefault(c, set()).update(faces[fi])
Ccen = np.zeros((N, 3))
for c, vs in cellverts.items(): Ccen[c] = pts[sorted(vs)].mean(axis=0)
owner = owner_all[:nIF]; nei = nei_all[:nIF]
d_vec = Ccen[nei] - Ccen[owner]
deltaCoeffs = magSf_all[:nIF]/np.einsum("ij,ij->i", Sf_all[:nIF], d_vec)
V = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")[:, 0]
w = np.loadtxt(MESH + "/b16mesh_weights.mtx")
sf = Sf_all[:nIF]
bnd = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
bc_cell = bnd[:, 0].astype(np.int64)
bSf = bnd[:, 1:4]; bmag = bnd[:, 4]; bdel = bnd[:, 5]
Ufixed_b = bnd[:, 7] > 0.5
nBnd = len(bc_cell)
nFaceTot = nIF + nBnd

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
import re
def read_boundary_vals(path):
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
            res[name] = np.array(vals); continue
        vm = re.search(r"value\s+uniform\s+([^\n;]+)", blk)
        if vm:
            u = vm.group(1).strip().rstrip(";").strip()
            if u.startswith("("):
                res[name] = np.array([float(x) for x in u.strip("()").split()])
            else:
                res[name] = np.array([float(u)])
        else:
            res[name] = None
    return res
U_i = read_of_internal(S + "/U").reshape(N, 3)
p_i = read_of_internal(S + "/p")
phi_i = read_of_internal(S + "/phi")
alpha_i = read_of_internal(S + "/alpha")
patch_nf = [("inlet", 84), ("outlet", 84), ("hotInlet", 196), ("hotOutlet", 196),
            ("solidEndWalls", 280), ("bottomWall", 1120), ("topWall", 1120),
            ("sideWalls", 4800)]
phi_bv = read_boundary_vals(S + "/phi")
phi_bvals = np.zeros(nBnd); o = 0
for pn, nf in patch_nf:
    v = phi_bv[pn]
    phi_bvals[o:o+nf] = (v[0] if v.size == 1 else v)
    o += nf
phi_all = np.concatenate([phi_i, phi_bvals])
log("mesh+fields ready (t=%.1fs)" % (time.time()-t0))

# ---------------- rebuild quantities ----------------
upw = phi_i >= 0.0
qp = upw.astype(float)          # pos0(phi): 1 if phi>=0
qn = 1.0 - qp
k_f = NU*magSf_all[:nIF]*deltaCoeffs
lower_f = -qp*phi_i + k_f
upper_f = qn*phi_i + k_f
diag_int = np.zeros(N)
np.add.at(diag_int, owner, qp*phi_i)
np.add.at(diag_int, nei, -qn*phi_i)
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
unclamped = (np.abs(Dbound) >= sumOff)
D_rel = np.maximum(np.abs(Dbound), sumOff)/ALPHA_REL
rAU_rel = 1.0/(D_rel/V)
rAU_u = rAU_rel/ALPHA_REL
mobF = w*rAU_rel[owner] + (1.0-w)*rAU_rel[nei]
kf = mobF*deltaCoeffs*magSf_all[:nIF]
jump = U_i[nei] - U_i[owner]
down_nei = upw
down_cells = np.where(down_nei, nei, owner)
# current H (relaxed basis) for HbyA feedback channels:
# H = (src - offdiag U)/V ; src: relax-source (D_rel - D0)*U with D0=diag_int
src = (D_rel - diag_int)[:, None]*U_i          # + boundary/dev ignored (small)
offU = np.zeros((N, 3))
np.add.at(offU, owner, qn[:, None]*phi_i[:, None]*U_i[nei])
np.add.at(offU, nei, -qp[:, None]*phi_i[:, None]*U_i[owner])
np.add.at(offU, owner, k_f[:, None]*U_i[nei])
np.add.at(offU, nei, k_f[:, None]*U_i[owner])
Hcur = (src - offU)/V[:, None]
HbyA_cur = rAU_rel[:, None]*Hcur
log("rebuild ready (t=%.1fs)" % (time.time()-t0))

# ---------------- Psi maps (faces_total x {U-block, p-block, phi-block}) ----
# PHI: interp contraction (internal faces)
rowsP = []; colsP = []; valsP = []
for k in range(3):
    rowsP.append(np.arange(nIF)); colsP.append(3*owner + k)
    valsP.append(w*sf[:, k])
    rowsP.append(np.arange(nIF)); colsP.append(3*nei + k)
    valsP.append((1.0-w)*sf[:, k])
PHI = sp.coo_matrix((np.concatenate(valsP),
                     (np.concatenate(rowsP), np.concatenate(colsP))),
                    shape=(nIF, 3*N)).tocsr()
# OFF = offdiag action /V (dH_op)
rO = np.concatenate([3*owner, 3*owner+1, 3*owner+2, 3*nei, 3*nei+1, 3*nei+2])
cO = np.concatenate([3*nei, 3*nei+1, 3*nei+2, 3*owner, 3*owner+1, 3*owner+2])
vO = np.concatenate([upper_f, upper_f, upper_f, lower_f, lower_f, lower_f])
vO = vO/np.concatenate([V[owner]]*3 + [V[nei]]*3)
OFF = sp.coo_matrix((vO, (rO, cO)), shape=(3*N, 3*N)).tocsr()
RAUdiag = sp.diags(np.tile(rAU_u, 3))
# SFB boundary flux contraction (uAssignable faces)
assign_b = ~Ufixed_b
bidx = np.where(assign_b)[0]
SFB = sp.lil_matrix((nBnd, 3*N))
for bf in bidx:
    c = bc_cell[bf]
    for k in range(3):
        SFB[bf, 3*c + k] = bSf[bf, k]
# Psi_U (relaxed production map; no dev/bsrc)
PsiU_int = (PHI @ (ALPHA_REL*(RAUdiag @ OFF)
                + (1.0-ALPHA_REL)*sp.identity(3*N))).tocsr()
# boundary (uAssignable): sf_b&(alpha*rAU*dH + (1-alpha)*dU_c)
PsiU_bnd = (SFB.tocsr() @ (ALPHA_REL*(RAUdiag @ OFF)
                + (1.0-ALPHA_REL)*sp.identity(3*N))).tocsr()
PsiU = sp.vstack([PsiU_int, PsiU_bnd]).tocsr()
# Psi_p: internal: -kf (dpn-dpo); boundary outlet: +kf_b dp_c
PsiP = sp.lil_matrix((nFaceTot, N))
PsiP[np.arange(nIF), nei] = -kf
PsiP[np.arange(nIF), owner] = kf
assign_b = ~Ufixed_b
bidx = np.where(assign_b)[0]
for bf in bidx:
    c = bc_cell[bf]
    if True:  # outlet is the only uAssignable patch (p fixedValue)
        PsiP[nIF + bf, c] = rAU_rel[c]*bdel[bf]*bmag[bf]
PsiP = PsiP.tocsr()

# ---- Psi_phi: coefficient feedback ----
# channel (a): dH/dphi (offdiag derivative + relax-source derivative)
#   per face f: own row: -qn*U_n/V_o + (1/a-1)*dDdphi_own*U_o/V_o
#               nei row: +qp*U_o/V_n + (1/a-1)*dDdphi_nei*U_n/V_n
#   dDdphi: own = -qn, nei = +qp
# channel (b): drAU_rel/dphi = -rAU_rel^2 * dDdphi /(a V)  (unclamped only)
# dHbyA/dphi per cell-face: drAU*H + rAU*(dH)
rowsH = []; colsH = []; valsH = []          # (3N x nIF): cell-vec dHbyA per dphi_f
for k in range(3):
    dH_o = (-qn*U_i[nei, k]/V[owner]
            + (1.0/ALPHA_REL - 1.0)*(-qn)*U_i[owner, k]/V[owner])
    dH_n = (+qp*U_i[owner, k]/V[nei]
            + (1.0/ALPHA_REL - 1.0)*(+qp)*U_i[nei, k]/V[nei])
    rowsH.append(3*owner + k); colsH.append(np.arange(nIF)); valsH.append(dH_o)
    rowsH.append(3*nei + k); colsH.append(np.arange(nIF)); valsH.append(dH_n)
DHPHI = sp.coo_matrix((np.concatenate(valsH),
                       (np.concatenate(rowsH), np.concatenate(colsH))),
                      shape=(3*N, nIF)).tocsr()
dDdphi_o = -qn/ (ALPHA_REL*V[owner]) * (-rAU_rel[owner]**2)
dDdphi_n = +qp/ (ALPHA_REL*V[nei]) * (-rAU_rel[nei]**2)
rowsR = []; colsR = []; valsR = []          # (N x nIF): drAU_rel per dphi_f
for k in range(3):
    # HbyA_k = rAU*H_k -> d = drAU*H_k + rAU*dH_k
    vals_ro = dDdphi_o*Hcur[owner, k] + rAU_rel[owner]*(
        -qn*U_i[nei, k]/V[owner] + (1.0/ALPHA_REL-1.0)*(-qn)*U_i[owner, k]/V[owner])
    vals_rn = dDdphi_n*Hcur[nei, k] + rAU_rel[nei]*(
        +qp*U_i[owner, k]/V[nei] + (1.0/ALPHA_REL-1.0)*(+qp)*U_i[nei, k]/V[nei])
    rowsR.append(3*owner + k); colsR.append(np.arange(nIF)); valsR.append(vals_ro)
    rowsR.append(3*nei + k); colsR.append(np.arange(nIF)); valsR.append(vals_rn)
DHB = sp.coo_matrix((np.concatenate(valsR),
                     (np.concatenate(rowsR), np.concatenate(colsR))),
                    shape=(3*N, nIF)).tocsr()
# clamped cells: zero the rAU-feedback rows
clampmask = np.repeat(~unclamped, 3)
DHB = sp.diags((~clampmask).astype(float)).dot(DHB)
# phiHbyA feedback: PHI @ DHB (internal faces), boundary: SFB @ DHB
PsiPHI_hA = (PHI @ DHB).tolil()
PsiPHI_bnd = (SFB.tocsr() @ DHB).tolil()
# kf feedback: dphi_f' affects mob of its two cells -> kf of ALL their faces
# d(kf_f)/d(phi_f') = (w_f*[dmob_o/df'] + (1-w_f)*[dmob_n/df']) *del_f*mag_f
# dmob_c/df' = -rAU_rel_c^2 * dDdphi(c,f')/(a V_c)   (unclamped)
rowsK = []; colsK = []; valsK = []
# build per-cell face incidence
cellfaces = [[] for _ in range(N)]
for fi in range(nIF):
    cellfaces[owner[fi]].append(fi)
    cellfaces[nei[fi]].append(fi)
dmob = {}   # (cell -> list of (face', coeff)) sparse via arrays
dmc_rows = []; dmc_cols = []; dmc_vals = []
for fi in range(nIF):
    if unclamped[owner[fi]]:
        dmc_rows.append(owner[fi]); dmc_cols.append(fi)
        dmc_vals.append(-rAU_rel[owner[fi]]**2*(-qn[fi])/(ALPHA_REL*V[owner[fi]]))
    if unclamped[nei[fi]]:
        dmc_rows.append(nei[fi]); dmc_cols.append(fi)
        dmc_vals.append(-rAU_rel[nei[fi]]**2*(+qp[fi])/(ALPHA_REL*V[nei[fi]]))
DMOB = sp.coo_matrix((dmc_vals, (dmc_rows, dmc_cols)), shape=(N, nIF)).tocsr()
# internal faces: kf_f = (w mob_o + (1-w) mob_n)*del*mag
rowsK = []; colsK = []; valsK = []
for k in range(nIF):
    for cell, wgt in ((owner[k], w[k]), (nei[k], 1.0-w[k])):
        row = DMOB.getrow(cell)
        for j in row.indices:
            rowsK.append(k); colsK.append(j)
            valsK.append(wgt*row.data[list(row.indices).index(j)]
                         *deltaCoeffs[k]*magSf_all[k])
KFB = sp.coo_matrix((valsK, (rowsK, colsK)), shape=(nIF, nIF)).tocsr()
# boundary outlet kf_b = mob_c*del_b*mag_b: cols = faces of c (internal) via DMOB
rowsK2 = []; colsK2 = []; valsK2 = []
for bf in bidx:
    c = bc_cell[bf]
    row = DMOB.getrow(c)
    for j, vv in zip(row.indices, row.data):
        rowsK2.append(bf); colsK2.append(j)
        valsK2.append(vv*bdel[bf]*bmag[bf]*p_i[c])
KFB2 = sp.coo_matrix((valsK2, (rowsK2, colsK2)), shape=(nBnd, nIF)).tocsr()
# (rows are global boundary-face ids 0..nBnd-1; mask to assignable later)
# assemble Psi_phi (nFaceTot x nFaceTot): internal block + boundary block
PsiPHI = sp.lil_matrix((nFaceTot, nFaceTot))
Pint = (sp.hstack([PsiPHI_hA.tocsr(), sp.csr_matrix((nIF, nBnd))])
       - sp.hstack([KFB, sp.csr_matrix((nIF, nBnd))])).tocsr()
Pbnd = (sp.hstack([PsiPHI_bnd.tocsr(), sp.csr_matrix((nBnd, nBnd))])
       - sp.hstack([KFB2, sp.csr_matrix((nBnd, nBnd))])).tocsr()
PsiPHI = sp.vstack([Pint, Pbnd]).tocsr()
log("Psi_phi built nnz=%d (t=%.1fs)" % (PsiPHI.nnz, time.time()-t0))

# ---- Psi_alpha (design direction): from dAlphaDxh*d (stageB6 dirs) ----
dA = np.loadtxt(E + "/stageB6_dalphadxh.mtx") if False else None
# use per-direction deltaAlpha*d from stageB6 exports
deltaAlpha = np.loadtxt(E + "/stageB6_deltaAlpha.mtx") if False else None
dirs = np.loadtxt(E + "/stageB6_dirs.mtx").reshape(3, N)
# direction contraction: da_c = dAlphaDxh_c . d  -> use exported dAlphaDxh? It is
# per-cell scalar field times direction? stageB6_dalphadxh.mtx is per-cell dAlpha/dxh
dAlphaDxh = read_of_internal(S + "/dAlphaDxh")
PsiA = {}
for di, nm in enumerate(["D1", "D2", "D3"]):
    da = dAlphaDxh*dirs[di]
    drAU = -rAU_rel**2/ALPHA_REL*da
    dHsrc = (1.0/ALPHA_REL - 1.0)*da[:, None]*U_i
    dHbyA = drAU[:, None]*Hcur + rAU_rel[:, None]*dHsrc
    v_int = np.einsum("fi,fi->f", sf, w[:, None]*dHbyA[owner]
                      + (1-w)[:, None]*dHbyA[nei])
    v_bnd = np.zeros(nBnd)
    for bf in bidx:
        c = bc_cell[bf]
        v_bnd[bf] = bSf[bf] @ dHbyA[c]
    # pEqnFlux alpha channel: d(kf)*(pn-po) + boundary d(kf_b)*p_c
    v_int += (drAU[owner]*w + drAU[nei]*(1-w))*deltaCoeffs*magSf_all[:nIF]*(
        p_i[nei] - p_i[owner])
    for bf in bidx:
        c = bc_cell[bf]
        v_bnd[bf] += drAU[c]*bdel[bf]*bmag[bf]*p_i[c]
    PsiA[nm] = np.concatenate([v_int, v_bnd])
log("Psi_alpha ready (t=%.1fs)" % (time.time()-t0))

# ---------------- w_true and checks ----------------
rxc = np.loadtxt(E + "/stageB6_rxc_analytic.mtx", dtype=np.float64)
def load_state(dname, tag, sign):
    Uu = np.loadtxt(STATE_CASE + "/stageB2/wstate_%s_h%s_%s_U.mtx" % (dname, tag, sign))
    pp = np.loadtxt(STATE_CASE + "/stageB2/wstate_%s_h%s_%s_p.mtx" % (dname, tag, sign))
    wst = np.zeros(NUNK)
    wst[0:NV:3] = Uu[:, 0]; wst[1:NV:3] = Uu[:, 1]; wst[2:NV:3] = Uu[:, 2]
    wst[P0:] = pp
    return wst
IMPP = spla.splu((sp.identity(nFaceTot) - PsiPHI).tocsc())
log("(I-Psi_phi) factorized (t=%.1fs)" % (time.time()-t0))
res = {}
for di, nm in enumerate(["D1", "D2", "D3"]):
    wt = load_state(nm, "0.001", "p") - load_state(nm, "0.001", "m")
    wt /= 2e-3
    dU = wt[:NV]; dp = wt[P0:]
    rhs = PsiU @ dU + PsiP @ dp + PsiA[nm]
    dphi = IMPP.solve(rhs)
    dphi_op = rhs - PsiA[nm]          # unclosed operator estimate (no closure)
    # P-row check: div(dphi) = 0
    divd = np.zeros(N)
    np.add.at(divd, owner, dphi[:nIF]); np.add.at(divd, nei, -dphi[:nIF])
    np.add.at(divd, bc_cell, dphi[nIF:])
    # U-row check: J_UU dU + J_UP dp + J_Uphi dphi + rxd_U = 0
    rxd = rxc[di*NUNK:(di+1)*NUNK]
    JUphi = np.zeros(NV)
    for k in range(3):
        np.add.at(JUphi, 3*down_cells + k, jump[:, k]*dphi[:nIF])
    # J_UU dU + J_UP dp = -(rxd_U + JUphi) ?  compute via residual identity:
    # rxd_U = V*U*da (momentum row); J_UU*dU + J_UP*dp computed from the rebuild:
    froz = (OFF @ dU)*0  # placeholder: full J_UU dU via matrix pieces below
    res[nm] = {
        "|dphi|": float(np.linalg.norm(dphi)),
        "|dphi_op(unclosed)|": float(np.linalg.norm(dphi_op)),
        "closure_rel": float(np.linalg.norm(dphi - dphi_op)
                             / max(np.linalg.norm(dphi), 1e-300)),
        "|div(dphi)|": float(np.linalg.norm(divd)),
        "|div(dphi)|/|dphi|": float(np.linalg.norm(divd)
                                    / max(np.linalg.norm(dphi), 1e-300)),
    }
    log("%s: closure_relL2(dphi vs unclosed)=%.4e ; |div dphi|/|dphi|=%.3e"
        % (nm, res[nm]["closure_rel"], res[nm]["|div(dphi)|/|dphi|"]))
json.dump(res, open(OUT + "/b20_augmented_check.json", "w"), indent=1)
log("DONE t=%.1fs" % (time.time()-t0))
