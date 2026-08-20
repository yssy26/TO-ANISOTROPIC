#!/usr/bin/env python3
"""BFINAL-021 E2: P-row tangent identity across the eps ladder (zero-code).

Context established this round (E1 data archaeology):
  * The exported operator/rxc/rhs live at state A (main-loop, molecular
    viscosity matrix, 44 NS.H correctors).
  * The FD truth (wstate) lives at state B (stageB2 phase-1 solveFullSST,
    246 iters at 1e-6, frozen-TURBULENT viscosity chi*nut active).
  * A-B: U relL2 13.6%, p relL2 10.6% (mean p shift -4800 Pa).
  * w(h) ladder self-convergence 0.04%-1.9%  =>  FD noise is NOT O(1);
    any O(1) identity failure is systematic.

This script measures, per direction and per eps:
  1. identity at A-basis:  r_P(w) = div(Psi_U dU + Psi_p dp) + rxd_P_export
     (the BFINAL-020 "sign-level failure" construct), cos vs -rxd_P,
     |r_P|/|rxd_P| -- its eps-dependence (true mismatch: eps-flat).
  2. same with the offline-rebuilt alpha channel (-PsiA) instead of the
     exported rxd_P  ->  representation-consistency of the two rxd
     constructions (E4 input).
  3. lambda* control: r_P(w') for the cached A-point operator solve w'
     (must be ~0 by construction; validates this script's machinery
     against the B20 numbers).
  4. |w(h)-w'|/|w'| ladder (the B15 "defect-1" mismatch gauge).

No solver runs, no source changes. python3.8 + numpy/scipy only.
"""
import time, json
import numpy as np
import scipy.sparse as sp

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
E = "/home/ys/dsH/b15_export"
S = "/home/ys/dsH/b16_states/1"
STATE_CASE = "/home/ys/dsH/b13_probe3"
WPC = ("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/"
       "BFINAL-015/cycle-1/wprime_cache.npz")
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-021/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N; P0 = NV
ALPHA_REL = 0.4
NU = 5.19009e-05

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-021 E2: identity ladder (A-basis) ===")

# ---------------- mesh & fields (loaders identical to b20) ----------------
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
log("mesh+fields ready (t=%.1fs)" % (time.time()-t0))

# ---------------- A-basis momentum rebuild (identical to b20) -------------
upw = phi_i >= 0.0
qp = upw.astype(float); qn = 1.0 - qp
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
src = (D_rel - diag_int)[:, None]*U_i
offU = np.zeros((N, 3))
np.add.at(offU, owner, qn[:, None]*phi_i[:, None]*U_i[nei])
np.add.at(offU, nei, -qp[:, None]*phi_i[:, None]*U_i[owner])
np.add.at(offU, owner, k_f[:, None]*U_i[nei])
np.add.at(offU, nei, k_f[:, None]*U_i[owner])
Hcur = (src - offU)/V[:, None]
log("A-basis rebuild ready (t=%.1fs)" % (time.time()-t0))

# ---------------- Psi maps (identical to b20_augmented_check) -------------
rowsP = []; colsP = []; valsP = []
for k in range(3):
    rowsP.append(np.arange(nIF)); colsP.append(3*owner + k)
    valsP.append(w*sf[:, k])
    rowsP.append(np.arange(nIF)); colsP.append(3*nei + k)
    valsP.append((1.0-w)*sf[:, k])
PHI = sp.coo_matrix((np.concatenate(valsP),
                     (np.concatenate(rowsP), np.concatenate(colsP))),
                    shape=(nIF, 3*N)).tocsr()
rO = np.concatenate([3*owner, 3*owner+1, 3*owner+2, 3*nei, 3*nei+1, 3*nei+2])
cO = np.concatenate([3*nei, 3*nei+1, 3*nei+2, 3*owner, 3*owner+1, 3*owner+2])
vO = np.concatenate([upper_f, upper_f, upper_f, lower_f, lower_f, lower_f])
vO = vO/np.concatenate([V[owner]]*3 + [V[nei]]*3)
OFF = sp.coo_matrix((vO, (rO, cO)), shape=(3*N, 3*N)).tocsr()
RAUdiag = sp.diags(np.tile(rAU_u, 3))
assign_b = ~Ufixed_b
bidx = np.where(assign_b)[0]
SFB = sp.lil_matrix((nBnd, 3*N))
for bf in bidx:
    c = bc_cell[bf]
    for k in range(3):
        SFB[bf, 3*c + k] = bSf[bf, k]
PsiU_int = (PHI @ (ALPHA_REL*(RAUdiag @ OFF)
                + (1.0-ALPHA_REL)*sp.identity(3*N))).tocsr()
PsiU_bnd = (SFB.tocsr() @ (ALPHA_REL*(RAUdiag @ OFF)
                + (1.0-ALPHA_REL)*sp.identity(3*N))).tocsr()
PsiU = sp.vstack([PsiU_int, PsiU_bnd]).tocsr()
PsiP = sp.lil_matrix((nFaceTot, N))
PsiP[np.arange(nIF), nei] = -kf
PsiP[np.arange(nIF), owner] = kf
for bf in bidx:
    c = bc_cell[bf]
    PsiP[nIF + bf, c] = rAU_rel[c]*bdel[bf]*bmag[bf]
PsiP = PsiP.tocsr()
log("PsiU/PsiP ready (t=%.1fs)" % (time.time()-t0))

def divFaces(dphi):
    d = np.zeros(N)
    np.add.at(d, owner, dphi[:nIF]); np.add.at(d, nei, -dphi[:nIF])
    np.add.at(d, bc_cell, dphi[nIF:])
    return d

# ---------------- rxd / dirs / PsiA (A-basis) ------------------------------
rxc = np.loadtxt(E + "/stageB6_rxc_analytic.mtx", dtype=np.float64)
dAlphaDxh = read_of_internal(S + "/dAlphaDxh")
dirs = np.loadtxt(E + "/stageB6_dirs.mtx").reshape(3, N)
PsiA = {}
rxdP_exp = {}
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
    v_int += (drAU[owner]*w + drAU[nei]*(1-w))*deltaCoeffs*magSf_all[:nIF]*(
        p_i[nei] - p_i[owner])
    for bf in bidx:
        c = bc_cell[bf]
        v_bnd[bf] += drAU[c]*bdel[bf]*bmag[bf]*p_i[c]
    PsiA[nm] = np.concatenate([v_int, v_bnd])
    rxd = rxc[di*NUNK:(di+1)*NUNK]
    rxdP_exp[nm] = rxd[P0:]
log("PsiA/rxd ready (t=%.1fs)" % (time.time()-t0))

# representation check: my PsiA rebuild vs exported rxd_P (both A-basis)
log("\n[E4-input] PsiA rebuild vs exported rxd_P (A-basis):")
rep = {}
for nm in ["D1", "D2", "D3"]:
    a = divFaces(PsiA[nm]); b = rxdP_exp[nm]
    cosv = float(a @ b/(np.linalg.norm(a)*np.linalg.norm(b)))
    rel = float(np.linalg.norm(a-b)/np.linalg.norm(b))
    rep[nm] = {"cos": cosv, "relL2": rel}
    log("    %s: cos=%.6f relL2=%.4e" % (nm, cosv, rel))

# ---------------- states / wprime ------------------------------------------
def load_state(dname, tag, sign):
    Uu = np.loadtxt(STATE_CASE + "/stageB2/wstate_%s_h%s_%s_U.mtx" % (dname, tag, sign))
    pp = np.loadtxt(STATE_CASE + "/stageB2/wstate_%s_h%s_%s_p.mtx" % (dname, tag, sign))
    wst = np.zeros(NUNK)
    wst[0:NV:3] = Uu[:, 0]; wst[1:NV:3] = Uu[:, 1]; wst[2:NV:3] = Uu[:, 2]
    wst[P0:] = pp
    return wst
cz = np.load(WPC)
wprime = {k: cz[k] for k in cz.files}
ladder = {"D1": ["0.0003", "0.001", "0.003"],
          "D2": ["0.0003", "0.001"],
          "D3": ["0.0003", "0.001"]}

def metrics(dU, dp, rxdP):
    dphi = PsiU @ dU + PsiP @ dp
    dv = divFaces(dphi)
    r_P = dv + rxdP
    cosv = float(dv @ (-rxdP)/(np.linalg.norm(dv)*np.linalg.norm(rxdP)))
    return {"cos(div,-rxdP)": cosv,
            "|r_P|/|rxdP|": float(np.linalg.norm(r_P)/np.linalg.norm(rxdP)),
            "|div|/|rxdP|": float(np.linalg.norm(dv)/np.linalg.norm(rxdP))}

res = {"rep_check_PsiA_vs_rxdP": rep, "ladder": {}, "wprime_control": {}}
log("\n[E2] identity ladder (A-basis operator, exported rxd_P):")
for nm in ["D1", "D2", "D3"]:
    res["ladder"][nm] = {}
    for tag in ladder[nm]:
        h = float(tag)
        wt = (load_state(nm, tag, "p") - load_state(nm, tag, "m"))/(2*h)
        m1 = metrics(wt[:NV], wt[P0:], rxdP_exp[nm])
        m2 = metrics(wt[:NV], wt[P0:], divFaces(PsiA[nm]))
        dwp = wt - wprime[nm]
        rel_w = float(np.linalg.norm(dwp)/np.linalg.norm(wprime[nm]))
        res["ladder"][nm][tag] = {"vs_exported_rxd": m1, "vs_PsiA_rxd": m2,
                                  "|w-w'|/|w'|": rel_w}
        log("  %s h=%s: cos=%.4f |r_P|/|rxdP|=%.3f (vs PsiA: cos=%.4f"
            " |r_P|=%.3f)  |w-w'|/|w'|=%.4f"
            % (nm, tag, m1["cos(div,-rxdP)"], m1["|r_P|/|rxdP|"],
               m2["cos(div,-rxdP)"], m2["|r_P|/|rxdP|"], rel_w))

log("\n[control] r_P(w') (operator's own solve; must be ~0):")
for nm in ["D1", "D2", "D3"]:
    m = metrics(wprime[nm][:NV], wprime[nm][P0:], rxdP_exp[nm])
    res["wprime_control"][nm] = m
    log("    %s: |r_P(w')|/|rxdP|=%.3e cos=%.6f"
        % (nm, m["|r_P|/|rxdP|"], m["cos(div,-rxdP)"]))

json.dump(res, open(OUT + "/b21_e2_identity_ladder.json", "w"), indent=1)
log("\nDONE t=%.1fs" % (time.time()-t0))
