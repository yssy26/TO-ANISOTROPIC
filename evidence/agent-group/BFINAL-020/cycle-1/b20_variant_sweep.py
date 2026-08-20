#!/usr/bin/env python3
"""BFINAL-020 phase-0c: variant sweep of the P-row flux-tangent basis.

Channels measured (b20_flux_decomp.json): the operator-vs-rebuild mismatch is
dominated by the relax BASIS (|basis|=1.73e-5, i.e. the BFINAL-003 relaxed-map
choice, alphaRel-scaled hA + (1-alphaRel) direct term), with the outlet
boundary-source channel (|bsrc|=3.2e-7) second and dev negligible.

Variants (each: M_v = M + DeltaM_v, splu, w'_v vs w_true + gDP/J factors):
  R3 : UNRELAXED basis + bsrc + dev      (drop alphaRel from hA/boundary hB,
        remove the direct relax term, kf-conv sign flip, + bsrc + dev channels)
  R2 : UNRELAXED basis only              (R3 minus bsrc/dev)
  S1 : relaxed basis + bsrc + dev        (control: basis kept)
Usage: python3.8 b20_variant_sweep.py R3 R2 S1
"""
import sys, time, json
import numpy as np
import scipy.io
import scipy.sparse as sp
import scipy.sparse.linalg as spla

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
E = "/home/ys/dsH/b15_export"
S = "/home/ys/dsH/b16_states/1"
EXPORT_JT = "/home/ys/dsH/b8_verify_diag/explicitJT.mtx"
STATE_CASE = "/home/ys/dsH/b13_probe3"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-020/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N; P0 = NV; PREF = NV
ALPHA_REL = 0.4
NU = 5.19009e-05
VARIANTS = sys.argv[1:] if len(sys.argv) > 1 else ["R3"]

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-020 phase-0c variant sweep: %s ===" % VARIANTS)

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
nIF = len(nei_all); nF = len(faces)
Sf_geo = np.zeros((nF, 3))
for fi in range(nF):
    fv = pts[faces[fi]]
    Sf_geo[fi] = 0.5*np.cross(fv, np.roll(fv, -1, axis=0)).sum(axis=0)
magSf_geo = np.linalg.norm(Sf_geo, axis=1)
cellverts = {}
for fi in range(nF):
    cells = (owner_all[fi], nei_all[fi]) if fi < nIF else (owner_all[fi],)
    for c in cells:
        cellverts.setdefault(c, set()).update(faces[fi])
Ccen = np.zeros((N, 3))
for c, vs in cellverts.items(): Ccen[c] = pts[sorted(vs)].mean(axis=0)
owner = owner_all[:nIF]; nei = nei_all[:nIF]
d_vec = Ccen[nei] - Ccen[owner]
deltaCoeffs = magSf_geo[:nIF]/np.einsum("ij,ij->i", Sf_geo[:nIF], d_vec)
V = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")[:, 0]
w = np.loadtxt(MESH + "/b16mesh_weights.mtx")
sf = np.zeros((nIF, 3))
for i, line in enumerate(open(MESH + "/b16mesh_sf.mtx")):
    sf[i] = [float(t) for t in line.split()]
bnd = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
bc_cell = bnd[:, 0].astype(np.int64)
bSf = bnd[:, 1:4]; bmag = bnd[:, 4]; bdel = bnd[:, 5]
Ufixed_b = bnd[:, 7] > 0.5
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
phi_i = read_of_internal(S + "/phi")
alpha_i = read_of_internal(S + "/alpha")
patch_nf = [("inlet", 84), ("outlet", 84), ("hotInlet", 196), ("hotOutlet", 196),
            ("solidEndWalls", 280), ("bottomWall", 1120), ("topWall", 1120),
            ("sideWalls", 4800)]
phi_bv = read_boundary_vals(S + "/phi")
phi_bvals = np.zeros(len(bc_cell)); o = 0
for pn, nf in patch_nf:
    v = phi_bv[pn]
    phi_bvals[o:o+nf] = (v[0] if v.size == 1 else v)
    o += nf
log("fields ready (t=%.1fs)" % (time.time()-t0))

# ---------------- momentum rebuild ----------------
upw = phi_i >= 0.0
k_f = NU*magSf_geo[:nIF]*deltaCoeffs
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
rAU_rel = 1.0/(D_rel/V)
rAU_u = rAU_rel/ALPHA_REL
mobF = w*rAU_rel[owner] + (1.0-w)*rAU_rel[nei]
kf = mobF*deltaCoeffs*magSf_geo[:nIF]
contRes = divphi.copy()
jump = U_i[nei] - U_i[owner]
down_nei = upw
log("rebuild ready (t=%.1fs)" % (time.time()-t0))

# ---------------- load M, rhs, truth ----------------
M = scipy.io.mmread(EXPORT_JT).tocsr(); M.sort_indices()
log("M loaded nnz=%d (t=%.1fs)" % (M.nnz, time.time()-t0))
def load_rhs_export(path):
    v = np.asarray(scipy.io.mmread(path)).ravel().reshape(N, 4)
    b = np.zeros(NUNK)
    b[0:NV:3] = v[:, 0]; b[1:NV:3] = v[:, 1]; b[2:NV:3] = v[:, 2]
    b[P0:] = v[:, 3]
    return b
bPD = load_rhs_export(E + "/explicitRhs_pressureDrop.mtx")
bTC = load_rhs_export(E + "/explicitRhs_thermalCoupling.mtx")
rxc = np.loadtxt(E + "/stageB6_rxc_analytic.mtx", dtype=np.float64)
def load_state(dname, tag, sign):
    Uu = np.loadtxt(STATE_CASE + "/stageB2/wstate_%s_h%s_%s_U.mtx" % (dname, tag, sign))
    pp = np.loadtxt(STATE_CASE + "/stageB2/wstate_%s_h%s_%s_p.mtx" % (dname, tag, sign))
    wst = np.zeros(NUNK)
    wst[0:NV:3] = Uu[:, 0]; wst[1:NV:3] = Uu[:, 1]; wst[2:NV:3] = Uu[:, 2]
    wst[P0:] = pp
    return wst
w_true = {}
for nm in ["D1", "D2", "D3"]:
    w_true[nm] = (load_state(nm, "0.001", "p") - load_state(nm, "0.001", "m"))/2e-3
cz = np.load("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/"
             "BFINAL-015/cycle-1/wprime_cache.npz")
wprime_old = {k: cz[k] for k in cz.files}
mask = np.ones(NUNK, dtype=bool); mask[PREF] = False
log("data ready (t=%.1fs)" % (time.time()-t0))

# ---------------- dev operator matrix K_dev (3N x 3N) ----------------
# devdiv(U)[c][j] = (1/V) sum_f s * sum_i Sf[i] * Tf[i][j],
# Tf = interp(tens), tens[c] = nu*(gU[c]^T - (2/3) I tr gU[c]),
# gU[c][i][j] = (1/V) sum_f s * Sf[j] * Uf[i], Uf = interp(U) (+bnd: U_c)
def build_dev():
    # G: gU rows (c*9 + i*3 + j) <- U cols (m*3 + k)
    rowsG = []; colsG = []; valsG = []
    for i in range(3):
        for j in range(3):
            r = np.repeat(9*np.arange(N) + i*3 + j, 2)
            # internal faces: own gets +w*sf[j]/V_own (col U_o+i), +(1-w)... 
            cO = 3*owner + i; cN = 3*nei + i
            rowsG.append(9*owner + i*3 + j); colsG.append(cO)
            valsG.append(w*sf[:, j]/V[owner])
            rowsG.append(9*owner + i*3 + j); colsG.append(cN)
            valsG.append((1.0-w)*sf[:, j]/V[owner])
            rowsG.append(9*nei + i*3 + j); colsG.append(cO)
            valsG.append(-w*sf[:, j]/V[nei])
            rowsG.append(9*nei + i*3 + j); colsG.append(cN)
            valsG.append(-(1.0-w)*sf[:, j]/V[nei])
    # boundary faces: gU[c] += Sf_b[j]*U_c[i]/V_c
    for i in range(3):
        for j in range(3):
            rowsG.append(9*bc_cell + i*3 + j); colsG.append(3*bc_cell + i)
            valsG.append(bSf[:, j]/V[bc_cell])
    G = sp.coo_matrix((np.concatenate(valsG),
                       (np.concatenate(rowsG), np.concatenate(colsG))),
                      shape=(9*N, 3*N)).tocsr()
    # TENS: tens[c][i][j] = nu*(gU[c][j][i] - (2/3) delta_ij sum_l gU[c][l][l])
    rowsT = []; colsT = []; valsT = []
    for c in range(N):
        for i in range(3):
            for j in range(3):
                r = 9*c + i*3 + j
                rowsT.append(r); colsT.append(9*c + j*3 + i); valsT.append(NU)
                if i == j:
                    for l in range(3):
                        rowsT.append(r); colsT.append(9*c + l*3 + l)
                        valsT.append(-2.0/3.0*NU)
    T = sp.coo_matrix((valsT, (rowsT, colsT)), shape=(9*N, 9*N)).tocsr()
    # DIV: devdiv[c][j] = (1/V_c) sum_f s sum_i Sf[i]*Tf[i][j], Tf=interp(tens)
    rowsD = []; colsD = []; valsD = []
    for i in range(3):
        for j in range(3):
            rowsD.append(3*owner + j); colsD.append(9*owner + i*3 + j)
            valsD.append(sf[:, i]/V[owner])
            rowsD.append(3*owner + j); colsD.append(9*nei + i*3 + j)
            valsD.append((1.0-w)*sf[:, i]/V[owner])
            rowsD.append(3*nei + j); colsD.append(9*owner + i*3 + j)
            valsD.append(-w*sf[:, i]/V[nei])
            rowsD.append(3*nei + j); colsD.append(9*nei + i*3 + j)
            valsD.append(-(1.0-w)*sf[:, i]/V[nei])
            # boundary: Tf_b = tens[c]
            rowsD.append(3*bc_cell + j); colsD.append(9*bc_cell + i*3 + j)
            valsD.append(bSf[:, i]/V[bc_cell])
    D = sp.coo_matrix((np.concatenate(valsD),
                       (np.concatenate(rowsD), np.concatenate(colsD))),
                      shape=(3*N, 9*N)).tocsr()
    return D @ T @ G
K_dev = build_dev().tocsr()
log("K_dev built nnz=%d (t=%.1fs)" % (K_dev.nnz, time.time()-t0))

# face-flux map from a cell vector field: PHI(X)_f = sf & (w X_o + (1-w) X_n)
rowsP_ = []; colsP_ = []; valsP_ = []
for k in range(3):
    rowsP_.append(np.arange(nIF)); colsP_.append(3*owner + k)
    valsP_.append(w*sf[:, k])
    rowsP_.append(np.arange(nIF)); colsP_.append(3*nei + k)
    valsP_.append((1.0-w)*sf[:, k])
PHI = sp.coo_matrix((np.concatenate(valsP_),
                     (np.concatenate(rowsP_), np.concatenate(colsP_))),
                    shape=(nIF, 3*N)).tocsr()
# div map: R_P(o) += phi_f, R_P(n) -= phi_f
DIVM = sp.coo_matrix(
    (np.concatenate([np.ones(nIF), -np.ones(nIF)]),
     (np.concatenate([P0 + owner, P0 + nei]),
      np.concatenate([np.arange(nIF), np.arange(nIF)]))),
    shape=(NUNK, nIF)).tocsr()
# boundary flux map (uAssignable faces only): row = boundary face, col = 3c+k
assign_b = ~Ufixed_b
brows = np.repeat(np.where(assign_b)[0], 3)
bcols = np.tile(3*bc_cell[assign_b][:, None], (1, 3)).ravel() + np.tile(np.arange(3), assign_b.sum())
bvals = (bSf[assign_b][:, None]*np.ones((1, 3))).ravel()
SFB = sp.coo_matrix((bvals, (brows, bcols)),
                    shape=(len(bc_cell), 3*N)).tocsr()
# boundary P-row scatter: P(c) += flux_b  (uAssignable faces only)
BSCAT = sp.coo_matrix(
    (np.ones(assign_b.sum()),
     (P0 + bc_cell[assign_b], np.where(assign_b)[0])),
    shape=(NUNK, len(bc_cell))).tocsr()
# jump map: U(down,c) += jump_c * phi_f
jrows = []; jcols = []; jvals = []
down_cells = np.where(down_nei, nei, owner)
for c in range(3):
    jrows.append(3*down_cells + c); jcols.append(np.arange(nIF))
    jvals.append(jump[:, c])
JUMPM = sp.coo_matrix((np.concatenate(jvals),
                       (np.concatenate(jrows), np.concatenate(jcols))),
                      shape=(3*N, nIF)).tocsr()
JUMPM_FULL = sp.vstack([JUMPM, sp.csr_matrix((N, nIF))]).tocsr()
# offdiag action /V (dH_op)
OFF = sp.coo_matrix(
    (np.concatenate([upper_f, upper_f, upper_f, lower_f, lower_f, lower_f]),
     (np.concatenate([3*owner+0, 3*owner+1, 3*owner+2,
                      3*nei+0, 3*nei+1, 3*nei+2]),
      np.concatenate([3*nei+0, 3*nei+1, 3*nei+2,
                      3*owner+0, 3*owner+1, 3*owner+2]))),
    shape=(3*N, 3*N)).tocsr()
rO = np.concatenate([3*owner, 3*owner+1, 3*owner+2, 3*nei, 3*nei+1, 3*nei+2])
cO = np.concatenate([3*nei, 3*nei+1, 3*nei+2, 3*owner, 3*owner+1, 3*owner+2])
vO = np.concatenate([upper_f, upper_f, upper_f, lower_f, lower_f, lower_f])
vO = vO/np.concatenate([V[owner]]*3 + [V[nei]]*3)
OFF = sp.coo_matrix((vO, (rO, cO)), shape=(3*N, 3*N)).tocsr()
RAUdiag = sp.diags(np.tile(rAU_u, 3))
# bsrc diagonal (outlet zeroGradient-U cells): dH += phi_b/V * dU
bsrcD = np.zeros(3*N)
_bcell = bc_cell[assign_b]
np.add.at(bsrcD, ((3*_bcell)[:, None] + np.arange(3)).ravel(),
          np.repeat(phi_bvals[assign_b]/V[_bcell], 3))
BSRC = sp.diags(bsrcD)
EXTRAS = (-K_dev + BSRC).tocsr()
aR = ALPHA_REL; oneMa = 1.0 - ALPHA_REL
# old flux maps
FLUX_old_P = PHI @ (aR*(RAUdiag @ OFF) + (1-aR)*sp.identity(3*N))
FLUX_old_jump = PHI @ (RAUdiag @ OFF)

def build_delta(variant):
    unrelaxed = variant in ("R2", "R3")
    with_extras = variant in ("R3", "S1")
    if unrelaxed:
        FLUX_new = PHI @ (RAUdiag @ (OFF + (EXTRAS if with_extras else 0)))
        FLUXB_old = SFB @ (aR*(RAUdiag @ OFF) + (1-aR)*sp.identity(3*N))
        FLUXB_new = SFB @ (RAUdiag @ (OFF + (EXTRAS if with_extras else 0)))
    else:
        FLUX_new = PHI @ (aR*(RAUdiag @ (OFF + (EXTRAS if with_extras else 0)))
                          + (1-aR)*sp.identity(3*N))
        FLUXB_old = SFB @ (aR*(RAUdiag @ OFF) + (1-aR)*sp.identity(3*N))
        FLUXB_new = SFB @ (aR*(RAUdiag @ (OFF + (EXTRAS if with_extras else 0)))
                           + (1-aR)*sp.identity(3*N))
    dP = (FLUX_new - FLUX_old_P).tocsc()      # faces x 3N (U-cols only)
    dJ = (FLUX_new - FLUX_old_jump).tocsc()   # jump basis difference
    dP_full = sp.hstack([dP, sp.csr_matrix((nIF, N))]).tocsr()
    dJ_full = sp.hstack([dJ, sp.csr_matrix((nIF, N))]).tocsr()
    # kf sign flip on the jump part only: old jump kf-part = +kf*(dpn-dpo)
    rowsK = []; colsK = []; valsK = []
    rowsK.append(np.arange(nIF)); colsK.append(P0 + nei)
    valsK.append(-2.0*kf)     # new(-kf) - old(+kf) = -2kf on P(n) col
    rowsK.append(np.arange(nIF)); colsK.append(P0 + owner)
    valsK.append(2.0*kf)
    KFLIP = sp.coo_matrix((np.concatenate(valsK),
                           (np.concatenate(rowsK), np.concatenate(colsK))),
                          shape=(nIF, NUNK)).tocsr()
    dB_full = sp.hstack([(FLUXB_new - FLUXB_old).tocsr(),
                         sp.csr_matrix((len(bc_cell), N))]).tocsr()
    DeltaJ = (DIVM @ dP_full + BSCAT @ dB_full
              + JUMPM_FULL @ (dJ_full + KFLIP)).tolil()
    # D4 diagonal: -contRes
    for c in range(3):
        for cell in range(N):
            v = -contRes[cell]
            if v != 0.0:
                DeltaJ[3*cell + c, 3*cell + c] += v
    return DeltaJ.tocsc()

for variant in VARIANTS:
    log("--- variant %s ---" % variant)
    DJ = build_delta(variant)
    DeltaM = DJ.T.tocsc()
    Mv = (M + DeltaM).tolil()
    Mv[:, PREF] = 0.0; Mv[PREF, PREF] = 1.0
    Mv = Mv.tocsc(); Mv.sort_indices()
    log("  factorizing (nnz=%d) ..." % Mv.nnz)
    lu = spla.splu(Mv)
    lam = lu.solve(bPD)
    out = {}
    truth = {"D1": -2.537459771850754, "D2": 0.5234094521578995,
             "D3": 6.477215547073193}
    for di, nm in enumerate(["D1", "D2", "D3"]):
        rxd = rxc[di*NUNK:(di+1)*NUNK]
        rhs = -rxd.copy(); rhs[PREF] = 0.0
        x = lu.solve(rhs, trans="T"); x[PREF] = 0.0
        wt = w_true[nm]; xo = wprime_old[nm]
        def mets(wp):
            d = wp - wt
            cosv = float(wp[mask] @ wt[mask] /
                         (np.linalg.norm(wp[mask])*np.linalg.norm(wt[mask])))
            rel = float(np.linalg.norm(d[mask])/np.linalg.norm(wt[mask]))
            relU = float(np.linalg.norm(d[:NV])/np.linalg.norm(wt[:NV]))
            relP = float(np.linalg.norm(d[P0:])/np.linalg.norm(wt[P0:]))
            return cosv, rel, relU, relP
        cN, rN, ruN, rpN = mets(x); cO, rO, ruO, rpO = mets(xo)
        gdp = float(bPD @ x)/float(bPD @ wt)
        jj = float(bTC @ x)/float(bTC @ wt)
        ld = -float(lam @ rxd)
        out[nm] = {"old": {"cos": cO, "relL2": rO, "relU": ruO, "relP": rpO},
                   "new": {"cos": cN, "relL2": rN, "relU": ruN, "relP": rpN},
                   "gdp_factor": gdp, "J_factor": jj, "lam_direct": ld,
                   "lam_direct_ratio": ld/truth[nm]}
        log("  %s: relL2 %.4f -> %.4f (relP %.3f -> %.3f) gDPf %.4f -> %.4f "
            "Jf %.4f lamdir-ratio -> %.4f" %
            (nm, rO, rN, rpO, rpN,
             float(bPD @ xo)/float(bPD @ wt), gdp, jj, ld/truth[nm]))
    json.dump(out, open(OUT + "/b20_variant_%s.json" % variant, "w"), indent=1)
log("SWEEP_DONE t=%.1fs" % (time.time()-t0))
