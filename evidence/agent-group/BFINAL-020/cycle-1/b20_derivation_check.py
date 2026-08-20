#!/usr/bin/env python3
"""BFINAL-020 phase-0: offline pre-validation of the derived operator fix.

Derivation (DERIVATION.md section 3-4): the U-row convective flux feedback of
`bounded Gauss upwind` nets to  jump*dphi*  on the DOWNWIND row only, with
dphi* = deltaPhiFacePU (the SAME relaxed closed-flux tangent the P rows use).
The current operator multiplies the jump by the OLD deltaPhiFace (no alphaRel
on the hA part, no relax-source part, kf sign flipped) and misses the Sp
diagonal -contRes*dU.

This script builds DeltaM = J_new^T - J_old^T from the exported state
(b16 mesh + b16_states) plus the exported explicitJT.mtx, then adjudicates:

  V1  slot-level rebuild self-checks against M (kf / P-row U-slots / hA slots)
  BL  baseline reproduction of b17 lambda*-direct (ratios 1.768/0.856/2.274)
  V2  w'_new vs w_true  (gate: cos>0.999, relL2<5e-2; old relL2 0.72-1.13)
  V3  gDP factor prediction  bPD^T w'_new / bPD^T w_true  (gate 1+-0.10)
  V4  lambda-direct_new / truth

Run: python3.8 b20_derivation_check.py   (~2 splu factorizations, ~40 min)
"""
import time, json, re
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

N = 33600; NV = 3*N; NUNK = 4*N; P0 = NV
ALPHA_REL = 0.4
NU = 5.19009e-05

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-020 phase-0 offline pre-validation ===")

# ---------------- mesh (b17-style geometry for deltaCoeffs) ----------------
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
assert np.allclose(sf, Sf_geo, rtol=0, atol=1e-15), "sf mismatch b16mesh vs geometry"
bnd = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
bc_cell = bnd[:, 0].astype(np.int64)
bSf = bnd[:, 1:4]; bmag = bnd[:, 4]; bdel = bnd[:, 5]
Ufixed_b = bnd[:, 7] > 0.5
log("mesh ready: nIF=%d nBnd=%d (t=%.1fs)" % (nIF, len(bnd), time.time()-t0))

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
def read_of_boundary_list(path, names):
    """boundary values per requested patch name, in the field's boundaryField"""
    txt = open(path).read()
    out = {}
    for name in names:
        m = re.search(r"\b%s\s*\n?\s*\{" % re.escape(name), txt)
        if m is None:
            raise RuntimeError("patch %s not found in %s" % (name, path))
        i = m.end()
        depth = 1; j = i
        while depth > 0:
            if txt[j] == "{": depth += 1
            elif txt[j] == "}": depth -= 1
            j += 1
        blk = txt[i:j]
        vm = re.search(r"value\s+nonuniform[^0-9\-]*\d+\s*\(", blk)
        if vm:
            k = blk.index("(", vm.end()-1)
            kk = blk.index(")", k)
            vals = []
            for tok in blk[k+1:kk].replace("(", " ").replace(")", " ").split():
                try: vals.append(float(tok))
                except ValueError: pass
            out[name] = np.array(vals)
            continue
        vm = re.search(r"value\s+uniform\s+([^\n;]+)", blk)
        if vm:
            u = vm.group(1).strip().rstrip(";").strip()
            if u.startswith("("):
                out[name] = np.array([float(x) for x in u.strip("()").split()])
            else:
                out[name] = np.array([float(u)])
            continue
        raise RuntimeError("no value entry for patch %s in %s" % (name, path))
    return out

# patch table from polyMesh/boundary (name, nFaces) in file order
btxt = open(POLY + "/boundary").read()
pat = re.compile(r"(\w+)\s*\{[^}]*?nFaces\s+(\d+)\s*;", re.S)
pairs = pat.findall(btxt)
patch_names = [p for p, _ in pairs]
patch_nfaces = [int(nf) for _, nf in pairs]

U_i = read_of_internal(S + "/U").reshape(N, 3)
p_i = read_of_internal(S + "/p")
phi_i = read_of_internal(S + "/phi")
alpha_i = read_of_internal(S + "/alpha")
phi_b = read_of_boundary_list(S + "/phi", patch_names)
phi_bvals = np.zeros(len(bc_cell))
o = 0
for pn, nf in zip(patch_names, patch_nfaces):
    v = phi_b[pn]
    if v.size == 1: phi_bvals[o:o+nf] = v[0]
    elif v.size == nf: phi_bvals[o:o+nf] = v
    else: raise RuntimeError("phi boundary size mismatch " + pn)
    o += nf
assert o == len(bc_cell)
log("fields ready (t=%.1fs)" % (time.time()-t0))

# ---------------- momentum matrix rebuild (unrelaxed, frozen phi) ----------
upw = phi_i >= 0.0
k_f = NU*magSf_geo*deltaCoeffs
lower_f = -np.where(upw, phi_i, 0.0) + k_f     # ldu lower  (conv + lap)
upper_f = np.where(upw, 0.0, phi_i) + k_f      # ldu upper
diag_int = np.zeros(N)
np.add.at(diag_int, owner, np.where(upw, phi_i, 0.0))
np.add.at(diag_int, nei, -np.where(upw, 0.0, phi_i))
# bounded Sp: diag += -(sum_int s*phi + sum_b phi_b)
divphi = np.zeros(N)
np.add.at(divphi, owner, phi_i)
np.add.at(divphi, nei, -phi_i)
np.add.at(divphi, bc_cell, phi_bvals)
diag_int -= divphi
# laplacian
np.add.at(diag_int, owner, -k_f)
np.add.at(diag_int, nei, -k_f)
# Sp(alpha)
diag_int += alpha_i*V
# boundary internalCoeffs (component-isotropic):
#   conv ic = phi_b*vic (vic=1 zeroGradient, 0 fixedValue)
#   lap  ic = -a_b for fixedValue U (gradientInternalCoeffs=-deltaCoeffs), 0 zeroGrad
a_b = NU*bmag*bdel
ic_b = np.where(Ufixed_b, -a_b, phi_bvals)
# rAUAdj basis: D() = diag + cmptAv(ic) signed
diag0 = diag_int.copy()
np.add.at(diag0, bc_cell, ic_b)
A_u = diag0/V
rAU_u = 1.0/A_u
# mob basis (relax): D_rel = max(|diag_int + sum_b cmptMax(|ic|)|, sumOff)/alpha
sumOff = np.zeros(N)
np.add.at(sumOff, owner, np.abs(upper_f))
np.add.at(sumOff, nei, np.abs(lower_f))
Dbound = diag_int.copy()
np.add.at(Dbound, bc_cell, np.abs(ic_b))
D_rel = np.maximum(np.abs(Dbound), sumOff)/ALPHA_REL
mob = 1.0/(D_rel/V)
mobF = w*mob[owner] + (1.0-w)*mob[nei]
kf_mine = mobF*deltaCoeffs*magSf_geo
contRes = divphi.copy()
jump = U_i[nei] - U_i[owner]
downwind_is_nei = upw
log("matrix rebuilt: rAU_u avg=%.6e |contRes|max=%.3e (t=%.1fs)"
    % (rAU_u.mean(), np.abs(contRes).max(), time.time()-t0))

# ---------------- load exported M ----------------
M = scipy.io.mmread(EXPORT_JT).tocsr()
M.sort_indices()
log("M loaded nnz=%d (t=%.1fs)" % (M.nnz, time.time()-t0))

# ---------------- V1: slot-level self-checks ----------------
def m_at(rows, cols):
    return np.asarray(M[rows, cols]).ravel()
log("[V1a] kf check: M[P(o),P(n)] vs -kf ...")
v = m_at(P0+owner, P0+nei)
err = np.abs(v + kf_mine)/np.maximum(np.abs(kf_mine), 1e-300)
log("      maxRel=%.3e  medRel=%.3e  (fails>1e-9: %d/%d)"
    % (err.max(), np.median(err), (err > 1e-9).sum(), err.size))
kf_fromM = -v
# re-derive mob from the exported kf for max fidelity below
mobF_M = kf_fromM/(deltaCoeffs*magSf_geo)

log("[V1b] P-row U-slot: M[P(o),U(down,c)] = -kf*jump_c -/+ w*sf_c ...")
maxerr = 0.0; mederr = 0.0
for c in range(3):
    col_down = np.where(downwind_is_nei, 3*nei + c, 3*owner + c)
    pred = -kf_fromM*jump[:, c] - np.where(downwind_is_nei, w, -w)*sf[:, c]
    got = m_at(P0+owner, col_down)
    e = np.abs(got - pred)/np.maximum(np.abs(pred), 1e-300)
    maxerr = max(maxerr, e.max()); mederr = max(mederr, np.median(e))
log("      maxRel=%.3e medRel=%.3e" % (maxerr, mederr))

# cell -> faces adjacency (needed by V1c and block C)
cellfaces = [[] for _ in range(N)]
for fi in range(nIF):
    cellfaces[owner[fi]].append(fi)
    cellfaces[nei[fi]].append(fi)

log("[V1c] hA conv-slot: M[U(n_f,0),U(o_f,0)] = upper_f + "
    "sum_{fp of o, down_fp=o}(-upper_f/V_o*convCo*sf_fp,0*jump_fp,0) ...")
rng = np.random.default_rng(20200820)
sample = rng.choice(nIF, size=min(3000, nIF), replace=False)
maxerr = 0.0; mederr = 0.0; ncmp = 0
for f in sample:
    o_, n_ = int(owner[f]), int(nei[f])
    pred = upper_f[f]     # A^T slot: M[U(nei),U(own)] = upper (CSR line 2484)
    for fp in cellfaces[o_]:
        # down_fp == o_  <=>  (phi_fp>=0 and o_ is nei) or (phi_fp<0 and o_ is own)
        if (upw[fp] and owner[fp] != o_) or ((not upw[fp]) and owner[fp] == o_):
            coef = w[fp] if owner[fp] == o_ else (1.0-w[fp])
            jf = U_i[nei[fp], 0] - U_i[owner[fp], 0]
            pred += -upper_f[f]/V[o_]*(coef*rAU_u[o_])*sf[fp, 0]*jf
    got = m_at(np.array([3*n_+0]), np.array([3*o_+0]))[0]
    if abs(pred) > 1e-300:
        e = abs(got-pred)/abs(pred)
        maxerr = max(maxerr, e); mederr = max(mederr, e); ncmp += 1
log("      maxRel=%.3e maxRel2=%.3e over %d faces (c=0)" % (maxerr, mederr, ncmp))

# ---------------- build DeltaM = J_new^T - J_old^T ----------------
log("building DeltaM ...")
rows = []; cols = []; vals = []
aR = ALPHA_REL; oneMa = 1.0 - ALPHA_REL
# (A) kf-conv sign flip: +2kf*jump on P(o)<-U(down), -2kf*jump on P(n)<-U(down)
for c in range(3):
    col_down = np.where(downwind_is_nei, 3*nei + c, 3*owner + c)
    rows.append(P0 + owner); cols.append(col_down); vals.append(2.0*kf_fromM*jump[:, c])
    rows.append(P0 + nei);   cols.append(col_down); vals.append(-2.0*kf_fromM*jump[:, c])
# (B) NEW direct conv part: (1-a)*convTerm*w*Sf on U rows
for c in range(3):
    for c2 in range(3):
        col_down = np.where(downwind_is_nei, 3*nei + c2, 3*owner + c2)
        rows.append(3*owner + c); cols.append(col_down)
        vals.append(oneMa*w*sf[:, c]*jump[:, c2])
        rows.append(3*nei + c); cols.append(col_down)
        vals.append(oneMa*(1.0-w)*sf[:, c]*jump[:, c2])
# (C) hA-conv alpha-rescale: +(1-a)*[upper_f/V_o*convCo]*sf_fp,c*jump_fp,c2
cr = []; cc = []; cv = []
for f in range(nIF):
    o_, n_ = int(owner[f]), int(nei[f])
    for fp in cellfaces[o_]:
        coef = w[fp] if owner[fp] == o_ else (1.0-w[fp])
        cco = coef*rAU_u[o_]
        sfj = np.outer(sf[fp], jump[fp])          # (c, c2)
        down_col = 3*(nei[fp] if upw[fp] else owner[fp])
        base = oneMa*upper_f[f]/V[o_]*cco
        for c in range(3):
            for c2 in range(3):
                vv = base*sfj[c, c2]
                if vv != 0.0:
                    cr.append(3*n_ + c); cc.append(down_col + c2); cv.append(vv)
    for fp in cellfaces[n_]:
        coef = w[fp] if owner[fp] == n_ else (1.0-w[fp])
        cco = coef*rAU_u[n_]
        sfj = np.outer(sf[fp], jump[fp])
        down_col = 3*(nei[fp] if upw[fp] else owner[fp])
        base = oneMa*lower_f[f]/V[n_]*cco
        for c in range(3):
            for c2 in range(3):
                vv = base*sfj[c, c2]
                if vv != 0.0:
                    cr.append(3*o_ + c); cc.append(down_col + c2); cv.append(vv)
rows.append(np.array(cr, dtype=np.int64)); cols.append(np.array(cc, dtype=np.int64))
vals.append(np.array(cv))
# (D) Sp diagonal: -contRes per component
for c in range(3):
    rows.append(3*np.arange(N) + c); cols.append(3*np.arange(N) + c)
    vals.append(-contRes)
rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
DeltaM = sp.coo_matrix((vals, (rows, cols)), shape=(NUNK, NUNK)).tocsr()
DeltaM.sum_duplicates()
log("DeltaM built nnz=%d |DeltaM|=%.6e (t=%.1fs)"
    % (DeltaM.nnz, np.sqrt(DeltaM.power(2).sum()), time.time()-t0))

# ---------------- rhs / sources ----------------
def load_rhs_export(path):
    v = np.asarray(scipy.io.mmread(path)).ravel()
    m = v.reshape(N, 4)
    b = np.zeros(NUNK)
    b[0:NV:3] = m[:, 0]; b[1:NV:3] = m[:, 1]; b[2:NV:3] = m[:, 2]
    b[P0:P0+N] = m[:, 3]
    return b
bPD = load_rhs_export(E + "/explicitRhs_pressureDrop.mtx")
bTC = load_rhs_export(E + "/explicitRhs_thermalCoupling.mtx")
rxc = np.loadtxt(E + "/stageB6_rxc_analytic.mtx", dtype=np.float64)

# w_true (h=0.001)
def load_state(dname, tag, sign):
    Uu = np.loadtxt(STATE_CASE + "/stageB2/wstate_%s_h%s_%s_U.mtx" % (dname, tag, sign))
    pp = np.loadtxt(STATE_CASE + "/stageB2/wstate_%s_h%s_%s_p.mtx" % (dname, tag, sign))
    wst = np.zeros(NUNK)
    wst[0:NV:3] = Uu[:, 0]; wst[1:NV:3] = Uu[:, 1]; wst[2:NV:3] = Uu[:, 2]
    wst[P0:] = pp
    return wst

truth_ref = {"D1": -2.537459771850754, "D2": 0.5234094521578995,
             "D3": 6.477215547073193}

# ---------------- BL + V4: lambda-direct old vs new ----------------
res = {}
log("[BL] factorizing M (baseline) ...")
lu_old = spla.splu(M.tocsc())
log("    splu done (t=%.1fs)" % (time.time()-t0))
lam_old = lu_old.solve(bPD)
log("[V4] factorizing M_new = M + DeltaM (pinned column PREF mirrors b15's J row pin) ...")
PREF = NV
Mnew = (M + DeltaM).tolil()
Mnew[:, PREF] = 0.0
Mnew[PREF, PREF] = 1.0
Mnew = Mnew.tocsc()
Mnew.sort_indices()
lu_new = spla.splu(Mnew)
log("    splu done (t=%.1fs)" % (time.time()-t0))
lam_new = lu_new.solve(bPD)
for di, nm in enumerate(["D1", "D2", "D3"]):
    rxd = rxc[di*NUNK:(di+1)*NUNK]
    ld_old = -float(lam_old @ rxd)
    ld_new = -float(lam_new @ rxd)
    tr = truth_ref[nm]
    res[nm] = {"lam_direct_old": ld_old, "lam_direct_new": ld_new,
               "truth": tr, "ratio_old": ld_old/tr, "ratio_new": ld_new/tr}
    log("    %s: lam-direct old=%.6f (ratio %.4f) -> new=%.6f (ratio %.4f)"
        % (nm, ld_old, ld_old/tr, ld_new, ld_new/tr))

# ---------------- V2/V3: w' new vs w_true ----------------
# mirror b15: J_new = M_new^T carries the pinned ROW PREF (from the pinned
# column of M_new); solve via trans='T'.
cz = np.load("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/"
             "BFINAL-015/cycle-1/wprime_cache.npz")
wprime_old = {k: cz[k] for k in cz.files}
log("[V2] w' solves (new, via M_new^T with pinned row) ...")
compare = {}
for di, nm in enumerate(["D1", "D2", "D3"]):
    rxd = rxc[di*NUNK:(di+1)*NUNK]
    rhs = -rxd.copy(); rhs[PREF] = 0.0
    x = lu_new.solve(rhs, trans="T")     # M_new^T x = rhs
    # pin row: enforce x[PREF] consistent with identity row: x[PREF]=rhs[PREF]=0
    x[PREF] = 0.0
    w_true = load_state(nm, "0.001", "p") - load_state(nm, "0.001", "m")
    w_true /= 2.0*0.001
    mask = np.ones(NUNK, dtype=bool); mask[PREF] = False
    def mets(wp):
        d = wp - w_true
        cosv = float(wp[mask] @ w_true[mask] /
                     (np.linalg.norm(wp[mask])*np.linalg.norm(w_true[mask])))
        rel = float(np.linalg.norm(d[mask])/np.linalg.norm(w_true[mask]))
        relU = float(np.linalg.norm(d[0:NV])/np.linalg.norm(w_true[0:NV]))
        relP = float(np.linalg.norm(d[P0:])/np.linalg.norm(w_true[P0:]))
        return cosv, rel, relU, relP
    cO, rO, ruO, rpO = mets(wprime_old[nm])
    cN, rN, ruN, rpN = mets(x)
    gdpN = float(bPD @ x); gdpT = float(bPD @ w_true)
    jN = float(bTC @ x); jT = float(bTC @ w_true)
    gdpO = float(bPD @ wprime_old[nm]); jO = float(bTC @ wprime_old[nm])
    compare[nm] = {
        "old": {"cos": cO, "relL2": rO, "relU": ruO, "relP": rpO},
        "new": {"cos": cN, "relL2": rN, "relU": ruN, "relP": rpN},
        "gdp": {"wp_old": gdpO, "wp_new": gdpN, "wt": gdpT,
                "ratio_old": gdpO/gdpT, "ratio_new": gdpN/gdpT},
        "J": {"wp_old": jO, "wp_new": jN, "wt": jT,
              "ratio_old": jO/jT, "ratio_new": jN/jT},
    }
    log("    %s: OLD cos=%.6f relL2=%.4e (relU=%.3e relP=%.3e)"
        % (nm, cO, rO, ruO, rpO))
    log("    %s: NEW cos=%.6f relL2=%.4e (relU=%.3e relP=%.3e)"
        % (nm, cN, rN, ruN, rpN))
    log("    %s: gDP factor old=%.4f -> new=%.4f | J factor old=%.4f -> new=%.4f"
        % (nm, gdpO/gdpT, gdpN/gdpT, jO/jT, jN/jT))

json.dump({"lam_direct": res, "compare": compare,
           "meta": {"DeltaM_nnz": int(DeltaM.nnz),
                    "max_kf_slot_rel": float(err.max())}},
          open(OUT + "/b20_derivation_check.json", "w"), indent=1)
log("B20_CHECK_DONE t=%.1fs" % (time.time()-t0))
