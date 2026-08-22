#!/usr/bin/env python3
"""BFINAL-027 M2/M3: thermal operator segment a/b adjudication.

M2 (judge segment a): assemble thermal operator A_T offline from state-B fields
  (bounded upwind conv + laplacian + boundary coeffs, b18 part-8 template with
  defects corrected), gate (a) A_T@T - srcT - Qv ~= 0, gate (b) A_T.T@Tb - dJdT ~= 0,
  then solve A_T^T lambda_ref = b_Q (splu) and compare vs production Tb.

M3 (judge segment b): recompute thermalC field from production Tb per
  AdjHeatTransfer.H L178-328 (with mixed-BC snGrad correction) and compare vs
  B19 L2=0.00111307150907 and M1 rhs-chain anchors.

Sign conventions (all confirmed by OpenFOAM-7 source reading):
  - conv (gaussConvectionScheme.C L95-107 + LduMatrixATmul.C Amul):
      M[owner,nei] = upper = (1-w)*phi = pnT ; M[nei,owner] = lower = -w*phi = -poT
      negSumDiag: diag[owner] += poT ; diag[nei] += -pnT
  - lap (gaussLaplacianScheme.C L63-64 + negSumDiag):
      upper = deltaCoeffs*gammaMagSf > 0 ; diag -= upper on both sides
      TEqn = div - lap -> diag += kfT both, off -= kfT both  (b18 code correct)
  - bounded Sp (boundedConvectionScheme.C L68-78): diag += -sum_all_faces(phi)
    (V cancels exactly; surfaceIntegrate includes boundary faces)
  - conv boundary: outflow(phi_b>=0, zeroGradient T) -> diag += phi_b ;
      inflow(phi_b<0, fixedValue T) -> src += -phi_b*T_b
      (addBoundarySource ADDS boundaryCoeffs: boundaryCoeffs = -patchFlux*valueBoundaryCoeffs)
  - lap boundary (pGamma = |Sf|_b * (e.D.e) with e = Sf_b/|Sf_b|):
      fixedValue T: gradientInternalCoeffs=-deltaCoeffs, gBou=deltaCoeffs*T_b
        -> TEqn diag += kf_b, src += kf_b*T_b
      mixed T: gInt=-vf*deltaCoeffs, gBou=vf*deltaCoeffs*refValue+(1-vf)*refGrad
        -> TEqn diag += vf*kf_b, src += vf*kf_b*refValue (+(1-vf)*kf_b/bdel*refGrad)
      zeroGradient: nothing
    DTEffective is zeroGradient on ALL patches -> boundary D = internal value (nonzero!),
    so lap boundary contributions through T BCs ARE present (b18's T_fixed handling kept,
    mixed bottomWall handling ADDED here - b18 defect).

b18 defects fixed here:
  D1 exact pyramid centroids (b26a recipe) instead of vertex-mean
  D2 internal kfT projection uses Sf/|Sf| direction (production SfGammaSn) and
     tensor interpolated THEN projected (linear in tensor), not per-cell projection
  D3 mixed bottomWall lap boundary diag/src contribution added
  D4 gate (a) relL2 referenced vs operator scale (Qv=0 makes b18's norm(Qv) div-by-0)
  D5 paths repointed to b25_qgate

Phases:
  python b27_m23_instrument.py assemble   -> load+assemble+gates+M3, save b27_at.npz
  python b27_m23_instrument.py solve      -> load b27_at.npz, splu solve, compare, tsv
"""
import time, re, sys, json
import numpy as np
import scipy.sparse as sp

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
SB = "/home/ys/dsH/b25_qgate/stageB2"
Q1 = "/home/ys/dsH/b25_qgate/1"
E26 = "/home/ys/dsH/b8_verify_diag"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-027/cycle-1"
N = 33600
SMALL = 1e-15
Tref = 600.0
B19_C_L2 = 1.11307150907e-03

t0 = time.time()
def log(m):
    print("[%.1fs] %s" % (time.time() - t0, m), flush=True)

# ---------------- field parsing (b26a read_field, validated) ----------------
def read_field(path):
    txt = open(path).read()
    i = txt.index("internalField")
    seg_line = txt[i:txt.index("\n", i)]
    if re.match(r"internalField\s+uniform\b", seg_line):
        u = seg_line[seg_line.index("uniform") + len("uniform"):].strip().rstrip(";").strip()
        toks = []
        for tok in u.replace("(", " ").replace(")", " ").split():
            try: toks.append(float(tok))
            except ValueError: break
        val = np.array(toks) if len(toks) > 1 else float(toks[0])
        internal = None
    else:
        j = txt.index("\n(", i)
        k = txt.index("\n)", j)
        vals = []
        for line in txt[j + 1:k].splitlines():
            for tok in line.replace(";", " ").split():
                for tt in tok.replace("(", " ").replace(")", " ").split():
                    vals.append(float(tt))
        internal = np.array(vals)
        val = None
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
        elif (p + 1 < len(lines) and lines[p + 1].strip() == "{"
              and " " not in s and not s.startswith(("type", "value", "internalField"))):
            name = s; hdr = p + 1
        else:
            p += 1; continue
        if name == "boundaryField":
            p += 1; continue
        depth = 1; q = hdr + 1; bval = None; btype = None
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
                    u = t[t.index("uniform") + len("uniform"):].rstrip(";").strip()
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
            for line in lines[qq:kk + 1]:
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

def get_internal(field, Nexp):
    internal, val, _ = field
    if internal is None:
        v = val if np.isscalar(val) else np.asarray(val)
        if np.isscalar(v): return np.full(Nexp, v)
        return np.broadcast_to(v, (Nexp,) + np.asarray(v).shape)
    return internal

def dir_deriv(D6, ex, ey, ez):
    return (ex * ex * D6[:, 0] + ex * ey * D6[:, 1] + ex * ez * D6[:, 2]
            + ey * ex * D6[:, 1] + ey * ey * D6[:, 3] + ey * ez * D6[:, 4]
            + ez * ex * D6[:, 2] + ez * ey * D6[:, 4] + ez * ez * D6[:, 5])

# ============================ mesh ============================
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
        faces.append([int(v) for v in s[s.index("(") + 1:-1].split()])
def rl(p):
    t = open(p).read(); i = t.index("\n("); j = t.index("\n)", i)
    return np.array([int(x) for x in t[i + 2:j].split()], dtype=np.int64)
oa = rl(POLY + "/owner"); na = rl(POLY + "/neighbour")
nIF = len(na); nF = len(faces)
Sf = np.zeros((nF, 3))
for fi in range(nF):
    fv = pts[faces[fi]]
    Sf[fi] = 0.5 * np.cross(fv, np.roll(fv, -1, axis=0)).sum(axis=0)
magSf = np.linalg.norm(Sf, axis=1)
fCtrs = np.array([pts[fv].mean(axis=0) for fv in faces])
owner = oa[:nIF]; nei = na[:nIF]
bc = oa[nIF:]
nB = nF - nIF
bfile = open(POLY + "/boundary").read()
patch_names = re.findall(r"^\s{4}(\w+)\s*$", bfile, re.M)
patch_sizes = [int(m) for m in re.findall(r"nFaces\s+(\d+);", bfile)]
patch_of_bface = np.concatenate(
    [np.full(s, i, dtype=np.int64) for i, s in enumerate(patch_sizes)])
bface_start = {}
o = nIF
for pn, sz in zip(patch_names, patch_sizes):
    bface_start[pn] = o; o += sz
log("mesh: N=%d nIF=%d nB=%d nF=%d patches=%s sizes=%s"
    % (N, nIF, nB, nF, patch_names, patch_sizes))

# exact pyramid centroids (b26a)
cEst = np.zeros((N, 3)); nFc = np.zeros(N, dtype=np.int64)
np.add.at(cEst, oa, fCtrs); np.add.at(nFc, oa, 1)
np.add.at(cEst, na, fCtrs[:nIF]); np.add.at(nFc, na, 1)
cEst /= nFc[:, None]
num = np.zeros((N, 3)); vol3 = np.zeros(N)
SfI = Sf[:nIF]; fCtrsI = fCtrs[:nIF]
p3o = np.einsum("ij,ij->i", SfI, fCtrsI - cEst[owner])
pc_o = 0.75 * fCtrsI + 0.25 * cEst[owner]
np.add.at(num, owner, pc_o * p3o[:, None]); np.add.at(vol3, owner, p3o)
p3n = np.einsum("ij,ij->i", -SfI, fCtrsI - cEst[nei])
pc_n = 0.75 * fCtrsI + 0.25 * cEst[nei]
np.add.at(num, nei, pc_n * p3n[:, None]); np.add.at(vol3, nei, p3n)
SfB = Sf[nIF:]; fCtrsB = fCtrs[nIF:]
p3b = np.einsum("ij,ij->i", SfB, fCtrsB - cEst[bc])
pc_b = 0.75 * fCtrsB + 0.25 * cEst[bc]
np.add.at(num, bc, pc_b * p3b[:, None]); np.add.at(vol3, bc, p3b)
C = num / vol3[:, None]
V_pyr = vol3 / 3.0
V = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")[:, 0]
log("V gate (pyramid) relL2=%.3e" % (np.linalg.norm(V_pyr - V) / np.linalg.norm(V)))
dvec = C[nei] - C[owner]
magSfI = magSf[:nIF]
dc = magSfI / np.einsum("ij,ij->i", SfI, dvec)          # internal deltaCoeffs
w = np.loadtxt(MESH + "/b16mesh_weights.mtx")
bnd = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
bc_cell = bnd[:, 0].astype(np.int64); bSf = bnd[:, 1:4]
bmag = bnd[:, 4]; bdel = bnd[:, 5]

# ============================ fields ============================
T_i, _, T_b = read_field(Q1 + "/T");         T = get_internal((T_i, None, None), N)
Tb_i, _, Tb_b = read_field(Q1 + "/Tb");      Tb = get_internal((Tb_i, None, None), N)
phTh_i, _, phTh_b2 = read_field(Q1 + "/phiThermal"); phiTh_int = get_internal((phTh_i, None, None), nIF)
Q_i, Q_u, _ = read_field(Q1 + "/Q");           Qv = get_internal((Q_i, Q_u, None), N)
DT_i, _, _ = read_field(Q1 + "/DTEffective"); DT6 = get_internal((DT_i, None, None), N).reshape(N, 6)
dtdxh_i, _, _ = read_field(Q1 + "/dDTDxh");  dDTDxh6 = get_internal((dtdxh_i, None, None), N).reshape(N, 6)
phTh_b = np.zeros(nB)
T_bvals = np.zeros(nB); T_btype = {}
for pn, sz in zip(patch_names, patch_sizes):
    o0 = bface_start[pn] - nIF
    sl = slice(o0, o0 + sz)
    btp, bvp = phTh_b2[pn]
    if isinstance(bvp, float): phTh_b[sl] = bvp
    elif bvp is not None: phTh_b[sl] = bvp
    btv, bvv = T_b[pn]; T_btype[pn] = btv
    if isinstance(bvv, float): T_bvals[sl] = bvv
    elif bvv is not None: T_bvals[sl] = bvv
    else: T_bvals[sl] = T[bc_cell[sl]]
log("phiTh_b per patch:")
for pn, sz in zip(patch_names, patch_sizes):
    o0 = bface_start[pn] - nIF
    sl = slice(o0, o0 + sz)
    log("  %-14s %4d faces  min=%+.4e max=%+.4e type=%s"
        % (pn, sz, phTh_b[sl].min(), phTh_b[sl].max(), T_btype[pn]))
# mixed params
Ttxt = open(Q1 + "/T").read()
mixed_params = {}
for pn in patch_names:
    m = re.search(r"^\s{4}%s\s*\{(.*?)\}" % pn, Ttxt[Ttxt.index("boundaryField"):], re.S | re.M)
    seg = m.group(1) if m else ""
    vf = re.search(r"valueFraction\s+uniform\s+([\d.eE+-]+)", seg)
    rv = re.search(r"refValue\s+uniform\s+([\d.eE+-]+)", seg)
    rg = re.search(r"refGradient\s+uniform\s+([\d.eE+-]+)", seg)
    mixed_params[pn] = (float(vf.group(1)) if vf else None,
                        float(rv.group(1)) if rv else None,
                        float(rg.group(1)) if rg else None)
log("T patch types: %s" % {pn: T_btype[pn] for pn in patch_names})
log("T mixed params: %s" % {pn: mixed_params[pn] for pn in patch_names if T_btype[pn] == "mixed"})

dJdT = np.loadtxt(SB + "/wstate_dJdT.mtx")
Mflow, Tref0 = np.loadtxt(SB + "/wstate_outletMeta.mtx")
log("dJdT: nz=%d maxabs=%.6e L2=%.6e  M_frozen=%.6e Tref=%.1f"
    % (int((dJdT != 0).sum()), np.abs(dJdT).max(), np.linalg.norm(dJdT), Mflow, Tref0))

# ============================ A_T assembly ============================
log("assembling A_T ...")
# internal face directions (production: e = Sf/|Sf| for gammaMagSf = SfGammaSn)
ei = SfI / magSfI[:, None]
# tensor interpolated THEN projected (production: interpolate(DTEffective) -> Sf & gamma -> Sn)
Dsurf = w[:, None] * DT6[owner] + (1.0 - w)[:, None] * DT6[nei]
kfT = magSfI * dir_deriv(Dsurf, ei[:, 0], ei[:, 1], ei[:, 2]) * dc
poT = np.where(phiTh_int >= 0.0, phiTh_int, 0.0)
pnT = np.where(phiTh_int >= 0.0, 0.0, phiTh_int)   # <= 0
diagT = np.zeros(N)
np.add.at(diagT, owner, poT)
np.add.at(diagT, nei, -pnT)
A_off_own = pnT.copy()    # A[owner,nei]
A_off_nei = -poT.copy()   # A[nei,owner]
np.add.at(diagT, owner, kfT)
np.add.at(diagT, nei, kfT)
A_off_own -= kfT
A_off_nei -= kfT
# bounded Sp: diag += -sum_all_faces(phi) (internal + boundary)
divphiT = np.zeros(N)
np.add.at(divphiT, owner, phiTh_int)
np.add.at(divphiT, nei, -phiTh_int)
np.add.at(divphiT, bc_cell, phTh_b)
diagT += -divphiT
# conv boundary: outflow diag += phi_b ; inflow src += -phi_b*T_b
obT = phTh_b >= 0.0
np.add.at(diagT, bc_cell[obT], phTh_b[obT])
srcT = np.zeros(N)
np.add.at(srcT, bc_cell[~obT], -phTh_b[~obT] * T_bvals[~obT])
# lap boundary: e = Sf_b/|Sf_b| (production SfGammaSn boundary)
eb = SfB / np.maximum(bmag, SMALL)[:, None]
DTb6 = DT6[bc_cell]
pGammaB = bmag * dir_deriv(DTb6, eb[:, 0], eb[:, 1], eb[:, 2])
kfTb = pGammaB * bdel
for pn, sz in zip(patch_names, patch_sizes):
    o0 = bface_start[pn] - nIF
    sl = slice(o0, o0 + sz)
    t_ = T_btype[pn]
    if t_ == "fixedValue":
        np.add.at(diagT, bc_cell[sl], kfTb[sl])
        np.add.at(srcT, bc_cell[sl], kfTb[sl] * T_bvals[sl])
    elif t_ == "mixed":
        vf, rv, rg = mixed_params[pn]
        vf = vf if vf is not None else 0.0
        rv = rv if rv is not None else 0.0
        rg = rg if rg is not None else 0.0
        np.add.at(diagT, bc_cell[sl], vf * kfTb[sl])
        np.add.at(srcT, bc_cell[sl], vf * kfTb[sl] * rv + (kfTb[sl] / bdel[sl]) * (1.0 - vf) * rg)
rows = np.concatenate([owner, nei, np.arange(N)])
cols = np.concatenate([nei, owner, np.arange(N)])
vals = np.concatenate([A_off_own, A_off_nei, diagT])
AT = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
log("A_T assembled: nnz=%d" % AT.nnz)

# ============================ gates ============================
# gate (a): A_T@T - srcT - Qv = 0   (Qv = 0 in state B)
resT = AT @ T - srcT - Qv
scale_a = np.linalg.norm(AT @ T) + np.linalg.norm(srcT)
rel_a = np.linalg.norm(resT) / max(scale_a, SMALL)
log("GATE (a) A_T*T - src - Q: maxabs=%.4e relL2(vs op scale)=%.3e  |Q|max=%.3e"
    % (np.abs(resT).max(), rel_a, np.abs(Qv).max()))
# gate (b): A_T^T Tb - dJdT = 0
resTb = (AT.T @ Tb) - dJdT
rel_b = np.linalg.norm(resTb) / max(np.linalg.norm(dJdT), SMALL)
log("GATE (b) A_T^T*Tb - dJdT: maxabs=%.4e relL2=%.3e" % (np.abs(resTb).max(), rel_b))

# ============================ M3: C_field ============================
log("M3: recompute thermalC from production Tb ...")
dvecC = C[nei] - C[owner]
efc = dvecC / np.maximum(np.linalg.norm(dvecC, axis=1), 1e-30)[:, None]
D6o = dDTDxh6[owner]; D6n = dDTDxh6[nei]
dd_o = dir_deriv(D6o, efc[:, 0], efc[:, 1], efc[:, 2])
dd_n = dir_deriv(D6n, efc[:, 0], efc[:, 1], efc[:, 2])
base = (Tb[owner] - Tb[nei]) * (T[owner] - T[nei]) * dc * magSfI
Cfield = np.zeros(N)
np.add.at(Cfield, owner, -w * dd_o * base)
np.add.at(Cfield, nei, -(1.0 - w) * dd_n * base)
# boundary: patch-aware snGrad (production T.boundaryField()[patchi].snGrad())
bDvec = fCtrsB - C[bc_cell]
bef = bDvec / np.maximum(np.linalg.norm(bDvec, axis=1), 1e-30)[:, None]
dd_b = dir_deriv(dDTDxh6[bc_cell], bef[:, 0], bef[:, 1], bef[:, 2])
snb = np.zeros(nB)
for pn, sz in zip(patch_names, patch_sizes):
    o0 = bface_start[pn] - nIF
    sl = slice(o0, o0 + sz)
    t_ = T_btype[pn]
    if t_ == "fixedValue":
        snb[sl] = (T_bvals[sl] - T[bc_cell[sl]]) * bdel[sl]
    elif t_ == "mixed":
        vf, rv, rg = mixed_params[pn]
        vf = vf if vf is not None else 0.0
        rv = rv if rv is not None else 0.0
        rg = rg if rg is not None else 0.0
        snb[sl] = vf * (rv - T[bc_cell[sl]]) * bdel[sl] + (1.0 - vf) * rg
    else:
        snb[sl] = 0.0
nbz = np.abs(snb) > SMALL
np.add.at(Cfield, bc_cell[nbz], dd_b[nbz] * Tb[bc_cell[nbz]] * snb[nbz] * bmag[nbz])
C_L2 = np.linalg.norm(Cfield)
log("M3 C_field L2=%.8e (B19 anchor %.8e) ratio=%.4f"
    % (C_L2, B19_C_L2, C_L2 / B19_C_L2))
z = np.loadtxt(E26 + "/stageB6_rxc_z_analytic.mtx").reshape(3, N)
projC = [float(np.dot(Cfield, zk)) for zk in z]

# M1 anchors from b26a_m1_identity.tsv (state-B measured)
M1 = {"D1": dict(C=-3.0341656277e-02, Gx=-1.5704339057e-03, FD=4.2524696236e-02, rhs=7.4436786419e-02),
      "D2": dict(C=-8.2108899666e-03, Gx=4.0431596817e-03, FD=-3.5289725180e-03, rhs=6.3875776689e-04),
      "D3": dict(C=-1.9791957452e-02, Gx=8.4067076831e-03, FD=-4.4495118952e-04, rhs=1.0940298579e-02)}

# ============================ save + outputs ============================
np.savez_compressed(OUT + "/b27_at.npz",
                    AT_data=AT.data, AT_indices=AT.indices, AT_indptr=AT.indptr,
                    srcT=srcT, dJdT=dJdT, Tb=Tb, T=T, Qv=Qv,
                    bc_cell=bc_cell, bface_start_keys=list(bface_start.keys()),
                    bface_start_vals=np.array([bface_start[k] for k in bface_start]),
                    patch_sizes=np.array(patch_sizes, dtype=np.int64),
                    Cfield=Cfield)
log("saved %s/b27_at.npz" % OUT)

gate_rows = [dict(gate="a", maxabs=float(np.abs(resT).max()), relL2=float(rel_a),
                  scale=float(scale_a), n_ref="op-scale"),
             dict(gate="b", maxabs=float(np.abs(resTb).max()), relL2=float(rel_b),
                  scale=float(np.linalg.norm(dJdT)), n_ref="dJdT")]
with open(OUT + "/b27_gates.tsv", "w") as f:
    f.write("gate\tmaxabs\trelL2\tscale\tnorm_ref\n")
    for r in gate_rows:
        f.write("%s\t%.6e\t%.6e\t%.6e\t%s\n" % (r["gate"], r["maxabs"], r["relL2"], r["scale"], r["n_ref"]))
log("wrote %s/b27_gates.tsv" % OUT)

m3_rows = [dict(qty="C_L2", val=C_L2, ref=B19_C_L2, ratio=C_L2 / B19_C_L2)]
for d, zi, v in zip(("D1", "D2", "D3"), z, projC):
    m3_rows.append(dict(qty="proj_%s" % d, val=v, ref=M1[d]["C"], ratio=v / M1[d]["C"]))
with open(OUT + "/b27_m3.tsv", "w") as f:
    f.write("qty\tval\tref\tratio\n")
    for r in m3_rows:
        f.write("%s\t%.8e\t%.8e\t%.6f\n" % (r["qty"], r["val"], r["ref"], r["ratio"]))
log("wrote %s/b27_m3.tsv" % OUT)

# ============================ M2: solve + compare ============================
log("M2: splu solve A_T^T lambda_ref = b_Q (=dJdT) ...")
from scipy.sparse.linalg import splu
ATcsc = AT.T.tocsc()
t1 = time.time()
lu = splu(ATcsc)
log("splu factorized (t=%.1fs)" % (time.time() - t1))
lam_ref = lu.solve(dJdT)
t2 = time.time()
log("lambda_ref solved (t=%.1fs)" % (time.time() - t2))
res_solve = AT.T @ lam_ref - dJdT
log("solve residual |AT.T*lam_ref - b_Q|/|b_Q| = %.3e"
    % (np.linalg.norm(res_solve) / max(np.linalg.norm(dJdT), 1e-30)))

# comparison vs production Tb
relL2 = np.linalg.norm(lam_ref - Tb) / np.linalg.norm(Tb)
cosT = float(np.dot(lam_ref, Tb) / (np.linalg.norm(lam_ref) * np.linalg.norm(Tb)))
sign_flips = int(((lam_ref > 0) != (Tb > 0)).sum())
log("M2 vs Tb: relL2=%.6e cos=%.8f sign_flips=%d (of %d)"
    % (relL2, cosT, sign_flips, N))
# per-patch boundary bands
bface_start2 = {k: int(v) for k, v in bface_start.items()}
bc_all = np.sort(np.unique(bc_cell))
internal_cells = np.setdiff1d(np.arange(N), bc_all)
rows2 = [dict(region="INTERNAL", n=int(len(internal_cells)),
              relL2=float(np.linalg.norm(lam_ref[internal_cells] - Tb[internal_cells])
                          / max(np.linalg.norm(Tb[internal_cells]), 1e-30)),
              maxabs=float(np.abs(lam_ref[internal_cells] - Tb[internal_cells]).max()))]
for pn, sz in zip(patch_names, patch_sizes):
    o0 = bface_start2[pn] - nIF
    sl = slice(o0, o0 + sz)
    cells = bc_cell[sl]
    rows2.append(dict(region=pn, n=int(sz),
                      relL2=float(np.linalg.norm(lam_ref[cells] - Tb[cells])
                                  / max(np.linalg.norm(Tb[cells]), 1e-30)),
                      maxabs=float(np.abs(lam_ref[cells] - Tb[cells]).max())))
with open(OUT + "/b27_m2_compare.tsv", "w") as f:
    f.write("region\tn\trelL2\tmaxabs\n")
    for r in rows2:
        f.write("%s\t%d\t%.6e\t%.6e\n" % (r["region"], r["n"], r["relL2"], r["maxabs"]))
log("wrote %s/b27_m2_compare.tsv" % OUT)

np.savez_compressed(OUT + "/b27_lamref.npz", lam_ref=lam_ref, Tb=Tb)
log("done (t=%.1fs)" % (time.time() - t0))
