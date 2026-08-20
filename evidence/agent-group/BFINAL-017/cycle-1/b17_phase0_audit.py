#!/usr/bin/env python3
"""
BFINAL-017 phase-0: assembly-contradiction audit (AD1 + AD2).  AUDIT ONLY.

AD1  chain adjoint point test on filter_chainrule.H's operator L:
       <L(g), d>  ==  <g, L^T(d)>
     with L^T implemented INDEPENDENTLY in numpy (transpose of the production
     adjoint formulas) and z implemented INDEPENDENTLY as the forward tangent
     of the production map filter_x.H + diff.c (implicit adaptive-eta).
     Plus linearity L(2g)==2L(g) (adaptive-eta rank-one check).
AD2  per-cell / per-term comparison of the production R_x contraction
     (sensitivity.H gsenshMomentum/-rxPressureRowT*dAlphaDxh) against the
     stageB6 oracle rxc (anRUa/assembleWeightedRPa), splitting the observed
     assembly-vs-lambda-direct mismatch into  chain  +  contraction  parts.

All production-side numbers come from exported data only (no production code
is reused on the transpose side).
"""
import json, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

MESH = "/home/ys/dsH/b16_mesh"
EXPORT_CASE = "/home/ys/dsH/b15_export"
STATES = "/home/ys/dsH/b16_states/1"
POLY = "/home/ys/dsH/b16_mesh/constant/polyMesh"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-017/cycle-1"

N = 33600
NV = 3 * N
NUNK = 4 * N
DEL = 8.0                      # validationProjectionBeta (mmaUpdateEnabled=false)
ETA_TOL = 1e-10                # projectionEtaTolerance
FILTER_R = 3.0
DESIGN_VOL = 5040 * 1.25e-10   # domainIntegrate(designMask); 5040 active cells
BFINAL012_ADJ = {"D1": -5.4050587339, "D2": 1.10592698054, "D3": 14.0447142088}
BFINAL016_LAMBDADIRECT = {"D1": -4.485902855631131, "D2": 0.44787866660803377,
                          "D3": 14.730268781930677}

t0 = time.time()
def log(m): print(m, flush=True)

res = {"_meta": {
    "N": N, "DEL": DEL, "FILTER_R": FILTER_R,
    "b15": EXPORT_CASE, "b16_states": STATES}}

# ----------------------------------------------------------------------
# loaders
# ----------------------------------------------------------------------
def read_of_scalar(path):
    vals = []; intxt = False
    for line in open(path):
        s = line.strip()
        if not intxt:
            if s == "(":
                intxt = True
            continue
        if s == ")":
            break
        for tok in s.replace(";", " ").split():
            vals.append(float(tok))
    return np.array(vals)

def read_of_vector(path):
    txt = open(path).read()
    i = txt.index("\n("); j = txt.index("\n)", i)
    return np.fromstring(txt[i+2:j].replace("(", " ").replace(")", " "),
                         sep=" ").reshape(-1, 3)

def read_polymesh_faces(path):
    txt = open(path).read()
    i = txt.index("\n(")
    j = txt.index("\n)", i)
    faces = []
    for line in txt[i+1:j].splitlines():
        s = line.strip().rstrip(")").split("(")
        if len(s) == 2 and s[0].strip().isdigit():
            faces.append([int(v) for v in s[1].split()])
    return faces

log("=== BFINAL-017 phase-0 audit ===")

# ----------------------------------------------------------------------
# 0. mesh: Sf / deltaCoeffs from polyMesh (independent reconstruction)
# ----------------------------------------------------------------------
pts = read_of_vector(POLY + "/points")
faces = read_polymesh_faces(POLY + "/faces")
def read_labels(path):
    txt = open(path).read()
    i = txt.index("\n(")
    j = txt.index("\n)", i)
    return np.array([int(t) for t in txt[i+2:j].split()], dtype=np.int64)
owner_all = read_labels(POLY + "/owner")
nei_all = read_labels(POLY + "/neighbour")
owner = owner_all[:len(nei_all)]
nei = nei_all.copy()
nIF = len(nei_all)
log("polymesh: points=%d faces=%d internalFaces=%d (t=%.1fs)"
    % (len(pts), len(faces), nIF, time.time()-t0))

# Sf = 0.5 * sum_p cross(p_i, p_{i+1}) over the face loop (OpenFOAM convention)
Sf = np.zeros((nIF, 3))
for fi in range(nIF):
    fv = pts[faces[fi]]
    Sf[fi] = 0.5*np.cross(fv, np.roll(fv, -1, axis=0)).sum(axis=0)
Sf_ref = np.loadtxt(MESH + "/b16mesh_sf.mtx")
log("[mesh] Sf reproduction vs b16mesh_sf: maxabs=%.3e"
    % np.max(np.abs(Sf - Sf_ref)))

# cell centres: average of the cell's unique vertices (exact for hexes)
nvtx = np.zeros(N); csum = np.zeros((N, 3))
seen = [np.zeros(0, dtype=np.int64)]*0
# accumulate with dedupe via per-cell sets is slow in pure python; use the
# face-vertex incidence: each hex cell has exactly 8 distinct vertices over
# its 6 faces.  Use a dict-of-sets once.
cellverts = {}
for fi in range(len(faces)):
    fv = faces[fi]
    cells = (owner[fi], nei[fi]) if fi < nIF else (owner_all[fi],)
    for c in cells:
        s = cellverts.get(c)
        if s is None:
            s = set(); cellverts[c] = s
        s.update(fv)
Ccen = np.zeros((N, 3))
for c, vs in cellverts.items():
    Ccen[c] = pts[sorted(vs)].mean(axis=0)
log("[mesh] cell centres built (%d cells, t=%.1fs)" % (len(cellverts), time.time()-t0))

# orthogonality + deltaCoeffs (=|Sf|/(Sf&d))
d_vec = Ccen[nei] - Ccen[owner]
Sfdot = np.einsum("ij,ij->i", Sf, d_vec)
magSf = np.linalg.norm(Sf, axis=1)
cosorth = Sfdot/(magSf*np.linalg.norm(d_vec, axis=1))
log("[mesh] orthogonality cos(Sf,d): min=%.15f max=%.15f"
    % (cosorth.min(), cosorth.max()))
deltaCoeffs = magSf/np.maximum(Sfdot, 1e-300)

# ----------------------------------------------------------------------
# 1. filter matrix A = laplacian(designFilterFaceMask,.) - Sp(b,.)
#    (zeroGradient BCs -> no boundary terms)
# ----------------------------------------------------------------------
maskS = read_of_scalar(STATES + "/designFilterFaceMask")   # 96860 internal
assert len(maskS) == nIF
len_c = (DESIGN_VOL/5040.0)**(1.0/3.0)
bcoef = 1.0/max((FILTER_R*len_c/3.464)**2, 1e-15)
V = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")[:, 0]
log("[filter] len=%.17g b=%.17g  V sum=%.17g" % (len_c, bcoef, V.sum()))

af = maskS*magSf*deltaCoeffs
rows = np.concatenate([owner, nei, owner, nei])
cols = np.concatenate([owner, nei, nei, owner])
vals = np.concatenate([-af, -af, af, af])
A = sp.coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsr()
A = A - sp.diags(bcoef*V)
log("[filter] A built: nnz=%d (t=%.1fs)" % (A.nnz, time.time()-t0))

# --- validation pair 1: forward  A @ xp == -b*x
x = read_of_scalar(STATES + "/x")
xp = read_of_scalar(STATES + "/xp")
r1 = A@xp + bcoef*V*x
log("[filter][VALIDATE fwd ] |A@xp + b*V*x|=%.3e rel=%.3e"
    % (np.linalg.norm(r1), np.linalg.norm(r1)/np.linalg.norm(bcoef*V*x)))
# --- validation pair 2: adjoint  A @ gsensPD == -b*gsenshPD (written fields)
gsensPD16 = read_of_scalar(STATES + "/gsensPressureDrop")
gsenshPD16 = read_of_scalar(STATES + "/gsenshPressureDrop")  # post-eta (in-place)
r2 = A@gsensPD16 + bcoef*V*gsenshPD16
log("[filter][VALIDATE adj ] |A@gsens + b*V*gsensh|=%.3e rel=%.3e"
    % (np.linalg.norm(r2), np.linalg.norm(r2)/np.linalg.norm(bcoef*V*gsenshPD16)))

lu = spla.splu(A.tocsc())
def Filt(u):        # production filter solve: A s = -b*V*u (V validated on the xp/x pair)
    return lu.solve(-bcoef*V*u)

# ----------------------------------------------------------------------
# 2. projection data: eta5 (bisection replica of filter_x.H/diff.c),
#    drho / dProjectionDeta  (replica of filter_x.H L62-86)
# ----------------------------------------------------------------------
mask = read_of_scalar(STATES + "/designMask").astype(np.float64)
maskb = mask > 0.5

def heaviside(g, eta):
    e = np.exp(-DEL)
    lo = g <= eta
    out = np.empty_like(g)
    gl = g[lo]/eta
    out[lo] = eta*(np.exp(-DEL*(1.0-gl)) - (1.0-gl)*e)
    gh = (g[~lo] - eta)/(1.0-eta)
    out[~lo] = eta + (1.0-eta)*(1.0 - np.exp(-DEL*gh) + gh*e)
    return out

def diffvol(g, eta):            # diff.c: sum_mask (xp - xh)*V
    return np.sum((g[maskb] - heaviside(g[maskb], eta))*V[maskb])

e0, e1 = 1e-4, 0.9999
y0, y1 = diffvol(xp, e0), diffvol(xp, e1)
assert y0*y1 <= 0.0, "eta bracket failed"
while (e1 - e0) > ETA_TOL:
    e5 = 0.5*(e0 + e1)
    y5 = diffvol(xp, e5)
    if y0*y5 < 0.0:
        e1, y1 = e5, y5
    else:
        e0, y0 = e5, y5
eta5 = 0.5*(e0 + e1)
log("[proj] eta5 mine=%.15f  (b15 log: 0.748005161603)" % eta5)

idx = np.where(maskb)[0]
xp_m = xp[idx]
e = np.exp(-DEL)
drho = np.zeros(N); dPdeta = np.zeros(N)
lo = xp_m <= eta5
pe = np.exp(-DEL*(1.0 - xp_m[lo]/eta5))
drho[idx[lo]] = DEL*pe + e
dPdeta[idx[lo]] = pe*(1.0 - DEL*xp_m[lo]/eta5) - e
pe2 = np.exp(-DEL*(xp_m[~lo] - eta5)/(1.0 - eta5))
drho[idx[~lo]] = DEL*pe2 + e
dPdeta[idx[~lo]] = pe2*(1.0 - DEL*(1.0 - xp_m[~lo])/(1.0 - eta5)) - e
Dden = np.sum(mask*dPdeta*V)
log("[proj] projectionEtaDenominator mine=%.12e (b15 log: -3.1385996818e-06)"
    % Dden)
drho_ref = np.loadtxt(EXPORT_CASE + "/stageB6_proj.mtx")
log("[proj] drho vs stageB6_proj.mtx: maxabs=%.3e"
    % np.max(np.abs(drho - drho_ref)))

# ----------------------------------------------------------------------
# 3. the chain operators (independent numpy implementations)
# ----------------------------------------------------------------------
# E(g):   production eta step  (filter_chainrule.H L74-161)
# E_T(v): INDEPENDENT transpose of E (diagonal + rank-one transposed)
# z(d):   INDEPENDENT forward tangent of filter_x.H + diff.c:
#         y = -b A^-1 (mask d);  deta = <mask V (1-drho) y>/D  from
#         d/dxp[ sum_mask (xp - xh(xp,eta))V ] = 0  (volume preservation)
def E(g):
    u = mask*g
    C = np.sum(mask*dPdeta*u)/Dden
    return drho*u + mask*V*(1.0-drho)*C

def E_T(v):
    return drho*mask*v + mask*dPdeta*np.sum(mask*V*(1.0-drho)*v)/Dden

def L(g):            # production chain: mask -> E -> filter solve -> mask
    return mask*Filt(E(g))

def L_T(d):          # independent transpose: mask -> filter -> E^T
    return E_T(Filt(mask*d))

def z_tangent(d):    # independent forward tangent (dxh/dx * d)
    y = Filt(mask*d)
    return drho*y + dPdeta*np.sum(mask*V*(1.0-drho)*y)/Dden

def filt_tangent(d): # dxp/dx * d  (forward filter tangent; A symmetric)
    return Filt(d)

# ----------------------------------------------------------------------
# 4. load the exported gauge data
# ----------------------------------------------------------------------
D_all = np.loadtxt(EXPORT_CASE + "/stageB6_dirs.mtx").reshape(3, N)
g15_mom = np.loadtxt(EXPORT_CASE + "/rxpr_prod_gsensh_momentum.mtx")
g15_pr = np.loadtxt(EXPORT_CASE + "/rxpr_prod_gsensh_pressurerow.mtx")
g15_tot = np.loadtxt(EXPORT_CASE + "/rxpr_prod_gsensh_total.mtx")
gsens15 = np.loadtxt(EXPORT_CASE + "/rxpr_prod_gsens.mtx")
z15 = np.loadtxt(EXPORT_CASE + "/stageB6_rxc_z_analytic.mtx").reshape(3, N)
zfd = np.loadtxt(EXPORT_CASE + "/stageB6_rxc_xh_FD_all_eps.mtx").reshape(3, 5, N)
rxc = np.loadtxt(EXPORT_CASE + "/stageB6_rxc_analytic.mtx")
lam15 = np.loadtxt(EXPORT_CASE + "/stageB6_lambda.mtx")
Uc16 = read_of_vector(STATES + "/Uc"); pc16 = read_of_scalar(STATES + "/pc")
lam16 = np.zeros(NUNK)
lam16[0:NV:3] = Uc16[:, 0]; lam16[1:NV:3] = Uc16[:, 1]; lam16[2:NV:3] = Uc16[:, 2]
lam16[NV:] = pc16
U16 = read_of_vector(STATES + "/U")
daDxh = read_of_scalar(STATES + "/dAlphaDxh")
log("[data] gauge data loaded (t=%.1fs)" % (time.time()-t0))

# formula cross-checks against b15 raw exports (garbage-lambda but same state)
g15_mom_mine = -daDxh*np.einsum("ij,ij->i", U16, read_of_vector(EXPORT_CASE + "/1/Uc"))*V
log("[xchk] -dAlphaDxh(U&Uc15)V vs rxpr_prod_gsensh_momentum: maxabs=%.3e"
    % np.max(np.abs(g15_mom_mine - g15_mom)))

# ----------------------------------------------------------------------
# AD1 -- chain adjoint point tests
# ----------------------------------------------------------------------
log("--- AD1: chain adjoint point tests ---")
rng = np.random.default_rng(20260820)
g_rand = rng.standard_normal(N)
d_rand = rng.standard_normal(N)

# (a) production pair reproduction: L(g15_tot) vs exported rxpr_prod_gsens
lg = L(g15_tot)
dn = np.linalg.norm(lg - gsens15)
log("AD1a L(g15_tot) vs rxpr_prod_gsens (production chain pair): "
    "|diff|=%.3e |gsens|=%.3e rel=%.3e maxabs=%.3e"
    % (dn, np.linalg.norm(gsens15), dn/np.linalg.norm(gsens15),
       np.max(np.abs(lg-gsens15))))
# also reproduce with the b15 pair through E: E(g15_tot) vs written? (b16 path below)

# (b) b16 written pair: mask*Filt(gsenshPD16) vs gsensPD16
lg2 = mask*Filt(gsenshPD16)
dn2 = np.linalg.norm(lg2 - gsensPD16)
log("AD1b filter(gsenshPD16_written) vs gsensPD16_written: rel=%.3e"
    % (dn2/np.linalg.norm(gsensPD16)))

# (c) adjoint closures <L(g),d> == <g,L^T(d)>
pairs = [("g15_tot", g15_tot, "D1", D_all[0]), ("g15_tot", g15_tot, "D2", D_all[1]),
         ("g15_tot", g15_tot, "D3", D_all[2]),
         ("g_rand", g_rand, "d_rand", d_rand),
         ("g_rand", g_rand, "D2", D_all[1])]
ad1 = {}
for gn, g, dn_, d in pairs:
    lhs = float(np.dot(L(g), d)); rhs = float(np.dot(g, L_T(d)))
    sc = max(abs(lhs), abs(rhs))
    rel = abs(lhs-rhs)/sc if sc > 0 else 0.0
    ad1["%s/%s" % (gn, dn_)] = {"lhs": lhs, "rhs": rhs, "rel": rel}
    log("AD1c <L(%s),%s>=%.12e   <g,L^T(%s)>=%.12e   rel=%.3e"
        % (gn, dn_, lhs, gn, rhs, rel))

# (d) linearity: L(2g) == 2 L(g)
for gn, g in [("g15_tot", g15_tot), ("g_rand", g_rand)]:
    lin = np.linalg.norm(L(2.0*g) - 2.0*L(g))
    ad1["linearity_"+gn] = float(lin)
    log("AD1d |L(2%s)-2L(%s)|=%.3e" % (gn, gn, lin))

# (e) z consistency: my independent tangent vs stageB6 z_analytic vs FD
ad1z = {}
for k, nm in enumerate(["D1", "D2", "D3"]):
    zm = z_tangent(D_all[k])
    ad1z[nm] = {
        "vs_z_analytic_rel": float(np.linalg.norm(zm-z15[k])/np.linalg.norm(z15[k])),
        "vs_zFD_eps1e-5_rel": float(np.linalg.norm(zm-zfd[k,4])/np.linalg.norm(zfd[k,4])),
        "zFD_eps-scan_rel": [float(np.linalg.norm(z15[k]-zfd[k,i])/np.linalg.norm(zfd[k,i]))
                             for i in range(5)]}
    log("AD1e z(%s): vs z_analytic rel=%.3e ; vs zFD(eps=1e-5) rel=%.3e ; "
        "z_analytic vs zFD eps=1e-3..1e-5 rel=%s"
        % (nm, ad1z[nm]["vs_z_analytic_rel"], ad1z[nm]["vs_zFD_eps1e-5_rel"],
           ["%.2e" % v for v in ad1z[nm]["zFD_eps-scan_rel"]]))
res["AD1"] = {"closures": ad1, "z": ad1z}

# ----------------------------------------------------------------------
# AD2 -- R_x contraction path comparison (good-lambda b16 side)
# ----------------------------------------------------------------------
log("--- AD2: R_x contraction path comparison ---")
# invert the eta step on the written b16 source field -> g_total_16
w = gsenshPD16.copy()
assert np.max(np.abs(w[~maskb])) == 0.0, "written gsensh nonzero outside mask"
rw = w[maskb]/drho[maskb]
num = np.sum((mask*dPdeta)[maskb]*rw)
den = Dden + np.sum((mask*dPdeta*V*(1.0-drho))[maskb]/drho[maskb])
Cinv = num/den
g_tot16 = np.zeros(N)
g_tot16[maskb] = rw - (V*(1.0-drho)*Cinv)[maskb]/drho[maskb]
rt = np.linalg.norm(E(g_tot16) - w)
log("AD2a eta-inversion round-trip |E(g_tot16)-written|=%.3e (|written|=%.3e)"
    % (rt, np.linalg.norm(w)))

# momentum term from the formula (production sensitivity.H:89-90), lambda16
g_mom16 = -daDxh*np.einsum("ij,ij->i", U16, Uc16)*V
g_pr16 = g_tot16 - g_mom16
log("AD2b |g_tot16|=%.6e  |g_mom16|=%.6e  |g_pr16|=%.6e"
    % (np.linalg.norm(g_tot16), np.linalg.norm(g_mom16), np.linalg.norm(g_pr16)))

ad2 = {}
for k, nm in enumerate(["D1", "D2", "D3"]):
    d = D_all[k]; z = z15[k]
    rxd = rxc[k*NUNK:(k+1)*NUNK]
    anRUd = rxd[:NV].reshape(N, 3)   # 3N cell-major
    anRPd = rxd[NV:]
    assembly = float(np.dot(gsensPD16, d))                  # production chain
    gz = float(np.dot(g_tot16, z))                          # contraction @ z
    ldirect = -float(np.dot(lam16, rxd))                    # -lambda^T rxd
    mom_pr = float(np.dot(g_mom16, z))
    mom_or = -float(np.sum(Uc16*anRUd))
    pr_pr = float(np.dot(g_pr16, z))
    pr_or = -float(np.dot(pc16, anRPd))
    ad2[nm] = {
        "assembly_b16field": assembly, "B12_ADJ_ref": BFINAL012_ADJ[nm],
        "contraction_at_z": gz, "lambda_direct": ldirect,
        "B16_lambda_direct_ref": BFINAL016_LAMBDADIRECT[nm],
        "delta_chain": assembly - gz, "delta_contraction": gz - ldirect,
        "mom_production": mom_pr, "mom_oracle": mom_or,
        "pr_production": pr_pr, "pr_oracle": pr_or,
        "delta_total": assembly - ldirect}
    log("AD2 %s: assembly=%.10f (B12 %.10f)\n"
        "      <g_tot,z>=%.10f   -lam^T rxd=%.10f (B16 %.10f)\n"
        "      delta_chain=%.3e   delta_contraction=%.3e   delta_total=%.3e\n"
        "      momentum: prod=%.10f oracle=%.10f (diff %.3e)\n"
        "      pressure: prod=%.10f oracle=%.10f (diff %.3e)"
        % (nm, assembly, BFINAL012_ADJ[nm], gz, ldirect,
           BFINAL016_LAMBDADIRECT[nm], assembly-gz, gz-ldirect, assembly-ldirect,
           mom_pr, mom_or, mom_pr-mom_or, pr_pr, pr_or, pr_pr-pr_or))

res["AD2"] = ad2

# per-cell localisation for the dominant differing direction
loc = {}
for k, nm in enumerate(["D1", "D2", "D3"]):
    z = z15[k]
    rxd = rxc[k*NUNK:(k+1)*NUNK]
    anRUd = rxd[:NV].reshape(N, 3)
    anRPd = rxd[NV:]
    cm_pr = g_mom16*z                       # production per-cell (momentum)
    cm_or = -np.sum(Uc16*anRUd, axis=1)     # oracle per-cell (momentum)
    cp_pr = g_pr16*z
    cp_or = -pc16*anRPd
    dm = cm_pr - cm_or; dp = cp_pr - cp_or
    loc[nm] = {
        "mom_maxabs": float(np.max(np.abs(dm))),
        "mom_l2": float(np.linalg.norm(dm)),
        "pr_maxabs": float(np.max(np.abs(dp))),
        "pr_l2": float(np.linalg.norm(dp)),
        "mom_argmax": int(np.argmax(np.abs(dm))),
        "pr_argmax": int(np.argmax(np.abs(dp))),
        "mom_top10_abs": float(np.sort(np.abs(dm))[::-1][:10].sum()),
        "pr_top10_abs": float(np.sort(np.abs(dp))[::-1][:10].sum())}
    log("AD2loc %s: |mom diff| L2=%.3e max=%.3e@%d |pr diff| L2=%.3e max=%.3e@%d"
        % (nm, loc[nm]["mom_l2"], loc[nm]["mom_maxabs"], loc[nm]["mom_argmax"],
           loc[nm]["pr_l2"], loc[nm]["pr_maxabs"], loc[nm]["pr_argmax"]))
res["AD2_localisation"] = loc

# ----------------------------------------------------------------------
# cross-identity with the b15 garbage lambda (algebra check only)
# ----------------------------------------------------------------------
lam_id = {}
for k, nm in enumerate(["D1", "D2", "D3"]):
    gz15 = float(np.dot(g15_tot, z15[k]))
    ld15 = -float(np.dot(lam15, rxc[k*NUNK:(k+1)*NUNK]))
    lam_id[nm] = {"g15_at_z": gz15, "lambda_direct_15": ld15,
                  "rel": abs(gz15-ld15)/max(abs(ld15), 1e-300)}
    log("XID %s (b15 in-run lambda): <g15_tot,z>=%.10f  -lam15^T rxd=%.10f  rel=%.3e"
        % (nm, gz15, ld15, lam_id[nm]["rel"]))
res["cross_identity_b15"] = lam_id

# z without eta-response (rank-one OFF) for context
noeta = {}
for k, nm in enumerate(["D1", "D2", "D3"]):
    y = Filt(D_all[k])
    z_no = drho*y
    gz_no = float(np.dot(g_tot16, z_no))
    noeta[nm] = gz_no
    log("CTX %s: <g_tot16, drho*y (eta-response OFF)>=%.10f" % (nm, gz_no))
res["context_no_eta"] = noeta

with open(OUT + "/b17_phase0_results.json", "w") as f:
    json.dump(res, f, indent=1)
log("AUDIT_DONE t=%.1fs" % (time.time()-t0))
