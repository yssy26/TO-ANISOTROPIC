#!/usr/bin/env python3
"""
BFINAL-008 S4 (Executor, offline) — independent recomputation of gates P1
(targeted null-vector actual-residual closure), P3 (general actual-residual
B-F2), P4 (J/J^T transpose closure, from run log + explicit-matrix structure),
P5 (explicit / matrix-free consistency), P7 (R_x re-audit).

Fresh code path; reads ONLY raw artifacts (stageB8_* exports from run-1,
corrected explicitJT.mtx, BFINAL-007 laplacian/boundary exports, BFINAL-006
null vector, stageB6 R_x/design artifacts) and the two run logs.  Read-only.
"""
import json, os, re, math
import numpy as np
import scipy.io, scipy.sparse as sp

CYCLE = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-008/cycle-1"
ART = CYCLE + "/artifacts"
S3 = ART + "/s3"
B7 = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-007/artifacts"
CASE = "/home/ys/dsH/b2_case_smoke"
LOG1 = "/home/ys/dsH/b8_scratch_r1/Log.stageB8.r1.txt"

N = 33600
NV = 3 * N
NUNK = 4 * N
P0 = NV          # 100800
PREF = P0

def read_vec(path, n=None):
    with open(path) as f:
        vals = [float(x) for x in f.read().split()]
    return np.array(vals[:n] if n else vals)

def read_celltype(path):
    with open(path) as f:
        return np.array([int(x) for x in f.read().split()])

def read_of(path, n):
    t = open(path).read()
    i = t.index('(')
    return np.array([float(x) for x in t[i + 1:t.rindex(')')].split()[:n]])

def block(a, b, idx):
    av = a[idx]; bv = b[idx]
    d = av - bv
    nA = np.dot(av, av); nB = np.dot(bv, bv); nD = np.dot(d, d)
    if nB == 0.0:
        return {"relL2": None, "cos": None, "normRatio": None,
                "|Jv|": math.sqrt(nA), "|FD|": 0.0}
    return {"relL2": math.sqrt(nD / nB),
            "cos": float(np.dot(av, bv) / math.sqrt(nA * nB)) if nA > 0 else None,
            "normRatio": math.sqrt(nA / nB),
            "|Jv|": math.sqrt(nA), "|FD|": math.sqrt(nB)}

out = {"HEAD": "ca8a772b8339a36e7906ea32a7125061c47aa280",
       "branch": "agent/dsH-stage-b-validation", "P1": {}, "P3": {}, "P4": {}, "P5": {}, "P7": {}}

# ---------------- cell classification ----------------
ct = read_celltype(ART + "/stageB8_celltype.mtx")
interior = [i for i, c in enumerate(ct) if c == 0]
other = [i for i, c in enumerate(ct) if c == 1]
outlet = [i for i, c in enumerate(ct) if c == 2]
offout = [i for i, c in enumerate(ct) if c != 2]
tot = list(range(N))
Pidx = [P0 + i for i in tot]
Pint = [P0 + i for i in interior]
Pout = [P0 + i for i in outlet]
Poth = [P0 + i for i in other]
Uidx = list(range(NV))

# ---------------- P1 + P3 dir1 (null-vector direction) ----------------
Jv1 = read_vec(ART + "/stageB8_Jv1.mtx", NUNK)
dir1 = read_vec(ART + "/stageB8_dir1.mtx", NUNK)
eps_list = ["1e-2", "1e-3", "1e-4", "1e-5"]
plateau = {}
for eps in eps_list:
    fd = read_vec(ART + ("/stageB8_FD1_eps_%s.mtx" % eps), NUNK)
    m = block(Jv1, fd, Pidx)
    plateau[eps] = m["relL2"]
out["P1"]["full_P"] = block(Jv1, read_vec(ART + "/stageB8_FD1_eps_1e-2.mtx", NUNK), Pidx)
out["P1"]["outlet_adjacent"] = block(Jv1, read_vec(ART + "/stageB8_FD1_eps_1e-2.mtx", NUNK), Pout)
out["P1"]["off_outlet"] = block(Jv1, read_vec(ART + "/stageB8_FD1_eps_1e-2.mtx", NUNK), [P0 + i for i in offout])
fd1 = read_vec(ART + "/stageB8_FD1_eps_1e-2.mtx", NUNK)
sgn = (Jv1[Pidx] * fd1[Pidx] > 0.0).mean()
out["P1"]["sameSignFraction"] = float(sgn)
out["P1"]["FD_mass_outlet_fraction"] = float(np.abs(fd1[Pout]).sum() / np.abs(fd1[Pidx]).sum())
out["P1"]["Jv_mass_outlet_fraction"] = float(np.abs(Jv1[Pout]).sum() / np.abs(Jv1[Pidx]).sum())
out["P1"]["full_P_relL2_eps_plateau"] = plateau
out["P1"]["||n_P||"] = float(np.linalg.norm(dir1[Pidx]))

# ---------------- P3 dir1/dir2/dir3 @ eps=1e-2 ----------------
for name, jvf, fdf in [("dir1", "stageB8_Jv1.mtx", "stageB8_FD1_eps_1e-2.mtx"),
                       ("dir2", "stageB8_Jv2.mtx", "stageB8_FD2_eps_1e-2.mtx"),
                       ("dir3", "stageB8_Jv3.mtx", "stageB8_FD3_eps_1e-2.mtx")]:
    Jv = read_vec(ART + "/" + jvf, NUNK)
    FD = read_vec(ART + "/" + fdf, NUNK)
    out["P3"][name] = {
        "U_rows": block(Jv, FD, Uidx),
        "P_internal": block(Jv, FD, Pint),
        "P_outlet_adjacent": block(Jv, FD, Pout),
        "P_other_boundary": block(Jv, FD, Poth),
        "P_total": block(Jv, FD, Pidx),
    }

# ---------------- P4: run-log BlockDot + reduced cold-flow dot ----------------
logtxt = open(LOG1).read()
p4 = {"BlockDot": [], "reduced_coldflow_dot": []}
for label in ["thermalCoupling", "pressureDrop"]:
    m1 = re.search(r"BlockDot \(%s\): PU=(\S+) PP=(\S+) UU=(\S+) UP=(\S+) boundaryU=(\S+)" % label, logtxt)
    m2 = re.search(r"Reduced cold-flow operator transpose dot-test \(%s\) max relative error=(\S+)" % label, logtxt)
    p4["BlockDot"].append({"label": label,
                           "PU": m1.group(1) if m1 else None, "PP": m1.group(2) if m1 else None,
                           "UU": m1.group(3) if m1 else None, "UP": m1.group(4) if m1 else None,
                           "boundaryU": m1.group(5) if m1 else None})
    p4["reduced_coldflow_dot"].append({"label": label,
                                       "maxRelErr": m2.group(1) if m2 else None})
out["P4"] = p4

# ---------------- P5: explicit / matrix-free consistency ----------------
# (a) corrected explicit P-P block == -L_prod (BFINAL-007 laplacian exports,
#     same design point/HEAD; primal untouched by the patch)
d = np.loadtxt(B7 + "/stageB7_laplacian_diag.mtx")          # L_int diag
db = np.loadtxt(B7 + "/stageB7_laplacian_diag_bnd.mtx")     # full diag (int+bnd)
u = np.loadtxt(B7 + "/stageB7_laplacian_upper.mtx")         # interior upper kf
bt = db - d                                                  # bnd internalCoeffs (84 cells)
Ldiag = -(d + bt)                                            # -L_prod diag (residual conv)
Lup = -u                                                     # -L_prod upper

M = scipy.io.mmread(S3 + "/explicitJT.mtx").tocsr()
B = M[P0:P0 + N, P0:P0 + N]                                  # P-P block of J^T (symmetric)
# explicit P-P diag vs -L_prod diag (exclude pRef cell 0, pinned to 1)
excl0 = [c for c in range(N) if c != 0]
diag_expl = np.array([B[c, c] for c in excl0])
diag_ref = Ldiag[excl0]
p5_pp_diag = float(np.linalg.norm(diag_expl - diag_ref) / np.linalg.norm(diag_ref))
# explicit P-P off-diag vs -L_prod upper (owner->nei pairs)
off_ref = {}
for fi in range(len(u)):
    # upper file: face order matches interior faces; we need owner/nei — not
    # exported here, so compare the SYMMETRIC difference instead: the P-P block
    # must be symmetric (laplacian + diagonal) -> (B - B^T) ~ 0.  The pRef
    # row/col (cell 0, pinned to identity in the export) is excluded.
    pass
sub = (B[1:, 1:] - B[1:, 1:].T).tocsr()
p5_pp_sym = float(np.abs(sub.data).max()) if sub.nnz else 0.0
# outlet-cell explicit diag vs -d - bt (corrected boundary internalCoeffs entry)
bnd_ref = -bt[outlet]
bnd_expl = np.array([B[c, c] for c in outlet]) - (-d[outlet])
p5_bnd = float(np.linalg.norm(bnd_expl - bnd_ref) / np.linalg.norm(bnd_ref))
# pRef row: row P0 of export must be identity (unchanged pin)
rrow = M.getrow(P0).tocoo()
p5_pref = {"nnz": int(rrow.nnz),
           "cols": rrow.col.tolist(),
           "vals": [float(v) for v in rrow.data]}
# (b) explicit J (= export^T) applied to dir_i vs matrix-free Jv_i (excl. pRef
#     row).  The export pins row 100800 of J^T, i.e. COLUMN 100800 of J = e, so
#     explicit J drops the physical p_cell0 coupling: directions with
#     v[PREF]=0 give a clean comparison; directions with v[PREF]!=0 differ
#     exactly by J_phys[i,PREF]*v[PREF] on the pRef-coupled rows (documented
#     BFINAL-003 G4 convention).
p5_jv = {}
for name, jvf, dirf in [("dir1", "stageB8_Jv1.mtx", "stageB8_dir1.mtx"),
                        ("dir2", "stageB8_Jv2.mtx", "stageB8_dir2.mtx"),
                        ("dir3", "stageB8_Jv3.mtx", "stageB8_dir3.mtx")]:
    v = read_vec(ART + "/" + dirf, NUNK)
    Jv_mf = read_vec(ART + "/" + jvf, NUNK)
    Jv_ex = (M.T @ v)
    keep = [i for i in range(NUNK) if i != PREF]
    d = Jv_ex - Jv_mf
    big = np.where(np.abs(d) > 1e-10)[0]
    rel = np.linalg.norm(d[keep]) / np.linalg.norm(Jv_mf[keep])
    relP = np.linalg.norm(d[P0 + 1:]) / np.linalg.norm(Jv_mf[P0 + 1:])
    relU = np.linalg.norm(d[:NV]) / np.linalg.norm(Jv_mf[:NV])
    p5_jv[name] = {"vPREF": float(v[PREF]),
                   "relL2_excl_pRef": float(rel), "relL2_P_excl_cell0": float(relP),
                   "relL2_U": float(relU),
                   "max_abs_diff": float(np.abs(d).max()),
                   "n_diff_rows_gt_1e-10": int(len(big)),
                   "diff_rows": big.tolist()[:12]}

out["P5"] = {"PP_diag_vs_minusLprod": p5_pp_diag,
             "PP_sym_max_abs_excl_cell0": p5_pp_sym,
             "outlet_bnd_vs_minusbt": p5_bnd,
             "pRef_row": p5_pref,
             "explicitJ_vs_matrixfree_Jv": p5_jv}

# ---------------- P7: R_x re-audit ----------------
op = np.loadtxt(B7 + "/stageB7_outlet_patch.mtx")   # cell delta magSf ic bc mobB mobC cand
ocells = op[:, 0].astype(int)
delta = op[:, 1]; magSf = op[:, 2]
rAU = np.loadtxt(B7 + "/stageB7_rebuilt_rAU.mtx")
p = read_of(CASE + "/1/p", N)
alphaRel = 0.4
drAU = -rAU ** 2 / alphaRel
dphi_b_da = drAU[ocells] * delta * magSf * p[ocells]
da = np.loadtxt(S3 + "/stageB6_deltaAlpha.mtx")
dm = read_of(CASE + "/1/designMask", N)
out["P7"] = {
    "Q5_dphi_b_da_min": float(dphi_b_da.min()),
    "Q5_dphi_b_da_max": float(dphi_b_da.max()),
    "Q5_dphi_b_da_L2": float(np.linalg.norm(dphi_b_da)),
    "Q5_deltaAlpha_outlet_nnz": int((np.abs(da[ocells]) > 1e-300).sum()),
    "Q5_designMask_outlet_max": float(dm[ocells].max()),
    "Q5_ref_Jp_w_L2": 6.22344977639e-10,
    "Q5_ratio_bnd_vs_Jp_w": float(np.linalg.norm(dphi_b_da) / 6.22344977639e-10),
}
# BFINAL-005 D1/D2/D3 contraction recompute from raw artifacts.
# Convention (BFINAL-005 post-reviewer): D = -lambda^T R_x d (Lagrangian minus,
# J^T lambda = g_w with the residual sign).  Oracle values must reproduce the
# locked BFINAL-004/005 anchors exactly.
lam = read_vec(S3 + "/stageB6_lambda.mtx", NUNK)
rxc = read_vec(S3 + "/stageB6_rxc_analytic.mtx", 3 * NUNK)
prod = read_vec(S3 + "/stageB6_prod_gsens.mtx", N)
prodVol = read_vec(S3 + "/stageB6_prod_gsensVol.mtx", N)
dirs = read_vec(S3 + "/stageB6_dirs.mtx", 3 * N)
locked_tot = [-0.0562356333787, 0.297672431753, -0.0363414298395]
locked_mom = [-0.0486867554277, 0.229079662362, -0.0301314016208]
locked_pre = [-0.00754887795106, 0.0685927693909, -0.00621002821874]
locked_vol = [0.1555186511144686, -0.08449667741658891, 0.06396237684432686]
cont = []
for di in range(3):
    rxd = rxc[di * NUNK:(di + 1) * NUNK]
    dvec = dirs[di * N:(di + 1) * N]
    lU = float(np.dot(lam[:NV], rxd[:NV]))
    lP = float(np.dot(lam[NV:], rxd[NV:]))
    D_mom = -lU
    D_pre = -lP
    D_tot = -(lU + lP)
    g3m = float(np.dot(prod, dvec))       # production momentum field contraction
    gv = float(np.dot(prodVol, dvec))     # production volume field contraction
    cont.append({"dir": "D%d" % (di + 1),
                 "D_total_oracle": D_tot, "locked_D_total": locked_tot[di],
                 "relErr_Dtot_vs_locked": abs(D_tot - locked_tot[di]) / abs(locked_tot[di]),
                 "D_momentum_oracle": D_mom, "locked_D_momentum": locked_mom[di],
                 "relErr_Dmom_vs_locked": abs(D_mom - locked_mom[di]) / abs(locked_mom[di]),
                 "D_pressure_oracle": D_pre, "locked_D_pressure": locked_pre[di],
                 "relErr_Dpre_vs_locked": abs(D_pre - locked_pre[di]) / abs(locked_pre[di]),
                 "prod_momentum_dot_d": g3m,
                 "relErr_G2_momentum": abs(g3m - D_mom) / abs(D_mom),
                 "prod_volume_dot_d": gv,
                 "relErr_G5_volume": abs(gv - locked_vol[di]) / abs(locked_vol[di])})
out["P7"]["BFINAL005_contraction"] = cont

with open(CYCLE + "/s4_p1_p3_p4_p5_p7.json", "w") as f:
    json.dump(out, f, indent=1)
print(json.dumps(out, indent=1))
print("S4_FAST_DONE")
