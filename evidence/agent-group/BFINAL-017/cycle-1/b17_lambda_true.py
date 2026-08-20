#!/usr/bin/env python3
"""BFINAL-017 phase-0 FINAL: exact-adjoint adjudication of the assembly-vs-
lambda-direct contradiction.

Established (b17_phase0 + followups):
  * chain is an exact transpose (AD1, 1e-15) -> assembly == <g_total(lambda), z>
  * the B2-baseline assembly pass used lambda_679 (gradProxy -219717.0666,
    B13 in-pass momentum -3.5036), while the WRITTEN Uc/pc are lambda_902
    (gradProxy -205194.1411, bitwise equal in b12/b13/b16 runs)
  * B16's lambda-direct gauge used the WRITTEN lambda_902
    -> the "147% assembly contradiction" compares two different adjoints.

Here: solve the exported unpinned system exactly (machine precision):
    lambda* = M^{-1} bPD        (M = explicitJT = J^T export)
    w*     = M^{-T} bTC         (primal truth for the operator gauge)
and re-evaluate every gauge:
    lambda-direct(lambda*), lambda-direct(lambda_902), assembly(lambda_679)
    operator ratios vs truth b^T w*
    |M^{-1} rxd| (ill-conditioning amplification of the lambda-direct gauge)
"""
import time, json
import numpy as np
import scipy.io
import scipy.sparse as sp
import scipy.sparse.linalg as spla

EXPORT_JT = "/home/ys/dsH/b8_verify_diag/explicitJT.mtx"
E = "/home/ys/dsH/b15_export"
S = "/home/ys/dsH/b16_states/1"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-017/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N; P0 = NV

t0 = time.time()
def log(m): print(m, flush=True)

def load_rhs_export(path):
    v = np.asarray(scipy.io.mmread(path)).ravel()
    m = v.reshape(N, 4)
    b = np.zeros(NUNK)
    b[0:NV:3] = m[:, 0]; b[1:NV:3] = m[:, 1]; b[2:NV:3] = m[:, 2]
    b[P0:P0+N] = m[:, 3]
    return b

def read_of_scalar(path):
    vals = []; intxt = False
    for line in open(path):
        s = line.strip()
        if not intxt:
            if s == "(": intxt = True
            continue
        if s == ")": break
        for tok in s.replace(";", " ").split():
            for tt in tok.replace("(", " ").replace(")", " ").split():
                vals.append(float(tt))
    return np.array(vals)

def read_of_vector(path):
    txt = open(path).read()
    i = txt.index("\n("); j = txt.index("\n)", i)
    return np.fromstring(txt[i+2:j].replace("(", " ").replace(")", " "),
                         sep=" ").reshape(-1, 3)

log("loading explicitJT (467MB) ...")
M = scipy.io.mmread(EXPORT_JT).tocsr()
log("matrix loaded nnz=%d (t=%.1fs)" % (M.nnz, time.time()-t0))
bPD = load_rhs_export(E + "/explicitRhs_pressureDrop.mtx")
bTC = load_rhs_export(E + "/explicitRhs_thermalCoupling.mtx")

lu = spla.splu(M.tocsc())
log("splu done (t=%.1fs)" % (time.time()-t0))

lam_star = lu.solve(bPD)
r = M @ lam_star - bPD
log("lambda*: |M@lam*-bPD|=%.3e  rel=%.3e"
    % (np.linalg.norm(r), np.linalg.norm(r)/np.linalg.norm(bPD)))
w_star = lu.solve(bTC, trans="T")   # M^T w = bTC  <=>  J w = bTC
r2 = M.T @ w_star - bTC
log("w*: |J@w*-bTC|=%.3e rel=%.3e"
    % (np.linalg.norm(r2), np.linalg.norm(r2)/np.linalg.norm(bTC)))

# written lambda_902
pc_w = read_of_scalar(S + "/pc"); Uc_w = read_of_vector(S + "/Uc")
lam902 = np.zeros(NUNK)
lam902[0:NV:3] = Uc_w[:, 0]; lam902[1:NV:3] = Uc_w[:, 1]; lam902[2:NV:3] = Uc_w[:, 2]
lam902[P0:] = pc_w
cosv = float(np.dot(lam902, lam_star)/(np.linalg.norm(lam902)*np.linalg.norm(lam_star)))
log("lambda_902(written) vs lambda*: cos=%.12f relL2=%.3e"
    % (cosv, np.linalg.norm(lam902-lam_star)/np.linalg.norm(lam_star)))
# gradProxy(lambda*) = sum (U & lamU) * D1
U = read_of_vector(S + "/U")
D = np.loadtxt(E + "/stageB6_dirs.mtx").reshape(3, N)
gp = float(np.sum(np.einsum("ij,ij->i", U, lam_star[:NV].reshape(N, 3))*D[0]))
log("gradProxy(lambda*)=%.10f   [round1=-205194.1411 round2=-219717.0666]" % gp)
gp902 = float(np.sum(np.einsum("ij,ij->i", U, lam902[:NV].reshape(N, 3))*D[0]))
log("gradProxy(lambda_902)=%.10f" % gp902)

# lambda-direct contractions
rxc = np.loadtxt(E + "/stageB6_rxc_analytic.mtx")
res = {}
for k, nm in enumerate(["D1", "D2", "D3"]):
    rxd = rxc[k*NUNK:(k+1)*NUNK]
    ld_star = -float(np.dot(lam_star, rxd))
    ld_902 = -float(np.dot(lam902, rxd))
    # ill-conditioning amplification: x = M^{-1} rxd ; functional err = r^T x
    x = lu.solve(rxd)
    amp = np.linalg.norm(x)/np.linalg.norm(rxd)
    res[nm] = {"lambda_direct_true": ld_star, "lambda_direct_902": ld_902,
               "assembly_B12": [-5.4050587339, 1.10592698054, 14.0447142088][k],
               "inv_norm_amplification": float(amp)}
    log("%s: -lam*^T rxd=%.10f | -lam902^T rxd=%.10f | assembly(679)=%.10f "
        "| |M^-1 rxd|/|rxd|=%.3e"
        % (nm, ld_star, ld_902, res[nm]["assembly_B12"], amp))

# operator-defect gauge with the exact lambda: ratio vs the B15 state-FD truth
# (gDP_wt = bPD^T w_true, per-direction FD of the perturbed primal solves;
#  h=0.001 column of BFINAL-015 t2t3_results.json)
log("--- operator gauge re-rating ---")
tt = json.load(open("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/"
                    "BFINAL-015/cycle-1/t2t3_results.json"))["contractions"]
for k, nm in enumerate(["D1", "D2", "D3"]):
    rxd = rxc[k*NUNK:(k+1)*NUNK]
    ld_star = -float(np.dot(lam_star, rxd))
    ld_902 = -float(np.dot(lam902, rxd))
    truth = tt[nm + "_h0.001"]["gDP_wt"]
    wp15 = tt[nm + "_h0.001"]["gDP_wp"]
    res[nm]["truth_gdp_wt"] = truth
    res[nm]["truth_gdp_wp_b15"] = wp15
    res[nm]["ratio_true"] = ld_star/truth
    res[nm]["ratio_902"] = ld_902/truth
    res[nm]["ratio_assembly_679"] = res[nm]["assembly_B12"]/truth
    log("%s: truth=%.6f | lam*-direct=%.6f (ratio %.4f) | "
        "lam902-direct=%.6f (ratio %.4f) | assembly679=%.6f (ratio %.4f) "
        "| b15 wprime-tangent=%.6f"
        % (nm, truth, ld_star, ld_star/truth, ld_902, ld_902/truth,
           res[nm]["assembly_B12"], res[nm]["assembly_B12"]/truth, wp15))

with open(OUT + "/b17_lambda_true.json", "w") as f:
    json.dump(res, f, indent=1)
log("LAMBDA_TRUE_DONE t=%.1fs" % (time.time()-t0))
