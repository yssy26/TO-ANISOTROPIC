#!/usr/bin/env python3
"""
BFINAL-007 POST-REVIEWER independent sigma_min rerun (candidate J).
Fresh code path: builds the diagnostic-only candidate J (P-P := -L_prod) and
estimates the smallest singular value via 3 random solves of a single splu
factorization (same method as the executor, but independent implementation).
Output -> post_review_sigma.json (does NOT overwrite executor's s3_q4a_sigma.json).
"""
import json
import time
import numpy as np
import scipy.io
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ART = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-007/artifacts"
CASE = "/home/ys/dsH/b2_case_smoke"
N = 33600
NV = 3 * N
NUNK = 4 * N
PREF = NV

t0 = time.time()
print("loading exported J^T ...", flush=True)
Jt = scipy.io.mmread(CASE + "/explicitJT.mtx").tocsr()
J = Jt.T.tocsr().tolil()
J[PREF, :] = 0.0
J[PREF, PREF] = 1.0
J = J.tocsr()
del Jt

diagB = np.loadtxt(ART + "/stageB7_laplacian_diag_bnd.mtx")
upper = np.loadtxt(ART + "/stageB7_laplacian_upper.mtx")
nIntF = upper.shape[0]
own = np.array([int(x) for x in open(CASE + "/constant/polyMesh/owner").read()[open(CASE + "/constant/polyMesh/owner").read().index("(")+1:open(CASE + "/constant/polyMesh/owner").read().rindex(")")].split()[:nIntF]], dtype=np.int64)
nei = np.array([int(x) for x in open(CASE + "/constant/polyMesh/neighbour").read()[open(CASE + "/constant/polyMesh/neighbour").read().index("(")+1:open(CASE + "/constant/polyMesh/neighbour").read().rindex(")")].split()[:nIntF]], dtype=np.int64)

# candidate P-P block = -L_prod (residual derivative of the production pEqn)
candPP = sp.dok_matrix((N, N), dtype=np.float64)
for c in range(N):
    candPP[c, c] = -diagB[c]
for f in range(nIntF):
    candPP[own[f], nei[f]] += -upper[f]
    candPP[nei[f], own[f]] += -upper[f]
candPP = candPP.tocsr()

# full candidate J: exported J with P-P replaced
Jc = J.tolil()
for c in range(N):
    r = NV + c
    keep = [(cc, v) for cc, v in zip(Jc.rows[r], Jc.data[r]) if not (NV <= cc < NUNK)]
    Jc.rows[r] = [cc for cc, _ in keep]
    Jc.data[r] = [v for _, v in keep]
for c in range(N):
    r = NV + c
    for j in range(candPP.indptr[c], candPP.indptr[c + 1]):
        Jc[r, NV + candPP.indices[j]] = candPP.data[j]
Jc = Jc.tocsr().tolil()
Jc[PREF, :] = 0.0; Jc[PREF, PREF] = 1.0
Jc = Jc.tocsr()
print("candidate J built nnz=%d t=%.1fs" % (Jc.nnz, time.time() - t0), flush=True)

lu = spla.splu(Jc.tocsr(), permc_spec="COLAMD")
print("splu done t=%.1fs" % (time.time() - t0), flush=True)
trials = []
for k in range(3):
    rng = np.random.default_rng(2000 + k)
    z = rng.standard_normal(NUNK)
    x = lu.solve(z)
    s = np.linalg.norm(z) / np.linalg.norm(x)
    trials.append(float(s))
    print("trial %d sigma_min~%.3e" % (k, s), flush=True)
res = {"sigma_min_trials": trials, "sigma_min_min": min(trials),
       "sigma_min_mean": float(np.mean(trials)), "_runtime_s": time.time() - t0}
with open(ART + "/post_review_sigma.json", "w") as f:
    json.dump(res, f, indent=1)
print("POST_REVIEW_SIGMA_DONE", json.dumps(res), flush=True)
