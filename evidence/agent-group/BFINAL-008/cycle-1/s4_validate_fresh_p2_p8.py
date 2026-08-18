#!/usr/bin/env python3
"""
BFINAL-008 S4 (Executor, offline) — Gate P2 (nullspace removal / invertibility)
and Gate P8 (direct tangent solvability smoke) on the CORRECTED production J.

Fresh code path (no reuse of BFINAL-006/007 scripts beyond file conventions):

  J := transpose(corrected explicitJT.mtx)  [corrected export = J^T]
  pin row 100800 (P(0)) -> identity          [BFINAL-006 S1 gate convention]
  splu (COLAMD) once; then:
    * sigma_min via ||z||/||J^-1 z|| random solves (3 trials, fresh seeds)
    * ||J n||/||n|| for the BFINAL-006/007 null vector n = wPrime_TAN_D1 (norm.)
    * inverse-iteration estimate of the smallest right-singular direction and
      its P-block outlet concentration (old mode: top-1% mass fraction 0.895,
      outlet-plane localization) -- must NOT be outlet-concentrated now
    * unique identity-row scan of the P block (expect exactly [100800])
    * P8: solve J wPrime = -R_x d for D1 (strong) + D2/D3 (supplementary)
      with iterative refinement (30 iters), true relative residual <= 1e-9

Read-only: no source/case file is modified.
"""
import os, time, json
import numpy as np
import scipy.io, scipy.sparse as sp, scipy.sparse.linalg as spla

CYCLE = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-008/cycle-1"
EXPORT = "/home/ys/dsH/b8_scratch_v1/explicitJT.mtx"
RXC    = "/home/ys/dsH/b8_scratch_v1/stageB6_rxc_analytic.mtx"
NPVEC  = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-006/cycle-2/artifacts/wPrime_TAN_D1.mtx"
OUTCELLS = CYCLE + "/artifacts/stageB8_outlet_cells.mtx"

N = 33600
NV = 3 * N
NUNK = 4 * N
PREF = NV            # 100800 = discretePIndex(0)
P0 = NV

t0 = time.time()
def log(m):
    print(m, flush=True)

log("=== BFINAL-008 S4 P2+P8 (corrected production J) ===")
log("HEAD ca8a772b branch agent/dsH-stage-b-validation (recorded in JSON)")

# ---------------- load corrected export, build J, pin ----------------
M = scipy.io.mmread(EXPORT).tocsr()
log("export loaded nnz=%d (t=%.1fs)" % (M.nnz, time.time() - t0))
J = M.T.tocsr()
del M
J = J.tolil()
J[PREF, :] = 0.0
J[PREF, PREF] = 1.0
J = J.tocsr()
log("re-pinned corrected J built nnz=%d (t=%.1fs)" % (J.nnz, time.time() - t0))

# ---------------- identity-row scan of the P block ----------------
rows = J[P0:P0 + N, :].tocsr()
idrows = []
for r in range(N):
    st, en = rows.indptr[r], rows.indptr[r + 1]
    if en - st == 1 and abs(rows.data[st]) > 1e-12 and abs(rows.data[st] - 1.0) < 1e-12:
        idrows.append(r)
log("identity rows in P block [100800,134400): %s" % idrows)
# row 106400 = discretePIndex(5600) physical check
r1064 = J.getrow(106400).tocoo()
log("row 106400: nnz=%d (physical, not identity) %s" % (r1064.nnz, "" if r1064.nnz > 1 else "!!"))

# ---------------- null vector n (full state) ----------------
n = np.loadtxt(NPVEC, dtype=np.float64)
nn = np.linalg.norm(n)
n = n / nn
Jn = J @ n
ratio = np.linalg.norm(Jn) / np.linalg.norm(n)
nP = n[P0:P0 + N]
log("null vector: ||wPrime_TAN_D1||=%.6e ||n_P||=%.6f n[PREF]=%.3e" %
    (nn, np.linalg.norm(nP), n[PREF]))
log("||J_corrected n||/||n|| = %.6e   (OLD exported J: 2.717e-19; candidate: 3.020e-09)"
    % ratio)

# ---------------- splu + sigma_min ----------------
log("[LU] splu(corrected J, COLAMD) ...")
lu = spla.splu(J, permc_spec="COLAMD")
log("splu done (t=%.1fs)" % (time.time() - t0))

sig = []
for k in range(3):
    rng = np.random.default_rng(20260818 + k)
    z = rng.standard_normal(NUNK)
    x = lu.solve(z)
    s = np.linalg.norm(z) / np.linalg.norm(x)
    sig.append(float(s))
    log("  trial %d: sigma_min~%.3e (||x||=%.3e)" % (k, s, np.linalg.norm(x)))
log("sigma_min(corrected J): min=%.3e mean=%.3e" % (min(sig), float(np.mean(sig))))

# ---------------- null-vector geometry of the corrected J ----------------
# estimate the smallest right-singular direction by inverse iteration
rng = np.random.default_rng(777)
v = rng.standard_normal(NUNK)
v = v / np.linalg.norm(v)
for it in range(3):
    v = lu.solve(v)
    v = v / np.linalg.norm(v)
vP = np.abs(v[P0:P0 + N])
order = np.argsort(vP)[::-1]
top1pct = order[: int(0.01 * N)]
outcells = set(int(x) for x in open(OUTCELLS).read().split())
outmass = vP[list(outcells)].sum() / vP.sum()
topmass = vP[top1pct].sum() / vP.sum()
nPabs = np.abs(nP)
ordern = np.argsort(nPabs)[::-1]
oldtop = nPabs[ordern[: int(0.01 * N)]].sum() / nPabs.sum()
oldout = nPabs[list(outcells)].sum() / nPabs.sum()
log("old null mode n_P: top-1%% mass=%.3f outlet-cell mass=%.3f (BFINAL-006: top-1%% 0.895)"
    % (oldtop, oldout))
log("corrected J v_min(P-block): top-1%% mass=%.3f outlet-cell mass=%.3f"
    % (topmass, outmass))
log("  max |v_P| at cell %d (outlet-cell? %s)" % (int(order[0]), str(int(order[0]) in outcells)))
log("  corrected J ||J n_old||/||n_old||=%.3e  (old exported J: 2.717e-19)"
    % ratio)

# ---------------- P8 tangent solves ----------------
rxc = np.loadtxt(RXC, dtype=np.float64)
log("[solves] J wPrime = -R_x d, D1/D2/D3, iterative refinement (30 iters)")
res8 = []
for di, name in enumerate(["D1", "D2", "D3"]):
    off = di * NUNK
    rxd = rxc[off:off + NUNK]
    rhs = -rxd.copy()
    rhs[PREF] = 0.0
    rnorm = np.linalg.norm(rhs)
    x = lu.solve(rhs)
    rr = np.linalg.norm(J @ x - rhs) / rnorm
    for it in range(1, 30):
        if rr < 1e-13:
            break
        x = x + lu.solve(rhs - J @ x)
        rr = np.linalg.norm(J @ x - rhs) / rnorm
    res = J @ x - rhs
    relU = np.linalg.norm(res[:NV]) / max(np.linalg.norm(rhs[:NV]), 1e-300)
    relP = np.linalg.norm(res[P0:]) / max(np.linalg.norm(rhs[P0:]), 1e-300)
    log("  %s: ||rhs||=%.6e trueRelRes=%.6e Urel=%.6e Prel=%.6e ||wPrime||=%.6e w[PREF]=%.3e"
        % (name, rnorm, rr, relU, relP, np.linalg.norm(x), x[PREF]))
    res8.append({"dir": name, "||rhs||": rnorm, "trueRelRes": rr,
                 "Urel": relU, "Prel": relP, "||wPrime||": float(np.linalg.norm(x)),
                 "wPREF": float(x[PREF])})

out = {
    "HEAD": "ca8a772b8339a36e7906ea32a7125061c47aa280",
    "branch": "agent/dsH-stage-b-validation",
    "P2": {
        "identity_rows_P_block": idrows,
        "row_106400_nnz": int(r1064.nnz),
        "||J n||/||n||": ratio,
        "sigma_min_trials": sig,
        "sigma_min_min": min(sig),
        "sigma_min_mean": float(np.mean(sig)),
        "old_mode_top1pct_mass": float(oldtop),
        "old_mode_outlet_mass": float(oldout),
        "corrected_vmin_top1pct_mass": float(topmass),
        "corrected_vmin_outlet_mass": float(outmass),
        "corrected_vmin_max_cell": int(order[0]),
        "corrected_vmin_max_is_outlet": bool(int(order[0]) in outcells),
    },
    "P8": {"solves": res8},
    "t_elapsed_s": time.time() - t0,
}
with open(CYCLE + "/s4_validate_fresh_p2_p8.json", "w") as f:
    json.dump(out, f, indent=1)
log("S4_P2_P8_DONE t=%.1fs -> s4_validate_fresh_p2_p8.json" % (time.time() - t0))
