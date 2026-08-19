#!/usr/bin/env python3
"""
Fresh verification (post-6320943) — Gate P2 (nullspace removal, UNPINNED
physical matrix) and Gate P8 (direct tangent solvability) plus the
identity-row audit, on the freshly exported corrected explicitJT.mtx.

Key difference vs the archived s4_p2_p8_nullspace_tangent.py: the old script
RE-PINNED row 100800 (P(0)) to an identity row following the BFINAL-006
legacy convention.  The current HEAD removes the pin for fixedValue-pressure
cases (discretePressureNeedsReference == p.needReference() == false), so this
script must:
  * verify the exported P block contains ZERO identity rows;
  * NOT re-pin anything;
  * compute sigma_min / null-vector ratios / tangent solves on the physical
    matrix as exported.

Read-only: no source/case file is modified.
"""
import os, time, json
import numpy as np
import scipy.io, scipy.sparse as sp, scipy.sparse.linalg as spla

CYCLE = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-008/cycle-1"
EXPORT = "/home/ys/dsH/b8_verify_diag/explicitJT.mtx"
RXC    = "/home/ys/dsH/b8_verify_diag/stageB6_rxc_analytic.mtx"
NPVEC  = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-006/cycle-2/artifacts/wPrime_TAN_D1.mtx"
OUTCELLS = CYCLE + "/artifacts/stageB8_outlet_cells.mtx"

N = 33600
NV = 3 * N
NUNK = 4 * N
P0 = NV            # 100800 = discretePIndex(0)

t0 = time.time()
def log(m):
    print(m, flush=True)

log("=== FRESH UNPINNED P2+P8 verification (HEAD 6320943) ===")

# ---------------- load fresh export, build J (NO re-pin) ----------------
M = scipy.io.mmread(EXPORT).tocsr()
log("fresh export loaded nnz=%d (t=%.1fs)" % (M.nnz, time.time() - t0))
J = M.T.tocsr()
del M
log("J built (physical, unpinned) nnz=%d (t=%.1fs)" % (J.nnz, time.time() - t0))

# ---------------- identity-row audit of the P block ----------------
rows = J[P0:P0 + N, :].tocsr()
idrows = []
for r in range(N):
    st, en = rows.indptr[r], rows.indptr[r + 1]
    if en - st == 1 and abs(rows.data[st]) > 1e-12 and abs(rows.data[st] - 1.0) < 1e-12:
        idrows.append(r)
log("identity rows in P block [100800,134400): %s  (MUST be [] for fixedValue outlet)"
    % idrows)
# row 100800 physical check
r1008 = J.getrow(100800).tocoo()
log("row 100800: nnz=%d (physical, NOT identity) %s"
    % (r1008.nnz, "" if r1008.nnz > 1 else "!!ARTIFICIAL IDENTITY ROW!!"))
r1064 = J.getrow(106400).tocoo()
log("row 106400: nnz=%d (physical, not identity) %s"
    % (r1064.nnz, "" if r1064.nnz > 1 else "!!"))

# ---------------- null vector n (full state) ----------------
n = np.loadtxt(NPVEC, dtype=np.float64)
nn = np.linalg.norm(n)
n = n / nn
Jn = J @ n
ratio = np.linalg.norm(Jn) / np.linalg.norm(n)
nP = n[P0:P0 + N]
log("null vector: ||wPrime_TAN_D1||=%.6e ||n_P||=%.6f n[PREF]=%.3e" %
    (nn, np.linalg.norm(nP), n[P0]))
log("||J_physical n||/||n|| = %.6e   (OLD exported J: 2.717e-19; BFINAL-007 candidate: 3.020e-09)"
    % ratio)

# ---------------- splu + sigma_min on the UNPINNED physical J ----------------
log("[LU] splu(physical unpinned J, COLAMD) ...")
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
log("sigma_min(physical unpinned J): min=%.3e mean=%.3e" % (min(sig), float(np.mean(sig))))

# ---------------- smallest right-singular direction geometry ----------------
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
log("physical J v_min(P-block): top-1%% mass=%.3f outlet-cell mass=%.3f"
    % (topmass, outmass))
log("  max |v_P| at cell %d (outlet-cell? %s)" % (int(order[0]), str(int(order[0]) in outcells)))

# ---------------- P8 tangent solves on the UNPINNED physical J ----------------
rxc = np.loadtxt(RXC, dtype=np.float64)
log("[solves] J_physical wPrime = -R_x d, D1/D2/D3, iterative refinement (30 iters)")
res8 = []
for di, name in enumerate(["D1", "D2", "D3"]):
    off = di * NUNK
    rxd = rxc[off:off + NUNK]
    rhs = -rxd.copy()
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
        % (name, rnorm, rr, relU, relP, np.linalg.norm(x), x[P0]))
    res8.append({"dir": name, "||rhs||": rnorm, "trueRelRes": rr,
                 "Urel": relU, "Prel": relP, "||wPrime||": float(np.linalg.norm(x)),
                 "wPREF": float(x[P0])})

out = {
    "verification": "fresh-unpinned post-6320943",
    "P2": {
        "identity_rows_P_block": idrows,
        "row_100800_nnz": int(r1008.nnz),
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
with open(CYCLE + "/s4_verify_fresh_unpinned.json", "w") as f:
    json.dump(out, f, indent=1)
log("FRESH_VERIFY_DONE t=%.1fs -> s4_verify_fresh_unpinned.json" % (time.time() - t0))
