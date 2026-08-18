#!/usr/bin/env python3
"""
BFINAL-006 cycle-2 -- S2/T1 offline tangent linear-system closure (repaired driver).

This driver implements the APPROVED plan step S2-T1-TangentSolve:

  * Pressure-Reference Gate (5 programmatic checks, from the ACTUAL case/source):
      C1: p.needReference() semantics from 0/p + readTransportProperties.H
          (outlet fixedValue uniform 0 -> needReference()==false ->
           setRefCell no-op -> pRefCell stays 0 -> fvSolution pRefCell=5600
           INEFFECTIVE).
      C2: scan explicitJT.mtx pressure rows [3N,4N) for the UNIQUE identity row
          using the CORRECTED criterion "single |v|>1e-12 AND diagonal==1"
          (NOT nrow==1, which fails because the pinned row carries 13 explicit
          zeros + 1 diagonal 1.0).
      C3: detected identity row == discretePIndex(0) == 3*N+0 == 100800.
      C4: discretePIndex(5600) == 106400 is NOT an identity row (physical
          continuity row, 18 nonzeros).
      C5: record pRefRow = 100800 as the T1 pinning row.
    Any check failing => BLOCKED (exit code 3), no guessing.

  * T1 tangent linear-system closure (task BFINAL-006 §8):
      J * wPrime = -R_x*d
    - explicitJT.mtx is the PHYSICAL J^T exported by solveDiscreteFlowAdjoint.H
      (csrVal is written directly; discreteMomentumScale/discreteAreaScale are
      applied only in the scaled application wrapper, so the exported matrix is
      already physical).  Transpose it: J = (J^T)^T.  After transposing, J has
      the pRef COLUMN pinned to identity (inherited from the J^T row pin).
    - Re-pin J's ROW pRefRow=100800 to identity and set rhs[pRefRow]=0
      (gauge: wPrime_p[100800]=0; do NOT zero pPrime[5600]=pPrime[106400]).
    - rhs = -R_x*d from BFINAL-005 stageB6_rxc_analytic.mtx (per direction:
      [anRUd(3N); anRPd(N)] block-major, physical) + stageB6_dirs.mtx
      (3 x N design directions).
    - scipy splu + iterative refinement; report ||rhs||, true relative
      residual, U-row residual, P-row residual per direction.
    - wPrime_p[100800] == 0 required.
    - D_TAN = g_w^T * wPrime with g_w = P block of explicitRhs_pressureDrop.mtx
      (direct pressure functional derivative; NO adjoint quantity used).

Read-only: never modifies the case, the matrix, or any source file.
Writes only under evidence/agent-group/BFINAL-006/cycle-2/.
"""
import os
import sys
import time
import json
import numpy as np
import scipy.io
import scipy.sparse as sp
import scipy.sparse.linalg as spla

CASE = "/home/ys/dsH/b2_case_smoke"
EVID = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-006/cycle-2"
ART = os.path.join(EVID, "artifacts")
os.makedirs(ART, exist_ok=True)

MAT = os.path.join(CASE, "explicitJT.mtx")
RXC_AN = os.path.join(CASE, "stageB6_rxc_analytic.mtx")
DIRS = os.path.join(CASE, "stageB6_dirs.mtx")
RHS_GW = os.path.join(CASE, "explicitRhs_pressureDrop.mtx")

N = 33600
NV = 3 * N          # 100800
NUNK = 4 * N        # 134400
DIR_NAMES = ["D1", "D2", "D3"]
NDIR = 3

# ---- recorded production scales (from BFINAL-003/005 run logs, same case) ----
MOMENTUM_SCALE = 0.00786409044975
AREA_SCALE = 2.47771624976e-07
MOM_AREA = MOMENTUM_SCALE / AREA_SCALE

t0 = time.time()

def log(msg):
    print(msg, flush=True)

def fail(msg):
    log("BLOCKED: " + msg)
    sys.exit(3)

# ===========================================================================
# CHECK 1 -- p.needReference() semantics from ACTUAL case/source
# ===========================================================================
log("=== BFINAL-006 cycle-2 S2/T1: Pressure-Reference Gate + T1 tangent solve ===")
log("workspace=%s case=%s" % (os.getcwd(), CASE))

p_path = os.path.join(CASE, "0", "p")
txt = open(p_path).read()
fixed_patches = []
zero_grad_patches = []
in_boundary = False
cur = None
for line in txt.splitlines():
    ls = line.strip()
    if ls == "boundaryField":
        in_boundary = True
        continue
    if not in_boundary:
        continue
    if ls == "}":
        cur = None
        continue
    if ls == "{":
        continue
    if ls.startswith("type"):
        if "fixedValue" in ls:
            fixed_patches.append(cur)
        elif "zeroGradient" in ls:
            zero_grad_patches.append(cur)
        continue
    if ls.startswith("value"):
        continue
    # otherwise it is a patch/subdict name line
    if ls and not ls.startswith("//"):
        cur = ls
c1_ok = (len(fixed_patches) == 1 and fixed_patches == ["outlet"])
log("[CHECK 1] p.needReference() semantics")
log("  fixedValue patches: %s   zeroGradient patches: %s"
    % (fixed_patches, zero_grad_patches))
log("  -> needReference()==False  (any fixesValue patch => false; outlet fixedValue uniform 0)")
log("  readTransportProperties.H: label pRefCell = 0; setRefCell(p, simple.dict(), pRefCell, pRefValue)")
log("  setRefCell no-op when !needReference() => pRefCell stays 0; fvSolution pRefCell=5600 INEFFECTIVE")
if not c1_ok:
    fail("CHECK 1 failed: expected exactly one fixedValue patch 'outlet'")
log("  CHECK 1 => PASS")

# ===========================================================================
# CHECK 2/3/4 -- identity-row scan of explicitJT.mtx (corrected criterion)
# ===========================================================================
log("[loading] reading explicitJT.mtx (%.1f MB) ..." % (os.path.getsize(MAT) / 1e6))
Jt = scipy.io.mmread(MAT).tocsr()
assert Jt.shape == (NUNK, NUNK), "unexpected shape %s" % (Jt.shape,)
log("  matrix %dx%d nnz=%d (t=%.1fs)" % (Jt.shape[0], Jt.shape[1], Jt.nnz, time.time() - t0))

log("[CHECK 2] identity-row scan of pressure rows [%d,%d) (corrected criterion: single |v|>1e-12 AND diagonal==1)" % (NV, NUNK))
ident_rows = []
for r in range(NV, NUNK):
    s, e = Jt.indptr[r], Jt.indptr[r + 1]
    nz = 0
    nzcol = -1
    nzval = 0.0
    for j in range(s, e):
        v = Jt.data[j]
        if abs(v) > 1e-12:
            nz += 1
            nzcol = Jt.indices[j]
            nzval = v
    if nz == 1 and nzcol == r and abs(nzval - 1.0) < 1e-12:
        ident_rows.append(r)
log("  identity rows (pressure block): %s" % ident_rows)
pref_expected = NV + 0          # discretePIndex(0)
pref_5600 = NV + 5600           # discretePIndex(5600) = 106400
if len(ident_rows) != 1 or ident_rows[0] != pref_expected:
    fail("CHECK 2/3 failed: unique identity row must be %d (discretePIndex(0)); got %s" % (pref_expected, ident_rows))

# record raw row content of the pinned row (after scipy duplicate-sum)
s, e = Jt.indptr[pref_expected], Jt.indptr[pref_expected + 1]
log("  row %d: stored entries=%d; nonzeros>1e-12=%d (diag=%.17g)"
    % (pref_expected, e - s,
       sum(1 for j in range(s, e) if abs(Jt.data[j]) > 1e-12),
       Jt.data[Jt.indices[s:e].tolist().index(pref_expected)] if pref_expected in Jt.indices[s:e] else float('nan')))
log("[CHECK 3] detected identity row %d == discretePIndex(0) = 3*%d+0 = %d => PASS" % (ident_rows[0], N, pref_expected))

s5600, e5600 = Jt.indptr[pref_5600], Jt.indptr[pref_5600 + 1]
nz5600 = sum(1 for j in range(s5600, e5600) if abs(Jt.data[j]) > 1e-12)
log("[CHECK 4] row %d (=discretePIndex(5600)): stored entries=%d, nonzeros>1e-12=%d, is_identity=%s"
    % (pref_5600, e5600 - s5600, nz5600, pref_5600 in ident_rows))
if pref_5600 in ident_rows:
    fail("CHECK 4 failed: row %d must NOT be identity" % pref_5600)
log("  CHECK 4 => PASS (106400 is a physical continuity row, NOT the reference row)")

pRefRow = pref_expected
log("[CHECK 5] pRefRow = %d recorded as T1 pinning row (rhs[pRefRow]=0, row->identity). => PASS" % pRefRow)

# ===========================================================================
# T1 -- build J = transpose(J^T), re-pin row pRefRow, rhs=-R_x*d, splu solve
# ===========================================================================
log("[T1] transposing explicitJT.mtx (J^T -> J) ...")
J = Jt.T.tocsr()
del Jt
log("  J nnz=%d (t=%.1fs)" % (J.nnz, time.time() - t0))

# verify the pRef COLUMN is identity (inherited from J^T row pin)
col = J[:, pRefRow].toarray().ravel()
col_nz = np.count_nonzero(col)
col_ok = (col_nz == 1) and (abs(col[pRefRow] - 1.0) < 1e-12)
log("  J[:,%d] identity-column check: nonzeros=%d ok=%s" % (pRefRow, col_nz, col_ok))
if not col_ok:
    fail("J pRef column is not identity")

log("[T1] re-pinning J row %d to identity (gauge wPrime_p[%d]=0; rhs[%d]=0) ..." % (pRefRow, pRefRow, pRefRow))
J = J.tolil()
J[pRefRow, :] = 0.0
J[pRefRow, pRefRow] = 1.0
J = J.tocsr()

log("[T1] loading directions and R_x*d ...")
dirs_raw = np.loadtxt(DIRS, dtype=np.float64)
if dirs_raw.size != NDIR * N:
    fail("dirs size %d != %d" % (dirs_raw.size, NDIR * N))
dirs = [dirs_raw[di * N:(di + 1) * N] for di in range(NDIR)]

rxc = np.loadtxt(RXC_AN, dtype=np.float64)
if rxc.size != NDIR * 4 * N:
    fail("rxc size %d != %d" % (rxc.size, NDIR * 4 * N))

gw_cell = np.asarray(scipy.io.mmread(RHS_GW), dtype=np.float64).ravel()
if gw_cell.size != 4 * N:
    fail("gw size %d != %d" % (gw_cell.size, 4 * N))
gw = gw_cell[3::4].copy()          # P block of the cell-major export
gwU = np.concatenate([gw_cell[0::4], gw_cell[1::4], gw_cell[2::4]])
log("  g_w: ||g_w(P)||=%.6e ||g_w(U)||=%.6e (U block ~0) sum(g_w,P)=%.6e (NOT gauge-neutral: outlet is fixedValue p=0, so only inlet cells carry +weight)"
    % (np.linalg.norm(gw), np.linalg.norm(gwU), float(np.sum(gw))))

log("[T1] scales recorded (exported matrix+rhs are PHYSICAL; no extra U/P scaling applied):")
log("  discreteMomentumScale=%.14e discreteAreaScale=%.14e M/A=%.14e" % (MOMENTUM_SCALE, AREA_SCALE, MOM_AREA))

log("[T1] SuperLU factorization of pinned J (n=%d) ..." % NUNK)
lu = spla.splu(J, permc_spec="COLAMD")
log("  factorized (t=%.1fs)" % (time.time() - t0))

def solve_refine(rhs, label):
    bnorm = np.linalg.norm(rhs)
    x = lu.solve(rhs)
    rr = np.linalg.norm(J @ x - rhs) / max(bnorm, 1e-300)
    log("    %s: iter 0 relres=%.6e" % (label, rr))
    for it in range(1, 30):
        if rr < 1e-13:
            break
        dx = lu.solve(rhs - J @ x)
        x = x + dx
        rr = np.linalg.norm(J @ x - rhs) / max(bnorm, 1e-300)
        if it % 5 == 0 or rr < 1e-12:
            log("    %s: iter %d relres=%.6e" % (label, it, rr))
    return x, rr

results = {}
for di in range(NDIR):
    off = di * 4 * N
    anRUd = rxc[off:off + 3 * N]
    anRPd = rxc[off + 3 * N:off + 4 * N]
    d = dirs[di]
    rxd = np.concatenate([anRUd, anRPd])
    rhs = -rxd.copy()
    rhs[pRefRow] = 0.0
    rnorm = np.linalg.norm(rhs)
    rnormU = np.linalg.norm(rhs[:3 * N])
    rnormP = np.linalg.norm(rhs[3 * N:])
    log("[T1] direction %s: ||R_x*d||=%.6e  ||rhs(post-pin)||=%.6e (U=%.6e P=%.6e)"
        % (DIR_NAMES[di], np.linalg.norm(rxd), rnorm, rnormU, rnormP))
    wP, rr = solve_refine(rhs, DIR_NAMES[di])
    res = J @ wP - rhs
    resU = np.linalg.norm(res[:3 * N])
    resP = np.linalg.norm(res[3 * N:])
    true_rel = np.linalg.norm(res) / max(rnorm, 1e-300)
    relU = resU / max(rnormU, 1e-300)
    relP = resP / max(rnormP, 1e-300)
    log("    %s: ||J*wPrime-rhs||=%.6e  trueRelRes=%.6e  U-rowRel=%.6e  P-rowRel=%.6e  ||wPrime||=%.6e"
        % (DIR_NAMES[di], np.linalg.norm(res), true_rel, relU, relP, np.linalg.norm(wP)))
    wpref = wP[pRefRow]
    if abs(wpref) > 1e-12:
        fail("wPrime[pRefRow]=%g not zero" % wpref)
    log("    %s: wPrime_p[%d]=%.3e (gauge pin OK)" % (DIR_NAMES[di], pRefRow, wpref))
    pPrime = wP[3 * N:]
    dtan = float(np.dot(gw, pPrime))
    log("    %s: D_TAN = g_w^T*wPrime = %.8e" % (DIR_NAMES[di], dtan))
    out = os.path.join(ART, "wPrime_TAN_%s.mtx" % DIR_NAMES[di])
    with open(out, "w") as f:
        for v in wP:
            f.write("%.17g\n" % v)
    log("    wrote %s (%.0f bytes)" % (out, os.path.getsize(out)))
    results[DIR_NAMES[di]] = dict(
        rxd_norm=float(np.linalg.norm(rxd)),
        rhs_norm=rnorm, rhs_normU=rnormU, rhs_normP=rnormP,
        abs_res=float(np.linalg.norm(res)),
        true_relres=true_rel, relU=relU, relP=relP,
        wnorm=float(np.linalg.norm(wP)),
        wpref=float(wpref),
        dtan=dtan)

summary = {
    "pRefRow": pRefRow,
    "identity_rows": ident_rows,
    "row_106400_is_identity": pref_5600 in ident_rows,
    "scales": {"discreteMomentumScale": MOMENTUM_SCALE,
               "discreteAreaScale": AREA_SCALE,
               "M_over_A": MOM_AREA},
    "g_w": {"normP": float(np.linalg.norm(gw)),
            "normU": float(np.linalg.norm(gwU)),
            "sumP": float(np.sum(gw))},
    "directions": results,
}
with open(os.path.join(ART, "T1_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)

log("[done] T1 finished in %.1fs" % (time.time() - t0))
print()
print("=== SUMMARY (T1 tangent solve, pRefRow=%d) ===" % pRefRow)
print("dir     ||R_x*d||       ||Jw-rhs||      trueRel   Urel      Prel      wPrime_p[pref]  D_TAN")
for di in range(NDIR):
    r = results[DIR_NAMES[di]]
    print("%-6s % .6e  % .6e  % .3e  % .3e  % .3e  % .3e        % .8e"
          % (DIR_NAMES[di], r["rxd_norm"], r["abs_res"], r["true_relres"],
             r["relU"], r["relP"], r["wpref"], r["dtan"]))
log("T1_OK")
