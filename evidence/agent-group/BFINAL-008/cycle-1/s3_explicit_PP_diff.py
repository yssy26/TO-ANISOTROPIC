#!/usr/bin/env python3
# BFINAL-008 S3 — verify the corrected explicitJT.mtx (run-1) vs the OLD
# explicitJT.mtx (real case, pre-patch export): P-P block must differ by
# (1) interior sign flip (NEW = -OLD on interior P-P entries),
# (2) +rAtU_c*deltaCoeffs_b*|Sf|_b on the 84 outlet-cell diagonals,
# (3) unchanged pRef identity row (row 100800 = P(0)).
# Layout: 134400 unknowns = 33600 cells x 4 (3 U + 1 P); P index = 100800 + c.
import numpy as np
import scipy.io
import time

OLD = "/home/ys/dsH/b2_case_smoke/explicitJT.mtx"
NEW = "/home/ys/dsH/b8_scratch_r1/explicitJT.mtx"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-008/cycle-1/artifacts/s3/explicitJT.mtx"

N = 134400
NP = 33600
P0 = 100800

t0 = time.time()
print("loading OLD ...", flush=True)
A_old = scipy.io.mmread(OLD).tocsr()
print(f"  OLD loaded {time.time()-t0:.1f}s shape={A_old.shape} nnz={A_old.nnz}", flush=True)
t0 = time.time()
print("loading NEW ...", flush=True)
A_new = scipy.io.mmread(NEW).tocsr()
print(f"  NEW loaded {time.time()-t0:.1f}s nnz={A_new.nnz}", flush=True)

# P-P block comparison: rows in [P0, P0+NP), cols in [P0, P0+NP)
# (NEW is the archived copy; use it — same file as b8_scratch_r1)
B_old = A_old[P0:P0+NP, P0:P0+NP]
B_new = A_new[P0:P0+NP, P0:P0+NP]
print(f"P-P block OLD nnz={B_old.nnz} NEW nnz={B_new.nnz}")

# structural set
coo_o = B_old.tocoo(); coo_n = B_new.tocoo()
set_o = set(zip(coo_o.row.tolist(), coo_o.col.tolist()))
set_n = set(zip(coo_n.row.tolist(), coo_n.col.tolist()))
only_o = set_o - set_n
only_n = set_n - set_o
print(f"only-in-OLD P-P entries: {len(only_o)}  only-in-NEW: {len(only_n)}")

# interior entries present in both: check sign flip NEW == -OLD
common = set_o & set_n
rows_o = coo_o.row; cols_o = coo_o.col; vals_o = coo_o.data
rows_n = coo_n.row; cols_n = coo_n.col; vals_n = coo_n.data
d_old = {(int(r), int(c)): float(v) for r, c, v in zip(rows_o, cols_o, vals_o)}
d_new = {(int(r), int(c)): float(v) for r, c, v in zip(rows_n, cols_n, vals_n)}

import random
random.seed(7)
sample = random.sample(sorted(common), 4000)
rel = []
for (r, c) in sample:
    vo, vn = d_old[(r, c)], d_new[(r, c)]
    if abs(vo) > 1e-30:
        rel.append((vn + vo)/abs(vo))  # NEW + OLD should be ~0 if flip (interior)
rel = np.array(rel)
print(f"interior common entries (sample 4000): mean((NEW+OLD)/|OLD|)={rel.mean():.3e} "
      f"max|.|={np.abs(rel).max():.3e}")

# outlet cells: NEW diag - (-OLD diag) should equal +rAtU*delta*|Sf| (boundary term)
outlet_cells = [int(x) for x in open(
    "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-008/cycle-1/artifacts/stageB8_outlet_cells.mtx").read().split()]
print(f"outlet cells: {len(outlet_cells)}")
bd = []
for c in outlet_cells:
    r = c  # P block row = cell index
    vo = d_old.get((r, r), 0.0)
    vn = d_new.get((r, r), 0.0)
    bd.append(vn + vo)  # NEW - (-OLD) = NEW + OLD = boundary term
bd = np.array(bd)
print(f"outlet-cell (NEW+OLD) diag: min={bd.min():.6e} max={bd.max():.6e} mean={bd.mean():.6e}  "
      f"(expect [+5.48e-10, +1.71e-9] per BFINAL-007 Q2)")

# pRef row: row 100800 (P cell 0) must be identity (unchanged)
row_old = A_old.getrow(P0).tocoo()
row_new = A_new.getrow(P0).tocoo()
print(f"pRef row OLD: {sorted(zip(row_old.col.tolist(), row_old.data.tolist()))}")
print(f"pRef row NEW: {sorted(zip(row_new.col.tolist(), row_new.data.tolist()))}")
print(f"pRef row identical: {row_old.col.tolist()==row_new.col.tolist() and np.allclose(row_old.data, row_new.data)}")

# off-outlet P-P diagonals (excl pRef): NEW == -OLD (flip only, no boundary term)
offout = [c for c in range(NP) if c not in set(outlet_cells) and c != 0]
sample2 = random.sample(offout, 2000)
rel2 = []
for c in sample2:
    vo = d_old.get((c, c), 0.0); vn = d_new.get((c, c), 0.0)
    if abs(vo) > 1e-30:
        rel2.append((vn + vo)/abs(vo))
rel2 = np.array(rel2)
print(f"off-outlet diag (sample 2000): mean((NEW+OLD)/|OLD|)={rel2.mean():.3e} max|.|={np.abs(rel2).max():.3e}")
