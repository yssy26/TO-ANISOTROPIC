#!/usr/bin/env python3
"""
BFINAL-016 A1c - clean follow-up (no factorization needed):
  1. r = M_unpinned @ lambda_prod - bPD  (block split, spatial structure)
     -> localizes the operator rows where the production adjoint fails the
        exported system (the pin artifact is removed).
  2. proper ADJ contractions: -lambda^T rxd for lambda_prod (per direction)
     vs the BFINAL-012 ADJ table.
"""
import time, json
import numpy as np
import scipy.io, scipy.sparse as sp

EXPORT_JT = "/home/ys/dsH/b8_verify_diag/explicitJT.mtx"
EXPORT_CASE = "/home/ys/dsH/b15_export"
FIELDS = "/home/ys/dsH/b16_states/1"
CYCLE = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-016/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N; P0 = NV

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-016 A1c: clean block differential ===")

def read_of_scalar(path):
    vals = []
    intxt = False
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

pc = read_of_scalar(FIELDS + "/pc"); Uc = read_of_vector(FIELDS + "/Uc")
lam = np.zeros(NUNK)
lam[0:NV:3] = Uc[:, 0]; lam[1:NV:3] = Uc[:, 1]; lam[2:NV:3] = Uc[:, 2]
lam[P0:] = pc

def load_rhs_export(path):
    v = np.asarray(scipy.io.mmread(path)).ravel()
    m = v.reshape(N, 4)
    b = np.zeros(NUNK)
    b[0:NV:3] = m[:, 0]; b[1:NV:3] = m[:, 1]; b[2:NV:3] = m[:, 2]
    b[P0:P0 + N] = m[:, 3]
    return b

bPD = load_rhs_export(EXPORT_CASE + "/explicitRhs_pressureDrop.mtx")
M = scipy.io.mmread(EXPORT_JT).tocsr()   # M = J^T as exported, unpinned
log("loaded (t=%.1fs)" % (time.time() - t0))

r = M @ lam - bPD
rU = r[:NV]; rP = r[P0:]
log("[1] M_unpinned@lam_prod - bPD: |r|=%.6e |rU|=%.6e |rP|=%.6e" %
    (np.linalg.norm(r), np.linalg.norm(rU), np.linalg.norm(rP)))
log("    relative: |rU|/|bU|=%.3e  |rP|/|bP|=%.3e" %
    (np.linalg.norm(rU)/max(np.linalg.norm(bPD[:NV]), 1e-300),
     np.linalg.norm(rP)/max(np.linalg.norm(bPD[P0:]), 1e-300)))
ru_c = np.sqrt(rU[0::3]**2 + rU[1::3]**2 + rU[2::3]**2)
log("    top U-res cells: %s" % list(map(int, np.argsort(ru_c)[::-1][:8])))
log("    top P-res cells: %s" % list(map(int, np.argsort(np.abs(rP))[::-1][:8])))

rxc = np.loadtxt(EXPORT_CASE + "/stageB6_rxc_analytic.mtx", dtype=np.float64)
log("[2] ADJ contractions -lam^T rxd:")
out = {"rU": float(np.linalg.norm(rU)), "rP": float(np.linalg.norm(rP))}
for di, name in enumerate(["D1", "D2", "D3"]):
    rxd = rxc[di*NUNK:(di+1)*NUNK]
    c = -float(lam @ rxd)
    log("    %s: -lam_prod^T rxd=%.10f" % (name, c))
    out[name] = c
with open(CYCLE + "/a1c_clean.json", "w") as f:
    json.dump(out, f, indent=1)
log("A1C_DONE t=%.1fs" % (time.time() - t0))
