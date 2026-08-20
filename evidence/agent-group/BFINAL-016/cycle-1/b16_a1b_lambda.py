#!/usr/bin/env python3
"""
BFINAL-016 A1b - localize the third production-vs-export operator difference.

  lambda_prod : from the state run's written adjoint fields (1/Uc, 1/pc)
                (production FGMRES solve at the same main-loop state as the
                 b15 export run - bit-identical lineage).
  lambda_exp  : lu.solve(bPD, trans='T') on the exported matrix M (= J^T),
                i.e. the exported-system adjoint.

Diagnostics:
  1. bPD^T lambda_prod (independent recomputation of the production ADJ
     contraction) and bPD^T lambda_exp.
  2. Block metrics lambda_prod vs lambda_exp (U/P split, cos, relL2).
  3. THE BLOCK DIFFERENTIAL: r = M @ lambda_prod - bPD.  If the production
     adjoint fails the EXPORTED system specifically on one block, that block
     is where the two operator assemblies differ.  Report block norms and
     the spatial structure of the largest residuals.
"""
import time, json
import numpy as np
import scipy.io, scipy.sparse as sp, scipy.sparse.linalg as spla

EXPORT_JT = "/home/ys/dsH/b8_verify_diag/explicitJT.mtx"
EXPORT_CASE = "/home/ys/dsH/b15_export"
FIELDS = "/home/ys/dsH/b16_states/1"
CYCLE = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-016/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N; PREF = NV; P0 = NV

t0 = time.time()
def log(m): print(m, flush=True)
log("=== BFINAL-016 A1b: production-vs-export lambda localization ===")

def read_of_scalar(path):
    vals = []
    with open(path) as f:
        intxt = False
        for line in f:
            s = line.strip()
            if not intxt:
                if s == "(":
                    intxt = True
                continue
            if s == ")":
                break
            # uniform list may be 'value;' lines
            for tok in s.replace(";", " ").split():
                vals.append(float(tok))
    return np.array(vals)

def read_of_vector(path):
    txt = open(path).read()
    i = txt.index("\n(")
    j = txt.index("\n)", i)
    body = txt[i+2:j]
    nums = np.fromstring(
        body.replace("(", " ").replace(")", " "), sep=" ")
    return nums.reshape(-1, 3)

pc = read_of_scalar(FIELDS + "/pc")
Uc = read_of_vector(FIELDS + "/Uc")
log("fields read: |Uc|=%.6e |pc|=%.6e (t=%.1fs)" %
    (np.linalg.norm(Uc), np.linalg.norm(pc), time.time() - t0))
lam_prod = np.zeros(NUNK)
lam_prod[0:NV:3] = Uc[:, 0]
lam_prod[1:NV:3] = Uc[:, 1]
lam_prod[2:NV:3] = Uc[:, 2]
lam_prod[P0:] = pc

def load_rhs_export(path):
    v = np.asarray(scipy.io.mmread(path)).ravel()
    m = v.reshape(N, 4)
    b = np.zeros(NUNK)
    b[0:NV:3] = m[:, 0]; b[1:NV:3] = m[:, 1]; b[2:NV:3] = m[:, 2]
    b[P0:P0 + N] = m[:, 3]
    return b

bPD = load_rhs_export(EXPORT_CASE + "/explicitRhs_pressureDrop.mtx")
bTC = load_rhs_export(EXPORT_CASE + "/explicitRhs_thermalCoupling.mtx")

# exported matrix M = J^T (as exported); build J (pinned) for the existing
# splu convention and use transposed solves for the adjoint system M x = b.
Mraw = scipy.io.mmread(EXPORT_JT).tocsr()
J = Mraw.T.tolil()
J[PREF, :] = 0.0; J[PREF, PREF] = 1.0
J = J.tocsr()
lu = spla.splu(J, permc_spec="COLAMD")
log("splu done (t=%.1fs)" % (time.time() - t0))
del Mraw

# lambda_exp: J^T x = b  (M x = b)
lam_exp = lu.solve(bPD, trans="T")
for it in range(30):
    rres = bPD - J.T @ lam_exp
    if np.linalg.norm(rres)/np.linalg.norm(bPD) < 1e-13:
        break
    lam_exp = lam_exp + lu.solve(rres, trans="T")

log("[1] contractions: bPD.lam_prod=%.10f  bPD.lam_exp=%.10f  "
    "bTC.lam_prod=%.10f" % (bPD @ lam_prod, bPD @ lam_exp, bTC @ lam_prod))

# [2] block metrics
d = lam_prod - lam_exp
cosv = lam_prod @ lam_exp/(np.linalg.norm(lam_prod)*np.linalg.norm(lam_exp))
log("[2] lam_prod vs lam_exp: cos=%.6f relL2=%.4e relU=%.4e relP=%.4e" %
    (cosv, np.linalg.norm(d)/np.linalg.norm(lam_exp),
     np.linalg.norm(d[:NV])/np.linalg.norm(lam_exp[:NV]),
     np.linalg.norm(d[P0:])/np.linalg.norm(lam_exp[P0:])))

# [3] block differential: exported operator applied to production lambda
r = J.T @ lam_prod - bPD
rU = r[:NV]; rP = r[P0:]
log("[3] M@lam_prod - bPD: |r|=%.6e (rel %.4e)  |rU|=%.6e  |rP|=%.6e" %
    (np.linalg.norm(r), np.linalg.norm(r)/np.linalg.norm(bPD),
     np.linalg.norm(rU), np.linalg.norm(rP)))
ru_c = np.sqrt(rU[0::3]**2 + rU[1::3]**2 + rU[2::3]**2)
topU = np.argsort(ru_c)[::-1][:10]
topP = np.argsort(np.abs(rP))[::-1][:10]
log("    top U-row residual cells: %s" % list(map(int, topU[:8])))
log("    |rU|/|bU|=%.4e  |rP|/|bP|=%.4e" %
    (np.linalg.norm(rU)/max(np.linalg.norm(bPD[:NV]), 1e-300),
     np.linalg.norm(rP)/max(np.linalg.norm(bPD[P0:]), 1e-300)))

out = {"bPD_lam_prod": float(bPD @ lam_prod),
       "bPD_lam_exp": float(bPD @ lam_exp),
       "bTC_lam_prod": float(bTC @ lam_prod),
       "lam_cos": float(cosv),
       "lam_relL2": float(np.linalg.norm(d)/np.linalg.norm(lam_exp)),
       "lam_relU": float(np.linalg.norm(d[:NV])/np.linalg.norm(lam_exp[:NV])),
       "lam_relP": float(np.linalg.norm(d[P0:])/np.linalg.norm(lam_exp[P0:])),
       "resid_rU": float(np.linalg.norm(rU)),
       "resid_rP": float(np.linalg.norm(rP)),
       "resid_rel": float(np.linalg.norm(r)/np.linalg.norm(bPD))}
with open(CYCLE + "/a1b_lambda_localization.json", "w") as f:
    json.dump(out, f, indent=1)
log("A1B_DONE t=%.1fs" % (time.time() - t0))
