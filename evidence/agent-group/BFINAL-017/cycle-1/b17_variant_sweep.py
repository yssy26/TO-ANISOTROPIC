#!/usr/bin/env python3
"""BFINAL-017 phase-0: variant sweep to identify what the b16/B2-module
context changed in the production rxPressureRowT rebuild.

Baseline rebuild (b17_rebuildT.py) reproduces the oracle side (V2 relL2~8e-3)
and the b15 raw export (V1 relL2~6e-3), but NOT the b16-implied production
field (rel 0.97).  Sweep plausible variants of the UEqn rebuild and see which
one reproduces the b16-implied T16.
"""
import time
import numpy as np

_src = open("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-017/cycle-1/b17_rebuildT.py").read()
_pre = _src.split("rAU, HbyA = build_ueqn(bounded_corr=True)")[0]
exec(_pre)

t0 = time.time()
def log(m): print(m, flush=True)

mask = read_field(S + "/designMask")[0]
xp = read_field(S + "/xp")[0]

# rebuild T16-implied + T15 export + oracle targets (same as b17_rebuildT)
lam15 = np.loadtxt(E + "/stageB6_lambda.mtx")
pc15 = lam15[NV:]
g15_pr = np.loadtxt(E + "/rxpr_prod_gsensh_pressurerow.mtx")
supp = daD != 0
T15_export = np.zeros(N); T15_export[supp] = -g15_pr[supp]/daD[supp]
maskb = mask > 0.5
DEL = 8.0
def heaviside(g, eta):
    e = np.exp(-DEL); lo = g <= eta; out = np.empty_like(g)
    gl = g[lo]/eta; out[lo] = eta*(np.exp(-DEL*(1.0-gl))-(1.0-gl)*e)
    gh = (g[~lo]-eta)/(1.0-eta); out[~lo] = eta+(1.0-eta)*(1.0-np.exp(-DEL*gh)+gh*e)
    return out
def diffvol(g, eta): return np.sum((g[maskb]-heaviside(g[maskb], eta))*V[maskb])
e0, e1 = 1e-4, 0.9999; y0, y1 = diffvol(xp, e0), diffvol(xp, e1)
while (e1-e0) > 1e-10:
    e5 = 0.5*(e0+e1); y5 = diffvol(xp, e5)
    if y0*y5 < 0: e1, y1 = e5, y5
    else: e0, y0 = e5, y5
eta5 = 0.5*(e0+e1)
idx = np.where(maskb)[0]; xpm = xp[idx]; ex = np.exp(-DEL); lo = xpm <= eta5
drho = np.zeros(N); dPde = np.zeros(N)
pe = np.exp(-DEL*(1.0-xpm[lo]/eta5)); drho[idx[lo]] = DEL*pe+ex; dPde[idx[lo]] = pe*(1.0-DEL*xpm[lo]/eta5)-ex
pe2 = np.exp(-DEL*(xpm[~lo]-eta5)/(1.0-eta5)); drho[idx[~lo]] = DEL*pe2+ex; dPde[idx[~lo]] = pe2*(1.0-DEL*(1.0-xpm[~lo])/(1.0-eta5))-ex
Dden = np.sum(mask*dPde*V)
w16 = read_field(S + "/gsenshPressureDrop")[0]
rw = w16[idx]/drho[idx]
Cinv = np.sum((mask*dPde)[idx]*rw)/(Dden+np.sum((mask*dPde*V*(1.0-drho))[idx]/drho[idx]))
gtot = np.zeros(N); gtot[idx] = rw-(V*(1.0-drho)*Cinv)[idx]/drho[idx]
gmom = -daD*np.einsum("ij,ij->i", U, Uc16)*V
gpr = gtot-gmom
T16_implied = np.zeros(N); T16_implied[supp] = -gpr[supp]/daD[supp]

z15 = np.loadtxt(E + "/stageB6_rxc_z_analytic.mtx").reshape(3, N)
rxc = np.loadtxt(E + "/stageB6_rxc_analytic.mtx")
oracle_pr = {}
for k, nm in enumerate(["D1", "D2", "D3"]):
    oracle_pr[nm] = -float(np.dot(pc16, rxc[k*NUNK:(k+1)*NUNK][NV:]))
implied_pr = {}
for k, nm in enumerate(["D1", "D2", "D3"]):
    implied_pr[nm] = float(np.dot(gpr, z15[k]))

def build_T_var(pc, alpharel=0.4, use_dev2=True, bounded=True,
                dhbya_mode="prod"):
    global ALPHAREL
    ALPHAREL = alpharel
    rAU, HbyA = build_ueqn(bounded_corr=bounded)
    if not use_dev2:
        # rebuild HbyA without the dev2 source term: approximate by
        # reconstructing from the same pipeline minus dev2 is complex; skip
        raise NotImplementedError
    drAU = -rAU*rAU/alpharel
    if dhbya_mode == "prod":
        dHbyA = (rAU/alpharel)[:, None]*((1.0-alpharel)*U - HbyA)
    elif dhbya_mode == "noalpha":   # forget the (1-alphaRel)U relaxation term
        dHbyA = (rAU/alpharel)[:, None]*(-HbyA)
    g0_int = magSf*deltaCoeffs*(p[nei] - p[owner])
    T1 = np.zeros(N); T2 = np.zeros(N)
    pcdiff = pc[owner] - pc[nei]
    SfdH_own = np.einsum("fi,fi->f", Sf, dHbyA[owner])
    SfdH_nei = np.einsum("fi,fi->f", Sf, dHbyA[nei])
    np.add.at(T1, owner, weights*SfdH_own*pcdiff)
    np.add.at(T1, nei, (1.0-weights)*SfdH_nei*pcdiff)
    np.add.at(T2, owner, -weights*g0_int*pcdiff*drAU[owner])
    np.add.at(T2, nei, -(1.0-weights)*g0_int*pcdiff*drAU[nei])
    am = U_assign
    np.add.at(T1, bc_cell[am],
              np.einsum("fi,fi->f", bSf[am], dHbyA[bc_cell[am]])*pc[bc_cell[am]])
    return T1+T2, T1, T2

log("=== variant sweep: which rebuild matches the b16-implied T16? ===")
variants = [
    ("baseline   aRel=0.4", dict(alpharel=0.4)),
    ("aRel=1.0           ", dict(alpharel=1.0)),
    ("aRel=0.7           ", dict(alpharel=0.7)),
    ("aRel=0.5           ", dict(alpharel=0.5)),
    ("aRel=0.3           ", dict(alpharel=0.3)),
    ("no (1-a)U term     ", dict(alpharel=0.4, dhbya_mode="noalpha")),
    ("no bounded corr    ", dict(alpharel=0.4, bounded=False)),
]
for nm, kw in variants:
    try:
        T16v, T1v, T2v = build_T_var(pc16, **kw)
        T15v, _, _ = build_T_var(pc15, **kw)
        d16 = np.linalg.norm((T16v-T16_implied)[supp])/np.linalg.norm(T16_implied[supp])
        d15 = np.linalg.norm((T15v-T15_export)[supp])/np.linalg.norm(T15_export[supp])
        # contraction against z for this variant
        cs = [float(np.dot(-T16v*daD, z15[k])) for k in range(3)]
        log("%s: vs T16_implied rel=%.3e | vs T15_export rel=%.3e | pr(d1,d2,d3)=%.4f %.4f %.4f (implied %.4f %.4f %.4f, oracle %.4f %.4f %.4f)"
            % (nm, d16, d15, cs[0], cs[1], cs[2],
               implied_pr["D1"], implied_pr["D2"], implied_pr["D3"],
               oracle_pr["D1"], oracle_pr["D2"], oracle_pr["D3"]))
    except Exception as ex:
        log("%s: FAILED %s" % (nm, ex))

# per-cell ratio structure of the baseline mismatch
T16b, T1b, T2b = build_T_var(pc16)
rat = np.zeros(N); ok = supp & (np.abs(T16b) > 1e-14)
rat[ok] = T16_implied[ok]/T16b[ok]
log("ratio implied/baseline on support: median=%.4f p10=%.4f p90=%.4f"
    % (np.median(rat[ok]), np.percentile(rat[ok], 10), np.percentile(rat[ok], 90)))
# how much of implied is explained by T1-only / T2-only / -T1 / -T2 / combos
for nm, fld in [("T1", T1b), ("T2", T2b), ("-T1", -T1b), ("-T2", -T2b)]:
    d = np.linalg.norm((fld-T16_implied)[supp])/np.linalg.norm(T16_implied[supp])
    log("  |%s - implied|/|implied| = %.3e" % (nm, d))
log("SWEEP_DONE t=%.1fs" % (time.time()-t0))
