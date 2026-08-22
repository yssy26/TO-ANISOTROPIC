#!/usr/bin/env python3
"""BFINAL-029 Phase 1 M1-direction state probe (decisive).

The b29 slot-table M1 preview (state-mix: B2 basis + /1 g-layer) gives
D1 ratio=+0.4665 / D3 ratio=+0.2049, while the export bexp gives
D1 ratio=-0.6561 / D3 ratio=+1.0573 (B26 b26a_m1_results.json).  The
state_test found all-B2 route cos=0.9877 with the export (vs 0.9594 for
state-mix) but did NOT compute M1 contractions.  This script computes
lhs = b^T wH7[d] and ratio vs RHS[d] for the export and every route state
(state-mix / all-/1 / all-B2) plus each single-field swap, to determine
whether ANY route state reproduces the export's D1 sign and D3 magnitude,
or whether the residual is structural (slot folding) rather than state.

w_true (wH7) and rhs = FD_J - C - Gx come from B26 artifacts; no FD
recompute (Phase-1 offline discipline).
"""
import time, json, re
import numpy as np
import importlib.util

MESH = "/home/ys/dsH/b16_mesh"
SB = "/home/ys/dsH/b25_qgate/stageB2"
Q1 = "/home/ys/dsH/b25_qgate/1"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-029/cycle-1"
B25 = "/home/ys/dsH/b25_qgate"
E26 = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-026/cycle-1"

t0 = time.time()
log = lambda m: print(m, flush=True)

# load b29_state_test machinery (re-runs its idempotent computations)
spec = importlib.util.spec_from_file_location("b29_state_test", OUT + "/b29_state_test.py")
st = importlib.util.module_from_spec(spec)
spec.loader.exec_module(st)

N = st.N; NV = st.NV; NUNK = st.NUNK

# ---------------- exported b_TC (same packing as st) ----------------
m_exp = np.loadtxt(B25 + "/b18rhs_thermalCoupling.mtx")
bexp = np.zeros(NUNK)
bexp[0:NV:3] = m_exp[:, 0]; bexp[1:NV:3] = m_exp[:, 1]
bexp[2:NV:3] = m_exp[:, 2]; bexp[NV:] = m_exp[:, 3]

# ---------------- w_true + rhs from B26 (no recompute) ----------------
wH7 = np.load(E26 + "/b26_wstar_h7.npz")
res26 = json.load(open(E26 + "/b26a_m1_results.json"))
RHS = {}
for r in res26["rows"]:
    RHS[r["dir"]] = r["rhs"]

def m1_row(name, b):
    out = {}
    for k in ("D1", "D2", "D3"):
        lhs = float(b @ wH7[k])
        rhs = RHS[k]
        out[k] = dict(lhs=lhs, rhs=rhs, ratio=lhs / rhs,
                      reldev=abs(lhs / rhs - 1.0),
                      sign_lhs=int(np.sign(lhs)), sign_rhs=int(np.sign(rhs)))
        log("[M1 %-10s] %s: lhs=%+.8e ratio=%+.5f  (bexp %+.5f)"
            % (name, k, lhs, lhs / rhs, bexp @ wH7[k] / rhs))
    return out

log("=== M1 state probe ===")
rows = {"export": m1_row("export", bexp)}
rows["state-mix"] = m1_row("state-mix", st.r_mix)
rows["all-1"] = m1_row("all-1", st.r_all1)
rows["all-B2"] = m1_row("all-B2", st.r_allB)

# single-field swaps from state-mix anchor (each swap -> route)
def swp(**kw):
    d = dict(st.fMix); d.update(kw); return st.route(d)

for name, kw in [("swap T->B2", dict(T=st.T_B)),
                 ("swap phiTh->B2", dict(phiTh=st.phTh_B)),
                 ("swap phibval->B2", dict(phi_bvals=st.phiB_B)),
                 ("swap Tb->(/1 proxy)", dict(Tb=st.Tb1))]:
    b = swp(**kw)
    rows[name] = m1_row(name, b)

# ---------------- summary table ----------------
print("\n=== M1 ratio matrix (lhs/rhs) ===")
print("%-20s %10s %10s %10s" % ("route", "D1", "D2", "D3"))
for nm, r in rows.items():
    print("%-20s %10.4f %10.4f %10.4f" % (nm, r["D1"]["ratio"], r["D2"]["ratio"], r["D3"]["ratio"]))
print("export(bexp) is the reference: D1=-0.6561 D2=-10.0836 D3=+1.0573")
print("\namplitude: |bexp U|L2=%.3e P=%.3e" % (np.linalg.norm(bexp[:NV]), np.linalg.norm(bexp[NV:])))
for nm in ("state-mix", "all-1", "all-B2"):
    b = {"state-mix": st.r_mix, "all-1": st.r_all1, "all-B2": st.r_allB}[nm]
    print("  %-10s |U|L2=%.3e P=%.3e  cos=%.4f"
          % (nm, np.linalg.norm(b[:NV]), np.linalg.norm(b[NV:]),
             float(b @ bexp) / (np.linalg.norm(b) * np.linalg.norm(bexp))))

json.dump({nm: {k: v for k, v in r.items()} for nm, r in rows.items()},
          open(OUT + "/b29_m1_state_probe.json", "w"), indent=1)
log("\nwrote %s/b29_m1_state_probe.json (t=%.1fs)" % (OUT, time.time() - t0))
