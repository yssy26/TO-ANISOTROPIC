#!/usr/bin/env python3
"""BFINAL-029 Phase 2: preregistered variant x identity matrix.

Preregistration (b29_PREREGISTRATION.md section 4, binding criteria 4.2):
  V0   = base full route (Phase 1 production mirror [g1_TIn, T8 ON])
  Z_S  = slot S zeroed  (S in {H7s2, T4} + Phase-1 divergent {T3, T6, T2, T1})
  F_S  = slot S sign-flipped
  X    = double cross: H7s2 & T4 simultaneously Z or simultaneously F

Each variant v computes b_TC_variant^T w_true(H7) for D1/D2/D3 and compares
against rhs = FD_J(Q) - C - Gx (B26 artifacts, frozen, no FD recompute).

Binding criteria (4.2):
  CONVICTED <=> D1 ratio in [0.85, 1.15] AND D3 ratio in [0.90, 1.20].
  D2 recorded only.  If no variant passes -> "candidate excluded", STOP.

The slot instrument module is re-executed via importlib (idempotent: mesh,
fields, basis, self-gate assert, slot table).  Its `route()` function is
reused for every variant.
"""
import time, json, importlib.util
import numpy as np

OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-029/cycle-1"
E26 = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-026/cycle-1"
B25 = "/home/ys/dsH/b25_qgate"

t0 = time.time()
log = lambda m: print(m, flush=True)
log("=== BFINAL-029 Phase 2 variant matrix ===")

# re-execute Phase 1 instrument (idempotent, defines route / FULL / ALL_SLOTS)
spec = importlib.util.spec_from_file_location(
    "b29_slot_instrument", OUT + "/b29_slot_instrument.py")
inst = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inst)

N = inst.N; NV = inst.NV; NUNK = inst.NUNK
route = inst.route
FULL = inst.FULL
ALL_SLOTS = inst.ALL_SLOTS

# exported b_TC (reference) + B26 frozen w_true / rhs
m_exp = np.loadtxt(B25 + "/b18rhs_thermalCoupling.mtx")
bexp = np.zeros(NUNK)
bexp[0:NV:3] = m_exp[:, 0]; bexp[1:NV:3] = m_exp[:, 1]
bexp[2:NV:3] = m_exp[:, 2]; bexp[NV:] = m_exp[:, 3]

wH7 = np.load(E26 + "/b26_wstar_h7.npz")
res26 = json.load(open(E26 + "/b26a_m1_results.json"))
RHS = {r["dir"]: r["rhs"] for r in res26["rows"]}
B26RATIO = {r["dir"]: r["ratio"] for r in res26["rows"]}

SUSPECTS = ["H7s2", "T4", "T3", "T6", "T2", "T1"]  # B26 + Phase-1 divergent

def relL2(b): return float(np.linalg.norm(b - bexp) / np.linalg.norm(bexp))
def cosv(b):
    return float(b @ bexp / (np.linalg.norm(b) * np.linalg.norm(bexp)))

def m1_row(b):
    out = {}
    for k in ("D1", "D2", "D3"):
        lhs = float(b @ wH7[k])
        rhs = RHS[k]
        out[k] = dict(lhs=lhs, rhs=rhs, ratio=lhs / rhs,
                      reldev=abs(lhs / rhs - 1.0))
    out["D1_ok"] = 0.85 <= out["D1"]["ratio"] <= 1.15
    out["D3_ok"] = 0.90 <= out["D3"]["ratio"] <= 1.20
    out["verdict"] = "CONVICTED" if (out["D1_ok"] and out["D3_ok"]) else ""
    return out

log("base route compute ...")
V0 = route("TIn", FULL)[0]

variants = {"V0_base": V0}
for S in SUSPECTS:
    variants["Z_" + S] = route("TIn", FULL - {S})[0]
    variants["F_" + S] = route("TIn", FULL, slots_flipped={S})[0]
variants["X_Z_H7s2_T4"] = route("TIn", FULL - {"H7s2", "T4"})[0]
variants["X_F_H7s2_T4"] = route("TIn", FULL, slots_flipped={"H7s2", "T4"})[0]

# per-slot M1 contribution (interpretability): delta = full - (full minus S)
log("per-slot M1 contributions ...")
contrib = {}
for S in ALL_SLOTS:
    b_m = route("TIn", FULL - {S})[0]
    d = V0 - b_m
    contrib[S] = {k: float(d @ wH7[k]) for k in ("D1", "D2", "D3")}

rows = {}
for name, b in variants.items():
    r = m1_row(b)
    r["relL2"] = relL2(b)
    r["cos"] = cosv(b)
    rows[name] = r
    log("[%s] D1=%+.5f D2=%+.5f D3=%+.5f  relL2=%.4f cos=%.4f  %s"
        % (name, r["D1"]["ratio"], r["D2"]["ratio"], r["D3"]["ratio"],
           r["relL2"], r["cos"], r["verdict"]))

# ---------------- summary table ----------------
print("\n=== PHASE 2 VARIANT x IDENTITY (ratio = lhs/rhs; B26 rhs frozen) ===")
print("export(bexp) reference: D1=%.4f D2=%.4f D3=%.4f"
      % (B26RATIO["D1"], B26RATIO["D2"], B26RATIO["D3"]))
print("%-14s %10s %10s %10s %6s %6s %s"
      % ("variant", "D1", "D2", "D3", "D1ok", "D3ok", "verdict"))
for name, r in rows.items():
    print("%-14s %10.4f %10.4f %10.4f %6s %6s %s"
          % (name, r["D1"]["ratio"], r["D2"]["ratio"], r["D3"]["ratio"],
             "Y" if r["D1_ok"] else "n", "Y" if r["D3_ok"] else "n",
             r["verdict"]))

print("\n=== per-slot M1 contribution (delta = full - (full minus S)) ===")
print("%-7s %12s %12s %12s" % ("slot", "D1", "D2", "D3"))
for S in ALL_SLOTS:
    c = contrib[S]
    print("%-7s %+12.4e %+12.4e %+12.4e" % (S, c["D1"], c["D2"], c["D3"]))

conv = [n for n, r in rows.items() if r["verdict"] == "CONVICTED"]
print("\nCONVICTED variants: %s" % (conv if conv else "NONE -> candidate excluded"))

json.dump({"reference": {"export_D1": B26RATIO["D1"],
                         "export_D2": B26RATIO["D2"],
                         "export_D3": B26RATIO["D3"],
                         "rhs": RHS},
           "criteria": {"D1": [0.85, 1.15], "D3": [0.90, 1.20]},
           "suspects": SUSPECTS,
           "rows": rows, "per_slot_contrib": contrib},
          open(OUT + "/b29_variant_matrix.json", "w"), indent=1)
log("\nwrote %s/b29_variant_matrix.json (t=%.1fs)" % (OUT, time.time() - t0))
