#!/usr/bin/env python3
"""BFINAL-029 Phase 2 completeness: full 2-toggle cross scan.

Task book (NEXT_TASK_BFINAL029.md) Phase 2: "另测双变体交叉（如 H7s2+T4）" —
"如" = such-as, so the double-cross is not limited to H7s2xT4.  To make the
"candidate excluded" verdict airtight we enumerate the COMPLETE 2-toggle space
over the suspect set {H7s2, T4, T3, T6, T2, T1} plus H7s2b (H7-family slot
that participates in the D3 near-cancellation H7s2+0.274 / H7s2b-0.260):

  Z_i+Z_j : route("TIn", FULL - {i,j})
  F_i+F_j : route("TIn", FULL, slots_flipped={i,j})

Same binding criteria as b29_PREREGISTRATION.md 4.2:
  CONVICTED <=> D1 ratio in [0.85,1.15] AND D3 ratio in [0.90,1.20].
w_true and rhs frozen from B26.  Complete enumeration (no selective hunting);
all 42 variants reported, verdicts per variant.
"""
import time, json, itertools, importlib.util
import numpy as np

OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-029/cycle-1"
E26 = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-026/cycle-1"
B25 = "/home/ys/dsH/b25_qgate"

t0 = time.time()
log = lambda m: print(m, flush=True)
log("=== BFINAL-029 Phase 2 complete 2-toggle cross scan ===")

spec = importlib.util.spec_from_file_location(
    "b29_slot_instrument", OUT + "/b29_slot_instrument.py")
inst = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inst)

N = inst.N; NV = inst.NV; NUNK = inst.NUNK
route = inst.route
FULL = inst.FULL

m_exp = np.loadtxt(B25 + "/b18rhs_thermalCoupling.mtx")
bexp = np.zeros(NUNK)
bexp[0:NV:3] = m_exp[:, 0]; bexp[1:NV:3] = m_exp[:, 1]
bexp[2:NV:3] = m_exp[:, 2]; bexp[NV:] = m_exp[:, 3]

wH7 = np.load(E26 + "/b26_wstar_h7.npz")
res26 = json.load(open(E26 + "/b26a_m1_results.json"))
RHS = {r["dir"]: r["rhs"] for r in res26["rows"]}

SCAN = ["H7s2", "T4", "T3", "T6", "T2", "T1", "H7s2b"]  # suspects + H7-family
PAIRS = list(itertools.combinations(SCAN, 2))

def m1_ok(b):
    ratio = {k: float(b @ wH7[k]) / RHS[k] for k in ("D1", "D2", "D3")}
    d1 = 0.85 <= ratio["D1"] <= 1.15
    d3 = 0.90 <= ratio["D3"] <= 1.20
    return ratio, d1, d3

rows = {}
convicted = []
log("scanning %d pairs x {Z, F} ..." % len(PAIRS))
for i, j in PAIRS:
    for form in ("Z", "F"):
        name = "%s_%s+%s" % (form, i, j)
        if form == "Z":
            b = route("TIn", FULL - {i, j})[0]
        else:
            b = route("TIn", FULL, slots_flipped={i, j})[0]
        ratio, d1, d3 = m1_ok(b)
        rows[name] = {"D1": ratio["D1"], "D2": ratio["D2"],
                      "D3": ratio["D3"], "D1_ok": d1, "D3_ok": d3,
                      "verdict": "CONVICTED" if (d1 and d3) else ""}
        if d1 and d3:
            convicted.append(name)
        log("[%s] D1=%+.4f D3=%+.4f  %s"
            % (name, ratio["D1"], ratio["D3"], rows[name]["verdict"]))

print("\n=== 2-TOGGLE CROSS SCAN (42 variants) ===")
print("%-16s %10s %10s %10s %6s %6s" % ("variant", "D1", "D2", "D3", "D1ok", "D3ok"))
for nm, r in rows.items():
    print("%-16s %10.4f %10.4f %10.4f %6s %6s %s"
          % (nm, r["D1"], r["D2"], r["D3"], "Y" if r["D1_ok"] else "n",
             "Y" if r["D3_ok"] else "n", r["verdict"]))
print("\nCONVICTED: %s" % (convicted if convicted else "NONE -> candidate excluded"))

json.dump({"criteria": {"D1": [0.85, 1.15], "D3": [0.90, 1.20]},
           "scan_set": SCAN, "rows": rows, "convicted": convicted},
          open(OUT + "/b29_crossscan.json", "w"), indent=1)
log("wrote %s/b29_crossscan.json (t=%.1fs)" % (OUT, time.time() - t0))
