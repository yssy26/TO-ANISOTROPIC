#!/usr/bin/env python3
"""
BFINAL-016 B1 - thermal-mediated quantification for the J objective.

Identity under test (task spec):
    FD_J = flow-mediated (+/- b_TC^T w_true) + thermal-mediated
where the thermal-mediated term is computed INDEPENDENTLY as
    thermal_D = sum_c (dJ/dT)[c] * (T+ - T-)/(2h)
from the exported T states and the P1b-verified dJ/dT source.

Criterion: (flow + thermal) explains FD_J to >= 90% on all directions.
Both sign conventions of the flow-mediated term are reported; the
pre-registration of this gauge is the task book itself.
"""
import os, json
import numpy as np

STATE = "/home/ys/dsH/b16_states/stageB2"
EXPORT_CASE = "/home/ys/dsH/b15_export"
CYCLE = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-016/cycle-1"

N = 33600; NV = 3*N; NUNK = 4*N; P0 = NV
FD_J = {"D1": 0.0428559144263, "D2": -0.00346892642156, "D3": -0.000282658160322}

import scipy.io

def load_rhs_export(path):
    v = np.asarray(scipy.io.mmread(path)).ravel()
    m = v.reshape(N, 4)
    b = np.zeros(NUNK)
    b[0:NV:3] = m[:, 0]; b[1:NV:3] = m[:, 1]; b[2:NV:3] = m[:, 2]
    b[P0:P0 + N] = m[:, 3]
    return b

bTC = load_rhs_export(EXPORT_CASE + "/explicitRhs_thermalCoupling.mtx")
dJdT = np.loadtxt(STATE + "/wstate_dJdT.mtx")
Mflow, Tref = np.loadtxt(STATE + "/wstate_outletMeta.mtx")
print("outlet M=%.6e Tref=%.6f  |dJdT|max=%.3e nz=%d" %
      (Mflow, Tref, np.abs(dJdT).max(), int((dJdT != 0).sum())))

def load_state(name, tag, sign, field):
    return np.loadtxt(STATE + "/wstate_%s_h%s_%s_%s.mtx" % (name, tag, sign, field))

hmap = {"0.0003": 3e-4, "0.001": 1e-3}
results = {}
for name in ["D1", "D2", "D3"]:
    tags = sorted(set(
        f.split("_h")[1].split("_")[0]
        for f in os.listdir(STATE) if f.startswith("wstate_%s_h" % name) and f.endswith("_T.mtx")))
    for tag in tags:
        h = hmap.get(tag)
        if h is None: continue
        Up = load_state(name, tag, "p", "U"); Um = load_state(name, tag, "m", "U")
        pp = load_state(name, tag, "p", "p"); pm = load_state(name, tag, "m", "p")
        Tp = load_state(name, tag, "p", "T"); Tm = load_state(name, tag, "m", "T")
        w_true = np.zeros(NUNK)
        w_true[0:NV:3] = ((Up - Um)/(2.0*h))[:, 0]
        w_true[1:NV:3] = ((Up - Um)/(2.0*h))[:, 1]
        w_true[2:NV:3] = ((Up - Um)/(2.0*h))[:, 2]
        w_true[P0:] = (pp - pm)/(2.0*h)
        flow = float(bTC @ w_true)
        thermal = float(dJdT @ ((Tp - Tm)/(2.0*h)))
        fd = FD_J[name]
        cand_plus = flow + thermal
        cand_minus = -flow + thermal
        print("B1 %s h=%s: flow(bTC^T w_true)=%.6e thermal=%.6e | "
              "+flow+thermal=%.6e (vs FD %.6e, err %.2e) | "
              "-flow+thermal=%.6e (err %.2e)" %
              (name, tag, flow, thermal,
               cand_plus, fd, abs(cand_plus-fd)/abs(fd),
               cand_minus, abs(cand_minus-fd)/abs(fd)))
        results["%s_h%s" % (name, tag)] = {
            "flow": flow, "thermal": thermal, "fd": fd,
            "err_plus": abs(cand_plus-fd)/abs(fd),
            "err_minus": abs(cand_minus-fd)/abs(fd)}

with open(CYCLE + "/b1_thermal.json", "w") as f:
    json.dump(results, f, indent=1)
print("B1_DONE")
