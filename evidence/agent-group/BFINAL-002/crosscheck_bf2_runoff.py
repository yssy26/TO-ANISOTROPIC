#!/usr/bin/env python3
"""Read-only Python B-F2 partial cross-check using exported artifacts.

DIAGNOSTIC_ONLY (cycle-2, D4). Reads ONLY existing evidence artifacts:
  /home/ys/b2_case_smoke/explicitJT.mtx             (complete physical J^T, 134400x134400, pRef row pinned)
  /home/ys/b2_case_smoke/stageB2_dUdirPhys.mtx      (direction dUdir, 3N = 100800)
  /home/ys/b2_case_smoke/stageB2_ruDudU.mtx         (frozen-phi full-residual momentum FD along dUdir)
  /home/ys/b2_case_smoke/stageB2_rpDudU_rebuilt.mtx (rebuilt-phi continuity FD along dUdir)
  /home/ys/b2_case_smoke/stageB2_Jdq_P_forward.mtx  (J's matrix-free P-row on the same direction)
Writes nothing into the case; prints a compact report.

Convention (verified from src/solveDiscreteFlowAdjoint.H L101-108 and the
export code L2896-2919): U index = 3*cell + cmpt, P index = 3N + cell,
N = 33600 cells, 4N = 134400 unknowns. explicitJT is the PHYSICAL J^T;
J*v is computed by scattering: for each J^T entry (r,c,val),
out[c] += val * v[r]   (because J[c,r] = J^T[r,c]).
"""
import sys, time

CASE = "/home/ys/dsH/b2_case_smoke/"
N = 33600
NU = 3 * N          # 100800
TOT = 4 * N         # 134400

def load_vec(name, count):
    with open(CASE + name, "r") as f:
        vals = [float(x) for x in f.read().split()]
    if len(vals) != count:
        raise SystemExit(f"{name}: expected {count} values, got {len(vals)}")
    return vals

def rel_l2(a, b, lo, hi, label):
    n = 0.0; na = 0.0; nb = 0.0; dot = 0.0; mx = 0.0; arg = -1
    for i in range(lo, hi):
        d = a[i] - b[i]
        n += d * d; na += a[i] * a[i]; nb += b[i] * b[i]; dot += a[i] * b[i]
        if abs(d) > mx: mx = abs(d); arg = i
    r = (n / nb) ** 0.5 if nb > 0 else float("nan")
    cos = dot / ((na * nb) ** 0.5) if (na > 0 and nb > 0) else float("nan")
    print(f"{label}: relL2={r:.6e} cos={cos:.9f} |ref|L2={nb**0.5:.6e} "
          f"|Jv|L2={na**0.5:.6e} maxAbsDiff={mx:.6e} @idx={arg}")

t0 = time.time()
print("loading direction + FD vectors ...", flush=True)
dU = load_vec("stageB2_dUdirPhys.mtx", NU)
ru = load_vec("stageB2_ruDudU.mtx", NU)
rp = load_vec("stageB2_rpDudU_rebuilt.mtx", N)
jq = load_vec("stageB2_Jdq_P_forward.mtx", N)
print(f"loaded in {time.time()-t0:.1f}s", flush=True)

# full 4N state vector: v = [dU; 0]
v = [0.0] * TOT
for i in range(NU):
    v[i] = dU[i]

out = [0.0] * TOT
t1 = time.time()
nrows = 0; nnz = 0
with open(CASE + "explicitJT.mtx", "r") as f:
    # skip 2 header lines
    f.readline(); f.readline()
    for line in f:
        p = line.split()
        if len(p) < 3:
            continue
        r = int(p[0]) - 1
        c = int(p[1]) - 1
        val = float(p[2])
        if val != 0.0:
            out[c] += val * v[r]
        nnz += 1
        nrows += 1
print(f"scatter matvec over {nnz} nnz in {time.time()-t1:.1f}s", flush=True)

outU = out[:NU]
outP = out[NU:]

print("\n--- B-F2 partial cross-check (direction dUdir = sin(0.271*(celli+1)), single hU) ---")
rel_l2(outU, ru, 0, NU, "momentum rows:   J*v  vs frozen-phi FD (ruDudU)")
rel_l2(outP, rp, 0, N,  "continuity rows: J*v  vs rebuilt-phi FD (rpDudU_rebuilt)")
rel_l2(outP, jq, 0, N,  "continuity rows: J*v  vs J matrix-free P-row (Jdq_P_forward)")

# P-row comparison excluding the pinned pRef row (explicitJT row P(pRef) -> identity)
print("\n--- continuity P-row, pRef row excluded ---")
rel_l2([outP[i] if i != 0 else 0.0 for i in range(N)],
       [rp[i] if i != 0 else 0.0 for i in range(N)], 0, N,
       "continuity rows: J*v  vs rebuilt-phi FD (excl pRef)")

# magnitude context: |outU|, |ru|, |outP|, |rp|, |jq|
def l2(x):
    return sum(a*a for a in x) ** 0.5
print(f"\n|Jv_U|L2={l2(outU):.6e}  |ruDudU|L2={l2(ru):.6e}  "
      f"|Jv_P|L2={l2(outP):.6e}  |rpDudU_rebuilt|L2={l2(rp):.6e}  "
      f"|Jdq_P_forward|L2={l2(jq):.6e}")
print(f"dUdir max|.|={max(abs(x) for x in dU):.6e}")
print(f"total wall {time.time()-t0:.1f}s")
