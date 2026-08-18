#!/usr/bin/env python3
"""
BFINAL-005 Fixer verification: production (stageB6RxDesignOracle=false)
pressure-row activation check.

Evidence chain:
  - switch-ON  run (official case dir, diagnostic switch ON):  rxpr_prod_*.mtx
    decomposition exports + dot-test + G1-G5 anchors.
  - switch-OFF run (prod copy, switch OFF): no .mtx exports (gated), no
    dot-test (gated); the ONLY production observables are:
      (a) the written xh-level field 1/gsenshPressureDrop (momentum+pressureRow)
      (b) the written raw-design field 1/gsensPressureDrop (post-chain)
      (c) the log line  "Active gradient ranges: ... dDP=[min, max]"
      (d) dgdx[1] / gsensPressureDrop written fields

Checks:
  C1: prod written 1/gsenshPressureDrop is BYTE-IDENTICAL to the switch-ON
      run's written 1/gsenshPressureDrop  (both contain momentum+pressureRow
      because the transpose is now computed unconditionally).
  C2: prod written 1/gsensPressureDrop (raw-design, post filter_chainrule.H)
      byte-identical to switch-ON run's.
  C3: prod log "Active gradient ranges ... dDP=[...]" equals the patched
      switch-ON dDP range [-0.0277320682043, 0.0259232358427], NOT the
      momentum-only BFINAL-004 range [-0.0200726738411, 0.0194203165319].
  C4: prod log has NO "RxPressureRowTranspose:" dot-test line (gated) and NO
      "RxPressureRow production xh fields" line (gated) — expected, since only
      diagnostics are gated.
  C5: prod written 1/gsenshPressureDrop != momentum-only field. Momentum-only
      reconstruction: gM_mtx from switch-ON rxpr_prod_gsensh_momentum.mtx is
      NOT written in prod (gated), so instead compare prod written total
      against the switch-ON xh-level total mtx (should match) and note the
      known |pressureRow|L2 > 0 from the switch-ON decomposition.
"""
import os, sys, math, re

ON_DIR  = sys.argv[1]   # switch-ON  case dir (official, artifacts present)
OFF_DIR = sys.argv[2]   # switch-OFF prod copy case dir
ON_LOG  = sys.argv[3]   # switch-ON  run log
OFF_LOG = sys.argv[4]   # switch-OFF prod run log

def read_field(path):
    """Read OpenFOAM ascii volScalarField internalField list."""
    with open(path) as f:
        txt = f.read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\((.*?)\)\s*;", txt, re.S)
    if not m:
        raise SystemExit(f"{path}: internalField not parsed")
    n = int(m.group(1))
    vals = [float(x) for x in m.group(2).split()]
    if len(vals) != n:
        raise SystemExit(f"{path}: expected {n} values, got {len(vals)}")
    return vals

def read_mtx(path, n):
    vals = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("%"):
                vals.append(float(line))
    if len(vals) != n:
        raise SystemExit(f"{path}: expected {n}, got {len(vals)}")
    return vals

def log_line(log, pattern):
    with open(log) as f:
        for line in f:
            if pattern in line:
                return line.strip()
    return None

print("="*100)
print("BFINAL-005 Fixer — production activation check (switch OFF vs switch ON)")
print("="*100)

# C1/C2: byte-compare written fields (exact bytes, not values)
for name in ("gsenshPressureDrop", "gsensPressureDrop"):
    p_on  = os.path.join(ON_DIR,  "1", name)
    p_off = os.path.join(OFF_DIR, "1", name)
    if not (os.path.exists(p_on) and os.path.exists(p_off)):
        print(f"C1/C2 {name}: MISSING on={p_on} off={p_off}")
        continue
    b_on  = open(p_on, "rb").read()
    b_off = open(p_off, "rb").read()
    ident = (b_on == b_off)
    print(f"C1/C2 {name}: byte-identical switch-ON vs switch-OFF = {ident}")
    if not ident:
        # numeric compare for diagnostics
        v_on  = read_field(p_on)
        v_off = read_field(p_off)
        maxd = max(abs(a-b) for a, b in zip(v_on, v_off))
        print(f"   numeric max|diff| = {maxd:.3e}  (byte mismatch but values may match)")

# C3: dDP range comparison
on_ddp  = log_line(ON_LOG,  "Active gradient ranges")
off_ddp = log_line(OFF_LOG, "Active gradient ranges")
print(f"C3 switch-ON  dDP line: {on_ddp}")
print(f"C3 switch-OFF dDP line: {off_ddp}")
MOM_ONLY = ("[-0.0200726738411, 0.0194203165319]")   # BFINAL-003/004 momentum-only
PATCHED  = ("[-0.0277320682043, 0.0259232358427]")   # BFINAL-005 cycle-2 patched
if off_ddp:
    off_range = off_ddp.split("dDP=")[1]
    print(f"C3 prod dDP range: {off_range}  == patched? {PATCHED in off_ddp}  == momentum-only? {MOM_ONLY in off_ddp}")

# C4: gated diagnostics absent in prod log
has_dot = log_line(OFF_LOG, "RxPressureRowTranspose:") is not None
has_xh  = log_line(OFF_LOG, "RxPressureRow production xh fields") is not None
print(f"C4 prod log dot-test line present (should be False): {has_dot}")
print(f"C4 prod log xh-fields line present (should be False): {has_xh}")
on_dot = log_line(ON_LOG, "RxPressureRowTranspose:")
print(f"C4 switch-ON dot-test line (should be present): {on_dot}")

# C5: production field semantics.
# NOTE: the WRITTEN 1/gsenshPressureDrop / 1/gsensPressureDrop are POST
# filter_chainrule.H fields (the chain mutates gsenshPressureDrop in place:
# designMask masking, drho scaling, eta correction — see filter_chainrule.H),
# so they must NOT be compared to the PRE-chain .mtx decomposition exports
# (which are only written when stageB6RxDesignOracle=true, i.e. switch-ON
# runs).  The production-activation facts are:
#   (1) C1/C2: written fields byte-identical between switch-ON and switch-OFF
#       => the production computation path is IDENTICAL (pressure-row active).
#   (2) C3: production dDP range == patched total range (NOT momentum-only).
#   (3) G3-in-production (below): raw-design projection of the production
#       written gsensPressureDrop onto D1/D2/D3 == oracle total to ~1e-9.
# Pre-chain decomposition self-check (from switch-ON .mtx, verified by
# verify_gates_fixer.py): max|total-(momentum+pressurerow)| = 5.204e-18 and
# |pressureRow|L2 = 0.0145680567705 > 0.
nC = 33600
gP_mtx = read_mtx(os.path.join(ON_DIR, "rxpr_prod_gsensh_pressurerow.mtx"), nC)
nP2 = math.sqrt(sum(v*v for v in gP_mtx))
print(f"C5 switch-ON |pressureRow|L2 (pre-chain) = {nP2:.12e}  (non-zero: {nP2>0})")
print("C5 written fields are POST-chain; byte-identity of C1/C2 + dDP range C3")
print("   + G3-in-production projection are the production-activation proofs.")

# G3-in-production: project the PRODUCTION written raw-design field onto D1/D2/D3
def read_field(path):
    with open(path) as f: txt = f.read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\((.*?)\)\s*;", txt, re.S)
    n = int(m.group(1)); return [float(x) for x in m.group(2).split()]
gRaw_prod = read_field(os.path.join(OFF_DIR, "1", "gsensPressureDrop"))
dirs = read_mtx(os.path.join(ON_DIR, "stageB6_dirs.mtx"), 3*nC)
oracles = {"D1": -0.0562356333787, "D2": 0.297672431753, "D3": -0.0363414298395}
for di, name in enumerate(["D1","D2","D3"]):
    d = dirs[di*nC:(di+1)*nC]
    p = sum(gRaw_prod[i]*d[i] for i in range(nC))
    rel = abs(p-oracles[name])/abs(oracles[name])
    print(f"C5 G3-in-production {name}: proj={p:.16e} oracle={oracles[name]} relErr={rel:.3e} signOK={math.copysign(1,p)==math.copysign(1,oracles[name])}")

print("="*100)
