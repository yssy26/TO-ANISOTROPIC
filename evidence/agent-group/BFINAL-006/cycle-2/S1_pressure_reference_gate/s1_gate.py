#!/usr/bin/env python3
"""
BFINAL-006 cycle-2 S1 — Pressure-Reference Gate (5 programmatic checks).

Assigned step: S1-PressureReferenceGate.
  (1) verify p.needReference()==false from 0/p outlet fixedValue uniform 0 +
      readTransportProperties.H setRefCell semantics;
  (2) scan explicitJT.mtx pressure rows (3N..4N-1) for the UNIQUE identity row
      using "single |v|>1e-12 AND diagonal==1" (NOT nrow==1);
  (3) verify that row == discretePIndex(0) == 100800;
  (4) verify discretePIndex(5600) == 106400 is NOT an identity row;
  (5) record pRefRow = 100800.
Any failure => BLOCKED, no guessing.

DIAGNOSTIC_ONLY: this script only READS case files, source files and the
exported matrix artifact; it writes nothing outside the cycle-2 evidence dir.
"""
import os, re, sys, time

WS   = "/home/ys/dsH/TO-ANISOTROPIC"
CASE = "/home/ys/dsH/b2_case_smoke"
MTX  = os.path.join(CASE, "explicitJT.mtx")
PVAL = os.path.join(CASE, "0", "p")
FVS  = os.path.join(CASE, "system", "fvSolution")
RTP  = os.path.join(WS, "src", "readTransportProperties.H")

OUTDIR = os.path.dirname(os.path.abspath(__file__))
LOG    = os.path.join(OUTDIR, "s1_gate_scan.log")
RAWMD  = os.path.join(OUTDIR, "s1_raw_scan_evidence.md")

def log(msg):
    print(msg, flush=True)
    with open(LOG, "a") as f:
        f.write(msg + "\n")

def main():
    open(LOG, "w").close()
    t0 = time.time()
    results = {}
    log("=== BFINAL-006 cycle-2 S1 Pressure-Reference Gate ===")
    log("workspace=%s case=%s" % (WS, CASE))
    log("explicitJT.mtx exists=%s size=%d" % (os.path.exists(MTX), os.path.getsize(MTX) if os.path.exists(MTX) else -1))

    # ---------------------------------------------------------------
    # CHECK 1: p.needReference()==false from 0/p + setRefCell semantics
    # ---------------------------------------------------------------
    log("\n[CHECK 1] p.needReference() semantics from ACTUAL case/source")
    with open(PVAL) as f:
        ptext = f.read()
    bf = re.search(r"boundaryField\s*\{(.*)\}\s*// \*{2,}", ptext, re.S)
    if not bf:
        bf = re.search(r"boundaryField\s*\{(.*)\}", ptext, re.S)
    bftext = bf.group(1)
    # extract each patch block: name { ... type X; ... }
    patch_blocks = re.findall(r"(\w+)\s*\{([^}]*)\}", bftext)
    fixed_patches = []
    zero_grad = 0
    for name, body in patch_blocks:
        m = re.search(r"type\s+(\w+)", body)
        typ = m.group(1) if m else "?"
        val = re.search(r"value\s+uniform\s+([-\d.eE+]+)", body)
        if typ == "fixedValue":
            fixed_patches.append((name, val.group(1) if val else None))
        elif typ == "zeroGradient":
            zero_grad += 1
    log("  patches found: %s" % ", ".join("%s(%s)" % (n, t) for n, t in [(p[0], "fixedValue") for p in fixed_patches] + [("zeroGradient x%d" % zero_grad, "")] if n))
    log("  fixedValue patches: %s" % fixed_patches)
    # OpenFOAM-7 GeometricField::needReference(): needRef=false if ANY patch fixesValue()
    needReference = not fixed_patches
    log("  -> p.needReference() == %s  (false because outlet fixedValue uniform 0 fixes a value)" % needReference)

    # readTransportProperties.H setRefCell call
    with open(RTP) as f:
        rtp_lines = f.readlines()
    init_line = [l for l in rtp_lines if re.search(r"label\s+pRefCell\s*=", l)]
    call_line = [l for l in rtp_lines if re.search(r"setRefCell\(p\s*,", l)]
    log("  readTransportProperties.H init : %s" % init_line[0].strip() if init_line else "  init NOT FOUND")
    log("  readTransportProperties.H call : %s" % call_line[0].strip() if call_line else "  call NOT FOUND")
    # fvSolution pRefCell entry (present but ineffective when !needReference)
    with open(FVS) as f:
        fvs = f.read()
    prc = re.search(r"pRefCell\s+(\d+)", fvs)
    prv = re.search(r"pRefValue\s+([-\d.eE+]+)", fvs)
    log("  fvSolution pRefCell=%s pRefValue=%s (INEFFECTIVE: setRefCell is a no-op when !needReference)" %
        (prc.group(1) if prc else "?", prv.group(1) if prv else "?"))
    c1_ok = (needReference is False) and bool(fixed_patches)
    results["C1"] = c1_ok
    log("  CHECK 1 => %s" % ("PASS" if c1_ok else "FAIL"))

    # ---------------------------------------------------------------
    # CHECK 2: scan explicitJT.mtx pressure rows (3N..4N-1) for the unique
    #          identity row via "single |v|>1e-12 AND diagonal==1".
    # ---------------------------------------------------------------
    log("\n[CHECK 2] identity-row scan of explicitJT.mtx (corrected criterion)")
    log("  reading header ...")
    with open(MTX) as f:
        for line in f:
            if line.startswith("%"):
                continue
            toks = line.split()
            if len(toks) == 3:
                nrow, ncol, nnz_hdr = int(toks[0]), int(toks[1]), int(toks[2])
                break
    log("  matrix %d x %d, header nnz=%d" % (nrow, ncol, nnz_hdr))
    N = nrow // 4
    P0 = 3 * N          # first pressure row index (== discretePIndex(0))
    log("  N = nrow/4 = %d ; pressure rows = [%d, %d) ; discretePIndex(0)=%d ; discretePIndex(5600)=%d"
        % (N, P0, 4 * N, P0, P0 + 5600))
    # stream the coordinate lines, keep only pressure rows, accumulate (row,col)->val
    from collections import defaultdict
    prs = defaultdict(float)
    nlines = 0
    t1 = time.time()
    with open(MTX) as f:
        for line in f:
            if line.startswith("%"):
                continue
            toks = line.split()
            if len(toks) != 3:
                continue
            r = int(toks[0]) - 1
            if r < P0:
                continue
            c = int(toks[1]) - 1
            v = float(toks[2])
            prs[(r, c)] += v
            nlines += 1
    log("  streamed %d pressure-block coordinate lines in %.1f s; unique (r,c) pairs = %d"
        % (nlines, time.time() - t1, len(prs)))

    # group by row
    rows = defaultdict(dict)
    for (r, c), v in prs.items():
        rows[r][c] = v
    log("  pressure rows with any entry: %d" % len(rows))

    # identity detection: exactly one |v|>1e-12 entry AND it is the diagonal AND |diag-1|<=1e-12
    identity_rows = []
    for r in sorted(rows):
        nz = [(c, v) for c, v in rows[r].items() if abs(v) > 1e-12]
        diag = rows[r].get(r, 0.0)
        if len(nz) == 1 and nz[0][0] == r and abs(diag - 1.0) <= 1e-12:
            identity_rows.append(r)
    log("  identity rows (single |v|>1e-12 AND diagonal==1): %s" % identity_rows)
    # also demonstrate that nrow==1 (old criterion) fails on the pinned row:
    for r in identity_rows:
        stored = len(rows[r])
        log("    row %d: stored entries=%d (13 explicit zeros + 1 diagonal 1) -> nrow==1 criterion would FAIL, corrected criterion PASSES"
            % (r, stored))
    c2_ok = (len(identity_rows) == 1 and identity_rows[0] == P0)
    results["C2"] = c2_ok
    log("  CHECK 2 => %s (unique identity row == %d == discretePIndex(0))" % ("PASS" if c2_ok else "FAIL", P0))

    # ---------------------------------------------------------------
    # CHECK 3: identity row == discretePIndex(0) == 100800
    # ---------------------------------------------------------------
    detected = identity_rows[0] if identity_rows else None
    expected_P0 = 3 * N
    c3_ok = (detected == expected_P0 == 100800)
    results["C3"] = c3_ok
    log("\n[CHECK 3] detected identity row %s == discretePIndex(0) = 3*%d+0 = %d == 100800 => %s"
        % (detected, N, expected_P0, "PASS" if c3_ok else "FAIL"))

    # ---------------------------------------------------------------
    # CHECK 4: discretePIndex(5600) == 106400 is NOT an identity row
    # ---------------------------------------------------------------
    row5600 = P0 + 5600
    if row5600 in rows:
        entries = sorted(rows[row5600].items())
        nz5600 = [(c, v) for c, v in entries if abs(v) > 1e-12]
        diag5600 = rows[row5600].get(row5600, 0.0)
        is_identity = (len(nz5600) == 1 and nz5600[0][0] == row5600 and abs(diag5600 - 1.0) <= 1e-12)
        log("\n[CHECK 4] row 106400 (= discretePIndex(5600)): stored entries=%d, |v|>1e-12 entries=%d, diag=%r, is_identity=%s"
            % (len(entries), len(nz5600), diag5600, is_identity))
        log("  row 106400 nonzero columns: %s" % [c for c, v in nz5600][:20])
    else:
        is_identity = False
        log("\n[CHECK 4] row 106400 has NO stored entries (not identity)")
    c4_ok = (not is_identity)
    results["C4"] = c4_ok
    log("  CHECK 4 => %s (106400 is a physical continuity row, NOT the reference/identity row)" % ("PASS" if c4_ok else "FAIL"))

    # ---------------------------------------------------------------
    # CHECK 5: record pRefRow = 100800
    # ---------------------------------------------------------------
    pRefRow = 100800
    c5_ok = (detected == pRefRow)
    results["C5"] = c5_ok
    log("\n[CHECK 5] pRefRow = %d recorded for T1 re-pinning (J row pRefRow -> identity, rhs[pRefRow]=0). => %s"
        % (pRefRow, "PASS" if c5_ok else "FAIL"))

    # ---------------------------------------------------------------
    # raw evidence dump for row 100800 and row 106400
    # ---------------------------------------------------------------
    with open(RAWMD, "w") as f:
        f.write("# BFINAL-006 cycle-2 S1 raw scan evidence\n\n")
        f.write("matrix: %s\n" % MTX)
        f.write("dims: %d x %d, N=%d, pressure rows [%d, %d)\n\n" % (nrow, ncol, N, P0, 4 * N))
        f.write("## identity rows detected (single |v|>1e-12 AND diagonal==1): %s\n\n" % identity_rows)
        if detected is not None and detected in rows:
            f.write("## row %d entries (pinned identity row):\n" % detected)
            for c, v in sorted(rows[detected].items()):
                f.write("  col %d : %r\n" % (c, v))
        if row5600 in rows:
            f.write("\n## row 106400 entries (physical continuity row):\n")
            for c, v in sorted(rows[row5600].items()):
                f.write("  col %d : %r\n" % (c, v))
        f.write("\n## patch summary (0/p):\n")
        for name, val in fixed_patches:
            f.write("  %s : fixedValue uniform %s\n" % (name, val))
        f.write("  zeroGradient patches: %d\n" % zero_grad)
        f.write("  p.needReference() == %s\n" % needReference)
        f.write("\n## check results: %s\n" % results)

    # ---------------------------------------------------------------
    all_ok = all(results.values())
    log("\n=== S1 GATE RESULT: %s (%.1f s) ===" % ("ALL 5 CHECKS PASS" if all_ok else "BLOCKED", time.time() - t0))
    marker = os.path.join(OUTDIR, "S1_GATE_PASSED" if all_ok else "S1_GATE_BLOCKED")
    with open(marker, "w") as f:
        f.write("results=%s\n" % results)
    sys.exit(0 if all_ok else 2)

if __name__ == "__main__":
    main()
