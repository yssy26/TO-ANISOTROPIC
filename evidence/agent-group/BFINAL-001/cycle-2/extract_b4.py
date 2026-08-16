#!/usr/bin/env python3
"""Read-only extractor for the Stage B4 diagnostic log (D1 re-anchor).
Reads a log file and prints the key B4 numbers with line numbers.
No modification of any project/case artifact.
"""
import re, sys

TARGETS = [
    ("g_w thermal/pressure", r"Discrete objective derivative errors: (.*)"),
    ("Rx-A", r"RxProbe Rx-A.*relL2=([0-9.eE+-]+)"),
    ("Rx-B raw", r"RxProbe Rx-B raw:.*\|dRp/dalpha\|L2=([0-9.eE+-]+)"),
    ("Rx-B final", r"RxProbe Rx-B \(fixed U,p.*\|dRp/dalpha\|L2=([0-9.eE+-]+)"),
    ("T1", r"Stage B4 T1 .*relL2=([0-9.eE+-]+) cos=([0-9.eE+-]+)"),
    ("T-cont", r"Stage B4 T-cont .*relL2 = ([0-9.eE+-]+)"),
    ("T2", r"Stage B4 T2 .*relL2\(J_P, div\(dphi_SIMPLE\)\)=([0-9.eE+-]+)"),
    ("T-fvm", r"Stage B4 T-fvm .*relL2 = ([0-9.eE+-]+)"),
    ("T-dev", r"Stage B4 T-dev .*relL2 = ([0-9.eE+-]+)"),
    ("T-conv", r"Stage B4 T-conv .*relL2 = ([0-9.eE+-]+)"),
    ("GateH1 same-basis dH", r"GateH1 same-basis dH:.*relL2=([0-9.eE+-]+) cos=([0-9.eE+-]+)"),
]

def main(path):
    lines = open(path, errors="replace").read().splitlines()
    found = {}
    for i, ln in enumerate(lines, 1):
        for name, pat in TARGETS:
            m = re.search(pat, ln)
            if m:
                found.setdefault(name, []).append((i, m.groups()))
    for name, hits in found.items():
        for i, g in hits:
            print(f"L{i:5d}  {name:24s}  {g}")
    # sanity: also count lines
    print(f"\n[total lines: {len(lines)}]")

if __name__ == "__main__":
    main(sys.argv[1])
