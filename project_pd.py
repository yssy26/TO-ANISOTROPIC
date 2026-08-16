#!/usr/bin/env python3
"""
Project the (exact) pressure-drop gradient gsensPressureDrop onto D1/D2/D3 and
compare with the FD values from the mini-FD baseline run.

Reads OpenFOAM ascii fields (gsensPressureDrop, designMask, x) and cell centres
C (written via `postProcess -func writeCellCentres`). Reconstructs D1/D2/D3 with
the exact formulas from validateStageB2GradientAmplitude.H, then projects over
active design cells (designMask > 0.5).
"""
import re
import sys
import numpy as np

BASE = r"\\wsl.localhost\Ubuntu-20.04\home\ys\b2_case_smoke"

def read_internal_field(path):
    """Read the internalField of an OpenFOAM ascii field file."""
    with open(path, "r") as f:
        txt = f.read()
    # locate "internalField"
    m = re.search(r'internalField\s+(?:nonuniform\s+)?List<[^>]+>\s*\n?\s*(\d+)\s*\n?\s*\((.*?)\)\s*;', txt, re.S)
    if not m:
        raise RuntimeError("cannot parse internalField in " + path)
    n = int(m.group(1))
    body = m.group(2)
    # scalars: whitespace/line separated numbers
    vals = re.findall(r'[-+0-9.eE]+', body)
    if len(vals) != n:
        raise RuntimeError("expected %d values, got %d in %s" % (n, len(vals), path))
    return np.array([float(v) for v in vals])

def read_vector_field(path):
    """Read a vector field's internalField (N vectors of 3 comps)."""
    with open(path, "r") as f:
        txt = f.read()
    m = re.search(r'internalField\s+(?:nonuniform\s+)?List<[^>]+>\s*\n?\s*(\d+)\s*\n?\s*\((.*?)\)\s*;', txt, re.S)
    if not m:
        raise RuntimeError("cannot parse internalField in " + path)
    n = int(m.group(1))
    body = m.group(2)
    # vectors: (a b c)
    triples = re.findall(r'\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)', body)
    if len(triples) != n:
        raise RuntimeError("expected %d vectors, got %d in %s" % (n, len(triples), path))
    return np.array([[float(a), float(b), float(c)] for a, b, c in triples])

gsens = read_internal_field(BASE + r"\1\gsensPressureDrop")
mask = read_internal_field(BASE + r"\1\designMask")
x = read_internal_field(BASE + r"\1\x")
C = read_vector_field(BASE + r"\1\C")

N = gsens.size
assert mask.size == N and x.size == N and C.shape[0] == N, "size mismatch"

Lx = C[:, 0].max() - C[:, 0].min()
Ly = C[:, 1].max() - C[:, 1].min()
Lz = C[:, 2].max() - C[:, 2].min()
xx = C[:, 0] / Lx
yy = C[:, 1] / Ly
zz = C[:, 2] / Lz

D = {}
D["D1"] = np.sin(2*np.pi*xx)*np.cos(np.pi*yy) + 0.5*np.sin(3*np.pi*zz)
D["D2"] = np.sin(4*np.pi*xx) + 0.3*np.cos(2*np.pi*xx)
D["D3"] = np.cos(2*np.pi*yy)*np.sin(2*np.pi*zz)

active = mask > 0.5
for k in D:
    d = D[k].copy()
    d[~active] = 0.0
    d[(x < 0.02) | (x > 0.98)] = 0.0
    mx = np.max(np.abs(d))
    if mx > 1e-30:
        d = d / mx
    proj = float(np.sum(gsens[active] * d[active]))
    print("%s projDP = %.12e" % (k, proj))

print("--- FD reference (from mini-FD baseline) ---")
print("D1 FD_gDP = -2.537")
print("D2 FD_gDP = +0.523")
print("D3 FD_gDP = +6.477")
