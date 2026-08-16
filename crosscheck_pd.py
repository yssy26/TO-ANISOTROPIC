#!/usr/bin/env python3
"""
Cross-check: project the VOLUME gradient (gsensVol, whose chain is known-correct
to relV~1e-6) onto D1 and compare with the FD value 0.155512 (from the mini-FD
baseline).  If this matches, the projection script is correct and the pressure
mismatch is a real operator bug.
"""
import re
import numpy as np

BASE = r"\\wsl.localhost\Ubuntu-20.04\home\ys\b2_case_smoke"

def read_internal_field(path):
    with open(path, "r") as f:
        txt = f.read()
    m = re.search(r'internalField\s+(?:nonuniform\s+)?List<[^>]+>\s*\n?\s*(\d+)\s*\n?\s*\((.*?)\)\s*;', txt, re.S)
    n = int(m.group(1))
    vals = re.findall(r'[-+0-9.eE]+', m.group(2))
    return np.array([float(v) for v in vals])

def read_vector_field(path):
    with open(path, "r") as f:
        txt = f.read()
    m = re.search(r'internalField\s+(?:nonuniform\s+)?List<[^>]+>\s*\n?\s*(\d+)\s*\n?\s*\((.*?)\)\s*;', txt, re.S)
    n = int(m.group(1))
    triples = re.findall(r'\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)', m.group(2))
    return np.array([[float(a), float(b), float(c)] for a, b, c in triples])

gsensV = read_internal_field(BASE + r"\1\gsensVol")
gsensPD = read_internal_field(BASE + r"\1\gsensPressureDrop")
mask = read_internal_field(BASE + r"\1\designMask")
x = read_internal_field(BASE + r"\1\x")
C = read_vector_field(BASE + r"\1\C")
N = gsensV.size

Lx = C[:,0].max() - C[:,0].min()
Ly = C[:,1].max() - C[:,1].min()
Lz = C[:,2].max() - C[:,2].min()
xx = C[:,0]/Lx; yy = C[:,1]/Ly; zz = C[:,2]/Lz

D1 = np.sin(2*np.pi*xx)*np.cos(np.pi*yy) + 0.5*np.sin(3*np.pi*zz)
active = mask > 0.5
D1[~active] = 0.0
D1[(x < 0.02) | (x > 0.98)] = 0.0
D1 /= np.max(np.abs(D1))

projV = np.sum(gsensV[active] * D1[active])
projPD = np.sum(gsensPD[active] * D1[active])
print("projV_D1  = %.12e  (FD_gV_D1 = 0.1555125114, expect match ~1e-6)" % projV)
print("projPD_D1 = %.12e  (FD_gDP_D1 = -2.53708557)" % projPD)
print("gsensV range:", gsensV.min(), gsensV.max())
print("gsensPD range:", gsensPD.min(), gsensPD.max())
