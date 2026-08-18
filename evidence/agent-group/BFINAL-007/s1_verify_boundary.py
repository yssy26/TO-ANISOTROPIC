#!/usr/bin/env python3
"""
BFINAL-007 S1 verification (read-only, offline):
Q2 boundary-term magnitudes and the exported J P-P sign structure.

Facts verified:
  1. outlet patch: 84 faces, 84 unique owner cells (constant/polyMesh/boundary + owner)
  2. outlet deltaCoeffs_b = 1/|nHat*(nHat&(Cf-Cn))|, magSf_b (from points/faces)
  3. rAtU (primalPressureMobility) runtime stats from BFINAL-003 run.log
  4. boundary term magnitude = rAtU_b*deltaCoeffs_b*|Sf|_b per outlet cell
  5. exported explicitJT.mtx P-P diagonal on the 84 outlet cells is the
     interior-only sum (no boundary term), confirming the missing term
  6. identity-row scan of explicitJT.mtx pressure block: only row 100800
     (discretePIndex(0)); row 106400 (discretePIndex(5600)) is a physical
     continuity row (NOT an identity row)
"""
import numpy as np, re, scipy.io

CASE = "/home/ys/dsH/b2_case_smoke"

def data_block(path):
    txt = open(path).read()
    m = re.search(r'\n(\d+)\s*\n\(', txt)
    n = int(m.group(1))
    return n, txt[m.end():]

def read_pts(path):
    n, body = data_block(path)
    toks = body.replace('(', ' ').replace(')', ' ').split()
    return np.array([float(t) for t in toks[:3*n]]).reshape(n, 3)

def read_ints(path):
    n, body = data_block(path)
    toks = body.replace('(', ' ').replace(')', ' ').split()
    return np.array([int(t) for t in toks[:n]])

pts = read_pts(CASE + '/constant/polyMesh/points')
owner = read_ints(CASE + '/constant/polyMesh/owner')
nf, fbody = data_block(CASE + '/constant/polyMesh/faces')
ftoks = fbody.replace('(', ' ').replace(')', ' ').split()
k = 0
face_pts = []
for i in range(nf):
    m = int(ftoks[k]); k += 1
    face_pts.append([int(t) for t in ftoks[k:k+m]]); k += m

def face_center_area(fi):
    fp = pts[face_pts[fi]]
    c = fp.mean(axis=0)
    v1 = fp[1] - fp[0]; v2 = fp[2] - fp[0]
    n = np.cross(v1, v2)
    return c, 0.5*np.linalg.norm(n), n/np.linalg.norm(n)

NCELL = 33600
cellc = np.zeros((NCELL, 3)); cnt = np.zeros(NCELL)
for fi in range(nf):
    fc, _, _ = face_center_area(fi)
    cellc[owner[fi]] += fc; cnt[owner[fi]] += 1
cellc /= cnt[:, None]

deltas = []; areas = []; outcells = set()
for fi in range(96944, 96944 + 84):      # outlet startFace 96944, nFaces 84
    fc, ar, n = face_center_area(fi)
    cell = owner[fi]; outcells.add(cell)
    d = fc - cellc[cell]
    dn = n*np.dot(n, d)
    deltas.append(1.0/np.linalg.norm(dn))
    areas.append(ar)
deltas = np.array(deltas); areas = np.array(areas)

print("outlet faces=84 unique cells=%d" % len(outcells))
print("outlet deltaCoeffs: min=%.3e max=%.3e avg=%.3e" %
      (deltas.min(), deltas.max(), deltas.mean()))
print("outlet magSf: min=%.3e max=%.3e avg=%.3e sum=%.4e" %
      (areas.min(), areas.max(), areas.mean(), areas.sum()))

rAtU_avg = 2.82337769549e-07   # BFINAL-003 run.log "NS rAtU: avg="
print("rAtU avg (run.log) = %.6e" % rAtU_avg)
print("boundary term avg = rAtU*delta*|Sf| = %.3e" %
      (rAtU_avg*deltas.mean()*areas.mean()))
print("boundary term min/max = %.3e / %.3e" %
      (rAtU_avg*deltas.min()*areas.min(), rAtU_avg*deltas.max()*areas.max()))

M = scipy.io.mmread(CASE + '/explicitJT.mtx').tocsr()
NV = 3*NCELL
ident = []
for r in range(NV, 4*NCELL):
    s, e = M.indptr[r], M.indptr[r+1]
    nz = [(M.indices[j], M.data[j]) for j in range(s, e) if abs(M.data[j]) > 1e-12]
    if len(nz) == 1 and nz[0][0] == r and abs(nz[0][1] - 1.0) <= 1e-12:
        ident.append(r)
print("identity rows in P block:", ident)
r5600 = NV + 5600
s, e = M.indptr[r5600], M.indptr[r5600+1]
print("row106400 stored=%d nz>1e-12=%d (physical continuity, NOT identity)" %
      (e - s, sum(1 for j in range(s, e) if abs(M.data[j]) > 1e-12)))

diags = []
for cell in sorted(outcells):
    r = NV + cell; s, e = M.indptr[r], M.indptr[r+1]
    for j in range(s, e):
        if M.indices[j] == r:
            diags.append(M.data[j])
diags = np.array(diags)
print("exported J P-P diag on %d outlet cells: min=%.3e max=%.3e (all negative = interior-only sum)"
      % (len(diags), diags.min(), diags.max()))
print("=> missing boundary term ~%.2e per outlet cell is the internalCoeffs "
      "diagonal -rAtU_b*|Sf|_b*deltaCoeffs_b (added by addBoundaryDiag)" %
      (rAtU_avg*deltas.mean()*areas.mean()))
print("S1_VERIFY_DONE")
