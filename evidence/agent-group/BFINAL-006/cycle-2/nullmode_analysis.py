#!/usr/bin/env python3
"""Correlate the T1 null-mode with alpha and cell geometry (read-only diagnostic)."""
import numpy as np, re

def _data_block(path):
    txt = open(path).read()
    m = re.search(r'\n(\d+)\s*\n\(', txt)
    n = int(m.group(1))
    return n, txt[m.end():]

def read_pts(path):
    n, body = _data_block(path)
    toks = body.replace('(', ' ').replace(')', ' ').split()
    return np.array([float(t) for t in toks[:3 * n]]).reshape(n, 3)

def read_ints(path):
    n, body = _data_block(path)
    toks = body.replace('(', ' ').replace(')', ' ').split()
    return np.array([int(t) for t in toks[:n]])

pts = read_pts('/home/ys/dsH/b2_case_smoke/constant/polyMesh/points')
owner = read_ints('/home/ys/dsH/b2_case_smoke/constant/polyMesh/owner')
neigh = read_ints('/home/ys/dsH/b2_case_smoke/constant/polyMesh/neighbour')
nf, fbody = _data_block('/home/ys/dsH/b2_case_smoke/constant/polyMesh/faces')
ftoks = fbody.replace('(', ' ').replace(')', ' ').split()
k = 0
face_pts = []
for i in range(nf):
    m = int(ftoks[k]); k += 1
    face_pts.append([int(t) for t in ftoks[k:k + m]]); k += m
NCELL = 33600
cent = np.zeros((NCELL, 3)); count = np.zeros(NCELL)
for fi in range(nf):
    fc = pts[face_pts[fi]].mean(axis=0)
    cent[owner[fi]] += fc; count[owner[fi]] += 1
    if fi < len(neigh):
        cent[neigh[fi]] += fc; count[neigh[fi]] += 1
cent /= count[:, None]
print('cell centers: x %.4f..%.4f y %.4f..%.4f z %.4f..%.4f' % (
    cent[:, 0].min(), cent[:, 0].max(), cent[:, 1].min(), cent[:, 1].max(),
    cent[:, 2].min(), cent[:, 2].max()))

def read_internal_scalar(path):
    txt = open(path).read()
    m = re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(([^)]*)\)', txt)
    return np.array([float(t) for t in m.group(2).replace('\n', ' ').split()])

alpha = read_internal_scalar('/home/ys/dsH/b2_case_smoke/0/alpha')
print('alpha: min %.4e max %.4e mean %.4e  #(alpha>1e-3)=%d #(alpha>0.5)=%d' % (
    alpha.min(), alpha.max(), alpha.mean(), (alpha > 1e-3).sum(), (alpha > 0.5).sum()))

w = np.loadtxt('/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-006/cycle-2/artifacts/wPrime_TAN_D1.mtx', dtype=np.float64)
u = w[:3 * NCELL]; p = w[3 * NCELL:]
ucm = np.linalg.norm(u.reshape(-1, 3), axis=1)
print('null-mode |U| per cell: max %.3e' % ucm.max())
print('null-mode P: max %.3e min %.3e' % (p.max(), p.min()))
corr = np.corrcoef(np.abs(p), alpha)[0, 1]
print('corr(|p_null|, alpha) = %.4f' % corr)
thr = np.percentile(np.abs(p), 99)
big = np.abs(p) > thr
print('|p|>p99: n=%d  mean alpha(big)=%.4e  mean alpha(rest)=%.4e' % (
    big.sum(), alpha[big].mean(), alpha[~big].mean()))
c = cent
print('big-P x %.3f..%.3f y %.3f..%.3f z %.3f..%.3f' % (
    c[big, 0].min(), c[big, 0].max(), c[big, 1].min(), c[big, 1].max(),
    c[big, 2].min(), c[big, 2].max()))
thrU = np.percentile(ucm, 99.9)
bigU = ucm > thrU
print('big-U(>p99.9): n=%d x %.3f..%.3f y %.3f..%.3f z %.3f..%.3f' % (
    bigU.sum(), c[bigU, 0].min(), c[bigU, 0].max(), c[bigU, 1].min(),
    c[bigU, 1].max(), c[bigU, 2].min(), c[bigU, 2].max()))
solid = alpha > 0.5
print('solid cells: %d  mean|p_null| solid=%.4e fluid=%.4e' % (
    solid.sum(), np.abs(p[solid]).mean(), np.abs(p[~solid]).mean()))
