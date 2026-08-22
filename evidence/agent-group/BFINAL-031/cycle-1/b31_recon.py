#!/usr/bin/env python3
"""BFINAL-031: assembly-chain unified adjudication instrument (offline only).

Stages:
  python3 b31_recon.py m0m1  -> mesh+fields+gates, M4 three-caliber gate (H1/H2/H3)
  python3 b31_recon.py m2    -> itemized reconciliation matrix (both labels x D1/D2/D3)

Conventions pinned in PREREGISTRATION.md / NOTEBOOK.md (B13 labels):
  Uc = pressureDrop adjoint velocity, Ub = thermalCoupling adjoint velocity,
  pc = pressureDrop adjoint pressure, pb = thermalCoupling adjoint pressure.
  gsensh_PD = mom(-dAlphaDxh*(U&Uc)*V) + prow(-T_pc*dAlphaDxh), T_pc = J_P^T pc.
  fsenshMeanT = mom(-dAlphaDxh*(U&Ub)*V) + prow(-T_pb*dAlphaDxh) + flux(+Gx field) + C.
  Chain: mask -> drho -> eta corr -> Helmholtz (L-I)post=-mid, b=1 uniform -> mask.
"""
import sys, time, json, re
import numpy as np
import scipy.sparse as sp

t0 = time.time()
def log(m): print("[%.1fs] %s" % (time.time() - t0, m), flush=True)

MESH = "/home/ys/dsH/b16_mesh"
POLY = MESH + "/constant/polyMesh"
SB = "/home/ys/dsH/b25_qgate/stageB2"
Q1 = "/home/ys/dsH/b25_qgate/1"
Q8 = "/home/ys/dsH/b8_verify_diag/1"
E26 = "/home/ys/dsH/b8_verify_diag"
OUT = "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-031/cycle-1"
JT_H7 = "/home/ys/dsH/b24_diag/explicitJT_H7.mtx"

N = 33600; NV = 3 * N; NUNK = 4 * N
ALPHA_REL = 0.4
NU = 5.19009e-05
SMALL = 1e-15
DEL = 8.0                    # validationProjectionBeta (b25 optProperties)

ADJ_J = {"D1": -0.0358350400485765, "D2": -0.02492599038981423, "D3": -0.01200301226991504}
ADJ_GDP = {"D1": -5.23029384626718, "D2": 1.0505249811003, "D3": 13.60342258150915}
FD_J = {"D1": 0.04252469623637622, "D2": -0.003528972517971574, "D3": -0.0004449511895182612}
FD_GDP = {"D1": -2.536260986929495, "D2": 0.5363662735724528, "D3": 6.480676419641229}

# ---------------- field parsing (b26a/b27 validated read_field) ----------------
def read_field(path):
    txt = open(path).read()
    i = txt.index("internalField")
    seg_line = txt[i:txt.index("\n", i)]
    if re.match(r"internalField\s+uniform\b", seg_line):
        u = seg_line[seg_line.index("uniform") + len("uniform"):].strip().rstrip(";").strip()
        toks = []
        for tok in u.replace("(", " ").replace(")", " ").split():
            try: toks.append(float(tok))
            except ValueError: break
        val = np.array(toks) if len(toks) > 1 else float(toks[0])
        internal = None
    else:
        j = txt.index("\n(", i)
        k = txt.index("\n)", j)
        vals = []
        for line in txt[j + 1:k].splitlines():
            for tok in line.replace(";", " ").split():
                for tt in tok.replace("(", " ").replace(")", " ").split():
                    vals.append(float(tt))
        internal = np.array(vals)
        val = None
    return internal, val

def get_internal(field, Nexp):
    internal, val = field
    if internal is None:
        v = val if np.isscalar(val) else np.asarray(val)
        if np.isscalar(v): return np.full(Nexp, v)
        return np.broadcast_to(v, (Nexp,) + np.asarray(v).shape).copy()
    return internal

def rl(p):
    t = open(p).read(); i = t.index("\n("); j = t.index("\n)", i)
    return np.array([int(x) for x in t[i + 2:j].split()], dtype=np.int64)

# ---------------- mesh (exact pyramid centroids, b26a recipe) ----------------
def load_mesh():
    pts = []
    intxt = False
    for line in open(POLY + "/points"):
        s = line.strip()
        if not intxt:
            if s == "(": intxt = True
            continue
        if s == ")": break
        pts.append([float(t) for t in s.strip("()").split()])
    pts = np.array(pts)
    faces = []
    intxt = False
    for line in open(POLY + "/faces"):
        s = line.strip()
        if not intxt:
            if s == "(": intxt = True; continue
        if s == ")": break
        if "(" in s and s.endswith(")"):
            faces.append([int(v) for v in s[s.index("(") + 1:-1].split()])
    oa = rl(POLY + "/owner"); na = rl(POLY + "/neighbour")
    nIF = len(na); nF = len(faces)
    Sf = np.zeros((nF, 3))
    for fi in range(nF):
        fv = pts[faces[fi]]
        Sf[fi] = 0.5 * np.cross(fv, np.roll(fv, -1, axis=0)).sum(axis=0)
    magSf = np.linalg.norm(Sf, axis=1)
    fCtrs = np.array([pts[fv].mean(axis=0) for fv in faces])
    owner = oa[:nIF]; nei = na[:nIF]; bc = oa[nIF:]
    nB = nF - nIF
    bfile = open(POLY + "/boundary").read()
    patch_names = re.findall(r"^\s{4}(\w+)\s*$", bfile, re.M)
    patch_sizes = [int(m) for m in re.findall(r"nFaces\s+(\d+);", bfile)]
    cEst = np.zeros((N, 3)); nFc = np.zeros(N, dtype=np.int64)
    np.add.at(cEst, oa, fCtrs); np.add.at(nFc, oa, 1)
    np.add.at(cEst, na, fCtrs[:nIF]); np.add.at(nFc, na, 1)
    cEst /= nFc[:, None]
    num = np.zeros((N, 3)); vol3 = np.zeros(N)
    SfI = Sf[:nIF]; fCtrsI = fCtrs[:nIF]
    p3o = np.einsum("ij,ij->i", SfI, fCtrsI - cEst[owner])
    pc_o = 0.75 * fCtrsI + 0.25 * cEst[owner]
    np.add.at(num, owner, pc_o * p3o[:, None]); np.add.at(vol3, owner, p3o)
    p3n = np.einsum("ij,ij->i", -SfI, fCtrsI - cEst[nei])
    pc_n = 0.75 * fCtrsI + 0.25 * cEst[nei]
    np.add.at(num, nei, pc_n * p3n[:, None]); np.add.at(vol3, nei, p3n)
    SfB = Sf[nIF:]; fCtrsB = fCtrs[nIF:]
    p3b = np.einsum("ij,ij->i", SfB, fCtrsB - cEst[bc])
    pc_b = 0.75 * fCtrsB + 0.25 * cEst[bc]
    np.add.at(num, bc, pc_b * p3b[:, None]); np.add.at(vol3, bc, p3b)
    C = num / vol3[:, None]
    V_pyr = vol3 / 3.0
    V = np.loadtxt(MESH + "/b16mesh_V_nueff.mtx")[:, 0]
    vgate = float(np.linalg.norm(V_pyr - V) / np.linalg.norm(V))
    dvec = C[nei] - C[owner]
    dc = magSf[:nIF] / np.einsum("ij,ij->i", SfI, dvec)
    w = np.loadtxt(MESH + "/b16mesh_weights.mtx")
    bnd = np.loadtxt(MESH + "/b16mesh_boundary.mtx")
    bc_cell = bnd[:, 0].astype(np.int64); bSf = bnd[:, 1:4]
    bmag = bnd[:, 4]; bdel = bnd[:, 5]
    Ufix = bnd[:, 7] > 0.5
    return dict(pts=pts, owner=owner, nei=nei, bc=bc, nIF=nIF, nB=nB,
                Sf=Sf, magSf=magSf, fCtrs=fCtrs, patch_names=patch_names,
                patch_sizes=patch_sizes, C=C, V=V, vgate=vgate, dc=dc, w=w,
                bc_cell=bc_cell, bSf=bSf, bmag=bmag, bdel=bdel, Ufix=Ufix,
                SfI=Sf[:nIF], SfB=Sf[nIF:], fCtrsB=fCtrsB)

def dir_deriv(D6, ex, ey, ez):
    return (ex * ex * D6[:, 0] + ex * ey * D6[:, 1] + ex * ez * D6[:, 2]
            + ey * ex * D6[:, 1] + ey * ey * D6[:, 3] + ey * ez * D6[:, 4]
            + ez * ex * D6[:, 2] + ez * ey * D6[:, 4] + ez * ez * D6[:, 5])

def rel(a, b):
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), SMALL))

def load_directions(mesh, designMask, x):
    Lx = mesh["C"][:, 0].max() - mesh["C"][:, 0].min()
    Ly = mesh["C"][:, 1].max() - mesh["C"][:, 1].min()
    Lz = mesh["C"][:, 2].max() - mesh["C"][:, 2].min()
    xx = mesh["C"][:, 0] / max(Lx, SMALL); yy = mesh["C"][:, 1] / max(Ly, SMALL)
    zz = mesh["C"][:, 2] / max(Lz, SMALL)
    act = designMask > 0.5
    pi_ = np.pi
    D = {}
    D["D1"] = np.where(act, np.sin(2*pi_*xx)*np.cos(pi_*yy) + 0.5*np.sin(3*pi_*zz), 0.0)
    D["D2"] = np.where(act, np.sin(4*pi_*xx) + 0.3*np.cos(2*pi_*xx), 0.0)
    D["D3"] = np.where(act, np.cos(2*pi_*yy)*np.sin(2*pi_*zz), 0.0)
    oob = (x < 0.02) | (x > 0.98)
    for k in D: D[k] = np.where(oob, 0.0, D[k])
    for k in D:
        mx = np.abs(D[k]).max()
        if mx > 1e-30: D[k] /= mx
    return D

# ============================ field helpers ============================
def sc(name, n=N, qdir=None):
    return get_internal(read_field((qdir or Q1) + "/" + name), n)

def vec(name, qdir=None):
    internal, val = read_field((qdir or Q1) + "/" + name)
    if internal is None:
        v = np.asarray(val, float).reshape(3)
        return np.broadcast_to(v, (N, 3)).copy()
    return internal.reshape(N, 3)

# ============================ M0M1 ============================
def m0m1():
    log("=== B31 M0M1: gates + M4 three-caliber gate ===")
    mesh = load_mesh()
    log("mesh: nIF=%d nB=%d patches=%s" % (mesh["nIF"], mesh["nB"], mesh["patch_names"]))
    log("V gate (pyramid) relL2=%.3e (need <1e-12)" % mesh["vgate"])
    assert mesh["vgate"] < 1e-12, "BLOCKED: V gate"

    U = vec("U")
    dAlphaDxh = sc("dAlphaDxh")
    designMask = sc("designMask")
    xp = sc("xp")
    alpha = sc("alpha")
    fsenshMeanT = sc("fsenshMeanT")
    fsensMeanT = sc("fsensMeanT")
    gsenshPD = sc("gsenshPressureDrop")
    gsensPD = sc("gsensPressureDrop")
    dfdx = sc("dfdx")
    Uc25 = vec("Uc")
    Ub25 = vec("Ub")
    pc25 = sc("pc"); pb25 = sc("pb")
    T = sc("T"); Tb = sc("Tb")
    dDTDxh = sc("dDTDxh").reshape(N, 6)
    phiTh_int = get_internal(read_field(Q1 + "/phiThermal"), mesh["nIF"])
    coldmask = get_internal(read_field(Q1 + "/coldFaceMask"), mesh["nIF"])
    phi_int = get_internal(read_field(Q1 + "/phi"), mesh["nIF"])
    maskF = get_internal(read_field(Q1 + "/designFilterFaceMask"), mesh["nIF"])
    log("b25 fields loaded")

    # b8 adjoints (export-era lambda; Ub/pb may be uniform zero - thermal label unsolved)
    Uc8 = vec("Uc", Q8)
    Ub8 = vec("Ub", Q8)
    _, pb8u = read_field(Q8 + "/pb")
    pb8 = np.full(N, float(pb8u) if np.isscalar(pb8u) else 0.0)
    pc8 = sc("pc", qdir=Q8)
    log("b8 adjoints loaded: |Uc8|L2=%.4e |pc8|L2=%.4e (pb8 uniform=%r)"
        % (np.linalg.norm(Uc8), np.linalg.norm(pc8), np.isscalar(pb8u)))

    rxprM = np.loadtxt(E26 + "/rxpr_prod_gsensh_momentum.mtx")
    rxprP = np.loadtxt(E26 + "/rxpr_prod_gsensh_pressurerow.mtx")
    rxprT = np.loadtxt(E26 + "/rxpr_prod_gsensh_total.mtx")
    rxprG = np.loadtxt(E26 + "/rxpr_prod_gsens.mtx")
    z = np.loadtxt(E26 + "/stageB6_rxc_z_analytic.mtx").reshape(3, N)
    D = load_directions(mesh, designMask, xp)
    log("directions rebuilt; |z| rows: %s" % [float(np.linalg.norm(zk)) for zk in z])

    # ---- segment recomputation, both lambdas ----
    mom_b8 = -dAlphaDxh * np.einsum("ci,ci->c", U, Uc8) * mesh["V"]
    mom_b25 = -dAlphaDxh * np.einsum("ci,ci->c", U, Uc25) * mesh["V"]
    rows = []
    def cmp(name, a, b_):
        r = rel(a, b_)
        mx = float(np.abs(a - b_).max())
        sf = int((((a > 0) != (b_ > 0)) & (np.abs(b_) > 1e-300)).sum())
        rows.append((name, r, mx, sf))
        log("CMP %-42s relL2=%.6e maxabs=%.4e signflips=%d" % (name, r, mx, sf))
        return r
    cmp("rxpr_momentum vs mom_recomp_b8", rxprM, mom_b8)
    cmp("rxpr_momentum vs mom_recomp_b25", rxprM, mom_b25)
    cmp("rxpr_total vs gsenshPressureDrop_disk_b25", rxprT, gsenshPD)
    cmp("rxpr_total vs mom+prow?? (mom_b8 only)", rxprT, mom_b8)
    cmp("rxpr_gsens(post,b8run) vs gsensPressureDrop_disk_b25", rxprG, gsensPD)

    # ---- projections: three calibers on the table ----
    res = {}
    for k, zk, Dk in zip(("D1", "D2", "D3"), z, (D["D1"], D["D2"], D["D3"])):
        res[k] = dict(
            dot_rxprT_z=float(np.dot(rxprT, zk)),
            dot_rxprM_z=float(np.dot(rxprM, zk)),
            dot_rxprP_z=float(np.dot(rxprP, zk)),
            dot_gsenshDisk_z=float(np.dot(gsenshPD, zk)),
            dot_gsensPDdisk_D=float(np.dot(gsensPD, Dk)),
            dot_dfdx_D=float(np.dot(dfdx, Dk)),
            dot_fsensMDisk_D=float(np.dot(fsensMeanT, Dk)),
            dot_z_D=float(np.dot(zk, Dk)),
            ADJ_gDP=ADJ_GDP[k], ADJ_J=ADJ_J[k])
        log("[%s] rxprT.z=%+.6e gsenshDisk.z=%+.6e | gsensPDisk.D=%+.6e(tsv %+.6e) "
            "dfdx.D=%+.6e(tsv %+.6e) | z.D=%.4f"
            % (k, res[k]["dot_rxprT_z"], res[k]["dot_gsenshDisk_z"],
               res[k]["dot_gsensPDdisk_D"], ADJ_GDP[k],
               res[k]["dot_dfdx_D"], ADJ_J[k], res[k]["dot_z_D"]))

    verdict = None
    r8 = rel(rxprM, mom_b8); r25 = rel(rxprM, mom_b25)
    if r8 < 1e-8 and r25 > 1e-3:
        verdict = "H1: rxpr exports are self-consistent with the b8-run lambda, NOT state-B"
    elif r25 < 1e-8:
        verdict = "H2/H3: rxpr exports ARE state-B; discrepancy sits in chain/projection layer"
    else:
        verdict = "INCONCLUSIVE per-cell; see ratios r8=%.3e r25=%.3e" % (r8, r25)
    log("M4 GATE VERDICT: %s" % verdict)

    with open(OUT + "/b31_m0m1_compare.tsv", "w") as f:
        f.write("cmp\trelL2\tmaxabs\tsignflips\n")
        for nm, r, mx, sf in rows:
            f.write("%s\t%.6e\t%.6e\t%d\n" % (nm, r, mx, sf))
    with open(OUT + "/b31_m4_projections.tsv", "w") as f:
        f.write("dir\tdot_rxprT_z\tdot_rxprM_z\tdot_rxprP_z\tdot_gsenshDisk_z\t"
                "dot_gsensPDdisk_D\tADJ_gDP_tsv\tdot_dfdx_D\tADJ_J_tsv\t"
                "dot_fsensMDisk_D\tdot_z_D\n")
        for k in ("D1", "D2", "D3"):
            r = res[k]
            f.write("%s\t%+.10e\t%+.10e\t%+.10e\t%+.10e\t%+.10e\t%+.10e\t%+.10e\t%+.10e\t%+.10e\t%.6f\n"
                    % (k, r["dot_rxprT_z"], r["dot_rxprM_z"], r["dot_rxprP_z"],
                       r["dot_gsenshDisk_z"], r["dot_gsensPDdisk_D"], r["ADJ_gDP"],
                       r["dot_dfdx_D"], r["ADJ_J"], r["dot_fsensMDisk_D"], r["dot_z_D"]))
    json.dump(dict(verdict=verdict, r8=r8, r25=r25, proj=res),
              open(OUT + "/b31_m0m1.json", "w"), indent=1)

    # persist for m2
    np.savez_compressed(OUT + "/b31_stage_cache.npz",
                        dAlphaDxh=dAlphaDxh, designMask=designMask, xp=xp,
                        U=U, Uc=Uc25, Ub=Ub25, pc=pc25, pb=pb25,
                        fsenshMeanT=fsenshMeanT, fsensMeanT=fsensMeanT,
                        gsenshPD=gsenshPD, gsensPD=gsensPD, dfdx=dfdx,
                        T=T, Tb=Tb, dDTDxh=dDTDxh.reshape(N, 6),
                        phiTh=phiTh_int, coldmask=coldmask, phi_int=phi_int,
                        maskF=maskF, z=z, mom_b8=mom_b8, mom_b25=mom_b25,
                        rxprM=rxprM, rxprP=rxprP, rxprT=rxprT, rxprG=rxprG,
                        D1=D["D1"], D2=D["D2"], D3=D["D3"])
    log("saved stage cache; M0M1 done t=%.1fs" % (time.time() - t0))

# ============================ M2 ============================
def rho_project(xpv, eta):
    """Heaviside projected density (diff.c / filter_x.H formulas), del=DEL."""
    pe = np.where(
        xpv <= eta,
        eta * (np.exp(-DEL * (1.0 - xpv / eta)) - (1.0 - xpv / eta) * np.exp(-DEL)),
        eta + (1.0 - eta) * (1.0 - np.exp(-DEL * (xpv - eta) / (1.0 - eta))
                             + (xpv - eta) * np.exp(-DEL) / (1.0 - eta)))
    return pe

def drho_deta(xpv, eta):
    pe1 = np.exp(-DEL * (1.0 - xpv / eta))
    pe2 = np.exp(-DEL * (xpv - eta) / (1.0 - eta))
    dr = np.where(xpv <= eta, DEL * pe1 + np.exp(-DEL), DEL * pe2 + np.exp(-DEL))
    de = np.where(xpv <= eta,
                  pe1 * (1.0 - DEL * xpv / eta) - np.exp(-DEL),
                  pe2 * (1.0 - DEL * (1.0 - xpv) / (1.0 - eta)) - np.exp(-DEL))
    return dr, de

def m2():
    log("=== B31 M2: chain replication + itemized reconciliation matrix ===")
    mesh = load_mesh()
    log("V gate relL2=%.3e" % mesh["vgate"])
    ca = np.load(OUT + "/b31_stage_cache.npz")
    dAlphaDxh = ca["dAlphaDxh"]; designMask = ca["designMask"]; xp = ca["xp"]
    U = ca["U"]; Uc25 = ca["Uc"]; Ub25 = ca["Ub"]; pc25 = ca["pc"]; pb25 = ca["pb"]
    fsensh = ca["fsenshMeanT"]; fsensM = ca["fsensMeanT"]
    gsenshPD = ca["gsenshPD"]; gsensPD = ca["gsensPD"]; dfdx = ca["dfdx"]
    T = ca["T"]; Tb = ca["Tb"]; dDTDxh = ca["dDTDxh"]
    phiTh = ca["phiTh"]; coldmask = ca["coldmask"]; phi_int = ca["phi_int"]
    maskF = ca["maskF"]; z = ca["z"]
    mom_b25 = ca["mom_b25"]
    D = {"D1": ca["D1"], "D2": ca["D2"], "D3": ca["D3"]}
    xh_disk = sc("xh")

    # ---- gates G1/G2: eta5 bisection + projection replication vs disk xh ----
    V = mesh["V"]; act = designMask > 0.5
    def diffvol(eta):
        return float(np.sum(np.where(act, (xp - rho_project(xp, eta)) * V, 0.0)))
    eta0, eta1 = 1e-4, 0.9999
    y0, y1 = diffvol(eta0), diffvol(eta1)
    if y0 * y1 > 0:
        eta5 = 0.5
        log("G1 WARNING: volume eq not bracketed y0=%.3e y1=%.3e -> eta=0.5" % (y0, y1))
    else:
        for _ in range(200):
            eta5 = 0.5 * (eta0 + eta1)
            y5 = diffvol(eta5)
            if y0 * y5 < 0.0: eta1 = eta5; y1 = y5
            else: eta0 = eta5; y0 = y5
            if (eta1 - eta0) <= 1e-14: break
    log("G1 eta5 bisected: eta5=%.16f  diffvol=%.3e" % (eta5, diffvol(eta5)))
    drho, dProjEta = drho_deta(xp, eta5)
    drho = np.where(act, drho, 0.0); dProjEta = np.where(act, dProjEta, 0.0)
    xh_recomp = np.where(act, rho_project(xp, eta5), 0.0)
    xh_disk_m = np.where(act, xh_disk, 0.0)
    log("G2 xh projection replication vs disk: relL2=%.3e maxabs=%.3e"
        % (rel(xh_recomp, xh_disk_m), float(np.abs(xh_recomp - xh_disk_m).max())))

    # ---- Helmholtz operator assembly + gate G0: forward filter on x ----
    owner, nei = mesh["owner"], mesh["nei"]
    dc, w, SfI, magSfI = mesh["dc"], mesh["w"], mesh["SfI"], mesh["magSf"][:mesh["nIF"]]
    kfv = maskF * magSfI * dc
    diag = np.zeros(N)
    np.add.at(diag, owner, -kfv); np.add.at(diag, nei, -kfv)
    rows_ = np.concatenate([owner, nei, np.arange(N)])
    cols_ = np.concatenate([nei, owner, np.arange(N)])
    vals_ = np.concatenate([kfv, kfv, diag - 1.0])       # (L - I)
    A = sp.csr_matrix((vals_, (rows_, cols_)), shape=(N, N))
    from scipy.sparse.linalg import splu
    lu = splu(A.tocsc())
    x_full = sc("x")                       # raw design field, ALL cells (filter_x.H RHS)
    x_fwd = lu.solve(-x_full)
    log("G0 forward filter replication: (L-I)xp=-x vs disk xp: relL2=%.3e maxabs=%.3e"
        % (rel(x_fwd, xp), float(np.abs(x_fwd - xp).max())))

    denom = float(np.sum(designMask * dProjEta * V))
    def chain_apply(pre):
        p = pre * designMask
        corr = float(np.sum(designMask * dProjEta * p)) / denom
        mid = drho * p + designMask * V * (1.0 - drho) * corr
        return lu.solve(-mid) * designMask
    log("chain replication ready (denom=%.6e)" % denom)

    # ---- H2/H3 adjudication: chain disk pre fields, project on Dk ----
    post_gdp_recomp = chain_apply(gsenshPD)
    post_J_recomp = chain_apply(fsensh)
    rows_c = []
    for k in ("D1", "D2", "D3"):
        Dk = D[k]
        r_g = rel(post_gdp_recomp, gsensPD)
        r_j = rel(post_J_recomp, fsensM)
        rows_c.append(dict(dir=k,
                           gdp_recomp=float(np.dot(post_gdp_recomp, Dk)),
                           gdp_disk=float(np.dot(gsensPD, Dk)),
                           gdp_tsv=ADJ_GDP[k],
                           J_recomp=float(np.dot(post_J_recomp, Dk)),
                           J_disk=float(np.dot(fsensM, Dk)),
                           J_tsv=ADJ_J[k],
                           fieldrel_gdp=r_g, fieldrel_J=r_j))
        log("[CHAIN] %s: gdp recomp=%+.6e disk=%+.6e tsv=%+.6e (fieldrel=%.2e) | "
            "J recomp=%+.6e disk=%+.6e tsv=%+.6e (fieldrel=%.2e)"
            % (k, rows_c[-1]["gdp_recomp"], rows_c[-1]["gdp_disk"], ADJ_GDP[k], r_g,
               rows_c[-1]["J_recomp"], rows_c[-1]["J_disk"], ADJ_J[k], r_j))
    chain_verdict = ("H3: production chain machinery reproduces post fields; z_file stale"
                     if max(r_g, r_j) < 1e-6 else
                     "H2: production chain itself deviates from replication")
    log("CHAIN VERDICT: %s" % chain_verdict)

    # ---- pressure-row segments via explicitJT_H7 (M = J^T stored) ----
    # T = J_P^T lambda_P = M[:, NV:] @ lambda_P  (pressure rows of J = cols NV: of M)
    log("loading %s (523MB MM coordinate) ..." % JT_H7)
    with open(JT_H7) as f:
        f.readline(); f.readline()
        d3 = np.loadtxt(f)
    Mrows = (d3[:, 0] - 1).astype(np.int64); Mcols = (d3[:, 1] - 1).astype(np.int64)
    Mvals = d3[:, 2]; del d3
    JT = sp.coo_matrix((Mvals, (Mrows, Mcols)), shape=(NUNK, NUNK)).tocsr()
    del Mrows, Mcols, Mvals
    log("explicitJT_H7 loaded: nnz=%d (duplicates summed by tocsr)" % JT.nnz)
    vpd = np.zeros(NUNK); vpd[NV:] = pc25
    T_pd = JT @ vpd                       # J_P^T pc  (pressureDrop label)
    vtb = np.zeros(NUNK); vtb[NV:] = pb25
    T_tb = JT @ vtb                       # J_P^T pb  (thermalCoupling label)
    prow_gdp = -T_pd * dAlphaDxh
    prow_J = -T_tb * dAlphaDxh
    mom_J = -dAlphaDxh * np.einsum("ci,ci->c", U, Ub25) * V

    # ---- C field (b27 M3 recipe, validated vs B19 anchors at 1e-12) ----
    Cc = mesh["C"]
    dvecC = Cc[nei] - Cc[owner]
    efc = dvecC / np.maximum(np.linalg.norm(dvecC, axis=1), 1e-30)[:, None]
    D6o = dDTDxh[owner]; D6n = dDTDxh[nei]
    dd_o = dir_deriv(D6o, efc[:, 0], efc[:, 1], efc[:, 2])
    dd_n = dir_deriv(D6n, efc[:, 0], efc[:, 1], efc[:, 2])
    base = (Tb[owner] - Tb[nei]) * (T[owner] - T[nei]) * dc * magSfI
    Cfield = np.zeros(N)
    np.add.at(Cfield, owner, -w * dd_o * base)
    np.add.at(Cfield, nei, -(1.0 - w) * dd_n * base)

    # ---- g face functional + outlet dJ/dphi (b26a recipe) ----
    nIF = mesh["nIF"]; bc_cell = mesh["bc_cell"]; bSf = mesh["bSf"]
    bmag = mesh["bmag"]; bdel = mesh["bdel"]; Ufix = mesh["Ufix"]
    patch_names = mesh["patch_names"]; patch_sizes = mesh["patch_sizes"]
    bface_start = {}; o = nIF
    for pn, sz in zip(patch_names, patch_sizes):
        bface_start[pn] = o; o += sz
    jumpT = T[nei] - T[owner]
    adjDown = np.where(phiTh >= 0.0, Tb[nei], Tb[owner])
    g_int = -coldmask * adjDown * jumpT
    Mflow, Tref = np.loadtxt(SB + "/wstate_outletMeta.mtx")
    o0 = bface_start["outlet"] - nIF
    no = patch_sizes[patch_names.index("outlet")]
    # outlet boundary flux values from phi file boundary block
    txt_phi = open(Q1 + "/phi").read()
    mb = re.search(r"outlet\s*\{.*?nonuniform.*?\n\(", txt_phi, re.S)
    ophi = None
    if mb is not None:
        j = txt_phi.index("\n(", mb.start()); k = txt_phi.index("\n)", j)
        vals = [float(x) for x in txt_phi[j+2:k].split()]
        ophi = np.array(vals)
    else:
        mu = re.search(r"outlet\s*\{[^}]*uniform\s+([^;]+);", txt_phi, re.S)
        ophi = np.full(no, float(mu.group(1).strip("()")))
    assert ophi is not None and len(ophi) == no, "outlet phi parse failed"
    Toc = T[bc_cell[o0:o0+no]]
    Tmix = float(np.sum(ophi * Toc) / Mflow)
    dJdphi_out = -(Toc - Tmix) / (Mflow * Tref)
    g_b = np.zeros(mesh["nB"])
    g_b[o0:o0+no] += dJdphi_out
    log("g functional: |int|max=%.4e Tmix=%.4f |dJdphi_out|max=%.4e"
        % (np.abs(g_int).max(), Tmix, np.abs(dJdphi_out).max()))

    # ---- mobility basis + dHbyA/drAU + offline Gx field (b26a recipe) ----
    def build_basis(phi_, phiB_, alpha_, nuEff_):
        qp = (phi_ >= 0).astype(float); qn = 1 - qp
        divp = np.zeros(N)
        np.add.at(divp, owner, phi_); np.add.at(divp, nei, -phi_)
        np.add.at(divp, bc_cell, phiB_)
        gam = w * nuEff_[owner] + (1 - w) * nuEff_[nei]
        kfv_ = gam * magSfI * dc
        lower = -qp * phi_ - kfv_; upper = qn * phi_ - kfv_
        diag_ = np.zeros(N)
        np.add.at(diag_, owner, -lower); np.add.at(diag_, nei, -upper)
        diag_ -= divp; diag_ += alpha_ * V
        sumOff = np.zeros(N)
        np.add.at(sumOff, owner, np.abs(upper)); np.add.at(sumOff, nei, np.abs(lower))
        ab = nuEff_[bc_cell] * bmag * bdel
        ic = np.where(Ufix, -ab, phiB_)
        Db = diag_.copy(); np.add.at(Db, bc_cell, np.abs(ic))
        Drel = np.maximum(np.abs(Db), sumOff) / ALPHA_REL
        return dict(lower=lower, upper=upper, diag=diag_, D_rel=Drel,
                    mob=V / Drel)
    U_B = np.loadtxt(SB + "/wstate_baseline_U.mtx")
    p_B = np.loadtxt(SB + "/wstate_baseline_p.mtx")
    phi_B = np.loadtxt(SB + "/wstate_baseline_phi.mtx")
    phiB_B = np.loadtxt(SB + "/wstate_baseline_phiB.mtx")
    nut = np.loadtxt(SB + "/wstate_baseline_nutFrozen.mtx")
    al_B = np.loadtxt(SB + "/wstate_baseline_alpha.mtx")
    nuE = NU + nut
    bs = build_basis(phi_B, phiB_B, al_B, nuE)
    mob_ex = np.loadtxt(SB + "/wstate_baseline_primalPressureMobility.mtx")
    mob_gate = float(np.linalg.norm(bs["mob"] - mob_ex) / np.linalg.norm(mob_ex))
    log("MOBILITY GATE relL2=%.3e (need <1e-12)" % mob_gate)
    assert mob_gate < 1e-12, "BLOCKED: mobility gate"
    mob = bs["mob"]

    def grad_p(p_):
        pf = w * p_[owner] + (1 - w) * p_[nei]
        g0 = np.zeros((N, 3))
        np.add.at(g0, owner, SfI * pf[:, None])
        np.add.at(g0, nei, -(SfI * pf[:, None]))
        pb_ = p_[bc_cell].copy()
        out_start = bface_start["outlet"] - nIF
        pb_[out_start:out_start + no] = 0.0
        np.add.at(g0, bc_cell, bSf * pb_[:, None])
        return g0 / V[:, None]

    src = (bs["D_rel"] - bs["diag"])[:, None] * U_B
    offU = np.zeros((N, 3))
    np.add.at(offU, owner, bs["upper"][:, None] * U_B[nei])
    np.add.at(offU, nei, bs["lower"][:, None] * U_B[owner])
    HbyA_base = mob[:, None] * (src - offU) / V[:, None] - mob[:, None] * grad_p(p_B)
    drAU = -mob * mob / ALPHA_REL
    dHbyA = (mob / ALPHA_REL)[:, None] * ((1.0 - ALPHA_REL) * U_B - HbyA_base)
    g0_int = magSfI * dc * (p_B[nei] - p_B[owner])
    am = np.where(~Ufix)[0]

    kf_int = (w * mob[owner] + (1 - w) * mob[nei]) * magSfI * dc

    def Gx_of(zdir):
        wc = dAlphaDxh * zdir
        dHw = dHbyA * wc[:, None]
        phx_int = np.einsum("fi,fi->f", SfI,
                            w[:, None] * dHw[owner] + (1 - w)[:, None] * dHw[nei])
        gf_int = w * drAU[owner] * wc[owner] + (1 - w) * drAU[nei] * wc[nei]
        phx_int -= gf_int * g0_int
        val = float(np.dot(g_int, phx_int))
        phx_b = np.einsum("fi,fi->f", bSf[am], dHw[bc_cell[am]])
        val += float(np.dot(g_b[am], phx_b))
        return val

    # offline Gx FIELD (linear in z): collect per-cell coefficients of wc=dAlphaDxh*z
    B1v = np.zeros((N, 3))
    np.add.at(B1v, owner, (w * g_int)[:, None] * SfI)
    np.add.at(B1v, nei, ((1.0 - w) * g_int)[:, None] * SfI)
    np.add.at(B1v, bc_cell[am], g_b[am][:, None] * bSf[am])
    B2v = np.zeros(N)
    np.add.at(B2v, owner, w * g0_int)
    np.add.at(B2v, nei, (1.0 - w) * g0_int)
    Fxfield = dAlphaDxh * (np.einsum("ci,ci->c", dHbyA, B1v) - drAU * B2v)
    # self-gate: field contraction must reproduce the direct formula AND b26a anchors
    GX_ANCHOR = {"D1": -1.5704339057e-03, "D2": 4.0431596817e-03, "D3": 8.4067076831e-03}
    fx_ok = True
    for k_i, k in (("D1", "D1"), ("D2", "D2"), ("D3", "D3")):
        ki = ("D1", "D2", "D3").index(k)
        vf = float(np.dot(Fxfield, z[ki]))
        vd = Gx_of(z[ki])
        anch = GX_ANCHOR[k]
        log("[GXGATE] %s field=%+.8e direct=%+.8e b26a_anchor=%+.8e rel(f/d)=%.2e ratio(f/a)=%.6f"
            % (k, vf, vd, anch, abs(vf - vd) / max(abs(vd), SMALL), vf / anch))
        if abs(vf - vd) > 1e-9 * max(abs(vd), 1e-30): fx_ok = False
    assert fx_ok, "BLOCKED: offline Gx field fails self-gate"

    # ---- source contractions ----
    wtrue = np.load("/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/"
                    "BFINAL-026/cycle-1/b26_wstar_h7.npz")
    def load_b(tag):
        m = np.loadtxt("/home/ys/dsH/b25_qgate/b18rhs_%s.mtx" % tag)
        bexp = np.zeros(NUNK)
        bexp[0:NV:3] = m[:, 0]; bexp[1:NV:3] = m[:, 1]
        bexp[2:NV:3] = m[:, 2]; bexp[NV:] = m[:, 3]
        return bexp
    b_TC = load_b("thermalCoupling")
    b_PD = load_b("pressureDrop")
    S_TC = {k: float(b_TC @ wtrue[k]) for k in ("D1", "D2", "D3")}
    S_PD = {k: float(b_PD @ wtrue[k]) for k in ("D1", "D2", "D3")}
    log("source contractions: b_TC.w = %s" % {k: "%.8e" % v for k, v in S_TC.items()})
    log("                      b_PD.w = %s" % {k: "%.8e" % v for k, v in S_PD.items()})

    # ---- implied production segments vs offline recomputation ----
    flux_implied = fsensh - Cfield - mom_J - prow_J      # production fluxDirect implied
    prow_gdp_implied = gsenshPD - mom_b25                # production prow(gDP) implied
    segcmp = [
        ("J_mom_off_vs_implied(fsensh-C-Fx-prow)", mom_J,
         fsensh - Cfield - Fxfield - prow_J),
        ("J_prow_off_vs_implied", prow_J, fsensh - Cfield - Fxfield - mom_J),
        ("J_C_off_vs_implied", Cfield, fsensh - Fxfield - mom_J - prow_J),
        ("J_Fx_off_vs_implied", Fxfield, fsensh - Cfield - mom_J - prow_J),
        ("GDP_mom_recomp_vs_(gsenshPD-prow_recomp)", mom_b25, gsenshPD - prow_gdp),
        ("GDP_prow_recomp_vs_implied(gsenshPD-mom)", prow_gdp, prow_gdp_implied),
    ]
    with open(OUT + "/b31_segcomp.tsv", "w") as f:
        f.write("cmp\trelL2\tmaxabs\tsignflips\n")
        for nm, a, b_ in segcmp:
            sf = int((((a > 0) != (b_ > 0)) & (np.abs(b_) > 1e-300)).sum())
            f.write("%s\t%.6e\t%.6e\t%d\n" % (nm, rel(a, b_), float(np.abs(a - b_).max()), sf))
            log("[SEG] %-46s relL2=%.6e maxabs=%.4e flips=%d"
                % (nm, rel(a, b_), float(np.abs(a - b_).max()), sf))

    # ---- identity B gates (flow-mediated source identity) ----
    for k_i, k in enumerate(("D1", "D2", "D3")):
        zk = z[k_i]
        lhsJ = float(np.dot(mom_J + prow_J, zk))
        lhsP = float(np.dot(mom_b25 + prow_gdp, zk))
        log("[IDB] %s: (mom+prow)_J.z = %+.8e vs b_TC.w = %+.8e (rel %.3e) | "
            "(mom+prow)_PD.z = %+.8e vs b_PD.w = %+.8e (rel %.3e)"
            % (k, lhsJ, S_TC[k], abs(lhsJ - S_TC[k]) / max(abs(S_TC[k]), SMALL),
               lhsP, S_PD[k], abs(lhsP - S_PD[k]) / max(abs(S_PD[k]), SMALL)))

    # ---- reconciliation matrix (chained per-segment contributions) ----
    segs_J = dict(mom=mom_J, prow=prow_J, flux=Fxfield, C=Cfield)
    segs_P = dict(mom=mom_b25, prow=prow_gdp)
    chain_seg_J = {n: chain_apply(f) for n, f in segs_J.items()}
    chain_seg_P = {n: chain_apply(f) for n, f in segs_P.items()}
    out_rows = []
    for k_i, k in enumerate(("D1", "D2", "D3")):
        Dk = D[k]; zk = z[k_i]
        cj = {n: float(np.dot(f, Dk)) for n, f in chain_seg_J.items()}
        cp = {n: float(np.dot(f, Dk)) for n, f in chain_seg_P.items()}
        sumJ = sum(cj.values()); sumP = sum(cp.values())
        preJ_z = {n: float(np.dot(f, zk)) for n, f in segs_J.items()}
        preP_z = {n: float(np.dot(f, zk)) for n, f in segs_P.items()}
        out_rows.append(dict(label="J", dir=k, projP=ADJ_J[k],
                             post_disk=float(np.dot(fsensM, Dk)),
                             **{"chain_" + n: v for n, v in cj.items()},
                             sum_chain=sumJ, pre_disk_z=float(np.dot(fsensh, zk)),
                             **{"prez_" + n: v for n, v in preJ_z.items()},
                             src=S_TC[k], FD=FD_J[k], gap_proj_FD=ADJ_J[k] - FD_J[k],
                             gap_proj_src=ADJ_J[k] - S_TC[k]))
        out_rows.append(dict(label="gDP", dir=k, projP=ADJ_GDP[k],
                             post_disk=float(np.dot(gsensPD, Dk)),
                             **{"chain_" + n: v for n, v in cp.items()},
                             sum_chain=sumP, pre_disk_z=float(np.dot(gsenshPD, zk)),
                             **{"prez_" + n: v for n, v in preP_z.items()},
                             src=S_PD[k], FD=FD_GDP[k], gap_proj_FD=ADJ_GDP[k] - FD_GDP[k],
                             gap_proj_src=ADJ_GDP[k] - S_PD[k]))
    keys = ["label", "dir", "projP", "post_disk",
            "chain_mom", "chain_prow", "chain_flux", "chain_C", "sum_chain",
            "pre_disk_z", "prez_mom", "prez_prow", "prez_flux", "prez_C",
            "src", "FD", "gap_proj_FD", "gap_proj_src"]
    with open(OUT + "/b31_matrix.tsv", "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in out_rows:
            f.write("\t".join((("%.10e" % r[kk]) if isinstance(r[kk], float) else str(r[kk]))
                              for kk in keys) + "\n")
    json.dump(dict(chain_verdict=chain_verdict, eta5=eta5,
                   xh_gate=float(rel(xh_recomp, xh_disk_m)),
                   fwd_gate=float(rel(x_fwd, xp)),
                   mob_gate=mob_gate, S_TC=S_TC, S_PD=S_PD,
                   rows=out_rows), open(OUT + "/b31_results.json", "w"), indent=1)
    log("matrix written; M2 done t=%.1fs" % (time.time() - t0))



if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "m0m1"
    if what == "m0m1":
        m0m1()
    elif what == "m2":
        m2()
    else:
        log("stage '%s' not implemented yet" % what)
