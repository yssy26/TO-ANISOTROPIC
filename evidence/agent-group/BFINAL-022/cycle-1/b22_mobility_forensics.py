#!/usr/bin/env python3
"""BFINAL-022 W4.3 addendum: forensics of the offline rAU-rebuild error.

The exact exported primalPressureMobility(B) (W1 capture) differs from the
offline D_rel rebuild (the b16-lineage instrument used by B15-B21) by
12.9% relL2.  This script pins the STRUCTURE of that deviation and tests
three OF7-semantics suspects:

  S1  boundary diagonal: relax() removes cmptMin(ic) after the clamp but
      D() re-adds addCmptAvBoundaryDiag -> for uniform-component ic these
      cancel and the OLD formula (add |ic| pre-clamp, no removal) is right.
      (fvMatrix.C relax tail + D().)
  S2  "bounded Gauss upwind": fvmDiv = upwind - fvm::Sp(surfaceIntegrate
      (phi)) -> the bounded correction enters the diag as -V*div(phi)
      (fvmSup.C: diag += mesh.V()*sp), NOT -div(phi).  Numerically
      negligible here (V=1.2e-10) -- both variants agree.
  S3  laplacian face gamma: true lower=upper=face-interpolated
      nuEff*magSf*deltaCoeffs (gaussLaplacianScheme.C); the b21 script
      uses one-sided owner/neighbour nuEff.  Effect: ~0.4% only.

Structure of the residual deviation (the part none of S1-S3 explains):
  solid cells exact, deviation concentrated in nut-active and
  clamp-bound (sumOff>|D|) transition/fluid cells.
"""
import numpy as np

MESH="/home/ys/dsH/b16_mesh"; POLY=MESH+"/constant/polyMesh"
SB="/home/ys/dsH/b22_gauge/stageB2"
N=33600; ALPHA_REL=0.4; NU=5.19009e-05

def read_labels(path):
    txt=open(path).read(); i=txt.index("\n("); j=txt.index("\n)",i)
    return np.array([int(t) for t in txt[i+2:j].split()],dtype=np.int64)
owner_all=read_labels(POLY+"/owner"); nei_all=read_labels(POLY+"/neighbour")
nIF=len(nei_all); owner=owner_all[:nIF]; nei=nei_all[:nIF]; nF=len(owner_all)
pts=[]; intxt=False
for line in open(POLY+"/points"):
    s=line.strip()
    if not intxt:
        if s=="(": intxt=True
        continue
    if s==")": break
    pts.append([float(t) for t in s.strip("()").split()])
pts=np.array(pts)
faces=[]; intxt=False
for line in open(POLY+"/faces"):
    s=line.strip()
    if not intxt:
        if s=="(": intxt=True; continue
    if s==")": break
    if "(" in s and s.endswith(")"):
        faces.append([int(v) for v in s[s.index("(")+1:-1].split()])
Sf_all=np.zeros((nF,3))
for fi in range(nF):
    fv=pts[faces[fi]]
    Sf_all[fi]=0.5*np.cross(fv,np.roll(fv,-1,axis=0)).sum(axis=0)
magSf_all=np.linalg.norm(Sf_all,axis=1)
cellverts={}
for fi in range(nF):
    cells=(owner_all[fi],nei_all[fi]) if fi<nIF else (owner_all[fi],)
    for c in cells: cellverts.setdefault(c,set()).update(faces[fi])
Ccen=np.zeros((N,3))
for c,vs in cellverts.items(): Ccen[c]=pts[sorted(vs)].mean(axis=0)
d_vec=Ccen[nei]-Ccen[owner]
deltaCoeffs=magSf_all[:nIF]/np.einsum("ij,ij->i",Sf_all[:nIF],d_vec)
V=np.loadtxt(MESH+"/b16mesh_V_nueff.mtx")[:,0]
w=np.loadtxt(MESH+"/b16mesh_weights.mtx")
bnd=np.loadtxt(MESH+"/b16mesh_boundary.mtx")
bc_cell=bnd[:,0].astype(np.int64); bmag=bnd[:,4]; bdel=bnd[:,5]
Ufixed_b=bnd[:,7]>0.5
bndcell=np.zeros(N,bool); bndcell[bc_cell]=True

mob=np.loadtxt(SB+"/wstate_baseline_primalPressureMobility.mtx")
alpha=np.loadtxt(SB+"/wstate_baseline_alpha.mtx")
phi_B=np.loadtxt(SB+"/wstate_baseline_phi.mtx"); phiB_B=np.loadtxt(SB+"/wstate_baseline_phiB.mtx")
nut=np.loadtxt(SB+"/wstate_baseline_nutFrozen.mtx")
nuEff=NU+nut
upw=phi_B>=0.0; qp=upw.astype(float); qn=1.0-qp

def rebuild(interp_gamma=False, divcoef=1.0):
    if interp_gamma:
        gam=(w*nuEff[owner]+(1-w)*nuEff[nei])*magSf_all[:nIF]*deltaCoeffs
        kfv_o=kfv_n=gam
    else:
        kfv_o=nuEff[owner]*magSf_all[:nIF]*deltaCoeffs
        kfv_n=nuEff[nei]*magSf_all[:nIF]*deltaCoeffs
    lower=-qp*phi_B+kfv_o; upper=qn*phi_B+kfv_n
    divp=np.zeros(N)
    np.add.at(divp,owner,phi_B); np.add.at(divp,nei,-phi_B)
    np.add.at(divp,bc_cell,phiB_B)
    diag=np.zeros(N)
    np.add.at(diag,owner,qp*phi_B); np.add.at(diag,nei,-qn*phi_B)
    diag-=divcoef*divp
    np.add.at(diag,owner,-kfv_o); np.add.at(diag,nei,-kfv_n)
    diag+=alpha*V
    sumOff=np.zeros(N)
    np.add.at(sumOff,owner,np.abs(upper)); np.add.at(sumOff,nei,np.abs(lower))
    ab=nuEff[bc_cell]*bmag*bdel
    ic=np.where(Ufixed_b,-ab,phiB_B)
    Db=diag.copy(); np.add.at(Db,bc_cell,np.abs(ic))
    D=np.maximum(np.abs(Db),sumOff)/ALPHA_REL
    return 1.0/(D/V), Db, sumOff

r0,Db,sumOff=rebuild()
def relerr(r):
    return np.abs(mob-r)/np.maximum(np.abs(mob),1e-300)
e0=relerr(r0)
print("[baseline rebuild] relL2=%.4e" % (np.linalg.norm(mob-r0)/np.linalg.norm(mob)))
r1,_,_=rebuild(interp_gamma=True)
print("[S3 interp gamma ] relL2=%.4e (delta %.2e)"
      % (np.linalg.norm(mob-r1)/np.linalg.norm(mob),
         np.linalg.norm(r1-r0)/np.linalg.norm(r0)))
r2,_,_=rebuild(divcoef=V[0])
print("[S2 V*divphi     ] relL2=%.4e (delta %.2e)"
      % (np.linalg.norm(mob-r2)/np.linalg.norm(mob),
         np.linalg.norm(r2-r0)/np.linalg.norm(r0)))
print()
print("[structure of the residual deviation]")
solid=alpha>1e6; trans=(alpha>1e2)&(~solid); fluid=alpha<=1e2
for nm,msk in (("solid(alpha>1e6)",solid),("transition(1e2..1e6)",trans),
               ("fluid(alpha<=1e2)",fluid),("bnd-adjacent",bndcell),
               ("interior",~bndcell)):
    print("  %-22s n=%5d meanRelErr=%.3e frac>1%%=%.3f"
          % (nm,msk.sum(),e0[msk].mean(),(e0[msk]>0.01).mean()))
clampbind=np.abs(Db)<sumOff
print("  clamp binds (sumOff>|D|): %.1f%% of cells; err mean=%.3e frac>1%%=%.3f"
      % (100*clampbind.mean(), e0[clampbind].mean(), (e0[clampbind]>0.01).mean()))
print("  diag-dominant            : err mean=%.3e frac>1%%=%.3f"
      % (e0[~clampbind].mean(), (e0[~clampbind]>0.01).mean()))
m=nut>1e-6
print("  nut>1e-6 (turbulent)     : n=%d err mean=%.3e frac>1%%=%.3f"
      % (m.sum(), e0[m].mean(), (e0[m]>0.01).mean()))
print("  nut<=1e-6 (molecular)    : n=%d err mean=%.3e frac>1%%=%.3f"
      % ((~m).sum(), e0[~m].mean(), (e0[~m]>0.01).mean()))
print()
print("verdict: S1 (boundary diag) is exactly represented by the old formula")
print("(relax's -cmptMin + D()'s +cmptAv cancel for uniform ic); S2/S3 change")
print("<0.5%%; the ~13%% deviation is concentrated in nut-active clamp-bound")
print("transition/fluid cells -- offline off-diagonal assembly semantics,")
print("not yet fully identified.")
