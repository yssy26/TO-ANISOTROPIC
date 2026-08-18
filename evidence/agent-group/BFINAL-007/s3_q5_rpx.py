#!/usr/bin/env python3
"""BFINAL-007 S3 Q5: does the outlet fixedValue boundary term affect R_P,x?
R_P,x (BFINAL-004/005 closed, rxPressureRowTranspose.H) = J_P*(dAlphaDxh*z)
with J_P*w = div(dphiHbyA_w - dflux_w), dflux_w = flux(fvm::laplacian(drAU*w,p)).
The drAU*w field is created with default calculated BC value 0
(rxPressureRowTranspose.H L36-39; stageB6 assembleWeightedRPa L453-527),
so the interpolated boundary gamma_b = 0 and dflux_b == 0 EXACTLY in the
operator as implemented -- the runtime dot test (2.58e-15) arbitrates.

Quantified here:
  (a) the production pEqn boundary-flux design derivative d(phi_b)/d(alpha)
      = +drAtU_b*deltaCoeffs_b*p_c*|Sf|_b  (drAtU = -rAU^2/alphaRel; the SAME
      boundary-face-flux mechanism as the state Jacobian term)
  (b) its magnitude vs the current dflux_b = 0 assumption and vs |J_P*w|.
"""
import numpy as np, json
ART="/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-007/artifacts"
CASE="/home/ys/dsH/b2_case_smoke"
N=33600
def read_of(path,n):
    t=open(path).read(); i=t.index('(')
    return np.array([float(x) for x in t[i+1:t.rindex(')')].split()[:n]])
op=np.loadtxt(ART+"/stageB7_outlet_patch.mtx")   # cell delta magSf ic bc mobB mobC cand
ocells=op[:,0].astype(int)
delta=op[:,1]; magSf=op[:,2]
rAU=np.loadtxt(ART+"/stageB7_rebuilt_rAU.mtx")
p=read_of(CASE+"/1/p",N)
alphaRel=0.4
drAU=-rAU**2/alphaRel          # BFINAL-004 locked per-cell design derivative
# (a) production boundary-flux design derivative per outlet cell
dphi_b_da = drAU[ocells]*delta*magSf*p[ocells]
res={"Q5_alphaRel":alphaRel,
     "Q5_outlet_ncells":int(len(ocells)),
     "Q5_drAU_outlet_range":[float(drAU[ocells].min()),float(drAU[ocells].max())],
     "Q5_pc_outlet_range":[float(p[ocells].min()),float(p[ocells].max())],
     "Q5_delta_b_range":[float(delta.min()),float(delta.max())],
     "Q5_magSf_b":float(magSf[0])}
res["Q5_dphi_b_da_min"]=float(dphi_b_da.min())
res["Q5_dphi_b_da_max"]=float(dphi_b_da.max())
res["Q5_dphi_b_da_L2"]=float(np.linalg.norm(dphi_b_da))
res["Q5_dphi_b_da_meanabs"]=float(np.abs(dphi_b_da).mean())
# (b) current dflux_b=0: boundary value of drAU*w field is exactly 0 (calculated BC)
res["Q5_current_dflux_b"]=0.0
res["Q5_current_dflux_b_reason"]=("drAU*w field created with dimensionedScalar 0 + "
    "calculated BC; correctBoundaryConditions() leaves value 0 (calculatedFvPatchField"
    "::updateCoeffs empty); interpolated boundary gamma_b=0 -> dflux_b==0 exactly;"
    " runtime dot test relErr=2.57680498181e-15 (Log.stageB7.on.txt L1900)")
# magnitudes for comparison (from probe run log)
res["Q5_ref_Jp_w_L2"]=6.22344977639e-10       # |J_P*deltaAlpha|L2 (RX-A weighted total_w)
res["Q5_ref_RxT_L2"]=0.000256656999292        # |T|L2 (rxPressureRowTranspose dot test)
res["Q5_ref_g0_bnd_max"]=0.0650620425279      # |g0| boundary max (design face gradient)
res["Q5_ratio_boundaryterm_vs_Jp_w"]=float(np.linalg.norm(dphi_b_da)/6.22344977639e-10)
# note: deltaAlpha (design direction) has ZERO support on the 84 outlet cells
da=np.loadtxt(CASE+"/stageB6_deltaAlpha.mtx")
dm=read_of(CASE+"/1/designMask",N)
res["Q5_deltaAlpha_outlet_nnz"]=int((np.abs(da[ocells])>1e-300).sum())
res["Q5_designMask_outlet_max"]=float(dm[ocells].max())
res["Q5_note_design_support"]=("deltaAlpha (RX-A weighted direction) and designMask are "
    "ZERO on all 84 outlet cells -> the design-weighted boundary term is IDENTICALLY 0 "
    "in the validated directions; the dflux_b=0 assumption is exact for this case.")
print(json.dumps(res,indent=1))
print("\nQ5: d(phi_b)/d(alpha) = +drAtU_b*delta_b*|Sf|_b*p_c per outlet cell:")
print("  range [%.4e, %.4e], L2=%.4e"%(res["Q5_dphi_b_da_min"],res["Q5_dphi_b_da_max"],res["Q5_dphi_b_da_L2"]))
print("  |term|L2 / |J_P*w|L2 = %.3e  (current operator: dflux_b = 0 exactly)"%res["Q5_ratio_boundaryterm_vs_Jp_w"])
print("  deltaAlpha support on outlet = %d cells; designMask max on outlet = %.1f"
      %(res["Q5_deltaAlpha_outlet_nnz"],res["Q5_designMask_outlet_max"]))
print("S3_Q5_DONE")
