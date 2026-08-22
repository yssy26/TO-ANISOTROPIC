# BFINAL-026 NOTEBOOK (executor working log)

Authoritative: PREREGISTRATION.md (immutable). Discipline: ORIENTATION_FOR_EXECUTORS.md v1.
HEAD: 055db8c (agent/dsH-stage-b-validation). Case: /home/ys/dsH/b25_qgate.

## Key pins (B25 tsv, never rerun)
- FD_J: D1 +0.04244654087520727 (h=1e-5, signJ=0) | D2 -0.003699503574594587 (h=1e-4, signJ=1) | D3 -0.0004211702528400529 (h=1e-4, signJ=1)
- ADJ_J: D1 -0.0358350400485765 | D2 -0.02492599038981423 | D3 -0.01200301226991504
- projV anchors (volume): D1 0.1555186511144686 | D2 -0.08449667741658891 | D3 0.06396237684432686 (rel 1e-6..1e-7 with FD gV at small h)
- baselineJ=-0.1674504032463462; relJ D1 2.18499 D2 0.85158 D3 0.96293; signJ D1=0 D2=1 D3=1
- adjoint_identity.tsv: thermalCoupling gradProxy=-510.7047296429894 L2=2.451199924458097; pressureDrop gradProxy=-315867.3295218884 L2=97.47897160928267
- stage_b2_summary projJ_D1/2/3 = ADJ_J values above
- M_frozen=0.00407336999996414, Tref=600.0, NU=5.19009e-05, alphaRel=0.4, NV=100800, NUNK=134400, N=33600
- b18rhs_thermalCoupling.mtx: 33600x4 plain whitespace, PHYSICAL (STALE note: NO descaling per stageB18RhsExport L870-898), 35840 nnz, max 1.56193e-05

## Work log
### 2026-08-21
- Read ORIENTATION_FOR_EXECUTORS.md (v1). Discipline: NOTEBOOK + extract-read logs + one-compile + reuse archived data + background long runs.
- Read b18_folding_lab.py L1-560 (template): read_field parser, mesh build, route_V0/V1, wtrue (FD central), C_field, C_direct, build_dHbyA_drAU, Gx_of, part-2 decomposition. DEFECTS to replace in B26: vertex-mean centroids (L164-165), B16 FD pins (L38), B16 z/dJdT/outletMeta deps (L373-374, L235), B16 SB.
- Read stageB6RxDesignOracle.H: dAlphaDxh (L153-162), D1/D2/D3 (L192-234, pyramid C = mesh.C()), anRUa (L395-403), drAU/dHbyA (L407-425), dphiHbyA + dflux + anRPa (L426-519).
- Read primitiveMeshCellCentresAndVols.C + primitiveMeshFaceCentresAndAreas.C: EXACT pyramid centroid formula captured for offline replication.
- D formula detail: Lx=Ly=Lz = gMax(C.comp)-gMin(C.comp) of cell centres (NOT centered); xx=Cx/Lx etc; active=designMask>0.5; x<0.02||x>0.98 -> 0; max-norm normalize (global max mag).

### Source pinning (this session, pre-instrument)
- dAlphaDxh (stageB6 L153-162): activeCell = designMask>0.5; unclipped = activeCell && alpha > alphamin*(1+1e-10); dAlphaDxh = -alphaMax*(1+qu)*qu/sqr(qu+xh+SMALL) on unclipped cells only.
- D1/D2/D3 (stageB6 L192-234): C=mesh.C() pyramid centroids; Lx=gMax(C.x)-gMin(C.x), same Ly/Lz (NOT centered); xx=C.x/max(Lx,SMALL); D1=sin(2pi*xx)cos(pi*yy)+0.5sin(3pi*zz); D2=sin(4pi*xx)+0.3cos(2pi*xx); D3=cos(2pi*yy)sin(2pi*zz); only activeCell (designMask>0.5); if x<0.02||x>0.98 -> 0 (physical x); global max-norm normalize.
- EXACT pyramid centroid (primitiveMeshCellCentresAndVols.C L95-140): cEst[cell] = mean of its face centres; internal faces: pyr3Vol = fAreas[facei] & (fCtrs[facei] - cEst[own]); pc = (3/4)*fCtrs[facei] + (1/4)*cEst[own]; cellCtrs[own] += pyr3Vol*pc; cellVols[own] += pyr3Vol. Neighbour: pyr3Vol = fAreas[facei] & (cEst[nei] - fCtrs[facei]); pc = (3/4)*fCtrs[facei] + (1/4)*cEst[nei]. FINAL cellCtrs /= cellVols.
- Face centres/areas (primitiveMeshFaceCentresAndAreas.C L55-110): triangle -> fCtrs=(1/3)sum p, fAreas=0.5*((p1-p0)^(p2-p0)); polygon -> fCentre=vertex mean; sumN += (next-pi)^(fCentre-pi); sumA += mag; sumAc += a*(pi+next+fCentre); fCtrs=(1/3)*sumAc/sumA; fAreas=0.5*sumN. Boundary faces: only own-cell side in pyramid loop.
- diff.c eta functional (src/diff.c, 47 lines): z = sum_{designMask>0.5}(gamma - projected)*V, two-branch exponential with eta, del.
- filterTangent (stageB6 L943-982): (M-B)y = -b*d via fvm::laplacian(designFilterFaceMask, y) - fvm::Sp(b, y) + dfield*b; PCG/DIC tol 1e-9.
- projectionTangent (stageB6 L986-1000): num = sum designMask*(1-drhoF)*V*y; etaResp = safeEtaDivide(num) (raw denom, mag<=SMALL->0); z = drhoF*y + dProjEta*etaResp.
- drhoF/dProjEta (stageB6 L890-920): xp<=eta5: pe=exp(-del*(1-xp/eta5)); drhoF=del*pe+exp(-del); dProjEta=pe*(1-del*xp/eta5)-exp(-del). else branch symmetric. projEtaDenom = sum designMask*dProjEta*V (raw, may be negative; old SMALL-clamp was WRONG ~3.14e9 amplification).
- analyticRxDirection (stageB6 L1057-1078): w = dAlphaDxh*z; anRUd = anRUa*w (= V*U*w per component); anRPd = assembleWeightedRPa(w) (weight inside face-flux assembly).
- anRUa = V*U (per cell component); drAU = -rAU_rel^2/ALPHAREL; dHbyA = (rAU_rel/ALPHAREL)*((1-ALPHAREL)*U - HbyA) (b18_folding_lab L462-513 validated).
- assembleWeightedRPa (stageB6 L426-519): dphiHbyA_w = face-interp of dHbyA*w; dflux_w = face-interp of drAU*w; anRPd = div(dphiHbyA_w - dflux_w) (P-row of R_x).
- Templates: b23_t9_abform_corrected.py (BFINAL-023, 205 lines) = corrected-symbols basis builder (build_basis L78-99, mobility gate L143, phi_map L114-132, grad_p L101-109, patch sizes via regex from polyMesh/boundary). DEFECT: SB points to b22_gauge/stageB2 -> MUST repoint to b25_qgate/stageB2.
  b18_folding_lab.py (BFINAL-018, 1199 lines) = M1/M3/M4 lab template. Reusable: read_field (L44-127), mesh build (L130-178), g face functional (L229-240), C_field (L376-408), build_dHbyA_drAU (L462-513), Gx_of (L522-533), route_V1_g (L1064-1087), thermal operator assembly (L661-764), rebuild_flux (L964-1010), block regression (L1141-1199). DEFECTS to NOT replicate: L164-165 vertex-mean centroids; L38 B16 FD pins; L373-374 uses E15=b15_export z3 analytic; L235 ST=b16_states/stageB2; L181-189 hardcoded B16 patch faces.
- w_true template: b21_wstar_clean.py (BFINAL-021, 72 lines): M = mmread(b8_verify_diag/explicitJT.mtx).tocsc(); sort_indices(); lu=spla.splu(M); w=lu.solve(-rxd, trans='T'); closure r_P=(M.T@w)[P0:]+rxd[P0:], rel-norm ~1e-11, BLOCKED if >1e-9; ~25 min.

### CRITICAL: b8_verify_diag == b25_qgate (state identity)
- /home/ys/dsH/b8_verify_diag/1/ bit-identical to /home/ys/dsH/b25_qgate/1/ for ALL compared List fields (alpha, U, p, phi, x, xh, xp, designMask, T, Tb, phiThermal, dAlphaDxh, dDTDxh: maxdiff=0.000e+00, relL2=0.000e+00); nuEffFrozen uniform 5.19009e-05 both.
- Therefore b8_verify_diag exports serve as production-side artifacts at the B25 state:
  - explicitJT.mtx (M=J^T, b21-validated closure ~1e-11)
  - stageB6_rxc_analytic.mtx: 403200 lines = 3x134400 (rxc = R_x*D for D1/D2/D3)
  - stageB6_rxc_z_analytic.mtx: 100800 = 3x33600 (reshape(3,N) = z per dir)
  - stageB6_rxc_xh_analytic.mtx: 100800 (=NV), stageB6_rxa_analytic.mtx: 134400, stageB6_rpa_terms.mtx: 100800
  - stageB6_lambda.mtx: 134400 (=NUNK), stageB6_dirs.mtx: 100800 (=3x33600 = D1/D2/D3)
  - stageB6_dalphadxh.mtx: 33600 (=N)
  - Plain-whitespace line format (first values "0\n0\n-0"). B25 itself has NO rxc/z/dirs exports -> b8_verify_diag fills the gap at identical state.
