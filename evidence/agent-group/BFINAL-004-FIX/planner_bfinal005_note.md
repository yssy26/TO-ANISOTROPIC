# Planner note — BFINAL-005 (PATCH_RX_PRESSURE_ROW) authorization check
- BFINAL-004 DIAGNOSTIC_ONLY PASSED (cycle-2). MISSING_PRESSURE_RX=CONFIRMED.
- task/BFINAL-004.md §12 Case 2 authorizes next task BFINAL-005 PATCH_RX_PRESSURE_ROW.
- Review feedback next_route = PATCH_RX_PRESSURE_ROW (authorize BFINAL-005).
- Verified: production sensitivity.H L67-68 gsenshPressureDrop = -dAlphaDxh*(U&Uc)*V (momentum-only).
- Probe chain (stageB6RxDesignOracle.H): drAU=-rAU^2/alphaRel; dHbyA=(rAU/alphaRel)((1-alphaRel)U-HbyA);
  assembleWeightedRPa(w)=div(dphiHbyA_w)-div(dflux_w); D_pressure=-pc^T R_P,x d.
- BFINAL-004 reference D_pressure (D1/D2/D3): -0.00754887795, +0.0685927694, -0.00621002822.
- B-final Δp_c gate = validateStageB2GradientAmplitude.H (stageB2Enabled=true), gDP: sign + bestRelDP<=10% + plateau.
- Current optProperties: stageB2Enabled=false, stageB4JacobianProbe=true, gradientValidated=true, mmaUpdateEnabled=false.
