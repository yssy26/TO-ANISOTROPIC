# BFINAL-005 cycle-2 — step4_verify: Gate table G1–G5 (independent)

All values independently recomputed by `verify_gates.py` from the raw
artifacts; G4/G5 anchors re-grepped directly from `run_stageb6_cycle2.log`
(raw log lines quoted below, not taken from the step3 summary).

## G1 — local pressure-row contraction: D_pressure_production vs D_pressure_oracle

| dir | D_pressure_production (Σ gsensh_pressurerow·z) | D_pressure_oracle (-Σ pc·R_P,xd) | relErr | sign | PASS |
|---|---|---|---|---|---|
| D1 | -0.00754887795106 | -0.00754887795106 | 4.37e-15 | same | PASS |
| D2 | +0.0685927693909 | +0.0685927693909 | 2.02e-15 | same | PASS |
| D3 | -0.00621002821874 | -0.00621002821874 | 4.55e-14 | same | PASS |

Threshold: relErr ≤ 1e-3 (preferred 1e-4) AND sign agreement. **PASS**
(preferred 1e-4 achieved by ~10 orders of magnitude).

## G2 — full R_x contraction: D_total_production vs D_total_oracle

| dir | D_total_production (Σ gsensh_total·z) | D_total_oracle (D_mom+D_press) | relErr | sign | \|D_pressure/D_total\| | PASS |
|---|---|---|---|---|---|---|
| D1 | -0.0562356333787 | -0.0562356333787 | 1.23e-15 | same | 0.134237 | PASS |
| D2 | +0.297672431753 | +0.297672431753 | 2.42e-15 | same | 0.230430 | PASS |
| D3 | -0.0363414298395 | -0.0363414298395 | 1.26e-13 | same | 0.170880 | PASS |

D_momentum component also closes: relErr 2.7e-15 / 4.8e-16 / 1.2e-13.
Anchors match BFINAL-004 §9/§10 exactly (D_momentum −0.0486867554277 /
+0.229079662362 / −0.0301314016208; D_pressure −0.00754887795106 /
+0.0685927693909 / −0.00621002821874; D_total −0.0562356333787 /
+0.297672431753 / −0.0363414298395). **PASS** (preferred 1e-4 achieved).

## G3 — raw-design projection/filter closure

| dir | D_total_raw (Σ rxpr_prod_gsens·d) | D_total_oracle | relErr | sign | PASS |
|---|---|---|---|---|---|
| D1 | -0.0562356332429 | -0.0562356333787 | 2.42e-09 | same | PASS |
| D2 | +0.297672431913 | +0.297672431753 | 5.38e-10 | same | PASS |
| D3 | -0.0363414298017 | -0.0363414298395 | 1.04e-09 | same | PASS |

The patched production raw-design field (single `gsenshPressureDrop` total fed
through the EXISTING filter_chainrule.H) projects onto the raw-design
directions to ~1e-9 of the oracle total. Volume projection UNCHANGED (see G5).
**PASS** (preferred 1e-4 achieved).

## G4 — R_w/J no regression (re-grepped from raw log)

| anchor | locked (BFINAL-003) | run log value (line) | PASS |
|---|---|---|---|
| momentum relL2 | 1.653e-4 | 1.65311363075e-4 (L1291 GateA eps=0.001, `relL2=0.000165311363075 cos=0.999999986345`) | PASS |
| P-total relL2 | 1.53e-4 | 1.52828911949e-4 (L1294 `relL2=0.000152828911949 cos=0.999999988323`) | PASS |
| J/J^T dot | 6.42e-14 | 6.42244770078e-14 (L1394 `Reduced cold-flow operator transpose dot-test (pressureDrop)`) | PASS |
| explicit oracle maxRelL_U | 3.53e-16 | 3.52693970502e-16 (L1266/L1656 `ExplicitJToracle: ... maxRelL_U=3.52693970502e-16`) | PASS |
| explicit oracle maxRelL_P | 3.98e-16 | 3.98062510668e-16 (L1266 `... maxRelL_P=3.98062510668e-16`) | PASS |
| GatePR | BFINAL-003 anchors | `GatePR h=0.001 |J_P|=4.33741818206e-06 |FD_P|=6.2437647062e-06 relL2=1.03546813403 cos=0.306822860167` (L1234) — identical across both occurrences (L1233-1249 and L1623-1639) | PASS |

Also: thermalCoupling transpose dot 6.42244770078e-14 (L1004) — same as
pressureDrop; BlockDot (thermalCoupling) PU=2.078e-13 PP=2.666e-14 UU=7.310e-14
UP=1.796e-12 boundaryU=1.804e-13 (L1201) — unchanged; Rx-A relL2=1.6135e-14
cos=1 (L147). **PASS** — R_w/J remains frozen from BFINAL-003.

## G5 — non-flow anchors (re-grepped from raw log)

| anchor | locked | run log value | PASS |
|---|---|---|---|
| volume proj D1 | 0.1555186511144686 | 0.1555186511144686 (L1868) — independent dot 0.1555186511144686 | PASS |
| volume proj D2 | -0.08449667741658891 | -0.08449667741658891 (L1869) — independent dot identical | PASS |
| volume proj D3 | 0.06396237684432686 | 0.06396237684432686 (L1870) — independent dot identical | PASS |
| filter/projection | stageB6 RX-C tangent closure | xh-tangent D1 eps=1e-4 relL2=2.142e-05 cos=0.999999999786 (L1789); R_U,xd D1 relL2=2.298e-05 cos=0.999999999752; R_P,xd D1 relL2=4.141e-05 cos=0.999999999157 (L1792) | PASS |
| Stage-A anisotropy | solidK=(21.6 0 0 21.6 0 16.6) | `Anisotropic conductivity: solidK=(21.6 0 0 21.6 0 16.6)` (L96) — unchanged | PASS |
| discrete objective derivative | thermal/pressure ~1e-11/1e-12 | `Discrete objective derivative errors: thermal=2.7025486001e-11, pressure=7.07793481308e-12` (L122) | PASS |
| projection eta | eta5=0.748005161603 del=8 | `StageB6 RX-C projection: eta5=0.748005161603 del=8 projectionEtaDenominator=-3.1385996818e-06` (L1773) | PASS |

**PASS** — volume/filter/projection/Stage-A/discrete-objective all unchanged.

## Summary

| Gate | verdict |
|---|---|
| G1 pressure-row contraction | **PASS** (4.4e-15 … 4.6e-14) |
| G2 total R_x contraction | **PASS** (1.2e-15 … 1.3e-13) |
| G3 raw-design projection | **PASS** (5.4e-10 … 2.4e-9) |
| G4 R_w/J no regression | **PASS** (all BFINAL-003 anchors reproduced) |
| G5 non-flow anchors | **PASS** (all unchanged, volume proj byte-identical) |

ALL GATES PASS. Production `gsenshPressureDrop` now contains both
momentum-row and pressure/continuity-row design derivatives
(gsenshPressureDrop = gsenshPressureDropMomentum + gsenshPressureDropPressureRow,
with gsenshPressureDropPressureRow = -rxPressureRowT*dAlphaDxh = -(J_P^T pc)*dAlphaDxh),
and the total reproduces the BFINAL-004 locked analytic R_x to machine precision
at the xh level and ~1e-9 through the production filter/chain path.
