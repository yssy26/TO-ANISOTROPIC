# BFINAL-006 cycle-2 S1 raw scan evidence

matrix: /home/ys/dsH/b2_case_smoke/explicitJT.mtx
dims: 134400 x 134400, N=33600, pressure rows [100800, 134400)

## identity rows detected (single |v|>1e-12 AND diagonal==1): [100800]

## row 100800 entries (pinned identity row):
  col 0 : 0.0
  col 1 : 0.0
  col 2 : 0.0
  col 3 : 0.0
  col 4 : 0.0
  col 5 : 0.0
  col 240 : 0.0
  col 241 : 0.0
  col 242 : 0.0
  col 3362 : 0.0
  col 100800 : 1.0
  col 100801 : 0.0
  col 100880 : 0.0
  col 101920 : 0.0

## row 106400 entries (physical continuity row):
  col 13440 : 6.970983468546996e-08
  col 13441 : 3.469018402950499e-11
  col 13442 : 1.2241258627846214e-07
  col 16800 : -1.2506929384968876e-07
  col 16801 : -1.2499987363493728e-07
  col 16802 : 2.373068031106467e-09
  col 16803 : -1.2484482570497253e-07
  col 16804 : -2.901084429996748e-11
  col 16805 : 5.570394721625742e-11
  col 17040 : -7.512582961638338e-11
  col 17041 : -1.250302931484936e-07
  col 17042 : -3.322265373829812e-11
  col 20162 : -1.250000000000001e-07
  col 105280 : 3.601577421815226e-10
  col 106400 : -1.899160629508871e-09
  col 106401 : 5.133736877383788e-10
  col 106480 : 5.127675200486909e-10
  col 107520 : 5.128616795402787e-10

## patch summary (0/p):
  outlet : fixedValue uniform 0
  zeroGradient patches: 7
  p.needReference() == False

## check results: {'C1': True, 'C2': True, 'C3': True, 'C4': True, 'C5': True}
