# BFINAL-004-FIX cycle-3 — step3-build record (executor)

Date: 2026-08-17 · Workspace: /home/ys/dsH/TO-ANISOTROPIC (git root verified) ·
HEAD ca8a772b8339a36e7906ea32a7125061c47aa280 (unchanged)

## Environment (exactly per approved plan validation_commands)

```
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin
unset FOAM_SIGFPE
cd /home/ys/dsH/TO-ANISOTROPIC/src && wmake
```

## Result

- **WMAKE_EXIT = 0** (START 2026-08-16T20:20:35Z → END 2026-08-16T20:43:55Z)
- Full log: `cycle-3/build_step3.log`
- Header-only incremental: only `MTO_HF.C` recompiled (MTO_HF.o 04:43) and
  relinked (binary 04:43); no other object rebuilt. The probe header
  `src/stageB6RxDesignOracle.H` was `touch`ed (content UNCHANGED, sha256 still
  `09832dd6ae5a3c7f37201b200bfb19dc0c48f42a3989f9593979ca73ca098baec` — identical
  to cycle-3/pre-state and cycle-2 post-fix) to force a genuine recompile against
  the audited header.

## Binary state

| item | value |
|---|---|
| binary | /home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF |
| size | 4747104 bytes |
| mtime | 2026-08-17 04:43 (fresh; object MTO_HF.o 04:43, 7964192 B) |
| fixed-header markers in binary | `J_P*deltaAlpha` / `degenerate(nB=0)` / `projectionEtaDenominator singular` present (grep count 3) |
| old pointwise-multiply string | `anRPa[celli]*da` absent (grep count 0) |

Warnings in the log are pre-existing unused-variable warnings from production
headers (createFields.H stageB2FDStepsMulti etc.); no errors.

## Scope integrity (post-build)

- `git status --porcelain=v1` identical to cycle-3/pre-state (only the two
  pre-existing BFINAL-003 production heads modified; hashes unchanged:
  solveDiscreteFlowAdjoint.H=f0c81497..., solveDiscreteFlowAdjointProduction.H=
  8c4901bf...). No production file touched. Probe header content unchanged by the
  touch (sha256 identical).

Proceeding to step4-run on the writable original /home/ys/dsH/b2_case_smoke
(stageB6RxDesignOracle=true, line 87 of constant/optProperties, verified).
