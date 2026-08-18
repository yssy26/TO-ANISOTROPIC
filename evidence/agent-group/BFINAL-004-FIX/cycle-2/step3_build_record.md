# BFINAL-004-FIX cycle-2 — step3-build record

Date: 2026-08-17
Workspace: /home/ys/dsH/TO-ANISOTROPIC (git root verified via `git rev-parse --show-toplevel` = /home/ys/dsH/TO-ANISOTROPIC)
HEAD: ca8a772b8339a36e7906ea32a7125061c47aa280 (unchanged)

## Environment (exactly per plan validation_commands)

```
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin
unset FOAM_SIGFPE
cd /home/ys/dsH/TO-ANISOTROPIC/src && wmake
```

## Result

- **WMAKE_EXIT = 0**
- Full log: `evidence/agent-group/BFINAL-004-FIX/cycle-2/build_step3.log`
- Only `MTO_HF.C` recompiled and relinked (header-only incremental; probe header
  `src/stageB6RxDesignOracle.H` is included by MTO_HF.C, tracked in
  `Make/linux64GccDPInt32Opt/MTO_HF.C.dep`, 2 references).
- No other object file rebuilt (find for `*.o ! MTO_HF.o` newer than binary → empty).

## Binary state (acceptance)

| item | value |
|---|---|
| binary path | /home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF |
| size | 4747104 bytes |
| mtime | 2026-08-17 03:23:23 (updated from 02:46 by this step's wmake) |
| object | src/Make/linux64GccDPInt32Opt/MTO_HF.o, 7964192 bytes, 03:23:23 |

## Notes

- First step3 wmake attempt (before touching the header) reported WMAKE_EXIT=0 with
  empty log: the step2 sanity build (02:46) had already compiled the fixed header
  (02:23), so the tree was up-to-date.
- To produce unambiguous evidence that wmake recompiles MTO_HF.C against the current
  fixed header and relinks, the probe header was `touch`ed (content unchanged:
  sha256 still `09832dd6ae5a3c7f37201b200bfb19dc0c48f42a3989f9593979c73ca098baec`,
  identical to step2 record) and wmake rerun. This forced a genuine recompile+relink
  with fresh mtime 03:23:23.
- Compilation warnings are pre-existing unused-variable warnings from production
  headers (NS.H, solveDiscreteFlowAdjoint.H, createFields.H); no errors.
- Probe header content unchanged by the touch; production headers untouched:
  - solveDiscreteFlowAdjoint.H        sha256 f0c814975cf8bc7f98b9f69efa2f0bc01bd251579920bbedda7d266b72927507 (== cycle-2/pre-state)
  - solveDiscreteFlowAdjointProduction.H sha256 8c4901bf8da4edbc7cc1efe73a7583963e4688bba93c8c58539c4c39cae1ad91 (== cycle-2/pre-state)

## Verification of step2 fix presence in current header (before build)

- `assembleWeightedRPa` defined at src/stageB6RxDesignOracle.H:557; called at
  L686 (RX-A w=deltaAlpha), L795 (RX-B w=dAlphaDxh*deltaXh), L1077 (RX-C w=dAlphaDxh*z).
- `Foam::max(projEtaDenom,SMALL)` occurrences: none (D4 fixed via safeEtaDivide).
- Post-build binary is the executable to be used by step4-run-report.
