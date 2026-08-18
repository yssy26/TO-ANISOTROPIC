#!/usr/bin/env bash
# BFINAL-005 step3-run: MTO_HF (fresh build with sensitivity.H PATCH_RX_PRESSURE_ROW)
# on the writable original /home/ys/dsH/b2_case_smoke with stageB6RxDesignOracle=true.
# Clean-PATH OpenFOAM-7 env (FOAM_USER_APPBIN exported AFTER sourcing bashrc),
# unset FOAM_SIGFPE, absolute-path binary (PATH would pick the stale
# /home/ys/OpenFOAM/ys-7/.../MTO_HF).
#
# NOTE: do NOT use `set -u` (OpenFOAM bashrc references unset vars -> shell aborts).
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc >/dev/null 2>&1
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin
unset FOAM_SIGFPE
cd /home/ys/dsH/b2_case_smoke || exit 90
BIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF
echo "RUN_CWD=$(pwd)"
echo "BIN=$BIN"
echo "BIN_SHA=$(sha256sum "$BIN" | cut -d' ' -f1)"
echo "START=$(date -u +%FT%TZ)"
"$BIN"
rc=$?
echo "MTO_RC=$rc"
echo "END=$(date -u +%FT%TZ)"
exit $rc
