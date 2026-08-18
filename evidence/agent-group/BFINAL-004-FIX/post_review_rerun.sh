#!/usr/bin/env bash
# Post-Reviewer INDEPENDENT rerun of the load-bearing diagnostic.
# Re-runs MTO_HF on the writable b2_case_smoke with stageB6RxDesignOracle=true.
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc >/dev/null 2>&1
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin
unset FOAM_SIGFPE
cd /home/ys/dsH/b2_case_smoke || exit 90
BIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF
echo "POSTREVIEW_RUN_CWD=$(pwd)"
echo "POSTREVIEW_BIN=$BIN"
echo "POSTREVIEW_BIN_SHA=$(sha256sum "$BIN" | cut -d' ' -f1)"
echo "POSTREVIEW_START=$(date -u +%FT%TZ)"
"$BIN"
rc=$?
echo "POSTREVIEW_MTO_RC=$rc"
echo "POSTREVIEW_END=$(date -u +%FT%TZ)"
exit $rc
