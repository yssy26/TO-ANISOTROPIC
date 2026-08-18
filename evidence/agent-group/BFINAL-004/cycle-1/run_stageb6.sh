#!/usr/bin/env bash
# BFINAL-004 step-3-run: MTO_HF with stageB4JacobianProbe+stageB5BoundaryRelaxOracle+stageB6RxDesignOracle
# on a writable copy of b2_case_smoke.  Clean-PATH OpenFOAM-7 env (FOAM_USER_APPBIN exported AFTER
# sourcing bashrc, see step-1 design note: settings.sh L163 overrides it), unset FOAM_SIGFPE.
#
# NOTE: do NOT use `set -u` (OpenFOAM bashrc references unset vars -> shell aborts),
# and invoke the fresh build bin by ABSOLUTE path: PATH resolution would pick the stale
# /home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin/MTO_HF (Aug 15, no stageB6).
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc >/dev/null 2>&1
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin
unset FOAM_SIGFPE
cd /home/ys/dsH/b2_case_smoke_b6 || exit 90
BIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF
echo "RUN_CWD=$(pwd)"
echo "BIN=$BIN"
echo "START=$(date -u +%FT%TZ)"
"$BIN"
rc=$?
echo "MTO_RC=$rc"
echo "END=$(date -u +%FT%TZ)"
exit $rc
