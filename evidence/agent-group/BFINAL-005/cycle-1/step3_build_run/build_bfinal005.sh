#!/usr/bin/env bash
# BFINAL-005 step3-build: header-only incremental wmake of src/ (sensitivity.H
# patch) into FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin.
# Exact env per approved plan:
#   export PATH=(clean); source /opt/openfoam7/etc/bashrc;
#   export PATH=/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH;
#   export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin;
#   unset FOAM_SIGFPE; cd src && wmake
# NOTE: do NOT use `set -u` (OpenFOAM bashrc references unset vars -> abort).
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc >/dev/null 2>&1
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin
unset FOAM_SIGFPE
cd /home/ys/dsH/TO-ANISOTROPIC/src || exit 90
echo "WMAKE_START=$(date -u +%FT%TZ)"
echo "WM_PROJECT_DIR=$WM_PROJECT_DIR"
echo "FOAM_USER_APPBIN=$FOAM_USER_APPBIN"
echo "WM_OPTIONS=$WM_OPTIONS"
echo "WM_COMPILER=$WM_COMPILER"
wmake
rc=$?
echo "WMAKE_EXIT=$rc"
echo "WMAKE_END=$(date -u +%FT%TZ)"
exit $rc
