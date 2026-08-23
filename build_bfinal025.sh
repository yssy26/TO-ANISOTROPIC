#!/bin/bash
# BFINAL-025 build: wclean; wmake in src with OpenFOAM-7 env.
# Matches prior rounds: clean PATH, source OF7 bashrc, custom FOAM_USER_APPBIN,
# unset FOAM_SIGFPE, sequential wmake (no WM_NCOMPPROCS), no tee.
set -o pipefail
export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
. /opt/openfoam7/etc/bashrc
export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin
unset FOAM_SIGFPE
cd /home/ys/dsH/TO-ANISOTROPIC/src || exit 90
wclean
wmake
rc=$?
echo "WMAKE_EXIT=$rc"
exit "$rc"
