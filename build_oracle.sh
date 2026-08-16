#!/bin/bash
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
unset FOAM_SIGFPE
cd /home/ys/TO-ANISOTROPIC/src || exit 1
wmake > /home/ys/TO-ANISOTROPIC/src/validation_stage_b2diag_build3.log 2>&1
rc=$?
echo "BUILD_EXIT=$rc"
tail -40 /home/ys/TO-ANISOTROPIC/src/validation_stage_b2diag_build3.log
exit $rc
