#!/bin/bash
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
unset FOAM_SIGFPE
cd /home/ys/b2_case_smoke || exit 1

echo "=== write cell centres ==="
postProcess -func writeCellCentres -time 1 > /dev/null 2>&1
echo "C field:"; ls -la 1/C 2>&1

echo "=== verify gsensPressureDrop / designMask / x written ==="
ls -la 1/gsensPressureDrop 1/designMask 1/x 2>&1
