#!/bin/bash
echo "=== binary ==="
stat -c '%y %n' /home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin/MTO_HF
echo "=== src file ==="
stat -c '%y %n' /home/ys/TO-ANISOTROPIC/src/solveDiscreteFlowAdjoint.H
echo "=== obj ==="
stat -c '%y %n' /home/ys/TO-ANISOTROPIC/src/Make/linux64GccDPInt32Opt/MTO_HF.o 2>/dev/null
echo "=== grep the fixed line in src (confirm fix present) ==="
grep -n '2.0/3.0)\*n\[i\]\*n\[j\]\*sf\[k\]' /home/ys/TO-ANISOTROPIC/src/solveDiscreteFlowAdjoint.H
grep -n 'twoThirdDelta\*sf\[k\]' /home/ys/TO-ANISOTROPIC/src/solveDiscreteFlowAdjoint.H
