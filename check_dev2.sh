#!/bin/bash
echo "=== dev2 in TensorI.H ==="
grep -n 'dev2' /opt/openfoam7/src/OpenFOAM/primitives/Tensor/TensorI.H | head -10
echo "=== dev / dev2 in Tensor.H ==="
grep -rn 'dev2\|dev(' /opt/openfoam7/src/OpenFOAM/primitives/Tensor/Tensor.H | head -20
echo "=== spherical twoThirds / oneThird ==="
grep -rn 'twoThirds\|oneThird' /opt/openfoam7/src/OpenFOAM/primitives/SphericalTensor/SphericalTensorI.H | head -10
