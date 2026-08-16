#!/bin/bash
WPY=/mnt/c/Users/admin/AppData/Local/Programs/Python/Python313/python.exe
echo "=== windows python version ==="
"$WPY" --version 2>&1
echo "=== numpy/scipy in windows python ==="
"$WPY" -c 'import numpy; print("numpy", numpy.__version__)' 2>&1 | head -2
"$WPY" -c 'import scipy; print("scipy", scipy.__version__)' 2>&1 | head -2
"$WPY" -c 'import scipy.sparse.linalg as s; print("splu/spsolve ok")' 2>&1 | head -2
echo "=== pip in windows python ==="
"$WPY" -m pip --version 2>&1 | head -2
