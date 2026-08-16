#!/bin/bash
echo "=== python3 location ==="
which python3
python3 --version 2>&1
echo "=== numpy / scipy in WSL python3 ==="
python3 -c 'import numpy; print("numpy", numpy.__version__)' 2>&1 | head -2
python3 -c 'import scipy; print("scipy", scipy.__version__)' 2>&1 | head -2
python3 -c 'import scipy.sparse.linalg as s; print("splu ok")' 2>&1 | head -2
echo "=== pip3 list (numpy/scipy) ==="
pip3 list 2>/dev/null | grep -iE 'numpy|scipy' || echo "(pip3 list failed or none)"
echo "=== pip available? ==="
which pip3 pip 2>&1
echo "=== Windows python? ==="
ls /mnt/c/Windows/py.exe 2>/dev/null && echo "py.exe exists" || echo "no py.exe"
ls /mnt/c/Python*/python.exe 2>/dev/null || echo "no C:\\Python"
ls /mnt/c/Users/admin/AppData/Local/Programs/Python 2>/dev/null || echo "no user python"
