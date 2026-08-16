#!/bin/bash
echo "=== b2_case explicit/py files ==="
ls -la /home/ys/b2_case/ | grep -iE 'explicit|\.py|\.mtx' || echo "(none)"
echo "=== b2_case_smoke explicit/py files ==="
ls -la /home/ys/b2_case_smoke/ | grep -iE 'explicit|\.py|\.mtx' || echo "(none)"
echo "=== py files in repo ==="
find /home/ys/TO-ANISOTROPIC -name '*.py' 2>/dev/null | head -30
echo "=== b4 log oracle lines ==="
grep -iE 'ExplicitJToracle|Loaded explicit|relRes|oracle|Explicit' /home/ys/b2_case/log.b4 2>/dev/null | head -40 || echo "(no log.b4)"
echo "=== python3 available? ==="
python3 -c 'import scipy; print("scipy", scipy.__version__)' 2>&1 | head -3
