#!/bin/bash
# Sequential experiment runner (parallel `conda run` calls race on activation).
C="$USERPROFILE/miniconda3/Scripts/conda.exe"
cd "$(dirname "$0")"
run() { "$C" run --no-capture-output -n class-and-track python -u "$@" 2>&1 | grep -v -E "^\s*[0-9]+r?\s+[-0-9.e+]+\s" ; }
for spec in "t0 15 -1 0.5 1e-5" "t0 15 -1 inf 1e-5" "t0 15 -1 0.5 1e-3"; do run run_asre_reactor_A.py $spec; done
run run_direct_reactor_A.py -1 0.5 1e-5
