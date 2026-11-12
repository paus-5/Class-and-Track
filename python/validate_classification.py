"""Solve the saved classification instances with free solvers and compare to MATLAB."""
import sys
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.io import loadmat

sys.path.insert(0, str(Path(__file__).parent))
from class_and_track.classification import classify
from class_and_track.data import REPO

for tag, objective in [("intlinprog", "l1"), ("gurobi", "l2")]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = loadmat(REPO / "Classification" / "MAT_files" / f"Classification_day_183_315_{tag}.mat",
                    squeeze_me=True)
    tic = time.time()
    out = classify(c["partial_weights"], c["xG1"], c["xG2"], float(c["yA_ref"]), float(c["yB_ref"]),
                   float(c["delta"]), objective=objective)
    dt = time.time() - tic
    same_a = np.array_equal(out["a"], np.round(c["a"]))
    same_b = np.array_equal(out["b"], np.round(c["b"]))
    matlab_obj = (np.sum(c["error_AOB"]) + np.sum(c["error_NOB"]) if objective == "l1"
                  else np.sum(c["error_AOB"] ** 2) + np.sum(c["error_NOB"] ** 2))
    print(f"{tag:10s} ({objective}, {'HiGHS' if objective == 'l1' else 'SCIP'}): {dt:5.1f}s  "
          f"status={out['status']!s:.40}  obj={out['objective']:.6g}  MATLAB obj={matlab_obj:.6g}  "
          f"nAOB={int(out['a'].sum())} nNOB={int(out['b'].sum())}  same class as MATLAB: {same_a and same_b}")
    if not (same_a and same_b):
        print("   OTUs that differ (index within considered OTUs):",
              np.flatnonzero((out["a"] != np.round(c["a"])) | (out["b"] != np.round(c["b"]))))
