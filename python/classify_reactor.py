"""Classify the OTUs of a reactor into AOB / NOB (port of Classification/main.m).

usage: python classify_reactor.py [A|B] [l1|l2]

Writes results/classification_<reactor>_<objective>.npz with full-length (44) vectors
class_AOB, class_NOB, yields_AOB, yields_NOB (same layout as the MATLAB .mat files).
For reactor A it also checks the inputs against the saved MATLAB workspace.
Times are absolute days (B.xls starts at day 5; MATLAB shifted B by those 5 days).
"""
import sys
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.io import loadmat

sys.path.insert(0, str(Path(__file__).parent))
from class_and_track.classification import classify
from class_and_track.data import REPO, classification_inputs

reactor = sys.argv[1] if len(sys.argv) > 1 else "B"
objective = sys.argv[2] if len(sys.argv) > 2 else "l2"     # l2 = MIP_gurobi.m (used for A)
yA_ref, yB_ref, delta = 0.2229, 0.0215, 0.01               # Classification/main.m

ci = classification_inputs(reactor)
print(f"reactor {reactor}: {ci['W'].shape[0]} OTUs considered, {ci['W'].shape[1]} time points")

if reactor == "A":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = loadmat(REPO / "Classification" / "MAT_files" / "Classification_day_183_315_gurobi.mat",
                    squeeze_me=True)
    for k in ("W", "xG1", "xG2"):
        ref = c["partial_weights"] if k == "W" else c[k]
        print(f"  check {k} vs MATLAB: max abs diff = {np.max(np.abs(ci[k] - ref)):.2e}")

tic = time.time()
out = classify(ci["W"], ci["xG1"], ci["xG2"], yA_ref, yB_ref, delta, objective=objective)
print(f"  {objective}: {time.time() - tic:.1f}s  status={out['status']}  obj={out['objective']:.6g}")

n, cons = ci["n_otu"], ci["considered"]
res = {k: np.zeros(n) for k in ("class_AOB", "class_NOB", "yields_AOB", "yields_NOB")}
res["class_AOB"][cons], res["class_NOB"][cons] = out["a"], out["b"]
res["yields_AOB"][cons], res["yields_NOB"][cons] = out["y_AOB"], out["y_NOB"]
print(f"  AOB ({int(res['class_AOB'].sum())}):", " ".join(f"E{i + 1}" for i in np.flatnonzero(res["class_AOB"])))
print(f"  NOB ({int(res['class_NOB'].sum())}):", " ".join(f"E{i + 1}" for i in np.flatnonzero(res["class_NOB"])))

outdir = Path(__file__).parent / "results"
outdir.mkdir(exist_ok=True)
np.savez(outdir / f"classification_{reactor}_{objective}.npz", **res,
         error_AOB=out["error_AOB"], error_NOB=out["error_NOB"], objective=out["objective"])
