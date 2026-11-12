"""Check the Python port against the saved MATLAB runs in Tracking/MAT_files.

1. open-loop simulation (iter 0) vs MATLAB 191218_Reactor_A_iter_0.mat
2. one ASRE step fed with MATLAB's iter-0 trajectory vs MATLAB iter 1
"""
import sys
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.io import loadmat

sys.path.insert(0, str(Path(__file__).parent))
from class_and_track.asre import asre_step, solve
from class_and_track.data import REPO, lin, problem_from_classification_mat


def mat(name):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return loadmat(REPO / "Tracking" / "MAT_files" / name, squeeze_me=True)


def report(label, t_ref, X_ref, fun, n):
    X = fun(t_ref).T
    scale = np.abs(X_ref).max(0)
    scale[scale == 0] = 1
    err = np.abs(X - X_ref) / scale
    print(f"{label}: max rel. error  biomass {err[:, :n].max():.2e}   substrates {err[:, n:].max():.2e}")


m0, m1 = mat("191218_Reactor_A_iter_0.mat"), mat("191218_Reactor_A_iter_1.mat")
prob = problem_from_classification_mat(
    REPO / "Classification" / "MAT_files" / "Classification_day_183_315_gurobi.mat",
    x0_from="row0", muA=float(m0["muARef"]), muB=float(m0["muBRef"]),
    kSA=float(m0["kSARef"]), kSB=float(m0["kSBRef"]))
M, n = prob.model, prob.model.n
assert np.allclose(prob.x0, m0["x0"]), "x0 differs from MATLAB"
print(f"n = {n} (nA={M.nA}, nB={M.nB}), horizon [{prob.t0}, {prob.tf}] days")

tic = time.time()
ol = solve(lambda t, X: M.rhs(t, X), (prob.t0, prob.tf), prob.x0, rtol=1e-8, atol=1e-11).sol
print(f"open loop solved in {time.time() - tic:.1f}s")
report("iter 0 (open loop)", m0["t_new"], m0["x_new"], ol, n)

prev_matlab = lin(m0["t_new"].astype(float), m0["x_new"])
r = np.broadcast_to(np.asarray(m1["lambda"], float), (n,))
q = np.broadcast_to(np.asarray(m1["Q"], float), (n,))
for vmax in (0.5, 0.2, 0.0):
    tic = time.time()
    step = asre_step(M, prev_matlab, prob.z, prob.x0, prob.t0, prob.tf, q, r, np.zeros(n),
                     v_bounds=(-1.0, vmax), rtol=1e-8, atol=1e-11)
    print(f"ASRE step (v_max={vmax}) solved in {time.time() - tic:.1f}s")
    report(f"iter 1, v in [-1,{vmax}]", m1["t_new"], m1["x_new"], step.traj, n)
