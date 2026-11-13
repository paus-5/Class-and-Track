"""Direct optimal control (CasADi + IPOPT) on Reactor A, same cost/bounds as the ASRE runs.
The optimal piecewise-constant control is re-simulated with scipy and J is evaluated with
the same routine used for ASRE, so the numbers are directly comparable.

usage: python run_direct_reactor_A.py [vmin] [vmax] [lambda]
"""
import pickle
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from class_and_track.asre import ASREStep, evaluate, solve
from class_and_track.data import REPO, problem_from_classification_mat
from class_and_track.direct import direct_tracking, piecewise_control

vmin = float(sys.argv[1]) if len(sys.argv) > 1 else -1.0
vmax = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
lam = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-5

prob = problem_from_classification_mat(
    REPO / "Classification" / "MAT_files" / "Classification_day_183_315_gurobi.mat", x0_from="t0")
M, n = prob.model, prob.model.n
q, r = np.ones(n), np.full(n, lam)

open_loop = solve(lambda t, X: M.rhs(t, X), (prob.t0, prob.tf), prob.x0).sol
tic = time.time()
res = direct_tracking(M, prob.z, prob.x0, prob.t0, prob.tf, q, r, (vmin, vmax), init_traj=open_loop)
print(f"IPOPT: {res['stats']['return_status']}  in {time.time() - tic:.1f}s, "
      f"discrete J = {res['J_discrete']:.5e}")

vfun = piecewise_control(res["t_nodes"], res["V"])
traj = solve(lambda t, X: M.rhs(t, X, vfun(t)), (prob.t0, prob.tf), prob.x0).sol
step = ASREStep(traj=traj, control=lambda t, X: vfun(t))
t, X, V, J = evaluate(step, M, prob.z, q, r, prob.t0, prob.tf)
print(f"re-simulated J (same evaluator as ASRE) = {J:.5e}")
sat = float(np.mean((res["V"] <= vmin + 1e-6) | (res["V"] >= vmax - 1e-6)))
print(f"fraction of control values at a bound = {sat:.2f}")

out = Path(__file__).parent / "results"
out.mkdir(exist_ok=True)
with open(out / f"direct_reactor_A_v[{vmin},{vmax}]_lam{lam:g}.pkl", "wb") as fh:
    pickle.dump(dict(res={k: v for k, v in res.items() if k != "stats"}, J=J, t=t, X=X, V=V), fh)
