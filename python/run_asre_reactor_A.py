"""Run the ASRE tracking on Reactor A (days 184-313) and log convergence and cost.

usage: python run_asre_reactor_A.py [row0|t0] [max_iter]
"""
import pickle
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from class_and_track.asre import asre_step, evaluate, solve
from class_and_track.data import REPO, problem_from_classification_mat

x0_from = sys.argv[1] if len(sys.argv) > 1 else "t0"
max_iter = int(sys.argv[2]) if len(sys.argv) > 2 else 25
vmin = float(sys.argv[3]) if len(sys.argv) > 3 else -1.0     # use -inf / inf to remove bounds
vmax = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
lam = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-5
tag = f"{x0_from}_v[{vmin},{vmax}]_lam{lam:g}"

prob = problem_from_classification_mat(
    REPO / "Classification" / "MAT_files" / "Classification_day_183_315_gurobi.mat", x0_from=x0_from)
M, n = prob.model, prob.model.n
q, r, F = np.ones(n), np.full(n, lam), np.zeros(n)   # lam = 1e-5 in the saved Reactor A runs
grid = np.linspace(prob.t0, prob.tf, 2000)
zs = np.abs(prob.z(grid)).max(0)

prev = solve(lambda t, X: M.rhs(t, X), (prob.t0, prob.tf), prob.x0).sol
from class_and_track.asre import ASREStep
_, _, _, J0 = evaluate(ASREStep(traj=prev), M, prob.z, q, r, prob.t0, prob.tf)
print(f"[{tag}] iter 0 (no control): J = {J0:.5e}", flush=True)
log = [dict(it=0, J=J0, diff=np.nan)]
best = None
for it in range(1, max_iter + 1):
    tic = time.time()
    step = asre_step(M, prev, prob.z, prob.x0, prob.t0, prob.tf, q, r, F, v_bounds=(vmin, vmax))
    step.diff = float(np.max(np.abs(step.traj(grid)[:n] - prev(grid)[:n]).T / zs))
    t, X, V, J = evaluate(step, M, prob.z, q, r, prob.t0, prob.tf)
    sat = float(np.mean((V <= vmin + 1e-9) | (V >= vmax - 1e-9)))
    print(f"[{tag}] iter {it}: J = {J:.5e}  max rel change = {step.diff:.3e}  "
          f"fraction of time-OTU points saturated = {sat:.2f}  ({time.time() - tic:.1f}s)", flush=True)
    log.append(dict(it=it, J=J, diff=step.diff, sat=sat))
    if best is None or J < best[1]:
        best = (it, J, t, X, V)
    if step.diff < 1e-3:
        print("converged")
        break
    prev = step.traj

out = Path(__file__).parent / "results"
out.mkdir(exist_ok=True)
with open(out / f"asre_reactor_A_{tag}.pkl", "wb") as fh:
    pickle.dump(dict(log=log, best=best, x0=prob.x0, t0=prob.t0, tf=prob.tf), fh)
print(f"best iterate {best[0]}  J = {best[1]:.5e}")
