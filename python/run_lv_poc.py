"""Proof of concept of the article (Article_Processes.tex, section "Proof of concept"):

 1. simulate the gLV nitrification model, case study 1 (1 AOB + 1 NOB), 300 days,
    s_in 1.25 -> 1.95 at day 150, D = 0.24; sample the biomasses daily (synthetic data),
 2. track the biomass data with the control model (u = 1 + v replaces 1 + A x), direct method,
    Q = I, R = diag(1e-4, 1e-5) as in the article,
 3. identify the interaction matrix A:
      a) article's two-step approach: least squares  u(t) - 1 = A x(t)
      b) direct parameter estimation of A in the ODE (collocation + IPOPT).

usage: python run_lv_poc.py [noise_rel] [seed]
    noise_rel: multiplicative log-normal noise on the biomass samples (default 0 = none)
"""
import sys
import time
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from class_and_track.data import lin
from class_and_track.direct import direct_tracking
from class_and_track.lv import case1_model, estimate_interactions, regress_interactions, simulate_glv

noise = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
tag = f"noise{noise:g}_seed{seed}"
rng = np.random.default_rng(seed)
HERE = Path(__file__).parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
np.set_printoptions(precision=4, suppress=True)

# ---- 1. synthetic data -----------------------------------------------------------------
M, A_true = case1_model()
n = M.n
t0, tf = 0.0, 300.0
x0 = np.array([0.003, 0.0004, 0.1, 0.1, 0.1])   # read off the article plots (MATLAB: random)
t_data = np.arange(t0, tf + 1)
sim = simulate_glv(M, A_true, x0, (t0, tf), t_eval=t_data)
Ytrue = sim.y
Z = Ytrue[:n] * np.exp(noise * rng.standard_normal((n, t_data.size))) if noise else Ytrue[:n].copy()
Z[:, 0] = x0[:n]                                              # x0 treated as known
print(f"true A =\n{A_true}")
print(f"x0 = {x0[:n]},  x at day 149: {Ytrue[:n, 149]},  x at day 300: {Ytrue[:n, -1]}")
print(f"s at day 149: {Ytrue[n:, 149]},  s at day 300: {Ytrue[n:, -1]}")

# ---- 2. tracking (direct method) -------------------------------------------------------
zfun = lin(t_data, Z.T)
z = lambda t: zfun(t)
q, r = np.ones(n), np.array([1e-4, 1e-5])
ol = simulate_glv(M, np.zeros((n, n)), x0, (t0, tf))           # no-interaction model
tic = time.time()
tr = direct_tracking(M, z, x0, t0, tf, q, r, v_bounds=(-1.0, 1.0), init_traj=ol.sol,
                     ipopt_opts={"ipopt.print_level": 0, "ipopt.sb": "yes"})
print(f"\ntracking: {tr['stats']['return_status']} in {time.time() - tic:.1f}s, J = {tr['J_discrete']:.4e}")
tk, Xtr, V = tr["t_nodes"], tr["X"], tr["V"]
tm = 0.5 * (tk[:-1] + tk[1:])                                 # interval midpoints
U = 1 + V
x_mid_tr = 0.5 * (Xtr[:n, :-1] + Xtr[:n, 1:])                 # tracked x on each interval
x_mid_true = sim.sol(tm)[:n]
U_true = 1 + A_true @ x_mid_true                              # the interaction to recover
rel_err = np.abs(Xtr[:n] - Ytrue[:n]) / np.abs(Ytrue[:n]).max(1, keepdims=True)
print(f"tracking max rel. error x: {rel_err.max(1)}")
print(f"control vs true 1 + A x (days 5-280): rms = "
      f"{np.sqrt(np.mean((U - U_true)[:, (tm > 5) & (tm < 280)] ** 2, 1))}")

# ---- 3a. two-step identification (article) ---------------------------------------------
win = (5.0, 280.0)                                             # drop start and terminal relaxation
A_reg, res_reg = regress_interactions(tm, U, x_mid_tr, win)
A_reg_true_u, _ = regress_interactions(tm, U_true, x_mid_true, win)   # sanity: exact u
print(f"\n(a) two-step, regress u - 1 = A x on the tracked trajectory, days {win}:\n{A_reg}")
print(f"    residual rms per row {res_reg};  rel. error |A_hat - A| / |A| = "
      f"{np.linalg.norm(A_reg - A_true) / np.linalg.norm(A_true):.3f}")
print(f"    (with the exact u = 1 + A x the regression returns\n{A_reg_true_u})")
Phi = x_mid_tr[:, (tm >= win[0]) & (tm <= win[1])].T
sv = np.linalg.svd(Phi / np.abs(Phi).max(0), compute_uv=False)
print(f"    singular values of the scaled regressors x(t): {sv}  (cond = {sv[0] / sv[-1]:.1f})")

# ---- figures ---------------------------------------------------------------------------
C_TRUE, C_TR, C_REG, C_EST, C_DATA, C_OL = "#0b0b0b", "#2a78d6", "#eb6834", "#1baf7a", "#52514e", "#898781"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": "#e1e0d9",
                     "grid.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#c3c2b7", "figure.facecolor": "#fcfcfb",
                     "axes.facecolor": "#fcfcfb", "legend.frameon": False, "savefig.dpi": 150})
names = ["x1 (AOB)", "x2 (NOB)", "s1 (NH4)", "s2 (NO2)", "s3 (NO3)"]
sim_reg = simulate_glv(M, A_reg, x0, (t0, tf), t_eval=t_data)
fig, axes = plt.subplots(2, 3, figsize=(12, 6), sharex=True)
for i, ax in enumerate(axes.flat[:5]):
    if i < n:
        ax.plot(t_data[::3], Z[i, ::3], "o", ms=2.5, color=C_DATA, label="synthetic data", zorder=5)
    else:
        ax.plot(t_data, Ytrue[i], color=C_TRUE, lw=1, ls=":", label="true (not used)")
    ax.plot(ol.t, ol.y[i], color=C_OL, lw=1, label="no interactions (A = 0)")
    ax.plot(tk, Xtr[i], color=C_TR, lw=1.6, label="tracking (direct)")
    ax.plot(t_data, sim_reg.y[i], color=C_REG, lw=1.2, ls="--", label="gLV with A from (a)")
    ax.set_title(names[i], loc="left")
    ax.axvline(150, color="#c3c2b7", lw=0.8)
ax = axes.flat[5]
for i, c in enumerate((C_TR, C_REG)):
    ax.plot(tm, U[i], color=c, lw=1.4, label=f"control u{i + 1}")
    ax.plot(tm, U_true[i], color=C_TRUE, lw=1, ls="--" if i else ":",
            label=f"true 1 + (A x)_{i + 1}")
ax.set_title("control vs true interaction term", loc="left")
for a in axes[1]:
    a.set_xlabel("day")
h1, l1 = axes.flat[0].get_legend_handles_labels()
h2, l2 = axes.flat[2].get_legend_handles_labels()
fig.legend(h1 + [h2[0]], l1 + [l2[0]], loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.0))
axes.flat[5].legend(fontsize=7, loc="lower left")
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(FIG / f"lv_poc_{tag}.png")

fig, axes = plt.subplots(1, n, figsize=(9, 3.2))
for i, ax in enumerate(axes):
    m = (tm >= win[0]) & (tm <= win[1])
    ax.plot(tm[m], U[i, m] - 1, color=C_TR, lw=1.2, label="u - 1 (tracking control)")
    ax.plot(tm[m], (A_true @ x_mid_true)[i, m], color=C_TRUE, lw=1, ls="--", label="true (A x)")
    ax.plot(tm[m], (A_reg @ x_mid_tr)[i, m], color=C_REG, lw=1.2, label="fitted (A_hat x)")
    ax.set_title(f"row {i + 1}: a = {A_reg[i].round(3)}  (true {A_true[i]})", loc="left", fontsize=8)
    ax.set_xlabel("day")
axes[0].legend(fontsize=7)
fig.tight_layout()
fig.savefig(FIG / f"lv_poc_regression_{tag}.png")
print(f"\nfigures: {FIG / f'lv_poc_{tag}.png'}, {FIG / f'lv_poc_regression_{tag}.png'}", flush=True)

# ---- 3b. direct estimation of A in the ODE (benchmark) ---------------------------------
print("\n(b) direct estimation of A in the gLV ODE from the biomass data (several starts):")
best = None
for A_init in (np.zeros((n, n)), A_reg, np.full((n, n), 1.0), rng.normal(0, 1, (n, n))):
    tic = time.time()
    try:
        e = estimate_interactions(M, t_data, Z, x0, t0, tf, A_init=A_init)
        print(f"    from A0 = {A_init.ravel().round(2)} -> A = {e['A'].ravel().round(4)}  "
              f"J = {e['J']:.3e}  ({e['status']}, {time.time() - tic:.0f}s)", flush=True)
        if best is None or e["J"] < best["J"]:
            best = e
    except RuntimeError as err:
        print(f"    from A0 = {A_init.ravel().round(2)} -> failed: {str(err).split(chr(39))[-2]} "
              f"({time.time() - tic:.0f}s)", flush=True)
if best is not None:
    print(f"    best: rel. error |A_hat - A| / |A| = "
          f"{np.linalg.norm(best['A'] - A_true) / np.linalg.norm(A_true):.3f}")

# Profile of the fit cost in a12: is a12 constrained by the data at all?
try:
    e = estimate_interactions(M, t_data, Z, x0, t0, tf, A_init=A_true,
                              fixed={(i, j): A_true[i, j] for i in range(n) for j in range(n)})
    print(f"\ncost at the true A (all entries fixed): J = {e['J']:.3e}  "
          f"(> 0 only through discretisation / noise)")
except RuntimeError as err:
    print(f"\ncost at the true A: failed: {str(err).split(chr(39))[-2]}")
print("\nprofile in a12 (other entries free, warm start from the true A):")
for a12 in (-0.5, -0.2, 0.0, 0.1, 0.2, 0.5):
    try:
        e = estimate_interactions(M, t_data, Z, x0, t0, tf, A_init=A_true, fixed={(0, 1): a12})
        print(f"    a12 = {a12:5.2f}: J = {e['J']:.3e}  A = {e['A'].ravel().round(4)}", flush=True)
    except RuntimeError as err:
        print(f"    a12 = {a12:5.2f}: failed: {str(err).split(chr(39))[-2]}", flush=True)
