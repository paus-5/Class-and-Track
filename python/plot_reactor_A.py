"""Plot Reactor A results: open loop vs ASRE (best iterate) vs direct method.

Reads the pickles written by run_asre_reactor_A.py / run_direct_reactor_A.py and writes
PNGs to python/figures/:
    cost_<tag>.png        J per ASRE iteration (all ASRE runs found) + direct and open loop
    otus_<tag>.png        per-OTU biomass trajectories vs data (small multiples)
    control_<tag>.png     v_i(t) heatmaps, ASRE vs direct, same colour scale
    substrates_<tag>.png  s1, s2, s3 vs measurements

usage: python plot_reactor_A.py [vmin] [vmax] [lambda] [--show]
"""
import pickle
import sys
from pathlib import Path

import matplotlib
import numpy as np

if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, str(Path(__file__).parent))
from class_and_track.asre import solve
from class_and_track.data import REPO, problem_from_classification_mat

args = [a for a in sys.argv[1:] if not a.startswith("--")]
vmin = float(args[0]) if len(args) > 0 else -1.0
vmax = float(args[1]) if len(args) > 1 else 0.5
lam = float(args[2]) if len(args) > 2 else 1e-5
tag = f"v[{vmin},{vmax}]_lam{lam:g}"

HERE = Path(__file__).parent
RES, FIG = HERE / "results", HERE / "figures"
FIG.mkdir(exist_ok=True)

# Palette: fixed roles (validated categorical slots 1-3; reference/open loop in neutral ink).
C_DIRECT, C_ASRE, C_ASRE2 = "#2a78d6", "#eb6834", "#1baf7a"
C_REF, C_OPEN, C_MUTED, C_GRID = "#0b0b0b", "#898781", "#52514e", "#e1e0d9"
plt.rcParams.update({
    "font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"], "font.size": 9,
    "axes.edgecolor": "#c3c2b7", "axes.labelcolor": C_MUTED, "axes.titlesize": 9,
    "axes.titlecolor": C_REF, "xtick.color": C_MUTED, "ytick.color": C_MUTED,
    "axes.grid": True, "grid.color": C_GRID, "grid.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.dpi": 150,
    "legend.frameon": False,
})


def load(name):
    p = RES / name
    if not p.exists():
        return None
    with open(p, "rb") as fh:
        return pickle.load(fh)


prob = problem_from_classification_mat(
    REPO / "Classification" / "MAT_files" / "Classification_day_183_315_gurobi.mat", x0_from="t0")
M, n, nA = prob.model, prob.model.n, prob.model.nA
labels = [f"AOB {i + 1}" for i in range(nA)] + [f"NOB {i + 1}" for i in range(n - nA)]

asre = load(f"asre_reactor_A_t0_{tag}.pkl")
direct = load(f"direct_reactor_A_{tag}.pkl")
if asre is None and direct is None:
    sys.exit(f"no results for {tag} in {RES}")

t = np.linspace(prob.t0, prob.tf, 2000)
X_open = solve(lambda tt, X: M.rhs(tt, X), (prob.t0, prob.tf), prob.x0).sol(t)
runs = [("open loop", t, X_open, None, C_OPEN)]
if asre is not None:
    it, J, ta, Xa, Va = asre["best"]
    runs.append((f"ASRE (iter {it}, J={J:.3g})", ta, Xa, Va, C_ASRE))
if direct is not None:
    runs.append((f"direct (J={direct['J']:.3g})", direct["t"], direct["X"], direct["V"], C_DIRECT))
in_win = (prob.t_otu >= prob.t0) & (prob.t_otu <= prob.tf)
t_data, z_data = prob.t_otu[in_win], prob.z_data[in_win]

# ---- 1. cost per ASRE iteration -------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 3.6))
asre_runs = sorted(RES.glob(f"asre_reactor_A_t0_*_lam{lam:g}.pkl"))
colors = {f"v[{vmin},{vmax}]": C_ASRE, "v[-1.0,inf]": C_ASRE2, "v[-inf,inf]": C_MUTED}
J_open = None
for p in asre_runs:
    d = load(p.name)
    bounds = p.stem.split("_t0_")[1].split("_lam")[0]
    its = [e["it"] for e in d["log"]]
    Js = [e["J"] for e in d["log"]]
    J_open = Js[0]
    c = colors.get(bounds, C_MUTED)
    ax.plot(its, Js, "-o", color=c, lw=2, ms=4, mec="#fcfcfb", mew=1, label=f"ASRE {bounds}")
if direct is not None:
    ax.axhline(direct["J"], color=C_DIRECT, lw=2, ls="--", label=f"direct v[{vmin},{vmax}]")
if J_open is not None:
    ax.axhline(J_open, color=C_OPEN, lw=1, ls=":", label="open loop (no control)")
ax.set_yscale("log")
ax.set_xlabel("ASRE iteration")
ax.set_ylabel("tracking cost J")
ax.set_title(f"Reactor A — cost, λ = {lam:g}", loc="left")
ax.legend(fontsize=8, loc="upper right")
fig.tight_layout()
fig.savefig(FIG / f"cost_{tag}.png")

# ---- 2. OTU trajectories (small multiples) --------------------------------------------
ncol = 6
nrow = int(np.ceil(n / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(13, 1.9 * nrow), sharex=True)
for i, ax in enumerate(axes.flat):
    if i >= n:
        ax.axis("off")
        continue
    ax.plot(t_data, z_data[:, i], "o", color=C_REF, ms=3, zorder=5)
    for name, tt, X, _, c in runs:
        ax.plot(tt, X[i], color=c, lw=1.0 if name == "open loop" else 1.5)
    ax.set_title(labels[i], loc="left", pad=2)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    ax.yaxis.get_offset_text().set_fontsize(7)
    ax.tick_params(labelsize=7)
handles = [plt.Line2D([], [], color=C_REF, marker="o", ls="", ms=4)] + \
          [plt.Line2D([], [], color=c, lw=2) for *_, c in runs]
fig.legend(handles, ["data"] + [r[0] for r in runs], loc="upper center", ncol=len(runs) + 1,
           bbox_to_anchor=(0.5, 1.0))
fig.supxlabel("day", color=C_MUTED)
fig.supylabel("biomass (g/L)", color=C_MUTED)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(FIG / f"otus_{tag}.png")

# ---- 3. control heatmaps --------------------------------------------------------------
ctl = [(name, tt, V) for name, tt, _, V, _ in runs if V is not None]
if ctl:
    lo = min(min(V.min() for *_, V in ctl), -1e-6)
    hi = max(max(V.max() for *_, V in ctl), 1e-6)
    norm = TwoSlopeNorm(vmin=lo, vcenter=0.0, vmax=hi)
    fig, axes = plt.subplots(len(ctl), 1, figsize=(10, 3.2 * len(ctl)), sharex=True, squeeze=False)
    for ax, (name, tt, V) in zip(axes[:, 0], ctl):
        im = ax.pcolormesh(tt, np.arange(n + 1) - 0.5, V[:, :-1] if V.shape[1] == tt.size else V,
                           cmap="RdBu_r", norm=norm, shading="flat")
        ax.set_yticks(range(n), labels, fontsize=6)
        ax.invert_yaxis()
        ax.axhline(nA - 0.5, color=C_REF, lw=0.8)
        ax.grid(False)
        ax.set_title(f"{name}: control v = u − 1   (bounds [{vmin:g}, {vmax:g}])", loc="left")
    axes[-1, 0].set_xlabel("day")
    cb = fig.colorbar(im, ax=axes[:, 0], shrink=0.8, pad=0.01)
    cb.set_label("v  (blue: inhibition, red: promotion)", color=C_MUTED)
    fig.savefig(FIG / f"control_{tag}.png", bbox_inches="tight")

# ---- 4. substrates --------------------------------------------------------------------
in_obs = (prob.t_obs >= prob.t0) & (prob.t_obs <= prob.tf)
fig, axes = plt.subplots(1, 3, figsize=(12, 3.2), sharex=True)
for j, (ax, nm) in enumerate(zip(axes, ["s1  (NH4)", "s2  (NO2)", "s3  (NO3)"])):
    ax.plot(prob.t_obs[in_obs], prob.S_data[in_obs, j], "o", color=C_REF, ms=3, zorder=5)
    for name, tt, X, _, c in runs:
        ax.plot(tt, X[n + j], color=c, lw=1.0 if name == "open loop" else 1.5)
    ax.set_title(nm, loc="left")
    ax.set_xlabel("day")
axes[0].set_ylabel("concentration")
fig.legend(handles, ["data"] + [r[0] for r in runs], loc="upper center", ncol=len(runs) + 1,
           bbox_to_anchor=(0.5, 1.02))
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(FIG / f"substrates_{tag}.png")

print(f"figures written to {FIG}")
if "--show" in sys.argv:
    plt.show()
