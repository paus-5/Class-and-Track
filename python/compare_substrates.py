"""Substrate fit on Reactor A for direct-method runs with different OTU classifications.

The tracked OTU sets differ between classifications, so J is not comparable; the substrates
(s1 = NH4, s2 = NO2, s3 = NO3) are measured independently of any classification and are
not in the cost, so how well they follow the data is a fair, out-of-objective check.

usage: python compare_substrates.py [vmin] [vmax] [lambda]
"""
import pickle
import sys
import warnings
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import loadmat

sys.path.insert(0, str(Path(__file__).parent))
from class_and_track.asre import solve
from class_and_track.data import REPO, problem_from_classification_mat

args = sys.argv[1:]
vmin = float(args[0]) if len(args) > 0 else -1.0
vmax = float(args[1]) if len(args) > 1 else 0.5
lam = float(args[2]) if len(args) > 2 else 1e-5
tag = f"v[{vmin},{vmax}]_lam{lam:g}"
HERE = Path(__file__).parent
RES, FIG = HERE / "results", HERE / "figures"
FIG.mkdir(exist_ok=True)
MAT = REPO / "Classification" / "MAT_files" / "Classification_day_183_315_gurobi.mat"

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    class_A = loadmat(MAT, squeeze_me=True)

# (label, classification, direct-result file, colour) -- validated categorical slots
cases = [("classification A", class_A, f"direct_reactor_A_{tag}.pkl", "#2a78d6"),
         ("classification B", dict(np.load(RES / "classification_B_l2.npz")),
          f"direct_reactor_A_{tag}_classB.pkl", "#eb6834")]

prob = problem_from_classification_mat(MAT, x0_from="t0")
win = (prob.t_obs >= prob.t0) & (prob.t_obs <= prob.tf)
t_d, S_d = prob.t_obs[win], prob.S_data[win]
names = ["s1 (NH4)", "s2 (NO2)", "s3 (NO3)"]

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": "#e1e0d9",
                     "grid.linewidth": 0.5, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": "#c3c2b7",
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
                     "legend.frameon": False, "savefig.dpi": 150})
fig, axes = plt.subplots(2, 3, figsize=(12, 5.6), sharex=True, sharey="col")
print(f"RMSE vs measured substrates on Reactor A, days {prob.t0:g}-{prob.tf:g}  ({tag})")
print(f"{'':28s}" + "".join(f"{nm:>12s}" for nm in names))
for row, (name, cl, fname, color) in enumerate(cases):
    p = problem_from_classification_mat(MAT, classification=cl, x0_from="t0")
    n = p.model.n
    ol = solve(lambda t, X: p.model.rhs(t, X), (p.t0, p.tf), p.x0).sol
    S_open = ol(t_d)[n:].T
    rm_open = np.sqrt(np.mean((S_open - S_d) ** 2, 0))
    print(f"{name + ', open loop':28s}" + "".join(f"{v:12.4f}" for v in rm_open))
    fp = RES / fname
    if fp.exists():
        d = pickle.load(open(fp, "rb"))
        S_dir = np.array([np.interp(t_d, d["t"], d["X"][n + j]) for j in range(3)]).T
        rm_dir = np.sqrt(np.mean((S_dir - S_d) ** 2, 0))
        print(f"{name + ', direct':28s}" + "".join(f"{v:12.4f}" for v in rm_dir)
              + f"     (J = {d['J']:.4g}, {p.model.nA} AOB / {p.model.nB} NOB)")
    else:
        d = None
        print(f"{name + ', direct':28s}  missing {fname}")
    tt = np.linspace(p.t0, p.tf, 1500)
    for j, ax in enumerate(axes[row]):
        ax.plot(t_d, S_d[:, j], "o", color="#0b0b0b", ms=3, zorder=5, label="measured")
        ax.plot(tt, ol(tt)[n + j], color="#898781", lw=1, label="open loop")
        if d is not None:
            ax.plot(d["t"], d["X"][n + j], color=color, lw=1.6, label="direct")
            ax.text(0.02, 0.95, f"RMSE {rm_dir[j]:.3f}", transform=ax.transAxes, va="top",
                    fontsize=8, color="#52514e")
        ax.set_title(f"{name}: {names[j]}", loc="left")
    axes[row, 0].set_ylabel("concentration")
for ax in axes[1]:
    ax.set_xlabel("day")
handles, labs = [], []
for ax in axes[:, 0]:
    for h, l in zip(*ax.get_legend_handles_labels()):
        if l == "direct":
            l = f"direct ({ax.get_title(loc='left').split(':')[0]})"
        if l not in labs:
            handles.append(h); labs.append(l)
fig.legend(handles, labs, loc="upper center", ncol=len(labs), bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(FIG / f"substrates_classA_vs_classB_{tag}.png")
print(f"figure: {FIG / f'substrates_classA_vs_classB_{tag}.png'}")
