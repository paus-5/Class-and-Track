"""Simulate the gLV proof-of-concept model (case study 1) and plot it in the layout of the
article's synthetic-data figure (tex/Images/proof_of_concept/250309_POC_*.png), for
D = 0.22 (synthetic_data.m, which made those images) and D = 0.24 (article text)."""
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from class_and_track.lv import case1_model, simulate_glv

FIG = Path(__file__).parent / "figures"
FIG.mkdir(exist_ok=True)
x0 = np.array([0.003, 0.0004, 0.1, 0.1, 0.1])   # read off the article plots (MATLAB: random)
t = np.linspace(0, 300, 3001)

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": "#e1e0d9",
                     "grid.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#c3c2b7", "figure.facecolor": "#fcfcfb",
                     "axes.facecolor": "#fcfcfb", "legend.frameon": False, "savefig.dpi": 150})
fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
S_COL = ["#e34948", "#008300", "#2a78d6"]         # s1, s2, s3 (red/green/blue as in the article)
for D, ls in ((0.22, "-"), (0.24, "--")):
    M, A = case1_model(D=D)
    Y = simulate_glv(M, A, x0, (0, 300), t_eval=t).y
    print(f"D = {D}: x1 {Y[0, 1490]:.3f} -> {Y[0, -1]:.3f}, x2 peak {Y[1].max():.4f} at day "
          f"{t[Y[1].argmax()]:.0f}, x2(300) = {Y[1, -1]:.4f}, s = {Y[2:, 1490].round(3)} -> {Y[2:, -1].round(3)}")
    axes[0].plot(t, Y[0], ls, color="#e34948", lw=1.6, label=f"D = {D}")
    axes[1].plot(t, Y[1], ls, color="#2a78d6", lw=1.6, label=f"D = {D}")
    for j in range(3):
        axes[2].plot(t, Y[2 + j], ls, color=S_COL[j], lw=1.4,
                     label=f"s{j + 1}" if D == 0.22 else None)
for ax, title in zip(axes, ["AOB (OTU 1)", "NOB (OTU 2)", "metabolites"]):
    ax.set_title(title, loc="left")
    ax.set_xlabel("time [days]")
    ax.axvline(150, color="#c3c2b7", lw=0.8)
axes[0].set_ylabel("concentration [g/L]")
axes[0].legend(title="solid / dashed")
axes[1].legend()
axes[2].legend(title="solid D=0.22, dashed 0.24", fontsize=8)
fig.suptitle("gLV case study 1: A = [[1, 0.1], [-1.5, 1]], s_in 1.25 -> 1.95 at day 150",
             x=0.01, ha="left", fontsize=10)
fig.tight_layout()
fig.savefig(FIG / "lv_dynamics_case1.png")
print(f"figure: {FIG / 'lv_dynamics_case1.png'}")
