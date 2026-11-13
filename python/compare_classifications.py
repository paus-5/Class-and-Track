"""Compare reactor A's saved classification with reactor B's (classify_reactor.py), and how
much of reactor A's community (days 183-315) each one assigns to AOB / NOB."""
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.io import loadmat

sys.path.insert(0, str(Path(__file__).parent))
from class_and_track.data import REPO

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    c = loadmat(REPO / "Classification" / "MAT_files" / "Classification_day_183_315_gurobi.mat",
                squeeze_me=True)
B = np.load(Path(__file__).parent / "results" / "classification_B_l2.npz")


def label(a, b):
    return np.where(np.round(a) == 1, "AOB", np.where(np.round(b) == 1, "NOB", "-"))


LA = label(c["class_AOB"], c["class_NOB"])
LB = label(B["class_AOB"], B["class_NOB"])
t = c["t_OTU"].astype(float)
win = (t > 183) & (t < 315)
X = c["OTU_interp"][win]                       # reactor A absolute OTU biomass
mean_share = X.mean(0) / X.sum(1).mean()       # mean share of reactor A biomass

print(f"{'OTU':>4} {'A':>4} {'B':>4}  share of A biomass")
for i in np.argsort(-mean_share):
    if LA[i] != "-" or LB[i] != "-":
        flag = "" if LA[i] == LB[i] else "  <- differs"
        print(f"E{i + 1:<3} {LA[i]:>4} {LB[i]:>4}  {100 * mean_share[i]:6.2f}%{flag}")
agree = np.sum((LA == LB) & ((LA != "-") | (LB != "-")))
print(f"\nagree on {agree} of {np.sum((LA != '-') | (LB != '-'))} OTUs classified by either")
for name, L in (("A", LA), ("B", LB)):
    print(f"classification {name}: AOB cover {100 * mean_share[L == 'AOB'].sum():.1f}%, "
          f"NOB cover {100 * mean_share[L == 'NOB'].sum():.1f}% of reactor A biomass")
both = (LA == "AOB") & (LB == "NOB") | (LA == "NOB") & (LB == "AOB")
print("swapped AOB<->NOB:", " ".join(f"E{i + 1}" for i in np.flatnonzero(both)) or "none")
