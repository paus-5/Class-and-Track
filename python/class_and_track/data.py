"""Data loading (port of 'General Functions/') and tracking-problem setup (port of the
preamble of Tracking/main.m)."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from .model import NitrificationModel

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "Data"
REACTOR_VOLUME_L = 6.5   # load_biomass.m: dilution_rate = flow / 6.5


def _numeric_block(path: Path) -> np.ndarray:
    """Mimic MATLAB xlsread: keep the numeric block, drop text header rows/columns."""
    df = pd.read_excel(path, header=None).apply(pd.to_numeric, errors="coerce")
    df = df.dropna(how="all").dropna(axis=1, how="all")
    return df.to_numpy(float)


def load_biomass(reactor: str, absolute_time: bool = False):
    """absolute_time=False reproduces MATLAB (t - t(1)); note B.xls starts at day 5, so
    MATLAB's Reactor B times are shifted by 5 days. Use absolute_time=True to compare reactors."""
    d = _numeric_block(DATA / f"{reactor}.xls")
    t_obs = d[:, 0] - (0.0 if absolute_time else d[0, 0])
    return dict(t_obs=t_obs, biomass=d[:, 1], s_in=d[:, 2], S1=d[:, 3], S2=d[:, 4],
                S3=d[:, 5], dilution_rate=d[:, 6] / REACTOR_VOLUME_L)


def load_relative_abundance(reactor: str, absolute_time: bool = False):
    d = _numeric_block(DATA / f"especes{reactor}.xls")
    return d[:, 0] - (0.0 if absolute_time else d[0, 0]), d[:, 1:]


def classification_inputs(reactor: str, t_change=183.0, delta_t=132.0, absolute_time=True):
    """Port of the preamble of Classification/main.m: returns the MIP inputs
    W (n_considered, T), xG1, xG2 (T,), plus considered OTU indices and OTU_interp."""
    from scipy.integrate import solve_ivp
    b = load_biomass(reactor, absolute_time)
    t_otu, rel = load_relative_abundance(reactor, absolute_time)
    rel = rel / rel.sum(1, keepdims=True)
    t_obs = b["t_obs"]
    idx = np.flatnonzero((t_obs > t_change) & (t_obs < t_change + delta_t))
    idx_otu = np.flatnonzero((t_otu > t_change) & (t_otu < t_change + delta_t))
    considered = np.flatnonzero(rel[idx_otu].max(0) > 0)
    otu_interp = interpolate_biomass(t_obs, b["biomass"], t_otu, rel)
    D, sin = lin(t_obs, b["dilution_rate"]), prev(t_obs, b["s_in"])
    t0, tF = t_obs[idx[0]], t_obs[idx[-1]]
    # Invariant Z1' = -D (Z1 - s_in), Z1(t0) = s_in(t0)
    sol = solve_ivp(lambda t, z: -D(t) * (z - sin(t)), (t0, tF), [b["s_in"][idx[0]]],
                    method="LSODA", t_eval=t_obs[idx], rtol=1e-8, atol=1e-10, max_step=1.0)
    Z1 = sol.y[0]
    W = lin(t_otu, otu_interp[:, considered])(t_obs[idx]).T
    xG1 = Z1 - b["S1"][idx]
    xG2 = Z1 - b["S1"][idx] - b["S2"][idx]
    return dict(W=W, xG1=xG1, xG2=xG2, considered=considered, n_otu=rel.shape[1],
                t_obs=t_obs, t_otu=t_otu, otu_interp=otu_interp, time_steps_index=idx)


def interpolate_biomass(t_obs, biomass, t_otu, otu_rel):
    """Absolute OTU abundance = total biomass (interpolated at OTU times) x relative abundance."""
    return np.interp(t_otu, t_obs, biomass)[:, None] * otu_rel


def lin(t, y):
    """interp1(t, y, ., 'linear') returning arrays (column-wise for 2-D y)."""
    return interp1d(t, y, axis=0, kind="linear", bounds_error=False, fill_value=np.nan,
                    assume_sorted=True)


def prev(t, y):
    """interp1(t, y, ., 'previous') (zero-order hold) -- used for s_in."""
    return interp1d(t, y, axis=0, kind="previous", bounds_error=False,
                    fill_value=(y[0], y[-1]), assume_sorted=True)


@dataclass
class TrackingProblem:
    model: NitrificationModel
    z: callable          # z(t): measured OTU biomasses (length n)
    x0: np.ndarray
    t0: float
    tf: float
    t_otu: np.ndarray
    z_data: np.ndarray   # OTU data at t_otu (rows) for plotting
    t_obs: np.ndarray
    S_data: np.ndarray   # (len(t_obs), 3) s1,s2,s3 measurements


def build_tracking_problem(t_obs, s_in, S, dilution_rate, t_otu, otu_abs,
                           class_aob, class_nob, yields_aob, yields_nob,
                           t_start_index: int, t_end=315.0,
                           muA=0.77, muB=1.07, kSA=0.07, kSB=0.14982,
                           x0_from="t0") -> TrackingProblem:
    """Port of the set-up in Tracking/main.m / simulate_no_control.m.

    x0_from: "t0"   -> OTU data interpolated at the start time t0 (corrected)
             "row0" -> first OTU sample of the whole series (what MATLAB did:
                       x0 = biomass_filtered(1,:), i.e. the day-0 community).
    Default kinetics are those of simulate_no_control.m / the saved Reactor A runs;
    Tracking/main.m used muA=1.97, muB=1.87, kSA=0.7, kSB=0.54 instead.
    """
    idx_aob = np.flatnonzero(np.round(class_aob))
    idx_nob = np.flatnonzero(np.round(class_nob))
    nA, nB = idx_aob.size, idx_nob.size
    bf = otu_abs[:, np.r_[idx_aob, idx_nob]].copy()
    bf[bf == 0] = bf[bf > 0].min()           # MATLAB: zeros replaced by global min
    z = lin(t_otu, bf)

    last = np.flatnonzero(t_obs < t_end)[-1]
    span = np.arange(t_start_index, last + 1)
    tS = t_obs[span]
    t0, tf = float(tS[0]), float(tS[-1])
    S_lin = lin(tS, S[span])
    if x0_from == "row0":
        xb0 = bf[0]
    else:
        xb0 = z(t0)
    x0 = np.concatenate([xb0, S_lin(t0)])

    D, sin = lin(t_obs, dilution_rate), prev(t_obs, s_in)
    model = NitrificationModel(
        nA=nA, nB=nB,
        kA=np.asarray(yields_aob)[idx_aob], kB=np.asarray(yields_nob)[idx_nob],
        muA=np.full(nA, muA), muB=np.full(nB, muB),
        kSA=np.full(nA, kSA), kSB=np.full(nB, kSB),
        D=lambda t: float(D(t)), s_in=lambda t: float(sin(t)))
    return TrackingProblem(model, z, x0, t0, tf, t_otu, bf, t_obs, S)


def problem_from_classification_mat(path, classification=None, **kw) -> TrackingProblem:
    """Build the problem directly from a MATLAB Classification_*.mat workspace.

    classification: optional mapping with class_AOB, class_NOB, yields_AOB, yields_NOB
    (e.g. another reactor's results/classification_*.npz) replacing the saved ones.
    """
    from scipy.io import loadmat
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = loadmat(path, squeeze_me=True)
    if classification is not None:
        for k in ("class_AOB", "class_NOB", "yields_AOB", "yields_NOB"):
            c[k] = np.asarray(classification[k], float)
    S = np.column_stack([c["S1"], c["S2"], c["S3"]])
    return build_tracking_problem(
        c["t_obs"].astype(float), c["s_in"].astype(float), S, c["dilution_rate"].astype(float),
        c["t_OTU"].astype(float), c["OTU_interp"].astype(float),
        c["class_AOB"], c["class_NOB"], c["yields_AOB"], c["yields_NOB"],
        t_start_index=int(c["time_steps_index"][0]) - 1, **kw)   # MATLAB 1-based
