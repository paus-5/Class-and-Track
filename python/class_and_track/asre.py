"""ASRE nonlinear optimal tracking (Cimen & Banks, Automatica 40 (2004) 1845-1863),
specialised to the nitrification model as in Tracking/main.m.

Thanks to the structure of A(X), B(X), C = [I 0] and diagonal Q, R, the (n+3)x(n+3)
Riccati matrix P and the feed-forward vector s_f reduce to n scalar ODEs each
(appendix "Tracking problem reformulation" of the article):

    p_j' = -q_j - 2 (f_j - D) p_j + p_j^2 f_j^2 x_j^2 / r_j ,      p_j(tf) = F_j
    s_j' = -q_j z_j - (f_j - D - f_j^2 x_j^2 p_j / r_j) s_j ,     s_j(tf) = F_j z_j(tf)
    v_j  = clip( -(f_j x_j^{prev} / r_j) (p_j x_j - s_j), v_min, v_max )

with f_j, x_j evaluated on the previous iterate X^{[i-1]} and x_j in the feedback on
the current state, exactly as in ricatti_diff.m, sF_diff.m and tracking_control.m.

Differences from the MATLAB code (deliberate; flagged with NOTE below):
  * q_j multiplies z_j in the s_f equation (the paper's -C^T Q z).  MATLAB omitted Q;
    identical results because MATLAB used Q = I.
  * Trajectories are dense ODE solutions (solve_ivp dense_output) instead of interp1
    closures over ode15s output.
  * Convergence is measured on a fixed time grid (MATLAB: Frobenius norm over the
    adaptive ode15s grid, whose length changes from one iteration to the next).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

from .model import NitrificationModel

TrajFun = Callable[[float], np.ndarray]


def solve(fun, t_span, y0, method="LSODA", rtol=1e-6, atol=1e-9):
    sol = solve_ivp(fun, t_span, y0, method=method, dense_output=True, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")
    return sol


@dataclass
class ASREStep:
    traj: TrajFun                      # X^{[i]}(t)
    P: TrajFun | None = None           # diag of P~^{[i]}(t)
    sf: TrajFun | None = None          # s~_f^{[i]}(t)
    control: Callable | None = None    # v(t, X)
    diff: float = np.nan


@dataclass
class ASREResult:
    iterations: list[ASREStep] = field(default_factory=list)
    converged: bool = False

    @property
    def last(self) -> ASREStep:
        return self.iterations[-1]


def asre_step(model: NitrificationModel, prev: TrajFun, z: TrajFun, x0, t0, tf,
              q, r, F, v_bounds=(-1.0, 0.5), **ode) -> ASREStep:
    """One ASRE iteration given the previous iterate X^{[i-1]} = prev(t)."""
    n = model.n
    vmin, vmax = v_bounds

    def riccati(t, p):
        X = prev(t)
        f, x = model.growth(X[n:]), X[:n]
        return -q - 2.0 * p * (f - model.D(t)) + p**2 * f**2 * x**2 / r

    P = solve(riccati, (tf, t0), F.copy(), **ode).sol       # backward in time

    def sf_rhs(t, s):
        X = prev(t)
        f, x = model.growth(X[n:]), X[:n]
        return -q * z(t) - (f - model.D(t) - f**2 * x**2 * P(t) / r) * s

    SF = solve(sf_rhs, (tf, t0), F * z(tf), **ode).sol

    def control(t, X):
        b1 = model.B1_diag(prev(t))                    # B on the previous iterate (paper eq. 9)
        v = -(b1 / r) * (P(t) * X[:n] - SF(t))
        return np.clip(v, vmin, vmax)                  # NOTE: ad-hoc saturation, not in the paper

    # NOTE: the paper integrates the LTV system A(X^{[i-1]})X + B(X^{[i-1]})u (eq. 12);
    # the MATLAB code (and this port) integrate the full nonlinear plant under the feedback.
    traj = solve(lambda t, X: model.rhs(t, X, control(t, X)), (t0, tf), x0, **ode).sol
    return ASREStep(traj=traj, P=P, sf=SF, control=control)


def asre_tracking(model: NitrificationModel, z: TrajFun, x0, t0: float, tf: float,
                  q, r, F=None, v_bounds=(-1.0, 0.5), tol=1e-3, max_iter=30,
                  initial: TrajFun | None = None, n_grid=2000, verbose=True,
                  method="LSODA", rtol=1e-6, atol=1e-9) -> ASREResult:
    """Run the ASRE iteration.

    q, r     : diagonals of Q and R (MATLAB main.m: q = 1, r = lambda).
    v_bounds : saturation of v = u - 1 (MATLAB tracking_control.m: [-1, 0.5]).
    tol      : stop when max_t max_j |x_j^[i] - x_j^[i-1]| / max_t |z_j| < tol.
    initial  : X^{[0]}; default = open-loop simulation with v = 0 (as in MATLAB).
               The paper uses X^{[-1]}(t) = x0 (Remark 2 allows a better guess).
    """
    n = model.n
    ode = dict(method=method, rtol=rtol, atol=atol)
    vec = lambda a: np.broadcast_to(np.asarray(a, float), (n,)).copy()
    q, r = vec(q), vec(r)
    F = np.zeros(n) if F is None else vec(F)
    x0 = np.asarray(x0, float)
    grid = np.linspace(t0, tf, n_grid)
    z_scale = np.max(np.abs(np.array([z(t) for t in grid])), axis=0)
    z_scale[z_scale == 0] = 1.0

    if initial is None:
        initial = solve(lambda t, X: model.rhs(t, X), (t0, tf), x0, **ode).sol
    res = ASREResult([ASREStep(traj=initial)])
    prev = initial
    for it in range(1, max_iter + 1):
        step = asre_step(model, prev, z, x0, t0, tf, q, r, F, v_bounds, **ode)
        step.diff = float(np.max(np.abs(step.traj(grid)[:n] - prev(grid)[:n]).T / z_scale))
        res.iterations.append(step)
        if verbose:
            print(f"ASRE iter {it}: max relative change in x = {step.diff:.3e}", flush=True)
        if step.diff < tol:
            res.converged = True
            break
        prev = step.traj
    return res


def evaluate(step: ASREStep, model: NitrificationModel, z: TrajFun, q, r, t0, tf, n_grid=4000):
    """Return (t, X, V, J) with J = int (x-z)^T Q (x-z) + v^T R v dt (trapezoid)."""
    n = model.n
    t = np.linspace(t0, tf, n_grid)
    X = step.traj(t)
    Z = np.array([z(tt) for tt in t]).T
    if step.control is None:
        V = np.zeros((n, t.size))
    else:
        V = np.array([step.control(tt, X[:, k]) for k, tt in enumerate(t)]).T
    e = X[:n] - Z
    integrand = (np.asarray(q)[:, None] * e**2).sum(0) + (np.asarray(r)[:, None] * V**2).sum(0)
    return t, X, V, float(np.trapezoid(integrand, t))
