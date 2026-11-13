"""Direct (first-discretise-then-optimise) solution of the same tracking problem.

    min  int_{t0}^{tf} (x - z)^T Q (x - z) + v^T R v dt
    s.t. X' = rhs(t, X, v),  X(t0) = x0,  v_min <= v <= v_max,  X >= 0

Direct collocation (Radau IIA, piecewise-constant v on a uniform grid), solved with IPOPT
through CasADi.  Unlike ASRE + clipping, the control bounds are part of the optimisation,
so the result is a (local) optimum of the constrained problem.

Why collocation and not explicit multiple shooting: the substrate equations are stiff.
d(s1')/ds1 ~ -kA muA KA/(s1+KA)^2 sum x_A u  is O(-30)/day and the s2 analogue O(-250)/day,
so explicit RK4 with h = 1/4 day (h|lambda| >> 2.8) is unstable; that made IPOPT fail
with "Restoration Failed" and dual infeasibility ~1e20.  Radau IIA is L-stable, and since
every collocation state is a decision variable, positivity can be imposed everywhere
(which also keeps Monod away from its pole at s = -K).
"""
from __future__ import annotations

import casadi as ca
import numpy as np

from .model import NitrificationModel


def _rhs_casadi(M: NitrificationModel, X, v, D, s_in):
    n, nA = M.n, M.nA
    x, s1, s2, s3 = X[:n], X[n], X[n + 1], X[n + 2]
    fA = ca.DM(M.muA) * s1 / (s1 + ca.DM(M.kSA))
    fB = ca.DM(M.muB) * s2 / (s2 + ca.DM(M.kSB))
    r = ca.vertcat(fA, fB) * (1 + v) * x
    rA, rB = r[:nA], r[nA:]
    dx = r - D * x
    ds1 = (s_in - s1) * D - ca.dot(ca.DM(M.kA), rA)
    ds2 = -s2 * D + ca.dot(ca.DM(M.kA), rA) - ca.dot(ca.DM(M.kB), rB)
    ds3 = -s3 * D + ca.dot(ca.DM(M.kB), rB)
    return ca.vertcat(dx, ds1, ds2, ds3)


def _radau(d: int):
    """Collocation matrices: C[r, j] = l_r'(tau_j), quadrature weights B[j] (j = 1..d)."""
    tau = np.r_[0.0, ca.collocation_points(d, "radau")]
    C, B = np.zeros((d + 1, d + 1)), np.zeros(d + 1)
    for r in range(d + 1):
        p = np.poly1d([1.0])
        for m in range(d + 1):
            if m != r:
                p *= np.poly1d([1.0, -tau[m]]) / (tau[r] - tau[m])
        C[r] = np.polyder(p)(tau)
        B[r] = np.polyint(p)(1.0)
    return tau, C, B


def direct_tracking(M: NitrificationModel, z, x0, t0, tf, q, r, v_bounds=(-1.0, 0.5),
                    n_intervals: int | None = None, degree: int = 3,
                    init_traj=None, init_control=None, ipopt_opts: dict | None = None):
    """Returns dict(t_nodes, X (nx, N+1), V (n, N), J_discrete, stats).

    init_traj(t) -> (nx, len(t)) initial state guess (e.g. open loop or an ASRE iterate);
    init_control (n, N) initial control guess.
    """
    n, nx, d = M.n, M.n + 3, degree
    N = n_intervals or int(round(tf - t0))           # default: one interval per day
    tk = np.linspace(t0, tf, N + 1)
    h = (tf - t0) / N
    q, r = np.broadcast_to(q, (n,)).astype(float), np.broadcast_to(r, (n,)).astype(float)
    tau, C, B = _radau(d)
    tc = tk[:-1, None] + h * tau[None, 1:]           # (N, d) collocation times

    # Initial guess and per-state scaling (OTU biomasses span 1e-5 .. 1e-1 g/L).
    x0 = np.asarray(x0, float)
    if init_traj is not None:
        ref_n, ref_c = init_traj(tk), init_traj(tc.ravel()).reshape(nx, N, d)
    else:
        ref_n, ref_c = np.repeat(x0[:, None], N + 1, 1), np.repeat(x0[:, None, None], N, 1).repeat(d, 2)
    ref_n, ref_c = np.maximum(ref_n, 0), np.maximum(ref_c, 0)
    Zc = np.array([z(t) for t in tc.ravel()]).T.reshape(n, N, d)
    scale = np.maximum(np.r_[np.maximum(np.abs(Zc).max((1, 2)), x0[:n]),
                             ref_n[n:].max(1)], 1e-6)

    # One collocation interval as a CasADi function (in scaled variables), mapped over N.
    Xa, Xcs = ca.SX.sym("Xa", nx), ca.SX.sym("Xc", nx, d)
    vs, Ds, Ss, Zs = ca.SX.sym("v", n), ca.SX.sym("D", d), ca.SX.sym("sin", d), ca.SX.sym("z", n, d)
    sc = ca.DM(scale)
    Xall = [Xa * sc] + [Xcs[:, j] * sc for j in range(d)]
    defects, L = [], 0
    for j in range(1, d + 1):
        xp = sum(C[rr, j] * Xall[rr] for rr in range(d + 1))
        fj = _rhs_casadi(M, Xall[j], vs, Ds[j - 1], Ss[j - 1])
        defects.append((h * fj - xp) / sc)
        e = Xall[j][:n] - Zs[:, j - 1]
        L += h * B[j] * (ca.dot(q * e, e) + ca.dot(r * vs, vs))
    interval = ca.Function("interval", [Xa, Xcs, vs, Ds, Ss, Zs],
                           [ca.vertcat(*defects), L]).map(N)

    opti = ca.Opti()
    Xn = opti.variable(nx, N + 1)                    # scaled node states
    Xc = opti.variable(nx, N * d)                    # scaled collocation states
    V = opti.variable(n, N)
    Dc = np.array([M.D(t) for t in tc.ravel()]).reshape(N, d).T
    Sc = np.array([M.s_in(t) for t in tc.ravel()]).reshape(N, d).T
    defect, Lk = interval(Xn[:, :N], Xc, V, Dc, Sc, Zc.reshape(n, N * d))
    opti.subject_to(ca.vec(defect) == 0)
    opti.subject_to(Xn[:, 1:] == Xc[:, d - 1::d])    # Radau: end state = last colloc. point
    opti.subject_to(Xn[:, 0] == x0 / scale)
    opti.subject_to(ca.vec(Xc) >= 0)
    lo, hi = v_bounds
    if np.isfinite(lo) or np.isfinite(hi):
        opti.subject_to(opti.bounded(lo if np.isfinite(lo) else -ca.inf, ca.vec(V),
                                     hi if np.isfinite(hi) else ca.inf))
    J = ca.sum2(Lk)

    # Objective scaled so IPOPT sees O(1) numbers (cost of the initial guess with v = 0).
    e = ref_c[:n] - Zc
    J_scale = float((q[:, None, None] * e**2 * (h * B[1:])[None, None, :]).sum()) or 1.0
    opti.minimize(J / J_scale)
    opti.set_initial(Xn, ref_n / scale[:, None])
    opti.set_initial(Xc, ref_c.reshape(nx, N * d) / scale[:, None])
    if init_control is not None:
        opti.set_initial(V, np.clip(init_control, lo, hi))
    opts = {"ipopt.print_level": 5, "print_time": False, "ipopt.max_iter": 3000,
            "ipopt.tol": 1e-8, "ipopt.mu_strategy": "adaptive",
            "ipopt.print_frequency_iter": 10}
    opts.update(ipopt_opts or {})
    opti.solver("ipopt", opts)
    sol = opti.solve()
    return dict(t_nodes=tk, X=sol.value(Xn) * scale[:, None], V=np.atleast_2d(sol.value(V)),
                J_discrete=float(sol.value(J)), stats=opti.stats())


def piecewise_control(t_nodes, V):
    """v(t) for re-simulating a direct solution with scipy (zero-order hold)."""
    def v(t):
        k = np.clip(np.searchsorted(t_nodes, t, side="right") - 1, 0, V.shape[1] - 1)
        return V[:, k]
    return v
