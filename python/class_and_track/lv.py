"""Generalised Lotka-Volterra (gLV) nitrification model: proof of concept of the article
(Article_Processes.tex, "Proof of concept"; MATLAB Tracking/synthetic_data.m).

    mu_i(s, x) = f_i(s) (1 + sum_j a_ij x_j),   f_i Monod in s1 (AOB) or s2 (NOB)

i.e. the NitrificationModel with v = A x.  Two ways to recover A from biomass data:
  * regress_interactions: the article's two-step approach, fit u(t) - 1 = A x(t) to the
    tracking control (linear least squares; MATLAB used fmincon + MultiStart for this).
  * estimate_interactions: fit A directly in the ODE (collocation + IPOPT), i.e. classical
    parameter estimation; tells whether A is identifiable from the data at all.
"""
from __future__ import annotations

import casadi as ca
import numpy as np
from scipy.integrate import solve_ivp

from .direct import _radau, _rhs_casadi
from .model import NitrificationModel

# Case study 1. A as in the article's figure (tex/Images/Stability/Interactions_case_study_1.png)
# and Tracking/synthetic_data.m; default_case.mat has 0.1 on the diagonal instead, and the
# article table gives K_2 = 0.3 while all code uses 1.3.  The article text says D = 0.24, but
# the article's POC images (250309_POC_*) were made by synthetic_data.m with D = 0.22.
CASE1 = dict(muA=0.77, muB=1.07, kSA=0.7, kSB=1.3, kA=1 / 0.251, kB=1 / 0.062,
             A=np.array([[1.0, 0.1], [-1.5, 1.0]]))


def case1_model(D=0.22, s_in_before=1.25, s_in_after=1.95, t_pert=150.0):
    c = CASE1
    M = NitrificationModel(
        nA=1, nB=1, kA=np.array([c["kA"]]), kB=np.array([c["kB"]]),
        muA=np.array([c["muA"]]), muB=np.array([c["muB"]]),
        kSA=np.array([c["kSA"]]), kSB=np.array([c["kSB"]]),
        D=lambda t: D, s_in=lambda t: s_in_before if t < t_pert else s_in_after)
    return M, c["A"].copy()


def simulate_glv(M: NitrificationModel, A, x0, t_span, t_eval=None):
    """gLV simulation (v = A x). s_in has a jump, so integrate the two pieces separately."""
    def rhs(t, X):
        return M.rhs(t, X, A @ np.maximum(X[:M.n], 0))
    return solve_ivp(rhs, t_span, x0, method="LSODA", t_eval=t_eval, rtol=1e-10,
                     atol=1e-12, max_step=0.5, dense_output=True)


def regress_interactions(t, U, X, t_window=None, weights=None):
    """Least squares  u_i(t) - 1 = sum_j a_ij x_j(t)  for each row i.
    t (T,), U (n, T) control u = 1 + v, X (n, T) biomasses. Returns (A_hat, residual rms)."""
    mask = np.ones_like(t, bool) if t_window is None else (t >= t_window[0]) & (t <= t_window[1])
    Phi = X[:, mask].T                                    # (T, n)
    Y = (U[:, mask] - 1).T                                # (T, n)
    w = np.ones(mask.sum()) if weights is None else weights[mask]
    sw = np.sqrt(w)[:, None]
    A_hat = np.linalg.lstsq(Phi * sw, Y * sw, rcond=None)[0].T
    res = np.sqrt(np.mean((Y - Phi @ A_hat.T) ** 2, 0))
    return A_hat, res


def estimate_interactions(M: NitrificationModel, t_data, Z, x0, t0, tf, q=None,
                          n_intervals=None, degree=3, S_data=None, s_weight=0.0,
                          A_init=None, fixed=None, ipopt_opts=None):
    """Fit A (n x n) in the gLV model to biomass data Z (n, len(t_data)) by collocation.
    Cost: sum_k q.(x(t_k) - z_k)^2 (+ s_weight * |s(t_k) - S_k|^2 if S_data given).
    The data times must lie on the grid t0 + k h (default h = 1 day). x0 is fixed.
    fixed: {(i, j): value} entries of A held fixed (for profiling)."""
    n, nx, d = M.n, M.n + 3, degree
    N = n_intervals or int(round(tf - t0))
    tk = np.linspace(t0, tf, N + 1)
    h = (tf - t0) / N
    tau, C, _ = _radau(d)
    tc = tk[:-1, None] + h * tau[None, 1:]
    q = np.ones(n) if q is None else np.asarray(q, float)
    k_data = np.rint((np.asarray(t_data) - t0) / h).astype(int)
    if not np.allclose(tk[k_data], t_data):
        raise ValueError("data times must lie on the collocation grid")

    # Initial guess: biomasses from the data (interpolated), substrates from a simulation
    # with A_init (default 0).  Starting the biomasses at the simulation instead made IPOPT
    # hit its iteration limit.
    A0 = np.zeros((n, n)) if A_init is None else np.asarray(A_init, float)
    sim = simulate_glv(M, A0, x0, (t0, tf))
    ref_n, ref_c = np.maximum(sim.sol(tk), 0), np.maximum(sim.sol(tc.ravel()), 0).reshape(nx, N, d)
    for i in range(n):
        ref_n[i] = np.interp(tk, t_data, Z[i])
        ref_c[i] = np.interp(tc.ravel(), t_data, Z[i]).reshape(N, d)
    scale = np.maximum(np.r_[np.abs(Z).max(1), np.abs(ref_n[n:]).max(1)], 1e-6)

    Xa, Xcs, As = ca.SX.sym("Xa", nx), ca.SX.sym("Xc", nx, d), ca.SX.sym("A", n, n)
    Ds, Ss = ca.SX.sym("D", d), ca.SX.sym("sin", d)
    sc = ca.DM(scale)
    Xall = [Xa * sc] + [Xcs[:, j] * sc for j in range(d)]
    defects = []
    for j in range(1, d + 1):
        xp = sum(C[r, j] * Xall[r] for r in range(d + 1))
        fj = _rhs_casadi(M, Xall[j], As @ Xall[j][:n], Ds[j - 1], Ss[j - 1])
        defects.append((h * fj - xp) / sc)
    interval = ca.Function("interval", [Xa, Xcs, As, Ds, Ss], [ca.vertcat(*defects)])
    interval = interval.map(N, [False, False, True, False, False], [False])   # A shared

    opti = ca.Opti()
    Xn, Xc = opti.variable(nx, N + 1), opti.variable(nx, N * d)
    # Fixed entries are substituted as constants (not equality constraints, which made
    # IPOPT's restoration phase fail); only the free entries are decision variables.
    fixed = fixed or {}
    free = [(i, j) for i in range(n) for j in range(n) if (i, j) not in fixed]
    a_free = opti.variable(len(free)) if free else None
    entries = {ij: a_free[k] for k, ij in enumerate(free)}
    entries.update({ij: ca.MX(float(v)) for ij, v in fixed.items()})
    A = ca.vertcat(*[ca.horzcat(*[entries[(i, j)] for j in range(n)]) for i in range(n)])
    Dc = np.array([M.D(t) for t in tc.ravel()]).reshape(N, d).T
    Sc = np.array([M.s_in(t) for t in tc.ravel()]).reshape(N, d).T
    opti.subject_to(ca.vec(interval(Xn[:, :N], Xc, A, Dc, Sc)) == 0)
    opti.subject_to(Xn[:, 1:] == Xc[:, d - 1::d])
    opti.subject_to(Xn[:, 0] == np.asarray(x0) / scale)
    opti.subject_to(ca.vec(Xc) >= 0)

    Xd =Xn[:n, k_data] * ca.repmat(ca.DM(scale[:n]), 1, len(k_data))
    E = (Xd - ca.DM(Z)) / ca.repmat(ca.DM(np.abs(Z).max(1)), 1, len(k_data))
    J = ca.sum2(ca.sum1(ca.repmat(ca.DM(q), 1, len(k_data)) * E**2))
    if S_data is not None and s_weight > 0:
        Sd = Xn[n:, k_data] * ca.repmat(ca.DM(scale[n:]), 1, len(k_data))
        Es = (Sd - ca.DM(S_data)) / ca.repmat(ca.DM(np.abs(S_data).max(1) + 1e-12), 1, len(k_data))
        J += s_weight * ca.sum2(ca.sum1(Es**2))
    opti.minimize(J)
    opti.set_initial(Xn, ref_n / scale[:, None])
    opti.set_initial(Xc, ref_c.reshape(nx, N * d) / scale[:, None])
    if free:
        opti.set_initial(a_free, [A0[i, j] for i, j in free])
    opts = {"ipopt.print_level": 0, "print_time": False, "ipopt.max_iter": 3000,
            "ipopt.tol": 1e-10, "ipopt.sb": "yes"}
    opts.update(ipopt_opts or {})
    opti.solver("ipopt", opts)
    sol = opti.solve()
    A_val = np.array([[float(sol.value(entries[(i, j)])) for j in range(n)] for i in range(n)])
    return dict(A=A_val, J=float(sol.value(J)),
                status=opti.stats()["return_status"], t_nodes=tk,
                X=sol.value(Xn) * scale[:, None])
