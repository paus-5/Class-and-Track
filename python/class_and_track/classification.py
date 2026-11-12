"""OTU -> functional group classification (port of Classification/MIP_*.m).

Decision variables (same order as MATLAB):
    eA (T), eB (T)  : absolute errors of the two nitrogen balances
    yA (n), yB (n)  : per-OTU inverse yields if the OTU is AOB / NOB
    a (n),  b (n)   : binaries, OTU i is AOB / NOB (a_i + b_i <= 1)

    |W^T yA - xG1| <= eA ,  |W^T yB - xG2| <= eB
    m1 a <= yA <= M1 a ,    m2 b <= yB <= M2 b      (big-M; yields within +-delta of reference)

objective  "l1": sum eA + sum eB                  -> MILP  (MIP_intlinprog.m)
           "l2": sum eA^2 + sum eB^2              -> MIQP  (MIP_gurobi.m)

Free solvers: HiGHS through scipy.optimize.milp for the MILP; SCIP (pyscipopt) for
the MIQP (HiGHS does not do mixed-integer QP).
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.optimize import Bounds, LinearConstraint, milp


def _structure(W, xG1, xG2, yA_ref, yB_ref, delta):
    n, T = W.shape
    M1, M2 = 1 / (yA_ref * (1 - delta)), 1 / (yB_ref * (1 - delta))
    m1, m2 = 1 / (yA_ref * (1 + delta)), 1 / (yB_ref * (1 + delta))
    I_T, I_n = sparse.eye(T), sparse.eye(n)
    Z = lambda r, c: sparse.csr_matrix((r, c))
    Wt = sparse.csr_matrix(W.T)
    A = sparse.vstack([
        sparse.hstack([-I_T, Z(T, T), Wt, Z(T, 3 * n)]),
        sparse.hstack([-I_T, Z(T, T), -Wt, Z(T, 3 * n)]),
        sparse.hstack([Z(T, T), -I_T, Z(T, n), Wt, Z(T, 2 * n)]),
        sparse.hstack([Z(T, T), -I_T, Z(T, n), -Wt, Z(T, 2 * n)]),
        sparse.hstack([Z(n, 2 * T), I_n, Z(n, n), -M1 * I_n, Z(n, n)]),
        sparse.hstack([Z(n, 2 * T), Z(n, n), I_n, Z(n, n), -M2 * I_n]),
        sparse.hstack([Z(n, 2 * T), -I_n, Z(n, n), m1 * I_n, Z(n, n)]),
        sparse.hstack([Z(n, 2 * T), Z(n, n), -I_n, Z(n, n), m2 * I_n]),
        sparse.hstack([Z(n, 2 * T + 2 * n), I_n, I_n]),
    ]).tocsr()
    ub = np.concatenate([xG1, -xG1, xG2, -xG2, np.zeros(4 * n), np.ones(n)])
    lb_x = np.zeros(2 * T + 4 * n)
    ub_x = np.concatenate([np.full(2 * T, np.inf), np.full(2 * n, M2), np.ones(2 * n)])
    integrality = np.r_[np.zeros(2 * T + 2 * n), np.ones(2 * n)]
    return A, ub, lb_x, ub_x, integrality, n, T


def _unpack(sol, n, T):
    return dict(error_AOB=sol[:T], error_NOB=sol[T:2 * T],
                y_AOB=sol[2 * T:2 * T + n], y_NOB=sol[2 * T + n:2 * T + 2 * n],
                a=np.round(sol[2 * T + 2 * n:2 * T + 3 * n]), b=np.round(sol[2 * T + 3 * n:]))


def classify(W, xG1, xG2, yA_ref, yB_ref, delta=0.01, objective="l1", time_limit=600.0):
    """W: (n_otu, T) OTU biomasses at the observation times; xG1, xG2: (T,) balances."""
    A, ub, lb_x, ub_x, integ, n, T = _structure(np.asarray(W, float), np.asarray(xG1, float),
                                                np.asarray(xG2, float), yA_ref, yB_ref, delta)
    if objective == "l1":
        c = np.r_[np.ones(2 * T), np.zeros(4 * n)]
        res = milp(c, constraints=LinearConstraint(A, -np.inf, ub), integrality=integ,
                   bounds=Bounds(lb_x, ub_x), options=dict(time_limit=time_limit))
        if res.x is None:
            raise RuntimeError(res.message)
        out = _unpack(res.x, n, T)
        out.update(objective=res.fun, status=res.message)
        return out
    if objective == "l2":
        from pyscipopt import Model, quicksum
        m = Model()
        m.hideOutput()
        m.setParam("limits/time", time_limit)
        N = A.shape[1]
        x = [m.addVar(lb=lb_x[j], ub=None if np.isinf(ub_x[j]) else ub_x[j],
                      vtype="B" if integ[j] else "C") for j in range(N)]
        A = A.tocsr()
        for i in range(A.shape[0]):
            row = A.getrow(i)
            m.addCons(quicksum(v * x[j] for j, v in zip(row.indices, row.data)) <= ub[i])
        # SCIP takes a linear objective: minimise an epigraph variable for the quadratic.
        tau = m.addVar(lb=0)
        m.addCons(quicksum(x[j] * x[j] for j in range(2 * T)) <= tau)
        m.setObjective(tau, "minimize")
        m.optimize()
        sol = np.array([m.getVal(v) for v in x])
        out = _unpack(sol, n, T)
        out.update(objective=m.getObjVal(), status=m.getStatus())
        return out
    raise ValueError(objective)
