"""Two-step nitrification chemostat model with OTU-level biomasses.

State X = (x_1..x_n, s1, s2, s3):  x_1..x_nA are AOB (group G1), x_nA+1..x_n are NOB (G2),
s1 = ammonium, s2 = nitrite, s3 = nitrate.

    x_i'  = (f_i(s) (1 + v_i) - D) x_i
    s1'   = (s_in - s1) D - sum_{G1} kA_i f_i (1+v_i) x_i
    s2'   = -s2 D + sum_{G1} kA_i f_i (1+v_i) x_i - sum_{G2} kB_i f_i (1+v_i) x_i
    s3'   = -s3 D + sum_{G2} kB_i f_i (1+v_i) x_i

v = u - 1 is the "interaction" control (u_i = I_i(x) in the article).  The MATLAB code
(Ax.m / Bx.m) evaluated this as A(X) X + B(X) v; here the right-hand side is written
directly, which is algebraically identical but avoids the (s_in/s1 - 1) * s1 term
that becomes 0 * inf when s1 hits 0.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


def monod(s, mu_max, ks):
    return mu_max * s / (s + ks)


def haldane(s, mu_max, ks, ki):
    return mu_max * s / (s + ks + s**2 / ki)


@dataclass
class NitrificationModel:
    nA: int
    nB: int
    kA: np.ndarray          # 1/yield of each AOB (g N / g biomass)
    kB: np.ndarray          # 1/yield of each NOB
    muA: np.ndarray
    muB: np.ndarray
    kSA: np.ndarray
    kSB: np.ndarray
    D: Callable[[float], float]      # dilution rate D(t)
    s_in: Callable[[float], float]   # inlet ammonium s_in(t)

    @property
    def n(self) -> int:
        return self.nA + self.nB

    def growth(self, s: np.ndarray) -> np.ndarray:
        """f(s): AOB grow on s1, NOB on s2 (Monod, as in the MATLAB code)."""
        return np.concatenate([monod(s[0], self.muA, self.kSA),
                               monod(s[1], self.muB, self.kSB)])

    def rhs(self, t: float, X: np.ndarray, v: np.ndarray | float = 0.0) -> np.ndarray:
        n, nA = self.n, self.nA
        X = np.maximum(X, 0.0)               # replaces odeset('NonNegative', ...)
        x, s = X[:n], X[n:]
        D = self.D(t)
        r = self.growth(s) * (1.0 + v) * x   # specific growth * biomass
        rA, rB = r[:nA], r[nA:]
        dx = r - D * x
        ds1 = (self.s_in(t) - s[0]) * D - self.kA @ rA
        ds2 = -s[1] * D + self.kA @ rA - self.kB @ rB
        ds3 = -s[2] * D + self.kB @ rB
        return np.concatenate([dx, [ds1, ds2, ds3]])

    def B1_diag(self, X: np.ndarray) -> np.ndarray:
        """Diagonal of the x-block of B(X) = diag(f(s)) diag(x)."""
        n = self.n
        return self.growth(X[n:]) * X[:n]
