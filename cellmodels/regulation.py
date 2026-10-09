"""Network dynamics: regulated expression of a pathogen gene during infection.

A pathogen gene is transcribed in proportion to pathogen biomass (colonization c(t)), with a
basal rate and an extra rate switched by a host signal s through a Hill function:

    dm/dt = c(t) * (alpha0 + alpha1 * H(s; K, n)) - delta * m,     H = s^n / (K^n + s^n)

Two hypotheses about aggressiveness determinants:

* fixed program: alpha1 = 0. Expression per unit pathogen is the same on every host.
* host-responsive: alpha1 > 0. Expression per unit pathogen depends on the host signal s_h.

Negative autoregulation, dm/dt = alpha / (1 + (m/K)^n) - delta m, is included as the standard
example of how network structure changes response time.
"""

from __future__ import annotations

from typing import Callable, Sequence

import numpy as np
from scipy.integrate import solve_ivp


def hill_activation(x, K, n):
    """x^n / (K^n + x^n)."""
    x = np.asarray(x, float)
    return x ** n / (K ** n + x ** n)


def hill_repression(x, K, n):
    """K^n / (K^n + x^n)."""
    return 1.0 - hill_activation(x, K, n)


def regulated_expression(t, alpha0, alpha1, delta, signal, K=1.0, n=2.0,
                         colonization: Callable[[float], float] = lambda s: 1.0, m0=0.0):
    """mRNA of a host-regulated pathogen gene at times t (host signal held constant)."""
    t = np.asarray(t, float)
    rate = alpha0 + alpha1 * float(hill_activation(signal, K, n))

    def rhs(s, y):
        return [colonization(s) * rate - delta * y[0]]

    sol = solve_ivp(rhs, (0.0, float(t[-1])), [m0], t_eval=t, rtol=1e-9, atol=1e-12)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[0]


def across_hosts(t, alpha0, alpha1, delta, signals: Sequence[float], **kw):
    """Stack regulated_expression over hosts with different signal levels: shape (hosts, times)."""
    return np.array([regulated_expression(t, alpha0, alpha1, delta, s, **kw) for s in signals])


def negative_autoregulation(t, alpha, delta, K, n=1.0, m0=0.0):
    """dm/dt = alpha / (1 + (m/K)^n) - delta m, solved numerically."""
    t = np.asarray(t, float)
    sol = solve_ivp(lambda s, y: [alpha / (1.0 + (y[0] / K) ** n) - delta * y[0]],
                    (0.0, float(t[-1])), [m0], t_eval=t, rtol=1e-9, atol=1e-12)
    return sol.y[0]


def rise_time(t, y, fraction=0.5):
    """First time y reaches `fraction` of its final value (linear interpolation)."""
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    target = fraction * y[-1]
    i = int(np.argmax(y >= target))
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (target - y[i - 1]) * (t[i] - t[i - 1]) / (y[i] - y[i - 1]))
