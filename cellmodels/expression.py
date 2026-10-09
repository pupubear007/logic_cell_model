"""Gene expression as a two-stage model.

    dm/dt = alpha_m * u(t) - delta_m * m        (mRNA; u is an input, 1 by default)
    dp/dt = alpha_p * m    - delta_p * p        (protein)

alpha_m: transcription rate (molecules/time), delta_m: mRNA decay rate (1/time),
alpha_p: translation rate per mRNA (1/time), delta_p: protein loss rate (dilution + decay).
RNA-seq measures m (in relative units), not p and not the rates themselves.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp


def steady_state(alpha_m, delta_m, alpha_p=None, delta_p=None):
    """Steady state (m*, p*) for constant input u = 1. Returns m* alone if protein rates are None."""
    m = alpha_m / delta_m
    if alpha_p is None or delta_p is None:
        return m
    return m, alpha_p * m / delta_p


def mrna(t, alpha_m, delta_m, m0=0.0):
    """mRNA after a step to constant transcription at t = 0."""
    t = np.asarray(t, float)
    ms = alpha_m / delta_m
    return ms + (m0 - ms) * np.exp(-delta_m * t)


def mrna_shutoff(t, m0, delta_m):
    """mRNA after transcription is blocked at t = 0; separates decay from synthesis."""
    return m0 * np.exp(-delta_m * np.asarray(t, float))


def protein(t, alpha_m, delta_m, alpha_p, delta_p, m0=0.0, p0=0.0):
    """Protein after a step to constant transcription at t = 0 (closed form)."""
    t = np.asarray(t, float)
    ms, ps = steady_state(alpha_m, delta_m, alpha_p, delta_p)
    A = m0 - ms
    if np.isclose(delta_p, delta_m):
        d = delta_m
        return ps + (p0 - ps) * np.exp(-d * t) + alpha_p * A * t * np.exp(-d * t)
    c = alpha_p * A / (delta_p - delta_m)
    B = p0 - ps - c
    return ps + B * np.exp(-delta_p * t) + c * np.exp(-delta_m * t)


def simulate(t, alpha_m, delta_m, alpha_p, delta_p,
             u: Callable[[float], float] = lambda s: 1.0, m0=0.0, p0=0.0):
    """Numerical solution for a time-varying input u(t), e.g. pathogen colonization.

    Returns arrays (m, p) at times t (t must be increasing and start at or after 0).
    """
    t = np.asarray(t, float)

    def rhs(s, y):
        return [alpha_m * u(s) - delta_m * y[0], alpha_p * y[0] - delta_p * y[1]]

    sol = solve_ivp(rhs, (0.0, float(t[-1])), [m0, p0], t_eval=t, rtol=1e-9, atol=1e-12)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[0], sol.y[1]
