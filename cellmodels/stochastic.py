"""Stochastic gene expression: few molecules per cell, simulated exactly (Gillespie).

Constitutive birth–death: production at rate alpha, decay at rate delta per molecule. The
stationary distribution is Poisson with mean alpha/delta, so the Fano factor
(variance / mean) is 1.

Bursty expression: bursts arrive at rate a, each adds a geometric number of molecules with mean
b. Then mean = a b / delta and Fano factor = 1 + b. Bulk (population-average) data see only
a b; single-cell variance separates burst frequency a from burst size b.
"""

from __future__ import annotations

import numpy as np


def poisson_moments(alpha, delta):
    """(mean, Fano factor) of the constitutive birth–death process."""
    return alpha / delta, 1.0


def bursty_moments(a, b, delta):
    """(mean, Fano factor) for bursts at rate a with geometric size of mean b."""
    return a * b / delta, 1.0 + b


def gillespie_bursty(a, b, delta, t_end, n_cells, rng=None, m0=0):
    """Exact simulation of bursty expression; returns molecule counts per cell at t_end.

    For constitutive production (one molecule at a time) use ``gillespie_birth_death``.
    """
    rng = np.random.default_rng(rng)
    p = 1.0 / (1.0 + b)  # geometric on {0, 1, 2, ...} with mean b
    out = np.empty(n_cells, dtype=int)
    for c in range(n_cells):
        t, m = 0.0, m0
        while True:
            total = a + delta * m
            t += rng.exponential(1.0 / total)
            if t > t_end:
                break
            if rng.random() < a / total:
                m += rng.geometric(p) - 1
            else:
                m -= 1
        out[c] = m
    return out


def gillespie_birth_death(alpha, delta, t_end, n_cells, rng=None, m0=0):
    """Exact simulation of constitutive production; returns counts per cell at t_end."""
    rng = np.random.default_rng(rng)
    out = np.empty(n_cells, dtype=int)
    for c in range(n_cells):
        t, m = 0.0, m0
        while True:
            total = alpha + delta * m
            t += rng.exponential(1.0 / total)
            if t > t_end:
                break
            m += 1 if rng.random() < alpha / total else -1
        out[c] = m
    return out


def fano(x):
    """Variance / mean of a sample."""
    x = np.asarray(x, float)
    return float(x.var(ddof=1) / x.mean())
