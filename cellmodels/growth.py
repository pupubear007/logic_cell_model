"""Growth of hyphae, colonies and lesions."""

from __future__ import annotations

import numpy as np


def radial_front(t, rate, lag=0.0, r0=0.0):
    """Colony or lesion radius growing at a constant front speed after a lag."""
    t = np.asarray(t, float)
    return r0 + rate * np.clip(t - lag, 0.0, None)


def fit_radial_rate(t, r):
    """Least-squares front speed and lag from radius measurements on the linear phase.

    Points with r <= 0 are treated as still in the lag phase and ignored. Returns (rate, lag).
    """
    t = np.asarray(t, float)
    r = np.asarray(r, float)
    keep = r > 0
    if keep.sum() < 2:
        raise ValueError("need at least two points with r > 0")
    rate, intercept = np.polyfit(t[keep], r[keep], 1)
    return float(rate), float(-intercept / rate)


def logistic(t, rate, capacity, n0):
    """Logistic growth N(t) = K / (1 + (K/N0 - 1) exp(-r t))."""
    t = np.asarray(t, float)
    return capacity / (1.0 + (capacity / n0 - 1.0) * np.exp(-rate * t))


def doubling_time(rate):
    """Doubling time of exponential growth at specific rate r."""
    return np.log(2.0) / rate


def lesion_area(radius):
    """Area of a circular lesion."""
    return np.pi * np.asarray(radius, float) ** 2
