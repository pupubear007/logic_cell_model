"""Diffusion with reaction, e.g. a secreted acid spreading into tissue while being consumed.

Steady diffusion with first-order consumption, D c'' = k c, has the length scale
lambda = sqrt(D / k): the distance over which the concentration falls by a factor e.
Units: D in m^2/s, k in 1/s, lengths in m.
"""

from __future__ import annotations

import numpy as np


def penetration_length(D, k):
    """Reaction–diffusion length sqrt(D / k)."""
    return np.sqrt(D / k)


def thiele_modulus(L, D, k):
    """phi = L sqrt(k / D): ratio of tissue thickness to penetration length."""
    return L * np.sqrt(k / D)


def effectiveness_factor(phi):
    """tanh(phi) / phi: mean concentration in the slab relative to the source concentration."""
    phi = np.asarray(phi, float)
    return np.where(phi > 0, np.tanh(phi) / np.where(phi > 0, phi, 1.0), 1.0)


def slab_profile(x, c0, L, D, k):
    """Steady profile in a slab 0 <= x <= L, source c0 at x = 0, no flux at x = L."""
    phi = thiele_modulus(L, D, k)
    x = np.asarray(x, float)
    return c0 * np.cosh(phi * (1.0 - x / L)) / np.cosh(phi)


def semi_infinite_profile(x, c0, D, k):
    """Steady profile into a deep tissue: c0 exp(-x / lambda)."""
    return c0 * np.exp(-np.asarray(x, float) / penetration_length(D, k))


def diffusion_time(L, D):
    """Characteristic time L^2 / D to diffuse a distance L."""
    return L * L / D
