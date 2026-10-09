"""Flow in microchannels, where the Navier–Stokes equations reduce to Stokes flow.

At small Reynolds number the inertial term of the Navier–Stokes equations is negligible, the
flow is laminar and linear in the pressure gradient, and mixing across a channel happens by
diffusion. Units: SI (Pa, Pa s, m, m/s, m^3/s).
"""

from __future__ import annotations

import numpy as np


def reynolds(rho, U, L, mu):
    """Re = rho U L / mu: inertia relative to viscosity."""
    return rho * U * L / mu


def peclet(U, L, D):
    """Pe = U L / D: transport by flow relative to diffusion."""
    return U * L / D


def damkohler(k, L, U):
    """Da = k L / U: reaction rate relative to transport by flow."""
    return k * L / U


def parallel_plate_velocity(y, h, dpdx, mu):
    """Pressure-driven (Poiseuille) velocity between plates at y = 0 and y = h.

    u(y) = (-dp/dx) / (2 mu) * y (h - y), an exact solution of the Stokes equations.
    """
    y = np.asarray(y, float)
    return (-dpdx) / (2.0 * mu) * y * (h - y)


def parallel_plate_mean_velocity(h, dpdx, mu):
    """Mean velocity (-dp/dx) h^2 / (12 mu)."""
    return (-dpdx) * h * h / (12.0 * mu)


def wall_shear_stress(Q, w, h, mu):
    """Wall shear stress 6 mu Q / (w h^2) in a wide, shallow channel (h << w)."""
    return 6.0 * mu * Q / (w * h * h)


def hydraulic_resistance_rect(w, h, L, mu):
    """Approximate hydraulic resistance of a rectangular channel (h < w), Delta p = R Q.

    R ~= 12 mu L / (w h^3 (1 - 0.63 h / w)), accurate to a few percent for h <= w.
    """
    if h > w:
        w, h = h, w
    return 12.0 * mu * L / (w * h ** 3 * (1.0 - 0.63 * h / w))


def mixing_length(U, w, D):
    """Downstream distance for diffusion to mix across a channel of width w: about U w^2 / D."""
    return U * w * w / D
