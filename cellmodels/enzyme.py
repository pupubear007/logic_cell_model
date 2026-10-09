"""Enzyme kinetics: Michaelis–Menten rates, inhibition, and progress curves.

v = Vmax S / (Km + S), with Vmax = kcat E_total. Units: concentrations in M, time in s.
Plant-pathology use: secreted cell-wall-degrading enzymes acting on their substrates.
"""

from __future__ import annotations

import numpy as np
from scipy.special import lambertw


def michaelis_menten(S, Vmax, Km):
    """Initial rate v = Vmax S / (Km + S)."""
    S = np.asarray(S, float)
    return Vmax * S / (Km + S)


def specificity_constant(kcat, Km):
    """kcat / Km: the second-order rate constant that governs the rate at S << Km."""
    return kcat / Km


def competitive(S, I, Vmax, Km, Ki):
    """Competitive inhibitor: apparent Km = Km (1 + I/Ki), Vmax unchanged."""
    return michaelis_menten(S, Vmax, Km * (1.0 + np.asarray(I, float) / Ki))


def noncompetitive(S, I, Vmax, Km, Ki):
    """Pure noncompetitive inhibitor: apparent Vmax = Vmax / (1 + I/Ki), Km unchanged."""
    return michaelis_menten(S, Vmax / (1.0 + np.asarray(I, float) / Ki), Km)


def uncompetitive(S, I, Vmax, Km, Ki):
    """Uncompetitive inhibitor: Vmax and Km both divided by (1 + I/Ki)."""
    f = 1.0 + np.asarray(I, float) / Ki
    return michaelis_menten(S, Vmax / f, Km / f)


def progress_curve(t, S0, Vmax, Km):
    """Substrate remaining at time t, from the integrated Michaelis–Menten equation.

    S(t) = Km W( (S0/Km) exp((S0 - Vmax t)/Km) ), with W the Lambert W function.
    """
    t = np.asarray(t, float)
    arg = (S0 / Km) * np.exp((S0 - Vmax * t) / Km)
    return Km * np.real(lambertw(arg))
