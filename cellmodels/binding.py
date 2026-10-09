"""Noncovalent binding: receptor R + ligand L <-> complex C.

Conventions: concentrations in mol/L (M), kon in 1/(M s), koff in 1/s, Kd = koff / kon in M.
"""

from __future__ import annotations

import numpy as np


def kd(kon: float, koff: float) -> float:
    """Equilibrium dissociation constant Kd = koff / kon."""
    return koff / kon


def fraction_bound(L, Kd):
    """Fraction of receptors bound at free ligand concentration L (ligand in excess)."""
    L = np.asarray(L, float)
    return L / (L + Kd)


def bound_with_depletion(R_tot, L_tot, Kd):
    """Equilibrium complex concentration when binding depletes the ligand.

    Solves C^2 - (R_tot + L_tot + Kd) C + R_tot L_tot = 0 and returns the physical root.
    """
    b = np.asarray(R_tot, float) + np.asarray(L_tot, float) + Kd
    return (b - np.sqrt(b * b - 4.0 * np.asarray(R_tot, float) * np.asarray(L_tot, float))) / 2.0


def observed_rate(L, kon, koff):
    """Pseudo-first-order relaxation rate kobs = kon L + koff (ligand in excess)."""
    return kon * np.asarray(L, float) + koff


def association(t, L, kon, koff, R_tot, C0=0.0):
    """Complex concentration C(t) after ligand L is added at t = 0 (ligand in excess)."""
    t = np.asarray(t, float)
    Ceq = R_tot * fraction_bound(L, koff / kon)
    return Ceq + (C0 - Ceq) * np.exp(-observed_rate(L, kon, koff) * t)


def dissociation(t, C0, koff):
    """Complex concentration C(t) after free ligand is removed at t = 0."""
    return C0 * np.exp(-koff * np.asarray(t, float))


def time_to_fraction(f, kobs):
    """Time for a first-order relaxation to cover fraction f of the way to equilibrium."""
    if not 0 < f < 1:
        raise ValueError("f must be in (0, 1)")
    return -np.log(1.0 - f) / kobs
