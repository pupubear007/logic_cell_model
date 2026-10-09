import numpy as np
import pytest

from cellmodels import enzyme, regulation, stochastic
from cellmodels.identifiability import local_identifiability, witness


# --- enzyme kinetics ---------------------------------------------------------------------------

def test_half_maximal_rate_at_km():
    assert enzyme.michaelis_menten(2e-3, 1.0, 2e-3) == pytest.approx(0.5)


def test_inhibition_limits():
    S, Vmax, Km, Ki = 1e-3, 1.0, 1e-3, 1e-6
    assert enzyme.competitive(S, 0.0, Vmax, Km, Ki) == pytest.approx(enzyme.michaelis_menten(S, Vmax, Km))
    # saturating substrate overcomes a competitive inhibitor but not a noncompetitive one
    big = 1e3
    assert enzyme.competitive(big, 1e-5, Vmax, Km, Ki) == pytest.approx(Vmax, rel=1e-3)
    assert enzyme.noncompetitive(big, 1e-5, Vmax, Km, Ki) == pytest.approx(Vmax / 11, rel=1e-3)


def test_progress_curve_solves_michaelis_menten():
    S0, Vmax, Km = 5e-3, 1e-4, 2e-3
    t = np.linspace(0, 200, 4001)
    S = enzyme.progress_curve(t, S0, Vmax, Km)
    assert S[0] == pytest.approx(S0)
    dSdt = np.gradient(S, t)
    assert np.allclose(dSdt[2:-2], -enzyme.michaelis_menten(S[2:-2], Vmax, Km), rtol=1e-3)


def _rates(S):
    return lambda th: enzyme.michaelis_menten(S, th[0], th[1])


def test_low_substrate_rates_identify_only_specificity_constant():
    low = np.array([1e-6, 2e-6, 4e-6])            # S << Km
    r = local_identifiability(_rates(low), [1.0, 1e-3], ("Vmax", "Km"), rel_tol=1e-2)
    assert r.rank == 1                            # only Vmax/Km
    wide = np.array([1e-4, 1e-3, 1e-2])           # around and above Km
    assert local_identifiability(_rates(wide), [1.0, 1e-3]).identifiable


# --- network dynamics: fixed program vs host-responsive ----------------------------------------

T = np.array([1.0, 2.0, 4.0, 8.0])


def test_fixed_program_gives_same_expression_on_every_host():
    y = regulation.across_hosts(T, 2.0, 0.0, 0.5, signals=[0.1, 1.0, 10.0])
    assert np.allclose(y, y[0])


def test_host_responsive_differs_across_hosts():
    y = regulation.across_hosts(T, 2.0, 3.0, 0.5, signals=[0.1, 10.0])
    assert not np.allclose(y[0], y[1])


def test_one_host_cannot_separate_basal_from_induced_transcription():
    one = lambda th: regulation.regulated_expression(T, th[0], th[1], th[2], signal=1.0)
    r = local_identifiability(one, [2.0, 3.0, 0.5], ("alpha0", "alpha1", "delta"))
    assert r.rank == 2
    two = lambda th: regulation.across_hosts(T, th[0], th[1], th[2], signals=[0.2, 5.0])
    assert local_identifiability(two, [2.0, 3.0, 0.5]).identifiable


def test_witness_pair_on_one_host():
    # alpha0 + alpha1 * H(1) is the same: H(1; K=1, n=2) = 1/2
    one = lambda th: regulation.regulated_expression(T, th[0], th[1], 0.5, signal=1.0)
    assert witness(one, [2.0, 2.0], [1.0, 4.0])


def test_negative_autoregulation_speeds_response():
    t = np.linspace(0, 20, 4001)
    nar = regulation.negative_autoregulation(t, alpha=10.0, delta=0.5, K=0.5, n=2.0)
    simple = nar[-1] * (1 - np.exp(-0.5 * t))  # same steady state and decay rate, no feedback
    assert regulation.rise_time(t, nar) < regulation.rise_time(t, simple)


# --- stochastic expression ---------------------------------------------------------------------

def test_birth_death_is_poisson():
    x = stochastic.gillespie_birth_death(20.0, 1.0, t_end=10.0, n_cells=3000, rng=1)
    assert x.mean() == pytest.approx(20.0, rel=0.03)
    assert stochastic.fano(x) == pytest.approx(1.0, abs=0.1)


def test_bursty_moments():
    a, b, delta = 4.0, 5.0, 1.0
    x = stochastic.gillespie_bursty(a, b, delta, t_end=10.0, n_cells=3000, rng=2)
    mean, fano = stochastic.bursty_moments(a, b, delta)
    assert x.mean() == pytest.approx(mean, rel=0.05)
    assert stochastic.fano(x) == pytest.approx(fano, rel=0.1)


def test_bulk_data_cannot_separate_burst_frequency_from_size():
    bulk = lambda th: np.array([stochastic.bursty_moments(th[0], th[1], 1.0)[0]])
    single_cell = lambda th: np.array(stochastic.bursty_moments(th[0], th[1], 1.0))
    assert witness(bulk, [4.0, 5.0], [10.0, 2.0])               # same mean
    assert local_identifiability(bulk, [4.0, 5.0]).rank == 1
    assert local_identifiability(single_cell, [4.0, 5.0]).identifiable
