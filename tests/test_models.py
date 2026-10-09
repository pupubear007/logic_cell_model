import numpy as np
import pytest

from cellmodels import binding, expression, growth, microfluidics, transport
from cellmodels.identifiability import local_identifiability, witness


# --- binding -----------------------------------------------------------------------------------

def test_half_bound_at_kd():
    assert binding.fraction_bound(1e-9, 1e-9) == pytest.approx(0.5)


def test_depletion_reduces_to_excess_ligand():
    Kd, R = 1e-9, 1e-12  # receptor far below ligand: no depletion
    L = np.array([1e-10, 1e-9, 1e-8])
    assert np.allclose(binding.bound_with_depletion(R, L, Kd), R * binding.fraction_bound(L, Kd),
                       rtol=1e-3)


def test_association_relaxes_to_equilibrium_at_kobs():
    kon, koff, L, R = 1e6, 1e-3, 1e-9, 1.0
    t = binding.time_to_fraction(0.5, binding.observed_rate(L, kon, koff))
    Ceq = R * binding.fraction_bound(L, koff / kon)
    assert binding.association(t, L, kon, koff, R) == pytest.approx(Ceq / 2)


# --- expression --------------------------------------------------------------------------------

@pytest.mark.parametrize("delta_p", [0.3, 1.0])  # distinct and equal decay rates
def test_protein_closed_form_matches_ode(delta_p):
    t = np.linspace(0, 20, 41)
    m, p = expression.simulate(t, 2.0, 1.0, 0.5, delta_p)
    assert np.allclose(m, expression.mrna(t, 2.0, 1.0), atol=1e-6)
    assert np.allclose(p, expression.protein(t, 2.0, 1.0, 0.5, delta_p), atol=1e-6)


def test_steady_state():
    assert expression.steady_state(2.0, 0.5, 3.0, 1.5) == pytest.approx((4.0, 8.0))


# --- identifiability (assay resolution) ---------------------------------------------------------

def snapshot(theta):
    """RNA-seq at steady state: one number, alpha / delta."""
    a, d = theta
    return np.array([expression.steady_state(a, d)])


def time_course(theta):
    """RNA-seq after induction at several times."""
    a, d = theta
    return expression.mrna(np.array([0.5, 1.0, 2.0, 4.0]), a, d)


def test_snapshot_cannot_separate_synthesis_from_decay():
    assert witness(snapshot, [2.0, 1.0], [4.0, 2.0])          # same ratio: unresolved
    assert not witness(time_course, [2.0, 1.0], [4.0, 2.0])   # time course separates them


def test_snapshot_rank_one_with_ratio_direction_free():
    r = local_identifiability(snapshot, [2.0, 1.0], names=("alpha", "delta"))
    assert r.rank == 1 and not r.identifiable
    v = r.null_directions[0]
    assert abs(abs(v[0]) - abs(v[1])) < 1e-6 and np.sign(v[0]) == np.sign(v[1])  # scale both


def test_time_course_identifies_both_rates():
    assert local_identifiability(time_course, [2.0, 1.0]).identifiable


def test_mrna_data_never_identify_protein_rates():
    t = np.array([0.5, 1.0, 2.0, 4.0])

    def rnaseq(theta):
        a_m, d_m, a_p, d_p = theta
        return expression.mrna(t, a_m, d_m)  # protein rates do not enter the observation

    r = local_identifiability(rnaseq, [2.0, 1.0, 0.5, 0.3])
    assert r.rank == 2


# --- growth ------------------------------------------------------------------------------------

def test_fit_radial_rate_recovers_parameters():
    t = np.arange(0, 10.0)
    rate, lag = growth.fit_radial_rate(t, growth.radial_front(t, 1.5, lag=2.0))
    assert rate == pytest.approx(1.5) and lag == pytest.approx(2.0)


def test_logistic_limits():
    assert growth.logistic(0.0, 1.0, 100.0, 1.0) == pytest.approx(1.0)
    assert growth.logistic(100.0, 1.0, 100.0, 1.0) == pytest.approx(100.0)


# --- transport ---------------------------------------------------------------------------------

def test_slab_profile_solves_reaction_diffusion():
    D, k, L, c0 = 1e-9, 1e-3, 1e-3, 1.0
    x = np.linspace(0, L, 2001)
    c = transport.slab_profile(x, c0, L, D, k)
    d2 = np.gradient(np.gradient(c, x), x)
    assert c[0] == pytest.approx(c0)
    assert np.allclose(D * d2[5:-5], k * c[5:-5], rtol=1e-3)
    assert np.mean(c) == pytest.approx(float(transport.effectiveness_factor(
        transport.thiele_modulus(L, D, k))), rel=1e-3)


def test_deep_slab_matches_semi_infinite():
    D, k = 1e-9, 1e-2
    lam = transport.penetration_length(D, k)
    x = np.linspace(0, 3 * lam, 5)
    assert np.allclose(transport.slab_profile(x, 1.0, 50 * lam, D, k),
                       transport.semi_infinite_profile(x, 1.0, D, k), rtol=1e-6)


# --- microfluidics -----------------------------------------------------------------------------

def test_poiseuille_mean_velocity():
    h, dpdx, mu = 50e-6, -1e4, 1e-3
    y = np.linspace(0, h, 20001)
    u = microfluidics.parallel_plate_velocity(y, h, dpdx, mu)
    assert np.trapezoid(u, y) / h == pytest.approx(
        microfluidics.parallel_plate_mean_velocity(h, dpdx, mu), rel=1e-6)


def test_typical_microchannel_is_stokes_flow():
    # water, 100 um/s, 50 um channel
    assert microfluidics.reynolds(1000.0, 100e-6, 50e-6, 1e-3) < 1e-2
