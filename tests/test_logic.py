"""The Python logic matches the Lean theorems in lean/CellLogic and WangLogic.

Each test names the Lean theorem it mirrors. The Lean proofs hold for all real parameter values;
these tests check the same statements on finite grids of worlds.
"""

import itertools

import numpy as np
import pytest

from cellmodels import enzyme, expression, stochastic
from cellmodels.logic import entails, resolves, resolves_via_decides, supported, witnesses

GRID = [0.5, 1.0, 2.0, 4.0]


def key(x):
    return tuple(np.round(np.atleast_1d(x), 9))


# Theorem 8.2 holds on every finite example ---------------------------------------------------

@pytest.mark.parametrize("seed", range(5))
def test_resolves_iff_decides(seed):
    rng = np.random.default_rng(seed)
    worlds = list(range(12))
    labels = rng.integers(0, 3, 12)
    truth = rng.integers(0, 2, 12).astype(bool)
    a = lambda w: int(labels[w])
    H = lambda w: bool(truth[w])
    assert resolves(a, H, worlds) == resolves_via_decides(a, H, worlds)


# CellLogic.Expression ----------------------------------------------------------------------------

TURNOVER = [(a, d) for a, d in itertools.product(GRID, GRID)]
fast = lambda w: w[1] > 1.0


def test_snapshot_not_resolves_fast():          # CellLogic.snapshot_not_resolves_fast
    snapshot = lambda w: key(expression.steady_state(*w))
    assert not resolves(snapshot, fast, TURNOVER)
    assert ((2.0, 1.0), (4.0, 2.0)) in witnesses(snapshot, fast, TURNOVER)


def test_snapshot_shutoff_resolves():           # CellLogic.snapshotShutoff_resolves
    design = lambda w: key([expression.steady_state(*w), expression.mrna_shutoff(1.0, 1.0, w[1])])
    assert resolves(design, fast, TURNOVER)


# CellLogic.Regulation ----------------------------------------------------------------------------

REG = [(b, i) for b, i in itertools.product([0.0, 1.0, 2.0], [0.0, 1.0, 2.0])]
responsive = lambda w: w[1] != 0
one_host = lambda s: (lambda w: key(w[0] + w[1] * s))


def test_one_host_not_resolves():               # CellLogic.oneHost_not_resolves
    assert not resolves(one_host(1.0), responsive, REG)


def test_two_hosts_resolve():                   # CellLogic.twoHosts_resolves
    two = lambda w: key([w[0] + w[1] * 0.2, w[0] + w[1] * 0.9])
    assert resolves(two, responsive, REG)


def test_fixed_program_predicts_equal_and_is_refuted():  # fixedProgram_predicts_equal / _refuted
    fixed = lambda w: w[1] == 0
    equal = lambda w: one_host(0.2)(w) == one_host(0.9)(w)
    assert entails(fixed, equal, REG)
    observed_difference = [lambda w: not equal(w)]
    assert not supported(observed_difference, fixed, REG)


# CellLogic.Enzyme --------------------------------------------------------------------------------

MM = list(itertools.product(GRID, GRID))         # (Vmax, Km)
high_affinity = lambda w: w[1] <= 1.0


def test_low_substrate_not_resolves():          # CellLogic.lowSubstrate_not_resolves
    slope = lambda w: key(w[0] / w[1])
    assert not resolves(slope, high_affinity, MM)


def test_two_substrates_resolve():              # CellLogic.twoSubstrates_resolves
    two = lambda w: key([enzyme.michaelis_menten(0.5, *w), enzyme.michaelis_menten(3.0, *w)])
    assert resolves(two, high_affinity, MM)


# CellLogic.Bursting ------------------------------------------------------------------------------

BURSTS = list(itertools.product(GRID, GRID))     # (frequency, size)
large = lambda w: w[1] > 1.0


def test_bulk_not_resolves():                   # CellLogic.bulk_not_resolves
    assert not resolves(lambda w: key(stochastic.bursty_moments(*w, 1.0)[0]), large, BURSTS)


def test_single_cell_resolves():                # CellLogic.singleCell_resolves
    assert resolves(lambda w: key(stochastic.bursty_moments(*w, 1.0)), large, BURSTS)
