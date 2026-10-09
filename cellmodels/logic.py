"""The logic of *Deduction and Induction in One Diagram*, for finite sets of worlds.

These functions mirror the Lean definitions in ``WangLogic`` (namespace ``WangLogic``) so the
same reasoning runs on real data, where the worlds are isolates, samples or candidate
mechanisms:

=========================  ================================  ============================
Python                     Lean (WangLogic)                  Paper
=========================  ================================  ============================
``entails(phi, psi, W)``   ``Entails``                       Definition 3.2
``supported(obs, t, W)``   ``Supported``                     Definition 3.4
``resolves(a, H, W)``      ``Resolves``                      Definition 8.1
``decides(a, r, H, W)``    ``Decides (AssayResult a r) H``   Theorem 8.2
``witnesses(a, H, W)``     ``undetermined``                  Corollary 8.3
=========================  ================================  ============================

A thought is any predicate ``world -> bool``; an assay is any function ``world -> result``
whose results can be compared with ``==`` (round floats first, or use a tolerance-aware key).
"""

from __future__ import annotations

from itertools import combinations
from typing import Callable, Hashable, Iterable, Sequence, TypeVar

World = TypeVar("World")
Thought = Callable[[World], bool]
Assay = Callable[[World], Hashable]


def entails(phi: Thought, psi: Thought, worlds: Iterable[World]) -> bool:
    """phi entails psi: every world where phi holds has psi (Definition 3.2)."""
    return all(psi(w) for w in worlds if phi(w))


def supported(observations: Sequence[Thought], theory: Thought, worlds: Iterable[World]) -> bool:
    """Some world satisfies the theory and every observation (Definition 3.4)."""
    return any(theory(w) and all(o(w) for o in observations) for w in worlds)


def witnesses(assay: Assay, hypothesis: Thought, worlds: Sequence[World]) -> list[tuple]:
    """Pairs of worlds with the same result but different verdicts (Corollary 8.3)."""
    return [(w, v) for w, v in combinations(worlds, 2)
            if assay(w) == assay(v) and hypothesis(w) != hypothesis(v)]


def resolves(assay: Assay, hypothesis: Thought, worlds: Sequence[World]) -> bool:
    """Worlds with the same result agree on the hypothesis (Definition 8.1)."""
    return not witnesses(assay, hypothesis, worlds)


def decides(assay: Assay, result: Hashable, hypothesis: Thought, worlds: Iterable[World]) -> bool:
    """The result entails the hypothesis or entails its negation."""
    matching = [w for w in worlds if assay(w) == result]
    return all(hypothesis(w) for w in matching) or not any(hypothesis(w) for w in matching)


def resolves_via_decides(assay: Assay, hypothesis: Thought, worlds: Sequence[World]) -> bool:
    """Theorem 8.2: resolution is equivalent to every possible result deciding the hypothesis."""
    return all(decides(assay, assay(w), hypothesis, worlds) for w in worlds)
