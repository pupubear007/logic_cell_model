"""What a measurement design can and cannot determine.

This is assay resolution (Definition 8.1 and Theorem 8.2 of *Deduction and Induction in One
Diagram*) applied to model parameters. Worlds are parameter vectors theta; the "assay" is the
measurement design, the map theta -> predicted observations. Two parameter vectors that give the
same predictions cannot be told apart by any amount of that data (Corollary 8.3), so a
conclusion that depends on which of them is true is not resolved.

Two tools:

* ``witness``: given two parameter vectors, decide whether the design separates them.
* ``local_identifiability``: near a parameter vector, find which parameter combinations the
  design determines (rank of the sensitivity matrix) and which it leaves free (null directions).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np

Predict = Callable[[np.ndarray], np.ndarray]


def witness(predict: Predict, theta_a: Sequence[float], theta_b: Sequence[float],
            rtol: float = 1e-6, atol: float = 1e-12) -> bool:
    """True if theta_a and theta_b give the same predictions under the design.

    A True result is a witness of underdetermination (Corollary 8.3): the design cannot resolve
    any conclusion on which theta_a and theta_b disagree.
    """
    ya = np.asarray(predict(np.asarray(theta_a, float)), float)
    yb = np.asarray(predict(np.asarray(theta_b, float)), float)
    return bool(np.allclose(ya, yb, rtol=rtol, atol=atol))


@dataclass
class LocalIdentifiability:
    rank: int
    n_params: int
    singular_values: np.ndarray
    null_directions: np.ndarray  # rows: unresolved directions in log-parameter space
    names: tuple[str, ...]

    @property
    def identifiable(self) -> bool:
        return self.rank == self.n_params

    def describe(self) -> str:
        lines = [f"rank {self.rank} of {self.n_params} parameters"
                 + (" (all identifiable)" if self.identifiable else "")]
        for v in self.null_directions:
            v = v / np.max(np.abs(v))
            terms = " ".join(f"{c:+.2f}*log({n})" for c, n in zip(v, self.names) if abs(c) > 1e-6)
            lines.append(f"  unresolved direction: {terms}")
        return "\n".join(lines)


def local_identifiability(predict: Predict, theta: Sequence[float],
                          names: Sequence[str] | None = None, rel_step: float = 1e-6,
                          rel_tol: float = 1e-6) -> LocalIdentifiability:
    """Sensitivity analysis in log-parameters at theta (all parameters must be positive).

    The columns of the sensitivity matrix are d(prediction)/d(log theta_i). A singular value
    below rel_tol times the largest one counts as zero; the matching right singular vectors are
    parameter combinations the design does not determine.
    """
    theta = np.asarray(theta, float)
    if np.any(theta <= 0):
        raise ValueError("parameters must be positive (sensitivities are taken in log space)")
    names = tuple(names) if names is not None else tuple(f"theta{i}" for i in range(len(theta)))
    base = np.asarray(predict(theta), float).ravel()
    S = np.empty((base.size, theta.size))
    for i in range(theta.size):
        up, dn = theta.copy(), theta.copy()
        up[i] *= np.exp(rel_step)
        dn[i] *= np.exp(-rel_step)
        S[:, i] = (np.asarray(predict(up), float).ravel()
                   - np.asarray(predict(dn), float).ravel()) / (2 * rel_step)
    _, s, Vt = np.linalg.svd(S, full_matrices=True)
    s_full = np.zeros(theta.size)
    s_full[:s.size] = s
    cutoff = rel_tol * (s_full.max() if s_full.max() > 0 else 1.0)
    rank = int(np.sum(s_full > cutoff))
    return LocalIdentifiability(rank, theta.size, s_full, Vt[rank:], names)
