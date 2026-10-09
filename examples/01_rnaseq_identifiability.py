"""What can RNA-seq tell about synthesis and decay rates?

Two genes with different rates but the same ratio look identical in a steady-state snapshot,
and different after induction. The snapshot does not resolve the rates (Corollary 8.3 of the
logic paper); the time course does. All numbers are illustrative.

Run:  python examples/01_rnaseq_identifiability.py  ->  examples/out/01_rnaseq_identifiability.png
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cellmodels import expression
from cellmodels.identifiability import local_identifiability, witness

slow = (2.0, 0.25)   # alpha (mRNA/h), delta (1/h): half-life ~2.8 h
fast = (8.0, 1.0)    # four times faster synthesis and decay, same ratio
sample_times = np.array([1.0, 2.0, 4.0, 8.0])


def snapshot(theta):
    return np.array([expression.steady_state(*theta)])


def time_course(theta):
    return expression.mrna(sample_times, *theta)


print("Snapshot gives the same value for both genes:", witness(snapshot, slow, fast))
print("Time course gives the same values:          ", witness(time_course, slow, fast))
print("\nSnapshot design:\n" + local_identifiability(snapshot, slow, ("alpha", "delta")).describe())
print("\nTime-course design:\n" + local_identifiability(time_course, slow, ("alpha", "delta")).describe())

t = np.linspace(0, 12, 300)
fig, ax = plt.subplots(figsize=(6, 3.6))
for (a, d), label, color in [(slow, "slow turnover", "#2a6f97"), (fast, "fast turnover", "#d17a22")]:
    ax.plot(t, expression.mrna(t, a, d), color=color, label=f"{label} (alpha={a}, delta={d})")
ax.plot(sample_times, expression.mrna(sample_times, *slow), "o", color="#2a6f97")
ax.plot(sample_times, expression.mrna(sample_times, *fast), "o", color="#d17a22")
ax.axhline(expression.steady_state(*slow), color="grey", ls=":", lw=1)
ax.text(12, expression.steady_state(*slow) + 0.2, "steady-state snapshot: identical for both",
        fontsize=8, color="grey", ha="right", va="bottom")
ax.set_ylim(top=9.2)
ax.set_xlabel("time after induction (h)")
ax.set_ylabel("mRNA (relative)")
ax.set_title("Same steady state, different dynamics", fontsize=10)
ax.legend(fontsize=8, frameon=False, loc="lower right")
fig.tight_layout()
out = Path(__file__).parent / "out"
out.mkdir(exist_ok=True)
fig.savefig(out / "01_rnaseq_identifiability.png", dpi=200)
print(f"\nfigure: {out / '01_rnaseq_identifiability.png'}")
