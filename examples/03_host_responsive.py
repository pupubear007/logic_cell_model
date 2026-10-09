"""Fixed program or host-responsive regulation: which experiment can tell them apart?

Two candidate explanations of a pathogen gene's expression:
  A. basal transcription only (alpha0 = 3, alpha1 = 0): a fixed program;
  B. lower basal transcription plus induction by a host signal (alpha0 = 1.5, alpha1 = 3).
On a host whose signal gives H(s) = 1/2 they predict the same expression at every time point, so
no amount of data from that host can separate them (Corollary 8.3 of the logic paper). A second
host with a different signal level does. All numbers are illustrative.

Run:  python examples/03_host_responsive.py  ->  examples/out/03_host_responsive.png
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cellmodels import regulation
from cellmodels.identifiability import local_identifiability, witness

T = np.linspace(0.0, 96.0, 300)        # hours post inoculation
samples = np.array([24.0, 48.0, 96.0])
delta = 0.05                           # 1/h, illustrative
hosts = {"host 1 (signal 1)": 1.0, "host 2 (signal 5)": 5.0}
A = dict(alpha0=3.0, alpha1=0.0)
B = dict(alpha0=1.5, alpha1=3.0)


def predict(hyp, signal, t):
    return regulation.regulated_expression(t, hyp["alpha0"], hyp["alpha1"], delta, signal)


one_host = lambda th: regulation.regulated_expression(samples, th[0], th[1], delta, 1.0)
both_hosts = lambda th: regulation.across_hosts(samples, th[0], th[1], delta, list(hosts.values()))

print("Host 1 alone separates A from B:", not witness(one_host, [3.0, 1e-9], [1.5, 3.0]))
print("Hosts 1 and 2 separate A from B:", not witness(both_hosts, [3.0, 1e-9], [1.5, 3.0]))
print("\nHost 1 only:\n" + local_identifiability(one_host, [1.5, 3.0], ("alpha0", "alpha1")).describe())
print("\nHosts 1 and 2:\n" + local_identifiability(both_hosts, [1.5, 3.0], ("alpha0", "alpha1")).describe())

fig, axes = plt.subplots(1, 2, figsize=(8, 3.4), sharey=True)
for ax, (name, s) in zip(axes, hosts.items()):
    for hyp, label, color, ls in [(A, "A: fixed program", "#2a6f97", "-"),
                                  (B, "B: host-responsive", "#d17a22", "--")]:
        ax.plot(T, predict(hyp, s, T), color=color, ls=ls, label=label)
        ax.plot(samples, predict(hyp, s, samples), "o", color=color, ms=4)
    ax.set_title(name, fontsize=10)
    ax.set_xlabel("hours post inoculation")
axes[0].set_ylabel("expression (relative)")
axes[0].legend(fontsize=8, frameon=False, loc="lower right")
fig.suptitle("Same data on host 1, different on host 2", fontsize=10)
fig.tight_layout()
out = Path(__file__).parent / "out"
out.mkdir(exist_ok=True)
fig.savefig(out / "03_host_responsive.png", dpi=200)
print(f"\nfigure: {out / '03_host_responsive.png'}")
