"""How far does a secreted acid reach into tissue before it is neutralized?

Steady diffusion with first-order consumption (buffering) into a leaf of thickness L. The
penetration length sqrt(D/k) sets the scale; the Thiele modulus L sqrt(k/D) says whether the
whole leaf is exposed (phi << 1) or only a surface layer (phi >> 1). D is the order of magnitude
for a small molecule in water; the consumption rates and thickness are illustrative.

Run:  python examples/02_acid_penetration.py  ->  examples/out/02_acid_penetration.png
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cellmodels import transport

D = 1e-9          # m^2/s, small molecule in water (order of magnitude)
L = 200e-6        # m, leaf thickness (illustrative)
x = np.linspace(0, L, 400)

fig, ax = plt.subplots(figsize=(6, 3.6))
for k, color in [(1e-3, "#2a6f97"), (1e-1, "#5a9e6f"), (10.0, "#d17a22")]:
    phi = transport.thiele_modulus(L, D, k)
    lam = transport.penetration_length(D, k)
    ax.plot(x * 1e6, transport.slab_profile(x, 1.0, L, D, k), color=color,
            label=f"k={k:g}/s: lambda={lam * 1e6:.0f} um, phi={phi:.2g}")
ax.set_xlabel("depth into leaf (um)")
ax.set_ylabel("concentration / source")
ax.set_title("Reaction–diffusion: penetration of a secreted acid", fontsize=10)
ax.legend(fontsize=8, frameon=False)
fig.tight_layout()
out = Path(__file__).parent / "out"
out.mkdir(exist_ok=True)
fig.savefig(out / "02_acid_penetration.png", dpi=200)
print(f"figure: {out / '02_acid_penetration.png'}")
