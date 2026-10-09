"""Quantitative cell and molecule models for plant–pathogen interactions.

Each module implements one topic of molecular and cellular bioengineering as small, tested
functions in SI units:

* ``binding``: receptor–ligand equilibrium and kinetics
* ``expression``: mRNA and protein dynamics
* ``enzyme``: Michaelis–Menten kinetics, inhibition, progress curves
* ``regulation``: network dynamics, host-regulated expression, autoregulation
* ``stochastic``: exact simulation of few-molecule expression
* ``identifiability``: what a measurement design can and cannot determine (assay resolution)
* ``growth``: hyphal and lesion growth
* ``transport``: diffusion with reaction (penetration length, Thiele modulus)
* ``microfluidics``: low-Reynolds-number channel flow (Stokes flow)
"""

__version__ = "0.1.0"
