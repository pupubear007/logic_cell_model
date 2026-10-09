# Conventions for this repository

- Public repository: **never commit experimental data, lab results, or unpublished thesis
  material.** Real data are read from local paths given in untracked config files.
- One topic per module in `cellmodels/`, in SI units unless a docstring says otherwise.
  Each function documents its equation.
- Every model gets tests against a closed form, a limit, or a conservation law
  (`tests/test_models.py`).
- Every new model gets an identifiability check (`cellmodels.identifiability`): which parameters
  the intended measurement design determines and which it leaves free.
- Examples use illustrative parameter values and say so. Do not present them as measurements.
- Cite the textbook by topic (binding, enzyme kinetics, gene expression and trafficking, network
  dynamics, growth, transport and reaction, stochastic processes), not by quoting it.
