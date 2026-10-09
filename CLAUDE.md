# Conventions for this repository

- Public repository: **never commit experimental data, lab results, or unpublished thesis
  material.** Real data are read from local paths given in untracked config files.
- One topic per module in `cellmodels/`, in SI units unless a docstring says otherwise.
  Each function documents its equation.
- Every model gets tests against a closed form, a limit, or a conservation law
  (`tests/test_models.py`).
- The logic design comes first. A hypothesis about a cell is a `Thought` over worlds
  (mechanisms), an experiment is an assay, and claims about what an experiment can decide are
  `Resolves` statements, proved in `lean/CellLogic` on top of `WangLogic` and mirrored in
  `tests/test_logic.py`. Every new model gets one: what the intended experiment cannot decide
  (witness pair) and what design does decide it.
- Lean: no `sorry`; only `propext`, `Classical.choice`, `Quot.sound` (`lean/AxiomCheck.lean`).
  Do not change a theorem statement to make a proof go through without saying so.
- Examples use illustrative parameter values and say so. Do not present them as measurements.
- Cite the textbook by topic (binding, enzyme kinetics, gene expression and trafficking, network
  dynamics, growth, transport and reaction, stochastic processes), not by quoting it.
