import WangLogic.Plant
import Mathlib.Basic.Real.Basic

/-!
# Cell models in the logic of *Deduction and Induction in One Diagram*

The reading of the paper's symbols used throughout this package:

* A **world** is a cell mechanism: a vector of rate constants, or a choice of regulatory
  structure.
* A **thought** (`Thought W`) is a hypothesis about the cell, such as "turnover is fast" or "the
  gene is host-responsive".
* An **assay** (`W → R`) is an experimental design: the map from mechanism to what is measured.
  Its result is the thought `AssayResult a r`.
* **Deduction** (`Entails`): a mechanism hypothesis predicts a measurement. **Falsification**
  (`modus_tollens`): a failed prediction eliminates the hypothesis.
* **Resolution** (`Resolves`, Theorem 8.2 `resolves_iff`): an experiment can decide a
  hypothesis exactly when every possible result entails it or rules it out. Two mechanisms that
  give the same result but disagree on the hypothesis (`undetermined`, Corollary 8.3) show that
  no amount of that experiment decides it. In modelling this is called (non-)identifiability.

This file adds two general tools. Each topic file then proves, for one design, that it does not
resolve a hypothesis (with an explicit witness pair) and, for a better design, that it does.
-/

namespace CellLogic

open WangLogic

/-- An experiment that tells all mechanisms apart resolves every hypothesis. -/
theorem resolves_of_injective {W R : Type} {a : W → R} (h : Function.Injective a)
    (H : Thought W) : Resolves a H :=
  fun w w' e => by rw [h e]

/-- A witness pair (same result, different verdict on `H`) shows that `a` does not resolve `H`.
This is Corollary 8.3 (`undetermined`) packaged as a non-resolution statement. -/
theorem not_resolves_of_witness {W R : Type} {a : W → R} {H : Thought W} {w w' : W}
    (hsame : a w = a w') (hH : H w) (hnH : ¬H w') : ¬Resolves a H :=
  fun hres => hnH ((hres w w' hsame).mp hH)

end CellLogic
