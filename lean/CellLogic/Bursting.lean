import CellLogic.Basic
import Mathlib.Tactic.Linarith

/-!
# Stochastic expression: bulk versus single-cell data

With bursts arriving at frequency `a` and geometric burst size of mean `b` (decay rate 1), the
stationary mean is `a b` and the Fano factor (variance / mean) is `1 + b`.

* `bulk_not_resolves`: bulk (population-average) data, which see only the mean, cannot decide
  whether bursts are large.
* `singleCell_resolves`: single-cell data, which give the mean and the Fano factor, decide every
  hypothesis about burst frequency and size.
-/

namespace CellLogic

open WangLogic

/-- A world: positive burst frequency and mean burst size. -/
structure Bursting where
  freq : ℝ
  size : ℝ
  freq_pos : 0 < freq
  size_pos : 0 < size

/-- Assay: the population mean (bulk RNA-seq). -/
def bulkMean (w : Bursting) : ℝ := w.freq * w.size

/-- Assay: mean and Fano factor (single-cell counts). -/
def singleCell (w : Bursting) : ℝ × ℝ := (w.freq * w.size, 1 + w.size)

/-- Hypothesis: expression comes in large bursts, mean size above `c`. -/
def LargeBursts (c : ℝ) : Thought Bursting := fun w => c < w.size

/-- Bulk data cannot decide whether expression is bursty. -/
theorem bulk_not_resolves : ¬Resolves bulkMean (LargeBursts 1) := by
  refine not_resolves_of_witness (w := ⟨1, 2, by norm_num, by norm_num⟩)
    (w' := ⟨2, 1, by norm_num, by norm_num⟩) ?_ ?_ ?_
  · norm_num [bulkMean]
  · norm_num [LargeBursts]
  · norm_num [LargeBursts]

theorem singleCell_injective : Function.Injective singleCell := by
  rintro ⟨a, b, ha, hb⟩ ⟨a', b', ha', hb'⟩ h
  simp only [singleCell, Prod.mk.injEq] at h
  obtain ⟨h1, h2⟩ := h
  have hbb : b = b' := by linarith
  subst hbb
  have haa : a = a' := mul_right_cancel₀ hb.ne' h1
  subst haa
  rfl

/-- Single-cell data decide every hypothesis about burst frequency and size. -/
theorem singleCell_resolves (H : Thought Bursting) : Resolves singleCell H :=
  resolves_of_injective singleCell_injective H

end CellLogic
