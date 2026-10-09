import CellLogic.Basic
import Mathlib.Analysis.SpecialFunctions.Exp

/-!
# Gene expression: what RNA-seq can decide about synthesis and decay

mRNA follows `dm/dt = α - δ m`. Its steady state is `α / δ`; after transcription is blocked it
decays as `m₀ exp(-δ t)`.

* `snapshot_not_resolves_fast`: a steady-state snapshot does not decide "turnover is fast",
  because doubling both rates leaves the snapshot unchanged.
* `snapshot_shutoff_resolves`: a snapshot together with a transcription-shutoff measurement
  separates every pair of mechanisms, so it decides every hypothesis about the rates.
-/

namespace CellLogic

open WangLogic

/-- A world: positive synthesis rate `α` and decay rate `δ` of one mRNA. -/
structure Turnover where
  synth : ℝ
  decay : ℝ
  synth_pos : 0 < synth
  decay_pos : 0 < decay

/-- Assay: steady-state abundance `α / δ` (a single RNA-seq snapshot). -/
noncomputable def snapshot (w : Turnover) : ℝ := w.synth / w.decay

/-- Assay: fraction of mRNA left a time `t` after transcription is blocked, `exp(-δ t)`. -/
noncomputable def shutoff (t : ℝ) (w : Turnover) : ℝ := Real.exp (-(w.decay * t))

/-- Assay: snapshot and shutoff together. -/
noncomputable def snapshotShutoff (t : ℝ) (w : Turnover) : ℝ × ℝ := (snapshot w, shutoff t w)

/-- Hypothesis: the mRNA turns over faster than rate `c`. -/
def FastTurnover (c : ℝ) : Thought Turnover := fun w => c < w.decay

/-- A steady-state snapshot cannot decide whether turnover is fast. -/
theorem snapshot_not_resolves_fast : ¬Resolves snapshot (FastTurnover 1) := by
  refine not_resolves_of_witness (w := ⟨4, 2, by norm_num, by norm_num⟩)
    (w' := ⟨2, 1, by norm_num, by norm_num⟩) ?_ ?_ ?_
  · norm_num [snapshot]
  · norm_num [FastTurnover]
  · norm_num [FastTurnover]

/-- Snapshot plus shutoff tells all mechanisms apart (for any shutoff time `t ≠ 0`). -/
theorem snapshotShutoff_injective {t : ℝ} (ht : t ≠ 0) :
    Function.Injective (snapshotShutoff t) := by
  rintro ⟨a, d, ha, hd⟩ ⟨a', d', ha', hd'⟩ h
  simp only [snapshotShutoff, snapshot, shutoff, Prod.mk.injEq] at h
  obtain ⟨h1, h2⟩ := h
  have hdd : d = d' := mul_right_cancel₀ ht (neg_inj.mp (Real.exp_injective h2))
  subst hdd
  have haa : a = a' := (div_left_inj' hd.ne').mp h1
  subst haa
  rfl

/-- Snapshot plus shutoff decides every hypothesis about synthesis and decay. -/
theorem snapshotShutoff_resolves {t : ℝ} (ht : t ≠ 0) (H : Thought Turnover) :
    Resolves (snapshotShutoff t) H :=
  resolves_of_injective (snapshotShutoff_injective ht) H

end CellLogic
