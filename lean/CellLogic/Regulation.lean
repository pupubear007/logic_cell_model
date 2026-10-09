import CellLogic.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.FieldSimp

/-!
# Fixed program or host-responsive regulation

Transcription of a pathogen gene per unit of pathogen is `basal + induced * s`, where `s ∈ [0,1]`
is the host signal seen through the gene's response (the Hill factor `H(signal)`). The two
hypotheses about aggressiveness determinants are:

* `FixedProgram`: `induced = 0`, the same expression on every host;
* `HostResponsive`: `induced ≠ 0`.

Results:

* `oneHost_not_resolves`: expression on a single host never decides between them.
* `fixedProgram_predicts_equal` and `fixedProgram_refuted`: a fixed program predicts equal
  expression on two hosts (deduction), so a measured difference eliminates it (modus tollens).
* `twoHosts_resolves`: two hosts with different signals decide every hypothesis about the
  mechanism, in particular fixed program versus host-responsive.
-/

namespace CellLogic

open WangLogic

/-- A world: basal and host-induced transcription rates (per unit pathogen). -/
structure Regulation where
  basal : ℝ
  induced : ℝ

/-- Assay: expression per unit pathogen on a host whose signal factor is `s`. -/
def oneHost (s : ℝ) (w : Regulation) : ℝ := w.basal + w.induced * s

/-- Assay: expression on two hosts with signal factors `s₁` and `s₂`. -/
def twoHosts (s₁ s₂ : ℝ) (w : Regulation) : ℝ × ℝ := (oneHost s₁ w, oneHost s₂ w)

def FixedProgram : Thought Regulation := fun w => w.induced = 0

def HostResponsive : Thought Regulation := fun w => w.induced ≠ 0

/-- On one host, a host-responsive gene can look exactly like a fixed program. -/
theorem oneHost_not_resolves (s : ℝ) : ¬Resolves (oneHost s) HostResponsive :=
  not_resolves_of_witness (w := ⟨0, 1⟩) (w' := ⟨s, 0⟩) (by simp [oneHost])
    (by simp [HostResponsive]) (by simp [HostResponsive])

/-- Deduction: a fixed program predicts the same expression on any two hosts. -/
theorem fixedProgram_predicts_equal (s₁ s₂ : ℝ) :
    Entails FixedProgram fun w => oneHost s₁ w = oneHost s₂ w := by
  intro w hw
  simp only [FixedProgram] at hw
  simp [oneHost, hw]

/-- Falsification: a measured difference between two hosts eliminates the fixed program. -/
theorem fixedProgram_refuted {s₁ s₂ : ℝ} {w : Regulation}
    (hdiff : Polarity.down.apply (oneHost s₁ w = oneHost s₂ w)) :
    Polarity.down.apply (FixedProgram w) :=
  modus_tollens (fixedProgram_predicts_equal s₁ s₂) hdiff

/-- Two hosts with different signal factors tell all mechanisms apart. -/
theorem twoHosts_injective {s₁ s₂ : ℝ} (hs : s₁ ≠ s₂) : Function.Injective (twoHosts s₁ s₂) := by
  rintro ⟨b, i⟩ ⟨b', i'⟩ h
  simp only [twoHosts, oneHost, Prod.mk.injEq] at h
  obtain ⟨h1, h2⟩ := h
  have hi : (i - i') * (s₁ - s₂) = 0 := by linarith
  have hii : i = i' := by
    rcases mul_eq_zero.mp hi with h | h
    · linarith
    · exact absurd (sub_eq_zero.mp h) hs
  subst hii
  have hbb : b = b' := by linarith
  subst hbb
  rfl

/-- Two hosts with different signal factors decide every hypothesis, including whether the gene
is a fixed program or host-responsive. -/
theorem twoHosts_resolves {s₁ s₂ : ℝ} (hs : s₁ ≠ s₂) (H : Thought Regulation) :
    Resolves (twoHosts s₁ s₂) H :=
  resolves_of_injective (twoHosts_injective hs) H

end CellLogic
