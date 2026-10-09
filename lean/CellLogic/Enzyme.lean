import CellLogic.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.LinearCombination

/-!
# Enzyme kinetics: what initial rates can decide

Michaelis–Menten rate `v(S) = Vmax S / (Km + S)`. At substrate concentrations far below `Km` the
rate is `(Vmax / Km) S`, so such measurements see only the specificity constant `Vmax / Km`.

* `lowSubstrate_not_resolves`: the low-substrate slope does not decide whether `Km ≤ c`.
* `twoSubstrates_resolves`: exact rates at two different substrate concentrations decide every
  hypothesis about `(Vmax, Km)`.
-/

namespace CellLogic

open WangLogic

/-- A world: positive `Vmax` and `Km` of one enzyme. -/
structure MichaelisMenten where
  vmax : ℝ
  km : ℝ
  vmax_pos : 0 < vmax
  km_pos : 0 < km

/-- Initial rate at substrate concentration `S`. -/
noncomputable def rate (S : ℝ) (w : MichaelisMenten) : ℝ := w.vmax * S / (w.km + S)

/-- Assay: the slope `Vmax / Km` of rate against substrate at low substrate. -/
noncomputable def lowSubstrateSlope (w : MichaelisMenten) : ℝ := w.vmax / w.km

/-- Assay: rates at two substrate concentrations. -/
noncomputable def twoSubstrates (S₁ S₂ : ℝ) (w : MichaelisMenten) : ℝ × ℝ :=
  (rate S₁ w, rate S₂ w)

/-- Hypothesis: the enzyme is half-saturated at or below concentration `c`. -/
def HighAffinity (c : ℝ) : Thought MichaelisMenten := fun w => w.km ≤ c

/-- Low-substrate rates cannot decide the enzyme's affinity. -/
theorem lowSubstrate_not_resolves : ¬Resolves lowSubstrateSlope (HighAffinity 1) := by
  refine not_resolves_of_witness (w := ⟨1, 1, by norm_num, by norm_num⟩)
    (w' := ⟨2, 2, by norm_num, by norm_num⟩) ?_ ?_ ?_
  · norm_num [lowSubstrateSlope]
  · norm_num [HighAffinity]
  · norm_num [HighAffinity]

/-- Equal rates at a positive substrate concentration, cross-multiplied. -/
private lemma rate_eq_iff {S : ℝ} (hS : 0 < S) (w w' : MichaelisMenten) :
    rate S w = rate S w' ↔ w.vmax * (w'.km + S) = w'.vmax * (w.km + S) := by
  have h1 : 0 < w.km + S := by linarith [w.km_pos]
  have h2 : 0 < w'.km + S := by linarith [w'.km_pos]
  unfold rate
  rw [div_eq_div_iff h1.ne' h2.ne']
  constructor
  · intro h
    have : S * (w.vmax * (w'.km + S)) = S * (w'.vmax * (w.km + S)) := by linarith
    exact mul_left_cancel₀ hS.ne' this
  · intro h
    linear_combination S * h

/-- Rates at two different positive substrate concentrations tell all enzymes apart. -/
theorem twoSubstrates_injective {S₁ S₂ : ℝ} (h1 : 0 < S₁) (h2 : 0 < S₂) (hS : S₁ ≠ S₂) :
    Function.Injective (twoSubstrates S₁ S₂) := by
  intro w w' h
  simp only [twoSubstrates, Prod.mk.injEq] at h
  have e1 := (rate_eq_iff h1 w w').mp h.1
  have e2 := (rate_eq_iff h2 w w').mp h.2
  have hv : (w.vmax - w'.vmax) * (S₁ - S₂) = 0 := by linear_combination e1 - e2
  have hvv : w.vmax = w'.vmax := by
    rcases mul_eq_zero.mp hv with h | h
    · linarith
    · exact absurd (sub_eq_zero.mp h) hS
  have hkk : w.km = w'.km := by
    have : w.vmax * (w'.km - w.km) = 0 := by rw [← hvv] at e1; linear_combination e1
    rcases mul_eq_zero.mp this with h | h
    · exact absurd h w.vmax_pos.ne'
    · linarith
  cases w; cases w'
  simp_all

/-- Rates at two different substrate concentrations decide every hypothesis about the enzyme. -/
theorem twoSubstrates_resolves {S₁ S₂ : ℝ} (h1 : 0 < S₁) (h2 : 0 < S₂) (hS : S₁ ≠ S₂)
    (H : Thought MichaelisMenten) : Resolves (twoSubstrates S₁ S₂) H :=
  resolves_of_injective (twoSubstrates_injective h1 h2 hS) H

end CellLogic
