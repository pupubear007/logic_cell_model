import Lake

open Lake DSL

package CellLogic where
  leanOptions := #[⟨`autoImplicit, false⟩]

-- The logic framework of "Deduction and Induction in One Diagram" (Thought, Entails, Resolves,
-- resolves_iff, modus_tollens, ...). Mathlib comes with it, at the same pinned commit.
require WangLogic from git
  "https://github.com/pupubear007/deductive-inductive-logic.git" @ "210ae35ac76b3c76dcd90e853273b68b29e32aa0" / "lean"

@[default_target] lean_lib CellLogic
