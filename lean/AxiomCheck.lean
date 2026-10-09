import CellLogic

/-! Axioms used by every theorem of `CellLogic`. CI fails on `sorryAx` or any axiom other than
`propext`, `Classical.choice` and `Quot.sound`. Run with `lake env lean AxiomCheck.lean`. -/

open CellLogic

#print axioms resolves_of_injective
#print axioms not_resolves_of_witness
#print axioms snapshot_not_resolves_fast
#print axioms snapshotShutoff_injective
#print axioms snapshotShutoff_resolves
#print axioms oneHost_not_resolves
#print axioms fixedProgram_predicts_equal
#print axioms fixedProgram_refuted
#print axioms twoHosts_injective
#print axioms twoHosts_resolves
#print axioms lowSubstrate_not_resolves
#print axioms twoSubstrates_injective
#print axioms twoSubstrates_resolves
#print axioms bulk_not_resolves
#print axioms singleCell_injective
#print axioms singleCell_resolves
