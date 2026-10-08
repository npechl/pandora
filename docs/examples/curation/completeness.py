from pandora.datasets import chain_completeness
from pandora.parsing import mmcif_to_structure

structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1aui.cif")
records, diagnostics = chain_completeness(structure)
for r in records:
    print(
        f"{r.chain_id}: SEQRES {r.seqres_length}, missing "
        f"{r.missing_n_term} at the N-terminus, {r.missing_c_term} at the "
        f"C-terminus, {r.missing_middle} in the middle"
    )

# 1p58 is a CA-only trace: every residue lacks most of its backbone.
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1p58.cif")
for definition in ["unobserved", "incomplete_backbone"]:
    records, _ = chain_completeness(structure, definition)
    a = records[0]
    print(
        f"1p58 chain A, {definition}: missing {a.missing_n_term} + "
        f"{a.missing_c_term} tails, {a.missing_middle} middle "
        f"of {a.seqres_length}"
    )
