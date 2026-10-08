# Benchmark quality filters

The quality filters of two published PPI benchmarks, each written as one curation policy in YAML. Both are built only from the rules on the [Datasets](../usage/datasets.md#curate-one-structure) page; nothing here is benchmark-specific code. The benchmark descriptions are in `examples/ppi/benchmarks/` in the repository, and neither benchmark is a preferred target: the same rules express the others too.

The test suite runs both the library and the CLI version below and checks they still print what is shown.

## The policies

Both files are in the repository, so you can pass them straight to `pandora curate --policy` or load them in Python.

=== "PPI v2.1"

    ```yaml
    # datasets/policies/ppi-v2.1-quality.yaml
    policy_id: ppi-v2.1-quality
    policy_name: PPI v2.1 quality filters
    policy_version: 1.0.0
    description: >
      Step 1 of the PPI v2.1 benchmark: sharp X-ray and cryo-EM structures,
      EPPIC's R-factor rules (X-ray only), no non-standard residues, and
      chains with few missing residues.
    quality_rules:
      max_resolution_by_method:
        X-RAY DIFFRACTION: 2.5
        ELECTRON MICROSCOPY: 2.0
      null_resolution_behavior: exclude
      rfactor_methods: [X-RAY DIFFRACTION]
      max_r_free: 0.35
      max_r_free_gap: 0.07
      max_r_sym: 0.10
      null_rfactor_behavior: include
      exclude_nonstandard_residues: true
      missing_residue_definition: incomplete_backbone
      max_missing_tail_fraction: 0.30
      max_missing_middle_fraction: 0.10
    ```

=== "ATOM3D PIP"

    ```yaml
    # datasets/policies/atom3d-pip-quality.yaml
    policy_id: atom3d-pip-quality
    policy_name: ATOM3D PIP filters
    policy_version: 1.0.0
    description: >
      The ATOM3D PIP corpus filter: X-ray and cryo-EM entries at 3.5 A or
      better with a chain of at least 50 residues.
    quality_rules:
      max_resolution: 3.5
      include_experimental_methods: [X-RAY DIFFRACTION, ELECTRON MICROSCOPY]
      min_chain_length: 50
    ```

## Run them

Both versions map MSE to MET during canonicalisation, then curate the same eight fixtures with each policy.

=== "`library`"

    ```python
    import yaml

    from pandora.canonicalisation import canonicalise_structure
    from pandora.datasets import curate_structure
    from pandora.metadata import collect_metadata
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.canonicalisation import (
        ModifiedResidueRules,
        canonicalisationPolicy,
    )
    from pandora.schemas.dataset import DatasetCurationPolicy

    # Map MSE to MET first, so the non-standard residue rule doesn't drop
    # every selenomethionine crystal.
    mse_policy = canonicalisationPolicy(
        policy_id="mse",
        policy_name="Map MSE",
        policy_version="1.0.0",
        modified_residue_rules=ModifiedResidueRules(strategy="map_to_parent"),
    )

    entries = ["104m", "1a08", "1aui", "1b6w", "10mv", "10tm", "22jy", "1p58"]
    structures = {}
    for entry_id in entries:
        structure, _, _ = mmcif_to_structure(f"datasets/dev/mmcif/{entry_id}.cif")
        canonical, _, _ = canonicalise_structure(structure, mse_policy)
        structures[entry_id] = (canonical, collect_metadata(canonical))

    for path in [
        "datasets/policies/atom3d-pip-quality.yaml",
        "datasets/policies/ppi-v2.1-quality.yaml",
    ]:
        with open(path) as handle:
            policy = DatasetCurationPolicy.model_validate(yaml.safe_load(handle))
        print(policy.policy_name)
        for entry_id, (structure, metadata) in structures.items():
            curated, exclusions, _ = curate_structure(structure, metadata, policy)
            line = f"  {entry_id}: {'kept' if curated else 'excluded'}"
            if exclusions:
                line += " - " + "; ".join(
                    f"{e.chain_id or 'entry'} {e.reason_code}" for e in exclusions
                )
            print(line)
    ```

    ```text
    ATOM3D PIP filters
      104m: kept
      1a08: kept
      1aui: kept
      1b6w: kept
      10mv: kept
      10tm: kept
      22jy: kept
      1p58: excluded - entry RESOLUTION_THRESHOLD
    PPI v2.1 quality filters
      104m: excluded - entry RFREE_GAP_THRESHOLD
      1a08: excluded - entry NONSTANDARD_RESIDUE
      1aui: kept - A MISSING_MIDDLE
      1b6w: excluded - entry RSYM_THRESHOLD
      10mv: kept
      10tm: excluded - entry RESOLUTION_THRESHOLD
      22jy: excluded - entry RESOLUTION_THRESHOLD
      1p58: excluded - entry RESOLUTION_THRESHOLD
    ```

=== "`cli`"

    ```bash
    mkdir -p raw
    for entry in 104m 1a08 1aui 1b6w 10mv 10tm 22jy 1p58; do
      cp "datasets/dev/mmcif/$entry.cif" raw/
    done
    cat > mse.yaml <<'EOF'
    policy_id: mse
    policy_name: Map MSE
    policy_version: 1.0.0
    modified_residue_rules:
      strategy: map_to_parent
    EOF
    pandora canonicalise --input-dir raw/ --policy mse.yaml --output-dir mapped/
    for name in atom3d-pip-quality ppi-v2.1-quality; do
      pandora curate --input-dir mapped/ \
        --policy "datasets/policies/$name.yaml" --output-dir "$name/"
      jq -r '.[] | "  \(.entry_id) \(.chain_id // "entry") \(.reason_code)"' \
        "$name/curation_exclusions.json"
    done
    ```

    ```text
    canonicalised 8 structures -> mapped
    curated: 7 retained, 1 excluded -> atom3d-pip-quality
      1P58 entry RESOLUTION_THRESHOLD
    curated: 2 retained, 6 excluded -> ppi-v2.1-quality
      104M entry RFREE_GAP_THRESHOLD
      10TM entry RESOLUTION_THRESHOLD
      1A08 entry NONSTANDARD_RESIDUE
      1AUI A MISSING_MIDDLE
      1B6W entry RSYM_THRESHOLD
      1P58 entry RESOLUTION_THRESHOLD
      22JY entry RESOLUTION_THRESHOLD
    ```

    Each run writes the structures that pass to its output directory, and its exclusion records to `curation_exclusions.json`.

22jy shows the difference between the two: a 2.2 Å cryo-EM structure passes ATOM3D's 3.5 Å limit but not PPI v2.1's 2.0 Å cryo-EM limit. 1aui passes PPI v2.1 as an entry, but loses chain A, which is missing 20% of its middle.

## How the policies answer the benchmarks' open questions

The PPI v2.1 description leaves some points open. A policy has to pick an answer, and these are the ones above:

| Open point | Answer in the policy |
|---|---|
| Do the Rfree and Rsym rules apply to cryo-EM? | No: `rfactor_methods: [X-RAY DIFFRACTION]`. |
| Which mmCIF field is Rsym, and what if it's missing? | `_reflns.pdbx_Rsym_value`, falling back to `pdbx_Rmerge_I_obs`; a missing value keeps the entry (`null_rfactor_behavior: include`). |
| What happens to entries with no resolution (NMR)? | Excluded (`null_resolution_behavior: exclude`). |
| Is `MSE` a non-standard residue? | No: it is mapped to `MET` before curation. |
| What is the 30% of missing tails a percentage of? | The whole SEQRES length; the middle's 10% is of SEQRES minus the missing tails. |

## What these policies don't cover

- **PPI v2.1:** the cryo-EM FSC ≤ 0.143 rule, which as written filters nothing; PDB-REDO models for X-ray entries; tag trimming and the FASTA comparison; and the interface-gap check, which needs chain pairs.
- **ATOM3D PIP:** "at least 50 amino acids" is applied to the longest chain (`min_chain_length`), not the whole entry, and "contains protein" has no rule of its own; the PDB snapshot date isn't a curation rule either.
