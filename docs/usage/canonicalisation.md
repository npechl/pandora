# Canonicalisation

`pandora.canonicalisation` normalizes a parsed `Structure` according to a policy. The steps always run in this order: modified residues, chain IDs, assemblies, missing data, altlocs, residue numbering, entities, ligands, then validation. Residue numbering runs after the missing-data checks so `renumber` can't hide sequence gaps from them. The [policy reference](../reference/policies.md#canonicalisation) lists every field and value; [Architecture](../contributing/architecture.md#canonicalisation-policy) shows how the rule groups relate, and [Functions](../reference/functions.md#pandora.canonicalisation) has full signatures.

Each section changes one rule group on fixtures from `datasets/dev/mmcif/`. The `library` tab is the Python code and what it prints; the `cli` tab does the same with `pandora canonicalise`, writing the policy as YAML and reading the results with [`jq`](https://jqlang.org/). Run the CLI examples from the repository root. The test suite runs every example on this page and checks it still prints what is shown.

The CLI writes three things to `--output-dir`: the canonical mmCIF files, `canonicalisation_provenance.json` (what the library returns as `provenance`) and `canonicalisation_mappings.json` (what it returns as `mappings`), each keyed by entry ID.

## Canonicalise with the default policy

An empty `canonicalisationPolicy` isn't a no-op: several rule groups default to something other than `"preserve"` (missing-atom/residue flagging, best-occupancy altloc selection), so even the default policy records transforms.

=== "`library`"

    ```python
    from pandora.canonicalisation import canonicalise_structure
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.canonicalisation import canonicalisationPolicy

    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")

    policy = canonicalisationPolicy(
        policy_id="p1", policy_name="Default", policy_version="1.0.0"
    )
    canonical, mappings, provenance = canonicalise_structure(structure, policy)

    print(provenance.transforms)
    ```

    ```text
    ['missing_atoms:annotate', 'missing_residues:annotate', 'altloc:select_best_occupancy']
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/104m.cif raw/
    cat > policy.yaml <<'EOF'
    policy_id: p1
    policy_name: Default
    policy_version: 1.0.0
    EOF
    pandora canonicalise --input-dir raw/ --policy policy.yaml --output-dir canonical/
    jq -c '.["104M"].transforms' canonical/canonicalisation_provenance.json
    ```

    ```text
    canonicalised 1 structures -> canonical
    ["missing_atoms:annotate","missing_residues:annotate","altloc:select_best_occupancy"]
    ```

## Choose chain IDs

`identifier_rules.chain_id` decides what each chain is called. `remap` relabels chains `A`, `B`, `C`, ... in the order they appear. `use_auth_chain_id` uses the author's chain IDs; when several asym units share one author ID (a protein and its ligands usually do), the later ones get a suffix.

=== "`library`"

    ```python
    from pandora.canonicalisation import canonicalise_structure
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.canonicalisation import (
        ChainIdRules,
        IdentifierRules,
        canonicalisationPolicy,
    )


    def chain_map(path: str, strategy: str) -> str:
        structure, _, _ = mmcif_to_structure(path)
        policy = canonicalisationPolicy(
            policy_id="p",
            policy_name="p",
            policy_version="1.0.0",
            identifier_rules=IdentifierRules(
                chain_id=ChainIdRules(strategy=strategy)
            ),
        )
        _, mappings, _ = canonicalise_structure(structure, policy)
        return ", ".join(
            f"{item.original_chain_id}->{item.canonical_chain_id}"
            for item in mappings.chain_id_mapping.items
        )


    # 13dg lists its chains as A, E, B, F, ...; remap relabels them in order.
    print("remap:", chain_map("datasets/dev/mmcif/13dg.cif", "remap"))
    # 1a08: label chain C is author chain B, and ligands share author ids.
    print(
        "use_auth_chain_id:",
        chain_map("datasets/dev/mmcif/1a08.cif", "use_auth_chain_id"),
    )
    ```

    ```text
    remap: A->A, E->B, B->C, F->D, C->E, G->F, D->G, H->H
    use_auth_chain_id: A->A, C->B, B->C, D->D, E->A_2, F->C_2, G->B_2, H->D_2
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/13dg.cif datasets/dev/mmcif/1a08.cif raw/
    cat > remap.yaml <<'EOF'
    policy_id: p
    policy_name: p
    policy_version: 1.0.0
    identifier_rules:
      chain_id:
        strategy: remap
    EOF
    cat > auth.yaml <<'EOF'
    policy_id: p
    policy_name: p
    policy_version: 1.0.0
    identifier_rules:
      chain_id:
        strategy: use_auth_chain_id
    EOF
    pandora canonicalise --input-dir raw/ --policy remap.yaml --output-dir remap/
    pandora canonicalise --input-dir raw/ --policy auth.yaml --output-dir auth/
    jq -r '.["13DG"].chain_id_mapping.items[]
      | "\(.original_chain_id)->\(.canonical_chain_id)"' \
      remap/canonicalisation_mappings.json | paste -sd, -
    jq -r '.["1A08"].chain_id_mapping.items[]
      | "\(.original_chain_id)->\(.canonical_chain_id)"' \
      auth/canonicalisation_mappings.json | paste -sd, -
    ```

    ```text
    canonicalised 2 structures -> remap
    canonicalised 2 structures -> auth
    A->A,E->B,B->C,F->D,C->E,G->F,D->G,H->H
    A->A,C->B,B->C,D->D,E->A_2,F->C_2,G->B_2,H->D_2
    ```

## Renumber residues

`identifier_rules.residue_numbering` decides each residue's `seq_id`. `preserve` keeps the SEQRES position, `use_auth_seq` uses the author's numbering, and `renumber` counts observed residues from 1, closing gaps. Every change is recorded in the residue-number mapping, so you can always get back to the original numbers.

=== "`library`"

    ```python
    from pandora.canonicalisation import canonicalise_structure
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.canonicalisation import (
        IdentifierRules,
        ResidueNumberingRules,
        canonicalisationPolicy,
    )

    # 10mv: chain A starts at SEQRES position 22 and has a gap in the middle.
    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/10mv.cif")

    for strategy in ["preserve", "use_auth_seq", "renumber"]:
        policy = canonicalisationPolicy(
            policy_id="p",
            policy_name="p",
            policy_version="1.0.0",
            identifier_rules=IdentifierRules(
                residue_numbering=ResidueNumberingRules(strategy=strategy)
            ),
        )
        canonical, mappings, _ = canonicalise_structure(structure, policy)
        seq_ids = sorted(
            {a.label_seq_id for a in canonical.atoms if a.label_asym_id == "A"}
        )
        print(f"{strategy}: first {seq_ids[:3]}, last {seq_ids[-1]}")

    item = mappings.residue_number_mapping.items[0]
    print(
        f"renumber mapping: chain {item.original_chain_id} "
        f"{item.original_seq_id} (auth {item.original_auth_seq_id}) "
        f"-> {item.canonical_seq_id}"
    )
    ```

    ```text
    preserve: first [22, 23, 24], last 257
    use_auth_seq: first [391, 392, 393], last 626
    renumber: first [1, 2, 3], last 231
    renumber mapping: chain A 22 (auth 391) -> 1
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/10mv.cif raw/
    cat > auth.yaml <<'EOF'
    policy_id: p
    policy_name: p
    policy_version: 1.0.0
    identifier_rules:
      residue_numbering:
        strategy: use_auth_seq
    EOF
    cat > renumber.yaml <<'EOF'
    policy_id: p
    policy_name: p
    policy_version: 1.0.0
    identifier_rules:
      residue_numbering:
        strategy: renumber
    EOF
    pandora canonicalise --input-dir raw/ --policy auth.yaml --output-dir auth/
    pandora canonicalise --input-dir raw/ --policy renumber.yaml --output-dir renumber/
    for run in auth renumber; do
      jq -r --arg run "$run" '.["10MV"].residue_number_mapping.items[0]
        | "\($run): chain \(.original_chain_id) \(.original_seq_id)"
          + " (auth \(.original_auth_seq_id)) -> \(.canonical_seq_id)"' \
        "$run/canonicalisation_mappings.json"
    done
    ```

    ```text
    canonicalised 1 structures -> auth
    canonicalised 1 structures -> renumber
    auth: chain A 22 (auth 391) -> 391
    renumber: chain A 22 (auth 391) -> 1
    ```

!!! warning "Renumbering breaks the link to SEQRES"
    After `use_auth_seq` or `renumber`, `seq_id` no longer indexes `_entity_poly_seq`. [`chain_completeness()`](datasets.md#inspect-chain-completeness) then reports such chains as `SEQRES_MISMATCH` instead of measuring them, and curation's missing-residue rules drop them. Keep `preserve` if you use those rules.

## Pick alternate locations

`altloc_rules` decides which conformer survives when a residue has several. `select_best_occupancy` keeps the most occupied one and breaks ties with `tie_breaker`; `select_user_defined` keeps the altloc you name. Each choice is recorded in the altloc selection mapping.

=== "`library`"

    ```python
    from pandora.canonicalisation import canonicalise_structure
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.canonicalisation import (
        AltlocRules,
        canonicalisationPolicy,
    )

    # 1b6w: residue 35 has two conformers, A and B, at occupancy 0.5 each.
    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1b6w.cif")

    rules = [
        AltlocRules(strategy="preserve"),
        AltlocRules(
            strategy="select_best_occupancy", tie_breaker="lowest_b_factor"
        ),
        AltlocRules(strategy="select_user_defined", user_defined_altloc="B"),
    ]
    for altloc_rules in rules:
        policy = canonicalisationPolicy(
            policy_id="p",
            policy_name="p",
            policy_version="1.0.0",
            altloc_rules=altloc_rules,
        )
        canonical, mappings, _ = canonicalise_structure(structure, policy)
        atoms = [a for a in canonical.atoms if a.label_seq_id == 35]
        print(f"{altloc_rules.strategy}: {len(atoms)} atoms in residue 35")
        for item in mappings.altloc_selection_mapping.items:
            print(
                f"  {item.residue_id}: kept {item.selected_altloc} of "
                f"{item.available_altlocs} ({item.selection_reason})"
            )
    ```

    ```text
    preserve: 9 atoms in residue 35
    select_best_occupancy: 7 atoms in residue 35
      35:MSE: kept A of ['A', 'B'] (best_occupancy)
    select_user_defined: 7 atoms in residue 35
      35:MSE: kept B of ['A', 'B'] (user_defined)
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/1b6w.cif raw/
    cat > best.yaml <<'EOF'
    policy_id: p
    policy_name: p
    policy_version: 1.0.0
    altloc_rules:
      strategy: select_best_occupancy
      tie_breaker: lowest_b_factor
    EOF
    cat > pick_b.yaml <<'EOF'
    policy_id: p
    policy_name: p
    policy_version: 1.0.0
    altloc_rules:
      strategy: select_user_defined
      user_defined_altloc: B
    EOF
    pandora canonicalise --input-dir raw/ --policy best.yaml --output-dir best/
    pandora canonicalise --input-dir raw/ --policy pick_b.yaml --output-dir pick_b/
    for run in best pick_b; do
      jq -r --arg run "$run" '.["1B6W"].altloc_selection_mapping.items[]
        | "\($run): \(.residue_id) kept \(.selected_altloc)"
          + " of \(.available_altlocs) (\(.selection_reason))"' \
        "$run/canonicalisation_mappings.json"
    done
    ```

    ```text
    canonicalised 1 structures -> best
    canonicalised 1 structures -> pick_b
    best: 35:MSE kept A of ["A","B"] (best_occupancy)
    pick_b: 35:MSE kept B of ["A","B"] (user_defined)
    ```

## Expand to the biological assembly

`assembly_rules.strategy: standardize_biological_assembly` picks one assembly (see `preferred_assembly_source`) and applies its symmetry operators, adding each copy as a new chain named `<chain>_<operator>`. The assembly mapping records which chain each copy came from.

=== "`library`"

    ```python
    from pandora.canonicalisation import canonicalise_structure
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.canonicalisation import (
        AssemblyRules,
        canonicalisationPolicy,
    )

    # 13dg: assembly 1 is chains A-D under operators 1, 2 and 3.
    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/13dg.cif")
    print("asymmetric unit:", [u.id for u in structure.asym_units])

    policy = canonicalisationPolicy(
        policy_id="p",
        policy_name="p",
        policy_version="1.0.0",
        assembly_rules=AssemblyRules(strategy="standardize_biological_assembly"),
    )
    canonical, mappings, _ = canonicalise_structure(structure, policy)
    print("after expansion:", [u.id for u in canonical.asym_units])
    for item in mappings.assembly_mapping.items:
        for copy in item.chain_copies[:2]:
            print(
                f"  {copy.canonical_chain_id} = chain {copy.source_chain_id} "
                f"under operator {copy.operator_id}"
            )
    ```

    ```text
    asymmetric unit: ['A', 'E', 'B', 'F', 'C', 'G', 'D', 'H']
    after expansion: ['A', 'E', 'B', 'F', 'C', 'G', 'D', 'H', 'A_2', 'A_3', 'B_2', 'B_3', 'C_2', 'C_3', 'D_2', 'D_3']
      A_2 = chain A under operator 2
      A_3 = chain A under operator 3
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/13dg.cif raw/
    cat > policy.yaml <<'EOF'
    policy_id: p
    policy_name: p
    policy_version: 1.0.0
    assembly_rules:
      strategy: standardize_biological_assembly
    EOF
    pandora canonicalise --input-dir raw/ --policy policy.yaml --output-dir canonical/
    pandora export --input canonical/13dg.cif --output 13dg.json
    jq -c '[.asym_units[].id]' 13dg.json
    jq -r '.["13DG"].assembly_mapping.items[].chain_copies[:2][]
      | "\(.canonical_chain_id) = chain \(.source_chain_id)"
        + " under operator \(.operator_id)"' \
      canonical/canonicalisation_mappings.json
    ```

    ```text
    canonicalised 1 structures -> canonical
    exported -> 13dg.json
    ["A","E","A_2","A_3","B","F","B_2","B_3","C","G","C_2","C_3","D","H","D_2","D_3"]
    A_2 = chain A under operator 2
    A_3 = chain A under operator 3
    ```

!!! bug "Known issue: chains outside the assembly are kept"
    Chains E–H above belong to assembly 2, not assembly 1, but they stay in the expanded structure. Until this is fixed, remove them yourself before computing contacts. `test_assembly_expansion_keeps_only_assembly_chains` pins the bug.

## Map modified residues

`modified_residue_rules.strategy: map_to_parent` renames modified polymer residues to their standard parent, for example selenomethionine (`MSE`) to methionine (`MET`), whose `SE` atom becomes `SD`. By default only `MSE` is mapped; list other residues in `comp_ids`, or leave it empty to map every modified residue whose parent the file names in `_pdbx_struct_mod_residue`. A residue with no known parent is left unchanged with a `MODIFIED_RESIDUE_UNMAPPED` warning; Pandora doesn't guess.

=== "`library`"

    ```python
    from pandora.canonicalisation import canonicalise_structure
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.canonicalisation import (
        ModifiedResidueRules,
        canonicalisationPolicy,
    )


    def mapped(path: str, rules: ModifiedResidueRules) -> None:
        structure, _, _ = mmcif_to_structure(path)
        policy = canonicalisationPolicy(
            policy_id="p",
            policy_name="p",
            policy_version="1.0.0",
            modified_residue_rules=rules,
        )
        _, mappings, provenance = canonicalise_structure(structure, policy)
        print(f"{structure.entry_id} {provenance.transforms[0]}")
        for item in mappings.modified_residue_mapping.items:
            print(
                f"  chain {item.chain_id} {item.seq_id}: "
                f"{item.original_comp_id} -> {item.parent_comp_id}"
            )


    # Selenomethionine, mapped by default once the strategy is on.
    mapped(
        "datasets/dev/mmcif/1b6w.cif",
        ModifiedResidueRules(strategy="map_to_parent"),
    )
    # Any residue named in comp_ids whose parent the file lists.
    mapped(
        "datasets/dev/mmcif/1a08.cif",
        ModifiedResidueRules(strategy="map_to_parent", comp_ids=["FTY"]),
    )
    ```

    ```text
    1B6W modified_residues:map_to_parent
      chain A 1: MSE -> MET
      chain A 35: MSE -> MET
    1A08 modified_residues:map_to_parent
      chain B 2: FTY -> TYR
      chain D 2: FTY -> TYR
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/1b6w.cif datasets/dev/mmcif/1a08.cif raw/
    cat > mse.yaml <<'EOF'
    policy_id: p
    policy_name: p
    policy_version: 1.0.0
    modified_residue_rules:
      strategy: map_to_parent
    EOF
    cat > fty.yaml <<'EOF'
    policy_id: p
    policy_name: p
    policy_version: 1.0.0
    modified_residue_rules:
      strategy: map_to_parent
      comp_ids: [FTY]
    EOF
    pandora canonicalise --input-dir raw/ --policy mse.yaml --output-dir mse/
    pandora canonicalise --input-dir raw/ --policy fty.yaml --output-dir fty/
    jq -r '.["1B6W"].modified_residue_mapping.items[]
      | "1B6W chain \(.chain_id) \(.seq_id): \(.original_comp_id) -> \(.parent_comp_id)"' \
      mse/canonicalisation_mappings.json
    jq -r '.["1A08"].modified_residue_mapping.items[]
      | "1A08 chain \(.chain_id) \(.seq_id): \(.original_comp_id) -> \(.parent_comp_id)"' \
      fty/canonicalisation_mappings.json
    ```

    ```text
    canonicalised 2 structures -> mse
    canonicalised 2 structures -> fty
    1B6W chain A 1: MSE -> MET
    1B6W chain A 35: MSE -> MET
    1A08 chain B 2: FTY -> TYR
    1A08 chain D 2: FTY -> TYR
    ```

## Handle missing atoms and residues

`missing_data_rules` has three parts. `missing_atoms` flags or drops residues missing a backbone atom. `missing_residues` flags or drops chains with gaps. `incomplete_chains` can also cut a gapped chain down to its longest unbroken stretch.

=== "`library`"

    ```python
    from pandora.canonicalisation import canonicalise_structure
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.canonicalisation import (
        IncompleteChainRules,
        MissingAtomsRules,
        MissingDataRules,
        canonicalisationPolicy,
    )


    def residues(structure, chain: str) -> int:
        return len(
            {a.label_seq_id for a in structure.atoms if a.label_asym_id == chain}
        )


    def canonicalise(structure, rules: MissingDataRules):
        policy = canonicalisationPolicy(
            policy_id="p",
            policy_name="p",
            policy_version="1.0.0",
            missing_data_rules=rules,
        )
        canonical, _, _ = canonicalise_structure(structure, policy)
        return canonical


    # 1a08 chain B: two residues lack part of their backbone.
    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1a08.cif")
    dropped = canonicalise(
        structure,
        MissingDataRules(
            missing_atoms=MissingAtomsRules(strategy="drop_partial_residue")
        ),
    )
    print(
        "drop_partial_residue, 1a08 chain B:",
        residues(structure, "B"),
        "->",
        residues(dropped, "B"),
        "residues",
    )

    # 10mv chain A has a gap; keep only its longest unbroken stretch.
    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/10mv.cif")
    truncated = canonicalise(
        structure,
        MissingDataRules(
            incomplete_chains=IncompleteChainRules(
                strategy="truncate_to_complete_regions"
            )
        ),
    )
    print(
        "truncate_to_complete_regions, 10mv chain A:",
        residues(structure, "A"),
        "->",
        residues(truncated, "A"),
        "residues",
    )
    ```

    ```text
    drop_partial_residue, 1a08 chain B: 4 -> 2 residues
    truncate_to_complete_regions, 10mv chain A: 231 -> 118 residues
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/1a08.cif raw/
    cat > drop.yaml <<'EOF'
    policy_id: p
    policy_name: p
    policy_version: 1.0.0
    missing_data_rules:
      missing_atoms:
        strategy: drop_partial_residue
    EOF
    pandora canonicalise --input-dir raw/ --policy drop.yaml --output-dir drop/
    pandora export --input raw/1a08.cif --output before.json
    pandora export --input drop/1a08.cif --output after.json
    # residues in 1a08 chain B, before and after
    jq '[.atoms[] | select(.label_asym_id == "B") | .label_seq_id] | unique | length' before.json after.json
    ```

    ```text
    canonicalised 1 structures -> drop
    exported -> before.json
    exported -> after.json
    4
    2
    ```

    `pandora export` writes a structure as JSON, which is the easiest way to count what is left from the shell.

!!! bug "Known issue: truncated structures can't be read back"
    The CLI example leaves out `truncate_to_complete_regions`: the mmCIF that `pandora canonicalise` writes after truncating 10mv can't be parsed again, because raw categories such as `_pdbx_sifts_xref_db` still refer to the removed residues. The library result is correct; only writing it to mmCIF and reading it back fails. `test_truncated_structure_round_trips_through_mmcif` pins the bug.

To measure how much of each chain is missing against SEQRES, without changing the structure, use [`chain_completeness()`](datasets.md#inspect-chain-completeness).

## Normalise entities

`entity_rules.strategy: standardize` renumbers entity IDs `1`, `2`, ...; `merge_equivalent_entities` gives polymer entities with the same canonical sequence one ID. PDB entries already number entities sequentially and give identical chains one entity, so on most files neither changes anything; they matter for files from other sources. See the [policy reference](../reference/policies.md#entity_rules).

## Strip waters, ions and ligands

`ligand_rules.strategy: filter` drops waters, ions and other ligands according to `keep_waters`, `keep_ions` and `keep_nonpolymer_ligands`.

=== "`library`"

    ```python
    from pandora.canonicalisation import canonicalise_structure, filter_ligands
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.canonicalisation import (
        LigandRules,
        canonicalisationPolicy,
    )
    from pandora.schemas.common import DiagnosticBundle

    # 104m: myoglobin, a sulfate ion, a heme, n-butyl isocyanide and waters.
    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")
    names = {e.id: e.pdbx_description for e in structure.entities}

    policy = canonicalisationPolicy(
        policy_id="p",
        policy_name="p",
        policy_version="1.0.0",
        ligand_rules=LigandRules(
            strategy="filter", keep_waters=False, keep_ions=False
        ),
    )
    canonical, _, _ = canonicalise_structure(structure, policy)
    print("in the policy:", len(structure.atoms), "->", len(canonical.atoms))
    kept = sorted({names[u.entity_id] for u in canonical.asym_units})
    print("  kept:", kept)

    # The same filter on its own, without the rest of canonicalisation.
    atoms, asym_units = filter_ligands(
        list(structure.atoms),
        list(structure.asym_units),
        structure.entities,
        LigandRules(strategy="filter", keep_waters=False, keep_ions=True),
        DiagnosticBundle(),
        structure.entry_id,
    )
    print("filter_ligands, keep ions:", len(structure.atoms), "->", len(atoms))
    ```

    ```text
    in the policy: 1450 -> 1266
      kept: ['MYOGLOBIN', 'N-BUTYL ISOCYANIDE', 'PROTOPORPHYRIN IX CONTAINING FE']
    filter_ligands, keep ions: 1450 -> 1271
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/104m.cif raw/
    cat > policy.yaml <<'EOF'
    policy_id: p
    policy_name: p
    policy_version: 1.0.0
    ligand_rules:
      strategy: filter
      keep_waters: false
      keep_ions: false
    EOF
    pandora canonicalise --input-dir raw/ --policy policy.yaml --output-dir canonical/
    pandora export --input raw/104m.cif --output before.json
    pandora export --input canonical/104m.cif --output after.json
    jq '.atoms | length' before.json after.json
    ```

    ```text
    canonicalised 1 structures -> canonical
    exported -> before.json
    exported -> after.json
    1450
    1266
    ```

### Filter ligands directly

`filter_ligands()`, used at the end of the library example above, is what `ligand_rules` calls internally. It's exported so you can filter ligands without the rest of canonicalisation; [curation's content rules](datasets.md#keep-only-the-polymer) reuse it the same way. No CLI command wraps it on its own.

## Trace what changed

`canonicalise_structure` returns three things. The structure is the result; `provenance.transforms` lists the rule groups that changed something, in the order they ran; and `mappings` holds one list per kind of change, so any canonical ID can be traced back to the deposited one. With `provenance_rules.emit_canonicalisation_report`, `provenance.report` also counts the warnings and errors the run raised.

=== "`library`"

    ```python
    from pandora.canonicalisation import canonicalise_structure
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.canonicalisation import (
        ChainIdRules,
        IdentifierRules,
        ModifiedResidueRules,
        ResidueNumberingRules,
        canonicalisationPolicy,
        canonicalisationProvenanceRules,
    )

    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1b6w.cif")

    policy = canonicalisationPolicy(
        policy_id="p",
        policy_name="p",
        policy_version="1.0.0",
        identifier_rules=IdentifierRules(
            chain_id=ChainIdRules(strategy="use_auth_chain_id"),
            residue_numbering=ResidueNumberingRules(strategy="renumber"),
        ),
        modified_residue_rules=ModifiedResidueRules(strategy="map_to_parent"),
        provenance_rules=canonicalisationProvenanceRules(
            emit_canonicalisation_report=True
        ),
    )
    canonical, mappings, provenance = canonicalise_structure(structure, policy)

    # Which rule groups changed something, in the order they ran.
    print("transforms:", provenance.transforms)
    # How many diagnostics the run raised (the diagnostics themselves are
    # not returned yet).
    print("report:", provenance.report)
    # One mapping list per kind of change, to trace any id back.
    for name, mapping in mappings:
        print(f"{name}: {len(mapping.items)} items")
    ```

    ```text
    transforms: ['modified_residues:map_to_parent', 'chain_id:use_auth_chain_id', 'missing_atoms:annotate', 'missing_residues:annotate', 'altloc:select_best_occupancy', 'residue_numbering:renumber']
    report: {'warnings': 0, 'errors': 0}
    chain_id_mapping: 2 items
    residue_number_mapping: 82 items
    assembly_mapping: 1 items
    entity_mapping: 2 items
    altloc_selection_mapping: 1 items
    modified_residue_mapping: 2 items
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/1b6w.cif raw/
    cat > policy.yaml <<'EOF'
    policy_id: p
    policy_name: p
    policy_version: 1.0.0
    identifier_rules:
      chain_id:
        strategy: use_auth_chain_id
      residue_numbering:
        strategy: renumber
    modified_residue_rules:
      strategy: map_to_parent
    provenance_rules:
      emit_canonicalisation_report: true
    EOF
    pandora canonicalise --input-dir raw/ --policy policy.yaml --output-dir canonical/
    jq -c '.["1B6W"] | .transforms, .report' canonical/canonicalisation_provenance.json
    jq -r '.["1B6W"] | to_entries[] | "\(.key): \(.value.items | length) items"' \
      canonical/canonicalisation_mappings.json
    ```

    ```text
    canonicalised 1 structures -> canonical
    ["modified_residues:map_to_parent","chain_id:use_auth_chain_id","missing_atoms:annotate","missing_residues:annotate","altloc:select_best_occupancy","residue_numbering:renumber"]
    {"warnings":0,"errors":0}
    chain_id_mapping: 2 items
    residue_number_mapping: 82 items
    assembly_mapping: 1 items
    entity_mapping: 2 items
    altloc_selection_mapping: 1 items
    modified_residue_mapping: 2 items
    ```

## Load a policy from YAML

Policies are plain data, so you can keep them in YAML next to your dataset. The library loads them with [`load_policy()`](ingestion.md#load-a-canonicalisation-policy); the CLI's `--policy` takes the same file. `datasets/canonicalisation.yaml` remaps chain IDs, renumbers residues, merges equivalent entities and drops waters and ions:

```yaml
# datasets/canonicalisation.yaml
policy_id: overview-remap
policy_name: Remap Chains And Renumber
policy_version: 1.0.0
description: >
  Remaps chain IDs and renumbers residues sequentially; selects altlocs by
  best occupancy (lowest B-factor tiebreak); merges equivalent entities;
  drops waters and ions. Mirrors the policy built in examples/overview.py.

identifier_rules:
  chain_id:
    strategy: remap
  residue_numbering:
    strategy: renumber

altloc_rules:
  strategy: select_best_occupancy
  tie_breaker: lowest_b_factor

entity_rules:
  strategy: merge_equivalent_entities

ligand_rules:
  strategy: filter
  keep_waters: false
  keep_ions: false
```

=== "`library`"

    ```python
    from pandora.canonicalisation import canonicalise_structure
    from pandora.ingestion.policy import load_policy
    from pandora.parsing import mmcif_to_structure

    policy = load_policy("datasets/canonicalisation.yaml")
    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")
    canonical, mappings, provenance = canonicalise_structure(structure, policy)

    print(policy.policy_id, policy.policy_version)
    print(provenance.transforms)
    ```

    ```text
    overview-remap 1.0.0
    ['chain_id:remap', 'missing_atoms:annotate', 'missing_residues:annotate', 'altloc:select_best_occupancy', 'residue_numbering:renumber', 'entity:merge_equivalent_entities', 'ligands:filter']
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/104m.cif raw/
    pandora canonicalise --input-dir raw/ --policy datasets/canonicalisation.yaml --output-dir canonical/
    jq -c '.["104M"] | [.policy_id, .policy_version], .transforms' \
      canonical/canonicalisation_provenance.json
    ```

    ```text
    canonicalised 1 structures -> canonical
    ["overview-remap","1.0.0"]
    ["chain_id:remap","missing_atoms:annotate","missing_residues:annotate","altloc:select_best_occupancy","residue_numbering:renumber","entity:merge_equivalent_entities","ligands:filter"]
    ```
