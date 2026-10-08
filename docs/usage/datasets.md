# Datasets

`pandora.datasets` filters structures for inclusion (`curate_structure`, `deduplicate_structures`), measures chain completeness (`chain_completeness`), and reshapes a canonical `Structure` into flat, ML-friendlier records (`extract_*_records`, `entry_sequences`). See [Functions](../reference/functions.md#pandora.datasets) for full signatures.

Each curation section applies one kind of rule to fixtures from `datasets/dev/mmcif/`. The `library` tab is the Python code and what it prints; the `cli` tab does the same with `pandora curate`, writing the policy as YAML and reading `curation_exclusions.json` with [`jq`](https://jqlang.org/). Run the CLI examples from the repository root. The test suite runs every example on this page and checks it still prints what is shown. The [curation policy reference](../reference/policies.md#curation) lists every field.

The sections after curation use this setup:

```python
from pandora.parsing import mmcif_to_structure
from pandora.canonicalisation import canonicalise_structure
from pandora.metadata import collect_metadata
from pandora.schemas.canonicalisation import canonicalisationPolicy

policy = canonicalisationPolicy(
    policy_id="p", policy_name="p", policy_version="1.0.0"
)
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")
canonical, _, _ = canonicalise_structure(structure, policy)
metadata = collect_metadata(canonical)
```

## Curate one structure

`curate_structure()` applies quality, organism, and content rules and returns `(curated_structure, exclusions, provenance)`. A kept entry has `exclusions` empty, or one record per chain the chain rules removed; an excluded entry has `curated_structure=None` and ends with an entry-level record (`chain_id=None`). Provenance is always populated, regardless of outcome. `metadata` can be `None`; rules that need missing data (such as resolution) then apply their configured default instead of being skipped.

=== "`library`"

    ```python
    from pandora.datasets import curate_structure
    from pandora.metadata import collect_metadata
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.dataset import DatasetCurationPolicy, QualityRules

    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")
    metadata = collect_metadata(structure)

    policy = DatasetCurationPolicy(
        policy_id="c1", policy_name="Default", policy_version="1.0.0"
    )
    curated, exclusions, provenance = curate_structure(structure, metadata, policy)
    print("default:", curated is not None, exclusions)

    strict = policy.model_copy(
        update={"quality_rules": QualityRules(max_resolution=1.0)}
    )
    curated, exclusions, provenance = curate_structure(structure, metadata, strict)
    print("strict:", curated is not None)
    for e in exclusions:
        print(f"  {e.reason_code} - {e.message}")
    ```

    ```text
    default: True []
    strict: False
      RESOLUTION_THRESHOLD - resolution 1.71 exceeds limit 1.0
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/104m.cif raw/
    cat > default.yaml <<'EOF'
    policy_id: c1
    policy_name: Default
    policy_version: 1.0.0
    EOF
    cat > strict.yaml <<'EOF'
    policy_id: c2
    policy_name: Strict
    policy_version: 1.0.0
    quality_rules:
      max_resolution: 1.0
    EOF
    pandora curate --input-dir raw/ --policy default.yaml --output-dir default/
    pandora curate --input-dir raw/ --policy strict.yaml --output-dir strict/
    jq -r '.[] | "\(.entry_id) \(.chain_id // "entry"): \(.reason_code) - \(.message)"' strict/curation_exclusions.json
    ```

    ```text
    curated: 1 retained, 0 excluded -> default
    curated: 0 retained, 1 excluded -> strict
    104M entry: RESOLUTION_THRESHOLD - resolution 1.71 exceeds limit 1.0
    ```

    The CLI always runs `collect_metadata()` on each structure, and writes `curation_exclusions.json` only when something was excluded or removed.

!!! note "Migrating from the single-exclusion return"
    Before chain-level rules, `curate_structure` returned `(structure, exclusion | None, provenance)`. It now returns a list. Replace `if curated is None: excluded.append(exclusion)` with `excluded.extend(exclusions)`, which also keeps the records of chains removed from entries that were kept.

Rules run in this order: entry rules (resolution, method, R-factors, organism, non-standard residues, `max_atoms`), then chain rules, then `min_chain_length` and `min_polymer_chains` on the chains that are left, then content rules.

## Set resolution limits per method

`max_resolution` applies one limit to every entry. `max_resolution_by_method` sets a limit per experimental method and overrides `max_resolution` for entries with that method; an entry with several methods gets the strictest. Cryo-EM resolution comes from `_em_3d_reconstruction`. An entry whose method has no limit, when there is no `max_resolution` either, isn't checked at all.

=== "`library`"

    ```python
    from pandora.datasets import curate_structure
    from pandora.metadata import collect_metadata
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.dataset import DatasetCurationPolicy, QualityRules

    # 22jy: a cryo-EM structure at 2.2 A.
    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/22jy.cif")
    metadata = collect_metadata(structure)
    print(metadata.quality.experimental_method, metadata.quality.resolution)

    rules = {
        "one limit, 3.5": QualityRules(max_resolution=3.5),
        "X-ray 2.5, cryo-EM 2.0": QualityRules(
            max_resolution_by_method={
                "X-RAY DIFFRACTION": 2.5,
                "ELECTRON MICROSCOPY": 2.0,
            }
        ),
    }
    for label, quality_rules in rules.items():
        policy = DatasetCurationPolicy(
            policy_id="c",
            policy_name="c",
            policy_version="1.0.0",
            quality_rules=quality_rules,
        )
        curated, exclusions, _ = curate_structure(structure, metadata, policy)
        print(f"{label}: kept={curated is not None}", end="")
        print("".join(f" ({e.reason_code}: {e.message})" for e in exclusions))
    ```

    ```text
    Electron Microscopy 2.2
    one limit, 3.5: kept=True
    X-ray 2.5, cryo-EM 2.0: kept=False (RESOLUTION_THRESHOLD: resolution 2.2 exceeds limit 2.0)
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/22jy.cif raw/
    cat > one_limit.yaml <<'EOF'
    policy_id: c
    policy_name: c
    policy_version: 1.0.0
    quality_rules:
      max_resolution: 3.5
    EOF
    cat > by_method.yaml <<'EOF'
    policy_id: c
    policy_name: c
    policy_version: 1.0.0
    quality_rules:
      max_resolution_by_method:
        X-RAY DIFFRACTION: 2.5
        ELECTRON MICROSCOPY: 2.0
    EOF
    pandora curate --input-dir raw/ --policy one_limit.yaml --output-dir one_limit/
    pandora curate --input-dir raw/ --policy by_method.yaml --output-dir by_method/
    jq -r '.[] | "\(.entry_id) \(.chain_id // "entry"): \(.reason_code) - \(.message)"' by_method/curation_exclusions.json
    ```

    ```text
    curated: 1 retained, 0 excluded -> one_limit
    curated: 0 retained, 1 excluded -> by_method
    22JY entry: RESOLUTION_THRESHOLD - resolution 2.2 exceeds limit 2.0
    ```

## Check X-ray R-factors

`max_r_free`, `max_r_free_gap` (\|Rfree − Rwork\|) and `max_r_sym` check the entries whose method is in `rfactor_methods` (X-ray by default), so cryo-EM and NMR entries never fail them. `max_r_sym` uses Rsym, or Rmerge when Rsym isn't reported. A missing value keeps the entry unless `null_rfactor_behavior` is `exclude`.

=== "`library`"

    ```python
    from pandora.datasets import curate_structure
    from pandora.metadata import collect_metadata
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.dataset import DatasetCurationPolicy, QualityRules


    def check(entry_id: str, quality_rules: QualityRules) -> str:
        structure, _, _ = mmcif_to_structure(f"datasets/dev/mmcif/{entry_id}.cif")
        policy = DatasetCurationPolicy(
            policy_id="c",
            policy_name="c",
            policy_version="1.0.0",
            quality_rules=quality_rules,
        )
        curated, exclusions, _ = curate_structure(
            structure, collect_metadata(structure), policy
        )
        if curated is not None:
            return "kept"
        return f"{exclusions[-1].reason_code} ({exclusions[-1].message})"


    # 1a08 is X-ray. It reports Rmerge but neither Rsym nor Rfree.
    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1a08.cif")
    quality = collect_metadata(structure).quality
    print(
        f"1a08: r_free={quality.r_free} r_work={quality.r_work} "
        f"r_sym={quality.r_sym} r_merge={quality.r_merge}"
    )

    # Rsym falls back to Rmerge.
    print("max_r_sym=0.05:", check("1a08", QualityRules(max_r_sym=0.05)))
    # A missing Rfree is kept by default...
    print("max_r_free=0.25:", check("1a08", QualityRules(max_r_free=0.25)))
    # ...or excluded if you ask for it.
    print(
        "max_r_free=0.25, exclude nulls:",
        check(
            "1a08",
            QualityRules(max_r_free=0.25, null_rfactor_behavior="exclude"),
        ),
    )
    # R-factor rules only check X-ray entries, so cryo-EM 22jy passes.
    print(
        "22jy, same rule:",
        check(
            "22jy",
            QualityRules(max_r_free=0.25, null_rfactor_behavior="exclude"),
        ),
    )
    ```

    ```text
    1a08: r_free=None r_work=0.192 r_sym=None r_merge=0.086
    max_r_sym=0.05: RSYM_THRESHOLD (r_sym 0.086 exceeds 0.05)
    max_r_free=0.25: kept
    max_r_free=0.25, exclude nulls: NULL_RFACTOR (r_free is null and null_rfactor_behavior='exclude')
    22jy, same rule: kept
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/1a08.cif datasets/dev/mmcif/22jy.cif raw/
    cat > rsym.yaml <<'EOF'
    policy_id: c
    policy_name: c
    policy_version: 1.0.0
    quality_rules:
      max_r_sym: 0.05
    EOF
    cat > rfree.yaml <<'EOF'
    policy_id: c
    policy_name: c
    policy_version: 1.0.0
    quality_rules:
      max_r_free: 0.25
      null_rfactor_behavior: exclude
    EOF
    pandora curate --input-dir raw/ --policy rsym.yaml --output-dir rsym/
    jq -r '.[] | "\(.entry_id) \(.chain_id // "entry"): \(.reason_code) - \(.message)"' rsym/curation_exclusions.json
    pandora curate --input-dir raw/ --policy rfree.yaml --output-dir rfree/
    jq -r '.[] | "\(.entry_id) \(.chain_id // "entry"): \(.reason_code) - \(.message)"' rfree/curation_exclusions.json
    ```

    ```text
    curated: 1 retained, 1 excluded -> rsym
    1A08 entry: RSYM_THRESHOLD - r_sym 0.086 exceeds 0.05
    curated: 1 retained, 1 excluded -> rfree
    1A08 entry: NULL_RFACTOR - r_free is null and null_rfactor_behavior='exclude'
    ```

## Reject non-standard residues

`exclude_nonstandard_residues` excludes an entry with any polymer residue that isn't one of the 20 amino acids, `UNK`, or a standard nucleotide. `allowed_nonstandard_residues` lets named residues through. Curation sees the canonical structure, so residues [mapped to their parent](canonicalisation.md#map-modified-residues) during canonicalisation, such as MSE, already count as standard.

=== "`library`"

    ```python
    from pandora.canonicalisation import canonicalise_structure
    from pandora.datasets import curate_structure
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.canonicalisation import (
        ModifiedResidueRules,
        canonicalisationPolicy,
    )
    from pandora.schemas.dataset import DatasetCurationPolicy, QualityRules


    def curate(structure, quality_rules: QualityRules) -> str:
        policy = DatasetCurationPolicy(
            policy_id="c",
            policy_name="c",
            policy_version="1.0.0",
            quality_rules=quality_rules,
        )
        curated, exclusions, _ = curate_structure(structure, None, policy)
        return "kept" if curated else exclusions[-1].message


    reject = QualityRules(exclude_nonstandard_residues=True)

    # 1a08 contains ACE and DIP caps and FTY, a modified tyrosine.
    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1a08.cif")
    print("1a08:", curate(structure, reject))
    allow = reject.model_copy(
        update={"allowed_nonstandard_residues": ["ACE", "DIP", "FTY"]}
    )
    print("1a08, allow-list:", curate(structure, allow))

    # 1b6w contains MSE. Mapping it to MET in canonicalisation makes it
    # standard by the time curation runs.
    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1b6w.cif")
    print("1b6w as deposited:", curate(structure, reject))
    canonical, _, _ = canonicalise_structure(
        structure,
        canonicalisationPolicy(
            policy_id="p",
            policy_name="p",
            policy_version="1.0.0",
            modified_residue_rules=ModifiedResidueRules(strategy="map_to_parent"),
        ),
    )
    print("1b6w after MSE -> MET:", curate(canonical, reject))
    ```

    ```text
    1a08: non-standard polymer residues ['ACE', 'DIP', 'FTY']
    1a08, allow-list: kept
    1b6w as deposited: non-standard polymer residues ['MSE']
    1b6w after MSE -> MET: kept
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/1a08.cif datasets/dev/mmcif/1b6w.cif raw/
    cat > mse.yaml <<'EOF'
    policy_id: p
    policy_name: Map MSE
    policy_version: 1.0.0
    modified_residue_rules:
      strategy: map_to_parent
    EOF
    cat > reject.yaml <<'EOF'
    policy_id: c
    policy_name: c
    policy_version: 1.0.0
    quality_rules:
      exclude_nonstandard_residues: true
    EOF
    cat > allow.yaml <<'EOF'
    policy_id: c
    policy_name: c
    policy_version: 1.0.0
    quality_rules:
      exclude_nonstandard_residues: true
      allowed_nonstandard_residues: [ACE, DIP, FTY]
    EOF
    pandora canonicalise --input-dir raw/ --policy mse.yaml --output-dir mapped/
    pandora curate --input-dir raw/ --policy reject.yaml --output-dir reject/
    jq -r '.[] | "\(.entry_id) \(.chain_id // "entry"): \(.reason_code) - \(.message)"' reject/curation_exclusions.json
    pandora curate --input-dir raw/ --policy allow.yaml --output-dir allow/
    jq -r '.[] | "\(.entry_id) \(.chain_id // "entry"): \(.reason_code) - \(.message)"' allow/curation_exclusions.json
    pandora curate --input-dir mapped/ --policy reject.yaml --output-dir mapped_reject/
    jq -r '.[] | "\(.entry_id) \(.chain_id // "entry"): \(.reason_code) - \(.message)"' mapped_reject/curation_exclusions.json
    ```

    ```text
    canonicalised 2 structures -> mapped
    curated: 0 retained, 2 excluded -> reject
    1A08 entry: NONSTANDARD_RESIDUE - non-standard polymer residues ['ACE', 'DIP', 'FTY']
    1B6W entry: NONSTANDARD_RESIDUE - non-standard polymer residues ['MSE']
    curated: 1 retained, 1 excluded -> allow
    1B6W entry: NONSTANDARD_RESIDUE - non-standard polymer residues ['MSE']
    curated: 1 retained, 1 excluded -> mapped_reject
    1A08 entry: NONSTANDARD_RESIDUE - non-standard polymer residues ['ACE', 'DIP', 'FTY']
    ```

## Drop chains with missing residues

`max_missing_tail_fraction` and `max_missing_middle_fraction` remove chains with too many missing residues, measured against SEQRES by [`chain_completeness()`](#inspect-chain-completeness). The tail fraction is the missing N- and C-terminal residues over the SEQRES length; the middle fraction is the missing residues between the first and last present ones over the SEQRES length minus the missing tails. `missing_residue_definition` decides whether a residue missing part of its backbone counts as missing (the default, `incomplete_backbone`) or only one with no atoms at all (`unobserved`).

Removing a chain keeps the rest of the entry. `min_polymer_chains` and `min_chain_length` then check what is left.

=== "`library`"

    ```python
    from pandora.datasets import curate_structure
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.dataset import DatasetCurationPolicy, QualityRules

    # 1aui chain A misses 95 residues in its middle; chain B is nearly whole.
    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1aui.cif")


    def curate(quality_rules: QualityRules) -> None:
        policy = DatasetCurationPolicy(
            policy_id="c",
            policy_name="c",
            policy_version="1.0.0",
            quality_rules=quality_rules,
        )
        curated, exclusions, _ = curate_structure(structure, None, policy)
        for e in exclusions:
            print(f"  {e.chain_id or 'entry'}: {e.reason_code} - {e.message}")
        if curated is not None:
            chains = sorted(
                {a.label_asym_id for a in curated.atoms if a.label_seq_id}
            )
            print(f"  kept polymer chains: {chains}")


    print("max 30% missing tails, 10% missing middle:")
    curate(
        QualityRules(max_missing_tail_fraction=0.3, max_missing_middle_fraction=0.1)
    )
    print("the same, but the entry needs two chains:")
    curate(
        QualityRules(
            max_missing_tail_fraction=0.3,
            max_missing_middle_fraction=0.1,
            min_polymer_chains=2,
        )
    )
    ```

    ```text
    max 30% missing tails, 10% missing middle:
      A: MISSING_MIDDLE - missing middle 0.201 exceeds max_missing_middle_fraction=0.1
      kept polymer chains: ['B']
    the same, but the entry needs two chains:
      A: MISSING_MIDDLE - missing middle 0.201 exceeds max_missing_middle_fraction=0.1
      entry: TOO_FEW_CHAINS - 1 polymer chain(s) found, fewer than min_polymer_chains=2
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/1aui.cif raw/
    cat > chains.yaml <<'EOF'
    policy_id: c
    policy_name: c
    policy_version: 1.0.0
    quality_rules:
      max_missing_tail_fraction: 0.3
      max_missing_middle_fraction: 0.1
    EOF
    cat > two_chains.yaml <<'EOF'
    policy_id: c
    policy_name: c
    policy_version: 1.0.0
    quality_rules:
      max_missing_tail_fraction: 0.3
      max_missing_middle_fraction: 0.1
      min_polymer_chains: 2
    EOF
    pandora curate --input-dir raw/ --policy chains.yaml --output-dir chains/
    jq -r '.[] | "\(.entry_id) \(.chain_id // "entry"): \(.reason_code) - \(.message)"' chains/curation_exclusions.json
    pandora curate --input-dir raw/ --policy two_chains.yaml --output-dir two_chains/
    jq -r '.[] | "\(.entry_id) \(.chain_id // "entry"): \(.reason_code) - \(.message)"' two_chains/curation_exclusions.json
    ```

    ```text
    curated: 1 retained, 0 excluded -> chains
    1AUI A: MISSING_MIDDLE - missing middle 0.201 exceeds max_missing_middle_fraction=0.1
    curated: 0 retained, 1 excluded -> two_chains
    1AUI A: MISSING_MIDDLE - missing middle 0.201 exceeds max_missing_middle_fraction=0.1
    1AUI entry: TOO_FEW_CHAINS - 1 polymer chain(s) found, fewer than min_polymer_chains=2
    ```

When a missing-fraction rule is set, a chain that can't be measured (no SEQRES, or a renumbered chain) is removed as `NO_SEQRES`.

## Cap the size

`max_atoms` excludes a whole entry above an atom count; `max_chain_length` removes single chains with more observed residues than the limit. `max_atoms` is checked after parsing, so it protects later steps such as contacts and export but not the parse itself.

=== "`library`"

    ```python
    from pandora.datasets import curate_structure
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.dataset import DatasetCurationPolicy, QualityRules


    def curate(entry_id: str, quality_rules: QualityRules) -> None:
        structure, _, _ = mmcif_to_structure(f"datasets/dev/mmcif/{entry_id}.cif")
        policy = DatasetCurationPolicy(
            policy_id="c",
            policy_name="c",
            policy_version="1.0.0",
            quality_rules=quality_rules,
        )
        curated, exclusions, _ = curate_structure(structure, None, policy)
        print(f"{entry_id}: kept={curated is not None}")
        for e in exclusions:
            print(f"  {e.chain_id or 'entry'}: {e.reason_code} - {e.message}")


    # Whole entries: 1aui has 4,832 atoms.
    curate("1aui", QualityRules(max_atoms=4000))
    # Single chains: 1aui chain A has 378 observed residues, chain B 165.
    curate("1aui", QualityRules(max_chain_length=300))
    ```

    ```text
    1aui: kept=False
      entry: TOO_MANY_ATOMS - 4832 atoms exceed max_atoms=4000
    1aui: kept=True
      A: CHAIN_TOO_LONG - 378 residues exceed max_chain_length=300
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/1aui.cif raw/
    cat > atoms.yaml <<'EOF'
    policy_id: c
    policy_name: c
    policy_version: 1.0.0
    quality_rules:
      max_atoms: 4000
    EOF
    cat > length.yaml <<'EOF'
    policy_id: c
    policy_name: c
    policy_version: 1.0.0
    quality_rules:
      max_chain_length: 300
    EOF
    pandora curate --input-dir raw/ --policy atoms.yaml --output-dir atoms/
    jq -r '.[] | "\(.entry_id) \(.chain_id // "entry"): \(.reason_code) - \(.message)"' atoms/curation_exclusions.json
    pandora curate --input-dir raw/ --policy length.yaml --output-dir length/
    jq -r '.[] | "\(.entry_id) \(.chain_id // "entry"): \(.reason_code) - \(.message)"' length/curation_exclusions.json
    ```

    ```text
    curated: 0 retained, 1 excluded -> atoms
    1AUI entry: TOO_MANY_ATOMS - 4832 atoms exceed max_atoms=4000
    curated: 1 retained, 0 excluded -> length
    1AUI A: CHAIN_TOO_LONG - 378 residues exceed max_chain_length=300
    ```

## Filter by organism

`organism_rules` keeps (`include_taxa`) or drops (`exclude_taxa`) entries by NCBI taxonomy ID, read from the metadata. Without taxonomy metadata, an active organism filter excludes the entry as `MISSING_TAXONOMY`.

=== "`library`"

    ```python
    from pandora.datasets import curate_structure
    from pandora.metadata import collect_metadata
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.dataset import DatasetCurationPolicy, OrganismRules

    policy = DatasetCurationPolicy(
        policy_id="c",
        policy_name="Human only",
        policy_version="1.0.0",
        organism_rules=OrganismRules(include_taxa=["9606"]),
    )
    for entry_id in ["1aui", "1ayi"]:
        structure, _, _ = mmcif_to_structure(f"datasets/dev/mmcif/{entry_id}.cif")
        metadata = collect_metadata(structure)
        curated, exclusions, _ = curate_structure(structure, metadata, policy)
        organisms = sorted({t.organism_scientific for t in metadata.taxonomies})
        outcome = "kept" if curated else exclusions[-1].reason_code
        print(f"{entry_id} {organisms}: {outcome}")
    ```

    ```text
    1aui ['Homo sapiens']: kept
    1ayi ['Escherichia coli']: ORGANISM_EXCLUDED
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/1aui.cif datasets/dev/mmcif/1ayi.cif raw/
    cat > human.yaml <<'EOF'
    policy_id: c
    policy_name: Human only
    policy_version: 1.0.0
    organism_rules:
      include_taxa: ["9606"]
    EOF
    pandora curate --input-dir raw/ --policy human.yaml --output-dir human/
    jq -r '.[] | "\(.entry_id) \(.chain_id // "entry"): \(.reason_code) - \(.message)"' human/curation_exclusions.json
    ```

    ```text
    curated: 1 retained, 1 excluded -> human
    1AYI entry: ORGANISM_EXCLUDED - taxa=['562'] excluded by organism_rules
    ```

## Keep only the polymer

`content_rules` never excludes anything. It strips ligands, waters or ions from structures that pass, the same way [`filter_ligands()`](canonicalisation.md#filter-ligands-directly) does.

=== "`library`"

    ```python
    from pandora.datasets import curate_structure
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.dataset import ContentRules, DatasetCurationPolicy

    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")

    policy = DatasetCurationPolicy(
        policy_id="c",
        policy_name="Protein only",
        policy_version="1.0.0",
        content_rules=ContentRules(
            keep_ligands=False, keep_waters=False, keep_ions=False
        ),
    )
    curated, exclusions, _ = curate_structure(structure, None, policy)
    print(len(structure.atoms), "->", len(curated.atoms), "atoms;", exclusions)
    ```

    ```text
    1450 -> 1217 atoms; []
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/104m.cif raw/
    cat > protein.yaml <<'EOF'
    policy_id: c
    policy_name: Protein only
    policy_version: 1.0.0
    content_rules:
      keep_ligands: false
      keep_waters: false
      keep_ions: false
    EOF
    pandora curate --input-dir raw/ --policy protein.yaml --output-dir protein/
    pandora export --input raw/104m.cif --output before.json
    pandora export --input protein/104m.cif --output after.json
    jq '.atoms | length' before.json after.json
    ```

    ```text
    curated: 1 retained, 0 excluded -> protein
    exported -> before.json
    exported -> after.json
    1450
    1217
    ```

## Collect exclusions across a batch

Curating many entries, extend one list with every call's records. It then holds entry records (`chain_id` is `None`) for entries left out, and chain records for chains removed from entries that were kept. Write it out with `write_records()` or pass it to [`build_dataset_manifest()`](provenance.md#assemble-a-dataset-manifest) as `excluded`. The CLI writes the same list to `curation_exclusions.json`.

=== "`library`"

    ```python
    from pandora.datasets import curate_structure
    from pandora.metadata import collect_metadata
    from pandora.parsing import mmcif_to_structure
    from pandora.schemas.dataset import DatasetCurationPolicy, QualityRules

    policy = DatasetCurationPolicy(
        policy_id="c",
        policy_name="Sharp, complete chains",
        policy_version="1.0.0",
        quality_rules=QualityRules(
            max_resolution=2.5,
            max_missing_middle_fraction=0.1,
        ),
    )

    kept = {}
    exclusions = []
    for entry_id in ["104m", "1aui", "1a08", "1p58"]:
        structure, _, _ = mmcif_to_structure(f"datasets/dev/mmcif/{entry_id}.cif")
        curated, records, _ = curate_structure(
            structure, collect_metadata(structure), policy
        )
        exclusions.extend(records)  # chain records from kept entries too
        if curated is not None:
            kept[entry_id] = curated

    entries_out = [e for e in exclusions if e.chain_id is None]
    chains_out = [e for e in exclusions if e.chain_id is not None]
    print("kept:", sorted(kept))
    print("entries excluded:", [(e.entry_id, e.reason_code) for e in entries_out])
    print(
        "chains removed:",
        [(e.entry_id, e.chain_id, e.reason_code) for e in chains_out],
    )
    ```

    ```text
    kept: ['104m', '1a08', '1aui']
    entries excluded: [('1P58', 'RESOLUTION_THRESHOLD')]
    chains removed: [('1AUI', 'A', 'MISSING_MIDDLE')]
    ```

=== "`cli`"

    ```bash
    mkdir -p raw && cp datasets/dev/mmcif/104m.cif datasets/dev/mmcif/1aui.cif datasets/dev/mmcif/1a08.cif datasets/dev/mmcif/1p58.cif raw/
    cat > policy.yaml <<'EOF'
    policy_id: c
    policy_name: Sharp, complete chains
    policy_version: 1.0.0
    quality_rules:
      max_resolution: 2.5
      max_missing_middle_fraction: 0.1
    EOF
    pandora curate --input-dir raw/ --policy policy.yaml --output-dir curated/
    jq -r '.[] | select(.chain_id == null)
      | "entry excluded: \(.entry_id) \(.reason_code)"' curated/curation_exclusions.json
    jq -r '.[] | select(.chain_id != null)
      | "chain removed: \(.entry_id) \(.chain_id) \(.reason_code)"' \
      curated/curation_exclusions.json
    ```

    ```text
    curated: 3 retained, 1 excluded -> curated
    entry excluded: 1P58 RESOLUTION_THRESHOLD
    chain removed: 1AUI A MISSING_MIDDLE
    ```

For whole benchmark policies built from these rules, see the [benchmark quality filters recipe](../recipes/benchmark-quality-filters.md).

## Inspect chain completeness

`chain_completeness()` counts each polymer chain's missing residues against SEQRES (`_entity_poly_seq`): at the N-terminus, at the C-terminus, and in the middle. The missing-residue rules use it; you can also call it yourself, for example to pick the chain with the fewest missing residues. A CA-only model shows why the definition matters: every residue is observed, but none has a full backbone. No CLI command wraps it.

```python
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
```

```text
A: SEQRES 521, missing 13 at the N-terminus, 35 at the C-terminus, 95 in the middle
B: SEQRES 169, missing 4 at the N-terminus, 0 at the C-terminus, 0 in the middle
1p58 chain A, unobserved: missing 0 + 0 tails, 5 middle of 495
1p58 chain A, incomplete_backbone: missing 495 + 0 tails, 0 middle of 495
```

Run it on a structure whose `label_seq_id`s still index SEQRES. After `residue_numbering: renumber` or `use_auth_seq` they don't, and the chain is reported as `SEQRES_MISMATCH` instead of measured.

## Deduplicate a batch

`deduplicate_structures()` drops every structure after the first with a
given `entry_id`:

=== "`library`"

    ```python
    from pandora.datasets import deduplicate_structures
    from pandora.schemas.dataset import DeduplicationRules

    duplicate, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")
    duplicate_canonical, _, _ = canonicalise_structure(duplicate, policy)

    retained, removed, provenance = deduplicate_structures(
        [canonical, duplicate_canonical], DeduplicationRules(enabled=True)
    )
    print(
        len(retained),
        "retained,",
        len(removed),
        "removed:",
        [r.reason_code for r in removed],
    )
    # 1 retained, 1 removed: ['DUPLICATE']
    ```

    `DeduplicationRules(enabled=False)` (the default) is a no-op — every
    structure passes through, `removed` is always `[]`.

=== "`cli`"

    ```bash
    pandora dedup --input-dir curated/ --output-dir deduped/
    # dedup: 5 retained, 0 removed -> deduped/
    ```

    `--disable` runs `DeduplicationRules(enabled=False)` instead — every
    structure passes through, but dedup provenance is still recorded.

## Reshape into flat records

Each `extract_*_records()` function turns one canonical `Structure`
into a list of flat, independently-serializable records — see
[Export](export.md) for writing them to disk. No CLI subcommand wraps
this reshaping step.

```python
from pandora.datasets import (
    entry_sequences,
    extract_chain_records,
    extract_residue_records,
    extract_interface_records,
)

# One representative sequence per entry (its longest polymer chain) —
# what compute_sequence_similarity()/entry-keyed clustering expect.
sequences = entry_sequences({canonical.entry_id: canonical})
print(sequences)
# {'104M': 'MVLSEGEWQLVLHVWAKVEAD...'}  (157 residues)

chains = extract_chain_records(canonical)
print([(c.chain_id, c.residue_count) for c in chains])
# [('A', 153)]

residues = extract_residue_records(canonical)
print(len(residues), residues[0].chain_id, residues[0].comp_id)
# 153 A VAL

structure2, _, _ = mmcif_to_structure("datasets/dev/mmcif/1a3n.cif")
canonical2, _, _ = canonicalise_structure(structure2, policy)
interfaces = extract_interface_records(canonical2, distance_cutoff=4.0)
print(
    len(interfaces),
    interfaces[0].chain_id_1,
    interfaces[0].chain_id_2,
    interfaces[0].contact_count,
)
# 5 A B 35
```

`extract_residue_records()` carries each residue's full atom list (coordinates, B-factor, ...) rather than a bare count — residue-level ML use cases (contact maps, solvent exposure) don't need to go back to the structure. `extract_interface_records()` reshapes [`annotate_chain_interfaces()`](annotation.md#chain-chain-interfaces)'s output; it doesn't compute contacts itself. It takes the same `atom_set` and `polymer_types` arguments, and each `InterfaceRecord` carries the `atom_set` used and the `residue_pairs`.
