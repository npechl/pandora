# Parsing and export

Two halves of one round-trip: `pandora.parsing.mmcif_to_structure()`
reads raw mmCIF into Pandora's typed `Structure`; `pandora.export`
writes a `Structure` (or the flat records from
[Datasets](datasets.md)) back out. See
[Functions](../reference/functions.md#pandora.parsing) and
[Functions](../reference/functions.md#pandora.export) for full
signatures.

## Parse an mmCIF file

`mmcif_to_structure()` returns a 3-tuple: the parsed `Structure` (or
`None` if parsing failed outright), a `DiagnosticBundle` of warnings
encountered, and a status string. No CLI subcommand exposes parsing on
its own — `pandora export` (below) and every other stage's `--input-dir`
call it internally.

```python
from pandora.parsing import mmcif_to_structure

structure, diagnostics, status = mmcif_to_structure(
    "datasets/dev/mmcif/104m.cif"
)
print(
    status,
    len(structure.atoms),
    "atoms,",
    len(diagnostics.warnings),
    "warnings",
)
# success 1450 atoms, 0 warnings
```

`status` is `"success"`, `"warning"` (parsed, but check
`diagnostics.warnings`), or `"failed"` (`structure` is `None` — check
`diagnostics.errors`). Every other mmCIF category not promoted to a
typed field is preserved verbatim in `structure.raw` — see
[Metadata](metadata.md#read-a-raw-category-directly) for reading it.

## Write a Structure back to mmCIF

`structure_to_mmcif()` rebuilds the standard mmCIF categories from a
`Structure`'s typed fields and re-emits everything in `.raw`. It's a
lossy round-trip (loop ordering, comments, and any field mmCIF
supports but Pandora's schema doesn't carry won't survive), but the
result reparses cleanly:

=== "`library`"

    ```python
    from pathlib import Path
    from pandora.canonicalisation import canonicalise_structure
    from pandora.export import structure_to_mmcif
    from pandora.schemas.canonicalisation import canonicalisationPolicy

    policy = canonicalisationPolicy(
        policy_id="p", policy_name="p", policy_version="1.0.0"
    )
    canonical, _, _ = canonicalise_structure(structure, policy)

    output_dir = Path("./datasets/output/export/")
    mmcif_path = structure_to_mmcif(canonical, output_dir / "104m.canonical.cif")

    reparsed, _, reparse_status = mmcif_to_structure(str(mmcif_path))
    print(reparse_status, len(reparsed.atoms), "atoms")
    # success 1450 atoms
    ```

=== "`cli`"

    ```bash
    pandora export --input canonical/104m.cif --output 104m.export.cif
    # exported -> 104m.export.cif
    ```

    The output format is inferred from `--output`'s suffix — `.cif`
    (or anything but `.json`) calls `structure_to_mmcif()`.

## Write any model to JSON

`write_json()` writes a single Pydantic model — a `Structure`,
provenance record, `AnnotationLayer`, anything in `pandora/schemas/` —
as pretty-printed JSON.

=== "`library`"

    ```python
    from pandora.export import write_json

    json_path = write_json(canonical, output_dir / "104m.json")
    print(json_path)
    # datasets/output/export/104m.json
    ```

=== "`cli`"

    ```bash
    pandora export --input canonical/104m.cif --output 104m.json
    # exported -> 104m.json
    ```

    `pandora export` only ever writes the parsed `Structure` itself as
    JSON this way — for provenance records, `AnnotationLayer`s, or
    anything else `write_json()` accepts, use the library form.

## Write a list of records

`write_records()` writes a list of record models — `ChainRecord`,
`ResidueRecord`, `InterfaceRecord` from [Datasets](datasets.md) — and
dispatches on the output path's suffix. `pandora export` only handles a
single `Structure`, so there's no CLI equivalent for record lists:

```python
from pandora.datasets import extract_chain_records
from pandora.export import write_records

chain_records = extract_chain_records(canonical)

write_records(
    chain_records, output_dir / "chains.jsonl"
)  # one JSON object per line
write_records(chain_records, output_dir / "chains.json")  # a JSON array
write_records(
    chain_records, output_dir / "chains.parquet"
)  # columnar, needs the `export` extra
```

`.parquet` needs `pandas`/`pyarrow` (`uv sync --extra export`) and
raises `RuntimeError` if they're missing, so it's safe to try/except
around when you don't know if the extra is installed. An empty list
always raises `ValueError("no records to write")` rather than writing
an empty file — a batch that curated down to nothing is a case to
handle explicitly, not silently produce an empty artifact.
