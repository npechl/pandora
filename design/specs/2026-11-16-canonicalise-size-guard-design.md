# Canonicalise size guard (short design)

## Why

`pandora canonicalise` ran out of memory on PDB entry 10ic in a full-PDB
run of the PPI v2.1 example. The parser builds one pydantic object per
atom (about 3× gemmi's memory; spike: 6.5 GB for 3j3q, 2.44M atoms), and
expanding assembly 1 multiplies the atom count again. Curation's
`max_atoms` exists to keep such entries out, but it runs after parsing
and canonicalisation, so it cannot prevent the crash. Entries this large
may stay out of the dataset for now.

## Goals / non-goals

- Goal: `canonicalise` skips an entry that would be too large *before*
  building the `Structure`, records why, and carries on.
- Non-goal: a leaner data model or streaming parser.
- Anti-goal: a hard-coded limit, or a silent skip.

## Approaches

- **A. Estimate before parsing, in `canonicalise` (picked).** Read the
  file with gemmi, count atoms in the asymmetric unit, and multiply each
  assembly generator's chains by its operator count for the assembly the
  policy will build. Over the limit → failure line, skip. The only option
  that prevents the crash in general.
- B. Check after parsing, before canonicalising. Exact count, but the
  pydantic build is what runs out of memory.
- C. Pre-filter the example's entry list (file size or an ID list).
  Approximate, and every new large entry needs a manual fix.

## Design

- **Policy field:** `max_atoms: int | None = None` on
  `canonicalisationPolicy` (top level; `None` = no check, today's
  behaviour). Documented in `docs/reference/policies.md` and
  `docs/policies/canonicalisation.yaml`.
- **Estimate:** a public function in `pandora/canonicalisation/` that
  takes a file path and the policy and returns the estimated atom count:
  model 1 of the asymmetric unit when the assembly strategy does not
  expand; otherwise, for the assembly the policy selects (same rule as
  `_select_biological_assembly`), the sum over generators of
  (atoms in the generator's chains × number of operators). Uses gemmi
  only, never builds pydantic records.
- **Meaning:** the count is taken before ligand/water filtering, so it
  is an upper bound and not the same number as curation's `max_atoms`.
  The docs say so.
- **CLI:** `canonicalise` runs the estimate when `max_atoms` is set.
  Over the limit: one `failures.jsonl` line with
  `error_type: "TooLarge"` and a message giving the estimate and the
  limit; counts gain `entries_too_large` (also counted in
  `entries_failed`).
- **PPI v2.1 example:** `dataset/policies/canonicalisation.yaml` sets
  `max_atoms: 1000000`, matching the curation limit.
- **Test:** a dev fixture with `max_atoms` set below its size lands in
  `failures.jsonl` as `TooLarge`; one under the limit is kept; the
  estimate equals the atom count of the canonicalised structure before
  ligand filtering on a fixture with a multi-operator assembly.

## Calls made (yours to override)

- Top-level policy field instead of a new rules group (one value, no
  group needed yet).
- Limit 1,000,000 for v2.1, same as curation. Because the estimate
  includes waters and ligands, an entry just over the line may be skipped
  that curation would have kept.
- Peak memory of the gemmi read itself is about a third of the pydantic
  parse; if that still fails on codon, a lighter `_atom_site` row count
  is the next step.
