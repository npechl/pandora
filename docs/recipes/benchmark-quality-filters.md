# Benchmark quality filters

The quality filters of two published PPI benchmarks, each written as one curation policy in YAML. Both are built only from the rules on the [Datasets](../usage/datasets.md#curate-one-structure) page; nothing here is benchmark-specific code. The benchmark descriptions are in `examples/ppi/benchmarks/` in the repository, and neither benchmark is a preferred target: the same rules express the others too.

The script and its output below are checked by the test suite, like the usage-page examples.

## The policies

=== "PPI v2.1"

    ```yaml
    --8<-- "recipes/ppi_v2_1_quality.yaml"
    ```

=== "ATOM3D PIP"

    ```yaml
    --8<-- "recipes/atom3d_pip_quality.yaml"
    ```

## Run them

The script maps MSE to MET during canonicalisation, then curates the same eight fixtures with each policy.

```python
--8 < --"recipes/benchmark_quality.py"
```

```text
--8<-- "recipes/benchmark_quality.out"
```

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
