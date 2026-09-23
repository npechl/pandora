# Step 1

## Filters (very similar to ProteinFlow):
- Remove biounits contained in another biounit (same seq = at least 90% seq id)
- Remove biounits containing more than 5’000 amino acids
- Keep only if resolution ≤ 3Å
- No more than 30% missing residues (= at least 1 missing backbone coord) in the tails and 10% in the middle (sequence left after removing the missing tails)
- Exclude unnatural atoms and amino acids
- Remove when missing fasta sequence or discrepancy between fasta and PDB sequence
- Keep biounit if at least one protein passes the filters in the biounit

## Pre-processing:
- Align fasta sequence to PDB sequence
- Remove poly-X tails (like poly-histidine)
- Compute DSSP
- Get contacts (with types) from COCɑDA and map them to the sequence
- Detect missing residues at the interface – missing residue at interface = missing residue in-between 2 non-missing residues at the interface
- Get pH (csv from COCɑDA, see code)

# Step 2

First detect chains that are repeated in a biounit and take one representative if each (90% seq id)
Then Mmseqs search to get all seq-ids above 90%, coverage 90% (see code)
Cluster with connected components
Select one representative sequence & structure per cluster:
Select sequence and structure with the least missing residues
If multiple sequences with the same number of missing residues, take the most connected sequence in the cluster as representative

# Step 3

- Foldseek search, threshold 0.5, coverage 80%  retrieve TM-scores from the search (see code)

# Step 4

- Based on COCɑDA contacts, find all interacting chains (use all chains, not just representatives from previous step)
- Remove redundancies. Redundancy = The 2 sequences in one dimer are the same (at 90% seq id) than the 2 sequences in  the other dimer
- The contacts are the at the same location in the 2 interfaces (Jaccard score >  0.5?)
- Here we need 2 sets of dimers (but keep the link between the 2):
- Complete set with everything after removing the redundancies
- Representatives set = complete set after removing ”broad” redundancies (same definition as above but taking criterion 1) into account only)

# Step 5

- Chain-level similarity: For each pair of dimer (i, j), mat[i, j] = 1 if TM(i1, j1) ≥ 0.5 or TM(i1, j2) ≥ 0.5 or TM(i2, j1) ≥ 0.5 or TM(i2, j2) ≥ 0.5 else 0
- Interface-level similarity: iDist ? Foldseek interface ?
- Get similarities between all the interfaces (complete set of dimers) and get the max similarity between each pair of representatives
- Select threshold to binarize for the similarity matrix
- Aggregate the 2 similarity matrices (union)

# Step 6

- Just connected components on the dimers similarity matrix from the previous step

# Step 7

- Problem: large over-representation of homomers compared to heteromers (1:10 ratio)
- Solution:
- For clusters with only homomers, select 1 representative (the most connected dimer in the cluster)
- For clusters with heteromers, discard homomers, then iteratively until there is no protein left in the cluster:
- Take the least connected heteromer
- Remove all heteromers connected to it

# Step 8

- Discard all chains present in the positive set from the 90% seq id representatives set
- For each positive pair: 
- Take one of the two chains
- Discard:
- All chains with TM-score ≥ 0.5 to the selected chain (set 1)
- All chains binding to any of the chains in set 1 (set 2)
- All chains with TM-score ≥ 0.5 to any of the chains in set 2
- Select a chain at random among the remaining chains
- For the next negative pairs, discard this newly selected chain

# Step 9

- Using the indexation between the complete set of dimers and the representative set of dimers, for each selected positive dimers, get the list of all distinct interfaces that can occur between the two chains and link each interface to its corresponding structure
