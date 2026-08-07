> **CORRECTION (Aug 7, 2026):** Two confirmed defects in the original pipeline
> (gap-stripping column misalignment; 8-bit QM wrap-around on 400-cell maps)
> invalidated the biological numbers in this file. See
> [CORRECTION_NOTICE.md](CORRECTION_NOTICE.md) for the verified corrected
> results (21 variable positions, 10 co-evolving pairs, 36 distinct rules (2 essential)) and the
> corrected pipeline (analysis/corrected_pipeline.py). Numbers in this file
> describe the original (buggy) run unless stated otherwise.

# 01. coevolution_shared.py - The Shared Toolbox

## What the Program Does

`coevolution_shared.py` is not an analysis by itself. It is the shared toolbox that 7 other scripts import (master_boolean, kmap_boolean, allseq_constraint, predictive_constraint, flipped_boolean, advanced, gpu_full, plus the report generators). It performs the low-level operations that every analysis needs:

1. Parse the FASTA alignment file into (header, sequence) pairs.
2. Convert each sequence into an array of integer codes (0 to 19 for the 20 amino acids, -1 for gaps or unknown characters) using the base-20 He 2012 encoding.
3. Compute Shannon entropy per position.
4. Compute mutual information (MI) between position pairs.
5. Compute coupling constants and the constraint function.
6. Find co-evolving position pairs using GPU acceleration with a CPU fallback.

The key design decision is that there is exactly one implementation of each operation. Every script uses the same MI formula, the same encoding, the same threshold. This guarantees that results from different scripts are comparable.

**Sequence length analyzed: 1,276 positions (full length), all 1,299 sequences.**

## What Data It Looks At

It loads all 1,299 sequences. For MI it looks at **position pairs within a sliding window**. Two positions i and j are compared if they are at most `max_gap` positions apart (default 30). It does not compare all 813,450 unordered pairs (1,276 choose 2) for the pair-finder; it compares only pairs with |i - j| <= 30, which is the biologically relevant range because co-evolution is usually local in the sequence.

The MI matrix computation (`mi_matrix_gpu`) can compute all pairs when asked, but the pair-finder uses the window.

```mermaid
flowchart TD
    A[FASTA file: 1,299 sequences] --> B[parse_fasta]
    B --> C[Base-20 He 2012 encoding]
    C --> D[position arrays: int 0-19, -1 for gap]
    D --> E[entropy per position]
    D --> F["MI for pairs |i-j| <= 30"]
    D --> G[majority reference per position]
    F --> H[co-evolving pairs list]
    G --> H
    H --> I[other scripts consume]
```

## The Encoding: Which Residue Gets Which Number

The He 2012 ordering maps each amino acid to a number 0 to 19:

| Index | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 |
|-------|---|---|---|---|---|---|---|---|---|---|----|----|----|----|----|----|----|----|----|----|
| AA | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |

This ordering is taken from He et al. 2012 (Physica A 391:93) and is the same ordering proved in the Lean theorem file `BaseNAminoEncoding.lean` (n-ary-kmap repo).

### Worked Example: Sequence to Position Array

Take the first sequence in the file, **WRU87367.1**, which begins `MFVFLVLLPLVSSQCVNLITRTQ...`. The first 10 residues and their He 2012 codes:

| Residue | M | F | V | F | L | V | L | L | P | L |
|---------|---|---|---|---|---|---|---|---|---|---|
| He index | 4 | 5 | 3 | 5 | 2 | 3 | 2 | 2 | 18 | 2 |
| 5-bit Gray | 6 | 7 | 2 | 7 | 3 | 2 | 3 | 3 | 27 | 3 |

The position array for this sequence contains 4 at index 0, 5 at index 1, 3 at index 2, and so on. When the script counts joint occurrences at (372, 401), it counts how many sequences have each (code_i, code_j) pair.

The 5-bit Gray code used by the binary K-map scripts is `gray(n) = n XOR (n >> 1)`: for M (index 4), gray(4) = 4 XOR 2 = 6; for F (index 5), gray(5) = 5 XOR 2 = 7.

## Formulas

### Shannon entropy (per position)

$$H(p) = -\sum_{a=1}^{20} P(a) \log_2 P(a)$$

Units: bits. H = 0 means fully conserved (one amino acid everywhere). H = log2(20) = 4.32 bits means uniform over all 20 amino acids.

**Worked example.** Position 372 has entropy H = 1.6328 bits and perplexity 2^H = 3.101. This means the effective number of amino acids at position 372 is about 3.1, even though the position varies across the dataset.

Position 401 has entropy H = 1.6715 bits and perplexity 3.185.

### Mutual information (between two positions)

$$MI(i, j) = \sum_{a=1}^{20} \sum_{b=1}^{20} P(a,b) \log_2 \frac{P(a,b)}{P(a) P(b)}$$

Units: bits. MI = 0 means the two positions are statistically independent. MI > 0 means knowing the residue at one position reduces uncertainty about the other.

**Worked example.** MI(372, 401) = 1.5917 bits. The joint distribution at these two positions over all 1,299 sequences is:

| (residue at 372, residue at 401) | Count | Frequency |
|----------------------------------|-------|-----------|
| (A, N) | 814 | 0.6266 |
| (T, D) | 260 | 0.2002 |
| (K, V) | 114 | 0.0878 |
| (F, E) | 45 | 0.0346 |
| (C, R) | 32 | 0.0246 |
| others | 34 | 0.0262 |

The marginal frequencies: P(372 = A) = 0.6266 + other occurrences of A with different partners. Because the joint is far from the product of marginals (A at 372 nearly always pairs with N at 401), the MI is large.

### Coupling constant / constraint function

$$J(i, j, a, b) = \ln \frac{P(a, b)}{P(a) P(b)}$$

J > 0 means the pair (a, b) is more common than expected under independence (co-evolutionary). J < 0 means less common than expected (anti-correlated). This is the DCA coupling sign convention (Morcos et al. 2011).

## Results Used by Other Scripts

- 1,249 variable positions (entropy > 0.3) out of 1,276.
- 36,918 co-evolving pairs (mutation-only MI > 0.1, window 30).
- Majority reference residues per position (the most common amino acid at that position).
- The pair-finder returns (pos_i, pos_j, MI, n_mutations, ref_i, ref_j) tuples sorted by MI descending; the top 15 of these feed the master Boolean rule extraction.

## Inference

The results of this toolbox are the inputs to every other script: position arrays, entropy, MI, references, and co-evolving pairs. The inference that carries through the whole project is that the dataset contains 1,249 variable positions and 36,918 co-evolving pairs: the Spike protein is densely coupled, and every downstream analysis builds on these shared numbers.

The most important downstream consequence: the co-evolving pairs (by mutation-only MI) become the K-maps whose minimization produces the 36-rule master Boolean function (script 11). All 10 pairs are used (see script 11).

## Scholar Questions and Answers

**Q: Why is the window 30 and not the whole sequence?**
A: Co-evolution between positions far apart in the sequence is usually mediated by folding (3D proximity), not direct sequence coupling. The window keeps the computation focused on locally coupled positions and is a common choice in co-evolution analysis. The full MI matrix (all 813,450 pairs) is still computed by `create_mi_heatmap.py` and `gpu_full_analysis.py`; only the pair-finder uses the window.

**Q: What does -1 mean in the position arrays?**
A: It marks gaps ('-') and ambiguous characters in the alignment. Positions with -1 are excluded from frequency and MI computations because they carry no amino acid information.

**Q: Why mutation-only MI?**
A: The reference (majority) pair at each position pair is the "wild-type-like" combination. Counting only mutations (pairs that differ from the reference) isolates the co-variation of changes, rather than the trivial fact that most sequences share the reference pair. This is why the master Boolean rules are about mutations.

**Q: How does the 5-bit Gray code differ from the He 2012 code?**
A: The He 2012 code is the direct 0-19 index used by the base-20 (n-ary) K-map. The 5-bit Gray code is `gray(n) = n XOR (n >> 1)`, applied to the physicochemical group order (script 03), and is used by the binary 32x32 K-map. Both are formally verified in the Lean proofs; they answer different questions (base-20 direct adjacency vs binary Gray adjacency).
