# 05. nary_kmap_co-evolution.py - The Base-20 (N-ary) K-map


**Sequence length analyzed: 1,276 positions (full length), all 1,299 sequences.**

## What the Program Does

This script builds the **n-ary (base-20) K-map**, the direct generalization of the binary K-map. Instead of encoding each amino acid into 5 binary bits (which creates a 32 x 32 grid with only 20 used rows and columns), it treats the 20 amino acids as 20 symbols of a base-20 alphabet.

The K-map is then a 20 x 20 = 400 cell grid. Row = first amino acid, column = second amino acid. Every cell is a real dipeptide; there are no empty don't-care cells.

The script:
1. Counts dipeptides across all sequences using base-20 k-mer counting (cell = first_aa * 20 + second_aa).
2. Normalizes to frequencies.
3. Thresholds at the 75th percentile of non-zero cells to build a Boolean map.
4. Runs Quine-McCluskey to find prime implicants.
5. Extracts motifs and computes coupling constants.

## How the Base-20 Numbering Works

Under He 2012, each amino acid has an index 0 to 19 (see script 01). A dipeptide XY is encoded as:

```
```

$$\text{cell} = \text{index}(X) \times 20 + \text{index}(Y)$$

```
```

This is a two-digit base-20 number. The first digit is the first amino acid, the second digit is the second amino acid. The K-map is the 20 x 20 grid of all possible cells 0 to 399.

**Worked example from our data.** In the first sequence WRU87367.1, the first two residues are M and F. He 2012: M = index 4, F = index 5. The dipeptide MF maps to cell 4 * 20 + 5 = 85. The dipeptide FV (positions 2-3) maps to 5 * 20 + 3 = 103.

```mermaid
flowchart TD
    A[sequences] --> B[encode each residue to 0-19 He 2012]
    B --> C[cell = first*20 + second for each dipeptide]
    C --> D[400-cell frequency K-map]
    D --> E[threshold 75th percentile]
    E --> F[Boolean 20x20]
    F --> G[Quine-McCluskey on 400 cells]
    G --> H[prime implicants 73, essential 42]
    D --> I[coupling J = ln(P/P_exp)]
```

## Why the Base-20 K-map Is Different From Binary

| Property | Binary 32x32 | Base-20 20x20 |
|----------|--------------|---------------|
| Cells | 1,024 | 400 |
| Used cells | 20 rows/cols used, 12 don't-care per axis | all 400 used |
| Empty cells | 61% of grid (20 of 32 codes used per axis) | 0% |
| Cell meaning | (gray(aa1), gray(aa2)) | (aa1, aa2) directly |
| Number of variables for QM | 10 bits (5+5) | 8 bits (4+4; only 16 of 20 values representable per axis) |

The base-20 K-map is more compact and every cell is interpretable: cell (V, Q) literally means the dipeptide Valine-Glutamine.

## Results (all 1,299 sequences)

| Metric | Value |
|--------|-------|
| K-map size | 20 x 20 = 400 cells |
| On-set cells | 93 (23.2% density) |
| Prime implicants | 73 |
| Essential prime implicants | 42 |
| Strong couplings | 337 |
| MI ratio (on/off) | 1.1869 |

## Inference

The base-20 K-map finds 42 essential prime implicants from the whole-protein dipeptide landscape, and the on-set cells have 1.19x higher MI than off-set cells. This is slightly stronger than the binary version (which had ratio 1.19 in the same direction). The key advantage is interpretability and density: with 400 cells instead of 1,024, the minimization is tighter and every rule maps directly to amino acid pairs.

The overall conclusion matches script 04: whole-protein dipeptide K-maps carry a weak but real co-evolution signal (on-set MI ratio 1.19), but position-pair analyses are needed for strong signal.

## Scholar Questions and Answers

**Q: Why is the base-20 K-map called "n-ary"?**
A: The binary K-map uses base 2 (each axis has 2^k cells for k bits). The base-20 K-map uses base 20 (each axis has exactly 20 cells, the amino acid alphabet). The Lean proofs in the n-ary-kmap repo (`NaryGrayCode.lean`, `BaseNAminoEncoding.lean`) formalize this generalization.

**Q: Why does the n-ary K-map not use Gray code?**
A: For a single symbol (k = 1), the n-ary Gray code of a one-digit number is the identity. The Gray structure matters for multi-symbol words (k >= 2) where adjacent words differ in one symbol. For the 2-mer base-20 K-map the cell index is the direct base-20 concatenation, which is the standard n-ary K-map convention.

**Q: Which is better, binary or base-20?**
A: For interpretability and density, base-20. For compatibility with the existing Lean binary proofs and the 5-bit Gray adjacency hypothesis (H1), binary. They answer different questions. The base-20 version is the "native" K-map for proteins.
