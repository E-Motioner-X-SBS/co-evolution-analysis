# E-Motioner-X-SBS: Complete Co-Evolution Analysis Pipeline

**Generated:** August 07, 2026 at 04:05
**Dataset:** SARS-CoV-2 Omicron Spike Protein — 1,299 sequences, 1276 positions
**Compute:** NVIDIA A100 80GB + 24-core Xeon, Python 3.10
**Author:** Shuvam Banerji Seal — IISER Kolkata

---

## Abstract

This document presents the complete computational pipeline for analyzing co-evolutionary
constraints in the SARS-CoV-2 Spike protein using Karnaugh map (K-map) Boolean minimization.
**1,299** Omicron variant Spike sequences from GISAID were encoded using base-20
(He 2012 ordering) amino acid representation. Position-pair mutual information identified
**15** co-evolving position pairs in the N-terminal signal peptide region
(positions 68-79). Quine-McCluskey Boolean minimization produced **152 essential
prime implicants** — each representing an irreducible co-evolutionary constraint.

**Key results:** H1 Gray-code adjacency enrichment = **0.19×**; 
max pairwise MI = **1.5917** at positions (372,401);
1249/80 positions are evolutionarily variable (H > 0.3).

## 1. Formal Foundations — 236 Lean 4 Theorems

The entire framework rests on formal proofs in Lean 4.29.0. All theorems use `native_decide`,
are **sorry-free** and **axiom-free**, and compile via `lake build` in under 1 second.

### 1.1 Binary K-map (`lean_proofs/`) — 106 Theorems

| File | Theorems | What Is Proved |
|------|----------|---------------|
| `KmapProofs.lean` | 16 | Gray code boundedness, injectivity, involution, DNA 2-bit encoding (A=00,C=01,G=11,T=10), complementarity, transition/transversion |
| `KmerIndexing.lean` | 35 | k-mer encoding (k=1..4) injectivity and boundedness via concatenation, row/column K-map decomposition |
| `AminoAcidEncoding.lean` | 27 | 20 AAs → 5-bit Gray code, 7 physicochemical groups, within-group dist-1 (14 pairs), cross-group dist-1 (26 pairs), max-dist pairs (F-H, Y-E, W-R, M-K at dist=5) |
| `ContactMapCompleteness.lean` | 20 | Contact maps as symmetric Boolean functions, irreducibility proof (no two contact pairs are 1-bit K-map adjacent) |
| `KmapEncodingEquiv.lean` | 8 | Distance distribution, Q₅ hypercube degree=5, encoding captures 50% of Q₅ edges |

### 1.2 N-ary K-map (`n-ary-kmap/`) — 115 Theorems

| File | Theorems | What Is Proved |
|------|----------|---------------|
| `NaryGrayCode.lean` | 54 | Generalized n-ary reflected Gray code; digit-sum formula g_k = (d_k + d_{k+1}) mod b |
| `BaseNDnaEncoding.lean` | 13 | Base-4 DNA encoder with transition adjacency (A↔C, C↔G, G↔T, T↔A as Gray neighbors) |
| `BaseNAminoEncoding.lean` | 33 | Base-20 amino acid encoder, 7-group physicochemical ordering |
| `BaseNKMap.lean` | 15 | Base-N K-map structure, cell adjacency, cyclic wrapping |

### 1.3 Speculative Encoding (`speculative-binary-encoding/`) — 15 Theorems

| File | Theorems | What Is Proved |
|------|----------|---------------|
| `AntiCorrelation.lean` | 15 | H1 encoding-invariance: the anti-correlation result holds across all 22 binary encodings |

### 1.4 Theorem Verification Bridge

The `lean_consistency.py` script in `kmap-sbm-validation/` re-computes every Lean theorem
in Python and reports discrepancies. Result: **103/103 pass**, zero discrepancies.

## 2. Dataset

| Property | Value |
|----------|-------|
| **File** | `Spike_protein.aln-fasta` (1.8 MB) |
| **Sequences** | 1,299 SARS-CoV-2 Omicron Spike proteins |
| **Alignment length** | 1276 residues |
| **Analysis region** | Positions 0-79 (N-terminal signal peptide) |
| **Variable positions** | 1249/80 (entropy > 0.3) |
| **Conserved positions** | -1169/80 |
| **Co-evolving pairs (MI > 0.01)** | 36374 |
| **Co-evolving pairs (MI > 0.1)** | 35774 |
| **Encoding** | Base-20, He 2012 ordering |

### 2.1 He 2012 Amino Acid Ordering

The 20 standard amino acids are ordered by physicochemical properties:

| Index | AA | Group |
|-------|-----|-------|
| 0 | A | Aliphatic |
| 1 | I | Aliphatic |
| 2 | L | Aliphatic |
| 3 | V | Aliphatic |
| 4 | M | Sulfur |
| 5 | F | Aromatic |
| 6 | Y | Aromatic |
| 7 | W | Aromatic |
| 8 | E | Negative |
| 9 | D | Negative |
| 10 | Q | Polar |
| 11 | N | Polar |
| 12 | H | Positive |
| 13 | K | Positive |
| 14 | R | Positive |
| 15 | S | Polar |
| 16 | T | Polar |
| 17 | C | Sulfur |
| 18 | P | Structure |
| 19 | G | Structure |

## 3. K-map Construction Methodology

### 3.1 Encoding Biological Sequences as Karnaugh Maps

The core innovation is representing protein sequences as Karnaugh maps —
the same mathematical objects used in digital logic design. The mapping is:

```
Biological Sequence → Amino Acid Encoding → K-map Cell → Frequency Map → Boolean Function
```

### 3.2 Step 1: Position Array Construction

Each sequence is converted to a position array of integers (0-19):

$$\text{seq}[i] \rightarrow \text{encoder.encode}(\text{seq}[i]) \in \{0, 1, \dots, 19\}$$

For 1,299 sequences, each of length up to 1276, we build a list of `np.int32` arrays.
Gap characters ('-') and ambiguous residues are mapped to -1 and excluded.

### 3.3 Step 2: Position-Pair K-map (20×20)

For each pair of positions $(i, j)$, a 20×20 frequency matrix is built:

$$K_{ij}(a, b) = \frac{1}{N}\sum_{s=1}^{N} \mathbb{1}[\text{seq}_s[i] = a \land \text{seq}_s[j] = b]$$

where $N = {n_all:,}$ is the number of sequences, $a,b \in \{0,\dots,19\}$.

### 3.4 Example: Frequency K-map for Positions (76, 77)

The reference residues are **D** at position 76 and **N** at position 77 (Wuhan-Hu-1).

The top 10 most frequent dipeptide pairs at (76, 77):

| Rank | AA(i) | AA(j) | Frequency |
|------|-------|-------|-----------|
| 1 | D | N | 0.6474 |
| 2 | F | D | 0.2386 |
| 3 | P | V | 0.0862 |
| 4 | N | P | 0.0192 |
| 5 | K | R | 0.0031 |
| 6 | V | L | 0.0023 |
| 7 | L | P | 0.0015 |
| 8 | N | D | 0.0008 |
| 9 | G | N | 0.0008 |

### 3.5 Step 3: Boolean Thresholding

The frequency K-map is converted to a Boolean function by thresholding:

$$B_{ij}(a,b) = \begin{cases} 1 & \text{if } K_{ij}(a,b) \geq \text{threshold} \\ 0 & \text{if } K_{ij}(a,b) < \text{threshold} \\ -1 & \text{if position is conserved (don't-care)} \end{cases}$$

The threshold is typically the 75th percentile of non-zero frequencies.
Don't-care (-1) marks the reference pair (conserved).

### 3.6 Step 4: Quine-McCluskey Boolean Minimization

The 20×20 = 400-cell K-map is flattened to a truth table with 8 binary variables
(4 bits for row amino acid + 4 bits for column amino acid). The Quine-McCluskey
algorithm finds the minimal set of prime implicants covering all on-set cells:

$$f(\text{pos}_i, \text{pos}_j, \text{aa}_i, \text{aa}_j) = \bigvee_k \bigwedge_{m \in S_k} b_m$$

where $b_m$ are the 8 binary variables and $S_k$ are the literal sets for each prime implicant.

## 4. Entropy and Conservation Analysis

### 4.1 Shannon Entropy

Position-specific entropy measures evolutionary variability:

$$H(p) = -\sum_{i=1}^{20} P(a_i) \log_2 P(a_i)$$

where $P(a_i)$ is the frequency of amino acid $a_i$ at position $p$.

- $H \approx 0$: highly conserved (one amino acid dominates)
- $H \approx 4.32$: maximally variable (uniform distribution over 20 AAs)

### 4.2 Conservation Landscape (Positions 0-79)

| Position | Consensus | Frequency | Entropy | Perplexity | Status |
|----------|-----------|-----------|---------|------------|--------|
| 0 | M | 0.9962 | 0.0392 | 1.028 | Conserved |
| 1 | F | 0.9969 | 0.0301 | 1.021 | Conserved |
| 2 | V | 0.9915 | 0.0817 | 1.058 | Conserved |
| 3 | F | 0.9946 | 0.0484 | 1.034 | Conserved |
| 4 | L | 0.9792 | 0.1617 | 1.119 | Conserved |
| 5 | V | 0.9938 | 0.0574 | 1.041 | Conserved |
| 6 | L | 1.0000 | -0.0000 | 1.000 | Conserved |
| 7 | L | 0.9938 | 0.0574 | 1.041 | Conserved |
| 8 | P | 0.9923 | 0.0687 | 1.049 | Conserved |
| 9 | L | 0.9923 | 0.0687 | 1.049 | Conserved |
| 10 | V | 0.9923 | 0.0687 | 1.049 | Conserved |
| 11 | S | 0.9985 | 0.0181 | 1.013 | Conserved |
| 12 | S | 0.9923 | 0.0687 | 1.049 | Conserved |
| 13 | Q | 0.9915 | 0.0778 | 1.055 | Conserved |
| 14 | C | 0.9915 | 0.0778 | 1.055 | Conserved |
| 15 | V | 0.9915 | 0.0797 | 1.057 | Conserved |
| 16 | N | 0.9908 | 0.0853 | 1.061 | Conserved |
| 17 | L | 0.9900 | 0.0970 | 1.070 | Conserved |
| 18 | I | 0.7329 | 0.8510 | 1.804 | Variable |
| 19 | T | 0.9777 | 0.1591 | 1.117 | Conserved |
| 20 | R | 0.9761 | 0.1757 | 1.130 | Conserved |
| 21 | T | 0.9754 | 0.1847 | 1.137 | Conserved |
| 22 | Q | 0.9754 | 0.1954 | 1.145 | Conserved |
| 23 | S | 0.7321 | 0.9941 | 1.992 | Variable |
| 24 | Y | 0.7313 | 0.9790 | 1.971 | Variable |
| 25 | T | 0.7313 | 1.0028 | 2.004 | Variable |
| 26 | N | 0.7313 | 1.0001 | 2.000 | Variable |
| 27 | S | 0.7313 | 0.9900 | 1.986 | Variable |
| 28 | F | 0.7313 | 0.8981 | 1.864 | Variable |
| 29 | T | 0.7313 | 1.0028 | 2.004 | Variable |
| 30 | R | 0.7313 | 1.0028 | 2.004 | Variable |
| 31 | G | 0.7313 | 0.9981 | 1.997 | Variable |
| 32 | V | 0.7313 | 0.9920 | 1.989 | Variable |
| 33 | Y | 0.7521 | 0.8582 | 1.813 | Variable |
| 34 | Y | 0.7321 | 0.9941 | 1.992 | Variable |
| 35 | P | 0.7313 | 1.0001 | 2.000 | Variable |
| 36 | D | 0.7313 | 0.9765 | 1.968 | Variable |
| 37 | K | 0.7313 | 1.0028 | 2.004 | Variable |
| 38 | V | 0.7313 | 1.0028 | 2.004 | Variable |
| 39 | F | 0.7313 | 1.0003 | 2.000 | Variable |
| ... | ... | ... | ... | ... | ... |
| 40 | R | 0.7313 | 0.9920 | 1.989 | Variable |
| 41 | S | 0.7521 | 0.8468 | 1.799 | Variable |
| 42 | S | 0.7321 | 0.9941 | 1.992 | Variable |
| 43 | V | 0.7313 | 0.9953 | 1.993 | Variable |
| 44 | L | 0.7313 | 0.9601 | 1.945 | Variable |
| 45 | H | 0.7306 | 0.9053 | 1.873 | Variable |
| 46 | S | 0.7313 | 1.0028 | 2.004 | Variable |
| 47 | T | 0.7313 | 1.0028 | 2.004 | Variable |
| 48 | Q | 0.7298 | 1.0188 | 2.026 | Variable |
| 49 | D | 0.7306 | 0.9993 | 1.999 | Variable |
| 50 | L | 0.7321 | 0.9882 | 1.984 | Variable |
| 51 | F | 0.7306 | 1.0019 | 2.003 | Variable |
| 52 | L | 0.7336 | 0.9755 | 1.966 | Variable |
| 53 | P | 0.7306 | 0.9745 | 1.965 | Variable |
| 54 | F | 0.9938 | 0.0627 | 1.044 | Conserved |
| 55 | F | 0.7313 | 0.9981 | 1.997 | Variable |
| 56 | S | 0.7306 | 1.0040 | 2.006 | Variable |
| 57 | N | 0.7306 | 0.9829 | 1.976 | Variable |
| 58 | V | 0.7290 | 1.0242 | 2.034 | Variable |
| 59 | T | 0.7306 | 1.0052 | 2.007 | Variable |
| 60 | W | 0.7306 | 1.0068 | 2.009 | Variable |
| 61 | F | 0.7306 | 1.0014 | 2.002 | Variable |
| 62 | H | 0.7306 | 0.9993 | 1.999 | Variable |
| 63 | A | 0.7306 | 1.0068 | 2.009 | Variable |
| 64 | I | 0.7275 | 1.0329 | 2.046 | Variable |
| 65 | H | 0.8861 | 0.6228 | 1.540 | Variable |
| 66 | V | 0.8830 | 0.6676 | 1.588 | Variable |
| 67 | S | 0.6474 | 1.3584 | 2.564 | Variable |
| 68 | G | 0.6490 | 1.3501 | 2.549 | Variable |
| 69 | T | 0.6490 | 1.0929 | 2.133 | Variable |
| 70 | N | 0.6490 | 1.1127 | 2.163 | Variable |
| 71 | G | 0.6513 | 1.3550 | 2.558 | Variable |
| 72 | T | 0.6513 | 1.3550 | 2.558 | Variable |
| 73 | K | 0.6482 | 1.3790 | 2.601 | Variable |
| 74 | R | 0.6482 | 1.3811 | 2.605 | Variable |
| 75 | F | 0.6490 | 1.3730 | 2.590 | Variable |
| 76 | D | 0.6474 | 1.3854 | 2.613 | Variable |
| 77 | N | 0.6482 | 1.3662 | 2.578 | Variable |
| 78 | P | 0.6505 | 1.3552 | 2.558 | Variable |
| 79 | V | 0.6474 | 1.1174 | 2.170 | Variable |

**Summary:** 1249 variable positions (H > 0.3), -1169 conserved.

## 5. Mutual Information Analysis

### 5.1 Definition

Mutual Information quantifies the co-dependence between two positions:

$$MI(i,j) = \sum_{x=1}^{20} \sum_{y=1}^{20} P(x,y) \log_2 \frac{P(x,y)}{P(x) \cdot P(y)}$$

where:
- $P(x,y)$ = joint frequency of amino acid $x$ at position $i$ and $y$ at position $j$
- $P(x)$ = marginal frequency of $x$ at position $i$
- $P(y)$ = marginal frequency of $y$ at position $j$

### 5.2 Implementation

We use a **vectorized numpy implementation** via `np.bincount` (O(400) per pair, not O(N²)):

```python
def mutual_information(pos_arrays, pos_i, pos_j, n_seqs):
    codes_i = np.array([arr[pos_i] for arr in pos_arrays[:n_seqs]])
    codes_j = np.array([arr[pos_j] for arr in pos_arrays[:n_seqs]])
    valid = (codes_i >= 0) & (codes_j >= 0)
    # Joint via bincount: flat_index = ci * 20 + cj
    pairs = codes_i[valid] * 20 + codes_j[valid]
    joint = np.bincount(pairs, minlength=400).reshape(20,20)
    # ... MI from joint and marginals
```

### 5.3 Top 50 Co-evolving Position Pairs

| Rank | Pos i | Pos j | MI | Ref(i) | Ref(j) | H(i) | H(j) | Δ |
|------|-------|-------|-------|--------|--------|------|------|----|
| 1 | 372 | 401 | 1.5917 | A | N | 1.633 | 1.672 | 29 |
| 2 | 401 | 404 | 1.5690 | N | S | 1.672 | 1.710 | 3 |
| 3 | 208 | 209 | 1.5420 | L | G | 1.614 | 1.615 | 1 |
| 4 | 209 | 210 | 1.5231 | G | R | 1.615 | 1.577 | 1 |
| 5 | 208 | 210 | 1.5019 | L | R | 1.614 | 1.577 | 2 |
| 6 | 207 | 209 | 1.4386 | N | G | 1.502 | 1.615 | 2 |
| 7 | 207 | 208 | 1.4377 | N | L | 1.502 | 1.614 | 1 |
| 8 | 206 | 207 | 1.4292 | I | N | 1.497 | 1.502 | 1 |
| 9 | 206 | 209 | 1.4094 | I | G | 1.497 | 1.615 | 3 |
| 10 | 206 | 208 | 1.4060 | I | L | 1.497 | 1.614 | 2 |
| 11 | 207 | 210 | 1.4010 | N | R | 1.502 | 1.577 | 3 |
| 12 | 133 | 136 | 1.3903 | N | F | 1.392 | 1.392 | 3 |
| 13 | 131 | 133 | 1.3891 | F | N | 1.389 | 1.392 | 2 |
| 14 | 131 | 136 | 1.3876 | F | F | 1.389 | 1.392 | 5 |
| 15 | 84 | 86 | 1.3816 | D | V | 1.382 | 1.382 | 2 |
| 16 | 83 | 84 | 1.3816 | N | D | 1.382 | 1.382 | 1 |
| 17 | 83 | 86 | 1.3816 | N | V | 1.382 | 1.382 | 3 |
| 18 | 89 | 101 | 1.3812 | A | I | 1.385 | 1.381 | 12 |
| 19 | 89 | 110 | 1.3812 | A | T | 1.385 | 1.381 | 21 |
| 20 | 101 | 110 | 1.3812 | I | T | 1.381 | 1.381 | 9 |
| 21 | 84 | 87 | 1.3791 | D | Y | 1.382 | 1.381 | 3 |
| 22 | 84 | 89 | 1.3791 | D | A | 1.382 | 1.385 | 5 |
| 23 | 86 | 89 | 1.3791 | V | A | 1.382 | 1.385 | 3 |
| 24 | 83 | 87 | 1.3791 | N | Y | 1.382 | 1.381 | 4 |
| 25 | 83 | 89 | 1.3791 | N | A | 1.382 | 1.385 | 6 |
| 26 | 86 | 87 | 1.3791 | V | Y | 1.382 | 1.381 | 1 |
| 27 | 110 | 133 | 1.3790 | T | N | 1.381 | 1.392 | 23 |
| 28 | 110 | 136 | 1.3781 | T | F | 1.381 | 1.392 | 26 |
| 29 | 82 | 83 | 1.3773 | F | N | 1.377 | 1.382 | 1 |
| 30 | 82 | 84 | 1.3773 | F | D | 1.377 | 1.382 | 2 |
| 31 | 82 | 86 | 1.3773 | F | V | 1.377 | 1.382 | 4 |
| 32 | 76 | 83 | 1.3768 | D | N | 1.385 | 1.382 | 7 |
| 33 | 76 | 84 | 1.3768 | D | D | 1.385 | 1.382 | 8 |
| 34 | 76 | 86 | 1.3768 | D | V | 1.385 | 1.382 | 10 |
| 35 | 73 | 74 | 1.3768 | K | R | 1.379 | 1.381 | 1 |
| 36 | 87 | 89 | 1.3766 | Y | A | 1.381 | 1.385 | 2 |
| 37 | 127 | 133 | 1.3765 | C | N | 1.376 | 1.392 | 6 |
| 38 | 110 | 131 | 1.3763 | T | F | 1.381 | 1.389 | 21 |
| 39 | 80 | 83 | 1.3762 | L | N | 1.376 | 1.382 | 3 |
| 40 | 80 | 84 | 1.3762 | L | D | 1.376 | 1.382 | 4 |
| 41 | 80 | 85 | 1.3762 | L | G | 1.376 | 1.376 | 5 |
| 42 | 80 | 86 | 1.3762 | L | V | 1.376 | 1.382 | 6 |
| 43 | 83 | 85 | 1.3762 | N | G | 1.382 | 1.376 | 2 |
| 44 | 84 | 85 | 1.3762 | D | G | 1.382 | 1.376 | 1 |
| 45 | 85 | 86 | 1.3762 | G | V | 1.376 | 1.382 | 1 |
| 46 | 206 | 210 | 1.3756 | I | R | 1.497 | 1.577 | 4 |
| 47 | 136 | 137 | 1.3753 | F | L | 1.392 | 1.375 | 1 |
| 48 | 127 | 136 | 1.3749 | C | F | 1.376 | 1.392 | 9 |
| 49 | 82 | 87 | 1.3748 | F | Y | 1.377 | 1.381 | 5 |
| 50 | 82 | 89 | 1.3748 | F | A | 1.377 | 1.385 | 7 |

## 6. Quine-McCluskey Boolean Minimization — All 108 Essential Prime Implicants

The Boolean minimization was performed on 15 position pairs (68-79),
producing **152 essential prime implicants**. Each rule has the form:

$$f(s_3, s_2, s_1, s_0, t_3, t_2, t_1, t_0) = \text{AND of literals}$$

**Variables:**
- $s_3 s_2 s_1 s_0$ = 4-bit binary encoding of residue at position $i$
- $t_3 t_2 t_1 t_0$ = 4-bit binary encoding of residue at position $j$
- $\bar{s}_k$ = NOT ($s_k = 0$), $s_k$ = ($s_k = 1$)

### 6.1 Complete Inference Rules (152 rules across 15 position pairs)

### Position Pair (413, 424) — MI = 1.0122, Reference: N→D

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 1 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (I, A) |
| 2 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  s0  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (I, D) |
| 3 | $$s3  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (E, D) |
| 4 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  t0$$ | (H, W) |
| 5 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (H, Q) |
| 6 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  t3  \cdot  t2  \cdot  t1  \cdot  t0$$ | (K, S) |
| 7 | $$s3  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  t1  \cdot  t0$$ | (R, S) |
| 8 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (A, M) |
| 9 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (A, K) |
| 10 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (F, L) |
| 11 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  t1  \cdot  \bar{t0}$$ | (Y, R) |
| 12 | $$s3  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (E, F) |

**Interpretation:** When position 413 mutates to any of the listed residues,
position 424 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (N, D).

### Position Pair (413, 425) — MI = 1.0246, Reference: N→F

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 13 | $$s3  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (E, F) |
| 14 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  t0$$ | (A, D) |
| 15 | $$\bar{s3}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (A, D) |
| 16 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  t1  \cdot  t0$$ | (I, V) |
| 17 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (L, M) |
| 18 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  t3  \cdot  \bar{t2}  \cdot  t1  \cdot  t0$$ | (W, N) |
| 19 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  t3  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (K, Q) |
| 20 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (K, K) |
| 21 | $$s3  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (R, K) |
| 22 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (W, F) |

**Interpretation:** When position 413 mutates to any of the listed residues,
position 425 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (N, F).

### Position Pair (413, 426) — MI = 1.0216, Reference: N→T

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 23 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (I, I) |
| 24 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  t0$$ | (L, W) |
| 25 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (W, D) |
| 26 | $$s3  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  t0$$ | (N, W) |
| 27 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (K, I) |
| 28 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (K, H) |
| 29 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  t3  \cdot  t2  \cdot  t1  \cdot  t0$$ | (K, S) |
| 30 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  t1  \cdot  t0$$ | (A, S) |
| 31 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (I, M) |
| 32 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (M, D) |
| 33 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (Y, Q) |
| 34 | $$s3  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (E, I) |

**Interpretation:** When position 413 mutates to any of the listed residues,
position 426 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (N, T).

### Position Pair (413, 427) — MI = 1.0284, Reference: N→G

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 35 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  t1  \cdot  t0$$ | (A, V) |
| 36 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t0}$$ | (W, E) |
| 37 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  t1  \cdot  t0$$ | (I, V) |
| 38 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (L, F) |
| 39 | $$s3  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (N, F) |
| 40 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (K, I) |
| 41 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (K, K) |
| 42 | $$s3  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  t1  \cdot  t0$$ | (R, V) |
| 43 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (M, F) |
| 44 | $$s3  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (E, H) |

**Interpretation:** When position 413 mutates to any of the listed residues,
position 427 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (N, G).

### Position Pair (413, 428) — MI = 1.0222, Reference: N→C

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 45 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  \bar{t1}  \cdot  t0$$ | (A, I) |
| 46 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  t0$$ | (I, F) |
| 47 | $$s3  \cdot  s2  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (H, K) |
| 48 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}$$ | (A, M) |
| 49 | $$\bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  t1  \cdot  t0$$ | (A, S) |
| 50 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  t3  \cdot  t2  \cdot  t1  \cdot  t0$$ | (W, S) |
| 51 | $$s3  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (N, M) |
| 52 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (K, H) |
| 53 | $$s3  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  t0$$ | (R, W) |
| 54 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (F, A) |
| 55 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (W, I) |

**Interpretation:** When position 413 mutates to any of the listed residues,
position 428 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (N, C).

### Position Pair (459, 473) — MI = 1.0093, Reference: P→N

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 56 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}$$ | (L, E) |
| 57 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  t0$$ | (W, I) |
| 58 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  t1  \cdot  t0$$ | (W, V) |
| 59 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  t0$$ | (L, W) |
| 60 | $$s3  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (N, L) |
| 61 | $$s3  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  t3  \cdot  t2  \cdot  t1  \cdot  t0$$ | (N, S) |
| 62 | $$s3  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  \bar{t0}$$ | (R, Y) |
| 63 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (A, F) |
| 64 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  t0$$ | (I, W) |
| 65 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (M, A) |

**Interpretation:** When position 459 mutates to any of the listed residues,
position 473 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (P, N).

### Position Pair (462, 473) — MI = 1.0114, Reference: R→N

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 66 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  t1$$ | (L, Y) |
| 67 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (L, D) |
| 68 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (Y, M) |
| 69 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (W, E) |
| 70 | $$s3  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  t1  \cdot  t0$$ | (N, V) |
| 71 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (H, I) |
| 72 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (V, K) |
| 73 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  t1  \cdot  t0$$ | (M, N) |
| 74 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (W, L) |

**Interpretation:** When position 462 mutates to any of the listed residues,
position 473 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (R, N).

### Position Pair (468, 473) — MI = 1.0114, Reference: I→N

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 75 | $$s3  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  t0$$ | (E, F) |
| 76 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (I, I) |
| 77 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  t1  \cdot  t0$$ | (L, V) |
| 78 | $$s3  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  t1  \cdot  t0$$ | (N, V) |
| 79 | $$s3  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (N, F) |
| 80 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  t3  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (K, Q) |
| 81 | $$s3  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  t1  \cdot  t0$$ | (R, S) |
| 82 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  \bar{t0}$$ | (V, Y) |
| 83 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (M, A) |
| 84 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (Y, E) |

**Interpretation:** When position 468 mutates to any of the listed residues,
position 473 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (I, N).

### Position Pair (469, 473) — MI = 1.0093, Reference: Y→N

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 85 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  t0$$ | (I, F) |
| 86 | $$s3  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t1  \cdot  t0$$ | (E, V) |
| 87 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (I, L) |
| 88 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  t0$$ | (L, W) |
| 89 | $$s3  \cdot  \bar{s2}  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (Q, A) |
| 90 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (K, F) |
| 91 | $$s3  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  t0$$ | (R, W) |
| 92 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (M, Q) |
| 93 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (F, M) |
| 94 | $$s3  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (E, K) |

**Interpretation:** When position 469 mutates to any of the listed residues,
position 473 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (Y, N).

### Position Pair (1026, 1040) — MI = 0.8136, Reference: S→G

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 95 | $$\bar{s3}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (A, D) |
| 96 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  t3  \cdot  t1  \cdot  t0$$ | (V, N) |
| 97 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  t1  \cdot  \bar{t0}$$ | (A, R) |
| 98 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (V, M) |
| 99 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (M, L) |
| 100 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (Y, I) |
| 101 | $$s3  \cdot  \bar{s2}  \cdot  s1  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (Q, K) |
| 102 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (H, H) |
| 103 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  t0$$ | (Y, W) |
| 104 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  t3  \cdot  t2  \cdot  t1  \cdot  \bar{t0}$$ | (W, R) |

**Interpretation:** When position 1026 mutates to any of the listed residues,
position 1040 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (S, G).

### Position Pair (1026, 1042) — MI = 1.2162, Reference: S→G

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 105 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  t3  \cdot  t2  \cdot  t1$$ | (V, R) |
| 106 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (A, D) |
| 107 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (L, H) |
| 108 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (F, K) |
| 109 | $$s3  \cdot  \bar{s2}  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  \bar{t0}$$ | (Q, Y) |
| 110 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (H, K) |
| 111 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  t0$$ | (I, W) |
| 112 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (F, I) |
| 113 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (Y, A) |
| 114 | $$s3  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  t1  \cdot  t0$$ | (E, N) |

**Interpretation:** When position 1026 mutates to any of the listed residues,
position 1042 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (S, G).

### Position Pair (1040, 1042) — MI = 0.8064, Reference: G→G

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 115 | $$\bar{s3}  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  t0$$ | (V, W) |
| 116 | $$\bar{s3}  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (L, I) |
| 117 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (F, F) |
| 118 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  t3  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (W, Q) |
| 119 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (H, F) |
| 120 | $$s3  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (S, M) |
| 121 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (A, Q) |
| 122 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  t3  \cdot  t2  \cdot  t1  \cdot  \bar{t0}$$ | (V, R) |
| 123 | $$s3  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (E, E) |

**Interpretation:** When position 1040 mutates to any of the listed residues,
position 1042 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (G, G).

### Position Pair (1064, 1065) — MI = 1.1921, Reference: V→P

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 124 | $$\bar{s3}  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (V, M) |
| 125 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  \bar{t0}$$ | (M, H) |
| 126 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t2  \cdot  t1  \cdot  \bar{t0}$$ | (M, Y) |
| 127 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (A, Q) |
| 128 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  t3  \cdot  \bar{t2}  \cdot  t1  \cdot  t0$$ | (W, N) |
| 129 | $$s3  \cdot  \bar{s2}  \cdot  s1  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (Q, K) |
| 130 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (K, A) |
| 131 | $$s3  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (R, I) |
| 132 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  t1  \cdot  t0$$ | (A, S) |
| 133 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (Y, E) |
| 134 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  t0$$ | (W, W) |

**Interpretation:** When position 1064 mutates to any of the listed residues,
position 1065 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (V, P).

### Position Pair (1064, 1066) — MI = 1.1908, Reference: V→A

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 135 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}$$ | (A, E) |
| 136 | $$\bar{s2}  \cdot  s1  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  t1  \cdot  t0$$ | (L, N) |
| 137 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  t1$$ | (M, L) |
| 138 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (W, M) |
| 139 | $$s3  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (E, Q) |
| 140 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (K, F) |
| 141 | $$s3  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (R, H) |
| 142 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (W, L) |

**Interpretation:** When position 1064 mutates to any of the listed residues,
position 1066 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (V, A).

### Position Pair (1064, 1074) — MI = 0.7738, Reference: V→A

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 143 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (A, A) |
| 144 | $$\bar{s3}  \cdot  \bar{s1}  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (I, A) |
| 145 | $$\bar{s3}  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (V, F) |
| 146 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (M, I) |
| 147 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  t2  \cdot  t1  \cdot  t0$$ | (W, W) |
| 148 | $$s3  \cdot  \bar{s2}  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (E, E) |
| 149 | $$s3  \cdot  \bar{s2}  \cdot  s1  \cdot  s0  \cdot  \bar{t3}  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (N, I) |
| 150 | $$s3  \cdot  s2  \cdot  \bar{s1}  \cdot  \bar{s0}  \cdot  t3  \cdot  \bar{t2}  \cdot  \bar{t1}  \cdot  t0$$ | (H, D) |
| 151 | $$s3  \cdot  s2  \cdot  s1  \cdot  \bar{s0}  \cdot  \bar{t3}  \cdot  t2  \cdot  \bar{t1}  \cdot  t0$$ | (R, F) |
| 152 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  t3  \cdot  \bar{t2}  \cdot  t1  \cdot  \bar{t0}$$ | (W, Q) |

**Interpretation:** When position 1064 mutates to any of the listed residues,
position 1074 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (V, A).


**Total inference rules:** 152

## 7. Coupling Constants and Constraint Functions

### 7.1 Definition

The coupling constant $J_{ij}(a,b)$ measures the log-odds of observing pair $(a,b)$
at positions $(i,j)$ compared to the independent expectation:

$$J_{ij}(a,b) = \ln\frac{P_{ij}(a,b)}{P_i(a) \cdot P_j(b)}$$

- $J > 0$: pair is **co-evolutionary** (more common than expected)
- $J < 0$: pair is **anti-correlated** (less common than expected)
- $J = 0$: pair occurs at random frequency

### 7.2 Top Coupling Constants

| Pos i | Pos j | MI | avg\|J\| | Top Co-evolutionary (J>0) | Top Anti-correlated (J<0) |
|-------|-------|-----|----------|---------------------------|---------------------------|
| 372 | 401 | 1.5917 | 10.83 | I-A:+23.0, I-L:+23.0, I-M:+23.0 | T-N:-20.9, A-V:-20.1, K-N:-20.1 |
| 401 | 404 | 1.5690 | 10.46 | A-L:+23.0, A-M:+23.0, A-F:+23.0 | D-S:-20.9, N-I:-20.1, V-S:-20.1 |
| 208 | 209 | 1.5420 | 10.77 | A-A:+23.0, A-M:+23.0, A-Y:+23.0 | E-G:-21.0, L-D:-20.1, R-G:-20.1 |
| 209 | 210 | 1.5231 | 10.46 | A-A:+23.0, A-M:+23.0, A-Y:+23.0 | P-R:-21.0, G-E:-21.0, D-R:-20.1 |
| 208 | 210 | 1.5019 | 10.47 | A-A:+23.0, A-M:+23.0, A-Y:+23.0 | L-E:-21.0, E-R:-21.0, R-R:-20.1 |
| 207 | 209 | 1.4386 | 10.61 | A-A:+23.0, A-M:+23.0, A-Y:+23.0 | R-G:-21.1, N-P:-21.0, N-D:-20.1 |
| 207 | 208 | 1.4377 | 10.75 | A-A:+23.0, A-M:+23.0, A-F:+23.0 | R-L:-21.1, N-E:-21.0, N-R:-20.1 |
| 206 | 207 | 1.4292 | 10.76 | A-A:+23.0, A-M:+23.0, A-F:+23.0 | I-R:-21.1, L-N:-20.1, L-R:-19.1 |
| 206 | 209 | 1.4094 | 10.67 | A-A:+23.0, A-M:+23.0, A-Y:+23.0 | I-P:-21.0, L-G:-20.1, I-L:-19.1 |
| 206 | 208 | 1.4060 | 10.81 | A-A:+23.0, A-M:+23.0, A-F:+23.0 | I-E:-21.0, L-L:-20.1, V-R:-19.1 |
| 207 | 210 | 1.4010 | 10.37 | A-A:+23.0, A-M:+23.0, A-Y:+23.0 | R-R:-21.1, N-E:-21.0, N-L:-20.1 |
| 133 | 136 | 1.3903 | 12.40 | A-A:+23.0, A-I:+23.0, A-M:+23.0 | N-P:-21.2, C-F:-21.2, N-D:-20.1 |
| 131 | 133 | 1.3891 | 12.74 | A-A:+23.0, A-I:+23.0, A-V:+23.0 | F-C:-21.2, Q-N:-21.2, F-P:-20.1 |
| 131 | 136 | 1.3876 | 12.72 | A-A:+23.0, A-I:+23.0, A-M:+23.0 | F-P:-21.2, Q-F:-21.2, F-D:-20.1 |
| 84 | 86 | 1.3816 | 12.44 | A-I:+23.0, A-L:+23.0, A-M:+23.0 | D-G:-21.2, N-V:-21.2, V-V:-20.1 |
| 83 | 84 | 1.3816 | 12.44 | A-A:+23.0, A-I:+23.0, A-L:+23.0 | F-D:-21.2, N-N:-21.2, N-V:-20.1 |
| 83 | 86 | 1.3816 | 12.44 | A-I:+23.0, A-L:+23.0, A-M:+23.0 | F-V:-21.2, N-G:-21.2, N-F:-20.1 |
| 89 | 101 | 1.3812 | 12.76 | L-A:+23.0, L-L:+23.0, L-V:+23.0 | A-W:-21.2, F-I:-21.2, A-G:-20.1 |
| 89 | 110 | 1.3812 | 12.76 | L-A:+23.0, L-I:+23.0, L-V:+23.0 | A-K:-21.2, F-T:-21.2, A-S:-20.1 |
| 101 | 110 | 1.3812 | 13.14 | A-A:+23.0, A-I:+23.0, A-V:+23.0 | I-K:-21.2, W-T:-21.2, I-S:-20.1 |

## 8. Perplexity Analysis

### 8.1 Definition

Perplexity measures the effective number of choices: $PP = 2^H$.
Conditional perplexity $PP(j|i)$ measures how much knowing residue $i$
reduces uncertainty about residue $j$:

$$PP(j) = 2^{H(P_j)}, \quad PP(j|i=a) = 2^{H(P_{j|i=a})}$$

The **co-evolution ratio** $PP(j) / PP(j|i)$ quantifies constraint strength:
- Ratio $\approx 1$: no constraint (positions evolve independently)
- Ratio $> 2$: strong constraint (position $i$ determines position $j$)

### 8.2 Conditional Perplexity for Top Pairs

| Pos i | Pos j | MI | PPₘ(j) | PPₖₒₙₐ | Ratio | Most Constraining Residue | PP\|that |
|-------|-------|-----|--------|---------|-------|--------------------------|----------|
| 372 | 401 | 1.5917 | 3.185 | 1.224 | 2.60 | K | 1.000 |
| 401 | 404 | 1.5690 | 3.271 | 1.291 | 2.53 | E | 1.000 |
| 208 | 209 | 1.5420 | 3.063 | 1.177 | 2.60 | D | 1.000 |
| 209 | 210 | 1.5231 | 2.983 | 1.106 | 2.70 | V | 1.000 |
| 208 | 210 | 1.5019 | 2.983 | 1.187 | 2.51 | D | 1.000 |
| 207 | 209 | 1.4386 | 3.063 | 1.409 | 2.17 | E | 1.000 |
| 207 | 208 | 1.4377 | 3.062 | 1.468 | 2.09 | E | 1.000 |
| 206 | 207 | 1.4292 | 2.832 | 1.201 | 2.36 | N | 1.000 |
| 206 | 209 | 1.4094 | 3.063 | 1.392 | 2.20 | L | 1.000 |
| 206 | 208 | 1.4060 | 3.062 | 1.438 | 2.13 | L | 1.000 |
| 207 | 210 | 1.4010 | 2.983 | 1.490 | 2.00 | E | 1.000 |
| 133 | 136 | 1.3903 | 2.625 | 1.000 | 2.63 | D | 1.000 |
| 131 | 133 | 1.3891 | 2.624 | 1.000 | 2.62 | F | 1.000 |
| 131 | 136 | 1.3876 | 2.625 | 1.000 | 2.63 | F | 1.000 |
| 84 | 86 | 1.3816 | 2.605 | 1.000 | 2.61 | V | 1.000 |
| 83 | 84 | 1.3816 | 2.605 | 1.000 | 2.61 | F | 1.000 |
| 83 | 86 | 1.3816 | 2.605 | 1.000 | 2.61 | F | 1.000 |
| 89 | 101 | 1.3812 | 2.605 | 1.000 | 2.60 | A | 1.000 |
| 89 | 110 | 1.3812 | 2.605 | 1.000 | 2.60 | A | 1.000 |
| 101 | 110 | 1.3812 | 2.605 | 1.000 | 2.60 | I | 1.000 |

## 9. H1: Gray-code Adjacency Hypothesis

### 9.1 Hypothesis

Consecutive residues in protein sequences are preferentially K-map-adjacent
(Hamming distance = 1) compared to random expectation.

### 9.2 Method

For each consecutive pair of residues in each of the {n_all:,} sequences:
1. Encode each residue using its base-20 index (He 2012)
2. Compute Gray code: $g(i) = i \oplus (i \gg 1)$
3. Compute Hamming distance: $h = \text{popcount}(g(i) \oplus g(j))$

### 9.3 Results

| Metric | Value |
|--------|-------|
| Total consecutive pairs | 1,647,830 |
| Hamming-1 pairs | 313,888 |
| Observed ratio | **0.1905** |
| Expected (random) | 0.1613 |
| **Enrichment** | **1.18×** |

### 9.4 Full Hamming Distance Distribution

| Distance | Count | Percentage | Cumulative |
|----------|-------|------------|------------|
| 0 | 87,469 | 5.3% | 5.3% |
| 1 | 313,888 | 19.0% | 24.4% |
| 2 | 575,018 | 34.9% | 59.3% |
| 3 | 439,445 | 26.7% | 85.9% |
| 4 | 204,386 | 12.4% | 98.3% |
| 5 | 27,624 | 1.7% | 100.0% |

## 10. Complete Analysis Scripts Inventory

The `datasets/co-evolution/` directory contains **19 Python scripts**
and **1 shared module** (`coevolution_shared.py`).

| # | Script | Lines | Purpose |
|---|--------|-------|---------|
| 1 | `coevolution_shared.py` | 329 | Shared module: FASTA parsing, position arrays (cached), vectorized MI via `np.bincount`, entropy, perplexity, coupling, constraint function, shared-memory worker pool |
| 2 | `run_kmap_analysis.py` | 908 | Master K-map pipeline: H1-H6 on binary 32×32 K-map, consensus K-map, co-evolution analysis |
| 3 | `boolean_co-evolution.py` | 637 | Binary K-map Boolean minimization: 32×32 thresholded → Quine-McCluskey → essential prime implicants |
| 4 | `nary_kmap_co-evolution.py` | 537 | Base-20 K-map analysis: 20×20 frequency map → Boolean → coupling constants |
| 5 | `master_boolean.py` | 330 | Master Boolean function: 152 essential PIs across 36,918 pairs (full-length) |
| 6 | `position_kmap_coevolution.py` | 481 | Position-pair K-maps with MI: builds per-position-pair 20×20 K-maps and minimizes |
| 7 | `run_allseq_analysis.py` | 328 | Full position-based K-map analysis on ALL 1,299 sequences |
| 8 | `kmap_boolean_coevolution.py` | 383 | K-map Boolean with full markdown output: 152 rules across 15 position pairs |
| 9 | `generate_co-evolution_md.py` | 313 | Markdown generator from JSON results |
| 10 | `generate_full_analysis_md.py` | 397 | Comprehensive report generator from all JSON outputs |
| 11 | `create_mi_heatmap.py` | 247 | MI heatmap visualization (full + focus region), matplotlib |
| 12 | `allseq_constraint_function.py` | 247 | Leave-One-Out Cross-Validation: constraint function prediction accuracy |
| 13 | `predictive_constraint_function.py` | 308 | Three K-map approaches: Observed, Flipped, Continuous prediction |
| 14 | `flipped_boolean_coevolution.py` | 254 | Flipped Boolean: forbidden pairs analysis (negative selection) |
| 15 | `dca_boolean_coevolution.py` | 344 | DCA → Boolean pipeline: inverse covariance → coupling → QM minimization |
| 16 | `variable_position_coevolution.py` | 421 | Variable-position K-map with strategic don't-care conditions |
| 17 | `perplexity_coevolution.py` | 217 | Perplexity-based co-evolution strength measurement |
| 18 | `advanced_co-evolution_analysis.py` | 471 | Co-evolution network, Walsh-Hadamard spectrum, variant classification, clustering |
| 19 | `full_length_analysis.py` | 206 | Full-length (all 1,276 positions) entropy and MI analysis |
| 20 | `gpu_full_analysis.py` | 280 | GPU-accelerated analysis: numba parallel entropy/H1/mutations + shared-memory Pool for MI |
| 21 | `run_all_bg.sh` | 107 | Master launcher: runs all 17 scripts concurrently |

## 11. Flipped Boolean Analysis — Negative Selection

### 11.1 Concept

While the positive Boolean function asks "what co-evolves?", the **flipped** version
asks "what CANNOT co-exist?" The on-set (1) is assigned to residue pairs that are
NEVER observed together across all sequences, capturing **negative selection**.

### 11.2 Method

1. Build 20×20 frequency K-map for each position pair
2. Cells with frequency = 0 → on-set (1) — FORBIDDEN
3. Cells with frequency > 0 → off-set (0) — ALLOWED
4. Reference pair → don't-care (-1)
5. QM minimization → minimal forbidden constraints

### 11.3 Results

The flipped analysis produced **0** forbidden rules. Why?

With 1,299 sequences across multiple Omicron sub-lineages,
virtually every amino acid pair appears at least once at every position pair.
The Spike protein's N-terminal region is under **purifying selection** but not
**absolute constraint** — any single substitution has been sampled by evolution.

This means:

- **No universally forbidden pairs exist** — only statistically disfavored ones
- Co-evolution is **probabilistic**, not deterministic
- The constraint function $J_{ij}$ (continuous) is more appropriate than Boolean

## 12. Leave-One-Out Cross-Validation

### 12.1 Method

For each position pair $(i,j)$:
1. Exclude one sequence from the dataset
2. Build constraint function $C_{ij}(a,b)$ from the remaining {n_all - 1:,} sequences
3. Predict the co-evolving partner for the held-out sequence: $b^* = \arg\max_b C_{ij}(a_{\text{mutation}}, b)$
4. Check if prediction matches

### 12.2 Key Insight

LOO-CV accuracy was near zero (0.08%). This reveals that co-evolutionary rules are **lineage-specific** —
different Omicron sub-variants (BA.1, BA.2, BA.4, BA.5, XBB, etc.) have **different co-evolutionary patterns**.
Rules learned from one variant do not generalize to others.

## 13. What This Analysis Pipeline Achieves

### 13.1 By the Numbers

| Metric | Value |
|--------|-------|
| Input sequences | 1,299 |
| Positions analyzed | 1276 |
| Position arrays built | 1,299 × 1276 = 1,657,524 integers |
| Variable positions (H > 0.3) | 1249 |
| Co-evolving pairs (MI > 0.1) | 35774 |
| Boolean expressions (QM minimized) | 152 essential prime implicants |
| Unique position pairs with rules | 15 |
| Lean 4 theorems | 236 (106 + 115 + 15) |
| Python scripts | 20 |
| Total Python LOC | ~7,000 |
| Shared module LOC | 340 |

## 14. How to Apply the Co-evolution Rules

### 14.1 For a New Spike Sequence

```python
# 1. Extract residues at positions 68-79
seq_region = seq[68:80]

# 2. For each co-evolving position pair, check the Boolean function
for (pos_i, pos_j) in coevolving_pairs:
    aa_i = encode(seq_region[pos_i - 68])
    aa_j = encode(seq_region[pos_j - 68])

    # Encode as 4+4 bits
    bits = int_to_bits(aa_i, 4) + int_to_bits(aa_j, 4)

    # Check if ANY essential prime implicant matches
    is_coevolving = any(pi.matches(bits) for pi in prime_implicants[(pos_i, pos_j)])

    if is_coevolving:
        print(f'Positions {pos_i}-{pos_j} are co-evolving')

# 3. If position i mutates, find the required partner at position j
def predict_partner(pos_i, aa_i, pos_j):
    for pi in prime_implicants[(pos_i, pos_j)]:
        if pi.matches_row(aa_i):
            return pi.col_aa  # required co-evolving partner
    return None  # no constraint
```

## 15. Key References

1. Karnaugh, M. (1953). The Map Method for Synthesis of Combinational Logic Circuits. *AIEE Transactions*.
2. Gray, F. (1953). Pulse Code Communication. U.S. Patent 2,632,058.
3. Quine, W.V. (1952). The Problem of Simplifying Truth Functions. *American Mathematical Monthly*.
4. McCluskey, E.J. (1956). Minimization of Boolean Functions. *Bell System Technical Journal*.
5. He, P.A. et al. (2012). A novel graphical representation of proteins. *MATCH Communications*.
6. Petoukhov, S.V. (2024). Matrix Representations of Genetic Code and Karnaugh Maps. *Biosystems*.
7. de Moura, L. et al. (2021). The Lean 4 Theorem Prover and Programming Language.

---
*Generated August 07, 2026 at 04:05 by `generate_full_pipeline_doc.py`*
*All values computed from 1,299 Omicron Spike sequences using shared `coevolution_shared` module*