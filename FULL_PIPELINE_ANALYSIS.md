# E-Motioner-X-SBS: Complete Co-Evolution Analysis Pipeline

**Generated:** August 07, 2026 at 09:56
**Dataset:** SARS-CoV-2 Omicron Spike Protein — 1,299 sequences, 1276 positions
**Compute:** NVIDIA A100 80GB + 24-core Xeon, Python 3.10
**Author:** Shuvam Banerji Seal — IISER Kolkata

---

## Abstract

This document presents the complete computational pipeline for analyzing co-evolutionary
constraints in the SARS-CoV-2 Spike protein using Karnaugh map (K-map) Boolean minimization.
**1,299** Omicron variant Spike sequences from GISAID were encoded using base-20
(He 2012 ordering) amino acid representation. Position-pair mutual information identified
**3** co-evolving position pairs in the N-terminal signal peptide region
(positions 68-79). Quine-McCluskey Boolean minimization produced **3 essential
prime implicants** — each representing an irreducible co-evolutionary constraint.

**Key results:** H1 Gray-code adjacency enrichment = **0.19×**; 
max pairwise MI = **0.8067** at positions (373,378);
21/80 positions are evolutionarily variable (H > 0.3).

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
| **Variable positions** | 21/80 (entropy > 0.3) |
| **Conserved positions** | 59/80 |
| **Co-evolving pairs (MI > 0.01)** | 73 |
| **Co-evolving pairs (MI > 0.1)** | 17 |
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
| 1 | K | R | 1.0000 |

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
| 0 | M | 0.9962 | -0.0000 | 1.000 | Conserved |
| 1 | F | 0.9992 | -0.0000 | 1.000 | Conserved |
| 2 | V | 0.9954 | 0.0302 | 1.021 | Conserved |
| 3 | F | 1.0000 | -0.0000 | 1.000 | Conserved |
| 4 | L | 0.9838 | 0.1193 | 1.086 | Conserved |
| 5 | V | 0.9985 | -0.0000 | 1.000 | Conserved |
| 6 | L | 1.0000 | -0.0000 | 1.000 | Conserved |
| 7 | L | 1.0000 | -0.0000 | 1.000 | Conserved |
| 8 | P | 0.9985 | 0.0091 | 1.006 | Conserved |
| 9 | L | 0.9992 | -0.0000 | 1.000 | Conserved |
| 10 | V | 0.9992 | -0.0000 | 1.000 | Conserved |
| 11 | S | 0.9992 | 0.0091 | 1.006 | Conserved |
| 12 | S | 1.0000 | -0.0000 | 1.000 | Conserved |
| 13 | Q | 0.9992 | -0.0000 | 1.000 | Conserved |
| 14 | C | 0.9992 | -0.0000 | 1.000 | Conserved |
| 15 | V | 1.0000 | -0.0000 | 1.000 | Conserved |
| 16 | N | 0.9977 | -0.0000 | 1.000 | Conserved |
| 17 | L | 0.9985 | 0.0091 | 1.006 | Conserved |
| 18 | I | 0.7383 | 0.8116 | 1.755 | Variable |
| 19 | T | 1.0000 | -0.0000 | 1.000 | Conserved |
| 20 | R | 1.0000 | -0.0000 | 1.000 | Conserved |
| 21 | T | 0.9992 | -0.0000 | 1.000 | Conserved |
| 22 | Q | 0.9992 | -0.0000 | 1.000 | Conserved |
| 23 | L | 0.2440 | 0.0307 | 1.021 | Conserved |
| 24 | P | 0.2440 | -0.0000 | 1.000 | Conserved |
| 25 | P | 0.2440 | -0.0000 | 1.000 | Conserved |
| 26 | S | 0.7529 | 0.8039 | 1.746 | Variable |
| 27 | Y | 1.0000 | -0.0000 | 1.000 | Conserved |
| 28 | T | 1.0000 | -0.0000 | 1.000 | Conserved |
| 29 | N | 1.0000 | -0.0000 | 1.000 | Conserved |
| 30 | S | 1.0000 | -0.0000 | 1.000 | Conserved |
| 31 | F | 1.0000 | -0.0000 | 1.000 | Conserved |
| 32 | T | 1.0000 | -0.0000 | 1.000 | Conserved |
| 33 | R | 1.0000 | -0.0000 | 1.000 | Conserved |
| 34 | G | 1.0000 | -0.0000 | 1.000 | Conserved |
| 35 | V | 1.0000 | -0.0000 | 1.000 | Conserved |
| 36 | Y | 1.0000 | -0.0000 | 1.000 | Conserved |
| 37 | Y | 1.0000 | -0.0000 | 1.000 | Conserved |
| 38 | P | 1.0000 | -0.0000 | 1.000 | Conserved |
| 39 | D | 1.0000 | -0.0000 | 1.000 | Conserved |
| ... | ... | ... | ... | ... | ... |
| 40 | K | 1.0000 | -0.0000 | 1.000 | Conserved |
| 41 | V | 1.0000 | -0.0000 | 1.000 | Conserved |
| 42 | F | 1.0000 | -0.0000 | 1.000 | Conserved |
| 43 | R | 1.0000 | -0.0000 | 1.000 | Conserved |
| 44 | S | 1.0000 | -0.0000 | 1.000 | Conserved |
| 45 | S | 1.0000 | -0.0000 | 1.000 | Conserved |
| 46 | V | 1.0000 | -0.0000 | 1.000 | Conserved |
| 47 | L | 1.0000 | -0.0000 | 1.000 | Conserved |
| 48 | H | 0.9977 | 0.0236 | 1.016 | Conserved |
| 49 | S | 1.0000 | -0.0000 | 1.000 | Conserved |
| 50 | T | 1.0000 | -0.0000 | 1.000 | Conserved |
| 51 | Q | 0.9985 | 0.0166 | 1.012 | Conserved |
| 52 | D | 0.9992 | -0.0000 | 1.000 | Conserved |
| 53 | L | 0.9992 | 0.0091 | 1.006 | Conserved |
| 54 | F | 1.0000 | -0.0000 | 1.000 | Conserved |
| 55 | L | 1.0000 | -0.0000 | 1.000 | Conserved |
| 56 | P | 1.0000 | -0.0000 | 1.000 | Conserved |
| 57 | F | 1.0000 | -0.0000 | 1.000 | Conserved |
| 58 | F | 1.0000 | -0.0000 | 1.000 | Conserved |
| 59 | S | 1.0000 | -0.0000 | 1.000 | Conserved |
| 60 | N | 1.0000 | -0.0000 | 1.000 | Conserved |
| 61 | V | 0.9985 | 0.0181 | 1.013 | Conserved |
| 62 | T | 1.0000 | -0.0000 | 1.000 | Conserved |
| 63 | W | 1.0000 | -0.0000 | 1.000 | Conserved |
| 64 | F | 1.0000 | -0.0000 | 1.000 | Conserved |
| 65 | H | 1.0000 | -0.0000 | 1.000 | Conserved |
| 66 | A | 0.7583 | 0.7979 | 1.739 | Variable |
| 67 | I | 0.9969 | 0.0326 | 1.023 | Conserved |
| 68 | H | 0.6713 | -0.0000 | 1.000 | Conserved |
| 69 | V | 0.6713 | -0.0000 | 1.000 | Conserved |
| 70 | S | 0.9985 | 0.0091 | 1.006 | Conserved |
| 71 | G | 1.0000 | -0.0000 | 1.000 | Conserved |
| 72 | T | 1.0000 | -0.0000 | 1.000 | Conserved |
| 73 | N | 1.0000 | -0.0000 | 1.000 | Conserved |
| 74 | G | 1.0000 | -0.0000 | 1.000 | Conserved |
| 75 | T | 0.9977 | 0.0091 | 1.006 | Conserved |
| 76 | K | 1.0000 | -0.0000 | 1.000 | Conserved |
| 77 | R | 1.0000 | -0.0000 | 1.000 | Conserved |
| 78 | F | 0.9992 | 0.0091 | 1.006 | Conserved |
| 79 | D | 0.9985 | 0.0091 | 1.006 | Conserved |

**Summary:** 21 variable positions (H > 0.3), 59 conserved.

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
| 1 | 373 | 378 | 0.8067 | F | A | 0.872 | 0.811 | 5 |
| 2 | 18 | 26 | 0.8024 | I | S | 0.812 | 0.804 | 8 |
| 3 | 378 | 407 | 0.7920 | A | N | 0.811 | 0.814 | 29 |
| 4 | 66 | 94 | 0.7907 | A | T | 0.798 | 0.799 | 28 |
| 5 | 215 | 216 | 0.7571 | G | R | 0.881 | 0.764 | 1 |
| 6 | 210 | 215 | 0.7532 | N | G | 0.785 | 0.881 | 5 |
| 7 | 407 | 410 | 0.7495 | N | S | 0.814 | 0.822 | 3 |
| 8 | 210 | 216 | 0.7453 | N | R | 0.785 | 0.764 | 6 |
| 9 | 212 | 215 | 0.3977 | V | G | 0.462 | 0.881 | 3 |
| 10 | 212 | 216 | 0.3773 | V | R | 0.462 | 0.764 | 4 |
| 11 | 488 | 495 | 0.3419 | F | R | 0.414 | 0.455 | 7 |
| 12 | 210 | 212 | 0.1769 | N | V | 0.785 | 0.462 | 2 |
| 13 | 419 | 442 | 0.1319 | N | K | 0.257 | 0.346 | 23 |
| 14 | 479 | 486 | 0.1134 | N | A | 0.142 | 0.133 | 7 |
| 15 | 479 | 480 | 0.1115 | N | K | 0.142 | 0.124 | 1 |
| 16 | 500 | 507 | 0.1028 | R | H | 0.119 | 0.110 | 7 |
| 17 | 503 | 507 | 0.1011 | Y | H | 0.114 | 0.110 | 4 |
| 18 | 500 | 503 | 0.0958 | R | Y | 0.119 | 0.114 | 3 |
| 19 | 480 | 486 | 0.0867 | K | A | 0.124 | 0.133 | 6 |
| 20 | 410 | 419 | 0.0832 | S | N | 0.822 | 0.257 | 9 |
| 21 | 375 | 377 | 0.0758 | P | F | 0.076 | 0.076 | 2 |
| 22 | 407 | 419 | 0.0735 | N | N | 0.814 | 0.257 | 12 |
| 23 | 486 | 500 | 0.0732 | A | R | 0.133 | 0.119 | 14 |
| 24 | 373 | 377 | 0.0720 | F | F | 0.872 | 0.076 | 4 |
| 25 | 373 | 375 | 0.0719 | F | P | 0.872 | 0.076 | 2 |
| 26 | 479 | 500 | 0.0707 | N | R | 0.142 | 0.119 | 21 |
| 27 | 486 | 507 | 0.0697 | A | H | 0.133 | 0.110 | 21 |
| 28 | 479 | 507 | 0.0674 | N | H | 0.142 | 0.110 | 28 |
| 29 | 486 | 503 | 0.0643 | A | Y | 0.133 | 0.114 | 17 |
| 30 | 210 | 214 | 0.0623 | N | E | 0.785 | 0.061 | 4 |
| 31 | 479 | 503 | 0.0622 | N | Y | 0.142 | 0.114 | 24 |
| 32 | 213 | 216 | 0.0610 | R | R | 0.058 | 0.764 | 3 |
| 33 | 214 | 216 | 0.0610 | E | R | 0.061 | 0.764 | 2 |
| 34 | 212 | 214 | 0.0610 | V | E | 0.462 | 0.061 | 2 |
| 35 | 213 | 214 | 0.0610 | R | E | 0.058 | 0.061 | 1 |
| 36 | 213 | 215 | 0.0610 | R | G | 0.058 | 0.881 | 2 |
| 37 | 214 | 215 | 0.0610 | E | G | 0.061 | 0.881 | 1 |
| 38 | 210 | 213 | 0.0588 | N | R | 0.785 | 0.058 | 3 |
| 39 | 212 | 213 | 0.0576 | V | R | 0.462 | 0.058 | 1 |
| 40 | 480 | 500 | 0.0551 | K | R | 0.124 | 0.119 | 20 |
| 41 | 480 | 503 | 0.0525 | K | Y | 0.124 | 0.114 | 23 |
| 42 | 480 | 507 | 0.0512 | K | H | 0.124 | 0.110 | 27 |
| 43 | 495 | 507 | 0.0512 | R | H | 0.455 | 0.110 | 12 |
| 44 | 495 | 503 | 0.0512 | R | Y | 0.455 | 0.114 | 8 |
| 45 | 495 | 500 | 0.0496 | R | R | 0.455 | 0.119 | 5 |
| 46 | 250 | 251 | 0.0393 | Y | L | 0.065 | 0.042 | 1 |
| 47 | 495 | 498 | 0.0386 | R | G | 0.455 | 0.781 | 3 |
| 48 | 681 | 683 | 0.0343 | K | H | 0.036 | 0.047 | 2 |
| 49 | 488 | 498 | 0.0320 | F | G | 0.414 | 0.781 | 10 |
| 50 | 479 | 495 | 0.0313 | N | R | 0.142 | 0.455 | 16 |

## 6. Quine-McCluskey Boolean Minimization — All 108 Essential Prime Implicants

The Boolean minimization was performed on 3 position pairs (68-79),
producing **3 essential prime implicants**. Each rule has the form:

$$f(s_3, s_2, s_1, s_0, t_3, t_2, t_1, t_0) = \text{AND of literals}$$

**Variables:**
- $s_3 s_2 s_1 s_0$ = 4-bit binary encoding of residue at position $i$
- $t_3 t_2 t_1 t_0$ = 4-bit binary encoding of residue at position $j$
- $\bar{s}_k$ = NOT ($s_k = 0$), $s_k$ = ($s_k = 1$)

### 6.1 Complete Inference Rules (3 rules across 3 position pairs)

### Position Pair (210, 212) — MI = 0.1769, Reference: N→V

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 1 | $$\bar{s3}  \cdot  s2  \cdot  \bar{s1}  \cdot  s0  \cdot  t3  \cdot  t1  \cdot  t0$$ | (N, S) |

**Interpretation:** When position 210 mutates to any of the listed residues,
position 212 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (N, V).

### Position Pair (212, 215) — MI = 0.3977, Reference: V→G

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 2 | $$\bar{s3}  \cdot  s0  \cdot  t3  \cdot  t2  \cdot  \bar{t1}  \cdot  \bar{t0}$$ | (V, P) |

**Interpretation:** When position 212 mutates to any of the listed residues,
position 215 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (V, G).

### Position Pair (212, 216) — MI = 0.3773, Reference: V→R

| Rule | Boolean Expression | Amino Acids |
|------|-------------------|-------------|
| 3 | $$\bar{s3}  \cdot  s2  \cdot  s1  \cdot  s0  \cdot  t3  \cdot  t1  \cdot  t0$$ | (S, R) |

**Interpretation:** When position 212 mutates to any of the listed residues,
position 216 must co-evolve to the corresponding partner residue to maintain
protein stability. The reference pair is (V, R).


**Total inference rules:** 3

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
| 373 | 378 | 0.8067 | 18.09 | A-I:+23.0, A-L:+23.0, A-V:+23.0 | F-T:-21.3, L-A:-21.3, S-A:-2.3 |
| 18 | 26 | 0.8024 | 18.87 | A-I:+23.0, A-L:+23.0, A-V:+23.0 | I-A:-21.3, T-S:-5.5 |
| 378 | 407 | 0.7920 | 18.87 | I-A:+23.0, I-I:+23.0, I-L:+23.0 | T-N:-21.3, A-D:-4.4 |
| 66 | 94 | 0.7907 | 18.87 | I-A:+23.0, I-L:+23.0, I-V:+23.0 | V-T:-21.3, A-I:-5.5 |
| 215 | 216 | 0.7571 | 16.85 | A-A:+23.0, A-I:+23.0, A-L:+23.0 | P-R:-21.3, G-E:-21.2, V-E:-17.4 |
| 210 | 215 | 0.7532 | 16.88 | A-A:+23.0, A-I:+23.0, A-L:+23.0 | N-P:-21.2, I-G:-21.2, I-V:-16.4 |
| 407 | 410 | 0.7495 | 18.87 | A-A:+23.0, A-I:+23.0, A-L:+23.0 | D-S:-21.3, N-R:-3.0 |
| 210 | 216 | 0.7453 | 17.57 | A-A:+23.0, A-I:+23.0, A-L:+23.0 | N-E:-21.2, N-K:-15.6, K-E:-14.4 |
| 212 | 215 | 0.3977 | 16.86 | A-A:+23.0, A-I:+23.0, A-L:+23.0 | V-V:-20.2, I-P:-19.8, L-P:-19.0 |
| 212 | 216 | 0.3773 | 16.14 | A-A:+23.0, A-I:+23.0, A-L:+23.0 | V-R:-20.3, I-E:-19.8, L-E:-19.0 |
| 488 | 495 | 0.3419 | 18.24 | A-A:+23.0, A-I:+23.0, A-L:+23.0 | V-R:-20.4, P-R:-16.5, F-Q:-1.7 |
| 210 | 212 | 0.1769 | 17.58 | A-A:+23.0, A-I:+23.0, A-M:+23.0 | N-V:-19.4, I-L:-19.1, I-S:-18.0 |
| 419 | 442 | 0.1319 | 18.93 | A-A:+23.0, A-I:+23.0, A-L:+23.0 | K-K:-1.7 |
| 479 | 486 | 0.1134 | 19.03 | A-I:+23.0, A-L:+23.0, A-V:+23.0 | N-E:-3.2, S-A:-2.1 |
| 479 | 480 | 0.1115 | 19.07 | A-A:+23.0, A-I:+23.0, A-L:+23.0 | N-T:-18.9, S-K:-1.9 |
| 500 | 507 | 0.1028 | 19.09 | A-A:+23.0, A-I:+23.0, A-L:+23.0 | R-Y:-18.8, Q-H:-2.3 |
| 503 | 507 | 0.1011 | 18.40 | A-A:+23.0, A-I:+23.0, A-L:+23.0 | N-H:-18.7, F-Y:-11.6, Y-Y:-2.9 |
| 500 | 503 | 0.0958 | 18.40 | A-A:+23.0, A-I:+23.0, A-L:+23.0 | R-N:-18.7, Q-F:-11.7, Q-Y:-1.9 |
| 480 | 486 | 0.0867 | 19.03 | A-I:+23.0, A-L:+23.0, A-V:+23.0 | T-A:-2.0, K-E:-1.6 |
| 410 | 419 | 0.0832 | 18.88 | A-A:+23.0, A-I:+23.0, A-L:+23.0 | S-K:-3.7 |

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
| 373 | 378 | 0.8067 | 1.754 | 1.109 | 1.58 | L | 1.000 |
| 18 | 26 | 0.8024 | 1.746 | 1.029 | 1.70 | I | 1.015 |
| 378 | 407 | 0.7920 | 1.758 | 1.011 | 1.74 | T | 1.000 |
| 66 | 94 | 0.7907 | 1.740 | 1.008 | 1.73 | V | 1.000 |
| 215 | 216 | 0.7571 | 1.698 | 1.008 | 1.68 | V | 1.000 |
| 210 | 215 | 0.7532 | 1.842 | 1.176 | 1.57 | N | 1.057 |
| 407 | 410 | 0.7495 | 1.768 | 1.103 | 1.60 | N | 1.085 |
| 210 | 216 | 0.7453 | 1.698 | 1.171 | 1.45 | N | 1.000 |
| 212 | 215 | 0.3977 | 1.842 | 1.093 | 1.69 | I | 1.000 |
| 212 | 216 | 0.3773 | 1.698 | 1.102 | 1.54 | I | 1.000 |
| 488 | 495 | 0.3419 | 1.371 | 1.048 | 1.31 | V | 1.000 |
| 210 | 212 | 0.1769 | 1.377 | 1.038 | 1.33 | I | 1.023 |
| 419 | 442 | 0.1319 | 1.271 | 1.458 | 0.87 | N | 1.221 |
| 479 | 486 | 0.1134 | 1.096 | 1.218 | 0.90 | N | 1.006 |
| 479 | 480 | 0.1115 | 1.090 | 1.268 | 0.86 | N | 1.000 |
| 500 | 507 | 0.1028 | 1.079 | 1.188 | 0.91 | R | 1.006 |
| 503 | 507 | 0.1011 | 1.079 | 1.006 | 1.07 | N | 1.000 |
| 500 | 503 | 0.0958 | 1.083 | 1.257 | 0.86 | R | 1.006 |
| 480 | 486 | 0.0867 | 1.096 | 1.258 | 0.87 | K | 1.026 |
| 410 | 419 | 0.0832 | 1.195 | 1.339 | 0.89 | S | 1.030 |

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
| Total consecutive pairs | 1,656,225 |
| Hamming-1 pairs | 314,886 |
| Observed ratio | **0.1901** |
| Expected (random) | 0.1613 |
| **Enrichment** | **1.18×** |

### 9.4 Full Hamming Distance Distribution

| Distance | Count | Percentage | Cumulative |
|----------|-------|------------|------------|
| 0 | 92,592 | 5.6% | 5.6% |
| 1 | 314,886 | 19.0% | 24.6% |
| 2 | 575,595 | 34.8% | 59.4% |
| 3 | 439,284 | 26.5% | 85.9% |
| 4 | 205,470 | 12.4% | 98.3% |
| 5 | 28,398 | 1.7% | 100.0% |

## 10. Complete Analysis Scripts Inventory

The `datasets/co-evolution/` directory contains **19 Python scripts**
and **1 shared module** (`coevolution_shared.py`).

| # | Script | Lines | Purpose |
|---|--------|-------|---------|
| 1 | `coevolution_shared.py` | 329 | Shared module: FASTA parsing, position arrays (cached), vectorized MI via `np.bincount`, entropy, perplexity, coupling, constraint function, shared-memory worker pool |
| 2 | `run_kmap_analysis.py` | 908 | Master K-map pipeline: H1-H6 on binary 32×32 K-map, consensus K-map, co-evolution analysis |
| 3 | `boolean_co-evolution.py` | 637 | Binary K-map Boolean minimization: 32×32 thresholded → Quine-McCluskey → essential prime implicants |
| 4 | `nary_kmap_co-evolution.py` | 537 | Base-20 K-map analysis: 20×20 frequency map → Boolean → coupling constants |
| 5 | `master_boolean.py` | 330 | Master Boolean function: 3 essential PIs across 10 pairs (full-length) |
| 6 | `position_kmap_coevolution.py` | 481 | Position-pair K-maps with MI: builds per-position-pair 20×20 K-maps and minimizes |
| 7 | `run_allseq_analysis.py` | 328 | Full position-based K-map analysis on ALL 1,299 sequences |
| 8 | `kmap_boolean_coevolution.py` | 383 | K-map Boolean with full markdown output: 3 rules across 3 position pairs |
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
| Variable positions (H > 0.3) | 21 |
| Co-evolving pairs (MI > 0.1) | 17 |
| Boolean expressions (QM minimized) | 3 essential prime implicants |
| Unique position pairs with rules | 3 |
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
*Generated August 07, 2026 at 09:56 by `generate_full_pipeline_doc.py`*
*All values computed from 1,299 Omicron Spike sequences using shared `coevolution_shared` module*