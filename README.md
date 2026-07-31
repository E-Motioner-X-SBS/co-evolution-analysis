# Co-evolution Analysis of SARS-CoV-2 Spike Protein

Karnaugh map (K-map) Boolean minimization applied to predict co-evolutionary constraints in the SARS-CoV-2 Spike protein.

## Overview

This pipeline analyzes **1,299 SARS-CoV-2 Omicron Spike protein sequences** using:
- Base-20 amino acid encoding (He 2012 ordering)
- Position-pair mutual information to identify co-evolving positions
- Quine-McCluskey Boolean minimization to extract minimal co-evolutionary rules
- Constraint function for predicting co-evolutionary pairs

## Quick Start

```bash
# Install dependencies
pip install numpy scipy matplotlib

# Run the master Boolean function analysis
python master_boolean.py

# Generate comprehensive markdown with all rules
python generate_co-evolution_md.py

# Create MI heatmap
python create_mi_heatmap.py

# Run LOO-CV prediction test
python allseq_constraint_function.py

# Run all scripts at once
bash run_all_bg.sh
```

## Key Results

| Metric | Value |
|--------|-------|
| Sequences analyzed | 1,299 |
| Variable positions | 1,249 (H > 0.3) |
| Co-evolutionary pairs (MI > 0.1) | 36,918 |
| Essential prime implicants | 152 |
| LOO-CV accuracy | 7.26% |
| Flipped Boolean forbidden rules | 345 |
| Sequence clusters | 10 |

## Scripts

| # | Script | Purpose | Output |
|---|--------|---------|--------|
| 1 | `coevolution_shared.py` | Shared module (FASTA, MI, entropy) | — |
| 2 | `master_boolean.py` | Master Boolean function (152 rules) | `master_boolean/` |
| 3 | `boolean_co-evolution.py` | Binary K-map Boolean minimization | `boolean_results/` |
| 4 | `nary_kmap_co-evolution.py` | Base-20 K-map analysis | `nary_kmap_results/` |
| 5 | `position_kmap_coevolution.py` | Position-pair K-maps with MI | `position_kmap_results/` |
| 6 | `run_allseq_analysis.py` | Full analysis on all sequences | `full_position_results/` |
| 7 | `run_kmap_analysis.py` | Master K-map pipeline (H1-H6) | `kmap_results/` |
| 8 | `flipped_boolean_coevolution.py` | Forbidden pairs (negative selection) | `flipped_boolean_results/` |
| 9 | `kmap_boolean_coevolution.py` | K-map Boolean with markdown output | `kmap_boolean_coevolution/` |
| 10 | `variable_position_coevolution.py` | Variable-position K-map with don't-care | `variable_position_results/` |
| 11 | `predictive_constraint_function.py` | Three K-map approaches (obs/flipped/cont) | `constraint_function_results/` |
| 12 | `allseq_constraint_function.py` | LOO-CV constraint function | `allseq_constraint_results/` |
| 13 | `dca_boolean_coevolution.py` | Local precision matrix → Boolean | `dca_boolean_results/` |
| 14 | `perplexity_coevolution.py` | Perplexity-based co-evolution | `perplexity_results/` |
| 15 | `advanced_co-evolution_analysis.py` | Network, Walsh-Hadamard, clustering | `advanced_analysis_results/` |
| 16 | `full_length_analysis.py` | Full-length (1276 positions) entropy/MI | `full_length_results/` |
| 17 | `gpu_full_analysis.py` | GPU-accelerated analysis (numba) | `full_gpu_results/` |
| 18 | `create_mi_heatmap.py` | MI heatmap visualization | `mi_heatmap/` |
| 19 | `generate_co-evolution_md.py` | Markdown from master_boolean JSON | `COEVOLUTION_KMAP_BOOLEAN.md` |
| 20 | `generate_full_analysis_md.py` | Comprehensive report from all JSONs | `FULL_COEVOLUTION_ANALYSIS.md` |
| 21 | `generate_full_pipeline_doc.py` | Full pipeline documentation | `FULL_PIPELINE_ANALYSIS.md` |

## Methodology

### 1. Position Identification
Variable positions identified by Shannon entropy: `H(p) > 0.3`

### 2. Mutual Information
For position pair (i, j):
```
MI(i,j) = Σ_{x,y} P(x,y) log2(P(x,y) / (P(x)·P(y)))
```

### 3. K-map Construction
For each co-evolving position pair (i, j), build a 20×20 frequency matrix:
```
K_{ij}(a, b) = (1/N) Σ_s 1[seq_s[i] = a AND seq_s[j] = b]
```

### 4. Boolean Minimization (Quine-McCluskey)
Threshold the frequency matrix to create a Boolean function, then minimize:
```
f(pos_i, pos_j, aa_i, aa_j) = OR of all essential prime implicants
```

### 5. Constraint Function
```
C(aa_i, aa_j) = ln(P(aa_i, aa_j) / (P(aa_i)·P(aa_j)))
```
- C > 0: co-evolutionary (more common than expected)
- C < 0: anti-correlated (less common than expected)

### 6. Prediction
```
P_co-evolution = σ(C) = 1/(1+e^{-C})
```

## Dependencies

- Python 3.10+
- NumPy
- SciPy
- Matplotlib
- [kmap-sbm-validation](https://github.com/E-Motioner-X-SBS/kmap-sbm-validation) (for `gray_amino.py` encoding)
- [n-ary-kmap](https://github.com/E-Motioner-X-SBS/n-ary-kmap) (for `bio_sequences.py` encoding)

## Repository Structure

```
co-evolution-analysis/
├── *.py                          # 21 analysis scripts
├── run_all_bg.sh                 # Master launcher
├── Spike_protein.aln-fasta       # Input: 1,299 Omicron sequences
├── *.md                          # Generated reports
├── master_boolean/               # Master Boolean function results
├── kmap_results/                 # K-map frequency analysis
├── boolean_results/              # Binary K-map Boolean results
├── nary_kmap_results/            # Base-20 K-map results
├── position_kmap_results/        # Position-pair K-map results
├── full_position_results/        # Full analysis results
├── flipped_boolean_results/      # Forbidden pairs results
├── constraint_function_results/  # Constraint function results
├── allseq_constraint_results/    # LOO-CV results
├── dca_boolean_results/          # Local precision matrix results
├── perplexity_results/           # Perplexity analysis results
├── advanced_analysis_results/    # Network, clustering, Walsh-Hadamard
├── full_length_results/          # Full-length analysis
├── full_gpu_results/             # GPU-accelerated results
├── kmap_boolean_coevolution/     # K-map Boolean with markdown
├── mi_heatmap/                   # MI heatmap visualizations
├── variable_position_results/    # Variable-position results
├── README.md
└── .gitignore
```

## Related Repositories

- [kmap-sbm-validation](https://github.com/E-Motioner-X-SBS/kmap-sbm-validation) — Binary K-map validation via SBM MD
- [n-ary-kmap](https://github.com/E-Motioner-X-SBS/n-ary-kmap) — Base-N K-map generalization
- [lean_proofs](https://github.com/E-Motioner-X-SBS/lean_proofs) — Lean 4 formal proofs (106 theorems)
- [datasets](https://github.com/E-Motioner-X-SBS/datasets) — Curated PDB structures

## Citation

```bibtex
@software{coevolution_analysis,
  title = {Co-evolution Analysis of SARS-CoV-2 Spike Protein via K-map Boolean Minimization},
  author = {Shuvam Banerji Seal},
  year = {2026},
  url = {https://github.com/E-Motioner-X-SBS/co-evolution-analysis}
}
```

## License

MIT

---

*"The pursuit of knowledge is the highest form of motion."*
