# Co-evolution Analysis of SARS-CoV-2 Spike Protein

> **CORRECTION (Aug 7, 2026):** two defects (A1 gap-stripping misalignment,
> A2 8-bit QM wrap-around) were fixed; all scripts re-run on aligned columns
> with padded 32x32 K-maps. Verified corrected results: 21 variable
> positions, 10 co-evolving pairs, 36 distinct Boolean rules (2 essential).
> See [My_Interpretation_of_the_Results/CORRECTION_NOTICE.md](My_Interpretation_of_the_Results/CORRECTION_NOTICE.md).

Karnaugh map (K-map) Boolean minimization applied to predict co-evolutionary constraints in the SARS-CoV-2 Spike protein.

**FULL-LENGTH + GPU:** All 20 scripts run on the complete sequence (1,276 positions) with all 1,299 sequences, GPU-accelerated via torch CUDA (A100). See [RESULTS_ANALYSIS.md](RESULTS_ANALYSIS.md) for the rigorous audit.

## Overview

This pipeline analyzes **1,299 SARS-CoV-2 Omicron Spike protein sequences** using:
- Base-20 amino acid encoding (He 2012 ordering)
- Position-pair mutual information to identify co-evolving positions
- Quine-McCluskey Boolean minimization to extract minimal co-evolutionary rules
- Constraint function for predicting co-evolutionary pairs
- **GPU acceleration** (`coevolution_gpu.py`): full 813K-pair MI matrix in 1.6s on A100

## Quick Start

```bash
# GPU environment (torch 2.12.1+cu130, numba 0.66.0)
PY=/store/shuvam/.venv/bin/python

# Run the master Boolean function analysis
$PY master_boolean.py

# Generate comprehensive markdown with all rules
$PY generate_co-evolution_md.py

# Create MI heatmap (GPU)
$PY create_mi_heatmap.py

# Run LOO-CV prediction test
$PY allseq_constraint_function.py

# GPU full analysis
nohup $PY -u gpu_full_analysis.py > logs/gpu_full.log 2>&1 &
```

## Key Results (FULL LENGTH, all 1,299 sequences)

| Metric | Value |
|--------|-------|
| Sequences analyzed | 1,299 |
| Positions analyzed | 1,276 (FULL length) |
| Variable positions | 21 (H > 0.3, corrected) |
| Co-evolutionary pairs (mutation-only MI > 0.1) | 10 |
| Max full-MI | 0.8067 @ (373, 378) |
| Max mutation-only MI | 0.8710 @ (495, 498) |
| Distinct Boolean rules (essential) | 36 (2) |
| Flipped Boolean forbidden rules | 490 |
| High-MI pairs (full length) | 5 (MI > 0.5) |
| LOO-CV accuracy | 9.24% (corrected) |
| DCA-style accuracy | 17.6% (local precision, NOT real DCA) |
| Variant signatures | 40 |
| GPU MI matrix speed | 813K pairs in 1.6 s (~800× vs CPU) |

## Scripts

| # | Script | Purpose | Output |
|---|--------|---------|--------|
| 1 | `coevolution_shared.py` | Shared module (FASTA, MI, entropy, GPU pair-finder) | — |
| 2 | `coevolution_gpu.py` | **GPU kernels** (torch CUDA: MI matrix, entropy, refs, coupling, H1) | — |
| 3 | `master_boolean.py` | Master Boolean function (36 rules, 2 essential) | `master_boolean/` |
| 4 | `boolean_co-evolution.py` | Binary K-map Boolean minimization (796,953 pairs) | `boolean_results/` |
| 5 | `nary_kmap_co-evolution.py` | Base-20 K-map analysis | `nary_kmap_results/` |
| 6 | `position_kmap_coevolution.py` | Position-pair K-maps with MI (25,199 pairs) | `position_kmap_results/` |
| 7 | `run_allseq_analysis.py` | Full analysis on all sequences | `full_position_results/` |
| 8 | `run_kmap_analysis.py` | Master K-map pipeline (H1-H6) | `kmap_results/` |
| 9 | `flipped_boolean_coevolution.py` | Forbidden pairs (negative selection, 345 rules) | `flipped_boolean_results/` |
| 10 | `kmap_boolean_coevolution.py` | K-map Boolean with markdown output | `kmap_boolean_coevolution/` |
| 11 | `variable_position_coevolution.py` | Variable-position K-map with don't-care | `variable_position_results/` |
| 12 | `predictive_constraint_function.py` | Constraint function train/test | `constraint_function_results/` |
| 13 | `allseq_constraint_function.py` | LOO-CV constraint function | `allseq_constraint_results/` |
| 14 | `dca_boolean_coevolution.py` | Local precision matrix → Boolean (NOT real DCA) | `dca_boolean_results/` |
| 15 | `perplexity_coevolution.py` | Perplexity-based co-evolution | `perplexity_results/` |
| 16 | `advanced_co-evolution_analysis.py` | Network, Walsh-Hadamard, clustering, signatures | `advanced_analysis_results/` |
| 17 | `full_length_analysis.py` | Full-length (1276 positions) entropy/MI | `full_length_results/` |
| 18 | `gpu_full_analysis.py` | GPU full analysis (torch CUDA, saves full MI matrix) | `full_gpu_results/` |
| 19 | `create_mi_heatmap.py` | MI heatmap visualization (GPU) | `mi_heatmap/` |
| 20 | `generate_co-evolution_md.py` | Markdown from master_boolean JSON | `COEVOLUTION_KMAP_BOOLEAN.md` |
| 21 | `generate_full_analysis_md.py` | Comprehensive report from all JSONs | `FULL_COEVOLUTION_ANALYSIS.md` |
| 22 | `generate_full_pipeline_doc.py` | Full pipeline documentation (full-length) | `FULL_PIPELINE_ANALYSIS.md` |

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
- NumPy, SciPy, Matplotlib
- **torch (CUDA)** — GPU kernels (`coevolution_gpu.py`); env: `/store/shuvam/.venv` (torch 2.12.1+cu130)
- numba (optional, for legacy parallel code)
- [kmap-sbm-validation](https://github.com/E-Motioner-X-SBS/kmap-sbm-validation) (for `gray_amino.py` encoding)
- [n-ary-kmap](https://github.com/E-Motioner-X-SBS/n-ary-kmap) (for `bio_sequences.py` encoding)

## Repository Structure

```
co-evolution-analysis/
├── *.py                          # 22 analysis scripts (incl. coevolution_gpu.py)
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
├── full_gpu_results/             # GPU results + full MI matrix (npy/csv)
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
