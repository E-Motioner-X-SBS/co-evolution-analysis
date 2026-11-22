# Co-Evolution Pipeline — Corrected Full-Length GPU Results & Rigorous Analysis

**Date:** Aug 7, 2026 (corrected) · **Last verified:** Aug 17, 2026
**Repo:** https://github.com/E-Motioner-X-SBS/co-evolution-analysis
**Dataset:** 1,299 SARS-CoV-2 Omicron Spike sequences, 1,276 positions (FULL LENGTH)
**Compute:** NVIDIA A100 80GB PCIe, torch 2.12.1+cu130 (CUDA), numba 0.66.0

> **CORRECTION NOTICE:** Two confirmed defects in the original pipeline were
> fixed on Aug 7, 2026 (A1: gap-stripping column misalignment; A2: 8-bit QM
> wrap-around). All numbers below are the **verified corrected** values. The
> pre-correction values (1,249 variable positions, MI 1.5917, 152 rules, LOO
> 2.93%) were artifacts of those defects and are documented in
> [My_Interpretation_of_the_Results/CORRECTION_NOTICE.md](My_Interpretation_of_the_Results/CORRECTION_NOTICE.md).
> This file supersedes any earlier draft that carried the stale numbers.
> **Corrections A3/A4 (Nov 2026)** — minimiser bookkeeping, constraint-function
> predictor, MIp check, mfDCA sign, 3D statistics — are applied below; see the
> CORRECTION_NOTICE.

---

## 1. GPU Acceleration (coevolution_gpu.py)

All heavy computations run on the A100 via torch CUDA:

| Kernel | Full-length scale | GPU time |
|--------|------------------|----------|
| Full MI matrix (all 813K pairs) | 1,299 seqs × 1,275 positions | **1.6 s** (497K pairs/s) |
| Entropy (all positions) | 1,275 positions | <0.1 s |
| Majority refs | 1,275 positions | <0.1 s |
| Coupling J = ln(P/P_exp) | per pair 20×20 | <0.1 s |
| H1 Gray adjacency | 1.65M consecutive pairs | <0.1 s |

CPU baseline for the same MI matrix: **20+ minutes** (timed out). GPU: **1.6 s** — ~800× speedup.

---

## 2. All 23 Scripts — FULL-LENGTH Corrected Results (no truncation)

| # | Script | Status | Key Result (corrected) |
|---|--------|--------|------------------------|
| 1 | `coevolution_shared.py` | ✅ | shared module: FASTA, MI, entropy, coupling, GPU pair-finder |
| 2 | `coevolution_gpu.py` | ✅ GPU | CUDA kernels: MI matrix, entropy, refs, H1, coupling |
| 3 | `master_boolean.py` | ✅ GPU | **36 prime implicants, 35 essential**, 10 pairs (A3) |
| 4 | `boolean_co-evolution.py` | ✅ GPU | whole-protein dipeptide Boolean minimization |
| 5 | `nary_kmap_co-evolution.py` | ✅ | base-20 (n-ary) K-map motifs (no position info) |
| 6 | `position_kmap_coevolution.py` | ✅ GPU | per-position-pair K-maps with MI |
| 7 | `run_allseq_analysis.py` | ✅ GPU | all-sequence MI matrix (window 30) |
| 8 | `run_kmap_analysis.py` | ✅ | H1 Gray adjacency 0.2167 (1.34× vs the uniform 5/31 null), Walsh–Hadamard (fixed transform: mode 0 + 4 modes carry 47% of the energy) |
| 9 | `flipped_boolean_coevolution.py` | ✅ GPU | **246 never-observed cubes** over the 10 pairs (52 essential) |
| 10 | `kmap_boolean_coevolution.py` | ✅ GPU | 35 essential rules across 10 pairs (literals = He-2012 code bits) |
| 11 | `variable_position_coevolution.py` | ✅ GPU | variable-positions-only K-maps with don't-cares |
| 12 | `predictive_constraint_function.py` | ✅ | train/test (file order) **25.9%**; conditional baseline 15.3%; 1,441 of 2,779 test mutations have a residue never seen in training |
| 13 | `allseq_constraint_function.py` | ✅ CPU | LOO-CV **99.4%** (3241/3259); conditional baseline 65.5%; 98.3% on the 300 unique sequences |
| 14 | `dca_boolean_coevolution.py` | ✅ | local-precision → Boolean, avg accuracy 0.0 (NOT real DCA) |
| 15 | `dca_mf_analysis.py` | ⚠ | mfDCA; on this alignment Meff = 1 (all sequences within 80% identity), so its pair ranking is not meaningful; DI sign fixed (A4) |
| 16 | `perplexity_coevolution.py` | ✅ | no position has perplexity > 3 on the corrected alignment, so no pairs |
| 17 | `advanced_co-evolution_analysis.py` | ✅ | network **21 nodes / 8 edges**, 24 variant signatures |
| 18 | `full_length_analysis.py` | ✅ GPU | 21 var, **5 high-MI pairs** (MI > 0.5) |
| 19 | `gpu_full_analysis.py` | ✅ GPU | max MI = **0.8067** at (373,378), full matrix saved |
| 20 | `create_mi_heatmap.py` | ✅ GPU | 813,450 pairs, max MI = 0.8067 |
| 21 | `generate_co-evolution_md.py` | ✅ | COEVOLUTION_KMAP_BOOLEAN.md (36 rules) |
| 22 | `generate_full_analysis_md.py` | ✅ | FULL_COEVOLUTION_ANALYSIS.md |
| 23 | `generate_full_pipeline_doc.py` | ✅ | FULL_PIPELINE_ANALYSIS.md |

---

## 3. Cross-Script Consistency (verified, corrected)

| Metric | Value | Scripts agreeing |
|--------|-------|------------------|
| Max MI | **0.8067** at (373,378) | gpu_full ✓, full_length ✓ |
| Variable positions | **21** / 1,276 | all scripts ✓ |
| Co-evolving pairs (mutation-only MI > 0.1) | **10** | master_boolean ✓, flipped ✓ |
| Master rules | **36 prime implicants, 35 essential** | master_boolean ✓, kmap_boolean ✓, exact_master_boolean ✓ |
| Never-observed cubes | **246** | flipped ✓ |
| H1 (5-bit group-order Gray) | 0.2167 (1.34× vs 5/31) | run_kmap |
| H1 (Gray of the He-2012 index, gaps excluded) | 0.1908 (1.18× vs 5/31) | gpu_full |
| Network | 21 nodes / 8 edges | advanced ✓ |

**NOTE on H1:** run_kmap and gpu_full use different residue orders before the
Gray map (group order vs He-2012 order), so they give different H1 values. The
5/31 null assumes uniform, distinct 5-bit codes; the composition-matched null on
PSICOV150 (3D_Structure_Analysis/11 §5) puts the effect at 1.18×.

---

## 4. LOO-CV (corrected, A4)

LOO-CV accuracy = **99.4% (3241/3259)**: for each sequence carrying a mutation at
one of the 10 pairs, the constraint function built from the other 1,298 sequences
predicts its residue at j from its residue at i. The earlier 9.24% came from the
eps-smoothed argmax, which picked a never-observed residue for 90.8% of the
held-out mutations. The LOO is exact and O(N) per pair. On the 300 unique
sequences it is 98.3% (conditional-frequency baseline 67.6%): the partner residues
are linked alleles of one clonal population.

---

## 5. Biological Interpretation (Full-Length, Corrected)

1. **Strongest co-evolution:** (373,378) MI = 0.8067 — but this pair is NOT
   3D-close (13–27 Å); its covariation is immune/lineage-driven.
2. **Pairs in 3D contact:** re-tested against residue pairs at the same sequence
   separation, with the 11 templates as units (`analysis/spike_3d_matched_check.py`):
   none of the 10 co-evolving pairs is enriched (the five 442–498 RBD pairs are in
   contact in 0/11 templates); only the i,i+3/4 pairs (407,410), (500,503), (503,507)
   are, which reflects the fold, not co-evolution.
3. **Never-observed combinations:** 246 cubes; absence of evidence in 299 unique
   sequences, not demonstrated incompatibility.
4. **Network:** 21 variable positions, 8 edges — small, interpretable.
5. **5 high-MI pairs** (MI > 0.5) full-length.
6. **Housekeeping confirmed:** conserved positions buried (52%), variable
   exposed (26%); Spearman(entropy, burial) = −0.17, p = 6e-9.
7. **Perplexity ratio** correlates with 3D contact (ρ = 0.67, p = 0.006, n = 15
   pairs; same template caveat).
8. **LOO-CV 99.4%:** the residue at one position of a pair nearly determines the
   other (linkage). Withholding relatives within 12 mismatches: 67.3% vs 65.0%
   frequency baseline (`analysis/spike_lineage_loo.py`).

---

## 6. The Two Corrected Defects (summary)

| Defect | Effect | Fix |
|--------|--------|-----|
| **A1** gap-stripping | "column j" mixed different raw positions → 1,249 fake variable positions, MI 1.59 artifacts | keep alignment, gap = state 20, exclude from counts (`aligned=True`) |
| **A2** 8-bit QM wrap-around | 20×20 (400-cell) map → 8 bits → cells 256–399 wrapped onto 0–143 → 143/152 phantom rules | pad to 32×32 (5 bits/axis, 10 bits), rows/cols 20–31 don't-care; `kmap_truth_table` raises on non-power-of-4 |

**Corrected headline numbers:** 21 variable positions · 10 co-evolving pairs ·
36 prime implicants (35 essential) · 246 never-observed cubes · max MI
0.8067 · LOO-CV 99.4% · 0 phantom rules.

---

## 7. Reproduce

```bash
G=/store/shuvam/.venv/bin/python   # torch 2.12.1+cu130, numba 0.66.0
export OMP_PROC_BIND=FALSE         # REQUIRED: torch otherwise pins to CPU 0
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
D=co-evolution-analysis   # this repository
cd $D

# Reproduce the corrected rule set (36 primes, 35 essential)
$G master_boolean.py

# Full-length GPU analyses
nohup $G -u $D/gpu_full_analysis.py > logs/gpu_full_analysis.log 2>&1 &
nohup $G -u $D/full_length_analysis.py > logs/full_length_gpu.log 2>&1 &
nohup $G -u $D/create_mi_heatmap.py > logs/mi_heatmap_gpu.log 2>&1 &
# monitor: nvidia-smi, tail -f logs/*.log
```

---

## 8. Multi-Dataset Validation (the campaign)

The 23-script suite was run on 6 external benchmark datasets (3,735 audited
runs; 3,729 ok, 6 documented-infeasible GPU-OOM giants on Pfam 222–486k seqs).
Coverage audit: 0 truncation flags. GPU/CPU consistency: 5/5 PASS. See
`My_Own_Interpretation_Across_New_Datasets/` and `plans/` for the full
campaign, and `agents.md` for the reproducibility commands.

---

*All numbers above are taken from the result JSONs in this repository
(master_boolean_summary, flipped_boolean_summary, full_gpu_results, etc.).*
