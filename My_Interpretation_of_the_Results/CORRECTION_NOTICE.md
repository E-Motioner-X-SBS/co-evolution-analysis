# CORRECTION NOTICE: Two Confirmed Defects in the Original Pipeline

**Date:** August 7, 2026
**Status:** Verified by direct computation and independent audit

This document records two systemic defects discovered in the original co-evolution pipeline, their verification, and the corrected results. It supersedes the biological claims in the earlier documentation where they conflict with the numbers below.

---

## Defect A1: Gap-stripping caused column misalignment

### The bug

The original pipeline built position arrays with:

```python
clean = "".join(aa for aa in seq if aa in encoder.encode)  # DELETES gaps
arr = np.array([encoder.encode.get(aa, -1) for aa in clean[:max_pos]])
```

This **deletes** gap characters (`-`, 8,038 total) and ambiguous characters (`X`, 356) instead of keeping them. Because different sequences have gaps at different positions, "column j" of the cleaned array corresponds to **different raw alignment positions** in different sequences.

### Verification

- Raw-aligned columns 372, 401, 413, 427, 852 are **conserved** (H = 0.00 to 0.05).
- Only **21 of 1,276** raw-aligned positions have H > 0.3 (gaps excluded from the count); 32 if gaps are counted as a state.
- The gap count before raw column 372 ranges from 1 to 9 across sequences, so the cleaned "column 372" mixes raw positions 373 to 381.
- Therefore the reported "1,249 variable positions", H(372) = 1.63, MI(372, 401) = 1.5917, the 36,918 co-evolving pairs, and the (413, 427) "271 mutations" are **artifacts of misalignment** - they are reproducible from the code, but they describe misaligned columns, not protein positions.

### The fix

Keep the full alignment; encode `-` and unknown characters as state 20 (21 states total):

```python
arr = np.array([encoder.encode.get(aa, 20) for aa in seq[:actual_max]])
```

All gap-aware functions (entropy, MI, majority refs, coupling) now exclude state 20 from counts. `coevolution_shared.load_position_arrays(..., aligned=True)` implements this; `aligned=False` reproduces the buggy behavior.

## Defect A2: 8-bit QM encoding corrupted rule labels

### The bug

The original fed a 400-cell (20 x 20) K-map to `boolean_minimize_kmap`, which derives `k_bits = int(log2(400))//2 = 4`, i.e., 8 bits = 256 minterms. Cells 256-399 **wrap onto cells 0-143**, creating phantom on-set minterms. The decoder then read 4+4 bits (a 16 x 16 grid) instead of the 20 x 20 grid, scrambling every rule label.

### Verification

- Direct test: cell 256 (row 13 = K, col 0 = A) wraps to minterm 0.
- Of the 152 "essential rules" in the original `boolean_functions.json`, only 9 have residue labels that were ever observed in the data. 143 are phantom (e.g., "(413=W, 427=E)" is never observed; the true on-set at (413, 427) is A-V, I-C, Y-A, D-I, Q-D, N-I, N-W, K-S, K-G, T-F, P-P, G-T).
- 4 of 12 on-set cells at (413, 427) involve residues T, C, P, G (He indices 16-19), which the 4-bit encoding cannot represent.

### The fix

Pad the 20 x 20 map to 32 x 32 (5 bits per axis = 10 bits), mark rows/cols 20-31 as don't-care, run QM on 1,024 cells (no wrap-around), and decode 5 bits per axis. `kmap_truth_table` now raises an error if given a non-power-of-4 cell count, preventing silent corruption.

## Corrected Results (verified, all scripts re-run with aligned columns + padded QM)

From `analysis/corrected_pipeline.py` and the re-run of every script (Aug 7, 2026):

| Metric | Original (buggy) | Corrected |
|--------|------------------|-----------|
| Variable positions (H > 0.3) | 1,249 | **21** |
| Co-evolving pairs (MI > 0.1, window 30) | 36,918 | **10-12** |
| Strongest full-MI pair | (372, 401) MI = 1.5917 | **(373, 378) MI = 0.8067** |
| Strongest mutation-only MI pair | (413, 427) MI = 2.352 | **(495, 498) MI = 0.8710** |
| Master Boolean prime implicants | 162 (143 phantom) | **36 distinct pairs** (all cover real cells) |
| Essential rules | 152 | **2** (minimal irredundant core) |
| Forbidden rules (flipped) | 345 | **490** |
| Mean mutations per sequence | 1,061 | **11.1** |
| High-MI pairs (MI > 0.5) | 35,858 | **5** |
| MI > 1.0 pairs | 106,626 | **0** |
| LOO-CV accuracy | 2.93% | **9.24%** (301/3259) |
| Train/test accuracy | 5.84% | **0.11%** |
| Network nodes/edges | 1,249 / 35,098 | **21 / 8** |

### Code corrections applied to ALL scripts (Aug 7, 2026)

Every analysis script was corrected and re-run:

1. **Aligned columns (FIX A1):** all scripts now encode the full alignment with
   gap = state 20 (was gap-stripped, misaligned). Changed in: coevolution_shared,
   coevolution_gpu, master_boolean, gpu_full_analysis, full_length_analysis,
   create_mi_heatmap, advanced_co-evolution_analysis, perplexity_coevolution,
   run_allseq_analysis, position_kmap_coevolution, variable_position_coevolution,
   flipped_boolean_coevolution, kmap_boolean_coevolution, dca_boolean_coevolution,
   allseq_constraint_function, predictive_constraint_function,
   generate_co-evolution_md. The binary-path scripts (run_kmap_analysis,
   boolean_co-evolution, nary_kmap_co-evolution) use the group-order binary
   encoding and are documented as historical.

2. **Padded QM (FIX A2):** every Quine-McCluskey input is now 32x32 padded
   (10 bits, 5 per axis) instead of raw 20x20 (8 bits, wrap-around). The QM
   merge was optimized (popcount grouping, 70x faster) and verified against
   a brute-force reference (50 random trials). kmap_truth_table now raises
   on non-power-of-4 inputs instead of silently corrupting.

3. **All result JSONs regenerated** with the corrected pipeline (see table above).

### Corrected top co-evolving pairs

| Pair | MI | Reference |
|------|-----|-----------|
| (373, 378) | 0.8067 | (F, A) |
| (18, 26) | 0.8024 | (I, S) |
| (378, 407) | 0.7920 | (A, N) |
| (66, 94) | 0.7907 | (A, T) |
| (215, 216) | 0.7571 | (G, R) |
| (210, 215) | 0.7532 | (N, G) |
| (407, 410) | 0.7495 | (N, S) |
| (210, 216) | 0.7453 | (N, R) |
| (212, 215) | 0.3977 | (V, G) |
| (212, 216) | 0.3773 | (V, R) |
| (488, 495) | 0.3419 | (F, R) |
| (210, 212) | 0.1769 | (N, V) |

### Corrected rules

The corrected master Boolean function has **36 distinct prime implicants** (residue-pair rules) across 10 pairs, of which **2 are essential** (the minimal irredundant core). Every rule was verified programmatically against the alignment: each listed residue pair is actually observed, so there are no phantom rules. The 2 essential rules:

```
IF pos 212 = G AND pos 216 = R THEN co-evolutionary (mutation-only MI = 0.377)
IF pos 210 = K AND pos 215 = G THEN co-evolutionary (mutation-only MI = 0.238)
```

Note on MI conventions: the top-pair table in this notice uses FULL MI (all pairs incl. reference); the master Boolean pipeline ranks by MUTATION-ONLY MI (reference pair excluded). Both are valid; they rank pairs differently - e.g. (495, 498) has full MI 0.0386 but mutation-only MI 0.8710, because its signal lives in the non-reference pairs.

### What this means

1. The Spike Omicron dataset (1,299 sequences) is **much more conserved** than the original analysis suggested: only 21 positions vary meaningfully, and only 10-12 position pairs genuinely co-evolve.

2. The strongest genuine co-evolution signals are (373, 378), (18, 26), (66, 94), and the 210-216 cluster - NOT the (372, 401) hub reported earlier.

3. The original 152 rules were mostly artifacts. The corrected rule set is 36 distinct real prime implicants with a 2-rule essential core.

4. The DCA analysis (dca_mf_analysis.py) was NOT affected by A1 (it uses the raw alignment with gap as state 20); its direct-coupling results stand, but the DCA-vs-MI comparison must be re-interpreted against the corrected MI matrix (max 0.807).

5. The LOO-CV accuracy IMPROVED to 9.24% on corrected data: the genuine pairs are more predictable than the artifact pairs.

## Reproducibility

```bash
cd /store/shuvam/E-motioner-X-SBS/datasets/co-evolution
python analysis/corrected_pipeline.py
# -> analysis/corrected_results.json
```

The corrected loader is `coevolution_shared.load_position_arrays(aligned=True)` (default). All scripts that consume position arrays should use it. The `My_Interpretation_of_the_Results` documentation has been updated to state which numbers are corrected and which are historical (buggy).
