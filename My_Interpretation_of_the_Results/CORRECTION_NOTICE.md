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
| Strongest MI pair | (372, 401) MI = 1.5917 | **(373, 378) MI = 0.8067** |
| Master Boolean prime implicants | 162 (143 phantom) | **255 (all cover real cells)** |
| Essential rules | 152 | **3** (minimal irredundant core) |
| Mean mutations per sequence | 1,061 | **11.1** |
| High-MI pairs (MI > 0.5) | 35,858 | **5** |
| MI > 1.0 pairs | 106,626 | **0** |
| LOO-CV accuracy | 2.93% | **9.24%** (301/3259) |
| Network nodes/edges | 1,249 / 35,098 | **21 / 8** |

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

The corrected master Boolean function has **255 total prime implicants** across 10 pairs, of which **3 are essential** (the minimal irredundant cover). Every one of the 255 PIs covers at least one real observed cell (verified programmatically); there are no phantom rules. Example essential rules:

```
IF pos 210 = N AND pos 212 = S THEN co-evolutionary
IF pos 212 = V AND pos 215 = P THEN co-evolutionary
IF pos 212 = S AND pos 216 = R THEN co-evolutionary
```

### What this means

1. The Spike Omicron dataset (1,299 sequences) is **much more conserved** than the original analysis suggested: only 21 positions vary meaningfully, and only 10-12 position pairs genuinely co-evolve.

2. The strongest genuine co-evolution signals are (373, 378), (18, 26), (66, 94), and the 210-216 cluster - NOT the (372, 401) hub reported earlier.

3. The original 152 rules were mostly artifacts. The corrected rule set is 255 real prime implicants with a 3-rule essential core.

4. The DCA analysis (dca_mf_analysis.py) was NOT affected by A1 (it uses the raw alignment with gap as state 20); its direct-coupling results stand, but the DCA-vs-MI comparison must be re-interpreted against the corrected MI matrix (max 0.807).

5. The LOO-CV accuracy IMPROVED to 9.24% on corrected data: the genuine pairs are more predictable than the artifact pairs.

## Reproducibility

```bash
cd /store/shuvam/E-motioner-X-SBS/datasets/co-evolution
python analysis/corrected_pipeline.py
# -> analysis/corrected_results.json
```

The corrected loader is `coevolution_shared.load_position_arrays(aligned=True)` (default). All scripts that consume position arrays should use it. The `My_Interpretation_of_the_Results` documentation has been updated to state which numbers are corrected and which are historical (buggy).
