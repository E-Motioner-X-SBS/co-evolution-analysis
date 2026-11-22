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
| Essential rules | 152 | **35** (correction A3; first reported as 2) |
| Forbidden rules (flipped) | 345 | **490** |
| Mean mutations per sequence | 1,061 | **11.1** |
| High-MI pairs (MI > 0.5) | 35,858 | **5** |
| MI > 1.0 pairs | 106,626 | **0** |
| LOO-CV accuracy | 2.93% | **99.4%** (3241/3259; first reported as 9.24%, correction A4) |
| Train/test accuracy | 5.84% | **25.9%** (first reported as 0.11%, correction A4) |
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

The corrected master Boolean function has **36 distinct prime implicants** (residue-pair rules) across 10 pairs, of which **35 are essential** (see correction A3 below; the "2 essential" first reported here came from a defect in the minimiser). Every rule was verified programmatically against the alignment: each listed residue pair is actually observed, so there are no phantom rules.

Note on MI conventions: the top-pair table in this notice uses FULL MI (all pairs incl. reference); the master Boolean pipeline ranks by MUTATION-ONLY MI (reference pair excluded). Both are valid; they rank pairs differently - e.g. (495, 498) has full MI 0.0386 but mutation-only MI 0.8710, because its signal lives in the non-reference pairs.

### What this means

1. The Spike Omicron dataset (1,299 sequences) is **much more conserved** than the original analysis suggested: only 21 positions vary meaningfully, and only 10-12 position pairs genuinely co-evolve.

2. The strongest genuine co-evolution signals are (373, 378), (18, 26), (66, 94), and the 210-216 cluster - NOT the (372, 401) hub reported earlier.

3. The original 152 rules were mostly artifacts. The corrected rule set is 36 distinct real prime implicants, 35 of them essential (correction A3).

4. The DCA analysis (dca_mf_analysis.py) was NOT affected by A1 (it uses the raw alignment with gap as state 20); its direct-coupling results stand, but the DCA-vs-MI comparison must be re-interpreted against the corrected MI matrix (max 0.807).

5. The LOO-CV accuracy of 9.24% reported here was itself an artifact of the predictor (correction A4); the corrected value is 99.4%.

## Reproducibility

```bash
cd co-evolution-analysis
python analysis/corrected_pipeline.py
# -> analysis/corrected_results.json
```

The corrected loader is `coevolution_shared.load_position_arrays(aligned=True)` (default). All scripts that consume position arrays should use it. The `My_Interpretation_of_the_Results` documentation has been updated to state which numbers are corrected and which are historical (buggy).

---

## Correction A3 (Nov 2026): the essential count — 35, not 2

The published "36 distinct prime implicants, 2 essential" came from a defect
in `kmap_sbm.analysis.prime_implicants.boolean_minimize_kmap`
(kmap-sbm-validation), not from the data:

* **Index mismatch.** Each implicant's coverage indexed one minterm order
  (on-set first, then don't-cares) while the essential test and the greedy
  cover used another (cells in table order, don't-cares interleaved).  On the
  padded Spike maps every on-cell after the first row's padding was shifted,
  so the reported essentials were arbitrary and the reported **cover reached
  only 2 of the 36 observed mutations** (8 of the 10 pairs had an empty
  cover).  `essential_prime_implicants` also included the greedy picks, so
  `n_essential` always equalled `covering_size`.
* **Not exact.** Only don't-cares adjacent to an on-cell took part in
  merging, so cubes could not grow through the padding and some reported
  cubes were not prime.

Both are fixed in kmap-sbm-validation (exact merging, coverage indexed by the
on-set, essentials and cover reported separately).  `analysis/exact_master_boolean.py`
rebuilds the same 10 padded maps and checks the minimiser against an
independent exact QM and brute-force enumeration (`analysis/exact_master_boolean.json`):

| | prime implicants | essential | cover |
|---|---|---|---|
| published (defective minimiser) | 36 distinct decoded rules | **2** | reached 2 of 36 on-cells |
| fixed minimiser / exact QM, padding = don't-care (as master_boolean.py) | **36** | **35** | 35, complete |
| exact QM, padding = off-set | 35 | **35** | 35, complete |

35 of the 36 observed co-mutations are isolated K-map cells, so almost every
rule is essential; the one merge is at pair (210, 212).  On these maps QM
therefore compresses 36 observations into 35 rules: the rules are, in
effect, the list of observed mutation pairs.  `master_boolean_summary.json`,
`kmap_boolean_coevolution/`, `flipped_boolean_results/` and the other
minimiser outputs have been regenerated with the fixed minimiser.  Read
"2 essential" anywhere in this repository as "35 essential".

## Correction A4 (Nov 2026): other defects found in a full read

| Defect | Effect | Fix / corrected value |
|---|---|---|
| Constraint-function predictor took argmax of the eps-smoothed C over all 20 residues | a never-observed partner scores ln(1/P(a)) > 0 and won 90.8% of predictions → "LOO-CV 9.24%", "train/test 0.11%" | argmax over partners that co-occur with the residue (ties by count); exact O(N) LOO: **99.4%** (3241/3259), conditional baseline 65.5%, **98.3%** on the 300 unique sequences; train/test 25.9%. Partners are linked lineage alleles: withholding every sequence within 12 mismatches of the test sequence drops accuracy to 67.3%, level with the 65.0% baseline; at 20 mismatches 45.8% vs 61.4% on the 54% of mutations still predictable (`analysis/spike_lineage_loo.py`). |
| `path_m` M4 "MIp" never accumulated the column sum | MIp = MI, so "APC adds nothing" | true MIp on the PSICOV holdout **0.176** vs circuit 0.179 (p = 0.55); 5-fold CV 0.189 vs 0.165. MIp is now a standing baseline in `run_contact_campaign.py`. |
| `dca_mf_analysis.py` used W = exp(−e_ij) | direct information with flipped couplings | sign fixed; on PSICOV150 (`analysis/mfdca_psicov_check.py`) DI prec@L/5 **0.566** vs 0.517 with the old sign (better on 93, worse on 25 targets). This also shows `path_h`'s APC-Frobenius mfDCA (0.169), used as "rawDCA" in the rank-fusion result, is far below a working mfDCA. |
| `path_h_dca_circuit.mfdca_di` (the rank-fusion DCA input) | pair frequencies weighted by W³ (Gram matrix of one-hot·W·√W) against W for single sites; Frobenius norm outside the zero-sum gauge → prec@L/5 0.169, no better than MIp (0.189, p = 0.24) | `analysis/dca_fusion_recheck.py`: the published fusion reproduces exactly (0.204 vs 0.169; `path_l` re-run unchanged on regenerated caches), but with a correct mfDCA the K-map prior **lowers** prec@L/5 — DI 0.566 → 0.471 (−0.095, CI −0.116…−0.073), zero-sum APC-Frobenius 0.692 → 0.565. The "orthogonal information" claim is withdrawn; the "264 rescued / 185 lost" decomposition has no surviving code and is withdrawn |
| Circuit vs MIp, intervals | — | `analysis/circuit_vs_mip_ci.py`: 150 CV proteins, prec@L/5 circuit − MIp −0.024 (CI −0.055…+0.007); MIp has the higher AUC (0.676 vs 0.651, p = 1.5e-5) |
| Spike mfDCA | Meff = 1.0 (all sequences within 80% identity) | Spike DCA pairs, incl. "(454,495) DCA-only 5.97×", are not meaningful |
| `validate_contacts.py` mapped columns via the first alignment row and its random control counted out-of-range draws | PSICOV first rows are gapped homologues | best-matching-row mapping, in-range control: top-20 MI pairs **2.38×**, precision **0.071** (`analysis/psicov_top_pairs_check.py`; first reported 2.49× / 0.066) |
| Omicron 3D statistics | the 6 "contact" pairs are 3–5 residues apart; control pooled over separations; binomial test over 299 models built on 11 templates (8d55: 203 models) | re-tested with same-separation controls (all of chain A, and ±30 residues) and templates as units (`analysis/spike_3d_matched_check.py`): none of the 10 co-evolving pairs is in contact more than its separation predicts (the five 442–498 RBD pairs: 0/11 templates); only the i,i+3/4 pairs (407,410), (500,503), (503,507) are enriched (2.9–6.6× vs local, 8–11/11 templates), a property of the fold rather than of co-evolution |
| Other script bugs | `path_l` had no entry point; `path_j` used wrong separation bins, identity mapping and trained on its test protein 1atzA; `stage_chem` crashed; Walsh–Hadamard collapsed to one value; residue decoders read cubes by their lowest corner or as Gray/8-bit; gaps counted as G / as matching cube bits; `justify_top15` gap-stripped; `run_all_bg.sh` ran the superseded snapshot; report generators hard-coded false claims | all fixed; affected Spike outputs regenerated |
