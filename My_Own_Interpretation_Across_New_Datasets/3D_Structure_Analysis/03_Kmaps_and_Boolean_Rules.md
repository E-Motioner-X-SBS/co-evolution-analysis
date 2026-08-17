# 03 — K-maps, Boolean Rules, and the 23-Script Pipeline

This document explains how the sequence-deduced rules were built: every K-map,
every pair statistic, the Boolean minimization, the complexity of each step,
and how the computations were optimized for GPU (A100, torch CUDA) and CPU.
All of this ran for days across multi-core/GPU machines; 3,687+ verified
script runs were produced.

---

## 1. The 23 scripts of the co-evolution pipeline

| # | Script | What it computes |
|---|---|---|
| 1 | `coevolution_shared.py` | shared module: FASTA parsing, entropy, MI, coupling, GPU pair-finder (single source of truth) |
| 2 | `coevolution_gpu.py` | GPU kernels (torch CUDA): MI matrix, entropy, majority refs, H1, coupling |
| 3 | `master_boolean.py` | master Boolean function: per-pair mutation K-maps → QM → **36 prime implicants, 2 essential** |
| 4 | `boolean_co-evolution.py` | whole-protein dipeptide Boolean minimization + coupling constants |
| 5 | `nary_kmap_co-evolution.py` | base-20 (n-ary) K-map analysis |
| 6 | `position_kmap_coevolution.py` | per-position-pair K-maps with MI |
| 7 | `run_allseq_analysis.py` | all-sequence MI matrix (window 30) |
| 8 | `run_kmap_analysis.py` | binary Gray K-map pipeline: H1 adjacency enrichment, signatures, Walsh–Hadamard, mutations |
| 9 | `flipped_boolean_coevolution.py` | flipped (forbidden) K-maps → **490 forbidden rules (117 unique)** |
| 10 | `kmap_boolean_coevolution.py` | K-map Boolean equations (5-bit Gray literals) |
| 11 | `variable_position_coevolution.py` | variable-positions-only K-maps with don't-cares |
| 12 | `predictive_constraint_function.py` | constraint function train/test prediction |
| 13 | `allseq_constraint_function.py` | leave-one-out CV with the constraint function |
| 14 | `dca_boolean_coevolution.py` | local precision matrix → Boolean (NOT real DCA) |
| 15 | `dca_mf_analysis.py` | proper mfDCA (Frobenius norm, APC, direct information) |
| 16 | `perplexity_coevolution.py` | perplexity ratio per pair |
| 17 | `advanced_co-evolution_analysis.py` | network, clustering, Walsh–Hadamard, variant signatures |
| 18 | `full_length_analysis.py` | full-length entropy + MI (window 30) |
| 19 | `gpu_full_analysis.py` | GPU full analysis: full MI matrix, entropy, refs, H1 |
| 20 | `create_mi_heatmap.py` | MI heatmap visualizations |
| 21 | `generate_co-evolution_md.py` | markdown of the Boolean functions |
| 22 | `generate_full_analysis_md.py` | comprehensive report |
| 23 | `generate_full_pipeline_doc.py` | full pipeline documentation |

## 2. The K-map construction (per position pair)

**Encoding.** Each amino acid is mapped to a 5-bit Gray code
(`kmap_sbm.encoding.gray_amino`, proven in Lean 4 — see `lean_proofs/`):
the 20 canonical residues occupy 20 of the 32 codes of the Q₅ hypercube such
that single bit-flips connect physicochemical neighbours.

**Frequency K-map.** For a position pair (i, j), from the aligned sequences
(gaps = state 20, excluded from counts):

```
K[a, b] = (1/N) * Σ_s 1[seq_s[i] = a AND seq_s[j] = b],   a, b in 1..20
```

**Mutation K-map** (used for Boolean minimization): threshold the frequency
map to a ternary function on the 20×20 grid:

```
f(a, b) = 1   if (a, b) observed with count > 0 and (a, b) != reference pair
        = -1  (don't care) if (a, b) is the majority/reference pair
        = 0   if (a, b) never observed
```

**32×32 padded encoding (CORRECTED — A2).** The original pipeline fed the
400-cell 20×20 map to a Quine–McCluskey implementation that derived
`k_bits = log2(400)//2 = 4` → 8 bits → cells 256–399 silently wrapped onto
cells 0–143, producing 143 phantom rules. The fix: pad to 32×32 (5 bits per
axis, 10 bits total = 1,024 cells), mark rows/cols 20–31 as don't-care, and
decode 5 bits per axis. `kmap_truth_table` now raises on non-power-of-4
inputs.

**Quine–McCluskey minimization.** Input: the 1,024-cell truth vector
(0/1/don't-care). Steps:
1. enumerate minterms (on-set cells) grouped by popcount;
2. iteratively merge adjacent minterms differing in one bit
   (the merge is GPU/vector-optimized: popcount grouping reduces the pairwise
   comparisons from O(m²) to O(Σ groups) — ~70× faster than the naive
   implementation; verified against a brute-force reference on 50 random
   trials);
3. collect **prime implicants**; select the **essential prime implicants**
   (covering on-set cells covered by no other PI) — the minimal irredundant
   core.

The result for the Omicron Spike: **36 distinct prime implicants** (residue
rules `pos_i = aa_i AND pos_j = aa_j`) with a **2-rule essential core**:
`(212=G, 216=R)` and `(210=K, 215=G)` — after the A2 correction (the
pre-correction "152 rules" were 143 phantom).

**n-ary K-map.** The base-20 variant: k-mer frequency maps
`cell = Σ_{window} ∏ 20^{position} code(aa)` for k = 2, normalized to
frequencies, minimized directly in base-20 (20×20 = 400 cells, no binary
projection). Produces protein-wide dipeptide motifs (row_aa, col_aa) — these
carry **no position information** and therefore cannot be mapped to 3D
positions (honest limitation documented).

## 3. The pair statistics

**Shannon entropy** per position (gaps excluded):

```
H(p) = − Σ_a P(a|p) log2 P(a|p)
```

**Mutual information** per position pair (gaps excluded, state-20 handled):

```
MI(i,j) = Σ_{a,b} P(a,b) log2( P(a,b) / (P(a)P(b)) )
```

**Mutation-only MI** (the master-Boolean convention): the reference pair
(majority residues) is excluded from the joint counts — isolates the
mutation-driven covariation.

**Perplexity ratio** (determinism):

```
r(i,j) = PP(j) / mean_a PP(j | i = a),     PP = 2^H
```

**Constraint function / coupling:**

```
C(a,b) = ln( P(a,b) / (P(a)P(b)) )
P_co-evolution(a,b) = σ(C) = 1 / (1 + e^{−C})
```

**Gray adjacency (H1).** For consecutive residues, the Hamming distance in
5-bit Gray space is counted; the observed fraction at distance 1 is compared
to the random expectation 5/31 ≈ 0.161.

## 4. Complexity and optimization

| Step | Complexity | Optimization |
|---|---|---|
| Entropy (all positions) | O(N·L) | vectorized numpy bincount; GPU < 0.1 s |
| MI for one pair | O(N) | numpy bincount over 400 cells |
| Full MI matrix (window 30, L positions) | O(L·N) | **GPU**: dense [N×L] int32 tensor; per-chunk `scatter_add` bincount; **adaptive chunking** bounds memory (flat int64 buffer ≤ 512 MB → chunk = max(256, 512MB/(N·8))) so 28k–486k-sequence families never OOM |
| QM minimization | O(m · 2^bits) worst case (m = on-set cells) | popcount-grouped merge; 32×32 padded maps; verified vs brute force |
| LOO-CV constraint function | O(10 · N²) | **inherently O(N²)** — for N = 28k ≈ 15–20 h; for N = 222–486k ≈ 27 days – 2.5 years of single-core compute → **documented as infeasible at full scale** (runs per no-timeout directive; or documented sampling is the only practical route) |
| H3 pairwise K-map distances | O(n² · d), d = 1,024 | GPU norm-trick `D = a² + b² − 2abᵀ` with **chunked statistics** (never materializes the n×n matrix; the CPU fallback is also chunked — the original scipy `pdist` fallback attempted an 883 GB allocation and was replaced) |
| Pairwise clustering (advanced) | O(n²) | capped at 2,000 sequences with explicit logging |

**Key GPU/CPU engineering facts** (all verified this session):
- full MI matrix 813K pairs × 1,299 seqs in 1.6 s (GPU) vs >20 min CPU (~800×).
- H3 on 152,041 sequences (11.6 billion pairs): GPU chunked, exit 0.
- GPU/CPU numerical consistency: 5/5 datasets PASS (max |Δ| < 1e-7;
  reference-residue tie-breaking unified to lowest-code between CPU and GPU).
- **torch CPU-pinning bug found and fixed**: `import torch` pins the process to
  CPU 0 unless `OMP_PROC_BIND=FALSE` (was using 1 of 24 cores; now all cores).
- memory-fraction cap (0.1) + `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`
  prevent allocator deadlocks on the shared GPU.

## 5. Reproducibility of the rule sets

The corrected rule sets (21 variable positions, 10 pairs, 36 PIs, 2
essential, 490 flipped) are reproduced by running the parameterized scripts
with defaults (no `COEVO_FASTA`/`COEVO_RESULTS` env) — regression-verified:
`master_boolean.py` re-run gives byte-identical summary values (36/2).

All run artifacts carry `run_log.json` (exit status, duration, stdout tail,
output file list) — every one of the 3,687 runs is auditable.
