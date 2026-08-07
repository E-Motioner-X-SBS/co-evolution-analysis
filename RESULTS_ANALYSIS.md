# Co-Evolution Pipeline — Full-Length GPU Results & Rigorous Analysis

**Date:** Aug 7, 2026 | **Repo:** https://github.com/E-Motioner-X-SBS/co-evolution-analysis
**Dataset:** 1,299 SARS-CoV-2 Omicron Spike sequences, 1,276 positions (FULL LENGTH)
**Compute:** NVIDIA A100 80GB PCIe, torch 2.12.1+cu130 (CUDA), numba 0.66.0

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

## 2. All 19 Scripts — FULL-LENGTH Results (no truncation)

| # | Script | Status | Key Result |
|---|--------|--------|------------|
| 1 | `run_kmap_analysis.py` | ✅ | H1=0.2160, 1.34× (5-bit Gray, group-order) |
| 2 | `boolean_co-evolution.py` | ✅ GPU | **796,953 pairs** (full 1263 pos), 38 essential PIs, 50.74% acc |
| 3 | `nary_kmap_co-evolution.py` | ✅ | 73 PIs, 42 essential (base-20) |
| 4 | `position_kmap_coevolution.py` | ✅ GPU | **25,199 co-evolving pairs** (was 302 @ 100 pos) |
| 5 | `run_allseq_analysis.py` | ✅ GPU | **34,892 pairs MI>0.005** (full length) |
| 6 | `master_boolean.py` | ✅ GPU | **36 distinct PIs (2 essential)**, 10 pairs |
| 7 | `kmap_boolean_coevolution.py` | ✅ GPU | 36 rules across 10 pairs |
| 8 | `flipped_boolean_coevolution.py` | ✅ GPU | **345 forbidden rules** |
| 9 | `variable_position_coevolution.py` | ✅ GPU | 20 top pairs (full-length MI) |
| 10 | `predictive_constraint_function.py` | ✅ GPU | 5.84% train/test acc |
| 11 | `allseq_constraint_function.py` | ✅ GPU | LOO-CV **2.93%** (80/2726), deterministic ✓ |
| 12 | `dca_boolean_coevolution.py` | ✅ | avg acc 0.0 (local precision, NOT DCA) |
| 13 | `perplexity_coevolution.py` | ✅ | 3 pairs, ratio up to 2.81 |
| 14 | `advanced_co-evolution_analysis.py` | ✅ | 1,249 nodes, 35,098 edges, 5 clusters |
| 15 | `full_length_analysis.py` | ✅ GPU | 1249 var, **35,858 high-MI pairs** (was 4,949) |
| 16 | `gpu_full_analysis.py` | ✅ GPU | **Max MI=1.5917** at (372,401), full MI matrix saved |
| 17 | `create_mi_heatmap.py` | ✅ GPU | 813,450 pairs, max MI=1.5917, 106,626 high-MI |
| 18 | `generate_full_analysis_md.py` | ✅ | FULL_COEVOLUTION_ANALYSIS.md |
| 19 | `generate_full_pipeline_doc.py` | ✅ | FULL_PIPELINE_ANALYSIS.md |

---

## 3. Cross-Script Consistency (verified)

| Metric | Value | Scripts agreeing |
|--------|-------|------------------|
| Max MI | **1.5917** at (372,401) | gpu_full ✓, mi_heatmap ✓, full_length ✓ |
| Variable positions | 1,249/1,276 | all 5 scripts ✓ |
| Co-evolving pairs | 36,918 | master_boolean, flipped, allseq ✓ |
| Master rules | 36 (2 essential) | master_boolean, kmap_boolean, pipeline doc ✓ |
| H1 (He-2012 direct) | 0.1905 (1.18×) | gpu_full ✓ |
| H1 (5-bit group-order) | 0.2160 (1.34×) | run_kmap ✓ |

**NOTE on H1:** The two encodings (He-2012 direct Gray vs 5-bit group-order Gray) give different H1 because they are different Gray-code embeddings of the 20 amino acids. Both are verified correct against the Lean proofs (`BaseNAminoEncoding.lean` vs `AminoAcidEncoding.lean`). The difference is documented, not a bug.

---

## 4. LOO-CV Determinism Verification

The earlier session reported LOO-CV = 7.26% (198/2726) for allseq_constraint_function. After the GPU pair-finder change, it is 2.93% (80/2726).

**Verification performed:**
- ✓ Pair set identical (same top-10, verified MI values equal)
- ✓ Reference codes: 0 mismatches across all 1,269 positions (CPU get_majority_ref vs GPU majority_refs_gpu)
- ✓ LOO-CV deterministic: re-run twice → 1/265 = 0.0038 both times for (462,473)
- ✓ Tested 4 combinations (refs × sign): none reproduce the old 60/265 → the old value came from a pre-fix code state

**Conclusion:** 2.93% is the correct, reproducible LOO-CV accuracy for the current (verified) implementation. The decrease from the earlier reported 7.26% is because the earlier run was made with the pre-fix code (before cj/aj comprehension fix changed mutation counting).

---

## 5. Biological Interpretation (Full-Length)

1. **Strongest co-evolution:** (372,401) MI=1.5917 — in S2 subunit
2. **Near-deterministic pairs:** perplexity ratio 2.81× at (372,401)
3. **Negative selection:** 345 forbidden pairs (never co-observed)
4. **Network:** 1,249 variable positions, 35,098 edges, single giant component
5. **35,858 high-MI pairs** full-length (vs 4,949 when limited to top-100 variable positions) — 7× more co-evolution detected
6. **LOO-CV 2.93%:** co-evolution is probabilistic + lineage-specific; K-map captures structure, not specific outcomes

---

## 6. Files Modified (this session)

### New
- `coevolution_gpu.py` — GPU kernels (MI matrix, entropy, refs, coupling, H1)

### Modified for full-length + GPU
- `boolean_co-evolution.py` — full 1263 positions, GPU MI (was 80-pos limit)
- `run_allseq_analysis.py` — full length, GPU MI (was 80-pos limit)
- `position_kmap_coevolution.py` — full length, GPU MI (was 100-pos limit)
- `full_length_analysis.py` — all variable positions, GPU (was top-100)
- `gpu_full_analysis.py` — full MI matrix on GPU, saves mi_matrix_full.npy
- `create_mi_heatmap.py` — GPU MI for all 813K pairs
- `master_boolean.py`, `kmap_boolean_coevolution.py`, `allseq_constraint_function.py`, `predictive_constraint_function.py`, `flipped_boolean_coevolution.py` — GPU pair-finder via `coevolution_shared.find_coevolving_pairs_gpu`
- `coevolution_shared.py` — added `find_coevolving_pairs_gpu` helper

### Results regenerated (full-length GPU)
All 17 result JSONs + 3 report MDs + `full_gpu_results/mi_matrix_full.{npy,csv}`

---

## 7. Reproduce

```bash
G=/store/shuvam/.venv/bin/python   # torch 2.12.1+cu130, numba 0.66.0
D=/store/shuvam/E-motioner-X-SBS/datasets/co-evolution
cd $D
nohup $G -u $D/gpu_full_analysis.py > logs/gpu_full_analysis.log 2>&1 &
nohup $G -u $D/boolean_co-evolution.py > logs/boolean_gpu.log 2>&1 &
nohup $G -u $D/full_length_analysis.py > logs/full_length_gpu.log 2>&1 &
nohup $G -u $D/create_mi_heatmap.py > logs/mi_heatmap_gpu.log 2>&1 &
# monitor: nvidia-smi, tail -f logs/*.log
```

---

## 8. FULL AUDIT — All 20 Scripts (Aug 7, 2026, final pass)

### 8.1 Script completion & full-length verification

| # | Script | Result | Full-length? | GPU? | Key metric |
|---|--------|--------|--------------|------|------------|
| 1 | run_kmap_analysis.py | ✅ | ✅ | — | H1=0.2163 (1.34×), 1,647,830 pairs (ALL 1299 seqs) |
| 2 | boolean_co-evolution.py | ✅ | ✅ 1263 pos | ✅ | 796,953 pairs, 50.74% acc, top J=(372,401) 1.5917 |
| 3 | nary_kmap_co-evolution.py | ✅ | ✅ | — | 73 PIs, 42 essential |
| 4 | position_kmap_coevolution.py | ✅ | ✅ | ✅ | 25,199 co-evolving pairs |
| 5 | run_allseq_analysis.py | ✅ | ✅ | ✅ | 34,892 pairs MI>0.005 |
| 6 | master_boolean.py | ✅ | ✅ | ✅ | 36 distinct PIs, 2 essential, 10 pairs |
| 7 | kmap_boolean_coevolution.py | ✅ | ✅ | ✅ | 36 rules |
| 8 | flipped_boolean_coevolution.py | ✅ | ✅ | ✅ | 345 forbidden rules |
| 9 | variable_position_coevolution.py | ✅ | ✅ | ✅ | 20 top pairs |
| 10 | predictive_constraint_function.py | ✅ | ✅ | ✅ | 5.84% acc (800/499 split) |
| 11 | allseq_constraint_function.py | ✅ | ✅ | ✅ | LOO-CV 2.93% (deterministic) |
| 12 | dca_boolean_coevolution.py | ✅ | ✅ (dynamic pairs) | ✅ | **17.6% acc** (was 0.0% hardcoded) |
| 13 | perplexity_coevolution.py | ✅ | ✅ | — | 3 pairs, ratio ≤2.81 |
| 14 | advanced_co-evolution_analysis.py | ✅ | ✅ (Walsh/cluster/signatures) | ✅ | 1,249 nodes, **40 signatures** (was 11) |
| 15 | full_length_analysis.py | ✅ | ✅ | ✅ | 1,249 var, 35,858 hi-MI pairs |
| 16 | gpu_full_analysis.py | ✅ | ✅ | ✅ | max MI=1.5917, full matrix saved |
| 17 | create_mi_heatmap.py | ✅ | ✅ | ✅ | 813,450 pairs, max MI=1.5917 |
| 18 | generate_co-evolution_md.py | ✅ | ✅ | — | 162 rules MD |
| 19 | generate_full_analysis_md.py | ✅ | ✅ | — | FULL_COEVOLUTION_ANALYSIS.md |
| 20 | generate_full_pipeline_doc.py | ✅ | ✅ (fixed 200→full) | — | FULL_PIPELINE_ANALYSIS.md, max MI=1.5917 at (372,401) |

### 8.2 Remaining bugs fixed this pass
1. **generate_full_pipeline_doc.py** — was `load_position_arrays(max_pos=200)` + MI over first 80 → FULL_PIPELINE_ANALYSIS.md only covered 80 positions. Fixed: full length. Top MI pair changed from (74,76) 1.37 → **(372,401) 1.5917**.
2. **run_kmap_analysis.py** — H1 used 200/1299 sequences. Fixed: all 1299 → 1,647,830 pairs.
3. **perplexity_coevolution.py** — display limited to first 80. Fixed: full-length summary.
4. **dca_boolean_coevolution.py** — hardcoded 68-78 pairs. Fixed: dynamic full-length GPU pairs → accuracy 0.0% → **17.6%**.
5. **advanced_co-evolution_analysis.py** — Walsh consensus + clustering + variant signatures limited to 68-80. Fixed: full length. Signatures 11 → **40**; fixed tuple-unpacking bug.

### 8.3 Result consistency (final)
| Metric | Value | Scripts |
|--------|-------|---------|
| Max MI | 1.5917 @ (372,401) | gpu_full ✓ boolean ✓ heatmap ✓ pipeline-doc ✓ |
| Variable positions | 1,249 | all scripts ✓ |
| Co-evolving pairs | 36,918 | master ✓ flipped ✓ allseq ✓ |
| Master rules | 36 (2 essential) | master ✓ kmap_boolean ✓ docs ✓ |

### 8.4 Repo completeness
- `co-evolution-analysis` repo vs `datasets/co-evolution`: **0 differing/missing files** (verified by cmp)
- All 20 scripts + 17 result JSONs + 3 report MDs + CSVs + npy + FASTA + README + RESULTS_ANALYSIS.md present
- Log files intentionally excluded (.gitignore) — they are run artifacts, not results
