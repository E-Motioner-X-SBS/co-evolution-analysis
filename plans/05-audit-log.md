
## Multi-dataset validation — verification loop (Aug 8, ~21:00-23:30)
16 bugs found & fixed (all in co-evolution-analysis scripts; canonical Spike numbers re-verified after fixes):
1. run_kmap_analysis: header-derived filename with '/' (DYR_ECOLI/1-159) -> nested path crash; sanitize accessions
2-3. boolean_co-evolution: Spike-specific len(seq)>100 filter in 2 functions -> relative >=50% filter
4-5. boolean_co-evolution: hardcoded n_seqs=1299 (x2) -> len(sequences); max_pos=min(200,..) & n_test=min(30,..) -> FULL length via shared GPU MI helper (coevolution_shared.compute_mi_matrix_gpu)
6-7. nary_kmap_co-evolution: same filter + n_seqs=1299 x4 -> fixed; full-length prediction MI
8. create_mi_heatmap: hardcoded focus 60-80 -> adaptive; numpy chained-comparison bug
9-10. position_kmap_coevolution: len>100 filter + min(1299, len) -> fixed
11-12. run_kmap_analysis: n_seqs=1299 x8 -> len; H5 min(200,max_len) -> full; H3 pair cap 10k->100k (logged)
13. advanced_co-evolution_analysis: clustering cap 200 -> 2000 + LOGGED (O(n²) — documented sampling)
14. predictive_constraint_function: hardcoded 800/499 split -> proportional 62/38 (works for any size)
15. run_allseq_analysis: numpy array chained comparison -> element-wise
16. majority_ref tie-break: CPU first-seen vs GPU lowest-index at tied columns (24 mismatches on b2ar) -> unified lowest-code tie-break in shared + 10 inline sites; GPU/CPU consistency now 5/5 PASS
Plus: runner report-generator race (two-wave), canonical dir-name map, GPU memory-fraction fail-fast (0.1), GPU semaphore (2), stale-process hygiene, --retry-failed flag.
Result: unit_tests/gpcrdb/evmutation/pfam/psicov150 re-run from clean state with frozen code (FINAL PASS started ~23:20).

## Resource-management lessons (Aug 9, 00:00-02:00)
- GPU sharing: another user's job (aqua, 71.7 GB) starved ours; torch allocator BLOCKS on
  OOM instead of raising -> jobs deadlocked. Fix: per-process memory fraction 0.1
  (coevolution_gpu.get_device) -> allocation fails fast -> CPU fallback (correct, slower).
- Two concurrent runners fought over the GPU (8 concurrent GPU jobs -> OOM -> CPU fallbacks
  -> 1 h timeouts on deep alignments; pfam 0/63 ok). Fix: ONE runner at a time (sequential
  datasets), GPU semaphore(2) over 8 GPU-heavy scripts, parallel=4.
- Timeout raised 3600 -> 7200 s for deep alignments (Pfam 486k seqs, EVmutation 28k).
- retry flags: --retry-failed (normal) / --retry-failed-cpu (forces CUDA_VISIBLE_DEVICES="")
  for memory-starved situations.
- CPU fallback for 28k-seq alignments is >1 h (timeout) — GPU required for deep datasets.
- Coverage audit: 105/105 runs clean (parser handles train/test totals now).
- Contact validation: b2ar 18.5x, 1whzA 8.1x, 1atzA 1.8x; DHFR insert-column caveat;
  5ht1a 7E2X fragment caveat. a2m reference-frame mapping implemented.
- Phase 6 rules consolidation live (unit_tests + gpcrdb done; nary minimization sub-dict fixed).
- FINAL SEQUENTIAL RUN started ~01:15 (evmutation -> pfam -> psicov150).

## Deep-alignment performance fixes (Aug 9, 15:00-15:40)
- Root cause of 1h timeouts on 28k-seq targets: the post-analysis "Combined MI +
  Perplexity" block (added to all scripts) was O(pairs x seqs) in pure Python:
  * perplexity_ratio used per-pair Counter loops over ALL sequences
  * combined_pair_scores extracted columns via per-pair 28k-element list
    comprehensions (mutual_information path)
- Fixes (coevolution_shared):
  * perplexity_ratio vectorized (numpy bincount) — EXACT equivalence verified
    on random pairs (0 mismatches)
  * combined_pair_scores dense fast path (dense_from_arrays + _mi_dense +
    _ratio_dense) — exact equivalence verified; 14,094 pairs x 28k seqs in
    209 s (was hours / timeouts)
- GPU/CPU note: GPU kernels were never the bottleneck; the CPU-side Python
  loops were. All GPU paths re-verified fast (MI 1.5s at 28k seqs).
- FINAL chain restarted 15:40 with these fixes; evmutation progressing (7 ok
  in first 10 min).

## More GPU hardening (Aug 9, 16:30)
- dca_mf/gpu_full OOM on 28k-seq targets: torch caching-allocator fragmentation
  (20.9 GB reserved-but-unusable vs 39.6 GB fraction cap). Fix: runner sets
  PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True (as the error message itself
  prescribes) + fraction 0.5 -> 0.7.
- evmutation at 56/420 (~1 job/min on deep alignments) — the honest full-scale cost.

## H3 GPU rework (Aug 9, ~17:30) — user: "use the gpu"
- run_kmap_analysis H3 (pairwise K-map distances) blew up CPU RAM (22 GB+ growing)
  and got OOM-SIGKILLed (3x) on deep alignments. CPU pdist O(n^2 x 1024) infeasible.
- Reworked: GPU norm-trick (D = a2+b2-2abT) with CHUNKED statistics — never
  materializes the n x n matrix; accumulates sum/sumsq/min/max/argmin/argmax
  per upper-triangle chunk. Verified on PABP_YEAST: 152,041 sequences,
  11,558,156,820 pairs, exit 0.
- Pfam-scale H3 (486k seqs ~ 4.8e14 FLOPs) will exceed the 2h timeout — honest
  documented limitation (or future chunked-multiprocess).

## User directive (Aug 9, ~23:30): NO 2-hour timeouts — let jobs run
- TIMEOUT_S raised 7200 -> 172800 (48 h ceiling; effectively unlimited for this run).
- All runners restarted with --retry-failed so previously-timed-out/killed jobs re-run
  with the new timeout and the forced-CPU LOO-CV fix.
- Expected: the O(n^2) LOO-CV jobs on 222-486k-seq families may take many hours each
  but will now COMPLETE (or run until they genuinely finish) — no artificial kills.

## Timeout ceiling removed entirely (Aug 10, ~01:30)
- TIMEOUT_S = None — the user's directive "let them run" applied fully: no kills at any
  duration. pfam allseq LOO-CV at 222-486k seqs may take DAYS at the current shared-CPU
  share (~3% each) — accepted as the honest full-scale cost.
- psicov150 core nearly done (1,402/1,650) — full 23-script suite will follow per
  "complete without missing anything".

## MAJOR FIX: single-core pinning discovered (Aug 10, ~14:30)
- User observed no multi-core operation. Root cause: `import torch` pins the process
  to CPU 0 (verified: affinity {0} after import). All analysis children were running
  on ONE core (32 processes on CPU 0; total usage 1.3 cores of 24).
- Fix: OMP_PROC_BIND=FALSE in the runner env (verified: all 24 CPUs retained).
  Relaunched all runners; children now show Cpus_allowed_list 0-23.
- Expected: ~10-30x wall-clock speedup on CPU-bound jobs (LOO-CV, n-ary, clustering).

## Script-quality audit (Aug 11, ~23:50) — user: "make sure scripts are properly written"
- py_compile: 35/35 pass (24 analysis + 11 scripts/ tools)
- pyflakes 3.4.0 (installed via pip fallback; uv cache broken): 0 undefined names,
  0 redefinitions across all 35
- Running processes verified = repo-tracked scripts
- run_psicov_direct.py moved from /store/shuvam/ into scripts/ (md5-identical to
  the running copy)
- psicov150: 3,148/3,150 (2 deep LOO-CVs finishing)

## Final-state documentation (Aug 12)
- ALL DATASETS COMPLETE except 12 pfam items, all giant-family specific:
  * 6 allseq LOO-CV giants (222k-486k seqs): RUNNING (weeks-scale per user's
    no-timeout directive; O(n^2) = 38-190 days of single-core compute)
  * 6 dca_mf/gpu_full giants: DOCUMENTED INFEASIBLE (GPU OOM at 8GB cap; CPU
    needs 89 GB RAM for the one-hot at 486k seqs)
  * 3 generators: running now (fixed SafeDict generator, exit-0 verified on PF00071)
- Key fixes in this final stretch: H3 chunked CPU fallback (was 883GB pdist),
  generator SafeDict (missing-data tolerant), all-sequences filters (no 50% cap),
  pipeline-doc label, nary filter, OMP_PROC_BIND=FALSE (torch CPU pinning).

## Rule→3D validation suite (Aug 13) — 8 proper scripts written & run
- extract_rules_inventory.py (8 pair sources, 36 PIs, 2 essential, 117 flipped, 20+20 motifs)
- validate_3d_complete.py: position-level; validate_residue_rules_3d.py: PI-level;
  validate_flipped_3d.py: steric hypothesis; gray_flip_analysis.py: bit-flip test;
  structural_profile.py: burial/entropy; trimer_interface.py: inter-chain;
  consistency_cross_script.py: cross-script agreement.
- Key results: structural pairs (407,410) 9.3x, (212,215) 3.6x, (210,215) 3.4x,
  (500,503)/(503,507) 10.4-10.5x; non-structural: 442-498 cluster, 373-378;
  Gray-adjacency NOT significant (1.18x p=0.35); housekeeping CONFIRMED
  (entropy-burial rho=-0.17 p=6e-9); trimer: zero inter-chain contacts;
  flipped: partial steric support (210/212-215). Honest synthesis in
  My_Own_Interpretation_Across_New_Datasets/01_3d_structural_validation.md

## 3D Structure Analysis documentation (Aug 14)
- Created My_Own_Interpretation_Across_New_Datasets/3D_Structure_Analysis/
  with 6 documents (~5,900 words): 00_Overview, 01_Datasets_and_Files,
  02_Contact_Computation, 03_Kmaps_and_Boolean_Rules, 04_Structural_Analyses,
  05_Reproducibility.
- ALL paper links HTTP-verified (journal DOIs resolve to publishers —
  OUP/PNAS/Nature 403 = bot-blocking, valid in browser; GitHub/EBI/RCSB/
  SwissModel/DeepMSA = 200).
- Reproducibility: fixed seeds + scipy exact tests; spot-check re-run of
  gray_flip_analysis reproduced identical numbers (38 pairs, 1.18x, p=0.345).

## Scar cleanup + 3D extension + interpretations (Aug 17, 2026)
- **A. Scar cleanup**: 6 stale pre-correction reports corrected in place
  (CO-EVOLUTION_BOOLEAN_FUNCTIONS.md, RESULTS_ANALYSIS.md, FULL_COEVOLUTION_ANALYSIS.md,
  FULL_PIPELINE_ANALYSIS.md/.html, COEVOLUTION_CONSTRAINTS.md). Generators
  (generate_full_pipeline_doc.py, generate_full_analysis_md.py) fixed of 25
  hardcoded stale tokens (236→221 theorems, 108→dynamic, 68-79→full, 4-bit→5-bit,
  17/19→23 scripts, 80→full_len). 3 reports regenerated from corrected JSONs;
  2 hand-edited (RESULTS_ANALYSIS, COEVOLUTION_CONSTRAINTS + sign fix -ln→ln);
  html rebuilt from corrected md. Grep confirms only bug-doc context remains.
- **B. Campaign verified**: coverage 3,735 runs 0 flags; completeness only 6
  documented-infeasible Pfam giants; GPU/CPU 5/5 PASS (diff ~4e-7); py_compile
  35/35. Preservation status written (02_preservation_status.md).
- **C. 3D extended to ALL 14 pair-emitting sources** (was 8): added dca_mf,
  dca_boolean, predictive_constraint, flipped_positions, kmap_boolean, network.
  New finding: mfDCA top-DI pair (454,495) is structural 5.97× (DCA-only).
  New script validate_dca_3d.py: MI 1.82× vs DCA 0.87× on Spike (reversal of
  PSICOV hierarchy); kmap_boolean essential pairs 2.02×; network 1.78×.
- **D. Interpretations written** (3D_Structure_Analysis/): 06_Per_Script (all
  23 scripts, minute detail), 07_DCA_vs_MI, 08_Statistical_Methodology,
  09_Algorithmic_Details. 00_Overview + 05_Reproducibility updated.
- **E. Final audit**: py_compile 35/35; coverage 0 flags; 14 sources confirmed;
  stale grep clean. All green.
