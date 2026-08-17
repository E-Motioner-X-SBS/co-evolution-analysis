# Understanding: Co-evolution Analysis Project (Aug 7, 2026)

## What the project is
Karnaugh-map (K-map) Boolean minimization applied to SARS-CoV-2 Spike co-evolution.
Core idea: encode amino acids via Gray codes (binary 5-bit, or base-20 He-2012),
build frequency K-maps over position pairs, threshold to Boolean functions,
Quine-McCluskey minimize -> co-evolutionary inference rules.

## Repositories (all under /store/shuvam/E-motioner-X-SBS/)
- `co-evolution-analysis/` — repo copy of the analysis (21 scripts + docs)
- `datasets/co-evolution/` — LIVE working copy; scripts run here; **final results here**
  (logs/ has final Aug 7 runs; markdowns regenerated Jul 31)
- `lean_proofs/` — Lean 4 formal proofs of encoding framework (106 theorems)
- `n-ary-kmap/` — base-N K-map generalization (115 Lean theorems) + src/nkmap encoder
- `kmap-sbm-validation/` — empirical SBM MD validation (H1-H6) + prime_implicants.py (QM)

## Data
- `Spike_protein.aln-fasta`: 1,299 Omicron Spike sequences, 1276 positions
- Gaps '-' stripped; He 2012 base-20 encoding via `Base20AminoEncoder(version=1)` from n-ary-kmap

## Pipeline (shared module coevolution_shared.py is single source of truth)
1. entropy per position (H > 0.3 = variable) -> 1249/1276 variable
2. mutual information per position pair (mutation-only MI variant) -> pairs
3. 20x20 frequency/mutation K-map per pair; threshold; DC for reference pair
4. Quine-McCluskey (kmap_sbm.analysis.prime_implicants.boolean_minimize_kmap)
5. essential prime implicants = inference rules; coupling J = ln(P/P_exp); sigmoid prediction

## Final results (datasets/co-evolution, verified Aug 7)
- variable_position_summary.json (fresh Aug 7 03:17): 1249 var pos, top 20 pairs
  top pair (413,427) MI=2.3524, 271 mutations, 12 on-set, 11 PIs, 10 essential
- master_boolean_summary.json (Jul 27): 36,918 pairs (MI>0.1), 162 PIs, 152 essential, 152 rules
- boolean results: binary K-map 93 on-set/1024 cells, 70 PIs, 38 essential, 50.7% accuracy
- nary: 93 on-set/400, 73 PIs, 42 essential
- LOO-CV: 7.26% (198/2726) — lineage-specific rules, doesn't generalize
- DCA Boolean: 0% (singular covariance)
- constraint fn: 5.8% prediction; coupling has BOTH signs (see discrepancy D2)
- perplexity: ratio up to 2.81x (372,401); conditional PP ~1.0 = near-deterministic
- network: 1249 nodes, 35098 edges, hub pos 85
- full-length: 4949 high-MI pairs; top entropy pos 852 (1.774)
- GPU (Jul 17, stale convention): 113.2s A100

## Discrepancies found (important!)
- D1: CO-EVOLUTION_BOOLEAN_FUNCTIONS.md claims MI 8.75-8.83 for pairs 68-79 — IMPOSSIBLE
  (max MI for 20x20 joint = log2(20) = 4.32). Early doc, stale/wrong; final MI max 2.35.
- D2: FULL_COEVOLUTION_ANALYSIS.md "ALL coupling constants C < 0" contradicts
  constraint_function_results (positive C up to 5.22) and my recomputation
  (pair 413,427: 302 J>0 cells, max 23.0). GPU top_co empty due to freq>0.001 filter +
  old sign convention. Docs disagree; recomputed J has both signs.
- D3: H1 enrichment differs per script: 1.34x (kmap_results), 1.14x (FULL_PIPELINE §9),
  0.37x (GPU), 1.05x corrected (SBM). Different pair sets/definitions.
- D4: lean theorem counts: 103 (index.md) vs 106 (README/agents.md) vs 236 total
  (106+115+15). index.md table stale.
- D5: kmap_boolean_coevolution/COEVOLUTION_KMAP_BOOLEAN.md shows 0 on-set/0 PIs —
  stale artifact from older script version; real rules in master_boolean + variable_position.
- D6: FULL_PIPELINE abstract says H1 0.18x vs §9 1.14x — internal contradiction.
- D7: "positions 68-79 in signal peptide" claim in CO-EVOLUTION_BOOLEAN_FUNCTIONS.md
  likely wrong region labeling (alignment-indexed, not Spike numbering).

## Repo copies vs datasets copies
- FULL_COEVOLUTION_ANALYSIS.md / FULL_PIPELINE_ANALYSIS.md differ only in generation date
  (Jul 27 repo vs Jul 31 datasets) — content same.
- coevolution_shared.py identical. boolean_co-evolution.py and
  variable_position_coevolution.py differ (datasets = optimized vectorized versions, Aug 7).
