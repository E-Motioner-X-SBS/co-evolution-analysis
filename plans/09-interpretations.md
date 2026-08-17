# 09 — Write minute detailed interpretations into My_Own_Interpretation folders

## Goal
For every one of the 23 scripts: document the algorithm, the exact formulas,
how its pairs/rules are derived, the 3D structural validation of those pairs
(enrichment, p, distances, controls), and the honest interpretation — in
minute statistical/methodological/algorithmic detail.

## Deliverables (in My_Own_Interpretation_Across_New_Datasets/3D_Structure_Analysis/)
- 06_Per_Script_3D_Validation.md — one section per script (all 23), each with:
  algorithm, formulas, pair-derivation, 3D results table, interpretation.
- 07_DCA_vs_MI_3D.md — the mfDCA direct-coupling vs MI vs perplexity 3D
  comparison (the strongest methodological contrast).
- 08_Statistical_Methodology.md — contact definition, mapping integrity,
  matched-separation controls, exact binomial test, Spearman, permutation,
  sample-size caveats, multiple-testing note.
- 09_Algorithmic_Details.md — QM with don't-cares, 32×32 padding, GPU kernels
  (scatter_add bincount, adaptive chunking, H1 popcount), encoding (He 2012,
  5-bit Gray), complexity table.
- Update 00_Overview.md, 04_Structural_Analyses.md with the expanded 14-source
  results; update 05_Reproducibility.md with the new scripts/commands.

## Standard per-script section (the "minute detail" template)
1. What it computes (one paragraph).
2. Algorithm (step list).
3. Formulas (LaTeX: MI, mutation-only MI, perplexity ratio, coupling C, DI,
   QM minterm merge, H1 Hamming, etc. — whichever apply).
4. Pair-derivation (how positions/rules are selected; thresholds; window).
5. 3D validation result (n pairs, contact fraction, enrichment, p, median Å;
   intra vs inter; pooled vs per-pair).
6. Honest interpretation (structural? phylogenetic? motif? limitation?).

## Test
- Every one of the 23 scripts has a section with all 6 template parts.
- Numbers match the JSONs (spot-check 3 scripts against results/contacts).
- Negative results stated as findings, not failures.
