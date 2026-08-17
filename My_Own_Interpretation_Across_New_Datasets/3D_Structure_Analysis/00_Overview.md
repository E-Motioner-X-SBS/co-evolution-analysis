# 3D Structure Analysis — Complete Documentation

**Location of this analysis:** `co-evolution-analysis/` (all scripts, data, results)
**Folder contents (this directory):**
| File | Contents |
|------|----------|
| `00_Overview.md` | Executive summary + the research questions asked |
| `01_Datasets_and_Files.md` | Every dataset, its paper (verified links), files, and formats |
| `02_Contact_Computation.md` | How 3D contacts are computed: PDB parsing, mapping, formulas |
| `03_Kmaps_and_Boolean_Rules.md` | K-map construction, Quine-McCluskey, prime implicants, complexity, GPU/CPU |
| `04_Structural_Analyses.md` | The 9 structural analyses: questions, formulas, results, scripts |
| `05_Reproducibility.md` | Exact commands to reproduce every number |
| `06_Per_Script_3D_Validation.md` | **All 23 scripts** in minute detail: algorithm, formulas, pair-derivation, 3D results |
| `07_DCA_vs_MI_3D.md` | DCA vs MI vs perplexity as 3D-contact predictors (the methodological contrast) |
| `08_Statistical_Methodology.md` | Contact definition, mapping, controls, binomial test, Spearman, caveats |
| `09_Algorithmic_Details.md` | Encoding, MI, GPU kernels, QM, mfDCA, H1, complexity table |
| `10_Boolean_Contact_Circuits.md` | **NEW**: Boolean circuits from FASTA → 3D contact prediction (3 phases; finds salt bridges/H-bonds but doesn't generalize on the Spike) |

---

## Executive summary

The E-Motioner-X-SBS project represents biological sequences as **Karnaugh maps**
(K-maps) via Gray-code encoding, and applies **Boolean minimization
(Quine-McCluskey)** to deduce co-evolution rules from multiple sequence
alignments. This folder documents the **structural validation** of those
sequence-deduced rules: do the mathematically-deduced co-evolutionary pairs
(prime implicants, Boolean functions, flipped rules, perplexity-ranked pairs)
carry real **3-dimensional structural information**?

### What was done, in one paragraph

1. **Sequences → rules**: 1,299 Omicron Spike sequences (GISAID) + 5 external
   benchmark datasets (PSICOV 150, Pfam, EVmutation, unit tests, DeepMSA,
   GPCRdb) were analysed by the 23-script K-map co-evolution pipeline
   (entropy, mutual information, K-map construction, Quine-McCluskey
   minimization, n-ary K-maps, flipped/forbidden K-maps, perplexity,
   constraint functions, DCA-style analysis) — multi-day, multi-core, GPU
   (A100, torch CUDA) runs; **3,687+ verified script runs**.
2. **Rules → 3D models**: 1,128 SWISS-MODEL-built PDB models of the 299
   unique Omicron Spike sequences (3–8 models each) were used to test every
   rule type in 3D. Native experimental structures (PSICOV 150) were used as
   an additional independent structural test bed.
3. **9 structural analyses** answered distinct research questions
   (position-level, residue-level, flipped, Gray bit-flip, housekeeping,
   trimer interface, cross-script consistency, perplexity, contact
   validation) — each with a proper script, matched random controls, and
   exact binomial/Spearman statistics.

### Headline results (all honest, including negatives)

**Extended (Aug 17, 2026): the 3D validation now covers ALL 14 pair-emitting
scripts** (up from 8). New findings: mfDCA's top-DI pair **(454,495) is
structural at 5.97×** (a DCA-only finding MI misses); the essential-rule pairs
are 2.02× enriched; the network edges are 1.78×. See `06_Per_Script_3D_Validation.md`
for the full 23-script breakdown and `07_DCA_vs_MI_3D.md` for the method
comparison.

| Question | Result |
|---|---|
| Do co-evolving positions sit close in 3D? | **Partly yes**: 7 pairs at 3.4–11.3× enrichment (p ≤ 1e-8): (500,503), (503,507), (407,410), (454,495), (210,214), (212,215), (210,215) |
| Do the 36 prime implicants fire on structural contacts? | The 210–215/212–215 rules: 2.8× (p=1e-6); pooled 0.58× — most PIs are non-structural |
| Are forbidden (flipped) pairs sterically constrained? | Partial: (210,215)/(212,215) 2.6–2.7×; overall rarity |
| Is Gray bit-flip adjacency predictive? | **No**: 1.18×, p=0.345 |
| Is the perplexity ratio structural? | **Yes — strongest single metric**: ρ=0.670, p=0.006 |
| Conservation ↔ burial ("housekeeping")? | **Yes**: ρ=−0.170, p=6.4e-09 (conserved buried 52%, variable 26%) |
| Trimer interface? | **No inter-chain contacts** — all signal is intra-protomer |
| Does cross-script agreement imply structure? | **No** — the 5–6-source-consensus pairs are the least structural |
| DCA vs MI on the Spike models? | **MI 1.82× vs DCA 0.87×** (reversal of the PSICOV hierarchy; DCA finds (454,495) MI misses) |

### The two honest limits
1. **The Spike's co-evolution is dominated by immune pressure / lineage
   phylogenetics** — the widely-agreed pairs (442–498 cluster, 373–378) are
   real sequence signal with no spatial meaning (13–27 Å apart).
2. **Plain MI lacks DCA's average-product correction** — raw contact precision
   0.066 vs published PSICOV 0.44 on the 150-protein benchmark, even though
   enrichment (2.49×) proves real structural signal.
