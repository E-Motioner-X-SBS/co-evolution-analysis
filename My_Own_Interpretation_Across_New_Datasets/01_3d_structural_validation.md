# 3D Structural Validation of the Sequence-Deduced Co-evolution Rules

*Question: the co-evolution rules deduced by mathematical logic from FASTA
sequences (K-map + Quine-McCluskey prime implicants, n-ary, flipped, perplexity,
across all 23 scripts) — do they carry 3D structural information?*

*Data: 1,128 SwissModel-built models of the 299 unique Omicron Spikes (best-GMQE
model each); mapping integrity verified per model (model target_sequence ==
gap-stripped aligned sequence, ambiguous columns excluded); contact = Cb-Cb < 8 A,
sep >= 3; random controls matched to the rules' separation distributions;
exact binomial p-values (scipy).*

---

## 1. Headline: the logic is PARTLY structural — and we can say exactly which pairs

### 1a. Position-level (48 unique pairs from 8 source scripts)

| Source (script) | Contact fraction | Enrichment | p-value |
|---|---|---|---|
| position_kmap | 0.199 | **2.08x** | 1e-103 |
| run_allseq | 0.140 | **1.46x** | 1e-31 |
| perplexity (combined) | 0.085 | 0.89x | 0.99 |
| master_boolean (mutation-only) | 0.032 | 0.33x | 1.0 |
| full_length / gpu_full | 0.000 | 0.0x | 1.0 |

**Pairs with significant structural support (per-pair, intra-chain):**

| Pair | Contact frac | Enrichment | p | Region |
|---|---|---|---|---|
| (500, 503) | 1.000 | **10.5x** | ~0 | RBD |
| (503, 507) | 0.997 | **10.4x** | ~0 | RBD |
| (407, 410) | 0.889 | **9.3x** | 1.5e-226 | RBD |
| (210, 214) | 0.384 | **4.0x** | 4.5e-11 | NTD |
| (212, 215) | 0.346 | **3.6x** | 1.4e-09 | NTD |
| (210, 215) | 0.329 | **3.4x** | 3.6e-08 | NTD |

**Pairs with NO structural support** (median 13-27 A): (18,26), (66,94),
(373,378) — the top full-MI pair —, (378,407), (442,448), (442,454), (448,454),
(488,495), (488,498), (495,498) and the whole 479-507 periphery cluster.

### 1b. Residue-level (the 36 prime implicants: IF pos_i=aa_i AND pos_j=aa_j)

Validated on the models where the rule actually FIRES (the residues occur):

| Rule | Fired | Contact frac | Enrichment | p |
|---|---|---|---|---|
| [210=I, 215=P] | 71 | 0.338 | **2.85x** | 1.1e-06 |
| [212=V, 215=P] | 76 | 0.329 | **2.77x** | 1.2e-06 |
| [212=S, 215=G] | 2 | 1.000 | 8.4x | 0.014 |
| all other PIs (442-454, 488-498, 495-498, 212-216...) | — | 0.000 | 0.0x | 1.0 |
| **Pooled 36 PIs** | 780 | 0.069 | 0.58x | 1.0 |

The 2 ESSENTIAL rules never fire in the 299 models (each fires in exactly 1 of
1,299 sequences; in that model the residues are sequentially ADJACENT — a
trivial backbone pair, not an independent contact).

### 1c. Flipped (forbidden) rules — steric hypothesis: PARTIALLY supported

| Pair | Contact frac | Enrichment | p |
|---|---|---|---|
| (210, 215) | 0.329 | **2.57x** | 7.3e-06 |
| (212, 215) | 0.346 | **2.71x** | 6.4e-07 |
| others (442-498 cluster) | 0.000 | 0.0x | 1.0 |
| **Pooled 117 unique forbidden rules** | 0.032 | 0.25x | 1.0 |

Verdict: RARITY overall (forbiddenness is mostly lineage/sequence-space
rarity), but at the two NTD pairs the forbidden combos sit at 3D-close
positions — steric incompatibility is a real constraint there.

### 1d. Gray bit-flip adjacency (the pure K-map prediction)

Observed co-evolving residue pairs are Gray-adjacent (Hamming=1) 23.7% vs
20.0% expected — enrichment 1.18x, **p = 0.345: NOT significant.** The
bit-flip adjacency property of the encoding does NOT predict the observed
co-evolution, and it does NOT correlate with 3D contact (the union analysis
finds Gray-adjacency and contact co-occurring only at (212,215)).

---

## 2. What we are looking at in 3D — the structural profile

- **Housekeeping hypothesis CONFIRMED**: Spearman(entropy, burial) = -0.170,
  p = 6.4e-09. Conserved positions are 52% buried; variable positions 26%
  buried. Variable positions = exposed surface (immune pressure), conserved
  = buried core.
- All co-evolving positions (210-216 NTD, 373-498 RBD/SD) are 1-10% buried —
  surface-exposed, consistent with epitope/immune-driven co-evolution.
- **Trimer interface: ZERO inter-chain contacts.** The co-evolution signal is
  entirely intra-protomer; no pair sits at the inter-protomer interface.
- dssp annotations were not present in the Omicron model payloads (honest
  limitation — secondary-structure context unavailable from that source).

---

## 3. Cross-script consistency — consistency does NOT imply structure

The most-agreed pairs (5-6 sources: the 442-498 cluster) are the LEAST
structural (0.0x). The structurally-strongest pairs (407,410; 210-215;
500,503; 503,507) are agreed by 2-3 sources. All scripts measure CORRELATION;
the dominant variant's lineage covariation is what most of them share — that
signal is phylogenetic, not spatial.

---

## 4. Why the theory is right where it is right, and wrong where it is wrong

1. **Where it works**: MI-based co-evolution genuinely finds structural pairs
   when the covariation is contact-driven — the RBD/NTD pairs (407,410),
   (210,214), (212,215), (210,215), (500,503), (503,507) are deduced purely
   from sequence logic AND are in 3D contact (3.4-10.5x, p <= 1e-8). The
   PSICOV benchmark confirms the mechanism at scale: 2.49x mean enrichment,
   69/150 proteins >2x.
2. **Where it fails**: (a) pairs whose covariation is driven by the dominant
   lineage/immune pressure (442-498 cluster, 373-378) — real sequence signal,
   no structure; (b) plain MI without DCA's APC correction underperforms
   state-of-the-art contact prediction (0.066 vs 0.44 precision on PSICOV);
   (c) Gray-adjacency is NOT a structural predictor (1.18x, p=0.35);
   (d) the n-ary/boolean motifs are protein-wide dipeptide preferences — they
   carry no positional information and cannot be mapped to 3D.
3. **How much**: roughly 6 of 48 position pairs (12%) are strongly structural
   (enrichment > 2.5x, p < 1e-6); ~2/3 of pairs show no structural signal.
   The essential 2-rule core is a rare-event rule (1 firing) — it is
   sequence-logically minimal but structurally untestable at this sample size.

## 5. Files
results/contacts/rules_inventory.json, validate_3d_complete.json,
validate_residue_rules_3d.json, validate_flipped_3d.json,
gray_flip_analysis.json, structural_profile.json, trimer_interface.json,
cross_script_consistency.json, rules_3d_omicron.json

## 6. Perplexity metric — the strongest single correlate of 3D contact

Deeper analysis (perplexity_3d_analysis.py): for every variable-position pair
(window 30, corrected alignment), the perplexity ratio (determinism) was
correlated with the pair's 3D contact fraction across the 299 models:

| Metric | Spearman rho vs contact | p |
|---|---|---|
| Perplexity ratio | **0.670** | **0.006** |
| Combined (MI + ratio) | 0.444 | 0.098 |
| MI alone | 0.266 | 0.338 |

- The DETERMINISM of a co-evolutionary relationship (ratio) is the best
  sequence-metric predictor of 3D proximity — near-deterministic pairs
  ((212,215), (210,215), (407,410)) are the structurally-real ones.
- Ratio quartiles: Q1-Q2 mean contact 0.00 (0.0x), Q3 0.232 (1.80x),
  Q4 0.170 (1.31x) — the top half of ratio space carries the structural
  signal.
- Honest limits: n=15 pairs (only 21 corrected variable positions);
  high ratio does NOT cleanly separate structural from non-structural
  (378,407 and 66,94 have ratio ~1.7 but 0.0-0.3x contact); the
  "ratio rescues" claim yields no NEW structural pairs beyond MI.
