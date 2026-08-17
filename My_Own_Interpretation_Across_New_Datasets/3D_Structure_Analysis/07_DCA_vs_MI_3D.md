# 07 — DCA vs MI vs Perplexity: 3D-Contact Prediction Comparison

**Question:** On the 299 Omicron Spike models, which sequence-derived method
best predicts 3D residue-residue contact — plain mutual information, mfDCA
direct information, or the perplexity ratio?

This is the strongest methodological contrast in the project: DCA is the
state-of-the-art contact predictor on diverse benchmarks (PSICOV 0.44
precision), while plain MI is known to be weaker (0.066). Does that hierarchy
hold on the Spike's own redundant, lineage-dominated alignment?

---

## 1. The three methods

### Plain MI (mutual information)
$$MI(i,j) = \sum_{a,b} P(a,b)\log_2\frac{P(a,b)}{P(a)P(b)}$$
Measures total statistical dependence between positions i and j — including
both direct (contact) and indirect (phylogenetic/transitive) coupling. No
correction for phylogeny or transitivity.

### mfDCA direct information (Morcos et al 2011 PNAS)
Mean-field Direct Coupling Analysis: solves the inverse Potts model
$P(\text{seq}) \propto \exp(\sum_{i<j} J_{ij}(a_i,a_j) + \sum_i h_i(a_i))$
via the mean-field approximation (J = covariance⁻¹), applies the **average-
product correction** (APC: subtract the rank-1 component $J_{ij} \approx
\text{row}_i \cdot \text{col}_j^T$ that captures shared phylogenetic bias),
then computes **direct information** $DI = \sum_{a,b} P_{\text{direct}}(a,b)
\cdot J_{ij}^{\text{APC}}(a,b)^2$. DI isolates direct coupling.

### Perplexity ratio (determinism)
$$r(i,j) = \frac{PP(j)}{\text{mean}_a\, PP(j \mid i=a)}, \quad PP = 2^{H}$$
Measures how much knowing residue i reduces the effective number of choices at
j. Ratio > 1 = determinism. Not a coupling decomposition — a conditional-
entropy summary.

---

## 2. Head-to-head on the Omicron models (top-20 pairs each)

**Script:** `scripts/validate_dca_3d.py` → `results/contacts/validate_dca_3d.json`
**Setup:** 299 models, Cβ–Cβ < 8 Å, sep ≥ 3, matched-separation control
(pooled rate 0.1088), exact binomial p.

| Method | n_pairs | contacts/n_models | fraction | enrichment | p | precision@L |
|--------|---------|-------------------|----------|------------|---|-------------|
| **MI** (top-20 full MI) | 20 | 923/4649 | 0.199 | **1.82×** | 1.4e-71 | **0.15** (3/20) |
| **mfDCA** (top-20 DI) | 20 | 162/1719 | 0.094 | 0.87× | 0.98 | 0.05 (1/20) |
| **Perplexity ratio** (ρ vs contact) | 17 | — | — | — | — | ρ=0.670, p=0.006 |

**MI structural pairs** (contact fraction > 0.5): (407,410), (503,507),
(500,503) — all RBD, 10–11× enrichment.
**mfDCA structural pair**: (454,495) — 5.97×, median 8.0 Å. The ONLY DCA pair
that is structural; the rest of the DI top-20 are non-structural (median
9–110 Å).

---

## 3. The key finding: DCA and MI are complementary, not hierarchical

On the diverse PSICOV benchmark (150 proteins, deep diverse MSAs):
- DCA top-L/5 precision = **0.44** (published)
- Plain MI precision = **0.066** (our measurement)
→ **DCA wins by 6.7×** on diverse data.

On the Spike's own 299 models (redundant, lineage-dominated alignment):
- MI top-20 enrichment = **1.82×** (3 structural pairs)
- mfDCA top-20 enrichment = **0.87×** (1 structural pair)
→ **MI wins by 2.1×** on the Spike.

**Why the reversal?**
1. **Redundancy**: the 1,299 entries contain only 299 unique sequences (587
   identical copies of one variant). DCA's covariance inversion is less
   effective on redundant data (low effective N).
2. **Lineage domination**: the Spike's dominant covariation signal is the
   BA.2-lineage phylogenetic drift (the 442–498 cluster, 373–378). DCA's APC
   removes some of this, but the DI ranking is still contaminated; MI's full-
   MI ranking happens to include the RBD loop pairs (500,503)/(503,507) that
   are near-deterministic in the reference frame.
3. **Complementarity**: DCA finds (454,495) — a structural pair that MI's top-20
   misses. MI finds (500,503)/(503,507)/(407,410) — structural pairs that DCA's
   top-20 misses. **Union of both methods' structural pairs = 4**, vs 3 (MI
   alone) or 1 (DCA alone).

---

## 4. The perplexity ratio: the best single correlate

While MI and DCA are pair-ranking methods (top-L precision), the perplexity
ratio is a per-pair scalar whose **value** correlates with 3D contact:

| Metric | Spearman ρ vs contact fraction | p |
|--------|-------------------------------|---|
| **Perplexity ratio** | **0.670** | **0.006** |
| Combined (MI + ratio, rank-normalized) | 0.444 | 0.098 |
| MI alone | 0.266 | 0.338 |
| mfDCA DI | 0.06 (uncorrelated with MI) | — |

The ratio's correlation is the strongest single sequence-metric signal of 3D
proximity. Near-deterministic pairs (ratio > 2: (212,215), (210,215),
(407,410)) are the structurally-real ones. Honest caveat: n=15 pairs (only 21
variable positions); high ratio does not cleanly separate (378,407) and (66,94)
have ratio ~1.7 but 0.0–0.3× contact.

---

## 5. Recommendation

For the Spike (lineage-dominated, redundant): **combine MI + DCA + perplexity
ratio**. The union of structural pairs across methods (4 pairs) exceeds any
single method (1–3). The perplexity ratio is the best single threshold
(ratio > 2 → likely contact), and DCA's (454,495) is a unique contribution.

For diverse benchmarks (PSICOV-class): DCA remains the state of the art; plain
MI is a weak contact predictor (0.066 vs 0.44). The K-map framework's value on
diverse data is enrichment (2.49× mean) and rule extraction, not raw precision.

---

## 6. Reproduce
```bash
$PY scripts/validate_dca_3d.py
# -> results/contacts/validate_dca_3d.json
# Expected: DCA 0.87x (1 structural pair (454,495) 5.97x);
#           MI 1.82x (3 structural pairs); perplexity ρ=0.670
```
