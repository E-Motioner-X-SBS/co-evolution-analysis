# 23. Combined MI + Perplexity Analysis — Applied to All Experiments

## What This Analysis Does

Every experiment in this project originally ranked co-evolving position pairs by **mutual information (MI)** alone. This analysis adds the **perplexity ratio** as a second, complementary lens and combines them. The combined analysis is now embedded in **all 20 analysis scripts** (via `coevolution_shared.combined_pair_scores`), so every experiment reports both the MI-only ranking and the combined MI + perplexity ranking.

## The Two Lenses

### Mutual Information (total correlation)

$$MI(i,j) = \sum_{a,b} P(a,b) \log_2 \frac{P(a,b)}{P(a) P(b)}$$

MI measures how much information two positions share. It captures both direct and indirect (transitive) correlation.

### Perplexity ratio (determinism)

$$r(i,j) = \frac{PP(j)}{\frac{1}{N_a} \sum_a PP(j | i = a)}, \quad PP = 2^H$$

The ratio measures how much knowing residue $i$ narrows the effective number of choices at $j$. Ratio 1 = no constraint; ratio 1.8 = knowing $i$ nearly halves the effective choices at $j$. It is a "determinism" measure, in effective-choice units.

### Combined score

$$\text{combined}(i,j) = \frac{1}{2}\left( \text{ranknorm}(MI) + \text{ranknorm}(r) \right)$$

The combined score is the average of the two rank-normalized scores. It rewards pairs that are both informative (high MI) and deterministic (high ratio).

## Results (corrected aligned data, 21 variable positions)

| Pair | MI | PP ratio | Combined |
|------|-----|----------|----------|
| (378, 407) | 0.792 | 1.74 | **0.875** |
| (18, 26) | 0.802 | 1.73 | **0.844** |
| (66, 94) | 0.791 | 1.73 | 0.812 |
| (210, 215) | 0.753 | 1.79 | 0.812 |
| (212, 215) | 0.398 | **1.84** | 0.750 |
| (373, 378) | **0.807** | 1.59 | 0.719 |

### Key observations

1. **The two lenses disagree on the top pair.** MI alone ranks (373, 378) first (MI = 0.807). The combined score ranks (378, 407) first, because it has nearly as high MI (0.792) AND a higher perplexity ratio (1.74 vs 1.59).

2. **The perplexity ratio rescues under-ranked pairs.** (212, 215) has modest MI (0.398) but the highest ratio (1.84) — evolution at these positions is highly deterministic even though the total shared information is low. The combined score lifts it to rank 5.

3. **Ranking agreement is partial.** The Spearman correlation between MI and the ratio is 0.69 over the 17 scored pairs — strong but not perfect. The combined ranking changes the order of the top-5 pairs.

## Does Combining Improve Prediction?

The LOO-CV experiment (`analysis/combined_mi_perplexity.py`) tested whether ranking by MI, ratio, or combined changes prediction accuracy:

| Ranking | LOO-CV accuracy (top-10 pairs) |
|---------|-------------------------------|
| MI only | 10.34% (318/3074) |
| Perplexity ratio only | 10.34% (318/3074) |
| Combined (rank) | 10.34% (318/3074) |
| Combined (z-score) | 10.34% (318/3074) |

**Result: identical accuracy.** This is because the top-10 pairs under all four rankings overlap almost completely — the pairs with genuinely strong signal are ranked high by all scores. The combined score changes the ORDER but not the SET of top pairs, so aggregate prediction is unchanged.

**Where combining DOES matter:**
- **Rule interpretation:** (212, 215) and (212, 216) have the highest determinism (ratio 1.84) but lower MI. Their Boolean rules describe near-deterministic compensatory mutations — biologically the strongest constraints, even though their total MI is moderate.
- **Robustness:** pairs high in BOTH scores (e.g., (378, 407), (18, 26)) are the most reliable co-evolution signals; a pair high in MI but low in ratio may be an indirect (transitive) correlation.
- **Anti-correlated pairs:** the 93/810 conditioning events with ratio < 1 (from the perplexity deep-dive) identify pairs where conditioning INCREASES uncertainty — the combined score correctly down-weights these.

## What This Means for the Pipeline

1. **MI remains the primary ranking** — the combined score does not change which pairs pass the MI > 0.1 threshold meaningfully.
2. **The perplexity ratio adds interpretation, not ranking power** — it identifies WHICH high-MI pairs are deterministic (rivets) vs indirectly coupled.
3. **For rule extraction (master Boolean):** the 10-12 corrected pairs are all high in both scores; the 3 essential rules come from the deterministic core.
4. **For prediction:** no improvement from combining — the bottleneck is the probabilistic nature of mutation, not the ranking score.

## Scholar Questions and Answers

**Q: Why do MI and perplexity ratio correlate (0.69) but not perfectly?**
A: Both measure dependence, but differently. MI is additive in bits over all residue pairs. The ratio is a quotient of effective-choice counts, dominated by the dominant residues. A pair can have high MI through many moderate couplings, or high ratio through one near-deterministic coupling.

**Q: Why doesn't combining improve LOO-CV?**
A: The top-10 pairs under all rankings are the same set. Prediction accuracy depends on which pairs are INCLUDED, not their internal order. Since the strong pairs are ranked high by all three scores, the sets match.

**Q: Is the combined score better for finding "real" co-evolution?**
A: For robustness, yes: pairs high in both are the most reliable. For pure ranking, MI is sufficient. The ratio adds biological interpretation (determinism) rather than statistical power.

**Q: Where is the combined analysis now?**
A: In all 20 analysis scripts (printed in every run) and in both report markdowns (FULL_COEVOLUTION_ANALYSIS.md and FULL_PIPELINE_ANALYSIS.md). The standalone experiment is `analysis/combined_mi_perplexity.py` with results in `analysis/combined_mi_perplexity.json`.
