# 23. Combined MI + Perplexity Analysis - Applied to All Experiments

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

2. **The perplexity ratio rescues under-ranked pairs.** (212, 215) has modest MI (0.398) but the highest ratio (1.84) - evolution at these positions is highly deterministic even though the total shared information is low. The combined score lifts it to rank 5.

3. **Ranking agreement is partial.** The Spearman correlation between MI and the ratio is 0.69 over the 17 scored pairs - strong but not perfect. The combined ranking changes the order of the top-5 pairs.

## Worked Example: What MI and Perplexity Actually See

Both scores are computed from the same joint distribution. The difference is
the question each one asks. This section shows the exact data behind the
three representative pairs.

### The MI top pair: (373, 378), MI = 0.807 bits

Position 373: entropy 0.872 bits, perplexity 1.83 (F or L dominates).
Position 378: entropy 0.811 bits, perplexity 1.75 (A or T dominates).

Joint distribution over 1,294 sequences:

| 373 | 378 | Count | Frequency |
|-----|-----|-------|-----------|
| F | A | 970 | 75.0% |
| L | T | 311 | 24.0% |
| S | T | 12 | 0.9% |
| S | A | 1 | 0.1% |

**What MI sees:** F always co-occurs with A, L always with T. If the two
positions were independent, we would also see F+T and L+A combinations, but
these are essentially absent. The observed joint is far from the product of
the marginals, so the information shared is large. The theoretical maximum
is min(H(373), H(378)) = 0.811 bits; the observed 0.807 is 99.5% of that
maximum. This is why (373, 378) is the MI top pair: the two positions are
locked together, almost one binary variable expressed twice.

### The perplexity-ratio top pair: (212, 215), ratio = 1.84

Position 215: marginal perplexity PP(215) = 1.842 (H = 0.881 bits).

Conditional perplexity of 215 given each residue at 212:

| Residue at 212 | Sequences | PP(215 \| 212) | Effective choices left |
|----------------|-----------|----------------|------------------------|
| V | 279 | 1.000 | 1 |
| I | 13 | 1.000 | 1 |
| L | 6 | 1.000 | 1 |
| S | 2 | 1.000 | 1 |

Average conditional perplexity = 1.000.
Ratio = 1.842 / 1.000 = 1.842.

**What the ratio sees:** once residue at 212 is known, position 215 has
exactly 1 effective choice left - it is fully determined. This is perfect
determinism. The ratio 1.84 means "knowing 212 collapses 215 from 1.84
effective choices to 1."

**Why MI ranks it lower (0.398):** position 212 is 75% gaps, so only 300
of 1,299 sequences carry a valid (non-gap) residue at 212. The total
information is diluted by the large gap fraction, so MI is modest even
though the determinism among valid sequences is perfect.

### The combined top pair: (378, 407), combined = 0.875

Position 407: marginal PP(407) = 1.758 (H = 0.814 bits).

Conditional perplexity of 407 given each residue at 378:

| Residue at 378 | Sequences | PP(407 \| 378) |
|----------------|-----------|----------------|
| A | 973 | 1.021 |
| T | 324 | 1.000 |

Average conditional perplexity = 1.011.
Ratio = 1.758 / 1.011 = 1.740.

**What the combined score sees:** (378, 407) has MI = 0.792 (nearly as
high as the MI top pair) AND ratio = 1.74 (nearly as deterministic as the
ratio top pair). It is strong on both axes, which is why it is the top
combined pair. Pairs like (212, 215) are rescued by the ratio but hurt by
low MI; pairs like (373, 378) are hurt slightly by lower determinism
(conditional PP leaves 1.10 choices on average, not 1.00).

### Summary of the three lenses

| Pair | MI | PP ratio | What it is |
|------|-----|----------|------------|
| (373, 378) | **0.807** | 1.59 | Most total shared information; near-perfect 2-state lock |
| (212, 215) | 0.398 | **1.84** | Perfect determinism among valid sequences; gap-diluted MI |
| (378, 407) | 0.792 | 1.74 | Strong on both axes: the combined top pair |

**Top by MI** = the pair sharing the most total information (bits).
**Top by ratio** = the pair where knowing one position most completely
determines the other (collapse to 1 effective choice).
**Top by combined** = the pair that is informative AND deterministic - the
most robust co-evolution constraint.

**Biological reading:** a top-MI pair says "these two positions evolve
together" (compensatory mutation, structural coupling). A top-ratio pair
says "these two positions are riveted - if one changes, the other is
forced to a specific residue." The combined top is the pair you would bet
on as the most reliable constraint.

**Caveat:** (212, 215) appears perfectly deterministic partly because 75%
of sequences have a gap at 212; the ratio is computed over the 300
non-gap sequences. Using both lenses together is the honest approach: MI
captures the total signal, the ratio captures determinism, the combined
score captures the robust core.

## Does Combining Improve Prediction?

The LOO-CV experiment (`analysis/combined_mi_perplexity.py`) tested whether ranking by MI, ratio, or combined changes prediction accuracy:

| Ranking | LOO-CV accuracy (top-10 pairs) |
|---------|-------------------------------|
| MI only | 10.34% (318/3074) |
| Perplexity ratio only | 10.34% (318/3074) |
| Combined (rank) | 10.34% (318/3074) |
| Combined (z-score) | 10.34% (318/3074) |

**Result: identical accuracy.** This is because the top-10 pairs under all four rankings overlap almost completely - the pairs with genuinely strong signal are ranked high by all scores. The combined score changes the ORDER but not the SET of top pairs, so aggregate prediction is unchanged.

**Where combining DOES matter:**
- **Rule interpretation:** (212, 215) and (212, 216) have the highest determinism (ratio 1.84) but lower MI. Their Boolean rules describe near-deterministic compensatory mutations - biologically the strongest constraints, even though their total MI is moderate.
- **Robustness:** pairs high in BOTH scores (e.g., (378, 407), (18, 26)) are the most reliable co-evolution signals; a pair high in MI but low in ratio may be an indirect (transitive) correlation.
- **Anti-correlated pairs:** the 93/810 conditioning events with ratio < 1 (from the perplexity deep-dive) identify pairs where conditioning INCREASES uncertainty - the combined score correctly down-weights these.

## What This Means for the Pipeline

1. **MI remains the primary ranking** - the combined score does not change which pairs pass the MI > 0.1 threshold meaningfully.
2. **The perplexity ratio adds interpretation, not ranking power** - it identifies WHICH high-MI pairs are deterministic (rivets) vs indirectly coupled.
3. **For rule extraction (master Boolean):** the 10-12 corrected pairs are all high in both scores; the 3 essential rules come from the deterministic core.
4. **For prediction:** no improvement from combining - the bottleneck is the probabilistic nature of mutation, not the ranking score.

## Scholar Questions and Answers

**Q: Why do MI and perplexity ratio correlate (0.69) but not perfectly?**
A: Both measure dependence, but differently. MI is additive in bits over all residue pairs. The ratio is a quotient of effective-choice counts, dominated by the dominant residues. A pair can have high MI through many moderate couplings, or high ratio through one near-deterministic coupling.

**Q: Why doesn't combining improve LOO-CV?**
A: The top-10 pairs under all rankings are the same set. Prediction accuracy depends on which pairs are INCLUDED, not their internal order. Since the strong pairs are ranked high by all three scores, the sets match.

**Q: Is the combined score better for finding "real" co-evolution?**
A: For robustness, yes: pairs high in both are the most reliable. For pure ranking, MI is sufficient. The ratio adds biological interpretation (determinism) rather than statistical power.

**Q: Where is the combined analysis now?**
A: In all 20 analysis scripts (printed in every run) and in both report markdowns (FULL_COEVOLUTION_ANALYSIS.md and FULL_PIPELINE_ANALYSIS.md). The standalone experiment is `analysis/combined_mi_perplexity.py` with results in `analysis/combined_mi_perplexity.json`.
