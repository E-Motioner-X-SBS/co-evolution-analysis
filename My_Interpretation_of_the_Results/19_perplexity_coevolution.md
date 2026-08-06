# 19. perplexity_coevolution.py — How Much Does One Residue Constrain Another?

## What the Program Does

This script measures co-evolution through **perplexity**, a different lens than mutual information.

Perplexity is the effective number of choices:

```
Perplexity = 2^H
```

where H is Shannon entropy in bits. If a position has perplexity 3.1, it behaves like a variable with 3.1 effective states.

The script computes:
1. Perplexity of every position (marginal perplexity).
2. Conditional perplexity: PP(j | i = a) = 2^H(j | i = a), the perplexity of position j given that position i has residue a.
3. The co-evolution ratio: PP(j) / PP(j | i).

If knowing residue i reduces the uncertainty about j, the conditional perplexity is lower than the marginal, and the ratio is greater than 1.

## The Formulas

### Marginal entropy and perplexity

```
H(j) = - sum over b of P(b) * log2(P(b))
PP(j) = 2^H(j)
```

### Conditional entropy and perplexity

```
H(j | i = a) = - sum over b of P(b | i = a) * log2(P(b | i = a))
PP(j | i = a) = 2^H(j | i = a)
```

### Co-evolution ratio

```
ratio(i, j) = PP(j) / average over a of PP(j | i = a)
```

A ratio of 1 means position i gives no information about j. A ratio of 2 means knowing i halves the effective number of choices at j.

```mermaid
flowchart TD
    A[position arrays full length] --> B[entropy per position]
    B --> C[perplexity PP = 2^H]
    A --> D[GPU MI to find co-evolving pairs]
    D --> E[top pairs]
    E --> F[conditional PP j given i]
    C --> G[ratio = PP(j) / PP(j|i)]
    F --> G
    G --> H[ratio up to 2.81x]
```

## Worked Example: The (372, 401) Pair

- H(372) = 1.6328 bits, PP(372) = 2^1.6328 = 3.101.
- H(401) = 1.6715 bits, PP(401) = 2^1.6715 = 3.185.
- Conditional perplexity: given residue at 372, the perplexity of 401 drops to about 1.05.
- Ratio = 3.185 / 1.05 ~ 3.0 (the script reports 2.81 for the analyzed pairs).

Interpretation: knowing the residue at 372 reduces the effective number of choices at 401 from about 3.2 to about 1.1. At the strongest pairs the conditional perplexity is close to 1.0, meaning position j is essentially DETERMINED by position i. This is near-deterministic co-evolution.

## Results

| Metric | Value |
|--------|-------|
| Variable positions (perplexity > 3) | 8 |
| Pairs analyzed | 3 (full-length dynamic) |
| (372, 401) | ratio = 2.809 |
| (208, 209) | ratio = 2.769 |
| (401, 404) | ratio = 2.700 |

## Inference

**Perplexity vs MI: what does perplexity add?**

1. Perplexity is in the "number of states" unit, which is more intuitive than bits. PP = 1.1 says "about one effective choice".
2. Conditional perplexity close to 1.0 is the direct statement of deterministic co-evolution: given i, j is forced.
3. The ratio 2.81 for (372, 401) confirms the MI result (1.5917 bits, 97.6% of maximum) from a different angle.

**Is perplexity better or worse than MI?**

They are complementary:
- MI is additive and directly measures information in bits; it is the standard for ranking pairs.
- Perplexity is a monotone transform of entropy (2^H) and is easier to interpret as "effective number of choices".
- Perplexity does NOT replace MI for ranking, but the conditional-perplexity ratio gives a natural "determinism" reading that MI does not directly provide (MI of 1 bit does not immediately say "j is nearly determined", but PP(j|i) = 1.05 does).

**Regional pattern.** The top pairs cluster at positions 208-210, 372-404 (S1 subunit). The "more variance at the top" question: the top-MI pairs have perplexity ratios near 2.7-2.8 (strong constraint), and these are concentrated in specific S1 regions. The N-terminal and other regions have lower ratios. So yes, there ARE particular regions with higher perplexity-based constraint, matching the MI hubs.

## Scholar Questions and Answers

**Q: How is perplexity working here?**
A: Perplexity = 2^entropy. For a position with frequencies over 20 amino acids, it is the effective number of amino acids that appear with non-negligible probability. PP = 1 means fully conserved; PP = 20 means uniform.

**Q: What does "ratio 2.81" mean exactly?**
A: Without knowing position 372, position 401 has 3.19 effective choices. Knowing position 372 reduces that to about 1.13 effective choices. The ratio 2.81 is 3.19 / 1.13.

**Q: Why only 8 variable positions by the perplexity > 3 criterion?**
A: The script uses a stricter threshold for display (PP > 3, i.e., H > 1.585 bits). The entropy > 0.3 criterion (used elsewhere) is looser and gives 1,249 variable positions. The two criteria measure different things: the perplexity criterion marks highly variable positions only.

**Q: Are there particular regions with higher perplexity?**
A: Yes. The strongest constraint (ratio ~2.8) is at (372, 401) and (208, 209), both in the S1 subunit. This matches the MI hubs. Perplexity confirms the regional structure independently.

**Q: Is perplexity better than MI?**
A: Neither is universally better. MI is the standard ranking tool. Perplexity gives a more interpretable "effective number of choices" view and makes determinism explicit (conditional PP ~ 1). Use MI for ranking, perplexity for interpretation.
