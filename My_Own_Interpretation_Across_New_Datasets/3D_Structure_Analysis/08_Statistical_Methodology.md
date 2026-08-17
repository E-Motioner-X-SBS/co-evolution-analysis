# 08 — Statistical Methodology (minute detail)

Every number in the 3D structural validation is produced by a documented
statistical procedure. This file specifies each one exactly.

---

## 1. Contact definition

Two residues $(r_a, r_b)$ in the **same chain** are in contact iff:
$$d(r_a, r_b) = \sqrt{(x_a-x_b)^2 + (y_a-y_b)^2 + (z_a-z_b)^2} < 8.0\;\text{Å}$$
$$|r_a - r_b| \geq 3 \quad \text{(rules analysis; excludes backbone-adjacent)}$$
$$|r_a - r_b| \geq 6 \quad \text{(DCA/benchmark convention; PSICOV)}$$

- **Cβ coordinates** used (Cα for glycine, which has no Cβ). Only the 20
  canonical amino acids are kept.
- The 8 Å Cβ–Cβ cutoff is the standard residue-contact definition in the
  co-evolution literature (Morcos 2011 PNAS; Jones 2012 Bioinformatics).
- The separation cutoff removes trivial backbone-neighbour contacts (residues
  at separation 1–2 are always close regardless of fold). The rules analysis
  uses ≥3 (so short-loop pairs like (212,216), separation 4, are testable);
  the benchmark uses ≥6 (published convention).

## 2. PDB parsing (fixed-width columns)

| Columns (0-indexed) | Field | Used for |
|---|---|---|
| 0–3 | record name | only `ATOM` records |
| 12–15 | atom name | select `CA`/`CB` |
| 16–19 | residue name | map 3-letter → 1-letter (20 canonical) |
| 21 | chain id | blank → chain A (PSICOV repaired models) |
| 22–26 | residue sequence number | residue identity within chain |
| 30–38, 38–46, 46–54 | x, y, z | distance computation |

Multi-chain parsing: homo-trimer models have chains A/B/C; PSICOV natives are
single-chain (blank chain id → A).

## 3. Sequence-position → model-residue mapping (the critical step)

Co-evolution rules are in **alignment columns** (0-based, 1,276 columns). Models
were built from the **gap-stripped target sequence** (SWISS-MODEL strips gaps
before modelling). The mapping:
$$\text{residue}(\text{col}) = \text{col} + 1 - (\text{non-canonical chars before col})$$

**Integrity guarantee** (`rules3d_common.col_to_residue_map`):
```python
if gap_strip(aligned_seq) != target_seq:
    return {}, False   # model SKIPPED, never silently mapped
```
For every one of the 299 unique sequences, the model's `target_sequence` (from
`full_details.json`) was verified equal to the gap-stripped aligned sequence.
Initially 96 skipped because `X` (ambiguous) was not stripped → fixed → **0
skipped** after the fix.

## 4. The random control (how significance is established)

For a rule set with position pairs $P$, per model:
1. Collect the rule-pair residue separations
   $S = \{|res(i) - res(j)| : (i,j) \in P\}$ (after mapping, sep ≥ cutoff).
2. Draw 2,000 random pairs $(a, b)$ per model: $a$ uniform over the model's
   residues, $b = a + s$, $s$ uniform over $S$ (**matched separation
   distribution** — controls for the fact that nearby-sequence pairs are more
   likely to be in contact regardless of co-evolution).
3. Count contacts among random draws → pooled control rate $p_0$.

**Fixed seeds**: `random.Random(42 + model_count)` per model — exactly
reproducible controls.

## 5. The enrichment and p-value

For a pair with $n$ mappable models and $c$ contacts:
$$\text{enrichment} = \frac{c/n}{p_0}, \quad
p\text{-value} = P(X \geq c) \text{ with } X \sim \text{Binomial}(n, p_0)$$

**Exact binomial test** (scipy `binomtest`, `alternative="greater"`):
one-sided — we test whether the pair's residues are significantly *closer* than
random. Enrichment > 1 with p < 0.05 = significant structural contact.

For pooled rule sets: $n = \sum n_k$, $c = \sum c_k$ over all pairs in the set;
same binomial test against the shared $p_0$.

## 6. The Spearman correlation (perplexity-as-predictor)

For the perplexity-ratio analysis (Analysis 8): for all variable-position
pairs with ≥50 mappable models, correlate the per-pair perplexity ratio (a
sequence metric) with the per-pair intra-chain contact fraction (a 3D metric)
via **Spearman rank correlation** (scipy `spearmanr`):
$$\rho = 1 - \frac{6\sum d_i^2}{n(n^2-1)}$$
where $d_i$ are rank differences. Spearman is non-parametric (no normality
assumption). n=15 pairs.

## 7. The PSICOV benchmark precision

For each of the 150 PSICOV proteins:
1. Compute the native contact map from the experimental PDB (Cβ < 8 Å, sep ≥ 6).
2. Take our top-20 MI pairs.
3. **Precision@L** = fraction of top-L pairs that are native contacts.
4. Compare to published PSICOV top-L/5 precision (0.44 mean).

**Enrichment** = (our contact fraction) / (random-pair contact fraction with
matched separation). Reported per-target and aggregated (mean, median, count >
2×/5×).

## 8. Sample-size and multiple-testing caveats (honest)

- **n = 299 models, 21 variable positions, 10–30 pairs per source.** The
  per-pair binomial tests have high power for the strong pairs (p ≤ 1e-8) but
  limited power for weak enrichment (1.2–1.5×). Weak pairs should not be
  over-interpreted.
- **Multiple testing**: 48 unique position pairs tested across sources. A
  Bonferroni correction at α=0.05/48 ≈ 0.001 would still keep the 6 strong
  pairs (p ≤ 1e-8) significant; the weak pairs (p ~ 0.01–0.05) would not
  survive. We report raw p-values and flag the strong pairs as robust.
- **Model vs native**: the Omicron 3D models are homology models (GMQE
  0.69–0.72), not experimental structures. The PSICOV 150 natives provide the
  independent experimental ground truth. Model error could shift a few Å; the
  strong pairs (median 5–9 Å, well inside 8 Å) are robust to this.
- **DSSP/solvent accessibility**: not available in the model payloads → burial
  approximated by local Cβ density (within 12 Å). Secondary-structure context
  unavailable from this source (documented; could be computed with DSSP as a
  follow-up).
