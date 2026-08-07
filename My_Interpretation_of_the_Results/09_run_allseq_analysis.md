> **CORRECTION (Aug 7, 2026):** Two confirmed defects in the original pipeline
> (gap-stripping column misalignment; 8-bit QM wrap-around on 400-cell maps)
> invalidated the biological numbers in this file. See
> [CORRECTION_NOTICE.md](CORRECTION_NOTICE.md) for the verified corrected
> results (21 variable positions, 10 co-evolving pairs, 36 distinct rules (2 essential)) and the
> corrected pipeline (analysis/corrected_pipeline.py). Numbers in this file
> describe the original (buggy) run unless stated otherwise.

# 09. run_allseq_analysis.py - Full Position-Pair K-maps

## What the Program Does

This script performs the position-pair analysis on ALL 1,299 sequences at full length:

1. Compute MI for every position pair within a window of 30 (full 1,263+ positions).
2. Find the top co-evolving pairs.
3. For the top 15 pairs, build a 20 x 20 frequency K-map (rows = residue at position i, columns = residue at position j).
4. Compute coupling constants for each K-map.
5. Run Quine-McCluskey minimization.

The comparison is **between positions, across all sequences**. For each pair (i, j), the joint distribution of residues at those two positions is built from all 1,299 sequences.

```mermaid
flowchart TD
    A[1,299 sequences] --> B[position arrays full length]
    B --> C[GPU MI for pairs within window 30]
    C --> D[34,892 pairs with MI > 0.005]
    D --> E[top 15 pairs]
    E --> F[20x20 frequency K-map per pair]
    F --> G[coupling J = ln(P/P_exp)]
    F --> H[Quine-McCluskey minimization]
```

## The Algorithm

For each pair (i, j) within the window:

1. Extract residues at i and j for every sequence.
2. Count the joint 20 x 20 distribution.
3. Normalize to frequencies.
4. Compute MI.
5. For the top pairs, compute coupling J(a, b) = ln(P(a,b) / (P(a) * P(b))).

## Worked Example: Pair (372, 401)

The top pair from the full-length MI is (372, 401) with MI = 1.5917, matching the full matrix value. Its K-map has a dominant cell (A, N) at frequency 0.6266.

The coupling for (A, N):

```
$$J(A, N) = \ln \frac{P(A,N)}{P(372{=}A) \cdot P(401{=}N)} = \ln \frac{0.6266}{0.63 \times 0.63} = \ln 1.57 \approx 0.45$$
         = ln( 0.6266 / (0.63 * 0.63) )
         = ln( 1.57 ) = 0.45
```

A positive coupling means the pair appears more often than chance.

## Results (full length (1,276 positions), all sequences)

| Metric | Value |
|--------|-------|
| Position pairs with MI > 0.005 | 34,892 |
| Top pair | (372, 401) MI = 1.5917 |
| K-maps minimized | 15 |
| Example coupling | J(372,401) region: +I-A (J=23.03), -T-N (J=-20.95) |

## Inference

The position-pair K-maps show that co-evolution is highly specific: the joint distribution at (372, 401) is dominated by (A, N), while other pairs have different dominant combinations. The positive couplings (like A-L) identify which residue pairs are favored, and the negative couplings (like D-S) identify disfavored combinations.

This script is the bridge between the whole-protein view (scripts 03-05) and the rule extraction (scripts 11-14): it identifies WHICH pairs matter before the Boolean rules are extracted.

## Scholar Questions and Answers

**Q: Why a window of 30?**
A: Co-evolution signal is strongest for nearby positions. The window reduces the pair count from 813,450 to about 37,000 and focuses on locally coupled positions. The full matrix is still available from scripts 06 and 08.

**Q: What does the coupling sign mean?**
A: J > 0: the residue pair is observed more often than expected under independence (co-evolutionary). J < 0: observed less often (anti-correlated). J = 0: independent.

**Q: Why does the top MI differ between scripts (1.5917 vs 1.5690)?**
A: The full MI matrix (scripts 06, 08, and this one) gives MI(372, 401) = 1.5917. The value 1.5690 belongs to the second-ranked pair (401, 404). Script 10 reports (401, 404) as its top pair because it uses a different window and pair-selection convention. The ranking of the top pairs is stable across scripts; the exact values depend on the convention.
