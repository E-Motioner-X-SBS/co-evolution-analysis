# 09. run_allseq_analysis.py — Full Position-Pair K-maps

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

The top pair from the full-length MI is (372, 401) with MI = 1.5690 in this script's computation (matching the 1.5917 from the full matrix; small differences come from the window and the mutation handling). Its K-map has a dominant cell (A, N) at frequency 0.6266.

The coupling for (A, N):

```
J(A, N) = ln( P(A,N) / (P(372=A) * P(401=N)) )
         = ln( 0.6266 / (0.63 * 0.63) )
         = ln( 1.57 ) = 0.45
```

A positive coupling means the pair appears more often than chance.

## Results (full length, all sequences)

| Metric | Value |
|--------|-------|
| Position pairs with MI > 0.005 | 34,892 |
| Top pair | (372, 401) MI = 1.5690 |
| K-maps minimized | 15 |
| Example coupling | J(372,401) region: +A-L, -D-S patterns |

## Inference

The position-pair K-maps show that co-evolution is highly specific: the joint distribution at (372, 401) is dominated by (A, N), while other pairs have different dominant combinations. The positive couplings (like A-L) identify which residue pairs are favored, and the negative couplings (like D-S) identify disfavored combinations.

This script is the bridge between the whole-protein view (scripts 03-05) and the rule extraction (scripts 11-14): it identifies WHICH pairs matter before the Boolean rules are extracted.

## Scholar Questions and Answers

**Q: Why a window of 30?**
A: Co-evolution signal is strongest for nearby positions. The window reduces the pair count from 813,450 to about 37,000 and focuses on locally coupled positions. The full matrix is still available from scripts 06 and 08.

**Q: What does the coupling sign mean?**
A: J > 0: the residue pair is observed more often than expected under independence (co-evolutionary). J < 0: observed less often (anti-correlated). J = 0: independent.

**Q: Why is the top MI slightly different between scripts (1.5690 vs 1.5917)?**
A: Different scripts use different conventions: full MI (all pairs counted) vs mutation-only MI (reference pair excluded), and different window parameters. The pair (372, 401) is the top pair in all of them; the exact value depends on the convention.
