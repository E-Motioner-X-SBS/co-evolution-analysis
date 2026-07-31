# Co-evolution Boolean Functions for SARS-CoV-2 Spike Protein

## Master Boolean Function

**Dataset:** 1,299 SARS-CoV-2 Spike protein sequences  
**Encoding:** Base-20 (He 2012 order)  
**Method:** Variable-position K-map with don't-care conditions → Quine-McCluskey minimization  
**Variable positions:** 57 out of 80 analyzed  
**Co-evolutionary pairs:** 1,161 (MI > 0.1)  
**Total prime implicants:** 108  
**Essential prime implicants:** 108 (all essential)  

---

## Master Boolean Equation

```
f(pos_i, pos_j, aa_i, aa_j) = Rule_1 OR Rule_2 OR ... OR Rule_108
```

Each rule is an AND of residue conditions at two positions. The function returns **1 (co-evolutionary)** when ANY prime implicant matches.

---

## Variables

- `s3, s2, s1, s0` = binary representation of residue at position `pos_i` (0-19)
- `t3, t2, t1, t0` = binary representation of residue at position `pos_j` (0-19)
- `~s3` = NOT s3 (bit is 0)
- `s3` = bit is 1

---

## Co-evolutionary Position Pairs

| Rank | Pos i | Pos j | MI | Mutations | Ref i | Ref j |
|------|-------|-------|-----|-----------|-------|-------|
| 1 | 76 | 77 | 8.83 | 458 | D | N |
| 2 | 74 | 79 | 8.82 | 459 | R | V |
| 3 | 71 | 75 | 8.82 | 457 | G | F |
| 4 | 72 | 75 | 8.82 | 457 | T | F |
| 5 | 78 | 79 | 8.82 | 459 | P | V |
| 6 | 71 | 79 | 8.82 | 459 | G | V |
| 7 | 72 | 79 | 8.82 | 459 | T | V |
| 8 | 73 | 79 | 8.82 | 459 | K | V |
| 9 | 74 | 75 | 8.82 | 457 | R | F |
| 10 | 73 | 75 | 8.81 | 457 | K | F |
| 11 | 77 | 79 | 8.80 | 459 | N | V |
| 12 | 69 | 79 | 8.78 | 459 | T | V |
| 13 | 69 | 70 | 8.77 | 457 | T | N |
| 14 | 68 | 79 | 8.75 | 459 | G | V |
| 15 | 68 | 75 | 8.75 | 457 | G | F |

---

## ALL 108 Inference Rules

### Position Pair (68, 75) — MI = 8.75

| Rule | Boolean Expression | Amino Acids | Type |
|------|-------------------|-------------|------|
| 1 | `~s3 & s2 & ~s1 & ~s0 & t3 & ~t2 & t0` | (M, V) | MUTATION |
| 2 | `~s3 & s2 & ~s1 & ~s0 & t3 & ~t2 & ~t1 & t0` | (M, D) | MUTATION |
| 3 | `s3 & s2 & s1 & ~s0 & ~t3 & t2 & t1 & t0` | (R, W) | MUTATION |
| 4 | `s3 & s2 & s1 & ~s0 & t3 & t2 & t1 & ~t0` | (R, R) | MUTATION |
| 5 | `~s3 & ~s2 & ~s1 & ~s0 & ~t3 & ~t2 & ~t1 & ~t0` | (A, A) | MUTATION |
| 6 | `~s3 & ~s2 & ~s1 & s0 & ~t3 & ~t1 & t0` | (V, I) | MUTATION |
| 7 | `~s3 & ~s2 & ~s1 & s0 & ~t3 & t2 & t1 & ~t0` | (V, Q) | MUTATION |
| 8 | `s3 & ~s2 & ~s1 & ~s0 & ~t3 & t0` | (E, I) | MUTATION |

**Interpretation:** When position 68 mutates to M, R, A, V, or E, position 75 must co-evolve to V, D, W, R, A, I, Q, or I respectively.

---

### Position Pair (68, 79) — MI = 8.75

| Rule | Boolean Expression | Amino Acids | Type |
|------|-------------------|-------------|------|
| 9 | `~s3 & s2 & ~s1 & ~s0 & ~t3 & ~t1 & t0` | (M, L) | MUTATION |
| 10 | `s3 & s2 & s1 & ~s0 & t3 & t2 & ~t1 & ~t0` | (R, I) | MUTATION |
| 11 | `s3 & s2 & s1 & ~s0 & t3 & t2 & t1 & ~t0` | (R, R) | MUTATION |
| 12 | `s3 & s2 & s1 & s0 & ~t3 & t1 & t0` | (S, D) | MUTATION |
| 13 | `~s3 & ~s2 & ~s1 & s0 & t3 & t2 & ~t0` | (L, S) | MUTATION |
| 14 | `~s3 & ~s2 & ~s1 & s0 & t3 & ~t2 & t1 & t0` | (V, R) | MUTATION |
| 15 | `~s3 & s2 & ~s1 & ~s0 & ~t3 & t2 & t1 & t0` | (M, N) | MUTATION |
| 16 | `s3 & ~s2 & s1 & s0 & t3 & t2 & t1 & t0` | (W, H) | MUTATION |
| 17 | `s3 & ~s2 & ~s1 & ~s0 & ~t3 & t0` | (E, I) | MUTATION |

---

### Position Pair (69, 70) — MI = 8.77

| Rule | Boolean Expression | Amino Acids | Type |
|------|-------------------|-------------|------|
| 18 | `~s3 & s2 & ~s1 & ~s0 & ~t3 & ~t2 & t1 & t0` | (M, N) | MUTATION |
| 19 | `~s3 & s2 & ~s1 & s0 & ~t3 & t2 & t1 & t0` | (N, S) | MUTATION |
| 20 | `s3 & s2 & s1 & ~s0 & ~t3 & t2 & t1 & t0` | (R, S) | MUTATION |
| 21 | `~s3 & ~s2 & ~s1 & s0 & ~t3 & t2 & t1 & ~t0` | (I, L) | MUTATION |
| 22 | `~s3 & s2 & ~s1 & ~s0 & t3 & ~t2 & t1 & t0` | (M, K) | MUTATION |
| 23 | `s3 & ~s2 & s1 & s0 & t3 & ~t2 & t1 & t0` | (W, K) | MUTATION |
| 24 | `s3 & ~s2 & ~s1 & ~s0 & t3 & t2 & ~t1 & t0` | (E, D) | MUTATION |
| 25 | `s3 & ~s2 & ~s1 & ~s0 & t3 & ~t2 & t1 & t0` | (E, H) | MUTATION |

---

### Position Pair (69, 79) — MI = 8.78

| Rule | Boolean Expression | Amino Acids | Type |
|------|-------------------|-------------|------|
| 26 | `~s3 & s2 & ~s1 & ~s0 & t3 & ~t1 & t0` | (M, F) | MUTATION |
| 27 | `~s3 & s2 & ~s1 & s0 & t3 & t2 & ~t0` | (N, W) | MUTATION |
| 28 | `s3 & ~s2 & ~s1 & s0 & t3 & t2 & t1 & ~t0` | (K, R) | MUTATION |
| 29 | `s3 & s2 & s1 & ~s0 & t3 & t2 & t1 & ~t0` | (R, R) | MUTATION |
| 30 | `~s3 & ~s2 & ~s1 & ~s0 & ~t3 & t2 & t1 & t0` | (A, S) | MUTATION |
| 31 | `~s3 & s2 & ~s1 & ~s0 & ~t3 & ~t2 & t0` | (M, A) | MUTATION |
| 32 | `s3 & ~s2 & ~s1 & ~s0 & ~t3 & t2 & t1 & t0` | (E, I) | MUTATION |
| 33 | `s3 & ~s2 & ~s1 & ~s0 & t3 & t1 & ~t0` | (E, R) | MUTATION |

---

### Position Pair (71, 75) — MI = 8.82

| Rule | Boolean Expression | Amino Acids | Type |
|------|-------------------|-------------|------|
| 34 | `~s3 & s2 & ~s1 & ~s0 & t3 & ~t2 & ~t1 & t0` | (M, D) | MUTATION |
| 35 | `r3 & r2 & ~r1 & r0 & c3 & c2 & c1 & ~c0` | (Y, W) | MUTATION |
| 36 | `s3 & s2 & s1 & ~s0 & t3 & ~t1 & t0` | (R, Q) | MUTATION |
| 37 | `~s3 & ~s2 & ~s1 & ~s0 & ~t3 & ~t2 & t1 & t0` | (A, S) | MUTATION |
| 38 | `~s3 & ~s2 & s1 & ~s0 & t3 & t2 & t1 & ~t0` | (L, Q) | MUTATION |
| 39 | `r3 & r2 & ~r1 & r0 & t3 & t2 & ~t1 & t0` | (Y, K) | MUTATION |
| 40 | `s3 & ~s2 & ~s1 & ~s0 & t3 & t2 & t1 & t0` | (E, H) | MUTATION |

---

### Position Pair (71, 79) — MI = 8.82

| Rule | Boolean Expression | Amino Acids | Type |
|------|-------------------|-------------|------|
| 41 | `r3 & r2 & ~r1 & r0 & ~t3 & ~t2 & t1 & t0` | (Y, N) | MUTATION |
| 42 | `~s3 & s2 & ~s1 & ~s0 & ~t3 & ~t1 & t0` | (M, L) | MUTATION |
| 43 | `s3 & s2 & s1 & ~s0 & t3 & t2 & t1 & ~t0` | (R, R) | MUTATION |
| 44 | `~s3 & ~s2 & ~s1 & s0 & t3 & ~t2 & t1 & ~t0` | (I, Y) | MUTATION |
| 45 | `~s3 & ~s2 & ~s1 & s0 & t3 & t2 & t1 & ~t0` | (I, K) | MUTATION |
| 46 | `s3 & ~s2 & s1 & s0 & t3 & t2 & t1 & t0` | (W, H) | MUTATION |
| 47 | `s3 & ~s2 & ~s1 & ~s0 & ~t3 & t2 & t1 & t0` | (E, F) | MUTATION |

---

### Position Pair (72, 75) — MI = 8.82

| Rule | Boolean Expression | Amino Acids | Type |
|------|-------------------|-------------|------|
| 48 | `~s3 & ~s2 & ~s1 & ~s0 & t3 & ~t2 & t1 & t0` | (A, K) | MUTATION |
| 49 | `~s3 & s2 & ~s1 & ~s0 & ~t3 & t2 & t0` | (M, I) | MUTATION |
| 50 | `s3 & ~s2 & s1 & s0 & ~t3 & t2 & t1 & t0` | (W, Y) | MUTATION |
| 51 | `~s3 & s2 & ~s1 & s0 & ~t3 & t2 & t1 & t0` | (N, W) | MUTATION |
| 52 | `~s3 & ~s2 & s1 & ~s0 & ~t3 & ~t2 & ~t1 & t0` | (L, V) | MUTATION |
| 53 | `~s3 & ~s2 & ~s1 & s0 & t3 & ~t1 & ~t0` | (F, A) | MUTATION |
| 54 | `s3 & ~s2 & ~s1 & ~s0 & t3 & ~t2 & t1 & t0` | (E, Q) | MUTATION |

---

### Position Pair (72, 79) — MI = 8.82

| Rule | Boolean Expression | Amino Acids | Type |
|------|-------------------|-------------|------|
| 55 | `~s3 & ~s2 & ~s1 & s0 & ~t3 & t2 & t0` | (V, S) | MUTATION |
| 56 | `~s3 & s2 & ~s1 & ~s0 & t3 & ~t2 & ~t1 & t0` | (M, D) | MUTATION |
| 57 | `~s3 & ~s2 & ~s1 & ~s0 & ~t3 & t2 & t1 & t0` | (A, Y) | MUTATION |
| 58 | `~s3 & ~s2 & s1 & ~s0 & t3 & ~t2 & t1 & ~t0` | (L, Q) | MUTATION |
| 59 | `~s3 & s2 & ~s1 & ~s0 & ~t3 & ~t2 & t0` | (M, A) | MUTATION |
| 60 | `s3 & ~s2 & ~s1 & ~s0 & t3 & t1 & ~t0` | (E, R) | MUTATION |

---

### Position Pair (73, 75) — MI = 8.81

| Rule | Boolean Expression | Amino Acids | Type |
|------|-------------------|-------------|------|
| 61 | `~s3 & ~s2 & s1 & s0 & ~t3 & ~t2 & t0` | (L, I) | MUTATION |
| 62 | `~s3 & ~s2 & s1 & s0 & ~t3 & t2 & t1 & ~t0` | (L, K) | MUTATION |
| 63 | `r3 & r2 & ~r1 & r0 & ~t3 & t2 & t1 & ~t0` | (Y, S) | MUTATION |
| 64 | `s3 & ~s2 & ~s1 & ~s0 & t3 & t2 & t1 & t0` | (H, Y) | MUTATION |
| 65 | `s3 & ~s2 & ~s1 & s0 & t3 & t2 & ~t1 & ~t0` | (K, S) | MUTATION |
| 66 | `s3 & s2 & s1 & ~s0 & ~t3 & t2 & t1 & t0` | (R, H) | MUTATION |
| 67 | `~s3 & s2 & ~s1 & ~s0 & ~t3 & t2 & ~t1 & t0` | (M, R) | MUTATION |

---

### Position Pair (73, 79) — MI = 8.82

| Rule | Boolean Expression | Amino Acids | Type |
|------|-------------------|-------------|------|
| 68 | `~s3 & ~s2 & s1 & s0 & t3 & ~t2 & t1 & ~t0` | (L, Q) | MUTATION |
| 69 | `s3 & s2 & s1 & ~s0 & t3 & ~t1 & t0` | (R, F) | MUTATION |
| 70 | `s3 & ~s2 & s1 & s0 & ~t3 & t2 & t1 & t0` | (W, Y) | MUTATION |
| 71 | `~s3 & s2 & ~s1 & s0 & ~t3 & t2 & ~t1 & t0` | (N, D) | MUTATION |
| 72 | `~s3 & ~s2 & ~s1 & ~s0 & ~t3 & ~t2 & t1 & t0` | (A, M) | MUTATION |
| 73 | `~s3 & ~s2 & ~s1 & s0 & ~t3 & ~t2 & t1 & ~t0` | (I, Q) | MUTATION |
| 74 | `~s3 & ~s2 & ~s1 & s0 & t3 & t2 & ~t1 & ~t0` | (F, L) | MUTATION |

---

### Position Pair (74, 75) — MI = 8.82

| Rule | Boolean Expression | Amino Acids | Type |
|------|-------------------|-------------|------|
| 75 | `r3 & r2 & ~r1 & r0 & t3 & ~t2 & t1 & t0` | (Y, K) | MUTATION |
| 76 | `r3 & r2 & ~r1 & r0 & t3 & t2 & t1 & ~t0` | (Y, N) | MUTATION |
| 77 | `~s3 & s2 & ~s1 & s0 & ~t3 & t2 & t1 & t0` | (N, S) | MUTATION |
| 78 | `s3 & s2 & s1 & ~s0 & t3 & t2 & t1 & ~t0` | (R, R) | MUTATION |
| 79 | `~s3 & ~s2 & ~s1 & s0 & ~t3 & ~t1 & t0` | (I, L) | MUTATION |
| 80 | `s3 & ~s2 & ~s1 & ~s0 & t3 & t2 & t1 & t0` | (E, H) | MUTATION |

---

### Position Pair (74, 79) — MI = 8.82

| Rule | Boolean Expression | Amino Acids | Type |
|------|-------------------|-------------|------|
| 81 | `r3 & r2 & ~r1 & r0 & t3 & t2 & t1 & t0` | (Y, Y) | MUTATION |
| 82 | `s3 & ~s2 & ~s1 & ~s0 & t3 & t2 & t1 & t0` | (H, Y) | MUTATION |
| 83 | `s3 & s2 & s1 & ~s0 & ~t3 & t2 & t0` | (R, I) | MUTATION |
| 84 | `~s3 & ~s2 & ~s1 & s0 & t3 & t2 & t1 & t0` | (I, Y) | MUTATION |
| 85 | `~s3 & ~s2 & ~s1 & s0 & t3 & ~t2 & t1 & ~t0` | (I, E) | MUTATION |
| 86 | `r3 & r2 & ~r1 & r0 & ~t3 & ~t2 & t1 & t0` | (Y, N) | MUTATION |
| 87 | `s3 & ~s2 & s1 & s0 & t3 & t2 & t1 & t0` | (W, V) | MUTATION |
| 88 | `s3 & ~s2 & ~s1 & ~s0 & ~t3 & t2 & t1 & t0` | (E, F) | MUTATION |

---

### Position Pair (76, 77) — MI = 8.83

| Rule | Boolean Expression | Amino Acids | Type |
|------|-------------------|-------------|------|
| 89 | `~s3 & ~s2 & ~s1 & s0 & t3 & ~t2 & t0` | (V, Q) | MUTATION |
| 90 | `r3 & r2 & ~r1 & r0 & t3 & ~t2 & t1 & t0` | (Y, K) | MUTATION |
| 91 | `s3 & s2 & s1 & ~s0 & t3 & ~t1 & t0` | (R, F) | MUTATION |
| 92 | `s3 & s2 & s1 & ~s0 & t3 & t2 & t1 & ~t0` | (R, R) | MUTATION |
| 93 | `~s3 & ~s2 & ~s1 & s0 & ~t3 & t2 & t1 & t0` | (I, L) | MUTATION |
| 94 | `r3 & r2 & ~r1 & r0 & t3 & t2 & t1 & ~t0` | (Y, N) | MUTATION |
| 95 | `s3 & ~s2 & ~s1 & ~s0 & ~t3 & t2 & t1 & t0` | (E, W) | MUTATION |

---

### Position Pair (77, 79) — MI = 8.80

| Rule | Boolean Expression | Amino Acids | Type |
|------|-------------------|-------------|------|
| 96 | `~s3 & ~s2 & s1 & s0 & t3 & ~t2 & t1 & ~t0` | (L, K) | MUTATION |
| 97 | `~s3 & s2 & ~s1 & ~s0 & t3 & t2 & t1 & ~t0` | (M, R) | MUTATION |
| 98 | `~s3 & s2 & ~s1 & s0 & ~t3 & t2 & ~t1 & t0` | (N, W) | MUTATION |
| 99 | `s3 & ~s2 & ~s1 & ~s0 & t3 & t2 & t1 & t0` | (H, Y) | MUTATION |
| 100 | `s3 & ~s2 & ~s1 & s0 & ~t3 & t2 & t1 & ~t0` | (K, H) | MUTATION |
| 101 | `~s3 & ~s2 & s1 & s0 & t3 & ~t2 & t0` | (L, I) | MUTATION |
| 102 | `r3 & r2 & ~r1 & r0 & ~t3 & t2 & t1 & t0` | (Y, Q) | MUTATION |
| 103 | `s3 & ~s2 & s1 & s0 & t3 & t2 & t1 & t0` | (W, V) | MUTATION |

---

### Position Pair (78, 79) — MI = 8.82

| Rule | Boolean Expression | Amino Acids | Type |
|------|-------------------|-------------|------|
| 104 | `~s3 & ~s2 & ~s1 & s0 & t3 & ~t2 & t0` | (V, Q) | MUTATION |
| 105 | `r3 & r2 & ~r1 & r0 & t3 & ~t2 & t1 & t0` | (Y, K) | MUTATION |
| 106 | `r3 & r2 & ~r1 & r0 & t3 & t2 & ~t1 & ~t0` | (Y, S) | MUTATION |
| 107 | `s3 & s2 & s1 & ~s0 & t3 & t2 & t1 & ~t0` | (R, R) | MUTATION |
| 108 | `r3 & r2 & ~r1 & r0 & ~t3 & t2 & t1 & t0` | (Y, E) | MUTATION |

---

## Position Pair Summary

| Position Pair | MI | Rules | Reference (i→j) |
|---------------|-----|-------|------------------|
| (68, 75) | 8.75 | 8 | G→F |
| (68, 79) | 8.75 | 9 | G→V |
| (69, 70) | 8.77 | 8 | T→N |
| (69, 79) | 8.78 | 8 | T→V |
| (71, 75) | 8.82 | 7 | G→F |
| (71, 79) | 8.82 | 7 | G→V |
| (72, 75) | 8.82 | 7 | T→F |
| (72, 79) | 8.82 | 6 | T→V |
| (73, 75) | 8.81 | 7 | K→F |
| (73, 79) | 8.82 | 7 | K→V |
| (74, 75) | 8.82 | 6 | R→F |
| (74, 79) | 8.82 | 8 | R→V |
| (76, 77) | 8.83 | 7 | D→N |
| (77, 79) | 8.80 | 8 | N→V |
| (78, 79) | 8.82 | 5 | P→V |

---

## How to Use the Boolean Functions

### For a new sequence:
1. Extract residues at positions 68-79
2. Compare to the reference (G, T, N, T, T, G, T, K, R, D, N, P, V)
3. For each position pair, check if the (aa_i, aa_j) combination matches any rule
4. If ANY rule matches → that position pair is co-evolutionary

### Example:
```
Sequence: ... G T N T G T K R D N P V ...
                        ↑ ↑ ↑ ↑ ↑ ↑
Reference:              G T N T G T K R D N P V

If position 74 mutates R→Y, check rules 75, 81-88.
Rule 81: IF pos 74=Y AND pos 79=Y THEN co-evolutionary
→ Position 79 must also mutate V→Y to satisfy co-evolution.
```

---

## Biological Context

These co-evolutionary rules capture **compensatory mutations** in the SARS-CoV-2 Spike protein. Positions 68-79 are in the **N-terminal signal peptide region** (cleavage site), where mutations are tightly constrained because the protein must maintain its structure and function across variants.

The high mutual information (MI > 8.7) indicates that mutations at these positions are **strongly coupled** — when one position mutates, the other must follow with a specific compensating mutation to maintain protein stability.
