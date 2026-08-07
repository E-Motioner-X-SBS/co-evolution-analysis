# Co-evolution Boolean Functions for SARS-CoV-2 Spike Protein

**Dataset:** 1299 sequences
**Encoding:** Base-20 (He 2012 order)
**Method:** Variable-position K-map with don't-care conditions
**Quine-McCluskey minimization**

---

## Variables

| Variable | Meaning |
|----------|---------|
| s4, s3, s2, s1, s0 | Binary code for residue at position i (0-19, 5 bits) |
| t4, t3, t2, t1, t0 | Binary code for residue at position j (0-19, 5 bits) |
| ~s4 | NOT s4 (bit is 0) |
| s4 | bit is 1 |
| s4.s3.s2.s1.s0 | AND of bits |

---

## Position Pair (495, 498)

| Property | Value |
|----------|-------|
| Mutual Information | 0.8710 |
| Total mutations | 425 |
| Reference pos 495 | R |
| Reference pos 498 | G |
| On-set cells | 2 |
| Off-set cells | 397 |
| Don't-care cells | 625 |
| Prime implicants | 2 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (495, 498): Reference = (R, G)

Co-evolutionary residue pairs (on-set):
  Q-G, R-S

Don't-care positions (conserved): 625 cells
Never-seen pairs (off-set): 397 cells
```

### Boolean Function

```
f(pos_495, pos_498) = 1 if ANY of these residue pairs appear:

    PI_1: pos_495=Q AND pos_498=G
    PI_2: pos_495=R AND pos_498=S

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| Q-G | 0.2641 | 0.0956 | co-evolutionary |
| R-S | 0.1005 | 0.2321 | co-evolutionary |
| R-G | -0.0325 | 0.6723 | anti-correlated |

### Inference Rules (Natural Language)

Rule 1: IF position 495 = Q AND position 498 = G THEN co-evolutionary (MI = 0.871)
Rule 2: IF position 495 = R AND position 498 = S THEN co-evolutionary (MI = 0.871)

---

## Position Pair (448, 454)

| Property | Value |
|----------|-------|
| Mutual Information | 0.8344 |
| Total mutations | 393 |
| Reference pos 448 | G |
| Reference pos 454 | L |
| On-set cells | 3 |
| Off-set cells | 396 |
| Don't-care cells | 625 |
| Prime implicants | 3 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (448, 454): Reference = (G, L)

Co-evolutionary residue pairs (on-set):
  S-L, S-R, G-R

Don't-care positions (conserved): 625 cells
Never-seen pairs (off-set): 396 cells
```

### Boolean Function

```
f(pos_448, pos_454) = 1 if ANY of these residue pairs appear:

    PI_1: pos_448=G AND pos_454=R
    PI_2: pos_448=S AND pos_454=L
    PI_3: pos_448=S AND pos_454=R

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| G-R | 0.2373 | 0.0849 | co-evolutionary |
| S-L | 0.0860 | 0.2176 | co-evolutionary |
| G-L | -0.0254 | 0.6968 | anti-correlated |

### Inference Rules (Natural Language)

Rule 3: IF position 448 = G AND position 454 = R THEN co-evolutionary (MI = 0.834)
Rule 4: IF position 448 = S AND position 454 = L THEN co-evolutionary (MI = 0.834)
Rule 5: IF position 448 = S AND position 454 = R THEN co-evolutionary (MI = 0.834)

---

## Position Pair (488, 498)

| Property | Value |
|----------|-------|
| Mutual Information | 0.8219 |
| Total mutations | 405 |
| Reference pos 488 | F |
| Reference pos 498 | G |
| On-set cells | 3 |
| Off-set cells | 396 |
| Don't-care cells | 625 |
| Prime implicants | 3 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (488, 498): Reference = (F, G)

Co-evolutionary residue pairs (on-set):
  V-G, F-S, P-G

Don't-care positions (conserved): 625 cells
Never-seen pairs (off-set): 396 cells
```

### Boolean Function

```
f(pos_488, pos_498) = 1 if ANY of these residue pairs appear:

    PI_1: pos_488=V AND pos_498=G
    PI_2: pos_488=P AND pos_498=G
    PI_3: pos_488=F AND pos_498=S

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| V-G | 0.2638 | 0.0786 | co-evolutionary |
| F-S | 0.0835 | 0.2319 | co-evolutionary |
| F-G | -0.0266 | 0.6880 | anti-correlated |

### Inference Rules (Natural Language)

Rule 6: IF position 488 = V AND position 498 = G THEN co-evolutionary (MI = 0.822)
Rule 7: IF position 488 = P AND position 498 = G THEN co-evolutionary (MI = 0.822)
Rule 8: IF position 488 = F AND position 498 = S THEN co-evolutionary (MI = 0.822)

---

## Position Pair (442, 454)

| Property | Value |
|----------|-------|
| Mutual Information | 0.8110 |
| Total mutations | 188 |
| Reference pos 442 | K |
| Reference pos 454 | L |
| On-set cells | 3 |
| Off-set cells | 396 |
| Don't-care cells | 625 |
| Prime implicants | 3 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (442, 454): Reference = (K, L)

Co-evolutionary residue pairs (on-set):
  N-L, N-R, K-R

Don't-care positions (conserved): 625 cells
Never-seen pairs (off-set): 396 cells
```

### Boolean Function

```
f(pos_442, pos_454) = 1 if ANY of these residue pairs appear:

    PI_1: pos_442=N AND pos_454=L
    PI_2: pos_442=N AND pos_454=R
    PI_3: pos_442=K AND pos_454=R

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| N-L | 0.0158 | 0.0602 | co-evolutionary |
| K-R | 0.0115 | 0.0821 | co-evolutionary |
| K-L | -0.0011 | 0.8530 | anti-correlated |

### Inference Rules (Natural Language)

Rule 9: IF position 442 = N AND position 454 = L THEN co-evolutionary (MI = 0.811)
Rule 10: IF position 442 = N AND position 454 = R THEN co-evolutionary (MI = 0.811)
Rule 11: IF position 442 = K AND position 454 = R THEN co-evolutionary (MI = 0.811)

---

## Position Pair (442, 448)

| Property | Value |
|----------|-------|
| Mutual Information | 0.7284 |
| Total mutations | 363 |
| Reference pos 442 | K |
| Reference pos 448 | G |
| On-set cells | 3 |
| Off-set cells | 396 |
| Don't-care cells | 625 |
| Prime implicants | 3 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (442, 448): Reference = (K, G)

Co-evolutionary residue pairs (on-set):
  N-S, N-G, K-S

Don't-care positions (conserved): 625 cells
Never-seen pairs (off-set): 396 cells
```

### Boolean Function

```
f(pos_442, pos_448) = 1 if ANY of these residue pairs appear:

    PI_1: pos_442=N AND pos_448=G
    PI_2: pos_442=N AND pos_448=S
    PI_3: pos_442=K AND pos_448=S

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| N-G | 0.2243 | 0.0632 | co-evolutionary |
| K-S | 0.0599 | 0.2186 | co-evolutionary |
| K-G | -0.0176 | 0.7166 | anti-correlated |

### Inference Rules (Natural Language)

Rule 12: IF position 442 = N AND position 448 = G THEN co-evolutionary (MI = 0.728)
Rule 13: IF position 442 = N AND position 448 = S THEN co-evolutionary (MI = 0.728)
Rule 14: IF position 442 = K AND position 448 = S THEN co-evolutionary (MI = 0.728)

---

## Position Pair (212, 215)

| Property | Value |
|----------|-------|
| Mutual Information | 0.3977 |
| Total mutations | 300 |
| Reference pos 212 | V |
| Reference pos 215 | G |
| On-set cells | 4 |
| Off-set cells | 396 |
| Don't-care cells | 624 |
| Prime implicants | 4 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (212, 215): Reference = (V, G)

Co-evolutionary residue pairs (on-set):
  I-V, L-V, V-P, S-G

Don't-care positions (conserved): 624 cells
Never-seen pairs (off-set): 396 cells
```

### Boolean Function

```
f(pos_212, pos_215) = 1 if ANY of these residue pairs appear:

    PI_1: pos_212=V AND pos_215=P
    PI_2: pos_212=S AND pos_215=G
    PI_3: pos_212=I AND pos_215=V
    PI_4: pos_212=L AND pos_215=V

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| S-G | 5.0106 | 0.0067 | co-evolutionary |
| I-V | 2.7593 | 0.0433 | co-evolutionary |
| L-V | 2.7593 | 0.0200 | co-evolutionary |
| V-P | 0.0726 | 0.9300 | co-evolutionary |

### Inference Rules (Natural Language)

Rule 15: IF position 212 = V AND position 215 = P THEN co-evolutionary (MI = 0.398)
Rule 16: IF position 212 = S AND position 215 = G THEN co-evolutionary (MI = 0.398)
Rule 17: IF position 212 = I AND position 215 = V THEN co-evolutionary (MI = 0.398)
Rule 18: IF position 212 = L AND position 215 = V THEN co-evolutionary (MI = 0.398)

---

## Position Pair (215, 216)

| Property | Value |
|----------|-------|
| Mutual Information | 0.3773 |
| Total mutations | 301 |
| Reference pos 215 | G |
| Reference pos 216 | R |
| On-set cells | 4 |
| Off-set cells | 395 |
| Don't-care cells | 625 |
| Prime implicants | 4 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (215, 216): Reference = (G, R)

Co-evolutionary residue pairs (on-set):
  V-R, E-R, P-E, P-K

Don't-care positions (conserved): 625 cells
Never-seen pairs (off-set): 395 cells
```

### Boolean Function

```
f(pos_215, pos_216) = 1 if ANY of these residue pairs appear:

    PI_1: pos_215=E AND pos_216=R
    PI_2: pos_215=V AND pos_216=R
    PI_3: pos_215=P AND pos_216=K
    PI_4: pos_215=P AND pos_216=E

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| P-E | 1.5218 | 0.2175 | co-evolutionary |
| V-R | 0.2463 | 0.0164 | co-evolutionary |
| G-R | 0.2463 | 0.7645 | co-evolutionary |

### Inference Rules (Natural Language)

Rule 19: IF position 215 = E AND position 216 = R THEN co-evolutionary (MI = 0.377)
Rule 20: IF position 215 = V AND position 216 = R THEN co-evolutionary (MI = 0.377)
Rule 21: IF position 215 = P AND position 216 = K THEN co-evolutionary (MI = 0.377)
Rule 22: IF position 215 = P AND position 216 = E THEN co-evolutionary (MI = 0.377)

---

## Position Pair (212, 216)

| Property | Value |
|----------|-------|
| Mutual Information | 0.3773 |
| Total mutations | 301 |
| Reference pos 212 | V |
| Reference pos 216 | R |
| On-set cells | 6 |
| Off-set cells | 394 |
| Don't-care cells | 624 |
| Prime implicants | 6 |
| Essential PIs | 1 |

### K-map (Compact View)

```
Position pair (212, 216): Reference = (V, R)

Co-evolutionary residue pairs (on-set):
  I-R, L-R, V-E, V-K, S-R, G-R

Don't-care positions (conserved): 624 cells
Never-seen pairs (off-set): 394 cells
```

### Boolean Function

```
f(pos_212, pos_216) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_212=G AND pos_216=R
    PI_2: pos_212=V AND pos_216=E
    PI_3: pos_212=I AND pos_216=R
    PI_4: pos_212=L AND pos_216=R
    PI_5: pos_212=V AND pos_216=K
    PI_6: pos_212=S AND pos_216=R

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| I-R | 2.6161 | 0.0432 | co-evolutionary |
| L-R | 2.6161 | 0.0199 | co-evolutionary |
| S-R | 2.6161 | 0.0066 | co-evolutionary |
| V-E | 0.0759 | 0.9236 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 23:** IF position 212 = **G** AND position 216 = **R** THEN co-evolutionary (MI = 0.377)
Rule 24: IF position 212 = V AND position 216 = E THEN co-evolutionary (MI = 0.377)
Rule 25: IF position 212 = I AND position 216 = R THEN co-evolutionary (MI = 0.377)
Rule 26: IF position 212 = L AND position 216 = R THEN co-evolutionary (MI = 0.377)
Rule 27: IF position 212 = V AND position 216 = K THEN co-evolutionary (MI = 0.377)
Rule 28: IF position 212 = S AND position 216 = R THEN co-evolutionary (MI = 0.377)

---

## Position Pair (210, 215)

| Property | Value |
|----------|-------|
| Mutual Information | 0.2377 |
| Total mutations | 282 |
| Reference pos 210 | N |
| Reference pos 215 | G |
| On-set cells | 4 |
| Off-set cells | 395 |
| Don't-care cells | 625 |
| Prime implicants | 4 |
| Essential PIs | 1 |

### K-map (Compact View)

```
Position pair (210, 215): Reference = (N, G)

Co-evolutionary residue pairs (on-set):
  I-P, N-V, N-E, K-G

Don't-care positions (conserved): 625 cells
Never-seen pairs (off-set): 395 cells
```

### Boolean Function

```
f(pos_210, pos_215) = 1 if ANY of these residue pairs appear:

    PI_1: pos_210=I AND pos_215=P
    PI_2: pos_210=N AND pos_215=E
    PI_3: pos_210=N AND pos_215=V
  * PI_4: pos_210=K AND pos_215=G

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| I-P | 1.5315 | 0.2162 | co-evolutionary |
| N-V | 0.2446 | 0.0064 | co-evolutionary |
| N-G | 0.2436 | 0.7758 | co-evolutionary |

### Inference Rules (Natural Language)

Rule 29: IF position 210 = I AND position 215 = P THEN co-evolutionary (MI = 0.238)
Rule 30: IF position 210 = N AND position 215 = E THEN co-evolutionary (MI = 0.238)
Rule 31: IF position 210 = N AND position 215 = V THEN co-evolutionary (MI = 0.238)
**Rule 32:** IF position 210 = **K** AND position 215 = **G** THEN co-evolutionary (MI = 0.238)

---

## Position Pair (210, 212)

| Property | Value |
|----------|-------|
| Mutual Information | 0.1769 |
| Total mutations | 301 |
| Reference pos 210 | N |
| Reference pos 212 | V |
| On-set cells | 4 |
| Off-set cells | 396 |
| Don't-care cells | 624 |
| Prime implicants | 4 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (210, 212): Reference = (N, V)

Co-evolutionary residue pairs (on-set):
  I-V, I-G, N-L, N-S

Don't-care positions (conserved): 624 cells
Never-seen pairs (off-set): 396 cells
```

### Boolean Function

```
f(pos_210, pos_212) = 1 if ANY of these residue pairs appear:

    PI_1: pos_210=I AND pos_212=V
    PI_2: pos_210=I AND pos_212=G
    PI_3: pos_210=N AND pos_212=L
    PI_4: pos_210=N AND pos_212=S

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| N-L | 3.6277 | 0.0199 | co-evolutionary |
| N-S | 3.6277 | 0.0066 | co-evolutionary |
| I-V | 0.0269 | 0.9701 | co-evolutionary |

### Inference Rules (Natural Language)

Rule 33: IF position 210 = I AND position 212 = V THEN co-evolutionary (MI = 0.177)
Rule 34: IF position 210 = I AND position 212 = G THEN co-evolutionary (MI = 0.177)
Rule 35: IF position 210 = N AND position 212 = L THEN co-evolutionary (MI = 0.177)
Rule 36: IF position 210 = N AND position 212 = S THEN co-evolutionary (MI = 0.177)

---

## Summary

| Metric | Value |
|--------|-------|
| Sequences | 1299 |
| Variable positions | 21 |
| Co-evolutionary pairs | 10 |
| Total inference rules | 36 |
| Position pairs with rules | 10 |

## How to Apply

1. Extract the residue pair at each co-evolving position pair from a new sequence
2. For each position pair, check if the residue pair matches any rule below
3. If YES: that position pair is co-evolutionary (consistent with the observed data)
4. If position i mutates: find which residue at position j satisfies the co-evolutionary constraint

**Example:** If position 212 mutates to G, check rules for position 212.
Essential rule: IF pos 212 = G AND pos 216 = R THEN co-evolutionary.
So position 216 must also show R to satisfy the co-evolutionary constraint.