# Co-evolution Boolean Functions for SARS-CoV-2 Spike Protein

**Dataset:** 1299 sequences
**Encoding:** Base-20 (He 2012 order)
**Method:** Variable-position K-map with don't-care conditions
**Quine-McCluskey minimization**

---

## Variables

| Variable | Meaning |
|----------|---------|
| s3, s2, s1, s0 | Binary code for residue at position i (0-19) |
| t3, t2, t1, t0 | Binary code for residue at position j (0-19) |
| ~s3 | NOT s3 (bit is 0) |
| s3 | bit is 1 |
| s3.s2.s1.s0 | AND of bits |

---

## Position Pair (413, 427)

| Property | Value |
|----------|-------|
| Mutual Information | 2.3524 |
| Total mutations | 271 |
| Reference pos 413 | N |
| Reference pos 427 | G |
| On-set cells | 12 |
| Off-set cells | 387 |
| Don't-care cells | 1 |
| Prime implicants | 11 |
| Essential PIs | 10 |

### K-map (Compact View)

```
Position pair (413, 427): Reference = (N, G)

Co-evolutionary residue pairs (on-set):
  A-V, I-C, Y-A, D-I, Q-D, N-I, N-W, K-S, K-G, T-F, P-P, G-T

Don't-care positions (conserved): 1 cells
Never-seen pairs (off-set): 387 cells
```

### Boolean Function

```
f(pos_413, pos_427) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_413=A AND pos_427=V
  * PI_2: pos_413=W AND pos_427=E
  * PI_3: pos_413=I AND pos_427=V
  * PI_4: pos_413=L AND pos_427=F
  * PI_5: pos_413=N AND pos_427=F
  * PI_6: pos_413=K AND pos_427=I
  * PI_7: pos_413=K AND pos_427=K
  * PI_8: pos_413=R AND pos_427=V
    PI_9: pos_413=R AND pos_427=S
  * PI_10: pos_413=M AND pos_427=F
  * PI_11: pos_413=E AND pos_427=H

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| Q-D | 5.2234 | 0.0054 | co-evolutionary |
| Y-A | 4.3361 | 0.0131 | co-evolutionary |
| D-I | 3.4558 | 0.0308 | co-evolutionary |
| I-C | 3.3407 | 0.0354 | co-evolutionary |
| A-V | 2.4332 | 0.0878 | co-evolutionary |
| N-G | 0.1957 | 0.7914 | co-evolutionary |
| K-G | 0.1720 | 0.0285 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 1:** IF position 413 = **A** AND position 427 = **V** THEN co-evolutionary (MI = 2.352)
**Rule 2:** IF position 413 = **W** AND position 427 = **E** THEN co-evolutionary (MI = 2.352)
**Rule 3:** IF position 413 = **I** AND position 427 = **V** THEN co-evolutionary (MI = 2.352)
**Rule 4:** IF position 413 = **L** AND position 427 = **F** THEN co-evolutionary (MI = 2.352)
**Rule 5:** IF position 413 = **N** AND position 427 = **F** THEN co-evolutionary (MI = 2.352)
**Rule 6:** IF position 413 = **K** AND position 427 = **I** THEN co-evolutionary (MI = 2.352)
**Rule 7:** IF position 413 = **K** AND position 427 = **K** THEN co-evolutionary (MI = 2.352)
**Rule 8:** IF position 413 = **R** AND position 427 = **V** THEN co-evolutionary (MI = 2.352)
Rule 9: IF position 413 = R AND position 427 = S THEN co-evolutionary (MI = 2.352)
**Rule 10:** IF position 413 = **M** AND position 427 = **F** THEN co-evolutionary (MI = 2.352)
**Rule 11:** IF position 413 = **E** AND position 427 = **H** THEN co-evolutionary (MI = 2.352)

---

## Position Pair (413, 425)

| Property | Value |
|----------|-------|
| Mutual Information | 2.3345 |
| Total mutations | 271 |
| Reference pos 413 | N |
| Reference pos 425 | F |
| On-set cells | 12 |
| Off-set cells | 387 |
| Don't-care cells | 1 |
| Prime implicants | 11 |
| Essential PIs | 10 |

### K-map (Compact View)

```
Position pair (413, 425): Reference = (N, F)

Co-evolutionary residue pairs (on-set):
  A-G, I-T, Y-V, D-C, Q-P, N-I, N-C, K-F, K-W, T-D, P-K, G-D

Don't-care positions (conserved): 1 cells
Never-seen pairs (off-set): 387 cells
```

### Boolean Function

```
f(pos_413, pos_425) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_413=E AND pos_425=F
  * PI_2: pos_413=A AND pos_425=D
  * PI_3: pos_413=A AND pos_425=D
  * PI_4: pos_413=I AND pos_425=V
  * PI_5: pos_413=L AND pos_425=M
  * PI_6: pos_413=W AND pos_425=N
  * PI_7: pos_413=K AND pos_425=Q
  * PI_8: pos_413=K AND pos_425=K
    PI_9: pos_413=R AND pos_425=I
  * PI_10: pos_413=R AND pos_425=K
  * PI_11: pos_413=W AND pos_425=F

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| Q-P | 5.2234 | 0.0054 | co-evolutionary |
| Y-V | 4.3361 | 0.0131 | co-evolutionary |
| D-C | 3.4558 | 0.0308 | co-evolutionary |
| I-T | 3.3407 | 0.0354 | co-evolutionary |
| A-G | 2.4332 | 0.0878 | co-evolutionary |
| N-F | 0.1957 | 0.7914 | co-evolutionary |
| K-F | 0.1720 | 0.0285 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 12:** IF position 413 = **E** AND position 425 = **F** THEN co-evolutionary (MI = 2.334)
**Rule 13:** IF position 413 = **A** AND position 425 = **D** THEN co-evolutionary (MI = 2.334)
**Rule 14:** IF position 413 = **A** AND position 425 = **D** THEN co-evolutionary (MI = 2.334)
**Rule 15:** IF position 413 = **I** AND position 425 = **V** THEN co-evolutionary (MI = 2.334)
**Rule 16:** IF position 413 = **L** AND position 425 = **M** THEN co-evolutionary (MI = 2.334)
**Rule 17:** IF position 413 = **W** AND position 425 = **N** THEN co-evolutionary (MI = 2.334)
**Rule 18:** IF position 413 = **K** AND position 425 = **Q** THEN co-evolutionary (MI = 2.334)
**Rule 19:** IF position 413 = **K** AND position 425 = **K** THEN co-evolutionary (MI = 2.334)
Rule 20: IF position 413 = R AND position 425 = I THEN co-evolutionary (MI = 2.334)
**Rule 21:** IF position 413 = **R** AND position 425 = **K** THEN co-evolutionary (MI = 2.334)
**Rule 22:** IF position 413 = **W** AND position 425 = **F** THEN co-evolutionary (MI = 2.334)

---

## Position Pair (413, 426)

| Property | Value |
|----------|-------|
| Mutual Information | 2.3199 |
| Total mutations | 271 |
| Reference pos 413 | N |
| Reference pos 426 | T |
| On-set cells | 12 |
| Off-set cells | 387 |
| Don't-care cells | 1 |
| Prime implicants | 13 |
| Essential PIs | 12 |

### K-map (Compact View)

```
Position pair (413, 426): Reference = (N, T)

Co-evolutionary residue pairs (on-set):
  A-C, I-G, Y-I, D-V, Q-D, N-A, N-V, K-N, K-T, T-D, P-L, G-F

Don't-care positions (conserved): 1 cells
Never-seen pairs (off-set): 387 cells
```

### Boolean Function

```
f(pos_413, pos_426) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_413=I AND pos_426=I
  * PI_2: pos_413=L AND pos_426=W
  * PI_3: pos_413=W AND pos_426=D
  * PI_4: pos_413=N AND pos_426=W
  * PI_5: pos_413=K AND pos_426=I
  * PI_6: pos_413=K AND pos_426=H
  * PI_7: pos_413=K AND pos_426=S
    PI_8: pos_413=R AND pos_426=H
  * PI_9: pos_413=A AND pos_426=S
  * PI_10: pos_413=I AND pos_426=M
  * PI_11: pos_413=M AND pos_426=D
  * PI_12: pos_413=Y AND pos_426=Q
  * PI_13: pos_413=E AND pos_426=I

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| Q-D | 4.8668 | 0.0054 | co-evolutionary |
| Y-I | 4.3361 | 0.0131 | co-evolutionary |
| D-V | 3.4558 | 0.0308 | co-evolutionary |
| I-G | 3.3407 | 0.0354 | co-evolutionary |
| A-C | 2.4332 | 0.0878 | co-evolutionary |
| N-T | 0.1957 | 0.7914 | co-evolutionary |
| K-T | 0.1720 | 0.0285 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 23:** IF position 413 = **I** AND position 426 = **I** THEN co-evolutionary (MI = 2.320)
**Rule 24:** IF position 413 = **L** AND position 426 = **W** THEN co-evolutionary (MI = 2.320)
**Rule 25:** IF position 413 = **W** AND position 426 = **D** THEN co-evolutionary (MI = 2.320)
**Rule 26:** IF position 413 = **N** AND position 426 = **W** THEN co-evolutionary (MI = 2.320)
**Rule 27:** IF position 413 = **K** AND position 426 = **I** THEN co-evolutionary (MI = 2.320)
**Rule 28:** IF position 413 = **K** AND position 426 = **H** THEN co-evolutionary (MI = 2.320)
**Rule 29:** IF position 413 = **K** AND position 426 = **S** THEN co-evolutionary (MI = 2.320)
Rule 30: IF position 413 = R AND position 426 = H THEN co-evolutionary (MI = 2.320)
**Rule 31:** IF position 413 = **A** AND position 426 = **S** THEN co-evolutionary (MI = 2.320)
**Rule 32:** IF position 413 = **I** AND position 426 = **M** THEN co-evolutionary (MI = 2.320)
**Rule 33:** IF position 413 = **M** AND position 426 = **D** THEN co-evolutionary (MI = 2.320)
**Rule 34:** IF position 413 = **Y** AND position 426 = **Q** THEN co-evolutionary (MI = 2.320)
**Rule 35:** IF position 413 = **E** AND position 426 = **I** THEN co-evolutionary (MI = 2.320)

---

## Position Pair (413, 424)

| Property | Value |
|----------|-------|
| Mutual Information | 2.3104 |
| Total mutations | 271 |
| Reference pos 413 | N |
| Reference pos 424 | D |
| On-set cells | 12 |
| Off-set cells | 387 |
| Don't-care cells | 1 |
| Prime implicants | 13 |
| Essential PIs | 12 |

### K-map (Compact View)

```
Position pair (413, 424): Reference = (N, D)

Co-evolutionary residue pairs (on-set):
  A-T, I-F, Y-C, D-G, Q-L, N-V, N-G, K-A, K-D, T-P, P-Y, G-D

Don't-care positions (conserved): 1 cells
Never-seen pairs (off-set): 387 cells
```

### Boolean Function

```
f(pos_413, pos_424) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_413=I AND pos_424=A
  * PI_2: pos_413=I AND pos_424=D
  * PI_3: pos_413=E AND pos_424=D
  * PI_4: pos_413=H AND pos_424=W
  * PI_5: pos_413=H AND pos_424=Q
  * PI_6: pos_413=K AND pos_424=S
    PI_7: pos_413=R AND pos_424=F
  * PI_8: pos_413=R AND pos_424=S
  * PI_9: pos_413=A AND pos_424=M
  * PI_10: pos_413=A AND pos_424=K
  * PI_11: pos_413=F AND pos_424=L
  * PI_12: pos_413=Y AND pos_424=R
  * PI_13: pos_413=E AND pos_424=F

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| Q-L | 5.2234 | 0.0054 | co-evolutionary |
| Y-C | 4.3361 | 0.0131 | co-evolutionary |
| D-G | 3.4558 | 0.0308 | co-evolutionary |
| I-F | 3.3407 | 0.0354 | co-evolutionary |
| A-T | 2.4332 | 0.0878 | co-evolutionary |
| N-D | 0.1938 | 0.7914 | co-evolutionary |
| K-D | 0.1701 | 0.0285 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 36:** IF position 413 = **I** AND position 424 = **A** THEN co-evolutionary (MI = 2.310)
**Rule 37:** IF position 413 = **I** AND position 424 = **D** THEN co-evolutionary (MI = 2.310)
**Rule 38:** IF position 413 = **E** AND position 424 = **D** THEN co-evolutionary (MI = 2.310)
**Rule 39:** IF position 413 = **H** AND position 424 = **W** THEN co-evolutionary (MI = 2.310)
**Rule 40:** IF position 413 = **H** AND position 424 = **Q** THEN co-evolutionary (MI = 2.310)
**Rule 41:** IF position 413 = **K** AND position 424 = **S** THEN co-evolutionary (MI = 2.310)
Rule 42: IF position 413 = R AND position 424 = F THEN co-evolutionary (MI = 2.310)
**Rule 43:** IF position 413 = **R** AND position 424 = **S** THEN co-evolutionary (MI = 2.310)
**Rule 44:** IF position 413 = **A** AND position 424 = **M** THEN co-evolutionary (MI = 2.310)
**Rule 45:** IF position 413 = **A** AND position 424 = **K** THEN co-evolutionary (MI = 2.310)
**Rule 46:** IF position 413 = **F** AND position 424 = **L** THEN co-evolutionary (MI = 2.310)
**Rule 47:** IF position 413 = **Y** AND position 424 = **R** THEN co-evolutionary (MI = 2.310)
**Rule 48:** IF position 413 = **E** AND position 424 = **F** THEN co-evolutionary (MI = 2.310)

---

## Position Pair (1026, 1040)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2989 |
| Total mutations | 278 |
| Reference pos 1026 | S |
| Reference pos 1040 | G |
| On-set cells | 11 |
| Off-set cells | 388 |
| Don't-care cells | 1 |
| Prime implicants | 10 |
| Essential PIs | 10 |

### K-map (Compact View)

```
Position pair (1026, 1040): Reference = (S, G)

Co-evolutionary residue pairs (on-set):
  A-R, L-H, V-Y, M-C, E-K, Q-M, K-F, S-S, T-D, C-G, G-L

Don't-care positions (conserved): 1 cells
Never-seen pairs (off-set): 388 cells
```

### Boolean Function

```
f(pos_1026, pos_1040) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_1026=A AND pos_1040=D
  * PI_2: pos_1026=V AND pos_1040=N
  * PI_3: pos_1026=A AND pos_1040=R
  * PI_4: pos_1026=V AND pos_1040=M
  * PI_5: pos_1026=M AND pos_1040=L
  * PI_6: pos_1026=Y AND pos_1040=I
  * PI_7: pos_1026=Q AND pos_1040=K
  * PI_8: pos_1026=H AND pos_1040=H
  * PI_9: pos_1026=Y AND pos_1040=W
  * PI_10: pos_1026=W AND pos_1040=R

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| T-D | 5.2234 | 0.0054 | co-evolutionary |
| G-L | 4.7715 | 0.0085 | co-evolutionary |
| L-H | 3.8735 | 0.0208 | co-evolutionary |
| V-Y | 3.4805 | 0.0308 | co-evolutionary |
| E-K | 2.9647 | 0.0516 | co-evolutionary |
| C-G | 0.1323 | 0.0901 | co-evolutionary |
| S-G | 0.1313 | 0.7860 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 49:** IF position 1026 = **A** AND position 1040 = **D** THEN co-evolutionary (MI = 2.299)
**Rule 50:** IF position 1026 = **V** AND position 1040 = **N** THEN co-evolutionary (MI = 2.299)
**Rule 51:** IF position 1026 = **A** AND position 1040 = **R** THEN co-evolutionary (MI = 2.299)
**Rule 52:** IF position 1026 = **V** AND position 1040 = **M** THEN co-evolutionary (MI = 2.299)
**Rule 53:** IF position 1026 = **M** AND position 1040 = **L** THEN co-evolutionary (MI = 2.299)
**Rule 54:** IF position 1026 = **Y** AND position 1040 = **I** THEN co-evolutionary (MI = 2.299)
**Rule 55:** IF position 1026 = **Q** AND position 1040 = **K** THEN co-evolutionary (MI = 2.299)
**Rule 56:** IF position 1026 = **H** AND position 1040 = **H** THEN co-evolutionary (MI = 2.299)
**Rule 57:** IF position 1026 = **Y** AND position 1040 = **W** THEN co-evolutionary (MI = 2.299)
**Rule 58:** IF position 1026 = **W** AND position 1040 = **R** THEN co-evolutionary (MI = 2.299)

---

## Position Pair (1026, 1042)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2989 |
| Total mutations | 278 |
| Reference pos 1026 | S |
| Reference pos 1042 | G |
| On-set cells | 11 |
| Off-set cells | 388 |
| Don't-care cells | 1 |
| Prime implicants | 10 |
| Essential PIs | 10 |

### K-map (Compact View)

```
Position pair (1026, 1042): Reference = (S, G)

Co-evolutionary residue pairs (on-set):
  A-D, L-M, V-L, M-K, E-Y, Q-F, K-G, S-P, T-C, C-H, G-S

Don't-care positions (conserved): 1 cells
Never-seen pairs (off-set): 388 cells
```

### Boolean Function

```
f(pos_1026, pos_1042) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_1026=V AND pos_1042=R
  * PI_2: pos_1026=A AND pos_1042=D
  * PI_3: pos_1026=L AND pos_1042=H
  * PI_4: pos_1026=F AND pos_1042=K
  * PI_5: pos_1026=Q AND pos_1042=Y
  * PI_6: pos_1026=H AND pos_1042=K
  * PI_7: pos_1026=I AND pos_1042=W
  * PI_8: pos_1026=F AND pos_1042=I
  * PI_9: pos_1026=Y AND pos_1042=A
  * PI_10: pos_1026=E AND pos_1042=N

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| T-C | 5.2234 | 0.0054 | co-evolutionary |
| G-S | 4.7715 | 0.0085 | co-evolutionary |
| L-M | 3.8735 | 0.0208 | co-evolutionary |
| V-L | 3.4805 | 0.0308 | co-evolutionary |
| E-Y | 2.9647 | 0.0516 | co-evolutionary |
| C-H | 2.4072 | 0.0901 | co-evolutionary |
| S-G | 0.2379 | 0.7860 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 59:** IF position 1026 = **V** AND position 1042 = **R** THEN co-evolutionary (MI = 2.299)
**Rule 60:** IF position 1026 = **A** AND position 1042 = **D** THEN co-evolutionary (MI = 2.299)
**Rule 61:** IF position 1026 = **L** AND position 1042 = **H** THEN co-evolutionary (MI = 2.299)
**Rule 62:** IF position 1026 = **F** AND position 1042 = **K** THEN co-evolutionary (MI = 2.299)
**Rule 63:** IF position 1026 = **Q** AND position 1042 = **Y** THEN co-evolutionary (MI = 2.299)
**Rule 64:** IF position 1026 = **H** AND position 1042 = **K** THEN co-evolutionary (MI = 2.299)
**Rule 65:** IF position 1026 = **I** AND position 1042 = **W** THEN co-evolutionary (MI = 2.299)
**Rule 66:** IF position 1026 = **F** AND position 1042 = **I** THEN co-evolutionary (MI = 2.299)
**Rule 67:** IF position 1026 = **Y** AND position 1042 = **A** THEN co-evolutionary (MI = 2.299)
**Rule 68:** IF position 1026 = **E** AND position 1042 = **N** THEN co-evolutionary (MI = 2.299)

---

## Position Pair (1040, 1042)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2989 |
| Total mutations | 278 |
| Reference pos 1040 | G |
| Reference pos 1042 | G |
| On-set cells | 11 |
| Off-set cells | 388 |
| Don't-care cells | 1 |
| Prime implicants | 10 |
| Essential PIs | 9 |

### K-map (Compact View)

```
Position pair (1040, 1042): Reference = (G, G)

Co-evolutionary residue pairs (on-set):
  L-S, M-F, F-G, Y-L, D-C, H-M, K-Y, R-D, S-P, C-K, G-H

Don't-care positions (conserved): 1 cells
Never-seen pairs (off-set): 388 cells
```

### Boolean Function

```
f(pos_1040, pos_1042) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_1040=V AND pos_1042=W
  * PI_2: pos_1040=L AND pos_1042=I
  * PI_3: pos_1040=F AND pos_1042=F
  * PI_4: pos_1040=W AND pos_1042=Q
  * PI_5: pos_1040=H AND pos_1042=F
  * PI_6: pos_1040=S AND pos_1042=M
  * PI_7: pos_1040=A AND pos_1042=Q
  * PI_8: pos_1040=V AND pos_1042=R
  * PI_9: pos_1040=E AND pos_1042=E
    PI_10: pos_1040=E AND pos_1042=S

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| D-C | 5.2234 | 0.0054 | co-evolutionary |
| L-S | 4.7715 | 0.0085 | co-evolutionary |
| H-M | 3.8735 | 0.0208 | co-evolutionary |
| Y-L | 3.4805 | 0.0308 | co-evolutionary |
| K-Y | 2.9647 | 0.0516 | co-evolutionary |
| G-H | 0.1323 | 0.0901 | co-evolutionary |
| G-G | 0.1304 | 0.7860 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 69:** IF position 1040 = **V** AND position 1042 = **W** THEN co-evolutionary (MI = 2.299)
**Rule 70:** IF position 1040 = **L** AND position 1042 = **I** THEN co-evolutionary (MI = 2.299)
**Rule 71:** IF position 1040 = **F** AND position 1042 = **F** THEN co-evolutionary (MI = 2.299)
**Rule 72:** IF position 1040 = **W** AND position 1042 = **Q** THEN co-evolutionary (MI = 2.299)
**Rule 73:** IF position 1040 = **H** AND position 1042 = **F** THEN co-evolutionary (MI = 2.299)
**Rule 74:** IF position 1040 = **S** AND position 1042 = **M** THEN co-evolutionary (MI = 2.299)
**Rule 75:** IF position 1040 = **A** AND position 1042 = **Q** THEN co-evolutionary (MI = 2.299)
**Rule 76:** IF position 1040 = **V** AND position 1042 = **R** THEN co-evolutionary (MI = 2.299)
**Rule 77:** IF position 1040 = **E** AND position 1042 = **E** THEN co-evolutionary (MI = 2.299)
Rule 78: IF position 1040 = E AND position 1042 = S THEN co-evolutionary (MI = 2.299)

---

## Position Pair (462, 473)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2946 |
| Total mutations | 265 |
| Reference pos 462 | R |
| Reference pos 473 | N |
| On-set cells | 10 |
| Off-set cells | 389 |
| Don't-care cells | 1 |
| Prime implicants | 10 |
| Essential PIs | 9 |

### K-map (Compact View)

```
Position pair (462, 473): Reference = (R, N)

Co-evolutionary residue pairs (on-set):
  I-P, L-I, F-A, Y-A, E-G, D-K, R-S, S-C, T-N, P-Q

Don't-care positions (conserved): 1 cells
Never-seen pairs (off-set): 389 cells
```

### Boolean Function

```
f(pos_462, pos_473) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_462=L AND pos_473=Y
    PI_2: pos_462=L AND pos_473=V
  * PI_3: pos_462=L AND pos_473=D
  * PI_4: pos_462=Y AND pos_473=M
  * PI_5: pos_462=W AND pos_473=E
  * PI_6: pos_462=N AND pos_473=V
  * PI_7: pos_462=H AND pos_473=I
  * PI_8: pos_462=V AND pos_473=K
  * PI_9: pos_462=M AND pos_473=N
  * PI_10: pos_462=W AND pos_473=L

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| P-Q | 5.2234 | 0.0054 | co-evolutionary |
| S-C | 3.3852 | 0.0339 | co-evolutionary |
| D-K | 3.0918 | 0.0454 | co-evolutionary |
| I-P | 2.4509 | 0.0862 | co-evolutionary |
| R-S | 0.2176 | 0.0085 | co-evolutionary |
| T-N | 0.2071 | 0.0169 | co-evolutionary |
| R-N | 0.1965 | 0.7960 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 79:** IF position 462 = **L** AND position 473 = **Y** THEN co-evolutionary (MI = 2.295)
Rule 80: IF position 462 = L AND position 473 = V THEN co-evolutionary (MI = 2.295)
**Rule 81:** IF position 462 = **L** AND position 473 = **D** THEN co-evolutionary (MI = 2.295)
**Rule 82:** IF position 462 = **Y** AND position 473 = **M** THEN co-evolutionary (MI = 2.295)
**Rule 83:** IF position 462 = **W** AND position 473 = **E** THEN co-evolutionary (MI = 2.295)
**Rule 84:** IF position 462 = **N** AND position 473 = **V** THEN co-evolutionary (MI = 2.295)
**Rule 85:** IF position 462 = **H** AND position 473 = **I** THEN co-evolutionary (MI = 2.295)
**Rule 86:** IF position 462 = **V** AND position 473 = **K** THEN co-evolutionary (MI = 2.295)
**Rule 87:** IF position 462 = **M** AND position 473 = **N** THEN co-evolutionary (MI = 2.295)
**Rule 88:** IF position 462 = **W** AND position 473 = **L** THEN co-evolutionary (MI = 2.295)

---

## Position Pair (468, 473)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2946 |
| Total mutations | 265 |
| Reference pos 468 | I |
| Reference pos 473 | N |
| On-set cells | 11 |
| Off-set cells | 388 |
| Don't-care cells | 1 |
| Prime implicants | 11 |
| Essential PIs | 10 |

### K-map (Compact View)

```
Position pair (468, 473): Reference = (I, N)

Co-evolutionary residue pairs (on-set):
  A-C, I-S, Y-K, E-G, D-I, Q-P, N-G, S-Q, T-A, P-A, G-N

Don't-care positions (conserved): 1 cells
Never-seen pairs (off-set): 388 cells
```

### Boolean Function

```
f(pos_468, pos_473) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_468=E AND pos_473=F
  * PI_2: pos_468=I AND pos_473=I
    PI_3: pos_468=I AND pos_473=S
  * PI_4: pos_468=L AND pos_473=V
  * PI_5: pos_468=N AND pos_473=V
  * PI_6: pos_468=N AND pos_473=F
  * PI_7: pos_468=K AND pos_473=Q
  * PI_8: pos_468=R AND pos_473=S
  * PI_9: pos_468=V AND pos_473=Y
  * PI_10: pos_468=M AND pos_473=A
  * PI_11: pos_468=Y AND pos_473=E

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| S-Q | 5.2234 | 0.0054 | co-evolutionary |
| A-C | 3.3852 | 0.0339 | co-evolutionary |
| Y-K | 3.0918 | 0.0454 | co-evolutionary |
| Q-P | 2.4509 | 0.0862 | co-evolutionary |
| I-S | 0.2176 | 0.0085 | co-evolutionary |
| G-N | 0.2071 | 0.0169 | co-evolutionary |
| I-N | 0.1965 | 0.7960 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 89:** IF position 468 = **E** AND position 473 = **F** THEN co-evolutionary (MI = 2.295)
**Rule 90:** IF position 468 = **I** AND position 473 = **I** THEN co-evolutionary (MI = 2.295)
Rule 91: IF position 468 = I AND position 473 = S THEN co-evolutionary (MI = 2.295)
**Rule 92:** IF position 468 = **L** AND position 473 = **V** THEN co-evolutionary (MI = 2.295)
**Rule 93:** IF position 468 = **N** AND position 473 = **V** THEN co-evolutionary (MI = 2.295)
**Rule 94:** IF position 468 = **N** AND position 473 = **F** THEN co-evolutionary (MI = 2.295)
**Rule 95:** IF position 468 = **K** AND position 473 = **Q** THEN co-evolutionary (MI = 2.295)
**Rule 96:** IF position 468 = **R** AND position 473 = **S** THEN co-evolutionary (MI = 2.295)
**Rule 97:** IF position 468 = **V** AND position 473 = **Y** THEN co-evolutionary (MI = 2.295)
**Rule 98:** IF position 468 = **M** AND position 473 = **A** THEN co-evolutionary (MI = 2.295)
**Rule 99:** IF position 468 = **Y** AND position 473 = **E** THEN co-evolutionary (MI = 2.295)

---

## Position Pair (1064, 1074)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2888 |
| Total mutations | 278 |
| Reference pos 1064 | V |
| Reference pos 1074 | A |
| On-set cells | 11 |
| Off-set cells | 388 |
| Don't-care cells | 1 |
| Prime implicants | 11 |
| Essential PIs | 10 |

### K-map (Compact View)

```
Position pair (1064, 1074): Reference = (V, A)

Co-evolutionary residue pairs (on-set):
  A-A, L-K, V-F, F-G, Y-T, E-C, Q-I, N-D, K-H, T-T, P-P

Don't-care positions (conserved): 1 cells
Never-seen pairs (off-set): 388 cells
```

### Boolean Function

```
f(pos_1064, pos_1074) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_1064=A AND pos_1074=A
  * PI_2: pos_1064=I AND pos_1074=A
  * PI_3: pos_1064=V AND pos_1074=F
    PI_4: pos_1064=V AND pos_1074=H
  * PI_5: pos_1064=M AND pos_1074=I
  * PI_6: pos_1064=W AND pos_1074=W
  * PI_7: pos_1064=E AND pos_1074=E
  * PI_8: pos_1064=N AND pos_1074=I
  * PI_9: pos_1064=H AND pos_1074=D
  * PI_10: pos_1064=R AND pos_1074=F
  * PI_11: pos_1064=W AND pos_1074=Q

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| K-H | 4.6844 | 0.0092 | co-evolutionary |
| E-C | 3.9113 | 0.0200 | co-evolutionary |
| Q-I | 3.4805 | 0.0308 | co-evolutionary |
| P-P | 2.9647 | 0.0516 | co-evolutionary |
| V-F | 0.2340 | 0.0054 | co-evolutionary |
| A-A | 0.1323 | 0.0901 | co-evolutionary |
| V-A | 0.1255 | 0.7860 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 100:** IF position 1064 = **A** AND position 1074 = **A** THEN co-evolutionary (MI = 2.289)
**Rule 101:** IF position 1064 = **I** AND position 1074 = **A** THEN co-evolutionary (MI = 2.289)
**Rule 102:** IF position 1064 = **V** AND position 1074 = **F** THEN co-evolutionary (MI = 2.289)
Rule 103: IF position 1064 = V AND position 1074 = H THEN co-evolutionary (MI = 2.289)
**Rule 104:** IF position 1064 = **M** AND position 1074 = **I** THEN co-evolutionary (MI = 2.289)
**Rule 105:** IF position 1064 = **W** AND position 1074 = **W** THEN co-evolutionary (MI = 2.289)
**Rule 106:** IF position 1064 = **E** AND position 1074 = **E** THEN co-evolutionary (MI = 2.289)
**Rule 107:** IF position 1064 = **N** AND position 1074 = **I** THEN co-evolutionary (MI = 2.289)
**Rule 108:** IF position 1064 = **H** AND position 1074 = **D** THEN co-evolutionary (MI = 2.289)
**Rule 109:** IF position 1064 = **R** AND position 1074 = **F** THEN co-evolutionary (MI = 2.289)
**Rule 110:** IF position 1064 = **W** AND position 1074 = **Q** THEN co-evolutionary (MI = 2.289)

---

## Position Pair (1064, 1065)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2876 |
| Total mutations | 278 |
| Reference pos 1064 | V |
| Reference pos 1065 | P |
| On-set cells | 12 |
| Off-set cells | 387 |
| Don't-care cells | 1 |
| Prime implicants | 11 |
| Essential PIs | 11 |

### K-map (Compact View)

```
Position pair (1064, 1065): Reference = (V, P)

Co-evolutionary residue pairs (on-set):
  A-Q, L-H, V-T, F-T, Y-V, E-K, Q-E, N-F, K-N, T-Y, P-A, P-S

Don't-care positions (conserved): 1 cells
Never-seen pairs (off-set): 387 cells
```

### Boolean Function

```
f(pos_1064, pos_1065) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_1064=V AND pos_1065=M
  * PI_2: pos_1064=M AND pos_1065=H
  * PI_3: pos_1064=M AND pos_1065=Y
  * PI_4: pos_1064=A AND pos_1065=Q
  * PI_5: pos_1064=W AND pos_1065=N
  * PI_6: pos_1064=Q AND pos_1065=K
  * PI_7: pos_1064=K AND pos_1065=A
  * PI_8: pos_1064=R AND pos_1065=I
  * PI_9: pos_1064=A AND pos_1065=S
  * PI_10: pos_1064=Y AND pos_1065=E
  * PI_11: pos_1064=W AND pos_1065=W

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| K-N | 4.6844 | 0.0092 | co-evolutionary |
| E-K | 3.9113 | 0.0200 | co-evolutionary |
| Q-E | 3.4805 | 0.0308 | co-evolutionary |
| P-A | 2.9647 | 0.0508 | co-evolutionary |
| A-Q | 2.4072 | 0.0901 | co-evolutionary |
| V-P | 0.2340 | 0.7860 | co-evolutionary |
| V-T | 0.1004 | 0.0054 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 111:** IF position 1064 = **V** AND position 1065 = **M** THEN co-evolutionary (MI = 2.288)
**Rule 112:** IF position 1064 = **M** AND position 1065 = **H** THEN co-evolutionary (MI = 2.288)
**Rule 113:** IF position 1064 = **M** AND position 1065 = **Y** THEN co-evolutionary (MI = 2.288)
**Rule 114:** IF position 1064 = **A** AND position 1065 = **Q** THEN co-evolutionary (MI = 2.288)
**Rule 115:** IF position 1064 = **W** AND position 1065 = **N** THEN co-evolutionary (MI = 2.288)
**Rule 116:** IF position 1064 = **Q** AND position 1065 = **K** THEN co-evolutionary (MI = 2.288)
**Rule 117:** IF position 1064 = **K** AND position 1065 = **A** THEN co-evolutionary (MI = 2.288)
**Rule 118:** IF position 1064 = **R** AND position 1065 = **I** THEN co-evolutionary (MI = 2.288)
**Rule 119:** IF position 1064 = **A** AND position 1065 = **S** THEN co-evolutionary (MI = 2.288)
**Rule 120:** IF position 1064 = **Y** AND position 1065 = **E** THEN co-evolutionary (MI = 2.288)
**Rule 121:** IF position 1064 = **W** AND position 1065 = **W** THEN co-evolutionary (MI = 2.288)

---

## Position Pair (413, 428)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2865 |
| Total mutations | 271 |
| Reference pos 413 | N |
| Reference pos 428 | C |
| On-set cells | 14 |
| Off-set cells | 385 |
| Don't-care cells | 1 |
| Prime implicants | 12 |
| Essential PIs | 11 |

### K-map (Compact View)

```
Position pair (413, 428): Reference = (N, C)

Co-evolutionary residue pairs (on-set):
  A-I, I-V, Y-W, D-A, Q-F, N-A, N-N, K-A, K-I, K-N, K-C, T-T, P-D, G-G

Don't-care positions (conserved): 1 cells
Never-seen pairs (off-set): 385 cells
```

### Boolean Function

```
f(pos_413, pos_428) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_413=A AND pos_428=I
  * PI_2: pos_413=I AND pos_428=F
  * PI_3: pos_413=H AND pos_428=K
  * PI_4: pos_413=A AND pos_428=M
    PI_5: pos_413=A AND pos_428=F
  * PI_6: pos_413=A AND pos_428=S
  * PI_7: pos_413=W AND pos_428=S
  * PI_8: pos_413=N AND pos_428=M
  * PI_9: pos_413=K AND pos_428=H
  * PI_10: pos_413=R AND pos_428=W
  * PI_11: pos_413=F AND pos_428=A
  * PI_12: pos_413=W AND pos_428=I

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| Q-F | 5.2234 | 0.0054 | co-evolutionary |
| Y-W | 4.3361 | 0.0131 | co-evolutionary |
| D-A | 3.4317 | 0.0308 | co-evolutionary |
| I-V | 3.3407 | 0.0354 | co-evolutionary |
| A-I | 2.4244 | 0.0878 | co-evolutionary |
| N-C | 0.1976 | 0.7914 | co-evolutionary |
| K-C | 0.1183 | 0.0269 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 122:** IF position 413 = **A** AND position 428 = **I** THEN co-evolutionary (MI = 2.287)
**Rule 123:** IF position 413 = **I** AND position 428 = **F** THEN co-evolutionary (MI = 2.287)
**Rule 124:** IF position 413 = **H** AND position 428 = **K** THEN co-evolutionary (MI = 2.287)
**Rule 125:** IF position 413 = **A** AND position 428 = **M** THEN co-evolutionary (MI = 2.287)
Rule 126: IF position 413 = A AND position 428 = F THEN co-evolutionary (MI = 2.287)
**Rule 127:** IF position 413 = **A** AND position 428 = **S** THEN co-evolutionary (MI = 2.287)
**Rule 128:** IF position 413 = **W** AND position 428 = **S** THEN co-evolutionary (MI = 2.287)
**Rule 129:** IF position 413 = **N** AND position 428 = **M** THEN co-evolutionary (MI = 2.287)
**Rule 130:** IF position 413 = **K** AND position 428 = **H** THEN co-evolutionary (MI = 2.287)
**Rule 131:** IF position 413 = **R** AND position 428 = **W** THEN co-evolutionary (MI = 2.287)
**Rule 132:** IF position 413 = **F** AND position 428 = **A** THEN co-evolutionary (MI = 2.287)
**Rule 133:** IF position 413 = **W** AND position 428 = **I** THEN co-evolutionary (MI = 2.287)

---

## Position Pair (459, 473)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2842 |
| Total mutations | 265 |
| Reference pos 459 | P |
| Reference pos 473 | N |
| On-set cells | 11 |
| Off-set cells | 388 |
| Don't-care cells | 1 |
| Prime implicants | 10 |
| Essential PIs | 10 |

### K-map (Compact View)

```
Position pair (459, 473): Reference = (P, N)

Co-evolutionary residue pairs (on-set):
  I-G, L-A, F-K, E-P, D-N, N-Q, K-I, K-G, R-C, T-A, P-S

Don't-care positions (conserved): 1 cells
Never-seen pairs (off-set): 388 cells
```

### Boolean Function

```
f(pos_459, pos_473) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_459=L AND pos_473=E
  * PI_2: pos_459=W AND pos_473=I
  * PI_3: pos_459=W AND pos_473=V
  * PI_4: pos_459=L AND pos_473=W
  * PI_5: pos_459=N AND pos_473=L
  * PI_6: pos_459=N AND pos_473=S
  * PI_7: pos_459=R AND pos_473=Y
  * PI_8: pos_459=A AND pos_473=F
  * PI_9: pos_459=I AND pos_473=W
  * PI_10: pos_459=M AND pos_473=A

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| N-Q | 5.2234 | 0.0054 | co-evolutionary |
| R-C | 3.3852 | 0.0339 | co-evolutionary |
| F-K | 3.0918 | 0.0454 | co-evolutionary |
| E-P | 2.4509 | 0.0862 | co-evolutionary |
| P-S | 0.2176 | 0.0085 | co-evolutionary |
| D-N | 0.2071 | 0.0169 | co-evolutionary |
| P-N | 0.1965 | 0.7960 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 134:** IF position 459 = **L** AND position 473 = **E** THEN co-evolutionary (MI = 2.284)
**Rule 135:** IF position 459 = **W** AND position 473 = **I** THEN co-evolutionary (MI = 2.284)
**Rule 136:** IF position 459 = **W** AND position 473 = **V** THEN co-evolutionary (MI = 2.284)
**Rule 137:** IF position 459 = **L** AND position 473 = **W** THEN co-evolutionary (MI = 2.284)
**Rule 138:** IF position 459 = **N** AND position 473 = **L** THEN co-evolutionary (MI = 2.284)
**Rule 139:** IF position 459 = **N** AND position 473 = **S** THEN co-evolutionary (MI = 2.284)
**Rule 140:** IF position 459 = **R** AND position 473 = **Y** THEN co-evolutionary (MI = 2.284)
**Rule 141:** IF position 459 = **A** AND position 473 = **F** THEN co-evolutionary (MI = 2.284)
**Rule 142:** IF position 459 = **I** AND position 473 = **W** THEN co-evolutionary (MI = 2.284)
**Rule 143:** IF position 459 = **M** AND position 473 = **A** THEN co-evolutionary (MI = 2.284)

---

## Position Pair (469, 473)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2842 |
| Total mutations | 265 |
| Reference pos 469 | Y |
| Reference pos 473 | N |
| On-set cells | 11 |
| Off-set cells | 388 |
| Don't-care cells | 1 |
| Prime implicants | 10 |
| Essential PIs | 10 |

### K-map (Compact View)

```
Position pair (469, 473): Reference = (Y, N)

Co-evolutionary residue pairs (on-set):
  A-P, I-I, I-G, Y-S, E-A, Q-K, N-N, K-G, T-Q, C-A, G-C

Don't-care positions (conserved): 1 cells
Never-seen pairs (off-set): 388 cells
```

### Boolean Function

```
f(pos_469, pos_473) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_469=I AND pos_473=F
  * PI_2: pos_469=E AND pos_473=V
  * PI_3: pos_469=I AND pos_473=L
  * PI_4: pos_469=L AND pos_473=W
  * PI_5: pos_469=Q AND pos_473=A
  * PI_6: pos_469=K AND pos_473=F
  * PI_7: pos_469=R AND pos_473=W
  * PI_8: pos_469=M AND pos_473=Q
  * PI_9: pos_469=F AND pos_473=M
  * PI_10: pos_469=E AND pos_473=K

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| T-Q | 5.2234 | 0.0054 | co-evolutionary |
| G-C | 3.3852 | 0.0339 | co-evolutionary |
| Q-K | 3.0918 | 0.0454 | co-evolutionary |
| A-P | 2.4509 | 0.0862 | co-evolutionary |
| Y-S | 0.2176 | 0.0085 | co-evolutionary |
| N-N | 0.2071 | 0.0169 | co-evolutionary |
| Y-N | 0.1965 | 0.7960 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 144:** IF position 469 = **I** AND position 473 = **F** THEN co-evolutionary (MI = 2.284)
**Rule 145:** IF position 469 = **E** AND position 473 = **V** THEN co-evolutionary (MI = 2.284)
**Rule 146:** IF position 469 = **I** AND position 473 = **L** THEN co-evolutionary (MI = 2.284)
**Rule 147:** IF position 469 = **L** AND position 473 = **W** THEN co-evolutionary (MI = 2.284)
**Rule 148:** IF position 469 = **Q** AND position 473 = **A** THEN co-evolutionary (MI = 2.284)
**Rule 149:** IF position 469 = **K** AND position 473 = **F** THEN co-evolutionary (MI = 2.284)
**Rule 150:** IF position 469 = **R** AND position 473 = **W** THEN co-evolutionary (MI = 2.284)
**Rule 151:** IF position 469 = **M** AND position 473 = **Q** THEN co-evolutionary (MI = 2.284)
**Rule 152:** IF position 469 = **F** AND position 473 = **M** THEN co-evolutionary (MI = 2.284)
**Rule 153:** IF position 469 = **E** AND position 473 = **K** THEN co-evolutionary (MI = 2.284)

---

## Position Pair (1064, 1066)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2816 |
| Total mutations | 278 |
| Reference pos 1064 | V |
| Reference pos 1066 | A |
| On-set cells | 11 |
| Off-set cells | 388 |
| Don't-care cells | 1 |
| Prime implicants | 9 |
| Essential PIs | 8 |

### K-map (Compact View)

```
Position pair (1064, 1066): Reference = (V, A)

Co-evolutionary residue pairs (on-set):
  A-E, L-V, V-Y, F-T, Y-P, E-N, Q-K, N-T, K-F, T-V, P-Q

Don't-care positions (conserved): 1 cells
Never-seen pairs (off-set): 388 cells
```

### Boolean Function

```
f(pos_1064, pos_1066) = 1 if ANY of these residue pairs appear:

  * PI_1: pos_1064=A AND pos_1066=E
  * PI_2: pos_1064=L AND pos_1066=N
  * PI_3: pos_1064=M AND pos_1066=L
    PI_4: pos_1064=V AND pos_1066=H
  * PI_5: pos_1064=W AND pos_1066=M
  * PI_6: pos_1064=E AND pos_1066=Q
  * PI_7: pos_1064=K AND pos_1066=F
  * PI_8: pos_1064=R AND pos_1066=H
  * PI_9: pos_1064=W AND pos_1066=L

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| K-F | 4.6844 | 0.0092 | co-evolutionary |
| E-N | 3.9113 | 0.0200 | co-evolutionary |
| Q-K | 3.4805 | 0.0308 | co-evolutionary |
| P-Q | 2.9647 | 0.0516 | co-evolutionary |
| A-E | 2.4072 | 0.0901 | co-evolutionary |
| V-A | 0.2340 | 0.7860 | co-evolutionary |
| V-Y | 0.2340 | 0.0054 | co-evolutionary |

### Inference Rules (Natural Language)

**Rule 154:** IF position 1064 = **A** AND position 1066 = **E** THEN co-evolutionary (MI = 2.282)
**Rule 155:** IF position 1064 = **L** AND position 1066 = **N** THEN co-evolutionary (MI = 2.282)
**Rule 156:** IF position 1064 = **M** AND position 1066 = **L** THEN co-evolutionary (MI = 2.282)
Rule 157: IF position 1064 = V AND position 1066 = H THEN co-evolutionary (MI = 2.282)
**Rule 158:** IF position 1064 = **W** AND position 1066 = **M** THEN co-evolutionary (MI = 2.282)
**Rule 159:** IF position 1064 = **E** AND position 1066 = **Q** THEN co-evolutionary (MI = 2.282)
**Rule 160:** IF position 1064 = **K** AND position 1066 = **F** THEN co-evolutionary (MI = 2.282)
**Rule 161:** IF position 1064 = **R** AND position 1066 = **H** THEN co-evolutionary (MI = 2.282)
**Rule 162:** IF position 1064 = **W** AND position 1066 = **L** THEN co-evolutionary (MI = 2.282)

---

## Summary

| Metric | Value |
|--------|-------|
| Sequences | 1299 |
| Variable positions | 57 |
| Co-evolutionary pairs | 36918 |
| Total inference rules | 162 |
| Position pairs with rules | 15 |

## How to Apply

1. Extract residues at positions 68-79 from a new sequence
2. For each position pair, check if the residue pair matches any rule
3. If YES: that position pair is co-evolutionary
4. If position i mutates: find which residue at position j satisfies the co-evolutionary constraint

**Example:** If position 76 mutates to Y, check rules for position 76.
Rule 2 says: IF pos 76 = Y AND pos 77 = K THEN co-evolutionary.
So position 77 must also mutate to K.