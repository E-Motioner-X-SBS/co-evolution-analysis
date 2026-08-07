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

## Position Pair (495, 498)

| Property | Value |
|----------|-------|
| Mutual Information | 0.8710 |
| Total mutations | 425 |
| Reference pos 495 | R |
| Reference pos 498 | G |
| On-set cells | 2 |
| Off-set cells | 0 |
| Don't-care cells | 1022 |
| Prime implicants | 20 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (495, 498): Reference = (R, G)

Co-evolutionary residue pairs (on-set):
  Q-G, R-S

Don't-care positions (conserved): 1022 cells
Never-seen pairs (off-set): 0 cells
```

### Boolean Function

```
f(pos_495, pos_498) = 1 if ANY of these residue pairs appear:

    PI_1: pos_495=I AND pos_498=M
    PI_2: pos_495=M AND pos_498=M
    PI_3: pos_495=F AND pos_498=A
    PI_4: pos_495=F AND pos_498=M
    PI_5: pos_495=F AND pos_498=M
    PI_6: pos_495=F AND pos_498=M
    PI_7: pos_495=F AND pos_498=M
    PI_8: pos_495=F AND pos_498=M
    PI_9: pos_495=F AND pos_498=M
    PI_10: pos_495=F AND pos_498=M
    PI_11: pos_495=V AND pos_498=V
    PI_12: pos_495=F AND pos_498=V
    PI_13: pos_495=Y AND pos_498=V
    PI_14: pos_495=W AND pos_498=I
    PI_15: pos_495=W AND pos_498=L
    PI_16: pos_495=W AND pos_498=V
    PI_17: pos_495=W AND pos_498=V
    PI_18: pos_495=W AND pos_498=V
    PI_19: pos_495=W AND pos_498=V
    PI_20: pos_495=W AND pos_498=V

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| Q-G | 0.2641 | 0.0956 | co-evolutionary |
| R-S | 0.1005 | 0.2321 | co-evolutionary |
| R-G | -0.0325 | 0.6723 | anti-correlated |

### Inference Rules (Natural Language)

Rule 1: IF position 495 = I AND position 498 = M THEN co-evolutionary (MI = 0.871)
Rule 2: IF position 495 = M AND position 498 = M THEN co-evolutionary (MI = 0.871)
Rule 3: IF position 495 = F AND position 498 = A THEN co-evolutionary (MI = 0.871)
Rule 4: IF position 495 = F AND position 498 = M THEN co-evolutionary (MI = 0.871)
Rule 5: IF position 495 = F AND position 498 = M THEN co-evolutionary (MI = 0.871)
Rule 6: IF position 495 = F AND position 498 = M THEN co-evolutionary (MI = 0.871)
Rule 7: IF position 495 = F AND position 498 = M THEN co-evolutionary (MI = 0.871)
Rule 8: IF position 495 = F AND position 498 = M THEN co-evolutionary (MI = 0.871)
Rule 9: IF position 495 = F AND position 498 = M THEN co-evolutionary (MI = 0.871)
Rule 10: IF position 495 = F AND position 498 = M THEN co-evolutionary (MI = 0.871)
Rule 11: IF position 495 = V AND position 498 = V THEN co-evolutionary (MI = 0.871)
Rule 12: IF position 495 = F AND position 498 = V THEN co-evolutionary (MI = 0.871)
Rule 13: IF position 495 = Y AND position 498 = V THEN co-evolutionary (MI = 0.871)
Rule 14: IF position 495 = W AND position 498 = I THEN co-evolutionary (MI = 0.871)
Rule 15: IF position 495 = W AND position 498 = L THEN co-evolutionary (MI = 0.871)
Rule 16: IF position 495 = W AND position 498 = V THEN co-evolutionary (MI = 0.871)
Rule 17: IF position 495 = W AND position 498 = V THEN co-evolutionary (MI = 0.871)
Rule 18: IF position 495 = W AND position 498 = V THEN co-evolutionary (MI = 0.871)
Rule 19: IF position 495 = W AND position 498 = V THEN co-evolutionary (MI = 0.871)
Rule 20: IF position 495 = W AND position 498 = V THEN co-evolutionary (MI = 0.871)

---

## Position Pair (448, 454)

| Property | Value |
|----------|-------|
| Mutual Information | 0.8344 |
| Total mutations | 393 |
| Reference pos 448 | G |
| Reference pos 454 | L |
| On-set cells | 3 |
| Off-set cells | 0 |
| Don't-care cells | 1021 |
| Prime implicants | 22 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (448, 454): Reference = (G, L)

Co-evolutionary residue pairs (on-set):
  S-L, S-R, G-R

Don't-care positions (conserved): 1021 cells
Never-seen pairs (off-set): 0 cells
```

### Boolean Function

```
f(pos_448, pos_454) = 1 if ANY of these residue pairs appear:

    PI_1: pos_448=I AND pos_454=N
    PI_2: pos_448=W AND pos_454=E
    PI_3: pos_448=W AND pos_454=E
    PI_4: pos_448=D AND pos_454=Q
    PI_5: pos_448=D AND pos_454=D
    PI_6: pos_448=W AND pos_454=E
    PI_7: pos_448=W AND pos_454=E
    PI_8: pos_448=E AND pos_454=N
    PI_9: pos_448=D AND pos_454=V
    PI_10: pos_448=D AND pos_454=N
    PI_11: pos_448=V AND pos_454=E
    PI_12: pos_448=F AND pos_454=E
    PI_13: pos_448=Y AND pos_454=E
    PI_14: pos_448=W AND pos_454=A
    PI_15: pos_448=W AND pos_454=E
    PI_16: pos_448=Y AND pos_454=N
    PI_17: pos_448=W AND pos_454=V
    PI_18: pos_448=W AND pos_454=N
    PI_19: pos_448=D AND pos_454=N
    PI_20: pos_448=D AND pos_454=N
    PI_21: pos_448=W AND pos_454=N
    PI_22: pos_448=W AND pos_454=N

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| G-R | 0.2373 | 0.0849 | co-evolutionary |
| S-L | 0.0860 | 0.2176 | co-evolutionary |
| G-L | -0.0254 | 0.6968 | anti-correlated |

### Inference Rules (Natural Language)

Rule 21: IF position 448 = I AND position 454 = N THEN co-evolutionary (MI = 0.834)
Rule 22: IF position 448 = W AND position 454 = E THEN co-evolutionary (MI = 0.834)
Rule 23: IF position 448 = W AND position 454 = E THEN co-evolutionary (MI = 0.834)
Rule 24: IF position 448 = D AND position 454 = Q THEN co-evolutionary (MI = 0.834)
Rule 25: IF position 448 = D AND position 454 = D THEN co-evolutionary (MI = 0.834)
Rule 26: IF position 448 = W AND position 454 = E THEN co-evolutionary (MI = 0.834)
Rule 27: IF position 448 = W AND position 454 = E THEN co-evolutionary (MI = 0.834)
Rule 28: IF position 448 = E AND position 454 = N THEN co-evolutionary (MI = 0.834)
Rule 29: IF position 448 = D AND position 454 = V THEN co-evolutionary (MI = 0.834)
Rule 30: IF position 448 = D AND position 454 = N THEN co-evolutionary (MI = 0.834)
Rule 31: IF position 448 = V AND position 454 = E THEN co-evolutionary (MI = 0.834)
Rule 32: IF position 448 = F AND position 454 = E THEN co-evolutionary (MI = 0.834)
Rule 33: IF position 448 = Y AND position 454 = E THEN co-evolutionary (MI = 0.834)
Rule 34: IF position 448 = W AND position 454 = A THEN co-evolutionary (MI = 0.834)
Rule 35: IF position 448 = W AND position 454 = E THEN co-evolutionary (MI = 0.834)
Rule 36: IF position 448 = Y AND position 454 = N THEN co-evolutionary (MI = 0.834)
Rule 37: IF position 448 = W AND position 454 = V THEN co-evolutionary (MI = 0.834)
Rule 38: IF position 448 = W AND position 454 = N THEN co-evolutionary (MI = 0.834)
Rule 39: IF position 448 = D AND position 454 = N THEN co-evolutionary (MI = 0.834)
Rule 40: IF position 448 = D AND position 454 = N THEN co-evolutionary (MI = 0.834)
Rule 41: IF position 448 = W AND position 454 = N THEN co-evolutionary (MI = 0.834)
Rule 42: IF position 448 = W AND position 454 = N THEN co-evolutionary (MI = 0.834)

---

## Position Pair (488, 498)

| Property | Value |
|----------|-------|
| Mutual Information | 0.8219 |
| Total mutations | 405 |
| Reference pos 488 | F |
| Reference pos 498 | G |
| On-set cells | 3 |
| Off-set cells | 0 |
| Don't-care cells | 1021 |
| Prime implicants | 27 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (488, 498): Reference = (F, G)

Co-evolutionary residue pairs (on-set):
  V-G, F-S, P-G

Don't-care positions (conserved): 1021 cells
Never-seen pairs (off-set): 0 cells
```

### Boolean Function

```
f(pos_488, pos_498) = 1 if ANY of these residue pairs appear:

    PI_1: pos_488=I AND pos_498=M
    PI_2: pos_488=I AND pos_498=H
    PI_3: pos_488=D AND pos_498=M
    PI_4: pos_488=L AND pos_498=Q
    PI_5: pos_488=I AND pos_498=H
    PI_6: pos_488=D AND pos_498=M
    PI_7: pos_488=L AND pos_498=D
    PI_8: pos_488=I AND pos_498=H
    PI_9: pos_488=D AND pos_498=M
    PI_10: pos_488=A AND pos_498=N
    PI_11: pos_488=I AND pos_498=H
    PI_12: pos_488=D AND pos_498=M
    PI_13: pos_488=L AND pos_498=V
    PI_14: pos_488=L AND pos_498=N
    PI_15: pos_488=L AND pos_498=N
    PI_16: pos_488=A AND pos_498=H
    PI_17: pos_488=E AND pos_498=M
    PI_18: pos_488=I AND pos_498=E
    PI_19: pos_488=D AND pos_498=A
    PI_20: pos_488=I AND pos_498=H
    PI_21: pos_488=D AND pos_498=M
    PI_22: pos_488=I AND pos_498=H
    PI_23: pos_488=D AND pos_498=M
    PI_24: pos_488=L AND pos_498=N
    PI_25: pos_488=L AND pos_498=N
    PI_26: pos_488=L AND pos_498=N
    PI_27: pos_488=L AND pos_498=N

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| V-G | 0.2638 | 0.0786 | co-evolutionary |
| F-S | 0.0835 | 0.2319 | co-evolutionary |
| F-G | -0.0266 | 0.6880 | anti-correlated |

### Inference Rules (Natural Language)

Rule 43: IF position 488 = I AND position 498 = M THEN co-evolutionary (MI = 0.822)
Rule 44: IF position 488 = I AND position 498 = H THEN co-evolutionary (MI = 0.822)
Rule 45: IF position 488 = D AND position 498 = M THEN co-evolutionary (MI = 0.822)
Rule 46: IF position 488 = L AND position 498 = Q THEN co-evolutionary (MI = 0.822)
Rule 47: IF position 488 = I AND position 498 = H THEN co-evolutionary (MI = 0.822)
Rule 48: IF position 488 = D AND position 498 = M THEN co-evolutionary (MI = 0.822)
Rule 49: IF position 488 = L AND position 498 = D THEN co-evolutionary (MI = 0.822)
Rule 50: IF position 488 = I AND position 498 = H THEN co-evolutionary (MI = 0.822)
Rule 51: IF position 488 = D AND position 498 = M THEN co-evolutionary (MI = 0.822)
Rule 52: IF position 488 = A AND position 498 = N THEN co-evolutionary (MI = 0.822)
Rule 53: IF position 488 = I AND position 498 = H THEN co-evolutionary (MI = 0.822)
Rule 54: IF position 488 = D AND position 498 = M THEN co-evolutionary (MI = 0.822)
Rule 55: IF position 488 = L AND position 498 = V THEN co-evolutionary (MI = 0.822)
Rule 56: IF position 488 = L AND position 498 = N THEN co-evolutionary (MI = 0.822)
Rule 57: IF position 488 = L AND position 498 = N THEN co-evolutionary (MI = 0.822)
Rule 58: IF position 488 = A AND position 498 = H THEN co-evolutionary (MI = 0.822)
Rule 59: IF position 488 = E AND position 498 = M THEN co-evolutionary (MI = 0.822)
Rule 60: IF position 488 = I AND position 498 = E THEN co-evolutionary (MI = 0.822)
Rule 61: IF position 488 = D AND position 498 = A THEN co-evolutionary (MI = 0.822)
Rule 62: IF position 488 = I AND position 498 = H THEN co-evolutionary (MI = 0.822)
Rule 63: IF position 488 = D AND position 498 = M THEN co-evolutionary (MI = 0.822)
Rule 64: IF position 488 = I AND position 498 = H THEN co-evolutionary (MI = 0.822)
Rule 65: IF position 488 = D AND position 498 = M THEN co-evolutionary (MI = 0.822)
Rule 66: IF position 488 = L AND position 498 = N THEN co-evolutionary (MI = 0.822)
Rule 67: IF position 488 = L AND position 498 = N THEN co-evolutionary (MI = 0.822)
Rule 68: IF position 488 = L AND position 498 = N THEN co-evolutionary (MI = 0.822)
Rule 69: IF position 488 = L AND position 498 = N THEN co-evolutionary (MI = 0.822)

---

## Position Pair (442, 454)

| Property | Value |
|----------|-------|
| Mutual Information | 0.8110 |
| Total mutations | 188 |
| Reference pos 442 | K |
| Reference pos 454 | L |
| On-set cells | 3 |
| Off-set cells | 0 |
| Don't-care cells | 1021 |
| Prime implicants | 24 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (442, 454): Reference = (K, L)

Co-evolutionary residue pairs (on-set):
  N-L, N-R, K-R

Don't-care positions (conserved): 1021 cells
Never-seen pairs (off-set): 0 cells
```

### Boolean Function

```
f(pos_442, pos_454) = 1 if ANY of these residue pairs appear:

    PI_1: pos_442=F AND pos_454=E
    PI_2: pos_442=M AND pos_454=N
    PI_3: pos_442=I AND pos_454=E
    PI_4: pos_442=M AND pos_454=E
    PI_5: pos_442=F AND pos_454=A
    PI_6: pos_442=F AND pos_454=E
    PI_7: pos_442=F AND pos_454=E
    PI_8: pos_442=F AND pos_454=E
    PI_9: pos_442=F AND pos_454=E
    PI_10: pos_442=F AND pos_454=E
    PI_11: pos_442=I AND pos_454=N
    PI_12: pos_442=L AND pos_454=N
    PI_13: pos_442=F AND pos_454=V
    PI_14: pos_442=Y AND pos_454=V
    PI_15: pos_442=Y AND pos_454=D
    PI_16: pos_442=Y AND pos_454=Q
    PI_17: pos_442=F AND pos_454=N
    PI_18: pos_442=Y AND pos_454=N
    PI_19: pos_442=F AND pos_454=N
    PI_20: pos_442=Y AND pos_454=N
    PI_21: pos_442=F AND pos_454=N
    PI_22: pos_442=Y AND pos_454=N
    PI_23: pos_442=F AND pos_454=N
    PI_24: pos_442=Y AND pos_454=N

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| N-L | 0.0158 | 0.0602 | co-evolutionary |
| K-R | 0.0115 | 0.0821 | co-evolutionary |
| K-L | -0.0011 | 0.8530 | anti-correlated |

### Inference Rules (Natural Language)

Rule 70: IF position 442 = F AND position 454 = E THEN co-evolutionary (MI = 0.811)
Rule 71: IF position 442 = M AND position 454 = N THEN co-evolutionary (MI = 0.811)
Rule 72: IF position 442 = I AND position 454 = E THEN co-evolutionary (MI = 0.811)
Rule 73: IF position 442 = M AND position 454 = E THEN co-evolutionary (MI = 0.811)
Rule 74: IF position 442 = F AND position 454 = A THEN co-evolutionary (MI = 0.811)
Rule 75: IF position 442 = F AND position 454 = E THEN co-evolutionary (MI = 0.811)
Rule 76: IF position 442 = F AND position 454 = E THEN co-evolutionary (MI = 0.811)
Rule 77: IF position 442 = F AND position 454 = E THEN co-evolutionary (MI = 0.811)
Rule 78: IF position 442 = F AND position 454 = E THEN co-evolutionary (MI = 0.811)
Rule 79: IF position 442 = F AND position 454 = E THEN co-evolutionary (MI = 0.811)
Rule 80: IF position 442 = I AND position 454 = N THEN co-evolutionary (MI = 0.811)
Rule 81: IF position 442 = L AND position 454 = N THEN co-evolutionary (MI = 0.811)
Rule 82: IF position 442 = F AND position 454 = V THEN co-evolutionary (MI = 0.811)
Rule 83: IF position 442 = Y AND position 454 = V THEN co-evolutionary (MI = 0.811)
Rule 84: IF position 442 = Y AND position 454 = D THEN co-evolutionary (MI = 0.811)
Rule 85: IF position 442 = Y AND position 454 = Q THEN co-evolutionary (MI = 0.811)
Rule 86: IF position 442 = F AND position 454 = N THEN co-evolutionary (MI = 0.811)
Rule 87: IF position 442 = Y AND position 454 = N THEN co-evolutionary (MI = 0.811)
Rule 88: IF position 442 = F AND position 454 = N THEN co-evolutionary (MI = 0.811)
Rule 89: IF position 442 = Y AND position 454 = N THEN co-evolutionary (MI = 0.811)
Rule 90: IF position 442 = F AND position 454 = N THEN co-evolutionary (MI = 0.811)
Rule 91: IF position 442 = Y AND position 454 = N THEN co-evolutionary (MI = 0.811)
Rule 92: IF position 442 = F AND position 454 = N THEN co-evolutionary (MI = 0.811)
Rule 93: IF position 442 = Y AND position 454 = N THEN co-evolutionary (MI = 0.811)

---

## Position Pair (442, 448)

| Property | Value |
|----------|-------|
| Mutual Information | 0.7284 |
| Total mutations | 363 |
| Reference pos 442 | K |
| Reference pos 448 | G |
| On-set cells | 3 |
| Off-set cells | 0 |
| Don't-care cells | 1021 |
| Prime implicants | 22 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (442, 448): Reference = (K, G)

Co-evolutionary residue pairs (on-set):
  N-S, N-G, K-S

Don't-care positions (conserved): 1021 cells
Never-seen pairs (off-set): 0 cells
```

### Boolean Function

```
f(pos_442, pos_448) = 1 if ANY of these residue pairs appear:

    PI_1: pos_442=F AND pos_448=E
    PI_2: pos_442=M AND pos_448=N
    PI_3: pos_442=I AND pos_448=N
    PI_4: pos_442=L AND pos_448=N
    PI_5: pos_442=F AND pos_448=V
    PI_6: pos_442=Y AND pos_448=V
    PI_7: pos_442=Y AND pos_448=D
    PI_8: pos_442=Y AND pos_448=Q
    PI_9: pos_442=F AND pos_448=N
    PI_10: pos_442=Y AND pos_448=N
    PI_11: pos_442=F AND pos_448=N
    PI_12: pos_442=Y AND pos_448=N
    PI_13: pos_442=F AND pos_448=H
    PI_14: pos_442=F AND pos_448=H
    PI_15: pos_442=Y AND pos_448=N
    PI_16: pos_442=F AND pos_448=N
    PI_17: pos_442=Y AND pos_448=N
    PI_18: pos_442=I AND pos_448=H
    PI_19: pos_442=M AND pos_448=H
    PI_20: pos_442=F AND pos_448=M
    PI_21: pos_442=F AND pos_448=H
    PI_22: pos_442=F AND pos_448=H

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| N-G | 0.2243 | 0.0632 | co-evolutionary |
| K-S | 0.0599 | 0.2186 | co-evolutionary |
| K-G | -0.0176 | 0.7166 | anti-correlated |

### Inference Rules (Natural Language)

Rule 94: IF position 442 = F AND position 448 = E THEN co-evolutionary (MI = 0.728)
Rule 95: IF position 442 = M AND position 448 = N THEN co-evolutionary (MI = 0.728)
Rule 96: IF position 442 = I AND position 448 = N THEN co-evolutionary (MI = 0.728)
Rule 97: IF position 442 = L AND position 448 = N THEN co-evolutionary (MI = 0.728)
Rule 98: IF position 442 = F AND position 448 = V THEN co-evolutionary (MI = 0.728)
Rule 99: IF position 442 = Y AND position 448 = V THEN co-evolutionary (MI = 0.728)
Rule 100: IF position 442 = Y AND position 448 = D THEN co-evolutionary (MI = 0.728)
Rule 101: IF position 442 = Y AND position 448 = Q THEN co-evolutionary (MI = 0.728)
Rule 102: IF position 442 = F AND position 448 = N THEN co-evolutionary (MI = 0.728)
Rule 103: IF position 442 = Y AND position 448 = N THEN co-evolutionary (MI = 0.728)
Rule 104: IF position 442 = F AND position 448 = N THEN co-evolutionary (MI = 0.728)
Rule 105: IF position 442 = Y AND position 448 = N THEN co-evolutionary (MI = 0.728)
Rule 106: IF position 442 = F AND position 448 = H THEN co-evolutionary (MI = 0.728)
Rule 107: IF position 442 = F AND position 448 = H THEN co-evolutionary (MI = 0.728)
Rule 108: IF position 442 = Y AND position 448 = N THEN co-evolutionary (MI = 0.728)
Rule 109: IF position 442 = F AND position 448 = N THEN co-evolutionary (MI = 0.728)
Rule 110: IF position 442 = Y AND position 448 = N THEN co-evolutionary (MI = 0.728)
Rule 111: IF position 442 = I AND position 448 = H THEN co-evolutionary (MI = 0.728)
Rule 112: IF position 442 = M AND position 448 = H THEN co-evolutionary (MI = 0.728)
Rule 113: IF position 442 = F AND position 448 = M THEN co-evolutionary (MI = 0.728)
Rule 114: IF position 442 = F AND position 448 = H THEN co-evolutionary (MI = 0.728)
Rule 115: IF position 442 = F AND position 448 = H THEN co-evolutionary (MI = 0.728)

---

## Position Pair (212, 215)

| Property | Value |
|----------|-------|
| Mutual Information | 0.3977 |
| Total mutations | 300 |
| Reference pos 212 | V |
| Reference pos 215 | G |
| On-set cells | 4 |
| Off-set cells | 0 |
| Don't-care cells | 1020 |
| Prime implicants | 27 |
| Essential PIs | 1 |

### K-map (Compact View)

```
Position pair (212, 215): Reference = (V, G)

Co-evolutionary residue pairs (on-set):
  I-V, L-V, V-P, S-G

Don't-care positions (conserved): 1020 cells
Never-seen pairs (off-set): 0 cells
```

### Boolean Function

```
f(pos_212, pos_215) = 1 if ANY of these residue pairs appear:

    PI_1: pos_212=A AND pos_215=E
    PI_2: pos_212=I AND pos_215=A
  * PI_3: pos_212=I AND pos_215=H
    PI_4: pos_212=A AND pos_215=A
    PI_5: pos_212=A AND pos_215=E
    PI_6: pos_212=I AND pos_215=A
    PI_7: pos_212=A AND pos_215=E
    PI_8: pos_212=I AND pos_215=A
    PI_9: pos_212=A AND pos_215=E
    PI_10: pos_212=I AND pos_215=A
    PI_11: pos_212=A AND pos_215=E
    PI_12: pos_212=I AND pos_215=A
    PI_13: pos_212=A AND pos_215=E
    PI_14: pos_212=I AND pos_215=A
    PI_15: pos_212=I AND pos_215=H
    PI_16: pos_212=A AND pos_215=E
    PI_17: pos_212=I AND pos_215=A
    PI_18: pos_212=I AND pos_215=H
    PI_19: pos_212=I AND pos_215=H
    PI_20: pos_212=I AND pos_215=H
    PI_21: pos_212=Y AND pos_215=H
    PI_22: pos_212=W AND pos_215=M
    PI_23: pos_212=W AND pos_215=E
    PI_24: pos_212=W AND pos_215=H
    PI_25: pos_212=W AND pos_215=H
    PI_26: pos_212=W AND pos_215=H
    PI_27: pos_212=W AND pos_215=H

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

Rule 116: IF position 212 = A AND position 215 = E THEN co-evolutionary (MI = 0.398)
Rule 117: IF position 212 = I AND position 215 = A THEN co-evolutionary (MI = 0.398)
**Rule 118:** IF position 212 = **I** AND position 215 = **H** THEN co-evolutionary (MI = 0.398)
Rule 119: IF position 212 = A AND position 215 = A THEN co-evolutionary (MI = 0.398)
Rule 120: IF position 212 = A AND position 215 = E THEN co-evolutionary (MI = 0.398)
Rule 121: IF position 212 = I AND position 215 = A THEN co-evolutionary (MI = 0.398)
Rule 122: IF position 212 = A AND position 215 = E THEN co-evolutionary (MI = 0.398)
Rule 123: IF position 212 = I AND position 215 = A THEN co-evolutionary (MI = 0.398)
Rule 124: IF position 212 = A AND position 215 = E THEN co-evolutionary (MI = 0.398)
Rule 125: IF position 212 = I AND position 215 = A THEN co-evolutionary (MI = 0.398)
Rule 126: IF position 212 = A AND position 215 = E THEN co-evolutionary (MI = 0.398)
Rule 127: IF position 212 = I AND position 215 = A THEN co-evolutionary (MI = 0.398)
Rule 128: IF position 212 = A AND position 215 = E THEN co-evolutionary (MI = 0.398)
Rule 129: IF position 212 = I AND position 215 = A THEN co-evolutionary (MI = 0.398)
Rule 130: IF position 212 = I AND position 215 = H THEN co-evolutionary (MI = 0.398)
Rule 131: IF position 212 = A AND position 215 = E THEN co-evolutionary (MI = 0.398)
Rule 132: IF position 212 = I AND position 215 = A THEN co-evolutionary (MI = 0.398)
Rule 133: IF position 212 = I AND position 215 = H THEN co-evolutionary (MI = 0.398)
Rule 134: IF position 212 = I AND position 215 = H THEN co-evolutionary (MI = 0.398)
Rule 135: IF position 212 = I AND position 215 = H THEN co-evolutionary (MI = 0.398)
Rule 136: IF position 212 = Y AND position 215 = H THEN co-evolutionary (MI = 0.398)
Rule 137: IF position 212 = W AND position 215 = M THEN co-evolutionary (MI = 0.398)
Rule 138: IF position 212 = W AND position 215 = E THEN co-evolutionary (MI = 0.398)
Rule 139: IF position 212 = W AND position 215 = H THEN co-evolutionary (MI = 0.398)
Rule 140: IF position 212 = W AND position 215 = H THEN co-evolutionary (MI = 0.398)
Rule 141: IF position 212 = W AND position 215 = H THEN co-evolutionary (MI = 0.398)
Rule 142: IF position 212 = W AND position 215 = H THEN co-evolutionary (MI = 0.398)

---

## Position Pair (215, 216)

| Property | Value |
|----------|-------|
| Mutual Information | 0.3773 |
| Total mutations | 301 |
| Reference pos 215 | G |
| Reference pos 216 | R |
| On-set cells | 4 |
| Off-set cells | 0 |
| Don't-care cells | 1020 |
| Prime implicants | 32 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (215, 216): Reference = (G, R)

Co-evolutionary residue pairs (on-set):
  V-R, E-R, P-E, P-K

Don't-care positions (conserved): 1020 cells
Never-seen pairs (off-set): 0 cells
```

### Boolean Function

```
f(pos_215, pos_216) = 1 if ANY of these residue pairs appear:

    PI_1: pos_215=A AND pos_216=V
    PI_2: pos_215=D AND pos_216=L
    PI_3: pos_215=E AND pos_216=V
    PI_4: pos_215=I AND pos_216=D
    PI_5: pos_215=D AND pos_216=I
    PI_6: pos_215=I AND pos_216=Q
    PI_7: pos_215=I AND pos_216=N
    PI_8: pos_215=I AND pos_216=V
    PI_9: pos_215=M AND pos_216=V
    PI_10: pos_215=M AND pos_216=V
    PI_11: pos_215=M AND pos_216=V
    PI_12: pos_215=M AND pos_216=V
    PI_13: pos_215=I AND pos_216=N
    PI_14: pos_215=I AND pos_216=N
    PI_15: pos_215=D AND pos_216=V
    PI_16: pos_215=I AND pos_216=N
    PI_17: pos_215=D AND pos_216=V
    PI_18: pos_215=D AND pos_216=V
    PI_19: pos_215=I AND pos_216=N
    PI_20: pos_215=D AND pos_216=V
    PI_21: pos_215=D AND pos_216=V
    PI_22: pos_215=D AND pos_216=L
    PI_23: pos_215=M AND pos_216=I
    PI_24: pos_215=M AND pos_216=L
    PI_25: pos_215=M AND pos_216=V
    PI_26: pos_215=D AND pos_216=L
    PI_27: pos_215=D AND pos_216=L
    PI_28: pos_215=D AND pos_216=L
    PI_29: pos_215=D AND pos_216=L
    PI_30: pos_215=I AND pos_216=L
    PI_31: pos_215=E AND pos_216=L
    PI_32: pos_215=D AND pos_216=A

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| P-E | 1.5218 | 0.2175 | co-evolutionary |
| V-R | 0.2463 | 0.0164 | co-evolutionary |
| G-R | 0.2463 | 0.7645 | co-evolutionary |

### Inference Rules (Natural Language)

Rule 143: IF position 215 = A AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 144: IF position 215 = D AND position 216 = L THEN co-evolutionary (MI = 0.377)
Rule 145: IF position 215 = E AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 146: IF position 215 = I AND position 216 = D THEN co-evolutionary (MI = 0.377)
Rule 147: IF position 215 = D AND position 216 = I THEN co-evolutionary (MI = 0.377)
Rule 148: IF position 215 = I AND position 216 = Q THEN co-evolutionary (MI = 0.377)
Rule 149: IF position 215 = I AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 150: IF position 215 = I AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 151: IF position 215 = M AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 152: IF position 215 = M AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 153: IF position 215 = M AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 154: IF position 215 = M AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 155: IF position 215 = I AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 156: IF position 215 = I AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 157: IF position 215 = D AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 158: IF position 215 = I AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 159: IF position 215 = D AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 160: IF position 215 = D AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 161: IF position 215 = I AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 162: IF position 215 = D AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 163: IF position 215 = D AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 164: IF position 215 = D AND position 216 = L THEN co-evolutionary (MI = 0.377)
Rule 165: IF position 215 = M AND position 216 = I THEN co-evolutionary (MI = 0.377)
Rule 166: IF position 215 = M AND position 216 = L THEN co-evolutionary (MI = 0.377)
Rule 167: IF position 215 = M AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 168: IF position 215 = D AND position 216 = L THEN co-evolutionary (MI = 0.377)
Rule 169: IF position 215 = D AND position 216 = L THEN co-evolutionary (MI = 0.377)
Rule 170: IF position 215 = D AND position 216 = L THEN co-evolutionary (MI = 0.377)
Rule 171: IF position 215 = D AND position 216 = L THEN co-evolutionary (MI = 0.377)
Rule 172: IF position 215 = I AND position 216 = L THEN co-evolutionary (MI = 0.377)
Rule 173: IF position 215 = E AND position 216 = L THEN co-evolutionary (MI = 0.377)
Rule 174: IF position 215 = D AND position 216 = A THEN co-evolutionary (MI = 0.377)

---

## Position Pair (212, 216)

| Property | Value |
|----------|-------|
| Mutual Information | 0.3773 |
| Total mutations | 301 |
| Reference pos 212 | V |
| Reference pos 216 | R |
| On-set cells | 6 |
| Off-set cells | 0 |
| Don't-care cells | 1018 |
| Prime implicants | 32 |
| Essential PIs | 1 |

### K-map (Compact View)

```
Position pair (212, 216): Reference = (V, R)

Co-evolutionary residue pairs (on-set):
  I-R, L-R, V-E, V-K, S-R, G-R

Don't-care positions (conserved): 1018 cells
Never-seen pairs (off-set): 0 cells
```

### Boolean Function

```
f(pos_212, pos_216) = 1 if ANY of these residue pairs appear:

    PI_1: pos_212=I AND pos_216=L
    PI_2: pos_212=I AND pos_216=V
    PI_3: pos_212=A AND pos_216=Q
    PI_4: pos_212=I AND pos_216=Q
    PI_5: pos_212=A AND pos_216=N
    PI_6: pos_212=I AND pos_216=N
    PI_7: pos_212=A AND pos_216=N
    PI_8: pos_212=I AND pos_216=V
    PI_9: pos_212=I AND pos_216=N
    PI_10: pos_212=A AND pos_216=V
    PI_11: pos_212=I AND pos_216=Q
    PI_12: pos_212=A AND pos_216=N
    PI_13: pos_212=I AND pos_216=V
    PI_14: pos_212=A AND pos_216=D
    PI_15: pos_212=I AND pos_216=I
    PI_16: pos_212=I AND pos_216=Q
    PI_17: pos_212=I AND pos_216=Q
    PI_18: pos_212=I AND pos_216=Q
    PI_19: pos_212=A AND pos_216=N
    PI_20: pos_212=I AND pos_216=V
    PI_21: pos_212=I AND pos_216=D
    PI_22: pos_212=I AND pos_216=E
    PI_23: pos_212=I AND pos_216=N
    PI_24: pos_212=I AND pos_216=N
    PI_25: pos_212=I AND pos_216=N
    PI_26: pos_212=D AND pos_216=D
    PI_27: pos_212=W AND pos_216=D
    PI_28: pos_212=W AND pos_216=Q
    PI_29: pos_212=W AND pos_216=N
    PI_30: pos_212=D AND pos_216=N
    PI_31: pos_212=W AND pos_216=N
  * PI_32: pos_212=W AND pos_216=N

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

Rule 175: IF position 212 = I AND position 216 = L THEN co-evolutionary (MI = 0.377)
Rule 176: IF position 212 = I AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 177: IF position 212 = A AND position 216 = Q THEN co-evolutionary (MI = 0.377)
Rule 178: IF position 212 = I AND position 216 = Q THEN co-evolutionary (MI = 0.377)
Rule 179: IF position 212 = A AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 180: IF position 212 = I AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 181: IF position 212 = A AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 182: IF position 212 = I AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 183: IF position 212 = I AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 184: IF position 212 = A AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 185: IF position 212 = I AND position 216 = Q THEN co-evolutionary (MI = 0.377)
Rule 186: IF position 212 = A AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 187: IF position 212 = I AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 188: IF position 212 = A AND position 216 = D THEN co-evolutionary (MI = 0.377)
Rule 189: IF position 212 = I AND position 216 = I THEN co-evolutionary (MI = 0.377)
Rule 190: IF position 212 = I AND position 216 = Q THEN co-evolutionary (MI = 0.377)
Rule 191: IF position 212 = I AND position 216 = Q THEN co-evolutionary (MI = 0.377)
Rule 192: IF position 212 = I AND position 216 = Q THEN co-evolutionary (MI = 0.377)
Rule 193: IF position 212 = A AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 194: IF position 212 = I AND position 216 = V THEN co-evolutionary (MI = 0.377)
Rule 195: IF position 212 = I AND position 216 = D THEN co-evolutionary (MI = 0.377)
Rule 196: IF position 212 = I AND position 216 = E THEN co-evolutionary (MI = 0.377)
Rule 197: IF position 212 = I AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 198: IF position 212 = I AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 199: IF position 212 = I AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 200: IF position 212 = D AND position 216 = D THEN co-evolutionary (MI = 0.377)
Rule 201: IF position 212 = W AND position 216 = D THEN co-evolutionary (MI = 0.377)
Rule 202: IF position 212 = W AND position 216 = Q THEN co-evolutionary (MI = 0.377)
Rule 203: IF position 212 = W AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 204: IF position 212 = D AND position 216 = N THEN co-evolutionary (MI = 0.377)
Rule 205: IF position 212 = W AND position 216 = N THEN co-evolutionary (MI = 0.377)
**Rule 206:** IF position 212 = **W** AND position 216 = **N** THEN co-evolutionary (MI = 0.377)

---

## Position Pair (210, 215)

| Property | Value |
|----------|-------|
| Mutual Information | 0.2377 |
| Total mutations | 282 |
| Reference pos 210 | N |
| Reference pos 215 | G |
| On-set cells | 4 |
| Off-set cells | 0 |
| Don't-care cells | 1020 |
| Prime implicants | 26 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (210, 215): Reference = (N, G)

Co-evolutionary residue pairs (on-set):
  I-P, N-V, N-E, K-G

Don't-care positions (conserved): 1020 cells
Never-seen pairs (off-set): 0 cells
```

### Boolean Function

```
f(pos_210, pos_215) = 1 if ANY of these residue pairs appear:

    PI_1: pos_210=A AND pos_215=H
    PI_2: pos_210=F AND pos_215=E
    PI_3: pos_210=M AND pos_215=E
    PI_4: pos_210=A AND pos_215=H
    PI_5: pos_210=A AND pos_215=H
    PI_6: pos_210=A AND pos_215=H
    PI_7: pos_210=M AND pos_215=Q
    PI_8: pos_210=I AND pos_215=Q
    PI_9: pos_210=A AND pos_215=H
    PI_10: pos_210=F AND pos_215=L
    PI_11: pos_210=A AND pos_215=M
    PI_12: pos_210=A AND pos_215=E
    PI_13: pos_210=A AND pos_215=H
    PI_14: pos_210=I AND pos_215=E
    PI_15: pos_210=F AND pos_215=A
    PI_16: pos_210=F AND pos_215=Q
    PI_17: pos_210=F AND pos_215=Q
    PI_18: pos_210=F AND pos_215=Q
    PI_19: pos_210=F AND pos_215=Q
    PI_20: pos_210=F AND pos_215=E
    PI_21: pos_210=F AND pos_215=E
    PI_22: pos_210=Y AND pos_215=H
    PI_23: pos_210=Y AND pos_215=M
    PI_24: pos_210=Y AND pos_215=H
    PI_25: pos_210=Y AND pos_215=H
    PI_26: pos_210=Y AND pos_215=H

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| I-P | 1.5315 | 0.2162 | co-evolutionary |
| N-V | 0.2446 | 0.0064 | co-evolutionary |
| N-G | 0.2436 | 0.7758 | co-evolutionary |

### Inference Rules (Natural Language)

Rule 207: IF position 210 = A AND position 215 = H THEN co-evolutionary (MI = 0.238)
Rule 208: IF position 210 = F AND position 215 = E THEN co-evolutionary (MI = 0.238)
Rule 209: IF position 210 = M AND position 215 = E THEN co-evolutionary (MI = 0.238)
Rule 210: IF position 210 = A AND position 215 = H THEN co-evolutionary (MI = 0.238)
Rule 211: IF position 210 = A AND position 215 = H THEN co-evolutionary (MI = 0.238)
Rule 212: IF position 210 = A AND position 215 = H THEN co-evolutionary (MI = 0.238)
Rule 213: IF position 210 = M AND position 215 = Q THEN co-evolutionary (MI = 0.238)
Rule 214: IF position 210 = I AND position 215 = Q THEN co-evolutionary (MI = 0.238)
Rule 215: IF position 210 = A AND position 215 = H THEN co-evolutionary (MI = 0.238)
Rule 216: IF position 210 = F AND position 215 = L THEN co-evolutionary (MI = 0.238)
Rule 217: IF position 210 = A AND position 215 = M THEN co-evolutionary (MI = 0.238)
Rule 218: IF position 210 = A AND position 215 = E THEN co-evolutionary (MI = 0.238)
Rule 219: IF position 210 = A AND position 215 = H THEN co-evolutionary (MI = 0.238)
Rule 220: IF position 210 = I AND position 215 = E THEN co-evolutionary (MI = 0.238)
Rule 221: IF position 210 = F AND position 215 = A THEN co-evolutionary (MI = 0.238)
Rule 222: IF position 210 = F AND position 215 = Q THEN co-evolutionary (MI = 0.238)
Rule 223: IF position 210 = F AND position 215 = Q THEN co-evolutionary (MI = 0.238)
Rule 224: IF position 210 = F AND position 215 = Q THEN co-evolutionary (MI = 0.238)
Rule 225: IF position 210 = F AND position 215 = Q THEN co-evolutionary (MI = 0.238)
Rule 226: IF position 210 = F AND position 215 = E THEN co-evolutionary (MI = 0.238)
Rule 227: IF position 210 = F AND position 215 = E THEN co-evolutionary (MI = 0.238)
Rule 228: IF position 210 = Y AND position 215 = H THEN co-evolutionary (MI = 0.238)
Rule 229: IF position 210 = Y AND position 215 = M THEN co-evolutionary (MI = 0.238)
Rule 230: IF position 210 = Y AND position 215 = H THEN co-evolutionary (MI = 0.238)
Rule 231: IF position 210 = Y AND position 215 = H THEN co-evolutionary (MI = 0.238)
Rule 232: IF position 210 = Y AND position 215 = H THEN co-evolutionary (MI = 0.238)

---

## Position Pair (210, 212)

| Property | Value |
|----------|-------|
| Mutual Information | 0.1769 |
| Total mutations | 301 |
| Reference pos 210 | N |
| Reference pos 212 | V |
| On-set cells | 4 |
| Off-set cells | 0 |
| Don't-care cells | 1020 |
| Prime implicants | 23 |
| Essential PIs | 1 |

### K-map (Compact View)

```
Position pair (210, 212): Reference = (N, V)

Co-evolutionary residue pairs (on-set):
  I-V, I-G, N-L, N-S

Don't-care positions (conserved): 1020 cells
Never-seen pairs (off-set): 0 cells
```

### Boolean Function

```
f(pos_210, pos_212) = 1 if ANY of these residue pairs appear:

    PI_1: pos_210=A AND pos_212=E
    PI_2: pos_210=F AND pos_212=E
    PI_3: pos_210=A AND pos_212=A
    PI_4: pos_210=A AND pos_212=E
    PI_5: pos_210=A AND pos_212=E
    PI_6: pos_210=A AND pos_212=E
    PI_7: pos_210=A AND pos_212=E
    PI_8: pos_210=A AND pos_212=E
    PI_9: pos_210=A AND pos_212=E
    PI_10: pos_210=A AND pos_212=E
    PI_11: pos_210=A AND pos_212=E
    PI_12: pos_210=F AND pos_212=A
    PI_13: pos_210=F AND pos_212=E
    PI_14: pos_210=F AND pos_212=E
    PI_15: pos_210=F AND pos_212=E
    PI_16: pos_210=F AND pos_212=E
    PI_17: pos_210=I AND pos_212=N
    PI_18: pos_210=M AND pos_212=N
    PI_19: pos_210=F AND pos_212=V
    PI_20: pos_210=F AND pos_212=N
  * PI_21: pos_210=F AND pos_212=N
    PI_22: pos_210=F AND pos_212=N
    PI_23: pos_210=F AND pos_212=N

(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|
| N-L | 3.6277 | 0.0199 | co-evolutionary |
| N-S | 3.6277 | 0.0066 | co-evolutionary |
| I-V | 0.0269 | 0.9701 | co-evolutionary |

### Inference Rules (Natural Language)

Rule 233: IF position 210 = A AND position 212 = E THEN co-evolutionary (MI = 0.177)
Rule 234: IF position 210 = F AND position 212 = E THEN co-evolutionary (MI = 0.177)
Rule 235: IF position 210 = A AND position 212 = A THEN co-evolutionary (MI = 0.177)
Rule 236: IF position 210 = A AND position 212 = E THEN co-evolutionary (MI = 0.177)
Rule 237: IF position 210 = A AND position 212 = E THEN co-evolutionary (MI = 0.177)
Rule 238: IF position 210 = A AND position 212 = E THEN co-evolutionary (MI = 0.177)
Rule 239: IF position 210 = A AND position 212 = E THEN co-evolutionary (MI = 0.177)
Rule 240: IF position 210 = A AND position 212 = E THEN co-evolutionary (MI = 0.177)
Rule 241: IF position 210 = A AND position 212 = E THEN co-evolutionary (MI = 0.177)
Rule 242: IF position 210 = A AND position 212 = E THEN co-evolutionary (MI = 0.177)
Rule 243: IF position 210 = A AND position 212 = E THEN co-evolutionary (MI = 0.177)
Rule 244: IF position 210 = F AND position 212 = A THEN co-evolutionary (MI = 0.177)
Rule 245: IF position 210 = F AND position 212 = E THEN co-evolutionary (MI = 0.177)
Rule 246: IF position 210 = F AND position 212 = E THEN co-evolutionary (MI = 0.177)
Rule 247: IF position 210 = F AND position 212 = E THEN co-evolutionary (MI = 0.177)
Rule 248: IF position 210 = F AND position 212 = E THEN co-evolutionary (MI = 0.177)
Rule 249: IF position 210 = I AND position 212 = N THEN co-evolutionary (MI = 0.177)
Rule 250: IF position 210 = M AND position 212 = N THEN co-evolutionary (MI = 0.177)
Rule 251: IF position 210 = F AND position 212 = V THEN co-evolutionary (MI = 0.177)
Rule 252: IF position 210 = F AND position 212 = N THEN co-evolutionary (MI = 0.177)
**Rule 253:** IF position 210 = **F** AND position 212 = **N** THEN co-evolutionary (MI = 0.177)
Rule 254: IF position 210 = F AND position 212 = N THEN co-evolutionary (MI = 0.177)
Rule 255: IF position 210 = F AND position 212 = N THEN co-evolutionary (MI = 0.177)

---

## Summary

| Metric | Value |
|--------|-------|
| Sequences | 1299 |
| Variable positions | 57 |
| Co-evolutionary pairs | 10 |
| Total inference rules | 255 |
| Position pairs with rules | 15 |

## How to Apply

1. Extract residues at positions 68-79 from a new sequence
2. For each position pair, check if the residue pair matches any rule
3. If YES: that position pair is co-evolutionary
4. If position i mutates: find which residue at position j satisfies the co-evolutionary constraint

**Example:** If position 76 mutates to Y, check rules for position 76.
Rule 2 says: IF pos 76 = Y AND pos 77 = K THEN co-evolutionary.
So position 77 must also mutate to K.