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
| On-set cells | 0 |
| Off-set cells | 400 |
| Don't-care cells | 0 |
| Prime implicants | 0 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (413, 427): Reference = (N, G)

Co-evolutionary residue pairs (on-set):
  

Don't-care positions (conserved): 0 cells
Never-seen pairs (off-set): 400 cells
```

### Boolean Function

```
f(pos_413, pos_427) = 1 if ANY of these residue pairs appear:


(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|

### Inference Rules (Natural Language)


---

## Position Pair (413, 425)

| Property | Value |
|----------|-------|
| Mutual Information | 2.3345 |
| Total mutations | 271 |
| Reference pos 413 | N |
| Reference pos 425 | F |
| On-set cells | 0 |
| Off-set cells | 400 |
| Don't-care cells | 0 |
| Prime implicants | 0 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (413, 425): Reference = (N, F)

Co-evolutionary residue pairs (on-set):
  

Don't-care positions (conserved): 0 cells
Never-seen pairs (off-set): 400 cells
```

### Boolean Function

```
f(pos_413, pos_425) = 1 if ANY of these residue pairs appear:


(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|

### Inference Rules (Natural Language)


---

## Position Pair (413, 426)

| Property | Value |
|----------|-------|
| Mutual Information | 2.3199 |
| Total mutations | 271 |
| Reference pos 413 | N |
| Reference pos 426 | T |
| On-set cells | 0 |
| Off-set cells | 400 |
| Don't-care cells | 0 |
| Prime implicants | 0 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (413, 426): Reference = (N, T)

Co-evolutionary residue pairs (on-set):
  

Don't-care positions (conserved): 0 cells
Never-seen pairs (off-set): 400 cells
```

### Boolean Function

```
f(pos_413, pos_426) = 1 if ANY of these residue pairs appear:


(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|

### Inference Rules (Natural Language)


---

## Position Pair (413, 424)

| Property | Value |
|----------|-------|
| Mutual Information | 2.3104 |
| Total mutations | 271 |
| Reference pos 413 | N |
| Reference pos 424 | D |
| On-set cells | 0 |
| Off-set cells | 400 |
| Don't-care cells | 0 |
| Prime implicants | 0 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (413, 424): Reference = (N, D)

Co-evolutionary residue pairs (on-set):
  

Don't-care positions (conserved): 0 cells
Never-seen pairs (off-set): 400 cells
```

### Boolean Function

```
f(pos_413, pos_424) = 1 if ANY of these residue pairs appear:


(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|

### Inference Rules (Natural Language)


---

## Position Pair (1026, 1040)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2989 |
| Total mutations | 278 |
| Reference pos 1026 | S |
| Reference pos 1040 | G |
| On-set cells | 0 |
| Off-set cells | 400 |
| Don't-care cells | 0 |
| Prime implicants | 0 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (1026, 1040): Reference = (S, G)

Co-evolutionary residue pairs (on-set):
  

Don't-care positions (conserved): 0 cells
Never-seen pairs (off-set): 400 cells
```

### Boolean Function

```
f(pos_1026, pos_1040) = 1 if ANY of these residue pairs appear:


(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|

### Inference Rules (Natural Language)


---

## Position Pair (1026, 1042)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2989 |
| Total mutations | 278 |
| Reference pos 1026 | S |
| Reference pos 1042 | G |
| On-set cells | 0 |
| Off-set cells | 400 |
| Don't-care cells | 0 |
| Prime implicants | 0 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (1026, 1042): Reference = (S, G)

Co-evolutionary residue pairs (on-set):
  

Don't-care positions (conserved): 0 cells
Never-seen pairs (off-set): 400 cells
```

### Boolean Function

```
f(pos_1026, pos_1042) = 1 if ANY of these residue pairs appear:


(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|

### Inference Rules (Natural Language)


---

## Position Pair (1040, 1042)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2989 |
| Total mutations | 278 |
| Reference pos 1040 | G |
| Reference pos 1042 | G |
| On-set cells | 0 |
| Off-set cells | 400 |
| Don't-care cells | 0 |
| Prime implicants | 0 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (1040, 1042): Reference = (G, G)

Co-evolutionary residue pairs (on-set):
  

Don't-care positions (conserved): 0 cells
Never-seen pairs (off-set): 400 cells
```

### Boolean Function

```
f(pos_1040, pos_1042) = 1 if ANY of these residue pairs appear:


(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|

### Inference Rules (Natural Language)


---

## Position Pair (462, 473)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2946 |
| Total mutations | 265 |
| Reference pos 462 | R |
| Reference pos 473 | N |
| On-set cells | 0 |
| Off-set cells | 400 |
| Don't-care cells | 0 |
| Prime implicants | 0 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (462, 473): Reference = (R, N)

Co-evolutionary residue pairs (on-set):
  

Don't-care positions (conserved): 0 cells
Never-seen pairs (off-set): 400 cells
```

### Boolean Function

```
f(pos_462, pos_473) = 1 if ANY of these residue pairs appear:


(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|

### Inference Rules (Natural Language)


---

## Position Pair (468, 473)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2946 |
| Total mutations | 265 |
| Reference pos 468 | I |
| Reference pos 473 | N |
| On-set cells | 0 |
| Off-set cells | 400 |
| Don't-care cells | 0 |
| Prime implicants | 0 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (468, 473): Reference = (I, N)

Co-evolutionary residue pairs (on-set):
  

Don't-care positions (conserved): 0 cells
Never-seen pairs (off-set): 400 cells
```

### Boolean Function

```
f(pos_468, pos_473) = 1 if ANY of these residue pairs appear:


(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|

### Inference Rules (Natural Language)


---

## Position Pair (1064, 1074)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2888 |
| Total mutations | 278 |
| Reference pos 1064 | V |
| Reference pos 1074 | A |
| On-set cells | 0 |
| Off-set cells | 400 |
| Don't-care cells | 0 |
| Prime implicants | 0 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (1064, 1074): Reference = (V, A)

Co-evolutionary residue pairs (on-set):
  

Don't-care positions (conserved): 0 cells
Never-seen pairs (off-set): 400 cells
```

### Boolean Function

```
f(pos_1064, pos_1074) = 1 if ANY of these residue pairs appear:


(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|

### Inference Rules (Natural Language)


---

## Position Pair (1064, 1065)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2876 |
| Total mutations | 278 |
| Reference pos 1064 | V |
| Reference pos 1065 | P |
| On-set cells | 0 |
| Off-set cells | 400 |
| Don't-care cells | 0 |
| Prime implicants | 0 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (1064, 1065): Reference = (V, P)

Co-evolutionary residue pairs (on-set):
  

Don't-care positions (conserved): 0 cells
Never-seen pairs (off-set): 400 cells
```

### Boolean Function

```
f(pos_1064, pos_1065) = 1 if ANY of these residue pairs appear:


(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|

### Inference Rules (Natural Language)


---

## Position Pair (413, 428)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2865 |
| Total mutations | 271 |
| Reference pos 413 | N |
| Reference pos 428 | C |
| On-set cells | 0 |
| Off-set cells | 400 |
| Don't-care cells | 0 |
| Prime implicants | 0 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (413, 428): Reference = (N, C)

Co-evolutionary residue pairs (on-set):
  

Don't-care positions (conserved): 0 cells
Never-seen pairs (off-set): 400 cells
```

### Boolean Function

```
f(pos_413, pos_428) = 1 if ANY of these residue pairs appear:


(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|

### Inference Rules (Natural Language)


---

## Position Pair (459, 473)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2842 |
| Total mutations | 265 |
| Reference pos 459 | P |
| Reference pos 473 | N |
| On-set cells | 0 |
| Off-set cells | 400 |
| Don't-care cells | 0 |
| Prime implicants | 0 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (459, 473): Reference = (P, N)

Co-evolutionary residue pairs (on-set):
  

Don't-care positions (conserved): 0 cells
Never-seen pairs (off-set): 400 cells
```

### Boolean Function

```
f(pos_459, pos_473) = 1 if ANY of these residue pairs appear:


(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|

### Inference Rules (Natural Language)


---

## Position Pair (469, 473)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2842 |
| Total mutations | 265 |
| Reference pos 469 | Y |
| Reference pos 473 | N |
| On-set cells | 0 |
| Off-set cells | 400 |
| Don't-care cells | 0 |
| Prime implicants | 0 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (469, 473): Reference = (Y, N)

Co-evolutionary residue pairs (on-set):
  

Don't-care positions (conserved): 0 cells
Never-seen pairs (off-set): 400 cells
```

### Boolean Function

```
f(pos_469, pos_473) = 1 if ANY of these residue pairs appear:


(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|

### Inference Rules (Natural Language)


---

## Position Pair (1064, 1066)

| Property | Value |
|----------|-------|
| Mutual Information | 2.2816 |
| Total mutations | 278 |
| Reference pos 1064 | V |
| Reference pos 1066 | A |
| On-set cells | 0 |
| Off-set cells | 400 |
| Don't-care cells | 0 |
| Prime implicants | 0 |
| Essential PIs | 0 |

### K-map (Compact View)

```
Position pair (1064, 1066): Reference = (V, A)

Co-evolutionary residue pairs (on-set):
  

Don't-care positions (conserved): 0 cells
Never-seen pairs (off-set): 400 cells
```

### Boolean Function

```
f(pos_1064, pos_1066) = 1 if ANY of these residue pairs appear:


(* = essential prime implicant)
```

### Coupling Constants (J_ij)

| Residue Pair | J_ij | Frequency | Type |
|-------------|------|-----------|------|

### Inference Rules (Natural Language)


---

## Summary

| Metric | Value |
|--------|-------|
| Sequences | 1299 |
| Variable positions | 57 |
| Co-evolutionary pairs | 36918 |
| Total inference rules | 0 |
| Position pairs with rules | 15 |

## How to Apply

1. Extract residues at positions 68-79 from a new sequence
2. For each position pair, check if the residue pair matches any rule
3. If YES: that position pair is co-evolutionary
4. If position i mutates: find which residue at position j satisfies the co-evolutionary constraint

**Example:** If position 76 mutates to Y, check rules for position 76.
Rule 2 says: IF pos 76 = Y AND pos 77 = K THEN co-evolutionary.
So position 77 must also mutate to K.