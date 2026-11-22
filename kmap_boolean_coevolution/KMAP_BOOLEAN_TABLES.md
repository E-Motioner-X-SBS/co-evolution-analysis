# Co-evolution Boolean Functions for SARS-CoV-2 Spike Protein

**Dataset:** 1299 sequences
**Encoding:** Base-20 (He 2012)
**Variable positions:** 21 / 1276
**Co-evolutionary pairs:** 10

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
| Essential PIs | 2 |

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| I | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| L | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| V | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| M | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| F | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Y | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| W | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| E | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| D | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Q | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** |
| N | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| H | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| K | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| R | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | DC |
| S | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| T | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| C | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| P | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| G | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

### Boolean Function

$$f(s_4,s_3,s_2,s_1,s_0,\; t_4,t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|
| 1 | $s3 \cdot s1 \cdot \bar{s0} \cdot t4 \cdot t1 \cdot t0$ | (QR, G) | 4 | 1 |
| 2 | $s3 \cdot s2 \cdot s1 \cdot \bar{s0} \cdot t3 \cdot t2 \cdot t1 \cdot t0$ | (R, S) | 2 | 1 |

### Inference Rules (Natural Language)

**Rule 1:** IF position 495 in **{QR}** AND position 498 in **{G}** THEN observed mutation (MI = 0.871)
**Rule 2:** IF position 495 in **{R}** AND position 498 in **{S}** THEN observed mutation (MI = 0.871)

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
| Essential PIs | 3 |

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| I | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| L | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| V | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| M | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| F | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Y | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| W | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| E | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| D | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Q | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| N | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| H | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| K | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| R | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| S | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 |
| T | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| C | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| P | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| G | 0 | 0 | DC | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 |

### Boolean Function

$$f(s_4,s_3,s_2,s_1,s_0,\; t_4,t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|
| 1 | $s3 \cdot s2 \cdot s1 \cdot s0 \cdot \bar{t4} \cdot \bar{t3} \cdot \bar{t2} \cdot t1 \cdot \bar{t0}$ | (S, L) | 1 | 1 |
| 2 | $s3 \cdot s2 \cdot s1 \cdot s0 \cdot t3 \cdot t2 \cdot t1 \cdot \bar{t0}$ | (S, R) | 2 | 1 |
| 3 | $s4 \cdot s1 \cdot s0 \cdot t3 \cdot t2 \cdot t1 \cdot \bar{t0}$ | (G, R) | 3 | 1 |

### Inference Rules (Natural Language)

**Rule 3:** IF position 448 in **{S}** AND position 454 in **{L}** THEN observed mutation (MI = 0.834)
**Rule 4:** IF position 448 in **{S}** AND position 454 in **{R}** THEN observed mutation (MI = 0.834)
**Rule 5:** IF position 448 in **{G}** AND position 454 in **{R}** THEN observed mutation (MI = 0.834)

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
| Essential PIs | 3 |

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| I | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| L | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| V | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** |
| M | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| F | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | DC |
| Y | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| W | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| E | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| D | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Q | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| N | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| H | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| K | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| R | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| S | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| T | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| C | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| P | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** |
| G | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

### Boolean Function

$$f(s_4,s_3,s_2,s_1,s_0,\; t_4,t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|
| 1 | $\bar{s4} \cdot \bar{s3} \cdot \bar{s2} \cdot s1 \cdot s0 \cdot t4 \cdot t1 \cdot t0$ | (V, G) | 2 | 1 |
| 2 | $\bar{s3} \cdot s2 \cdot \bar{s1} \cdot s0 \cdot t3 \cdot t2 \cdot t1 \cdot t0$ | (F, S) | 2 | 1 |
| 3 | $s4 \cdot s1 \cdot \bar{s0} \cdot t4 \cdot t1 \cdot t0$ | (P, G) | 4 | 1 |

### Inference Rules (Natural Language)

**Rule 6:** IF position 488 in **{V}** AND position 498 in **{G}** THEN observed mutation (MI = 0.822)
**Rule 7:** IF position 488 in **{F}** AND position 498 in **{S}** THEN observed mutation (MI = 0.822)
**Rule 8:** IF position 488 in **{P}** AND position 498 in **{G}** THEN observed mutation (MI = 0.822)

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
| Essential PIs | 3 |

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| I | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| L | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| V | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| M | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| F | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Y | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| W | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| E | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| D | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Q | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| N | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 |
| H | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| K | 0 | 0 | DC | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 |
| R | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| S | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| T | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| C | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| P | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| G | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

### Boolean Function

$$f(s_4,s_3,s_2,s_1,s_0,\; t_4,t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|
| 1 | $s3 \cdot \bar{s2} \cdot s1 \cdot s0 \cdot \bar{t4} \cdot \bar{t3} \cdot \bar{t2} \cdot t1 \cdot \bar{t0}$ | (N, L) | 1 | 1 |
| 2 | $s3 \cdot \bar{s2} \cdot s1 \cdot s0 \cdot t3 \cdot t2 \cdot t1 \cdot \bar{t0}$ | (N, R) | 2 | 1 |
| 3 | $s3 \cdot s2 \cdot \bar{s1} \cdot s0 \cdot t3 \cdot t2 \cdot t1 \cdot \bar{t0}$ | (K, R) | 2 | 1 |

### Inference Rules (Natural Language)

**Rule 9:** IF position 442 in **{N}** AND position 454 in **{L}** THEN observed mutation (MI = 0.811)
**Rule 10:** IF position 442 in **{N}** AND position 454 in **{R}** THEN observed mutation (MI = 0.811)
**Rule 11:** IF position 442 in **{K}** AND position 454 in **{R}** THEN observed mutation (MI = 0.811)

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
| Essential PIs | 3 |

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| I | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| L | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| V | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| M | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| F | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Y | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| W | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| E | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| D | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Q | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| N | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | **1** |
| H | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| K | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | DC |
| R | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| S | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| T | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| C | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| P | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| G | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

### Boolean Function

$$f(s_4,s_3,s_2,s_1,s_0,\; t_4,t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|
| 1 | $s3 \cdot \bar{s2} \cdot s1 \cdot s0 \cdot t3 \cdot t2 \cdot t1 \cdot t0$ | (N, S) | 2 | 1 |
| 2 | $s3 \cdot \bar{s2} \cdot s1 \cdot s0 \cdot t4 \cdot t1 \cdot t0$ | (N, G) | 3 | 1 |
| 3 | $s3 \cdot s2 \cdot \bar{s1} \cdot s0 \cdot t3 \cdot t2 \cdot t1 \cdot t0$ | (K, S) | 2 | 1 |

### Inference Rules (Natural Language)

**Rule 12:** IF position 442 in **{N}** AND position 448 in **{S}** THEN observed mutation (MI = 0.728)
**Rule 13:** IF position 442 in **{N}** AND position 448 in **{G}** THEN observed mutation (MI = 0.728)
**Rule 14:** IF position 442 in **{K}** AND position 448 in **{S}** THEN observed mutation (MI = 0.728)

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
| Essential PIs | 4 |

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| I | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| L | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| V | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 |
| M | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| F | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Y | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| W | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| E | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| D | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Q | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| N | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| H | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| K | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| R | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| S | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** |
| T | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| C | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| P | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| G | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

### Boolean Function

$$f(s_4,s_3,s_2,s_1,s_0,\; t_4,t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|
| 1 | $\bar{s4} \cdot \bar{s3} \cdot \bar{s2} \cdot \bar{s1} \cdot s0 \cdot \bar{t4} \cdot \bar{t3} \cdot \bar{t2} \cdot t1 \cdot t0$ | (I, V) | 0 | 1 |
| 2 | $\bar{s4} \cdot \bar{s3} \cdot \bar{s2} \cdot s1 \cdot \bar{s0} \cdot \bar{t4} \cdot \bar{t3} \cdot \bar{t2} \cdot t1 \cdot t0$ | (L, V) | 0 | 1 |
| 3 | $\bar{s4} \cdot \bar{s3} \cdot \bar{s2} \cdot s1 \cdot s0 \cdot t4 \cdot t1 \cdot \bar{t0}$ | (V, P) | 2 | 1 |
| 4 | $s3 \cdot s2 \cdot s1 \cdot s0 \cdot t4 \cdot t1 \cdot t0$ | (S, G) | 3 | 1 |

### Inference Rules (Natural Language)

**Rule 15:** IF position 212 in **{I}** AND position 215 in **{V}** THEN observed mutation (MI = 0.398)
**Rule 16:** IF position 212 in **{L}** AND position 215 in **{V}** THEN observed mutation (MI = 0.398)
**Rule 17:** IF position 212 in **{V}** AND position 215 in **{P}** THEN observed mutation (MI = 0.398)
**Rule 18:** IF position 212 in **{S}** AND position 215 in **{G}** THEN observed mutation (MI = 0.398)

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
| Essential PIs | 4 |

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| I | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| L | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| V | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 |
| M | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| F | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Y | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| W | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| E | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 |
| D | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Q | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| N | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| H | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| K | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| R | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| S | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| T | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| C | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| P | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 | 0 |
| G | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | DC | 0 | 0 | 0 | 0 | 0 |

### Boolean Function

$$f(s_4,s_3,s_2,s_1,s_0,\; t_4,t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|
| 1 | $\bar{s3} \cdot \bar{s2} \cdot s1 \cdot s0 \cdot t3 \cdot t2 \cdot t1 \cdot \bar{t0}$ | (VG, R) | 2 | 1 |
| 2 | $s3 \cdot \bar{s2} \cdot \bar{s1} \cdot \bar{s0} \cdot t3 \cdot t2 \cdot t1 \cdot \bar{t0}$ | (E, R) | 2 | 1 |
| 3 | $s4 \cdot s1 \cdot \bar{s0} \cdot t3 \cdot \bar{t2} \cdot \bar{t1} \cdot \bar{t0}$ | (P, E) | 3 | 1 |
| 4 | $s4 \cdot s1 \cdot \bar{s0} \cdot t3 \cdot t2 \cdot \bar{t1} \cdot t0$ | (P, K) | 3 | 1 |

### Inference Rules (Natural Language)

**Rule 19:** IF position 215 in **{VG}** AND position 216 in **{R}** THEN observed mutation (MI = 0.377)
**Rule 20:** IF position 215 in **{E}** AND position 216 in **{R}** THEN observed mutation (MI = 0.377)
**Rule 21:** IF position 215 in **{P}** AND position 216 in **{E}** THEN observed mutation (MI = 0.377)
**Rule 22:** IF position 215 in **{P}** AND position 216 in **{K}** THEN observed mutation (MI = 0.377)

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
| Essential PIs | 6 |

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| I | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 |
| L | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 |
| V | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 | 0 |
| M | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| F | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Y | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| W | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| E | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| D | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Q | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| N | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| H | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| K | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| R | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| S | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 |
| T | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| C | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| P | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| G | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 |

### Boolean Function

$$f(s_4,s_3,s_2,s_1,s_0,\; t_4,t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|
| 1 | $\bar{s4} \cdot \bar{s3} \cdot \bar{s2} \cdot \bar{s1} \cdot s0 \cdot t3 \cdot t2 \cdot t1 \cdot \bar{t0}$ | (I, R) | 1 | 1 |
| 2 | $\bar{s4} \cdot \bar{s3} \cdot \bar{s2} \cdot s1 \cdot \bar{s0} \cdot t3 \cdot t2 \cdot t1 \cdot \bar{t0}$ | (L, R) | 1 | 1 |
| 3 | $\bar{s4} \cdot \bar{s3} \cdot \bar{s2} \cdot s1 \cdot s0 \cdot t3 \cdot \bar{t2} \cdot \bar{t1} \cdot \bar{t0}$ | (V, E) | 1 | 1 |
| 4 | $\bar{s4} \cdot \bar{s3} \cdot \bar{s2} \cdot s1 \cdot s0 \cdot t3 \cdot t2 \cdot \bar{t1} \cdot t0$ | (V, K) | 1 | 1 |
| 5 | $s3 \cdot s2 \cdot s1 \cdot s0 \cdot t3 \cdot t2 \cdot t1 \cdot \bar{t0}$ | (S, R) | 2 | 1 |
| 6 | $s4 \cdot s1 \cdot s0 \cdot t3 \cdot t2 \cdot t1 \cdot \bar{t0}$ | (G, R) | 3 | 1 |

### Inference Rules (Natural Language)

**Rule 23:** IF position 212 in **{I}** AND position 216 in **{R}** THEN observed mutation (MI = 0.377)
**Rule 24:** IF position 212 in **{L}** AND position 216 in **{R}** THEN observed mutation (MI = 0.377)
**Rule 25:** IF position 212 in **{V}** AND position 216 in **{E}** THEN observed mutation (MI = 0.377)
**Rule 26:** IF position 212 in **{V}** AND position 216 in **{K}** THEN observed mutation (MI = 0.377)
**Rule 27:** IF position 212 in **{S}** AND position 216 in **{R}** THEN observed mutation (MI = 0.377)
**Rule 28:** IF position 212 in **{G}** AND position 216 in **{R}** THEN observed mutation (MI = 0.377)

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
| Essential PIs | 4 |

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| I | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 |
| L | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| V | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| M | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| F | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Y | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| W | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| E | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| D | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Q | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| N | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | DC |
| H | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| K | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** |
| R | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| S | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| T | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| C | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| P | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| G | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

### Boolean Function

$$f(s_4,s_3,s_2,s_1,s_0,\; t_4,t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|
| 1 | $\bar{s4} \cdot \bar{s3} \cdot \bar{s2} \cdot \bar{s1} \cdot s0 \cdot t4 \cdot t1 \cdot \bar{t0}$ | (I, P) | 2 | 1 |
| 2 | $s3 \cdot \bar{s2} \cdot s1 \cdot s0 \cdot \bar{t3} \cdot \bar{t2} \cdot t1 \cdot t0$ | (N, VG) | 2 | 1 |
| 3 | $s3 \cdot \bar{s2} \cdot s1 \cdot s0 \cdot t3 \cdot \bar{t2} \cdot \bar{t1} \cdot \bar{t0}$ | (N, E) | 2 | 1 |
| 4 | $s3 \cdot s2 \cdot \bar{s1} \cdot s0 \cdot t4 \cdot t1 \cdot t0$ | (K, G) | 3 | 1 |

### Inference Rules (Natural Language)

**Rule 29:** IF position 210 in **{I}** AND position 215 in **{P}** THEN observed mutation (MI = 0.238)
**Rule 30:** IF position 210 in **{N}** AND position 215 in **{VG}** THEN observed mutation (MI = 0.238)
**Rule 31:** IF position 210 in **{N}** AND position 215 in **{E}** THEN observed mutation (MI = 0.238)
**Rule 32:** IF position 210 in **{K}** AND position 215 in **{G}** THEN observed mutation (MI = 0.238)

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
| Essential PIs | 3 |

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| I | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** |
| L | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| V | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| M | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| F | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Y | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| W | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| E | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| D | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Q | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| N | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **1** | 0 | 0 | 0 | 0 |
| H | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| K | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| R | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| S | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| T | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| C | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| P | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| G | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

### Boolean Function

$$f(s_4,s_3,s_2,s_1,s_0,\; t_4,t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|
| 1 | $\bar{s4} \cdot \bar{s3} \cdot \bar{s2} \cdot \bar{s1} \cdot s0 \cdot \bar{t3} \cdot \bar{t2} \cdot t1 \cdot t0$ | (I, VG) | 1 | 2 |
| 2 | $s3 \cdot \bar{s2} \cdot s1 \cdot s0 \cdot \bar{t4} \cdot \bar{t3} \cdot \bar{t2} \cdot t1 \cdot \bar{t0}$ | (N, L) | 1 | 1 |
| 3 | $s3 \cdot \bar{s2} \cdot s1 \cdot s0 \cdot t3 \cdot t2 \cdot t1 \cdot t0$ | (N, S) | 2 | 1 |

### Inference Rules (Natural Language)

**Rule 33:** IF position 210 in **{I}** AND position 212 in **{VG}** THEN observed mutation (MI = 0.177)
**Rule 34:** IF position 210 in **{N}** AND position 212 in **{L}** THEN observed mutation (MI = 0.177)
**Rule 35:** IF position 210 in **{N}** AND position 212 in **{S}** THEN observed mutation (MI = 0.177)

---

## Summary

| Metric | Value |
|--------|-------|
| Sequences | 1299 |
| Variable positions | 21 |
| Co-evolutionary pairs | 10 |
| Total inference rules | 35 |
| Position pairs with rules | 10 |

## How to Apply

For a new sequence, extract the residues at each co-evolving pair and check:
```
For each (pos_i, pos_j) pair:
  1. Get residues aa_i at pos_i and aa_j at pos_j
  2. Check if (aa_i, aa_j) matches any essential PI
  3. If YES → position pair is co-evolutionary
  4. If position i mutates → find which aa_j satisfies the PI
```