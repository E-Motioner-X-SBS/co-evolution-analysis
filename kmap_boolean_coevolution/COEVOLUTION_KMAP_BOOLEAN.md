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
| Off-set cells | 0 |
| Don't-care cells | 1022 |
| Prime implicants | 20 |
| Essential PIs | 0 |

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| I | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| L | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| V | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| M | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| F | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Y | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| W | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| E | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| D | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Q | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** |
| N | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| H | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| K | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| R | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC |
| S | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| T | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| C | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| P | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| G | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |

### Boolean Function

$$f(s_3,s_2,s_1,s_0,\; t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|

### Inference Rules (Natural Language)


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

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| I | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| L | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| V | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| M | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| F | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Y | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| W | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| E | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| D | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Q | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| N | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| H | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| K | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| R | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| S | DC | DC | **1** | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC | DC |
| T | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| C | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| P | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| G | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC | DC |

### Boolean Function

$$f(s_3,s_2,s_1,s_0,\; t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|

### Inference Rules (Natural Language)


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

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| I | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| L | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| V | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** |
| M | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| F | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC |
| Y | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| W | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| E | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| D | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Q | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| N | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| H | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| K | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| R | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| S | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| T | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| C | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| P | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** |
| G | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |

### Boolean Function

$$f(s_3,s_2,s_1,s_0,\; t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|

### Inference Rules (Natural Language)


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

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| I | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| L | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| V | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| M | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| F | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Y | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| W | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| E | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| D | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Q | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| N | DC | DC | **1** | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC | DC |
| H | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| K | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC | DC |
| R | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| S | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| T | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| C | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| P | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| G | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |

### Boolean Function

$$f(s_3,s_2,s_1,s_0,\; t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|

### Inference Rules (Natural Language)


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

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| I | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| L | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| V | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| M | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| F | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Y | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| W | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| E | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| D | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Q | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| N | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | **1** |
| H | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| K | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC |
| R | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| S | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| T | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| C | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| P | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| G | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |

### Boolean Function

$$f(s_3,s_2,s_1,s_0,\; t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|

### Inference Rules (Natural Language)


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

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| I | DC | DC | DC | **1** | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| L | DC | DC | DC | **1** | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| V | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC |
| M | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| F | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Y | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| W | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| E | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| D | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Q | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| N | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| H | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| K | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| R | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| S | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** |
| T | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| C | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| P | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| G | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |

### Boolean Function

$$f(s_3,s_2,s_1,s_0,\; t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|
| 1 | $\bar{s3} \cdot s0 \cdot t3 \cdot t2 \cdot \bar{t1} \cdot \bar{t0}$ | (V, P) | 3 | 8 |

### Inference Rules (Natural Language)

**Rule 1:** IF position 212 = **V** AND position 215 = **P** THEN co-evolutionary (MI = 0.398)

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

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| I | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| L | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| V | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC | DC |
| M | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| F | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Y | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| W | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| E | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC | DC |
| D | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Q | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| N | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| H | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| K | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| R | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| S | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| T | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| C | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| P | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC | **1** | DC | DC | DC | DC | DC | DC |
| G | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |

### Boolean Function

$$f(s_3,s_2,s_1,s_0,\; t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|

### Inference Rules (Natural Language)


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

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| I | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC | DC |
| L | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC | DC |
| V | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC | **1** | DC | DC | DC | DC | DC | DC |
| M | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| F | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Y | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| W | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| E | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| D | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Q | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| N | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| H | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| K | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| R | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| S | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC | DC |
| T | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| C | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| P | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| G | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC | DC |

### Boolean Function

$$f(s_3,s_2,s_1,s_0,\; t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|
| 1 | $\bar{s3} \cdot s2 \cdot s1 \cdot s0 \cdot t3 \cdot t1 \cdot t0$ | (S, R) | 1 | 2 |

### Inference Rules (Natural Language)

**Rule 2:** IF position 212 = **S** AND position 216 = **R** THEN co-evolutionary (MI = 0.377)

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

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| I | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC |
| L | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| V | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| M | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| F | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Y | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| W | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| E | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| D | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Q | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| N | DC | DC | DC | **1** | DC | DC | DC | DC | **1** | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| H | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| K | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** |
| R | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| S | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| T | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| C | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| P | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| G | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |

### Boolean Function

$$f(s_3,s_2,s_1,s_0,\; t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|

### Inference Rules (Natural Language)


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

### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)

| AA_i \ AA_j | A | I | L | V | M | F | Y | W | E | D | Q | N | H | K | R | S | T | C | P | G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| I | DC | DC | DC | **1** | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** |
| L | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| V | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| M | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| F | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Y | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| W | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| E | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| D | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| Q | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| N | DC | DC | **1** | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | **1** | DC | DC | DC | DC |
| H | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| K | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| R | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| S | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| T | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| C | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| P | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |
| G | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC | DC |

### Boolean Function

$$f(s_3,s_2,s_1,s_0,\; t_3,t_2,t_1,t_0) = 
\text{OR of all prime implicants below}$$

### Essential Prime Implicants (minimum covering set)

| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |
|---|-------------------|-------------|-------------|----------|
| 1 | $\bar{s3} \cdot s2 \cdot \bar{s1} \cdot s0 \cdot t3 \cdot t1 \cdot t0$ | (N, S) | 1 | 2 |

### Inference Rules (Natural Language)

**Rule 3:** IF position 210 = **N** AND position 212 = **S** THEN co-evolutionary (MI = 0.177)

---

## Summary

| Metric | Value |
|--------|-------|
| Sequences | 1299 |
| Variable positions | 21 |
| Co-evolutionary pairs | 10 |
| Total inference rules | 3 |
| Position pairs with rules | 3 |

## How to Apply

For a new sequence, extract residues at positions 68-79 and check:
```
For each (pos_i, pos_j) pair:
  1. Get residues aa_i at pos_i and aa_j at pos_j
  2. Check if (aa_i, aa_j) matches any essential PI
  3. If YES → position pair is co-evolutionary
  4. If position i mutates → find which aa_j satisfies the PI
```