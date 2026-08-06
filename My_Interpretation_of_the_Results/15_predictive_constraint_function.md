# 15. predictive_constraint_function.py - Can We Predict Mutations?


**Sequence length analyzed: 1,276 positions (full length), all 1,299 sequences (800 train, 499 test).**

## What the Program Does

This script tests whether the constraint function can predict the co-evolving partner of a mutation. It implements the continuous (frequency-based) approach:

1. Find co-evolving pairs (MI > 0.1, window 30).
2. For the top 10 pairs, build the frequency K-map.
3. Compute the constraint function C(a, b) = ln(P(a,b) / (P(a) * P(b))).
4. Convert to probability with the sigmoid: P_co-evolution = 1 / (1 + e^-C).
5. Train on the first 800 sequences, test on the remaining 499.
6. For each test mutation at position i, predict the partner at j as the argmax of C.

## The Constraint Function

```
C(a, b) = ln( P(a, b) / (P(a) * P(b)) )
```

- C > 0: pair more common than expected (co-evolutionary).
- C < 0: pair less common than expected (anti-correlated).
- C = 0: independent.

## The Sigmoid Prediction

```
P_co-evolution(a, b) = sigma(C) = 1 / (1 + e^-C)
```

This maps C to a probability between 0 and 1. C > 0 gives P > 0.5.

## The Prediction Protocol

For each test sequence (indices 800 to 1299):

1. Take the residue at position i (a mutation relative to the reference).
2. Look up the constraint row C(a, :).
3. Predict the partner as the argmax: best_j = argmax over b of C(a, b).
4. If best_j equals the actual residue at position j in the test sequence, the prediction is correct.

## Worked Example: Pair (413, 427)

The constraint row for the mutation at 413 = W gives the highest C value at 427 = E (the co-evolving partner). If the test sequence actually has E at 427, the prediction is correct.

The top co-evolutionary pairs found include (413, 427) with C values up to 4.68 (e.g., K-H at (1064, 1074)), corresponding to P = 0.99.

## Results

| Metric | Value |
|--------|-------|
| Training sequences | 800 |
| Test sequences | 499 |
| Test mutations | 1,062 |
| Correct predictions | 62 |
| **Prediction accuracy** | **5.84%** |
| Strongest predicted pairs | P = 0.99 (e.g., K-H), P = 0.98 (e.g., E-C) |

## Inference

5.84% accuracy is low, and this is an honest and important result. The constraint function captures the STRUCTURE of co-evolution (which pairs are favored) but cannot predict specific mutations in held-out sequences. The reasons:

1. **Probabilistic co-evolution.** Even strongly coupled positions have several possible partners; argmax picks one, but evolution picks from a distribution.
2. **Lineage specificity.** The 1,299 sequences span many Omicron sub-lineages with different reference residues. Rules learned from one lineage do not transfer to another.
3. **The argmax is a point estimate.** The full distribution matters; the mode is not always the outcome.

The value of the constraint function is descriptive (which pairs are coupled, how strongly), not predictive.

## Scholar Questions and Answers

**Q: Why is the accuracy so low?**
A: Because the task is fundamentally hard: predicting WHICH specific mutation occurs next. The model captures the coupling structure but not the stochastic mutation process. Even a perfect model of co-evolution would not predict individual mutations with high accuracy, because mutation is partly random.

**Q: What is the difference from LOO-CV (script 16)?**
A: This script uses an 800/499 train/test split. Script 16 uses leave-one-out (train on 1,298, test on 1). The LOO-CV accuracy is 2.93%, lower because it is a harder test (less training data per test... actually more, but the pair selection differs). Both are far below useful prediction.

**Q: Is the low accuracy a failure?**
A: No. It is a correct measurement of what K-map co-evolution analysis can and cannot do. The analysis is descriptive: it maps the constraint landscape. The low prediction accuracy is itself a finding about the nature of co-evolution (probabilistic, lineage-specific).

## Mermaid Diagram

```mermaid
flowchart TD
    A[co-evolving pairs top 10] --> B[frequency K-map]
    B --> C[C = ln(P / P_exp)]
    C --> D[sigmoid P = 1/(1+e^-C)]
    B --> E[train 800 sequences]
    C --> F[test 499 sequences]
    F --> G[argmax prediction]
    G --> H[accuracy 5.84%]
```
