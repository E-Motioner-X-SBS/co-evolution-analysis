# 10 — Boolean Contact Circuits from FASTA Data

## The research question
The existing co-evolution circuits (36 PIs from Quine–McCluskey) describe
*what co-evolved* in the Spike alignment. The 3D validation showed they are
**poor contact predictors** (0.36× pooled enrichment — worse than random).

**New question:** can we build Boolean circuits, trained on 3D contact labels,
that predict which residue pairs are in 3D contact — and do they generalize?

## Approach (3 phases)

### Phase 1: Per-pair contact circuits (residue identity → contact)
For each variable position pair (i,j), across the 299 models:
- Collect (aa_i, aa_j, is_contact) per model
- Build a 32×32 contact K-map: cell(a,b) = 1 if majority of models with
  (aa_i=a, aa_j=b) are in contact; 0 if not; -1 if unobserved
- QM minimize → Boolean expression (circuit)
- 5-fold CV: does the circuit predict contact on held-out models?
- **Key question:** does residue identity at borderline pairs determine contact?

### Phase 2: Cross-pair circuit (physicochemical group + separation → contact)
Pool ALL variable pairs' data. Features:
- group_i (3 bits, 8 physicochemical groups)
- group_j (3 bits)
- separation bin (2 bits: 3-5, 6-10, 11-20, 21+)
- Total: 8 bits = 256-cell K-map
- QM minimize → global circuit
- Train/test split by position pair (test on held-out pairs)
- **Key question:** is there a GENERAL rule (group + sep → contact)?

### Phase 3: Co-evolution transfer circuit (MI + perplexity + sep → contact)
Features per pair (from FASTA only, no 3D needed):
- MI bin (2 bits: <0.1, 0.1-0.3, 0.3-0.6, >0.6)
- perplexity ratio bin (2 bits: <1.0, 1.0-1.5, 1.5-2.0, >2.0)
- separation bin (2 bits)
- Total: 6 bits = 64-cell K-map
- Label: pair is "contact" if contact_frac > 0.3
- LOO-CV (leave out one pair)
- **Key question:** can co-evolution statistics predict contact?

## What "works" means
- Circuit beats the majority-class baseline (always predict "no contact")
- Precision > 0.20, recall > 0.30, F1 > 0.25, MCC > 0
- The Boolean expression is non-trivial (not constant 0 or 1)
- CV metrics ≈ training metrics (not overfit)

## Honest expectation
- Phase 1: should work for borderline pairs (contact varies by residue)
- Phase 2: may find general principles (hydrophobic-hydrophobic at sep 3-5)
- Phase 3: hardest test — co-evolution statistics are weak contact predictors
  (MI 1.82× enrichment); the circuit may not beat the baseline by much
