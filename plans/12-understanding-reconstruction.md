# 12 — Understanding: Structure Reconstruction from Entropic Tension + K-map Circuits

**Date:** Aug 22, 2026 (campaign reopened; INFINITY_DONE deleted)
**Directive (user, restated):**
1. Challenge on the baseline: "random chance should be way higher, like 0.5" —
   needs explicit clarification of imbalanced-precision baselines.
2. DCA reaches 0.7–0.8 because it encodes far more coupling structure. **Can we
   do coupling INSIDE the K-map?** Proposal: a *tension mechanism* — for ONE
   sequence element, look at its entropic relations with EVERYTHING else
   (full-length profiles, not just windowed pairs), then deduce structure.
3. **Grand question:** can the minimized Boolean circuits REPLACE / reconstruct
   the protein's 3D structure? How much can be reconstructed — including
   conserved positions? Quantify honestly.
4. Study across datasets; downloads allowed; compare against LITERATURE values
   only; verify 3D structures correctly; find where theory is wrong and fix it.

## Baseline clarification (must be written up first)

precision@L/5 answers: "of our top-L/5 ranked candidate pairs, what fraction is
a native contact?" The candidate universe is all pairs with sep ≥ 6 (~L²/2
pairs per protein). Only ~2.8% of those are true contacts ⇒ a RANDOM ranker
scores precision ≈ base rate ≈ 0.03, not 0.5. 0.5 would be the accuracy floor
of a balanced coin-flip problem; contact prediction is a ~1:30 imbalanced
retrieval problem. Our 0.179 = 5.1× enrichment over random ranking.
DCA's 0.7+ is the literature state of the art on this exact benchmark.

## Why DCA wins (theory to test)

mfDCA inverts the GLOBAL covariance (19L×19L) — it disentangles direct from
indirect (transitive/phylogenetic) correlations and applies APC. Our circuits
so far use only local consensus identities + separation: no global coupling
information at all. Hypothesis H-NEW: adding coupling-derived features
(mfDCA DI bins, APC-corrected MI bins) as K-map input dimensions will close
part of the gap while keeping the interpretable SOP form.

## New paths

### PATH G — Tension profiles (user's mechanism)
For position i, tension toward j := normalized entropic dependence
T(i,j) = MI(i,j) / max_k MI(i,k) over cached partners k of i (and symmetric).
Per-pair features: T(i,j), T(j,i), ranks. Questions:
 G1 Do high-tension positions correspond to contact-map hubs?
 G2 Does adding tension bins to the universal circuit lift precision?
 G3 Are tension profiles reproducible across homologous folds (literature)?

### PATH H — Coupling inside the K-map (the direct ask)
Compute proper mfDCA DI per PSICOV protein (Morcos 2011; we already have an
implementation pattern in dca_mf_analysis.py), bin DI into 2 bits, extend the
universal circuit to (g_i,g_j,sep,DI-bin) → QM → circuit. Compare lifts:
  circuit(sep) → circuit(+groups) → circuit(+DI) → raw mfDCA ranking.
Literature anchor: PSICOV paper's own published numbers (0.73 L/5).

### PATH I — Reconstruction quantification (the grand question)
Literature anchors: Jones 2012 ("118/150 targets ≥0.5 L/5 long-range precision
— sufficient to benefit structure prediction"); Weigt 2009 (top-L contacts fold
TPR); CONFOLD/Cheng work (contact precision vs TM-score curves). Deliverables:
 I1 Recall/coverage curves: fraction of native contacts recovered at K ∈ {L/10,
    L/5, L/2, L} per method (circuit-L2, MI, mfDCA-DI, PSICOV-con).
 I2 Reconstruction-readiness verdict vs literature thresholds: at OUR precision,
    how many correct contacts do we supply, and what TM-score-class does the
    literature predict from that? (Honest: likely below folding threshold.)
 I3 Conserved-core coverage: do circuit rules additionally COVER the conserved
    (non-coevolving) structural core via burial/group rules?

### PATH J — Second dataset (downloads allowed)
EVmutation MSAs already local (20 proteins, up to 28k seqs). Fetch their PDB
structures from RCSB (Hopf 2017 benchmark set), map via A2M reference-frame
(validate_contacts.py already implements this), rerun PATH B/G/H evaluation.
Also fetch ORIGINAL RCSB entries for PSICOV targets to cross-check the
repaired natives (GT verification, user's "make sure looking at 3D structures
correctly").

## Success criteria
- [ ] Baseline clarification written into the campaign doc (0.03 not 0.5).
- [ ] PATH G: tension features computed; hub test + circuit lift measured.
- [ ] PATH H: mfDCA-in-K-map circuit trained/evaluated with protein-level CV.
- [ ] PATH I: recall curves + reconstruction-readiness verdict vs literature.
- [ ] PATH J: ≥1 new dataset evaluated end-to-end (or documented blocker).
- [ ] All comparisons anchored to literature numbers (Jones 2012 Table 1;
      Morcos 2011; Hopf 2017) — never self-referential.
- [ ] Lean: any new discrete object formalized if load-bearing.
- [ ] Docs updated; determinism checks; double audit.

## Open questions
- Q1: mfDCA runtime per PSICOV protein (reweighting O(N²·L)) — measure first.
- Q2: EVmutation→PDB ID mapping provenance (web-search Hopf 2017 supplement).
- Q3: does DI binning add beyond MI binning (they correlate ρ≈0.06 on Spike)?
