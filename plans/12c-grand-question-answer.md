# 12c — Answering the Grand Questions (baseline, coupling, reconstruction)
**Date:** Aug 22, 2026

## Q1: "Random chance should be way higher, like 0.5 — how are you getting 0.179?"

precision@L/5 answers a *retrieval* question: "of the top-L/5 pairs our method
ranks highest, what fraction is a true native contact?" The candidate pool is
all residue pairs with sequence separation ≥ 6 — about L²/2 pairs per protein.
Only ~2.8% of those pairs are true contacts (contacts are rare!). So:

- A random ranker scores precision ≈ **0.028–0.034**, not 0.5.
- 0.5 would be coin-flip *accuracy* on a BALANCED problem; this is a ~1:30
  imbalanced retrieval problem.
- Our 0.179 means 17.9% of top predictions are real = **5.1× enrichment over
  random ranking**. DCA's 0.727 is 21× random — that is why DCA is the
  state of the art.

## Q2: "DCA encodes more coupling constants — can we do coupling in the K-map?"

Yes mechanically; no gain in precision. Three variants tested:

| Variant | prec@L/5 |
|---|---|
| groups × sep (L1) | 0.134 |
| groups × sep × **mfDCA-DI quartile** (circH) | 0.124 |
| identities × sep (L2) | 0.165 |
| raw mfDCA ranking (continuous) | 0.169 |
| PSICOV published (sparse inverse covariance) | 0.727 |

Why no lift: (a) binning continuous DI into 4 quartiles destroys exactly the
fine-grained ranking that makes DCA work; (b) at consensus level, group
identity cells already absorb most of what DI would add. The K-map's strength
is interpretable discrete rules; DCA's strength is continuous ranking fidelity.

**But one genuine positive emerged**: circH AUC 0.658 > raw-DCA AUC 0.583.
The K-map cell-rate smoothing acts as a denoising prior that improves GLOBAL
ranking even while losing top-K sharpness. Open direction: rank-fusion hybrid
(DCA provides physics; K-map prior re-ranks).

Implementation honesty: my first mfDCA had three real bugs (scrambled block
matrix from transpose+reshape; Frobenius over wrong axes; fixed-point
numerator broadcast on the wrong position axis). All caught by validating
against native contacts and against MI correlation; after fixes our mfDCA
prec@L/5 = 0.169 matches PUBLISHED mean-field benchmarks (~0.15–0.20,
Ekeberg/Morcos-generation methods on PSICOV150). The gap to 0.73 is PSICOV's
sparse-inverse-covariance method itself, not our code.

## Q3: "Can we reconstruct the complete 3D structure from the circuits?"

Quantified (PATH I), against literature thresholds only:

| Method | median LR prec@L/5 | recall@L | Reconstruction-ready? |
|---|---|---|---|
| circuit-L2 | 0.000 | 0.028 | NO |
| plain MI | ~0.06 | ~0.03 | NO |
| PSICOV published | **0.683** | **0.218** | **YES** |

Anchors: Jones2012 (≥0.5 L/5 long-range "sufficient"), Weigt2009 (~50%
folds TPRs), CONFOLD-line (~30% ⇒ useful folds).

**Verdict**: with current signal strength, NO — the circuits recover ~2.8% of
the contact graph at top-L versus DCA's 21.8%; folding needs orders of
magnitude more correct constraints than we supply. What circuits CAN do, and
now demonstrably do across held-out proteins:

1. **Positive channel**: sparse interpretable co-evolutionary rules
   (prec 2.6–5× base).
2. **Negative channel** (flipped circuit): forbidden-combination calls
   transfer across proteins — specificity 98.3%, held-out forbidden cells show
   exactly the predicted ≤1.75% contact rate.
3. **Conservation channel**: fully-conserved column pairs are **3.2× enriched**
   for native contacts (buried core signal) — the "genes staying constant"
   carry structure through conservation, not covariation.

So the honest reframing: circuits are a **compressed annotation layer** on the
structure (which pairs co-evolve / which are forbidden / which core positions
are locked), not a coordinate-generating replacement for it.

## Q4: "Look at both maps and complete the circuit"

Done (PATH K). The complete two-sided circuit (positive ∨ flipped, abstain
otherwise) achieves definitive-call MCC +0.169 at 22% coverage with the flipped
side transferring perfectly. This completes the COEVOLUTION_CONSTRAINTS.md
framework: f_obs AND NOT f_forb, now validated cross-protein.

## Reproduce
```bash
$PY scripts/path_g_tension.py                       # G1 hub + G2 lift
$PY scripts/path_h_dca_circuit.py --stage di        # mfDCA caches (~80 min)
$PY scripts/path_h_dca_circuit.py --stage cv        # circH vs rawDCA
$PY scripts/path_i_reconstruction.py                # recall curves + verdict
$PY scripts/path_k_complete_circuit.py              # complete two-sided circuit
```
All deterministic (fixed seeds); md5-stable on rerun (verified for G/H/I).

## Q5: Can the K-map layer ADD predictive power on top of DCA?

**YES — demonstrated for the first time** (PATH L rank-fusion):

| Method | prec@L/5 |
|---|---|
| rawDCA alone | 0.169 |
| circH (K-map prior) alone | 0.124 |
| **ranksum (DCA + K-map prior)** | **0.204 ± 0.104** |

Wilcoxon p≈0 vs both parents. The K-map cell-rate prior carries orthogonal
structural information not present in continuous coupling scores. This is
the strongest result of the campaign: the Boolean/K-map framework is not
merely an interpretability layer — it adds real predictive signal when fused
with standard co-evolution methods.
