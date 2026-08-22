# 11b — Strategy & Steps: K-map ↔ 3D-Structure Encoding Campaign

## Core design decision (binding)

**Unit of analysis = (protein, column-pair).** Unlike the Spike work (one protein, 299 models),
PSICOV-150 gives 149 usable proteins × ~10⁴ pairs each ≈ **1.5M labeled observations**
(label = native contact, Cβ<8Å, sep≥6 — literature convention, Jones 2012 / Morcos 2011).

Features are **per-protein-intrinsic** (no cross-protein alignment needed):
consensus He-code residue at i and j → 8 physicochemical groups; separation bin.
This is literally "the binary-encoded sequence predicts the 3D structure".

## Ground-truth sources (user directive: literature, not self-computed)
1. Native experimental PDBs (`raw/pdb/`) — physical truth. Identity mapping col↔residue
   verified 149/150; exception 1vjkA = one extra C-terminal residue (drop contacts w/ res 88).
2. PSICOV published predictions (`raw/con/*.out`) — the literature's own contact maps,
   format `i j 0 8 score` (cols 3-4 documented dummies), used ONLY as comparison baseline.

## Paths (each independently verifiable)

### PATH B — Universal contact circuits (headline)
- Level-1 coarse: K-map (g_i[3] × g_j[3] × sep_bin[2]) = 256 cells.
- Level-2 fine: per sep-bin, 32×32 padded (aa_i[5] × aa_j[5]) = 1024 cells ×4.
- Cell labeling: 1 if rate ≥ T·base AND n≥n_min; DC if n<n_min; else 0. Defaults T=2, n_min=30;
  sensitivity sweep {1.5,2,3}×{10,30}.
- Split: protein-level 100 train / 49 test (seed 42) + 5-fold protein CV.
- Eval: binary MCC/F1/enrichment at circuit membership; ranked precision@L/5,@L/2,@L using
  empirical cell-rate scores; baselines: random, sep-only ablation, **PSICOV published top-L/5**,
  per-protein plain-MI (recomputed via shared module).
- Interpretability: SOP rules ↔ biophysical priors (salt bridges acidic∧basic, H-bonds
  amide∧hydroxyl — rediscovered at scale?).

### PATH C — Three-way literature comparison
Per protein: precision@K of (a) our circuit scores, (b) PSICOV published scores,
(c) plain MI, vs native contacts. K ∈ {L/5, L/2, L}. Paired stats (Wilcoxon).

### PATH D — Inverse predictability
Pooled AUC (and per-protein mean AUC) of MI, perplexity-ratio, group-cell-rate for the contact
label; how much of the contact Boolean function is visible from sequence alone?

### PATH E — Gray-adjacency at proper power
h = Hamming(Gray(aa_i), Gray(aa_j)) of consensus residues over ALL native-contact pairs
(~10⁵ pairs vs Spike's 38). Background: resample (aa_i,aa_j) from the empirical joint
distribution of NON-contact pairs in the same sep-bin (controls composition+correlation).
Test h=1 enrichment; report effect size + CI.

### PATH F — Contact map AS a Boolean function (exploratory)
For L≤62 proteins: native contact matrix (padded to 64×64, 12-bit) → exact QM → cover size.
Compare compression (#PIs/#minterms) vs randomized maps matched on density AND separation
distribution. Question: do real contact maps minimize to fewer/larger implicants than chance?
(Runtime risk: measure first; fall back to smaller L or grouped merge.)

## Engineering
- New module `scripts/kmap_structure.py` (single source of truth for this campaign);
  runner scripts import it. Reuses `coevolution_shared`, `kmap_sbm.analysis.prime_implicants`,
  `nkmap.encoding.bio_sequences`.
- Dense MSA storage int8 [N,L]; max 74,836×266 ≈ 19.9 MB — trivial.
- GPU optional (MI); everything else numpy/scipy.
- Every run writes JSON under `results/kmap_structure/`; logs to stdout; fixed seeds.

## Verification gates (loop rule)
- G1 prototype: 3 proteins end-to-end, hand-checkable numbers.
- G2 adversarial tests green (empty columns, all-gap, tiny MSAs, degenerate cells, 1vjkA).
- G3 identity-mapping assertion per protein at load time (hard fail otherwise).
- G4 sensitivity sweep shows conclusions stable under labeling thresholds.
- G5 clean-state rerun byte-identical (fixed seeds).
- G6 Lean `lake build` clean with new theorems.
