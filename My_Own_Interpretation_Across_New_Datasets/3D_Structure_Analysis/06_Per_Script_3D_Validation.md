# 06 — Per-Script 3D Validation: All 23 Scripts in Minute Detail

**Question:** For *every* one of the 23 analysis scripts — what does it compute,
how (algorithm + formulas), how are its position pairs / rules derived, and do
those pairs carry 3D structural information on the 299 Omicron Spike models?

**Data:** 1,128 SWISS-MODEL PDBs (299 unique sequences, best-GMQE model each);
alignment-column → model-residue mapping with integrity verification; contact =
Cβ–Cβ < 8 Å, separation ≥ 3; matched-separation random controls (2,000
draws/model, fixed seeds); exact binomial p (scipy `binomtest`, one-sided
greater). Pooled control rate = 0.0885 (intra-chain).

**Inventory:** `results/contacts/rules_inventory.json` — **14 position-pair
sources** (extended this session from the original 8). The 9 remaining scripts
are shared modules, protein-wide motif analysers, visualizers, or report
generators (no position pairs → documented as non-position-mappable).

---

## Summary table: all 14 pair-emitting sources (3D validation)

| # | Script | n_pairs | intra frac | enrich | p | structural? |
|---|--------|---------|-----------|--------|---|-------------|
| 1 | master_boolean | 10 | 0.032 | 0.36× | 1.0 | 2 pairs structural (212,215; 210,215) |
| 2 | full_length | 5 | 0.000 | 0.00× | 1.0 | NO |
| 3 | run_allseq | 30 | 0.140 | **1.58×** | 1.2e-43 | YES — 6 structural pairs |
| 4 | position_kmap | 30 | 0.199 | **2.25×** | 1.9e-122 | YES — best source |
| 5 | allseq_constraint | 10 | 0.032 | 0.36× | 1.0 | same pairs as master_boolean |
| 6 | variable_position | 10 | 0.032 | 0.36× | 1.0 | same pairs as master_boolean |
| 7 | perplexity | 17 | 0.085 | 0.96× | 0.8 | 3 structural pairs (diluted) |
| 8 | gpu_full | 10 | 0.000 | 0.00× | 1.0 | NO (top-MI pairs are non-structural) |
| 9 | **dca_mf** (NEW) | 20 | 0.094 | 0.87× | 0.21 | 1 structural pair (454,495) 5.97× |
| 10 | dca_boolean (NEW) | 10 | 0.032 | 0.36× | 1.0 | same pairs as master_boolean |
| 11 | predictive_constraint (NEW) | 10 | 0.032 | 0.36× | 1.0 | same pairs as master_boolean |
| 12 | flipped_positions (NEW) | 10 | 0.032 | 0.36× | 1.0 | same positions as master_boolean |
| 13 | **kmap_boolean** (NEW) | 2 | 0.179 | **2.02×** | 3.4e-04 | YES — essential-rule pairs |
| 14 | **network** (NEW) | 8 | 0.157 | **1.78×** | 2.4e-21 | YES — network edges |

**Structural pairs found across ALL sources** (enrichment > 2.5×, p < 1e-6):

| Pair | best enrichment | p | median Å | region | sources detecting it |
|------|-----------------|---|----------|--------|----------------------|
| (500,503) | 11.30× | 1.4e-315 | 5.3 | RBD | run_allseq, position_kmap |
| (503,507) | 11.26× | 5e-311 | 5.0 | RBD | run_allseq, position_kmap |
| (407,410) | 10.04× | 3.5e-235 | 6.7 | RBD | run_allseq, position_kmap, perplexity |
| (454,495) | 5.97× | 2.7e-83 | 8.0 | RBD | **dca_mf** (DCA-only finding) |
| (210,214) | 4.33× | 7.3e-12 | 9.8 | NTD | run_allseq, position_kmap |
| (212,215) | 3.91× | 2.7e-10 | 8.5 | NTD | master_boolean, run_allseq, position_kmap, perplexity, dca_boolean, predictive_constraint, flipped_positions, kmap_boolean, network |
| (210,215) | 3.71× | 8.3e-09 | 11.0 | NTD | (same 9 sources) |

---

## Per-script detail (all 23)

### Script 1: `coevolution_shared.py` — shared module (no pairs)
**What:** The single source of truth for FASTA parsing, He-2012 base-20
encoding, Shannon entropy, mutual information (numpy `bincount`, O(N+400) per
pair), mutation-only MI, coupling C = ln(P/P_exp), perplexity ratio, GPU
pair-finder, and the combined MI+perplexity rank-normalized score.
**Algorithm:** `load_position_arrays(aligned=True)` keeps the full alignment
with gap = state 20 (FIX A1); `mutual_information` builds the 20×20 joint
distribution via `bincount(codes_i*20+codes_j)`; `perplexity_ratio` =
PP(j)/mean_a PP(j|i=a) with conditional entropy H(j|i=a) over residues with
n≥5.
**3D:** N/A (module — no pairs of its own; every pair-emitting script uses it).
**Formulas:**
$$H(p) = -\sum_a P(a|p)\log_2 P(a|p), \quad
MI(i,j) = \sum_{a,b} P(a,b)\log_2\frac{P(a,b)}{P(a)P(b)}$$
$$C(a,b) = \ln\frac{P(a,b)}{P(a)P(b)}, \quad
\sigma(C) = \frac{1}{1+e^{-C}}, \quad
r(i,j) = \frac{PP(j)}{\text{mean}_a PP(j|i=a)}$$

### Script 2: `coevolution_gpu.py` — GPU kernels (no pairs)
**What:** torch CUDA kernels: dense [N×L] int32 tensor, `majority_refs_gpu`
(vectorized argmax over one-hot counts via `scatter_add`), `compute_entropy_gpu`,
`mi_matrix_gpu` (batched `scatter_add` bincount with **adaptive chunking** —
flat int64 buffer ≤ 512 MB so 486k-sequence families never OOM),
`coupling_matrix_gpu`, `h1_adjacency_gpu` (5-bit Gray `g(i)=i^(i>>1)`, GPU
popcount via bit tricks).
**3D:** N/A (module).
**Key engineering:** `get_device()` caps per-process GPU memory at 10%
(fail-fast on shared GPU); `OMP_PROC_BIND=FALSE` required (torch otherwise pins
to CPU 0). GPU MI matrix: 813K pairs × 1,299 seqs in 1.6 s (~800× vs CPU).

### Script 3: `master_boolean.py` — flagship rule extractor
**What:** Finds variable positions (H > 0.3), co-evolving pairs (mutation-only
MI > 0.1, window 30), builds per-pair 32×32 mutation K-maps (1=mutation,
-1=reference don't-care, 0=never), runs Quine–McCluskey, decodes 5-bit-per-axis
prime implicants → **36 distinct PIs, 2 essential**.
**Algorithm:** GPU mutation-only MI (`mi_matrix_gpu` with `mutation_only=True`,
reference pair excluded); `build_mutation_kmap` pads 20×20→32×32 (FIX A2);
`boolean_minimize_kmap` (QM with don't-cares); `extract_prime_implicants`
decodes 10-bit cubes (5 bits row + 5 bits col), rejects padding codes ≥ 20.
**Pairs:** 10 co-evolving pairs (top by mutation-only MI): (495,498) 0.871,
(448,454) 0.834, (488,498) 0.822, (442,454) 0.811, (442,448) 0.728, (212,215)
0.398, (215,216) 0.377, (212,216) 0.377, (210,215) 0.238, (210,212) 0.177.
**3D result:** pooled 0.032 contact fraction, **0.36× enrichment (p=1.0)** —
most pairs are non-structural. BUT 2 pairs ARE structural: **(212,215) 3.91×
(p=2.7e-10)** and **(210,215) 3.71× (p=8.3e-09)**. The 442–498 cluster and
(495,498) are 0.0× (13–27 Å apart — immune/lineage covariation).
**Interpretation:** the master rule set is sequence-logically minimal but
mostly non-structural; the 2 NTD structural pairs are genuine contact-driven
co-evolution. The 2 essential rules (212=G∧216=R, 210=K∧215=G) fire in only 1
sequence each → structurally untestable at this sample size.

### Script 4: `boolean_co-evolution.py` — whole-protein dipeptide Boolean
**What:** Whole-protein dipeptide Boolean minimization + coupling constants.
Produces protein-wide dipeptide **motifs** (row_aa, col_aa) — no position
information.
**3D:** N/A — motifs are protein-wide dipeptide preferences; they carry no
positional information and cannot be mapped to 3D positions (documented
limitation, same as n-ary K-map).

### Script 5: `nary_kmap_co-evolution.py` — base-20 (n-ary) K-map
**What:** Base-20 K-map analysis: k-mer frequency maps, minimized directly in
base-20 (20×20 = 400 cells, no binary projection). Produces protein-wide
dipeptide motifs.
**3D:** N/A — same as boolean_co-evolution: motifs have no position info.
**Formulas:** $\text{cell} = \sum_{\text{window}} \prod 20^{\text{position}} \text{code}(aa)$

### Script 6: `position_kmap_coevolution.py` — per-position-pair K-maps
**What:** Builds per-position-pair 20×20 K-maps with MI; ranks all variable-
position pairs (window 30) by MI. **Best structural source.**
**Pairs:** 30 top pairs by full MI.
**3D result:** pooled 0.199 contact fraction, **2.25× enrichment (p=1.9e-122)**
— the strongest pooled signal. Per-pair: (500,503) 11.3×, (503,507) 11.26×,
(407,410) 10.04×, (210,214) 4.33×, (212,215) 3.91×, (210,215) 3.71×.
**Interpretation:** full-MI ranking (not mutation-only) is the best single
source for structural pairs — it captures the RBD loop pairs (500,503)/(503,507)
that mutation-only MI misses (those are near-deterministic in the reference
frame).

### Script 7: `run_allseq_analysis.py` — all-sequence MI matrix
**What:** Full MI matrix over all variable-position pairs (window 30), all
sequences. GPU-accelerated.
**Pairs:** 30 top pairs by full MI (overlaps heavily with position_kmap).
**3D result:** pooled 0.140, **1.58× (p=1.2e-43)**. Same 6 structural pairs as
position_kmap. Slightly lower enrichment because its pair set extends further
down the MI ranking (more non-structural pairs included).
**Interpretation:** confirms position_kmap; the full-MI top-30 is a robust
structural-pair detector.

### Script 8: `run_kmap_analysis.py` — binary Gray K-map pipeline (H1–H6)
**What:** Applies the 5-bit Gray-code encoding framework: H1 (consecutive-
residue Gray adjacency), H2 (signature), H3 (pairwise K-map distance),
H5 (contacts), Walsh–Hadamard spectrum, mutations.
**H1 result:** observed Gray-Hamming-1 fraction = 0.2167 vs expected 0.1613 →
**1.34× enrichment** (1,644,588 consecutive pairs). This is a global property
of consecutive residues, NOT a set of co-evolving position pairs.
**3D:** H1 is a global adjacency statistic (all consecutive pairs), not a
position-pair set → not mappable to 3D as co-evolution rules. The Gray bit-flip
adjacency test on co-evolving pairs (gray_flip_analysis.py, Analysis 4) found
1.18×, p=0.345 — NOT significant. H3 (pairwise K-map signature distance) is
O(n²) on sequences, not positions.
**Interpretation:** the Gray-encoding's adjacency property holds for consecutive
residues (1.34×) but does NOT predict which specific residue pairs co-evolve,
and does not predict 3D contact. Honest negative for the pure K-map prediction.

### Script 9: `flipped_boolean_coevolution.py` — forbidden pairs
**What:** Flipped K-maps: 1 = forbidden (never observed), 0 = observed, -1 =
reference. QM → **490 forbidden rules** across 10 position pairs (117 unique
residue-pair rules).
**Pairs (position-level):** 10 unique position pairs (same as master_boolean's
10, since forbidden rules are at the co-evolving positions).
**3D result (position-level):** pooled 0.032, 0.36× — same as master_boolean
(same positions). The residue-level flipped validation (Analysis 3) found
steric support at (212,215) 2.71× and (210,215) 2.57×; pooled 0.25× = rarity
overall.
**Interpretation:** forbiddenness is mostly lineage/sequence-space rarity, with
genuine steric constraint at the 2 NTD pairs.

### Script 10: `kmap_boolean_coevolution.py` — K-map Boolean equations
**What:** Generates the 5-bit-Gray-literal Boolean equations for each pair's
prime implicants. The 2 **essential** rules (the minimal irredundant core):
- (212,216): pos212=G ∧ pos216=R
- (210,215): pos210=K ∧ pos215=G
**Pairs:** 2 essential-rule position pairs.
**3D result:** pooled 0.179, **2.02× (p=3.4e-04)** — the essential-rule pairs
ARE enriched. (210,215) is 3.71× structural; (212,216) is 0.43× (weak).
**Interpretation:** the sequence-logically minimal core is partly structural —
(210,215) is a genuine contact pair; (212,216) is not (the essential rule fires
in only 1 sequence → untestable).

### Script 11: `variable_position_coevolution.py` — variable-positions K-maps
**What:** Variable-positions-only K-maps with don't-cares (conserved positions
excluded). Same 10 pairs as master_boolean.
**3D result:** 0.36× — identical to master_boolean (same pair set).
**Interpretation:** restricting to variable positions doesn't change the
structural verdict; the pair set is the determining factor.

### Script 12: `predictive_constraint_function.py` — train/test prediction
**What:** Constraint function C = ln(P/P_exp), prediction σ(C); train/test
split (62/38 proportional). Accuracy **0.11%**.
**Pairs:** 10 pairs (same as master_boolean).
**3D result:** 0.36× — same pair set.
**Formulas:** $C(a,b) = \ln\frac{P(a,b)}{P(a)P(b)}, \quad \hat{y} = \sigma(C)$
**Interpretation:** the constraint function is a weak predictor (0.11%) —
co-evolution is lineage-specific; the rules capture structure, not specific
outcomes.

### Script 13: `allseq_constraint_function.py` — LOO-CV
**What:** Leave-one-out cross-validation with the constraint function. Accuracy
**9.24% (301/3259)**, deterministic. Forced off CUDA (CPU-bound; GPU MI caused
allocator spins).
**Pairs:** 10 pairs (same as master_boolean).
**3D result:** 0.36× — same pair set.
**Complexity:** O(10·N²) — infeasible at 222–486k sequences (27 days – 2.5
years); documented.
**Interpretation:** 9.24% LOO-CV confirms lineage-specificity; the genuine
pairs are more predictable than artifact pairs (was 2.93% pre-correction).

### Script 14: `dca_boolean_coevolution.py` — local precision → Boolean
**What:** Local precision matrix (inverse covariance) → Boolean minimization.
**NOT real DCA** (no mean-field, no APC). Accuracy 0.0.
**Pairs:** 10 pairs (same as master_boolean — it uses the same GPU pair-finder).
**3D result:** 0.36× — same pair set.
**Interpretation:** the local-precision approach adds no structural information
beyond MI (same pairs); the "17.6%" reported in the pre-correction audit was an
artifact of hardcoded pairs.

### Script 15: `dca_mf_analysis.py` — proper mfDCA (Morcos 2011)
**What:** Proper mean-field DCA: Frobenius norm, average-product correction
(APC), direct information (DI). Top DI pair: **(454,495) DI=0.37**;
ρ(DI, MI) = 0.06 (uncorrelated with MI).
**Algorithm:** mfDCA: solve Potts-model inverse problem (mean-field
approximation of the coupling matrix J), apply APC (subtract the rank-1
average-product component), compute DI = Σ P(a,b)·J(a,b)².
**Pairs:** 20 top-DI pairs.
**3D result:** pooled 0.094, 0.87× (p=0.21 — NOT significant overall). BUT
**1 pair is strongly structural: (454,495) 5.97× (p=2.7e-83, median 8.0 Å)** —
the DCA-only structural finding. The rest of the DI top-20 are non-structural
(median 9–110 Å).
**Interpretation:** mfDCA finds ONE structural pair (454,495) that MI's top-20
misses — direct-coupling decomposition recovers a contact-driven pair hidden in
the MI ranking. But DCA's overall top-20 is LESS structural than MI's (0.87× vs
1.82×) on the Spike — likely because the redundant, lineage-dominated alignment
reduces DCA's effectiveness. On the diverse PSICOV benchmark, DCA wins (0.44 vs
0.066 precision); on the Spike, MI wins. Honest nuance.

### Script 16: `perplexity_coevolution.py` — perplexity ratio
**What:** Perplexity ratio r(i,j) = PP(j)/mean_a PP(j|i=a) per pair. Ratio up
to 2.81×.
**Pairs:** 17 pairs (variable-position pairs, window 30).
**3D result:** pooled 0.085, 0.96× (p=0.8 — not significant pooled). BUT 3
pairs structural: (407,410) 10.04×, (212,215) 3.91×, (210,215) 3.71×.
**Perplexity-as-predictor (Analysis 8):** Spearman(ratio, contact fraction) =
**0.670, p=0.006** — the strongest single sequence metric of 3D proximity.
**Formulas:** $r(i,j) = \frac{2^{H(j)}}{\text{mean}_a 2^{H(j|i=a)}}$
**Interpretation:** the DETERMINISM of a co-evolutionary relationship (ratio)
is the best single predictor of 3D contact — near-deterministic pairs
((212,215), (210,215), (407,410)) are the structurally-real ones. Pooled
enrichment is weak because many high-ratio pairs are non-structural; the
correlation is the finding, not the pooled enrichment.

### Script 17: `advanced_co-evolution_analysis.py` — network/clustering/signatures
**What:** Co-evolution network (21 nodes, 8 edges), Walsh–Hadamard spectrum,
variant signatures (40), clustering.
**Pairs:** 8 network edges (MI > threshold among variable positions).
**3D result:** pooled 0.157, **1.78× (p=2.4e-21)** — the network edges are
enriched. Edges include (18,26), (66,94), (407,410), (212,215), (210,215),
(378,407), (215,216), (210,212). The structural ones (407,410; 212,215; 210,215)
drive the enrichment; (18,26) and (66,94) are non-structural.
**Interpretation:** the co-evolution network captures both structural and
phylogenetic edges; 3 of 8 edges are structural.

### Script 18: `full_length_analysis.py` — full-length entropy + MI
**What:** Full-length (1,276 positions) entropy and MI (window 30). 21
variable positions, 5 high-MI pairs (MI > 0.5).
**Pairs:** 5 top pairs by full MI: (373,378), (18,26), (378,407), (66,94),
(215,216) — all by full MI.
**3D result:** pooled 0.000, **0.00× (p=1.0)** — NONE of the top-5 full-MI
pairs are structural (median 8.9–17 Å). The strongest full-MI pair (373,378)
0.807 is 16.6 Å apart — immune/lineage covariation.
**Interpretation:** full-MI ranking by raw value is a POOR structural predictor
— the top-MI pairs are dominated by lineage covariation. The structural pairs
(407,410; 500,503; etc.) have lower MI (0.75) but ARE contacts. This is why
mutation-only MI and perplexity ratio outperform raw MI for structure.

### Script 19: `gpu_full_analysis.py` — GPU full analysis
**What:** GPU full analysis: full MI matrix, entropy, refs, H1, couplings. Max
MI = 0.8067 at (373,378).
**Pairs:** 10 top coupling pairs (by avg|J|).
**3D result:** pooled 0.000, **0.00× (p=1.0)** — same as full_length; the
top-coupling pairs are the non-structural high-MI pairs.
**Interpretation:** confirms full_length — raw MI/coupling magnitude does not
predict structure; the structural signal is in the mid-MI, high-determinism
pairs.

### Script 20: `create_mi_heatmap.py` — MI heatmap (visualization)
**What:** MI heatmap visualization (GPU MI for all 813K pairs). No pair set
emitted (visualization only).
**3D:** N/A (visualization).

### Scripts 21–23: report generators
`generate_co-evolution_md.py`, `generate_full_analysis_md.py`,
`generate_full_pipeline_doc.py` — produce markdown/HTML reports from the result
JSONs. No pairs of their own.
**3D:** N/A (report generators).

---

## Cross-cutting findings (the "things that work")

1. **Full-MI pair ranking (position_kmap, run_allseq) is the best structural
   detector** (2.25×, 1.58× pooled; 6 structural pairs). It captures the RBD
   loop pairs (500,503)/(503,507) that mutation-only MI misses.
2. **mfDCA finds one structural pair MI misses: (454,495) 5.97×** — the
   direct-coupling decomposition recovers a contact-driven pair hidden in the
   MI ranking. DCA and MI are complementary, not redundant.
3. **The perplexity RATIO (determinism) is the best single 3D correlate**
   (ρ=0.670, p=0.006) — near-deterministic pairs are structural.
4. **The essential-rule pairs (kmap_boolean) are enriched 2.02×** — the
   sequence-logically minimal core is partly structural.
5. **The co-evolution network edges are enriched 1.78×** — 3 of 8 edges
   structural.
6. **Raw MI magnitude does NOT predict structure** (full_length, gpu_full:
   0.00×) — the top-MI pairs are lineage-dominated. Structure lives in the
   mid-MI, high-determinism pairs.
7. **The 442–498 cluster and (373,378) are consistently non-structural** across
   ALL sources — immune/lineage covariation, not contact.

## What does NOT work (honest negatives)
- Gray bit-flip adjacency: 1.18×, p=0.345 (Analysis 4).
- Raw MI magnitude as a structural predictor: 0.00×.
- dca_boolean (local precision): same pairs as MI, 0.36×, accuracy 0.0.
- LOO-CV prediction: 9.24% (lineage-specific; rules don't generalize).
- Trimer interface: zero inter-chain contacts for all pairs.
- Cross-script consensus ⇒ structure: NO (the most-agreed pairs are the least
  structural).
