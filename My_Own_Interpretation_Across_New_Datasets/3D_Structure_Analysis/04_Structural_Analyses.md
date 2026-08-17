# 04 — The 9 Structural Analyses: questions, formulas, results, scripts

Nine distinct analyses were run on the 3D models. For each: the research
question, the exact method (formulas), the statistical values obtained, and
the script that produces them.

**Common setup** (all analyses):
- 299 unique Omicron sequences → best-GMQE SWISS-MODEL model each (1,128 PDBs
  available; best used; robustness across models implicitly via per-model
  agreement).
- Alignment column → model residue mapping with integrity verification
  (see `02_Contact_Computation.md`).
- Contact: Cβ–Cβ < 8 Å, separation ≥ 3; matched-separation random controls
  (2,000 draws/model); exact binomial p (scipy `binomtest`, one-sided
  greater).

---

## Analysis 1 — Position-level rule→3D validation
**Script:** `scripts/validate_3d_complete.py` → `results/contacts/validate_3d_complete.json`
**Question:** do the co-evolving position pairs from EACH of the 8 pair-emitting
scripts sit close in 3D?

**Method:** for each source's pair set, pooled contact fraction and enrichment
vs the shared matched control (pooled rate p₀ = 0.0956), plus per-pair stats
and the trimer-interface (inter-chain) fraction.

**Results:**
| Source | contact fraction | enrichment | p |
|---|---|---|---|
| position_kmap (30 pairs) | 0.199 | **2.08×** | 1.3e-103 |
| run_allseq (30 pairs) | 0.140 | **1.46×** | 1.1e-31 |
| perplexity (17 pairs) | 0.085 | 0.89× | 0.99 |
| master_boolean (10) / allseq_constraint / variable_position | 0.032 | 0.33× | 1.0 |
| full_length / gpu_full | 0.000 | 0.0× | 1.0 |

Per-pair significant (intra-chain): **(407,410)** 0.889 / 9.3× / p=1.5e-226;
**(210,214)** 0.384 / 4.0× / p=4.5e-11; **(212,215)** 0.346 / 3.6× / p=1.4e-09;
**(210,215)** 0.329 / 3.4× / p=3.6e-08; **(500,503)** 1.000 / 10.5×;
**(503,507)** 0.997 / 10.4×. Non-structural: (18,26), (66,94), (373,378),
(378,407), (442,448), (442,454), (448,454), (488,495), (488,498), (495,498),
479–507 periphery (0 contacts, 13–73 Å).

## Analysis 2 — Residue-level: the 36 prime implicants
**Script:** `scripts/validate_residue_rules_3d.py` → `validate_residue_rules_3d.json`
**Question:** when a Boolean rule `IF pos_i=aa_i AND pos_j=aa_j` FIRES (both
residues occur in a model), are those residues in contact?

**Method:** per rule, count models where it fires; measure contact; control as
in Analysis 1 (pooled p₀ = 0.1187).

**Results:**
| Rule | fired | contact frac | enrichment | p |
|---|---|---|---|---|
| [210=I, 215=P] | 71 | 0.338 | **2.85×** | 1.1e-06 |
| [212=V, 215=P] | 76 | 0.329 | **2.77×** | 1.2e-06 |
| [212=S, 215=G] | 2 | 1.000 | 8.4× | 0.014 |
| all others (442–454, 488–498, 495–498, 212–216…) | — | 0.000 | 0.0× | 1.0 |
| **Pooled 36 PIs** | 780 | 0.069 | 0.58× | 1.0 |

The 2 essential rules fire in exactly 1 of 1,299 sequences each (212=G occurs
once); in that model the residues are sequentially adjacent (separation 1) —
a trivial backbone pair, excluded as non-independent. Honest limit: the
essential core is sequence-minimal but structurally untestable at this sample
size.

## Analysis 3 — Flipped (forbidden) rules: the steric hypothesis
**Script:** `scripts/validate_flipped_3d.py` → `validate_flipped_3d.json`
**Question:** are forbidden (never-observed) residue-pair POSITIONS 3D-close
(steric incompatibility) or distant (sequence-space rarity)?

**Results (117 unique forbidden rules on 10 position pairs):**
| Pair | contact frac | enrichment | p |
|---|---|---|---|
| (212, 215) | 0.346 | **2.71×** | 6.4e-07 |
| (210, 215) | 0.329 | **2.57×** | 7.3e-06 |
| others (442–498 cluster) | 0.000 | 0.0× | 1.0 |
| **Pooled** | 0.032 | 0.25× | 1.0 |

**Verdict: RARITY overall, with genuine steric constraint at the two NTD
pairs** — forbidden combos there would clash in 3D.

## Analysis 4 — Gray bit-flip adjacency (the pure K-map prediction)
**Script:** `scripts/gray_flip_analysis.py` → `gray_flip_analysis.json`
**Question:** are observed co-evolving residue pairs Gray-adjacent
(Hamming distance 1 in 5-bit Gray space) more than random — the "consistency
when you flip the bits" prediction?

**Method:** Hamming distance of the 38 observed on-set residue pairs vs the
null (all 20×20 pairs; P(H=1) = 5/31 = 0.200 under the 5-bit Gray embedding).

**Result: H=1 observed 23.7% vs 20.0% expected → 1.18×, p = 0.345 — NOT
significant.** Gray-adjacency does not predict observed co-evolution, and it
does not co-occur with 3D contact except at (212,215).

## Analysis 5 — Structural profile: the "housekeeping" hypothesis
**Script:** `scripts/structural_profile.py` → `structural_profile.json`
**Question:** are conserved positions buried (structural core) and variable
positions exposed (immune surface)?

**Method:** burial proxy = number of Cβ within 12 Å of a residue; per alignment
column, fraction of models where the mapped residue is more buried than the
model median; Spearman correlation with per-column Shannon entropy.

**Results:**
- **Spearman(entropy, burial) = −0.170, p = 6.4e-09** — significant negative:
  variable positions are exposed.
- Conserved positions (H = 0): mean buried fraction 0.519; variable
  (H > 80th percentile): 0.260.
- All co-evolving positions (210–216 NTD, 373–498 RBD/SD) are 0–10% buried —
  surface-exposed, consistent with epitope-driven co-evolution.
- (66, 94), (18) are 79–100% buried (conserved core) yet their co-evolution is
  NOT contact-based.

## Analysis 6 — Trimer interface
**Script:** `scripts/trimer_interface.py` → `trimer_interface.json`
**Question:** do any co-evolving pairs sit at the INTER-PROTOMER interface?

**Method:** inter-chain contact = min over chain pairs (X≠Y) of
d(residue_i in X, residue_j in Y) < 8 Å; matched control.

**Result: ZERO inter-chain contacts for all 48 pairs** (median inter-distances
18–69 Å). The co-evolution signal is entirely intra-protomer.

## Analysis 7 — Cross-script consistency
**Script:** `scripts/consistency_cross_script.py` → `cross_script_consistency.json`
**Question:** do the 23 scripts agree on which positions co-evolve, and does
agreement imply structure?

**Results:** overlap matrix (shared pairs / Jaccard) across 8 sources;
consensus pairs (≥2 sources) cross-checked against 3D. Key finding:
**(212,215) — 6 sources AND structural (3.6×)**; but the 5-source consensus
pairs (442,448), (442,454), (448,454), (488,498), (495,498), (215,216) are
ALL non-structural. **Consistency ≠ structure** — the shared signal is
lineage/immune covariation.

## Analysis 8 — Perplexity as a structural predictor
**Script:** `scripts/perplexity_3d_analysis.py` → `perplexity_3d_analysis.json`
**Question:** does the perplexity RATIO (determinism) correlate with 3D
contact?

**Method:** for all 17 variable-position pairs (window 30): ratio + MI +
combined (shared module), correlated (Spearman) with the pair's intra-chain
contact fraction (n = 15 pairs with ≥50 models; pooled control p₀ = 0.129).

**Results:**
| Metric | Spearman ρ | p |
|---|---|---|
| **Perplexity ratio** | **0.670** | **0.006** |
| Combined (MI+ratio) | 0.444 | 0.098 |
| MI | 0.266 | 0.338 |

Ratio quartiles: Q1–Q2 mean contact 0.000 (0.0×); Q3 0.232 (1.80×); Q4 0.170
(1.31×). Top-ratio pairs: (212,215) 2.68×, (210,215) 2.55×, (407,410) 6.88×.
**The determinism measure is the strongest single sequence metric of 3D
proximity** (honest caveats: n=15; high ratio does not cleanly separate —
(378,407), (66,94) at ratio ~1.7 are non-structural; no "rescued" pairs beyond
MI).

## Analysis 9 — PSICOV 150 benchmark (cross-dataset validation)
**Script:** `scripts/validate_contacts.py` → `results/contacts/contact_validation.json` + `psicov_aggregate.json`
**Question:** on 150 proteins with NATIVE structures, are our top co-evolving
pairs on native contacts, and how does precision compare to published DCA?

**Results:** 150/150 targets validated (blank-chain PDB handling added):
- mean enrichment **2.49×** (69/150 > 2×, 25/150 > 5×; max 1dbxA 16.7×);
- mean precision 0.066 vs published PSICOV 0.44 (top-L/5) — **the honest
  finding that plain MI (no APC) detects real signal but underperforms
  state-of-the-art DCA**.

---

## Consolidated statistical summary

| Analysis | Key statistic | p |
|---|---|---|
| 1 position-level (best source) | 2.08× (position_kmap) | 1.3e-103 |
| 1 per-pair best | (407,410) 9.3× | 1.5e-226 |
| 2 residue-level best | 2.85× | 1.1e-06 |
| 2 pooled 36 PIs | 0.58× | 1.0 |
| 3 flipped pooled | 0.25× | 1.0 |
| 3 flipped best pair | 2.71× | 6.4e-07 |
| 4 Gray H=1 | 1.18× | 0.345 |
| 5 entropy–burial | ρ = −0.170 | 6.4e-09 |
| 6 trimer inter-chain | 0 contacts | — |
| 7 consensus vs structure | (212,215) 6 sources, 3.6× | 1.4e-09 |
| 8 perplexity–contact | ρ = 0.670 | 0.006 |
| 9 PSICOV enrichment | 2.49× mean | — |
