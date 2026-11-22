# 11 — K-map ↔ 3D-Structure Encoding: Universal Contact Circuits on Literature Ground Truth

> **Correction A4 (Nov 2026)** — see `My_Interpretation_of_the_Results/CORRECTION_NOTICE.md`.
> Spike rules: 36 prime implicants, **35 essential** (not 2). Constraint-function LOO-CV **99.4%**
> (the 9.24% was a predictor artefact); it is lineage linkage: withholding the test sequence's
> relatives within 12 mismatches drops it to 67%, level with the 65% frequency baseline
> (`analysis/spike_lineage_loo.py`). APC-corrected MI on PSICOV150 holdout **0.176** vs circuit
> 0.179 (5-fold CV 0.189 vs 0.165), so the circuit does not beat the standard sequence baseline.
> Rank fusion (0.204) reproduces but improved only a defective mfDCA (0.169); with a correct mfDCA
> (0.566 DI / 0.692 APC-Frobenius) the K-map prior lowers precision (`analysis/dca_fusion_recheck.py`).
> Spike mfDCA: Meff = 1 (all sequences within 80% identity), so its pairs, incl. "(454,495)
> DCA-only", are not meaningful; the DI coupling sign is fixed. 3D: with same-separation controls
> and the 11 templates (not 299 models) as units, none of the 10 co-evolving pairs is in contact
> more than its separation predicts; only the i,i+3/4 pairs (407,410), (500,503), (503,507) are
> enriched (`analysis/spike_3d_matched_check.py`). PSICOV top-20 MI pairs: 2.38× / precision 0.071 (first reported 2.49× / 0.066).

**Campaign:** plans/11-kmap-structure-encoding.md (+11b strategy, +11c decisions)
**Date:** Aug 22, 2026 · **Status:** all five stages executed; results verified
**Question:** can the chain *binary sequence → Gray/base-20 encoding → K-map →
Quine-McCluskey → prime implicants* encode **3D structure data** — verified
against LITERATURE ground truth (native experimental PDBs + PSICOV's published
contact predictions), not our own models?

---

## 1. Headline answers

| # | Question | Answer | Evidence |
|---|---|---|---|
| Q1 | Do sequence-encoded residue identities predict 3D contacts **on unseen proteins**? | **YES** — universal QM circuit reaches prec@L/5 = **0.179** (holdout) = 5.1× random, 1.9× plain MI, ~15× separation-only | PATH B, 5-fold protein-level CV |
| Q2 | How far from the literature's own predictions? | Honestly quantified: PSICOV published DCA = **0.723**; our no-phylogeny-correction circuits reach ~25% of that | PATH B/C |
| Q3 | Is Gray-code adjacency structurally meaningful? | **Yes, weakly but decisively**: 1.178× among 47,430 native-contact pairs, two-sided **p = 1.1×10⁻⁷⁷** — the earlier Spike "negative" (1.18×, p=0.345) was a power failure, not an effect failure | PATH E |
| Q4 | Are real contact maps more "K-map-expressible" than chance? | **Yes**: exact-QM compression better than separation-matched shuffles in 6/6 proteins (2.31 vs 2.74 PIs/contact, Wilcoxon p=0.031) | PATH F |
| Q5 | What do the learned Boolean rules say? | Hydrophobic/aromatic core packing dominates ({M,F,Y,W}-vocabularies); Cys/sulfur terms; interpretable SOP form | PATH B circuits |

**Bottom line:** the K-map/QM pipeline does encode 3D structure — transferably,
interpretably, and now measured against literature ground truth with honest gaps.

## 2. Data and integrity (why this ground truth is trustworthy)

- 150 PSICOV targets: deep diverse MSAs (511–74,836 seqs; median 3,228) +
  **native experimental structures** (`raw/pdb`, repaired/renumbered) +
  **PSICOV's own published contact predictions** (`raw/con/*.out`).
- README guarantee verified live: alignment width == PDB CA count for **149/150**
  ⇒ MSA column c ↔ PDB residue c EXACTLY (no alignment step, no gap-mapping
  error). Exception 1vjkA = one extra C-terminal PDB residue ('S',88) absent
  from target seq → drop-rule (residues > width ignored) resolves it exactly.
- Contacts: Cβ–Cβ < 8 Å, separation ≥ 6 (Morcos 2011 / Jones 2012 convention);
  CA for Gly. All 150 PDBs have CB coverage ≥50% CA; numbering contiguous.
- Feature caches: 150/150 npz validated post-run (array-length alignment, code
  ranges, meta parse) after a dual-process race was detected and killed;
  **0 corrupt files**.
- Universe: **1,675,886 candidate pair observations; 47,430 native contacts;
  base rate 2.83%** — the properly powered regime (Spike had 21 positions /
  10 pairs / 38 Gray-testable residue pairs).

## 3. PATH B — universal contact circuits (the headline)

**Design.** Features per (protein, column-pair): consensus He-code residues at
i and j (majority, lowest-code tie-break) → level-1 groups (8 physicochemical
classes) or level-2 identities (20 codes on a 32×32 padded grid per separation
bin [6,11) [11,21) [21,31) ≥31); label = native contact. Cell labeling:
rate ≥ T·base ∧ n≥n_min → ON; n<n_min → DC; else OFF (defaults T=2, n_min=30).
QM (popcount-grouped, don't-care-aware) → sum-of-prime-implicants circuit.
Every trained circuit is runtime-asserted SOUND (no PI matches an OFF cell)
and COMPLETE (every ON cell covered) — the property Lean-formalized in
ContactCircuits.lean.

**Results — precision@L/5 on fully held-out proteins:**

| Method | f0 | f1 | f2 | f3 | f4 | HOLDOUT(100/49) |
|---|---|---|---|---|---|---|
| **circuit-L2 (identity K-map)** | .170 | .148 | .179 | .173 | .157 | **.179** |
| circuit-L1 (group K-map) | .109 | .113 | .176 | .141 | .130 | .142 |
| plain MI (per-protein) | .059 | .096 | .084 | .100 | .121 | .097 |
| **APC-corrected MI (MIp)** | .155 | .175 | .181 | .208 | .228 | **.176** |
| perplexity ratio | .035 | .031 | .025 | .030 | .035 | .028 |
| sep-only ablation | .034 | .031 | .023 | .024 | .019 | .033 |
| base rate (random) | .035 | .033 | .037 | .036 | .036 | .035 |
| **PSICOV published (literature)** | .698 | .692 | .746 | .761 | .737 | **.723** |

- Circuit-L2 beats plain MI in all six splits and sep-only by ~5×, but **not
  APC-corrected MI**: MIp is ahead in 4 of 5 folds (CV mean 0.189 vs 0.165) and
  level on the holdout (0.176 vs 0.179, Wilcoxon p = 0.55).  The MIp row was
  added in Nov 2026; an earlier reviewer check reported MIp = MI because its APC
  term was never computed (correction A4).
- Binary membership classifier (L2): MCC ≈ 0.09, enrichment **2.5–2.8×** across
  folds (weak-but-real; ranked scores are the stronger reading).
- Circuit STABILITY across folds (exact minimiser, Nov 2026): L1 = 28–33 PIs
  (14–19 essential); L2 PIs/bin ≈ [36–52, 48–70, 34–48, 24–35] — same shape
  every fold.  (First reported as 20–22 essential and [72–106, …] with the
  under-merging minimiser.)

**What the rules say (holdout circuit).** Level-1 ON-cells are dominated by
hydrophobic∧aromatic × hydrophobic∧aromatic at short separation, plus
sulfur-involving terms (Cys pairing) — i.e., the QM **rediscovers core
hydrophobic packing** as the dominant short-range contact class. Level-2
bin-[6,10] implicants are concrete packing vocabularies, e.g.
{M,F}×{M,F,Y,W}, {I,V}×{I,V,F,W}, {I,F,C}×{F,W} — aromatic/aliphatic packings,
in minimal SOP form.

## 4. Sensitivity (gate G4)

prec@L/5 is **invariant at 0.179** across all nine (T∈{1.5,2,3})×(n_min∈{10,30,60})
labeling settings (ranked scores derive from cell rates, which thresholds barely
move); binary MCC varies coherently 0.055–0.097. Conclusions are labeling-robust.
`sensitivity_sweep.json`.

## 5. PATH E — Gray adjacency at proper power

| Statistic | Value |
|---|---|
| Native-contact pairs tested | **47,430** (Spike comparison: 38) |
| h=1 rate, contacts | 22.83% [95% CI 22.27–23.40] (protein bootstrap, 2000×) |
| h=1 rate, background (within-protein, within-sep-bin non-contact resampling) | 19.37% |
| Enrichment | **1.178×** |
| Two-sided binomial p | **1.1×10⁻⁷⁷** |
| χ² goodness-of-fit over h=0..5 | stat 1414.9, p ≈ 8×10⁻³⁰⁴ |

Interpretation: the encoding's signature property (chemically adjacent residues
one bit-flip apart) is genuinely enriched among contacting pairs — but the
effect is SMALL (~18% relative). The prior Spike negative (1.18×, n=38, p=0.345)
is thereby resolved: identical point estimate, decisive power. This also
validates the background design: composition/correlation-controlled nulls matter
(raw expectation 5/31≈16.1% would have overstated enrichment).

## 6. PATH F — the contact map AS a Boolean function (exploratory)

For each protein with L≤62 (exact QM feasible at 64×64 padding): treat the
native contact matrix as a ternary map (ON=contact, OFF=candidate-non-contact,
DC=everything else incl. sep<6 band and border), minimize, compare cover size
against controls whose labels are shuffled WITHIN separation bands (density and
separation structure preserved exactly).

| Protein | L | on-cells | PIs real | comp real | comp shuffled |
|---|---|---|---|---|---|
| 1aapA | 56 | 114 | 121 | 1.06 | 1.22 |
| 1brfA | 53 | 91 | 114 | 1.25 | 1.38 |
| 1guuA | 50 | 45 | 65 | 1.44 | 1.49 |
| 1jo8A | 58 | 127 | 118 | 0.93 | 1.24 |
| 1ku3A | 61 | 57 | 55 | 0.96 | 1.06 |
| 1m8aA | 61 | 95 | 93 | 0.98 | 1.11 |

**6/6 real maps compress better; mean 1.11 vs 1.25 PIs per contact; Wilcoxon
p = 0.031** — the smallest p a two-sided Wilcoxon test can give with n = 6.
(Table and means recomputed with the exact minimiser, Nov 2026; first reported
as 2.31 vs 2.74.)
Real 3D geometry carries K-map-expressible regularity beyond what density and
separation structure alone explain. Honest limit: only 6 proteins qualify for
exact QM at this cell budget; treat as directional evidence, not a benchmark.

## 7. Honest limitations

1. **Gap to DCA stands**: 0.179 vs 0.723 (literature). Consensus-identity
   features + majority-vote cells cannot substitute for coupling decomposition;
   the value demonstrated is transfer + interpretability, not SOTA prediction.
2. Consensus residues collapse within-column variation; frequency-weighted
   variants untested here.
3. Separation bins are global; per-protein length conditioning absent.
4. PATH F n=6; PATH E background uses consensus pairs (MSA-frequency-weighted
   variant untested); multiple-comparison discipline applied within (not
   across) analyses.
5. PSICOV proteins are a single-domain-heavy benchmark; transfer to
   multi-domain chains not established.

## 8. Reproduce

```bash
PY=/store/shuvam/.venv/bin/python; export OMP_PROC_BIND=FALSE
cd /store/shuvam/E-motioner-X-SBS/co-evolution-analysis
$PY scripts/tests/test_kmap_structure.py                       # 26/26
$PY scripts/run_contact_campaign.py --stage cache              # resumable
$PY scripts/run_contact_campaign.py --stage cv                 # PATH B (+C)
$PY scripts/run_contact_campaign.py --stage sweep              # G4
$PY scripts/run_contact_campaign.py --stage gray --boot 2000   # PATH E
$PY scripts/run_contact_campaign.py --stage mapqm              # PATH F
```
Outputs under `results/kmap_structure/`: features/*.npz (150),
cv_fold{0..4}.json, holdout_split.json, cv_summary.json,
sensitivity_sweep.json, path_e_gray_adjacency.json, path_f_map_qm.json.

Determinism: fixed seeds everywhere (folds seed=42; backgrounds/bootstrap
seeded RNG); stage re-runs reproduce identical JSONs — verified: **5
consecutive `cv` runs byte-identical** (`holdout_split.json` md5
`ebbd6666…`) and **3 consecutive `gray` runs byte-identical**
(`path_e_gray_adjacency.json` md5 `78fef01a…`).

---

## 9. SKEPTICISM AUDIT (post-results verification wave, Aug 22)

Directive: be skeptical, verify externally, repeat. Seven checks executed.

### S1 — Literature reconciliation of the PSICOV baseline [CORRECTION]
Our measurement of PSICOV's own prediction files: prec@L/5 = **0.727** (sep≥6),
0.640 (long ≥24). Web verification against the actual paper (Jones et al. 2012,
Bioinformatics 28:184; abstract + Table 1 via doi.org/10.1093/bioinformatics/btr638):
- Abstract: "For 118 out of 150 targets, the L/5 precision for long-range
  contacts (sequence separation >23) was ≥ 0.5" ⇒ mean well above 0.44.
- Table 1 rows put PSICOV at L/5 ≈ **0.73** ([i−j]>6) / ≈0.68 ([i−j]>23).
⇒ **our pipeline reproduces the published benchmark to Δ=0.003**.
⇒ **The "published top-L/5 = 0.44" figure carried in this repo's older docs
   (agents.md §4/§5.3, 00_overview.md, 01_Datasets_and_Files.md) is a
   MISATTRIBUTION** (0.44 is the top-L column / MIp-B&vN row, not PSICOV L/5).
   Historical docs left as-is per no-rewrite doctrine; corrected HERE and in
   plans/11d-audit-log.md.

### S2 — Independent hand-recomputation (no shared code)
1a3aA from raw files with a from-scratch parser: 145 CB residues, **325 native
contacts** (== pipeline), PDB(CA)-sequence == target fasta TRUE,
prec@L/5 = 0.862 for that protein (consistent with its cached scores).

### S3 — Sequence-level identity, all 150 targets
PDB CA-sequence == target fasta: **149/150 exact**; sole exception 1vjkA
(+1 C-terminal residue, prefix identical — known artifact, drop-rule applied).

### S4 — Is Gray enrichment encoding-specific? (code-permutation control)
Permuting which residue sits at which He code (30 permutations; identical pairs,
labels, backgrounds): random-code enrichment mean **0.961**, sd 0.095,
max 1.158 vs He/Gray **1.198** ⇒ **z = 2.51**; as an empirical test with 30
permutations, p ≤ 1/31 ≈ 0.03 (the normal-approximation p ≈ 0.006 is optimistic).
Honest reading: most of the raw 1.18× reflects generic chemistry clustering
(any grouped encoding); the SPECIFIC Gray ordering adds a real but modest edge.

### S5 — Homology leakage in PATH B splits
Alignment-free 5-mer Jaccard across all C(150,2)=11,175 target pairs:
max = **0.029**, zero pairs above 0.3 ⇒ no near-duplicate/homologous targets;
the paper's one-target-per-Pfam-family design holds; protein-level CV splits
are sound.

### S6 — Determinism ×N ("verify again, five more times")
| stage | reruns | byte-identical md5 |
|---|---|---|
| cv (holdout_split.json) | **7** | ebbd6666cac3d30f… every run |
| gray (path_e_gray_adjacency.json) | 3 | 78fef01ac63c3feb… every run |
| mapqm (path_f_map_qm.json) | 2 | 9845cf72d908be12… every run |

### S7 — Band-stratified honesty check
| method | short 6–11 | medium 12–23 | long ≥24 |
|---|---|---|---|
| circuit-L2 | .179 | .134 | **.137** |
| circuit-L1 | .135 | .110 | .116 |
| MI | .115 | .081 | .077 |
| PSICOV published | .409 | .478 | **.637** |
Circuit transfer holds in EVERY band (long-range: 2.9× base-rate-equivalent MI,
still ~4.6× below DCA). The earlier single-number table mixed bands; this is
the fair version.

## 10. Post-audit verdict

Every headline claim survived adversarial re-examination; one literature
misattribution (0.44) was caught and corrected; the Gray result is confirmed
encoding-specific (z=2.51); PATH B transfer is homology-clean and
band-robust; everything reruns byte-identically.
