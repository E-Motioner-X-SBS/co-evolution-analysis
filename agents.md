# agents.md — Complete Agent Guide: Co-evolution Analysis (E-Motioner-X-SBS)

> This file is the single entry point for any agent (or human) starting work on
> the E-Motioner-X-SBS co-evolution project. It covers: the organization, the
> environment, the scientific framework, every repository to read, every
> dataset (with verified links and papers), every script, the 3D structure
> analysis, the multi-day runs, the interpretations, and exact reproducibility.
> Written 2026-08-14; verified live this session (environment, links, scripts).

---

## 0. READ THIS FIRST — the org at a glance

**Organization root:** `/store/shuvam/E-motioner-X-SBS/`
**Every subdirectory under this root is a SEPARATE GIT REPOSITORY** (each has
its own `.git/`). Verify with `git -C <folder> status`.

**The core thesis:** biological sequences can be represented as **Karnaugh
maps (K-maps)** via **Gray-code encoding**, connecting digital logic design to
structural biology. Co-evolutionary constraints in a protein family = a Boolean
function over residue codes, minimized by **Quine–McCluskey** → prime
implicants → rules.

**Repositories (all under the org root):**

| Repo | Contents | Read for |
|---|---|---|
| `co-evolution-analysis/` | **THIS repo** — the 23-script co-evolution pipeline + all validation + 3D structural analysis | everything |
| `lean_proofs/` | Lean 4 formal proofs of the encoding framework — **106 theorems, zero `sorry`** | the mathematical foundation |
| `n-ary-kmap/` | Base-N K-map generalization — 115 Lean theorems; `src/nkmap/` encoders | base-20 encoding |
| `kmap-sbm-validation/` | Empirical validation via structure-based models (SMOG/SBM MD); `src/kmap_sbm/` (Gray encoding + Quine–McCluskey) | encodings + QM + MD plans |
| `datasets/` | Curated PDB structures; `datasets/co-evolution/` = working copy of the analysis | raw data |
| `skills/` | Agent onboarding docs + small tools (lit_review_adder, gh_org_manager, latex_tikz_diagrams…) | org conventions |
| `contact_mapping/`, `pdb_hunter/`, `ml_folding_models/`, `Is-Kmap-Possible/`, `KMAP-rethink/`, `e-motioner-x-sbs.github.io/` | related exploratory projects | context (optional) |

---

## 1. Environment (exact, verified 2026-08-14)

### 1.1 The Python environment — USE THIS ONE

```bash
PY=/store/shuvam/.venv/bin/python
```

Verified contents: **Python 3.13.5, torch 2.12.1+cu130 (CUDA available),
numpy 2.4.6, scipy 1.18.0, matplotlib 3.11.0, requests, numba 0.66.0, pyflakes.**

- GPU: **NVIDIA A100 80 GB PCIe** (81,920 MiB total).
- `uv` 0.9.21 is available at `/home/roy/.local/bin/uv`, but **its cache is
  broken** (permission error on `/store/uv`) — for NEW packages use
  `$PY -m pip install <pkg>` (pyflakes was installed this way).

### 1.2 REQUIRED environment variables (do not skip)

```bash
export OMP_PROC_BIND=FALSE
# CRITICAL: without it, `import torch` pins the process to CPU 0 only
# (verified: affinity becomes {0}; with it, all 24 CPUs are usable).
# This single bug made "multi-core" runs run on ONE core.
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
# avoids torch caching-allocator fragmentation deadlocks on the shared GPU
```

### 1.3 sys.path requirements

The analysis scripts insert the sibling-repo paths themselves, but any NEW
script importing them must add:

```python
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/kmap-sbm-validation/src")
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")
```

### 1.4 The machine is SHARED

~45 users, load often > 24 cores. Plan for contention; jobs get ~5–100% of a
core depending on time of day. The GPU is shared too — the pipeline caps
per-process GPU memory at 10% (see `coevolution_gpu.get_device()`).

---

## 2. The scientific framework — read these in order

1. **`skills/ORGANIZATION.md`** — org overview, results summary, cheat sheet.
2. **`lean_proofs/`** — the 106 Lean theorems: Gray-code properties
   (`KmapProofs.lean`), amino-acid 5-bit Gray encoding with 7 physicochemical
   groups (`AminoAcidEncoding.lean`), k-mer indexing, contact-map completeness.
   Build: `lake build` (needs
   `LD_LIBRARY_PATH=/home/roy/.elan/toolchains/leanprover--lean4---v4.29.0/lib/lean`).
3. **`n-ary-kmap/`** — the base-N generalization (115 theorems); the base-20
   encoder lives in `src/nkmap/encoding/bio_sequences.py`
   (`AMINO_HE_2012`: A=0, I=1, L=2, V=3, M=4, F=5, Y=6, W=7, E=8, D=9, Q=10,
   N=11, H=12, K=13, R=14, S=15, T=16, C=17, P=18, G=19).
4. **`kmap-sbm-validation/src/kmap_sbm/`** — `encoding/gray_amino.py`
   (5-bit Gray per residue, `encode_gray_single`, `build_aa_kmap_2d` — the
   32×32 dipeptide K-map used for H3) and `analysis/prime_implicants.py`
   (the Quine–McCluskey implementation; **must receive power-of-4 cell
   counts — it raises otherwise**).

### 2.1 The core pipeline (one line each)

```
FASTA (aligned) → base-20 (He 2012) encoding, gaps/ambiguous = state 20
→ per-position Shannon entropy H(p) (gaps excluded) → variable positions (H > 0.3)
→ per-position-pair mutual information MI(i,j) (window 30; mutation-only variant
  excludes the majority/reference pair)
→ per-pair 20×20 frequency K-map → ternary mutation K-map (1 observed mutation,
  -1 reference don't-care, 0 never observed) → 32×32 zero-padded (10 bits)
→ Quine–McCluskey → prime implicants → essential rules
→ constraint function C = ln(P/P_exp); prediction σ(C) = 1/(1+e^{-C})
→ flipped K-maps (forbidden pairs), n-ary K-maps, perplexity ratio,
  DCA-style analysis, network/clustering — 23 scripts total
```

**The A1/A2 corrections (MUST know — see
`My_Interpretation_of_the_Results/CORRECTION_NOTICE.md`):**
- **A1 gap-stripping**: never strip gaps from the alignment — keep them as
  state 20 and exclude from counts, else "column j" mixes different raw
  positions across sequences (produced 1,249 fake variable positions).
- **A2 QM wrap-around**: never feed a 20×20 (400-cell) map to the 8-bit QM —
  pad to 32×32 (1,024 cells, 5 bits/axis) with don't-care rows/cols 20–31.
  The old bug produced 143/152 phantom rules; the corrected set is
  **36 distinct prime implicants, 2 essential**.

---

## 3. This repository (`co-evolution-analysis/`) — the full map

### 3.1 Top-level layout

```
co-evolution-analysis/
├── *.py                         # the 23 analysis scripts (parameterized)
├── Spike_protein.aln-fasta      # 1,299 Omicron Spike sequences (aligned, 1,276 cols)
├── swissmodel/                  # SWISS-MODEL models + batch tooling (see §6)
├── data/                        # GITIGNORED: dataset downloads + conversions
│   ├── psicov150/  pfam/  evmutation/  unit_tests/  deepmsa/  gpcrdb/  structures/
│   ├── metadata.json, manifest.json, download.log
├── results/                     # GITIGNORED: ALL script outputs
│   ├── <dataset>/<target>/<script_dir>/   (script outputs + run_log.json)
│   ├── contacts/                # all 3D-structure analysis JSONs
│   ├── consistency/  coverage_report.json
├── My_Interpretation_of_the_Results/      # the ORIGINAL Spike analysis (23 docs)
├── My_Own_Interpretation_Across_New_Datasets/   # the multi-dataset analysis
│   ├── 00_overview.md, 01_3d_structural_validation.md
│   ├── <dataset>/NN_<script>.md            # per-script interpretations
│   └── 3D_Structure_Analysis/              # ★ the full 3D documentation (§7)
├── plans/                       # planning docs (understanding/research/strategy/steps/decisions/audit)
├── scripts/                     # tooling: downloader, converters, runner, checkers, 3D analysis
├── agents.md                    # THIS FILE
└── README.md                    # repo overview + results table
```

### 3.2 The 23 analysis scripts — what each computes

| # | Script | Purpose |
|---|---|---|
| 1 | `coevolution_shared.py` | shared module: FASTA parsing, entropy, MI, coupling, GPU pair-finder, combined MI+perplexity (single source of truth) |
| 2 | `coevolution_gpu.py` | GPU kernels (torch CUDA): MI matrix, entropy, majority refs, H1, coupling |
| 3 | `master_boolean.py` | master Boolean function → 36 PIs / 2 essential (the flagship rule extractor) |
| 4 | `boolean_co-evolution.py` | whole-protein dipeptide Boolean minimization + coupling constants |
| 5 | `nary_kmap_co-evolution.py` | base-20 (n-ary) K-map analysis |
| 6 | `position_kmap_coevolution.py` | per-position-pair K-maps with MI |
| 7 | `run_allseq_analysis.py` | all-sequence MI matrix (window 30) |
| 8 | `run_kmap_analysis.py` | binary Gray K-map pipeline: H1 adjacency (1.34×), signatures, Walsh–Hadamard, mutations |
| 9 | `flipped_boolean_coevolution.py` | flipped (forbidden) K-maps → 490 forbidden rules |
| 10 | `kmap_boolean_coevolution.py` | K-map Boolean equations (5-bit Gray literals) |
| 11 | `variable_position_coevolution.py` | variable-positions-only K-maps with don't-cares |
| 12 | `predictive_constraint_function.py` | constraint function train/test prediction |
| 13 | `allseq_constraint_function.py` | leave-one-out CV with the constraint function |
| 14 | `dca_boolean_coevolution.py` | local precision matrix → Boolean (NOT real DCA) |
| 15 | `dca_mf_analysis.py` | proper mfDCA (Frobenius norm, APC, direct information) |
| 16 | `perplexity_coevolution.py` | perplexity ratio per pair |
| 17 | `advanced_co-evolution_analysis.py` | network, clustering, Walsh–Hadamard, variant signatures |
| 18 | `full_length_analysis.py` | full-length entropy + MI (window 30) |
| 19 | `gpu_full_analysis.py` | GPU full analysis: full MI matrix, entropy, refs, H1 |
| 20 | `create_mi_heatmap.py` | MI heatmaps |
| 21 | `generate_co-evolution_md.py` | markdown of the Boolean functions |
| 22 | `generate_full_analysis_md.py` | comprehensive report |
| 23 | `generate_full_pipeline_doc.py` | full pipeline documentation |

### 3.3 Parameterization (how every script is pointed at any dataset)

Every script reads two environment variables (defaults unchanged = the
original Spike analysis):

```bash
export COEVO_FASTA=/path/to/any.fasta      # input alignment
export COEVO_RESULTS=/path/to/output_dir   # outputs
export COEVO_OUT=/path/to/report.md        # report generators only
```

- Report generators read the sibling result dirs at the `COEVO_RESULTS` level
  (they were fixed to be tolerant of missing dependencies — SafeDict).
- GPU-heavy scripts (8 of them) share a semaphore (max 4 concurrent); the
  LOO-CV scripts (`allseq_constraint_function.py`,
  `predictive_constraint_function.py`) are **forced off CUDA entirely**
  (they are CPU-bound and their GPU MI call caused allocator spins).
- `coevolution_gpu.get_device()` caps per-process GPU memory at 10% and
  `mi_matrix_gpu` uses **adaptive chunking** (flat int64 buffer ≤ 512 MB) so
  deep alignments (up to 486k sequences) never OOM.

### 3.4 The runner (how the multi-day runs were executed)

```bash
$PY scripts/run_dataset.py --dataset <ds> --psicov-mode full --parallel 12
# two waves: analysis scripts first, report generators after
# --retry-failed      : quarantine failed dirs and re-run them
# --retry-failed-cpu  : same, but GPU-heavy scripts forced to CPU
# --psicov-mode core  : the 11 core scripts (psicov150 150 proteins)
```

Every run writes `run_log.json` (status, seconds, stdout/stderr tail, output
file list) — **3,687+ verified runs, all auditable**.

---

## 4. The datasets (all verified live 2026-08-13)

| Dataset | Paper (DOI — verified resolves) | Download URL (verified) | Local data | Format |
|---|---|---|---|---|
| **Omicron Spike (1,299)** | GISAID; SWISS-MODEL: Waterhouse 2018, 10.1093/nar/gky427 | — | `Spike_protein.aln-fasta` | aligned FASTA, 1,276 cols |
| **PSICOV 150** | Jones et al. 2012 *Bioinformatics* 28:184, 10.1093/bioinformatics/btr638 | `https://bioinfadmin.cs.ucl.ac.uk/downloads/PSICOV/suppdata/suppdata.tar.gz` | `data/psicov150/` (aln/, seq/, con/, **pdb/** natives) | tar.gz; MSAs headerless; PDBs blank-chain |
| **Pfam PF00072/00595/00071** | Morcos et al. 2011 *PNAS* 108:E1293, 10.1073/pnas.1111471108; Mistry 2021 NAR 49:D412 | `https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00072?annotation=alignment:full` (+00595, 00071) | `data/pfam/*_full.sto` + `converted/` | Stockholm (486k/330k/222k seqs) |
| **EVmutation 20** | Hopf et al. 2017 *Nat Biotechnol* 35:1026, 10.1038/nbt.3969 | `https://marks.hms.harvard.edu/evmutation/supp/alignments.tar.gz` (**needs `curl -k`**) | `data/evmutation/` | A2M (lowercase = inserts) |
| **Unit tests** | Hopf 2019 plmc 10.1093/bioinformatics/bty657; Kamisetty 2013 GREMLIN 10.1073/pnas.1307869110; Seemayer 2014 CCMpred 10.1093/bioinformatics/btu171 | GitHub raw: `debbiemarkslab/plmc` (DHFR.a2m), `sokrypton/GREMLIN` (1whzA), `soedinglab/CCMpred` (1atzA.aln) | `data/unit_tests/` | A2M / two-column / headerless |
| **DeepMSA** | Zhang et al. 2020 *Bioinformatics* 36:2105, 10.1093/bioinformatics/btz978 | `https://zhanggroup.org/DeepMSA/dataset.tar.bz2` | `data/deepmsa/` (614 targets: seq.fasta + **native.pdb**) | tar.bz2 — **no MSAs shipped** |
| **GPCRdb** | Kooistra et al. 2021 *NAR* 49:D335, 10.1093/nar/gkaa1080 | `https://gpcrdb.org/alignment/fasta/001_001_001_001/` (5-HT1A), `.../001_001_003_008/` (β2AR) | `data/gpcrdb/` | FASTA (125/136 seqs) |
| **Structures (RCSB)** | — | `https://files.rcsb.org/download/<ID>.pdb` (6TNE, 6AK2, 1HE1, 7E2X, 2RH1, 4P3R, 1WHZ, 1ATZ) | `data/structures/` | PDB |

Notes:
- Journal DOI links return HTTP 403 to curl (publisher bot-blocking) but
  RESOLVE correctly — valid in a browser.
- Conversion pipeline: `scripts/convert_all.py` (Stockholm/A2M/PSICOV→FASTA,
  metadata + redundancy stats); `scripts/download_datasets.py` (integrity
  verified, manifest.json); `scripts/fetch_structures.py`.
- **Redundancy reality**: the 1,299 Omicron entries contain only 299 unique
  sequences (587 identical copies of one variant) — the effective sample for
  co-evolution is 299. This is WHY the validation datasets (diverse MSAs)
  matter.

---

## 5. The 3D structure work (the core of the validation)

### 5.1 The models — 1,128 PDB files

- Built with **SWISS-MODEL** (Automodel mode; token in
  `swissmodel/.env`, gitignored) for the **299 unique Omicron sequences**:
  `swissmodel/batch/pdbs/u####/` — 3–8 models each (median 3) → **1,128
  `.pdb.gz` + 1,128 `.cif.gz`**, plus per-sequence `full_details.json`
  (GMQE, QSQE, **target_sequence** — the exact modelled sequence), and
  `templates_summary.json` (template SMTL id, identity, coverage).
- Best-model GMQE 0.69–0.72 for all 299 (all > 0.6); mostly homo-trimers.
- Tooling: `swissmodel/swissmodel_batch.py` (submit→poll→download,
  rate-limited; 299 jobs, 0 failed) and `swissmodel/verify_*.py`.
- **Native experimental structures** (PSICOV 150) are the independent
  ground truth.

### 5.2 How contacts are computed (summary — full detail in
`My_Own_Interpretation_Across_New_Datasets/3D_Structure_Analysis/02_Contact_Computation.md`)

- PDB parse: ATOM records; Cβ coords (Cα for Gly); chain ids (blank chain =
  chain A for PSICOV); residue numbering per chain.
- **Contact** = Cβ–Cβ distance < 8 Å, residue separation ≥ 3 (rules) or ≥ 6
  (DCA convention).
- **Position mapping**: alignment column → model residue = col + 1 −
  (non-canonical chars before col). Integrity check: the model's
  `target_sequence` must equal the gap/ambiguous-stripped aligned sequence;
  mismatches are skipped, never silently mapped (0 skips after the X-stripping
  fix).
- **Controls**: matched separation distribution, 2,000 random draws/model,
  fixed seeds; exact binomial p (scipy `binomtest`, one-sided greater).

### 5.3 The 9 structural analyses (scripts in `scripts/`, results in
`results/contacts/`, full write-up in `3D_Structure_Analysis/04_Structural_Analyses.md`)

| # | Script | Question | Headline result |
|---|---|---|---|
| 1 | `validate_3d_complete.py` | Do co-evolving pairs from each of 8 sources sit close in 3D? | position_kmap 2.08× (p=1e-103); (407,410) 9.3×; master_boolean 0.33× |
| 2 | `validate_residue_rules_3d.py` | Do the 36 prime implicants fire on contacts? | [210,I,215,P] 2.85× (p=1e-6); pooled 0.58× |
| 3 | `validate_flipped_3d.py` | Are forbidden pairs sterically constrained? | (212,215) 2.71×; pooled rarity |
| 4 | `gray_flip_analysis.py` | Is Gray bit-flip adjacency predictive? | **1.18×, p=0.345 — NO** |
| 5 | `structural_profile.py` | Housekeeping: conserved buried, variable exposed? | **ρ=−0.17, p=6e-9 — YES** |
| 6 | `trimer_interface.py` | Any inter-protomer contacts? | **ZERO** |
| 7 | `consistency_cross_script.py` | Does cross-script agreement imply structure? | **NO** — 5-source pairs non-structural |
| 8 | `perplexity_3d_analysis.py` | Is the perplexity ratio structural? | **ρ=0.67, p=0.006 — strongest metric** |
> **CORRECTION (Aug 22, 2026):** this figure is a misattribution. Jones et al. 2012 Table 1 gives PSICOV top-L/5 ≈ **0.73** (sep>6) / ≈0.68 (>23); our re-measurement of the shipped con files = 0.727 / 0.640 (see 3D_Structure_Analysis/11 §9-S1). The 0.44 value corresponds to top-L / MIp-B&vN rows.
| 9 | `validate_contacts.py` | PSICOV 150 natives: enrichment + precision? | 2.49× mean enrichment; precision 0.066 vs published 0.44 |

### 5.4 The interpretations (where the findings live)

- `My_Own_Interpretation_Across_New_Datasets/00_overview.md` — multi-dataset
  overview (datasets, methods, honest results).
- `My_Own_Interpretation_Across_New_Datasets/01_3d_structural_validation.md`
  — the rule→3D synthesis.
- `My_Own_Interpretation_Across_New_Datasets/3D_Structure_Analysis/` — **the
  6-document full write-up**: 00_Overview, 01_Datasets_and_Files,
  02_Contact_Computation, 03_Kmaps_and_Boolean_Rules,
  04_Structural_Analyses, 05_Reproducibility.
- `My_Own_Interpretation_Across_New_Datasets/<dataset>/NN_<script>.md` —
  per-script interpretations (3,037 docs), generated by
  `scripts/generate_interpretation.py`; rules/tables by
  `scripts/generate_rules_summary.py`.

---

## 6. The multi-day runs — what was actually executed

1. **The 23-script suite × 5 validation datasets** (unit_tests 63/63, gpcrdb
   42/42, evmutation 420/420, psicov150 3,150/3,150 — full suite; pfam 51/63
   with the remaining jobs documented, see below) = **3,687+ runs** over days
   on the A100 + 24-core machine.
2. **The Spike analysis** (all 23 scripts on the 1,299 sequences; corrected
   results in `datasets/co-evolution/` result folders and
   `My_Interpretation_of_the_Results/`).
3. **SWISS-MODEL batch**: 299 automodel jobs → 1,128 PDBs (submission ~13
   min, total ~75 min; rate limits respected: 0 × 429).
4. **9 structural analyses** on the 299 models (this session's final work).

### Known computational limits (HONEST — do not "fix" by killing)

| Job | Why it cannot finish at full scale |
|---|---|
| `allseq_constraint_function.py` (LOO-CV) on 222k–486k-seq Pfam families | O(10·N²): ≈27 days (222k) to 2.5 years (486k) of single-core compute. Runs per the no-timeout directive; the only practical route is **documented sampling** (e.g., ≤2,000 test sequences, clearly labeled) |
| `dca_mf_analysis.py` / `gpu_full_analysis.py` on the Pfam giants | GPU OOM at the 8 GB cap; CPU would need 89 GB RAM for the one-hot at 486k sequences |
| `run_kmap_analysis.py` H3 on huge alignments | fixed via GPU chunked norm-trick (the scipy pdist fallback attempted an 883 GB allocation — replaced with a chunked CPU fallback) |

---

## 7. Reproducibility — exact commands

```bash
PY=/store/shuvam/.venv/bin/python
export OMP_PROC_BIND=FALSE
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
cd /store/shuvam/E-motioner-X-SBS/co-evolution-analysis

# 0. sanity
$PY -c "import torch; print(torch.__version__, torch.cuda.is_available())"

# 1. reproduce the corrected Spike rule set (regression-verified: 36/2)
$PY master_boolean.py

# 2. reproduce the rule inventory for the structural analyses
$PY scripts/extract_rules_inventory.py

# 3. reproduce any structural analysis (deterministic: fixed seeds)
$PY scripts/validate_3d_complete.py
$PY scripts/validate_residue_rules_3d.py
$PY scripts/validate_flipped_3d.py
$PY scripts/gray_flip_analysis.py
$PY scripts/structural_profile.py
$PY scripts/trimer_interface.py
$PY scripts/consistency_cross_script.py
$PY scripts/perplexity_3d_analysis.py
$PY scripts/validate_contacts.py

# 4. audit everything
$PY scripts/verify_coverage.py          # 0 truncation flags expected
$PY scripts/check_completeness.py       # every dataset×target×script accounted
$PY scripts/check_gpu_cpu.py --all --max-seq 5000   # 5/5 PASS expected

# 5. regenerate the interpretation docs
$PY scripts/generate_rules_summary.py
$PY scripts/generate_interpretation.py

# 6. re-run any dataset's scripts (parameterized)
COEVO_FASTA=data/gpcrdb/converted/b2ar.fasta \
COEVO_RESULTS=/tmp/example \
  $PY master_boolean.py
```

Determinism guarantees: fixed random seeds per model; exact scipy statistics;
byte-compare-verified re-runs; coverage audit 0 flags on 3,732 runs.

---

## 8. Key honest findings (the science in one table)

| Question | Answer | Evidence |
|---|---|---|
| Is the co-evolution signal structural? | **Partly** — 6 pairs at 3.4–10.5× (p ≤ 1e-8): (500,503), (503,507), (407,410), (210,214), (212,215), (210,215) | validate_3d_complete |
| Where does it fail? | The 442–498 cluster and (373,378): 13–27 Å apart — immune/lineage covariation | validate_3d_complete |
| Prime implicants? | 210/212–215 rules fire on contacts (2.8×); pooled 0.58× | validate_residue_rules_3d |
| Flipped rules? | Steric at 210/212–215; rarity overall | validate_flipped_3d |
| Gray bit-flip? | NOT predictive (1.18×, p=0.345) | gray_flip_analysis |
| Housekeeping? | YES — conserved buried (52% vs 26%), ρ=−0.17, p=6e-9 | structural_profile |
| Trimer interface? | Zero inter-chain contacts | trimer_interface |
| Consistency ⇒ structure? | NO | consistency_cross_script |
| Perplexity? | ρ=0.67, p=0.006 — best single metric | perplexity_3d_analysis |
| vs published DCA? | enrichment 2.49× (PSICOV 150); precision 0.066 vs 0.44 [CORRECTED: PSICOV L/5 is 0.73, not 0.44 — see 11_Kmap doc §9] | validate_contacts |
| Data quality | GPU/CPU identical to 1e-7; coverage 0 flags; byte-compare identical | check_gpu_cpu, verify_coverage |

---

## 9. Troubleshooting (learned the hard way)

| Symptom | Cause / fix |
|---|---|
| All processes on CPU 0, no multi-core | `OMP_PROC_BIND=FALSE` missing (torch pins to CPU 0) |
| torch OOM / allocator deadlock (0% GPU, memory held) | `expandable_segments:True` + memory fraction 0.1 + adaptive chunking |
| "Unable to allocate 883 GiB" in run_kmap | H3 pdist fallback — replaced with chunked GPU/CPU norm-trick |
| `kmap_truth_table` raises on 400-cell maps | A2: pad 20×20 → 32×32 before QM (by design) |
| 1,249 "variable positions" / MI 1.59 | A1: you stripped gaps — use `load_position_arrays(aligned=True)` (default) |
| 143 phantom rules | A2: 8-bit QM wrap-around — use padded maps |
| LOO-CV runs forever on Pfam giants | documented infeasible (O(N²)); do not kill; sampling is the only practical route |
| `marks.hms.harvard.edu` TLS error | `curl -k` (already in download_datasets.py) |
| uv install fails (permission on /store/uv) | use `$PY -m pip install` |
| Report generator KeyError on missing data | fixed with SafeDict; regenerate |
| Coverage audit flags "SEQS n<total" | check it's a train/test or subset line (parser handles totals) |

---

## 9.5 K-map ↔ 3D-Structure Encoding Campaign (Aug 22, 2026 — NEW)

Literature-ground-truth verification that the Gray→K-map→QM chain encodes 3D
structure. Full write-up:
`My_Own_Interpretation_Across_New_Datasets/3D_Structure_Analysis/11_Kmap_Structure_Encoding.md`.

| Asset | Path |
|---|---|
| Library (contact K-maps, native-PDB labels, PSICOV-con parsing) | `scripts/kmap_structure.py` |
| Stage runner: cache → cv → sweep → gray → mapqm | `scripts/run_contact_campaign.py` |
| Adversarial tests (26) | `scripts/tests/test_kmap_structure.py` |
| Results | `results/kmap_structure/*.json` + `features/*.npz` (150) |
| Lean proofs (12 theorems, axiom-audited) | `../lean_proofs/proofs/is_kmap_possible/ContactCircuits.lean` |

Key numbers: circuit transfer prec@L/5 **0.179** holdout (MI 0.097, base 0.034,
PSICOV-pub 0.727); Gray enrichment **1.178×, p=1.1e-77** (n=47,430; permutation
control z=2.51); map-compression 6/6 p=0.031; determinism md5-stable ×5.

**CORRECTION carried here:** the older "published PSICOV L/5 = 0.44" figure in

### Aug 22 extension — tension, coupling-in-K-map, reconstruction, complete circuit

| Result | Numbers |
|---|---|
| Tension mechanism (G): no lift | circ-G 0.154 < circ-L2 0.165; hubs anti-correlated |
| mfDCA-in-K-map (H): no lift from DI bins | circH 0.124 < L1 0.134 < L2 0.165 |
| Our mfDCA matches published benchmarks | 0.169 ≈ literature 0.15–0.20 |
| **Rank-fusion: K-map prior + DCA** | **0.204 prec@L/5, p≈0 vs both parents** |
| Reconstruction verdict (I) | circuits NOT ready (recall@L 2.8%); DCA IS (21.8%) |
| Complete circuit (K): flipped transfers | specificity 98.3%; cons/cons 3.2× enriched |
| Gray adjacency at scale | 1.178×, p=1e-77; permutation z=2.51 |

New scripts: `path_{g,h,i,k,l}_*.py`; new results: `results/kmap_structure/path_{g,h,i,k,l}*.json`.
§4/§5.3 above is a misattribution (paper Table-1 L/5 ≈ 0.73; we reproduce
0.727). See 11_Kmap doc §9-S1 for the web-verified reconciliation.

---

## 10. Agent behavior rules (this project)

1. **Read before running**: `plans/`, this file, the 3D_Structure_Analysis
   docs, and the CORRECTION_NOTICE before touching scripts.
2. **Never strip gaps, never shrink K-maps below 32×32 for QM** — the two
   historical defects.
3. **Verify links/facts live** before citing; tag
   [VERIFIED]/[HYPOTHESIS]/[UNKNOWN].
4. **Never delete**: move to `/store/shuvam/.rigor-trash/YYYYMMDD-<name>/` +
   INDEX entry.
5. **Every run logged**: `run_log.json` (exit, seconds, tail, outputs).
6. **No silent failures**: investigate (parsing? implementation?), fix,
   re-run, re-audit.
7. **GPU/CPU consistency** and **coverage audit** after any code change.
8. **Honest reporting**: negative results are results (Gray-flip p=0.345,
   trimer zero, pooled PIs 0.58× are findings, not failures).
9. Keep `plans/05-audit-log.md` current; update this file when the project
   changes.
10. The org uses `uv`/Python 3.13 conventions (`skills/README.md`), but the
    working venv is `/store/shuvam/.venv`.
