# 00 — Understanding: Multi-Dataset Validation of the K-map / Co-evolution Pipeline

## The task (user's words, restated)
1. Download ALL 6 validation datasets COMPLETELY (no segments): PSICOV150, Pfam (PF00072/PDZ/Ras), EVmutation20, unit-test sets, DeepMSA CASP12/13, GPCRdb.
2. Run the FULL co-evolution analysis (all 23 Python scripts from the co-evolution-analysis repo) on each dataset, in an INDEPENDENT new repository (NOT datasets/), with `data/` (per-dataset folders) and `results/` (per-dataset per-script outputs).
3. Generate EVERYTHING the scripts produce: K-map tables, prime implicants, Boolean minimized equations — for main K-map, n-ary, flipped, perplexity.
4. Validate against the structural ground truth the datasets carry (native PDBs, contact maps): do our co-evolution rules/pairs translate to 3D contacts?
5. GPU/CPU consistency checks (sample per dataset).
6. RIGID plan + LOOPING: download → run → verify results → re-run → investigate anomalies (parsing? implementation?) → fix → re-run. Never stop after one pass.
7. `agents.md` in the co-evolution repository: env setup (uv), how to run, dataset provenance, agent behavior.
8. `My_Own_Interpretation_Across_New_Datasets/` folder with all results + interpretations.
9. HONEST reporting: where the algorithms work and where they don't.

## Why validation is needed
The current Omicron dataset (1,299 seqs) is unfit: 587 identical copies → effective N ≈ 299,
one variant dominates → MI/entropy biased, 9 variable positions hidden. Validation datasets
must be diverse, deep, structurally characterized, and benchmarked.

## The 6 datasets (all URLs verified live Aug 8, 2026 — see 01-research.md)
| # | Dataset | Paper | Contents | Ground truth |
|---|---------|-------|----------|--------------|
| 1 | PSICOV 150 | Jones 2012 Bioinformatics | 150 MSAs (median 3,228 seqs) + native PDBs + PSICOV contact predictions | **native PDB + published top-L/5 precision (0.44)** |
| 2 | Pfam PF00072/PF00595/PF00071 | Morcos 2011 PNAS (DCA canon) | 486k/330k/222k seq Stockholm alignments | representative PDBs (6TNE/6AK2/1HE1, to fetch) |
| 3 | EVmutation 20 | Hopf 2017 Nat Biotech | 20 A2M alignments (up to 28k seqs) | structures + published mutation-effect benchmarks |
| 4 | Unit tests | plmc/GREMLIN/CCMpred | DHFR.a2m (3,629×80), 1whzA (662×78), 1atzA (1,534×75) | PDBs 4P3R/1WHZ/1ATZ |
| 5 | DeepMSA | Zhang 2020 Bioinformatics | CASP12/13 FM target MSAs (deep, remote homology) | CASP native structures |
| 6 | GPCRdb | Kooistra 2021 NAR | 5-HT1A (125 seqs), β2AR (136 seqs), curated structure-based | PDBs 7E2X/2RH1 |

## Success criteria
- [ ] All 6 datasets fully downloaded, verified (sizes/magic/tar members/seq counts) → manifest.json
- [ ] All converted to canonical FASTA + per-dataset metadata (n_seq, length, redundancy stats, PDBs)
- [ ] All 23 scripts parameterized (COEVO_FASTA / COEVO_RESULTS) with regression proof (defaults unchanged)
- [ ] Every script run per dataset (or documented subset for PSICOV 150), outputs present, no silent failures
- [ ] K-map tables + prime implicants + Boolean equations generated for main/n-ary/flipped/perplexity per dataset
- [ ] GPU/CPU consistency verified per dataset (sample)
- [ ] Contact-map validation: co-evolving pairs/rules vs native contacts; precision reported vs published
- [ ] Second verification pass (re-run) → byte-identical results on sample
- [ ] agents.md complete; My_Own_Interpretation_Across_New_Datasets populated
- [ ] Honest final report: what works, what doesn't, anomalies investigated

## Open questions
- PSICOV 150 × 23 scripts = 3,450 runs → full-scale definition (see 02-strategy; adaptive)
- A2M lowercase-insert handling (uppercase-all vs strip) → decide in Phase 2, document
- torch env: fresh uv venv vs existing /store/shuvam/.venv (cu130) → test in Phase 0
