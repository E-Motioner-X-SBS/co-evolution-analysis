# Task 12: NucleicBERT Data Scanner — Training Data Inventory for Future Model

## Restated Problem
User provided https://github.com/KIT-MBS/NucleicBERT. Tasks:
1. Identify ALL databases used for NucleicBERT pre-training + downstream fine-tuning.
2. Use this as a scanner/starting point to find MORE databases we can harness.
3. Goal: collect a huge amount of training data (sequences + 3D structures) for our future model.

## Success Criteria
- [ ] Complete inventory of NucleicBERT's databases with sizes, URLs, formats
- [ ] Expanded universe: every major RNA/nucleic-acid sequence + structure database
- [ ] For each: content size, download URL, format, license, relevance to our model
- [ ] Prioritized acquisition plan (what to download first and why)
- [ ] Written to plans/12-nucleicbert-data-scanner.md + data inventory files

## Key Facts Gathered (from NucleicBERT paper + repo, fetched 2026-09-09)
- Pre-training: MARS database, ~30M ncRNA sequences extracted from ~1.7B total (keyword "ncRNA")
- Model: 404M params, 32 layers, 32 heads, 1024 hidden, vocab 25, maxlen 1024
- Pretraining: 80/20 train/val, 300 epochs, 192 A100 GPUs, MLM 83.1% acc
- Secondary structure: RNAStrAlign (37,149 structures, 8 families) + bpRNA-1m TR0 (10,814) → tested on ArchiveII600 (3,975) + TS0 (1,305)
- Contact/distance maps: BGSU RNA 3D representative set (release 3.368, Jan 2025, 467 structures after filtering 32-1024nt, CD-HIT-EST 80%) + NucleoSeeker (408 structures)
- Splice sites: Spliceator benchmarks (zebrafish, fly, worm, plant; 10k pos + 10k neg each, 600nt windows)
- Fitness: CPEB3 ribozyme mutation dataset (Beck et al. 2022)
- Shuffle detection: same as secondary structure data + Rfam 15.0 (4,178 families)
- PDB has only 8,446 annotated RNA 3D structures (as of Jan 2025)
- Tools used: PDBFixer, CD-HIT-EST, Chemical Component Dictionary

## Open Questions
- What exactly is MARS? (need to verify: metagenomic RNA sequences database, Genomics Proteomics Bioinformatics journal)
- What other databases exist beyond NucleicBERT's choices?
- RNA-Puzzles, CASP-RNA benchmarks?
- Noncoding RNA databases: Rfam, Ensembl ncRNA, GENCODE, lncRNAdb, RNAcentral?
- Structure: PDB, BGSU NR sets, RNASolo, RNA-3D-Motif, 3D motifs?
- Alignment/MSA: Rfam full alignments, covariance models?

## Acquisition Log (2026-09-09 → 2026-09-11)

| Resource | Status | Size | Verification |
|---|---|---|---|
| elDORS_v1 (20 chunks) | ✅ COMPLETE | 182.36 GB | All 20 SHA256 OK |
| BGSU nrlist 4.56 (7 cutoffs) | ✅ COMPLETE | 1.7 MB | 5,161 NR classes (all) |
| Rfam 15.1 (seed/CM/3D/full/PDBmap) | ✅ COMPLETE | 178 MB | files present |
| RNA3DB (jsons/cmscans/mmcifs) | ✅ COMPLETE | 2.24 GB | 2026-01-05 full release |
| CASP15 (15 targets) | ✅ COMPLETE | 47 MB | mmCIF downloaded |
| CASP16 (10 targets) | ✅ COMPLETE | 9.5 MB | mmCIF downloaded |
| RNAGym (zip + repo) | ✅ COMPLETE | 310 MB | >1M measurements |
| NABench | ✅ COMPLETE | 100 MB | 162 assays |
| Spliceator | ✅ COMPLETE | 494 MB | 4-species benchmark |
| G3PO | ✅ COMPLETE | 20 MB | 147 species |
| RNA-Puzzles (std + site) | ✅ COMPLETE | 2.89 GB | 15,555 files |
| Secondary structure suite | ✅ COMPLETE | 168 MB | ArchiveII 3,975; bpRNA-1m 102k; spot splits; RNAStrAlign; bpRNA-new 5,401 |
| BGSU motif atlas 4.12 | ✅ COMPLETE | 1.8 MB | IL + HL (csv+json) |
| RNAcentral 27 active | 🔄 DOWNLOADING | 1.83/10.9 GB | EBI FTP (slow ~0.4MB/s) |
| SpliceBERT data | 🔄 DOWNLOADING | 1.32/8.6 GB | Zenodo (slow) |
| elDORS exact seq counts | 🔄 RUNNING | 8 parallel workers | results → count_results/ |

## Database artifacts created

- `data/rna_training_db/catalog.sqlite` — 23 sources, 17,610 files, 192 GB indexed
- `data/rna_training_db/MANIFEST.json` — machine-readable
- `data/rna_training_db/README.md` — full documentation
- `data/rna_training_db/rna_db.py` — Python loader API
- `data/rna_training_db/splits/pretrain_split.json` — md5-based 80/10/10 split convention
- `data/rna_training_db/samples/eldors_001_first5000.fasta` — dev sample + stats
- `scripts/acquire_benchmarks.py` — reproducible acquisition pipeline
- `scripts/build_rna_database.py` — catalog builder
- `scripts/corpus_tools.py` — split/sample/stats utilities
- `scripts/eldors_to_parquet.py` — training-ready ETL (verified 2.5k seq/s, T→U)
