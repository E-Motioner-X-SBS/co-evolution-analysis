# RNA/Nucleic-Acid Training Data Master Inventory
**For: future RNA foundation model development (NucleicBERT-inspired + beyond)**
**Compiled: 2026-09-09 | All URLs verified live this session unless tagged [UNVERIFIED]**

## Part 1: What NucleicBERT Used (the scanner seed)

| Purpose | Database | Size used | Source URL |
|---|---|---|---|
| **Pre-training (MLM)** | MARS (ncRNA-tagged subset) | ~30M sequences (from 1.7B total) | https://ngdc.cncb.ac.cn/omix/release/OMIX003037 |
| Secondary structure (train) | RNAStrAlign | 37,149 structures | rna.urmc.rochester.edu / HF: rouskinlab/RNAstralign |
| Secondary structure (train) | bpRNA-1m TR0 | 10,814 structures | sparks-lab.org/server/spot-rna/ |
| Secondary structure (test) | ArchiveII600 | 3,975 (≤600nt cutoff) | rna.urmc.rochester.edu |
| Secondary structure (test) | bpRNA TS0 | 1,305 structures | sparks-lab.org |
| Contact/distance maps (train) | BGSU NR 3.368 (Jan 2025) | 467 structures (32-1024nt, CD-HIT-EST 80%) | rna.bgsu.edu/rna3dhub/nrlist |
| Contact/distance maps (train) | NucleoSeeker set | 408 structures | github.com/theuutkarsh/nucleoseeker |
| Splice sites | Spliceator benchmarks | 4 species × (10k pos + 10k neg), 600nt | bigest-icube.fr/spliceator |
| Fitness | CPEB3 ribozyme DMS | ~10^5 variants (Beck 2022) | Front. Mol. Biosci. 9:893864 |
| Shuffle detection | Rfam 15.0 | 4,178 families | ftp.ebi.ac.uk/pub/databases/Rfam/ |

NucleicBERT model: 404M params (32L × 32H × 1024d), vocab 25, maxlen 1024. Weights: https://zenodo.org/records/16989562 (18GB, 5 checkpoints).

## Part 2: The Expanded Universe — SEQUENCE DATABASES

### Tier 1: Pre-training scale (billions)

| Database | Size (verified) | Download | Notes |
|---|---|---|---|
| **elDORS_v1 (clustered)** | **1.32B seqs, 170GB** (80% ID, 10-4096nt) | AWS S3: registry.opendata.aws/eldors_v1 → s3://eldors_v1/ | **NEWEST (Jul 2026), MARS successor, explicitly built for RNA LM pre-training. CC BY 4.0. THIS IS THE TOP ACQUISITION TARGET.** |
| elDORS_raw_v1 | 1.91B seqs, 1.2TB | same registry (elDORS_v1_raw/) | Raw tier; disk-constrained for us |
| **MARS 1.0** | 1.73B seqs, 1.57TB (30 tgz parts) | https://download.cncb.ac.cn/OMIX/OMIX003037/ | What NucleicBERT used; elDORS supersedes it |
| MARS composition | RNAcentral + MG-RAST + GWH + MGnify + NCBI nt/env_nt/tsa_nt/pat_nt, deduped | — | 20× NCBI nt, 60× RNAcentral |
| NCBI nt + subsets | nt ~96M seqs; core_nt <½ of nt | ftp.ncbi.nlm.nih.gov/blast/db/ | Extract FASTA via blastdbcmd |
| MGnify | 5.7B+ proteins; RNA via mgnify_genomes | ftp.ebi.ac.uk/pub/databases/metagenomics/ | Metagenomic RNA source |
| MG-RAST | ~260k+ public metagenomes | mg-rast.org | MARS component |
| Genome Warehouse (GWH) | >1PB assemblies | ftp://download.big.ac.cn/gwh | MARS component |

### Tier 2: Curated ncRNA (millions)

| Database | Size | Download | Notes |
|---|---|---|---|
| **RNAcentral** | **Rel 27 (Jul 2026): >55M seqs**, 52 expert DBs | ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/ | THE ncRNA hub; used by RiNALMo (36M) |
| Rfam | **15.1 (Jan 2026): 4,227 families**, 26,106 genomes in Rfamseq | ftp.ebi.ac.uk/pub/databases/Rfam/ | Seed+full alignments, CMs. CC0 |
| Ensembl ncRNA | Rel 116 (Jun 2026), ~8,185 genomes | ftp.ensembl.org/pub/release-116/fasta/*/ncrna/ | Used by RNA-FM (23.7M) |
| RefSeq RNA | Rel 237 (Aug 2026): 86M RNAs | ftp.ncbi.nlm.nih.gov/refseq/release/ | Curated reference |
| GENCODE v50 | 35,885 human lncRNA + 7,608 small-ncRNA genes | ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/ | Human annotation gold standard |

### Tier 3: Specialized ncRNA families

| Database | Size | Download |
|---|---|---|
| miRBase v23 (Aug 2026) | 38,589 hairpins / 48,860 mature / 271 species | mirbase.org/download/ |
| tRNADB-CE | 19.6M tRNA genes | trna.ie.niigata-u.ac.jp |
| GtRNAdb rel 22 | 431k tRNAs from 4,867 genomes | gtrnadb.ucsc.edu |
| SILVA 138.2 | SSU 9.4M+ / SSU NR99 510k rRNAs | arb-silva.de/download/ |
| GTDB R10-RS226 | 732,475 MAG genomes (ncRNA mining) | data.gtdb.ecogenomic.org |
| IMG/VR v4 | 15M viral genomes incl. RNA viruses | img.jgi.doe.gov |

## Part 3: The Expanded Universe — 3D STRUCTURE DATABASES

| Database | Size (verified 2026-09) | Download | Notes |
|---|---|---|---|
| **PDB/RCSB (RNA-containing)** | **22,185 entries** w/ nucleic acid polymers | files.rcsb.org | Weekly; CC0. 259,747 total structures |
| **NAKB** (NDB successor) | 22,344 nucleic-acid 3D structures | nakb.org | Weekly annotations, RNAEQ sets |
| **BGSU nrlist (current 4.56)** | **5,161 NR equivalence classes** ("all" cutoff) | rna.bgsu.edu/rna3dhub/nrlist/download/rna/4.56/all/csv | Weekly since 2011; resolution cutoffs 1.5Å-20Å; CC BY 4.0. THE standard NR training source |
| **RNAsolo2** | **22,412 cleaned RNA structures / 5,107 classes / 2,688 precompiled benchmark sets** | rnasolo.cs.put.poznan.pl | Cleaned PDB-derived; torsion angles; ML-ready |
| RNA 3D Motif Atlas (4.12) | 2,972 IL loops (413 motifs) + 2,020 HL loops (254 motifs) | rna.bgsu.edu/rna3dhub/motifs | Recurrent motifs for motif-level supervision |
| NucleoSeeker | Filtering pipeline (pip install nucleoseeker) | github.com/theuutkarsh/nucleoseeker | PDB curation tool, not a DB |
| RNA-Puzzles | ≥79 puzzles (PZ1-PZ79 through Apr 2026) | rnapuzzles.org/puzzles/list/ | Blind benchmark; Round V = 23 targets |
| CASP15 RNA (2022) | 13 RNA targets (R1107-R1190) | predictioncenter.org/casp15 | Verified PDBs: 7qr4, 7qr3, 8s95, 8fza, 8tvz, 8btz, 7zj4, 7ptk, 7ptl, 8uys, 8uye, 8uyg, 8uyj, 7yr7, 7yr6 |
| CASP16 RNA (2024) | 42 NA targets, 65 groups | predictioncenter.org/casp16 | Verified PDBs: 8uo6, 9cfn, 9c2k, 9dcf, 9b0l, 9ely, 9bzc, 9bz1, 9cbu, 9cbx. No TM>0.8 on unseen natural RNAs |
| RNA3DB | All PDB RNA chains, Rfam-labeled, NR train/test splits | github.com/marcellszi/rna3db | DL-ready splits, periodically updated |
| Rfam 3D mapping | 65+ families w/ 3D seeds | ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.3d.seed.gz | Family↔structure join table: rna.bgsu.edu/data/pdb_chain_to_best_rfam.txt |
| NuFold DB | Computed RNA 3D models | nufold.kiharalab.org | Augmentation only, not ground truth |

**Key scarcity fact**: only ~22k experimental RNA 3D structures exist (vs ~1M+ protein). This is THE bottleneck for RNA 3D supervision.

## Part 4: BENCHMARK / EVALUATION DATASETS

### Secondary structure (2D)
| Dataset | Size | Download |
|---|---|---|
| ArchiveII | 3,975 original / 2,975 filtered | HF: multimolecule/archiveii |
| bpRNA-1m | 102,318 structures | bprna.cgrb.oregonstate.edu |
| bpRNA spot splits | TR0=10,814 / VL0=1,300 / TS0=1,305 | sparks-lab.org; HF: multimolecule/bprna-spot |
| RNAStrAlign | 37,149 raw / 27,125 valid after dedup | HF: rouskinlab/RNAstralign |
| bpRNA-new (successor; NO "bpRNA-2m" exists) | Rfam 14.2 distinct families | HF: multimolecule/bprna-new |

### Tertiary / contacts
| Dataset | Size | Download |
|---|---|---|
| RNA-Puzzles standardized | PZ1-PZ21 core + Round V | github.com/mmagnus/RNA-Puzzles-Standardized-Submissions |
| CASP15/16 RNA | 13 + 42 targets | predictioncenter.org |
| BGSU NR (any pinned release) | per-release CSV | rna.bgsu.edu/rna3dhub/nrlist |
| RNA-Benchmark (Poznań) | monomer/RNA-RNA/RNA-protein | github.com/iammarcol/RNA-Benchmark |

### Fitness / function
| Dataset | Size | Download |
|---|---|---|
| **RNAGym (2025)** | **>30 DMS assays, >1M measurements** (tRNA, aptamers, ribozymes, splicing mRNAs) | github.com/MarksLab-DasLab/RNAGym; rnagym.org; HF: Marks-lab/RNAgym |
| **NABench (Nov 2025)** | 162 assays, 2.6M mutated sequences, 29 foundation models | github.com/mrzzmrzz/NABench; arXiv:2511.02888 |
| CPEB3 DMS (Beck 2022) | ~10^5 variants | In RNAGym bundle |
| RNAcmap3 benchmark | 30 Rfam + 105 non-Rfam RNAs | via GPB paper |

### Splicing
| Dataset | Size | Download |
|---|---|---|
| Spliceator | 4-species test (zebrafish/fly/worm/plant), 10k+10k each | bigest-icube.fr/spliceator/static/data/ |
| SpliceBERT data | 2M+ pre-mRNAs from 72 vertebrates | zenodo.org/records/7995778 |
| G3PO+ | 147 eukaryotes, >20k genes | github.com/BiGEst-ICube/g3po |

### MSA / coevolution
| Dataset | Size | Download |
|---|---|---|
| Rfam full alignments | 4,227 families | ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/full_alignments/ |
| rMSA pipeline | nt + RNAcentral merged | github.com/pylelab/rMSA; zhanggroup.org/rMSA/ |
| RNAcmap3-elDORS | pre-indexed BLAST/infernal builds (861GB) | S3: elDORS registry (RNAcmap3_optimized_elDORS/) |

## Part 5: Cross-Reference — What Other RNA LMs Used

| Model | Params | Corpus | Key source |
|---|---|---|---|
| NucleicBERT | 404M | MARS ncRNA ~30M | OMIX003037 |
| RiNALMo | 650M (largest) | RNAcentral 36M | github.com/lbcb-sci/RiNALMo |
| RNA-FM | 100M | Ensembl ncRNA 23.7M | github.com/ml4bio/RNA-FM |
| RNAErnie | 86M | RNAcentral ~23M (motif-aware) | Nat. Mach. Intell. 2024 |
| SpliceBERT | 19M | 72-vertebrate pre-mRNA 2M+ | zenodo 7995778 |
| RNA-MSM | MSA-based | RNAcmap3 MSAs vs MARS | github.com/yikunpku/RNA-MSM |
| 3UTRBERT | BERT-base | human 3'UTRs | figshare 26082916 |
| UTR-LM | ESM2-style | 255,795 5'UTRs (5 species) | github.com/a96123155/UTR-LM |
| Nucleotide Transformer | 50M-2.5B | 3,202 genomes (DNA) | instadeepai |
| ERNIE-RNA | BERT-base | ncRNA, structure-aware MLM | Nat. Comm. 2025 |

## Part 6: Prioritized Acquisition Plan

**Disk constraint: 896GB free on /store. Full MARS (1.57TB) and elDORS_raw (1.2TB) DO NOT FIT. elDORS_v1 clustered (170GB) FITS.**

### Priority 1 (immediate, fits disk, highest value)
1. **elDORS_v1 clustered** (170GB, 1.32B seqs) — the single best pre-training corpus available; MARS successor; CC BY 4.0. S3: `aws s3 sync s3://eldors_v1/ ...` (check exact bucket path from registry page)
2. **RNAcentral active** (~10-20GB compressed, 55M seqs) — curated ncRNA tier
3. **BGSU nrlist 4.56 CSVs** (MBs) — NR structure classes for 3D training
4. **RNAsolo2 archive** (GBs) — 22,412 cleaned structures + 2,688 benchmark sets
5. **Rfam 15.1** (GBs) — seed/full alignments + CMs + 3D seeds

### Priority 2 (benchmarks + eval, all small)
6. ArchiveII, bpRNA spot splits, RNAStrAlign, bpRNA-new (HF multimolecule/rouskinlab)
7. RNA-Puzzles standardized submissions + CASP15/16 targets (PDBs)
8. RNAGym + NABench (fitness benchmarks)
9. Spliceator + SpliceBERT data + G3PO+ (splicing)
10. RNA3DB splits, RNA 3D Motif Atlas 4.12

### Priority 3 (scale-up, needs disk expansion or selective download)
11. MARS ncRNA-tagged subset only (~30M seqs, the NucleicBERT slice — re-extract from MARS tgz parts rather than downloading all 1.57TB)
12. Ensembl ncRNA (all divisions), RefSeq RNA, GENCODE
13. MGnify/GWH/MG-RAST metagenomic RNA (huge; sample selectively)
14. elDORS_raw / full MARS (only after disk expansion)

### Priority 4 (specialized families)
15. miRBase, tRNADB-CE, GtRNAdb, SILVA, GTDB MAGs, IMG/VR

## Open Questions / Next Steps
- [ ] Verify exact S3 bucket path for elDORS_v1 download (registry page → bucket ARN)
- [ ] Decide pre-training corpus composition (elDORS_v1 alone vs + RNAcentral + Ensembl)
- [ ] Decide 3D supervision strategy given ~22k structure scarcity (contact maps from BGSU NR like NucleicBERT? RNAsolo2 multimodal?)
- [ ] Plan MSA branch (RNAcmap3-elDORS pre-indexed builds, 861GB — needs disk)
- [ ] Model architecture decisions come AFTER data is secured (per user: data first)
