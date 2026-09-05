# 01 — Research: Validation Dataset Sources (all URLs verified live, Aug 8 2026)

> **CORRECTION (Aug 22, 2026):** the 'PSICOV published L/5 = 0.44' quoted below is a misattribution; the paper's PSICOV L/5 ≈ 0.73 (we reproduce 0.727). See 11_Kmap_Structure_Encoding.md §9-S1.

Two independent research agents + my own spot-checks (curl HEAD/partial + full fetch of small files).

| # | Dataset | Paper | Verified URL(s) | Format | Depth | Length | Structures | Benchmark ground truth |
|---|---------|-------|-----------------|--------|-------|--------|-----------|------------------------|
| 1 | PSICOV 150 | Jones et al. 2012, Bioinformatics 28:184 (also metaPSICOV 2015, CCMpred 2014 set) | `https://bioinfadmin.cs.ucl.ac.uk/downloads/PSICOV/suppdata/suppdata.tar.gz` (302→**200**, 54.6 MB) | tar.gz: `aln/` (150 MSAs, 1 gapped seq/line, no headers), `seq/` targets, `con/` PSICOV predictions, `pdb/` native structures | 511–74,836 (median 3,228; 140/150 ≥1,000) | 56–235 | **native PDBs in the tar** | **published top-L/5 precision (PSICOV 0.44 mean; metaPSICOV 0.54)** |
| 2 | Pfam full | Morcos et al. 2011 PNAS (DCA canon); InterPro (Pfam 38.2) | `https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00072?annotation=alignment:full` (**200**, gz; ~90 MB) — same pattern PF00595, PF00071 | Stockholm (gz) | 486,732 / 330,330 / 221,818 | ~1150 cols (domain ~120) | rep. PDBs 6TNE (PF00072), 6AK2 (PF00595), 1HE1 (PF00071) — to fetch from RCSB | DCA results on these families in Morcos 2011 + many papers |
| 3 | EVmutation 20 | Hopf et al. 2017, Nat Biotechnol 35:1026 | `https://marks.hms.harvard.edu/evmutation/supp/alignments.tar.gz` (**200**, 37,596,525 B; **needs curl -k**: broken TLS chain) | A2M (aligned FASTA, lowercase = inserts) | up to 28,048 (BRCA1) | 18–716 | per-protein (e.g., 5GJW BRCA1) | published mutation-effect benchmarks |
| 4 | Unit tests | plmc (Hopf 2019); GREMLIN (Kamisetty 2013); CCMpred (Seemayer 2014) | `raw.githubusercontent.com/debbiemarkslab/plmc/master/example/protein/DHFR.a2m` (**200**); `.../sokrypton/GREMLIN/master/1whzA.id90cov82.cut.msa` (**200**); `.../soedinglab/CCMpred/master/example/1atzA.aln` (**200**) | A2M / FASTA | 3,629 / 662 (id90-clustered) / 1,534 | 80 / 78 / 75 | 4P3R, 1WHZ, 1ATZ | tool example sets |
| 5 | DeepMSA | Zhang et al. 2020, Bioinformatics 36:2105 | `https://zhanggroup.org/DeepMSA/dataset.tar.bz2` (**200**, 10.9 MB) | tar.bz2 (CASP12/13 FM + membrane targets; a3m-style) | deep | domain-level | CASP targets (native structures at CASP/PDB) | contact/TM results in paper |
| 6 | GPCRdb | Kooistra et al. 2021, NAR 49:D335 | `https://gpcrdb.org/alignment/fasta/001_001_001_001/` (5-HT1A, **200**) and `.../001_001_003_008/` (β2AR, **200**) | FASTA (curated, gapped) | 125 / 136 | ~968 gapped | 7E2X, 2RH1 | GPCR contact preds (Hopf 2012) |

## Verified dead (agents checked — do NOT use)
`bioinf.cs.ucl.ac.uk/downloads/psicov` (404 — use bioinfadmin), `gremlin.bakerlab.org`, `evcomplex.org`, `pfam.xfam.org`, old EVcouplings precomputed API, `marks.hms.harvard.edu` without `-k`, PSICOV "original" path.

## Key format facts (for Phase 3 converters)
- InterPro returns **Stockholm gzipped**; `format=fasta` → 404. Must parse Stockholm (#=GS header lines).
- PSICOV `aln/` = one gapped sequence per line, NO headers; each filename = PDB id (e.g., 1ctfA.aln).
- A2M: lowercase = insert relative to reference; decision: **uppercase-all** (treat as residues) — document.
- DeepMSA: a3m-style inside tar.bz2 (inspect first, then convert).
- GPCRdb: already FASTA.

## Redundancy situation per dataset (for honest reporting)
- PSICOV: duplicate rows removed by construction; HHblits/UniProt20 builds → diverse.
- Pfam full: REDUNDANT by design (all homologs) → we will ALSO produce an id90-clustered copy for analysis (standard EVcouplings protocol) and report both.
- EVmutation: jackhmmer builds, unclustered → report redundancy stats, cluster copy.
- Unit tests: 1whzA already id90; others small (redundancy stats reported).
- DeepMSA: deep multi-source (UniClust30+UniRef90+Metaclust) → diverse.
- GPCRdb: curated one-seq-per-receptor → low redundancy.
