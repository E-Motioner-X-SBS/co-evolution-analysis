# 01 — Datasets, Papers, and Files

Every dataset used in the 3D structural validation, its source paper, the
verified link, the files we hold, and their exact formats. Links were
HTTP-verified on 2026-08-13 (200 = live; journal DOIs resolve to the
publisher — 403 is publisher bot-blocking, valid in a browser).

---

## 1. The Omicron Spike dataset (primary analysis target)

| Property | Value |
|---|---|
| Sequences | **1,299** SARS-CoV-2 Omicron Spike proteins |
| Source | GISAID submissions (aligned FASTA) |
| File | `co-evolution-analysis/Spike_protein.aln-fasta` |
| Format | FASTA alignment, 1,276 columns, gaps `-` and ambiguous `X` (356 ×) present |
| Unique sequences | **299** (deduplicated by gap-stripped sequence; max group 587 identical) |
| Corrected variable positions | 21 (Shannon entropy H > 0.3, gaps excluded) |
| Corrected co-evolving pairs | 10 (mutation-only MI > 0.1, window 30) |
| Prime implicants | **36 distinct** (2 essential) — after correcting the 8-bit QM wrap-around defect (the pre-correction "152 rules" contained 143 phantom labels) |

**3D models of this dataset** — built with SWISS-MODEL (Automodel mode):
- **1,128 model PDB files** (`.pdb.gz`) + 1,128 ModelCIF (`.cif.gz`),
  299 unique sequences × 3–8 models each (median 3).
- Folder: `co-evolution-analysis/swissmodel/batch/pdbs/u####/`
  (`u0000`…`u0298`), where `u####` is the unique-sequence id.
- Per-sequence metadata in the same folder:
  - `full_details.json` — the full API payload: model_id, GMQE, QSQE,
    QMEAN, oligo_state (homo-trimer for most), **target_sequence** (the
    exact sequence the model was built from), pdb_release, smtl_version,
    built_with (ProMod3 3.6.0).
  - `templates_summary.json` — per-model: template SMTL id (e.g., 8d55.1),
    sequence identity, coverage, method, GMQE.
  - `model_XX.json` — per-model metadata.
- Model quality: best-model GMQE 0.69–0.72 for all 299 sequences (all > 0.6).
- Building infrastructure: `swissmodel_batch.py` (submit → poll → download,
  rate-limited), documented in `swissmodel/README.md`.

Paper (SWISS-MODEL server & pipeline):
- Waterhouse, A. et al. "SWISS-MODEL: homology modelling of protein
  structures and complexes." *Nucleic Acids Res.* 46(W1), W296–W303 (2018).
  https://doi.org/10.1093/nar/gky427 — **verified (resolves to OUP)**
- Bienert, S. et al. "The SWISS-MODEL Repository—new features and
  functionality." *Nucleic Acids Res.* 45, D313–D319 (2017).
  https://doi.org/10.1093/nar/gkw1132
- Studer, G. et al. "ProMod3—A versatile homology modelling toolbox."
  *PLoS Comput. Biol.* 17(1), e1008667 (2021).
  https://doi.org/10.1371/journal.pcbi.1008667 — **verified**

## 2. PSICOV 150 benchmark (structural ground truth #1)

| Property | Value |
|---|---|
| What | 150 protein domains with MSAs + **native experimental PDB structures** + PSICOV contact predictions |
| Paper | Jones, D.T., Buchan, D.W.A., Cozzetto, D., Pontil, M. "PSICOV: precise structural contact prediction using sparse inverse covariance estimation on large multiple sequence alignments." *Bioinformatics* 28(2), 184–190 (2012). https://doi.org/10.1093/bioinformatics/btr638 — **verified (resolves to OUP)** |
| Download | `https://bioinfadmin.cs.ucl.ac.uk/downloads/PSICOV/suppdata/suppdata.tar.gz` (54.6 MB) — **verified 200** |
| Local files | `data/psicov150/raw/{aln,seq,con,pdb}/` + `data/psicov150/converted/*.fasta` |
| MSA format | one gapped sequence per line, no headers; filename = PDB id (e.g., `1ctfA.aln`) |
| PDB format | Modeller-"repaired" native structures; **blank chain column** (single-chain); our parser handles it |
| Benchmark number | published mean top-L/5 long-range precision **0.44** (PSICOV); metaPSICOV 0.54 |

## 3. Pfam full alignments (PF00072/00595/00071)

| Property | Value |
|---|---|
| What | Deep, diverse family alignments — the canonical DCA test families |
| Paper (DCA method) | Morcos, F. et al. "Direct-coupling analysis of residue coevolution captures native contacts across many protein families." *PNAS* 108(49), E1293–E1301 (2011). https://doi.org/10.1073/pnas.1111471108 — **verified (resolves to PNAS)** |
| Paper (database) | Mistry, J. et al. "Pfam: The protein families database in 2021." *Nucleic Acids Res.* 49, D412–D419 (2021). https://doi.org/10.1093/nar/gkaa913 |
| Download | InterPro API: `https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00072?annotation=alignment:full` (same pattern PF00595, PF00071) — **verified 200** |
| Local files | `data/pfam/PF00072_full.sto` (486,732 seqs), `PF00595_full.sto` (330,330), `PF00071_full.sto` (221,818); converted FASTA in `data/pfam/converted/` |
| Format | Stockholm (gz-served, requests auto-decompresses) |
| Representative structures | RCSB: 6TNE (PF00072), 6AK2 (PF00595), 1HE1 (PF00071) — `data/structures/*.pdb`, fetched from `https://files.rcsb.org/download/<ID>.pdb` — **verified 200** |

## 4. EVmutation 20 alignments (Hopf 2017)

| Property | Value |
|---|---|
| What | 20 protein alignments used to benchmark mutation-effect prediction |
| Paper | Hopf, T.A. et al. "Mutation effects predicted from sequence co-variation." *Nat. Biotechnol.* 35(11), 1026–1031 (2017). https://doi.org/10.1038/nbt.3969 — **verified 200 (Nature)** |
| Download | `https://marks.hms.harvard.edu/evmutation/supp/alignments.tar.gz` (37.6 MB; **needs `curl -k`** — broken TLS chain) — **verified 200** |
| Local files | `data/evmutation/raw/*.a2m`, `data/evmutation/converted/*.fasta` |
| Format | A2M (aligned FASTA; lowercase = insert columns relative to the reference) |
| Depth | up to 28,048 sequences (BG_STRSQ); 20 proteins total |
| Note | per-protein structures are in the paper's supplement only — sequence-level analysis only; documented |

## 5. Unit test sets (plmc / GREMLIN / CCMpred)

| Dataset | Paper | URL (all **verified 200**) | Local file | Content |
|---|---|---|---|---|
| DHFR.a2m | Hopf, T.A. et al. "The EVcouplings Python framework for coevolutionary sequence analysis." *Bioinformatics* 35(9), 1582–1583 (2019). https://doi.org/10.1093/bioinformatics/bty657 | `https://raw.githubusercontent.com/debbiemarkslab/plmc/master/example/protein/DHFR.a2m` | `data/unit_tests/DHFR.a2m` | 3,629 seqs × 171 cols |
| 1whzA.msa | Kamisetty, H., Ovchinnikov, S., Baker, D. "Assessing the utility of coevolution-based residue–residue contact predictions in a sequence- and structure-rich era." *PNAS* 110(39), 15674–15679 (2013). https://doi.org/10.1073/pnas.1307869110 | `https://raw.githubusercontent.com/sokrypton/GREMLIN/master/1whzA.id90cov82.cut.msa` | `data/unit_tests/1whzA.msa` | 662 seqs × 65 cols |
| 1atzA.aln | Seemayer, S., Gruber, M., Söding, J. "CCMpred—fast and precise prediction of protein residue–residue contacts from correlated mutations." *Bioinformatics* 30(21), 3128–3129 (2014). https://doi.org/10.1093/bioinformatics/btu171 | `https://raw.githubusercontent.com/soedinglab/CCMpred/master/example/1atzA.aln` | `data/unit_tests/1atzA.aln` | 3,068 seqs × 75 cols |
| Structures | — | `https://files.rcsb.org/download/4p3r.pdb` etc. | `data/structures/{4p3r,1whz,1atz}.pdb` | DHFR / 1WHZ / 1ATZ |

## 6. DeepMSA (CASP12/13 free-modeling targets)

| Property | Value |
|---|---|
| What | CASP12/13 free-modeling targets with native structures (ground truth) |
| Paper | Zhang, C., Zheng, W., Mortuza, S.M., Li, Y., Zhang, Y. "DeepMSA: constructing deep multiple sequence alignment to improve contact prediction and fold-recognition for distant-homology proteins." *Bioinformatics* 36(7), 2105–2113 (2020). https://doi.org/10.1093/bioinformatics/btz978 — **verified (resolves to OUP)** |
| Download | `https://zhanggroup.org/DeepMSA/dataset.tar.bz2` (10.9 MB) — **verified 200** |
| Local files | `data/deepmsa/raw/data/<target>/{seq.fasta, native.pdb}` — 614 targets |
| Note | The tar ships **sequences + native structures only** (no MSAs — those come from the DeepMSA pipeline). Structural ground truth used; MSA-based analysis requires external MSA building (no hhblits/jackhmmer on this system) — documented |

## 7. GPCRdb (GPCR family alignments)

| Property | Value |
|---|---|
| What | Curated structure-based alignments of GPCR families (one sequence per receptor) |
| Paper | Kooistra, A.J. et al. "GPCRdb in 2021: integrating GPCR sequence, structure and function." *Nucleic Acids Res.* 49, D335–D343 (2021). https://doi.org/10.1093/nar/gkaa1080 |
| Download | `https://gpcrdb.org/alignment/fasta/001_001_001_001/` (5-HT1A, 125 seqs) and `.../001_001_003_008/` (β2AR, 136 seqs) — **verified 200** |
| Local files | `data/gpcrdb/{5ht1a,b2ar}.fasta`, `data/structures/{7e2x,2rh1}.pdb` |
| Note | 7E2X covers only the 219-residue TM core (fragment) — contact validation limited to that region (documented) |

---

## Summary table

| Dataset | Paper link (verified) | Structures | Role here |
|---|---|---|---|
| Omicron Spike (1,299) | GISAID; SWISS-MODEL docs (Waterhouse 2018) | **1,128 models** (299 unique) | main 3D test bed |
| PSICOV 150 | Jones 2012 (bioinformatics/btr638) | 150 natives in tar | **structural ground truth** |
| Pfam ×3 | Morcos 2011 (pnas.1111471108) | 6TNE/6AK2/1HE1 | deep families |
| EVmutation 20 | Hopf 2017 (nbt.3969) | — | deep alignments |
| Units (DHFR/1whzA/1atzA) | Hopf 2019 / Kamisetty 2013 / Seemayer 2014 | 4P3R/1WHZ/1ATZ | smoke tests |
| DeepMSA 614 | Zhang 2020 (btz978) | 614 natives | structural targets |
| GPCRdb | Kooistra 2021 (gkaa1080) | 7E2X/2RH1 | membrane family |
