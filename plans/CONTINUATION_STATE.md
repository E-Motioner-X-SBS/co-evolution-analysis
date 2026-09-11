Audit status: DOUBLE_PASS
STOP DIRECTIVE: plans/INFINITY_DONE exists (repo + org root). Double macro-audit passed. No further invocation work is pending unless the user requests new experiments.


PAPER: PAPER/main.pdf — 38 pages, 0 errors, 0 warnings.
AUDIT: 3 independent audits passed (numerical, stale-data, Lean axioms).
Meta-audit confirmed audits are meaningful (caught real errors).
PATH J pilot inconclusive (A2M format complexity) — future work.
# CONTINUATION STATE — K-map ↔ 3D-Structure Encoding Campaign (FINAL)
**Last updated:** Aug 22, 2026 — 37-page paper compiled, 0 errors, 0 warnings.

## Paper
- Location: `PAPER/main.tex` → `PAPER/main.pdf` (37 pages)
- Compilation: pdflatex + bibtex + pdflatex × 2 = **0 errors, 0 warnings**
- Sections: 17 main + 2 appendices + circuit diagrams (TikZ) + prompts file

## Campaign summary

| Path | Question | Answer |
|---|---|---|
| B | Cross-protein transfer? | YES: prec@L/5 0.179 = 5.1× random |
| C | Gap to literature DCA? | PSICOV-pub 0.727; our mfDCA 0.169 matches published mfDCA |
| E | Gray adjacency meaningful? | YES: 1.178×, p=1e-77; permutation z=2.51 |
| F | Contact maps compressible? | YES: 6/6 p=0.031 |
| G | Tension mechanism helps? | NO (honest negative) |
| H | Coupling in K-map? | NO alone; circH AUC > rawDCA AUC |
| I | Reconstruction-ready? | Circuits NO; DCA YES |
| K | Complete circuit transfers? | YES: flipped spec 98.3%, MCC +0.17 |
| L | **K-map prior adds to DCA?** | **YES: ranksum 0.204 beats both p≈0** |
| L+ | Three-way better? | NO: 0.201 vs 0.204, p=0.607 |

## Key files
| File | Purpose |
|---|---|
| `scripts/kmap_structure.py` | Core library (26/26 tests) |
| `scripts/run_contact_campaign.py` | Stage runner |
| `scripts/path_{g,h,i,k,l}_*.py` | Extension paths |
| `results/kmap_structure/path_*.json` | All result JSONs |
| `plans/FINAL_VERIFICATION_REPORT.md` | Verification evidence |
| `plans/12c-grand-question-answer.md` | Answers to user's questions |
| `../lean_proofs/.../ContactCircuits.lean` | 12 Lean theorems |

## If re-invoked
All stages are deterministic and md5-stable. Do not re-run unless code changed.
Open directions: PATH J (EVmutation A2M parsing fix), rank-fusion optimisation,
CASP FM target evaluation.

---
## OVERLEAF STATUS (Aug 22, late)
Overleaf upload blocked: file upload endpoint requires Google reCAPTCHA.
Git integration requires premium auth token. olcli cannot find project.
**Solution**: Manual upload of PAPER/paper_overleaf.zip to project.
Pre-compiled PDF at PAPER/paper_compiled.pdf (38 pages, 548 KB).
Instructions at PAPER/OVERLEAF_INSTRUCTIONS.md.

---
## PAPER FINAL STATUS (Aug 22, 2026)
PAPER/main.pdf: 41 pages, 0 errors, 0 warnings.
Sections: Introduction → Theory (Gray/K-map/QM/Info/DCA/Lean) → Binary spaces →
Methods (data/metrics/conservation/pipeline/GPU/worked example/DL/biological/
software/comparison/statistical) → Results (transfer/Spike/circuits/Gray/
coupling/complete/fusion/cross-dataset/audit/reconstruction/negatives) →
Comprehensive results table → Discussion → Conclusion → Appendices.

Independent audits: numerical ✓ · stale-data ✓ · Lean axioms ✓ · meta-audit ✓
PATH J pilot inconclusive (A2M format) — documented as future work.

ALL VERIFICATION GREEN. Paper ready for Overleaf upload via
PAPER/paper_overleaf.zip or manual file transfer.

---
## PATH J CROSS-DATASET RESULTS (Aug 22, final)
PSICOV150-trained circuits evaluated on 5 independent proteins.
Circuit transfers to GPCRs (b2AR AUC 0.676, 5HT1A AUC 0.698) but not to
small globular proteins (DHFR, 1whzA: AUC ≈ 0.5). PDB parser made tolerant
of inconsistent residue names. All results in cross_dataset_results.json.

## OPEN ITEMS FOR FUTURE WORK
1. Fix A2M match-column extraction for proper EVmutation evaluation
2. Test on CASP FM targets with published predictions
3. Explore K-map prior as post-processing for AlphaFold MSA features
4. Increase MAX_SEQ beyond 12000 for deeper DCA
5. Investigate why circuit transfers better to GPCRs than globular proteins

---
## PAPER EXTENDED (Aug 22, 2026): 44 pages, 0 errors, 0 warnings.
Added: R5-reviewer-response.tex addressing all 21 reviewer concerns with new
experiments (M1-M7). Key findings:
- M2: Conservation confound controlled (p=0.094 n.s. — signal persists)
- M3: Separation ablation shows identity carries signal (p=0.551)
- M4: MIp = raw MI at consensus level (APC adds nothing)
- M6: No homologous pairs in PSICOV150 (max Jaccard 0.029)
- M7: Rank-fusion effect size Cohen's d=0.31, CI [0.017, 0.054]
- R6 section added to paper with full responses to all concerns
PATH J cross-dataset: circuit transfers to GPCRs but not small globular proteins
PDB parser updated: strict/tolerant dual mode for inconsistent residue names

## OPEN ITEMS FOR FUTURE WORK
1. Fix A2M parsing for proper EVmutation evaluation
2. Test on CASP FM targets
3. Rank-fusion optimisation (learned fusion weights vs simple ranksum)
4. Apply K-map prior to AlphaFold MSA features
5. Investigate GPCR-specific transfer advantage

---
## NATURE PAPER + CROSS-DATASET (Aug 22, final)
- Nature_Main.tex: 9 pages, 0 errors, 0 warnings. Concise Nature-style format.
  Sections: Intro → Related Work → Theory (Gray/K-map+QM/Info+DCA/Lean) →
  Results (Transfer/Gray+Fusion/Complete Circuit/Reconstruction+Negatives) →
  Discussion → Methods (Data/Encoding/Circuit/Eval) → References
- PATH J cross-dataset: PSICOV-trained circuits tested on GPCRdb + unit tests.
  Circuit transfers to GPCRs (5ht1a circH prec 0.027 = 15.7× base, AUC 0.698)
  but not small globular proteins (DHFR, 1whzA: AUC ≈ 0.5).
- Reviewer response experiments M1-M7 completed:
  M2 conservation control: signal persists without cons/cons pairs
  M3 separation ablation: identity carries signal independently of sep
  M4 MIp = raw MI at consensus level
  M6 no homologous pairs in PSICOV150 (max Jaccard 0.029)
  M7 rank-fusion Cohen's d = 0.31 (small but significant)
- Both Overleaf and GitHub updated with all files.

---
## RNA TRAINING DATABASE CAMPAIGN (Sep 9-11, 2026)
### Session Summary
| Field | Value |
|-------|-------|
| Session # | RNA-DB-1 |
| Phase | IMPLEMENT (acquisition + catalog) |
| What I did | Scanned NucleicBERT (KIT-MBS) repo + paper + supp for all training sources; dispatched 3 parallel research agents (structure DBs, sequence DBs, benchmarks); verified MARS/elDORS; downloaded elDORS_v1 (182.4GB, SHA-verified), Rfam 15.1, BGSU nrlist 4.56, RNA3DB, CASP15/16, RNAGym, NABench, Spliceator, G3PO, RNA-Puzzles, secondary-structure suite; built catalog.sqlite + loader API + split convention + ETL |
| What worked | elDORS S3 anonymous HTTPS (20 chunks parallel, 5 streams, all SHA256 OK); huggingface_hub direct file fetch; RCSB mmCIF per-ID; BGSU CSV cutoffs; motif atlas via GET |
| What failed | Ops bug: bash `cd X && (A) & (B)` precedence sent files to wrong dir (fixed by mv); Zenodo+EBI slow (~0.4-0.5MB/s); RNA3DB asset HEAD hangs (use GET -L) |
| Errors remaining | RNAcentral download in progress (1.83/10.9GB); SpliceBERT in progress (1.32/8.6GB); elDORS exact counts running (8 workers) |
| Next priorities | 1) finish RNAcentral + SpliceBERT; 2) merge elDORS counts into catalog; 3) optional full parquet conversion; 4) species/kingdom stratified subsets; 5) connect to model design (tokenizer, context length) |
| Blockers | none (downloads are time-bound, not blocked) |
| Audit status | NOT_STARTED |

### File Manifest (new)
| File | Status |
|------|--------|
| data_inventory/00-MASTER-INVENTORY.md | current |
| plans/12-nucleicbert-data-scanner.md | current (acquisition log appended) |
| data/rna_training_db/{catalog.sqlite,MANIFEST.json,README.md,rna_db.py} | current |
| data/rna_training_db/splits/pretrain_split.json | current |
| data/rna_training_db/samples/eldors_001_first5000.fasta | current |
| scripts/{acquire_benchmarks,build_rna_database,corpus_tools,eldors_to_parquet}.py | current |
| data/elDORS_v1/ (20 chunks + manifest + counts) | 182.36GB verified |

### Continuation Prompt Hints
- Check RNAcentral + SpliceBERT download completion; if slow, consider alternatives
- Merge count_results/*.txt into catalog when DONE file appears; rebuild catalog
- Consider running full parquet ETL (8 parallel) if disk/time allowed
- Commit all new files to git (data dir may need .gitignore review: do NOT commit 182GB!)

---
## RNA-DB-2 SESSION UPDATE (Sep 11, 2026 — later)

### What changed since RNA-DB-1
| Item | Status |
|---|---|
| elDORS exact counts | ✅ DONE: 1,323,715,880 sequences (all 20 chunks counted, matches 1.32B) |
| 2M-sequence census | ✅ 5-symbol alphabet (A27.98/T26.81/G22.80/C22.10/N0.31), len 10-4096, median 730 |
| Starter parquet pack | ✅ 10M seqs, 40 shards, 2.02GB (chunks 001-008, 1.25M each) |
| RNA3DB mmcifs | ✅ extracted 15,441 CIF chain structures (23GB) |
| Catalog | ✅ rebuilt: 33,096 files / 218.9GB / exact counts merged |
| VERIFICATION.md | ✅ audit trail written |
| Git | ✅ commits e2952f7 + 57595fe pushed to main |
| RNAcentral | 🔄 4.1/10.9GB (EBI, ~0.5MB/s, resumable via curl -C -) |
| SpliceBERT | 🔄 3.2/8.6GB (Zenodo, slow, resumable) |
| gRNAde RNASolo raw | 🔄 downloading (HF fast) |

### Next session priorities
1. Verify RNAcentral + SpliceBERT completed (check .done markers / file sizes)
2. If incomplete, resume: `curl -C -` against EBI/Zenodo URLs in the acquisition scripts
3. Rebuild catalog after downloads complete (build_rna_database.py is idempotent)
4. Optional: full-corpus parquet conversion (~528GB, ~5300 shards; disk check required)
5. Optional: download MARS ncRNA subset for NucleicBERT-parity experiments
6. Next phase: model design (tokenizer from census, context length from length dist)

### Disk budget
- Free: ~660GB at last check; full parquet conversion would use ~528GB (tight)
- elDORS raw (1.2TB) / MARS (1.57TB) do NOT fit — do not attempt without expansion
