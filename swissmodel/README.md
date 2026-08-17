# SwissModel — Template Matching for the Omicron Spike Dataset

SWISS-MODEL (https://swissmodel.expasy.org) integration for the 1,299-sequence
Omicron Spike co-evolution dataset. Answers: *which experimental/AFDB structures
match each sequence, and how well.*

## Contents

| File | Purpose |
|------|---------|
| `swissmodel_template_search.py` | Reusable tool: template search for any sequence (guest or official API mode) |
| `swissmodel_seq1_templates.json` | 50 templates for sequence #1 (WRU87367.1) — run 1 (guest, project 8XmGU6) |
| `swissmodel_seq1_templates_e2e.json` | Same sequence, run 2 (guest, project ZHyAMb) — identical results |
| `.env` | **SECRET — gitignored.** SWISS-MODEL API token (do not commit) |
| `plans/` | Research archive, API contract, decisions, audit log for this work |

## Two ways to get template matches

### 1. Guest mode (no account) — works out of the box
```bash
python swissmodel/swissmodel_template_search.py \
    --fasta Spike_protein.aln-fasta --seq-index 0 --mode guest --out results/seq0_templates.json
```
Submits to the web form (`POST /interactive`), polls the project, parses the
embedded `globalTemplateData` JSON. Runtime for a 1,269-aa Spike: **~7–15 min**.

### 2. Official Modelling API (token) — full automodel incl. 3D model building
```bash
python swissmodel/swissmodel_template_search.py \
    --fasta Spike_protein.aln-fasta --seq-index 0 --mode coreapi --out results/seq0_full.json
```
Token is read automatically from `swissmodel/.env` (or `SWISSMODEL_API_TOKEN` env).
This runs the real pipeline (template search + ProMod3 modelling + QMEAN), so
allow **30–60+ min** per sequence. API contract: `POST /automodel/` → 202,
poll `GET /project/{id}/models/full-details/` until COMPLETED/FAILED.

## Verified results (sequence #1: WRU87367.1, 1,269 aa)

- **50 templates** matched (41 HHblits + 9 BLAST; all cryo-EM)
- Top: `8d55.1` *Closed state of SARS-CoV-2 BA.2 variant spike* — **98.3% identity, 99.9% coverage**
- Identity range 76–99% (median 94%); coverage 93–100%
- Two independent runs → identical template sets (deterministic)
- Coreapi job `2870bf`: submitted 202, status RUNNING (model building)

## Template fields (per hit)

`rank, qname (SMTL id), pdb_id, chain, method, seq_identity (seq_id), coverage,
seq_similarity (seq_sim), found_by, resolution_A, oligo_state, qsqe,
gmqe_afdb_plddt (pred_lddt), title` — plus full alignments `trg_seq`/`tpl_seq`
in the raw JSON.

## Security

`.env` is gitignored (verified with `git check-ignore`). The token was shared in
plaintext chat once — consider regenerating it at https://swissmodel.expasy.org/account
after heavy use.

## Batch run (Aug 8, 2026) — all 1,299 sequences modelled ✅

`swissmodel_batch.py` shipped the **299 unique sequences** (dedup of the 1,299 —
the server caches by sequence, verified: resubmit → 200 + same project_id),
polled all projects, and downloaded every model. **All 299 COMPLETED, 0 FAILED.**

| Metric | Value |
|--------|-------|
| Original sequences covered | **1,299 / 1,299** |
| Unique sequences (jobs) | 299 |
| Model PDB + CIF files | 1,128 + 1,128 (850 MB) |
| Models per sequence | 3–8 (median 3) |
| Best-model GMQE | 0.69–0.72 (**all > 0.6**) |
| Top templates | 8hxk.1.B (218), 8d55.1.C (204), 8d56.1.C (204), 6zgf.1.C (93), 7cn8.1.C (81) |
| Runtime | ~75 min total (submit ~13 min; ~6 jobs/min peak) |
| Rate limit incidents | 0 × 429 (paced 30 sub/min, ≤60 poll/min) |
| Integrity | all 1,128 PDBs pass gzip; 60/60 sample have ATOM (median 26,814); target-sequence exact match on 5/5 spot check |

Structure:
```
batch/
├── manifest.json          ← every original index (0..1298) → uid, status, PDB path
├── results_index.json     ← per-uid summary: best model, GMQE, templates used
├── state/                 ← checkpoints (submitted/done/failed/uniques/orig_map)
├── pdbs/u####/            ← per unique sequence: *.pdb.gz, *.cif.gz, model_*.json,
│                             full_details.json, templates_summary.json
└── orchestrator.log       ← full run log
```
Re-run / resume: `python swissmodel_batch.py run` (idempotent; `verify` regenerates the manifest).

## Per-entry view (1,299 folders)
The 1,299 FASTA entries contain 299 unique sequences (identical sequences share a model —
the server itself caches by sequence). For pipelines that want one folder per FASTA entry:
`materialize_entries.py` creates `batch/pdbs_by_entry/entry_####_<accession>/` with hardlinks
to that entry's model files (+ `templates_summary.json`, `entry.json`). Hardlinks = 4.0 GB
worth of "copies" in only 828 MB of real storage.
