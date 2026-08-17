# 05 — Audit Log: SwissModel CoreAPI learning + single-sequence template test

## Micro-audits
| Step | Check | Result |
|------|-------|--------|
| API docs | coreapi page + schema.js + OpenAPI + help + wrapper all fetched live | ✅ agree |
| Auth requirement | anonymous POST /automodel/ | ✅ 401 verified |
| Sequence prep | seq1 WRU87367.1 gap-stripped 1269 aa (7 gaps) | ✅ |
| Guest submit | POST /interactive, CSRF flow | ✅ project 8XmGU6 |
| Poll | automodel_poll JSON + page status | ✅ |
| Results parse | globalTemplateData, 50 fields/template | ✅ |
| Biolog. sanity | top template = BA.2 Spike closed state 8d55.1 (98.3% id) | ✅ |
| Reproducibility | 2nd run ZHyAMb: identical (pdb,chain,rank) set | ✅ |
| Script | read_fasta + parse unit-tested; E2E guest run completed | ✅ |

## Files produced
- `co-evolution-analysis/swissmodel_template_search.py` — reusable tool (guest + coreapi modes)
- `co-evolution-analysis/swissmodel_seq1_templates.json` — 50 templates (run 1, 8XmGU6)
- `co-evolution-analysis/swissmodel_seq1_templates_e2e.json` — 50 templates (run 2, ZHyAMb)
- `/tmp/opencode/seq1.fasta`, `crambin_final.html`, `templates_final.html`, `e2e_final.html`, `spike_seq1_templates.json` — raw artifacts
- `plans/00..05` — this workspace record

## Remaining (flagged)
- coreapi automodel E2E: BLOCKED on token (O2 in 04-decisions.md)
- No files deleted; nothing quarantined

## Update (19:07) — token setup + folder reorganization
- `.env` created with SWISSMODEL_API_TOKEN; gitignored in BOTH repos:
  - co-evolution-analysis/.gitignore:47 → swissmodel/.env ✅ (git check-ignore verified)
  - datasets/.gitignore:37 → co-evolution/swissmodel/.env ✅ (added rule — was NOT ignored)
- Proper folder: co-evolution-analysis/swissmodel/ (script, 2 result JSONs, .env, README.md, plans/)
- Mirrored to datasets/co-evolution/swissmodel/ (working-copy convention)
- Script updated: default FASTA = repo root; token auto-loaded from swissmodel/.env; validated
  (syntax ✓, token load ✓, fasta read ✓, parse 50 templates ✓)
- Coreapi automodel with token: POST /automodel/ → **202**, project 2870bf, status RUNNING
  (poller PID 800391 → /tmp/opencode/coreapi_full_details.json when COMPLETED/FAILED)
- Found pre-existing plans/00-understanding-coevolution.md (earlier session, 03:20) —
  carried into folder, contains PRE-correction numbers (stale, flagged)

## Batch run COMPLETE (Aug 8, 06:20) — 299/299 unique, 1299/1299 originals
- 299 unique sequences (dedup of 1299; server caches by sequence — verified 200+same project_id)
- All 299 automodel jobs COMPLETED, 0 FAILED (1 transient retry mechanism armed, unused)
- 1,128 PDB + 1,128 CIF files, 850.4 MB, in batch/pdbs/u####/ (3–8 models per uid)
- Integrity: all 1,128 PDBs pass gzip read; 60/60 sample have ATOM records (median 26,814 atoms)
- Quality: best-model GMQE 0.69–0.72 for all 299 (all > 0.6); templates: 8hxk.1.B (218), 8d55.1.C (204), 8d56.1.C (204), 6zgf.1.C (93), 7cn8.1.C (81)…
- Timing: submissions ~13 min; total runtime ≈ 75 min (server processed ~6 jobs/min peak)
- Rate limits respected: 30 sub/min, ≤60 poll/min; zero 429s encountered
- Bugs found & fixed: verify() keyed uid vs md5 (manifest said pending for all); .pdb.gz suffix match; pending formula
- Files: batch/manifest.json (all 1299 → status/file), batch/results_index.json (299 summaries), batch/state/*.json (checkpoints)

## Independent re-verification (Aug 8, 06:4x) — user asked "make sure all 1299 are done"
- Built verify_all.py: re-parses FASTA from scratch (1299 seqs), md5 → uid → disk, validates
  every PDB (non-empty, full gzip read, >1000 lines), orphan-folder check. Trusts nothing but disk.
- Canonical repo: 1299/1299 valid PDB models; 1,128 unique PDB + 1,128 CIF; 299/299 folders OK; exit 0
- datasets mirror: identical result; exit 0
- Both FASTA copies byte-identical (md5 b11f4bdbd0f155febc6c1686b7def80d)
- Verifier bug fixed during this pass: total-file counter counted per original sequence
  (5,401 inflated) → now counts unique files (1,128 correct); per-sequence checks were always right
- verify_all.py copied to datasets mirror

## Per-entry materialization (Aug 8, ~07:00) — user asked "why only 299, shouldn't it be 1299?"
- Answer: 1,299 FASTA entries = 299 UNIQUE sequences; server caches by sequence (verified:
  identical resubmit -> 200 + same project_id), so 299 jobs produce the complete model set.
- Added materialize_entries.py: one folder per FASTA ENTRY (1,299), hardlinked to its uid's
  model files. Zero extra disk (4.0 GB would-be -> 828 MB actual).
- Verified: 1,299/1,299 entry folders, 0 without PDB, 0 invalid; duplicate entries correctly
  share their uid (spot-checked); naming entry_####_<accession>, all accessions unique+safe.
- Mirrored to datasets/co-evolution/swissmodel (1,299 folders, same sizes).
- Structure: batch/pdbs_by_entry/entry_####_ACC/{*.pdb.gz, *.cif.gz, templates_summary.json, entry.json}

## Redundancy impact analysis (Aug 8) — "is the co-evolution analysis wrong?"
- Empirically compared corrected (aligned, gap=20) methodology on 1,299 vs 299 unique:
  * variable positions: 21 vs 30 (9 extra found on uniques: 348,419,479,480,486,500,503,507,766)
  * co-evolving pairs MI>0.1: 12 vs 33
  * top-12 pairs: IDENTICAL sets, Spearman rho=0.993 on shared values
  * top pair (373,378): 0.8067 -> 0.8976 (MI attenuated by dominant-variant dilution in 1,299)
  * per-pair QM essential counts shift (e.g. (212,216): 6 -> 0) — rule essentiality is
    redundancy-sensitive; global essential core remains 2 rules but the specific rules shift
- NOTE: my first comparison script re-stripped gaps -> reproduced the buggy 1,247/36,875/1.59
  numbers, accidentally re-confirming the A1 diagnosis. Redone properly on aligned columns.
- Conclusion: analysis NOT invalidated — same clusters/pairs/rules qualitatively; but the
  1,299 numbers understate variation (587x dominant variant dilutes MI + hides 9 variable
  positions). Unique-set numbers are the cleaner estimate.
- Saved: batch/redundancy_analysis_1299_vs_299.json
