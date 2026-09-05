# 03 — Steps: The Rigid, Looping Execution Plan

> **CORRECTION (Aug 22, 2026):** the 'PSICOV published L/5 = 0.44' quoted below is a misattribution; the paper's PSICOV L/5 ≈ 0.73 (we reproduce 0.727). See 11_Kmap_Structure_Encoding.md §9-S1.

Master to-do lives in the session tool (todowrite); this file is the ground truth of WHAT and HOW.
Phases are executed in order; the VERIFICATION LOOP (Phase 8) may send us back to any phase.

## PHASE 0 — Environment  [todo #4]
- uv venv `.venv` python 3.13 in co-evolution-analysis; install numpy scipy matplotlib requests (+torch cu129 attempt)
- Verify: `python -c "import numpy, scipy, matplotlib; print(ok)"`, GPU check `torch.cuda.is_available()`
- Fallback documented: /store/shuvam/.venv (torch 2.12.1+cu130, numba 0.66)
- PASS: venv importable + GPU available (or documented fallback used)

## PHASE 1 — Downloads  [todo #5]
- scripts/download_datasets.py (written; resumable; integrity per file: size, magic, gzip -t, tar members, seq headers)
- All 6 datasets → data/<dataset>/ ; manifest.json (file, bytes, ok)
- PASS: manifest all ok; spot re-check of tar members + sequence counts (2nd independent check)

## PHASE 2 — Parameterization  [todo #6]
- 2a. Patch coevolution_shared.load_position_arrays: COEVO_FASTA env override
- 2b. Patch 20 scripts: fasta_file + results_dir env-aware (uniform, `__import__("os")` inline)
- 2c. Patch specials: dca_mf (FASTA, OUT), gpu_full (FASTA, OUT), generate_co-evolution_md (fasta + OUT), generate_full_analysis_md (OUT), generate_full_pipeline_doc (OUT)
- 2d. REGRESSION: no-env runs of master_boolean + full_length_analysis + kmap_boolean → outputs identical to git HEAD (diff); with-env run on unit-test DHFR → outputs in COEVO_RESULTS
- PASS: both regression directions green; git diff shows ONLY the parameterization lines

## PHASE 3 — Conversion + metadata  [todo #7]
- scripts/convert_all.py: per dataset → data/<dataset>/converted/<name>.fasta (+ per-protein files for PSICOV)
- Stockholm parser (strip # lines), PSICOV aln parser (header=filename), A2M (uppercase), a3m (uppercase, drop reference-insert lines), FASTA passthrough
- Pfam: full + id90-clustered copy (greedy; report stats)
- metadata.json per dataset: n_raw, n_converted, unique seqs, redundancy ratio, length range, source, paper, PDB ids
- Structures: extract PSICOV pdb/; fetch Pfam reps (6TNE, 6AK2, 1HE1), GPCR (7E2X, 2RH1), unit (4P3R, 1WHZ, 1ATZ), DeepMSA CASP natives (by target) from RCSB
- PASS: converted sequence counts == source counts (per-file logs); every sequence parses; metadata written

## PHASE 4 — Analysis runs  [todo #8]
- scripts/run_dataset.py: for each dataset × script: env COEVO_FASTA=..., COEVO_RESULTS=results/<dataset>/<script>/, cwd=that dir, subprocess run, capture exit+stdout+stderr → run_log.json + status
- PSICOV: per-protein loop for core scripts (all 150); full suite per strategy S2
- PASS: every script exit=0 or logged FAILED-with-reason; output files exist per script (parse script outputs to confirm non-empty JSONs)

## PHASE 5 — GPU/CPU consistency  [todo #9]
- scripts/check_gpu_cpu.py: entropy, refs, MI top-pairs GPU vs CPU per dataset sample
- PASS: max abs diff < 1e-6 on samples; table in results/consistency/

## PHASE 6 — Rules/equations/K-map tables  [todo #10]
- Verify per dataset: master_boolean JSONs (PIs, essential), kmap_boolean (equations), flipped (forbidden rules), nary, perplexity outputs
- Generate consolidated MD per dataset: K-map tables (CSV/MD), prime implicants, minimized Boolean equations (main/n-ary/flipped/perplexity) → My_Own_Interpretation_Across_New_Datasets/<dataset>/rules_and_tables/
- PASS: every dataset has rules + equations files; spot-verify equation ↔ truth table (evaluate equation on on-set cells)

## PHASE 7 — Contact-map validation  [todo #11]
- scripts/validate_contacts.py: contacts from native PDBs (Cβ<8Å, sep≥6); score top-MI pairs + rule pairs; precision@L/5; enrichment permutation test; compare vs published (PSICOV 0.44 etc.)
- PASS: per-dataset table + honest comparison (works / doesn't work, with reasons)

## PHASE 8 — Verification loop (the core of rigor)  [todo #12]
- 8a. Re-run pass: 1 protein/dataset (or whole small datasets) → byte-compare all output files (md5) with first run
- 8b. Sanity audit: no NaN/Inf in JSONs; MI in [0, ~5]; entropy in [0, log2(20)]; rule labels map to real residues
- 8c. Anomaly protocol: any mismatch/NaN/unexpected → reproduce, inspect parsing (converted FASTA vs source), inspect script logic, fix (either converter or script), re-run affected dataset, re-audit
- PASS: all re-runs byte-identical; zero uninvestigated anomalies; log of every anomaly + resolution

## PHASE 9 — Interpretation  [todo #13]
- My_Own_Interpretation_Across_New_Datasets/<dataset>/NN_<script>.md following the existing 01..23 style (what it does, algorithm, formulas, worked example, results, inference, honest notes)
- Cross-dataset summary: which algorithms validate, which fail, why
- PASS: interpretation files for every dataset×script run; honest failure notes present where applicable

## PHASE 10 — Final audit + agents.md  [todo #14]
- Macro-audit vs 00-understanding success criteria; update plans (04-decisions, 05-audit-log)
- agents.md updated with final env/run/provenance/behavior facts
- Final honest report to user

## Loop rule
Any phase that fails its PASS gate → fix at the correct level (converter, script, parameterization, data) → re-run from that phase → re-audit. Never proceed on an uninvestigated anomaly.
