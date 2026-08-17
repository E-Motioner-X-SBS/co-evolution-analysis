# 04 — Decisions & Discoveries: SwissModel test

## Decisions
| ID | Decision | Rationale | Rejected |
|----|----------|-----------|----------|
| D1 | Test the guest `/interactive` web path first (no account) | No token available; coreapi requires auth (401 verified) | Waiting for user token before any test |
| D2 | Parse `const globalTemplateData` JSON from the results page | Server-embedded, complete, structured; fields verified (50 fields/template) | Parsing HTML tables (lossy, no trg/tpl seq) |
| D3 | Deliver `swissmodel_template_search.py` in `co-evolution-analysis/` | Reusable for any of the 1,299 sequences; FASTA lives there | Keeping test in /tmp only |
| D4 | coreapi mode implemented but not executed (no token) | Contract fully verified from OpenAPI + official wrapper; runs once token provided | — |

## Discoveries
| ID | Finding | Impact |
|----|---------|--------|
| F1 | coreapi `/automodel/` → 401 without token (empirically verified) | Account needed for official API |
| F2 | Old anonymous `POST /interactive/template_identification` removed (404 × 4 URLs) | Must use current paths |
| F3 | Guest path: `POST /interactive` (fields target/project_title/email/automodel/is_alignment + CSRF) creates project; empty `multiTargetValidate` triggers "Target sequences must be unique." | Omit that field |
| F4 | `GET /interactive/{pid}/models/automodel_poll` returns JSON: tpl_search_status, tpl_count, progress steps, model_ids | Structured progress monitoring |
| F5 | Results page embeds `const globalTemplateData = [...]` — full template objects (50 fields: pdb_id, qname, chain, method, seq_id, coverage, seq_sim, found_by, resolution, oligo_state, qsqe, pred_lddt, trg_seq, tpl_seq, ...) | Canonical parse source |
| F6 | Spike seq1 (WRU87367.1): **50 templates**, top 8d55.1 (BA.2 Spike, 98.3% id, 99.9% cov); identity 76–99%, coverage 93–100%; 41 HHblits + 9 BLAST; all EM | Template matching works for Omicron data |
| F7 | Same sequence re-submitted → new project, search NOT cached across projects; runtime ~7–15 min for 1,269 aa | Budget 15 min/seq for batch runs; rate limits 200/min, 10000/6h |
| F8 | Crambin validation: 9 templates (1 AFDB rank0 + 8 PDB NMR/X-ray via HHblits) | Parser + flow correct on small case |

## Reversals
| ID | What changed | From → To | Why |
|----|--------------|-----------|-----|
| R1 | Initial plan was coreapi-only | → guest path as primary test | Auth blocked it; guest path returns same template search |

## Open
| ID | Q | Status |
|----|---|--------|
| O1 | Full automodel (model building) results for a Spike seq via coreapi | Needs token (O2) |
| O2 | SWISS-MODEL account/token availability | Ask user |
| O3 | Batch run for all 1,299 sequences (~15 min each → ~325 h serial; needs parallelization + polite rate limits) | User decision |

## Batch run (Aug 8) — decisions
| ID | Decision | Rationale |
|----|----------|-----------|
| D5 | Dedupe 1,299 seqs → 299 unique before submitting | Verified server caches by sequence (resubmit seq0 → 200, same project_id 2870bf); 299 jobs ≪ 1299 |
| D6 | Pacing: 30 submissions/min, ≤90 status checks/min, 429 → Retry-After/60s backoff | Docs: 200/min rapid, 10000/6h prolonged |
| D7 | Crash-safe checkpoints in batch/state/*.json (atomic tmp+rename) | Runs for hours; resumable after any interruption |
| D8 | One `run` orchestrator: submit → poll → auto-download PDB/CIF → manifest | Single background process (PID 4059787) |
| D9 | Store per-unique-sequence folder batch/pdbs/u####/ (PDB, CIF, model JSON, full_details.json) | Proper structure; manifest maps all 1,299 originals → uid |
| D10 | Seed submitted.json with existing 2870bf (seq0) | No duplicate job |
| D11 | `run` capped at 168 h; verify phase reports done/failed/pending per original sequence | Honest completion accounting |
