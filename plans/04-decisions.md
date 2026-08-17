# 04 — Decisions & Discoveries (validation project)

## Decisions
| ID | Decision | Rationale | Rejected |
|----|----------|-----------|----------|
| D1 | Run everything INSIDE co-evolution-analysis (data/, results/, interpretation folder, agents.md, scripts/) | Explicit user instruction | Separate repo |
| D2 | Parameterize scripts via COEVO_FASTA / COEVO_RESULTS env vars, defaults unchanged | Uniform pattern across 20 scripts; backward compatible; git-diffable | Monkeypatching, script copies |
| D3 | Use `__import__("os").environ.get()` inline in patches | Avoids touching import blocks (most scripts lack `import os`) | Adding imports everywhere |
| D4 | gitignore data/ + results/ (derived, reproducible) | Multi-GB artifacts not for git; scripts regenerate them | Tracking binaries |
| D5 | Pfam: analyze full + id90-clustered copies, report both | Pfam is redundant by design; cluster is the standard DCA protocol | Raw only |
| D6 | A2M lowercase → uppercase (keep inserts) | Inserts are real sequence; consistent with gap=20 encoding | Stripping inserts |
| D7 | Contact def: Cβ–Cβ < 8 Å, sep ≥ 6 | Standard DCA convention | Other thresholds (will note sensitivity if time) |
| D8 | PSICOV full scale: core pipeline on all 150; full 23-script suite per strategy S2 (adaptive, measured) | 150×23 ≈ 3,450 runs; runtime measured in Phase 4 | Blind 3,450 runs |
| D9 | Results structure: results/<dataset>/<script>/ with cwd there | Scripts write relative dirs; isolates runs | Shared dir |
| D10 | Quarantine (not delete) replaced scaffold dir → .rigor-trash/ | A9 no-delete rule | rm -rf |

## Discoveries (survey, Aug 8)
- All 20 main scripts share identical `base_dir = Path("/store/.../datasets/co-evolution")`, `fasta_file`, `results_dir` pattern → single uniform patch possible.
- coevolution_shared.load_position_arrays(fasta_path=None) is the single loader → 1-line env override covers all loader-based scripts.
- Most scripts lack `import os` (only 6 have it) → inline `__import__("os")` chosen.
- Scripts already sys.path-insert kmap-sbm-validation/src and n-ary-kmap/src (absolute) → runner only needs repo root.
- Report generators read result JSONs from base_dir — need COEVO_RESULTS-aware reads (patched in Phase 2c individually).

## Open
| Q | Status |
|---|--------|
| torch in fresh uv venv vs /store/shuvam/.venv fallback | test Phase 0 |
| DeepMSA tar internal structure (a3m vs aln) | inspect Phase 3 |
| CASP native structure fetch for DeepMSA targets | try RCSB Phase 3/7; else document |
| Runtime of full 23-script suite on PSICOV sample | measure Phase 4, adapt D8 |

## Full-scale confirmation (Aug 9, ~21:50)
- User confirmed: use FULL scale (all sequences, all positions) — no sampling.
  Accept the wall time (est. 14-18 h total). O(n²) jobs (LOO-CV at 222-486k seqs)
  that cannot finish will be logged as TIMEOUT and DOCUMENTED as infeasible at
  this scale — an honest result, not a skip.

## Incident log
- / (root fs) filled to 100% by 206 GB of FOREIGN leftover files in /tmp/opencode
  (legacy-indices, runs-first, smoke, uwater, underwater-redesign — dated Aug 2-3,
  not from this session). Moved to /store/.rigor-trash/20260809-tmp-opencode-foreign/
  (rsync --remove-source-files, source intact, verified). / now 57% used.
- Retry quarantine bug: --retry-failed quarantined non-ok dirs across ALL datasets
  (missing --dataset filter) — benign (runners re-ran them) but fixed.
