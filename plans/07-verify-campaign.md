# 07 — Verify the multi-dataset validation campaign is preserved & runnable

## Already verified (this session)
- scripts/ present: 24 tools (runner, downloader, converters, 3D analyses, audits).
- py_compile: 35/35 pass (exit 0).
- results/ structure: 178 targets (unit 3, gpcrdb 2, evmutation 20, pfam 3, psicov150 150).
- run_log.json: 3,735 runs (3,729 ok, 6 exit1). The 6 exit1 = pfam PF00071/595/072 × {gpu_full, dca_mf} = documented infeasible GPU-OOM giants. Honest, preserved.
- coverage_report.json: 3,732 runs, 0 truncation flags.
- consistency_report.json: 5/5 PASS (entropy/refs/MI max diff ~4e-7).
- data/: all 6 datasets + metadata.json + manifest.json (4.7 GB).
- swissmodel/batch: 1,128 PDBs (956 MB).

## Steps (B1–B5)
- B1. Re-run verify_coverage.py → expect 0 flags. Save summary.
- B2. Re-run check_completeness.py → expect every dataset×target×script accounted.
- B3. Re-run check_gpu_cpu.py --all --max-seq 5000 → expect 5/5 PASS.
- B4. Confirm the 6 exit1 are documented infeasible (not silent failures): cross-check vs plans/04-decisions + agents.md limits table.
- B5. Write a preservation-status summary into My_Own_Interpretation_Across_New_Datasets/ (a new 02_preservation_status.md).

## Test
- All 3 audit scripts exit 0; summaries match the numbers above.
- No new flags; no regressions vs the audit-log entries.
