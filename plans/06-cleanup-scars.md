# 06 — Cleanup: remove pre-correction "scars" from the repo

## Problem
Six top-level report files still carry pre-correction (buggy) numbers and
conflict with CORRECTION_NOTICE.md (the verified truth: 21 var / 10 pairs /
36 PIs / 2 essential / 490 forbidden / LOO-CV 9.24% / max MI 0.807).

## Stale files (verified by grep + reading)
| File | Stale content | Disposition |
|---|---|---|
| CO-EVOLUTION_BOOLEAN_FUNCTIONS.md | 57 var, 1161 pairs, 108 rules | quarantine; regenerate via generate_co-evolution_md.py |
| RESULTS_ANALYSIS.md | 1249 var, 36918 pairs, MI 1.5917, 152 rules, LOO 2.93/7.26 | quarantine; rewrite as corrected audit |
| FULL_COEVOLUTION_ANALYSIS.md | "18 scripts", "3 rules across 3 pairs" | quarantine; regenerate via generate_full_analysis_md.py |
| FULL_PIPELINE_ANALYSIS.md | "3 essential", "19/20/17 scripts", LOO 0.08%, MI 1.59 | quarantine; regenerate via generate_full_pipeline_doc.py |
| FULL_PIPELINE_ANALYSIS.html | older snapshot (Jul 17) | quarantine; regenerate if generator emits html, else drop |
| COEVOLUTION_CONSTRAINTS.md | old example positions (76/77) | quarantine; rewrite with corrected examples |

## Steps (A1–A7)
- A1. mkdir .rigor-trash/20260817-stale-reports/ + INDEX entry (A9 no-delete).
- A2. Move the 6 files into the trash dir (keep relative paths).
- A3. Regenerate FULL_COEVOLUTION_ANALYSIS.md:
  `COEVO_RESULTS=datasets/co-evolution COEVO_OUT=$PWD/FULL_COEVOLUTION_ANALYSIS.md $PY generate_full_analysis_md.py`
- A4. Regenerate FULL_PIPELINE_ANALYSIS.md:
  `COEVO_RESULTS=datasets/co-evolution COEVO_OUT=$PWD/FULL_PIPELINE_ANALYSIS.md $PY generate_full_pipeline_doc.py`
- A5. Regenerate Boolean-functions report:
  `COEVO_RESULTS=datasets/co-evolution $PY generate_co-evolution_md.py`
  → produces datasets/co-evolution/kmap_boolean_coevolution/COEVOLUTION_KMAP_BOOLEAN.md (corrected).
  Write a corrected top-level CO-EVOLUTION_BOOLEAN_FUNCTIONS.md that mirrors it.
- A6. Rewrite RESULTS_ANALYSIS.md as the corrected audit (21/10/36/2/490/9.24%/0.807).
- A7. Rewrite COEVOLUTION_CONSTRAINTS.md with corrected example positions.
- A8. Verify: grep stale tokens across repo root *.md → only agents.md (legit bug-doc context).

## Test
- `grep -rn "1,249\|36,918\|1.5917\|152 rules\|108 rules\|1,161\|57 of 80\|7.26%" *.md` → clean (except agents.md troubleshooting rows).
- Corrected numbers present in each regenerated file.
