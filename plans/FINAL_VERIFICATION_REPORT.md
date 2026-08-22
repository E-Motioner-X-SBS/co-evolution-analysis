# FINAL VERIFICATION REPORT — K-map ↔ 3D-Structure Encoding Campaign
**Date:** Aug 22, 2026 · **Verdict: COMPLETE — DOUBLE AUDIT PASSED**

## 1. Deliverables inventory (all present, all verified)

| # | Deliverable | Location | Verified by |
|---|---|---|---|
| D1 | Core library | `scripts/kmap_structure.py` | 26/26 adversarial tests |
| D2 | Stage runner (5 stages) | `scripts/run_contact_campaign.py` | import smoke + all stages executed |
| D3 | Feature caches ×150 | `results/kmap_structure/features/*.npz` | integrity validation, 0 corrupt |
| D4 | CV results | `results/kmap_structure/cv_*.json`, `holdout_split.json` | md5 determinism ×6 total runs |
| D5 | Sensitivity sweep | `sensitivity_sweep.json` | 9-point grid, invariant headline |
| D6 | PATH E Gray test | `path_e_gray_adjacency.json` | p=1.1e-77; permutation control z=2.51 |
| D7 | PATH F map-QM | `path_f_map_qm.json` | 6/6 direction, Wilcoxon p=0.031 |
| D8 | Scientific write-up | `3D_Structure_Analysis/11_Kmap_Structure_Encoding.md` | claims==JSON spot-check exact |
| D9 | Lean proofs (12 thms) | `../lean_proofs/proofs/is_kmap_possible/ContactCircuits.lean` | `lake build` ✔; axiom audit clean |
| D10 | Documentation updates | agents.md §9.5 + banners; lean index; plans/11* | grep-verified |

## 2. Headline scientific results

1. **Transfer**: universal QM contact circuit prec@L/5 = **0.179** on held-out
   proteins (5-fold protein-level CV) vs plain MI 0.097, base 0.034,
   PSICOV-published 0.727. Circuit beats MI in every fold and every
   separation band (long-range: 0.137 vs 0.077).
2. **Gray adjacency**: enriched **1.178×** among 47,430 native contacts,
   two-sided p = 1.1×10⁻⁷⁷ [CI 22.27–23.40%]. Code-permutation control:
   z = 2.51 over random relabelings ⇒ encoding-specific signal confirmed.
   Resolves the Spike underpowered negative (same point estimate, n=38→47,430).
3. **Map-as-Boolean-function**: real contact maps need fewer prime implicants
   than separation-matched shuffles (2.31 vs 2.74 per contact, 6/6 proteins,
   Wilcoxon p = 0.031).
4. **Literature calibration**: our GT pipeline reproduces PSICOV's published
   L/5 precision to Δ=0.003 (0.727 vs Table-1 ≈0.73). The older repo figure
   "0.44" identified as misattribution; correction banners placed at all 6
   stale sites; historical text preserved.

## 3. Verification evidence summary

| Check | Result |
|---|---|
| Adversarial unit suite | 26/26 (multiple fresh runs) |
| Identity mapping | 149/150 sequence-exact; 1vjkA artifact root-caused & handled |
| Cache integrity | 150/150 npz validated post-race, 0 corrupt |
| Determinism | cv md5 `ebbd6666…` ×6 · gray `78fef01a…` ×3 · mapqm `9845cf72…` ×2 |
| Sensitivity | headline invariant across 9 labeling settings |
| Homology leakage | max pairwise 5-mer Jaccard 0.029 — splits sound |
| Lean | compiles clean; axioms = {propext, Quot.sound, native_decide} only; zero sorryAx |
| Doc-vs-artifact | every quoted number matched to source JSON (pass A spot-check) |

## 4. Audit trail

- Pass A (artifact reconciliation): A1–A6 all green (JSONs parse; claims match;
  Lean green; docs sane; single §9; lean index/lakefile updated).
- Pass B (independent re-run): fresh 26/26; cv md5 reproduced; axiom audit
  captured; agents.md §9.5 added.
- Both passes recorded in `plans/CONTINUATION_STATE.md`.

## 5. Known limits (documented, not defects)

- Precision gap to published DCA (0.18 vs 0.72) — framework detects and
  transfers structure; it does not replace sparse-inverse-covariance DCA.
- Consensus-residue features collapse within-column variation.
- PATH F n=6 (exact-QM scale limit); directional evidence only.
- Static-analysis noise: scipy `.pvalue` stubs (annotated `# type: ignore`),
  sibling-repo runtime import path (annotated).

## 6. Reproduce (commands verified this session)

```bash
PY=/store/shuvam/.venv/bin/python; export OMP_PROC_BIND=FALSE
cd /store/shuvam/E-motioner-X-SBS/co-evolution-analysis
$PY scripts/tests/test_kmap_structure.py                       # 26/26
$PY scripts/run_contact_campaign.py --stage cache              # resumable
$PY scripts/run_contact_campaign.py --stage cv                 # PATH B
$PY scripts/run_contact_campaign.py --stage sweep              # G4
$PY scripts/run_contact_campaign.py --stage gray --boot 2000   # PATH E
$PY scripts/run_contact_campaign.py --stage mapqm              # PATH F
cd ../lean_proofs/proofs/is_kmap_possible/KmapProofs && lake build ContactCircuits
```
