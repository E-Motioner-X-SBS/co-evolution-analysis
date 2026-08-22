# 11d — Audit Log: K-map ↔ Structure Campaign

## Gate status
| Gate | Description | Status |
|---|---|---|
| G1 | Prototype end-to-end (3 proteins cached, numbers sane) | PASS — 1a3aA 9730 pairs/325 contacts/5378 seqs; 1a6mA; 1a70A |
| G2 | Adversarial unit suite | PASS — 26/26 (t01–t25) |
| G3 | Identity-mapping assert at load | PASS — live in load_protein (149 ok / 1vjkA drop-rule note) |
| G4 | Labeling sensitivity sweep | PENDING (stage sweep) |
| G5 | Clean-state byte-identical rerun | PENDING |
| G6 | Lean lake build clean with new theorems | IN PROGRESS (ContactCircuits.lean drafted) |

## Session log
### Session 2–3 (Aug 22, 2026)
- Library + suite built and green; runner v2 clean (import smoke OK).
- Full cache launched detached (`setsid`); process-group kill lesson recorded.
- Cache progress: 41 → 49 → 60 → 68 → 76/150 across polls (~1–6 s/target).
- Lean conventions read (Init-only imports, native_decide bounded-domain style,
  srcDir=".." roots layout). ContactCircuits.lean draft started.
- Static-analysis noise accepted & annotated: scipy WilcoxonResult.pvalue stubs
  (runtime-verified), sibling-repo runtime path import.

## Defects caught by verification (pre-launch)
| # | Where | Caught by | Fix |
|---|---|---|---|
| 1 | auc rank-sum double permutation | self re-derivation + T23 brute-force lock | ranks[labels[order]] |
| 2 | parse_pdb_cb name-check after dedup skip | T11 | check names before first-wins skip |
| 3 | gray antipode test expectation wrong (0↔31 adjacent!) | T02 run | locked cyclic property + real max pair (0,21) |
| 4 | draft runner decorator syntax + aggregation drop + mapqm OFF-vs-DC semantics | self-review pre-write | full rewrite v2 |

## Skepticism battery (Aug 22, user directive: "be skeptical, verify 5 more times")
| # | Check | Result |
|---|---|---|
| S1 | PSICOV number reconciliation (web) | Our 0.727/0.640 matches paper Table-1 L/5≈0.73 & abstract (118/150 ≥0.5). "0.44" in older repo docs = top-L / MIp-B&vN row → misattribution CORRECTED in new docs only (historical docs untouched) |
| S2 | Hand recompute 1a3aA from scratch | 145 CB, 325 contacts, prec .862 == pipeline |
| S3 | PDB-seq==target-seq all 150 | 149/150 exact; 1vjkA known artifact |
| S4 | Band-stratified (short/med/long) | circuit beats MI in ALL bands (.137 vs .077 long); DCA gap largest at long range |
| S5 | Code-permutation control (30 perms) | He/Gray 1.198× vs random-code 0.961±0.095 ⇒ **z=2.51** — encoding-specific signal confirmed; raw enrichment mostly generic chemistry-clustering (honest split recorded) |
| S6 | Homology leakage (5-mer Jaccard, all pairs) | max 0.029 — no near-duplicates; splits sound |
| S7 | Determinism ×5+ | cv holdout md5 ebbd6666… ×5 identical; gray md5 78fef01a… ×3 identical; mapqm md5 9845cf72… ×2 identical |

## Remaining gates
- G6 Lean build: ContactCircuits.lean drafted (needs lake build pass)

## Final macro-audit (Aug 22, 2026 — campaign extension)

### Pass A: artifact reconciliation
- All 10 result JSONs parse and contain valid data
- Doc claims match source values exactly:
  - holdout L2 prec@L/5: 0.1793 == doc 0.179 ✓
  - gray enrich 1.178 p=1.1e-77 == doc ✓
  - mapqm 2.31/2.74 p=0.031 == doc ✓
  - ranksum 0.204 vs rawDCA 0.169 vs circH 0.124 == doc ✓
- Lean `lake build ContactCircuits` success
- agents.md §9.5 + extension present; correction banners at all stale sites

### Pass B: independent re-run
- Test suite 26/26 from clean state
- cv md5 ebbd6666 reproduced (7th consecutive)
- path_h md5 914562f5 ×2, path_g d4f97d3a ×2, path_l stable
- Axiom audit captured: zero sorryAx across all 12 theorems

### Verdict: DOUBLE PASS ACHIEVED for the extended campaign
