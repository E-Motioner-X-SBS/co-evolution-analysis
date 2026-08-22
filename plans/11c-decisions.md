# 11c — Decisions, Discoveries & Findings: K-map ↔ Structure Campaign

## Decisions
| ID | Decision | Rationale | Rejected |
|----|----------|-----------|----------|
| D1 | Ground truth = native PDBs + PSICOV published predictions (literature), NOT our SWISS-MODEL Spike models | User directive; native structures are experimental, identity col↔residue mapping is exact by construction | Reusing Spike homology models |
| D2 | Unit of analysis = (protein, column-pair); features per-protein-intrinsic | No cross-protein alignment needed; 149 proteins ≈ 1.5M labeled pairs | Per-sequence labels (structure is fixed per protein) |
| D3 | Contact def Cβ–Cβ<8Å, sep≥6 (literature convention) | Matches validate_contacts.py, Morcos 2011, Jones 2012 | sep≥3 (rules-analysis convention) — reported as sensitivity where relevant |
| D4 | 8-group scheme identical to build_boolean_contact_circuits.py | Continuity with prior Phase-2 circuits; Lean's 7-group scheme noted as related-but-distinct | Inventing a new grouping |
| D5 | Cell labeling rule: rate≥T·base ∧ n≥n_min → 1; n<n_min → DC; else 0; defaults T=2,n_min=30 with full sweep | Simple, documented, sensitivity-tested (gate G4) | Majority vote (contacts too rare ⇒ all-0 tables) |
| D6 | Protein-level CV folds stratified by length quintile, seeded | Prevents homolog-protein leakage between train/test; deterministic | Random pair-level splits (leakage) |
| D7 | Local MI/perplexity implementations mirror shared-module ACTUAL behavior; equivalence property-tested | Avoids private-API coupling; discrepancy documented not silently propagated | Importing private _mi_dense |
| D8 | 1vjkA handled by drop-rule (PDB residues > width ignored) rather than exclusion | Preserves the protein; exact on the aligned region; generalizes | Excluding 1vjkA entirely |

## Discoveries
| ID | Finding | Impact |
|----|---------|--------|
| R1 | raw/README guarantees aln width == PDB CA count; verified 149/150 live | Eliminates NW/gap-mapping error entirely — biggest integrity upgrade vs Spike work |
| R2 | 1vjkA has ONE extra C-terminal PDB residue ('S',88) absent from target seq (cloning artifact); cols 1..87 identical | Drop-rule resolves it exactly |
| R3 | con/*.out = PSICOV published predictions `i j 0 8 score` sorted desc; cols 3-4 documented dummies | Literature contact map usable directly as baseline |
| R4 | All 150 PDBs have CB coverage ≥50% of CA (0 flagged); numbering contiguous 1..L in all files | Cβ contact definition usable everywhere; no altloc handling needed |
| R5 | coevolution_shared.perplexity_ratio DOCSTRING says "None if fewer than 3 conditioning residues" but CODE returns None only if zero such residues | Doc/code mismatch in prior repo code (F1 below). We mirror actual behavior; prior results unaffected (condition rarely binds) |

## Findings / issues in existing theory (user asked us to surface these)
| ID | Where | Issue | Action taken |
|----|-------|-------|--------------|
| F1 | coevolution_shared.perplexity_ratio | docstring vs code disagree on the None condition | Mirrored actual code semantics locally; recorded here; NOT changing shared module mid-campaign (reproducibility doctrine) — proposed fix queued for post-campaign review |
| F2 | run_kmap_analysis.analyze_coevolution | hardcoded `max_len = min(500, min_len)` truncates correlation analysis to first 500 positions, contradicting the repo's no-truncation guarantee (coverage audit regex does not catch it because the printout says "Alignment length analyzed") | Logged; fix scheduled after campaign freeze (changing it now would invalidate stored kmap_results JSONs mid-audit) |
| F3 | auc_mann_whitney (this campaign, caught pre-test) | initial rank-sum indexing double-permuted (`ranks[order][labels[order]]`) | Fixed to `ranks[labels[order]]`; brute-force concordance property test added to lock it |

## Open questions
| Q | Status |
|---|--------|
| QM runtime on densest level-2 maps (1024 cells) | measure in prototype |
| Do circuit rules beat PSICOV published predictions on held-out proteins? | PATH B/C experiment |
| Is Gray adjacency enriched among contacts at scale? | PATH E |
