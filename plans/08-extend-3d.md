# 08 — Extend the 3D structural validation to ALL 23 scripts' pair outputs

## Current state
rules_inventory.json position_pairs = 8 sources:
master_boolean(10), full_length(5), run_allseq(30), position_kmap(30),
allseq_constraint(10), variable_position(10), perplexity(17), gpu_full(10).

## Missing pair-emitting scripts (JSON structure verified)
| Script | Result JSON | Pair field | Add as |
|---|---|---|---|
| dca_mf_analysis | dca_results/dca_mf_summary.json | top_20_di (pos_i,pos_j,di) | dca_mf |
| dca_boolean_coevolution | dca_boolean_results/dca_boolean_summary.json | results[].pos_i/pos_j | dca_boolean |
| advanced_co-evolution_analysis | advanced_analysis_results/coevolution_network.json | edges[].source/target | network |
| predictive_constraint_function | constraint_function_results/constraint_function_summary.json | results[].pos_i/pos_j | predictive_constraint |
| flipped_boolean_coevolution | flipped_boolean_results/flipped_boolean_summary.json | rules[].pos_i/pos_j (unique) | flipped_positions |
| kmap_boolean_coevolution | kmap_boolean_coevolution/boolean_functions.json | rules[].pos_i/pos_j | kmap_boolean |
| run_kmap_analysis | kmap_results/h1_adjacency_results.json | consecutive Gray-adjacent (1.6M) | h1_adjacency (distinct category) |

Non-position scripts (documented, cannot map to 3D): coevolution_shared,
coevolution_gpu (modules); boolean_co-evolution, nary_kmap (motifs, no
positions); create_mi_heatmap (viz); generate_* (report generators).

## Steps (C1–C6)
- C1. Extend scripts/extract_rules_inventory.py to load the 6 new position-pair
  sources (+ h1 as a documented non-position category). Output updated
  rules_inventory.json with 14 position-pair sources.
- C2. Re-run extract_rules_inventory.py → verify 14 sources, pair counts.
- C3. Re-run scripts/validate_3d_complete.py on the expanded inventory →
  results/contacts/validate_3d_complete.json with all 14 sources scored
  against the 299 models (intra+inter, per-pair + pooled, binomial p).
- C4. Add a new script scripts/validate_dca_3d.py for the mfDCA top-DI pairs
  specifically (the strongest external-method comparison): top-20 DI pairs →
  3D contact enrichment + precision, mirroring validate_contacts.py but on
  the Omicron models.
- C5. Run C4; save results/contacts/validate_dca_3d.json.
- C6. Update results/contacts cross-script consistency with the expanded set.

## Test
- rules_inventory.json has 14 position-pair sources; counts match JSONs.
- validate_3d_complete.json reports all 14 sources; per-pair table complete.
- DCA 3D validation JSON present with enrichment + precision.
