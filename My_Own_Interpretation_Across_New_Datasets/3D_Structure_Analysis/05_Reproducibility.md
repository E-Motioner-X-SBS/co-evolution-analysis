# 05 — Reproducibility: exact commands for every number

Every number in this analysis is produced by a **plain Python script** run from
`co-evolution-analysis/`. No interactive REPL, no one-off inline commands —
each analysis has a single entry point that writes a JSON + prints the table.
This document gives the exact commands and the expected outputs.

## 0. Environment

```bash
PY=/store/shuvam/.venv/bin/python     # torch 2.12.1+cu130, numpy, scipy, matplotlib
export OMP_PROC_BIND=FALSE            # REQUIRED: torch otherwise pins to CPU 0
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
cd /store/shuvam/E-motioner-X-SBS/co-evolution-analysis
```

## 1. Rule extraction (the input to all structural analyses)

```bash
$PY scripts/extract_rules_inventory.py
# -> results/contacts/rules_inventory.json
# Expected summary (14 position-pair sources, extended Aug 17):
#   master_boolean 10, full_length 5, run_allseq 30, position_kmap 30,
#     allseq_constraint 10, variable_position 10, perplexity 17, gpu_full 10,
#     dca_mf 20, dca_boolean 10, predictive_constraint 10,
#     flipped_positions 10, kmap_boolean 2, network 8
#   n_residue_rules 36 | n_essential 2 | n_flipped 117
#   n_nary_motifs 20 | n_boolean_motifs 20
```
(Perplexity pairs are computed from the alignment via the shared module —
`combined_pair_scores` — because the stored perplexity summary was empty;
this is the same function every script's combined section uses.)

## 2. The 9 structural analyses

```bash
# 1. Position-level (8 sources × 3D, intra+inter, per-pair + pooled)
$PY scripts/validate_3d_complete.py
# -> results/contacts/validate_3d_complete.json
# Expected: models used 299; pooled control ~0.0956;
#   position_kmap enrich 2.08x p=1.3e-103; per-pair (407,410) 9.3x p=1.5e-226 ...

# 2. Residue-level (36 PIs; fires-on-rule validation)
$PY scripts/validate_residue_rules_3d.py
# -> results/contacts/validate_residue_rules_3d.json
# Expected: control ~0.1187; [210,I,215,P] 2.85x p=1.1e-06; pooled 0.58x

# 3. Flipped / forbidden (steric hypothesis)
$PY scripts/validate_flipped_3d.py
# -> results/contacts/validate_flipped_3d.json
# Expected: control ~0.1278; (212,215) 2.71x p=6.4e-07; pooled 0.25x -> RARITY

# 4. Gray bit-flip adjacency
$PY scripts/gray_flip_analysis.py
# -> results/contacts/gray_flip_analysis.json
# Expected: H=1 observed 9/38 (0.237) vs null 0.200 -> 1.18x p=0.345

# 5. Structural profile (burial / housekeeping)
$PY scripts/structural_profile.py
# -> results/contacts/structural_profile.json
# Expected: Spearman(entropy, burial) rho=-0.170 p=6.4e-09; conserved 0.519 buried, variable 0.260

# 6. Trimer interface (inter-chain)
$PY scripts/trimer_interface.py
# -> results/contacts/trimer_interface.json
# Expected: zero inter-chain contacts for all pairs

# 7. Cross-script consistency
$PY scripts/consistency_cross_script.py
# -> results/contacts/cross_script_consistency.json
# Expected: (212,215) 6 sources 3.6x; 5-source cluster non-structural

# 8. Perplexity as structural predictor
$PY scripts/perplexity_3d_analysis.py
# -> results/contacts/perplexity_3d_analysis.json
# Expected: Spearman(ratio, contact) rho=0.670 p=0.006; MI rho=0.266 p=0.338

# 9. PSICOV 150 cross-dataset contact validation
$PY scripts/validate_contacts.py
# -> results/contacts/contact_validation.json (+ per-target files)
# Expected: 150 targets; mean enrichment 2.49x; mean precision 0.066

# 10. DCA vs MI 3D comparison (NEW — Aug 17)
$PY scripts/validate_dca_3d.py
# -> results/contacts/validate_dca_3d.json
# Expected: DCA 0.87x (1 structural pair (454,495) 5.97x);
#           MI 1.82x (3 structural pairs); perplexity ρ=0.670
```

All scripts import the shared core `scripts/rules3d_common.py` (model
enumeration, mapping integrity, PDB parsing, distances, controls, binomial
statistics) — so the definitions of "contact", "control", and "p-value" are
single-sourced and identical across analyses.

## 3. Determinism of the runs

1. **Fixed seeds**: random controls use `random.Random(42 + model_count)` /
   `Random(7 + model_count)` / `Random(11 + model_count)` per model — exact
   reproducible controls.
2. **Exact statistics**: scipy `binomtest` (exact binomial, one-sided
   greater) and `spearmanr` — no approximations.
3. **Deterministic core**: re-running `master_boolean.py` (no env) reproduces
   the corrected rule set exactly (36 PIs / 2 essential; byte-identical
   summary) — regression-verified.
4. **Byte-compare verification**: re-running master_boolean/kmap_boolean/
   flipped on gpcrdb/b2ar produced byte-identical JSON outputs.
5. **Coverage audit**: `$PY scripts/verify_coverage.py` — 0 truncation flags
   across 3,732 runs (all sequences, full length, empirically proven).

## 4. The multi-day pipeline runs (context)

The 23-script suite was run on 6 datasets with the parameterized runner:
```bash
$PY scripts/run_dataset.py --dataset <ds> --psicov-mode full --parallel 12
# per-dataset results under results/<ds>/<target>/<script>/ with run_log.json
# resume: --retry-failed (quarantines failed dirs and re-runs them)
```
3,687+ verified runs; heavy jobs (LOO-CV at 222k–486k sequences) run per the
no-timeout directive and are documented as infeasible at full scale (O(N²);
27 days – 2.5 years of compute per job — see `03_Kmaps_and_Boolean_Rules.md`
§4). Every run is reproducible from `results/**/run_log.json`.

## 5. Expected-output checksums

The JSONs are deterministic; re-running any analysis script must reproduce the
JSON byte-for-byte (same seeds, same inputs). A quick verification:
```bash
$PY scripts/validate_3d_complete.py > /tmp/a.json.log 2>&1   # rerun
# compare printed "models used: 299" and the per-pair table to 04_Structural_Analyses.md
```
