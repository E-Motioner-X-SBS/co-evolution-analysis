# 02 — Multi-Dataset Campaign Preservation Status

**Verified:** Aug 17, 2026 (re-run of all three audit scripts).

## Campaign scope
The 23-script co-evolution pipeline was run on 6 external benchmark datasets,
fully parameterized via `COEVO_FASTA` / `COEVO_RESULTS` / `COEVO_OUT` env vars
(defaults unchanged → regression-verified on the original Spike analysis).

| Dataset | Targets | Role |
|---------|---------|------|
| unit_tests | 3 (DHFR, 1whzA, 1atzA) | pipeline smoke tests |
| gpcrdb | 2 (5ht1a, b2ar) | curated membrane proteins |
| evmutation | 20 | mutation-effect benchmark (up to 28k seqs) |
| pfam | 3 (PF00072/00595/00071) | deep diverse families (486k/330k/222k seqs) |
| psicov150 | 150 | **contact-prediction benchmark (native PDBs)** |
| deepmsa | 614 | structural ground truth (seqs + natives only; no MSAs shipped) |

## Audit results (all green)

### B1 — Coverage audit (`scripts/verify_coverage.py`)
- **3,735 runs** audited · **0 truncation flags**
- Every run used ALL sequences and FULL alignment length (no hardcoded caps).
- Exit 0.

### B2 — Completeness audit (`scripts/check_completeness.py`)
- Every dataset × target × script is accounted for.
- The only non-`ok` runs are the **6 documented-infeasible Pfam giants**:
  - `pfam/PF00072|00595|00071` × `gpu_full_analysis.py` → exit1 (GPU OOM at the
    8 GB per-process cap; CPU would need 89 GB RAM for the one-hot at 486k seqs).
  - `pfam/PF00072|00595|00071` × `dca_mf_analysis.py` → exit1 (same GPU OOM).
  - `pfam/PF00595/allseq_constraint_function.py` → NO RUN LOG (the O(10·N²)
    LOO-CV at 330k seqs ≈ 27 days – 2.5 years of single-core compute; runs per
    the no-timeout directive, or documented sampling is the only practical route).
- These are **honest limits, not silent failures** — documented in
  `plans/04-decisions.md`, `agents.md` §6, and the 3D_Structure_Analysis docs.

### B3 — GPU/CPU consistency (`scripts/check_gpu_cpu.py --all --max-seq 5000`)
- **5/5 PASS** (one sample per dataset).
- Max entropy diff ≈ 4.4e-7 · 0 reference mismatches · max MI diff ≈ 4.1e-7 ·
  top-10 rank agreement 10/10 on every dataset.
- Reference-residue tie-breaking unified to lowest code between CPU and GPU.

## Scripts present & runnable
- `scripts/`: 24 tools (runner, downloader, converters, 3D analyses, audits).
- `py_compile`: **35/35 pass** (24 analysis + 11 tooling scripts).
- Reproducibility: every run writes `run_log.json` (status, seconds, stdout/stderr
  tail, output file list) — all 3,735 runs are auditable.

## Data preserved
- `data/`: 4.7 GB — all 6 datasets + `metadata.json` + `manifest.json` (integrity
  verified: sizes, magic bytes, gzip -t, tar members, sequence counts).
- `results/`: 1.7 GB — `results/<dataset>/<target>/<script>/` with outputs +
  `run_log.json`; `results/contacts/` (all 3D-structure JSONs);
  `results/consistency/`; `results/coverage_report.json`.
- `swissmodel/batch/`: 956 MB — 1,128 PDB + 1,128 CIF models (299 unique Omicron
  sequences, best-GMQE 0.69–0.72, all > 0.6; 0 failures).

## How to re-run any dataset
```bash
PY=/store/shuvam/.venv/bin/python
export OMP_PROC_BIND=FALSE
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
cd /store/shuvam/E-motioner-X-SBS/co-evolution-analysis
$PY scripts/run_dataset.py --dataset <ds> --psicov-mode full --parallel 12
# --retry-failed quarantines failed dirs and re-runs; --retry-failed-cpu forces CPU
```

## Conclusion
The multi-dataset validation campaign is **fully preserved and reproducible**:
scripts present and compile, results stored with per-run logs, all three audits
green, and the only incomplete runs are the explicitly-documented infeasible
giants (a quantified limitation, not a gap).
