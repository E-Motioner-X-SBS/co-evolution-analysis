# 02 — Strategy: Multi-Dataset Validation inside co-evolution-analysis

## User decision (binding)
Everything lives INSIDE `/store/shuvam/E-motioner-X-SBS/co-evolution-analysis/`:
- `data/<dataset>/` — downloads + converted FASTA + metadata + structures
- `results/<dataset>/<script>/` — all script outputs
- `My_Own_Interpretation_Across_New_Datasets/` — interpretations
- `agents.md`, `plans/`, `scripts/` (converters/runners/checkers)

## Core strategy decisions

### S1. Parameterization (needed to run the 23 scripts on new FASTA)
All 20 main scripts share the identical pattern (surveyed):
```
base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
fasta_file = base_dir / "Spike_protein.aln-fasta"
results_dir = base_dir / "<script_results_name>"
```
Patch uniformly (backward-compatible; defaults UNCHANGED when env absent):
```
fasta_file  = Path(os.environ.get("COEVO_FASTA")   or (base_dir / "Spike_protein.aln-fasta"))
results_dir = Path(os.environ.get("COEVO_RESULTS") or (base_dir / "<name>"))
```
- Use `__import__("os").environ.get(...)` inline → no import-line changes needed.
- Special cases: dca_mf_analysis.py & gpu_full_analysis.py (`FASTA`/`OUT`), generate_full_analysis_md.py (`BASE = Path(__file__).parent` + `OUT`), generate_full_pipeline_doc.py (`OUT`), generate_co-evolution_md.py (`fasta_file` + output path).
- Shared loader: coevolution_shared.load_position_arrays(fasta_path=None) → env override (1 line; `import os` already present).
- **Regression gate**: with no env vars, re-run 3 scripts → outputs byte-identical to current files (git diff clean). With env → outputs land in COEVO_RESULTS.

### S2. "Full scale" definition for PSICOV 150 (150 proteins × 23 scripts ≈ 3,450 runs)
Adaptive, honest, recorded:
- **Core pipeline (per protein, ALL 150)**: entropy/MI (full_length_analysis), co-evolving pairs + master_boolean rules, kmap_boolean, flipped, nary, perplexity → results/<dataset>/<protein>/.
- **Full 23-script suite**: on all 150 if wall-clock allows (measure in Phase 4; expected ~minutes/protein on A100) — otherwise on a stratified representative subset (10 proteins by depth×length), documented.
- Report generators: per dataset once.

### S3. Conversion choices (Phase 3)
- Stockholm → FASTA (parse all sequences; keep full width; gaps `-`; strip `#` lines).
- PSICOV aln → FASTA with header = filename (PDB id).
- A2M → FASTA with lowercase uppercased (inserts kept; documented).
- DeepMSA → inspect tar; a3m → FASTA (uppercase+lowercase→uppercase; remove insert lines relative to reference if a3m format).
- Pfam: produce BOTH full and id90-clustered copies (cluster by CD-HIT if available, else greedy 90% by sequence identity — implement simple greedy in numpy; report cluster stats).
- Metadata JSON per dataset: source, paper, n_seqs (raw/converted/clustered), unique counts, length, redundancy ratio, PDB ids.

### S4. Contact-map validation (the scientific core)
- Contact definition: Cβ–Cβ (Cα for Gly) distance < 8 Å, residue separation ≥ 6 (standard DCA convention).
- Per protein: compute contact map from native PDB; score our top-MI pairs and Boolean-rule pairs by precision (fraction on contacts); compare to published top-L/5 precision (PSICOV: 0.44 mean).
- Enrichment test: are co-evolving pairs significantly more often contacts than random pairs of same length distribution? (hypergeometric / permutation).
- Rule→contact translation: does a Boolean rule (pos_i=aa_i ∧ pos_j=aa_j) mark positions that are in contact? Report per dataset.
- Sources: PSICOV native PDBs (in tar); Pfam reps (fetch 6TNE/6AK2/1HE1 from RCSB); GPCRdb (7E2X/2RH1); DeepMSA CASP natives (fetch by target id); EVmutation/unit: where structure available (4P3R/1WHZ/1ATZ; EVmutation per-protein), else sequence-only (documented).

### S5. GPU/CPU consistency
- Script `scripts/check_gpu_cpu.py`: on a sample per dataset (3 proteins or the family), compare GPU vs CPU: entropy vector, majority refs, MI matrix (top pairs) → max abs diff; assert < 1e-6; report table.

### S6. Honesty protocol
- Every script run records: exit code, stdout/stderr tail, output file list, duration → run_log.json per run.
- No silent failures: any exit≠0 or empty output → flagged, investigated, fixed or documented as a known limitation.
- Interpretation docs report failures honestly.

### S7. Environment
- uv venv `.venv` (python 3.13); `uv pip install numpy scipy matplotlib requests`; torch: try `--index-url https://download.pytorch.org/whl/cu129`; fallback: use `/store/shuvam/.venv/bin/python` (torch 2.12.1+cu130, numba) — documented in agents.md.
- sys.path needs: repo root (co-evolution-analysis), kmap-sbm-validation/src, n-ary-kmap/src (scripts already insert the latter two; runner adds repo root).

## Rejected alternatives
- Separate validation repo (user overrode: keep inside co-evolution-analysis).
- Monkeypatching scripts at import time (fragile) vs editing scripts (reviewable, git-diffable, backward-compatible) → edit.
- Copying scripts per dataset (duplication) → single parameterized set.
