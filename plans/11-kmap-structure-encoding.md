# 11 — Understanding: The K-map ↔ 3D-Structure Encoding Campaign

**Date:** 2026-08-22
**Directive (user, restated):**
1. Use ALL prior binary-level machinery (Gray encoding → K-map → Quine-McCluskey → prime
   implicants) and all co-evolution information + Lean proofs + 3D structures.
2. **Do NOT be dependent on our own computations**: the literature already has contact maps,
   connections, structures. Use literature ground truth itself.
3. Build a K-map giving an *intrinsic view* of how the chain
   `binary sequence → Gray code → QM → prime implicants → K-map` can **encode 3D structure data**.
4. Try as many different paths as possible. "This must be possible, or maybe some way it is
   possible. You need to verify it."
5. Free to download new datasets/things.
6. Rigorous verification of everything; if theory has issues or better ways exist, fix in THIS repo.
7. Document every line; verify mathematically with Lean.
8. Long-running: days are fine. "We will be finding some pattern."

## Why this can now succeed where the Spike circuits failed

The Boolean contact circuits (`scripts/build_boolean_contact_circuits.py`, doc 10) failed on the
Spike for three documented reasons: sparsity (21 variable positions × 2–3 residues each),
redundancy (587 identical copies), lineage domination. The fix is **data**, and the literature
provides it locally:

| Asset | What it gives | Status |
|---|---|---|
| `data/psicov150/raw/aln/` | 150 deep diverse MSAs (511–74,836 seqs) | local ✓ |
| `data/psicov150/raw/pdb/` | 150 **native experimental** PDBs (repaired, renumbered) | local ✓ |
| `data/psicov150/raw/con/` | **PSICOV's published contact predictions** (Jones 2012) — the literature contact map | local ✓ |
| `data/psicov150/raw/seq/` | target sequences | local ✓ |
| README (raw/) | format spec + **the identity-mapping guarantee** | read ✓ |

### THE CRITICAL STRUCTURAL FACT [VERIFIED from raw/README]

> "The alignments correspond to the sequence represented in each PDB file … So there are exactly
> the same number of columns in each alignment as C-alpha ATOM records in each PDB file."

⇒ **MSA column c ↔ PDB residue c (1-based), exactly, by construction.** No Needleman-Wunsch,
no gap-mapping error, no model/target mismatch. This eliminates the largest error source of the
Spike validation (which needed gap-count mapping through homology models).

### con/ format [VERIFIED from raw/README]
`i j min_dist(dummy=0) max_dist(dummy=8) score`, one row per candidate pair, sorted by score.
Columns 3–4 are explicitly documented dummies ("ignore!").

## Success criteria
- [ ] A verified pipeline that, per protein, builds **contact-labeled K-maps** from MSA × native
      contacts and QM-minimizes them into Boolean circuits (per-pair + universal).
- [ ] Universal cross-protein circuit (residue-group × group × separation → contact) evaluated
      with **protein-level holdout CV** across 150 proteins (the test Phase 2 failed with n_test=3).
- [ ] Three-way comparison: our circuit pairs vs **literature PSICOV predictions** vs native contacts.
- [ ] Properly powered redo of Gray-adjacency↔contact test (was n=38, p=0.345 on Spike).
- [ ] Inverse predictability: how much of the contact Boolean function do sequence stats recover?
- [ ] New Lean theorems covering the new objects; `lake build` clean.
- [ ] Every script fully commented; results + honest interpretation written to
      `3D_Structure_Analysis/11_*`; agents.md updated; audits green.

## Open questions
- Q1: Do repaired PSICOV PDBs contain CB atoms or CA-only? (affects contact def) → check Phase R.
- Q2: Are all 150 aln widths == PDB CA counts? (identity mapping must hold per protein) → assert.
- Q3: QM runtime on densest contact K-maps (1024 cells) → measure on prototype.
