# 02 — Contact Computation: what we look at in the 3D structure

This document defines precisely **what data the 3D structure contributes**, how
it is read from the PDB files, how sequence positions are mapped to model
residues, and how "contact" is defined and validated.

---

## 1. What data lives in the PDB files

PDB (Protein Data Bank) format, fixed-width columns (1-indexed):

```
ATOM      1  N   GLU     1      18.222  18.496 -16.203  1.00  0.00           N
01234567890123456789012345678901234567890123456789012345678901234567890123456789
```

| Columns (0-indexed) | Field | Used for |
|---|---|---|
| 0–3 | record name (`ATOM`) | only `ATOM` records |
| 12–15 | atom name (`CA`, `CB`, `N`…) | select `CA`/`CB` |
| 16–19 | residue name (`GLU`…) | map to 1-letter code (20 canonical) |
| 20 | chain id (blank for PSICOV's repaired models) | chain assignment (trimer A/B/C) |
| 21–25 | residue sequence number | residue identity within a chain |
| 29–37, 37–45, 45–53 | x, y, z coordinates | distance computation |

We keep only **Cβ coordinates** (Cα for glycine, which has no Cβ) of the 20
canonical amino acids. The result per model:

```
coords: {chain: {residue_number: (x, y, z)}}
```

- Omicron models: 3 chains (A, B, C) × ~1,270 residues each (homo-trimer).
- PSICOV repaired natives: 1 chain with **blank** chain id → treated as chain A.

## 2. Contact definition

Two residues `(r_a, r_b)` in the same chain are **in contact** iff

```
d(r_a, r_b) = sqrt( (x_a-x_b)^2 + (y_a-y_b)^2 + (z_a-z_b)^2 ) < 8.0 Å
|r_a − r_b| >= 3          (rules analysis; backbone-adjacent excluded)
|r_a − r_b| >= 6          (DCA-convention controls; PSICOV benchmark)
```

- The 8 Å Cβ–Cβ cutoff is the standard residue-contact definition used in the
  co-evolution literature (e.g., Morcos et al. 2011 PNAS; Jones et al. 2012).
- The sequence-separation cutoffs remove trivial backbone-neighbour contacts:
  residues at separation 1–2 are always "close" regardless of structure.
  The rules analysis uses ≥3 (we only exclude truly adjacent residues so that
  short-loop pairs like (212,216) — separation 4 — are testable); the
  benchmark comparison uses ≥6 (the published convention).

## 3. Mapping sequence positions → model residues (the critical step)

The co-evolution rules are expressed in **alignment columns** (0-based,
1,276 columns of the Omicron alignment). The models were built from the
**gap-stripped target sequence** (SWISS-MODEL strips gaps before modelling).
The mapping:

```
residue(col) = col + 1 − (number of non-canonical characters in columns < col)
```

i.e. every gap `-` and ambiguous character `X` before column `col` shifts the
residue number down by one. Implemented as:

```python
def col_to_residue_map(aligned_seq, target_seq):
    # integrity: the model's target sequence MUST equal the aligned sequence
    # with all non-canonical characters removed
    if gap_strip(aligned_seq) != target_seq:
        return {}, False            # model skipped, never silently mapped
    m = {}
    res = 1
    for col, ch in enumerate(aligned_seq):
        if ch not in CANONICAL:     # gaps AND ambiguous chars
            continue
        m[col] = res
        res += 1
    return m, True
```

**Integrity guarantee**: for every one of the 299 unique sequences, the model's
`target_sequence` (from `full_details.json`) was verified equal to the
gap-stripped aligned sequence; mismatching models were skipped and logged
(initially 96 skipped because `X` was not stripped — fixed, then 0 skipped).

## 4. The random control (how significance is established)

For a given rule set with position pairs `P`, the control is:

1. per model, collect the rule-pair residue separations
   `S = {|residue(i) − residue(j)| : (i,j) in P}` (after mapping, sep ≥ cutoff);
2. draw 2,000 random pairs `(a, b)` per model with `a` uniform over the model's
   residues and `b = a + s`, `s` uniform over `S` (matched separation
   distribution);
3. count contacts among the random draws → pooled control rate `p0`.

Then for a pair with `n` mappable models and `c` contacts:

```
enrichment  = (c / n) / p0
p_value     = P(X >= c)  with X ~ Binomial(n, p0)      (exact, scipy binomtest,
              alternative="greater")
```

Enrichment > 1 with p < 0.05 = the pair's residues are significantly closer in
3D than random residue pairs of the same separation distribution.

## 5. Chain handling (trimer models)

- **Intra-chain contact**: residues in the SAME chain (the natural monomer
  view) — chain A used as the reference chain.
- **Inter-chain contact**: residue `i` in chain X vs residue `j` in chain Y,
  X ≠ Y; the minimum over all chain pairs. This tests the trimer interface —
  a distinct structural dimension (the Spike's RBD-up/down conformational
  machinery involves inter-protomer contacts).

## 6. What is NOT used (honest limitations)

- **dssp secondary-structure annotations** were requested but are not present
  in the Omicron `full_details.json` payloads → secondary-structure context
  was not available from this source (documented; could be computed from
  geometry with DSSP as a follow-up).
- **Solvent accessibility** is not stored → approximated by **local Cβ
  density** (number of Cβ within 12 Å) as a burial proxy (see
  `04_Structural_Analyses.md`, analysis 5).
- **MD dynamics** (SMOG/SBM) are NOT used in this analysis — contacts are
  static-geometry based. SMOG-based dynamic validation is a planned follow-up
  via the `kmap-sbm-validation` repository.

## 7. Validation of the structural test bed itself

- Model quality: all 299 unique sequences have best-model GMQE 0.69–0.72
  (all > 0.6), QSQE up to 0.947 (homo-trimer), built with ProMod3 3.6.0 from
  templates such as 8d55.1 (BA.2 closed-state Spike, 98.3% identity).
- Target-sequence integrity verified per model (§3).
- GPU/CPU numerical consistency: entropy/refs/MI identical to 1e-7 across all
  5 datasets (5/5 PASS).
- The PSICOV 150 native experimental structures provide an independent
  structural ground truth (native, not modelled).
