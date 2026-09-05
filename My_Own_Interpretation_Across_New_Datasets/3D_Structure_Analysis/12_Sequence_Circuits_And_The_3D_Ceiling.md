# 12 — Sequence-as-Circuit, Contact-Map Blocks, and the Information Ceiling

**Date:** Sep 5, 2026 · **Status:** all four stages executed, gates green
**Engine:** `scripts/sequence_circuits.py` · **Gates:** `scripts/tests/test_sequence_circuits.py` (22/22)
**Driver:** `scripts/run_sequence_circuit_campaign.py` · **Outputs:** `results/sequence_circuits/*.json`
**Equations:** [`BOOLEAN_CIRCUIT_EQUATIONS.md`](../../BOOLEAN_CIRCUIT_EQUATIONS.md)
**Formal layer:** `lean_proofs/proofs/is_kmap_possible/SequenceCircuits.lean` (compiles clean, sorry-free)

---

## The three questions, and the three answers

| # | Question | Answer |
|---|---|---|
| Q1 | Can a protein sequence be written as an exact minimal Boolean circuit? | **Yes, provably and losslessly** — verified by round-trip on all 1,299 Spike sequences |
| Q2 | Does a 3D contact map carry real Boolean/K-map structure? | **Yes, decisively** — segment×segment blocks, 86× over separation-matched shuffles, p = 1.5e-25, 150/150 proteins |
| Q3 | Can a circuit over sequence features *predict* 3D contacts well enough to fold? | **No — and it is capped, not merely underperforming.** See §4 |

The three answers are not in tension. The sequence circuit is *faithful*; the
contact circuit is *structured*; the map between them is *information-starved*.
That last point is the new result, and it is quantitative.

---

## 1. The sequence as a circuit (R-B) — exact and lossless

Define, for a sequence `s` of length `L` over the 20-letter alphabet:

```
f_s : {0,1}^p × {0,1}^b → {0,1}        p = ⌈log₂ L⌉ = 11,  b = 5
f_s(i, a) = 1   ⟺   i < L  ∧  a = Gray(He(s_i))
```

**Proposition 1 (faithfulness).** `s ↦ f_s` is injective, and `s` is recoverable:
`s_i = He⁻¹(Gray⁻¹(the unique a with f_s(i,a) = 1))`.
Formalised as `sc_cell_injective` (Lean, exhaustive over 64×32 cells).
Verified by execution on **1,299/1,299 Spike sequences: `all_lossless = true`.**

**Proposition 2 (canonicity).** The prime-implicant set of a Boolean function is
unique, so `PI(f_s)` is a well-defined invariant *of the sequence*.

**Proposition 3 (structure theorem).** Every prime implicant of `f_s` has a
**fixed residue field** — no cube is free in a residue bit.
*Proof.* Cells `(i,a)` and `(i',a')` with `a ≠ a'` and `i ≠ i'` differ in at
least one position bit and at least one residue bit, hence Hamming distance
≥ 2, hence never merge. Since each position carries exactly one residue,
`a ≠ a'` forces `i ≠ i'`. ∎
Formalised as `sc_no_cross_residue_merge` (Lean, 2²⁰ cases).
**Verified empirically: 0 of 64,122 prime implicants across 60 Spike sequences
has a free residue bit.**

**Corollary 3.1.** `PI(f_s) = ⋃_a {maximal dyadic subcubes of P_a}` where
`P_a = {i : s_i = a}`. The minimal circuit is exactly a *per-residue dyadic
block cover of the sequence*. Measured on Spike sequence 1: 1,025 cubes across
20 residues (L: 84, S: 78, T: 73, V: 70, …).

### The honest negative

Faithful is not the same as informative. Circuit **size** carries essentially
no biological signal:

| quantity | Spike (n=1299) | composition-matched shuffle | uniform random |
|---|---|---|---|
| cubes in minimal cover | 1013.2 | 1020.2 | 1042.6 |
| compression (residues/cubes) | **1.2531** | **1.2429** | **1.2183** |
| dipeptide-circuit cubes | 41.1 | 41.7 | 25.3 |

A real Spike sequence is **0.8 % more compressible than its own shuffle**. The
minimal circuit is a canonical re-encoding of the sequence, not a discovery
about it. Anyone hoping that "sequence → K-map → QM" would expose folding
information *in the compression statistic* should stop here: it does not.

---

## 2. The contact map as a circuit (R-D) — where the structure is

Now encode the *structure* instead. With plain-binary position indexing:

```
C_X : {0,1}^p × {0,1}^p → {0,1}      C_X(i,j) = 1 ⟺ ‖r_i − r_j‖ < 8 Å ∧ |i−j| ≥ 6
```

**Proposition 4 (segment lemma).** A cube whose free bits are exactly the low
`u` bits of the i-field and the low `v` bits of the j-field asserts

```
[i₀, i₀+2^u) × [j₀, j₀+2^v)  ⊆  contacts
```

— every residue of one chain segment touching every residue of another.
Formalised as `sc_low_free_is_interval`, `sc_interval_span`, and
`sc_contact_cube_is_block` (Lean, exhaustive).

**Proposition 5 (encoding invariance).** The same holds under Gray-coded
positions: because `high bits of g(i) = g(high bits of i)`, the preimage is
again an aligned dyadic interval, merely relabelled. Formalised as
`sc_gray_cube_is_also_interval`; the relabelling is `sc_gray_relabels_segments`.

> **Correction to an earlier draft of this work.** The first version of
> `sequence_circuits.py` asserted that Gray coding *destroys* the interval
> structure and that plain binary was therefore required. That is false, and
> Lean refused the proof. Binary is preferred only because its decoding is the
> identity. This is the same kind of encoding-invariance already recorded in
> `AntiCorrelation.lean`.

### Result — 150 PSICOV targets, 5 separation-preserving shuffles each

Controls hold the contact count and the entire `|i−j|` distribution fixed and
randomise only *which* residues touch.

| statistic | real maps | shuffled | real > shuffled | Wilcoxon p |
|---|---|---|---|---|
| compression (contacts / cubes) | **1.7450** | 1.2246 | **150/150** | 2.3e-26 |
| fraction of contacts inside segment blocks | **0.2059** | 0.0024 | 145/150 | 1.5e-25 |

**Real contact maps put 20.6 % of their contacts inside contiguous
segment×segment blocks; separation-matched random maps manage 0.24 % — an 86×
enrichment**, with a mean of 32 block cubes per protein. Every circuit is
runtime-verified sound and complete.

This supersedes PATH F, which reached the same qualitative conclusion on **6**
proteins (`L ≤ 62`, Wilcoxon p = 0.031) because the dense reference minimiser
could not go further. The sparse engine lifts the limit to the whole benchmark
(L = 50–266) and the p-value by 24 orders of magnitude.

**Interpretation.** This is the Boolean-algebraic signature of secondary
structure. A β-sheet pairing or a helix–helix packing *is* a fully-connected
bipartite block between two chain segments, which is precisely what a large
prime implicant of `C_X` encodes. The K-map formalism sees real folds.

---

## 3. Co-evolution circuits (R-C) — three polarities, and a trap

Full equations in [`BOOLEAN_CIRCUIT_EQUATIONS.md`](../../BOOLEAN_CIRCUIT_EQUATIONS.md) §1.
Reproduces the published conventions exactly: **21 variable positions**
(H > 0.3) — matching the corrected repo figure — and 11 co-evolving pairs at
mutation-only MI > 0.1 (the published set is 10; the extra pair sits at the
threshold).

The `standard` polarity is clean and biologically legible. For the top pair
(373, 378), MI = 0.803, the entire circuit is four essential terms:

```
f₃₇₃,₃₇₈ = (L·T) ∨ (F·A) ∨ (S·T) ∨ (S·A)
```

### The don't-care trap (this is a methodological finding, not a result)

A Spike column carries only 2–4 distinct residues, so ~99 % of the 32×32 grid
is unobservable. If those cells are marked **don't-care**, Quine–McCluskey
absorbs the whole region and returns rules like

```
forbidden = {D,E,F,M,N,Q,W,Y} × {C,D,E,G,H,K,N,P,Q,R,S,T}      (8 free bits)
```

which is formally sound and scientifically vacuous — it is a statement about
the don't-care mass, not about the protein. Marking them **off** instead
confines cubes to the observed alphabet and yields testable rules. Both
policies are reported for every pair so the sensitivity is on the record.

This is the same degeneracy that sank Phase 1 of
[`10_Boolean_Contact_Circuits.md`](10_Boolean_Contact_Circuits.md); it is a
general hazard of applying QM to sparse biological K-maps, and it is worth
stating in the paper.

---

## 4. The information ceiling — why better minimisation cannot help

Here is the part that answers Q3 as a *bound* rather than a benchmark.

A circuit whose inputs are (residue at i, residue at j, separation bin) is **a
function of its cell**. Every candidate pair landing in the same cell receives
the same score, so no minimisation algorithm, threshold, don't-care policy, or
cover-selection rule can rank one above another. The achievable precision is
therefore capped by the *purity of the cells*, which is a property of the
feature map alone — computable by giving each cell its own pooled empirical
contact rate, i.e. handing the circuit the test labels.

Because a cell holds many pairs, the quantity that matters is the expectation
over within-cell orderings: a function of the cell cannot tell those pairs
apart, so it cannot break that symmetry. That is `achievable` below. (The
`tie-break max` column resolves ties in favour of contacts; it is vacuous for
coarse partitions — `sep_only` hits 0.976 with four cells purely by sorting on
the label — and is shown only to indicate how much the tie-break matters.)

150 PSICOV targets, Cβ–Cβ < 8 Å, separation ≥ 6:

| feature set | cells | **achievable** | tie-break max |
|---|---|---|---|
| identity + separation + MI quartile | 6390 | **0.2412** | 0.2448 |
| identity + separation (= PATH B level-2) | 1600 | **0.1933** | 0.2084 |
| 8-group chemistry + separation (= level-1) | 256 | **0.1418** | 0.2349 |
| identity only | 400 | 0.1383 | 0.2684 |
| Gray Hamming distance + separation | 24 | 0.0789 | 0.1397 |
| separation only | 4 | 0.0598 | 0.9760 |
| MI quartile only | 4 | 0.0430 | 1.0000 |

Against the same benchmark: **PATH B circuit-L2 holdout = 0.179**, plain
MI = 0.097, **PSICOV published DCA = 0.723**.

### The three numbers that matter

- **PATH B reaches 0.179 against an achievable 0.1933 — 93 % of its feature
  class's ceiling.** The Boolean machinery is essentially saturated. The
  remaining headroom from *any* improvement in minimisation, thresholding,
  don't-care policy or cover selection is at most ~8 % relative.
- **The ceiling itself is 3.7× below DCA** (0.1933 vs 0.723).
- **Adding an MI quartile lifts the ceiling to 0.2412** — the only feature that
  moves it appreciably, and still 3× short of DCA.

**What this settles.** The gap between the K-map circuits and published DCA is
**not a minimisation failure and cannot be closed by better Boolean
machinery.** It is the feature class running out of information. Exact QM,
Espresso, Petrick-optimal covers, smarter don't-care policies, n-ary K-maps —
none of them can move a predictor above its cell-purity ceiling.

---

## 5. What would actually raise the ceiling

The ceiling is a property of the *partition*, so the only lever is a finer
partition that separates contacting from non-contacting pairs **within** the
current cells. Ranked by expected value:

1. **Column statistics, not residue identity.** Adding an MI quartile already
   moves the ceiling (§4). Real DCA couplings would move it much further,
   because the coupling is computed from the whole alignment rather than from
   a consensus letter. This is why PATH L's rank fusion (circuit + mfDCA →
   prec@L/5 0.204, p = 1.8e-5 over raw DCA) is the one direction in this repo
   that beat its inputs.
2. **Frequency-weighted rather than consensus residues.** Collapsing a column
   to its majority letter discards the distribution the coupling lives in —
   already flagged as limitation #2 in
   [`11_Kmap_Structure_Encoding.md`](11_Kmap_Structure_Encoding.md).
3. **Chemistry is not the lever.** The 8-group physicochemical partition has a
   *lower* ceiling than plain residue identity (§4): grouping merges cells and
   can only lose information. Charged/aromatic/polar classes make rules
   readable, not stronger.

**On folding.** Even at the top of this range the pipeline stays below the
`path_i_reconstruction.json` verdict rule (median long-range prec@L/5 ≥ 0.30
**and** recall@L ≥ 0.10). Contact ranking is not structure prediction: there is
still no distance-geometry or folding stage anywhere in the organisation
(`ml_folding_models/src/` and `Is-Kmap-Possible/src/contact_kmap/` are empty).

---

## 6. Reproduce

```bash
PY=/store/shuvam/.venv/bin/python; export OMP_PROC_BIND=FALSE
cd /store/shuvam/E-motioner-X-SBS/co-evolution-analysis

$PY scripts/tests/test_sequence_circuits.py                      # 22/22 gates
$PY scripts/run_sequence_circuit_campaign.py --stage seq         # 1,299 sequences
$PY scripts/run_sequence_circuit_campaign.py --stage coevo       # circuit equations
$PY scripts/run_sequence_circuit_campaign.py --stage contact     # 150 PSICOV, 5 shuffles
$PY scripts/run_sequence_circuit_campaign.py --stage ceiling --with-mi
$PY scripts/generate_circuit_equations_md.py                     # -> BOOLEAN_CIRCUIT_EQUATIONS.md

cd ../lean_proofs/proofs/is_kmap_possible
lean SequenceCircuits.lean                                       # clean, sorry-free
```

Seeds: folds 42, controls/shuffles 20260905. All stages deterministic.

---

## 7. Two corrections to existing repo claims

### 7.1 The reference minimiser is not exact QM

`kmap_sbm.analysis.prime_implicants` keeps only don't-cares at Hamming distance
1 from a real on-minterm — a performance choice stated in its own docstring.
With clustered don't-cares (the contact-campaign regime, where DC = "cell has
n < n_min") it **under-merges**. Measured on a 10-variable table, 60 on-cells,
400 don't-cares:

| | exact QM | reference |
|---|---|---|
| prime implicants | **522** | 187 |

and every one of the 74 cubes the reference reports that exact QM does not is
**strictly inside** a true prime (`test_sequence_circuits.py::t06d`). The rules
stay sound; the *counts* and the *"essential"* designations are approximations.
Published PI counts derived from it — including the headline "36 distinct prime
implicants (2 essential)" — should be read as lower bounds under an
approximation, not as exact QM quantities.

### 7.2 "Axiom-free" is incorrect for most theorems

Verified with `#print axioms`:

| theorem | proof method | axioms |
|---|---|---|
| `cc_sep_bin_mono` | `omega`/induction | `propext, Quot.sound` |
| `cc_fixed_match_unique` | `native_decide` | `propext, …native_decide.ax_1_1` |
| `cc_cover_complete` | `native_decide` | `propext, …native_decide.ax_1_1` |
| `cc_padding_safety` | `omega` over a `native_decide` lemma | `propext, Classical.choice, Quot.sound, …native_decide.ax_1_1` |

`native_decide` discharges goals through the compiler and introduces a trust
axiom. **90 of the 118 theorems in `is_kmap_possible` use it** (KmapProofs
16/16, AminoAcidEncoding 25/27, ContactMapCompleteness 16/20, KmapEncodingEquiv
7/8, KmerIndexing 20/35, ContactCircuits 6/12), as do all the new ones here.

"Sorry-free" is accurate. "Axiom-free" is not, and appears in at least five
places: `kmap-sbm-validation/RESULTS.md:15` and `:178`, `lean_proofs/agents.md:11`,
`co-evolution-analysis/FULL_PIPELINE_ANALYSIS.md:27` (and its copy in
`datasets/co-evolution/`), `skills/kmap_sbm_validation.md:52`. Suggested
wording: *"no user-declared axioms; `native_decide` theorems additionally
depend on `ofReduceBool`/compiler trust."*
