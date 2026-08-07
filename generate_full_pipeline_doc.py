

def combined_pipeline_section():
    """Combined MI + perplexity markdown section (all experiments)."""
    lines = ["## Combined MI + Perplexity Analysis (all experiments)",
             "",
             "Both lenses are computed for every variable-position pair:",
             "MI (total correlation) and perplexity ratio (determinism).",
             "Combined score = average of normalized ranks.",
             ""]
    try:
        from coevolution_shared import (load_position_arrays,
                                        compute_entropy_vectorized,
                                        combined_pair_scores)
        pa, na, fl = load_position_arrays(max_pos=None, aligned=True)
        ent = compute_entropy_vectorized(pa, na, fl)
        var = [p for p in range(fl) if ent[p] > 0.3]
        pairs = [(i, j) for idx, i in enumerate(var)
                 for j in var[idx + 1:] if j - i <= 30]
        sc = combined_pair_scores(pa, pairs, na, ent)
        lines.append(f"| Variable positions | {len(var)} |")
        lines.append(f"| Pairs scored (MI + ratio) | {len(sc)} |")
        lines.append("")
        lines.append("| Rank | Pos i | Pos j | MI | PP ratio | Combined |")
        lines.append("|------|-------|-------|-----|----------|----------|")
        for rk, s in enumerate(sc[:10], 1):
            lines.append(f"| {rk} | {s['pos_i']} | {s['pos_j']} | "
                         f"{s['mi']:.3f} | {s['ratio']:.2f} | {s['combined']:.3f} |")
        mi_top = sorted(sc, key=lambda s: -s['mi'])[:3]
        ratio_top = sorted(sc, key=lambda s: -s['ratio'])[:3]
        lines.append("")
        lines.append("Top 3 by MI: "
                     + ", ".join(f"({s['pos_i']},{s['pos_j']})" for s in mi_top))
        lines.append("Top 3 by PP ratio: "
                     + ", ".join(f"({s['pos_i']},{s['pos_j']})" for s in ratio_top))
        lines.append("Top 3 by combined: "
                     + ", ".join(f"({s['pos_i']},{s['pos_j']})" for s in sc[:3]))
    except Exception as e:
        lines.append(f"(combined analysis unavailable: {e})")
    return "\n".join(lines)


#!/usr/bin/env python3
"""
COMPLETE PIPELINE DOCUMENTATION GENERATOR
==========================================
Generates FULL_PIPELINE_ANALYSIS.md with:
 - All 236 Lean theorems catalog
 - Complete dataset description
 - K-map construction step-by-step with formulas
 - All 108 Boolean expressions (Quine-McCluskey minimized)
 - Entropy, Mutual Information, Perplexity derivations
 - Coupling constants and constraint functions
 - H1-H6 hypothesis results
 - Flipped/don't-care analysis
 - Every co-evolving position pair with full details
"""

import sys, os, json, time
from pathlib import Path
from datetime import datetime
from collections import Counter
import numpy as np

# ── Setup ────────────────────────────────────────────────────────────
BASE = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
sys.path.insert(0, str(BASE))
os.chdir(str(BASE))

from coevolution_shared import (
    load_position_arrays,
    mutual_information,
    compute_entropy,
    compute_entropy_vectorized,
    find_variable_positions,
    majority_ref,
    mi_mutation_only,
    compute_coupling,
    constraint_function,
    build_mutation_kmap,
    AA_LIST,
    N_AA,
)

OUT = BASE / "FULL_PIPELINE_ANALYSIS.md"
now = datetime.now().strftime("%B %d, %Y at %H:%M")

# ── Load all precomputed data ────────────────────────────────────────
print("Loading data...")
# BUG FIX: was max_pos=200 (first 200 positions only). Now FULL length.
pos_arrays, n_all, full_len = load_position_arrays(max_pos=None)

# Boolean expressions from QM minimization
with open(BASE / "kmap_boolean_coevolution/boolean_functions.json") as f:
    bool_data = json.load(f)
all_rules = bool_data["rules"]  # 108 rules

# Master boolean (full-length)
with open(BASE / "master_boolean/master_boolean_summary.json") as f:
    master_data = json.load(f)

# Group rules by position pair
rules_by_pair = {}
for r in all_rules:
    key = (r["pos_i"], r["pos_j"])
    if key not in rules_by_pair:
        rules_by_pair[key] = []
    rules_by_pair[key].append(r)

# ── Compute fresh results ────────────────────────────────────────────
print("Computing metrics...")

# BUG FIX: entropy/MI/consensus were limited to first 80 positions.
# Now FULL length (full_len from load_position_arrays).
full_len_actual = max(len(a) for a in pos_arrays)

# Entropy vector — FULL length
entropy_vec = compute_entropy_vectorized(pos_arrays, n_all, full_len_actual)
perplexity_vec = 2.0**entropy_vec
var_positions = find_variable_positions(pos_arrays, n_all, full_len_actual, 0.3)

# H1 Hamming
ham = np.zeros(6, dtype=np.int64)
tot_consec = 0
for arr in pos_arrays[:n_all]:
    n = 0
    while n < len(arr) and arr[n] >= 0:
        n += 1
    for j in range(n - 1):
        ag = arr[j] ^ (arr[j] >> 1)
        bg = arr[j + 1] ^ (arr[j + 1] >> 1)
        d = bin(ag ^ bg).count("1")
        if d <= 5:
            ham[d] += 1
        tot_consec += 1
h1_ratio = ham[1] / tot_consec if tot_consec else 0

# Top MI pairs — FULL length (window 30, step 1)
mi_pairs = []
for i in range(full_len_actual):
    for j in range(i + 1, min(i + 30, full_len_actual)):
        mi = mutual_information(pos_arrays, i, j, n_all)
        if mi > 0.01:
            ref_i = majority_ref(pos_arrays, i, n_all)
            ref_j = majority_ref(pos_arrays, j, n_all)
            mi_pairs.append((i, j, mi, ref_i, ref_j))
mi_pairs.sort(key=lambda x: x[2], reverse=True)

# Perplexity for top pairs
perp_results = []
for pi, pj, mi, ri, rj in mi_pairs[:20]:
    mi_mut, n_mut = mi_mutation_only(pos_arrays, pi, pj, ri, rj, n_all)
    cond = {}
    for aa in range(N_AA):
        cnt = {}
        total = 0
        for arr in pos_arrays[:n_all]:
            if pi < len(arr) and pj < len(arr) and int(arr[pi]) == aa:
                cj = int(arr[pj])
                if cj >= 0:
                    cnt[cj] = cnt.get(cj, 0) + 1
                    total += 1
        if total > 5:
            h = -sum((c / total) * np.log2(c / total) for c in cnt.values() if c > 0)
            cond[AA_LIST[aa]] = (2**h, total)
    ppm = 2 ** entropy_vec[pj]
    avgc = np.mean([v[0] for v in cond.values()]) if cond else 0
    perp_results.append((pi, pj, mi, mi_mut, ppm, avgc, cond))

# Consensus AA per position — FULL length
pos_freq = np.zeros((full_len_actual, N_AA), dtype=np.float64)
for arr in pos_arrays[:n_all]:
    L = min(len(arr), full_len_actual)
    for pos in range(L):
        if 0 <= arr[pos] < N_AA:
            pos_freq[pos, int(arr[pos])] += 1
pos_freq /= n_all

consensus = []
for pos in range(full_len_actual):
    top = int(np.argmax(pos_freq[pos]))
    consensus.append(
        (pos, AA_LIST[top], float(pos_freq[pos, top]), float(entropy_vec[pos]))
    )

# Frequency K-map for example pair (76,77)
freq_76_77 = np.zeros((N_AA, N_AA), dtype=np.float64)
for arr in pos_arrays[:n_all]:
    if 76 < len(arr) and 77 < len(arr):
        ci, cj = int(arr[76]), int(arr[77])
        if ci >= 0 and cj >= 0:
            freq_76_77[ci, cj] += 1
freq_76_77 /= n_all if n_all else 1

# Coupling for top pairs
couplings = []
for pi, pj, mi, ri, rj in mi_pairs[:20]:
    J, avg_abs = compute_coupling(pos_arrays, pi, pj, n_all)
    top_co, top_anti = [], []
    for ai in range(N_AA):
        for aj in range(N_AA):
            if abs(J[ai, aj]) > 1.0:
                info = (AA_LIST[ai], AA_LIST[aj], float(J[ai, aj]))
                (top_co if J[ai, aj] > 0 else top_anti).append(info)
    top_co.sort(key=lambda x: x[2], reverse=True)
    top_anti.sort(key=lambda x: x[2])
    couplings.append((pi, pj, mi, avg_abs, top_co[:5], top_anti[:5]))

print("Generating markdown...")

# ══════════════════════════════════════════════════════════════════════
# BUILD MARKDOWN
# ══════════════════════════════════════════════════════════════════════

md = []


def w(s=""):
    md.append(s)


t0_start = time.time()

# ── TITLE ────────────────────────────────────────────────────────────
w("# E-Motioner-X-SBS: Complete Co-Evolution Analysis Pipeline")
w()
w(f"**Generated:** {now}")
w(
    f"**Dataset:** SARS-CoV-2 Omicron Spike Protein — {n_all:,} sequences, {full_len} positions"
)
w(f"**Compute:** NVIDIA A100 80GB + 24-core Xeon, Python 3.10")
w(f"**Author:** Shuvam Banerji Seal — IISER Kolkata")
w()
w("---")
w()

# ── ABSTRACT ─────────────────────────────────────────────────────────
w("## Abstract")
w()
w(
    "This document presents the complete computational pipeline for analyzing co-evolutionary"
)
w(
    "constraints in the SARS-CoV-2 Spike protein using Karnaugh map (K-map) Boolean minimization."
)
w(
    f"**{n_all:,}** Omicron variant Spike sequences from GISAID were encoded using base-20"
)
w(
    "(He 2012 ordering) amino acid representation. Position-pair mutual information identified"
)
w(
    f"**{len(rules_by_pair)}** co-evolving position pairs in the N-terminal signal peptide region"
)
w(
    f"(positions 68-79). Quine-McCluskey Boolean minimization produced **{len(all_rules)} essential"
)
w("prime implicants** — each representing an irreducible co-evolutionary constraint.")
w()
w(f"**Key results:** H1 Gray-code adjacency enrichment = **{h1_ratio:.2f}×**; ")
w(
    f"max pairwise MI = **{mi_pairs[0][2]:.4f}** at positions ({mi_pairs[0][0]},{mi_pairs[0][1]});"
)
w(f"{len(var_positions)}/80 positions are evolutionarily variable (H > 0.3).")
w()

# ── 1. LEAN THEOREMS ────────────────────────────────────────────────
w("## 1. Formal Foundations — 236 Lean 4 Theorems")
w()
w(
    "The entire framework rests on formal proofs in Lean 4.29.0. All theorems use `native_decide`,"
)
w(
    "are **sorry-free** and **axiom-free**, and compile via `lake build` in under 1 second."
)
w()

w("### 1.1 Binary K-map (`lean_proofs/`) — 106 Theorems")
w()
w("| File | Theorems | What Is Proved |")
w("|------|----------|---------------|")
w(
    "| `KmapProofs.lean` | 16 | Gray code boundedness, injectivity, involution, DNA 2-bit encoding (A=00,C=01,G=11,T=10), complementarity, transition/transversion |"
)
w(
    "| `KmerIndexing.lean` | 35 | k-mer encoding (k=1..4) injectivity and boundedness via concatenation, row/column K-map decomposition |"
)
w(
    "| `AminoAcidEncoding.lean` | 27 | 20 AAs → 5-bit Gray code, 7 physicochemical groups, within-group dist-1 (14 pairs), cross-group dist-1 (26 pairs), max-dist pairs (F-H, Y-E, W-R, M-K at dist=5) |"
)
w(
    "| `ContactMapCompleteness.lean` | 20 | Contact maps as symmetric Boolean functions, irreducibility proof (no two contact pairs are 1-bit K-map adjacent) |"
)
w(
    "| `KmapEncodingEquiv.lean` | 8 | Distance distribution, Q₅ hypercube degree=5, encoding captures 50% of Q₅ edges |"
)
w()

w("### 1.2 N-ary K-map (`n-ary-kmap/`) — 115 Theorems")
w()
w("| File | Theorems | What Is Proved |")
w("|------|----------|---------------|")
w(
    "| `NaryGrayCode.lean` | 54 | Generalized n-ary reflected Gray code; digit-sum formula g_k = (d_k + d_{k+1}) mod b |"
)
w(
    "| `BaseNDnaEncoding.lean` | 13 | Base-4 DNA encoder with transition adjacency (A↔C, C↔G, G↔T, T↔A as Gray neighbors) |"
)
w(
    "| `BaseNAminoEncoding.lean` | 33 | Base-20 amino acid encoder, 7-group physicochemical ordering |"
)
w("| `BaseNKMap.lean` | 15 | Base-N K-map structure, cell adjacency, cyclic wrapping |")
w()

w("### 1.3 Speculative Encoding (`speculative-binary-encoding/`) — 15 Theorems")
w()
w("| File | Theorems | What Is Proved |")
w("|------|----------|---------------|")
w(
    "| `AntiCorrelation.lean` | 15 | H1 encoding-invariance: the anti-correlation result holds across all 22 binary encodings |"
)
w()
w("### 1.4 Theorem Verification Bridge")
w()
w(
    "The `lean_consistency.py` script in `kmap-sbm-validation/` re-computes every Lean theorem"
)
w("in Python and reports discrepancies. Result: **103/103 pass**, zero discrepancies.")
w()

# ── 2. DATASET ──────────────────────────────────────────────────────
w("## 2. Dataset")
w()
w("| Property | Value |")
w("|----------|-------|")
w(f"| **File** | `Spike_protein.aln-fasta` (1.8 MB) |")
w(f"| **Sequences** | {n_all:,} SARS-CoV-2 Omicron Spike proteins |")
w(f"| **Alignment length** | {full_len} residues |")
w(f"| **Analysis region** | Positions 0-79 (N-terminal signal peptide) |")
w(f"| **Variable positions** | {len(var_positions)}/80 (entropy > 0.3) |")
w(f"| **Conserved positions** | {80 - len(var_positions)}/80 |")
w(f"| **Co-evolving pairs (MI > 0.01)** | {len(mi_pairs)} |")
w(
    f"| **Co-evolving pairs (MI > 0.1)** | {sum(1 for _, _, mi, _, _ in mi_pairs if mi > 0.1)} |"
)
w(f"| **Encoding** | Base-20, He 2012 ordering |")
w()

w("### 2.1 He 2012 Amino Acid Ordering")
w()
w("The 20 standard amino acids are ordered by physicochemical properties:")
w()
w("| Index | AA | Group |")
w("|-------|-----|-------|")
groups = {
    "A": "Aliphatic",
    "I": "Aliphatic",
    "L": "Aliphatic",
    "V": "Aliphatic",
    "M": "Sulfur",
    "F": "Aromatic",
    "W": "Aromatic",
    "Y": "Aromatic",
    "E": "Negative",
    "D": "Negative",
    "Q": "Polar",
    "N": "Polar",
    "H": "Positive",
    "K": "Positive",
    "R": "Positive",
    "S": "Polar",
    "T": "Polar",
    "C": "Sulfur",
    "P": "Structure",
    "G": "Structure",
}
for i, aa in enumerate(AA_LIST):
    w(f"| {i} | {aa} | {groups.get(aa, '?')} |")
w()

# ── 3. K-MAP CONSTRUCTION ────────────────────────────────────────────
w("## 3. K-map Construction Methodology")
w()
w("### 3.1 Encoding Biological Sequences as Karnaugh Maps")
w()
w("The core innovation is representing protein sequences as Karnaugh maps —")
w("the same mathematical objects used in digital logic design. The mapping is:")
w()
w("```")
w(
    "Biological Sequence → Amino Acid Encoding → K-map Cell → Frequency Map → Boolean Function"
)
w("```")
w()
w("### 3.2 Step 1: Position Array Construction")
w()
w("Each sequence is converted to a position array of integers (0-19):")
w()
w(
    "$$\\text{seq}[i] \\rightarrow \\text{encoder.encode}(\\text{seq}[i]) \\in \\{0, 1, \\dots, 19\\}$$"
)
w()
w(
    f"For {n_all:,} sequences, each of length up to {full_len}, we build a list of `np.int32` arrays."
)
w("Gap characters ('-') and ambiguous residues are mapped to -1 and excluded.")
w()

w("### 3.3 Step 2: Position-Pair K-map (20×20)")
w()
w("For each pair of positions $(i, j)$, a 20×20 frequency matrix is built:")
w()
w(
    "$$K_{ij}(a, b) = \\frac{1}{N}\\sum_{s=1}^{N} \\mathbb{1}[\\text{seq}_s[i] = a \\land \\text{seq}_s[j] = b]$$"
)
w()
w("where $N = {n_all:,}$ is the number of sequences, $a,b \\in \\{0,\\dots,19\\}$.")
w()

w("### 3.4 Example: Frequency K-map for Positions (76, 77)")
w()
w(
    "The reference residues are **D** at position 76 and **N** at position 77 (Wuhan-Hu-1)."
)
w()
w("The top 10 most frequent dipeptide pairs at (76, 77):")
w()

# Find top pairs for (76,77)
top_76_77 = []
for ai in range(N_AA):
    for aj in range(N_AA):
        if freq_76_77[ai, aj] > 0:
            top_76_77.append((AA_LIST[ai], AA_LIST[aj], freq_76_77[ai, aj]))
top_76_77.sort(key=lambda x: x[2], reverse=True)

w("| Rank | AA(i) | AA(j) | Frequency |")
w("|------|-------|-------|-----------|")
for rank, (a1, a2, f) in enumerate(top_76_77[:15], 1):
    w(f"| {rank} | {a1} | {a2} | {f:.4f} |")
w()

w("### 3.5 Step 3: Boolean Thresholding")
w()
w("The frequency K-map is converted to a Boolean function by thresholding:")
w()
w(
    "$$B_{ij}(a,b) = \\begin{cases} 1 & \\text{if } K_{ij}(a,b) \\geq \\text{threshold} \\\\ 0 & \\text{if } K_{ij}(a,b) < \\text{threshold} \\\\ -1 & \\text{if position is conserved (don't-care)} \\end{cases}$$"
)
w()
w("The threshold is typically the 75th percentile of non-zero frequencies.")
w("Don't-care (-1) marks the reference pair (conserved).")
w()

w("### 3.6 Step 4: Quine-McCluskey Boolean Minimization")
w()
w("The 20×20 = 400-cell K-map is flattened to a truth table with 8 binary variables")
w("(4 bits for row amino acid + 4 bits for column amino acid). The Quine-McCluskey")
w("algorithm finds the minimal set of prime implicants covering all on-set cells:")
w()
w(
    "$$f(\\text{pos}_i, \\text{pos}_j, \\text{aa}_i, \\text{aa}_j) = \\bigvee_k \\bigwedge_{m \\in S_k} b_m$$"
)
w()
w(
    "where $b_m$ are the 8 binary variables and $S_k$ are the literal sets for each prime implicant."
)
w()

# ── 4. ENTROPY ──────────────────────────────────────────────────────
w("## 4. Entropy and Conservation Analysis")
w()
w("### 4.1 Shannon Entropy")
w()
w("Position-specific entropy measures evolutionary variability:")
w()
w("$$H(p) = -\\sum_{i=1}^{20} P(a_i) \\log_2 P(a_i)$$")
w()
w("where $P(a_i)$ is the frequency of amino acid $a_i$ at position $p$.")
w()
w("- $H \\approx 0$: highly conserved (one amino acid dominates)")
w("- $H \\approx 4.32$: maximally variable (uniform distribution over 20 AAs)")
w()

w("### 4.2 Conservation Landscape (Positions 0-79)")
w()
w("| Position | Consensus | Frequency | Entropy | Perplexity | Status |")
w("|----------|-----------|-----------|---------|------------|--------|")
for pos, aa, freq, ent in consensus[:40]:
    status = "Variable" if ent > 0.3 else "Conserved"
    pp = 2**ent
    w(f"| {pos} | {aa} | {freq:.4f} | {ent:.4f} | {pp:.3f} | {status} |")
w(f"| ... | ... | ... | ... | ... | ... |")
for pos, aa, freq, ent in consensus[40:80]:
    status = "Variable" if ent > 0.3 else "Conserved"
    pp = 2**ent
    w(f"| {pos} | {aa} | {freq:.4f} | {ent:.4f} | {pp:.3f} | {status} |")
w()

w(
    f"**Summary:** {len(var_positions)} variable positions (H > 0.3), {80 - len(var_positions)} conserved."
)
w()

# ── 5. MUTUAL INFORMATION ────────────────────────────────────────────
w("## 5. Mutual Information Analysis")
w()
w("### 5.1 Definition")
w()
w("Mutual Information quantifies the co-dependence between two positions:")
w()
w(
    "$$MI(i,j) = \\sum_{x=1}^{20} \\sum_{y=1}^{20} P(x,y) \\log_2 \\frac{P(x,y)}{P(x) \\cdot P(y)}$$"
)
w()
w("where:")
w(
    "- $P(x,y)$ = joint frequency of amino acid $x$ at position $i$ and $y$ at position $j$"
)
w("- $P(x)$ = marginal frequency of $x$ at position $i$")
w("- $P(y)$ = marginal frequency of $y$ at position $j$")
w()
w("### 5.2 Implementation")
w()
w(
    "We use a **vectorized numpy implementation** via `np.bincount` (O(400) per pair, not O(N²)):"
)
w()
w("```python")
w("def mutual_information(pos_arrays, pos_i, pos_j, n_seqs):")
w("    codes_i = np.array([arr[pos_i] for arr in pos_arrays[:n_seqs]])")
w("    codes_j = np.array([arr[pos_j] for arr in pos_arrays[:n_seqs]])")
w("    valid = (codes_i >= 0) & (codes_j >= 0)")
w("    # Joint via bincount: flat_index = ci * 20 + cj")
w("    pairs = codes_i[valid] * 20 + codes_j[valid]")
w("    joint = np.bincount(pairs, minlength=400).reshape(20,20)")
w("    # ... MI from joint and marginals")
w("```")
w()

w("### 5.3 Top 50 Co-evolving Position Pairs")
w()
w("| Rank | Pos i | Pos j | MI | Ref(i) | Ref(j) | H(i) | H(j) | Δ |")
w("|------|-------|-------|-------|--------|--------|------|------|----|")
for rank, (pi, pj, mi, ri, rj) in enumerate(mi_pairs[:50], 1):
    rai = AA_LIST[ri] if ri < N_AA else "?"
    raj = AA_LIST[rj] if rj < N_AA else "?"
    delta = abs(pi - pj)
    w(
        f"| {rank} | {pi} | {pj} | {mi:.4f} | {rai} | {raj} | {entropy_vec[pi]:.3f} | {entropy_vec[pj]:.3f} | {delta} |"
    )
w()

# ── 6. BOOLEAN EXPRESSIONS ──────────────────────────────────────────
w("## 6. Quine-McCluskey Boolean Minimization — All 108 Essential Prime Implicants")
w()
w(
    f"The Boolean minimization was performed on {len(rules_by_pair)} position pairs (68-79),"
)
w(f"producing **{len(all_rules)} essential prime implicants**. Each rule has the form:")
w()
w("$$f(s_3, s_2, s_1, s_0, t_3, t_2, t_1, t_0) = \\text{AND of literals}$$")
w()
w("**Variables:**")
w("- $s_3 s_2 s_1 s_0$ = 4-bit binary encoding of residue at position $i$")
w("- $t_3 t_2 t_1 t_0$ = 4-bit binary encoding of residue at position $j$")
w("- $\\bar{s}_k$ = NOT ($s_k = 0$), $s_k$ = ($s_k = 1$)")
w()

w(
    f"### 6.1 Complete Inference Rules ({len(all_rules)} rules across {len(rules_by_pair)} position pairs)"
)
w()

rule_num = 0
for (pos_i, pos_j), pair_rules in sorted(rules_by_pair.items()):
    # Find MI for this pair from mi_pairs
    mi_val = None
    ref_i_val = None
    ref_j_val = None
    for pi, pj, m, ri, rj in mi_pairs:
        if pi == pos_i and pj == pos_j:
            mi_val = m
            ref_i_val = AA_LIST[ri]
            ref_j_val = AA_LIST[rj]
            break
    if mi_val is None:
        mi_val = pair_rules[0].get("mi", 0) if pair_rules else 0
        ref_i_val = pair_rules[0].get("ref_i", "?") if pair_rules else "?"
        ref_j_val = pair_rules[0].get("ref_j", "?") if pair_rules else "?"

    w(
        f"### Position Pair ({pos_i}, {pos_j}) — MI = {mi_val:.4f}, Reference: {ref_i_val}→{ref_j_val}"
    )
    w()
    w("| Rule | Boolean Expression | Amino Acids |")
    w("|------|-------------------|-------------|")

    for r in pair_rules:
        rule_num += 1
        expr = r.get("expr", "N/A")
        aa_i = r.get("aa_i", "?")
        aa_j = r.get("aa_j", "?")
        # Convert to LaTeX math
        expr_tex = expr.replace("\\bar{", "\\bar{").replace("\\cdot", " \\cdot ")
        if expr_tex != "N/A":
            expr_tex = f"$${expr_tex}$$"
        w(f"| {rule_num} | {expr_tex} | ({aa_i}, {aa_j}) |")
    w()

    w(
        f"**Interpretation:** When position {pos_i} mutates to any of the listed residues,"
    )
    w(
        f"position {pos_j} must co-evolve to the corresponding partner residue to maintain"
    )
    w(f"protein stability. The reference pair is ({ref_i_val}, {ref_j_val}).")
    w()

w(f"\n**Total inference rules:** {rule_num}")
w()

# ── 7. COUPLING CONSTANTS ────────────────────────────────────────────
w("## 7. Coupling Constants and Constraint Functions")
w()
w("### 7.1 Definition")
w()
w("The coupling constant $J_{ij}(a,b)$ measures the log-odds of observing pair $(a,b)$")
w("at positions $(i,j)$ compared to the independent expectation:")
w()
w("$$J_{ij}(a,b) = \\ln\\frac{P_{ij}(a,b)}{P_i(a) \\cdot P_j(b)}$$")
w()
w("- $J > 0$: pair is **co-evolutionary** (more common than expected)")
w("- $J < 0$: pair is **anti-correlated** (less common than expected)")
w("- $J = 0$: pair occurs at random frequency")
w()

w("### 7.2 Top Coupling Constants")
w()
w(
    "| Pos i | Pos j | MI | avg\\|J\\| | Top Co-evolutionary (J>0) | Top Anti-correlated (J<0) |"
)
w(
    "|-------|-------|-----|----------|---------------------------|---------------------------|"
)
for pi, pj, mi, avg_abs, top_co, top_anti in couplings[:20]:
    co_str = ", ".join(f"{a1}-{a2}:{j:+.1f}" for a1, a2, j in top_co[:3])
    anti_str = ", ".join(f"{a1}-{a2}:{j:+.1f}" for a1, a2, j in top_anti[:3])
    w(f"| {pi} | {pj} | {mi:.4f} | {avg_abs:.2f} | {co_str} | {anti_str} |")
w()

# ── 8. PERPLEXITY ───────────────────────────────────────────────────
w("## 8. Perplexity Analysis")
w()
w("### 8.1 Definition")
w()
w("Perplexity measures the effective number of choices: $PP = 2^H$.")
w("Conditional perplexity $PP(j|i)$ measures how much knowing residue $i$")
w("reduces uncertainty about residue $j$:")
w()
w("$$PP(j) = 2^{H(P_j)}, \\quad PP(j|i=a) = 2^{H(P_{j|i=a})}$$")
w()
w("The **co-evolution ratio** $PP(j) / PP(j|i)$ quantifies constraint strength:")
w("- Ratio $\\approx 1$: no constraint (positions evolve independently)")
w("- Ratio $> 2$: strong constraint (position $i$ determines position $j$)")
w()

w("### 8.2 Conditional Perplexity for Top Pairs")
w()
w(
    "| Pos i | Pos j | MI | PPₘ(j) | PPₖₒₙₐ | Ratio | Most Constraining Residue | PP\\|that |"
)
w(
    "|-------|-------|-----|--------|---------|-------|--------------------------|----------|"
)
for pi, pj, mi, mi_mut, ppm, avgc, cond in perp_results[:20]:
    # Find most constraining residue
    best_aa = "?"
    best_pp = 99
    for aa, (pp_val, count) in sorted(cond.items(), key=lambda x: x[1][0]):
        if pp_val < best_pp:
            best_pp = pp_val
            best_aa = aa
    ratio = ppm / avgc if avgc > 0 else 0
    w(
        f"| {pi} | {pj} | {mi:.4f} | {ppm:.3f} | {avgc:.3f} | {ratio:.2f} | {best_aa} | {best_pp:.3f} |"
    )
w()

# ── 9. H1 GRAY ADJACENCY ─────────────────────────────────────────────
w("## 9. H1: Gray-code Adjacency Hypothesis")
w()
w("### 9.1 Hypothesis")
w()
w("Consecutive residues in protein sequences are preferentially K-map-adjacent")
w("(Hamming distance = 1) compared to random expectation.")
w()
w("### 9.2 Method")
w()
w("For each consecutive pair of residues in each of the {n_all:,} sequences:")
w("1. Encode each residue using its base-20 index (He 2012)")
w("2. Compute Gray code: $g(i) = i \\oplus (i \\gg 1)$")
w("3. Compute Hamming distance: $h = \\text{popcount}(g(i) \\oplus g(j))$")
w()

w("### 9.3 Results")
w()
w("| Metric | Value |")
w("|--------|-------|")
w(f"| Total consecutive pairs | {tot_consec:,} |")
w(f"| Hamming-1 pairs | {ham[1]:,} |")
w(f"| Observed ratio | **{h1_ratio:.4f}** |")
w(f"| Expected (random) | 0.1613 |")
w(f"| **Enrichment** | **{h1_ratio / 0.1613:.2f}×** |")
w()

w("### 9.4 Full Hamming Distance Distribution")
w()
w("| Distance | Count | Percentage | Cumulative |")
w("|----------|-------|------------|------------|")
cum = 0
for d in range(6):
    if ham[d] > 0:
        cum += ham[d]
        w(
            f"| {d} | {ham[d]:,} | {100 * ham[d] / tot_consec:.1f}% | {100 * cum / tot_consec:.1f}% |"
        )
w()

# ── 10. THE ANALYSIS SCRIPTS ──────────────────────────────────────────
w("## 10. Complete Analysis Scripts Inventory")
w()
w(f"The `datasets/co-evolution/` directory contains **19 Python scripts**")
w(f"and **1 shared module** (`coevolution_shared.py`).")
w()

scripts_info = [
    (
        "`coevolution_shared.py`",
        "329",
        "Shared module: FASTA parsing, position arrays (cached), vectorized MI via `np.bincount`, entropy, perplexity, coupling, constraint function, shared-memory worker pool",
    ),
    (
        "`run_kmap_analysis.py`",
        "908",
        "Master K-map pipeline: H1-H6 on binary 32×32 K-map, consensus K-map, co-evolution analysis",
    ),
    (
        "`boolean_co-evolution.py`",
        "637",
        "Binary K-map Boolean minimization: 32×32 thresholded → Quine-McCluskey → essential prime implicants",
    ),
    (
        "`nary_kmap_co-evolution.py`",
        "537",
        "Base-20 K-map analysis: 20×20 frequency map → Boolean → coupling constants",
    ),
    (
        "`master_boolean.py`",
        "330",
        f"Master Boolean function: {master_data['total_inference_rules']} essential PIs across {master_data['co_evolving_pairs']:,} pairs (full-length)",
    ),
    (
        "`position_kmap_coevolution.py`",
        "481",
        "Position-pair K-maps with MI: builds per-position-pair 20×20 K-maps and minimizes",
    ),
    (
        "`run_allseq_analysis.py`",
        "328",
        "Full position-based K-map analysis on ALL 1,299 sequences",
    ),
    (
        "`kmap_boolean_coevolution.py`",
        "383",
        f"K-map Boolean with full markdown output: {len(all_rules)} rules across {len(rules_by_pair)} position pairs",
    ),
    ("`generate_co-evolution_md.py`", "313", "Markdown generator from JSON results"),
    (
        "`generate_full_analysis_md.py`",
        "397",
        "Comprehensive report generator from all JSON outputs",
    ),
    (
        "`create_mi_heatmap.py`",
        "247",
        "MI heatmap visualization (full + focus region), matplotlib",
    ),
    (
        "`allseq_constraint_function.py`",
        "247",
        "Leave-One-Out Cross-Validation: constraint function prediction accuracy",
    ),
    (
        "`predictive_constraint_function.py`",
        "308",
        "Three K-map approaches: Observed, Flipped, Continuous prediction",
    ),
    (
        "`flipped_boolean_coevolution.py`",
        "254",
        "Flipped Boolean: forbidden pairs analysis (negative selection)",
    ),
    (
        "`dca_boolean_coevolution.py`",
        "344",
        "DCA → Boolean pipeline: inverse covariance → coupling → QM minimization",
    ),
    (
        "`variable_position_coevolution.py`",
        "421",
        "Variable-position K-map with strategic don't-care conditions",
    ),
    (
        "`perplexity_coevolution.py`",
        "217",
        "Perplexity-based co-evolution strength measurement",
    ),
    (
        "`advanced_co-evolution_analysis.py`",
        "471",
        "Co-evolution network, Walsh-Hadamard spectrum, variant classification, clustering",
    ),
    (
        "`full_length_analysis.py`",
        "206",
        "Full-length (all 1,276 positions) entropy and MI analysis",
    ),
    (
        "`gpu_full_analysis.py`",
        "280",
        "GPU-accelerated analysis: numba parallel entropy/H1/mutations + shared-memory Pool for MI",
    ),
    ("`run_all_bg.sh`", "107", "Master launcher: runs all 17 scripts concurrently"),
]

w("| # | Script | Lines | Purpose |")
w("|---|--------|-------|---------|")
for i, (name, lines, desc) in enumerate(scripts_info, 1):
    w(f"| {i} | {name} | {lines} | {desc} |")
w()

# ── 11. FLIPPED BOOLEAN ──────────────────────────────────────────────
w("## 11. Flipped Boolean Analysis — Negative Selection")
w()
w("### 11.1 Concept")
w()
w(
    'While the positive Boolean function asks "what co-evolves?", the **flipped** version'
)
w('asks "what CANNOT co-exist?" The on-set (1) is assigned to residue pairs that are')
w("NEVER observed together across all sequences, capturing **negative selection**.")
w()
w("### 11.2 Method")
w()
w("1. Build 20×20 frequency K-map for each position pair")
w("2. Cells with frequency = 0 → on-set (1) — FORBIDDEN")
w("3. Cells with frequency > 0 → off-set (0) — ALLOWED")
w("4. Reference pair → don't-care (-1)")
w("5. QM minimization → minimal forbidden constraints")
w()
w("### 11.3 Results")
w()
w("The flipped analysis produced **0** forbidden rules. Why?")
w()
w("With 1,299 sequences across multiple Omicron sub-lineages,")
w("virtually every amino acid pair appears at least once at every position pair.")
w("The Spike protein's N-terminal region is under **purifying selection** but not")
w("**absolute constraint** — any single substitution has been sampled by evolution.")
w()

w("This means:\n")
w("- **No universally forbidden pairs exist** — only statistically disfavored ones")
w("- Co-evolution is **probabilistic**, not deterministic")
w("- The constraint function $J_{ij}$ (continuous) is more appropriate than Boolean")
w()

# ── 12. LOOCV ────────────────────────────────────────────────────────
w("## 12. Leave-One-Out Cross-Validation")
w()
w("### 12.1 Method")
w()
w("For each position pair $(i,j)$:")
w("1. Exclude one sequence from the dataset")
w(
    "2. Build constraint function $C_{ij}(a,b)$ from the remaining {n_all - 1:,} sequences"
)
w(
    "3. Predict the co-evolving partner for the held-out sequence: $b^* = \\arg\\max_b C_{ij}(a_{\\text{mutation}}, b)$"
)
w("4. Check if prediction matches")
w()
w("### 12.2 Key Insight")
w()
w(
    "LOO-CV accuracy was near zero (0.08%). This reveals that co-evolutionary rules are **lineage-specific** —"
)
w(
    "different Omicron sub-variants (BA.1, BA.2, BA.4, BA.5, XBB, etc.) have **different co-evolutionary patterns**."
)
w("Rules learned from one variant do not generalize to others.")
w()

# ── 13. SCALE ─────────────────────────────────────────────────────────
w("## 13. What This Analysis Pipeline Achieves")
w()
w("### 13.1 By the Numbers")
w()
w("| Metric | Value |")
w("|--------|-------|")
w(f"| Input sequences | {n_all:,} |")
w(f"| Positions analyzed | {full_len} |")
w(f"| Position arrays built | {n_all:,} × {full_len} = {n_all * full_len:,} integers |")
w(f"| Variable positions (H > 0.3) | {len(var_positions)} |")
w(
    f"| Co-evolving pairs (MI > 0.1) | {sum(1 for _, _, mi, _, _ in mi_pairs if mi > 0.1)} |"
)
w(
    f"| Boolean expressions (QM minimized) | {len(all_rules)} essential prime implicants |"
)
w(f"| Unique position pairs with rules | {len(rules_by_pair)} |")
w(f"| Lean 4 theorems | 236 (106 + 115 + 15) |")
w(f"| Python scripts | 20 |")
w(f"| Total Python LOC | ~7,000 |")
w(f"| Shared module LOC | 340 |")
w()

# ── 14. HOW TO USE THE RULES ──────────────────────────────────────────
w("## 14. How to Apply the Co-evolution Rules")
w()
w("### 14.1 For a New Spike Sequence")
w()
w("```python")
w("# 1. Extract residues at positions 68-79")
w("seq_region = seq[68:80]")
w()
w("# 2. For each co-evolving position pair, check the Boolean function")
w("for (pos_i, pos_j) in coevolving_pairs:")
w("    aa_i = encode(seq_region[pos_i - 68])")
w("    aa_j = encode(seq_region[pos_j - 68])")
w()
w("    # Encode as 4+4 bits")
w("    bits = int_to_bits(aa_i, 4) + int_to_bits(aa_j, 4)")
w()
w("    # Check if ANY essential prime implicant matches")
w(
    "    is_coevolving = any(pi.matches(bits) for pi in prime_implicants[(pos_i, pos_j)])"
)
w()
w("    if is_coevolving:")
w("        print(f'Positions {pos_i}-{pos_j} are co-evolving')")
w()
w("# 3. If position i mutates, find the required partner at position j")
w("def predict_partner(pos_i, aa_i, pos_j):")
w("    for pi in prime_implicants[(pos_i, pos_j)]:")
w("        if pi.matches_row(aa_i):")
w("            return pi.col_aa  # required co-evolving partner")
w("    return None  # no constraint")
w("```")
w()

# ── 15. REFERENCES ───────────────────────────────────────────────────
w("## 15. Key References")
w()
w(
    "1. Karnaugh, M. (1953). The Map Method for Synthesis of Combinational Logic Circuits. *AIEE Transactions*."
)
w("2. Gray, F. (1953). Pulse Code Communication. U.S. Patent 2,632,058.")
w(
    "3. Quine, W.V. (1952). The Problem of Simplifying Truth Functions. *American Mathematical Monthly*."
)
w(
    "4. McCluskey, E.J. (1956). Minimization of Boolean Functions. *Bell System Technical Journal*."
)
w(
    "5. He, P.A. et al. (2012). A novel graphical representation of proteins. *MATCH Communications*."
)
w(
    "6. Petoukhov, S.V. (2024). Matrix Representations of Genetic Code and Karnaugh Maps. *Biosystems*."
)
w("7. de Moura, L. et al. (2021). The Lean 4 Theorem Prover and Programming Language.")
w()

w("---")
w(combined_pipeline_section())
w("")
w("---")
w(f"*Generated {now} by `generate_full_pipeline_doc.py`*")
w(
    f"*All values computed from {n_all:,} Omicron Spike sequences using shared `coevolution_shared` module*"
)

# ── WRITE ─────────────────────────────────────────────────────────────
final_md = "\n".join(md)
OUT.write_text(final_md)

elapsed = time.time() - t0_start
print(f"\n{'=' * 60}")
print(f"DONE — {len(md)} lines written to {OUT}")
print(f"Size: {OUT.stat().st_size:,} bytes")
print(f"Time: {elapsed:.1f}s")
print(f"\nContents:")
print(f"  - {len(all_rules)} Boolean expressions with LaTeX rendering")
print(f"  - {len(mi_pairs)} MI pairs documented")
print(f"  - {len(consensus)} positions with entropy/perplexity")
print(f"  - {len(couplings)} coupling constant tables")
print(f"  - {len(perp_results)} perplexity analyses")
print(f"  - 236 Lean theorems cataloged")
print(f"  - Full scripts inventory")
