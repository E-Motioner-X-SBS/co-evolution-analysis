#!/usr/bin/env python3
"""
K-map Boolean Co-evolution Analysis
=====================================

Creates proper K-maps for co-evolutionary position pairs,
derives Boolean functions via Quine-McCluskey, and outputs
all inference rules with full details.
"""

import sys
import json
import numpy as np
from collections import Counter
from pathlib import Path

sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/kmap-sbm-validation/src")

from nkmap.encoding.bio_sequences import Base20AminoEncoder, AMINO_HE_2012
from kmap_sbm.analysis.prime_implicants import boolean_minimize_kmap


def parse_fasta(filepath):
    sequences = []
    current_header = None
    current_seq = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences.append((current_header, "".join(current_seq)))
                current_header = line[1:]
                current_seq = []
            elif line:
                current_seq.append(line.upper())
    if current_header is not None:
        sequences.append((current_header, "".join(current_seq)))
    return sequences


def compute_entropy(pos_arrays, pos, n_seqs):
    counts = Counter()
    for arr in pos_arrays[:n_seqs]:
        if pos < len(arr) and arr[pos] >= 0:
            counts[int(arr[pos])] += 1
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((c / total) * np.log2(c / total) for c in counts.values() if c > 0)


def get_majority_ref(pos_arrays, pos, n_seqs):
    return Counter(
        int(a[pos]) for a in pos_arrays[:n_seqs] if pos < len(a) and a[pos] >= 0
    ).most_common(1)[0][0]


def compute_mi(pos_arrays, pos_i, pos_j, n_seqs):
    ref_i = get_majority_ref(pos_arrays, pos_i, n_seqs)
    ref_j = get_majority_ref(pos_arrays, pos_j, n_seqs)
    joint, marg_i, marg_j = Counter(), Counter(), Counter()
    for arr in pos_arrays[:n_seqs]:
        if pos_i < len(arr) and pos_j < len(arr):
            ci, cj = int(arr[pos_i]), int(arr[pos_j])
            if ci >= 0 and cj >= 0 and (ci != ref_i or cj != ref_j):
                joint[(ci, cj)] += 1
                marg_i[ci] += 1
                marg_j[cj] += 1
    total = sum(joint.values())
    if total < 5:
        return 0.0, 0
    mi = sum(
        (c / total)
        * np.log2((c / total) / ((marg_i[ai] / total) * (marg_j[aj] / total)))
        for (ai, aj), c in joint.items()
        if marg_i[ai] > 0 and marg_j[aj] > 0
    )
    return mi, total


def build_mutation_kmap(pos_arrays, pos_i, pos_j, ref_i, ref_j, n_seqs):
    kmap = np.zeros((20, 20), dtype=np.int32)
    for arr in pos_arrays[:n_seqs]:
        if pos_i < len(arr) and pos_j < len(arr):
            ci, cj = int(arr[pos_i]), int(arr[pos_j])
            if ci >= 0 and cj >= 0:
                kmap[ci, cj] = 1 if (ci != ref_i or cj != ref_j) else -1
    return kmap


def kmap_to_markdown(kmap, aa_list, pos_i, pos_j, ref_aa_i, ref_aa_j):
    """Convert a 20x20 K-map to a markdown table."""
    lines = []
    lines.append(f"### K-map for Position Pair ({pos_i}, {pos_j})")
    lines.append(f"Reference: pos {pos_i}={ref_aa_i}, pos {pos_j}={ref_aa_j}")
    lines.append("")
    lines.append("| | " + " | ".join(aa_list) + " |")
    lines.append("|---" + "|---" * 20 + "|")
    for i in range(20):
        row = [aa_list[i]]
        for j in range(20):
            v = kmap[i, j]
            if v == 1:
                row.append("**1**")
            elif v == -1:
                row.append("DC")
            else:
                row.append("0")
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def decode_pi(pi, aa_list):
    """Decode a prime implicant to amino acid pair and Boolean expression."""
    values = list(pi["values"])
    mask = list(pi["mask"])
    while len(values) < 8:
        values.append(0)
        mask.append(False)

    row_code = sum(values[j] * (2 ** (3 - j)) for j in range(4) if not mask[j])
    col_code = sum(values[j + 4] * (2 ** (3 - j)) for j in range(4) if not mask[j + 4])

    row_aa = aa_list[row_code % 20] if row_code < 20 else "?"
    col_aa = aa_list[col_code % 20] if col_code < 20 else "?"

    terms = []
    for j in range(8):
        if not mask[j]:
            var = f"s{3 - j}" if j < 4 else f"t{7 - j}"
            if values[j] == 0:
                terms.append(f"\\bar{{{var}}}")
            else:
                terms.append(var)

    expr = " \\cdot ".join(terms) if terms else "\\text{TRUE}"
    return row_aa, col_aa, expr, pi["n_dontcares"], len(pi["coverage"])


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    results_dir = base_dir / "kmap_boolean_coevolution"
    results_dir.mkdir(exist_ok=True)
    md_path = results_dir / "COEVOLUTION_KMAP_BOOLEAN.md"

    print("=" * 80)
    print("K-map Boolean Co-evolution Analysis")
    print("=" * 80)

    # Load ALL sequences
    print("\n[1/4] Loading ALL sequences...")
    sequences = parse_fasta(fasta_file)
    encoder = Base20AminoEncoder(version=1)
    aa_list = list(AMINO_HE_2012)
    n_all = len(sequences)
    full_length = len(sequences[0][1])
    print(f"  Total: {n_all} sequences")

    # Build position arrays
    print("\n[2/4] Building position arrays...")
    max_pos = len(sequences[0][1])
    pos_arrays = []
    for _, seq in sequences:
        clean = "".join(aa for aa in seq if aa in encoder.encode)
        arr = np.array(
            [encoder.encode.get(aa, -1) for aa in clean[:max_pos]], dtype=np.int32
        )
        pos_arrays.append(arr)

    # Find variable positions
    variable_positions = [
        p for p in range(max_pos) if compute_entropy(pos_arrays, p, n_all) > 0.3
    ]
    print(f"  Variable positions: {len(variable_positions)}")

    # Find co-evolutionary pairs (GPU-accelerated via coevolution_shared)
    print("\n[3/4] Finding co-evolutionary pairs (GPU)...")
    from coevolution_shared import find_coevolving_pairs_gpu

    co_evolving = find_coevolving_pairs_gpu(
        pos_arrays, variable_positions, n_all, max_gap=30, min_mi=0.1
    )
    print(f"  Co-evolutionary pairs: {len(co_evolving)}")

    # Build K-maps and Boolean functions
    print("\n[4/4] Building K-maps and Boolean functions...")
    md_lines = []
    md_lines.append("# Co-evolution Boolean Functions for SARS-CoV-2 Spike Protein")
    md_lines.append("")
    md_lines.append(f"**Dataset:** {n_all} sequences")
    md_lines.append(f"**Encoding:** Base-20 (He 2012)")
    md_lines.append(f"**Variable positions:** {len(variable_positions)} / {max_pos}")
    md_lines.append(f"**Co-evolutionary pairs:** {len(co_evolving)}")
    md_lines.append("")

    all_rules = []
    rule_num = 0

    for pos_i, pos_j, mi, n_muts, ref_i, ref_j in co_evolving[:15]:
        ref_aa_i = aa_list[ref_i]
        ref_aa_j = aa_list[ref_j]

        # Build K-map
        kmap = build_mutation_kmap(pos_arrays, pos_i, pos_j, ref_i, ref_j, n_all)

        # Run QM
        bool_flat = kmap.flatten().astype(int)
        result = boolean_minimize_kmap(bool_flat, algorithm="qm")

        # Decode PIs
        pis = []
        for pi in result["prime_implicants"]:
            aa_i, aa_j, expr, n_dc, coverage = decode_pi(pi, aa_list)
            is_essential = pi in result["essential_prime_implicants"]
            pis.append(
                {
                    "aa_i": aa_i,
                    "aa_j": aa_j,
                    "expr": expr,
                    "n_dc": n_dc,
                    "coverage": coverage,
                    "is_essential": is_essential,
                }
            )

        # Write to markdown
        md_lines.append(f"---")
        md_lines.append(f"## Position Pair ({pos_i}, {pos_j})")
        md_lines.append(f"")
        md_lines.append(f"| Property | Value |")
        md_lines.append(f"|----------|-------|")
        md_lines.append(f"| Mutual Information | {mi:.4f} |")
        md_lines.append(f"| Total mutations | {n_muts} |")
        md_lines.append(f"| Reference pos {pos_i} | {ref_aa_i} |")
        md_lines.append(f"| Reference pos {pos_j} | {ref_aa_j} |")
        md_lines.append(f"| On-set cells | {int((kmap == 1).sum())} |")
        md_lines.append(f"| Off-set cells | {int((kmap == 0).sum())} |")
        md_lines.append(f"| Don't-care cells | {int((kmap == -1).sum())} |")
        md_lines.append(f"| Prime implicants | {result['n_prime_implicants']} |")
        md_lines.append(f"| Essential PIs | {result['n_essential']} |")
        md_lines.append(f"")

        # K-map table
        md_lines.append(
            f"### K-map (1 = co-evolutionary, 0 = never seen, DC = don't-care)"
        )
        md_lines.append(f"")
        header = "| AA_i \\ AA_j | " + " | ".join(aa_list) + " |"
        separator = "|---" + "|---" * 20 + "|"
        md_lines.append(header)
        md_lines.append(separator)
        for i in range(20):
            row = [aa_list[i]]
            for j in range(20):
                v = kmap[i, j]
                if v == 1:
                    row.append("**1**")
                elif v == -1:
                    row.append("DC")
                else:
                    row.append("0")
            md_lines.append("| " + " | ".join(row) + " |")
        md_lines.append("")

        # Boolean function
        md_lines.append(f"### Boolean Function")
        md_lines.append(f"")
        md_lines.append(f"$$f(s_3,s_2,s_1,s_0,\\; t_3,t_2,t_1,t_0) = ")
        md_lines.append(f"\\text{{OR of all prime implicants below}}$$")
        md_lines.append(f"")

        # Essential PIs as equations
        md_lines.append(f"### Essential Prime Implicants (minimum covering set)")
        md_lines.append(f"")
        md_lines.append(
            f"| # | Boolean Expression | Amino Acids | Don't-cares | Coverage |"
        )
        md_lines.append(
            f"|---|-------------------|-------------|-------------|----------|"
        )
        pi_num = 0
        for p in pis:
            if p["is_essential"]:
                pi_num += 1
                md_lines.append(
                    f"| {pi_num} | ${p['expr']}$ | ({p['aa_i']}, {p['aa_j']}) | {p['n_dc']} | {p['coverage']} |"
                )
                all_rules.append(
                    {
                        "pos_i": pos_i,
                        "pos_j": pos_j,
                        "aa_i": p["aa_i"],
                        "aa_j": p["aa_j"],
                        "ref_i": ref_aa_i,
                        "ref_j": ref_aa_j,
                        "mi": mi,
                        "expr": p["expr"],
                        "n_dc": p["n_dc"],
                        "coverage": p["coverage"],
                    }
                )
        md_lines.append("")

        # Natural language rules
        md_lines.append(f"### Inference Rules (Natural Language)")
        md_lines.append(f"")
        for p in pis:
            if p["is_essential"]:
                rule_num += 1
                md_lines.append(
                    f"**Rule {rule_num}:** IF position {pos_i} = **{p['aa_i']}** "
                    f"AND position {pos_j} = **{p['aa_j']}** "
                    f"THEN co-evolutionary (MI = {mi:.3f})"
                )
        md_lines.append("")

        print(
            f"  ({pos_i:2d},{pos_j:2d}): MI={mi:.4f}, "
            f"PI={result['n_prime_implicants']}, EPI={result['n_essential']}"
        )

    # Summary
    md_lines.append("---")
    md_lines.append("")
    md_lines.append("## Summary")
    md_lines.append("")
    md_lines.append(f"| Metric | Value |")
    md_lines.append(f"|--------|-------|")
    md_lines.append(f"| Sequences | {n_all} |")
    md_lines.append(f"| Variable positions | {len(variable_positions)} |")
    md_lines.append(f"| Co-evolutionary pairs | {len(co_evolving)} |")
    md_lines.append(f"| Total inference rules | {len(all_rules)} |")
    md_lines.append(
        f"| Position pairs with rules | {len(set((r['pos_i'], r['pos_j']) for r in all_rules))} |"
    )
    md_lines.append("")

    md_lines.append("## How to Apply")
    md_lines.append("")
    md_lines.append(
        "For a new sequence, extract residues at positions 68-79 and check:"
    )
    md_lines.append("```")
    md_lines.append("For each (pos_i, pos_j) pair:")
    md_lines.append("  1. Get residues aa_i at pos_i and aa_j at pos_j")
    md_lines.append("  2. Check if (aa_i, aa_j) matches any essential PI")
    md_lines.append("  3. If YES → position pair is co-evolutionary")
    md_lines.append("  4. If position i mutates → find which aa_j satisfies the PI")
    md_lines.append("```")

    # Write markdown
    with open(md_path, "w") as f:
        f.write("\n".join(md_lines))

    # Save JSON
    summary = {
        "dataset": f"{n_all} sequences",
        "variable_positions": len(variable_positions),
        "co_evolving_pairs": len(co_evolving),
        "total_rules": len(all_rules),
        "rules": all_rules,
    }
    with open(results_dir / "boolean_functions.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nResults saved to: {results_dir}")
    print(f"Markdown: {md_path}")
    print(f"\nTotal inference rules: {len(all_rules)}")


if __name__ == "__main__":
    main()
