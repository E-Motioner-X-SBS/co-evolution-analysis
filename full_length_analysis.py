#!/usr/bin/env python3
"""
Full-Length Co-evolution Analysis
=================================

Analyzes ALL 1,276 positions of the Spike protein,
not just the first 80.
"""

import sys
import json
import numpy as np
from collections import Counter
from pathlib import Path

sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")

from nkmap.encoding.bio_sequences import Base20AminoEncoder, AMINO_HE_2012


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


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    results_dir = base_dir / "full_length_results"
    results_dir.mkdir(exist_ok=True)

    print("=" * 80)
    print("Full-Length Co-evolution Analysis (ALL 1,276 positions)")
    print("=" * 80)

    # Load ALL sequences
    print("\n[1/4] Loading ALL sequences...")
    sequences = parse_fasta(fasta_file)
    encoder = Base20AminoEncoder(version=1)
    aa_list = list(AMINO_HE_2012)
    n_all = len(sequences)
    full_length = len(sequences[0][1])
    print(f"  Total: {n_all} sequences")
    print(f"  Full length: {full_length} residues")

    # Build position arrays for FULL length
    print("\n[2/4] Building position arrays (full length)...")
    pos_arrays = []
    for _, seq in sequences:
        clean = "".join(aa for aa in seq if aa in encoder.encode)
        arr = np.array(
            [encoder.encode.get(aa, -1) for aa in clean[:full_length]], dtype=np.int32
        )
        pos_arrays.append(arr)

    # Compute entropy for ALL positions
    print("\n[3/4] Computing entropy for ALL positions...")
    entropies = []
    for pos in range(full_length):
        ent = compute_entropy(pos_arrays, pos, n_all)
        entropies.append(ent)

    entropies = np.array(entropies)

    # Find variable positions
    variable_threshold = 0.3
    variable_positions = np.where(entropies > variable_threshold)[0]

    print(f"  Total positions: {full_length}")
    print(
        f"  Variable positions (entropy > {variable_threshold}): {len(variable_positions)}"
    )

    # Show top 20 most variable positions
    print(f"\n  Top 20 most variable positions:")
    sorted_indices = np.argsort(entropies)[::-1]
    print(f"  {'Position':>8s} {'Entropy':>8s} {'Perplexity':>12s}")
    print(f"  {'-' * 8} {'-' * 8} {'-' * 12}")
    for idx in sorted_indices[:20]:
        print(f"  {idx:8d} {entropies[idx]:8.3f} {2 ** entropies[idx]:12.3f}")

    # Compute MI for top variable position pairs
    print("\n[4/4] Computing MI for top variable position pairs...")
    print("  (Limiting to top 50 variable positions for efficiency)")

    top_var = sorted_indices[:100]  # Use top 100 most variable positions
    n_var = len(top_var)
    mi_matrix = np.zeros((n_var, n_var), dtype=np.float64)

    for idx_i in range(n_var):
        for idx_j in range(idx_i + 1, n_var):
            pos_i = top_var[idx_i]
            pos_j = top_var[idx_j]

            # Compute MI
            ref_i = Counter(
                int(a[pos_i])
                for a in pos_arrays[:n_all]
                if pos_i < len(a) and a[pos_i] >= 0
            ).most_common(1)[0][0]
            ref_j = Counter(
                int(a[pos_j])
                for a in pos_arrays[:n_all]
                if pos_j < len(a) and a[pos_j] >= 0
            ).most_common(1)[0][0]

            joint, marg_i, marg_j = Counter(), Counter(), Counter()
            for arr in pos_arrays[:n_all]:
                if pos_i < len(arr) and pos_j < len(arr):
                    ci, cj = int(arr[pos_i]), int(arr[pos_j])
                    if ci >= 0 and cj >= 0 and (ci != ref_i or cj != ref_j):
                        joint[(ci, cj)] += 1
                        marg_i[ci] += 1
                        marg_j[cj] += 1

            total = sum(joint.values())
            if total >= 5:
                mi = sum(
                    (c / total)
                    * np.log2(
                        (c / total) / ((marg_i[ai] / total) * (marg_j[aj] / total))
                    )
                    for (ai, aj), c in joint.items()
                    if marg_i[ai] > 0 and marg_j[aj] > 0
                )
                mi_matrix[idx_i, idx_j] = mi
                mi_matrix[idx_j, idx_i] = mi

    # Find top MI pairs
    high_mi = []
    for idx_i in range(n_var):
        for idx_j in range(idx_i + 1, n_var):
            if mi_matrix[idx_i, idx_j] > 0.5:
                high_mi.append(
                    (
                        int(top_var[idx_i]),
                        int(top_var[idx_j]),
                        float(mi_matrix[idx_i, idx_j]),
                    )
                )

    high_mi.sort(key=lambda x: x[2], reverse=True)

    print(f"\n  Top 20 co-evolutionary pairs (full length):")
    print(f"  {'Pos i':>8s} {'Pos j':>8s} {'MI':>8s}")
    print(f"  {'-' * 8} {'-' * 8} {'-' * 8}")
    for pos_i, pos_j, mi in high_mi[:20]:
        print(f"  {pos_i:8d} {pos_j:8d} {mi:8.4f}")

    # Save results
    summary = {
        "n_sequences": n_all,
        "full_length": full_length,
        "n_variable_positions": int(len(variable_positions)),
        "variable_threshold": variable_threshold,
        "top_20_variable_positions": [
            {
                "position": int(i),
                "entropy": float(entropies[i]),
                "perplexity": float(2 ** entropies[i]),
            }
            for i in sorted_indices[:20]
        ],
        "top_20_co_evolving_pairs": [
            {"pos_i": pi, "pos_j": pj, "mi": mi} for pi, pj, mi in high_mi[:20]
        ],
        "n_high_mi_pairs": len(high_mi),
    }

    with open(results_dir / "full_length_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nResults saved to: {results_dir}")
    print(f"\nSummary:")
    print(f"  Full length: {full_length} residues")
    print(f"  Variable positions: {len(variable_positions)}")
    print(f"  High MI pairs: {len(high_mi)}")


if __name__ == "__main__":
    main()
