#!/usr/bin/env python3
"""
Perplexity-Based Co-evolution Analysis
=======================================

Uses perplexity (from information theory) to measure predictability
of co-evolutionary patterns.

Key insight:
- Low perplexity = predictable (constrained by co-evolution)
- High perplexity = unpredictable (many possible outcomes)
- Difference between conditional and marginal perplexity = co-evolution strength
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


def compute_entropy(dist):
    """Compute Shannon entropy from a distribution."""
    total = sum(dist.values())
    if total == 0:
        return 0.0
    return -sum((c / total) * np.log2(c / total) for c in dist.values() if c > 0)


def compute_perplexity(entropy):
    """Compute perplexity from entropy: PP = 2^H."""
    return 2**entropy


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    results_dir = base_dir / "perplexity_results"
    results_dir.mkdir(exist_ok=True)

    print("=" * 80)
    print("Perplexity-Based Co-evolution Analysis")
    print("=" * 80)
    print("Perplexity = 2^H where H is Shannon entropy")
    print("Low perplexity = predictable (constrained)")
    print("High perplexity = unpredictable (many outcomes)")
    print("")

    # Load ALL sequences
    print("[1/4] Loading ALL sequences...")
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
        arr = np.array(
            [encoder.encode.get(aa, 20) for aa in seq[:max_pos]], dtype=np.int32
        )
        pos_arrays.append(arr)

    # Compute position-specific perplexity
    print("\n[3/4] Computing position-specific perplexity...")
    pos_entropy = np.zeros(max_pos)
    pos_perplexity = np.zeros(max_pos)
    pos_counts = np.zeros(max_pos, dtype=np.int64)

    for pos in range(max_pos):
        counts = Counter()
        for arr in pos_arrays[:n_all]:
            if pos < len(arr) and 0 <= arr[pos] < 20:
                counts[int(arr[pos])] += 1
        pos_counts[pos] = sum(counts.values())
        pos_entropy[pos] = compute_entropy(counts)
        pos_perplexity[pos] = compute_perplexity(pos_entropy[pos])

    print("\nPosition perplexity (lower = more constrained):")
    print(
        f"{'Position':>8s} {'Entropy':>8s} {'Perplexity':>12s} {'Count':>6s} {'Interpretation'}"
    )
    print(f"{'-' * 8} {'-' * 8} {'-' * 12} {'-' * 6} {'-' * 20}")
    # BUG FIX: was range(80) — only displayed first 80. Now FULL length.
    # Print summary stats for all positions + top 20 most/least variable.
    all_pp = pos_perplexity
    print(
        f"  FULL LENGTH: {max_pos} positions, "
        f"mean PP={all_pp.mean():.3f}, "
        f"conserved (<1.5): {(all_pp < 1.5).sum()}, "
        f"variable (>5): {(all_pp > 5).sum()}"
    )
    print("  Top 20 most variable positions:")
    for pos in np.argsort(all_pp)[::-1][:20]:
        interp = (
            "conserved"
            if all_pp[pos] < 1.5
            else "variable"
            if all_pp[pos] > 5
            else "moderate"
        )
        print(
            f"{int(pos):8d} {pos_entropy[int(pos)]:8.3f} {all_pp[int(pos)]:12.3f} {pos_counts[int(pos)]:6d} {interp}"
        )
    print("  Top 10 most conserved positions:")
    for pos in np.argsort(all_pp)[:10]:
        interp = (
            "conserved"
            if all_pp[pos] < 1.5
            else "variable"
            if all_pp[pos] > 5
            else "moderate"
        )
        print(
            f"{int(pos):8d} {pos_entropy[int(pos)]:8.3f} {all_pp[int(pos)]:12.3f} {pos_counts[int(pos)]:6d} {interp}"
        )

    # Compute pairwise conditional perplexity
    print("\n[4/4] Computing pairwise conditional perplexity...")
    print("Conditional perplexity = PP(j | i) = perplexity of j given residue at i")
    print("Co-evolution strength = PP(j) / PP(j | i)")
    print("  High ratio = position i strongly constrains position j")
    print("")

    # Find variable positions
    variable_positions = [p for p in range(max_pos) if pos_perplexity[p] > 3.0]
    print(f"Variable positions (perplexity > 3): {len(variable_positions)}")

    # Compute conditional perplexity for top co-evolutionary pairs
    # BUG FIX: Original code had hardcoded co_evolving_pairs. Now computes
    # them dynamically using MI (mutual information) on variable positions.
    print("  Computing MI to find co-evolving pairs dynamically...")
    co_evolving_pairs = []
    for idx_i, pos_i in enumerate(variable_positions):
        for idx_j in range(idx_i + 1, len(variable_positions)):
            pos_j = variable_positions[idx_j]
            if abs(pos_i - pos_j) > 30:
                continue
            # Compute MI for this pair
            joint = Counter()
            marg_i = Counter()
            marg_j = Counter()
            for arr in pos_arrays[:n_all]:
                if pos_i < len(arr) and pos_j < len(arr):
                    ci, cj = int(arr[pos_i]), int(arr[pos_j])
                    if 0 <= ci < 20 and 0 <= cj < 20:
                        joint[(ci, cj)] += 1
                        marg_i[ci] += 1
                        marg_j[cj] += 1
            total = sum(joint.values())
            if total < 5:
                continue
            mi = sum(
                (c / total)
                * np.log2((c / total) / ((marg_i[ai] / total) * (marg_j[aj] / total)))
                for (ai, aj), c in joint.items()
                if marg_i[ai] > 0 and marg_j[aj] > 0
            )
            if mi > 0.1:
                co_evolving_pairs.append((pos_i, pos_j, mi))

    co_evolving_pairs.sort(key=lambda x: x[2], reverse=True)
    co_evolving_pairs = [(p[0], p[1]) for p in co_evolving_pairs[:20]]
    print(f"  Found {len(co_evolving_pairs)} co-evolving pairs (MI > 0.1)")

    results = []
    for pos_i, pos_j in co_evolving_pairs:
        if pos_i >= max_pos or pos_j >= max_pos:
            continue

        # Marginal perplexity at pos_j
        pp_marginal = pos_perplexity[pos_j]

        # Conditional perplexity: PP(j | i = x) for each x
        conditional_entropies = {}
        for aa_code in range(20):
            counts = Counter()
            for arr in pos_arrays[:n_all]:
                if pos_i < len(arr) and pos_j < len(arr):
                    if int(arr[pos_i]) == aa_code:
                        cj = int(arr[pos_j])
                        if cj >= 0:
                            counts[cj] += 1
            if sum(counts.values()) > 0:
                h_cond = compute_entropy(counts)
                pp_cond = compute_perplexity(h_cond)
                conditional_entropies[aa_code] = (h_cond, pp_cond, sum(counts.values()))

        # Average conditional perplexity
        avg_pp_cond = np.mean(
            [v[1] for v in conditional_entropies.values() if v[2] > 0]
        )

        # Co-evolution ratio
        ratio = pp_marginal / avg_pp_cond if avg_pp_cond > 0 else 0

        print(f"Position pair ({pos_i}, {pos_j}):")
        print(f"  Marginal PP(j) = {pp_marginal:.3f}")
        print(f"  Average conditional PP(j|i) = {avg_pp_cond:.3f}")
        print(f"  Co-evolution ratio = {ratio:.3f}")
        print(f"  (Ratio > 1 means position i constrains position j)")

        # Show which residues at pos_i give lowest conditional PP
        sorted_residues = sorted(conditional_entropies.items(), key=lambda x: x[1][1])
        print(f"  Top 5 residues at pos {pos_i} that constrain pos {pos_j}:")
        for aa_code, (h, pp, n) in sorted_residues[:5]:
            if n > 0:
                print(f"    {aa_list[aa_code]}: PP={pp:.3f}, n={n}")
        print("")

        results.append(
            {
                "pos_i": pos_i,
                "pos_j": pos_j,
                "pp_marginal": float(pp_marginal),
                "pp_conditional": float(avg_pp_cond),
                "ratio": float(ratio),
            }
        )

    # Save results
    summary = {
        "n_sequences": n_all,
        "variable_positions": len(variable_positions),
        "results": results,
    }
    with open(results_dir / "perplexity_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nResults saved to: {results_dir}")


if __name__ == "__main__":
    main()
