#!/usr/bin/env python3
"""
All-Sequence Constraint Function for Co-evolution
==================================================

Builds the constraint function on ALL 1299 Omicron sequences.
Uses the Wuhan reference as baseline and tests prediction
within the same dataset via leave-one-out cross-validation.
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
        if pos < len(arr) and 0 <= arr[pos] < 20:
            counts[int(arr[pos])] += 1
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((c / total) * np.log2(c / total) for c in counts.values() if c > 0)


def get_majority_ref(pos_arrays, pos, n_seqs):
    return Counter(
        int(a[pos]) for a in pos_arrays[:n_seqs] if pos < len(a) and 0 <= a[pos] < 20
    ).most_common(1)[0][0]


def compute_frequency_kmap(pos_arrays, pos_i, pos_j, n_seqs, exclude_idx=None):
    """Build frequency K-map, optionally excluding one sequence."""
    kmap = np.zeros((20, 20), dtype=np.float64)
    for idx, arr in enumerate(pos_arrays[:n_seqs]):
        if exclude_idx is not None and idx == exclude_idx:
            continue
        if pos_i < len(arr) and pos_j < len(arr):
            ci, cj = int(arr[pos_i]), int(arr[pos_j])
            if 0 <= ci < 20 and 0 <= cj < 20:
                kmap[ci, cj] += 1
    total = kmap.sum()
    if total > 0:
        kmap /= total
    return kmap


def compute_constraint_function(kmap_freq):
    """C(aa_i, aa_j) = ln(P / P_expected) — constraint / coupling function.

    SIGN CONVENTION (verified against DCA literature, Morcos et al 2011 PNAS):
      C > 0 → pair is MORE common than expected (co-evolutionary)
      C < 0 → pair is LESS common than expected (anti-correlated)
      C = 0 → pair occurs at random frequency (independent)

    Prediction: P_co-evolution = σ(C) = 1/(1+e^{-C})
    """
    marg_i = kmap_freq.sum(axis=1)
    marg_j = kmap_freq.sum(axis=0)
    epsilon = 1e-10
    with np.errstate(divide="ignore", invalid="ignore"):
        C = np.log(
            (kmap_freq + epsilon)
            / ((marg_i[:, None] + epsilon) * (marg_j[None, :] + epsilon))
        )
    return C


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    results_dir = base_dir / "allseq_constraint_results"
    results_dir.mkdir(exist_ok=True)

    print("=" * 80)
    print("All-Sequence Constraint Function (ALL 1299 Omicron sequences)")
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
        arr = np.array(
            [encoder.encode.get(aa, 20) for aa in seq[:max_pos]], dtype=np.int32
        )
        pos_arrays.append(arr)

    # Find variable positions
    variable_positions = [
        p for p in range(max_pos) if compute_entropy(pos_arrays, p, n_all) > 0.3
    ]
    print(f"  Variable positions: {len(variable_positions)}")

    # Find co-evolutionary pairs (GPU-accelerated via coevolution_shared)
    from coevolution_shared import find_coevolving_pairs_gpu
    co_evolving = find_coevolving_pairs_gpu(
        pos_arrays, variable_positions, n_all, max_gap=30, min_mi=0.1
    )
    print(f"  Co-evolutionary pairs: {len(co_evolving)}")

    # Build constraint function on ALL sequences
    print("\n[4/4] Building constraint function on ALL sequences...")

    # Use leave-one-out cross-validation
    print("\nLeave-One-Out Cross-Validation:")
    print("For each sequence, build constraint function on ALL OTHER sequences,")
    print("then predict mutations for the held-out sequence.")

    all_results = []

    for pos_i, pos_j, mi, n_muts, ref_i, ref_j in co_evolving[:10]:
        # BUG FIX: ref_i and ref_j are ALREADY integer codes (0-19) returned
        # by get_majority_ref(). The previous code did aa_list.index(ref_i)
        # which looked up an INT in a list of STRINGS → always returned 0.
        # This meant ALL pairs were tested against Alanine (code 0) as reference,
        # invalidating the LOO-CV results.
        ref_i_code = int(ref_i)  # Already an int (0-19) from get_majority_ref
        ref_j_code = int(ref_j)  # Already an int (0-19) from get_majority_ref
        # For display: convert code back to amino acid letter
        ref_aa_i = aa_list[ref_i_code] if 0 <= ref_i_code < 20 else "?"
        ref_aa_j = aa_list[ref_j_code] if 0 <= ref_j_code < 20 else "?"

        correct = 0
        total = 0

        for holdout_idx in range(n_all):
            # Build constraint function on ALL OTHER sequences
            kmap_freq = compute_frequency_kmap(
                pos_arrays, pos_i, pos_j, n_all, exclude_idx=holdout_idx
            )

            if kmap_freq.sum() == 0:
                continue

            C = compute_constraint_function(kmap_freq)

            # Test on held-out sequence
            arr = pos_arrays[holdout_idx]
            if pos_i < len(arr) and pos_j < len(arr):
                ci, cj = int(arr[pos_i]), int(arr[pos_j])
                if 0 <= ci < 20 and 0 <= cj < 20:
                    if ci != ref_i_code or cj != ref_j_code:
                        # This is a mutation
                        # Predict: what's the best cj for this ci?
                        best_cj = np.argmax(C[ci, :])
                        if best_cj == cj:
                            correct += 1
                        total += 1

        accuracy = correct / total if total > 0 else 0
        all_results.append(
            {
                "pos_i": pos_i,
                "pos_j": pos_j,
                "mi": mi,
                "ref_i": ref_aa_i,
                "ref_j": ref_aa_j,
                "n_mutations": total,
                "correct": correct,
                "accuracy": accuracy,
            }
        )

        print(
            f"  ({pos_i:2d},{pos_j:2d}): MI={mi:.4f}, accuracy={accuracy:.4f} ({correct}/{total})"
        )

    # Overall accuracy
    total_correct = sum(r["correct"] for r in all_results)
    total_test = sum(r["n_mutations"] for r in all_results)
    overall_acc = total_correct / total_test if total_test > 0 else 0

    print(
        f"\nOverall LOO-CV accuracy: {overall_acc:.4f} ({total_correct}/{total_test})"
    )

    # Save results
    summary = {
        "dataset": f"{n_all} Omicron sequences",
        "method": "Leave-One-Out Cross-Validation",
        "variable_positions": len(variable_positions),
        "co_evolving_pairs": len(co_evolving),
        "overall_accuracy": overall_acc,
        "total_correct": total_correct,
        "total_test": total_test,
        "results": all_results,
    }

    with open(results_dir / "allseq_constraint_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nResults saved to: {results_dir}")


if __name__ == "__main__":
    main()
