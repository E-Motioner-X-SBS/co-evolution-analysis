#!/usr/bin/env python3
"""
Predictive Constraint Function for Co-evolution
================================================

Implements the three K-map approaches:
A. Original: observed pairs (co-evolutionary motifs)
B. Flipped: forbidden pairs (destabilization constraints)
C. Continuous: frequency-based prediction function

The constraint function predicts whether a mutation will be
destabilizing based on the K-map structure.
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


def compute_frequency_kmap(pos_arrays, pos_i, pos_j, n_seqs):
    """Build 20x20 frequency K-map."""
    kmap = np.zeros((20, 20), dtype=np.float64)
    for arr in pos_arrays[:n_seqs]:
        if pos_i < len(arr) and pos_j < len(arr):
            ci, cj = int(arr[pos_i]), int(arr[pos_j])
            if ci >= 0 and cj >= 0:
                kmap[ci, cj] += 1
    total = kmap.sum()
    if total > 0:
        kmap /= total
    return kmap


def compute_constraint_function(kmap_freq):
    """
    Compute constraint function C(aa_i, aa_j) = ln(P/P_expected).

    SIGN CONVENTION (verified against DCA literature, Morcos et al 2011 PNAS):
      C > 0 → pair is MORE common than expected (co-evolutionary)
      C < 0 → pair is LESS common than expected (anti-correlated)
      C = 0 → pair occurs at random frequency (independent)

    NOTE: Previous formula used -ln(P/P_exp) giving OPPOSITE sign — fixed.
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


def sigmoid(x):
    """Sigmoid function for probability."""
    return 1.0 / (1.0 + np.exp(-np.clip(x, -10, 10)))


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    results_dir = base_dir / "constraint_function_results"
    results_dir.mkdir(exist_ok=True)

    print("=" * 80)
    print("Predictive Constraint Function for Co-evolution")
    print("=" * 80)

    # Load ALL sequences
    print("\n[1/5] Loading ALL sequences...")
    sequences = parse_fasta(fasta_file)
    encoder = Base20AminoEncoder(version=1)
    aa_list = list(AMINO_HE_2012)
    n_all = len(sequences)
    full_length = len(sequences[0][1])
    print(f"  Total: {n_all} sequences")

    # Build position arrays
    print("\n[2/5] Building position arrays...")
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

    # Find co-evolutionary pairs
    print("\n[3/5] Finding co-evolutionary pairs...")
    co_evolving = []
    for idx_i, pos_i in enumerate(variable_positions):
        for idx_j in range(idx_i + 1, len(variable_positions)):
            pos_j = variable_positions[idx_j]
            if abs(pos_i - pos_j) > 30:
                continue

            ref_i = get_majority_ref(pos_arrays, pos_i, n_all)
            ref_j = get_majority_ref(pos_arrays, pos_j, n_all)
            joint, marg_i, marg_j = Counter(), Counter(), Counter()
            for arr in pos_arrays[:n_all]:
                if pos_i < len(arr) and pos_j < len(arr):
                    ci, cj = int(arr[pos_i]), int(arr[pos_j])
                    if ci >= 0 and cj >= 0 and (ci != ref_i or cj != ref_j):
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
                co_evolving.append((pos_i, pos_j, mi, total, ref_i, ref_j))

    co_evolving.sort(key=lambda x: x[2], reverse=True)
    print(f"  Co-evolutionary pairs: {len(co_evolving)}")

    # Build constraint functions for top pairs
    print("\n[4/5] Building constraint functions...")
    all_results = []

    for pos_i, pos_j, mi, n_muts, ref_i, ref_j in co_evolving[:10]:
        # BUG FIX: ref_i, ref_j are already int codes (0-19) from get_majority_ref
        ref_i_code = int(ref_i)
        ref_j_code = int(ref_j)
        ref_aa_i = aa_list[ref_i_code] if 0 <= ref_i_code < 20 else "?"
        ref_aa_j = aa_list[ref_j_code] if 0 <= ref_j_code < 20 else "?"

        # Frequency K-map
        kmap_freq = compute_frequency_kmap(pos_arrays, pos_i, pos_j, n_all)

        # Constraint function
        C = compute_constraint_function(kmap_freq)

        # Prediction function (sigmoid)
        P = sigmoid(C)

        print(f"\n  Position pair ({pos_i}, {pos_j}):")
        print(f"    MI: {mi:.4f}")

        # Top co-evolutionary pairs (C > 0)
        top_co = []
        for ai in range(20):
            for aj in range(20):
                if kmap_freq[ai, aj] > 0.005 and C[ai, aj] > 0:
                    top_co.append(
                        (
                            aa_list[ai],
                            aa_list[aj],
                            C[ai, aj],
                            P[ai, aj],
                            kmap_freq[ai, aj],
                        )
                    )
        top_co.sort(key=lambda x: x[2], reverse=True)

        print(f"    Top co-evolutionary (C > 0):")
        for aa1, aa2, c, p, freq in top_co[:5]:
            print(f"      {aa1}-{aa2}: C={c:.4f}, P={p:.4f}, freq={freq:.4f}")

        # Top anti-correlated (C < 0)
        top_anti = []
        for ai in range(20):
            for aj in range(20):
                if kmap_freq[ai, aj] > 0.001 and C[ai, aj] < 0:
                    top_anti.append(
                        (
                            aa_list[ai],
                            aa_list[aj],
                            C[ai, aj],
                            P[ai, aj],
                            kmap_freq[ai, aj],
                        )
                    )
        top_anti.sort(key=lambda x: x[2])

        print(f"    Top anti-correlated (C < 0):")
        for aa1, aa2, c, p, freq in top_anti[:5]:
            print(f"      {aa1}-{aa2}: C={c:.4f}, P={p:.4f}, freq={freq:.4f}")

        all_results.append(
            {
                "pos_i": pos_i,
                "pos_j": pos_j,
                "mi": mi,
                "ref_i": ref_aa_i,
                "ref_j": ref_aa_j,
                "top_co": [
                    (a1, a2, float(c), float(p), float(f))
                    for a1, a2, c, p, f in top_co[:10]
                ],
                "top_anti": [
                    (a1, a2, float(c), float(p), float(f))
                    for a1, a2, c, p, f in top_anti[:10]
                ],
            }
        )

    # Test prediction
    print("\n[5/5] Testing prediction...")
    # For each position pair, predict: given mutation at i, what's the best j?
    correct = 0
    total = 0

    for pos_i, pos_j, mi, n_muts, ref_i, ref_j in co_evolving[:10]:
        # BUG FIX: ref_i, ref_j are already int codes (0-19) from get_majority_ref
        ref_i_code = int(ref_i)
        ref_j_code = int(ref_j)
        # Build constraint function from first 800 sequences
        kmap_freq = np.zeros((20, 20), dtype=np.float64)
        for arr in pos_arrays[:800]:
            if pos_i < len(arr) and pos_j < len(arr):
                ci, cj = int(arr[pos_i]), int(arr[pos_j])
                if ci >= 0 and cj >= 0:
                    kmap_freq[ci, cj] += 1
        total_freq = kmap_freq.sum()
        if total_freq > 0:
            kmap_freq /= total_freq

        C = compute_constraint_function(kmap_freq)

        # Test on sequences 800-1299
        for arr in pos_arrays[800:]:
            if pos_i < len(arr) and pos_j < len(arr):
                ci, cj = int(arr[pos_i]), int(arr[pos_j])
                if ci >= 0 and cj >= 0:
                    # BUG FIX: ref_i, ref_j are already int codes — no lookup needed
                    if ci != ref_i_code or cj != ref_j_code:
                        # This is a mutation pair
                        # Predict: what's the best cj for this ci?
                        best_cj = np.argmax(C[ci, :])
                        if best_cj == cj:
                            correct += 1
                        total += 1

    accuracy = correct / total if total > 0 else 0
    print(f"  Prediction accuracy: {accuracy:.4f} ({correct}/{total})")

    # Save results
    print("\n" + "=" * 80)
    print(f"Results saved to: {results_dir}")

    summary = {
        "dataset": f"{n_all} sequences",
        "variable_positions": len(variable_positions),
        "co_evolving_pairs": len(co_evolving),
        "prediction_accuracy": accuracy,
        "results": all_results,
    }
    with open(results_dir / "constraint_function_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nPrediction accuracy: {accuracy:.4f}")
    print(f"The constraint function C(aa_i, aa_j) = ln(P/P_expected)")
    print(f"predicts co-evolutionary pairs with {accuracy * 100:.1f}% accuracy.")


if __name__ == "__main__":
    main()
