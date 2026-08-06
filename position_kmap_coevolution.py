#!/usr/bin/env python3
"""
Position-Based K-map Analysis for Co-evolution Prediction
=========================================================

Key insight: Co-evolution happens BETWEEN positions, not between residues.
We build K-maps from position-pair statistics:
- Position (i,j) has a K-map showing which residue pairs appear at those positions
- The structure of this K-map captures co-evolutionary constraints

This approach directly identifies which residue pairs co-evolve at which positions.
"""

import sys
import json
import numpy as np
from collections import Counter, defaultdict
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


def build_position_kmaps(sequences, encoder, max_positions=None, n_seqs=200):
    """
    Build position-pair K-maps.

    For each pair of positions (i, j), build a 20x20 frequency matrix
    showing which amino acid pairs appear at those positions across
    all sequences.

    This captures co-evolutionary constraints directly.

    max_positions=None → FULL sequence length (all positions).
    """
    print("\n=== Building Position-Pair K-maps ===")

    n = min(n_seqs, len(sequences))
    min_len = min(len(seq) for _, seq in sequences[:n])
    if max_positions is None:
        max_pos = min_len  # FULL length
    else:
        max_pos = min(max_positions, min_len)

    # Build frequency matrices for each position pair
    # We'll use a sliding window of position pairs
    position_kmaps = {}

    # First, build per-position frequency vectors
    pos_freq = np.zeros((max_pos, 20), dtype=np.float64)
    pos_counts = np.zeros(max_pos, dtype=np.int64)

    for i in range(n):
        _, seq = sequences[i]
        clean_seq = "".join(aa for aa in seq if aa in encoder.encode)
        for j in range(min(max_pos, len(clean_seq))):
            aa = clean_seq[j]
            if aa in encoder.encode:
                code = encoder.encode[aa]
                pos_freq[j, code] += 1
                pos_counts[j] += 1

    # Normalize per-position
    for j in range(max_pos):
        if pos_counts[j] > 0:
            pos_freq[j] /= pos_counts[j]

    # Build position-pair K-maps (20x20 for each pair)
    # For efficiency, only build for nearby positions (within window)
    window = 20  # Check co-evolution within 20-residue window

    print(f"  Sequences: {n}")
    print(f"  Positions analyzed: {max_pos}")
    print(f"  Co-evolution window: {window}")

    return pos_freq, pos_counts, max_pos, window


def compute_position_pair_kmap(sequences, encoder, pos_i, pos_j, n_seqs=200):
    """
    Build a 20x20 K-map for a specific position pair (i, j).

    Each cell (aa1, aa2) counts how many times amino acid aa1 appears
    at position i and aa2 at position j across all sequences.
    """
    n = min(n_seqs, len(sequences))
    kmap = np.zeros((20, 20), dtype=np.float64)

    for i in range(n):
        _, seq = sequences[i]
        clean_seq = "".join(aa for aa in seq if aa in encoder.encode)

        if pos_i < len(clean_seq) and pos_j < len(clean_seq):
            aa_i = clean_seq[pos_i]
            aa_j = clean_seq[pos_j]
            if aa_i in encoder.encode and aa_j in encoder.encode:
                code_i = encoder.encode[aa_i]
                code_j = encoder.encode[aa_j]
                kmap[code_i, code_j] += 1

    # Normalize
    total = kmap.sum()
    if total > 0:
        kmap /= total

    return kmap


def compute_mutual_information(sequences, pos_i, pos_j, encoder, n_seqs=200):
    """
    Compute mutual information between positions i and j.
    """
    n = min(n_seqs, len(sequences))
    joint = Counter()
    marg_i = Counter()
    marg_j = Counter()

    for i in range(n):
        _, seq = sequences[i]
        clean_seq = "".join(aa for aa in seq if aa in encoder.encode)

        if pos_i < len(clean_seq) and pos_j < len(clean_seq):
            aa_i = clean_seq[pos_i]
            aa_j = clean_seq[pos_j]
            if aa_i in encoder.encode and aa_j in encoder.encode:
                joint[(aa_i, aa_j)] += 1
                marg_i[aa_i] += 1
                marg_j[aa_j] += 1

    total = sum(joint.values())
    if total == 0:
        return 0.0

    mi = 0.0
    for (ai, aj), count in joint.items():
        p_joint = count / total
        p_i = marg_i[ai] / total
        p_j = marg_j[aj] / total
        if p_joint > 0 and p_i > 0 and p_j > 0:
            mi += p_joint * np.log2(p_joint / (p_i * p_j))

    return mi


def minimize_position_kmap(kmap_2d):
    """
    Run Quine-McCluskey on a position-pair K-map.
    """
    # Threshold to Boolean
    nonzero = kmap_2d[kmap_2d > 0]
    if len(nonzero) == 0:
        return None

    threshold = np.percentile(nonzero, 70)
    kmap_bool = (kmap_2d >= threshold).astype(int)

    # Flatten and minimize
    bool_flat = kmap_bool.flatten().astype(int)
    result = boolean_minimize_kmap(bool_flat, algorithm="qm")

    return {
        "result": result,
        "kmap_bool": kmap_bool,
        "threshold": threshold,
        "n_on_set": int(kmap_bool.sum()),
    }


def compute_coupling_from_kmap(kmap_2d, pos_i, pos_j, encoder):
    """
    Compute coupling constant J_ij from position-pair K-map.

    J_ij = log(P(aa_i, aa_j) / (P(aa_i) * P(aa_j)))

    High |J| = strong co-evolutionary constraint.
    """
    # Marginal frequencies
    row_marg = kmap_2d.sum(axis=1)
    col_marg = kmap_2d.sum(axis=0)

    # Coupling matrix
    epsilon = 1e-10
    with np.errstate(divide="ignore", invalid="ignore"):
        coupling = np.log(
            (kmap_2d + epsilon)
            / ((row_marg[:, None] + epsilon) * (col_marg[None, :] + epsilon))
        )

    # Average coupling strength
    avg_coupling = np.mean(np.abs(coupling))
    max_coupling = np.max(np.abs(coupling))

    # Find strongest positive and negative couplings
    aa_list = list(AMINO_HE_2012)

    max_pos_idx = np.unravel_index(np.argmax(coupling), coupling.shape)
    max_neg_idx = np.unravel_index(np.argmin(coupling), coupling.shape)

    return {
        "coupling_matrix": coupling.tolist(),
        "avg_coupling": float(avg_coupling),
        "max_coupling": float(max_coupling),
        "strongest_positive": {
            "aa_i": aa_list[max_pos_idx[0]],
            "aa_j": aa_list[max_pos_idx[1]],
            "J": float(coupling[max_pos_idx]),
        },
        "strongest_negative": {
            "aa_i": aa_list[max_neg_idx[0]],
            "aa_j": aa_list[max_neg_idx[1]],
            "J": float(coupling[max_neg_idx]),
        },
    }


def predict_coevolution_from_position_kmaps(
    position_kmaps, sequences, encoder, max_positions, n_seqs=100
):
    """
    Test if position-pair K-maps predict co-evolution.

    For each position pair (i,j):
    1. Build the K-map
    2. Compute the Boolean minimized function
    3. Check if "on-set" residues correlate with high MI
    """
    print("\n=== Predicting Co-evolution from Position K-maps ===")

    n = min(n_seqs, len(sequences))

    # Build clean sequences
    clean_seqs = []
    for i in range(n):
        _, seq = sequences[i]
        clean = "".join(aa for aa in seq if aa in encoder.encode)
        if len(clean) > 100:
            clean_seqs.append(clean)

    if len(clean_seqs) < 5:
        return {}

    # Test prediction for several position pairs
    test_pairs = []
    window = 20

    for i in range(0, min(50, max_positions), 5):
        for j in range(i + 1, min(i + window, max_positions)):
            mi = compute_mutual_information(sequences, i, j, encoder, n_seqs)
            if mi > 0.01:  # Only test pairs with non-trivial MI
                test_pairs.append((i, j, mi))

    test_pairs.sort(key=lambda x: x[2], reverse=True)

    print(f"  Position pairs with MI > 0.01: {len(test_pairs)}")
    print(f"  Top 10 co-evolving position pairs:")
    for i, j, mi in test_pairs[:10]:
        print(f"    Position ({i}, {j}): MI = {mi:.4f}")

    # For top pairs, build K-map and check prediction
    correct = 0
    total = 0

    for pos_i, pos_j, mi in test_pairs[:20]:
        kmap = compute_position_pair_kmap(sequences, encoder, pos_i, pos_j, n_seqs)

        # Threshold to Boolean
        nonzero = kmap[kmap > 0]
        if len(nonzero) == 0:
            continue

        threshold = np.percentile(nonzero, 70)
        on_set_mask = kmap >= threshold

        # Compute average MI for on-set vs off-set residue pairs
        on_set_mi = []
        off_set_mi = []

        for aa_i_code in range(20):
            for aa_j_code in range(20):
                if kmap[aa_i_code, aa_j_code] > 0:
                    # Get MI contribution for this residue pair
                    pair_count = 0
                    for seq in clean_seqs:
                        if pos_i < len(seq) and pos_j < len(seq):
                            if (
                                encoder.encode.get(seq[pos_i]) == aa_i_code
                                and encoder.encode.get(seq[pos_j]) == aa_j_code
                            ):
                                pair_count += 1

                    # Weight by frequency
                    pair_freq = pair_count / len(clean_seqs)

                    if on_set_mask[aa_i_code, aa_j_code]:
                        on_set_mi.append(pair_freq)
                    else:
                        off_set_mi.append(pair_freq)

        avg_on = np.mean(on_set_mi) if on_set_mi else 0
        avg_off = np.mean(off_set_mi) if off_set_mi else 0

        # Prediction: on-set residue pairs should be more frequent
        if avg_on > avg_off:
            correct += 1
        total += 1

    accuracy = correct / total if total > 0 else 0

    print(f"\n  Prediction accuracy: {accuracy:.4f} ({correct}/{total})")

    return {
        "n_test_pairs": len(test_pairs),
        "top_pairs": [(i, j, mi) for i, j, mi in test_pairs[:10]],
        "prediction_accuracy": accuracy,
        "n_correct": correct,
        "n_total": total,
    }


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    results_dir = base_dir / "position_kmap_results"
    results_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("Position-Based K-map Co-evolution Analysis")
    print("=" * 70)

    # Parse
    print("\n[1/5] Parsing FASTA...")
    sequences = parse_fasta(fasta_file)
    full_length = len(sequences[0][1])
    encoder = Base20AminoEncoder(version=1)
    print(f"  Loaded {len(sequences)} sequences")

    # Build position frequency vectors — FULL LENGTH
    print("\n[2/5] Building position frequency vectors...")
    pos_freq, pos_counts, max_pos, window = build_position_kmaps(
        sequences, encoder, max_positions=None, n_seqs=1299
    )

    # Find co-evolving positions
    print("\n[3/5] Finding co-evolving position pairs (FULL length, GPU)...")
    n_seqs = min(1299, len(sequences))

    # Compute MI for all nearby position pairs — GPU accelerated
    try:
        import coevolution_gpu as cg

        # Build position arrays
        pos_arrays = []
        for _, seq in sequences[:n_seqs]:
            clean = "".join(aa for aa in seq if aa in encoder.encode)
            arr = np.array([encoder.encode.get(aa, -1) for aa in clean], dtype=np.int32)
            pos_arrays.append(arr)
        dense = cg.dense_to_gpu(pos_arrays)
        pairs = cg.all_pairs(dense.shape[1], max_gap=window)
        mi_dict, _ = cg.mi_matrix_gpu(dense, pairs, min_total=10, chunk=16384)
        mi_results = [(i, j, mi) for (i, j), mi in mi_dict.items() if mi > 0.005]
        print(f"  GPU MI computed for {len(mi_dict)} pairs")
    except Exception as e:
        print(f"  GPU failed ({e}), using CPU...")
        mi_results = []
        for i in range(0, max_pos):
            for j in range(i + 1, min(i + window, max_pos)):
                mi = compute_mutual_information(sequences, i, j, encoder, n_seqs)
                if mi > 0.005:
                    mi_results.append((i, j, mi))

    mi_results.sort(key=lambda x: x[2], reverse=True)

    print(f"  Position pairs with MI > 0.005: {len(mi_results)}")
    print(f"  Top 15 co-evolving pairs:")
    for i, j, mi in mi_results[:15]:
        print(f"    ({i:3d}, {j:3d}): MI = {mi:.4f}")

    # Build K-maps for top co-evolving pairs
    print("\n[4/5] Building position-pair K-maps...")
    aa_list = list(AMINO_HE_2012)

    position_results = []
    for idx, (pos_i, pos_j, mi) in enumerate(mi_results[:10]):
        print(f"\n  --- Position pair ({pos_i}, {pos_j}), MI = {mi:.4f} ---")

        # Build K-map
        kmap = compute_position_pair_kmap(sequences, encoder, pos_i, pos_j, n_seqs)

        # Minimize
        minimization = minimize_position_kmap(kmap)
        if minimization is None:
            continue

        result = minimization["result"]

        # Compute coupling
        coupling = compute_coupling_from_kmap(kmap, pos_i, pos_j, encoder)

        print(f"    K-map: {minimization['n_on_set']} on-set cells")
        print(f"    Prime implicants: {result['n_prime_implicants']}")
        print(f"    Essential PIs: {result['n_essential']}")
        print(f"    Avg coupling: {coupling['avg_coupling']:.4f}")
        print(
            f"    Strongest +: {coupling['strongest_positive']['aa_i']}-"
            f"{coupling['strongest_positive']['aa_j']} "
            f"(J={coupling['strongest_positive']['J']:.4f})"
        )
        print(
            f"    Strongest -: {coupling['strongest_negative']['aa_i']}-"
            f"{coupling['strongest_negative']['aa_j']} "
            f"(J={coupling['strongest_negative']['J']:.4f})"
        )

        position_results.append(
            {
                "position_pair": (pos_i, pos_j),
                "mi": mi,
                "n_on_set": minimization["n_on_set"],
                "n_prime_implicants": result["n_prime_implicants"],
                "n_essential": result["n_essential"],
                "coupling": coupling,
                "kmap_sum": float(kmap.sum()),
            }
        )

    # Predict co-evolution
    print("\n[5/5] Predicting co-evolution...")
    prediction = predict_coevolution_from_position_kmaps(
        None, sequences, encoder, max_pos, n_seqs=100
    )

    # Save results
    print("\n" + "=" * 70)
    print("Saving Results...")
    print("=" * 70)

    summary = {
        "dataset": "SARS-CoV-2 Spike Protein Position-Based K-map",
        "encoding": "Base-20 (He 2012)",
        "num_sequences": len(sequences),
        "max_positions": max_pos,
        "mi_results": [
            {"pos_i": i, "pos_j": j, "mi": mi} for i, j, mi in mi_results[:30]
        ],
        "position_results": position_results,
        "prediction": prediction,
    }

    with open(results_dir / "position_kmap_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nResults saved to: {results_dir}")

    # Final summary
    print("\n" + "=" * 70)
    print("POSITION-BASED K-MAP ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nKey Results:")
    print(f"  Co-evolving position pairs found: {len(mi_results)}")
    print(
        f"  Top pair: ({mi_results[0][0]}, {mi_results[0][1]}) with MI = {mi_results[0][2]:.4f}"
    )
    print(f"  Position K-maps analyzed: {len(position_results)}")
    print(f"  Prediction accuracy: {prediction.get('prediction_accuracy', 0):.4f}")

    if position_results:
        print(f"\n  Coupling constants from top position pairs:")
        for pr in position_results[:5]:
            i, j = pr["position_pair"]
            print(
                f"    J({i},{j}) = {pr['coupling']['avg_coupling']:.4f} "
                f"[+{pr['coupling']['strongest_positive']['aa_i']}-"
                f"{pr['coupling']['strongest_positive']['aa_j']}, "
                f"-{pr['coupling']['strongest_negative']['aa_i']}-"
                f"{pr['coupling']['strongest_negative']['aa_j']}]"
            )


if __name__ == "__main__":
    main()
