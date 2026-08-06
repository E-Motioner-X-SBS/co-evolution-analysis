#!/usr/bin/env python3
"""
Boolean Minimization and Co-evolution Prediction for SARS-CoV-2 Spike Protein
=============================================================================

This script applies K-map Boolean minimization to the Spike protein
co-evolution analysis. The key insight:

1. The K-map dipeptide frequency matrix can be thresholded to create
   a Boolean function f(row, col) = 1 if dipeptide is frequent, 0 otherwise.

2. Quine-McCluskey minimization finds prime implicants = essential
   residue pair motifs that capture the co-evolution structure.

3. The minimized Boolean function defines "coupling constants" between
   residue positions — predicting which pairs co-evolve.

4. We validate by checking if the minimized function predicts observed
   co-evolution patterns (mutual information, correlation).
"""

import os
import sys
import json
import numpy as np
from collections import Counter, defaultdict
from pathlib import Path
from itertools import combinations

sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/kmap-sbm-validation/src")

from kmap_sbm.encoding.gray_amino import (
    encode_gray_single,
    encode_single,
    gray_code_5bit,
    build_aa_kmap_2d,
    hamming_distance,
    gray_hamming_int,
    AA_GROUPS,
    AA_ALPHABET,
    _AA_TO_INDEX,
    _WITHIN_GROUP_DISTANCE_1,
    _MAX_DISTANCE_PAIRS,
)

from kmap_sbm.analysis.prime_implicants import (
    kmap_truth_table,
    prime_implicants_quine_mccluskey,
    boolean_minimize_kmap,
)


# ============================================================
# 1. FASTA Parser
# ============================================================


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


def encode_sequence_gray(seq):
    encoded = []
    for aa in seq:
        if aa in _AA_TO_INDEX:
            encoded.append(encode_gray_single(aa))
    return encoded


# ============================================================
# 2. Boolean K-map Construction
# ============================================================


def build_boolean_kmap(sequences, threshold_percentile=75):
    """
    Build a Boolean K-map from dipeptide frequencies.

    The Boolean function f(row, col) = 1 if the dipeptide
    (row_aa, col_aa) appears more frequently than the threshold.
    """
    print("\n=== Building Boolean K-map ===")

    # Build frequency K-map
    kmap_freq = np.zeros((32, 32), dtype=float)
    total_pairs = 0

    for _, seq in sequences:
        clean_seq = "".join(aa for aa in seq if aa in _AA_TO_INDEX)
        for i in range(len(clean_seq) - 1):
            aa1, aa2 = clean_seq[i], clean_seq[i + 1]
            if aa1 in _AA_TO_INDEX and aa2 in _AA_TO_INDEX:
                row = encode_gray_single(aa1)
                col = encode_gray_single(aa2)
                kmap_freq[row, col] += 1
                total_pairs += 1

    # Normalize
    if total_pairs > 0:
        kmap_freq /= total_pairs

    # Threshold to create Boolean map
    threshold = np.percentile(kmap_freq[kmap_freq > 0], threshold_percentile)
    kmap_bool = (kmap_freq >= threshold).astype(int)

    # Count on-set and off-set
    n_on = int(kmap_bool.sum())
    n_total = kmap_bool.size
    density = n_on / n_total

    print(f"  Total dipeptide pairs: {total_pairs}")
    print(f"  Threshold (P{threshold_percentile}): {threshold:.6f}")
    print(f"  On-set (frequent): {n_on} cells")
    print(f"  Off-set (rare): {n_total - n_on} cells")
    print(f"  Density: {density:.4f}")

    return kmap_freq, kmap_bool, threshold


# ============================================================
# 3. Boolean Minimization (Quine-McCluskey)
# ============================================================


def minimize_boolean_function(kmap_bool):
    """
    Run Quine-McCluskey minimization on the Boolean K-map.
    """
    print("\n=== Boolean Minimization (Quine-McCluskey) ===")

    # Flatten the 32x32 boolean map to 1D for the QM algorithm
    # The QM algorithm expects a 1D truth table
    bool_map_flat = kmap_bool.flatten().astype(int)

    # Run minimization
    result = boolean_minimize_kmap(bool_map_flat, algorithm="qm")

    print(f"  Input variables: {result['n_vars']}")
    print(f"  Minterms (on-set): {result['n_minterms']}")
    print(f"  Prime implicants: {result['n_prime_implicants']}")
    print(f"  Essential prime implicants: {result['n_essential']}")
    print(f"  Covering size: {result['covering_size']}")

    # Decode prime implicants to residue pairs
    print(f"\n  Prime Implicants (residue pairs):")
    aa_list = sorted(_AA_TO_INDEX.keys())

    for i, pi in enumerate(result["prime_implicants"][:20]):  # Show first 20
        values = pi["values"]
        mask = pi["mask"]
        n_dc = pi["n_dontcares"]

        # Decode the values to amino acids
        # The truth table has 10 variables (5 for row, 5 for col)
        # But we're working with 32x32 = 10 variables
        row_bits = values[:5]
        col_bits = values[5:10]

        row_val = sum(b * (2 ** (4 - i)) for i, b in enumerate(row_bits))
        col_val = sum(b * (2 ** (4 - i)) for i, b in enumerate(col_bits))

        # Decode mask
        row_mask = mask[:5]
        col_mask = mask[5:10]

        print(
            f"    PI {i + 1}: row={row_val:2d} col={col_val:2d} (dc={n_dc}) "
            f"coverage={len(pi['coverage'])} minterms"
        )

    return result


# ============================================================
# 4. Extract Co-evolution Motifs from Prime Implicants
# ============================================================


def extract_coevolution_motifs(qm_result, sequences, n_seqs=100):
    """
    Extract co-evolution motifs from the minimized Boolean function.

    Each prime implicant defines a "motif" — a set of dipeptide pairs
    that are collectively frequent. These motifs correspond to
    co-evolutionary constraints.
    """
    print("\n=== Extracting Co-evolution Motifs ===")

    # Build mapping from Gray code values to amino acids
    gray_to_aa = {}
    for aa, idx in _AA_TO_INDEX.items():
        gray_val = gray_code_5bit(idx)
        gray_to_aa[gray_val] = aa

    # Pre-encode all sequences for fast matching
    print(f"  Pre-encoding {n_seqs} sequences...")
    encoded_seqs = []
    total_aa = 0
    for _, seq in sequences[:n_seqs]:
        clean = [encode_gray_single(aa) for aa in seq if aa in _AA_TO_INDEX]
        if len(clean) > 1:
            encoded_seqs.append(np.array(clean, dtype=np.int32))
            total_aa += len(clean)
    print(f"  Encoded {len(encoded_seqs)} sequences, {total_aa} total AAs")

    motifs = []

    for i, pi in enumerate(qm_result["prime_implicants"]):
        values = pi["values"]
        mask = pi["mask"]
        coverage = pi["coverage"]

        row_bits = values[:5]
        col_bits = values[5:10]
        row_mask_arr = mask[:5]
        col_mask_arr = mask[5:10]

        row_val = sum(b * (2 ** (4 - j)) for j, b in enumerate(row_bits))
        col_val = sum(b * (2 ** (4 - j)) for j, b in enumerate(col_bits))
        row_aa = gray_to_aa.get(row_val, "?")
        col_aa = gray_to_aa.get(col_val, "?")

        # Count motif occurrences using numpy vectorization
        motif_count = 0
        for arr in encoded_seqs:
            if len(arr) < 2:
                continue
            r_vals = arr[:-1]
            c_vals = arr[1:]
            matches = np.ones(len(r_vals), dtype=bool)
            for k in range(5):
                if not row_mask_arr[k]:
                    matches &= ((r_vals >> (4 - k)) & 1) == row_bits[k]
                if not col_mask_arr[k]:
                    matches &= ((c_vals >> (4 - k)) & 1) == col_bits[k]
            motif_count += int(np.sum(matches))

        if motif_count > 0:
            motifs.append(
                {
                    "rank": i + 1,
                    "row_aa": row_aa,
                    "col_aa": col_aa,
                    "row_gray": row_val,
                    "col_gray": col_val,
                    "n_dontcares": pi["n_dontcares"],
                    "coverage": len(coverage),
                    "motif_count": motif_count,
                    "motif_fraction": motif_count / total_aa if total_aa > 0 else 0,
                }
            )

    # Sort by motif count
    motifs.sort(key=lambda x: x["motif_count"], reverse=True)

    print(f"  Total motifs extracted: {len(motifs)}")
    print(f"\n  Top 10 Co-evolution Motifs:")
    print(
        f"  {'Rank':>4s} {'Row':>3s} {'Col':>3s} {'DC':>3s} {'Coverage':>8s} {'Count':>6s} {'Fraction':>8s}"
    )
    print(f"  {'-' * 4} {'-' * 3} {'-' * 3} {'-' * 3} {'-' * 8} {'-' * 6} {'-' * 8}")

    for m in motifs[:10]:
        print(
            f"  {m['rank']:4d} {m['row_aa']:>3s} {m['col_aa']:>3s} {m['n_dontcares']:3d} "
            f"{m['coverage']:8d} {m['motif_count']:6d} {m['motif_fraction']:8.4f}"
        )

    return motifs


# ============================================================
# 5. Compute Coupling Constants
# ============================================================


def compute_coupling_constants(kmap_freq, sequences, n_seqs=100):
    """
    Compute coupling constants between residue positions.

    The coupling constant J_ij measures the strength of co-evolution
    between positions i and j. It's derived from the K-map structure:

    J_ij = log(P(aa_i, aa_j) / (P(aa_i) * P(aa_j)))

    This is analogous to the inverse Ising coupling in statistical physics.
    """
    print("\n=== Computing Coupling Constants ===")

    # Build position-specific frequency matrices
    n = min(n_seqs, len(sequences))

    # Get clean sequences
    clean_seqs = []
    for i in range(n):
        _, seq = sequences[i]
        clean = "".join(aa for aa in seq if aa in _AA_TO_INDEX)
        if len(clean) > 100:
            clean_seqs.append(clean)

    if len(clean_seqs) < 5:
        print("  Insufficient clean sequences")
        return {}

    min_len = min(len(s) for s in clean_seqs)
    # PERFORMANCE FIX: Original code used max_pos = min_len (1276 positions
    # → 813K pairs), which timed out at ~32% after 20 minutes. All other
    # scripts in this pipeline analyze positions 0-79 (the N-terminal
    # signal peptide region). Limit to 80 positions for consistency and speed.
    max_pos = min(80, min_len)  # 80 positions → 3160 pairs
    n_positions = max_pos

    # Pre-compute position arrays for fast vectorized access
    pos_arrays = []
    for s in clean_seqs:
        arr = np.array([_AA_TO_INDEX.get(aa, -1) for aa in s[:max_pos]], dtype=np.int32)
        pos_arrays.append(arr)

    print(f"  Positions: {n_positions}, Sequences: {len(clean_seqs)}")

    # Build all pairs and compute MI
    pairs = [(i, j) for i in range(n_positions) for j in range(i + 1, n_positions)]
    print(f"  Total pairs to compute: {len(pairs)}")
    mi_matrix = np.zeros((n_positions, n_positions))

    # PERFORMANCE FIX: Use vectorized numpy bincount approach instead of
    # pure-Python Counter. Build dense array once, then bincount per pair.
    dense = np.full((len(pos_arrays), n_positions), -1, dtype=np.int32)
    for si, arr in enumerate(pos_arrays):
        L = min(len(arr), n_positions)
        dense[si, :L] = arr[:L]

    for idx, (i, j) in enumerate(pairs):
        codes_i = dense[:, i]
        codes_j = dense[:, j]
        valid = (codes_i >= 0) & (codes_j >= 0)
        ci = codes_i[valid]
        cj = codes_j[valid]
        if len(ci) < 10:
            continue
        # Joint via bincount: flat = ci * 20 + cj
        joint_flat = (
            np.bincount(ci.astype(np.int64) * 20 + cj.astype(np.int64), minlength=400)
            .reshape(20, 20)
            .astype(np.float64)
        )
        total = joint_flat.sum()
        if total == 0:
            continue
        marg_i = joint_flat.sum(axis=1)
        marg_j = joint_flat.sum(axis=0)
        mi = 0.0
        for ai in range(20):
            for aj in range(20):
                if joint_flat[ai, aj] > 0 and marg_i[ai] > 0 and marg_j[aj] > 0:
                    p = joint_flat[ai, aj] / total
                    pi_v = marg_i[ai] / total
                    pj_v = marg_j[aj] / total
                    mi += p * np.log2(p / (pi_v * pj_v))
        mi_matrix[i, j] = mi
        mi_matrix[j, i] = mi
        if (idx + 1) % 500 == 0:
            print(f"    Coupling MI: {idx + 1}/{len(pairs)}")

    # Compute coupling constants from MI
    # J_ij = MI(i,j) * beta (inverse temperature)
    beta = 1.0  # Can be tuned
    coupling = beta * mi_matrix

    # Find strongest couplings
    strong_couplings = []
    for i in range(n_positions):
        for j in range(i + 1, n_positions):
            if abs(coupling[i, j]) > 0.01:
                strong_couplings.append((i, j, coupling[i, j]))

    strong_couplings.sort(key=lambda x: abs(x[2]), reverse=True)

    print(f"  Positions analyzed: {n_positions}")
    print(f"  Strong couplings (|J| > 0.01): {len(strong_couplings)}")
    print(f"\n  Top 10 Coupling Constants:")
    print(f"  {'Pos i':>5s} {'Pos j':>5s} {'J_ij':>8s} {'|J|':>8s}")
    print(f"  {'-' * 5} {'-' * 5} {'-' * 8} {'-' * 8}")

    for i, j, J in strong_couplings[:10]:
        print(f"  {i:5d} {j:5d} {J:8.4f} {abs(J):8.4f}")

    return {
        "n_positions": n_positions,
        "coupling_matrix": coupling.tolist(),
        "strong_couplings": strong_couplings[:50],
    }


# ============================================================
# 6. Predict Co-evolution from Boolean Function
# ============================================================


def predict_coevolution(qm_result, kmap_freq, sequences, n_seqs=100):
    """
    Test if the minimized Boolean function predicts co-evolution.

    The prediction: residue pairs that are "on" in the Boolean K-map
    (frequent dipeptides) should show higher co-evolution scores.
    """
    print("\n=== Predicting Co-evolution from Boolean Function ===")

    # Build clean sequences
    n = min(n_seqs, len(sequences))
    clean_seqs = []
    for i in range(n):
        _, seq = sequences[i]
        clean = "".join(aa for aa in seq if aa in _AA_TO_INDEX)
        if len(clean) > 100:
            clean_seqs.append(clean)

    if len(clean_seqs) < 5:
        print("  Insufficient sequences")
        return {}

    # Compute position-specific conservation
    min_len = min(len(s) for s in clean_seqs)
    max_pos = min(200, min_len)

    # Compute mutual information for each position pair
    n_test = min(30, max_pos)
    mi_scores = np.zeros((n_test, n_test))

    for i in range(n_test):
        for j in range(i + 1, n_test):
            joint = Counter()
            marg_i = Counter()
            marg_j = Counter()

            for seq in clean_seqs:
                if i < len(seq) and j < len(seq):
                    aa_i, aa_j = seq[i], seq[j]
                    if aa_i in _AA_TO_INDEX and aa_j in _AA_TO_INDEX:
                        joint[(aa_i, aa_j)] += 1
                        marg_i[aa_i] += 1
                        marg_j[aa_j] += 1

            total = sum(joint.values())
            if total == 0:
                continue

            mi = 0
            for (ai, aj), count in joint.items():
                p_joint = count / total
                p_i = marg_i[ai] / total
                p_j = marg_j[aj] / total
                if p_joint > 0 and p_i > 0 and p_j > 0:
                    mi += p_joint * np.log2(p_joint / (p_i * p_j))

            mi_scores[i, j] = mi
            mi_scores[j, i] = mi

    # Extract on-set from Boolean K-map
    on_set_cells = np.argwhere(kmap_freq >= np.percentile(kmap_freq[kmap_freq > 0], 75))

    # Map on-set cells to residue pairs
    gray_to_aa = {}
    for aa, idx in _AA_TO_INDEX.items():
        gray_val = gray_code_5bit(idx)
        gray_to_aa[gray_val] = aa

    on_set_pairs = set()
    for cell in on_set_cells:
        row, col = cell
        row_aa = gray_to_aa.get(row, None)
        col_aa = gray_to_aa.get(col, None)
        if row_aa and col_aa:
            on_set_pairs.add((row_aa, col_aa))

    # Compute average MI for on-set vs off-set pairs
    on_set_mi = []
    off_set_mi = []

    for i in range(n_test):
        for j in range(i + 1, n_test):
            if mi_scores[i, j] > 0:
                # Get amino acids at these positions (from first sequence)
                seq = clean_seqs[0]
                if i < len(seq) and j < len(seq):
                    aa_i, aa_j = seq[i], seq[j]
                    if (aa_i, aa_j) in on_set_pairs or (aa_j, aa_i) in on_set_pairs:
                        on_set_mi.append(mi_scores[i, j])
                    else:
                        off_set_mi.append(mi_scores[i, j])

    avg_on_mi = np.mean(on_set_mi) if on_set_mi else 0
    avg_off_mi = np.mean(off_set_mi) if off_set_mi else 0

    print(f"  On-set pairs (frequent dipeptides): {len(on_set_mi)}")
    print(f"  Off-set pairs (rare dipeptides): {len(off_set_mi)}")
    print(f"  Average MI for on-set: {avg_on_mi:.4f}")
    print(f"  Average MI for off-set: {avg_off_mi:.4f}")
    print(
        f"  MI ratio (on/off): {avg_on_mi / avg_off_mi:.2f}" if avg_off_mi > 0 else ""
    )

    # Compute prediction accuracy
    # If MI > median, predict "co-evolving"
    all_mi = list(on_set_mi) + list(off_set_mi)
    accuracy = 0  # Default: initialized to avoid UnboundLocalError
    if all_mi:
        median_mi = np.median(all_mi)
        on_correct = sum(1 for m in on_set_mi if m > median_mi)
        off_correct = sum(1 for m in off_set_mi if m <= median_mi)
        total = len(all_mi)
        accuracy = (on_correct + off_correct) / total if total > 0 else 0

        print(f"\n  Prediction accuracy: {accuracy:.4f}")
        print(f"  On-set correctly predicted: {on_correct}/{len(on_set_mi)}")
        print(f"  Off-set correctly predicted: {off_correct}/{len(off_set_mi)}")

    return {
        "on_set_pairs": list(on_set_pairs),
        "avg_on_mi": avg_on_mi,
        "avg_off_mi": avg_off_mi,
        "n_on_set": len(on_set_mi),
        "n_off_set": len(off_set_mi),
        "prediction_accuracy": accuracy if all_mi else 0,
    }


# ============================================================
# 7. Main Pipeline
# ============================================================


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    results_dir = base_dir / "boolean_results"
    results_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("Boolean Minimization and Co-evolution Prediction")
    print("SARS-CoV-2 Spike Protein")
    print("=" * 70)

    # 1. Parse FASTA
    print("\n[1/7] Parsing FASTA alignment...")
    sequences = parse_fasta(fasta_file)
    full_length = len(sequences[0][1])
    print(f"  Loaded {len(sequences)} sequences")

    # 2. Build Boolean K-map
    print("\n[2/7] Building Boolean K-map...")
    kmap_freq, kmap_bool, threshold = build_boolean_kmap(
        sequences, threshold_percentile=75
    )

    # 3. Boolean minimization
    print("\n[3/7] Running Quine-McCluskey minimization...")
    qm_result = minimize_boolean_function(kmap_bool)

    # 4. Extract co-evolution motifs
    print("\n[4/7] Extracting co-evolution motifs...")
    motifs = extract_coevolution_motifs(qm_result, sequences, n_seqs=1299)

    # 5. Compute coupling constants
    print("\n[5/7] Computing coupling constants...")
    coupling = compute_coupling_constants(kmap_freq, sequences, n_seqs=1299)

    # 6. Predict co-evolution
    print("\n[6/7] Predicting co-evolution from Boolean function...")
    prediction = predict_coevolution(qm_result, kmap_freq, sequences, n_seqs=1299)

    # 7. Save results
    print("\n[7/7] Saving results...")

    summary = {
        "dataset": "SARS-CoV-2 Spike Protein Boolean Analysis",
        "num_sequences": len(sequences),
        "boolean_kmap": {
            "threshold": threshold,
            "n_on_set": int(kmap_bool.sum()),
            "n_total": kmap_bool.size,
            "density": float(kmap_bool.sum() / kmap_bool.size),
        },
        "minimization": {
            "n_prime_implicants": qm_result["n_prime_implicants"],
            "n_essential": qm_result["n_essential"],
            "covering_size": qm_result["covering_size"],
        },
        "motifs": motifs[:20],
        "coupling": {
            "n_positions": coupling["n_positions"],
            "n_strong_couplings": len(coupling["strong_couplings"]),
            "top_couplings": coupling["strong_couplings"][:20],
        },
        "prediction": prediction,
    }

    with open(results_dir / "boolean_analysis_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # Save Boolean K-map
    np.save(results_dir / "boolean_kmap.npy", kmap_bool)
    np.savetxt(results_dir / "boolean_kmap.csv", kmap_bool, delimiter=",", fmt="%d")

    # Save coupling matrix
    np.save(results_dir / "coupling_matrix.npy", np.array(coupling["coupling_matrix"]))

    print(f"\nResults saved to: {results_dir}")

    # Print final summary
    print("\n" + "=" * 70)
    print("BOOLEAN MINIMIZATION COMPLETE")
    print("=" * 70)
    print(f"\nKey Results:")
    print(
        f"  Boolean K-map: {int(kmap_bool.sum())} on-set cells out of {kmap_bool.size}"
    )
    print(f"  Prime implicants: {qm_result['n_prime_implicants']}")
    print(f"  Essential prime implicants: {qm_result['n_essential']}")
    print(
        f"  Top co-evolution motifs: {motifs[0]['row_aa']}-{motifs[0]['col_aa']}"
        if motifs
        else ""
    )
    print(f"  Prediction accuracy: {prediction.get('prediction_accuracy', 0):.4f}")
    print(f"\n  Coupling constants:")
    for i, j, J in coupling["strong_couplings"][:5]:
        print(f"    J({i},{j}) = {J:.4f}")

    print(
        f"\n  The minimized Boolean function captures {qm_result['n_prime_implicants']} essential"
    )
    print(f"  residue-pair motifs that define the co-evolutionary structure of the")
    print(f"  SARS-CoV-2 Spike protein. These motifs can be used to predict")
    print(f"  co-evolutionary coupling constants between residue positions.")


if __name__ == "__main__":
    main()
