#!/usr/bin/env python3
"""
Local Precision Matrix → Boolean Minimization Pipeline
========================================================

WARNING: This script was originally named "DCA" but does NOT implement
real Direct Coupling Analysis (DCA). Real DCA (Morcos et al 2011, PNAS)
inverts the GLOBAL 19L×19L covariance matrix to disentangle direct from
indirect correlations. This script instead inverts per-pair 20×20 covariance
matrices, which computes the LOCAL PRECISION MATRIX — NOT DCA couplings.

The local precision matrix does NOT remove transitive correlations, which
is the defining feature of DCA. Results from this script should NOT be
interpreted as DCA results.

Pipeline:
1. Build per-pair 20×20 covariance matrices C_ij
2. Compute local precision J_ij = C_ij^{-1} (per-pair pseudoinverse)
3. Threshold J to create Boolean function
4. Run Quine-McCluskey minimization
5. Test prediction accuracy

For proper DCA, use pydca (https://github.com/KIT-MBS/pydca) or
EVcouplings (https://github.com/debbiemarkslab/EVcouplings).
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


def build_position_arrays(sequences, encoder, max_pos=200):
    """Pre-compute position arrays."""
    pos_arrays = []
    for _, seq in sequences:
        clean = "".join(aa for aa in seq if aa in encoder.encode)
        arr = np.array(
            [encoder.encode.get(aa, -1) for aa in clean[:max_pos]], dtype=np.int32
        )
        pos_arrays.append(arr)
    return pos_arrays


def compute_covariance_matrix(pos_arrays, pos_i, pos_j, n_seqs):
    """
    Compute covariance matrix C_ij(sigma_i, sigma_j).

    C_ij = P(sigma_i, sigma_j) - P(sigma_i) * P(sigma_j)

    This removes the conservation signal and captures
    only the pairwise correlations.
    """
    # Joint distribution
    joint = np.zeros((20, 20), dtype=np.float64)
    for arr in pos_arrays[:n_seqs]:
        if pos_i < len(arr) and pos_j < len(arr):
            ci, cj = int(arr[pos_i]), int(arr[pos_j])
            if ci >= 0 and cj >= 0:
                joint[ci, cj] += 1

    total = joint.sum()
    if total == 0:
        return None, None, None

    joint /= total

    # Marginal distributions
    marg_i = joint.sum(axis=1)
    marg_j = joint.sum(axis=0)

    # Covariance matrix: C = P(x,y) - P(x)*P(y)
    C = joint - marg_i[:, None] * marg_j[None, :]

    return C, joint, marg_i


def compute_coupling_constants(C):
    """
    Compute coupling constants J = C^{-1} (pseudoinverse).

    The inverse covariance matrix captures the DIRECT couplings
    between residues that are NOT explained by single-position conservation.
    """
    # Add small regularization to make matrix invertible
    epsilon = 0.01 * np.eye(20)
    C_reg = C + epsilon

    # Compute pseudoinverse (more robust than direct inverse)
    try:
        J = np.linalg.pinv(C_reg)
    except np.linalg.LinAlgError:
        J = np.zeros_like(C)

    return J


def threshold_to_boolean(J, threshold_percentile=75):
    """
    Threshold coupling matrix to Boolean function.

    f(sigma_i, sigma_j) = 1 if |J_ij| > threshold
    """
    abs_J = np.abs(J)
    nonzero = abs_J[abs_J > 0]

    if len(nonzero) == 0:
        return np.zeros((20, 20), dtype=int), 0

    threshold = np.percentile(nonzero, threshold_percentile)
    J_bool = (abs_J >= threshold).astype(int)

    return J_bool, threshold


def minimize_boolean(J_bool):
    """Run Quine-McCluskey on Boolean coupling matrix."""
    bool_flat = J_bool.flatten().astype(int)
    result = boolean_minimize_kmap(bool_flat, algorithm="qm")
    return result


def predict_from_J(J, pos_i, pos_j, pos_arrays, n_seqs, train_end=800):
    """
    Predict co-evolution using J matrix.

    Given a residue at position i, predict the most likely
    residue at position j based on coupling constants.
    """
    # Training: build J from first train_end sequences
    C_train, joint_train, marg_train = compute_covariance_matrix(
        pos_arrays, pos_i, pos_j, train_end
    )

    if C_train is None:
        return 0, 0

    J_train = compute_coupling_constants(C_train)

    # Testing: predict on remaining sequences
    test_arrays = pos_arrays[train_end:]
    correct = 0
    total = 0

    for arr in test_arrays:
        if pos_i < len(arr) and pos_j < len(arr):
            ci, cj = int(arr[pos_i]), int(arr[pos_j])
            if ci >= 0 and cj >= 0:
                # Given residue ci at pos_i, predict best cj at pos_j
                predicted_j = np.argmax(J_train[ci, :])
                if predicted_j == cj:
                    correct += 1
                total += 1

    return correct, total


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    results_dir = base_dir / "dca_boolean_results"
    results_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("DCA → Boolean Minimization Pipeline")
    print("=" * 70)

    # Load ALL sequences
    print("\n[1/6] Loading ALL sequences...")
    sequences = parse_fasta(fasta_file)
    full_length = len(sequences[0][1])
    encoder = Base20AminoEncoder(version=1)
    aa_list = list(AMINO_HE_2012)
    n_all = len(sequences)
    print(f"  Total: {n_all} sequences")

    # Build position arrays
    print("\n[2/6] Building position arrays...")
    pos_arrays = build_position_arrays(sequences, encoder, max_pos=full_length)
    print(f"  Max position: {full_length} (FULL LENGTH)")

    # Top co-evolving pairs — dynamically computed from full-length MI
    # BUG FIX: was hardcoded pairs in positions 68-78. Now computed via
    # GPU mutation-only MI over ALL variable positions (full length).
    print("\n  Finding co-evolving pairs (GPU, full length)...")
    try:
        from coevolution_shared import find_coevolving_pairs_gpu
        from collections import Counter as _C

        # variable positions (entropy > 0.3)
        vp = []
        for p in range(full_length):
            cnt = _C(int(a[p]) for a in pos_arrays if p < len(a) and a[p] >= 0)
            t = sum(cnt.values())
            if t == 0:
                continue
            h = -sum((c / t) * np.log2(c / t) for c in cnt.values() if c > 0)
            if h > 0.3:
                vp.append(p)
        co_evolving_pairs = find_coevolving_pairs_gpu(
            pos_arrays, vp, n_all, max_gap=30, min_mi=0.1
        )
        co_evolving_pairs = [(p[0], p[1]) for p in co_evolving_pairs[:10]]
        print(f"  Top 10 co-evolving pairs (dynamic): {co_evolving_pairs}")
    except Exception as e:
        print(f"  GPU pair finding failed ({e}), using fallback pairs")
        co_evolving_pairs = [
            (74, 76),
            (74, 77),
            (72, 74),
            (76, 78),
            (74, 78),
            (68, 74),
            (72, 76),
            (68, 73),
            (68, 77),
            (72, 75),
        ]

    # ============================================================
    # STEP 1-2: Compute J_ij via DCA for all pairs
    # ============================================================
    print("\n[3/6] Computing J_ij via DCA (inverse covariance)...")

    all_results = []

    for pos_i, pos_j in co_evolving_pairs:
        print(f"\n  --- Position pair ({pos_i}, {pos_j}) ---")

        # Step 1: Covariance matrix
        C, joint, marg = compute_covariance_matrix(pos_arrays, pos_i, pos_j, n_all)
        if C is None:
            continue

        print(f"  Covariance matrix: 20x20")
        print(f"  Max |C_ij|: {np.max(np.abs(C)):.6f}")
        print(f"  Mean |C_ij|: {np.mean(np.abs(C)):.6f}")

        # Step 2: Coupling constants via inverse
        J = compute_coupling_constants(C)

        print(f"  Coupling matrix J: 20x20")
        print(f"  Max |J_ij|: {np.max(np.abs(J)):.6f}")
        print(f"  Mean |J_ij|: {np.mean(np.abs(J)):.6f}")

        # Find strongest couplings
        abs_J = np.abs(J)
        np.fill_diagonal(abs_J, 0)  # Ignore self-coupling
        max_idx = np.unravel_index(np.argmax(abs_J), abs_J.shape)

        print(
            f"  Strongest coupling: {aa_list[max_idx[0]]}-{aa_list[max_idx[1]]} J={J[max_idx]:.4f}"
        )

        # Step 3: Threshold to Boolean
        J_bool, threshold = threshold_to_boolean(J, threshold_percentile=75)
        print(f"  Boolean threshold: {threshold:.4f}")
        print(f"  On-set: {J_bool.sum()} pairs")
        print(f"  Off-set: {400 - J_bool.sum()} pairs")

        # Step 4: QM minimization
        result = minimize_boolean(J_bool)
        print(f"  Prime implicants: {result['n_prime_implicants']}")
        print(f"  Essential PIs: {result['n_essential']}")

        # Decode essential PIs
        print(f"  Essential prime implicants:")
        for pi_idx, pi in enumerate(result["essential_prime_implicants"]):
            values = list(pi["values"])
            mask = list(pi["mask"])
            while len(values) < 8:
                values.append(0)
                mask.append(False)

            row_code = sum(values[j] * (2 ** (3 - j)) for j in range(4) if not mask[j])
            col_code = sum(
                values[j + 4] * (2 ** (3 - j)) for j in range(4) if not mask[j + 4]
            )

            row_aa = aa_list[row_code % 20] if row_code < 20 else "?"
            col_aa = aa_list[col_code % 20] if col_code < 20 else "?"

            terms = []
            for j in range(8):
                if not mask[j]:
                    var = f"s{3 - j}" if j < 4 else f"t{7 - j}"
                    if values[j] == 0:
                        terms.append(f"~{var}")
                    else:
                        terms.append(var)
            term_str = " & ".join(terms) if terms else "TRUE"

            print(f"    PI_{pi_idx + 1}: {term_str} = ({row_aa}, {col_aa})")

        # Step 5: Test prediction
        print(f"\n  Testing prediction...")
        correct, total = predict_from_J(
            J, pos_i, pos_j, pos_arrays, n_all, train_end=800
        )
        accuracy = correct / total if total > 0 else 0

        print(f"  Prediction accuracy: {accuracy:.4f} ({correct}/{total})")

        all_results.append(
            {
                "pos_i": pos_i,
                "pos_j": pos_j,
                "max_coupling": float(np.max(np.abs(J))),
                "mean_coupling": float(np.mean(np.abs(J))),
                "n_on_set": int(J_bool.sum()),
                "n_prime_implicants": result["n_prime_implicants"],
                "n_essential": result["n_essential"],
                "accuracy": accuracy,
            }
        )

    # ============================================================
    # STEP 6: Summary
    # ============================================================
    print("\n[4/6] Summary...")

    avg_accuracy = np.mean([r["accuracy"] for r in all_results])

    print(f"\n{'=' * 70}")
    print(f"DCA → BOOLEAN PIPELINE RESULTS")
    print(f"{'=' * 70}")
    print(f"\nSequences: {n_all}")
    print(f"Position pairs analyzed: {len(all_results)}")
    print(f"Average prediction accuracy: {avg_accuracy:.4f}")
    print(f"\nPer-pair results:")
    for r in all_results:
        print(
            f"  ({r['pos_i']:2d},{r['pos_j']:2d}): acc={r['accuracy']:.4f}, "
            f"PI={r['n_prime_implicants']}, EPI={r['n_essential']}"
        )

    # Save
    summary = {
        "dataset": "SARS-CoV-2 Spike Protein (ALL sequences)",
        "num_sequences": n_all,
        "method": "DCA (inverse covariance) → Boolean minimization",
        "results": all_results,
        "avg_accuracy": float(avg_accuracy),
    }

    with open(results_dir / "dca_boolean_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nResults saved to: {results_dir}")


if __name__ == "__main__":
    main()
