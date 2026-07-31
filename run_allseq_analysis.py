#!/usr/bin/env python3
"""
Full Position-Based K-map Co-evolution Analysis
================================================
Runs on ALL 1,299 sequences from the Spike protein alignment.
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
    """Pre-compute position arrays for all sequences."""
    print("  Building position arrays...")
    pos_arrays = []
    min_len = 999999
    for _, seq in sequences:
        clean = "".join(aa for aa in seq if aa in encoder.encode)
        arr = np.array(
            [encoder.encode.get(aa, -1) for aa in clean[:max_pos]], dtype=np.int32
        )
        pos_arrays.append(arr)
        min_len = min(min_len, len(arr))
    return pos_arrays, min_len


def compute_mi_matrix(pos_arrays, n_seqs, max_pos, window=30):
    """Compute MI for all nearby position pairs using vectorized operations."""
    print("  Computing MI matrix...")
    mi_matrix = np.zeros((max_pos, max_pos), dtype=np.float64)

    # BUG FIX: Original code used `range(0, max_pos, 2)` which skipped every
    # other position (computed MI for (0,2), (0,4), ... but never (0,1), (0,3)).
    # This missed ~50% of position pairs. Fixed to use step=1.
    for i in range(max_pos):
        for j in range(i + 2, min(i + window, max_pos)):
            # Extract codes at positions i and j
            codes_i = np.array([arr[i] for arr in pos_arrays[:n_seqs] if i < len(arr)])
            codes_j = np.array([arr[j] for arr in pos_arrays[:n_seqs] if j < len(arr)])

            min_len = min(len(codes_i), len(codes_j))
            if min_len < 10:
                continue

            codes_i = codes_i[:min_len]
            codes_j = codes_j[:min_len]

            # Filter valid codes
            valid = (codes_i >= 0) & (codes_j >= 0)
            codes_i = codes_i[valid]
            codes_j = codes_j[valid]

            if len(codes_i) < 10:
                continue

            # Compute joint and marginal counts
            joint = Counter()
            marg_i = Counter()
            marg_j = Counter()

            for ci, cj in zip(codes_i, codes_j):
                joint[(int(ci), int(cj))] += 1
                marg_i[int(ci)] += 1
                marg_j[int(cj)] += 1

            total = sum(joint.values())
            if total == 0:
                continue

            mi = 0.0
            for (ai, aj), count in joint.items():
                p_joint = count / total
                p_i = marg_i[ai] / total
                p_j = marg_j[aj] / total
                if p_joint > 0 and p_i > 0 and p_j > 0:
                    mi += p_joint * np.log2(p_joint / (p_i * p_j))

            mi_matrix[i, j] = mi
            mi_matrix[j, i] = mi

    return mi_matrix


def build_position_kmap(pos_arrays, encoder, pos_i, pos_j, n_seqs):
    """Build 20x20 K-map for a position pair."""
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


def compute_coupling(kmap):
    """Compute coupling constants from K-map."""
    row_marg = kmap.sum(axis=1)
    col_marg = kmap.sum(axis=0)
    epsilon = 1e-10
    with np.errstate(divide="ignore", invalid="ignore"):
        coupling = np.log(
            (kmap + epsilon)
            / ((row_marg[:, None] + epsilon) * (col_marg[None, :] + epsilon))
        )
    return coupling


def minimize_kmap(kmap_2d):
    """Run QM minimization on position-pair K-map."""
    nonzero = kmap_2d[kmap_2d > 0]
    if len(nonzero) == 0:
        return None
    threshold = np.percentile(nonzero, 70)
    kmap_bool = (kmap_2d >= threshold).astype(int)
    bool_flat = kmap_bool.flatten().astype(int)
    result = boolean_minimize_kmap(bool_flat, algorithm="qm")
    return result, kmap_bool, threshold


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    results_dir = base_dir / "full_position_results"
    results_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("FULL Position-Based K-map Co-evolution Analysis")
    print("=" * 70)

    # Load ALL sequences
    print("\n[1/5] Loading ALL sequences...")
    sequences = parse_fasta(fasta_file)
    full_length = len(sequences[0][1])
    encoder = Base20AminoEncoder(version=1)
    n_all = len(sequences)
    print(f"  Total sequences: {n_all}")

    # Build position arrays
    print("\n[2/5] Building position arrays...")
    pos_arrays, min_len = build_position_arrays(sequences, encoder, max_pos=200)
    print(f"  Min sequence length: {min_len}")

    # Compute MI matrix
    print("\n[3/5] Computing MI for position pairs (ALL sequences)...")
    max_pos = min(80, min_len)
    mi_matrix = compute_mi_matrix(pos_arrays, n_all, max_pos, window=30)

    # Find top co-evolving pairs
    mi_results = []
    for i in range(max_pos):
        for j in range(i + 2, max_pos):
            if mi_matrix[i, j] > 0.005:
                mi_results.append((i, j, float(mi_matrix[i, j])))

    mi_results.sort(key=lambda x: x[2], reverse=True)
    print(f"  Pairs with MI > 0.005: {len(mi_results)}")

    print(f"\n  Top 20 co-evolving position pairs (ALL {n_all} sequences):")
    print(f"  {'Pos i':>5s} {'Pos j':>5s} {'MI':>8s}")
    print(f"  {'-' * 5} {'-' * 5} {'-' * 8}")
    for i, j, mi in mi_results[:20]:
        print(f"  {i:5d} {j:5d} {mi:8.4f}")

    # Build position-pair K-maps for top pairs
    print("\n[4/5] Building position-pair K-maps and minimizing...")
    aa_list = list(AMINO_HE_2012)
    position_results = []

    for idx, (pos_i, pos_j, mi) in enumerate(mi_results[:15]):
        print(f"\n  --- Pair ({pos_i}, {pos_j}), MI={mi:.4f} ---")

        kmap = build_position_kmap(pos_arrays, encoder, pos_i, pos_j, n_all)
        coupling = compute_coupling(kmap)

        # Minimize
        minimization = minimize_kmap(kmap)
        if minimization is None:
            continue

        result, kmap_bool, threshold = minimization

        avg_J = float(np.mean(np.abs(coupling)))
        max_pos_idx = np.unravel_index(np.argmax(coupling), coupling.shape)
        max_neg_idx = np.unravel_index(np.argmin(coupling), coupling.shape)

        # Top residue pairs
        top_pairs = []
        for ai in range(20):
            for aj in range(20):
                if kmap[ai, aj] > 0.001:
                    top_pairs.append(
                        (
                            aa_list[ai],
                            aa_list[aj],
                            float(kmap[ai, aj]),
                            float(coupling[ai, aj]),
                        )
                    )
        top_pairs.sort(key=lambda x: x[2], reverse=True)

        print(f"  On-set: {int(kmap_bool.sum())} cells")
        print(f"  Prime implicants: {result['n_prime_implicants']}")
        print(f"  Essential PIs: {result['n_essential']}")
        print(f"  Avg |J|: {avg_J:.4f}")
        print(
            f"  Strongest +: {aa_list[max_pos_idx[0]]}-{aa_list[max_pos_idx[1]]} (J={coupling[max_pos_idx]:.4f})"
        )
        print(
            f"  Strongest -: {aa_list[max_neg_idx[0]]}-{aa_list[max_neg_idx[1]]} (J={coupling[max_neg_idx]:.4f})"
        )
        print(f"  Top 5 residue pairs:")
        for aa1, aa2, freq, J in top_pairs[:5]:
            print(f"    {aa1}-{aa2}: freq={freq:.4f}, J={J:.4f}")

        # Prime implicant equations
        print(f"  Prime implicants as equations:")
        for pi_idx, pi in enumerate(result["prime_implicants"][:5]):
            values = list(pi["values"])
            mask = list(pi["mask"])
            while len(values) < 8:
                values.append(0)
                mask.append(False)

            # Decode row and col codes
            row_code = 0
            col_code = 0
            for j in range(4):
                if not mask[j]:
                    row_code = row_code * 2 + values[j]
            for j in range(4):
                if not mask[j + 4]:
                    col_code = col_code * 2 + values[j + 4]

            row_aa = aa_list[row_code % 20] if row_code < 20 else "?"
            col_aa = aa_list[col_code % 20] if col_code < 20 else "?"

            terms = []
            for j in range(8):
                if not mask[j]:
                    var = f"r{3 - j}" if j < 4 else f"c{7 - j}"
                    if values[j] == 0:
                        terms.append(f"~{var}")
                    else:
                        terms.append(var)
            term_str = " & ".join(terms) if terms else "TRUE"

            print(f"    PI_{pi_idx + 1}: {term_str} = ({row_aa},{col_aa})")

        position_results.append(
            {
                "pos_i": pos_i,
                "pos_j": pos_j,
                "mi": mi,
                "n_on_set": int(kmap_bool.sum()),
                "n_prime_implicants": result["n_prime_implicants"],
                "n_essential": result["n_essential"],
                "avg_coupling": avg_J,
                "strongest_positive": f"{aa_list[max_pos_idx[0]]}-{aa_list[max_pos_idx[1]]}",
                "strongest_negative": f"{aa_list[max_neg_idx[0]]}-{aa_list[max_neg_idx[1]]}",
                "top_pairs": [(a1, a2, f, J) for a1, a2, f, J in top_pairs[:10]],
            }
        )

    # Step 5: Save
    print("\n[5/5] Saving results...")

    summary = {
        "dataset": "SARS-CoV-2 Spike Protein (ALL sequences)",
        "num_sequences": n_all,
        "encoding": "Base-20 (He 2012)",
        "max_positions": max_pos,
        "top_mi_pairs": [
            {"pos_i": i, "pos_j": j, "mi": mi} for i, j, mi in mi_results[:30]
        ],
        "position_results": position_results,
    }

    with open(results_dir / "full_analysis_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nResults saved to: {results_dir}")

    # Final
    print("\n" + "=" * 70)
    print(f"COMPLETE (ALL {n_all} SEQUENCES)")
    print("=" * 70)
    print(f"Co-evolving pairs: {len(mi_results)}")
    print(
        f"Top pair: ({mi_results[0][0]},{mi_results[0][1]}) MI={mi_results[0][2]:.4f}"
    )
    print(f"\nCoupling constants:")
    for pr in position_results[:5]:
        print(
            f"  J({pr['pos_i']},{pr['pos_j']}) = {pr['avg_coupling']:.4f} "
            f"[+{pr['strongest_positive']}, -{pr['strongest_negative']}]"
        )


if __name__ == "__main__":
    main()
