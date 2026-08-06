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


def build_position_arrays(sequences, encoder, max_pos=None):
    """Pre-compute position arrays for all sequences.

    max_pos=None → FULL sequence length (all positions).
    """
    print("  Building position arrays...")
    pos_arrays = []
    min_len = 999999
    for _, seq in sequences:
        clean = "".join(aa for aa in seq if aa in encoder.encode)
        if max_pos is None:
            arr = np.array([encoder.encode.get(aa, -1) for aa in clean], dtype=np.int32)
        else:
            arr = np.array(
                [encoder.encode.get(aa, -1) for aa in clean[:max_pos]], dtype=np.int32
            )
        pos_arrays.append(arr)
        min_len = min(min_len, len(arr))
    return pos_arrays, min_len


def compute_mi_matrix(pos_arrays, n_seqs, max_pos, window=30):
    """Compute MI for all nearby position pairs using GPU (torch CUDA).

    Falls back to vectorized numpy if CUDA unavailable.
    """
    print("  Computing MI matrix...")
    mi_matrix = np.zeros((max_pos, max_pos), dtype=np.float64)

    # Build all pairs within window (i, i+1 .. i+window) — step=1 (no skipping)
    pairs = [
        (i, j) for i in range(max_pos) for j in range(i + 1, min(i + window, max_pos))
    ]

    # GPU path (torch CUDA on A100)
    try:
        import coevolution_gpu as cg

        dense = cg.dense_to_gpu(pos_arrays)
        mi_dict, _ = cg.mi_matrix_gpu(dense, pairs, min_total=10, chunk=16384)
        for (i, j), mi in mi_dict.items():
            mi_matrix[i, j] = mi
            mi_matrix[j, i] = mi
        print(f"  GPU MI computed for {len(mi_dict)} pairs (torch CUDA)")
        return mi_matrix
    except Exception as e:
        print(f"  GPU failed ({e}), using numpy fallback...")

    # NumPy fallback (vectorized, no Counter)
    dense = np.full((n_seqs, max_pos), -1, dtype=np.int32)
    for si, arr in enumerate(pos_arrays[:n_seqs]):
        L = min(len(arr), max_pos)
        dense[si, :L] = arr[:L]

    for idx, (i, j) in enumerate(pairs):
        codes_i = dense[:, i]
        codes_j = dense[:, j]
        valid = (codes_i >= 0) & (codes_j >= 0)
        ci = codes_i[valid]
        cj = codes_j[valid]
        if len(ci) < 10:
            continue
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

    # Build position arrays — FULL LENGTH (all 1276 positions, all sequences)
    print("\n[2/5] Building position arrays...")
    pos_arrays, min_len = build_position_arrays(sequences, encoder, max_pos=None)
    print(f"  Min sequence length: {min_len}")

    # Compute MI matrix — FULL LENGTH on GPU
    print("\n[3/5] Computing MI for position pairs (ALL sequences, FULL length)...")
    max_pos = min_len  # full length — NO truncation
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
