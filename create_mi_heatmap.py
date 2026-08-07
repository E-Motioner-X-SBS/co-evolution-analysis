#!/usr/bin/env python3
"""
Create MI Heatmap for Co-evolution Analysis
============================================
Generates a heatmap showing mutual information between all position pairs.
"""

import sys
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


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    output_dir = base_dir / "mi_heatmap"
    output_dir.mkdir(exist_ok=True)

    print("Loading ALL sequences...")
    sequences = parse_fasta(fasta_file)
    encoder = Base20AminoEncoder(version=1)
    n_all = len(sequences)
    full_length = len(sequences[0][1])
    print(f"Total: {n_all} sequences")

    # Build position arrays
    print("Building position arrays...")
    max_pos = len(sequences[0][1])
    pos_arrays_list = []
    for _, seq in sequences:
        arr = np.array(
            [encoder.encode.get(aa, 20) for aa in seq[:max_pos]], dtype=np.int32
        )
        pos_arrays_list.append(arr)
    pos_arrays = np.array(pos_arrays_list)

    # Compute MI matrix for ALL position pairs (FULL length) — GPU accelerated
    print("Computing MI matrix for all position pairs (FULL length, GPU)...")
    mi_matrix = np.zeros((max_pos, max_pos), dtype=np.float64)

    try:
        import coevolution_gpu as cg

        pos_arrays_list = [arr for arr in pos_arrays]  # list of 1D arrays
        dense = cg.dense_to_gpu(pos_arrays_list)
        pairs = cg.all_pairs(max_pos)  # ALL pairs, full matrix (~813K)
        print(f"  GPU: {len(pairs)} pairs on {dense.device}")
        mi_dict, _ = cg.mi_matrix_gpu(dense, pairs, min_total=10, chunk=32768)
        for (i, j), mi in mi_dict.items():
            mi_matrix[i, j] = mi
            mi_matrix[j, i] = mi
        print(f"  Completed: {len(mi_dict)} pairs (GPU)")
    except Exception as e:
        print(f"  GPU failed ({e}), using vectorized numpy...")

        def mi_vectorized(pos_arrays, pos_i, pos_j, n_seqs):
            codes_i = pos_arrays[:n_seqs, pos_i]
            codes_j = pos_arrays[:n_seqs, pos_j]
            valid = (0 <= codes_i < 20) & (0 <= codes_j < 20)
            codes_i = codes_i[valid]
            codes_j = codes_j[valid]
            if len(codes_i) < 10:
                return 0.0
            pairs = codes_i * 20 + codes_j
            joint_flat = (
                np.bincount(pairs, minlength=400).reshape(20, 20).astype(np.float64)
            )
            total = joint_flat.sum()
            if total == 0:
                return 0.0
            marg_i = joint_flat.sum(axis=1)
            marg_j = joint_flat.sum(axis=0)
            mi = 0.0
            for ai in range(20):
                for aj in range(20):
                    if joint_flat[ai, aj] > 0 and marg_i[ai] > 0 and marg_j[aj] > 0:
                        p = joint_flat[ai, aj] / total
                        p_i = marg_i[ai] / total
                        p_j = marg_j[aj] / total
                        if p > 0 and p_i > 0 and p_j > 0:
                            mi += p * np.log2(p / (p_i * p_j))
            return mi

        total_pairs = 0
        for i in range(max_pos):
            for j in range(i + 1, max_pos):
                mi = mi_vectorized(pos_arrays, i, j, n_all)
                mi_matrix[i, j] = mi
                mi_matrix[j, i] = mi
                total_pairs += 1
                if total_pairs % 10000 == 0:
                    print(f"    Progress: {total_pairs} pairs")
        print(f"  Completed: {total_pairs} pairs")

    # Save MI matrix
    np.save(output_dir / "mi_matrix.npy", mi_matrix)
    np.savetxt(output_dir / "mi_matrix.csv", mi_matrix, delimiter=",", fmt="%.4f")

    # Create heatmap using matplotlib
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap

        print("Creating MI heatmap...")

        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))

        # Custom colormap: white for low MI, blue for medium, red for high
        colors = [
            "#f7fbff",
            "#deebf7",
            "#c6dbef",
            "#9ecae1",
            "#6baed6",
            "#4292c6",
            "#2171b5",
            "#084594",
            "#08306b",
        ]
        cmap = LinearSegmentedColormap.from_list("mi_heatmap", colors)

        # Plot heatmap
        im = ax.imshow(mi_matrix, cmap=cmap, aspect="auto", vmin=0, vmax=2)

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label("Mutual Information (bits)", fontsize=12)

        # Labels
        ax.set_xlabel("Position j", fontsize=12)
        ax.set_ylabel("Position i", fontsize=12)
        ax.set_title(
            f"Mutual Information Heatmap - SARS-CoV-2 Spike Protein\n({n_all} sequences, positions 0-{max_pos - 1})",
            fontsize=14,
        )

        # Tick marks
        tick_positions = np.arange(0, max_pos, 5)
        ax.set_xticks(tick_positions)
        ax.set_yticks(tick_positions)
        ax.set_xticklabels([str(x) for x in tick_positions], fontsize=8)
        ax.set_yticklabels([str(x) for x in tick_positions], fontsize=8)

        # Mark high-MI pairs
        high_mi_pairs = []
        for i in range(max_pos):
            for j in range(i + 1, max_pos):
                if mi_matrix[i, j] > 1.0:
                    high_mi_pairs.append((i, j, mi_matrix[i, j]))

        high_mi_pairs.sort(key=lambda x: x[2], reverse=True)

        for i, j, mi in high_mi_pairs[:20]:
            ax.plot(j, i, "r*", markersize=8)
            ax.plot(i, j, "r*", markersize=8)

        # Add legend for stars
        ax.plot([], [], "r*", markersize=8, label="Top MI pairs")
        ax.legend(loc="upper right", fontsize=10)

        plt.tight_layout()
        plt.savefig(output_dir / "mi_heatmap.png", dpi=150, bbox_inches="tight")
        plt.savefig(output_dir / "mi_heatmap.pdf", bbox_inches="tight")
        print(f"Heatmap saved to {output_dir / 'mi_heatmap.png'}")

        # Also create a focused heatmap for the co-evolutionary region (positions 60-80)
        fig2, ax2 = plt.subplots(figsize=(10, 8))
        focus_start = 60
        focus_end = 80
        mi_focus = mi_matrix[focus_start:focus_end, focus_start:focus_end]

        im2 = ax2.imshow(mi_focus, cmap="hot", aspect="auto", vmin=0, vmax=10)
        cbar2 = plt.colorbar(im2, ax=ax2, shrink=0.8)
        cbar2.set_label("Mutual Information (bits)", fontsize=12)

        ax2.set_xlabel("Position j", fontsize=12)
        ax2.set_ylabel("Position i", fontsize=12)
        ax2.set_title(
            f"MI Heatmap - Co-evolutionary Region (positions {focus_start}-{focus_end - 1})\n({n_all} sequences)",
            fontsize=14,
        )

        tick_pos = np.arange(0, focus_end - focus_start, 2)
        ax2.set_xticks(tick_pos)
        ax2.set_yticks(tick_pos)
        ax2.set_xticklabels([str(focus_start + x) for x in tick_pos], fontsize=9)
        ax2.set_yticklabels([str(focus_start + x) for x in tick_pos], fontsize=9)

        # Mark high-MI pairs in this region
        for i in range(focus_end - focus_start):
            for j in range(i + 1, focus_end - focus_start):
                if mi_focus[i, j] > 2.0:
                    ax2.plot(j, i, "w*", markersize=6)
                    ax2.plot(i, j, "w*", markersize=6)

        ax2.plot([], [], "w*", markersize=6, label="MI > 2.0")
        ax2.legend(loc="upper right", fontsize=9)

        plt.tight_layout()
        plt.savefig(output_dir / "mi_heatmap_focus.png", dpi=150, bbox_inches="tight")
        plt.savefig(output_dir / "mi_heatmap_focus.pdf", bbox_inches="tight")
        print(f"Focused heatmap saved to {output_dir / 'mi_heatmap_focus.png'}")

        plt.close("all")

    except ImportError:
        print("matplotlib not available - saving raw data only")

    # Save summary
    summary = {
        "n_sequences": n_all,
        "max_position": max_pos,
        "mi_matrix_shape": list(mi_matrix.shape),
        "max_mi": float(np.max(mi_matrix)),
        "mean_mi": float(np.mean(mi_matrix[mi_matrix > 0])),
        "high_mi_pairs_count": int(np.sum(mi_matrix > 1.0) // 2),
    }

    with open(output_dir / "mi_heatmap_summary.json", "w") as f:
        import json

        json.dump(summary, f, indent=2)

    print(f"\nSummary:")
    print(f"  Max MI: {summary['max_mi']:.4f}")
    print(f"  Mean MI (non-zero): {summary['mean_mi']:.4f}")
    print(f"  Pairs with MI > 1.0: {summary['high_mi_pairs_count']}")
    print(f"\nFiles saved to: {output_dir}")


if __name__ == "__main__":
    main()
