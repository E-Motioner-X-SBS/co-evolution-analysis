#!/usr/bin/env python3
"""
Advanced Co-evolution Analysis
==============================

Implements:
1. Co-evolution network (graph of co-evolving positions)
2. Walsh-Hadamard spectrum of co-evolution signals
3. Variant classification from co-evolution signatures
4. Clustering of sequences by co-evolution patterns
"""

import sys
import json
import numpy as np
from collections import Counter
from pathlib import Path

sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")

from nkmap.encoding.bio_sequences import Base20AminoEncoder, AMINO_HE_2012

# Shared co-evolution module
sys.path.insert(0, str(Path(__file__).resolve().parent))
from coevolution_shared import mutual_information, compute_entropy, majority_ref


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


# (Removed duplicate mi_vectorized and compute_mi_vectorized — now in coevolution_shared.mutual_information)
# Use shared coevolution_shared.mutual_information() for MI computation


def compute_mi_for_pairs(pos_arrays, variable_positions, n_seqs, max_gap=30):
    """Compute MI for all variable position pairs using shared mutual_information.
    Serially loops — correct for the 24-core system where parallelism comes from
    run_all_bg.sh distributing scripts, not within-script thread pools.
    """
    pairs = []
    for idx_i, pos_i in enumerate(variable_positions):
        for idx_j in range(idx_i + 1, len(variable_positions)):
            pos_j = variable_positions[idx_j]
            if abs(pos_i - pos_j) <= max_gap:
                pairs.append((pos_i, pos_j))

    print(
        f"  Computing MI for {len(pairs)} pairs (numpy vectorized via coevolution_shared)..."
    )

    mi_results = []
    for pi, pj in pairs:
        mi = mutual_information(pos_arrays, pi, pj, n_seqs)
        if mi > 0.01:
            mi_results.append((pi, pj, mi, 0))

    mi_results.sort(key=lambda x: x[2], reverse=True)
    return mi_results


def fwht(a):
    """Fast Walsh-Hadamard Transform (in-place, iterative).

    Operates on arrays whose length is a power of 2.
    Uses the standard butterfly algorithm: O(n log n) time, O(1) extra space.
    """
    a = np.asarray(a, dtype=np.float64).copy()
    n = len(a)
    if n & (n - 1) != 0:
        # Pad to next power of 2
        n2 = 1 << (n - 1).bit_length()
        a = np.concatenate([a, np.zeros(n2 - n)])
        n = n2
    h = 1
    while h < n:
        for i in range(0, n, h * 2):
            for j in range(i, i + h):
                x = a[j]
                y = a[j + h]
                a[j] = x + y
                a[j + h] = x - y
        h *= 2
    return a


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    results_dir = base_dir / "advanced_analysis_results"
    results_dir.mkdir(exist_ok=True)

    print("=" * 80)
    print("Advanced Co-evolution Analysis")
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

    # ============================================================
    # ANALYSIS 1: Co-evolution Network
    # ============================================================
    print("\n" + "=" * 80)
    print("ANALYSIS 1: Co-evolution Network")
    print("=" * 80)

    # Build MI matrix for all variable position pairs (PARALLEL)
    print("\n  Computing MI matrix for variable positions (PARALLEL)...")
    n_var = len(variable_positions)
    mi_matrix = np.zeros((n_var, n_var), dtype=np.float64)

    mi_results = compute_mi_for_pairs(pos_arrays, variable_positions, n_all, max_gap=30)

    for pos_i, pos_j, mi, total in mi_results:
        idx_i = variable_positions.index(pos_i)
        idx_j = variable_positions.index(pos_j)
        mi_matrix[idx_i, idx_j] = mi
        mi_matrix[idx_j, idx_i] = mi

    print(f"  Completed: {len(mi_results)} pairs with MI > 0.01")

    # Build adjacency matrix (edges with MI > threshold)
    threshold = 0.5
    adj_matrix = (mi_matrix > threshold).astype(int)
    np.fill_diagonal(adj_matrix, 0)

    # Network properties
    degree = adj_matrix.sum(axis=1)
    n_edges = int(adj_matrix.sum() // 2)
    avg_degree = np.mean(degree)
    max_degree_pos = variable_positions[np.argmax(degree)]

    print(f"  Nodes (variable positions): {n_var}")
    print(f"  Edges (MI > {threshold}): {n_edges}")
    print(f"  Average degree: {avg_degree:.2f}")
    print(f"  Highest degree position: {max_degree_pos} (degree={int(np.max(degree))})")

    # Find clusters (positions that co-evolve together)
    print("\n  Finding co-evolution clusters...")
    # Simple clustering: positions with MI > 1.0 form clusters
    high_mi_adj = (mi_matrix > 1.0).astype(int)
    np.fill_diagonal(high_mi_adj, 0)

    # Find connected components using BFS
    visited = set()
    clusters = []
    for start in range(n_var):
        if start in visited:
            continue
        cluster = []
        queue = [start]
        while queue:
            node = queue.pop(0)
            if node in visited:
                continue
            visited.add(node)
            cluster.append(variable_positions[node])
            for neighbor in range(n_var):
                if high_mi_adj[node, neighbor] and neighbor not in visited:
                    queue.append(neighbor)
        if len(cluster) > 1:
            clusters.append(sorted(cluster))

    print(f"  Co-evolution clusters: {len(clusters)}")
    for i, cluster in enumerate(clusters):
        print(f"    Cluster {i + 1}: positions {cluster}")

    # Save network
    network_data = {
        "nodes": [
            {
                "position": int(p),
                "degree": int(degree[i]),
                "entropy": float(compute_entropy(pos_arrays, p, n_all)),
                "perplexity": float(2 ** compute_entropy(pos_arrays, p, n_all)),
            }
            for i, p in enumerate(variable_positions)
        ],
        "edges": [
            {
                "source": int(variable_positions[i]),
                "target": int(variable_positions[j]),
                "mi": float(mi_matrix[i, j]),
            }
            for i in range(n_var)
            for j in range(i + 1, n_var)
            if mi_matrix[i, j] > threshold
        ],
        "clusters": clusters,
        "properties": {
            "n_nodes": n_var,
            "n_edges": n_edges,
            "avg_degree": float(avg_degree),
            "threshold": threshold,
        },
    }
    with open(results_dir / "coevolution_network.json", "w") as f:
        json.dump(network_data, f, indent=2)

    # ============================================================
    # ANALYSIS 2: Walsh-Hadamard Spectrum
    # ============================================================
    print("\n" + "=" * 80)
    print("ANALYSIS 2: Walsh-Hadamard Spectrum")
    print("=" * 80)

    # Build consensus K-map for co-evolutionary region
    print("\n  Building consensus K-map for positions 68-79...")
    # BUG FIX: The original code built a 32×32 K-map but then re-encoded
    # aa_list[ci] back to its code via encoder.encode.get() — this was a
    # no-op (returned ci). The K-map was effectively a 20×20 matrix padded
    # with zeros in a 32×32 grid.
    # FIX: Build a proper 20×20 base-20 K-map (He 2012 encoding, no Gray).
    consensus_kmap = np.zeros((20, 20), dtype=np.float64)
    for arr in pos_arrays[:n_all]:
        for pos in range(68, min(79, max_pos)):
            if pos < len(arr) and pos + 1 < len(arr):
                ci, cj = int(arr[pos]), int(arr[pos + 1])
                if 0 <= ci < 20 and 0 <= cj < 20:
                    consensus_kmap[ci, cj] += 1

    total = consensus_kmap.sum()
    if total > 0:
        consensus_kmap /= total

    # Walsh-Hadamard Transform
    # BUG FIX: The original code computed FWHT on only the first 32 elements
    # of a flattened 1024-element K-map — meaningless (only first row).
    # FIX: Compute 2D Walsh-Hadamard Transform on the full 20×20 K-map.
    # Since 20 is not a power of 2, we pad to 32×32 (next power of 2) with
    # zeros, apply 1D FWHT to each row, then to each column (separable 2D WHT).
    print("  Computing 2D Walsh-Hadamard Transform on full K-map...")
    padded = np.zeros((32, 32), dtype=np.float64)
    padded[:20, :20] = consensus_kmap
    # Apply FWHT to rows then columns (2D separable transform)
    spectrum_2d = np.zeros_like(padded)
    for i in range(32):
        spectrum_2d[i, :] = fwht(padded[i, :])
    for j in range(32):
        spectrum_2d[:, j] = fwht(spectrum_2d[:, j])
    spectrum_flat = np.abs(spectrum_2d).flatten()
    spectrum = spectrum_flat / (spectrum_flat.sum() if spectrum_flat.sum() > 0 else 1)

    # Find dominant modes
    dominant_modes = np.argsort(spectrum)[::-1][:5]
    explained_var = spectrum[dominant_modes[:3]].sum() * 100

    print(f"  Top 5 Walsh-Hadamard modes:")
    for mode in dominant_modes:
        print(f"    Mode {mode}: {spectrum[mode]:.4f}")
    print(f"  Top 3 modes explain: {explained_var:.1f}% of variance")

    # Save
    wh_data = {
        "spectrum": spectrum.tolist(),
        "dominant_modes": dominant_modes.tolist(),
        "explained_variance_top3": float(explained_var),
    }
    with open(results_dir / "walsh_hadamard_spectrum.json", "w") as f:
        json.dump(wh_data, f, indent=2)

    # ============================================================
    # ANALYSIS 3: Variant Classification
    # ============================================================
    print("\n" + "=" * 80)
    print("ANALYSIS 3: Variant Classification from Co-evolution Signatures")
    print("=" * 80)

    # Build co-evolution signature for each sequence
    print("\n  Building co-evolution signatures...")
    signatures = []
    for arr in pos_arrays[:n_all]:
        sig = []
        for pos_i, pos_j in [(76, 77), (74, 79), (71, 75), (72, 75), (78, 79)]:
            if pos_i < len(arr) and pos_j < len(arr):
                ci, cj = int(arr[pos_i]), int(arr[pos_j])
                sig.append(ci * 20 + cj)  # Encode as single integer
            else:
                sig.append(0)
        signatures.append(sig)

    signatures = np.array(signatures)

    # Cluster sequences by signature
    print("  Clustering sequences by co-evolution signature...")
    unique_signatures = {}
    for i, sig in enumerate(signatures):
        key = tuple(sig)
        if key not in unique_signatures:
            unique_signatures[key] = []
        unique_signatures[key].append(i)

    print(f"  Unique co-evolution signatures: {len(unique_signatures)}")
    print(f"  Top 5 signature clusters:")
    sorted_clusters = sorted(
        unique_signatures.items(), key=lambda x: len(x[1]), reverse=True
    )
    for sig, indices in sorted_clusters[:5]:
        print(f"    Signature {sig}: {len(indices)} sequences")

    # Save classification
    classification = {
        "n_unique_signatures": len(unique_signatures),
        "clusters": [
            {
                "signature": [int(x) for x in sig],
                "count": len(indices),
                "indices": [int(x) for x in indices[:10]],
            }
            for sig, indices in sorted_clusters[:10]
        ],
    }
    with open(results_dir / "variant_classification.json", "w") as f:
        json.dump(classification, f, indent=2)

    # ============================================================
    # ANALYSIS 4: Sequence Clustering by Co-evolution Patterns
    # ============================================================
    print("\n" + "=" * 80)
    print("ANALYSIS 4: Sequence Clustering")
    print("=" * 80)

    # Build feature vectors from co-evolution patterns
    print("\n  Building feature vectors...")
    features = []
    for arr in pos_arrays[:n_all]:
        feat = []
        for pos in range(68, min(80, max_pos)):
            if pos < len(arr):
                feat.append(int(arr[pos]))
            else:
                feat.append(-1)
        features.append(feat)

    features = np.array(features)

    # Compute pairwise Hamming distances
    print("  Computing pairwise distances...")
    n_seq = min(200, n_all)  # Limit for efficiency
    dist_matrix = np.zeros((n_seq, n_seq), dtype=np.float64)

    for i in range(n_seq):
        for j in range(i + 1, n_seq):
            valid = (features[i] >= 0) & (features[j] >= 0)
            if valid.sum() > 0:
                dist = np.sum(features[i, valid] != features[j, valid]) / valid.sum()
            else:
                dist = 1.0
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist

    # Simple clustering using distance threshold
    threshold_dist = 0.1
    clusters_seq = []
    assigned = set()
    for i in range(n_seq):
        if i in assigned:
            continue
        cluster = [i]
        assigned.add(i)
        for j in range(i + 1, n_seq):
            if j not in assigned and dist_matrix[i, j] < threshold_dist:
                cluster.append(j)
                assigned.add(j)
        clusters_seq.append(cluster)

    print(f"  Sequence clusters (distance < {threshold_dist}): {len(clusters_seq)}")
    for i, cluster in enumerate(clusters_seq[:5]):
        print(f"    Cluster {i + 1}: {len(cluster)} sequences")

    # Save
    cluster_data = {
        "n_clusters": len(clusters_seq),
        "threshold": threshold_dist,
        "clusters": [{"size": len(c), "indices": c[:20]} for c in clusters_seq[:10]],
    }
    with open(results_dir / "sequence_clusters.json", "w") as f:
        json.dump(cluster_data, f, indent=2)

    # Summary
    print("\n" + "=" * 80)
    print("ADVANCED ANALYSIS SUMMARY")
    print("=" * 80)
    print(f"Sequences: {n_all}")
    print(f"Variable positions: {len(variable_positions)}")
    print(f"Co-evolution network: {n_var} nodes, {n_edges} edges")
    print(f"Co-evolution clusters: {len(clusters)}")
    print(f"Unique co-evolution signatures: {len(unique_signatures)}")
    print(f"Sequence clusters: {len(clusters_seq)}")
    print(f"Top WH mode explains: {explained_var:.1f}% variance")
    print(f"\nResults saved to: {results_dir}")


if __name__ == "__main__":
    main()
