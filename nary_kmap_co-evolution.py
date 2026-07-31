#!/usr/bin/env python3
"""
N-ary K-map Analysis of SARS-CoV-2 Spike Protein Co-evolution
=============================================================

Uses base-20 encoding (direct amino acid mapping) instead of binary
5-bit Gray code. This gives a 20x20 K-map for dipeptides.

Key difference from binary approach:
- Binary: 32x32 = 1024 cells (12 don't-care cells)
- N-ary: 20x20 = 400 cells (0 don't-care cells, all used)

The n-ary approach is more compact and directly captures the
biochemical relationships between amino acids.
"""

import os
import sys
import json
import numpy as np
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")

from nkmap.encoding.bio_sequences import (
    Base20AminoEncoder,
    AMINO_HE_2012,
    AMINO_BY_GROUP,
)
from nkmap.analysis.kmap_builder import (
    count_kmers_cpu,
    normalize_kmap,
    build_kmap_2d,
)
from nkmap.encoding.n_ary_gray import n_ary_hamming


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


# ============================================================
# 2. N-ary K-map Construction
# ============================================================


def build_nary_kmap(sequences, k=2, n_seqs=None):
    """
    Build a base-20 n-ary K-map for the alignment.

    For k=2 (dipeptides): 20x20 = 400 cells
    For k=3 (tripeptides): 20^3 = 8000 cells
    """
    encoder = Base20AminoEncoder(version=1)  # He 2012 ordering

    if n_seqs is None:
        n_seqs = len(sequences)
    n_seqs = min(n_seqs, len(sequences))

    # Aggregate k-mer counts across all sequences
    cell_size = 20**k
    total_counts = np.zeros(cell_size, dtype=np.int64)
    total_kmers = 0

    for i in range(n_seqs):
        _, seq = sequences[i]
        # Remove gaps and non-canonical residues
        clean_seq = "".join(aa for aa in seq if aa in encoder.encode)

        counts = count_kmers_cpu(clean_seq, k, encoder)
        total_counts += counts
        total_kmers += max(0, len(clean_seq) - k + 1)

    # Normalize to frequencies
    freq = normalize_kmap(total_counts, mode="frequency")

    print(f"\n=== N-ary K-map (base-20, k={k}) ===")
    print(f"  Sequences analyzed: {n_seqs}")
    print(f"  Total k-mers: {total_kmers}")
    print(f"  K-map size: {cell_size} cells ({20**k} = 20^{k})")
    print(f"  Non-zero cells: {np.count_nonzero(freq)}")
    print(f"  Density: {np.count_nonzero(freq) / cell_size:.4f}")

    return freq, total_counts, encoder


def build_nary_kmap_2d(freq_1d, k=2, encoder=None):
    """Reshape 1D K-map to 2D for visualization."""
    if k == 2:
        return freq_1d.reshape(20, 20)
    else:
        return build_kmap_2d(freq_1d, k, encoder)


# ============================================================
# 3. Boolean K-map from N-ary Frequencies
# ============================================================


def build_boolean_kmap_nary(freq_1d, threshold_percentile=75):
    """
    Build a Boolean K-map from n-ary frequencies.

    Threshold the frequency K-map to create a Boolean function
    f(cell) = 1 if frequency >= threshold.
    """
    # Threshold
    nonzero = freq_1d[freq_1d > 0]
    threshold = np.percentile(nonzero, threshold_percentile) if len(nonzero) > 0 else 0

    kmap_bool = (freq_1d >= threshold).astype(int)

    n_on = int(kmap_bool.sum())
    n_total = kmap_bool.size

    print(f"\n=== Boolean K-map (n-ary) ===")
    print(f"  Threshold (P{threshold_percentile}): {threshold:.6f}")
    print(f"  On-set (frequent): {n_on} cells")
    print(f"  Off-set (rare): {n_total - n_on} cells")
    print(f"  Density: {n_on / n_total:.4f}")

    return kmap_bool, threshold


# ============================================================
# 4. Quine-McCluskey Minimization on N-ary K-map
# ============================================================


def minimize_nary_boolean(kmap_bool):
    """
    Run Quine-McCluskey on the n-ary Boolean K-map.

    For a 20x20 map, we have 2*log2(20) ≈ 8.6 variables.
    The QM algorithm works on the flattened truth table.
    """
    print("\n=== Quine-McCluskey on N-ary Boolean K-map ===")

    # Import from kmap-sbm-validation
    sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/kmap-sbm-validation/src")
    from kmap_sbm.analysis.prime_implicants import (
        kmap_truth_table,
        prime_implicants_quine_mccluskey,
        boolean_minimize_kmap,
    )

    # Flatten for QM
    bool_flat = kmap_bool.flatten().astype(int)

    # Run minimization
    result = boolean_minimize_kmap(bool_flat, algorithm="qm")

    print(f"  Input variables: {result['n_vars']}")
    print(f"  Minterms (on-set): {result['n_minterms']}")
    print(f"  Prime implicants: {result['n_prime_implicants']}")
    print(f"  Essential prime implicants: {result['n_essential']}")
    print(f"  Covering size: {result['covering_size']}")

    return result


# ============================================================
# 5. Co-evolution Motifs from N-ary Prime Implicants
# ============================================================


def extract_nary_motifs(qm_result, encoder, sequences, n_seqs=100):
    """
    Extract co-evolution motifs from the minimized n-ary Boolean function.
    Optimized: pre-compute encoded arrays, use numpy for fast matching.
    """
    print("\n=== Extracting N-ary Co-evolution Motifs ===")

    decode = encoder.decode
    n = min(n_seqs, len(sequences))
    
    # Pre-encode all sequences into flat arrays for fast access
    print(f"  Pre-encoding {n} sequences...")
    encoded_seqs = []
    for _, seq in sequences[:n]:
        clean = [encoder.encode[aa] for aa in seq if aa in encoder.encode]
        if len(clean) > 1:
            encoded_seqs.append(np.array(clean, dtype=np.int32))
    print(f"  Encoded {len(encoded_seqs)} sequences")
    
    motifs = []
    for i, pi in enumerate(qm_result["prime_implicants"]):
        values = pi["values"]
        mask = pi["mask"]
        coverage = pi["coverage"]
        n_vars = len(values)
        n_digits = n_vars // 2

        # Extract row and column codes from prime implicant
        row_code = sum(values[j] * (2 ** (min(n_digits, 5) - 1 - j)) 
                       for j in range(min(n_digits, 5)) if not mask[j])
        col_code = sum(values[j + n_digits] * (2 ** (min(n_digits, 5) - 1 - j))
                       for j in range(min(n_digits, 5)) if not mask[j + n_digits])

        row_aa = decode.get(row_code % 20, "?")
        col_aa = decode.get(col_code % 20, "?")

        # Count occurrences using numpy vectorization
        motif_count = 0
        for arr in encoded_seqs:
            if len(arr) < 2:
                continue
            # Find all positions where consecutive pair matches (row_aa, col_aa)
            matches_first = (arr[:-1] % 20 == row_code % 20)
            matches_second = (arr[1:] % 20 == col_code % 20)
            motif_count += int(np.sum(matches_first & matches_second))

        if motif_count > 0:
            motifs.append({
                "rank": i + 1, "row_aa": row_aa, "col_aa": col_aa,
                "row_code": row_code % 20, "col_code": col_code % 20,
                "n_dontcares": pi["n_dontcares"], "coverage": len(coverage),
                "motif_count": motif_count,
            })

    motifs.sort(key=lambda x: x["motif_count"], reverse=True)
    print(f"  Total motifs: {len(motifs)}")
    if motifs:
        print(f"\n  Top 5 Co-evolution Motifs:")
        for m in motifs[:5]:
            print(f"    {m['rank']:4d} {m['row_aa']:>3s}-{m['col_aa']:>3s} count={m['motif_count']:6d}")
    return motifs


# ============================================================
# 6. Coupling Constants from N-ary K-map
# ============================================================


def compute_nary_couplings(freq_2d, sequences, encoder, n_seqs=100):
    """
    Compute coupling constants from the n-ary K-map structure.

    J_ij = log(P(aa_i, aa_j) / (P(aa_i) * P(aa_j)))

    This is the inverse Ising coupling, directly from the 20x20 frequency matrix.
    """
    print("\n=== Computing N-ary Coupling Constants ===")

    # Marginal frequencies
    row_marginal = freq_2d.sum(axis=1)  # P(aa_i)
    col_marginal = freq_2d.sum(axis=0)  # P(aa_j)

    # Compute coupling constants
    # J_ij = log(P(aa_i, aa_j) / (P(aa_i) * P(aa_j)))
    # Avoid log(0) by adding small epsilon
    epsilon = 1e-10

    with np.errstate(divide="ignore", invalid="ignore"):
        coupling = np.log(
            (freq_2d + epsilon)
            / ((row_marginal[:, None] + epsilon) * (col_marginal[None, :] + epsilon))
        )

    # Find strongest couplings
    strong_couplings = []
    aa_list = list(AMINO_HE_2012)

    for i in range(20):
        for j in range(20):
            if i != j and abs(coupling[i, j]) > 0.1:
                strong_couplings.append(
                    {
                        "aa_i": aa_list[i],
                        "aa_j": aa_list[j],
                        "code_i": i,
                        "code_j": j,
                        "J": float(coupling[i, j]),
                        "freq": float(freq_2d[i, j]),
                    }
                )

    strong_couplings.sort(key=lambda x: abs(x["J"]), reverse=True)

    print(f"  K-map size: 20x20 = 400 cells")
    print(f"  Strong couplings (|J| > 0.1): {len(strong_couplings)}")
    print(f"\n  Top 15 Coupling Constants:")
    print(f"  {'AA i':>5s} {'AA j':>5s} {'Code':>8s} {'J_ij':>8s} {'Freq':>8s}")
    print(f"  {'-' * 5} {'-' * 5} {'-' * 8} {'-' * 8} {'-' * 8}")

    for c in strong_couplings[:15]:
        print(
            f"  {c['aa_i']:>5s} {c['aa_j']:>5s} {c['code_i']:2d},{c['code_j']:2d}  "
            f"{c['J']:8.4f} {c['freq']:8.6f}"
        )

    return {
        "coupling_matrix": coupling.tolist(),
        "strong_couplings": strong_couplings,
        "row_marginal": row_marginal.tolist(),
        "col_marginal": col_marginal.tolist(),
    }


# ============================================================
# 7. Co-evolution Prediction Test
# ============================================================


def predict_coevolution_nary(freq_2d, encoder, sequences, n_seqs=100):
    """
    Test if the n-ary K-map predicts co-evolution.

    The prediction: dipeptide pairs that are "on" in the Boolean K-map
    should show higher mutual information between their positions.
    """
    print("\n=== Predicting Co-evolution from N-ary K-map ===")

    # Build clean sequences
    n = min(n_seqs, len(sequences))
    clean_seqs = []
    for i in range(n):
        _, seq = sequences[i]
        clean = "".join(aa for aa in seq if aa in encoder.encode)
        if len(clean) > 100:
            clean_seqs.append(clean)

    if len(clean_seqs) < 5:
        print("  Insufficient sequences")
        return {}

    # Compute pairwise mutual information
    min_len = min(len(s) for s in clean_seqs)
    n_test = min(30, min_len)

    mi_matrix = np.zeros((n_test, n_test))

    for i in range(n_test):
        for j in range(i + 1, n_test):
            joint = Counter()
            marg_i = Counter()
            marg_j = Counter()

            for seq in clean_seqs:
                if i < len(seq) and j < len(seq):
                    aa_i, aa_j = seq[i], seq[j]
                    if aa_i in encoder.encode and aa_j in encoder.encode:
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

            mi_matrix[i, j] = mi
            mi_matrix[j, i] = mi

    # Classify pairs as "on-set" (frequent dipeptide) vs "off-set"
    # Use the n-ary frequency as threshold
    threshold = np.percentile(freq_2d[freq_2d > 0], 75)

    on_set_pairs = set()
    for i in range(20):
        for j in range(20):
            if freq_2d[i, j] >= threshold:
                aa_i = AMINO_HE_2012[i]
                aa_j = AMINO_HE_2012[j]
                on_set_pairs.add((aa_i, aa_j))

    # Compute MI for on-set vs off-set
    on_set_mi = []
    off_set_mi = []

    for i in range(n_test):
        for j in range(i + 1, n_test):
            if mi_matrix[i, j] > 0:
                seq = clean_seqs[0]
                if i < len(seq) and j < len(seq):
                    aa_i, aa_j = seq[i], seq[j]
                    if (aa_i, aa_j) in on_set_pairs or (aa_j, aa_i) in on_set_pairs:
                        on_set_mi.append(mi_matrix[i, j])
                    else:
                        off_set_mi.append(mi_matrix[i, j])

    avg_on_mi = np.mean(on_set_mi) if on_set_mi else 0
    avg_off_mi = np.mean(off_set_mi) if off_set_mi else 0

    print(f"  On-set pairs (frequent): {len(on_set_mi)}")
    print(f"  Off-set pairs (rare): {len(off_set_mi)}")
    print(f"  Average MI for on-set: {avg_on_mi:.4f}")
    print(f"  Average MI for off-set: {avg_off_mi:.4f}")
    if avg_off_mi > 0:
        print(f"  MI ratio (on/off): {avg_on_mi / avg_off_mi:.2f}")

    return {
        "n_on_set": len(on_set_mi),
        "n_off_set": len(off_set_mi),
        "avg_on_mi": avg_on_mi,
        "avg_off_mi": avg_off_mi,
        "mi_ratio": avg_on_mi / avg_off_mi if avg_off_mi > 0 else 0,
    }


# ============================================================
# 8. Main Pipeline
# ============================================================


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    results_dir = base_dir / "nary_kmap_results"
    results_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("N-ary K-map Analysis of SARS-CoV-2 Spike Protein Co-evolution")
    print("=" * 70)

    # 1. Parse FASTA
    print("\n[1/6] Parsing FASTA alignment...")
    sequences = parse_fasta(fasta_file)
    print(f"  Loaded {len(sequences)} sequences")

    # 2. Build n-ary K-map (k=2)
    print("\n[2/6] Building base-20 n-ary K-map (k=2)...")
    freq_1d, counts, encoder = build_nary_kmap(sequences, k=2, n_seqs=1299)
    freq_2d = build_nary_kmap_2d(freq_1d, k=2, encoder=encoder)

    # 3. Boolean minimization
    print("\n[3/6] Building Boolean K-map and minimizing...")
    kmap_bool, threshold = build_boolean_kmap_nary(freq_1d, threshold_percentile=75)
    qm_result = minimize_nary_boolean(kmap_bool)

    # 4. Extract motifs
    print("\n[4/6] Extracting co-evolution motifs...")
    motifs = extract_nary_motifs(qm_result, encoder, sequences, n_seqs=1299)

    # 5. Compute couplings
    print("\n[5/6] Computing coupling constants...")
    coupling = compute_nary_couplings(freq_2d, sequences, encoder, n_seqs=1299)

    # 6. Predict co-evolution
    print("\n[6/6] Predicting co-evolution...")
    prediction = predict_coevolution_nary(freq_2d, encoder, sequences, n_seqs=1299)

    # Save results
    print("\n" + "=" * 70)
    print("Saving Results...")
    print("=" * 70)

    summary = {
        "dataset": "SARS-CoV-2 Spike Protein N-ary K-map Analysis",
        "encoding": "Base-20 (He 2012 order)",
        "kmap_size": "20x20 = 400 cells",
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
            "n_strong_couplings": len(coupling["strong_couplings"]),
            "top_couplings": coupling["strong_couplings"][:20],
        },
        "prediction": prediction,
    }

    with open(results_dir / "nary_analysis_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # Save K-maps
    np.save(results_dir / "nary_kmap_freq.npy", freq_1d)
    np.save(results_dir / "nary_kmap_2d.npy", freq_2d)
    np.savetxt(results_dir / "nary_kmap_2d.csv", freq_2d, delimiter=",", fmt="%.6f")
    np.save(results_dir / "nary_boolean_kmap.npy", kmap_bool)
    np.save(
        results_dir / "nary_coupling_matrix.npy", np.array(coupling["coupling_matrix"])
    )

    print(f"\nResults saved to: {results_dir}")

    # Final summary
    print("\n" + "=" * 70)
    print("N-ARY K-MAP ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nKey Results:")
    print(f"  Encoding: Base-20 (He 2012 order)")
    print(f"  K-map size: 20x20 = 400 cells (no don't-care)")
    print(f"  Boolean on-set: {int(kmap_bool.sum())} cells")
    print(f"  Prime implicants: {qm_result['n_prime_implicants']}")
    print(f"  Essential prime implicants: {qm_result['n_essential']}")
    print(f"  Top motif: {motifs[0]['row_aa']}-{motifs[0]['col_aa']}" if motifs else "")
    print(f"  Strong couplings: {len(coupling['strong_couplings'])}")
    print(f"  MI ratio (on/off): {prediction.get('mi_ratio', 0):.2f}")

    print(f"\n  The n-ary K-map captures the co-evolutionary structure of the")
    print(f"  Spike protein using direct base-20 encoding. The 400-cell K-map")
    print(f"  has no don't-care cells, providing a complete picture of all")
    print(f"  dipeptide relationships. The prime implicants identify the essential")
    print(f"  co-evolutionary motifs, and the coupling constants quantify the")
    print(f"  strength of co-evolution between residue pairs.")


if __name__ == "__main__":
    main()
