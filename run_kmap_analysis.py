#!/usr/bin/env python3
"""
K-map Analysis of SARS-CoV-2 Spike Protein Co-evolution
=======================================================

Applies the complete K-map biological encoding framework to the
Spike protein alignment dataset. Uses the 5-bit Gray code encoding
from AminoAcidEncoding.lean (Lean 4 verified).

Analysis pipeline:
1. Parse FASTA alignment
2. Encode sequences using 5-bit Gray code
3. Build K-maps (dipeptide frequency)
4. Run H1-H6 hypothesis tests
5. Generate co-evolution analysis
"""

import os
import sys
import json
import numpy as np
from collections import Counter, defaultdict
from pathlib import Path

# Add the kmap-sbm-validation src to path
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

# ============================================================
# 1. FASTA Parser
# ============================================================


def parse_fasta(filepath):
    """Parse a FASTA file and return list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_seq = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences.append((current_header, "".join(current_seq)))
                current_header = line[1:]  # Remove '>'
                current_seq = []
            elif line:
                current_seq.append(line.upper())

    if current_header is not None:
        sequences.append((current_header, "".join(current_seq)))

    return sequences


def extract_accession(header):
    """Extract accession number from FASTA header."""
    return header.split()[0]


def extract_description(header):
    """Extract description from FASTA header."""
    parts = header.split(" | ", 1)
    if len(parts) > 1:
        return parts[1]
    return header


# ============================================================
# 2. Sequence Statistics
# ============================================================


def compute_sequence_stats(sequences):
    """Compute basic statistics for the alignment."""
    seq_lengths = [len(seq) for _, seq in sequences]

    # Count amino acid frequencies across all sequences
    aa_counts = Counter()
    gap_counts = 0
    unknown_counts = 0

    for _, seq in sequences:
        for aa in seq:
            if aa in _AA_TO_INDEX:
                aa_counts[aa] += 1
            elif aa == "-":
                gap_counts += 1
            else:
                unknown_counts += 1

    total_residues = sum(aa_counts.values()) + gap_counts + unknown_counts

    return {
        "num_sequences": len(sequences),
        "alignment_length": seq_lengths[0] if seq_lengths else 0,
        "min_seq_length": min(seq_lengths) if seq_lengths else 0,
        "max_seq_length": max(seq_lengths) if seq_lengths else 0,
        "mean_seq_length": np.mean(seq_lengths) if seq_lengths else 0,
        "total_residues": total_residues,
        "aa_frequencies": {
            aa: count / sum(aa_counts.values()) for aa, count in aa_counts.items()
        },
        "gap_count": gap_counts,
        "unknown_count": unknown_counts,
    }


# ============================================================
# 3. Gray Code Encoding
# ============================================================


def encode_sequence_gray(seq):
    """Encode a sequence to 5-bit Gray code values.

    CORRECTED (FIX A1): gaps/unknown are encoded as -1 (kept in the array
    to preserve alignment positions), not stripped. Consecutive-pair
    analyses must skip pairs where either code is -1.
    """
    encoded = []
    for aa in seq:
        if aa in _AA_TO_INDEX:
            encoded.append(encode_gray_single(aa))
        else:
            encoded.append(-1)
    return encoded


def compute_hamming_distances(seq1_encoded, seq2_encoded):
    """Compute pairwise Hamming distances between two encoded sequences."""
    min_len = min(len(seq1_encoded), len(seq2_encoded))
    distances = []
    for i in range(min_len):
        d = gray_hamming_int(seq1_encoded[i], seq2_encoded[i])
        distances.append(d)
    return distances


def compute_hamming_matrix(encoded_seqs, n_seqs=None):
    """Compute pairwise Hamming distance matrix for a subset of sequences."""
    if n_seqs is None:
        n_seqs = len(encoded_seqs)
    n_seqs = min(n_seqs, len(encoded_seqs))

    ham_matrix = np.zeros((n_seqs, n_seqs), dtype=float)

    for i in range(n_seqs):
        for j in range(i + 1, n_seqs):
            min_len = min(len(encoded_seqs[i]), len(encoded_seqs[j]))
            total_dist = 0
            for k in range(min_len):
                total_dist += gray_hamming_int(encoded_seqs[i][k], encoded_seqs[j][k])
            avg_dist = total_dist / min_len if min_len > 0 else 0
            ham_matrix[i, j] = avg_dist
            ham_matrix[j, i] = avg_dist

    return ham_matrix


# ============================================================
# 4. K-map Construction
# ============================================================


def build_kmap_for_sequence(seq):
    """Build a 32x32 dipeptide frequency K-map for a single sequence."""
    return build_aa_kmap_2d(seq)


def build_kmap_consensus(sequences, n_seqs=None):
    """Build a consensus K-map averaged over multiple sequences."""
    if n_seqs is None:
        n_seqs = len(sequences)
    n_seqs = min(n_seqs, len(sequences))

    kmap_sum = np.zeros(1024, dtype=float)
    count = 0

    for i in range(n_seqs):
        _, seq = sequences[i]
        # Remove gaps
        clean_seq = "".join(aa for aa in seq if aa in _AA_TO_INDEX)
        if len(clean_seq) > 1:
            kmap = build_aa_kmap_2d(clean_seq)
            kmap_sum += kmap
            count += 1

    if count > 0:
        kmap_sum /= count

    return kmap_sum


# ============================================================
# 5. H1: Gray-code Adjacency Analysis
# ============================================================


def analyze_h1_adjacency(sequences, n_seqs=100):
    """
    H1: Gray-code adjacency enrichment.

    Tests whether consecutive residues in sequences tend to have
    Hamming distance 1 in Gray code space (adjacent on K-map).
    """
    print("\n=== H1: Gray-code Adjacency Analysis ===")

    total_pairs = 0
    hamming_1_count = 0
    hamming_dist_counts = Counter()

    n = min(n_seqs, len(sequences))

    for i in range(n):
        _, seq = sequences[i]
        encoded = encode_sequence_gray(seq)

        for j in range(len(encoded) - 1):
            # CORRECTED (FIX A1): only count pairs where BOTH residues are
            # canonical and truly consecutive in the alignment (no gap).
            if encoded[j] < 0 or encoded[j + 1] < 0:
                continue
            d = gray_hamming_int(encoded[j], encoded[j + 1])
            hamming_dist_counts[d] += 1
            total_pairs += 1
            if d == 1:
                hamming_1_count += 1

    ratio = hamming_1_count / total_pairs if total_pairs > 0 else 0

    # Expected ratio for random 5-bit encoding
    # In 5-bit space, probability of hamming distance 1 = 5/31 ≈ 0.161
    expected_ratio = 5 / 31

    print(f"  Total consecutive pairs: {total_pairs}")
    print(f"  Hamming distance 1 pairs: {hamming_1_count}")
    print(f"  Observed ratio: {ratio:.4f}")
    print(f"  Expected ratio (random): {expected_ratio:.4f}")
    print(f"  Enrichment ratio: {ratio / expected_ratio:.4f}")
    print(f"  Hamming distance distribution:")
    for d in sorted(hamming_dist_counts.keys()):
        print(
            f"    d={d}: {hamming_dist_counts[d]} ({100 * hamming_dist_counts[d] / total_pairs:.1f}%)"
        )

    return {
        "total_pairs": total_pairs,
        "hamming_1_count": hamming_1_count,
        "observed_ratio": ratio,
        "expected_ratio": expected_ratio,
        "enrichment_ratio": ratio / expected_ratio if expected_ratio > 0 else 0,
        "distribution": dict(hamming_dist_counts),
    }


# ============================================================
# 6. H2: K-map Signature vs Sequence Properties
# ============================================================


def analyze_h2_signature(sequences, n_seqs=100):
    """
    H2: K-map signature analysis.

    Tests whether the K-map signature correlates with sequence properties
    like length, GC content, and amino acid composition.
    """
    print("\n=== H2: K-map Signature Analysis ===")

    n = min(n_seqs, len(sequences))

    seq_lengths = []
    gc_contents = []
    entropy_values = []
    kmap_norms = []

    for i in range(n):
        _, seq = sequences[i]
        clean_seq = "".join(aa for aa in seq if aa in _AA_TO_INDEX)

        seq_lengths.append(len(clean_seq))

        # GC content (using nucleotide analogy - count G and C amino acids)
        gc = (
            (clean_seq.count("G") + clean_seq.count("C")) / len(clean_seq)
            if len(clean_seq) > 0
            else 0
        )
        gc_contents.append(gc)

        # Shannon entropy of amino acid composition
        aa_freq = Counter(clean_seq)
        total = sum(aa_freq.values())
        entropy = -sum(
            (c / total) * np.log2(c / total) for c in aa_freq.values() if c > 0
        )
        entropy_values.append(entropy)

        # K-map norm
        kmap = build_aa_kmap_2d(clean_seq)
        kmap_norms.append(np.linalg.norm(kmap))

    print(f"  Sequences analyzed: {n}")
    print(f"  Sequence length range: {min(seq_lengths)} - {max(seq_lengths)}")
    print(f"  Mean GC content: {np.mean(gc_contents):.4f}")
    print(f"  Mean entropy: {np.mean(entropy_values):.4f} bits")
    print(f"  Mean K-map norm: {np.mean(kmap_norms):.6f}")

    # Compute correlations
    corr_len_gc = np.corrcoef(seq_lengths, gc_contents)[0, 1] if n > 1 else 0
    corr_len_entropy = np.corrcoef(seq_lengths, entropy_values)[0, 1] if n > 1 else 0

    print(f"  Correlation (length vs GC): {corr_len_gc:.4f}")
    print(f"  Correlation (length vs entropy): {corr_len_entropy:.4f}")

    return {
        "num_sequences": n,
        "seq_lengths": seq_lengths,
        "gc_contents": gc_contents,
        "entropy_values": entropy_values,
        "kmap_norms": kmap_norms,
        "corr_length_gc": corr_len_gc,
        "corr_length_entropy": corr_len_entropy,
    }


# ============================================================
# 7. H3: Structural Classification from K-map Distances
# ============================================================


def analyze_h3_structure(sequences, n_seqs=50):
    """
    H3: K-map distance-based structural classification.

    Tests whether K-map distances can distinguish different
    sequence clusters (analogous to structural classification).
    """
    print("\n=== H3: Structural Classification from K-map Distances ===")

    n = min(n_seqs, len(sequences))

    # Build K-maps for subset
    kmaps = []
    seq_ids = []
    for i in range(n):
        header, seq = sequences[i]
        clean_seq = "".join(aa for aa in seq if aa in _AA_TO_INDEX)
        if len(clean_seq) > 10:
            kmap = build_aa_kmap_2d(clean_seq)
            kmaps.append(kmap)
            seq_ids.append(extract_accession(header))

    if len(kmaps) < 2:
        print("  Insufficient sequences for analysis")
        return {}

    # Compute pairwise K-map distances
    kmap_array = np.array(kmaps)
    n_kmaps = len(kmap_array)

    # Euclidean distance matrix
    from scipy.spatial.distance import pdist, squareform

    dist_matrix = squareform(pdist(kmap_array, metric="euclidean"))

    # Average distance
    upper_tri = dist_matrix[np.triu_indices(n_kmaps, k=1)]
    avg_dist = np.mean(upper_tri)
    std_dist = np.std(upper_tri)

    print(f"  Sequences analyzed: {n_kmaps}")
    print(f"  Average K-map distance: {avg_dist:.6f}")
    print(f"  Std K-map distance: {std_dist:.6f}")
    print(f"  Min distance: {np.min(upper_tri):.6f}")
    print(f"  Max distance: {np.max(upper_tri):.6f}")

    # Find most similar and most different pairs
    min_idx = np.unravel_index(
        np.argmin(dist_matrix + np.eye(n_kmaps) * 1e10), dist_matrix.shape
    )
    max_idx = np.unravel_index(np.argmax(dist_matrix), dist_matrix.shape)

    print(
        f"  Most similar pair: {seq_ids[min_idx[0]]} vs {seq_ids[min_idx[1]]} (dist={dist_matrix[min_idx]:.6f})"
    )
    print(
        f"  Most different pair: {seq_ids[max_idx[0]]} vs {seq_ids[max_idx[1]]} (dist={dist_matrix[max_idx]:.6f})"
    )

    return {
        "num_sequences": n_kmaps,
        "avg_distance": avg_dist,
        "std_distance": std_dist,
        "min_distance": float(np.min(upper_tri)),
        "max_distance": float(np.max(upper_tri)),
        "most_similar": (
            seq_ids[min_idx[0]],
            seq_ids[min_idx[1]],
            float(dist_matrix[min_idx]),
        ),
        "most_different": (
            seq_ids[max_idx[0]],
            seq_ids[max_idx[1]],
            float(dist_matrix[max_idx]),
        ),
    }


# ============================================================
# 8. H5: Contact Map Invariants
# ============================================================


def analyze_h5_contacts(sequences, n_seqs=100):
    """
    H5: Contact map invariants.

    Tests whether the K-map encoding satisfies contact map properties:
    - Symmetry: K-map should show symmetric dipeptide frequencies
    - Distance distribution: Hamming distances between encoded residues
    """
    print("\n=== H5: Contact Map Invariants ===")

    n = min(n_seqs, len(sequences))

    # Build consensus K-map
    kmap_consensus = build_kmap_consensus(sequences, n)

    # Check symmetry (K-map should be approximately symmetric for
    # dipeptide frequencies if the encoding is well-balanced)
    kmap_2d = kmap_consensus.reshape(32, 32)
    symmetry_error = np.mean(np.abs(kmap_2d - kmap_2d.T))

    # Compute Hamming distance distribution for all encoded AAs
    all_encoded = []
    for i in range(n):
        _, seq = sequences[i]
        encoded = encode_sequence_gray(seq)
        all_encoded.extend(encoded)

    # Hamming distance distribution
    ham_dist_counts = Counter()
    n_pairs = min(10000, len(all_encoded) * (len(all_encoded) - 1) // 2)

    # Sample pairs for efficiency
    if len(all_encoded) > 100:
        indices = np.random.choice(len(all_encoded), size=(n_pairs, 2), replace=True)
        for i, j in indices:
            if i != j:
                d = gray_hamming_int(all_encoded[i], all_encoded[j])
                ham_dist_counts[d] += 1
    else:
        for i in range(len(all_encoded)):
            for j in range(i + 1, len(all_encoded)):
                d = gray_hamming_int(all_encoded[i], all_encoded[j])
                ham_dist_counts[d] += 1

    total = sum(ham_dist_counts.values())
    ham_dist_ratio = {d: c / total for d, c in ham_dist_counts.items()}

    print(f"  Consensus K-map symmetry error: {symmetry_error:.6f}")
    print(f"  Hamming distance distribution (sampled {total} pairs):")
    for d in sorted(ham_dist_ratio.keys()):
        print(f"    d={d}: {ham_dist_ratio[d]:.4f}")

    # Check if K-map has long-range order
    # (non-zero correlations at large distances in the K-map)
    autocorr = np.correlate(kmap_consensus[:64], kmap_consensus[:64], mode="full")
    autocorr = autocorr[63:]  # Take positive lags
    autocorr /= autocorr[0] if autocorr[0] != 0 else 1

    print(f"  K-map autocorrelation at lag 1: {autocorr[1]:.4f}")
    print(
        f"  K-map autocorrelation at lag 5: {autocorr[5]:.4f}"
        if len(autocorr) > 5
        else ""
    )

    return {
        "symmetry_error": symmetry_error,
        "hamming_distribution": dict(ham_dist_ratio),
        "autocorrelation_lag1": float(autocorr[1]) if len(autocorr) > 1 else 0,
    }


# ============================================================
# 9. Co-evolution Analysis
# ============================================================


def analyze_coevolution(sequences, n_seqs=50):
    """
    Analyze co-evolution patterns in the Spike protein alignment.

    Uses mutual information and correlation analysis between
    residue positions.
    """
    print("\n=== Co-evolution Analysis ===")

    n = min(n_seqs, len(sequences))

    # Extract clean sequences (no gaps)
    clean_seqs = []
    for i in range(n):
        header, seq = sequences[i]
        clean = "".join(aa for aa in seq if aa in _AA_TO_INDEX)
        if len(clean) > 100:
            clean_seqs.append((extract_accession(header), clean))

    if len(clean_seqs) < 5:
        print("  Insufficient clean sequences")
        return {}

    print(f"  Clean sequences: {len(clean_seqs)}")

    # Compute position-specific amino acid frequencies
    min_len = min(len(seq) for _, seq in clean_seqs)
    max_len = min(500, min_len)  # Limit to first 500 positions

    # Count amino acids at each position
    pos_counts = [Counter() for _ in range(max_len)]
    for _, seq in clean_seqs:
        for i in range(min(max_len, len(seq))):
            pos_counts[i][seq[i]] += 1

    # Compute conservation at each position
    conservation = []
    for i in range(max_len):
        total = sum(pos_counts[i].values())
        if total > 0:
            max_freq = max(pos_counts[i].values()) / total
            conservation.append(max_freq)
        else:
            conservation.append(0)

    print(f"  Alignment length analyzed: {max_len}")
    print(f"  Mean conservation: {np.mean(conservation):.4f}")
    print(
        f"  Highly conserved positions (>90%): {sum(1 for c in conservation if c > 0.9)}"
    )
    print(
        f"  Variable positions (<70% conservation): {sum(1 for c in conservation if c < 0.7)}"
    )

    # Compute pairwise position correlations
    n_positions = min(200, max_len)  # Limit computation
    encoding_matrix = np.zeros((len(clean_seqs), n_positions), dtype=int)
    for i, (_, seq) in enumerate(clean_seqs):
        for j in range(n_positions):
            if j < len(seq) and seq[j] in _AA_TO_INDEX:
                encoding_matrix[i, j] = encode_gray_single(seq[j])

    # Compute correlation matrix between positions
    corr_matrix = np.corrcoef(encoding_matrix.T)

    # Find strongly correlated position pairs
    strong_correlations = []
    for i in range(n_positions):
        for j in range(i + 1, n_positions):
            if abs(corr_matrix[i, j]) > 0.7:
                strong_correlations.append((i, j, corr_matrix[i, j]))

    strong_correlations.sort(key=lambda x: abs(x[2]), reverse=True)

    print(f"  Position pairs with |correlation| > 0.7: {len(strong_correlations)}")
    if strong_correlations:
        print(f"  Top 5 correlated pairs:")
        for i, j, corr in strong_correlations[:5]:
            print(f"    Position {i} vs {j}: r = {corr:.4f}")

    return {
        "num_clean_sequences": len(clean_seqs),
        "alignment_length": max_len,
        "conservation": conservation,
        "mean_conservation": np.mean(conservation),
        "num_highly_conserved": sum(1 for c in conservation if c > 0.9),
        "num_variable": sum(1 for c in conservation if c < 0.7),
        "strong_correlations": strong_correlations[:20],
    }


# ============================================================
# 10. Walsh-Hadamard Analysis
# ============================================================


def analyze_walsh_hadamard(sequences, n_seqs=100):
    """
    Walsh-Hadamard transform analysis of K-map.

    Computes the Walsh-Hadamard spectrum of the consensus K-map
    to identify dominant frequency modes.
    """
    print("\n=== Walsh-Hadamard Transform Analysis ===")

    # Build consensus K-map
    kmap = build_kmap_consensus(sequences, n_seqs)

    # Reshape to 32x32
    kmap_2d = kmap.reshape(32, 32)

    # Fast Walsh-Hadamard Transform (simplified)
    # For a 32x32 map, we compute the 1D FWHT on flattened vector
    def fwht(a):
        """Fast Walsh-Hadamard Transform."""
        n = len(a)
        if n <= 1:
            return a
        h = n // 2
        result = np.zeros(n)
        for i in range(h):
            result[i] = a[2 * i] + a[2 * i + 1]
            result[h + i] = a[2 * i] - a[2 * i + 1]
        return fwht(result[:h]) + fwht(result[h:])

    # Compute FWHT on first 32 values (power of 2)
    spectrum = fwht(kmap[:32])
    spectrum = np.abs(spectrum)

    # Normalize
    spectrum /= spectrum.sum() if spectrum.sum() > 0 else 1

    # Find dominant modes
    dominant_modes = np.argsort(spectrum)[::-1][:5]

    print(f"  K-map size: 32x32 = 1024 cells")
    print(f"  Non-zero cells: {np.count_nonzero(kmap)}")
    print(f"  Walsh-Hadamard spectrum (top 5 modes):")
    for mode in dominant_modes:
        print(f"    Mode {mode}: {spectrum[mode]:.4f}")

    # Compute explained variance
    top_modes = spectrum[dominant_modes[:3]]
    explained_var = top_modes.sum() * 100

    print(f"  Top 3 modes explain: {explained_var:.1f}% of variance")

    return {
        "nonzero_cells": int(np.count_nonzero(kmap)),
        "dominant_modes": dominant_modes.tolist(),
        "spectrum_values": spectrum.tolist(),
        "explained_variance_top3": explained_var,
    }


# ============================================================
# 11. Mutation Analysis (H6 analog)
# ============================================================


def analyze_mutations(sequences, n_seqs=100):
    """
    Analyze mutation patterns using Gray code distances.

    Compares sequences to the first (reference) sequence and
    computes mutation profiles.
    """
    print("\n=== Mutation Analysis (H6 analog) ===")

    n = min(n_seqs, len(sequences))

    # Use first sequence as reference
    ref_header, ref_seq = sequences[0]
    ref_encoded = encode_sequence_gray(ref_seq)
    ref_clean = "".join(aa for aa in ref_seq if aa in _AA_TO_INDEX)

    print(f"  Reference: {extract_accession(ref_header)}")
    print(f"  Reference length: {len(ref_clean)}")

    # Compute mutations for each sequence
    mutation_counts = []
    mutation_positions = Counter()

    for i in range(1, n):
        header, seq = sequences[i]
        encoded = encode_sequence_gray(seq)

        min_len = min(len(ref_encoded), len(encoded))
        mutations = 0
        for j in range(min_len):
            if ref_encoded[j] != encoded[j]:
                mutations += 1
                mutation_positions[j] += 1

        mutation_counts.append(mutations)

    print(f"  Sequences compared: {n - 1}")
    print(f"  Mean mutations per sequence: {np.mean(mutation_counts):.1f}")
    print(f"  Max mutations: {max(mutation_counts)}")
    print(f"  Min mutations: {min(mutation_counts)}")

    # Find most variable positions
    top_variable = mutation_positions.most_common(10)
    print(f"  Top 10 most variable positions:")
    for pos, count in top_variable:
        pct = count / (n - 1) * 100
        print(f"    Position {pos}: {count} mutations ({pct:.1f}%)")

    # Compute Gray code distance distribution of mutations
    mutation_distances = []
    dist_counter = Counter()
    for i in range(1, n):
        header, seq = sequences[i]
        encoded = encode_sequence_gray(seq)
        min_len = min(len(ref_encoded), len(encoded))
        for j in range(min_len):
            if ref_encoded[j] != encoded[j]:
                d = gray_hamming_int(ref_encoded[j], encoded[j])
                mutation_distances.append(d)

    if mutation_distances:
        dist_counter = Counter(mutation_distances)
        total_muts = sum(dist_counter.values())
        print(f"  Mutation Hamming distance distribution:")
        for d in sorted(dist_counter.keys()):
            pct = dist_counter[d] / total_muts * 100
            print(f"    d={d}: {dist_counter[d]} ({pct:.1f}%)")

    return {
        "num_sequences": n - 1,
        "mean_mutations": float(np.mean(mutation_counts)),
        "max_mutations": max(mutation_counts),
        "min_mutations": min(mutation_counts),
        "mutation_counts": mutation_counts,
        "top_variable_positions": top_variable,
        "mutation_distance_distribution": dict(dist_counter),
    }


# ============================================================
# 12. Main Pipeline
# ============================================================


def main():
    """Run the complete K-map analysis pipeline."""

    # Paths
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    results_dir = base_dir / "kmap_results"
    results_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("K-map Analysis of SARS-CoV-2 Spike Protein Co-evolution")
    print("=" * 70)

    # 1. Parse FASTA
    print("\n[1/10] Parsing FASTA alignment...")
    sequences = parse_fasta(fasta_file)
    print(f"  Loaded {len(sequences)} sequences")

    # 2. Sequence statistics
    print("\n[2/10] Computing sequence statistics...")
    stats = compute_sequence_stats(sequences)
    print(f"  Alignment length: {stats['alignment_length']}")
    print(f"  Gap count: {stats['gap_count']}")
    print(f"  Unknown residues: {stats['unknown_count']}")

    # 3. Gray code encoding
    print("\n[3/10] Encoding sequences with 5-bit Gray code...")
    encoded_seqs = []
    for header, seq in sequences:
        encoded = encode_sequence_gray(seq)
        encoded_seqs.append(encoded)
    print(f"  Encoded {len(encoded_seqs)} sequences")

    # 4. H1: Adjacency analysis
    print("\n[4/10] Running H1: Gray-code adjacency analysis...")
    # BUG FIX: was n_seqs=200 (only 200 of 1299 sequences). Now ALL sequences.
    h1_results = analyze_h1_adjacency(sequences, n_seqs=1299)

    # 5. H2: Signature analysis
    print("\n[5/10] Running H2: K-map signature analysis...")
    h2_results = analyze_h2_signature(sequences, n_seqs=1299)

    # 6. H3: Structural classification
    print("\n[6/10] Running H3: Structural classification...")
    h3_results = analyze_h3_structure(sequences, n_seqs=1299)

    # 7. H5: Contact map invariants
    print("\n[7/10] Running H5: Contact map invariants...")
    h5_results = analyze_h5_contacts(sequences, n_seqs=1299)

    # 8. Co-evolution analysis
    print("\n[8/10] Running co-evolution analysis...")
    coevo_results = analyze_coevolution(sequences, n_seqs=1299)

    # 9. Walsh-Hadamard analysis
    print("\n[9/10] Running Walsh-Hadamard transform analysis...")
    walsh_results = analyze_walsh_hadamard(sequences, n_seqs=1299)

    # 10. Mutation analysis
    print("\n[10/10] Running mutation analysis...")
    mutation_results = analyze_mutations(sequences, n_seqs=1299)

    # ============================================================
    # Save Results
    # ============================================================

    print("\n" + "=" * 70)
    print("Saving Results...")
    print("=" * 70)

    # Save summary JSON
    summary = {
        "dataset": "SARS-CoV-2 Spike Protein Alignment",
        "num_sequences": len(sequences),
        "alignment_length": stats["alignment_length"],
        "analysis_results": {
            "h1_adjacency": h1_results,
            "h2_signature": h2_results,
            "h3_structure": h3_results,
            "h5_contacts": h5_results,
            "coevolution": coevo_results,
            "walsh_hadamard": walsh_results,
            "mutations": mutation_results,
        },
    }

    # Remove non-serializable items
    for key in ["conservation", "mutation_counts"]:
        if key in summary["analysis_results"].get("coevolution", {}):
            summary["analysis_results"]["coevolution"][key] = summary[
                "analysis_results"
            ]["coevolution"][key][:50]
        if key in summary["analysis_results"].get("mutations", {}):
            summary["analysis_results"]["mutations"][key] = summary["analysis_results"][
                "mutations"
            ][key][:50]

    with open(results_dir / "analysis_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # Save consensus K-map
    consensus_kmap = build_kmap_consensus(sequences, n_seqs=1299)
    np.save(results_dir / "consensus_kmap.npy", consensus_kmap)
    np.savetxt(
        results_dir / "consensus_kmap.csv",
        consensus_kmap.reshape(32, 32),
        delimiter=",",
    )

    # Save individual K-maps for first 10 sequences
    for i in range(min(10, len(sequences))):
        header, seq = sequences[i]
        acc = extract_accession(header)
        clean_seq = "".join(aa for aa in seq if aa in _AA_TO_INDEX)
        if len(clean_seq) > 1:
            kmap = build_aa_kmap_2d(clean_seq)
            np.save(results_dir / f"kmap_{acc}.npy", kmap)

    # Save encoded sequences (first 50)
    with open(results_dir / "encoded_sequences.json", "w") as f:
        encoded_data = []
        for i in range(min(50, len(sequences))):
            header, seq = sequences[i]
            encoded = encode_sequence_gray(seq)
            encoded_data.append(
                {
                    "accession": extract_accession(header),
                    "length": len(seq),
                    "encoded_length": len(encoded),
                    "gray_values": encoded[:100],  # First 100 values
                }
            )
        json.dump(encoded_data, f, indent=2)

    # Save H1 detailed results
    with open(results_dir / "h1_adjacency_results.json", "w") as f:
        json.dump(h1_results, f, indent=2)

    print(f"\nResults saved to: {results_dir}")
    print(f"  - analysis_summary.json")
    print(f"  - consensus_kmap.npy / .csv")
    print(f"  - kmap_*.npy (individual K-maps)")
    print(f"  - encoded_sequences.json")
    print(f"  - h1_adjacency_results.json")

    # Print final summary
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nKey Results:")
    print(
        f"  H1 Adjacency: Ratio = {h1_results['observed_ratio']:.4f} (enrichment = {h1_results['enrichment_ratio']:.2f}x)"
    )
    print(
        f"  H3 Structure: Avg K-map distance = {h3_results.get('avg_distance', 'N/A')}"
    )
    print(f"  H5 Contacts: Symmetry error = {h5_results.get('symmetry_error', 'N/A')}")
    print(
        f"  Co-evolution: {coevo_results.get('num_clean_sequences', 0)} clean sequences analyzed"
    )
    print(
        f"  Walsh-Hadamard: Top 3 modes explain {walsh_results.get('explained_variance_top3', 0):.1f}% variance"
    )
    print(
        f"  Mutations: Mean = {mutation_results.get('mean_mutations', 0):.1f} per sequence"
    )


if __name__ == "__main__":
    main()
