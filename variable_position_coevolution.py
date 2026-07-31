#!/usr/bin/env python3
"""
Variable-Position Co-evolution via K-map with Don't-Care Conditions
===================================================================

Key insight: Co-evolution happens at VARIABLE positions, not conserved ones.
We use the don't-care condition strategically:
- Conserved positions → don't-care in Boolean function
- Variable positions → on-set/off-set based on co-evolutionary coupling

This extracts MUTATIONS as prime implicants, not conservation patterns.
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


def compute_position_entropy(pos_arrays, pos, n_seqs):
    """Compute entropy at a position (lower = more conserved)."""
    counts = Counter()
    for arr in pos_arrays[:n_seqs]:
        if pos < len(arr) and arr[pos] >= 0:
            counts[int(arr[pos])] += 1

    total = sum(counts.values())
    if total == 0:
        return 0.0

    entropy = 0.0
    for c in counts.values():
        p = c / total
        if p > 0:
            entropy -= p * np.log2(p)
    return entropy


def get_majority_ref(pos_arrays, pos, n_seqs):
    """Return the most common residue code at a position (int 0-19)."""
    return Counter(
        int(a[pos]) for a in pos_arrays[:n_seqs] if pos < len(a) and a[pos] >= 0
    ).most_common(1)[0][0]


def find_variable_positions(pos_arrays, n_seqs, entropy_threshold=0.5, max_pos=None):
    if max_pos is None:
        max_pos = (
            pos_arrays.shape[1] if hasattr(pos_arrays, "shape") else len(pos_arrays[0])
        )
    """Find positions with high entropy (variable positions)."""
    print("  Computing position entropies...")
    entropies = []
    for pos in range(max_pos):
        ent = compute_position_entropy(pos_arrays, pos, n_seqs)
        entropies.append(ent)

    entropies = np.array(entropies)

    # Variable positions: entropy > threshold
    variable_mask = entropies > entropy_threshold
    variable_positions = np.where(variable_mask)[0]

    print(f"  Entropy threshold: {entropy_threshold}")
    print(f"  Variable positions: {len(variable_positions)} / {max_pos}")
    print(f"  Conserved positions: {max_pos - len(variable_positions)} / {max_pos}")

    return variable_positions, entropies


def build_mutation_kmap(pos_arrays, pos_i, pos_j, variable_positions, n_seqs):
    """
    Build a K-map that captures MUTATIONS at variable positions.

    For each sequence:
    - Extract residues at positions pos_i and pos_j
    - If BOTH positions are variable → this is a co-evolutionary observation
    - If either is conserved → don't-care

    The K-map cell (aa_i, aa_j) = 1 if this mutation pair appears.
    """
    kmap = np.zeros((20, 20), dtype=np.float64)
    n_observations = 0

    for arr in pos_arrays[:n_seqs]:
        if pos_i < len(arr) and pos_j < len(arr):
            ci, cj = int(arr[pos_i]), int(arr[pos_j])
            if ci >= 0 and cj >= 0:
                kmap[ci, cj] += 1
                n_observations += 1

    if n_observations > 0:
        kmap /= n_observations

    return kmap, n_observations


def build_mutation_kmap_with_dontcare(
    pos_arrays,
    pos_i,
    pos_j,
    variable_positions_i,
    variable_positions_j,
    reference_seq,
    n_seqs,
):
    """
    Build K-map with don't-care conditions.

    - If position i is VARIABLE: use the actual residue
    - If position i is CONSERVED: mark as don't-care (-1)
    - Same for position j

    The Boolean function becomes:
    f(aa_i, aa_j) = 1 if (aa_i, aa_j) appears at variable positions
                    = -1 (don't-care) if either position is conserved
                    = 0 otherwise
    """
    kmap = np.zeros((20, 20), dtype=np.int32)  # -1=don't-care, 0=off, 1=on

    # First, determine which observations are "co-evolutionary"
    # (both positions variable) vs "conservation" (either position conserved)

    for arr in pos_arrays[:n_seqs]:
        if pos_i < len(arr) and pos_j < len(arr):
            ci, cj = int(arr[pos_i]), int(arr[pos_j])
            if ci >= 0 and cj >= 0:
                # Check if this is a variable observation
                is_variable = ci != reference_seq[pos_i] or cj != reference_seq[pos_j]

                if is_variable:
                    # This is a mutation — mark as on-set
                    kmap[ci, cj] = 1
                else:
                    # This is the reference — don't-care
                    kmap[ci, cj] = -1

    return kmap


def find_coevolutionary_pairs(pos_arrays, variable_positions, n_seqs, top_n=20):
    """
    Find position pairs where mutations at one position correlate with
    mutations at another position.

    Uses conditional entropy: H(j | i) = H(i,j) - H(i)
    Low conditional entropy = strong co-evolution
    """
    print("\n  Computing conditional entropies...")

    co_evolving = []
    var_set = set(variable_positions)

    # Only consider pairs where both positions are variable
    for idx_i, pos_i in enumerate(variable_positions):
        for idx_j in range(idx_i + 1, len(variable_positions)):
            pos_j = variable_positions[idx_j]

            # Skip if too far apart
            if abs(pos_i - pos_j) > 30:
                continue

            # Joint distribution
            joint = Counter()
            marg_i = Counter()
            marg_j = Counter()

            for arr in pos_arrays[:n_seqs]:
                if pos_i < len(arr) and pos_j < len(arr):
                    ci, cj = int(arr[pos_i]), int(arr[pos_j])
                    if ci >= 0 and cj >= 0:
                        # BUG FIX: Original code used pos_arrays[0] (first sequence)
                        # as reference instead of the majority reference. This is
                        # inconsistent with all other scripts. Fixed to use majority.
                        ref_i = get_majority_ref(pos_arrays, pos_i, n_seqs)
                        ref_j = get_majority_ref(pos_arrays, pos_j, n_seqs)

                        if ci != ref_i or cj != ref_j:
                            joint[(ci, cj)] += 1
                            marg_i[ci] += 1
                            marg_j[cj] += 1

            total = sum(joint.values())
            if total < 5:  # Need enough mutations
                continue

            # Mutual information
            mi = 0.0
            for (ai, aj), count in joint.items():
                p_joint = count / total
                p_i = marg_i[ai] / total
                p_j = marg_j[aj] / total
                if p_joint > 0 and p_i > 0 and p_j > 0:
                    mi += p_joint * np.log2(p_joint / (p_i * p_j))

            if mi > 0.1:  # Only significant co-evolution
                co_evolving.append((pos_i, pos_j, mi, total))

    co_evolving.sort(key=lambda x: x[2], reverse=True)
    return co_evolving[:top_n]


def minimize_with_dontcare(kmap_with_dc):
    """
    Run QM on a K-map with don't-care conditions.

    -1 = don't-care (can be 0 or 1)
     0 = off-set (must be 0)
     1 = on-set (must be 1)
    """
    # Convert to format expected by QM
    # QM expects: 0 = off, 1 = on, -1 = don't-care
    bool_map = kmap_with_dc.copy()

    # Flatten and run QM
    bool_flat = bool_map.flatten().astype(int)
    result = boolean_minimize_kmap(bool_flat, algorithm="qm")

    return result


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    results_dir = base_dir / "variable_position_results"
    results_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("Variable-Position Co-evolution via K-map with Don't-Care")
    print("=" * 70)

    # Load ALL sequences
    print("\n[1/5] Loading ALL sequences...")
    sequences = parse_fasta(fasta_file)
    full_length = len(sequences[0][1])
    encoder = Base20AminoEncoder(version=1)
    aa_list = list(AMINO_HE_2012)
    n_all = len(sequences)
    print(f"  Total: {n_all} sequences")

    # Build position arrays
    print("\n[2/5] Building position arrays...")
    pos_arrays = build_position_arrays(sequences, encoder, max_pos=full_length)

    # Find variable positions
    print("\n[3/5] Finding variable positions...")
    variable_positions, entropies = find_variable_positions(
        pos_arrays, n_all, entropy_threshold=0.3, max_pos=full_length
    )

    print(f"\n  Variable positions (entropy > 0.3):")
    for pos in variable_positions[:20]:
        print(f"    Position {pos:3d}: entropy = {entropies[pos]:.3f}")

    # Find co-evolutionary pairs
    print("\n[4/5] Finding co-evolutionary position pairs...")
    co_evolving = find_coevolutionary_pairs(
        pos_arrays, variable_positions, n_all, top_n=20
    )

    print(f"\n  Top co-evolutionary pairs (mutations only):")
    print(f"  {'Pos i':>5s} {'Pos j':>5s} {'MI':>8s} {'N muts':>8s}")
    print(f"  {'-' * 5} {'-' * 5} {'-' * 8} {'-' * 8}")
    for pos_i, pos_j, mi, n in co_evolving:
        print(f"  {pos_i:5d} {pos_j:5d} {mi:8.4f} {n:8d}")

    # Build mutation K-maps with don't-care conditions
    print("\n[5/5] Building mutation K-maps with don't-care conditions...")

    # BUG FIX: Original code used pos_arrays[0] (first sequence) as reference.
    # This is inconsistent with find_coevolutionary_pairs() which uses majority_ref.
    # Fixed: build a synthetic majority reference array.
    ref_codes = []
    for pos in range(len(pos_arrays[0])):
        ref_codes.append(get_majority_ref(pos_arrays, pos, n_all))
    reference = np.array(ref_codes, dtype=np.int32)

    all_results = []
    for idx, (pos_i, pos_j, mi, n_muts) in enumerate(co_evolving[:10]):
        print(f"\n  --- Pair ({pos_i}, {pos_j}), MI={mi:.4f}, {n_muts} mutations ---")

        # Build K-map with don't-care conditions
        kmap_dc = build_mutation_kmap_with_dontcare(
            pos_arrays,
            pos_i,
            pos_j,
            variable_positions,
            variable_positions,
            reference,
            n_all,
        )

        # Count on-set, off-set, don't-care
        n_on = int((kmap_dc == 1).sum())
        n_off = int((kmap_dc == 0).sum())
        n_dc = int((kmap_dc == -1).sum())

        print(f"  On-set (mutations): {n_on} cells")
        print(f"  Off-set (never seen): {n_off} cells")
        print(f"  Don't-care (conserved): {n_dc} cells")

        # Run QM with don't-care
        result = minimize_with_dontcare(kmap_dc)

        print(f"  Prime implicants: {result['n_prime_implicants']}")
        print(f"  Essential PIs: {result['n_essential']}")

        # Decode essential PIs
        print(f"  Essential prime implicants (mutation motifs):")
        for pi_idx, pi in enumerate(result["essential_prime_implicants"][:10]):
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

            # Check if this is a mutation (different from reference)
            ref_i_val = int(reference[pos_i]) if pos_i < len(reference) else -1
            ref_j_val = int(reference[pos_j]) if pos_j < len(reference) else -1

            is_mutation = row_code % 20 != ref_i_val or col_code % 20 != ref_j_val
            mut_str = "MUTATION" if is_mutation else "reference"

            print(f"    PI_{pi_idx + 1}: {term_str} = ({row_aa},{col_aa}) [{mut_str}]")

        all_results.append(
            {
                "pos_i": pos_i,
                "pos_j": pos_j,
                "mi": mi,
                "n_mutations": n_muts,
                "n_on_set": n_on,
                "n_off_set": n_off,
                "n_dontcare": n_dc,
                "n_prime_implicants": result["n_prime_implicants"],
                "n_essential": result["n_essential"],
            }
        )

    # Save results
    print("\n" + "=" * 70)
    print("Saving results...")
    summary = {
        "dataset": "SARS-CoV-2 Spike Protein (ALL sequences)",
        "num_sequences": n_all,
        "method": "Variable-position K-map with don't-care conditions",
        "variable_positions": variable_positions.tolist(),
        "co_evolving_pairs": co_evolving,
        "results": all_results,
    }

    with open(results_dir / "variable_position_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nResults saved to: {results_dir}")

    # Summary
    print("\n" + "=" * 70)
    print("VARIABLE-POSITION CO-EVOLUTION ANALYSIS")
    print("=" * 70)
    print(f"\nSequences: {n_all}")
    print(f"Variable positions: {len(variable_positions)}")
    print(f"Co-evolutionary pairs: {len(co_evolving)}")
    print(f"\nThe don't-care condition allows us to:")
    print(f"  1. Focus on MUTATIONS (variable positions)")
    print(f"  2. Ignore CONSERVATION (conserved positions)")
    print(f"  3. Extract prime implicants that capture co-evolutionary MOTIFS")
    print(f"     (not conservation patterns)")


if __name__ == "__main__":
    main()
