#!/usr/bin/env python3
"""
Flipped Boolean Co-evolution Analysis
======================================

Instead of: "what co-evolves?" (on-set = observed mutations)
We ask:     "what CANNOT co-evolve?" (on-set = forbidden pairs)

The forbidden pairs are more restrictive and give us
negative constraints: "IF pos i = X THEN pos j CANNOT be Y"
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


def compute_entropy(pos_arrays, pos, n_seqs):
    counts = Counter()
    for arr in pos_arrays[:n_seqs]:
        if pos < len(arr) and arr[pos] >= 0:
            counts[int(arr[pos])] += 1
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((c / total) * np.log2(c / total) for c in counts.values() if c > 0)


def get_majority_ref(pos_arrays, pos, n_seqs):
    return Counter(
        int(a[pos]) for a in pos_arrays[:n_seqs] if pos < len(a) and a[pos] >= 0
    ).most_common(1)[0][0]


def compute_mi(pos_arrays, pos_i, pos_j, n_seqs):
    ref_i = get_majority_ref(pos_arrays, pos_i, n_seqs)
    ref_j = get_majority_ref(pos_arrays, pos_j, n_seqs)
    joint, marg_i, marg_j = Counter(), Counter(), Counter()
    for arr in pos_arrays[:n_seqs]:
        if pos_i < len(arr) and pos_j < len(arr):
            ci, cj = int(arr[pos_i]), int(arr[pos_j])
            if ci >= 0 and cj >= 0 and (ci != ref_i or cj != ref_j):
                joint[(ci, cj)] += 1
                marg_i[ci] += 1
                marg_j[cj] += 1
    total = sum(joint.values())
    if total < 5:
        return 0.0, 0
    mi = sum(
        (c / total)
        * np.log2((c / total) / ((marg_i[ai] / total) * (marg_j[aj] / total)))
        for (ai, aj), c in joint.items()
        if marg_i[ai] > 0 and marg_j[aj] > 0
    )
    return mi, total


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    results_dir = base_dir / "flipped_boolean_results"
    results_dir.mkdir(exist_ok=True)

    print("=" * 80)
    print("FLIPPED Boolean Co-evolution Analysis")
    print("=" * 80)
    print("On-set (1) = FORBIDDEN pairs (never observed together)")
    print("Off-set (0) = ALLOWED pairs (observed together)")
    print("DC = conserved positions (ignore)")
    print("")

    # Load ALL sequences
    print("[1/4] Loading ALL sequences...")
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

    # Find co-evolutionary pairs
    print("\n[3/4] Finding co-evolutionary pairs...")
    co_evolving = []
    for idx_i, pos_i in enumerate(variable_positions):
        for idx_j in range(idx_i + 1, len(variable_positions)):
            pos_j = variable_positions[idx_j]
            if abs(pos_i - pos_j) > 30:
                continue
            mi, n_muts = compute_mi(pos_arrays, pos_i, pos_j, n_all)
            if mi > 0.1:
                ref_i = get_majority_ref(pos_arrays, pos_i, n_all)
                ref_j = get_majority_ref(pos_arrays, pos_j, n_all)
                co_evolving.append((pos_i, pos_j, mi, n_muts, ref_i, ref_j))

    co_evolving.sort(key=lambda x: x[2], reverse=True)
    print(f"  Co-evolutionary pairs: {len(co_evolving)}")

    # Build FLIPPED K-maps and Boolean functions
    print("\n[4/4] Building flipped K-maps...")
    all_rules = []
    rule_num = 0

    for pos_i, pos_j, mi, n_muts, ref_i, ref_j in co_evolving[:15]:
        # BUG FIX: ref_i, ref_j are already int codes (0-19) from get_majority_ref.
        # Previous code did aa_list.index(ref_i) on a list of strings → always 0.
        ref_i_code = int(ref_i)
        ref_j_code = int(ref_j)
        ref_aa_i = aa_list[ref_i_code] if 0 <= ref_i_code < 20 else "?"
        ref_aa_j = aa_list[ref_j_code] if 0 <= ref_j_code < 20 else "?"

        # Build K-map with FLIPPED logic:
        # 1 = FORBIDDEN (never observed together at variable positions)
        # 0 = ALLOWED (observed together)
        # DC = conserved (reference pair)
        #
        # BUG FIX: The original code had two bugs:
        #   (a) Line 156 reset kmap to all 0 immediately after line 155 set it to all 1,
        #       defeating the purpose. No cell was ever set to 1 (forbidden).
        #   (b) The second loop (lines 172-179) was a duplicate of the first (159-168),
        #       doing nothing new.
        # The fix: start with all 1 (forbidden), mark observed pairs as 0 (allowed),
        # mark reference as -1 (DC). Cells that remain 1 are truly forbidden.
        kmap = np.full((20, 20), 1, dtype=np.int32)  # All FORBIDDEN by default

        # Mark observed pairs as ALLOWED (0), reference pair as DC (-1)
        for arr in pos_arrays[:n_all]:
            if pos_i < len(arr) and pos_j < len(arr):
                ci, cj = int(arr[pos_i]), int(arr[pos_j])
                if ci >= 0 and cj >= 0:
                    if ci == ref_i_code and cj == ref_j_code:
                        kmap[ci, cj] = -1  # DC (conserved reference)
                    else:
                        kmap[ci, cj] = 0  # ALLOWED (observed mutation)

        # Count cells
        n_forbidden = int((kmap == 1).sum())
        n_allowed = int((kmap == 0).sum())
        n_dc = int((kmap == -1).sum())

        print(f"\n  Position pair ({pos_i}, {pos_j}):")
        print(f"    Forbidden (never observed): {n_forbidden} cells")
        print(f"    Allowed (observed): {n_allowed} cells")
        print(f"    Don't-care (conserved): {n_dc} cells")

        # Run QM on FLIPPED map
        bool_flat = kmap.flatten().astype(int)
        result = boolean_minimize_kmap(bool_flat, algorithm="qm")

        print(f"    Prime implicants: {result['n_prime_implicants']}")
        print(f"    Essential PIs: {result['n_essential']}")

        # Extract rules
        print(f"    Forbidden pairs (co-evolution constraints):")
        for pi in result["prime_implicants"]:
            values = list(pi["values"])
            mask = list(pi["mask"])
            while len(values) < 8:
                values.append(0)
                mask.append(False)

            row_code = sum(values[j] * (2 ** (3 - j)) for j in range(4) if not mask[j])
            col_code = sum(
                values[j + 4] * (2 ** (3 - j)) for j in range(4) if not mask[j + 4]
            )

            if row_code < 20 and col_code < 20:
                aa_i = aa_list[row_code]
                aa_j = aa_list[col_code]
                is_essential = pi in result["essential_prime_implicants"]
                marker = "*" if is_essential else " "
                rule_num += 1
                print(
                    f"    {marker} Rule {rule_num}: IF pos {pos_i}={aa_i} AND pos {pos_j}={aa_j} THEN FORBIDDEN"
                )

                all_rules.append(
                    {
                        "pos_i": pos_i,
                        "pos_j": pos_j,
                        "aa_i": aa_i,
                        "aa_j": aa_j,
                        "ref_i": ref_aa_i,
                        "ref_j": ref_aa_j,
                        "mi": mi,
                        "is_essential": is_essential,
                        "type": "forbidden",
                    }
                )

    # Save results
    print("\n" + "=" * 80)
    print(f"Total forbidden rules: {rule_num}")
    print(f"Results saved to: {results_dir}")

    summary = {
        "dataset": f"{n_all} sequences",
        "variable_positions": len(variable_positions),
        "co_evolving_pairs": len(co_evolving),
        "total_forbidden_rules": rule_num,
        "rules": all_rules,
    }
    with open(results_dir / "flipped_boolean_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)


if __name__ == "__main__":
    main()
