#!/usr/bin/env python3
"""
Master Boolean Function for Co-evolution Prediction
====================================================

Builds the complete set of prime implicants and essential PIs
over ALL sequences, then constructs the master Boolean function.

The master function is: f(pos_i, pos_j, aa_i, aa_j) = 1
if ANY essential prime implicant matches.

This gives us a complete set of co-evolutionary inference rules.
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
    """Parse FASTA alignment file."""
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
    """Compute Shannon entropy at a position (gaps excluded)."""
    counts = Counter()
    for arr in pos_arrays[:n_seqs]:
        if pos < len(arr):
            c = int(arr[pos])
            if 0 <= c < 20:  # exclude gaps (20) and legacy -1
                counts[c] += 1
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((c / total) * np.log2(c / total) for c in counts.values() if c > 0)


def find_variable_positions(pos_arrays, n_seqs, max_pos=80, threshold=0.3):
    """Find positions with entropy above threshold."""
    return [
        p for p in range(max_pos) if compute_entropy(pos_arrays, p, n_seqs) > threshold
    ]


def find_coevolutionary_pairs(
    pos_arrays, variable_positions, n_seqs, max_gap=30, min_mi=0.1
):
    """Find position pairs with significant co-evolution (MI > min_mi).

    GPU-accelerated via torch CUDA (coevolution_gpu) — mutation-only MI.
    Falls back to CPU Counter implementation if CUDA unavailable.
    """
    co_evolving = []

    # Build pair list (within max_gap)
    pairs = [
        (pos_i, pos_j)
        for idx_i, pos_i in enumerate(variable_positions)
        for pos_j in variable_positions[idx_i + 1 :]
        if abs(pos_i - pos_j) <= max_gap
    ]

    # GPU path
    try:
        import coevolution_gpu as cg

        dense = cg.dense_to_gpu(pos_arrays)
        refs = cg.majority_refs_gpu(dense)
        mi_dict, cnt_dict = cg.mi_matrix_gpu(
            dense,
            pairs,
            refs=refs,
            mutation_only=True,
            min_total=5,
            chunk=16384,
        )
        for (pos_i, pos_j), mi in mi_dict.items():
            n_mut = cnt_dict[(pos_i, pos_j)]
            if mi > min_mi:
                co_evolving.append(
                    (pos_i, pos_j, mi, n_mut, int(refs[pos_i]), int(refs[pos_j]))
                )
        co_evolving.sort(key=lambda x: x[2], reverse=True)
        return co_evolving
    except Exception as e:
        print(f"  GPU MI failed ({e}), using CPU...")

    for idx_i, pos_i in enumerate(variable_positions):
        for idx_j in range(idx_i + 1, len(variable_positions)):
            pos_j = variable_positions[idx_j]
            if abs(pos_i - pos_j) > max_gap:
                continue

            # Majority reference (exclude gap state 20)
            ref_i = Counter(
                int(a[pos_i])
                for a in pos_arrays[:n_seqs]
                if pos_i < len(a) and 0 <= int(a[pos_i]) < 20
            ).most_common(1)[0][0]
            ref_j = Counter(
                int(a[pos_j])
                for a in pos_arrays[:n_seqs]
                if pos_j < len(a) and 0 <= int(a[pos_j]) < 20
            ).most_common(1)[0][0]

            # Joint distribution of mutations
            joint = Counter()
            marg_i = Counter()
            marg_j = Counter()

            for arr in pos_arrays[:n_seqs]:
                if pos_i < len(arr) and pos_j < len(arr):
                    ci, cj = int(arr[pos_i]), int(arr[pos_j])
                    if 0 <= ci < 20 and 0 <= cj < 20 and (ci != ref_i or cj != ref_j):
                        joint[(ci, cj)] += 1
                        marg_i[ci] += 1
                        marg_j[cj] += 1

            total = sum(joint.values())
            if total < 5:
                continue

            mi = sum(
                (c / total)
                * np.log2((c / total) / ((marg_i[ai] / total) * (marg_j[aj] / total)))
                for (ai, aj), c in joint.items()
                if marg_i[ai] > 0 and marg_j[aj] > 0
            )

            if mi > min_mi:
                co_evolving.append((pos_i, pos_j, mi, total, ref_i, ref_j))

    co_evolving.sort(key=lambda x: x[2], reverse=True)
    return co_evolving


def build_mutation_kmap(pos_arrays, pos_i, pos_j, ref_i, ref_j, n_seqs):
    """Build 20x20 mutation K-map, padded to 32x32 for correct QM encoding.

    CORRECTED (FIX A2): The original passed a 400-cell map to a QM that
    derives k_bits = int(log2(400))//2 = 4 -> 8 bits, wrapping cells
    256-399 onto 0-143 and producing phantom rules (verified: 143/152
    rules had labels never observed). Fix: build the 20x20 mutation map,
    pad to 32x32 (5 bits per axis = 10 bits total), rows/cols 20-31 as
    don't-care (-1). All 20 residues are now representable.

    CORRECTED (FIX A1): callers must pass ALIGNED position arrays
    (gap = 20 state) so column j is the same raw position in every
    sequence. Cells with gaps are excluded (gap-gap pairs are not
    counted as mutations; positions with gap in a sequence are skipped).
    """
    kmap = np.full((32, 32), -1, dtype=np.int32)  # padding rows/cols 20-31
    kmap[:20, :20] = 0  # 20x20 region: default OFF-SET (never observed)
    for arr in pos_arrays[:n_seqs]:
        if pos_i < len(arr) and pos_j < len(arr):
            ci, cj = int(arr[pos_i]), int(arr[pos_j])
            # skip gaps (state 20) and legacy -1
            if ci < 0 or cj < 0 or ci >= 20 or cj >= 20:
                continue
            if ci == ref_i and cj == ref_j:
                kmap[ci, cj] = -1  # reference: don't-care
            else:
                kmap[ci, cj] = 1  # observed mutation
    return kmap


def extract_prime_implicants(result, aa_list):
    """Extract all prime implicants from QM result.

    CORRECTED (FIX A2): the QM now runs on a 32x32 (10-bit) map, so each
    implicant has 10 bits: 5 for row, 5 for column. Values 20-31 on either
    axis are padding (don't-care) and are rejected.
    """
    pis = []
    for pi in result["prime_implicants"]:
        values = list(pi["values"])
        mask = list(pi["mask"])
        while len(values) < 10:
            values.append(0)
            mask.append(False)

        row_code = sum(values[j] * (2 ** (4 - j)) for j in range(5) if not mask[j])
        col_code = sum(
            values[j + 5] * (2 ** (4 - j)) for j in range(5) if not mask[j + 5]
        )

        if row_code < 20 and col_code < 20:
            pis.append(
                {
                    "aa_i": aa_list[row_code],
                    "aa_j": aa_list[col_code],
                    "row_code": row_code,
                    "col_code": col_code,
                    "n_dontcares": pi["n_dontcares"],
                    "coverage": len(pi["coverage"]),
                    "is_essential": pi in result["essential_prime_implicants"],
                }
            )
    return pis


def main():
    base_dir = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
    fasta_file = base_dir / "Spike_protein.aln-fasta"
    results_dir = base_dir / "master_boolean"
    results_dir.mkdir(exist_ok=True)

    print("=" * 80)
    print("MASTER BOOLEAN FUNCTION FOR CO-EVOLUTION")
    print("=" * 80)

    # 1. Load ALL sequences
    print("\n[1/5] Loading ALL sequences...")
    sequences = parse_fasta(fasta_file)
    encoder = Base20AminoEncoder(version=1)
    aa_list = list(AMINO_HE_2012)
    n_all = len(sequences)
    full_length = len(sequences[0][1])
    print(f"  Total: {n_all} sequences")

    # 2. Build position arrays — CORRECTED (FIX A1): aligned columns, gap = 20
    print("\n[2/5] Building position arrays (aligned, gap=20)...")
    max_pos = len(sequences[0][1])
    pos_arrays = []
    for _, seq in sequences:
        arr = np.array(
            [encoder.encode.get(aa, 20) for aa in seq[:max_pos]], dtype=np.int32
        )
        pos_arrays.append(arr)

    # 3. Find variable positions and co-evolutionary pairs
    print("\n[3/5] Finding variable positions and co-evolutionary pairs...")
    variable_positions = find_variable_positions(
        pos_arrays, n_all, max_pos, threshold=0.3
    )
    print(f"  Variable positions: {len(variable_positions)} / {max_pos}")

    co_evolving = find_coevolutionary_pairs(pos_arrays, variable_positions, n_all)
    print(f"  Co-evolutionary pairs (MI > 0.1): {len(co_evolving)}")

    print(f"\n  Top 20 co-evolutionary pairs:")
    print(
        f"  {'Pos i':>5s} {'Pos j':>5s} {'MI':>8s} {'N muts':>8s} {'Ref i':>5s} {'Ref j':>5s}"
    )
    print(f"  {'-' * 5} {'-' * 5} {'-' * 8} {'-' * 8} {'-' * 5} {'-' * 5}")
    for pos_i, pos_j, mi, n, ref_i, ref_j in co_evolving[:20]:
        ri = aa_list[ref_i] if 0 <= ref_i < 20 else "?"
        rj = aa_list[ref_j] if 0 <= ref_j < 20 else "?"
        print(f"  {pos_i:5d} {pos_j:5d} {mi:8.4f} {n:8d} {ri:>5s} {rj:>5s}")

    # 4. Build prime implicants for all pairs
    print("\n[4/5] Building prime implicants...")
    all_pis = []
    all_essential = []
    all_inferences = []

    for pos_i, pos_j, mi, n_muts, ref_i, ref_j in co_evolving[:15]:
        kmap = build_mutation_kmap(pos_arrays, pos_i, pos_j, ref_i, ref_j, n_all)
        bool_flat = kmap.flatten().astype(int)
        result = boolean_minimize_kmap(bool_flat, algorithm="qm")

        pis = extract_prime_implicants(result, aa_list)

        for pi in pis:
            pi["pos_i"] = pos_i
            pi["pos_j"] = pos_j
            pi["ref_i"] = aa_list[ref_i]
            pi["ref_j"] = aa_list[ref_j]
            pi["mi"] = mi
            all_pis.append(pi)

        essential = [p for p in pis if p["is_essential"]]
        all_essential.extend(essential)

        # Build inference rules
        for pi in essential:
            all_inferences.append(
                {
                    "pos_i": pos_i,
                    "pos_j": pos_j,
                    "aa_i": pi["aa_i"],
                    "aa_j": pi["aa_j"],
                    "ref_i": aa_list[ref_i],
                    "ref_j": aa_list[ref_j],
                    "mi": mi,
                }
            )

    # Deduplicate: multiple QM cubes can decode to the same residue pair
    # (differing only in which don't-care cells they cover). Report the
    # distinct-pair count so it matches the markdown generator.
    seen_pis = {}
    deduped_pis = []
    for pi in all_pis:
        key = (pi["pos_i"], pi["pos_j"], pi["aa_i"], pi["aa_j"])
        if key not in seen_pis:
            seen_pis[key] = len(deduped_pis)
            deduped_pis.append(dict(pi))
        else:
            # OR the essential flag: a later cube may be the essential one
            idx = seen_pis[key]
            deduped_pis[idx]["is_essential"] = (
                deduped_pis[idx]["is_essential"] or pi["is_essential"]
            )
    all_pis = deduped_pis
    all_essential = [p for p in all_pis if p["is_essential"]]

    print(f"  Total prime implicants (distinct pairs): {len(all_pis)}")
    print(f"  Essential prime implicants: {len(all_essential)}")

    # 5. Write all inference rules
    print("\n[5/5] Writing inference rules...")

    # Group by position pair
    pair_groups = {}
    for inf in all_inferences:
        key = (inf["pos_i"], inf["pos_j"])
        if key not in pair_groups:
            pair_groups[key] = []
        pair_groups[key].append(inf)

    print("\n" + "=" * 80)
    print("ALL CO-EVOLUTIONARY INFERENCE RULES")
    print("=" * 80)

    rule_count = 0
    for (pos_i, pos_j), rules in sorted(pair_groups.items()):
        print(f"\nPosition pair ({pos_i}, {pos_j}):")
        print(
            f"  Reference: pos {pos_i}={rules[0]['ref_i']}, pos {pos_j}={rules[0]['ref_j']}"
        )
        print(f"  MI: {rules[0]['mi']:.4f}")
        print(f"  Rules ({len(rules)}):")
        for r in rules:
            rule_count += 1
            print(
                f"    Rule {rule_count:3d}: IF pos {pos_i}={r['aa_i']} AND pos {pos_j}={r['aa_j']} "
                f"THEN co-evolutionary (MI={r['mi']:.3f})"
            )

    print(f"\nTotal inference rules: {rule_count}")

    # Save results
    summary = {
        "dataset": f"SARS-CoV-2 Spike Protein ({n_all} sequences)",
        "variable_positions": len(variable_positions),
        "co_evolving_pairs": len(co_evolving),
        "total_prime_implicants": len(all_pis),
        "essential_prime_implicants": len(all_essential),
        "total_inference_rules": rule_count,
        "top_co_evolving_pairs": [
            {
                "pos_i": pi,
                "pos_j": pj,
                "mi": mi,
                "n": n,
                "ref_i": aa_list[ri],
                "ref_j": aa_list[rj],
            }
            for pi, pj, mi, n, ri, rj in co_evolving[:20]
        ],
        "inferences": all_inferences,
    }

    with open(results_dir / "master_boolean_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nResults saved to: {results_dir}")

    # Final summary
    print("\n" + "=" * 80)
    print("MASTER BOOLEAN FUNCTION SUMMARY")
    print("=" * 80)
    print(f"Sequences: {n_all}")
    print(f"Variable positions: {len(variable_positions)}")
    print(f"Co-evolutionary pairs: {len(co_evolving)}")
    print(f"Total prime implicants: {len(all_pis)}")
    print(f"Essential prime implicants: {len(all_essential)}")
    print(f"Total inference rules: {rule_count}")
    print(f"\nThe master Boolean function:")
    print(f"  f(pos_i, pos_j, aa_i, aa_j) = OR of all essential prime implicants")
    print(f"  Each PI is an AND of residue conditions at two positions")
    print(f"  The function returns 1 (co-evolutionary) when ANY PI matches")




    # ============================================================
    # COMBINED MI + PERPLEXITY ANALYSIS (all experiments)
    # ============================================================
    print("\n=== Combined MI + Perplexity Analysis ===")
    try:
        from coevolution_shared import (
            combined_pair_scores, compute_entropy_vectorized,
            load_position_arrays as _lpa,
        )
        _pa, _na, _fl = _lpa(max_pos=None, aligned=True)
        _ent = compute_entropy_vectorized(_pa, _na, _fl)
        _var = [p for p in range(_fl) if _ent[p] > 0.3]
        _pairs = [(i, j) for idx, i in enumerate(_var)
                  for j in _var[idx + 1:] if j - i <= 30]
        _scored = combined_pair_scores(_pa, _pairs, _na, _ent)
        print(f"  Variable positions: {len(_var)}")
        print(f"  Pairs scored (MI + perplexity ratio): {len(_scored)}")
        print(f"  Top 5 combined (MI + ratio):")
        for _s in _scored[:5]:
            print(f"    ({_s['pos_i']},{_s['pos_j']}): MI={_s['mi']:.3f} "
                  f"ratio={_s['ratio']:.2f} combined={_s['combined']:.3f}")
        _mi_top = sorted(_scored, key=lambda s: -s['mi'])[:5]
        print(f"  Top 5 by MI alone:")
        for _s in _mi_top:
            print(f"    ({_s['pos_i']},{_s['pos_j']}): MI={_s['mi']:.3f} "
                  f"ratio={_s['ratio']:.2f}")
        # ranking agreement
        _r_mi = {(_s['pos_i'], _s['pos_j']): idx
                 for idx, _s in enumerate(sorted(_scored, key=lambda s: -s['mi']))}
        _r_cb = {(_s['pos_i'], _s['pos_j']): idx
                 for idx, _s in enumerate(_scored)}
        _same = sum(1 for k in _r_mi if _r_mi[k] == _r_cb[k])
        print(f"  Ranking agreement (MI vs combined, top-5 same): "
              f"{len([k for k in _r_mi if k in _r_cb and _r_mi[k] < 5 and _r_cb[k] < 5])}/5")
    except Exception as _e:
        print(f"  Combined analysis skipped: {_e}")


if __name__ == "__main__":
    main()
