#!/usr/bin/env python3
"""
Flipped Boolean Co-evolution Analysis
======================================

Instead of: "what co-evolves?" (on-set = observed mutations)
We ask:     "which combinations are never observed?" (on-set = unobserved pairs)

With 1,299 sequences (299 unique) "never observed" is absence of evidence,
not proof of incompatibility; each rule is a cube of residue combinations
(sets at i and j) that no sequence carries.
"""

from kmap_paths import KMAP_SBM_SRC, NARY_SRC, REPO_ROOT  # noqa: E402

import sys
import json
import numpy as np
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(NARY_SRC))
sys.path.insert(0, str(KMAP_SBM_SRC))

from nkmap.encoding.bio_sequences import Base20AminoEncoder, AMINO_HE_2012
from kmap_sbm.analysis.prime_implicants import boolean_minimize_kmap


from coevolution_shared import parse_fasta, decode_pair_cube  # single source of truth


from coevolution_shared import compute_entropy  # single source of truth


from coevolution_shared import majority_ref as get_majority_ref  # single source of truth


def main():
    base_dir = REPO_ROOT
    fasta_file = Path(__import__("os").environ.get("COEVO_FASTA") or (base_dir / "Spike_protein.aln-fasta"))
    results_dir = Path(__import__("os").environ.get("COEVO_RESULTS") or (base_dir / "flipped_boolean_results"))
    results_dir.mkdir(parents=True, exist_ok=True)

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
        arr = np.array(
            [encoder.encode.get(aa, 20) for aa in seq[:max_pos]], dtype=np.int32
        )
        pos_arrays.append(arr)

    # Find variable positions
    variable_positions = [
        p for p in range(max_pos) if compute_entropy(pos_arrays, p, n_all) > 0.3
    ]
    print(f"  Variable positions: {len(variable_positions)}")

    # Find co-evolutionary pairs (GPU-accelerated via coevolution_shared)
    print("\n[3/4] Finding co-evolutionary pairs (GPU)...")
    from coevolution_shared import find_coevolving_pairs_gpu

    co_evolving = find_coevolving_pairs_gpu(
        pos_arrays, variable_positions, n_all, max_gap=30, min_mi=0.1
    )
    print(f"  Co-evolutionary pairs: {len(co_evolving)}")

    # Build FLIPPED K-maps and Boolean functions
    print("\n[4/4] Building flipped K-maps...")
    all_rules = []
    rule_num = 0

    for pos_i, pos_j, mi, n_muts, ref_i, ref_j in co_evolving[:15]:
        # BUG FIX: ref_i, ref_j are already int codes (0-19) from majority ref.
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
        # CORRECTED (FIX A2): 32x32 padded, rows/cols 20-31 don't-care.
        # 20x20 region starts ALL FORBIDDEN (1); observed pairs -> 0 (allowed);
        # reference -> -1 (DC). Cells that remain 1 are truly forbidden.
        kmap = np.full((32, 32), -1, dtype=np.int32)  # padding = don't-care
        kmap[:20, :20] = 1  # 20x20 region: all forbidden by default

        # Mark observed pairs as ALLOWED (0), reference pair as DC (-1)
        for arr in pos_arrays[:n_all]:
            if pos_i < len(arr) and pos_j < len(arr):
                ci, cj = int(arr[pos_i]), int(arr[pos_j])
                if 0 <= ci < 20 and 0 <= cj < 20:
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
        essential_keys = {(tuple(q["values"]), tuple(q["mask"]))
                          for q in result["essential_prime_implicants"]}
        for pi in result["prime_implicants"]:
            rows, cols = decode_pair_cube(pi)
            if not rows or not cols:
                continue  # cube lies wholly in the padding
            aa_i = "".join(aa_list[c] for c in rows)
            aa_j = "".join(aa_list[c] for c in cols)
            is_essential = (tuple(pi["values"]), tuple(pi["mask"])) in essential_keys
            marker = "*" if is_essential else " "
            rule_num += 1
            print(
                f"    {marker} Rule {rule_num}: IF pos {pos_i} in {{{aa_i}}} AND pos {pos_j} in {{{aa_j}}} THEN never observed"
            )

            all_rules.append(
                {
                    "pos_i": pos_i,
                    "pos_j": pos_j,
                    "aa_i": aa_i,
                    "aa_j": aa_j,
                    "n_cells": len(rows) * len(cols),
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
