#!/usr/bin/env python3
"""
All-Sequence Constraint Function for Co-evolution
==================================================

Builds the constraint function C(a, b) = ln(P(a,b) / (P(a) P(b))) on all
sequences at full length and tests it by leave-one-out cross-validation: for
each sequence carrying a mutation at a co-evolving pair (relative to the
majority-residue reference), predict its residue at j from its residue at i
using the other N-1 sequences.

Prediction = argmax of C over the partners that co-occur with the residue at
i (coevolution_shared.predict_partner).  The earlier argmax over all 20
residues was dominated by the eps smoothing: a never-observed partner scores
ln(1/P(a)) > 0 and won almost every argmax, which is where the old 9.24%
came from.  The LOO is exact and O(N) per pair (one count is removed from the
pooled table instead of rebuilding it N times).

Also reported: the conditional-frequency baseline argmax P(b | a), and both
accuracies on the de-duplicated alignment (identical sequences make plain LOO
optimistic: a held-out copy is still in the training data).
"""

from kmap_paths import NARY_SRC, REPO_ROOT  # noqa: E402

import sys
import json
import numpy as np
from pathlib import Path

sys.path.insert(0, str(NARY_SRC))

from nkmap.encoding.bio_sequences import Base20AminoEncoder, AMINO_HE_2012


from coevolution_shared import parse_fasta  # single source of truth


from coevolution_shared import compute_entropy  # single source of truth



def main():
    base_dir = REPO_ROOT
    fasta_file = Path(__import__("os").environ.get("COEVO_FASTA") or (base_dir / "Spike_protein.aln-fasta"))
    results_dir = Path(__import__("os").environ.get("COEVO_RESULTS") or (base_dir / "allseq_constraint_results"))
    results_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("All-Sequence Constraint Function (ALL sequences, full length)")
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
    from coevolution_shared import find_coevolving_pairs_gpu
    co_evolving = find_coevolving_pairs_gpu(
        pos_arrays, variable_positions, n_all, max_gap=30, min_mi=0.1
    )
    print(f"  Co-evolutionary pairs: {len(co_evolving)}")

    # Leave-one-out cross-validation (exact, O(N) per pair)
    print("\n[4/4] Leave-one-out cross-validation...")
    from coevolution_shared import loo_partner_accuracy

    dense = np.stack(pos_arrays)
    unique = np.unique(dense, axis=0)
    all_results = []
    unique_totals = {"n_mutations": 0, "correct": 0, "correct_conditional": 0}

    for pos_i, pos_j, mi, n_muts, ref_i, ref_j in co_evolving[:10]:
        ref_i_code = int(ref_i)  # integer He-2012 codes from the majority ref
        ref_j_code = int(ref_j)
        ref_aa_i = aa_list[ref_i_code] if 0 <= ref_i_code < 20 else "?"
        ref_aa_j = aa_list[ref_j_code] if 0 <= ref_j_code < 20 else "?"

        r = loo_partner_accuracy(dense[:, pos_i], dense[:, pos_j], ref_i_code, ref_j_code)
        u = loo_partner_accuracy(unique[:, pos_i], unique[:, pos_j], ref_i_code, ref_j_code)
        for k in unique_totals:
            unique_totals[k] += u[k]
        total = r["n_mutations"]
        accuracy = r["correct"] / total if total > 0 else 0
        all_results.append(
            {
                "pos_i": pos_i,
                "pos_j": pos_j,
                "mi": mi,
                "ref_i": ref_aa_i,
                "ref_j": ref_aa_j,
                "n_mutations": total,
                "correct": r["correct"],
                "accuracy": accuracy,
                "correct_conditional": r["correct_conditional"],
                "abstained": r["abstained"],
                "unique_sequences": u,
            }
        )
        print(
            f"  ({pos_i:2d},{pos_j:2d}): MI={mi:.4f}, accuracy={accuracy:.4f} "
            f"({r['correct']}/{total}), conditional baseline {r['correct_conditional']}/{total}"
        )

    # Overall accuracy
    total_correct = sum(r["correct"] for r in all_results)
    total_test = sum(r["n_mutations"] for r in all_results)
    overall_acc = total_correct / total_test if total_test > 0 else 0
    cond_correct = sum(r["correct_conditional"] for r in all_results)
    cond_acc = cond_correct / total_test if total_test > 0 else 0
    un = unique_totals["n_mutations"]

    print(
        f"\nOverall LOO-CV accuracy: {overall_acc:.4f} ({total_correct}/{total_test}); "
        f"conditional baseline {cond_acc:.4f}"
    )
    print(
        f"De-duplicated ({len(unique)} unique sequences): constraint "
        f"{unique_totals['correct'] / un if un else 0:.4f}, conditional "
        f"{unique_totals['correct_conditional'] / un if un else 0:.4f} (n={un})"
    )

    # Save results
    summary = {
        "dataset": f"{n_all} Omicron sequences",
        "method": "Leave-One-Out Cross-Validation",
        "variable_positions": len(variable_positions),
        "co_evolving_pairs": len(co_evolving),
        "predictor": "argmax_b C(a,b) over partners b co-occurring with a (LOO, exact)",
        "overall_accuracy": overall_acc,
        "total_correct": total_correct,
        "total_test": total_test,
        "conditional_baseline_accuracy": cond_acc,
        "n_unique_sequences": int(len(unique)),
        "unique_sequences": {
            **unique_totals,
            "accuracy": unique_totals["correct"] / un if un else 0,
            "conditional_accuracy": unique_totals["correct_conditional"] / un if un else 0,
        },
        "results": all_results,
    }

    with open(results_dir / "allseq_constraint_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nResults saved to: {results_dir}")




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
