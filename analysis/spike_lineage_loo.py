#!/usr/bin/env python3
"""Spike partner prediction when the held-out sequence's relatives are held out too.

allseq_constraint_function.py predicts, for every sequence carrying a mutation
at one of the 10 co-evolving pairs, its residue at j from its residue at i,
using the pooled 20x20 table of the other N-1 sequences (99.4%; 98.3% on the
300 unique sequences).  Plain leave-one-out keeps the held-out sequence's
near-identical relatives in the training table, so this measures how well a
lineage predicts itself.

Here every training sequence within Hamming distance h of the held-out one is
removed as well (h = -1 is plain LOO; h = 0 also drops identical copies).
Distances are counted over all alignment columns where both sequences have a
residue, excluding the columns of the 10 pairs themselves.  For each h the
script reports

  coverage     fraction of held-out mutations whose residue at i still
               occurs in the training table (otherwise no prediction)
  accuracy     constraint-function prediction (predict_partner), on the
               covered mutations
  conditional  argmax P(b | a) baseline on the same mutations
  majority     always predicting the training table's most common
               partner of residue a's column j (residue-blind baseline)

Output: analysis/spike_lineage_loo.json
Run:    python analysis/spike_lineage_loo.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO)]

import numpy as np  # noqa: E402

from kmap_paths import NARY_SRC  # noqa: E402

sys.path.insert(0, str(NARY_SRC))
from nkmap.encoding.bio_sequences import Base20AminoEncoder  # noqa: E402

from coevolution_shared import (N_AA, compute_entropy, find_coevolving_pairs_gpu,  # noqa: E402
                                parse_fasta, predict_partner)

RADII = (-1, 0, 1, 2, 3, 5, 8, 12, 20)


def hamming_matrix(dense: np.ndarray, skip_cols) -> np.ndarray:
    keep = np.setdiff1d(np.arange(dense.shape[1]), np.asarray(sorted(skip_cols)))
    X = dense[:, keep]
    res = X < N_AA
    D = np.zeros((len(X), len(X)), dtype=np.int32)
    for s in range(0, len(X), 128):
        blk, rb = X[s:s + 128], res[s:s + 128]
        both = rb[:, None, :] & res[None, :, :]
        D[s:s + 128] = ((blk[:, None, :] != X[None, :, :]) & both).sum(axis=2)
    return D


def evaluate_pair(ci, cj, ref_i, ref_j, D, h):
    ok = (ci < N_AA) & (cj < N_AA)
    full = np.zeros((N_AA, N_AA))
    np.add.at(full, (ci[ok], cj[ok]), 1)
    out = {"n_mutations": 0, "covered": 0, "correct": 0, "correct_conditional": 0,
           "correct_majority": 0, "mean_removed": 0.0}
    for s in np.flatnonzero(ok):
        a, b = int(ci[s]), int(cj[s])
        if a == ref_i and b == ref_j:
            continue
        out["n_mutations"] += 1
        drop = np.flatnonzero(ok & (D[s] <= h)) if h >= 0 else np.array([], int)
        drop = np.union1d(drop, [s])
        counts = full.copy()
        np.add.at(counts, (ci[drop], cj[drop]), -1)
        out["mean_removed"] += len(drop)
        pick, cond = predict_partner(counts, a)
        if pick is None:
            continue
        out["covered"] += 1
        out["correct"] += int(pick == b)
        out["correct_conditional"] += int(cond == b)
        out["correct_majority"] += int(int(np.argmax(counts.sum(axis=0))) == b)
    return out


def main() -> int:
    seqs = parse_fasta(REPO / "Spike_protein.aln-fasta")
    enc = Base20AminoEncoder(version=1)
    L = len(seqs[0][1])
    pos = [np.array([enc.encode.get(c, 20) for c in s[:L]], dtype=np.int32) for _, s in seqs]
    dense = np.stack(pos)
    var = [p for p in range(L) if compute_entropy(pos, p, len(pos)) > 0.3]
    pairs = find_coevolving_pairs_gpu(pos, var, len(pos), max_gap=30, min_mi=0.1)[:10]
    cols = {c for p in pairs for c in p[:2]}
    D = hamming_matrix(dense, cols)
    iu = np.triu_indices(len(D), 1)
    dist_summary = {"median": float(np.median(D[iu])), "p90": float(np.percentile(D[iu], 90)),
                    "max": int(D[iu].max()),
                    "nearest_neighbour_median": float(np.median(
                        np.where(np.eye(len(D), dtype=bool), 10**9, D).min(axis=1)))}
    print("pairwise distance (mismatches):", dist_summary, flush=True)

    by_radius = {}
    for h in RADII:
        tot = {"n_mutations": 0, "covered": 0, "correct": 0, "correct_conditional": 0,
               "correct_majority": 0, "mean_removed": 0.0}
        per_pair = []
        for pos_i, pos_j, mi, n_muts, ref_i, ref_j in pairs:
            r = evaluate_pair(dense[:, pos_i], dense[:, pos_j], int(ref_i), int(ref_j), D, h)
            per_pair.append({"pair": [int(pos_i), int(pos_j)], **r})
            for k in tot:
                tot[k] += r[k]
        n, cov = tot["n_mutations"], tot["covered"]
        by_radius[str(h)] = {
            "n_mutations": n, "coverage": cov / n if n else 0.0,
            "accuracy": tot["correct"] / cov if cov else 0.0,
            "conditional_accuracy": tot["correct_conditional"] / cov if cov else 0.0,
            "majority_accuracy": tot["correct_majority"] / cov if cov else 0.0,
            "accuracy_over_all": tot["correct"] / n if n else 0.0,
            "mean_sequences_removed": tot["mean_removed"] / n if n else 0.0,
            "per_pair": per_pair}
        b = by_radius[str(h)]
        print(f"h={h:3d}: coverage {b['coverage']:.3f}  accuracy {b['accuracy']:.3f}  "
              f"conditional {b['conditional_accuracy']:.3f}  majority {b['majority_accuracy']:.3f}  "
              f"(removed {b['mean_sequences_removed']:.0f})", flush=True)
    out = {"n_sequences": len(dense), "pairs": [[int(p[0]), int(p[1])] for p in pairs],
           "distance": dist_summary, "radii": list(RADII), "by_radius": by_radius}
    (REPO / "analysis" / "spike_lineage_loo.json").write_text(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
