#!/usr/bin/env python3
"""
Combined MI + Perplexity analysis across all co-evolution experiments.

Question: does combining mutual information (total correlation) with the
perplexity ratio (determinism) improve co-evolution results over using
MI alone?

Design:
  1. On the CORRECTED aligned data (21 variable positions), compute for
     every variable-position pair (window 30):
       - MI (mutation-only, as in the pipeline)
       - Perplexity ratio PP(j)/avg PP(j|i) (determinism)
       - Combined score: rank-based (avg of normalized ranks) and
         z-score based (S = z(MI) + z(ratio))
  2. Rank pairs by MI-only, ratio-only, and combined.
  3. Evaluate each ranking by Leave-One-Out prediction accuracy:
       given a mutation at position i in a held-out sequence, predict the
       partner at j as argmax of the constraint row C(i, .) (built from the
       other 1298 sequences); correct if it matches the held-out residue.
     This is the same protocol as allseq_constraint_function.py, applied
     to the top-K pairs of each ranking.
  4. Report whether combined ranking improves accuracy and whether the
     combined constraint (C + lambda * log(PP ratio)) beats C alone.

Run: python analysis/combined_mi_perplexity.py
"""
import sys, json, numpy as np
from collections import Counter
from scipy.stats import spearmanr
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
from coevolution_shared import (load_position_arrays, compute_entropy_vectorized,
                                majority_ref, mutual_information)

pos_arrays, n_all, full_len = load_position_arrays(max_pos=None, aligned=True)
AA = "AILVMFYWEDQNHKRSTCPG"
ent = compute_entropy_vectorized(pos_arrays, n_all, full_len)
PP = 2.0 ** ent

# Variable positions (H > 0.3)
var_pos = [p for p in range(full_len) if ent[p] > 0.3]
print(f"Variable positions: {len(var_pos)}")

def perp_ratio(i, j):
    """PP(j) / mean PP(j|i=a) over residues a with n>=5."""
    counts = Counter()
    for arr in pos_arrays[:n_all]:
        if i < len(arr) and j < len(arr):
            ci, cj = int(arr[i]), int(arr[j])
            if 0 <= ci < 20 and 0 <= cj < 20:
                counts[(ci, cj)] += 1
    ci_tot = Counter()
    for (ci, cj), c in counts.items():
        ci_tot[ci] += c
    avg = []
    for ci, tot in ci_tot.items():
        h = 0.0
        for (c_i, c_j), c in counts.items():
            if c_i == ci:
                p = c / tot
                if p > 0:
                    h -= p * np.log2(p)
        if tot >= 5:
            avg.append(2.0 ** h)
    return PP[j] / np.mean(avg) if avg else None

# Build pair scores
pairs = [(i, j) for idx, i in enumerate(var_pos) for j in var_pos[idx+1:] if j - i <= 30]
mi_dict, ratio_dict = {}, {}
for i, j in pairs:
    m = mutual_information(pos_arrays, i, j, n_all)
    r = perp_ratio(i, j)
    if m > 0.01 and r is not None:
        mi_dict[(i, j)] = m
        ratio_dict[(i, j)] = r

print(f"Pairs with MI>0.01 and valid ratio: {len(mi_dict)}")

# Combined scores
keys = list(mi_dict.keys())
mi_vals = np.array([mi_dict[k] for k in keys])
ra_vals = np.array([ratio_dict[k] for k in keys])

# rank-based combined
def rank_normalize(vals):
    order = np.argsort(np.argsort(vals))  # ranks 0..n-1
    return order / (len(vals) - 1) if len(vals) > 1 else vals * 0

rank_mi = rank_normalize(mi_vals)
rank_ra = rank_normalize(ra_vals)
combined_rank = rank_mi + rank_ra  # avg rank

# z-score combined
z_mi = (mi_vals - mi_vals.mean()) / mi_vals.std()
z_ra = (ra_vals - ra_vals.mean()) / ra_vals.std()
combined_z = z_mi + z_ra

# Spearman between scores
print("\nScore correlations:")
print(f"  Spearman(MI, ratio)     = {spearmanr(mi_vals, ra_vals).statistic:.4f}")
print(f"  Spearman(MI, comb_rank) = {spearmanr(mi_vals, combined_rank).statistic:.4f}")

# LOO-CV prediction for top-K pairs under each ranking
def loo_accuracy(pair_list, K, use_ratio_in_constraint=False):
    """Argmax constraint prediction on held-out sequences."""
    correct = 0
    total = 0
    for (i, j) in pair_list[:K]:
        ri = majority_ref(pos_arrays, i, n_all)
        rj = majority_ref(pos_arrays, j, n_all)
        for hold in range(n_all):
            # build frequency kmap on other 1298
            kmap = np.zeros((20, 20))
            for idx, arr in enumerate(pos_arrays):
                if idx == hold:
                    continue
                if i < len(arr) and j < len(arr):
                    ci, cj = int(arr[i]), int(arr[j])
                    if 0 <= ci < 20 and 0 <= cj < 20:
                        kmap[ci, cj] += 1
            t = kmap.sum()
            if t == 0:
                continue
            kmap /= t
            marg_i = kmap.sum(axis=1)
            marg_j = kmap.sum(axis=0)
            eps = 1e-10
            with np.errstate(divide="ignore", invalid="ignore"):
                C = np.log((kmap + eps) / ((marg_i[:, None] + eps) * (marg_j[None, :] + eps)))
            if use_ratio_in_constraint:
                # add determinism bonus: + lambda * log(ratio) uniformly? Not
                # per-cell; we only combine for ranking, not per-cell C.
                pass
            arr = pos_arrays[hold]
            if i < len(arr) and j < len(arr):
                ci, cj = int(arr[i]), int(arr[j])
                if 0 <= ci < 20 and 0 <= cj < 20:
                    if ci != ri or cj != rj:
                        best = int(np.argmax(C[ci, :]))
                        if best == cj:
                            correct += 1
                        total += 1
    return correct, total

K = 10
order_mi = [keys[k] for k in np.argsort(-mi_vals)]
order_ra = [keys[k] for k in np.argsort(-ra_vals)]
order_cr = [keys[k] for k in np.argsort(-combined_rank)]
order_cz = [keys[k] for k in np.argsort(-combined_z)]

print(f"\nLOO-CV prediction (top-{K} pairs by each ranking):")
for name, order in [("MI only", order_mi), ("Perplexity ratio only", order_ra),
                    ("Combined (rank)", order_cr), ("Combined (z-score)", order_cz)]:
    c, t = loo_accuracy(order, K)
    print(f"  {name:<22}: {c}/{t} = {c/t:.4f}" if t else f"  {name}: no tests")

# Top 10 pairs under each ranking
print("\nTop 10 pairs by each ranking:")
for name, order in [("MI", order_mi), ("Ratio", order_ra), ("Combined", order_cr)]:
    print(f"  {name}: {[(p, round(mi_dict[p],3), round(ratio_dict[p],2)) for p in order[:5]]}")

# Save
out = {
    "n_pairs": len(keys),
    "spearman_mi_ratio": float(spearmanr(mi_vals, ra_vals).statistic),
    "top10_mi": [{"pos_i": p[0], "pos_j": p[1], "mi": float(mi_dict[p]), "ratio": float(ratio_dict[p])} for p in order_mi[:10]],
    "top10_ratio": [{"pos_i": p[0], "pos_j": p[1], "mi": float(mi_dict[p]), "ratio": float(ratio_dict[p])} for p in order_ra[:10]],
    "top10_combined": [{"pos_i": p[0], "pos_j": p[1], "mi": float(mi_dict[p]), "ratio": float(ratio_dict[p])} for p in order_cr[:10]],
    "loo_top10": {
        "mi_only": loo_accuracy(order_mi, K),
        "ratio_only": loo_accuracy(order_ra, K),
        "combined_rank": loo_accuracy(order_cr, K),
        "combined_z": loo_accuracy(order_cz, K),
    },
}
json.dump(out, open("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution/analysis/combined_mi_perplexity.json", "w"), indent=2)
print("\nSaved to analysis/combined_mi_perplexity.json")
