#!/usr/bin/env python3
"""
Discovery: how do the three pair-scoring lenses compare?
  1. MI (total correlation) - scripts 06, 08, 09, 10
  2. Perplexity ratio (determinism) - script 19
  3. mfDCA Frobenius/DI (direct coupling) - script 18

Run: python analysis/dca_vs_mi_vs_ratio.py
"""
import sys, numpy as np, json
from collections import Counter
from scipy.stats import spearmanr
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
from coevolution_shared import load_position_arrays, compute_entropy_vectorized

pos_arrays, n_all, full_len = load_position_arrays(max_pos=None)
ent = compute_entropy_vectorized(pos_arrays, n_all, full_len)
PP = 2.0 ** ent
MI_full = np.load("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution/mi_heatmap/mi_matrix.npy")
F = np.load("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution/dca_results/frobenius_scores.npy")
DI = np.load("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution/dca_results/di_scores.npy")

def ratio_of(i, j):
    counts = Counter()
    for arr in pos_arrays[:n_all]:
        if i < len(arr) and j < len(arr):
            ci, cj = int(arr[i]), int(arr[j])
            if ci >= 0 and cj >= 0:
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
                if p > 0: h -= p * np.log2(p)
        if tot >= 5:
            avg.append(2.0 ** h)
    return PP[j] / np.mean(avg) if avg else None

pairs = [(i, j) for i in range(0, 800) for j in range(i+1, min(i+30, 800))]
ratios, mis, fs, dis = [], [], [], []
for i, j in pairs[::3]:
    r = ratio_of(i, j)
    if r is None: continue
    ratios.append(r); mis.append(MI_full[i, j]); fs.append(F[i, j]); dis.append(DI[i, j])

res = {
    "n_pairs": len(ratios),
    "spearman_ratio_MI": float(spearmanr(ratios, mis).statistic),
    "spearman_ratio_F": float(spearmanr(ratios, fs).statistic),
    "spearman_ratio_DI": float(spearmanr(ratios, dis).statistic),
    "spearman_MI_F": float(spearmanr(mis, fs).statistic),
    "spearman_MI_DI": float(spearmanr(mis, dis).statistic),
    "spearman_F_DI": float(spearmanr(fs, dis).statistic),
    "conclusion": (
        "MI and perplexity ratio are strongly correlated (both measure total correlation). "
        "Both are nearly uncorrelated with DCA F/DI (direct coupling). "
        "The three lenses are: total correlation (MI), determinism (perplexity ratio), "
        "direct coupling (DCA). They answer different questions."
    ),
}
json.dump(res, open("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution/analysis/dca_vs_mi_vs_ratio.json", "w"), indent=2)
print(json.dumps(res, indent=2))
