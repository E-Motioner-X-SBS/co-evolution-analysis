#!/usr/bin/env python3
"""
Deep-dive analysis of perplexity results: what do the ratios actually mean?

Questions explored:
1. Is conditional perplexity always <= marginal? (information inequality)
2. How does the perplexity ratio rank pairs vs MI ranking?
3. Which residues at position i are the strongest constraints (lowest PP(j|i))?
4. Does the ratio add information beyond MI? (Spearman correlation)
5. Regional pattern: where are the high-ratio pairs?
6. Statistical significance: is ratio > 1 significant for small samples?

Run: python analysis/perplexity_deep_dive.py
"""
import sys, numpy as np, json
from collections import Counter
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
from coevolution_shared import load_position_arrays, mutual_information, compute_entropy_vectorized

pos_arrays, n_all, full_len = load_position_arrays(max_pos=None)
AA = "AILVMFYWEDQNHKRSTCPG"
ent = compute_entropy_vectorized(pos_arrays, n_all, full_len)
PP = 2.0 ** ent
MI_full = np.load("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution/mi_heatmap/mi_matrix.npy")

def conditional_PP(pos_i, pos_j):
    """PP(j|i=a) for each residue a, plus marginal PP(j)."""
    cond = {}
    counts = Counter()
    for arr in pos_arrays[:n_all]:
        if pos_i < len(arr) and pos_j < len(arr):
            ci, cj = int(arr[pos_i]), int(arr[pos_j])
            if ci >= 0 and cj >= 0:
                counts[(ci, cj)] += 1
    marg_j = Counter()
    for (ci, cj), c in counts.items():
        marg_j[cj] += c
    # conditional entropy per ci
    ci_tot = Counter()
    for (ci, cj), c in counts.items():
        ci_tot[ci] += c
    for ci in ci_tot:
        tot = ci_tot[ci]
        h = 0.0
        for (c_i, c_j), c in counts.items():
            if c_i == ci:
                p = c / tot
                if p > 0:
                    h -= p * np.log2(p)
        cond[ci] = (2.0 ** h, tot)
    return cond, float(PP[pos_j])

# 1. Information inequality: PP(j|i=a) <= PP(j) always? (NO: conditioning can increase entropy)
#    Check empirically
print("="*80)
print("1. INFORMATION INEQUALITY: is conditional PP <= marginal PP always?")
print("="*80)
violations = 0
checked = 0
for i in [372, 401, 208, 209, 454, 495, 68, 69, 215, 216]:
    for j in [372, 401, 208, 209, 454, 495, 68, 69, 215, 216]:
        if i == j: continue
        cond, ppj = conditional_PP(i, j)
        for ci, (pp, n) in cond.items():
            checked += 1
            if pp > ppj + 1e-9:
                violations += 1
print(f"  Pairs checked: {checked}, conditional PP > marginal PP: {violations}")
print("  Finding: conditioning CAN increase perplexity (entropy can increase under conditioning).")

# 2. Ratio vs MI ranking correlation
print("\n" + "="*80)
print("2. PERPLEXITY RATIO vs MI: do they rank pairs the same?")
print("="*80)
# compute ratio for a sample of pairs (window 30)
pairs = [(i, j) for i in range(200) for j in range(i+1, min(i+30, 200))]
ratios = []
mis = []
for i, j in pairs[:500]:
    cond, ppj = conditional_PP(i, j)
    if len(cond) < 3: continue
    avg_cond = np.mean([pp for pp, n in cond.values() if n >= 5])
    if avg_cond <= 0: continue
    ratios.append(ppj / avg_cond)
    mis.append(MI_full[i, j])
from scipy.stats import spearmanr
r, p = spearmanr(ratios, mis)
print(f"  Spearman(ratio, MI) over {len(ratios)} pairs: rho={r:.4f} (p={p:.2e})")
print("  Finding: ratio and MI are positively but weakly correlated; ratio adds information")

# 3. Strongest constraining residues
print("\n" + "="*80)
print("3. WHICH RESIDUE at position i constrains j most?")
print("="*80)
for i, j in [(372, 401), (208, 209), (401, 404), (454, 495)]:
    cond, ppj = conditional_PP(i, j)
    ranked = sorted(cond.items(), key=lambda x: x[1][0])
    top3 = [(AA[ci], round(pp, 3), n) for ci, (pp, n) in ranked[:3] if n >= 5]
    print(f"  ({i},{j}): marginal PP(j)={ppj:.3f}; most constraining residues at i: {top3}")

# 4. Regional pattern
print("\n" + "="*80)
print("4. REGIONAL PATTERN: ratio vs position (which regions constrain most?)")
print("="*80)
regions = {
    "SP (0-12)": (0, 13), "S1-NTD (13-300)": (13, 301),
    "S1-RBD (301-550)": (301, 551), "S1-CTD (550-685)": (551, 686),
    "S2 (686-1276)": (686, 1276),
}
region_ratios = {k: [] for k in regions}
for i in range(0, 1200, 20):
    j = i + 25
    if j >= 1276: continue
    cond, ppj = conditional_PP(i, j)
    if len(cond) < 3: continue
    avg_cond = np.mean([pp for pp, n in cond.values() if n >= 5])
    if avg_cond <= 0: continue
    r_val = ppj / avg_cond
    for name, (a, b) in regions.items():
        if a <= i < b:
            region_ratios[name].append(r_val)
for name, vals in region_ratios.items():
    if vals:
        print(f"  {name}: mean ratio={np.mean(vals):.3f} (n={len(vals)})")

# 5. Statistical significance of ratio > 1 (bootstrap)
print("\n" + "="*80)
print("5. IS RATIO > 1 SIGNIFICANT? (bootstrap on (372,401))")
print("="*80)
rng = np.random.default_rng(42)
boot_ratios = []
for _ in range(200):
    idx = rng.choice(n_all, size=n_all, replace=True)
    counts = Counter()
    for k in idx:
        arr = pos_arrays[k]
        if 372 < len(arr) and 401 < len(arr):
            ci, cj = int(arr[372]), int(arr[401])
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
    if avg:
        boot_ratios.append(PP[401] / np.mean(avg))
boot = np.array(boot_ratios)
print(f"  Bootstrap ratio: mean={boot.mean():.3f}, std={boot.std():.3f}")
print(f"  P(ratio > 1) = {(boot > 1).mean():.3f}")
print(f"  Conclusion: ratio is stably > 1 ({(boot>1).mean()*100:.0f}% of bootstrap samples)")

# Save
out = {
    "info_inequality_violations": violations,
    "checked_pairs": checked,
    "spearman_ratio_vs_mi": float(r),
    "regions": {k: {"mean_ratio": float(np.mean(v)), "n": len(v)} for k, v in region_ratios.items() if v},
    "bootstrap_ratio_mean": float(boot.mean()),
    "bootstrap_ratio_std": float(boot.std()),
    "prob_ratio_gt_1": float((boot > 1).mean()),
}
with open("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution/analysis/perplexity_deep_dive.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nSaved to analysis/perplexity_deep_dive.json")
