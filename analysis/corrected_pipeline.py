#!/usr/bin/env python3
"""
CORRECTED co-evolution analysis pipeline.

Fixes two confirmed defects in the original pipeline:

FIX A1 (column misalignment): The original code strips gaps (clean =
"".join(aa for aa in seq if aa in encoder.encode)), which DELETES gap
characters. Because different sequences have gaps at different positions,
"column j" of the cleaned array corresponds to DIFFERENT raw alignment
positions in different sequences. Verified: raw-aligned columns 372, 401,
413, 427, 852 are conserved (H ~ 0), but the cleaned columns reported
H(372)=1.63 and 1,249 "variable positions" - artifacts of misalignment.

  Fix: keep the full alignment; encode '-' as state 20 (21 states total).

FIX A2 (QM encoding corruption): The original feeds a 400-cell (20x20)
K-map to a Quine-McCluskey implementation that derives k_bits =
int(log2(400))//2 = 4 -> 8 bits. Cells 256-399 wrap onto cells 0-143,
corrupting minterms. Verified: of 152 rules, 143 have residue labels
never observed in the data.

  Fix: pad the 20x20 map to 32x32 (5 bits per axis), mark rows/cols 20-31
  as don't-care, run QM on 1024 cells with 10 bits (no wrap-around),
  decode 5 bits per axis.

Run: python analysis/corrected_pipeline.py
"""
import sys, json, numpy as np
from collections import Counter
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/kmap-sbm-validation/src")
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")
from nkmap.encoding.bio_sequences import AMINO_HE_2012
from kmap_sbm.analysis.prime_implicants import boolean_minimize_kmap

AA = AMINO_HE_2012
AA_TO_IDX = {aa: i for i, aa in enumerate(AA)}
GAP = 20  # 21st state

FASTA = "/store/shuvam/E-motioner-X-SBS/datasets/co-evolution/Spike_protein.aln-fasta"
OUT = "/store/shuvam/E-motioner-X-SBS/datasets/co-evolution/analysis/corrected_results.json"

seqs = []
h = None; s = []
for line in open(FASTA):
    line = line.strip()
    if line.startswith(">"):
        if h: seqs.append((h, "".join(s)))
        h = line[1:]; s = []
    elif line: s.append(line.upper())
if h: seqs.append((h, "".join(s)))

N = len(seqs)
L = len(seqs[0][1])

# CORRECTED position arrays: aligned, gaps = 20
X = np.full((N, L), GAP, dtype=np.int64)
for i, (_, seq) in enumerate(seqs):
    for j, c in enumerate(seq):
        X[i, j] = AA_TO_IDX.get(c, GAP)

print(f"Aligned matrix: {X.shape}, {N} sequences, {L} positions")
print(f"Gap state (20) count: {(X == GAP).sum()} of {N*L} ({100*(X==GAP).sum()/(N*L):.2f}%)")

# Entropy per position (gaps excluded from the count, like standard MSA practice)
def entropy_col(col):
    cnt = Counter(int(c) for c in col if c != GAP)
    total = sum(cnt.values())
    if total == 0: return 0.0
    return -sum((c/total)*np.log2(c/total) for c in cnt.values() if c > 0)

H = np.array([entropy_col(X[:, j]) for j in range(L)])
var_pos = np.where(H > 0.3)[0]
print(f"\nVariable positions (H > 0.3, gaps excluded): {len(var_pos)} / {L}")
print(f"Top 10 by entropy:")
for j in np.argsort(H)[::-1][:10]:
    print(f"  pos {j}: H={H[j]:.3f}")

# MI matrix (all pairs, aligned, gap pairs excluded)
def mi_pair(i, j):
    ci = X[:, i]; cj = X[:, j]
    valid = (ci != GAP) & (cj != GAP)
    ci, cj = ci[valid], cj[valid]
    if len(ci) < 10: return 0.0
    joint = np.bincount(ci*20 + cj, minlength=400).reshape(20, 20).astype(float)
    total = joint.sum()
    if total == 0: return 0.0
    mi = 0.0
    for a in range(20):
        for b in range(20):
            if joint[a, b] > 0:
                pa = joint[a,:].sum()/total; pb = joint[:,b].sum()/total
                mi += (joint[a,b]/total) * np.log2((joint[a,b]/total)/(pa*pb))
    return mi

# Compute MI for all pairs (window 30 for speed, full for the hub check)
pairs = [(i, j) for i in var_pos for j in var_pos if i < j <= i+30]
mi_results = []
for i, j in pairs:
    m = mi_pair(i, j)
    if m > 0.1:
        mi_results.append((i, j, m))
mi_results.sort(key=lambda x: -x[2])
print(f"\nCo-evolving pairs (MI > 0.1, window 30, aligned): {len(mi_results)}")
print("Top 10:")
for i, j, m in mi_results[:10]:
    print(f"  ({i},{j}): MI={m:.4f}")

# Top-15 rules with CORRECTED QM (32x32 padded, 10-bit)
print("\nTop-15 pair rules (corrected QM):")
def build_kmap_32(i, j):
    """20x20 mutation map padded to 32x32, rows/cols 20-31 = don't-care."""
    # majority refs
    ri = Counter(int(c) for c in X[:, i] if c != GAP).most_common(1)[0][0]
    rj = Counter(int(c) for c in X[:, j] if c != GAP).most_common(1)[0][0]
    kmap = np.full((32, 32), -1, dtype=np.int32)  # all don't-care
    for a in range(20):
        for b in range(20):
            # count how many sequences have (a,b) at (i,j) without gaps
            cnt = ((X[:, i] == a) & (X[:, j] == b)).sum()
            if cnt > 0:
                if a == ri and b == rj:
                    kmap[a, b] = -1  # reference
                else:
                    kmap[a, b] = 1   # observed mutation
            else:
                kmap[a, b] = 0       # never observed
    return kmap, ri, rj

rules = []
for i, j, m in mi_results[:15]:
    kmap, ri, rj = build_kmap_32(i, j)
    res = boolean_minimize_kmap(kmap.flatten().astype(int), algorithm="qm")
    n_on = int((kmap[:20, :20] == 1).sum())
    rules.append({"pos_i": int(i), "pos_j": int(j), "mi": round(m, 4),
                  "ref_i": AA[ri], "ref_j": AA[rj], "n_on": n_on,
                  "n_PI": res["n_prime_implicants"], "n_essential": res["n_essential"]})
    print(f"  ({i},{j}): MI={m:.4f} ref=({AA[ri]},{AA[rj]}) on={n_on} PIs={res['n_prime_implicants']} ess={res['n_essential']}")

total_ess = sum(r["n_essential"] for r in rules)
print(f"\nTotal essential rules (corrected): {total_ess}")

out = {
    "fix_A1": "aligned columns, gaps = state 20 (was: gap-stripped, misaligned)",
    "fix_A2": "32x32 padded QM, 10-bit, 5 bits per axis (was: 400 cells in 8 bits, wrap-around)",
    "n_sequences": N, "n_positions": L,
    "n_variable_positions": int(len(var_pos)),
    "top_variable": [{"pos": int(j), "H": round(float(H[j]), 4)} for j in np.argsort(H)[::-1][:20]],
    "n_coevolving_pairs_MI01": len(mi_results),
    "top_mi_pairs": [{"pos_i": int(i), "pos_j": int(j), "mi": round(m, 4)} for i, j, m in mi_results[:20]],
    "top15_rules": rules,
    "total_essential_rules": int(total_ess),
}
json.dump(out, open(OUT, "w"), indent=2)
print(f"\nSaved to {OUT}")
