#!/usr/bin/env python3
"""
Justification for the top-15 position-pair cutoff in master_boolean.py.

Question: does the top-15 cutoff in master_boolean.py drop anything?  This
script counts the co-evolving pairs (mutation-only MI > 0.1, window 30) and
the prime implicants / essential rules as a function of K, on the aligned
alignment (gaps kept as state 20) with padded 32x32 maps, exactly as
master_boolean.py builds them.  (An earlier version stripped gaps and fed
unpadded 20x20 maps, i.e. both corrected defects A1 and A2.)

Run: python analysis/justify_top15.py
"""

import sys as _sys_kp
from pathlib import Path as _Path_kp
_sys_kp.path.insert(0, str(_Path_kp(__file__).resolve().parents[1]))
from kmap_paths import KMAP_SBM_SRC, NARY_SRC, REPO_ROOT  # noqa: E402
import sys, numpy as np, json
sys.path.insert(0, str(NARY_SRC))
sys.path.insert(0, str(KMAP_SBM_SRC))
sys.path.insert(0, str(REPO_ROOT))
from nkmap.encoding.bio_sequences import Base20AminoEncoder, AMINO_HE_2012
from kmap_sbm.analysis.prime_implicants import boolean_minimize_kmap
from coevolution_shared import find_coevolving_pairs_gpu
from collections import Counter

enc = Base20AminoEncoder(version=1)
FASTA = str(REPO_ROOT / "Spike_protein.aln-fasta")
seqs = []
h = None; s = []
for line in open(FASTA):
    line = line.strip()
    if line.startswith(">"):
        if h: seqs.append((h, "".join(s)))
        h = line[1:]; s = []
    elif line: s.append(line.upper())
if h: seqs.append((h, "".join(s)))

max_pos = len(seqs[0][1])
pos_arrays = [np.array([enc.encode.get(c, 20) for c in seq[:max_pos]], dtype=np.int32)
              for _, seq in seqs]   # aligned: gap/unknown = 20
n_all = len(seqs)
aa_list = list(AMINO_HE_2012)

vp = []
for p in range(max_pos):
    cnt = Counter(int(a[p]) for a in pos_arrays if 0 <= a[p] < 20)
    t = sum(cnt.values())
    if t == 0: continue
    he = -sum((c/t)*np.log2(c/t) for c in cnt.values() if c > 0)
    if he > 0.3: vp.append(p)

co_evolving = find_coevolving_pairs_gpu(pos_arrays, vp, n_all, max_gap=30, min_mi=0.1)

import master_boolean as mb  # noqa: E402

results = {}
for K in [5, 10, 15, 20]:
    total_pis = 0; total_ess = 0
    for pos_i, pos_j, mi, n_mut, ref_i, ref_j in co_evolving[:K]:
        kmap = mb.build_mutation_kmap(pos_arrays, pos_i, pos_j, int(ref_i), int(ref_j), n_all)
        res = boolean_minimize_kmap(kmap.flatten().astype(int), algorithm="qm")
        total_pis += res["n_prime_implicants"]
        total_ess += res["n_essential"]
    results[K] = {"pairs_used": min(K, len(co_evolving)), "total_pis": total_pis, "essential": total_ess}

pairs_beyond_15 = [
    {"rank": i, "pos_i": int(p[0]), "pos_j": int(p[1]), "mi": round(float(p[2]), 4), "n_mut": int(p[3])}
    for i, p in enumerate(co_evolving[15:], start=16)
]

out = {
    "question": "Does the top-15 cutoff in master_boolean.py drop co-evolving pairs?",
    "n_coevolving_pairs": len(co_evolving),
    "answer": (
        f"{len(co_evolving)} pairs pass mutation-only MI > 0.1 within window 30, so the "
        f"top-15 cutoff keeps {'all of them' if len(co_evolving) <= 15 else 'the 15 strongest'}."
    ),
    "rules_by_K": {str(k): v for k, v in results.items()},
    "pairs_beyond_rank_15": pairs_beyond_15,
}
with open(str(REPO_ROOT / "analysis/top15_justification.json"), "w") as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))
