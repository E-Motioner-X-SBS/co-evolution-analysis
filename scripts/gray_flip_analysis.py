#!/usr/bin/env python3
"""
gray_flip_analysis.py — the K-map "bit-flip" consistency test.

The K-map theory places amino acids on a Gray-code hypercube: residues that
co-evolve should be GRAY-ADJACENT (Hamming distance 1) in the encoded space —
a single bit-flip connects them, exactly the "consistency when you flip the
bits" intuition.

Test, for every co-evolving position pair:
  1. collect the OBSERVED residue pairs (aa_i, aa_j) that co-occur (the
     on-set cells of the mutation K-map)
  2. compute their 5-bit Gray Hamming distance (kmap_sbm encoding)
  3. compare the Hamming-distance distribution to the null: random residue
     pairs (uniform over the 20 AA) — expected P(H=1) = 5/31 ≈ 0.161
     (5 neighbors of a 5-bit code among 31 others)
  4. ALSO: correlate Gray-adjacency with 3D contact (the union test:
     is a co-evolving pair BOTH Gray-adjacent AND 3D-close?)

Output: results/contacts/gray_flip_analysis.json
"""
import json
import sys
from collections import Counter
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/kmap-sbm-validation/src")
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")
from kmap_sbm.encoding.gray_amino import encode_gray_single  # noqa: E402
from nkmap.encoding.bio_sequences import AMINO_HE_2012  # noqa: E402

INV = REPO / "results" / "contacts" / "rules_inventory.json"
OUT = REPO / "results" / "contacts" / "gray_flip_analysis.json"

AA = list(AMINO_HE_2012)
AA_IDX = {a: i for i, a in enumerate(AA)}
GRAY = {a: encode_gray_single(a) for a in AA}


def hamming(a, b):
    return bin(a ^ b).count("1")


def main():
    inv = json.load(open(INV))
    residue_rules = inv["residue_rules"]  # 36 on-set residue pairs
    essential = inv["essential_rules"]
    # observed residue pairs per position pair (the on-set cells)
    per_pair = {}
    for r in residue_rules:
        key = (r["pos_i"], r["pos_j"])
        per_pair.setdefault(key, []).append((r["aa_i"], r["aa_j"]))
    for r in essential:
        key = (r["pos_i"], r["pos_j"])
        per_pair.setdefault(key, []).append((r["aa_i"], r["aa_j"]))

    # 1) Gray Hamming distribution of observed co-evolving residue pairs
    obs_h = Counter()
    obs_n = 0
    for pairs in per_pair.values():
        for (ai, aj) in pairs:
            if ai in GRAY and aj in GRAY:
                obs_h[hamming(GRAY[ai], GRAY[aj])] += 1
                obs_n += 1
    # null: all 20x20 residue pairs
    null_h = Counter()
    null_n = 0
    for a in AA:
        for b in AA:
            null_h[hamming(GRAY[a], GRAY[b])] += 1
            null_n += 1
    exp_h1 = null_h[1] / null_n  # 5/31
    obs_h1 = obs_h[1] / obs_n if obs_n else 0.0
    enrich_h1 = obs_h1 / exp_h1 if exp_h1 else float("nan")
    # binomial p for the H=1 count
    from scipy.stats import binomtest
    p_h1 = binomtest(obs_h[1], obs_n, exp_h1, alternative="greater").pvalue if obs_n else 1.0

    print("=== Gray bit-flip consistency of observed co-evolving pairs ===")
    print(f"observed pairs: {obs_n} | H=1 count: {obs_h[1]} "
          f"({obs_h1:.3f}, expected {exp_h1:.3f}, enrich {enrich_h1:.2f}x, p={p_h1:.3g})")
    print("Hamming distribution (observed vs null):")
    for d in sorted(set(obs_h) | set(null_h)):
        print(f"  H={d}: obs {obs_h.get(d, 0)} ({100*obs_h.get(d,0)/obs_n:.1f}%) | "
              f"null {null_h.get(d, 0)} ({100*null_h.get(d,0)/null_n:.1f}%)")

    # 2) per-pair Gray adjacency
    print("\nper position pair (H=1 fraction of its observed residue pairs):")
    per_pair_rows = []
    for key, pairs in sorted(per_pair.items()):
        h1 = sum(1 for (ai, aj) in pairs
                 if ai in GRAY and aj in GRAY and hamming(GRAY[ai], GRAY[aj]) == 1)
        n = sum(1 for (ai, aj) in pairs if ai in GRAY and aj in GRAY)
        per_pair_rows.append({"pair": list(key), "n_observed": n,
                              "gray_adjacent": h1,
                              "frac_h1": round(h1 / n, 3) if n else None})
        print(f"  {str(key):>10}: {h1}/{n} Gray-adjacent ({h1/n if n else 0:.3f})")

    # 3) UNION test: Gray-adjacent AND 3D-close (load the position-level 3D result)
    v3d = REPO / "results" / "contacts" / "validate_3d_complete.json"
    union = None
    if v3d.exists():
        d3 = json.load(open(v3d))
        pair3d = {}
        for src, info in d3["sources"].items():
            for pp in info.get("per_pair", []):
                key = tuple(pp["pair"])
                if key not in pair3d or pp["models"] > pair3d[key]["models"]:
                    pair3d[key] = pp
        rows = []
        for key, pairs in per_pair.items():
            g = sum(1 for (ai, aj) in pairs
                    if ai in GRAY and aj in GRAY and hamming(GRAY[ai], GRAY[aj]) == 1)
            n = sum(1 for (ai, aj) in pairs if ai in GRAY and aj in GRAY)
            c3 = pair3d.get(tuple(key), {})
            rows.append({"pair": list(key), "gray_adjacent": g, "n": n,
                         "frac_h1": round(g / n, 3) if n else None,
                         "contact_frac": c3.get("fraction"),
                         "contact_enrichment": c3.get("enrichment")})
        union = rows
        print("\n=== UNION: Gray-adjacency x 3D contact ===")
        for r in rows:
            print(f"  {str(r['pair']):>10}: GrayH1={r['frac_h1']} "
                  f"contact={r['contact_frac']} enrich={r['contact_enrichment']}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump({"n_observed_pairs": obs_n, "h1_observed": obs_h[1],
               "h1_expected_frac": exp_h1, "h1_observed_frac": obs_h1,
               "h1_enrichment": enrich_h1, "h1_p_value": p_h1,
               "h_distribution_observed": dict(obs_h),
               "h_distribution_null": dict(null_h),
               "per_pair": per_pair_rows, "union": union},
              open(OUT, "w"), indent=1)
    print(f"\nsaved -> {OUT}")


if __name__ == "__main__":
    main()
