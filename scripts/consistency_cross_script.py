#!/usr/bin/env python3
"""
consistency_cross_script.py — do the 23 scripts AGREE on which positions
co-evolve?

For each source script's position-pair list (from the rules inventory):
  - pair set, overlap with every other source (Jaccard + shared count)
  - agreement matrix + a consensus pair list (pairs reported by >= k sources)
  - the consensus pairs are then cross-checked against the 3D validation
    (are the most-agreed pairs the structurally-supported ones?)

Output: results/contacts/cross_script_consistency.json
"""
import json
import sys
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
INV = REPO / "results" / "contacts" / "rules_inventory.json"
V3D = REPO / "results" / "contacts" / "validate_3d_complete.json"
OUT = REPO / "results" / "contacts" / "cross_script_consistency.json"


def main():
    inv = json.load(open(INV))
    sources = inv["position_pairs"]
    sets = {name: {(p["pos_i"], p["pos_j"]) for p in ps}
            for name, ps in sources.items() if ps}
    names = sorted(sets)

    print("=== pair counts per source ===")
    for n in names:
        print(f"  {n}: {len(sets[n])}")

    print("\n=== overlap matrix (shared pairs / Jaccard) ===")
    print(f"{'':>16}", end="")
    for n in names:
        print(f"{n[:9]:>10}", end="")
    print()
    matrix = {}
    for a in names:
        row = []
        print(f"{a[:16]:>16}", end="")
        for b in names:
            inter = len(sets[a] & sets[b])
            union = len(sets[a] | sets[b])
            jac = inter / union if union else 0
            row.append({"shared": inter, "jaccard": round(jac, 2)})
            print(f"{inter:>10}", end="")
        print()
        matrix[a] = row

    # consensus: pairs reported by >= 2 sources
    from collections import Counter
    vote = Counter()
    for s in sets.values():
        for p in s:
            vote[p] += 1
    consensus = {p: v for p, v in vote.items() if v >= 2}
    print(f"\nconsensus pairs (>=2 sources): {len(consensus)}")
    for p, v in sorted(consensus.items(), key=lambda kv: -kv[1]):
        print(f"  {p}: {v} sources")

    # cross-check consensus vs 3D
    d3 = json.load(open(V3D)) if V3D.exists() else None
    rows = []
    if d3:
        pair3d = {}
        for src, info in d3["sources"].items():
            for pp in info.get("per_pair", []):
                key = tuple(pp["pair"])
                if key not in pair3d or pp["models"] > pair3d[key]["models"]:
                    pair3d[key] = pp
        print("\n=== consensus pairs vs 3D contact ===")
        for p, v in sorted(consensus.items(), key=lambda kv: -kv[1]):
            c = pair3d.get(p, {})
            rows.append({"pair": list(p), "sources": v,
                         "contact_frac": c.get("fraction"),
                         "enrichment": c.get("enrichment")})
            print(f"  {p}: {v} sources | contact {c.get('fraction')} "
                  f"enrich {c.get('enrichment')}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump({"sets": {n: sorted(sets[n]) for n in names},
               "overlap": matrix, "consensus": {str(k): v
                                                for k, v in consensus.items()},
               "consensus_vs_3d": rows},
              open(OUT, "w"), indent=1)
    print(f"\nsaved -> {OUT}")


if __name__ == "__main__":
    main()
