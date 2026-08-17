#!/usr/bin/env python3
"""
trimer_interface.py — INTER-CHAIN (across-protomer) contact analysis.

The Spike models are homo-trimers (chains A, B, C). Co-evolving pairs may be
in contact either WITHIN one protomer (intra-chain) or ACROSS protomers at
the trimer interface (inter-chain) — e.g., the RBD-up conformation brings the
RBD of one protomer near the NTD of the next.

For every co-evolving position pair (all sources):
  - intra-chain contact fraction (already in validate_3d_complete)
  - INTER-chain contact fraction: residue i in chain A vs residue j in chain B
    (min over chain pairs), and the reverse assignment
  - enrichment vs matched random control (inter-chain)
  - report pairs whose STRUCTURAL signal is inter-chain only (missed by the
    monomer view)

Output: results/contacts/trimer_interface.json
"""
import json
import math
import sys
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
sys.path.insert(0, str(REPO / "scripts"))
from rules3d_common import (OmicronModels, CONTACT, MIN_SEP, N_RAND,
                            dist)  # noqa: E402

INV = REPO / "results" / "contacts" / "rules_inventory.json"
OUT = REPO / "results" / "contacts" / "trimer_interface.json"


def main():
    inv = json.load(open(INV))
    all_pairs = sorted({(p["pos_i"], p["pos_j"])
                        for ps in inv["position_pairs"].values() for p in ps})
    models = OmicronModels()
    rows = {p: {"models": 0, "inter_contacts": 0, "dists": []} for p in all_pairs}
    rnd_c = 0
    rnd_n = 0
    n_models = 0
    import random
    for uid, aligned, cmap, coords in models.iter_models():
        n_models += 1
        chains = sorted(coords)
        if len(chains) < 2:
            continue
        mapped = {}
        for (i, j) in all_pairs:
            ri, rj = cmap.get(i), cmap.get(j)
            if ri is None or rj is None or abs(ri - rj) < MIN_SEP:
                continue
            mapped[(i, j)] = (ri, rj)
        # inter-chain: min distance between residue ri in chain X and rj in
        # chain Y (X != Y), over all ordered chain pairs
        for (i, j), (ri, rj) in mapped.items():
            best = None
            for x in chains:
                for y in chains:
                    if x == y:
                        continue
                    if ri in coords[x] and rj in coords[y]:
                        dd = dist(coords[x][ri], coords[y][rj])
                        best = dd if best is None else min(best, dd)
            if best is not None:
                rows[(i, j)]["models"] += 1
                rows[(i, j)]["dists"].append(best)
                if best < CONTACT:
                    rows[(i, j)]["inter_contacts"] += 1
        # control
        chain = chains[0]
        residues = sorted(coords[chain])
        seps = [abs(ri - rj) for (ri, rj) in mapped.values()]
        if not seps or len(residues) < 10:
            continue
        rng = random.Random(11 + n_models)
        for _ in range(N_RAND):
            a = rng.choice(residues)
            s = rng.choice(seps)
            b = a + s
            if b not in coords[chain]:
                continue
            # inter-chain control: same residue pair across two chains
            if len(chains) >= 2:
                y = chains[1] if chains[1] != chain else chains[0]
                if a in coords[chain] and b in coords[y]:
                    if dist(coords[chain][a], coords[y][b]) < CONTACT:
                        rnd_c += 1
            rnd_n += 1
    control = rnd_c / rnd_n if rnd_n else 0.0
    print(f"models with >=2 chains: {n_models} | inter-chain random control: {control:.4f}")

    from scipy.stats import binomtest
    print(f"\n{'pair':>10} {'n':>4} {'ct':>4} {'frac':>7} {'enr':>7} {'p':>9} {'medA':>6}")
    out = []
    for p in all_pairs:
        r = rows[p]
        if r["models"] == 0:
            continue
        frac = r["inter_contacts"] / r["models"]
        enr = frac / control if control else float("nan")
        pv = binomtest(r["inter_contacts"], r["models"], control,
                       alternative="greater").pvalue if control else 1.0
        md = sorted(r["dists"])[len(r["dists"]) // 2]
        out.append({"pair": list(p), "models": r["models"],
                    "inter_contacts": r["inter_contacts"],
                    "fraction": round(frac, 3), "enrichment": round(enr, 2),
                    "p_value": pv, "median_inter_dist": round(md, 2)})
        print(f"{str(p):>10} {r['models']:>4} {r['inter_contacts']:>4} "
              f"{frac:>7.3f} {enr:>7.2f} {pv:>9.2g} {md:>6.1f}")

    # pairs with strong INTER-chain signal (missed by monomer view)
    strong = [o for o in out if o["enrichment"] and o["enrichment"] > 2
              and o["p_value"] < 0.05]
    print(f"\ninter-chain-significant pairs (enrich>2x, p<0.05): {len(strong)}")
    for o in strong:
        print("  ", o["pair"], o["fraction"], o["enrichment"], "x")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump({"models_with_trimers": n_models, "control_rate": control,
               "per_pair": out, "inter_chain_significant": strong},
              open(OUT, "w"), indent=1)
    print(f"saved -> {OUT}")


if __name__ == "__main__":
    main()
