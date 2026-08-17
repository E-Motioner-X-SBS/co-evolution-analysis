#!/usr/bin/env python3
"""
validate_flipped_3d.py — FLIPPED (forbidden) rules vs 3D structure.

The flipped K-map marks residue pairs (pos_i=aa_i, pos_j=aa_j) that are NEVER
observed together in the alignment — negative selection. Two hypotheses:

  H_steric : forbidden pairs sit at 3D-close positions, so the combination is
             structurally incompatible (would clash) -> forbidden positions
             are ENRICHED for contacts.
  H_rarity : forbiddenness reflects sequence-space rarity / lineage, with no
             structural signature -> forbidden positions are NOT enriched.

Test (position level, per unique forbidden position pair):
  - contact fraction of forbidden position pairs across the 299 models
  - vs matched random control -> enrichment + binomial p
Also reports the residue-level detail: for the residues that DO occur at those
positions in the models (the observed combos), their distance distribution.

Output: results/contacts/validate_flipped_3d.json
"""
import json
import sys
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
sys.path.insert(0, str(REPO / "scripts"))
from rules3d_common import (OmicronModels, evaluate_pairs, row_summary,
                            aggregate_set, CONTACT)  # noqa: E402

INV = REPO / "results" / "contacts" / "rules_inventory.json"
OUT = REPO / "results" / "contacts" / "validate_flipped_3d.json"


def main():
    inv = json.load(open(INV))
    flipped = inv["flipped_rules"]
    pos_pairs = sorted({(r["pos_i"], r["pos_j"]) for r in flipped})
    print(f"unique forbidden position pairs: {len(pos_pairs)} "
          f"(from {len(flipped)} residue rules)")
    models = OmicronModels()
    rows, control, model_count = evaluate_pairs(models.iter_models(), pos_pairs)
    print(f"models used: {model_count} | random control: {control:.4f}")

    print("\n=== forbidden position pairs (intra-chain) ===")
    print(f"{'pair':>10} {'n':>4} {'ct':>4} {'frac':>7} {'enr':>7} {'p':>9} {'medA':>6}")
    table = []
    for p in pos_pairs:
        s = row_summary({p: rows[p]}, control)
        if not s:
            continue
        s = s[0]
        table.append(s)
        print(f"{str(p):>10} {s['models']:>4} {s['contacts']:>4} "
              f"{s['fraction']:>7.3f} {s['enrichment']:>7.2f} "
              f"{s['p_value']:>9.2g} {s['median_dist'] or 0:>6.1f}")

    agg = aggregate_set(rows, control, pos_pairs)
    print(f"\nPOOLED forbidden: {agg['contacts']}/{agg['n']} "
          f"({agg['fraction']}, enrich {agg['enrichment']}x, p={agg['p_value']:.3g})")
    n_sig = sum(1 for s in table if s["p_value"] < 0.05)
    n_enr = sum(1 for s in table if s["enrichment"] > 1.5)
    print(f"pairs with p<0.05: {n_sig}/{len(table)} | enrichment>1.5x: {n_enr}/{len(table)}")
    verdict = ("STERIC (forbidden positions are structurally close)"
               if agg["enrichment"] and agg["enrichment"] > 1.5
               else "RARITY (no structural signature)")
    print(f"VERDICT: {verdict}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump({"models_used": model_count, "control_rate": control,
               "n_forbidden_rules": len(flipped),
               "n_position_pairs": len(pos_pairs),
               "pooled": agg, "per_pair": table, "verdict": verdict},
              open(OUT, "w"), indent=1)
    print(f"saved -> {OUT}")


if __name__ == "__main__":
    main()
