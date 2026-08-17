#!/usr/bin/env python3
"""
validate_3d_complete.py — POSITION-LEVEL rule->3D validation.

Every co-evolving position pair from EVERY source script (8 sources, from the
rules_inventory) is checked against the 299 Omicron 3D models:
  - intra-chain contact fraction (Cb-Cb < 8 A, sep >= 3)
  - inter-chain contact fraction (trimer interface)
  - enrichment vs matched random control
  - binomial p-value per pair and per rule set

Output: results/contacts/validate_3d_complete.json + printed table.
"""
import json
import sys
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
sys.path.insert(0, str(REPO / "scripts"))
from rules3d_common import (OmicronModels, evaluate_pairs, row_summary,
                            aggregate_set)  # noqa: E402

INV = REPO / "results" / "contacts" / "rules_inventory.json"
OUT = REPO / "results" / "contacts" / "validate_3d_complete.json"


def main():
    inv = json.load(open(INV))
    sources = inv["position_pairs"]
    models = OmicronModels()

    results = {}
    # evaluate each source's pairs on a SHARED random control (pooled over all
    # rule pairs, so enrichments are comparable across sources)
    all_pairs = sorted({(p["pos_i"], p["pos_j"])
                        for ps in sources.values() for p in ps})
    print(f"unique position pairs across all sources: {len(all_pairs)}")
    rows, control, model_count = evaluate_pairs(models.iter_models(), all_pairs)
    print(f"models used: {model_count} | pooled random contact rate: {control:.4f}")

    for name, ps in sources.items():
        pair_list = [(p["pos_i"], p["pos_j"]) for p in ps]
        rrows = row_summary(rows, control, key_intra=True)
        agg = aggregate_set(rows, control, pair_list)
        # inter-chain aggregate
        n_i = sum(rows[p]["models"] for p in pair_list)
        c_i = sum(rows[p]["inter_contacts"] for p in pair_list)
        frac_i = c_i / n_i if n_i else 0.0
        results[name] = {
            "n_pairs": len(pair_list),
            "intra": agg,
            "inter": {"n": n_i, "contacts": c_i,
                      "fraction": round(frac_i, 3),
                      "enrichment": round(frac_i / control, 2) if control else None},
            "per_pair": [r for r in rrows if r["pair"] in set(pair_list)],
        }
        print(f"\n=== {name} ({len(pair_list)} pairs) ===")
        print(f"  INTRA: {agg['contacts']}/{agg['n']} contacts "
              f"({agg['fraction']}, enrich {agg['enrichment']}x, p={agg['p_value']:.3g})")
        print(f"  INTER: {c_i}/{n_i} contacts ({frac_i:.3f}, "
              f"enrich {frac_i / control:.2f}x)" if control else "  INTER: n/a")

    print("\n=== per-pair detail (intra-chain) ===")
    print(f"{'pair':>10} {'n':>4} {'ct':>4} {'frac':>6} {'enr':>6} {'p':>8} {'medA':>6}")
    for r in sorted(rows.keys()):
        s = row_summary({r: rows[r]}, control)
        if not s:
            print(f"{str(r):>10} {0:>4} {0:>4} {'n/a':>6} {'n/a':>6} {'-':>8} {'-':>6}")
            continue
        s = s[0]
        print(f"{str(r):>10} {s['models']:>4} {s['contacts']:>4} "
              f"{s['fraction']:>6.3f} {s['enrichment']:>6.2f} "
              f"{s['p_value']:>8.2g} {s['median_dist'] or 0:>6.1f}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump({"models_used": model_count, "control_rate": control,
               "sources": results,
               "skipped_models": models.skipped},
              open(OUT, "w"), indent=1)
    print(f"\nsaved -> {OUT}")


if __name__ == "__main__":
    main()
