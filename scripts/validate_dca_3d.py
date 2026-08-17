#!/usr/bin/env python3
"""
validate_dca_3d.py — DCA (mfDCA direct information) vs MI as a 3D-contact
predictor on the 299 Omicron Spike models.

Question: which method's top pairs are more 3D-close — the mfDCA direct-
coupling ranking (dca_mf_analysis, Morcos 2011) or the plain-MI ranking
(gpu_full / run_allseq)?

Method: for each method's top-20 pairs, evaluate intra-chain contact fraction
+ enrichment vs the shared matched-separation random control (rules3d_common),
and report per-pair stats. This is the Omicron-model analogue of the PSICOV
benchmark (validate_contacts.py), restricted to the Spike's own 3D models.

Output: results/contacts/validate_dca_3d.json + printed comparison table.
"""
import json
import sys
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
sys.path.insert(0, str(REPO / "scripts"))
from rules3d_common import OmicronModels, evaluate_pairs, row_summary  # noqa: E402

BASE = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
OUT = REPO / "results" / "contacts" / "validate_dca_3d.json"


def _pairs_from(path, key_i, key_j, score_key, limit=20):
    try:
        d = json.load(open(BASE / path))
    except Exception:
        return []
    items = d
    for k in (["results"], ["top_co_evolving_pairs"], ["top_20_di"],
              ["top_mi_pairs"], ["couplings"]):
        if isinstance(d, dict) and k[0] in d:
            items = d[k[0]]
            break
    if not isinstance(items, list):
        return []
    out = []
    for r in items[:limit]:
        if isinstance(r, dict) and r.get(key_i) is not None:
            out.append((int(r[key_i]), int(r[key_j])))
    return out


def main():
    # DCA direct-information top-20 (mfDCA, Morcos 2011)
    dca_pairs = _pairs_from("dca_results/dca_mf_summary.json",
                            "pos_i", "pos_j", "score", 20)
    # Plain-MI top-20 (gpu_full couplings / run_allseq top_mi_pairs)
    mi_pairs = _pairs_from("full_position_results/full_analysis_summary.json",
                           "pos_i", "pos_j", "mi", 20)
    if not mi_pairs:
        mi_pairs = _pairs_from("full_gpu_results/gpu_summary.json",
                               "pos_i", "pos_j", "mi", 20)

    all_pairs = sorted(set(dca_pairs) | set(mi_pairs))
    print(f"DCA top-20: {len(dca_pairs)} pairs | MI top-20: {len(mi_pairs)} pairs | "
          f"union: {len(all_pairs)}")

    models = OmicronModels()
    rows, control, model_count = evaluate_pairs(models.iter_models(), all_pairs)
    print(f"models used: {model_count} | pooled random contact rate: {control:.4f}")

    def summarize(name, pair_list):
        rrows = row_summary(rows, control)
        sel = set(pair_list)
        per = [r for r in rrows if r["pair"] in sel]
        n = sum(rows[p]["models"] for p in pair_list if p in rows)
        c = sum(rows[p]["intra_contacts"] for p in pair_list if p in rows)
        frac = c / n if n else 0.0
        enrich = frac / control if control else float("nan")
        from scipy.stats import binomtest
        pval = binomtest(c, n, control, alternative="greater").pvalue if (control and n) else 1.0
        # precision@L: fraction of the method's top-L pairs that are contacts
        contact_pairs = [p for p in pair_list if p in rows and rows[p]["models"] > 0
                         and rows[p]["intra_contacts"] / rows[p]["models"] > 0.5]
        prec = len(contact_pairs) / len(pair_list) if pair_list else 0.0
        return {"n_pairs": len(pair_list), "n_models": n, "contacts": c,
                "fraction": round(frac, 3), "enrichment": round(enrich, 2),
                "p_value": float(pval), "precision@L": round(prec, 3),
                "contact_pairs": [list(p) for p in contact_pairs],
                "per_pair": [{"pair": list(r["pair"]), "models": r["models"],
                              "contacts": r["contacts"], "fraction": r["fraction"],
                              "enrichment": r["enrichment"], "p_value": float(r["p_value"]),
                              "median_dist": r["median_dist"]} for r in per]}

    res = {
        "models_used": model_count, "control_rate": control,
        "dca_mf": summarize("dca_mf", dca_pairs),
        "mi": summarize("mi", mi_pairs),
    }

    print("\n=== DCA (mfDCA direct information) top-20 vs MI top-20 ===")
    for name in ("dca_mf", "mi"):
        s = res[name]
        print(f"  {name:>8}: {s['contacts']}/{s['n_models']} contacts "
              f"(frac {s['fraction']}, enrich {s['enrichment']}x, "
              f"p={s['p_value']:.2g}, precision@L={s['precision@L']})")
    print("\n  DCA contact pairs:", res["dca_mf"]["contact_pairs"])
    print("  MI  contact pairs:", res["mi"]["contact_pairs"])

    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(OUT, "w"), indent=1)
    print(f"\nsaved -> {OUT}")


if __name__ == "__main__":
    main()
