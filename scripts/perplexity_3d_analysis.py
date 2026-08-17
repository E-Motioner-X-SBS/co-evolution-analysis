#!/usr/bin/env python3
"""
perplexity_3d_analysis.py — does the PERPLEXITY metric carry 3D information?

Beyond "are the top-combined pairs structural", the deeper question: is the
perplexity RATIO (determinism: PP(j)/mean PP(j|i)) itself predictive of 3D
contact?

Method:
  1. all variable-position pairs (window 30) from the corrected alignment
  2. per pair: perplexity ratio + MI + combined score (shared module,
     vectorized — the exact functions every script uses)
  3. per pair: intra-chain contact fraction over the 299 models
     (rules3d_common, matched controls)
  4. tests:
     a. Spearman(ratio, contact_fraction) over pairs with >= 100 models
     b. ratio quartiles -> mean contact fraction + enrichment vs random
     c. Spearman(MI, contact_fraction) for comparison (is ratio better?)
     d. combined score vs contact (same)
  5. the "ratio rescues deterministic pairs" claim: are the pairs where
     ratio >> MI-rank the structural ones?

Output: results/contacts/perplexity_3d_analysis.json
"""
import json
import math
import sys
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/kmap-sbm-validation/src")
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")
from coevolution_shared import (load_position_arrays, compute_entropy_vectorized,
                                combined_pair_scores)  # noqa: E402
from rules3d_common import OmicronModels, evaluate_pairs  # noqa: E402

FASTA = REPO / "Spike_protein.aln-fasta"
OUT = REPO / "results" / "contacts" / "perplexity_3d_analysis.json"


def main():
    # 1. pairs + metrics from the corrected alignment
    pos_arrays, n_all, full_len = load_position_arrays(
        fasta_path=str(FASTA), max_pos=None, aligned=True, clear_cache=True)
    ent = compute_entropy_vectorized(pos_arrays, n_all, full_len)
    var = [p for p in range(full_len) if ent[p] > 0.3]
    pairs = [(i, j) for a, i in enumerate(var) for j in var[a + 1:]
             if j - i <= 30]
    scored = combined_pair_scores(pos_arrays, pairs, n_all, ent)
    metrics = {(s["pos_i"], s["pos_j"]): s for s in scored}
    print(f"variable positions: {len(var)} | pairs scored: {len(metrics)}")

    # 2. 3D contact fractions
    models = OmicronModels()
    rows, control, model_count = evaluate_pairs(models.iter_models(),
                                                 sorted(metrics))
    print(f"models: {model_count} | random contact rate: {control:.4f}")

    # 3. join + tests
    joined = []
    for p in sorted(metrics):
        m = metrics[p]
        r = rows.get(p)
        if not r or r["models"] < 50:
            continue
        frac = r["intra_contacts"] / r["models"]
        joined.append({"pair": list(p), "mi": m["mi"], "ratio": m["ratio"],
                       "combined": m["combined"],
                       "n_models": r["models"],
                       "contact_frac": round(frac, 4),
                       "enrichment": round(frac / control, 2) if control else None})
    print(f"pairs with >=50 models: {len(joined)}")

    import numpy as np
    from scipy.stats import spearmanr
    xs = [j["ratio"] for j in joined]
    xm = [j["mi"] for j in joined]
    xc = [j["combined"] for j in joined]
    y = [j["contact_frac"] for j in joined]
    r_ratio = spearmanr(xs, y)
    r_mi = spearmanr(xm, y)
    r_comb = spearmanr(xc, y)
    print(f"\nSpearman(ratio, contact):  rho={r_ratio.statistic:.3f} p={r_ratio.pvalue:.3g}")
    print(f"Spearman(MI, contact):      rho={r_mi.statistic:.3f} p={r_mi.pvalue:.3g}")
    print(f"Spearman(combined, contact):rho={r_comb.statistic:.3f} p={r_comb.pvalue:.3g}")

    # ratio quartiles -> mean contact fraction
    rs = sorted(j["ratio"] for j in joined)
    q = [rs[int(len(rs) * k / 4)] for k in range(1, 4)]
    print("\nratio quartile -> mean contact fraction:")
    quartiles = []
    for k in range(4):
        lo = 0 if k == 0 else q[k - 1]
        hi = q[k] if k < 3 else float("inf")
        grp = [j for j in joined if lo <= j["ratio"] < hi]
        if grp:
            mcf = sum(j["contact_frac"] for j in grp) / len(grp)
            enr = mcf / control
            quartiles.append({"range": [round(lo, 2), round(hi, 2)],
                              "n": len(grp), "mean_contact": round(mcf, 3),
                              "enrichment": round(enr, 2)})
            print(f"  Q{k+1} ratio[{lo:.2f},{hi:.2f}): n={len(grp)} "
                  f"mean_contact={mcf:.3f} enrich={enr:.2f}x")

    # 4. the "ratio rescues deterministic pairs" claim: top pairs by ratio
    #    that MI ranks low — are they the structural ones?
    print("\ntop 10 by ratio (with MI rank and contact):")
    mi_rank = {tuple(j["pair"]): rank for rank, j in
               enumerate(sorted(joined, key=lambda x: -x["mi"]))}
    top_ratio = sorted(joined, key=lambda x: -x["ratio"])[:10]
    rescued = []
    for j in top_ratio:
        rank = mi_rank.get(tuple(j["pair"]))
        print(f"  {j['pair']}: ratio={j['ratio']:.2f} (MI rank {rank + 1}/"
              f"{len(joined)}) contact={j['contact_frac']:.3f} "
              f"enrich={j['enrichment']}x")
        if rank is not None and rank > len(joined) // 2 and \
                (j["enrichment"] or 0) > 2:
            rescued.append(j["pair"])
    print(f"\npairs with LOW MI rank but HIGH ratio AND structural: {rescued}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump({"n_pairs": len(joined), "control_rate": control,
               "spearman_ratio": {"rho": r_ratio.statistic, "p": r_ratio.pvalue},
               "spearman_mi": {"rho": r_mi.statistic, "p": r_mi.pvalue},
               "spearman_combined": {"rho": r_comb.statistic, "p": r_comb.pvalue},
               "ratio_quartiles": quartiles,
               "pairs": joined},
              open(OUT, "w"), indent=1)
    print(f"\nsaved -> {OUT}")


if __name__ == "__main__":
    main()
