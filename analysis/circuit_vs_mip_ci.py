#!/usr/bin/env python3
"""Circuit-L2 vs APC-corrected MI (MIp): paired per-protein differences with CIs.

Reads the per-protein results of scripts/run_contact_campaign.py (holdout
split and the five CV folds, results/kmap_structure/) and reports, for
prec@L/5, prec@L and AUC, the mean paired difference circuit - MIp with a
95% protein-bootstrap interval (10,000 resamples), the Wilcoxon signed-rank
p-value and the win/loss count.

Output: analysis/circuit_vs_mip_ci.json
Run:    python analysis/circuit_vs_mip_ci.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

REPO = Path(__file__).resolve().parents[1]
RES = REPO / "results" / "kmap_structure"
METRICS = ("prec@L/5", "prec@L", "auc")


def paired(rows, metric):
    x = np.array([r["L2"]["ranked"][metric] - r["mip"]["ranked"][metric] for r in rows])
    rng = np.random.default_rng(0)
    boot = x[rng.integers(0, len(x), size=(10000, len(x)))].mean(axis=1)
    return {"n": len(x), "circuit": float(np.mean([r["L2"]["ranked"][metric] for r in rows])),
            "mip": float(np.mean([r["mip"]["ranked"][metric] for r in rows])),
            "mean_diff": float(x.mean()),
            "ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
            "wilcoxon_p": float(wilcoxon(x).pvalue) if np.any(x != 0) else 1.0,
            "circuit_better": int((x > 0).sum()), "mip_better": int((x < 0).sum())}


def main() -> int:
    hold = json.load(open(RES / "holdout_split.json"))["per_protein"]
    cv = [r for f in range(5) for r in json.load(open(RES / f"cv_fold{f}.json"))["per_protein"]]
    out = {"holdout": {m: paired(hold, m) for m in METRICS},
           "cv_pooled": {m: paired(cv, m) for m in METRICS}}
    (REPO / "analysis" / "circuit_vs_mip_ci.json").write_text(json.dumps(out, indent=1))
    for split, res in out.items():
        for m, v in res.items():
            print(f"{split:9s} {m:9s} circuit {v['circuit']:.3f} MIp {v['mip']:.3f} "
                  f"diff {v['mean_diff']:+.3f} [{v['ci95'][0]:+.3f}, {v['ci95'][1]:+.3f}] "
                  f"p={v['wilcoxon_p']:.2g} ({v['circuit_better']}/{v['mip_better']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
