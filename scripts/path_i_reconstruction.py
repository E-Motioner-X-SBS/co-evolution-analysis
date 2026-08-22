#!/usr/bin/env python3
"""
path_i_reconstruction.py — PATH I: how much 3D structure can the circuits
actually reconstruct?  Quantified against LITERATURE thresholds only.

Metrics
-------
For each method m ∈ {circuit-L2, plain-MI, perplexity, PSICOV-published} and
each protein:
    hits(m, K)          = #native contacts among top-K ranked candidates
    recall(m, K)        = hits / n_native_contacts
    coverage(m, K)      = hits / L            ("one correct contact per residue"
                                                  heuristic of Jones 2012)
Aggregates over 150 proteins with protein-level 100/49 split for circuit
scores (train-derived cell rates), literature baselines label-free.

Literature anchors (all quoted, none self-computed):
  [J12] Jones et al. 2012, abstract: "118/150 targets ≥ 0.5 precision at
        L/5 long-range ... sufficient to be of significant benefit in
        protein structure prediction."
  [W09] Weigt et al. 2009 PNAS: top-L coevolutionary contacts at ~50%
        precision were sufficient to fold TPR domains.
  [C18] Adhikari/Cheng CONFOLD line of work: model quality (TM-score)
        increases steeply with contact precision; ~30% precision at L/5–L/2
        already yields useful folds; 50%+ approaches native-like cores.
Verdict rule (documented): a method is *reconstruction-ready* iff median
long-range prec@L/5 ≥ 0.30 AND median recall@L ≥ 0.10 (conservative reading
of [C18]+[J12]); otherwise we report the gap explicitly.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import kmap_structure as ks  # noqa: E402

OUT = ks.RESULTS


def recall_at(score, lab, K, n_contacts):
    if len(score) == 0 or n_contacts == 0:
        return 0.0
    top = np.argsort(-score, kind="stable")[:min(K, len(score))]
    return float(lab[top].sum()) / n_contacts


def main():
    from run_contact_campaign import load_all_caches, train_counts, sep_only_counts

    caches = load_all_caches()
    targets = sorted(caches)
    fl3 = ks.protein_folds(targets, n_folds=3, seed=42)
    train, test = fl3[0] + fl3[1], fl3[2]

    _, c2 = train_counts(train, caches)
    cs = sep_only_counts(train, caches)

    KS = (("L/10", 1 / 10), ("L/5", 1 / 5), ("L/2", 1 / 2), ("L", 1))
    methods = ("circ_L2", "MI", "PP", "sep_only", "PSICOV_pub")
    rec = {m: {k: [] for k, _ in KS} for m in methods}
    cov = {m: {k: [] for k, _ in KS} for m in methods}
    prec_lr = {m: [] for m in methods}

    for name in test:
        tc = caches[name]
        sc = {
            "circ_L2": np.zeros(len(tc.pairs)),
            "MI": np.nan_to_num(tc.mi, nan=-1e9),
            "PP": np.nan_to_num(tc.pp, nan=-1e9),
            "sep_only": ks.cell_rate_score(cs, tc.sbin.astype(np.int64)),
            "PSICOV_pub": None,
        }
        for b in range(4):
            m = tc.sbin == b
            if m.any():
                sc["circ_L2"][m] = ks.cell_rate_score(c2[b], tc.level2_cells(b))
        con = ks.parse_con(ks.RAW / "con" / f"{name}.out")
        lut = {(i - 1, j - 1): s for i, j, s in con}
        sc["PSICOV_pub"] = np.array(
            [lut.get((int(a), int(b)), -np.inf) for a, b in tc.pairs])

        n_cont = int(tc.lab.sum())
        L = tc.n_cols
        sep = tc.pairs[:, 1] - tc.pairs[:, 0]
        lr = sep >= 24                       # long-range band (CASP/Jones >23)
        for meth, s in sc.items():
            for kn, frac in KS:
                K = max(1, round(frac * L))
                rec[meth][kn].append(recall_at(s, tc.lab, K, n_cont))
                hits = s[np.argsort(-s, kind="stable")[:min(K, len(s))]]
                cov[meth][kn].append(float(len(hits)) / L)
            if lr.sum() >= 20:
                k5 = max(1, round(L / 5))
                prec_lr[meth].append(
                    ks.precision_at_k(sc[meth][lr], tc.lab[lr], min(k5, int(lr.sum()))))

    agg: dict = {}
    verdict: dict = {}
    for m in methods:
        agg[m] = {"recall_" + k: float(np.mean(v)) for k, v in rec[m].items()}
        agg[m].update({"coverage_" + k: float(np.mean(v))
                       for k, v in cov[m].items()})
        agg[m]["median_longrange_prec_L5"] = float(np.median(prec_lr[m]))
        ready = (agg[m]["median_longrange_prec_L5"] >= 0.30
                 and agg[m]["recall_L"] >= 0.10)
        verdict[m] = {"reconstruction_ready": bool(ready)}
    agg["verdict"] = verdict
    agg["literature_anchors"] = {
        "Jones2012": "118/150 targets with L/5 long-range precision >= 0.5 "
                     "'sufficient to be of significant benefit in protein "
                     "structure prediction'",
        "Weigt2009": "top-L contacts near 50% precision folded TPR domains",
        "CONFOLD-line": "~30% precision at L/5-L/2 yields useful folds; "
                        "TM-score rises steeply with precision",
        "verdict_rule": "ready iff median LR prec@L/5 >= 0.30 AND recall@L >= 0.10",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "path_i_reconstruction.json").write_text(json.dumps(agg, indent=1))
    print(json.dumps({k: v for k, v in agg.items() if k != "verdict"}, indent=1))
    print(json.dumps(verdict, indent=1))


if __name__ == "__main__":
    main()
