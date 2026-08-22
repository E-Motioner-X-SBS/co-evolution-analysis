#!/usr/bin/env python3
"""
path_l_rank_fusion.py — Rank-fusion hybrid: DCA continuous score × K-map prior.

Motivation (plans/12b): circH AUC 0.658 > rawDCA AUC 0.583, but circH
prec@L/5 0.124 < rawDCA 0.169. The K-map cell-rate smoothing improves GLOBAL
ranking while losing top-K sharpness. Hypothesis: fusing the two signals
captures both properties — DCA provides continuous physics, the K-map prior
denoises phylogenetic noise.

Fusion rules tested (all label-free on test proteins):
  F1  rank-sum:      fused = rank(DCA) + rank(cell-rate)
  F2  z-score sum:   fused = z(DCA) + z(cell-rate)
  F3  product:       fused = z(DCA) * cell-rate   (gated by prior)
  F4  geometric:     fused = sqrt(rank(DCA) * rank(cell-rate))
  F5  DCA-gated:     fused = DCA where cell-rate > median else DCA - |z(prior)|

Evaluation: 5-fold protein-level CV identical to PATH B/H.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import kmap_structure as ks  # noqa: E402

OUT = ks.RESULTS
DI_DIR = OUT / "di"


class CountsH:
    N_CELLS = 1024

    def __init__(self):
        self.pos = np.zeros(self.N_CELLS, np.int64)
        self.tot = np.zeros(self.N_CELLS, np.int64)

    def add(self, idx, lab):
        np.add.at(self.tot, idx, 1)
        np.add.at(self.pos, idx, lab.astype(np.int64))


def di_bin(di_scores):
    r = np.argsort(np.argsort(di_scores, kind="stable"))
    return (4 * r // max(len(r), 1)).clip(0, 3).astype(np.int8)


def zscore(s):
    mu = np.mean(s)
    sd = np.std(s)
    return (s - mu) / max(sd, 1e-12)


def rank01(s):
    """Normalized rank in [0,1]."""
    r = np.argsort(np.argsort(s, kind="stable"))
    return r / max(len(s) - 1, 1)


def main():
    from run_contact_campaign import load_all_caches
    caches = load_all_caches()
    targets = sorted(caches)
    folds = ks.protein_folds(targets, n_folds=5, seed=42)

    methods = ("rawDCA", "circH", "ranksum", "zsum", "product", "geometric",
               "dcagated")
    agg = {m: {"prec@L/5": [], "auc": []} for m in methods}
    per_protein = []

    for fold in range(5):
        test = folds[fold]
        train = [t for g, grp in enumerate(folds) if g != fold for t in grp]
        tr = CountsH()
        for name in train:
            tct = caches[name]
            dca_t = np.load(DI_DIR/f"{name}.npz")["dca_pairs"]
            dib_t = di_bin(dca_t)
            idx_t = (((tct.g_i.astype(np.int64)*8 + tct.g_j)*4 + tct.sbin)*4 + dib_t)
            tr.add(idx_t, tct.lab)
        truth = ks.label_truth_table(
            ks.CellCounts(CountsH.N_CELLS, tr.pos.copy(), tr.tot.copy()), 2.0, 30)
        from run_contact_campaign import qm_extract
        qm_extract(truth)

        for name in test:
            tc = caches[name]
            dca = np.load(DI_DIR/f"{name}.npz")["dca_pairs"]
            dib = di_bin(dca)
            idx = (((tc.g_i.astype(np.int64)*8 + tc.g_j)*4 + tc.sbin)*4 + dib)
            prior = ks.cell_rate_score(ks.CellCounts(CountsH.N_CELLS, tr.pos, tr.tot), idx)
            L = tc.n_cols; K = max(1, round(L / 5))
            lab = tc.lab

            scores = {
                "rawDCA": dca,
                "circH": prior,
                "ranksum": rank01(dca) + rank01(prior),
                "zsum": zscore(dca) + zscore(prior),
                "product": zscore(dca) * prior,
                "geometric": np.sqrt(rank01(dca) * rank01(prior)),
                "dcagated": dca + np.where(prior > np.median(prior),
                                            np.abs(zscore(prior)),
                                            -np.abs(zscore(prior))),
            }
            row = {"target": name, "fold": fold}
            for m, s in scores.items():
                p5 = ks.precision_at_k(s, lab, K)
                auc = ks.auc_mann_whitney(s, lab)
                agg[m]["prec@L/5"].append(p5)
                agg[m]["auc"].append(auc)
                row[f"{m}_prec"] = round(p5, 4)
                row[f"{m}_auc"] = round(auc, 4)
            per_protein.append(row)

        print(f"fold {fold}: " + " ".join(
            f"{m}={np.mean(agg[m]['prec@L/5'][-len(test):]):.3f}"
            for m in methods), flush=True)

    summary = {}
    for m in methods:
        p5 = np.array(agg[m]["prec@L/5"])
        a = np.array(agg[m]["auc"])
        summary[m] = {"mean_prec_L5": float(p5.mean()),
                      "sd_prec_L5": float(p5.std()),
                      "mean_auc": float(a.mean())}

    # paired Wilcoxon: each fusion vs both parents
    stats = {}
    for fusion in ("ranksum", "zsum", "product", "geometric", "dcagated"):
        for parent in ("rawDCA", "circH"):
            x = np.array(agg[fusion]["prec@L/5"])
            y = np.array(agg[parent]["prec@L/5"])
            w = wilcoxon(x, y)
            stats[f"{fusion}_vs_{parent}"] = {
                "p_value": float(w.pvalue),
                "mean_diff": float(np.mean(x - y))}
    summary["paired_stats"] = stats

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/"path_l_rank_fusion.json").write_text(json.dumps(
        {"summary": summary, "per_fold_sizes":
         {m: len(v["prec@L/5"]) for m, v in agg.items()}}, indent=1))
    print(json.dumps({m: summary[m] for m in methods if m != "paired_stats"},
                     indent=1))
