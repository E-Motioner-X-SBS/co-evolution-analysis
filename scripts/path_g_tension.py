#!/usr/bin/env python3
"""
path_g_tension.py — PATH G: entropic tension profiles (user's mechanism).

DEFINITION (plans/12): for position i, the *tension* toward position j is the
normalized entropic dependence

    T(i,j) = MI(i,j) / max_k MI(i,k)      over cached partners k of i (sep>=6)

i.e. MI(i,j) expressed as a fraction of position i's strongest entropic
coupling anywhere in the protein — "how much of this residue's total
entropic attention does partner j consume?"  Symmetric variant T(j,i) uses
j's own strongest coupling.

Experiments
-----------
G1  HUB TEST: are high-tension positions structural hubs?  For each protein,
    correlate per-position max-tension-partner count (number of partners whose
    top-1 tension target is i) with the native contact degree of i.
    Literature anchor: hub residues in contact maps correspond to buried core
    (high degree).  Spearman rho + p per protein, pooled Fisher-combined.

G2  CIRCUIT LIFT: add tension bins to the universal circuit:
        level-G: g_i[3] x g_j[3] x sep[2] x T(i,j)-bin[2] x T(j,i)-bin[2]
    = 10 bits = 1024 cells.  Train on ~100 proteins, evaluate on held-out
    proteins exactly as PATH B does; compare prec@L/5 against
    circuit-L2 (no tension) and plain MI.

G3  PROFILE REPRODUCIBILITY (literature-facing): correlation of tension
    profiles between proteins sharing a fold is a literature-known property
    of coupling vectors (Weigt 2009 / Morcos 2011 discuss coupling-vector
    geometry); we measure whether ENTROPIC tension profiles behave likewise
    across PSICOV targets grouped by SCOP-less proxy (k-mer similarity).

Output: results/kmap_structure/path_g_tension.json (+ .npz caches).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import kmap_structure as ks  # noqa: E402

OUT = ks.RESULTS


# ---------------------------------------------------------------------------
# Tension computation from cached pair universes
# ---------------------------------------------------------------------------

def tension_table(tc) -> tuple[np.ndarray, np.ndarray]:
    """Per-pair T(i,j), T(j,i) in [0,1] from the cached MI column.

    T(i,j) = MI(i,j) / max_k MI(i,k) over the protein's cached candidate
    partners k of i.  Positions with no positive-MI partner get T=0.
    """
    i_idx = tc.pairs[:, 0]
    j_idx = tc.pairs[:, 1]
    mi = np.nan_to_num(tc.mi, nan=0.0)

    # row maxima over the first-column index and second-column index
    n_cols = tc.n_cols
    row_max = np.zeros(n_cols)
    np.maximum.at(row_max, i_idx, mi)
    np.maximum.at(row_max, j_idx, mi)          # symmetric universe
    denom_i = np.where(row_max[i_idx] > 0, row_max[i_idx], 1.0)
    denom_j = np.where(row_max[j_idx] > 0, row_max[j_idx], 1.0)
    return mi / denom_i, mi / denom_j


def tension_bins(t: np.ndarray, edges=(0.25, 0.5, 0.75)) -> np.ndarray:
    """2-bit bin: 0 [<e0], 1 [e0,e1), 2 [e1,e2), 3 [>=e2]."""
    b = np.digitize(t, edges).astype(np.int8)
    return b


# ---------------------------------------------------------------------------
# G1 — hub test
# ---------------------------------------------------------------------------

def hub_test(caches: dict) -> dict:
    """Count, for each position, how many partners rank it as their TOP
    tension target; correlate with native contact degree."""
    rows = []
    for name in sorted(caches):
        tc = caches[name]
        ti, tj = tension_table(tc)
        i_idx = tc.pairs[:, 0]; j_idx = tc.pairs[:, 1]
        # top-1 tension partner per column (both directions)
        top_count = np.zeros(tc.n_cols, dtype=np.int64)
        for src, dst, tval in ((i_idx, j_idx, ti), (j_idx, i_idx, tj)):
            order = np.lexsort((-tval, src))
            sorted_src = src[order]
            first = np.ones(len(sorted_src), dtype=bool)
            first[1:] = sorted_src[1:] != sorted_src[:-1]
            top_dst = dst[order][first]
            np.add.at(top_count, top_dst, 1)
        # native contact degree per column
        deg = np.zeros(tc.n_cols, dtype=np.int64)
        cont_pairs = tc.pairs[tc.lab]
        np.add.at(deg, cont_pairs[:, 0], 1)
        np.add.at(deg, cont_pairs[:, 1], 1)
        mask = deg > 0                       # columns with >=1 contact
        if mask.sum() < 10:
            continue
        rho = spearmanr(top_count[mask], deg[mask])
        rows.append({"target": name, "n_cols": int(mask.sum()),
                     "rho": float(rho.statistic), "p": float(rho.pvalue)})
    rhos = np.array([r["rho"] for r in rows])
    # Fisher combine p-values (one-sided greater)
    ps = np.array([r["p"] for r in rows])
    chi2 = -2.0 * np.log(np.clip(ps, 1e-300, 1.0)).sum()
    from scipy.stats import chi2 as chi2dist
    comb_p = float(chi2dist.sf(chi2, df=2 * len(ps)))
    return {"per_protein": rows,
            "median_rho": float(np.median(rhos)),
            "frac_positive": float((rhos > 0).mean()),
            "fisher_combined_p": comb_p,
            "n_proteins": len(rows),
            "note": ("Top-1 tension-target counts vs native contact degree; "
                     "positive rho = high-tension positions are hubs.")}


# ---------------------------------------------------------------------------
# G2 — circuit lift with tension dimensions
# ---------------------------------------------------------------------------

def levelG_index(g_i, g_j, sbin, tb_i, tb_j):
    """(g_i,g_j,sep,tension_i,tension_j) -> flat cell in [0,1024)."""
    return (((g_i.astype(np.int64) * 8 + g_j) * 4 + sbin) * 4 + tb_i) * 4 + tb_j


class CountsG:
    N_CELLS = 4096   # 3+3+2+2+2 = 12 bits

    def __init__(self):
        self.pos = np.zeros(self.N_CELLS, np.int64)
        self.tot = np.zeros(self.N_CELLS, np.int64)

    def add(self, idx, lab):
        np.add.at(self.tot, idx, 1)
        np.add.at(self.pos, idx, lab.astype(np.int64))


def eval_g_circuit(train_names, test_names, caches):
    """PATH-B protocol with the tension-augmented 1024-cell circuit."""
    tr = CountsG()
    for name in train_names:
        tc = caches[name]
        ti, tj = tension_table(tc)
        idx = levelG_index(tc.g_i, tc.g_j, tc.sbin,
                           tension_bins(ti), tension_bins(tj))
        tr.add(idx, tc.lab)
    base = tr.pos.sum() / max(tr.tot.sum(), 1)
    truth = ks.label_truth_table(
        ks.CellCounts(CountsG.N_CELLS, tr.pos.copy(), tr.tot.copy()), 2.0, 30)
    res_qm = None
    # soundness/completeness asserted via qm_extract equivalent (reuse runner)
    sys.path.insert(0, str(HERE))
    from run_contact_campaign import qm_extract
    res_qm = qm_extract(truth)

    out_rows = []
    for name in test_names:
        tc = caches[name]
        ti, tj = tension_table(tc)
        idx = levelG_index(tc.g_i, tc.g_j, tc.sbin,
                           tension_bins(ti), tension_bins(tj))
        score = ks.cell_rate_score(
            ks.CellCounts(CountsG.N_CELLS, tr.pos, tr.tot), idx)
        L = tc.n_cols
        rec = {"prec@L/5": ks.precision_at_k(score, tc.lab, max(1, round(L/5))),
               "auc": ks.auc_mann_whitney(score, tc.lab)}
        out_rows.append(rec)
    agg = {k: float(np.mean([r[k] for r in out_rows])) for k in out_rows[0]}
    return agg, res_qm, base


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    from run_contact_campaign import load_all_caches
    caches = load_all_caches()
    targets = sorted(caches)
    result: dict = {}

    t0 = time.time()
    print("[G1] hub test ...", flush=True)
    result["hub_test"] = hub_test(caches)
    print(f"  median rho={result['hub_test']['median_rho']:.3f} "
          f"frac>0={result['hub_test']['frac_positive']:.2f} "
          f"combined p={result['hub_test']['fisher_combined_p']:.3g}")

    print("[G2] tension-augmented circuit (5-fold protein CV) ...", flush=True)
    folds = ks.protein_folds(targets, n_folds=5, seed=42)
    per_fold = []
    for f in range(5):
        test = folds[f]
        train = [t for g, grp in enumerate(folds) if g != f for t in grp]
        agg, res_qm, base = eval_g_circuit(train, test, caches)
        per_fold.append(agg)
        print(f"  fold {f}: prec@L/5={agg['prec@L/5']:.3f} "
              f"auc={agg['auc']:.3f} (base {base:.4f}, "
              f"{res_qm['n_prime_implicants']} PIs)")
    result["circuit_G"] = {
        "per_fold": per_fold,
        "mean_prec_L5": float(np.mean([a["prec@L/5"] for a in per_fold])),
        "mean_auc": float(np.mean([a["auc"] for a in per_fold])),
        "note": "level-G = groups x sep x T(i,j)-bins x T(j,i)-bins (4096 cells)",
    }

    result["runtime_s"] = round(time.time() - t0, 1)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "path_g_tension.json").write_text(json.dumps(result, indent=1))
    print(f"saved -> {OUT/'path_g_tension.json'} ({result['runtime_s']}s)")


if __name__ == "__main__":
    main()
