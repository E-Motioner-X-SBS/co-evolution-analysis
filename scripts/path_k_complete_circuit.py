#!/usr/bin/env python3
"""
path_k_complete_circuit.py — PATH K: the COMPLETE circuit.

The user's observation: our circuits encode the VARIABLE-position subgraph
(positive mutation K-maps).  The FLIPPED map (never-observed combinations =
negative selection) plus CONSERVED positions complete the picture.  This
script builds both polarities over one shared feature space and evaluates
the combined two-sided classifier.

Feature space (10 bits = 1024 cells):
    g_i[3] x g_j[3] x sep_bin[2] x cons_type[2]
where cons_type ∈ {var/var, var/cons, cons/var, cons/cons} and a column is
CONSERVED iff its Shannon entropy H < 0.3 (repo-wide convention; gaps
excluded).  Conserved columns carry the burial/core signal documented by
structural_profile.py (rho(entropy,burial) = -0.17).

Two ternary tables are labeled from the SAME pooled counts:
    POSITIVE table: cell ON  if rate >= T*base   and n >= n_min   (promotes)
                    DC       if n < n_min                        (untestable)
    FLIPPED  table: cell ON  if rate <= base/T   and n >= n_min   (forbids)
                    DC       if n < n_min
Combined decision (documented rule):
    call CONTACT    iff positive-cell ON
    call FORBIDDEN  iff flipped-cell ON
    abstain         otherwise
Metrics: positive precision/recall; flipped specificity + forbidden-hit rate;
definitive-call coverage; MCC over definitive calls.  Protein-level CV
identical to PATH B (folds seed=42).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import kmap_structure as ks  # noqa: E402

OUT = ks.RESULTS
CONS_H = 0.3          # conserved-column entropy threshold (repo convention)
N_CELLS = 1024        # 3+3+2+2 bits


def levelK_index(g_i, g_j, sbin, ctype) -> np.ndarray:
    """((g_i*8+g_j)*4+sbin)*4+ctype -> [0,1024)."""
    return (((g_i.astype(np.int64) * 8 + g_j) * 4 + sbin) * 4 + ctype)


def cons_type(pd_dense: np.ndarray) -> np.ndarray:
    """Per-column conservation flag (True = conserved, H < CONS_H)."""
    n, L = pd_dense.shape
    ent = np.zeros(L)
    for c in range(L):
        col = pd_dense[:, c]
        v = col[(col >= 0) & (col < ks.GAP_STATE)]
        if v.size == 0:
            continue
        cnt = np.bincount(v, minlength=ks.GAP_STATE).astype(np.float64)
        p = cnt[cnt > 0] / cnt.sum()
        ent[c] = -float(np.sum(p * np.log2(p)))
    return ent < CONS_H


def ctype_pairs(cons: np.ndarray, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    """2-bit pair type: 0 var/var, 1 var/cons, 2 cons/var, 3 cons/cons."""
    ci = cons[i].astype(np.int64)
    cj = cons[j].astype(np.int64)
    return (ci * 2 + cj)


class CountsK:
    def __init__(self):
        self.pos = np.zeros(N_CELLS, np.int64)
        self.tot = np.zeros(N_CELLS, np.int64)

    def add(self, idx, lab):
        np.add.at(self.tot, idx, 1)
        np.add.at(self.pos, idx, lab.astype(np.int64))


def label_two_sided(counts: CountsK, t_mult: float, n_min: int):
    base = counts.pos.sum() / max(counts.tot.sum(), 1)
    n = counts.tot
    ok = n >= n_min
    rate = counts.pos / np.maximum(n, 1)
    pos = np.where(ok & (rate >= t_mult * base), 1, 0).astype(np.int8)
    flip = np.where(ok & (rate <= base / t_mult), 1, 0).astype(np.int8)
    return pos, flip, float(base)


def eval_fold(train_names, test_names, caches, t_mult=2.0, n_min=30):
    tr = CountsK()
    cons_cache = {}
    for name in train_names:
        tc = caches[name]
        if name not in cons_cache:
            pd = ks.load_protein(name)
            assert pd.dense is not None
            cons_cache[name] = cons_type(pd.dense)
        cons = cons_cache[name]
        idx = levelK_index(tc.g_i, tc.g_j, tc.sbin,
                           ctype_pairs(cons, tc.pairs[:, 0], tc.pairs[:, 1]))
        tr.add(idx, tc.lab)

    pos_tab, flip_tab, base = label_two_sided(tr, t_mult, n_min)

    rows = []
    for name in test_names:
        tc = caches[name]
        if name not in cons_cache:
            pd = ks.load_protein(name)
            assert pd.dense is not None
            cons_cache[name] = cons_type(pd.dense)
        cons = cons_cache[name]
        idx = levelK_index(tc.g_i, tc.g_j, tc.sbin,
                           ctype_pairs(cons, tc.pairs[:, 0], tc.pairs[:, 1]))
        pcall = pos_tab[idx] == 1          # predicted CONTACT
        fcall = flip_tab[idx] == 1         # predicted FORBIDDEN (no-contact)
        lab = tc.lab
        tp = int((pcall & lab).sum()); fp = int((pcall & ~lab).sum())
        fn = int((~pcall & lab).sum()); tn_def = int((fcall & ~lab).sum())
        ffor = int((fcall & lab).sum())    # forbidden-call but actually contact
        definitive = int((pcall | fcall).sum())
        prec = tp / (tp + fp) if tp + fp else 0.0
        spec = tn_def / (tn_def + ffor) if tn_def + ffor else 0.0
        # MCC over DEFINITIVE calls only (abstentions truly excluded):
        # negatives = flipped-called pairs that are indeed non-contacts
        mcc_v = ks.mcc(tp, fp, tn_def, ffor)
        rows.append({
            "target": name,
            "pos_precision": prec,
            "pos_recall": tp / max(tp + fn, 1),
            "flip_specificity": spec,
            "flip_forbidden_hit_rate": ffor / max(int(fcall.sum()), 1),
            "coverage": definitive / max(len(lab), 1),
            "mcc_definitive": mcc_v,
            "base_rate": float(lab.mean()),
        })
    agg = {k: float(np.mean([r[k] for r in rows]))
           for k in rows[0] if isinstance(rows[0][k], (int, float))}
    return agg, pos_tab, flip_tab, base, rows


def main():
    from run_contact_campaign import load_all_caches
    caches = load_all_caches()
    targets = sorted(caches)
    folds = ks.protein_folds(targets, n_folds=5, seed=42)
    per_fold, all_rows = [], []
    for f in range(5):
        test = folds[f]
        train = [t for g, grp in enumerate(folds) if g != f for t in grp]
        t0 = time.time()
        agg, pos_tab, flip_tab, base, rows = eval_fold(train, test, caches)
        per_fold.append(agg)
        print(f"fold {f} ({time.time()-t0:.0f}s): "
              f"pos-prec={agg['pos_precision']:.3f} "
              f"pos-rec={agg['pos_recall']:.3f} | "
              f"flip-spec={agg['flip_specificity']:.3f} "
              f"flip-hits-contact={agg['flip_forbidden_hit_rate']:.3f} | "
              f"coverage={agg['coverage']:.3f} "
              f"MCC={agg['mcc_definitive']:.3f}", flush=True)
        for r in rows:
            r["fold"] = f
        all_rows.extend(rows)
    mean = {k: float(np.mean([a[k] for a in per_fold])) for k in per_fold[0]}
    summary = {"per_fold": per_fold, "mean": mean,
               "note": ("Complete circuit: positive (promotes contact) + "
                        "flipped (forbids) tables over groups x sep x "
                        "conservation-type; two-sided decisions with "
                        "abstention.")}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "path_k_complete_circuit.json").write_text(json.dumps(
        {"summary": summary, "per_protein": all_rows}, indent=1))
    print(json.dumps(mean, indent=1))


if __name__ == "__main__":
    main()
