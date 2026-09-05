#!/usr/bin/env python3
"""
run_contact_campaign.py — The K-map ↔ 3D-structure encoding campaign runner.
=============================================================================
Campaign: plans/11-kmap-structure-encoding.md + plans/11b-strategy.md.

Stages (each idempotent, resumable, writes JSON under results/kmap_structure/):

  cache   Stage A: per-protein feature caches (pairs, group/identity features,
                   native-contact labels from LITERATURE PDBs, plain-MI scores,
                   perplexity-ratio scores).  The expensive part; run once.
  cv      Stage B / PATH B: universal contact circuits.
                   Protein-level CV: pooled K-map cell counts on train proteins
                   -> documented ternary labeling -> Quine-McCluskey -> Boolean
                   circuits; evaluated on held-out proteins against baselines:
                   random(base rate), separation-only ablation, plain MI,
                   perplexity ratio, and PSICOV's PUBLISHED predictions
                   (raw/con/*.out — the literature contact map).
  sweep   Stage C: labeling-rule sensitivity sweep (t_mult x n_min grid) on the
                   fixed ~100/~49 split — conclusions must be stable (gate G4).
  gray    Stage D / PATH E: Gray-adjacency vs native contacts at proper power;
                   background = within-protein, within-sep-bin resampling of
                   non-contact consensus residue pairs; bootstrap CI.
  mapqm   Stage E / PATH F (exploratory): treat a protein's native CONTACT MAP
                   itself as a Boolean function (64x64 padded, don't-care
                   outside the candidate universe), exact-QM it, compare cover
                   complexity vs separation-band label shuffles.

Usage:
  PY scripts/run_contact_campaign.py --stage cache [--limit N] [--force-cache]
  PY scripts/run_contact_campaign.py --stage cv [--folds 5]
  PY scripts/run_contact_campaign.py --stage sweep
  PY scripts/run_contact_campaign.py --stage gray [--boot 2000]
  PY scripts/run_contact_campaign.py --stage mapqm

Environment: OMP_PROC_BIND=FALSE recommended (repo doctrine).
Every trained circuit is soundness/completeness asserted at extraction time
(the property Lean-formalized as ContactCircuits.qm_cover_sound).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon, binomtest, chisquare

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import kmap_structure as ks  # noqa: E402

OUT = ks.RESULTS
FEAT_DIR = OUT / "features"

#: Ranked-metric definitions (label, fraction-of-L).  K = round(frac * L).
RANKED_KS = (("prec@L/5", 1 / 5), ("prec@L/2", 1 / 2), ("prec@L", 1.0))

#: Method namespaces that appear in every per-protein record.
METHOD_KEYS = ("L1", "L2", "sep_only", "mi", "pp", "psicov_published")


# ===========================================================================
# STAGE A — feature caches
# ===========================================================================

def _entropy_vector(dense: np.ndarray) -> np.ndarray:
    """Per-column Shannon entropy (bits), gaps excluded — mirrors shared module."""
    ent = np.zeros(dense.shape[1])
    for c in range(dense.shape[1]):
        col = dense[:, c]
        valid = col[(col >= 0) & (col < ks.GAP_STATE)]
        if valid.size == 0:
            continue
        cnt = np.bincount(valid, minlength=ks.GAP_STATE).astype(np.float64)
        p = cnt[cnt > 0] / cnt.sum()
        ent[c] = -float(np.sum(p * np.log2(p)))
    return ent


def stage_cache(limit: int = 0, force: bool = False) -> None:
    """Build one .npz per target with everything downstream needs. Resumable."""
    FEAT_DIR.mkdir(parents=True, exist_ok=True)
    targets = ks.list_targets()
    if limit:
        targets = targets[:limit]
    manifest_path = OUT / "cache_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    t_start = time.time()
    for n, target in enumerate(targets, 1):
        dest = FEAT_DIR / f"{target}.npz"
        if dest.exists() and not force:
            manifest.setdefault(target, "cached")
            continue
        t0 = time.time()
        pd = ks.load_protein(target)
        assert pd.dense is not None, f"{target}: dense matrix required"
        pairs = ks.pair_universe(pd)
        feats = ks.pair_features(pd, pairs)
        mi = ks.mi_scores(pd, pairs)
        ent = _entropy_vector(pd.dense)
        pp = np.empty(len(pairs), dtype=np.float64)
        for k, (i, j) in enumerate(pairs):
            r = ks.perplexity_ratio_columns(pd.dense, int(i), int(j), ent)
            pp[k] = np.nan if r is None else r
        np.savez_compressed(
            dest,
            pairs=pairs, sbin=feats["sbin"],
            aa_i=feats["aa_i"], aa_j=feats["aa_j"],
            g_i=feats["g_i"], g_j=feats["g_j"],
            is_contact=feats["is_contact"], gray_h=feats["gray_h"],
            mi=mi, pp=pp,
            meta=json.dumps({"target": target, "n_cols": pd.n_cols,
                             "n_seqs": pd.n_seqs, "identity_ok": pd.identity_ok,
                             "notes": pd.notes}),
        )
        manifest[target] = (f"ok pairs={len(pairs)} "
                            f"contacts={int(feats['is_contact'].sum())} "
                            f"seqs={pd.n_seqs} width={pd.n_cols} "
                            f"{time.time()-t0:.1f}s")
        print(f"[{n:3d}/{len(targets)}] {target}: {manifest[target]}", flush=True)
    manifest_path.write_text(json.dumps(manifest, indent=1))
    print(f"cache complete in {time.time()-t_start:.1f}s -> {FEAT_DIR}")


# ===========================================================================
# Cache loading
# ===========================================================================

class TargetCache:
    """In-memory view of one cached target (all arrays aligned to tc.pairs)."""

    __slots__ = ("target", "pairs", "sbin", "aa_i", "aa_j", "g_i", "g_j",
                 "lab", "gray_h", "mi", "pp", "meta")

    def __init__(self, path: Path):
        z = np.load(path, allow_pickle=False)
        self.target = path.stem
        self.pairs = z["pairs"]
        self.sbin = z["sbin"]; self.aa_i = z["aa_i"]; self.aa_j = z["aa_j"]
        self.g_i = z["g_i"]; self.g_j = z["g_j"]
        self.lab = z["is_contact"].astype(bool)
        self.gray_h = z["gray_h"]
        self.mi = z["mi"]; self.pp = z["pp"]
        self.meta = json.loads(str(z["meta"]))

    @property
    def n_cols(self) -> int:
        return int(self.meta["n_cols"])

    def level1_cells(self) -> np.ndarray:
        return ks.level1_cell_index(self.g_i, self.g_j, self.sbin)

    def level2_cells(self, b: int) -> np.ndarray:
        m = self.sbin == b
        return ks.level2_cell_index(self.aa_i[m], self.aa_j[m])


def load_all_caches() -> dict:
    out = {}
    for p in sorted(FEAT_DIR.glob("*.npz")):
        out[p.stem] = TargetCache(p)
    if not out:
        raise SystemExit("no caches found — run --stage cache first")
    return out


# ===========================================================================
# Circuit training / evaluation machinery (PATH B)
# ===========================================================================

def train_counts(train_names, caches):
    """Pooled K-map cell counters over TRAIN proteins only (no leakage)."""
    c1 = ks.new_level1_counts()
    c2 = ks.new_level2_counts()
    for name in train_names:
        tc = caches[name]
        c1.add(tc.level1_cells(), tc.lab)
        for b in range(ks.N_SEP_BINS):
            m = tc.sbin == b
            if m.any():
                c2[b].add(tc.level2_cells(b), tc.lab[m])
    return c1, c2


def sep_only_counts(train_names, caches):
    """Ablation counter: cells = separation bin ONLY (no residue information)."""
    c = ks.CellCounts(ks.N_SEP_BINS, np.zeros(ks.N_SEP_BINS, np.int64),
                      np.zeros(ks.N_SEP_BINS, np.int64))
    for name in train_names:
        tc = caches[name]
        c.add(tc.sbin.astype(np.int64), tc.lab)
    return c


def qm_extract(truth: np.ndarray) -> dict:
    """QM-minimize a labeled table, asserting SOUNDNESS + COMPLETENESS:

      completeness — every ON cell is covered by >= 1 prime implicant;
      soundness    — no prime implicant matches any OFF cell.

    (Don't-care cells may or may not be covered — standard QM semantics.)
    This is exactly the property formalized in Lean as
    ContactCircuits.qm_cover_complete / qm_implicants_avoid_offset.
    """
    res = ks.qm_minimize(truth)
    n_bits = int(np.log2(len(truth)))
    preds = []
    for pi in res["prime_implicants"]:
        vals, mask = pi["values"], pi["mask"]

        def pred(cell, _vals=tuple(vals), _mask=tuple(mask), _nb=n_bits):
            for b in range(_nb):
                if not _mask[b] and ((cell >> (_nb - 1 - b)) & 1) != _vals[b]:
                    return False
            return True
        preds.append(pred)
    for c in np.flatnonzero(truth == 1):
        assert any(p(int(c)) for p in preds), f"ON cell {c} uncovered"
    for p in preds:
        for c in np.flatnonzero(truth == 0):
            assert not p(int(c)), f"PI matches OFF cell {c}"
    return res


def matched_codes(values, mask, width: int):
    """All `width`-bit field values consistent with an implicant (None=all DC).

    Used for human-readable level-2 rules; padding codes >=20 cannot appear
    (soundness check t19 locks this for real tables).
    """
    if all(mask):
        return None
    out = []
    for x in range(1 << width):
        ok = True
        for k in range(width):
            if not mask[k] and ((x >> (width - 1 - k)) & 1) != values[k]:
                ok = False
                break
        if ok:
            out.append(x)
    return out


def psicov_published_scores(tc: TargetCache) -> np.ndarray:
    """Literature baseline: PSICOV published score per candidate pair.

    raw/con/<target>.out enumerates pairs sorted by prediction score; we map
    each candidate pair to its published score (1-based->0-based columns).
    Missing entries (none observed) would sort last via -inf.
    """
    con = ks.parse_con(ks.RAW / "con" / f"{tc.target}.out")
    lut = {(i - 1, j - 1): s for (i, j, s) in con}
    return np.array([lut.get((int(a), int(b)), -np.inf) for a, b in tc.pairs],
                    dtype=np.float64)


def eval_ranked_on_target(score: np.ndarray, tc: TargetCache) -> dict:
    """precision@K (K = L/5, L/2, L over the protein's candidate universe)
    plus full-universe AUC.  Deterministic; NaN-safe inputs not expected."""
    out = {}
    for name, frac in RANKED_KS:
        k = max(1, int(round(frac * tc.n_cols)))
        out[name] = ks.precision_at_k(score, tc.lab, k)
    out["auc"] = ks.auc_mann_whitney(score, tc.lab)
    return out


def evaluate_fold(train_names, test_names, caches, t_mult: float, n_min: int) -> dict:
    """Train circuits on train proteins; evaluate all methods on test proteins."""
    c1, c2 = train_counts(train_names, caches)
    cs = sep_only_counts(train_names, caches)

    tt1 = ks.label_truth_table(c1, t_mult=t_mult, n_min=n_min)
    tt2 = [ks.label_truth_table(c2[b], t_mult=t_mult, n_min=n_min)
           for b in range(ks.N_SEP_BINS)]
    res1 = qm_extract(tt1)
    res2 = [qm_extract(tt2[b]) for b in range(ks.N_SEP_BINS)]

    rows = []
    for name in test_names:
        tc = caches[name]
        idx1 = tc.level1_cells()
        pred1 = (tt1[idx1] == 1).astype(int)
        score_l1 = ks.cell_rate_score(c1, idx1)

        pred_l2 = np.zeros(len(tc.pairs), dtype=int)
        score_l2 = np.zeros(len(tc.pairs))
        for b in range(ks.N_SEP_BINS):
            m = tc.sbin == b
            if not m.any():
                continue
            cells_b = tc.level2_cells(b)
            score_l2[m] = ks.cell_rate_score(c2[b], cells_b)
            pred_l2[m] = (tt2[b][cells_b] == 1)

        rec = {
            "target": name, "n_pairs": len(tc.pairs),
            "n_contacts": int(tc.lab.sum()),
            "L1": {"binary": ks.binary_metrics(pred1, tc.lab),
                   "ranked": eval_ranked_on_target(score_l1, tc)},
            "L2": {"binary": ks.binary_metrics(pred_l2, tc.lab),
                   "ranked": eval_ranked_on_target(score_l2, tc)},
            "sep_only": {"ranked": eval_ranked_on_target(
                ks.cell_rate_score(cs, tc.sbin.astype(np.int64)), tc)},
            "mi": {"ranked": eval_ranked_on_target(
                np.nan_to_num(tc.mi, nan=-1e9), tc)},
            "pp": {"ranked": eval_ranked_on_target(
                np.nan_to_num(tc.pp, nan=-1e9), tc)},
            "psicov_published": {"ranked": eval_ranked_on_target(
                psicov_published_scores(tc), tc)},
        }
        rows.append(rec)

    agg = _aggregate(rows)
    agg["circuits"] = {
        "level1": _circuit_summary_level1(tt1, res1, c1),
        "level2_per_bin": [_circuit_summary_level2(tt2[b], res2[b], c2[b], b)
                           for b in range(ks.N_SEP_BINS)],
    }
    agg["train_n"] = len(train_names)
    agg["test_n"] = len(test_names)
    agg["params"] = {"t_mult": t_mult, "n_min": n_min}
    agg["per_protein"] = rows
    return agg


def _circuit_summary_level1(tt, res, counts) -> dict:
    rules = []
    for pi in res["prime_implicants"]:
        d = ks.decode_level1_implicant(pi["values"], pi["mask"])
        rules.append({
            "groups_i": None if d["groups_i"] is None else
                        [ks.GROUP_NAMES[g] for g in d["groups_i"]],
            "groups_j": None if d["groups_j"] is None else
                        [ks.GROUP_NAMES[g] for g in d["groups_j"]],
            "sep_bins": d["seps"],
            "n_dontcare": pi["n_dontcares"],
            "essential": pi in res["essential_prime_implicants"],
        })
    return {"cells": len(tt), "on_cells": int((tt == 1).sum()),
            "dontcare_cells": int((tt == -1).sum()),
            "base_rate": counts.base_rate(),
            "n_prime_implicants": res["n_prime_implicants"],
            "n_essential": res["n_essential"], "rules": rules}


def _circuit_summary_level2(tt, res, counts, b: int) -> dict:
    rules = []
    for pi in res["prime_implicants"]:
        rows_m = matched_codes(pi["values"][:5], pi["mask"][:5], 5)
        cols_m = matched_codes(pi["values"][5:], pi["mask"][5:], 5)
        rules.append({
            "row_codes": None if rows_m is None else [c for c in rows_m if c < 20],
            "col_codes": None if cols_m is None else [c for c in cols_m if c < 20],
            "row_aa": None if not rows_m else [ks.CODE_TO_AA[c] for c in rows_m if c < 20],
            "col_aa": None if not cols_m else [ks.CODE_TO_AA[c] for c in cols_m if c < 20],
            "n_dontcare": pi["n_dontcares"],
            "essential": pi in res["essential_prime_implicants"],
        })
    return {"sep_bin": b, "sep_range": _bin_range(b),
            "cells": len(tt), "on_cells": int((tt == 1).sum()),
            "dontcare_cells": int((tt == -1).sum()),
            "base_rate": counts.base_rate(),
            "n_prime_implicants": res["n_prime_implicants"],
            "n_essential": res["n_essential"], "rules": rules}


def _bin_range(b: int) -> str:
    lo = ks.MIN_SEPARATION if b == 0 else ks.SEP_BIN_EDGES[b - 1]
    hi = ks.SEP_BIN_EDGES[b] if b < len(ks.SEP_BIN_EDGES) else "inf"
    return f"[{lo},{hi})"


# ---------------------------------------------------------------------------
# aggregation
# ---------------------------------------------------------------------------

def _agg_metric_block(blocks: list) -> dict:
    """Mean/sd/median/n(+n_nan) over per-protein scalar metric values."""
    out = {}
    names = sorted({k for blk in blocks for k in blk})
    for name in names:
        vals = [blk[name] for blk in blocks if name in blk]
        arr = np.array([v for v in vals if isinstance(v, (int, float))
                        and not isinstance(v, bool)], dtype=np.float64)
        if arr.size == 0:
            continue
        finite = arr[np.isfinite(arr)]
        if finite.size == 0:
            out[name] = {"mean": None, "sd": None, "median": None,
                         "n": 0, "n_nonfinite": int(arr.size)}
            continue
        out[name] = {"mean": float(finite.mean()), "sd": float(finite.std()),
                     "median": float(np.median(finite)), "n": int(finite.size),
                     "n_nonfinite": int(arr.size - finite.size)}
    return out


def _aggregate(rows: list) -> dict:
    """Aggregate per-protein records -> {method: {binary:{...}, ranked:{...}}}."""
    agg = {}
    for key in METHOD_KEYS:
        entry = {}
        for sub in ("binary", "ranked"):
            blocks = [r[key][sub] for r in rows if sub in r.get(key, {})]
            if blocks:
                entry[sub] = _agg_metric_block(blocks)
        agg[key] = entry

    def col(key):     # per-protein ranked precision@L/5 (NaN where absent)
        return np.array([r[key].get("ranked", {}).get("prec@L/5", np.nan)
                         for r in rows], dtype=np.float64)

    l2 = col("L2")
    stats = {}
    for other in ("psicov_published", "mi", "pp", "sep_only"):
        o = col(other)
        ok = ~(np.isnan(l2) | np.isnan(o))
        if ok.sum() >= 5 and np.any(l2[ok] != o[ok]):
            w = wilcoxon(l2[ok], o[ok])
            stats[f"L2_vs_{other}"] = {
                "n": int(ok.sum()),
                "mean_diff": float(np.mean(l2[ok] - o[ok])),
                "p_value": float(w.pvalue),  # type: ignore[attr-defined]
            }
    agg["paired_stats_prec_L5"] = stats
    return agg


def _print_fold_summary(res: dict) -> None:
    """Compact console table: mean precision@L/5 per method + circuit sizes."""

    def m(key):
        e = res.get(key, {}).get("ranked", {}).get("prec@L/5")
        return f"{e['mean']:.3f}" if e and e.get("mean") is not None else "n/a"

    b2 = res.get("L2", {}).get("binary", {}).get("mcc")
    e2 = res.get("L2", {}).get("binary", {}).get("enrichment")
    base = res.get("L2", {}).get("binary", {}).get("base_rate")
    print(f"  prec@L/5  circuit-L1={m('L1')}  circuit-L2={m('L2')}  "
          f"sep-only={m('sep_only')}  MI={m('mi')}  PP={m('pp')}  "
          f"PSICOV-published={m('psicov_published')}")
    mcc_s = f"{b2['mean']:.3f}" if b2 and b2.get("mean") is not None else "n/a"
    enr_s = f"{e2['mean']:.2f}x" if e2 and e2.get("mean") is not None else "n/a"
    base_s = f"{base['mean']:.3f}" if base and base.get("mean") is not None else "n/a"
    print(f"  L2 binary: MCC={mcc_s} enrichment={enr_s} (base {base_s})")
    c = res.get("circuits", {})
    if c:
        print(f"  circuits: L1 {c['level1']['n_prime_implicants']} PIs "
              f"({c['level1']['n_essential']} ess, {c['level1']['on_cells']} on) "
              f"| L2 PIs/bin={[x['n_prime_implicants'] for x in c['level2_per_bin']]}")


# ===========================================================================
# STAGE B driver — CV
# ===========================================================================

def stage_cv(folds: int = 5) -> None:
    caches = load_all_caches()
    targets = sorted(caches)
    print(f"cv: {len(targets)} cached targets, {folds}-fold protein-level CV",
          flush=True)
    fold_summaries = []
    fl = ks.protein_folds(targets, n_folds=folds, seed=42)
    for f in range(folds):
        test = fl[f]
        train = [t for g, grp in enumerate(fl) if g != f for t in grp]
        t0 = time.time()
        res = evaluate_fold(train, test, caches, t_mult=2.0, n_min=30)
        res["fold"] = f
        (OUT / f"cv_fold{f}.json").write_text(json.dumps(res, indent=1))
        print(f"fold {f}: train={len(train)} test={len(test)} "
              f"({time.time()-t0:.1f}s)", flush=True)
        _print_fold_summary(res)
        fold_summaries.append({k: v for k, v in res.items()
                               if k not in ("per_protein",)})
    # fixed ~100/~49 holdout: 3-fold scheme, folds 0+1 train, fold 2 test
    fl3 = ks.protein_folds(targets, n_folds=3, seed=42)
    split_res = evaluate_fold(fl3[0] + fl3[1], fl3[2], caches,
                              t_mult=2.0, n_min=30)
    split_res["split"] = "holdout: 3-fold scheme, folds 0+1 train, 2 test"
    (OUT / "holdout_split.json").write_text(json.dumps(split_res, indent=1))
    print("\n=== FIXED HOLDOUT (~100 train / ~49 test) ===")
    _print_fold_summary(split_res)
    slim = [{k: v for k, v in r.items() if k != "per_protein"}
            for r in fold_summaries]
    (OUT / "cv_summary.json").write_text(json.dumps(
        {"folds": slim,
         "holdout": {k: v for k, v in split_res.items()
                     if k != "per_protein"}}, indent=1))


# ===========================================================================
# STAGE C — sensitivity sweep (gate G4)
# ===========================================================================

def stage_sweep() -> None:
    caches = load_all_caches()
    targets = sorted(caches)
    fl3 = ks.protein_folds(targets, n_folds=3, seed=42)
    train, test = fl3[0] + fl3[1], fl3[2]
    grid = []
    for t_mult in (1.5, 2.0, 3.0):
        for n_min in (10, 30, 60):
            res = evaluate_fold(train, test, caches, t_mult=t_mult, n_min=n_min)
            mcc_e = res["L2"]["binary"].get("mcc", {})
            pl5 = res["L2"]["ranked"].get("prec@L/5", {})
            row = {
                "t_mult": t_mult, "n_min": n_min,
                "L2_binary": {k: v for k, v in res["L2"]["binary"].items()},
                "L2_prec_L5_mean": pl5.get("mean"),
                "L1_prec_L5_mean": res["L1"]["ranked"].get("prec@L/5", {}).get("mean"),
                "level1_on_cells": res["circuits"]["level1"]["on_cells"],
                "level1_PIs": res["circuits"]["level1"]["n_prime_implicants"],
                "level2_on_cells": [x["on_cells"]
                                    for x in res["circuits"]["level2_per_bin"]],
                "level2_PIs": [x["n_prime_implicants"]
                               for x in res["circuits"]["level2_per_bin"]],
            }
            grid.append(row)
            mcc_s = (f"{mcc_e['mean']:.3f}"
                     if mcc_e.get("mean") is not None else "n/a")
            pl5_s = (f"{pl5['mean']:.3f}"
                     if pl5.get("mean") is not None else "n/a")
            print(f"T={t_mult} nmin={n_min}: L2 MCC={mcc_s} "
                  f"prec@L/5={pl5_s} L1 on={row['level1_on_cells']} "
                  f"L2 on={row['level2_on_cells']}")
    (OUT / "sensitivity_sweep.json").write_text(json.dumps(grid, indent=1))
    print("sensitivity sweep saved")


# ===========================================================================
# STAGE D — PATH E: Gray adjacency at proper power
# ===========================================================================

def stage_gray(boot: int = 2000) -> None:
    """Is Gray-adjacency (h=1) enriched among native-contact residue pairs?

    Background controls composition AND correlation: for each contact pair we
    draw a replacement NON-contact pair from the SAME protein and SAME sep-bin
    (falling back to whole-protein non-contacts if a bin is empty there), then
    compute its Gray distance.  Two-sided binomial on h==1 counts; chi-square
    goodness-of-fit over h=0..5; protein-level bootstrap CI on the h==1 rate.
    """
    caches = load_all_caches()
    targets = sorted(caches)
    rng = np.random.default_rng(42)
    obs_h, bg_h = [], []
    for name in targets:
        tc = caches[name]
        cont_idx = np.flatnonzero(tc.lab)
        obs_h.append(tc.gray_h[cont_idx].astype(np.int16))
        bg = np.empty(len(cont_idx), dtype=np.int16)
        non = np.flatnonzero(~tc.lab)
        bins = tc.sbin[cont_idx]
        for b in range(ks.N_SEP_BINS):
            sel = bins == b
            if not sel.any():
                continue
            cand = non[tc.sbin[non] == b]
            cand = cand if cand.size else non
            draws = cand[rng.integers(0, cand.size, size=int(sel.sum()))]
            bg[sel] = [ks.gray_hamming(int(tc.aa_i[d]), int(tc.aa_j[d]))
                       for d in draws]
        bg_h.append(bg)
    obs_all = np.concatenate(obs_h)
    bg_all = np.concatenate(bg_h)
    n_contacts = len(obs_all)
    assert len(bg_all) == n_contacts
    print(f"gray test: {n_contacts} native-contact pairs "
          f"across {len(targets)} proteins")

    o1 = float((obs_all == 1).mean())
    b1 = float((bg_all == 1).mean()) if n_contacts else float("nan")
    enrich = o1 / b1 if b1 else float("nan")
    k1 = int((obs_all == 1).sum())
    pval = binomtest(k1, n_contacts, b1, alternative="two-sided").pvalue \
        if n_contacts else float("nan")

    obs_hist = np.array([(obs_all == d).sum() for d in range(6)], float)
    bg_hist = np.array([(bg_all == d).sum() for d in range(6)], float)
    exp_full = bg_hist / bg_hist.sum() * obs_hist.sum()
    exp_m, obs_m = exp_full.copy(), obs_hist.copy()
    while (exp_m < 5).any() and len(exp_m) > 2:      # valid-chi-square merging
        i = int(np.argmin(exp_m))
        j = i - 1 if i > 0 else 1
        exp_m[j] += exp_m[i]; obs_m[j] += obs_m[i]
        exp_m = np.delete(exp_m, i); obs_m = np.delete(obs_m, i)
    stat, chi_p = chisquare(f_obs=obs_m,
                            f_exp=exp_m * (obs_m.sum() / exp_m.sum()))

    sizes = [len(h) for h in obs_h]
    starts = np.cumsum([0] + sizes[:-1])
    boots = []
    for _ in range(boot):
        pick = rng.integers(0, len(targets), size=len(targets))
        sel = np.concatenate([np.arange(starts[i], starts[i] + sizes[i])
                              for i in pick]) if sizes else np.zeros(0, int)
        boots.append(float((obs_all[sel] == 1).mean()))
    ci_lo, ci_hi = np.percentile(boots, [2.5, 97.5])

    result = {
        "n_contacts": int(n_contacts),
        "n_proteins": len(targets),
        "observed_h1_rate": o1,
        "background_h1_rate": b1,
        "enrichment_h1": enrich,
        "binomial_two_sided_p": float(pval),
        "observed_h_distribution": obs_hist.tolist(),
        "background_h_distribution": bg_hist.tolist(),
        "chi_square_gof": {"stat": float(stat), "p": float(chi_p),
                           "bins_after_merge": int(len(obs_m))},
        "bootstrap_ci95_h1_rate": [float(ci_lo), float(ci_hi)],
        "note": ("Background = within-protein, within-sep-bin resampling of "
                 "NON-contact consensus residue pairs (composition+correlation "
                 "controlled). Spike comparison: 38 pairs, 1.18x, p=0.345."),
    }
    (OUT / "path_e_gray_adjacency.json").write_text(json.dumps(result, indent=1))
    print(json.dumps(result, indent=1))


# ===========================================================================
# STAGE E — PATH F: contact map AS a Boolean function
# ===========================================================================

def stage_mapqm(max_width: int = 62, shuffles: int = 5) -> None:
    """Exact-QM each small protein's native CONTACT MATRIX (64x64 padded).

    Truth-table semantics (locked here):
      ON  (1)  : candidate pair (sep>=6, both columns mapped) that IS a contact
      OFF (0)  : candidate pair that is NOT a contact
      DC  (-1) : everything else — diagonal, sep<6 band, padding border

    Control: permute contact labels WITHIN each separation bin of the same
    protein (density AND separation structure preserved exactly), re-minimize.
    If real maps need systematically fewer/larger implicants than their
    sep-preserving shadows, 3D geometry carries K-map-expressible regularity
    beyond density alone.
    """
    caches = load_all_caches()
    targets = [t for t in sorted(caches) if caches[t].n_cols <= max_width]
    print(f"mapqm: {len(targets)} proteins with L<={max_width}", flush=True)
    rng = np.random.default_rng(7)
    rows = []

    def build_truth(tc, labels):
        g = np.full((64, 64), -1, dtype=np.int8)          # DC default
        g[tc.pairs[:, 0], tc.pairs[:, 1]] = np.where(labels, 1, 0).astype(np.int8)
        return g.reshape(-1)

    for name in targets:
        tc = caches[name]
        truth_real = build_truth(tc, tc.lab)
        res_real = ks.qm_minimize(truth_real)
        on_real = int((truth_real == 1).sum())
        pis_shuf = []
        for _ in range(shuffles):
            labs = tc.lab.copy()
            for b in range(ks.N_SEP_BINS):
                idx = np.flatnonzero(tc.sbin == b)
                labs[idx] = labs[idx][rng.permutation(len(idx))]
            pis_shuf.append(ks.qm_minimize(build_truth(tc, labs))["n_prime_implicants"])
        comp_real = res_real["n_prime_implicants"] / on_real if on_real else float("nan")
        comp_shuf = [p / on_real for p in pis_shuf] if on_real else [np.nan] * shuffles
        rows.append({
            "target": name, "L": tc.n_cols, "on_cells": on_real,
            "PIs_real": res_real["n_prime_implicants"],
            "PIs_essential_real": res_real["n_essential"],
            "compression_real": comp_real,
            "PIs_shuffle": pis_shuf,
            "compression_shuffle_mean": float(np.nanmean(comp_shuf)),
            "compression_shuffle_sd": float(np.nanstd(comp_shuf)),
        })
        print(f"  {name}: L={tc.n_cols} on={on_real} "
              f"PIs={rows[-1]['PIs_real']} comp={comp_real:.3f} "
              f"shuffle={rows[-1]['compression_shuffle_mean']:.3f}"
              f"+-{rows[-1]['compression_shuffle_sd']:.3f}", flush=True)
    real = np.array([r["compression_real"] for r in rows])
    shuf = np.array([r["compression_shuffle_mean"] for r in rows])
    ok = np.isfinite(real) & np.isfinite(shuf)
    w = wilcoxon(real[ok], shuf[ok]) if ok.sum() >= 5 else None
    summary = {
        "n_proteins": len(rows),
        "max_width": max_width,
        "mean_compression_real": float(np.nanmean(real)),
        "mean_compression_shuffled": float(np.nanmean(shuf)),
        "wilcoxon_p": float(w.pvalue) if w else None,  # type: ignore[attr-defined]
        "per_protein": rows,
        "note": ("Compression = #prime implicants / #contact cells. Controls "
                 "shuffle labels within separation bands (density+sep kept)."),
    }
    (OUT / "path_f_map_qm.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k != "per_protein"},
                     indent=1))


# ===========================================================================
# main
# ===========================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", required=True,
                    choices=["cache", "cv", "sweep", "gray", "mapqm"])
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--limit", type=int, default=0,
                    help="cache stage only: process first N targets")
    ap.add_argument("--force-cache", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    FEAT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if args.stage == "cache":
        stage_cache(limit=args.limit, force=args.force_cache)
    elif args.stage == "cv":
        stage_cv(folds=args.folds)
    elif args.stage == "sweep":
        stage_sweep()
    elif args.stage == "gray":
        stage_gray(boot=args.boot)
    elif args.stage == "mapqm":
        stage_mapqm()
    print(f"\nstage '{args.stage}' done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
