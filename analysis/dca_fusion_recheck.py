#!/usr/bin/env python3
"""Does the K-map prior add anything to a working mfDCA?  (PSICOV150, 5-fold CV)

The paper's rank-fusion result (scripts/path_l_rank_fusion.py: ranksum 0.204
vs rawDCA 0.169 prec@L/5) fused the K-map cell-rate prior with
path_h_dca_circuit.mfdca_di's APC-corrected Frobenius score.  That score is
far below what mean-field DCA reaches on PSICOV150 (corrected DI: 0.566,
analysis/mfdca_psicov_check.py).  Two reasons, both in mfdca_di:

  * pair frequencies are weighted by W^3: the Gram matrix is built from
    A = (one-hot * W) * sqrt(W), so A A^T sums W_n^3, while the single-site
    frequencies use W_n.  The covariance mixes two different weightings.
  * the Frobenius norm is taken in the gauge left by the inversion (one amino
    acid as the zero reference), not the zero-sum gauge that the APC-Frobenius
    score assumes (Ekeberg et al. 2013).

This script scores every cached candidate pair of every target with

  pathh_fn_apc  mfdca_di's score, same arithmetic (torch), i.e. the published
                "rawDCA" input; checked against mfdca_di itself (--check)
  di            direct information, dca_mf_analysis with the corrected sign
  fn_apc_zs     APC-corrected Frobenius norm in the zero-sum gauge

and repeats path_l's evaluation for each one: the DCA score alone, the
K-map prior built from it (cell = residue groups x separation x DCA
quartile, trained on the other folds), and their rank-sum fusion.  MIp
(APC-corrected MI) is the sequence-only baseline.  Statistics: paired
Wilcoxon over the 150 proteins and 95% protein-bootstrap intervals of the
mean differences.

It also writes results/kmap_structure/di/<target>.npz (path_h's cache format,
pathh_fn_apc) so path_h_dca_circuit.py --stage cv and path_l_rank_fusion.py
can be re-run unchanged.

Output: analysis/dca_fusion_recheck.json
Run:    python analysis/dca_fusion_recheck.py --stage scores   (about 30 min, CPU)
        python analysis/dca_fusion_recheck.py --stage eval
        python analysis/dca_fusion_recheck.py --check          (mfdca_di parity)
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
os.environ.setdefault("COEVO_RESULTS", str(REPO / "results" / "mfdca_psicov_check"))

import numpy as np  # noqa: E402
import torch  # noqa: E402
from scipy.stats import spearmanr, wilcoxon  # noqa: E402

import dca_mf_analysis as d  # noqa: E402
import kmap_structure as ks  # noqa: E402
import path_h_dca_circuit as ph  # noqa: E402
import path_l_rank_fusion as pl  # noqa: E402
from run_contact_campaign import load_all_caches  # noqa: E402

CACHE = ks.RESULTS / "dca_recheck"
VARIANTS = ("pathh_fn_apc", "di", "fn_apc_zs")
MAX_SEQ = 5000          # subsample for the corrected scores (as mfdca_psicov_check)
OUT = REPO / "analysis" / "dca_fusion_recheck.json"


# ---------------------------------------------------------------------------
# scores
# ---------------------------------------------------------------------------

def pathh_fn_apc(dense: np.ndarray) -> tuple[np.ndarray, float]:
    """path_h_dca_circuit.mfdca_di's returned score, same arithmetic, in torch.

    Kept bug-for-bug (W^3 pair weighting, non-zero-sum gauge, +-30 clip,
    APC over a Frobenius matrix that includes its diagonal) because the
    point is to reproduce the published rawDCA input.
    """
    N, L = dense.shape
    if N > ph.MAX_SEQ:
        keep = np.sort(np.random.default_rng(42).choice(N, ph.MAX_SEQ, replace=False))
        dense = dense[keep]
        N = ph.MAX_SEQ
    Q = ph.Q
    X = torch.from_numpy(dense.astype(np.int64))
    valid = X < Q
    W = torch.ones(N, dtype=torch.float64)
    if N > 1:
        counts = torch.ones(N, dtype=torch.float64)
        thr = int(0.2 * L)
        for s in range(0, N, 256):
            blk = X[s:s + 256]
            same = (blk[:, None, :] == X[None, :, :]) & valid[None, :, :]
            counts[s:s + 256] = ((~same).sum(dim=2) <= thr).sum(dim=1).double() + 1
        W = 1.0 / counts
    Meff = W.sum()
    oh = torch.zeros(N, L, Q, dtype=torch.float64)
    oh.scatter_(2, X.clamp(0, Q - 1).unsqueeze(-1), valid.double().unsqueeze(-1))
    Vw = oh * W[:, None, None]
    Pi = 0.5 * (Vw.sum(dim=0) / Meff) + 0.5 / Q
    A = (Vw * W.sqrt()[:, None, None]).reshape(N, L * Q).T
    G = (A @ A.T / Meff).reshape(L, Q, L, Q)
    Pij = 0.5 * G + 0.5 * G.permute(2, 3, 0, 1)
    Pij = 0.5 * Pij + 0.5 / (Q * Q)
    eye = torch.eye(Q, dtype=torch.float64)
    for l in range(L):
        Pij[l, :, l, :] = Pi[l][:, None] * eye
    Pic = Pi[:, :Q - 1]
    C = Pij[:, :Q - 1, :, :Q - 1] - Pic[:, :, None, None] * Pic[None, None, :, :]
    Cm = C.reshape(L * (Q - 1), L * (Q - 1))
    Cm = 0.5 * (Cm + Cm.T)
    invC = torch.linalg.inv(Cm + ph.REG * torch.eye(Cm.shape[0], dtype=torch.float64))
    E = (-invC).reshape(L, Q - 1, L, Q - 1).permute(0, 2, 1, 3).clamp(-30.0, 30.0)
    frob = E.pow(2).sum(dim=(2, 3)).sqrt()
    apc = frob - frob.mean(dim=1)[:, None] * frob.mean(dim=0)[None, :] / frob.mean()
    apc.fill_diagonal_(0.0)
    return apc.numpy(), float(Meff)


def corrected_scores(dense: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """(DI, zero-sum-gauge APC-Frobenius, Meff) from dca_mf_analysis."""
    dense = dense.astype(np.int64)
    if dense.shape[0] > MAX_SEQ:
        keep = np.sort(np.random.default_rng(42).choice(dense.shape[0], MAX_SEQ, replace=False))
        dense = dense[keep]
    X = np.where(dense == ks.GAP_STATE, d.Q - 1, dense)
    W = d.reweight(X)
    Pi, Pij = d.compute_frequencies(X, W)
    C = d.build_covariance(Pi, Pij)
    invC = torch.linalg.inv(C + 1e-4 * torch.eye(C.shape[0], dtype=torch.float64))
    L = X.shape[1]
    J = torch.zeros(L, L, d.Q, d.Q, dtype=torch.float64)
    J[:, :, : d.Q - 1, : d.Q - 1] = -invC.reshape(L, d.Q - 1, L, d.Q - 1).permute(0, 2, 1, 3)
    di = d.direct_information(Pi, J)
    Jzs = (J - J.mean(dim=2, keepdim=True) - J.mean(dim=3, keepdim=True)
           + J.mean(dim=(2, 3), keepdim=True))
    F = Jzs[:, :, : d.Q - 1, : d.Q - 1].pow(2).sum(dim=(2, 3)).sqrt().numpy()
    np.fill_diagonal(F, 0.0)
    return di, d.apc_correction(F), float(W.sum())


def stage_scores() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    (ks.RESULTS / "di").mkdir(parents=True, exist_ok=True)
    caches = load_all_caches()
    t0 = time.time()
    for n, name in enumerate(sorted(caches), 1):
        dest = CACHE / f"{name}.npz"
        if dest.exists():
            continue
        tc = caches[name]
        dense = ks.load_protein(name).dense
        i, j = tc.pairs[:, 0], tc.pairs[:, 1]
        old, meff_h = pathh_fn_apc(dense)
        di, fn_zs, meff = corrected_scores(dense)
        meta = {"Meff_pathh": meff_h, "Meff": meff, "max_seq": MAX_SEQ,
                "max_seq_pathh": ph.MAX_SEQ}
        np.savez_compressed(dest, pathh_fn_apc=old[i, j], di=di[i, j],
                            fn_apc_zs=fn_zs[i, j], meta=json.dumps(meta))
        np.savez_compressed(ks.RESULTS / "di" / f"{name}.npz", dca_pairs=old[i, j],
                            meta=json.dumps({"N_used": int(min(dense.shape[0], ph.MAX_SEQ)),
                                             "Meff": meff_h, "L": int(dense.shape[1]),
                                             "source": "analysis/dca_fusion_recheck.py"}))
        print(f"[{n:3d}/150] {name}: L={dense.shape[1]} N={dense.shape[0]} "
              f"Meff={meff:.0f} ({time.time() - t0:.0f}s)", flush=True)


def check_parity(names=("1jo8A", "1brfA", "1guuA")) -> dict:
    """pathh_fn_apc vs path_h_dca_circuit.mfdca_di on a few small targets."""
    caches = load_all_caches()
    out = {}
    for name in names:
        tc = caches[name]
        dense = ks.load_protein(name).dense
        ref, _ = ph.mfdca_di(dense)
        mine, _ = pathh_fn_apc(dense)
        i, j = tc.pairs[:, 0], tc.pairs[:, 1]
        K = max(1, round(tc.n_cols / 5))
        out[name] = {"max_abs_diff": float(np.abs(ref[i, j] - mine[i, j]).max()),
                     "spearman": float(spearmanr(ref[i, j], mine[i, j])[0]),
                     "prec_L5_ref": ks.precision_at_k(ref[i, j], tc.lab, K),
                     "prec_L5_torch": ks.precision_at_k(mine[i, j], tc.lab, K)}
        print(name, out[name], flush=True)
    return out


# ---------------------------------------------------------------------------
# evaluation (path_l's protocol)
# ---------------------------------------------------------------------------

def _cells(tc, dca):
    return (((tc.g_i.astype(np.int64) * 8 + tc.g_j) * 4 + tc.sbin) * 4 + pl.di_bin(dca))


def _boot_ci(x, n=10000, seed=0):
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), size=(n, len(x)))
    m = x[idx].mean(axis=1)
    return [float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))]


def stage_eval() -> dict:
    caches = load_all_caches()
    targets = sorted(caches)
    scores = {t: np.load(CACHE / f"{t}.npz") for t in targets}
    folds = ks.protein_folds(targets, n_folds=5, seed=42)
    rows = []
    for fold in range(5):
        test = folds[fold]
        train = [t for g, grp in enumerate(folds) if g != fold for t in grp]
        counts = {}
        for v in VARIANTS:
            tr = pl.CountsH()
            for name in train:
                tr.add(_cells(caches[name], scores[name][v]), caches[name].lab)
            counts[v] = ks.CellCounts(pl.CountsH.N_CELLS, tr.pos, tr.tot)
        for name in test:
            tc = caches[name]
            K = max(1, round(tc.n_cols / 5))
            row = {"target": name, "fold": fold, "L": tc.n_cols}
            mip = ks.apc_corrected(tc.mi, tc.pairs, tc.n_cols)
            row["mip"] = ks.precision_at_k(mip, tc.lab, K)
            row["mip_auc"] = ks.auc_mann_whitney(mip, tc.lab)
            for v in VARIANTS:
                s = scores[name][v]
                prior = ks.cell_rate_score(counts[v], _cells(tc, s))
                fused = pl.rank01(s) + pl.rank01(prior)
                for key, sc in ((v, s), (f"{v}+prior", prior), (f"{v}+ranksum", fused)):
                    row[key] = ks.precision_at_k(sc, tc.lab, K)
                    row[f"{key}_auc"] = ks.auc_mann_whitney(sc, tc.lab)
            row["di+mip_ranksum"] = ks.precision_at_k(
                pl.rank01(scores[name]["di"]) + pl.rank01(mip), tc.lab, K)
            rows.append(row)

    keys = [k for k in rows[0] if k not in ("target", "fold", "L") and not k.endswith("_auc")]
    summary = {k: {"mean_prec_L5": float(np.mean([r[k] for r in rows])),
                   "mean_auc": (float(np.mean([r[f"{k}_auc"] for r in rows]))
                                if f"{k}_auc" in rows[0] else None)} for k in keys}
    comparisons = [("pathh_fn_apc+ranksum", "pathh_fn_apc"),
                   ("di", "pathh_fn_apc"), ("di", "mip"), ("fn_apc_zs", "mip"),
                   ("di+ranksum", "di"), ("fn_apc_zs+ranksum", "fn_apc_zs"),
                   ("di+prior", "di"), ("di+mip_ranksum", "di"),
                   ("pathh_fn_apc+prior", "mip"), ("mip", "pathh_fn_apc")]
    paired = {}
    for a, b in comparisons:
        x = np.array([r[a] - r[b] for r in rows])
        paired[f"{a} vs {b}"] = {
            "mean_diff": float(x.mean()), "ci95_bootstrap": _boot_ci(x),
            "wilcoxon_p": float(wilcoxon(x).pvalue) if np.any(x != 0) else 1.0,
            "n_better": int((x > 0).sum()), "n_worse": int((x < 0).sum())}
    out = {"protocol": "5-fold protein CV (ks.protein_folds seed 42), prec@L/5, sep >= 6 universe",
           "n_targets": len(rows), "summary": summary, "paired": paired, "per_target": rows}
    OUT.write_text(json.dumps(out, indent=1))
    for k, v in summary.items():
        print(f"  {k:24s} prec@L/5 {v['mean_prec_L5']:.3f}")
    for k, v in paired.items():
        print(f"  {k:40s} {v['mean_diff']:+.3f} CI {v['ci95_bootstrap'][0]:+.3f}..{v['ci95_bootstrap'][1]:+.3f} "
              f"p={v['wilcoxon_p']:.2g} ({v['n_better']}/{v['n_worse']})")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["scores", "eval"], default="eval")
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    if args.check:
        res = check_parity()
        (REPO / "analysis" / "dca_fusion_parity.json").write_text(json.dumps(res, indent=1))
        return 0
    if args.stage == "scores":
        stage_scores()
    else:
        stage_eval()
    return 0


if __name__ == "__main__":
    sys.exit(main())
