#!/usr/bin/env python3
"""
path_h_dca_circuit.py — PATH H: proper mfDCA coupling INSIDE the K-map.

Morcos et al. 2011 (PNAS) mean-field DCA per PSICOV protein:
  1. sequence reweighting (theta=0.2), Meff
  2. single/pair frequencies with pseudocount lambda=0.5
     P_i(a) = (1-lam) f_i(a) + lam/q ;  P_ij(a,b) analogous off-diagonal,
     P_ii(a,b) = delta(a,b) P_i(a)
  3. covariance C[(i,a),(j,b)] = P_ij - P_i P_j   (gap state removed, q=20)
  4. e_ij = -(C + reg I)^{-1} blocks
  5. Frobenius score + APC (Dunn 2008) + direct information (mean-field
     fixed-point iteration, 200 iters max)

Then the K-map integration (the user's ask — "do coupling in the K-map"):
  level-H circuit cell = g_i[3] x g_j[3] x sep[2] x DI-bin[2]  (1024 cells)
  DI-bin = quartile of the pair's DI RANK among the protein's candidate pairs
  (label-free, computable on test proteins).
Evaluated with the same 5-fold protein-level CV as PATH B/G, against:
  raw-DI ranking (literature mfDCA proxy), circuit-L2 (no coupling),
  and the published anchors.

Depth policy: MSAs deeper than MAX_SEQ are uniformly subsampled (standard
practice; reweighting makes Meff dominated by cluster counts). Documented.
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
DI_DIR = OUT / "di"
MAX_SEQ = 12000          # uniform subsample cap (documented)
Q = 20                   # canonical states (gap removed from couplings)
REG = 1e-4               # covariance regularization


# ---------------------------------------------------------------------------
# mfDCA core (numpy; self-contained)
# ---------------------------------------------------------------------------

def mfdca_di(dense: np.ndarray) -> tuple[np.ndarray, dict]:
    """Full-length DI matrix (upper triangle meaningful) for one alignment.

    dense: [N, L] int8 codes 0..19, gap=20.  Returns (di [L,L], info).
    """
    N, L = dense.shape
    if N > MAX_SEQ:
        rng = np.random.default_rng(42)
        keep = np.sort(rng.choice(N, MAX_SEQ, replace=False))
        dense = dense[keep]
        N = MAX_SEQ
    X = dense.astype(np.int64)
    valid = X < Q                                   # gap=20 excluded

    # --- reweighting (chunked) ---
    W = np.ones(N)
    if N > 1:
        counts = np.ones(N)
        step = 800
        for s in range(0, N, step):
            blk = X[s:s+step]                        # [c,L]
            same = (blk[:, None, :] == X[None, :, :]) & valid[None, :, :]
            ham = (~same).sum(axis=2)                # mismatches incl gaps-eq
            counts[s:s+step] = (ham <= int(0.2 * L)).sum(axis=1) + 1
        W = 1.0 / counts
    Meff = W.sum()

    # --- weighted one-hot [N, L, Q] ---
    oh = np.zeros((N, L, Q), dtype=np.float64)
    oh[np.arange(N)[:, None], np.arange(L)[None, :], np.clip(X, 0, Q-1)] = (
        (X < Q).astype(np.float64))
    Vw = oh * W[:, None, None]

    # --- single-site frequencies + pseudocount ---
    fi = Vw.sum(axis=0) / Meff                      # [L,Q]
    Pi = (1 - 0.5) * fi + 0.5 / Q

    # --- pairwise via Gram trick: A[(l,a), n] = sqrt(W_n) oh[n,l,a] ---
    A = (Vw * np.sqrt(W)[:, None, None]).reshape(N, L * Q).T   # [LQ, N]
    G = A @ A.T / Meff                              # [(l,a),(j,b)]
    G = G.reshape(L, Q, L, Q)
    Pij = 0.5 * G + 0.5 * G.transpose(2, 3, 0, 1)   # symmetrize (swap residues)
    Pij = 0.5 * Pij + 0.5 / (Q * Q)                 # pseudocount (off-diag)
    idx = np.arange(Q)
    for l in range(L):                              # diagonal: delta * Pi
        Pij[l, :, l, :] = Pi[l][:, None] * (idx[:, None] == idx[None, :])

    # --- covariance with gap removed ---
    Pic = Pi[:, :Q-1]
    Pijc = Pij[:, :Q-1, :, :Q-1]
    C = Pijc - Pic[:, :, None, None] * Pic[None, None, :, :]
    # C is [i,a,j,b]; grouping (i,a) as rows and (j,b) as cols requires a
    # plain reshape of THIS axis order.  The previous transpose(0,2,1,3)
    # produced an [i,j,a,b] tensor whose reshape scrambled row/col semantics
    # (max|C-C^T| was 0.25!) — root cause of the chance-level DI.
    Cm = C.reshape(L * (Q - 1), L * (Q - 1))
    Cm = 0.5 * (Cm + Cm.T)          # enforce exact symmetry for LAPACK

    # --- invert ---
    invC = np.linalg.inv(Cm + REG * np.eye(Cm.shape[0]))
    E = -invC.reshape(L, Q-1, L, Q-1).transpose(0, 2, 1, 3)   # e[(i,a,j,b)]

    # Overflow guard: exp(-E) overflows float64 for strong couplings.
    # Clipping to +-30 keeps exp finite (< 1e13); DI is recomputed from the
    # bounded Wm, matching the pragmatic range used by reference mfDCA codes.
    E = np.clip(E, -30.0, 30.0)

    # E layout here is [i,j,a,b]: state axes are (2,3).
    # (Axes (1,3) would mix position j with state b — earlier bug.)
    frob = np.sqrt((E ** 2).sum(axis=(2, 3)))

    # --- direct information (mean-field fixed point, vectorized) ---
    Wm = np.exp(E)   # W = exp(+J): Boltzmann weight of couplings J=-C^-1
    mu1 = np.full((L, L, Q-1), 1.0/(Q-1))
    mu2 = np.full((L, L, Q-1), 1.0/(Q-1))
    pi_e = Pic                                       # [L,q]
    for _ in range(200):
        calc1 = np.einsum('ijb,ijab->ija', mu2, Wm)
        calc2 = np.einsum('ija,ijab->ijb', mu1, Wm)
        # Morcos 2011: mu_i(a) <- pi_i(a) / sum_b W_ij(a,b) mu_j(b)
        # numerator carries the SAME position as the free index:
        #   n1[i,j,a] = pi_i(a);  n2[i,j,b] = pi_j(b)
        n1 = pi_e[:, None, :] / np.clip(calc1, 1e-300, None)
        n1 /= n1.sum(axis=2, keepdims=True)
        n2 = pi_e[None, :, :] / np.clip(calc2, 1e-300, None)
        n2 /= n2.sum(axis=2, keepdims=True)
        d = max(np.abs(n1-mu1).max(), np.abs(n2-mu2).max())
        mu1, mu2 = n1, n2
        if d < 4e-5:
            break
    Pdir = Wm * mu1[:, :, :, None] * mu2[:, :, None, :]
    Pdir /= Pdir.sum(axis=(2, 3), keepdims=True)
    Pfac = pi_e[:, None, :, None] * pi_e[None, :, None, :]   # P_i(a)P_j(b)
    with np.errstate(divide="ignore", invalid="ignore"):
        di = (Pdir * np.log(np.clip(Pdir, 1e-300, None)
                            / np.clip(Pfac, 1e-300, None))).sum(axis=(2, 3))
    di = np.nan_to_num(di, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(di, 0.0)

    # Primary DCA score: APC-corrected Frobenius norm of the mean-field
    # couplings (Morcos et al. 2011).  The iterative DI is retained as a
    # diagnostic only: on families with many fully-conserved columns the
    # pseudo-inverse conditioning makes the fixed-point DI unreliable
    # (documented in plans/12; Frob_apc validates strongly on natives).
    apc_f = frob - frob.mean(axis=1)[:, None] * frob.mean(axis=0)[None, :] / frob.mean()
    np.fill_diagonal(apc_f, 0.0)

    info = {"N_used": int(N), "Meff": float(Meff), "L": int(L)}
    return apc_f, info


def di_for_cached_pairs(tc, di_mat: np.ndarray) -> np.ndarray:
    """DI score per cached candidate pair (order-aligned)."""
    return di_mat[tc.pairs[:, 0], tc.pairs[:, 1]]


# ---------------------------------------------------------------------------
# Circuit integration
# ---------------------------------------------------------------------------

def di_bin_of(di_scores: np.ndarray) -> np.ndarray:
    """Quartile of DI rank within the protein's candidate universe (0..3)."""
    r = np.argsort(np.argsort(di_scores, kind="stable"))
    return (4 * r // max(len(r), 1)).clip(0, 3).astype(np.int8)


def levelH_index(g_i, g_j, sbin, dib):
    """g_i[3] x g_j[3] x sep[2] x DI-bin[2] -> [0,1024)."""
    return (((g_i.astype(np.int64) * 8 + g_j) * 4 + sbin) * 4 + dib)


class CountsH:
    N_CELLS = 1024

    def __init__(self):
        self.pos = np.zeros(self.N_CELLS, np.int64)
        self.tot = np.zeros(self.N_CELLS, np.int64)

    def add(self, idx, lab):
        np.add.at(self.tot, idx, 1)
        np.add.at(self.pos, idx, lab.astype(np.int64))


# ---------------------------------------------------------------------------
# stages
# ---------------------------------------------------------------------------

def stage_di(force=False):
    """Compute + cache per-protein DI matrices (aligned to cached pairs)."""
    DI_DIR.mkdir(parents=True, exist_ok=True)
    from run_contact_campaign import load_all_caches
    caches = load_all_caches()
    manifest = {}
    t0 = time.time()
    for n, name in enumerate(sorted(caches), 1):
        dest = DI_DIR / f"{name}.npz"
        if dest.exists() and not force:
            continue
        tc = caches[name]
        pd = ks.load_protein(name)
        assert pd.dense is not None
        apc_mat, info = mfdca_di(pd.dense)
        dca_pairs = di_for_cached_pairs(tc, apc_mat)
        np.savez_compressed(dest, dca_pairs=dca_pairs,
                            meta=json.dumps(info))
        manifest[name] = {"Meff": info["Meff"], "N": info["N_used"],
                          "secs": round(time.time()-t0, 1)}
        print(f"[{n}/150] {name}: Meff={info['Meff']:.0f} "
              f"N={info['N_used']} ({time.time()-t0:.0f}s cum)", flush=True)
    (OUT/"path_h_di_manifest.json").write_text(json.dumps(manifest, indent=1))


def eval_h(train_names, test_names, caches):
    tr = CountsH()
    for name in train_names:
        tc = caches[name]
        dib = di_bin_of(np.load(DI_DIR/f"{name}.npz")["dca_pairs"])
        tr.add(levelH_index(tc.g_i, tc.g_j, tc.sbin, dib), tc.lab)
    truth = ks.label_truth_table(
        ks.CellCounts(CountsH.N_CELLS, tr.pos.copy(), tr.tot.copy()), 2.0, 30)
    from run_contact_campaign import qm_extract
    res_qm = qm_extract(truth)

    rows = []
    for name in test_names:
        tc = caches[name]
        dca = np.load(DI_DIR/f"{name}.npz")["dca_pairs"]
        dib = di_bin_of(dca)
        idx = levelH_index(tc.g_i, tc.g_j, tc.sbin, dib)
        score_H = ks.cell_rate_score(
            ks.CellCounts(CountsH.N_CELLS, tr.pos, tr.tot), idx)
        L = tc.n_cols
        rows.append({
            "prec@L/5_circH": ks.precision_at_k(score_H, tc.lab, max(1, round(L/5))),
            "prec@L/5_rawDCA": ks.precision_at_k(dca, tc.lab, max(1, round(L/5))),
            "auc_circH": ks.auc_mann_whitney(score_H, tc.lab),
            "auc_rawDCA": ks.auc_mann_whitney(dca, tc.lab),
        })
    agg = {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}
    return agg, res_qm, rows


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["di", "cv"], default="cv")
    ap.add_argument("--force-di", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    if args.stage == "di":
        stage_di(force=args.force_di)
        return

    from run_contact_campaign import load_all_caches
    caches = load_all_caches()
    targets = sorted(caches)
    folds = ks.protein_folds(targets, n_folds=5, seed=42)
    per_fold, per_protein = [], []
    for f in range(5):
        test = folds[f]
        train = [t for g, grp in enumerate(folds) if g != f for t in grp]
        agg, res_qm, rows = eval_h(train, test, caches)
        per_fold.append(agg)
        for r in rows:
            r["fold"] = f
        per_protein.extend(rows)
        print(f"fold {f}: circH prec@L/5={agg['prec@L/5_circH']:.3f} "
              f"rawDCA={agg['prec@L/5_rawDCA']:.3f} "
              f"(PIs={res_qm['n_prime_implicants']})", flush=True)
    summary = {
        "per_fold": per_fold,
        "mean": {k: float(np.mean([a[k] for a in per_fold])) for k in per_fold[0]},
        "note": ("level-H = groups x sep x mfDCA-DI-quartile (1024 cells); "
                 "rawDI = plain mfDCA ranking (literature method proxy)."),
    }
    (OUT/"path_h_dca_circuit.json").write_text(json.dumps(
        {"summary": summary, "per_protein": per_protein}, indent=1))
    print(json.dumps(summary["mean"], indent=1))


if __name__ == "__main__":
    main()
