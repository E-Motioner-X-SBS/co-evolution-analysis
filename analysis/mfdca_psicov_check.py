#!/usr/bin/env python3
"""mfDCA (dca_mf_analysis.py, corrected coupling sign) on all 150 PSICOV targets.

For every target: encode the PSICOV alignment with the dca_mf_analysis
alphabet (20 He-2012 codes + gap as state 20), subsample to at most MAX_SEQ
sequences (seeded), reweight (theta = 0.2), build the pseudocounted
covariance, invert, and score every candidate pair (sep >= 6, the
kmap_structure universe) by

  di        direct information with W = exp(e_ij), e_ij = -C^{-1}
  di_old    the same with the earlier, sign-flipped W = exp(-e_ij)
  fn_apc    APC-corrected Frobenius norm of e_ij (sign-invariant)

and report precision@L/5 against the native contacts.  This is the check
behind correction A4 in My_Interpretation_of_the_Results/CORRECTION_NOTICE.md.

Output: analysis/mfdca_psicov_check.json
Run:    python analysis/mfdca_psicov_check.py [--max-seq 5000]
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

import dca_mf_analysis as d  # noqa: E402
import kmap_structure as ks  # noqa: E402


def score_target(target: str, max_seq: int) -> dict:
    pdat = ks.load_protein(target)
    dense = pdat.dense.astype(np.int64)
    n_all = dense.shape[0]
    if n_all > max_seq:
        keep = np.sort(np.random.default_rng(42).choice(n_all, max_seq, replace=False))
        dense = dense[keep]
    X = np.where(dense == ks.GAP_STATE, d.Q - 1, dense)
    W = d.reweight(X)
    Pi, Pij = d.compute_frequencies(X, W)
    C = d.build_covariance(Pi, Pij)
    invC = torch.linalg.inv(C + 1e-4 * torch.eye(C.shape[0], dtype=torch.float64))
    L = X.shape[1]
    J4d = torch.zeros(L, L, d.Q, d.Q, dtype=torch.float64)
    J4d[:, :, : d.Q - 1, : d.Q - 1] = -invC.reshape(L, d.Q - 1, L, d.Q - 1).permute(0, 2, 1, 3)
    pairs = ks.pair_universe(pdat)
    lab = np.array([(int(a), int(b)) in pdat.contacts for a, b in pairs])
    K = max(1, round(L / 5))
    i, j = pairs[:, 0], pairs[:, 1]
    out = {"target": target, "L": int(L), "n_seqs": int(n_all), "n_used": int(X.shape[0]),
           "Meff": float(W.sum()), "n_contacts": int(lab.sum())}
    for name, J in (("di", J4d), ("di_old", -J4d)):
        out[f"prec_L5_{name}"] = ks.precision_at_k(d.direct_information(Pi, J)[i, j], lab, K)
    fn = d.apc_correction(d.frobenius_scores(J4d))
    out["prec_L5_fn_apc"] = ks.precision_at_k(fn[i, j], lab, K)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-seq", type=int, default=5000)
    args = ap.parse_args()
    out_path = REPO / "analysis" / "mfdca_psicov_check.json"
    done = {}
    if out_path.exists():
        prev = json.loads(out_path.read_text())
        if prev.get("max_seq") == args.max_seq:
            done = {r["target"]: r for r in prev["per_target"]}
    rows = []
    t0 = time.time()
    for n, t in enumerate(ks.list_targets(), 1):
        if t in done:
            rows.append(done[t])
            continue
        r = score_target(t, args.max_seq)
        rows.append(r)
        print(f"[{n:3d}/150] {t}: L={r['L']} Meff={r['Meff']:.0f} "
              f"DI={r['prec_L5_di']:.3f} DI_old={r['prec_L5_di_old']:.3f} "
              f"FN-APC={r['prec_L5_fn_apc']:.3f} ({time.time() - t0:.0f}s)", flush=True)
        if n % 10 == 0:
            _write(out_path, rows, args.max_seq)
    _write(out_path, rows, args.max_seq)
    s = json.loads(out_path.read_text())["summary"]
    print(json.dumps(s, indent=1))
    return 0


def _write(path: Path, rows: list, max_seq: int) -> None:
    keys = ("prec_L5_di", "prec_L5_di_old", "prec_L5_fn_apc")
    summary = {k: float(np.mean([r[k] for r in rows])) for k in keys}
    summary["n_targets"] = len(rows)
    summary["n_di_better_than_old"] = int(sum(r["prec_L5_di"] > r["prec_L5_di_old"] for r in rows))
    summary["n_di_worse_than_old"] = int(sum(r["prec_L5_di"] < r["prec_L5_di_old"] for r in rows))
    path.write_text(json.dumps({"max_seq": max_seq, "theta": d.THETA, "lambda": d.LAMBDA,
                                "summary": summary, "per_target": rows}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
