#!/usr/bin/env python3
"""
path_j_crossdataset.py — Test PSICOV150-trained circuits on independent datasets.

Trains the universal contact circuit on ALL 150 PSICOV proteins, then evaluates
on completely independent protein families from GPCRdb and unit_tests.
No protein in the test set shares evolutionary history with any training protein.

Also computes plain MI ranking as a per-protein baseline for comparison.

Usage:
    python scripts/path_j_crossdataset.py
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

OUT = ks.RESULTS / "cross_dataset"
OUT.mkdir(parents=True, exist_ok=True)

# ── Target definitions ──
TARGETS = [
    # (dataset, target_name, fasta_path, pdb_path, description)
    ("unit_tests", "DHFR",   "data/unit_tests/converted/DHFR.fasta",
     "data/structures/4p3r.pdb", "E. coli dihydrofolate reductase"),
    ("unit_tests", "1whzA",  "data/unit_tests/converted/1whzA.fasta",
     "data/structures/1whz.pdb", "Protein G B1 domain"),
    ("unit_tests", "1atzA",  "data/unit_tests/converted/1atzA.fasta",
     "data/structures/1atz.pdb", "Azurin"),
    ("gpcrdb", "b2ar",   "data/gpcrdb/converted/b2ar.fasta",
     "data/structures/2rh1.pdb", "β2-adrenergic receptor"),
    ("gpcrdb", "5ht1a",  "data/gpcrdb/converted/5ht1a.fasta",
     "data/structures/7e2x.pdb", "5-HT1A serotonin receptor"),
]


def mi_pair(dense, i, j):
    ci = dense[:, i]; cj = dense[:, j]
    v = (ci >= 0) & (ci < 20) & (cj >= 0) & (cj < 20)
    cv, cjv = ci[v].astype(np.int64), cj[v].astype(np.int64)
    if len(cv) < 10:
        return 0.0
    jnt = np.bincount(cv * 20 + cjv, minlength=400).reshape(20, 20).astype(float)
    t = jnt.sum()
    if t == 0:
        return 0.0
    pj = jnt / t; pim = jnt.sum(1) / t; pjm = jnt.sum(0) / t
    val = sum(pj[a, b] * np.log2(pj[a, b] / (pim[a] * pjm[b]))
              for a in range(20) for b in range(20)
              if jnt[a, b] > 0 and pim[a] > 0 and pjm[b] > 0)
    return max(val, 0.0)


def rank01(s):
    r = np.argsort(np.argsort(s, kind="stable"))
    return r / max(len(s) - 1, 1)


def _load_ragged(path):
    """Load FASTA padding ragged rows to max width."""
    seqs = []
    for ln in open(path):
        ln = ln.strip()
        if not ln or ln.startswith(">"):
            continue
        seqs.append(ln.upper())
    width = max(len(s) for s in seqs)
    lut = np.full(256, ks.GAP_STATE, dtype=np.int8)
    for ch, code in ks.AA_TO_CODE.items():
        lut[ord(ch)] = code; lut[ord(ch.lower())] = code
    buf = np.frombuffer("".join(s.ljust(width, '-')[:width] for s in seqs).encode(), dtype=np.uint8)
    return lut[buf].reshape(len(seqs), width)


def main():
    from run_contact_campaign import load_all_caches

    print("=" * 70)
    print("PATH J: Cross-Dataset Validation of K-Map Contact Circuits")
    print("Training: PSICOV150 (all 150 proteins)")
    print("Testing: GPCRdb + Unit tests (completely independent families)")
    print("=" * 70)

    # ── Step 1: Train universal circuit on ALL PSICOV150 ──
    print("\n[1] Loading PSICOV150 caches...")
    caches = load_all_caches()
    psicov_targets = sorted(caches)

    tr_pos = np.zeros(1024, np.int64)
    tr_tot = np.zeros(1024, np.int64)

    for name in psicov_targets:
        tc = caches[name]
        idx = tc.level1_cells()
        np.add.at(tr_tot, idx, 1)
        np.add.at(tr_pos, idx, tc.lab.astype(np.int64))

    prior_counts = ks.CellCounts(1024, tr_pos, tr_tot)
    base_psicov = float(tr_pos.sum() / max(tr_tot.sum(), 1))
    n_on = int((ks.label_truth_table(
        ks.CellCounts(1024, tr_pos.copy(), tr_tot.copy()), 2.0, 30) == 1).sum())
    print(f"  Trained on {len(psicov_targets)} proteins")
    print(f"  Base rate: {base_psicov:.4f} | ON cells: {n_on}/256")

    # Group map: He code → group index
    group_map = {}
    for aa, code in ks.AA_TO_CODE.items():
        group_map[code] = ks.AA_GROUP8.get(aa, 7)

    # ── Step 2: Evaluate each target ──
    all_results = []

    for ds_name, target_name, fasta_rel, pdb_rel, desc in TARGETS:
        print(f"\n{'='*60}")
        print(f"[{ds_name}] {target_name}: {desc}")

        fasta_path = ks.REPO / fasta_rel
        pdb_path = ks.REPO / pdb_rel
        if not fasta_path.exists() or not pdb_path.exists():
            print(f"  SKIP: missing files"); continue

        # Load MSA
        dense = _load_ragged(fasta_path)
        N, L = dense.shape
        cons = ks.consensus_codes(dense)
        valid = cons >= 0
        vc = np.flatnonzero(valid)

        # Candidate pairs
        ii, jj = np.triu_indices(vc.size, k=1)
        pi_arr, pj_arr = vc[ii], vc[jj]
        sep_arr = pj_arr - pi_arr
        kp = sep_arr >= 6
        pi_arr, pj_arr, sep_arr = pi_arr[kp], pj_arr[kp], sep_arr[kp]

        if len(pi_arr) < 50 or L < 30:
            print(f"  SKIP: only {len(pi_arr)} pairs, L={L}")
            continue

        # Native contacts
        contacts = ks.native_contacts(pdb_path, n_cols=L)
        lab = np.array([(int(a), int(b)) in contacts
                        for a, b in zip(pi_arr, pj_arr)])
        n_contacts = int(lab.sum())
        base_rate = float(lab.mean())

        if n_contacts == 0:
            print(f"  SKIP: zero native contacts"); continue

        # Features
        g_i = np.array([group_map.get(int(c), 7) for c in cons[pi_arr]], dtype=np.int8)
        g_j = np.array([group_map.get(int(c), 7) for c in cons[pj_arr]], dtype=np.int8)
        s_bins = np.array([min((s - 6) // 10, 3) for s in sep_arr], dtype=np.int8)

        # MI scores (per-protein baseline)
        print(f"  Computing MI ({len(pi_arr)} pairs, N={N})...")
        t0 = time.time()
        mi_vals = np.array([mi_pair(dense, int(i), int(j))
                            for i, j in zip(pi_arr, pj_arr)])
        mi_time = time.time() - t0

        # Circuit scores (using PSICOV-trained prior)
        idx_l1 = level1_index(g_i, g_j, s_bins)
        circ_scores = ks.cell_rate_score(prior_counts, idx_l1)

        # Rank fusion
        fused = rank01(mi_vals) + rank01(circ_scores)

        # Evaluate
        K = max(1, round(L / 5))
        results_row = {
            "dataset": ds_name, "target": target_name,
            "description": desc,
            "L": L, "N": N, "n_pairs": len(pi_arr),
            "n_contacts": n_contacts,
            "base_rate": round(base_rate, 4),
            "K": K,
        }

        for method, scores in [("MI", mi_vals), ("circH", circ_scores),
                               ("ranksum", fused)]:
            p5 = ks.precision_at_k(scores, lab, K)
            auc_v = ks.auc_mann_whitney(scores, lab)
            enr = p5 / base_rate if base_rate > 0 else 0
            results_row[method] = {"prec_L5": round(p5, 4),
                                   "auc": round(auc_v, 3),
                                   "enrichment": round(enr, 2)}

        all_results.append(results_row)

        print(f"  L={L} pairs={len(pi_arr)} contacts={n_contacts} "
              f"base={base_rate:.4f} K={K}")
        print(f"  MI:      prec@L/5={results_row['MI']['prec_L5']:.3f} "
              f"AUC={results_row['MI']['auc']:.3f} "
              f"({results_row['MI']['enrichment']:.1f}x)")
        print(f"  circH:   prec@L/5={results_row['circH']['prec_L5']:.3f} "
              f"AUC={results_row['circH']['auc']:.3f} "
              f"({results_row['circH']['enrichment']:.1f}x)")
        print(f"  ranksum: prec@L/5={results_row['ranksum']['prec_L5']:.3f} "
              f"AUC={results_row['ranksum']['auc']:.3f} "
              f"({results_row['ranksum']['enrichment']:.1f}x)")

    # ── Summary ──
    print("\n" + "=" * 70)
    print("CROSS-DATASET SUMMARY")
    print(f"\n{'Dataset':<12} {'Target':<10} {'L':>4} {'N':>6} "
          f"{'MI':>6} {'circH':>7} {'ranksum':>8} "
          f"{'MI_enr':>7} {'fus_enr':>7}")

    for r in all_results:
        mi_p = r.get("MI", {}).get("prec_L5", 0)
        ch_p = r.get("circH", {}).get("prec_L5", 0)
        rs_p = r.get("ranksum", {}).get("prec_L5", 0)
        mi_e = r.get("MI", {}).get("enrichment", 0)
        rs_e = r.get("ranksum", {}).get("enrichment", 0)
        print(f"{r['dataset']:<12} {r['target']:<10} {r['L']:>4} {r['N']:>6} "
              f"{mi_p:>6.3f} {ch_p:>7.3f} {rs_p:>8.3f} "
              f"{mi_e:>7.2f} {rs_e:>7.2f}")

    means = {}
    for m in ("MI", "circH", "ranksum"):
        vals = [r[m]["prec_L5"] for r in all_results if m in r]
        means[m] = float(np.mean(vals)) if vals else 0
    print(f"\nMean prec@L/5: MI={means.get('MI',0):.3f} "
          f"circH={means.get('circH',0):.3f} ranksum={means.get('ranksum',0):.3f}")

    (OUT / "cross_dataset_results.json").write_text(json.dumps(all_results, indent=1))
    print(f"\nsaved -> {OUT/'cross_dataset_results.json'}")


def level1_index(g_i, g_j, s_bin):
    return ((g_i.astype(np.int64) * 8 + g_j) * 4 + s_bin)


if __name__ == "__main__":
    main()
