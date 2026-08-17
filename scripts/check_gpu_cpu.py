#!/usr/bin/env python3
"""
PHASE 5: GPU/CPU consistency check per dataset (sample targets).

Compares GPU kernels (coevolution_gpu) vs CPU implementations (coevolution_shared)
on: entropy per position, majority reference residues, and MI for all variable
position pairs (window 30). PASS: max|diff| < 1e-6 (entropy/MI), 0 ref mismatches.

Usage: python scripts/check_gpu_cpu.py [--target NAME] [--all]
Writes results/consistency/consistency_report.json
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

REPO = Path("/store/shuvam/E-motioner-X-SBS/co-evolution-analysis")
sys.path.insert(0, str(REPO))
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/kmap-sbm-validation/src")
sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")

import numpy as np  # noqa: E402

SAMPLES = {
    "unit_tests": [("DHFR", "data/unit_tests/converted/DHFR.fasta")],
    "gpcrdb": [("b2ar", "data/gpcrdb/converted/b2ar.fasta")],
    "evmutation": [
        (
            "DYR_ECOLI_hmmerbit_plmc_n5_m30_f50_t0.2_r1-159_id100_b80",
            "data/evmutation/converted/DYR_ECOLI_hmmerbit_plmc_n5_m30_f50_t0.2_r1-159_id100_b80.fasta",
        )
    ],
    "pfam": [("PF00071", "data/pfam/converted/PF00071.fasta")],
    "psicov150": [("1ctfA", "data/psicov150/converted/1ctfA.fasta")],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", help="check one (name)")
    ap.add_argument("--all", action="store_true")
    ap.add_argument(
        "--max-seq",
        type=int,
        default=0,
        help="cap sequences (large Pfam families); 0 = all",
    )
    args = ap.parse_args()

    import coevolution_shared as cs
    import coevolution_gpu as cg

    if not cg.check_cuda():
        print("CUDA NOT AVAILABLE — consistency check impossible")
        return 1

    out_dir = REPO / "results" / "consistency"
    out_dir.mkdir(parents=True, exist_ok=True)
    report = {"timestamp": time.strftime("%Y-%m-%d %H:%M:%S"), "checks": []}

    for ds, items in SAMPLES.items():
        for name, fasta_rel in items:
            if args.target and name != args.target:
                continue
            fasta = REPO / fasta_rel
            if not fasta.exists():
                continue
            print(f"=== {ds}/{name} ===", flush=True)
            pos_arrays, n_all, full_len = cs.load_position_arrays(
                fasta_path=str(fasta), max_pos=None, aligned=True, clear_cache=True
            )
            if args.max_seq:
                pos_arrays = pos_arrays[: args.max_seq]
                n_all = args.max_seq
            L = len(pos_arrays[0])
            print(f"  {n_all} seqs x {L} positions", flush=True)

            # CPU
            ent_cpu = cs.compute_entropy_vectorized(pos_arrays, n_all, L)
            refs_cpu = np.array(
                [cs.majority_ref(pos_arrays, p, n_all) for p in range(L)],
                dtype=np.int32,
            )
            var = [p for p in range(L) if ent_cpu[p] > 0.3]
            pairs = [
                (i, j) for a, i in enumerate(var) for j in var[a + 1 :] if j <= i + 30
            ]
            mi_cpu = {}
            for i, j in pairs:
                m = cs.mutual_information(pos_arrays, i, j, n_all)
                if m > 0:
                    mi_cpu[(i, j)] = m

            # GPU
            dense = cg.dense_to_gpu(pos_arrays)
            ent_gpu = cg.compute_entropy_gpu(dense).cpu().numpy()
            refs_gpu = cg.majority_refs_gpu(dense).cpu().numpy()
            refs_g = refs_gpu.astype(np.int32)
            mi_gpu = {}
            if pairs:
                import torch

                refs_t = torch.tensor(refs_g, device=dense.device)
                md, cd = cg.mi_matrix_gpu(dense, pairs, refs=refs_t, chunk=4096)
                for i, j in pairs:
                    if (i, j) in md and md[(i, j)] > 0:
                        mi_gpu[(i, j)] = float(md[(i, j)])

            d_ent = float(np.max(np.abs(ent_cpu - ent_gpu)))
            n_ref_mm = int((refs_cpu != refs_g).sum())
            common = set(mi_cpu) & set(mi_gpu)
            d_mi = float(max((abs(mi_cpu[k] - mi_gpu[k]) for k in common), default=0.0))
            rank_cpu = sorted(list(mi_cpu), key=lambda k: mi_cpu[k], reverse=True)[:10]
            rank_gpu = sorted(list(mi_gpu), key=lambda k: mi_gpu[k], reverse=True)[:10]
            agree = len(set(rank_cpu) & set(rank_gpu))

            row = {
                "dataset": ds,
                "target": name,
                "n_seqs": n_all,
                "L": L,
                "pairs": len(pairs),
                "max_entropy_diff": d_ent,
                "ref_mismatches": n_ref_mm,
                "max_mi_diff": d_mi,
                "top10_rank_agree": agree,
                "pass": d_ent < 1e-6 and n_ref_mm == 0 and d_mi < 1e-6,
            }
            report["checks"].append(row)
            print(
                f"  entropy diff={d_ent:.2e} ref_mm={n_ref_mm} "
                f"mi_diff={d_mi:.2e} top10_agree={agree}/10 -> "
                f"{'PASS' if row['pass'] else 'FAIL'}",
                flush=True,
            )

    json.dump(report, open(out_dir / "consistency_report.json", "w"), indent=1)
    n_pass = sum(1 for c in report["checks"] if c["pass"])
    print(
        f"consistency: {n_pass}/{len(report['checks'])} PASS -> "
        f"{out_dir / 'consistency_report.json'}"
    )
    return 0 if n_pass == len(report["checks"]) else 1


if __name__ == "__main__":
    sys.exit(main())
