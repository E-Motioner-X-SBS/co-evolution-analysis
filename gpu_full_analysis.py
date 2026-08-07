#!/usr/bin/env python3
"""GPU-Accelerated Full Co-evolution Analysis — ALL 1,299 Sequences, FULL Length.

Uses torch CUDA (A100 80GB) for ALL computations:
  - Entropy (vectorized one-hot on GPU)
  - H1 Gray adjacency (tensor popcount on GPU)
  - Mutations (GPU comparison vs reference)
  - MI matrix: ALL position pairs (full 1276×1276, ~813K pairs) on GPU
  - Coupling constants on GPU

No position truncation, no sequence subsampling.
"""

import sys, json, time
import numpy as np
from pathlib import Path

sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")
from nkmap.encoding.bio_sequences import Base20AminoEncoder, AMINO_HE_2012

sys.path.insert(0, str(Path(__file__).resolve().parent))
import coevolution_gpu as cg

BASE = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
FASTA = BASE / "Spike_protein.aln-fasta"
OUT = BASE / "full_gpu_results"
OUT.mkdir(exist_ok=True)
AA = list(AMINO_HE_2012)

import torch

print(
    f"torch {torch.__version__}, CUDA: {torch.cuda.is_available()}, "
    f"device: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU'}"
)


def parse_fasta(fp):
    seqs, c_h, c_s = [], None, []
    with open(fp) as f:
        for l in f:
            l = l.strip()
            if l.startswith(">"):
                if c_h:
                    seqs.append((c_h, "".join(c_s)))
                c_h = l[1:]
                c_s = []
            elif l:
                c_s.append(l.upper())
    if c_h:
        seqs.append((c_h, "".join(c_s)))
    return seqs


def main():
    t0 = time.time()
    print("=" * 80)
    print("GPU CO-EVOLUTION ANALYSIS — ALL 1,299 SEQUENCES, FULL LENGTH")
    print("=" * 80)

    print("\n[1/3] Loading & encoding ALL sequences...")
    seqs = parse_fasta(FASTA)
    enc = Base20AminoEncoder(version=1)
    n_all = len(seqs)
    fl = len(seqs[0][1])
    print(f"  {n_all} sequences, {fl} positions (FULL LENGTH)")

    # Build position arrays (He 2012 codes 0-19, 20 = gap, ALIGNED)
    # CORRECTED (FIX A1): was gap-stripped (-1), which misaligned columns.
    # Now keeps the alignment: '-' and unknown -> 20.
    pos_arrays = []
    for _, s in seqs:
        arr = np.array([enc.encode.get(c, 20) for c in s], dtype=np.int32)
        pos_arrays.append(arr)

    # Move to GPU once
    dense = cg.dense_to_gpu(pos_arrays)
    print(
        f"  Dense tensor: {tuple(dense.shape)}, "
        f"{dense.nbytes / 1e9:.2f} GB on {dense.device}"
    )

    print("\n[2/3] Computing ALL metrics on GPU...")

    # ── Entropy (GPU) ──────────────────────────────────────────────
    ent = cg.compute_entropy_gpu(dense).cpu().numpy()
    pp = 2.0**ent
    vp = np.where(ent > 0.3)[0].tolist()
    se = np.argsort(ent)[::-1]
    print(f"  Variable positions: {len(vp)}/{fl} ({100 * len(vp) / fl:.1f}%)")

    # ── H1 Gray adjacency (GPU, He-2012 direct Gray) ───────────────
    h1_ratio, h1_count, h1_total = cg.h1_adjacency_gpu(dense)
    er = 0.1613
    enr = h1_ratio / er if er > 0 else 0
    print(f"  H1: {h1_count}/{h1_total} = {h1_ratio:.4f} (enrichment {enr:.2f}x)")

    # ── Mutations vs majority reference (GPU) ──────────────────────
    refs = cg.majority_refs_gpu(dense)  # [max_pos] majority refs
    ref_row = refs[None, :].expand_as(dense)  # [n, max_pos]
    valid = dense >= 0
    is_mut = valid & (dense != ref_row)
    mut_counts = is_mut.sum(dim=0).cpu().numpy()  # mutations per position
    tps = is_mut.sum(dim=1).cpu().numpy()  # mutations per sequence
    mc = mut_counts
    tv = np.argsort(mc)[::-1][:20]
    print(f"  Mutations: mean={tps.mean():.1f}/seq")

    # ── Full MI matrix: ALL pairs (FULL length, ~813K pairs) on GPU ─
    # BUG FIX: dense width (1275 after gap filtering) may be 1 less than raw
    # sequence length (1276). Use dense.shape[1] to avoid out-of-bounds GPU index.
    dense_fl = dense.shape[1]
    print("\n  Computing FULL MI matrix on GPU (all pairs, full length)...")
    all_pairs = cg.all_pairs(dense_fl)  # full matrix, ~813K pairs
    print(f"  Pairs: {len(all_pairs)}")
    mi_dict, cnt_dict = cg.mi_matrix_gpu(dense, all_pairs, min_total=5, chunk=32768)

    # Build MI matrix + sorted pair list
    mi_mat = np.zeros((dense_fl, dense_fl), dtype=np.float64)
    for (i, j), mi in mi_dict.items():
        mi_mat[i, j] = mi
        mi_mat[j, i] = mi
    mr = sorted(
        [(i, j, mi, cnt_dict[(i, j)]) for (i, j), mi in mi_dict.items() if mi > 0.01],
        key=lambda x: x[2],
        reverse=True,
    )
    hm = sum(1 for r in mr if r[2] > 0.5)
    print(f"  MI pairs (MI>0.01): {len(mr)} | MI>0.5: {hm}")
    print(f"  Max MI: {mr[0][2]:.4f} at ({mr[0][0]},{mr[0][1]})")

    # ── Couplings for top 10 pairs (GPU) ───────────────────────────
    print("  Computing couplings on GPU...")
    cr = []
    for pi, pj, mi, n in mr[:10]:
        J = cg.coupling_matrix_gpu(dense, pi, pj)
        avg_abs = float(np.mean(np.abs(J)))
        top, anti = [], []
        km = np.zeros((20, 20), dtype=np.float64)
        for a in pos_arrays:
            if pi < len(a) and pj < len(a):
                ci, cj = int(a[pi]), int(a[pj])
                if 0 <= ci < 20 and 0 <= cj < 20:
                    km[ci, cj] += 1
        t = km.sum()
        if t > 0:
            km /= t
            for ai in range(20):
                for aj in range(20):
                    if km[ai, aj] > 0.001 and J[ai, aj] > 0:
                        top.append(
                            (AA[ai], AA[aj], float(J[ai, aj]), float(km[ai, aj]))
                        )
                    elif km[ai, aj] > 0.001 and J[ai, aj] < 0:
                        anti.append(
                            (AA[ai], AA[aj], float(J[ai, aj]), float(km[ai, aj]))
                        )
        top.sort(key=lambda x: x[2], reverse=True)
        anti.sort(key=lambda x: x[2])
        ri, rj = int(refs[pi]), int(refs[pj])
        cr.append(
            {
                "pos_i": pi,
                "pos_j": pj,
                "mi": mi,
                "ref_i": AA[ri],
                "ref_j": AA[rj],
                "avg_coupling": avg_abs,
                "top_co": top[:5],
                "top_anti": anti[:5],
            }
        )
        print(f"    ({pi},{pj}): MI={mi:.4f}, avg|J|={avg_abs:.4f}")

    # Save
    print("\n[3/3] Saving results...")
    summary = {
        "dataset": f"{n_all} Omicron Spike sequences",
        "full_length": fl,
        "compute_time_s": round(time.time() - t0, 1),
        "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "CPU",
        "entropy": {
            "n_variable": len(vp),
            "n_conserved": fl - len(vp),
            "top_20": [
                {"pos": int(i), "ent": float(ent[i]), "pp": float(pp[i])}
                for i in se[:20]
            ],
        },
        "h1": {
            "total_pairs": int(h1_total),
            "h1_count": int(h1_count),
            "observed": float(h1_ratio),
            "expected": float(er),
            "enrichment": float(enr),
        },
        "mutations": {
            "n": int(n_all - 1),
            "mean": float(tps.mean()),
            "max": int(tps.max()),
            "min": int(tps.min()),
            "top_pos": [{"pos": int(p), "n": int(mc[p])} for p in tv],
        },
        "mi": {
            "n_pairs": len(mr),
            "n_high_mi_05": hm,
            "max_mi": float(mr[0][2]) if mr else 0,
            "max_pair": [int(mr[0][0]), int(mr[0][1])] if mr else [],
            "top_30": [{"pi": p, "pj": j, "mi": m} for p, j, m, _ in mr[:30]],
        },
        "couplings": cr,
    }
    with open(OUT / "gpu_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    np.save(OUT / "mi_matrix_full.npy", mi_mat)
    np.savetxt(OUT / "mi_matrix_full.csv", mi_mat, delimiter=",")
    with open(OUT / "entropies.csv", "w") as f:
        f.write("pos,entropy,perplexity\n")
        for i in range(len(ent)):
            f.write(f"{i},{ent[i]:.6f},{pp[i]:.4f}\n")

    print(f"\nSAVED to {OUT}")
    print(f"Total time: {time.time() - t0:.1f}s")
    print(f"\nKEY RESULTS:")
    print(f"  Variable positions: {len(vp)}/{fl} ({100 * len(vp) / fl:.1f}%)")
    print(f"  H1 enrichment: {enr:.2f}x")
    print(f"  Max MI (full length): {mr[0][2]:.4f} at ({mr[0][0]},{mr[0][1]})")
    print(f"  Mean mutations/seq: {tps.mean():.1f}")




    # ============================================================
    # COMBINED MI + PERPLEXITY ANALYSIS (all experiments)
    # ============================================================
    print("\n=== Combined MI + Perplexity Analysis ===")
    try:
        from coevolution_shared import (
            combined_pair_scores, compute_entropy_vectorized,
            load_position_arrays as _lpa,
        )
        _pa, _na, _fl = _lpa(max_pos=None, aligned=True)
        _ent = compute_entropy_vectorized(_pa, _na, _fl)
        _var = [p for p in range(_fl) if _ent[p] > 0.3]
        _pairs = [(i, j) for idx, i in enumerate(_var)
                  for j in _var[idx + 1:] if j - i <= 30]
        _scored = combined_pair_scores(_pa, _pairs, _na, _ent)
        print(f"  Variable positions: {len(_var)}")
        print(f"  Pairs scored (MI + perplexity ratio): {len(_scored)}")
        print(f"  Top 5 combined (MI + ratio):")
        for _s in _scored[:5]:
            print(f"    ({_s['pos_i']},{_s['pos_j']}): MI={_s['mi']:.3f} "
                  f"ratio={_s['ratio']:.2f} combined={_s['combined']:.3f}")
        _mi_top = sorted(_scored, key=lambda s: -s['mi'])[:5]
        print(f"  Top 5 by MI alone:")
        for _s in _mi_top:
            print(f"    ({_s['pos_i']},{_s['pos_j']}): MI={_s['mi']:.3f} "
                  f"ratio={_s['ratio']:.2f}")
        # ranking agreement
        _r_mi = {(_s['pos_i'], _s['pos_j']): idx
                 for idx, _s in enumerate(sorted(_scored, key=lambda s: -s['mi']))}
        _r_cb = {(_s['pos_i'], _s['pos_j']): idx
                 for idx, _s in enumerate(_scored)}
        _same = sum(1 for k in _r_mi if _r_mi[k] == _r_cb[k])
        print(f"  Ranking agreement (MI vs combined, top-5 same): "
              f"{len([k for k in _r_mi if k in _r_cb and _r_mi[k] < 5 and _r_cb[k] < 5])}/5")
    except Exception as _e:
        print(f"  Combined analysis skipped: {_e}")


if __name__ == "__main__":
    main()
