#!/usr/bin/env python3
"""
Shared Co-evolution Data Loading and Computation Module
======================================================
SINGLE source of truth for all 17 co-evolution scripts:
  - FASTA parsing → position arrays (He 2012 encoding)
  - Vectorized Mutual Information (bincount, O(400) per pair)
  - Entropy / Perplexity / Variable position detection
  - Coupling J = ln(P/P_expected) and constraint function
  - Multi-core worker pool with shared-memory initializer

This ELIMINATES:
  - 19 independent FASTA parsers    →  1 implementation
  -  6 independent MI routines      →  1 numpy-vectorized
  - Repeated 26 MB data pickling    →  data loaded once per worker
"""

from __future__ import annotations

import os
import sys
import json
import numpy as np
from collections import Counter
from pathlib import Path
from typing import Optional

sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")
from nkmap.encoding.bio_sequences import Base20AminoEncoder, AMINO_HE_2012

AA_LIST = list(AMINO_HE_2012)
N_AA = 20
_CACHE: dict = {}


# ══════════════════════════════════════════════════════════════════════
#  1. FASTA parser
# ══════════════════════════════════════════════════════════════════════


def parse_fasta(filepath):
    """Parse FASTA alignment -> list[(header, sequence)]."""
    sequences = []
    cur_h = None
    cur_s = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur_h is not None:
                    sequences.append((cur_h, "".join(cur_s)))
                cur_h = line[1:]
                cur_s = []
            elif line:
                cur_s.append(line.upper())
    if cur_h is not None:
        sequences.append((cur_h, "".join(cur_s)))
    return sequences


# ══════════════════════════════════════════════════════════════════════
#  2. Position arrays (lazy, cached)
# ══════════════════════════════════════════════════════════════════════


def load_position_arrays(
    fasta_path=None, max_pos=None, n_seqs=None, clear_cache=False, aligned=True
):
    """Load FASTA, He 2012 encode, return (pos_arrays, n_all, full_len).

    max_pos=None -> FULL sequence length (all positions).
    Cached on first call. Pass clear_cache=True to force reload.

    aligned=True (CORRECTED, default): keep the full alignment; encode '-'
    and unknown characters as 20 (21 states). Every sequence has the same
    length, so column j is the same raw alignment position in every sequence.

    aligned=False (LEGACY, BUGGY): strip gaps, so column j mixes different
    raw positions across sequences. This was the original behavior and is
    kept only for reproducing the buggy results; it produces artifacts
    (1,249 fake 'variable' positions, fake MI hubs). See analysis/
    corrected_pipeline.py and the correction document.
    """
    if fasta_path is None:
        fasta_path = os.environ.get("COEVO_FASTA") or (Path(__file__).resolve().parent / "Spike_protein.aln-fasta")
    fasta_path = Path(fasta_path)

    cache_key = ("pos_arrays", str(fasta_path), max_pos, n_seqs, aligned)
    if not clear_cache and cache_key in _CACHE:
        return _CACHE[cache_key]

    encoder = Base20AminoEncoder(version=1)
    seqs_raw = parse_fasta(fasta_path)
    n_all = min(n_seqs, len(seqs_raw)) if n_seqs else len(seqs_raw)
    full_len = len(seqs_raw[0][1]) if seqs_raw else 0
    if max_pos is None:
        actual_max = full_len  # FULL length
    else:
        actual_max = min(max_pos, full_len)

    pos_arrays = []
    for _, seq in seqs_raw[:n_all]:
        if aligned:
            # CORRECTED: keep alignment, gap/unknown = 20
            arr = np.array(
                [encoder.encode.get(aa, 20) for aa in seq[:actual_max]],
                dtype=np.int32,
            )
        else:
            # LEGACY (buggy): strip gaps -> column misalignment
            clean = "".join(aa for aa in seq if aa in encoder.encode)
            arr = np.array(
                [encoder.encode.get(aa, -1) for aa in clean[:actual_max]],
                dtype=np.int32,
            )
        pos_arrays.append(arr)

    result = (pos_arrays, n_all, full_len)
    _CACHE[cache_key] = result
    return result


# ══════════════════════════════════════════════════════════════════════
#  3. Entropy / Perplexity
# ══════════════════════════════════════════════════════════════════════


def compute_entropy(pos_arrays, pos, n_seqs):
    """Shannon entropy at a single position."""
    counts = Counter()
    for arr in pos_arrays[:n_seqs]:
        if pos < len(arr):
            c = int(arr[pos])
            if 0 <= c < N_AA:  # exclude gap (20)
                counts[c] += 1
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((c / total) * np.log2(c / total) for c in counts.values() if c > 0)


def compute_entropy_vectorized(pos_arrays, n_seqs, max_pos):
    """Entropy for all positions using numpy (fast)."""
    n = min(n_seqs, len(pos_arrays))
    dense = np.full((n, max_pos), -1, dtype=np.int32)
    for i, arr in enumerate(pos_arrays[:n]):
        L = min(len(arr), max_pos)
        dense[i, :L] = arr[:L]

    ent = np.zeros(max_pos, dtype=np.float64)
    for pos in range(max_pos):
        col = dense[:, pos]
        valid = col[(col >= 0) & (col < N_AA)]
        if len(valid) == 0:
            continue
        cnt = np.bincount(valid, minlength=N_AA).astype(np.float64)
        total = cnt.sum()
        if total == 0:
            continue
        p = cnt / total
        ent[pos] = -np.sum(p[p > 0] * np.log2(p[p > 0]))
    return ent


def perplexity(entropy):
    """Perplexity = 2^H."""
    return np.power(2.0, entropy)


def find_variable_positions(pos_arrays, n_seqs, max_pos=80, threshold=0.3):
    """List positions with entropy > threshold."""
    ent = compute_entropy_vectorized(pos_arrays, n_seqs, max_pos)
    return [p for p in range(max_pos) if ent[p] > threshold]


def majority_ref(pos_arrays, pos, n_seqs):
    """Most common residue code at a position (gaps excluded).

    Ties are broken by LOWEST code index — this matches the GPU kernel
    (argmax over one-hot counts), keeping CPU/GPU references consistent.
    """
    cnt = Counter(
        int(a[pos])
        for a in pos_arrays[:n_seqs]
        if pos < len(a) and 0 <= int(a[pos]) < N_AA
    )
    if not cnt:
        return 0
    best = max(cnt.values())
    return min(code for code, c in cnt.items() if c == best)


# ══════════════════════════════════════════════════════════════════════
#  4. Mutual Information — SINGLE vectorized implementation
# ══════════════════════════════════════════════════════════════════════


def mutual_information(pos_arrays, pos_i, pos_j, n_seqs):
    """MI(pos_i, pos_j) via numpy bincount — O(n_seqs + 400), NOT O(n²).

    This is the ONLY MI implementation. All scripts should use this.
    It replaces 6 independent, slower Counter-based implementations.
    """
    # Extract two columns
    codes_i = np.array(
        [arr[pos_i] for arr in pos_arrays[:n_seqs] if pos_i < len(arr)],
        dtype=np.int32,
    )
    codes_j = np.array(
        [arr[pos_j] for arr in pos_arrays[:n_seqs] if pos_j < len(arr)],
        dtype=np.int32,
    )
    min_len = min(len(codes_i), len(codes_j))
    if min_len < 10:
        return 0.0
    codes_i = codes_i[:min_len]
    codes_j = codes_j[:min_len]
    valid = (codes_i >= 0) & (codes_i < N_AA) & (codes_j >= 0) & (codes_j < N_AA)
    codes_i = codes_i[valid]
    codes_j = codes_j[valid]
    if len(codes_i) < 10:
        return 0.0

    # Joint distribution: flat index = ci * 20 + cj
    pairs = codes_i.astype(np.int64) * N_AA + codes_j.astype(np.int64)
    joint_flat = (
        np.bincount(pairs, minlength=N_AA * N_AA).reshape(N_AA, N_AA).astype(np.float64)
    )
    total = joint_flat.sum()
    if total == 0:
        return 0.0
    marg_i = joint_flat.sum(axis=1)
    marg_j = joint_flat.sum(axis=0)

    mi = 0.0
    for ai in range(N_AA):
        for aj in range(N_AA):
            if joint_flat[ai, aj] > 0 and marg_i[ai] > 0 and marg_j[aj] > 0:
                p = joint_flat[ai, aj] / total
                pi_val = marg_i[ai] / total
                pj_val = marg_j[aj] / total
                mi += p * np.log2(p / (pi_val * pj_val))
    return float(mi)


def compute_mi_matrix(pos_arrays, n_seqs, max_pos, window=30):
    """Full MI matrix for positions 0..max_pos-1 within sliding window."""
    mi_mat = np.zeros((max_pos, max_pos), dtype=np.float64)
    for i in range(max_pos):
        for j in range(i + 1, min(i + window, max_pos)):
            m = mutual_information(pos_arrays, i, j, n_seqs)
            mi_mat[i, j] = m
            mi_mat[j, i] = m
    return mi_mat


def mi_mutation_only(pos_arrays, pos_i, pos_j, ref_i, ref_j, n_seqs):
    """MI over mutation pairs only (excludes reference)."""
    joint = Counter()
    marg_i = Counter()
    marg_j = Counter()
    for arr in pos_arrays[:n_seqs]:
        if pos_i < len(arr) and pos_j < len(arr):
            ci, cj = int(arr[pos_i]), int(arr[pos_j])
            if 0 <= ci < N_AA and 0 <= cj < N_AA and (ci != ref_i or cj != ref_j):
                joint[(ci, cj)] += 1
                marg_i[ci] += 1
                marg_j[cj] += 1
    total = sum(joint.values())
    if total < 5:
        return 0.0, 0
    mi_val = sum(
        (c / total)
        * np.log2((c / total) / ((marg_i[ai] / total) * (marg_j[aj] / total)))
        for (ai, aj), c in joint.items()
        if marg_i[ai] > 0 and marg_j[aj] > 0
    )
    return float(mi_val), total

def compute_mi_matrix_gpu(pos_arrays, pairs, n_positions, min_total=10, chunk=4096):
    """Full MI matrix over ALL given position pairs, GPU (torch CUDA) with
    vectorized numpy fallback. Uses ALL sequences and ALL positions — no caps.

    Returns (mi_matrix, n_computed)."""
    mi_matrix = np.zeros((n_positions, n_positions), dtype=np.float64)
    try:
        import coevolution_gpu as cg
        dense = cg.dense_to_gpu(pos_arrays)
        mi_dict, _ = cg.mi_matrix_gpu(dense, pairs, min_total=min_total, chunk=chunk)
        for (i, j), mi in mi_dict.items():
            mi_matrix[i, j] = mi
            mi_matrix[j, i] = mi
        return mi_matrix, len(mi_dict)
    except Exception as e:
        print(f"    GPU MI failed ({e}); vectorized numpy fallback")
    dense = np.full((len(pos_arrays), n_positions), -1, dtype=np.int32)
    for si, arr in enumerate(pos_arrays):
        L = min(len(arr), n_positions)
        dense[si, :L] = arr[:L]
    n = 0
    for (i, j) in pairs:
        codes_i = dense[:, i]
        codes_j = dense[:, j]
        valid = (codes_i >= 0) & (codes_j >= 0) & (codes_i < N_AA) & (codes_j < N_AA)
        ci, cj = codes_i[valid], codes_j[valid]
        if len(ci) < min_total:
            continue
        joint = np.bincount(ci.astype(np.int64) * N_AA + cj.astype(np.int64),
                            minlength=N_AA * N_AA).reshape(N_AA, N_AA).astype(np.float64)
        total = joint.sum()
        if total == 0:
            continue
        marg_i = joint.sum(axis=1)
        marg_j = joint.sum(axis=0)
        mi = 0.0
        for ai in range(N_AA):
            for aj in range(N_AA):
                if joint[ai, aj] > 0 and marg_i[ai] > 0 and marg_j[aj] > 0:
                    p = joint[ai, aj] / total
                    pi_v = marg_i[ai] / total
                    pj_v = marg_j[aj] / total
                    mi += p * np.log2(p / (pi_v * pj_v))
        mi_matrix[i, j] = mi
        mi_matrix[j, i] = mi
        n += 1
    return mi_matrix, n


# ══════════════════════════════════════════════════════════════════════
#  5. Coupling constants, constraint function, K-map builders
# ══════════════════════════════════════════════════════════════════════


def compute_coupling(pos_arrays, pos_i, pos_j, n_seqs):
    """J = ln(P / P_expected) — DCA coupling convention.

    SIGN CONVENTION (verified against DCA literature, Morcos et al 2011 PNAS):
      J > 0 → pair is MORE common than expected (co-evolutionary / positively coupled)
      J < 0 → pair is LESS common than expected (anti-correlated / negatively selected)
      J = 0 → pair occurs at random (independent)

    NOTE: The previous formula used -ln(P/P_exp) which gave the OPPOSITE sign
    (negative for co-evolutionary). This was inconsistent with the docstrings
    and broke the sigmoid prediction σ(J) = 1/(1+e^{-J}).
    """
    kmap = np.zeros((N_AA, N_AA), dtype=np.float64)
    for arr in pos_arrays[:n_seqs]:
        if pos_i < len(arr) and pos_j < len(arr):
            ci, cj = int(arr[pos_i]), int(arr[pos_j])
            if 0 <= ci < N_AA and 0 <= cj < N_AA:
                kmap[ci, cj] += 1

    total = kmap.sum()
    if total == 0:
        return np.zeros((N_AA, N_AA)), 0.0
    kmap /= total
    mi_marg = kmap.sum(axis=1)
    mj_marg = kmap.sum(axis=0)
    eps = 1e-10
    with np.errstate(divide="ignore", invalid="ignore"):
        J = np.log((kmap + eps) / ((mi_marg[:, None] + eps) * (mj_marg[None, :] + eps)))
    return J, float(np.mean(np.abs(J)))


def constraint_function(kmap_freq):
    """C(aa_i, aa_j) = ln(P / P_expected) — constraint / coupling function.

    SIGN CONVENTION (verified against DCA literature):
      C > 0 → pair is MORE common than expected (co-evolutionary)
      C < 0 → pair is LESS common than expected (anti-correlated)
      C = 0 → pair occurs at random frequency (independent)

    Prediction: P_co-evolution = σ(C) = 1/(1+e^{-C})
    """
    mi_marg = kmap_freq.sum(axis=1)
    mj_marg = kmap_freq.sum(axis=0)
    eps = 1e-10
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(
            (kmap_freq + eps) / ((mi_marg[:, None] + eps) * (mj_marg[None, :] + eps))
        )


def build_mutation_kmap(pos_arrays, pos_i, pos_j, ref_i, ref_j, n_seqs):
    """20x20 K-map: 1=mutation, -1=reference."""
    kmap = np.zeros((N_AA, N_AA), dtype=np.int32)
    for arr in pos_arrays[:n_seqs]:
        if pos_i < len(arr) and pos_j < len(arr):
            ci, cj = int(arr[pos_i]), int(arr[pos_j])
            if 0 <= ci < N_AA and 0 <= cj < N_AA:
                kmap[ci, cj] = 1 if (ci != ref_i or cj != ref_j) else -1
    return kmap


def build_frequency_kmap(pos_arrays, pos_i, pos_j, n_seqs, exclude_idx=None):
    """20x20 frequency K-map, optionally excluding one sequence (for LOO-CV)."""
    kmap = np.zeros((N_AA, N_AA), dtype=np.float64)
    for idx, arr in enumerate(pos_arrays[:n_seqs]):
        if exclude_idx is not None and idx == exclude_idx:
            continue
        if pos_i < len(arr) and pos_j < len(arr):
            ci, cj = int(arr[pos_i]), int(arr[pos_j])
            if 0 <= ci < N_AA and 0 <= cj < N_AA:
                kmap[ci, cj] += 1
    total = kmap.sum()
    if total > 0:
        kmap /= total
    return kmap


# ══════════════════════════════════════════════════════════════════════
#  6. Multi-core workers with SHARED MEMORY (no pickling overhead)
# ══════════════════════════════════════════════════════════════════════

# Module-level globals set by pool initializer — data loaded ONCE per worker
_WORKER_DATA: list = []
_WORKER_N: int = 0


def _init_worker(pos_arrays, n_seqs):
    """Pool initializer: store shared data in worker global scope."""
    global _WORKER_DATA, _WORKER_N
    _WORKER_DATA = pos_arrays
    _WORKER_N = n_seqs


def _mi_worker(pair):
    """Compute MI for one position pair using shared worker data.
    No pickling of 26MB arrays — data is already in worker memory.
    """
    global _WORKER_DATA, _WORKER_N
    pi, pj = pair
    mi = mutual_information(_WORKER_DATA, pi, pj, _WORKER_N)
    ref_i = majority_ref(_WORKER_DATA, pi, _WORKER_N)
    ref_j = majority_ref(_WORKER_DATA, pj, _WORKER_N)
    mi_mut, n_mut = mi_mutation_only(_WORKER_DATA, pi, pj, ref_i, ref_j, _WORKER_N)
    return (int(pi), int(pj), float(mi), int(n_mut))


def get_worker_count():
    """Number of CPU cores, capped at 24."""
    n = os.cpu_count() or 4
    return min(n, 24)


# ══════════════════════════════════════════════════════════════════════
#  7. GPU-accelerated pair finder (torch CUDA, A100)
# ══════════════════════════════════════════════════════════════════════


def find_coevolving_pairs_gpu(
    pos_arrays,
    variable_positions,
    n_seqs,
    max_gap=30,
    min_mi=0.1,
    min_muts=5,
    mutation_only=True,
):
    """Find co-evolving position pairs with MI computed on GPU.

    Returns list of (pos_i, pos_j, mi, n_muts, ref_i, ref_j) sorted by MI desc.
    Falls back to CPU (mutual_information) if CUDA unavailable.
    """
    try:
        import coevolution_gpu as cg

        dense = cg.dense_to_gpu(pos_arrays[:n_seqs])
        refs = cg.majority_refs_gpu(dense)
        pairs = [
            (pi, pj)
            for idx_i, pi in enumerate(variable_positions)
            for pj in variable_positions[idx_i + 1 :]
            if abs(pi - pj) <= max_gap
        ]
        mi_dict, cnt_dict = cg.mi_matrix_gpu(
            dense,
            pairs,
            refs=refs,
            mutation_only=mutation_only,
            min_total=min_muts,
            chunk=4096,
        )
        results = []
        for (pi, pj), mi in mi_dict.items():
            if mi > min_mi:
                results.append(
                    (pi, pj, mi, cnt_dict[(pi, pj)], int(refs[pi]), int(refs[pj]))
                )
        results.sort(key=lambda x: x[2], reverse=True)
        return results
    except Exception:
        # CPU fallback: mutation-only MI via shared mutual_information
        results = []
        for idx_i, pi in enumerate(variable_positions):
            for pj in variable_positions[idx_i + 1 :]:
                if abs(pi - pj) > max_gap:
                    continue
                ri = majority_ref(pos_arrays, pi, n_seqs)
                rj = majority_ref(pos_arrays, pj, n_seqs)
                if mutation_only:
                    mi, n_mut = mi_mutation_only(pos_arrays, pi, pj, ri, rj, n_seqs)
                else:
                    mi = mutual_information(pos_arrays, pi, pj, n_seqs)
                    n_mut = 0
                if mi > min_mi:
                    results.append((pi, pj, mi, n_mut, ri, rj))
        results.sort(key=lambda x: x[2], reverse=True)
        return results


# ══════════════════════════════════════════════════════════════════════
#  8. Combined MI + Perplexity analysis (all experiments)
# ══════════════════════════════════════════════════════════════════════

def perplexity_ratio(pos_arrays, pos_i, pos_j, n_seqs, entropy_vec=None):
    """PP(j) / mean PP(j|i=a) over residues a with n>=5.

    Marginal PP uses 2^H(j); conditional PP uses 2^H(j|i=a). Ratio > 1
    means knowing residue i reduces the effective number of choices at j
    (determinism). Returns None if fewer than 3 conditioning residues.

    VECTORIZED (numpy bincount) — identical semantics to the original
    Counter version but ~100x faster on deep alignments (verified equal
    on random samples).
    """
    codes_i = np.array(
        [int(a[pos_i]) for a in pos_arrays[:n_seqs] if pos_i < len(a)],
        dtype=np.int32,
    )
    codes_j = np.array(
        [int(a[pos_j]) for a in pos_arrays[:n_seqs] if pos_j < len(a)],
        dtype=np.int32,
    )
    m = min(len(codes_i), len(codes_j))
    codes_i, codes_j = codes_i[:m], codes_j[:m]
    valid = ((codes_i >= 0) & (codes_i < N_AA)
             & (codes_j >= 0) & (codes_j < N_AA))
    codes_i, codes_j = codes_i[valid], codes_j[valid]
    if len(codes_i) < 5:
        return None
    joint = np.bincount(codes_i.astype(np.int64) * N_AA + codes_j.astype(np.int64),
                        minlength=N_AA * N_AA).reshape(N_AA, N_AA).astype(np.float64)
    tot = joint.sum()
    if tot == 0:
        return None
    ci_tot = joint.sum(axis=1)
    cond = []
    for ci in range(N_AA):
        if ci_tot[ci] >= 5:
            p = joint[ci] / ci_tot[ci]
            p = p[p > 0]
            if len(p):
                h = -float(np.sum(p * np.log2(p)))
                cond.append(2.0 ** h)
    if not cond:
        return None
    if entropy_vec is not None:
        ppj = 2.0 ** float(entropy_vec[pos_j])
    else:
        marg_j = joint.sum(axis=0)
        pj = marg_j[marg_j > 0] / tot
        h_j = -float(np.sum(pj * np.log2(pj)))
        ppj = 2.0 ** h_j
    avg_cond = float(np.mean(cond))
    if avg_cond <= 0:
        return None
    return ppj / avg_cond


def dense_from_arrays(pos_arrays, n_seqs=None):
    """Dense [n_seqs, L] int32 matrix (gaps=-1) — one build, many slices."""
    n = min(n_seqs, len(pos_arrays)) if n_seqs else len(pos_arrays)
    L = max(len(a) for a in pos_arrays[:n]) if n else 0
    dense = np.full((n, L), -1, dtype=np.int32)
    for i, arr in enumerate(pos_arrays[:n]):
        dense[i, :len(arr)] = arr[:L]
    return dense


def _mi_dense(dense, i, j):
    """MI from dense column slices (no per-pair list comprehensions)."""
    ci, cj = dense[:, i], dense[:, j]
    valid = (ci >= 0) & (ci < N_AA) & (cj >= 0) & (cj < N_AA)
    ci, cj = ci[valid], cj[valid]
    if len(ci) < 10:
        return 0.0
    joint = np.bincount(ci.astype(np.int64) * N_AA + cj.astype(np.int64),
                        minlength=N_AA * N_AA).reshape(N_AA, N_AA).astype(np.float64)
    total = joint.sum()
    if total == 0:
        return 0.0
    mi = 0.0
    for a in range(N_AA):
        for b in range(N_AA):
            if joint[a, b] > 0:
                pa = joint[a, :].sum() / total
                pb = joint[:, b].sum() / total
                mi += (joint[a, b] / total) * np.log2((joint[a, b] / total) / (pa * pb))
    return mi


def _ratio_dense(dense, i, j, entropy_vec=None):
    """Perplexity ratio from dense column slices (vectorized)."""
    ci, cj = dense[:, i], dense[:, j]
    valid = (ci >= 0) & (ci < N_AA) & (cj >= 0) & (cj < N_AA)
    ci, cj = ci[valid], cj[valid]
    if len(ci) < 5:
        return None
    joint = np.bincount(ci.astype(np.int64) * N_AA + cj.astype(np.int64),
                        minlength=N_AA * N_AA).reshape(N_AA, N_AA).astype(np.float64)
    tot = joint.sum()
    if tot == 0:
        return None
    ci_tot = joint.sum(axis=1)
    cond = []
    for a in range(N_AA):
        if ci_tot[a] >= 5:
            p = joint[a] / ci_tot[a]
            p = p[p > 0]
            if len(p):
                h = -float(np.sum(p * np.log2(p)))
                cond.append(2.0 ** h)
    if not cond:
        return None
    if entropy_vec is not None:
        ppj = 2.0 ** float(entropy_vec[j])
    else:
        marg_j = joint.sum(axis=0)
        pj = marg_j[marg_j > 0] / tot
        ppj = 2.0 ** (-float(np.sum(pj * np.log2(pj))))
    avg_cond = float(np.mean(cond))
    if avg_cond <= 0:
        return None
    return ppj / avg_cond


def combined_pair_scores(pos_arrays, pairs, n_seqs, entropy_vec=None,
                         mi_fn=None):
    """Combined MI + perplexity-ratio score for a list of pairs.

    Returns list of dicts: {pos_i, pos_j, mi, ratio, combined}
      combined = rank-normalized MI + rank-normalized ratio (averaged).
    mi_fn defaults to mutual_information; pass a mutation-only MI function
    to match the pipeline convention. Dense-matrix fast path when mi_fn is
    None (no per-pair Python column extraction — ~100x faster on deep MSAs).
    """
    if mi_fn is None:
        dense = dense_from_arrays(pos_arrays, n_seqs)
        scored = []
        for (i, j) in pairs:
            mi = _mi_dense(dense, i, j)
            ratio = _ratio_dense(dense, i, j, entropy_vec)
            if mi is not None and ratio is not None:
                scored.append({"pos_i": int(i), "pos_j": int(j),
                               "mi": float(mi), "ratio": float(ratio)})
        if not scored:
            return scored
        n = len(scored)
        mi_vals = np.array([s["mi"] for s in scored])
        ra_vals = np.array([s["ratio"] for s in scored])
        def rnorm(v):
            order = np.argsort(np.argsort(v))
            return order / (n - 1) if n > 1 else np.zeros(n)
        rm = rnorm(mi_vals)
        rr = rnorm(ra_vals)
        for idx, s in enumerate(scored):
            s["combined"] = float(0.5 * (rm[idx] + rr[idx]))
        scored.sort(key=lambda s: -s["combined"])
        return scored
    scored = []
    for (i, j) in pairs:
        mi = mi_fn(pos_arrays, i, j, n_seqs)
        ratio = perplexity_ratio(pos_arrays, i, j, n_seqs, entropy_vec)
        if mi is not None and ratio is not None:
            scored.append({"pos_i": int(i), "pos_j": int(j),
                           "mi": float(mi), "ratio": float(ratio)})
    if not scored:
        return scored
    n = len(scored)
    mi_vals = np.array([s["mi"] for s in scored])
    ra_vals = np.array([s["ratio"] for s in scored])
    def rnorm(v):
        order = np.argsort(np.argsort(v))
        return order / (n - 1) if n > 1 else np.zeros(n)
    rm = rnorm(mi_vals)
    rr = rnorm(ra_vals)
    for idx, s in enumerate(scored):
        s["combined"] = float(0.5 * (rm[idx] + rr[idx]))
    scored.sort(key=lambda s: -s["combined"])
    return scored
