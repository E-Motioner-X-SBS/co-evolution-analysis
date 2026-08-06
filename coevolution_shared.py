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


def load_position_arrays(fasta_path=None, max_pos=None, n_seqs=None, clear_cache=False):
    """Load FASTA, He 2012 encode, return (pos_arrays, n_all, full_len).

    max_pos=None → FULL sequence length (all positions).
    Cached on first call. Pass clear_cache=True to force reload.
    """
    if fasta_path is None:
        fasta_path = Path(__file__).resolve().parent / "Spike_protein.aln-fasta"
    fasta_path = Path(fasta_path)

    cache_key = ("pos_arrays", str(fasta_path), max_pos, n_seqs)
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
        if pos < len(arr) and arr[pos] >= 0:
            counts[int(arr[pos])] += 1
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
        valid = col[col >= 0]
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
    """Most common residue code at a position."""
    cnt = Counter(
        int(a[pos]) for a in pos_arrays[:n_seqs] if pos < len(a) and a[pos] >= 0
    )
    if not cnt:
        return 0
    return cnt.most_common(1)[0][0]


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
    valid = (codes_i >= 0) & (codes_j >= 0)
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
            if ci >= 0 and cj >= 0 and (ci != ref_i or cj != ref_j):
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
            if ci >= 0 and cj >= 0:
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
            if ci >= 0 and cj >= 0:
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
            if ci >= 0 and cj >= 0:
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
            chunk=16384,
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
