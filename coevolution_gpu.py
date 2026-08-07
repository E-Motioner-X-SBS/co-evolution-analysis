#!/usr/bin/env python3
"""
GPU-Accelerated Co-evolution Kernels (torch CUDA on A100)
==========================================================
Full-length, all-sequence computation kernels:
  - Dense position tensor on GPU (n_seqs × max_pos)
  - Majority references per position (vectorized argmax over counts)
  - Shannon entropy per position (vectorized one-hot + log2)
  - Mutual Information for ANY pair set (batched scatter_add bincount)
  - Mutation-only MI (excludes reference-reference pairs)
  - Coupling constants J = ln(P/P_exp) for arbitrary pairs
  - H1 Gray-adjacency ratio for consecutive residues

All kernels operate on the FULL sequence length and ALL sequences.
No position truncation, no sequence subsampling.
"""

from __future__ import annotations

import numpy as np
import torch

N_AA = 20
_log2 = torch.log2


def check_cuda() -> bool:
    """Return True if CUDA is available."""
    return torch.cuda.is_available()


def get_device() -> torch.device:
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def dense_to_gpu(pos_arrays, device=None) -> torch.Tensor:
    """Build dense [n_seqs, max_pos] int32 tensor on GPU.

    pos_arrays: list of np.int32 arrays (codes 0-19, -1 for gap/unknown).
    Returns tensor on CUDA (or CPU fallback).
    """
    if device is None:
        device = get_device()
    max_pos = max(len(a) for a in pos_arrays) if pos_arrays else 0
    n = len(pos_arrays)
    dense = torch.full((n, max_pos), -1, dtype=torch.int32, device=device)
    for i, arr in enumerate(pos_arrays):
        L = len(arr)
        dense[i, :L] = torch.from_numpy(arr[:L].astype(np.int32)).to(device)
    return dense


def majority_refs_gpu(dense: torch.Tensor) -> torch.Tensor:
    """Most common valid residue code per position (gaps excluded) → int32 [max_pos]."""
    max_pos = dense.shape[1]
    valid = ((dense >= 0) & (dense < N_AA)).float()  # exclude gap (20) and -1
    idx = dense.clamp(min=0, max=N_AA - 1).long()  # clamp so scatter never OOB
    onehot = torch.zeros(max_pos, N_AA, dtype=torch.float32, device=dense.device)
    onehot.scatter_add_(1, idx.T, valid.T)  # [max_pos, 20] counts
    return onehot.argmax(dim=1).to(torch.int32)


def compute_entropy_gpu(dense: torch.Tensor) -> torch.Tensor:
    """Shannon entropy (bits) per position, gaps excluded → float [max_pos]."""
    max_pos = dense.shape[1]
    valid = ((dense >= 0) & (dense < N_AA)).float()
    idx = dense.clamp(min=0, max=N_AA - 1).long()
    onehot = torch.zeros(max_pos, N_AA, dtype=torch.float32, device=dense.device)
    onehot.scatter_add_(1, idx.T, valid.T)  # [max_pos, 20]
    totals = onehot.sum(dim=1, keepdim=True).clamp(min=1.0)
    p = onehot / totals
    mask = p > 0
    ent = -(p * _log2(p.clamp(min=1e-30)) * mask).sum(dim=1)
    return ent


def mi_matrix_gpu(
    dense: torch.Tensor,
    pairs,
    refs: torch.Tensor | None = None,
    mutation_only: bool = False,
    min_total: int = 5,
    chunk: int = 4096,
):
    """Mutual Information for an arbitrary list of (i, j) pairs.

    Args:
        dense: [n_seqs, max_pos] int32 GPU tensor.
        pairs: iterable of (i, j) position pairs (0-based).
        refs: [max_pos] int32 ref codes; required if mutation_only=True.
        mutation_only: exclude reference-reference pairs from joint counts.
        min_total: minimum valid observations; below → MI = 0.
        chunk: pairs per batch.

    Returns:
        (mi_dict, count_dict) — {(i,j): float MI}, {(i,j): n_observations}.
    """
    device = dense.device
    n = dense.shape[0]
    pairs = list(pairs)
    mi_dict: dict = {}
    cnt_dict: dict = {}

    for start in range(0, len(pairs), chunk):
        cp = pairs[start : start + chunk]
        P = len(cp)
        ii = torch.tensor([p[0] for p in cp], dtype=torch.long, device=device)
        jj = torch.tensor([p[1] for p in cp], dtype=torch.long, device=device)

        ci = dense[:, ii]  # [n, P]
        cj = dense[:, jj]
        # CORRECTED (FIX A1): aligned arrays use gap = 20; exclude gaps (-1
        # legacy, 20 = gap) so MI is computed only over observed residue pairs.
        valid = (ci >= 0) & (ci < N_AA) & (cj >= 0) & (cj < N_AA)  # [n, P]

        if mutation_only and refs is not None:
            ref_i = refs[ii]  # [P]
            ref_j = refs[jj]
            is_ref = (ci == ref_i[None, :]) & (cj == ref_j[None, :])
            valid = valid & ~is_ref

        # CORRECTED (FIX A1): valid = residues only (0-19); gaps/legacy excluded.
        # Clamp flat to [0, 399] so scatter_add never sees an out-of-range index.
        flat = (
            ci.clamp(min=0, max=N_AA - 1).long() * N_AA
            + cj.clamp(min=0, max=N_AA - 1).long()
        )
        flat = flat.masked_fill(~valid, 0)

        joint = torch.zeros(P, N_AA * N_AA, dtype=torch.float32, device=device)
        joint.scatter_add_(1, flat.T, valid.float().T)  # [P, 400]
        joint = joint.reshape(P, N_AA, N_AA)

        totals = joint.sum(dim=(1, 2))  # [P]
        valid_tot = totals >= min_total
        safe_t = totals.clamp(min=1.0)

        marg_i = joint.sum(dim=2)  # [P, 20]
        marg_j = joint.sum(dim=1)
        p = joint / safe_t[:, None, None]
        pi = marg_i / safe_t[:, None]
        pj = marg_j / safe_t[:, None]
        ratio = p / (pi[:, :, None] * pj[:, None, :]).clamp(min=1e-30)
        contrib = torch.where(
            p > 0, p * _log2(ratio.clamp(min=1e-30)), torch.zeros_like(p)
        )
        mi = contrib.sum(dim=(1, 2))  # [P]

        # Zero out pairs with insufficient observations
        mi = torch.where(valid_tot, mi, torch.zeros_like(mi))

        mi_cpu = mi.cpu().numpy()
        tot_cpu = totals.cpu().numpy().astype(np.int64)
        for k in range(P):
            i, j = cp[k]
            mi_dict[(i, j)] = float(mi_cpu[k])
            cnt_dict[(i, j)] = int(tot_cpu[k])

    return mi_dict, cnt_dict


def coupling_matrix_gpu(
    dense: torch.Tensor,
    pos_i: int,
    pos_j: int,
    refs: torch.Tensor | None = None,
    mutation_only: bool = False,
    min_total: int = 5,
) -> np.ndarray:
    """Coupling J = ln(P(aa_i,aa_j) / (P(aa_i)·P(aa_j))) for one pair.

    Returns 20×20 float numpy array. J > 0 → co-evolutionary.
    """
    device = dense.device
    ci = dense[:, pos_i].long()
    cj = dense[:, pos_j].long()
    valid = (ci >= 0) & (ci < N_AA) & (cj >= 0) & (cj < N_AA)
    if mutation_only and refs is not None:
        is_ref = (ci == int(refs[pos_i])) & (cj == int(refs[pos_j]))
        valid = valid & ~is_ref
    flat = ci.clamp(min=0, max=N_AA - 1) * N_AA + cj.clamp(min=0, max=N_AA - 1)
    flat = flat.masked_fill(~valid, 0)
    joint = torch.zeros(N_AA * N_AA, dtype=torch.float32, device=device)
    joint.scatter_add_(0, flat, valid.float())
    total = joint.sum()
    if total < min_total:
        return np.zeros((N_AA, N_AA), dtype=np.float64)
    joint = joint.reshape(N_AA, N_AA) / total
    marg_i = joint.sum(dim=1)
    marg_j = joint.sum(dim=0)
    eps = 1e-10
    J = torch.log((joint + eps) / (marg_i[:, None] * marg_j[None, :] + eps))
    return J.cpu().numpy().astype(np.float64)


def h1_adjacency_gpu(dense: torch.Tensor) -> tuple[float, int, int]:
    """H1: fraction of consecutive residue pairs at Gray Hamming distance 1.

    Uses 5-bit Gray code g(i) = i ^ (i >> 1) on the He-2012 code directly.
    Returns (ratio, n_dist1, n_total).
    """
    device = dense.device
    n, max_pos = dense.shape
    a = dense[:, :-1].long()  # [n, max_pos-1]
    b = dense[:, 1:].long()
    valid = (a >= 0) & (b >= 0)
    ga = a.clamp(min=0) ^ (a.clamp(min=0) >> 1)
    gb = b.clamp(min=0) ^ (b.clamp(min=0) >> 1)
    xor = ga ^ gb

    # popcount via bit tricks on tensor
    def popcount(x: torch.Tensor) -> torch.Tensor:
        x = x - ((x >> 1) & 0x55555555)
        x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
        x = (x + (x >> 4)) & 0x0F0F0F0F
        return (x * 0x01010101) >> 24

    d = popcount(xor)
    valid = valid & (d <= 5)
    n_total = int(valid.sum().item())
    n_dist1 = int((valid & (d == 1)).sum().item())
    ratio = n_dist1 / n_total if n_total else 0.0
    return ratio, n_dist1, n_total


def all_pairs(max_pos: int, max_gap: int | None = None):
    """All (i, j) pairs with i < j. max_gap=None → all pairs (full matrix)."""
    if max_gap is None:
        return [(i, j) for i in range(max_pos) for j in range(i + 1, max_pos)]
    return [
        (i, j)
        for i in range(max_pos)
        for j in range(i + 1, min(i + max_gap + 1, max_pos))
    ]


if __name__ == "__main__":
    # Self-test on small synthetic data
    rng = np.random.default_rng(0)
    pos_arrays = [rng.integers(-1, 20, size=50).astype(np.int32) for _ in range(30)]
    d = dense_to_gpu(pos_arrays)
    refs = majority_refs_gpu(d)
    ent = compute_entropy_gpu(d)
    print(
        "dense:",
        tuple(d.shape),
        "refs:",
        refs[:5].cpu().numpy(),
        "entropy[0]:",
        float(ent[0]),
    )
    pairs = all_pairs(50, max_gap=5)
    mi, cnt = mi_matrix_gpu(d, pairs)
    print("pairs:", len(pairs), "MI(0,1):", mi.get((0, 1)))
    J = coupling_matrix_gpu(d, 0, 1)
    print("coupling shape:", J.shape)
    h1, n1, nt = h1_adjacency_gpu(d)
    print("H1:", round(h1, 4), f"({n1}/{nt})")
    print("SELF-TEST OK")
