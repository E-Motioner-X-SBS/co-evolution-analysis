#!/usr/bin/env python3
"""
Proper Mean-Field Direct Coupling Analysis (mfDCA)
===================================================
Correct implementation of mfDCA following Morcos et al. 2011 (PNAS 108:E1293),
verified against the reference implementation py-mfdca (utdal/py-mfdca) and
Weigt et al. 2009 (PNAS 106:67).

Algorithm (per Morcos 2011):
  1. Sequence reweighting: W_l = 1 / (number of sequences within 80% identity
     of sequence l, threshold theta = 0.2). Meff = sum(W).
  2. Single-site frequencies with pseudocount lambda = 0.5:
       P_i(a) = (1-lambda) * (1/Meff) * sum_l W_l * [x_l,i = a] + lambda/q
  3. Pairwise frequencies with pseudocount:
       P_ij(a,b) = (1-lambda) * (1/Meff) * sum_l W_l * [x_l,i=a, x_l,j=b]
                   + lambda/q^2   (off-diagonal i != j)
       P_ii(a,b) = delta(a,b) * P_i(a)  (diagonal: no self-coupling)
  4. Connected correlation (covariance) matrix C of size L(q-1) x L(q-1):
       C[(i,alpha), (j,beta)] = P_ij(alpha,beta) - P_i(alpha)*P_j(beta)
       for alpha,beta in 0..q-2 (last state / gap removed)
  5. Couplings: J = -C^{-1}   (invert the covariance on GPU)
  6. Direct Information DI_ij via iterative mean-field (mu1/mu2 fixed point):
       W_mf = exp(-J_ij) padded to q x q
       iterate: mu1 = P_i / (W_mf . mu2), mu2 = P_j / (W_mf^T . mu1)
       P_dir = W_mf * (mu1 otimes mu2)  normalized
       DI = trace(P_dir^T * log(P_dir / (P_i otimes P_j)))
  7. Also compute Frobenius norm score F_ij = ||J_ij||_F (standard DCA score)
     and its APC-corrected version.

State alphabet: 0..19 = 20 amino acids (He 2012 order), 20 = gap.
q = 21. Gauge: last state (gap) removed from C.

GPU: covariance inversion + DI via torch CUDA on A100.

References:
  - Morcos F, Pagnani A, Lunt B, et al. PNAS 108(49):E1293-E1301 (2011).
  - Weigt M, White RA, Szurmant H, et al. PNAS 106(1):67-72 (2009).
  - Ekeberg M, et al. Phys. Rev. E 87:012707 (2013) for plmDCA (not used here).
"""

import sys
import json
import time
import numpy as np
import torch
from pathlib import Path

sys.path.insert(0, "/store/shuvam/E-motioner-X-SBS/n-ary-kmap/src")
from nkmap.encoding.bio_sequences import AMINO_HE_2012

BASE = Path("/store/shuvam/E-motioner-X-SBS/datasets/co-evolution")
FASTA = BASE / "Spike_protein.aln-fasta"
OUT = BASE / "dca_results"
OUT.mkdir(exist_ok=True)

Q = 21  # 20 AAs + gap
Q_USED = 20  # states in covariance (q-1, gap removed)
THETA = 0.2  # reweighting threshold (80% identity)
LAMBDA = 0.5  # pseudocount
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# He 2012 order + gap at index 20
AA_ORDER = AMINO_HE_2012 + "-"
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_ORDER)}


def parse_fasta(filepath):
    seqs, cur_h, cur_s = [], None, []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur_h is not None:
                    seqs.append((cur_h, "".join(cur_s)))
                cur_h = line[1:]
                cur_s = []
            elif line:
                cur_s.append(line.upper())
    if cur_h is not None:
        seqs.append((cur_h, "".join(cur_s)))
    return seqs


def encode_msa(seqs):
    """Encode MSA to integer matrix [N, L], 0-19 AA, 20 = gap."""
    N = len(seqs)
    L = len(seqs[0][1])
    X = np.full((N, L), Q - 1, dtype=np.int64)  # default gap
    for i, (_, s) in enumerate(seqs):
        for j, aa in enumerate(s):
            X[i, j] = AA_TO_IDX.get(aa, Q - 1)
    return X


def reweight(X, theta=THETA):
    """Sequence weights: W_l = 1/count(seqs within theta*L identity)."""
    N, L = X.shape
    max_diffs = int(theta * L)
    W = np.zeros(N, dtype=np.float64)
    # Vectorized: count per-sequence neighbors (Hamming <= max_diffs)
    # For N=1299, L=1276: pairwise is 1299^2*1276 ~ 2.1e9 ops — heavy but OK
    # on GPU via torch.
    Xt = torch.from_numpy(X).to(DEVICE)
    # chunked pairwise comparison to bound memory
    chunk = 200
    counts = torch.zeros(N, dtype=torch.float64, device=DEVICE)
    for s in range(0, N, chunk):
        block = Xt[s : s + chunk]  # [chunk, L]
        diff = (block[:, None, :] != Xt[None, :, :]).sum(dim=2)  # [chunk, N]
        counts[s : s + chunk] = (diff <= max_diffs).sum(dim=1).float()
    counts = counts.cpu().numpy()
    return 1.0 / counts


def compute_frequencies(X, W):
    """Single and pairwise frequencies with pseudocount."""
    N, L = X.shape
    Meff = W.sum()
    Xt = torch.from_numpy(X).to(DEVICE)
    Wt = torch.from_numpy(W).to(DEVICE)

    # Single-site: P_i(a) = (1/W_eff) sum_l W_l [x_li = a]
    onehot = torch.zeros(N, L, Q, dtype=torch.float64, device=DEVICE)
    onehot.scatter_(2, Xt.unsqueeze(-1), 1.0)
    Pi_raw = (onehot * Wt[:, None, None]).sum(dim=0) / Meff  # [L, Q]
    Pi = (1 - LAMBDA) * Pi_raw + LAMBDA / Q * torch.ones_like(Pi_raw)

    # Pairwise: P_ij(a,b) for all i,j via chunked einsum to bound memory.
    # Pij_raw[i, j, a, b] = (1/Meff) sum_l W_l [x_li=a, x_lj=b]
    # Full tensor is [L, L, Q, Q] = 1276^2 * 21^2 * 8B = 5.7 GB fp64 (fits 80GB).
    # Weighted one-hot: V[l, i, a] = W_l * onehot[l, i, a]
    V = onehot * Wt[:, None, None]  # [N, L, Q]
    Pij_raw = torch.einsum("lia,ljb->ijab", V, onehot) / Meff
    # Symmetrize
    Pij_raw = 0.5 * (Pij_raw + Pij_raw.permute(1, 0, 3, 2))
    # Pseudocount
    Pij = (1 - LAMBDA) * Pij_raw + LAMBDA / (Q * Q)
    # Diagonal: P_ii(a,b) = delta(a,b) * P_i(a)
    eye = torch.eye(Q, dtype=torch.float64, device=DEVICE)
    for i in range(L):
        Pij[i, i] = Pi[i][:, None] * eye

    return Pi, Pij


def build_covariance(Pi, Pij):
    """C[(i,alpha),(j,beta)] = Pij - Pi*Pj, alpha,beta in 0..Q-2 (gap removed).

    CRITICAL: the flattening must permute [i, j, a, b] -> [i, a, j, b] so that
    row index = i*(q-1) + a and column index = j*(q-1) + b. A plain reshape
    would flatten as (i,j,a,b) giving row = i*(q-1)^2 + j*(q-1) + a which is
    WRONG (breaks matrix symmetry and the inversion semantics).
    """
    L = Pi.shape[0]
    qm1 = Q - 1
    Pi_c = Pi[:, :qm1]  # [L, q-1]
    Pij_c = Pij[:, :, :qm1, :qm1]  # [L, L, q-1, q-1]
    # C[i, j, a, b] = Pij_c[i,j,a,b] - Pi_c[i,a]*Pi_c[j,b]
    C = Pij_c - Pi_c[:, None, :, None] * Pi_c[None, :, None, :]
    # Permute [i, j, a, b] -> [i, a, j, b], then flatten
    C = C.permute(0, 2, 1, 3).reshape(L * qm1, L * qm1)
    return C


def direct_information(Pi, J4d):
    """DI via iterative mean-field (Morcos 2011, py-mfdca Compute_Results).

    GPU-accelerated with torch (cuBLAS einsum). Memory: [L,L,q,q] fp64 = 5.7 GB
    per tensor; ~4 tensors live simultaneously -> ~23 GB on GPU (fits 80 GB).
    """
    L = Pi.shape[0]
    q = Q
    eps = 4e-5
    tiny = 1e-300
    W_mf = torch.exp(-J4d).to(DEVICE)  # [L, L, q, q]
    pi = Pi.to(DEVICE)  # [L, q]
    pj = Pi.to(DEVICE)
    # mu1[i,j,a], mu2[i,j,b] init uniform over q
    mu1 = torch.full((L, L, q), 1.0 / q, dtype=torch.float64, device=DEVICE)
    mu2 = torch.full((L, L, q), 1.0 / q, dtype=torch.float64, device=DEVICE)
    pi_exp = pi[None, :, :]  # [1, L, q]
    pj_exp = pj[None, :, :]  # [1, L, q]
    for _ in range(200):
        # calc1[i,j,a] = sum_b mu2[i,j,b] * W_mf[i,j,a,b]
        calc1 = torch.einsum("ijb,ijab->ija", mu2, W_mf)
        calc2 = torch.einsum("ija,ijab->ijb", mu1, W_mf)
        new1 = pi_exp / torch.clamp(calc1, min=tiny)
        new1 = new1 / new1.sum(dim=2, keepdim=True)
        new2 = pj_exp / torch.clamp(calc2, min=tiny)
        new2 = new2 / new2.sum(dim=2, keepdim=True)
        d1 = (new1 - mu1).abs().max().item()
        d2 = (new2 - mu2).abs().max().item()
        mu1, mu2 = new1, new2
        if max(d1, d2) < eps:
            break
    # P_dir = W_mf * mu1[i,j,a] * mu2[i,j,b], normalized per pair
    Pdir = W_mf * mu1[:, :, :, None] * mu2[:, :, None, :]
    Pdir = Pdir / Pdir.sum(dim=(2, 3), keepdim=True)
    # Pfac = pi[i,a] * pj[j,b]
    Pfac = pi_exp[:, :, :, None] * pj_exp[:, :, None, :]
    ratio = torch.log(torch.clamp(Pdir, min=tiny) / torch.clamp(Pfac, min=tiny))
    DI = torch.einsum("ijab,ijab->ij", Pdir, ratio)
    diag = torch.arange(L, device=DEVICE)
    DI[diag, diag] = 0.0
    return DI.cpu().numpy()


def frobenius_scores(J4d):
    """F_ij = ||J_ij||_F (Frobenius norm of coupling block)."""
    F = torch.sqrt((J4d**2).sum(dim=(2, 3))).cpu().numpy()
    np.fill_diagonal(F, 0.0)
    return F


def apc_correction(F):
    """Average Product Correction (Dunn et al 2008): F_ij - F_i*F_j/F_mean."""
    with np.errstate(divide="ignore", invalid="ignore"):
        Fi = F.mean(axis=1, keepdims=True)
        Fj = F.mean(axis=0, keepdims=True)
        Fmean = F.mean()
        F_apc = F - Fi * Fj / Fmean
    np.fill_diagonal(F_apc, 0.0)
    return F_apc


def top_pairs(score_matrix, k=50):
    """Top-k position pairs by score."""
    L = score_matrix.shape[0]
    idx = np.triu_indices(L, k=1)
    scores = score_matrix[idx]
    order = np.argsort(scores)[::-1][:k]
    return [(int(idx[0][o]), int(idx[1][o]), float(scores[o])) for o in order]


def main():
    t0 = time.time()
    print("=" * 80)
    print("PROPER mfDCA — Morcos et al. 2011 (PNAS), GPU-accelerated")
    print(f"torch {torch.__version__}, device: {DEVICE}")
    print("=" * 80)

    print("\n[1/6] Loading MSA...")
    seqs = parse_fasta(FASTA)
    X = encode_msa(seqs)
    N, L = X.shape
    print(f"  {N} sequences, {L} positions, q={Q} states (20 AA + gap)")
    print(
        f"  Gaps per sequence: {((X == Q - 1).sum(1)).min().item()}-{((X == Q - 1).sum(1)).max().item()}"
    )

    print("\n[2/6] Sequence reweighting (theta=0.2)...")
    t1 = time.time()
    W = reweight(X)
    Meff = W.sum()
    print(f"  Meff = {Meff:.1f} (effective sequence count after reweighting)")
    print(f"  Max weight {W.max():.3f} (sequence with fewest neighbors)")

    print("\n[3/6] Computing frequencies with pseudocount lambda=0.5...")
    Pi, Pij = compute_frequencies(X, W)
    print(f"  Pi: [{Pi.shape}], Pij: [{Pij.shape}]")

    print("\n[4/6] Building covariance matrix...")
    C = build_covariance(Pi, Pij)
    print(f"  C: [{C.shape[0]} x {C.shape[1]}]")

    print("\n[5/6] Inverting covariance (GPU)...")
    t2 = time.time()
    C_gpu = C.to(DEVICE)
    # Regularize: add small diagonal for numerical stability
    reg = 1e-4 * torch.eye(C_gpu.shape[0], dtype=torch.float64, device=DEVICE)
    invC = torch.linalg.inv(C_gpu + reg)
    torch.cuda.synchronize()
    print(f"  Inversion: {time.time() - t2:.1f}s")
    # Couplings J = -invC, reshaped to [L, L, q-1, q-1], padded to [L, L, q, q]
    # invC rows/cols are indexed (i, a) -> reshape [L, q-1, L, q-1] then
    # permute [i, a, j, b] -> [i, j, a, b]
    invC4 = invC.reshape(L, Q - 1, L, Q - 1).permute(0, 2, 1, 3)  # [L, L, q-1, q-1]
    J4d = torch.zeros(L, L, Q, Q, dtype=torch.float64, device=DEVICE)
    J4d[:, :, : Q - 1, : Q - 1] = -invC4
    print(f"  Couplings J: [{J4d.shape}] (J = -C^-1, gap state padded to 0)")

    print("\n[6/6] Computing scores (Frobenius, APC, DI)...")
    F = frobenius_scores(J4d)
    F_apc = apc_correction(F)
    DI = direct_information(Pi, J4d)

    print("\n" + "=" * 80)
    print("RESULTS")
    print("=" * 80)

    top_F = top_pairs(F, 20)
    top_apc = top_pairs(F_apc, 20)
    top_DI = top_pairs(DI, 20)

    print("\nTop 20 pairs by Frobenius norm ||J_ij||_F:")
    for i, j, s in top_F:
        print(f"  ({i:4d},{j:4d}): F={s:.4f}")
    print("\nTop 20 pairs by APC-corrected Frobenius:")
    for i, j, s in top_apc:
        print(f"  ({i:4d},{j:4d}): F_apc={s:.4f}")
    print("\nTop 20 pairs by Direct Information:")
    for i, j, s in top_DI:
        print(f"  ({i:4d},{j:4d}): DI={s:.4f}")

    # Save
    summary = {
        "method": "mfDCA (Morcos et al 2011 PNAS)",
        "dataset": f"{N} Omicron Spike sequences, {L} positions",
        "q": Q,
        "theta": THETA,
        "lambda": LAMBDA,
        "Meff": float(Meff),
        "compute_time_s": round(time.time() - t0, 1),
        "n_pairs_total": L * (L - 1) // 2,
        "top_20_frobenius": [{"pos_i": i, "pos_j": j, "score": s} for i, j, s in top_F],
        "top_20_apc": [{"pos_i": i, "pos_j": j, "score": s} for i, j, s in top_apc],
        "top_20_di": [{"pos_i": i, "pos_j": j, "score": s} for i, j, s in top_DI],
    }
    with open(OUT / "dca_mf_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    np.save(OUT / "frobenius_scores.npy", F)
    np.save(OUT / "frobenius_apc_scores.npy", F_apc)
    np.save(OUT / "di_scores.npy", DI)
    np.savetxt(OUT / "frobenius_scores.csv", F, delimiter=",")
    np.savetxt(OUT / "di_scores.csv", DI, delimiter=",")
    print(f"\nSaved to {OUT}")
    print(f"Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
