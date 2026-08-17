# 09 — Algorithmic Details (minute detail)

The full algorithmic specification of every computational step in the pipeline.

---

## 1. Encoding

### Base-20 (He 2012) encoding
`AMINO_HE_2012`: A=0, I=1, L=2, V=3, M=4, F=5, Y=6, W=7, E=8, D=9, Q=10, N=11,
H=12, K=13, R=14, S=15, T=16, C=17, P=18, G=19. Gap `-` and ambiguous `X` =
state 20 (21 states; excluded from all counts). Proven in `n-ary-kmap/` (115
Lean theorems).

### 5-bit Gray encoding
`g(i) = i ^ (i >> 1)` on the He-2012 code. Single bit-flips connect
physicochemical neighbours (7 groups). The 20 canonical residues occupy 20 of
the 32 codes of Q₅. Proven in `lean_proofs/AminoAcidEncoding.lean` (106
theorems). Used by `run_kmap_analysis.py` (H1 adjacency) and the K-map Boolean
literals (5 bits per axis).

## 2. Shannon entropy (vectorized)
$$H(p) = -\sum_{a=0}^{19} P(a|p)\log_2 P(a|p)$$
Implementation (`compute_entropy_vectorized`): dense [N×L] int32 matrix;
per-column `np.bincount(valid, minlength=20)`; $p = \text{cnt}/\text{total}$;
$H = -\sum p\log_2 p$ over $p > 0$. Gaps (state 20) excluded. O(N·L).

## 3. Mutual information (vectorized, single source of truth)
$$MI(i,j) = \sum_{a,b} P(a,b)\log_2\frac{P(a,b)}{P(a)P(b)}$$
Implementation (`mutual_information`): extract columns i, j; mask valid
(0–19); `pairs = codes_i*20 + codes_j`; `joint = np.bincount(pairs,
minlength=400).reshape(20,20)`; marginals by axis-sum; MI loop over 400 cells.
O(N + 400) per pair.

### Mutation-only MI
Excludes the reference (majority) pair from joint counts — isolates mutation-
driven covariation. Reference = most common residue (lowest-code tie-break,
unified CPU/GPU).

## 4. GPU MI matrix (`mi_matrix_gpu`)
- Dense [N×L] int32 tensor on CUDA.
- Per chunk of pairs: gather columns `ci = dense[:, ii]`, `cj = dense[:, jj]`
  → [N, P].
- `valid = (ci≥0)&(ci<20)&(cj≥0)&(cj<20)`; for mutation-only, also
  `& ~is_ref`.
- `flat = ci.clamp(0,19)*20 + cj.clamp(0,19)`, masked_fill invalid → 0.
- `joint = torch.zeros(P, 400).scatter_add_(1, flat.T, valid.float().T)` →
  [P, 20, 20].
- MI = $\sum p \log_2(p / (p_i \cdot p_j))$ per pair, vectorized.
- **Adaptive chunking**: `chunk = max(256, 512MB / (N*8))` so the flat int64
  buffer never exceeds 512 MB → 486k-sequence families never OOM.
- Speed: 813K pairs × 1,299 seqs in 1.6 s (~800× vs CPU).

## 5. K-map construction (per position pair)
**Frequency K-map:** $K[a,b] = \frac{1}{N}\sum_s \mathbf{1}[s_i=a \wedge s_j=b]$
**Mutation K-map** (ternary, for QM):
$$f(a,b) = \begin{cases} 1 & \text{observed mutation} \\ -1 & \text{reference (don't-care)} \\ 0 & \text{never observed} \end{cases}$$
**32×32 padding (FIX A2):** pad 20×20 → 32×32 (5 bits/axis, 10 bits total =
1,024 cells); rows/cols 20–31 = don't-care. `kmap_truth_table` raises on
non-power-of-4 input.

## 6. Quine–McCluskey with don't-cares (`prime_implicants.py`)
1. **Truth table**: 1,024 cells → 10-bit minterms; on-set (1) + don't-care (-1)
   included for PI generation.
2. **Merge phase** (iterative): group minterms by popcount; merge pairs
   differing in exactly one bit → larger implicant with that bit as don't-care.
   Repeat until no merges (prime implicants). O(m · 2^bits) worst case.
3. **Essential PI**: a PI covering ≥1 on-set minterm covered by no other PI.
4. **Cover**: essential PIs + greedy cover of remaining on-set minterms.
- Verified against brute-force reference (50 random trials).
- Result: **36 distinct PIs, 2 essential** (after dedup of cubes decoding to
  the same residue pair).

## 7. Coupling / constraint function
$$C(a,b) = J(a,b) = \ln\frac{P(a,b)}{P(a)P(b)}, \quad \hat{y} = \sigma(C) = \frac{1}{1+e^{-C}}$$
Sign: C > 0 = co-evolutionary (more common than expected). (FIX: the original
used -ln, giving the opposite sign — corrected to +ln to match DCA convention.)

## 8. Perplexity ratio (vectorized)
$$r(i,j) = \frac{PP(j)}{\text{mean}_a\, PP(j \mid i=a)}, \quad PP = 2^{H}$$
Implementation (`perplexity_ratio` / `_ratio_dense`): joint 20×20 via bincount;
conditional entropy $H(j|i=a)$ over residues with $n \geq 5$; $PP(j|i=a) =
2^{H(j|i=a)}$; ratio = marginal PP / mean conditional PP. Vectorized (numpy
bincount) — ~100× faster than the original Counter version (verified equal).

## 9. mfDCA (`dca_mf_analysis.py`)
1. Compute the covariance matrix $C_{ij}(a,b) = P(a,b) - P_i(a)P_j(b)$.
2. Mean-field: $J = C^{-1}$ (regularized, $\lambda = 0.5$).
3. **APC**: $J^{\text{APC}}_{ij}(a,b) = J_{ij}(a,b) - \frac{\text{row}_i \cdot
   \text{col}_j^T}{\text{grand mean}}$ (removes the rank-1 phylogenetic bias).
4. **Direct information**: $DI = \sum_{a,b} P_{\text{direct}}(a,b) \cdot
   J^{\text{APC}}(a,b)^2$.
5. Frobenius norm, APC, DI reported. ρ(DI, MI) = 0.06 (uncorrelated).

## 10. H1 Gray adjacency (`h1_adjacency_gpu`)
- For consecutive residues: $g_a = a \oplus (a \gg 1)$, $g_b = b \oplus (b
  \gg 1)$; $d = \text{popcount}(g_a \oplus g_b)$.
- GPU popcount via bit tricks: `x - ((x>>1)&0x5555...)` → `(x&0x3333...) +
  ((x>>2)&0x3333...)` → `(x + (x>>4)) & 0x0F0F...` → `(x * 0x01010101) >> 24`.
- Observed H=1 fraction vs expected 5/31 = 0.1613. Result: 0.2167 → **1.34×**.
- 1,644,588 consecutive pairs on 1,299 sequences.

## 11. Complexity table

| Step | Complexity | Optimization | Scale achieved |
|---|---|---|---|
| Entropy (all positions) | O(N·L) | vectorized bincount; GPU < 0.1 s | 1,299 × 1,276 |
| MI one pair | O(N + 400) | numpy bincount | — |
| Full MI matrix (window 30) | O(L·N) | GPU scatter_add; adaptive chunk | 813K pairs in 1.6 s |
| QM minimization | O(m·2^bits) | popcount-grouped merge; 32×32 padded | 1,024 cells, 36 PIs |
| LOO-CV constraint | O(10·N²) | CPU-only (forced off CUDA) | infeasible > 222k seqs |
| H3 pairwise K-map distance | O(n²·1024) | GPU norm-trick, chunked stats | 152k seqs, 11.6B pairs |
| Clustering | O(n²) | capped 2,000 seqs (logged) | — |

## 12. The two corrected defects (algorithmic)

| Defect | Root cause | Fix |
|---|---|---|
| A1 gap-stripping | `clean = "".join(aa for aa in seq if aa in encode)` deleted gaps → column j mixed raw positions | keep alignment; gap = state 20; `aligned=True` default |
| A2 8-bit QM wrap | 400-cell map → `k_bits = log2(400)//2 = 4` → 8 bits → cells 256–399 wrapped onto 0–143 | pad to 32×32 (5 bits/axis); `kmap_truth_table` raises on non-power-of-4 |
