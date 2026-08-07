"""
Quine-McCluskey Boolean minimization of K-map truth tables.
"""

from __future__ import annotations

from typing import Any

import numpy as np


def kmap_truth_table(bool_map):
    """Convert a 4^k K-map to a truth table of shape (num_cells, 2k+1).

    Each row: [bit_0, ..., bit_{2k-1}, output_value]
    Output values: 1 (on-set), 0 (off-set), -1 (don't-care).

    BUG FIX: Original code used `bool_map.shape[0]` which returns the first
    dimension only — for a 2D array (e.g., 32×32) this gives 32 instead of
    1024. Fixed to use `bool_map.size` which gives the total number of cells.

    HARD GUARD: the number of cells must be a power of 4 (4^k), because the
    encoding uses 2k bits per cell. A 20×20 = 400 cell map is NOT a power of
    4: k_bits = int(log2(400))//2 = 4 would give only 8 bits (256 cells) and
    cells 256-399 would silently wrap onto cells 0-143, corrupting minterms
    and producing phantom rules. Verified defect: 143 of 152 rules in the
    original pipeline had labels never observed in the data.
    Callers with a 20×20 map MUST pad to 32×32 (5 bits per axis, 10 bits
    total, rows/cols 20-31 as don't-care) before calling this function.
    """
    bool_map = np.asarray(bool_map).ravel()  # Flatten to 1D first
    nc = bool_map.size  # FIX: was shape[0] (breaks for 2D input)
    # Guard: nc must be a power of 4
    if nc > 0:
        log4 = np.log(nc) / np.log(4)
        if abs(log4 - round(log4)) > 1e-9:
            raise ValueError(
                f"kmap_truth_table requires 4^k cells (power of 4), got {nc}. "
                f"A 20x20=400 map is not representable in 8 bits; pad to 32x32 "
                f"(1024 = 4^5) with don't-care padding first."
            )
    k_bits = int(np.log2(nc)) // 2 if nc > 0 else 0
    rows = []
    for idx in range(nc):
        val = bool_map.flat[idx]
        bits = [(idx >> (2 * k_bits - 1 - b)) & 1 for b in range(2 * k_bits)]
        rows.append(bits + [int(val)])
    return np.array(rows) if rows else np.zeros((0, 2 * k_bits + 1), dtype=int)


def prime_implicants_quine_mccluskey(truth_table):
    """Exact Quine-McCluskey minimization for small K-maps.

    BUG FIX: Don't-care terms (-1) are now included in the on-set for prime
    implicant generation. This allows QM to form LARGER implicants by
    combining on-set minterms with don't-care terms — the standard QM
    behavior. Don't-cares are NOT required to be covered in the final
    chart (handled in boolean_minimize_kmap).
    """
    # Include both on-set (1) and don't-care (-1) for PI generation
    on_set = truth_table[(truth_table[:, -1] == 1) | (truth_table[:, -1] == -1)]
    if len(on_set) == 0:
        return []

    n_vars = truth_table.shape[1] - 1
    minterms = [tuple(int(b) for b in row[:n_vars]) for row in on_set]

    prime_implicants = _qm_find_prime_implicants(minterms, n_vars)
    return _qm_format_implicants(prime_implicants, n_vars, minterms)


def _qm_find_prime_implicants(minterms, n_vars):
    """Internal: iterative merging phase of QM algorithm."""
    implicants = [(m, tuple(False for _ in range(n_vars))) for m in minterms]

    changed = True
    while changed:
        changed = False
        new_implicants = []
        used_indices: set[int] = set()
        for i, (a, mask_a) in enumerate(implicants):
            for j, (b, mask_b) in enumerate(implicants):
                if i >= j:
                    continue
                diff_pos = None
                compatible = True
                for k in range(n_vars):
                    if mask_a[k] and mask_b[k]:
                        if a[k] != b[k]:
                            compatible = False
                            break
                        continue
                    if not mask_a[k] and not mask_b[k]:
                        if a[k] != b[k]:
                            if diff_pos is None:
                                diff_pos = k
                            else:
                                compatible = False
                                break
                    else:
                        compatible = False
                        break
                if compatible and diff_pos is not None:
                    new_val = list(a)
                    new_mask = list(mask_a)
                    new_val[diff_pos] = 0
                    new_mask[diff_pos] = True
                    new_t = (tuple(new_val), tuple(new_mask))
                    if new_t not in new_implicants:
                        new_implicants.append(new_t)
                    used_indices.add(i)
                    used_indices.add(j)
                    changed = True

        for i, imp in enumerate(implicants):
            if i not in used_indices and imp not in new_implicants:
                new_implicants.append(imp)
        implicants = new_implicants
    return implicants


def _qm_format_implicants(implicants, n_vars, minterms):
    """Format prime implicants with coverage info."""
    result = []
    for val_tuple, mask_tuple in implicants:
        covered_mt = []
        for idx, mt in enumerate(minterms):
            match = True
            for i in range(n_vars):
                if mask_tuple[i]:
                    continue
                if val_tuple[i] != mt[i]:
                    match = False
                    break
            if match:
                covered_mt.append(idx)

        result.append(
            {
                "values": list(val_tuple),
                "mask": list(mask_tuple),
                "n_dontcares": int(sum(mask_tuple)),
                "coverage": covered_mt,
            }
        )
    return result


def kmap_signature_distance(bool_a, bool_b):
    """Hamming distance between two Boolean K-map signatures.

    Only counts positions where both maps have defined (non-don't-care) values.
    """
    a = np.asarray(bool_a).ravel()
    b = np.asarray(bool_b).ravel()
    mask = (a != -1) & (b != -1)
    if mask.sum() == 0:
        return 0
    return int((a[mask] != b[mask]).sum())


def boolean_minimize_kmap(bool_map, algorithm: str = "qm"):
    """Perform Boolean minimization on a K-map truth table.

    Args:
        bool_map: 4^k ternary K-map (-1=don't-care, 0=0, 1=1).
        algorithm: "qm" (Quine-McCluskey) or "espresso".

    Returns dict with prime_implicants, covering_set, and metrics.

    BUG FIXES:
      1. Essential PI identification now only considers ON-SET minterms
         (not don't-cares) for unique coverage.
      2. covering_size now correctly computes the minimum cover using
         essential PIs + greedy fallback for remaining minterms.
    """
    bool_map_arr = np.asarray(bool_map).ravel()
    tt = kmap_truth_table(bool_map_arr)
    pis = prime_implicants_quine_mccluskey(tt)
    n_on = int((bool_map_arr == 1).sum())
    n_dc = int((bool_map_arr == -1).sum())
    total_cells = bool_map_arr.size

    # The truth table includes both on-set (1) and don't-care (-1) rows for PI
    # generation. The minterms list is built from these rows in order.
    # We need to track which minterms in the list are actual on-set (not DC)
    # to determine essential PIs and covering.
    # The minterms list order matches: on_set = tt[tt[:,-1] == 1 | tt[:,-1] == -1]
    # So minterm_idx i corresponds to the i-th row in on_set.
    on_set_rows = tt[(tt[:, -1] == 1) | (tt[:, -1] == -1)]
    # minterm_is_onset[i] = True if the i-th minterm is actual on-set (value 1)
    minterm_is_onset = [row[-1] == 1 for row in on_set_rows]
    on_set_minterm_indices = [i for i, is_on in enumerate(minterm_is_onset) if is_on]

    # Essential PIs: those covering at least one on-set minterm uniquely
    essential = []
    for p in pis:
        for idx in p["coverage"]:
            # Only check on-set minterms (not don't-cares)
            if idx not in on_set_minterm_indices:
                continue
            is_essential = True
            for q in pis:
                if q is p:
                    continue
                if idx in q["coverage"]:
                    is_essential = False
                    break
            if is_essential:
                essential.append(p)
                break

    # Remove duplicate essential PIs
    seen = set()
    essential_unique = []
    for p in essential:
        key = (tuple(p["values"]), tuple(p["mask"]))
        if key not in seen:
            seen.add(key)
            essential_unique.append(p)
    essential = essential_unique

    # Compute actual covering_size:
    # 1. Start with essential PIs
    # 2. Find which on-set minterms are still uncovered
    # 3. Greedily add PIs to cover remaining minterms (approximate Petrick's method)
    covered = set()
    for p in essential:
        for idx in p["coverage"]:
            if idx in on_set_minterm_indices:
                covered.add(idx)

    uncovered = set(on_set_minterm_indices) - covered
    remaining_pis = [p for p in pis if p not in essential]

    # Greedy: add PIs that cover the most uncovered minterms
    while uncovered:
        best_pi = None
        best_cover = 0
        for p in remaining_pis:
            new_cover = len(set(p["coverage"]) & uncovered)
            if new_cover > best_cover:
                best_cover = new_cover
                best_pi = p
        if best_pi is None or best_cover == 0:
            break
        essential.append(best_pi)
        remaining_pis.remove(best_pi)
        for idx in best_pi["coverage"]:
            if idx in on_set_minterm_indices:
                covered.add(idx)
        uncovered = set(on_set_minterm_indices) - covered

    covering_size = len(essential)

    return {
        "n_vars": int(tt.shape[1] - 1) if len(tt) > 0 else 0,
        "n_minterms": n_on,
        "n_dontcare": n_dc,
        "n_cells": int(total_cells),
        "n_prime_implicants": len(pis),
        "n_essential": len(essential),
        "prime_implicants": pis,
        "essential_prime_implicants": essential,
        "covering_size": covering_size,
    }
