#!/usr/bin/env python3
"""Exact prime implicants for the master Boolean function (corrections A3/A4).

master_boolean.py minimises each co-evolving pair's padded 32x32 mutation K-map
with `kmap_sbm.analysis.prime_implicants.boolean_minimize_kmap`.  Before the
kmap-sbm-validation fix of 2026-11 that function (a) indexed implicant
coverage and the on-set in two different minterm orders, so its essentials
and cover were wrong -- the published cover reached 2 of the 36 observed
mutations, which is where "2 essential" came from -- and (b) merged only
don't-cares adjacent to an on-cell.

This script rebuilds exactly the same maps (same functions, same pairs) and
reports, per pair and in total:

  reference   the kmap_sbm minimiser as installed (exact after the fix)
  exact       an independent exact Quine-McCluskey
              (scripts/sequence_circuits.prime_implicants), cross-checked
              against brute-force enumeration of all prime cubes, with
              essential primes and a verified cover

under the two don't-care policies discussed in doc 12 §3:

  padding-DC  rows/cols 20-31 and the reference cell are don't-care
              (master_boolean.py's choice)
  padding-OFF only the reference cell is don't-care; padding is off-set, so
              cubes stay inside the 20x20 residue block and decode to rules

The pre-fix published values are recorded under "published_before_fix".

Output: analysis/exact_master_boolean.json
Run:    python analysis/exact_master_boolean.py
"""

from __future__ import annotations

import json
import os
import sys
from itertools import product
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from kmap_paths import KMAP_SBM_SRC, NARY_SRC  # noqa: E402

sys.path.insert(0, str(NARY_SRC))
sys.path.insert(0, str(KMAP_SBM_SRC))

import numpy as np  # noqa: E402

import master_boolean as mb  # noqa: E402
import sequence_circuits as sc  # noqa: E402
from coevolution_shared import parse_fasta  # noqa: E402
from kmap_sbm.analysis.prime_implicants import boolean_minimize_kmap  # noqa: E402
from nkmap.encoding.bio_sequences import AMINO_HE_2012, Base20AminoEncoder  # noqa: E402

N_VARS = 10  # 5 row bits (MSB) + 5 column bits


def brute_force_primes(on: set[int], dc: set[int]) -> set[tuple[int, int]]:
    """All prime cubes by direct enumeration: a cube (val, mask) is an implicant
    if every cell it covers is on or don't-care and it covers an on-cell; it is
    prime if freeing any further fixed bit breaks that."""
    allowed = on | dc
    def cells(val, mask):
        free = [b for b in range(N_VARS) if (mask >> b) & 1]
        for bits in product((0, 1), repeat=len(free)):
            yield val | sum(1 << b for b, x in zip(free, bits) if x)
    def implicant(val, mask):
        cs = list(cells(val, mask))
        return all(c in allowed for c in cs) and any(c in on for c in cs)
    # grow from every on-cell; every prime contains an on-cell
    primes, seen = set(), set()
    stack = [(m, 0) for m in on]
    while stack:
        val, mask = stack.pop()
        if (val, mask) in seen:
            continue
        seen.add((val, mask))
        grown = False
        for b in range(N_VARS):
            if (mask >> b) & 1:
                continue
            nv, nm = val & ~(1 << b), mask | (1 << b)
            if implicant(nv, nm):
                grown = True
                stack.append((nv, nm))
        if not grown:
            primes.add((val, mask))
    return primes


def decode(pi: sc.Implicant) -> dict:
    """Row/column residue sets covered by a cube (inside the 20x20 block)."""
    rows = sorted({c >> 5 for c in sc.cube_minterms(pi) if (c >> 5) < 20 and (c & 31) < 20})
    cols = sorted({c & 31 for c in sc.cube_minterms(pi) if (c >> 5) < 20 and (c & 31) < 20})
    return {"aa_i": "".join(AMINO_HE_2012[r] for r in rows), "aa_j": "".join(AMINO_HE_2012[c] for c in cols),
            "free_bits": pi.n_free, "in_block": all((c >> 5) < 20 and (c & 31) < 20 for c in sc.cube_minterms(pi))}


def main() -> int:
    fasta = Path(os.environ.get("COEVO_FASTA") or REPO / "Spike_protein.aln-fasta")
    sequences = parse_fasta(fasta)
    encoder = Base20AminoEncoder(version=1)
    pos_arrays = [np.array([encoder.encode.get(a, 20) for a in s], dtype=np.int32) for _, s in sequences]
    n_all, max_pos = len(pos_arrays), len(sequences[0][1])
    variable = mb.find_variable_positions(pos_arrays, n_all, max_pos, threshold=0.3)
    pairs = mb.find_coevolutionary_pairs(pos_arrays, variable, n_all)[:15]   # as master_boolean.py
    aa_list = list(AMINO_HE_2012)

    rows, totals = [], {"reference_distinct_rules": 0, "reference_essential": 0, "reference_primes_covering_on": 0}
    ref_keys = set()
    for policy in ("padding-DC", "padding-OFF"):
        totals[policy] = {"exact_primes": 0, "exact_essential": 0, "exact_cover": 0, "bruteforce_agree": True}
    for pos_i, pos_j, mi, n_muts, ref_i, ref_j in pairs:
        kmap = mb.build_mutation_kmap(pos_arrays, pos_i, pos_j, ref_i, ref_j, n_all)
        flat = kmap.flatten()
        on = {int(c) for c in np.flatnonzero(flat == 1)}
        dc_all = {int(c) for c in np.flatnonzero(flat == -1)}
        ref_cell = {ref_i * 32 + ref_j} & dc_all

        # reference minimiser, decoded exactly as master_boolean.py does
        ref = boolean_minimize_kmap(flat.astype(int), algorithm="qm")
        ref_pis = mb.extract_prime_implicants(ref, aa_list)
        for p in ref_pis:
            ref_keys.add((pos_i, pos_j, p["aa_i"], p["aa_j"], p["is_essential"]))
        ref_on = [q for q in ref["prime_implicants"] if any(
            sum(int(v) << (N_VARS - 1 - k) for k, (v, m) in enumerate(zip(q["values"], q["mask"])) if not m)
            == (c & ~sum(1 << (N_VARS - 1 - k) for k, m in enumerate(q["mask"]) if m)) for c in on)]
        row = {"pair": [int(pos_i), int(pos_j)], "mi": float(mi), "n_on": len(on),
               "reference": {"raw_primes": len(ref["prime_implicants"]), "primes_covering_on": len(ref_on),
                             "raw_essential": len(ref["essential_prime_implicants"])}}
        for policy, dc in (("padding-DC", dc_all), ("padding-OFF", ref_cell)):
            # QM over on u dc also yields cubes that cover only don't-cares; they are
            # primes of the extended function but never needed in a cover, and the
            # standard definition of a prime implicant of f excludes them.
            raw = sc.prime_implicants(on, dc, N_VARS)
            pis = [p for p in raw if any(p.covers(m) for m in on)]
            cover, essential = sc._essential_and_greedy_cover(pis, sorted(on))
            complete, sound = sc.verify_cover(cover, on, allowed=on | dc)
            brute = brute_force_primes(on, dc)
            agree = brute == {(p.val, p.mask) for p in pis}
            row[policy] = {"exact_primes": len(pis), "dc_only_cubes_dropped": len(raw) - len(pis), "exact_essential": len(essential), "exact_cover": len(cover),
                           "cover_complete": complete, "cover_sound": sound, "bruteforce_agree": agree,
                           "essential_rules": [decode(p) for p in essential]}
            t = totals[policy]
            t["exact_primes"] += len(pis); t["exact_essential"] += len(essential); t["exact_cover"] += len(cover)
            t["bruteforce_agree"] &= agree and complete and sound
        totals["reference_primes_covering_on"] += len(ref_on)
        rows.append(row)

    distinct = {k[:4] for k in ref_keys}
    totals["reference_distinct_rules"] = len(distinct)
    totals["reference_essential"] = len({k[:4] for k in ref_keys if k[4]})
    out = {"pairs": rows, "totals": totals, "n_pairs": len(rows),
           "published_before_fix": {"distinct_rules": 36, "essential": 2,
                                    "on_cells_covered_by_published_cover": 2,
                                    "on_cells": sum(r["n_on"] for r in rows)}}
    (REPO / "analysis" / "exact_master_boolean.json").write_text(json.dumps(out, indent=2))

    print(f"{len(rows)} co-evolving pairs (MI > 0.1), padded 32x32 maps, 10 variables")
    print(f"reference minimiser (kmap_sbm, as installed): {totals['reference_distinct_rules']} distinct decoded rules, "
          f"{totals['reference_essential']} essential; {totals['reference_primes_covering_on']} cubes cover an on-cell")
    for policy in ("padding-DC", "padding-OFF"):
        t = totals[policy]
        print(f"exact QM, {policy:11s}: {t['exact_primes']:4d} prime implicants, {t['exact_essential']:3d} essential, "
              f"cover {t['exact_cover']:3d} cubes; brute-force agreement + verified cover: {t['bruteforce_agree']}")
    print(f"{'pair':>12s} {'on':>3s} | {'ref PIs':>7s} {'ref ess':>7s} | {'DC: PIs':>7s} {'ess':>4s} | {'OFF: PIs':>8s} {'ess':>4s}")
    for r in rows:
        print(f"{str(tuple(r['pair'])):>12s} {r['n_on']:3d} | {r['reference']['primes_covering_on']:7d} {r['reference']['raw_essential']:7d} | "
              f"{r['padding-DC']['exact_primes']:7d} {r['padding-DC']['exact_essential']:4d} | "
              f"{r['padding-OFF']['exact_primes']:8d} {r['padding-OFF']['exact_essential']:4d}")
    return 0 if all(totals[p]["bruteforce_agree"] for p in ("padding-DC", "padding-OFF")) else 1


if __name__ == "__main__":
    sys.exit(main())
