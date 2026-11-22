#!/usr/bin/env python3
"""Top co-evolving pairs vs native contacts on PSICOV150, with correct mapping.

scripts/validate_contacts.py scored each PSICOV target's top co-evolving pairs
(master_boolean's top_co_evolving_pairs: mutation-only MI > 0.1, window 30,
entropy > 0.3, at most 20) after mapping alignment columns to PDB residues by
globally aligning the alignment's FIRST row to the PDB chain.  In the PSICOV
alignments the first row is a gapped homologue, not the target, so that
mapping was partly wrong; for PSICOV the correct map is the identity (column
c = residue c+1, verified for 149/150 targets in kmap_structure.py).

For every target this script recomputes the pairs exactly as master_boolean.py
does and reports, under both mappings:
  precision   fraction of the pairs that are native contacts (Cb < 8 A, sep >= 6)
  random      contact rate of random pairs with the same separations
              (draws always inside the chain)
  enrichment  precision / random

Output: analysis/psicov_top_pairs_check.json
Run:    python analysis/psicov_top_pairs_check.py
"""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]

import numpy as np  # noqa: E402

import kmap_structure as ks  # noqa: E402
import master_boolean as mb  # noqa: E402
from validate_contacts import contact_map, nw_align  # noqa: E402


def score(pairs, contacts, n_res, mapping, seed=42):
    mapped = [(mapping[i], mapping[j]) for i, j in pairs if i in mapping and j in mapping]
    if not mapped:
        return None
    mapped = [(min(a, b), max(a, b)) for a, b in mapped]
    prec = sum(p in contacts for p in mapped) / len(mapped)
    seps = [b - a for a, b in mapped]
    rng = random.Random(seed)
    hits = drawn = 0
    for _ in range(2000):
        s = rng.choice(seps)
        if s <= 0 or s >= n_res:
            continue
        a = rng.randrange(n_res - s)
        drawn += 1
        hits += (a, a + s) in contacts
    rand = hits / drawn if drawn else 0.0
    return {"n_mapped": len(mapped), "precision": prec, "random": rand,
            "enrichment": prec / rand if rand > 0 else None}


def main() -> int:
    rows = []
    for t in ks.list_targets():
        pdat = ks.load_protein(t)
        pa = [r.astype(np.int32) for r in pdat.dense]
        L = pdat.n_cols
        var = mb.find_variable_positions(pa, len(pa), L, threshold=0.3)
        pairs = [(int(p[0]), int(p[1])) for p in mb.find_coevolutionary_pairs(pa, var, len(pa))[:20]]
        seq_pdb, contacts = contact_map(str(ks.RAW / "pdb" / f"{t}.pdb"))
        identity = {c: c for c in range(min(L, len(seq_pdb)))}
        first = open(ks.RAW / "aln" / f"{t}.aln").readline().strip().upper()
        old = nw_align(first, seq_pdb)   # validate_contacts' former first-row mapping
        row = {"target": t, "n_pairs": len(pairs),
               "identity": score(pairs, contacts, len(seq_pdb), identity),
               "first_row_nw": score(pairs, contacts, len(seq_pdb), old)}
        rows.append(row)
        if len(rows) % 25 == 0:
            print(f"  {len(rows)}/150", flush=True)

    def agg(key):
        vals = [r[key] for r in rows if r[key] is not None]
        enr = [v["enrichment"] for v in vals if v["enrichment"] is not None]
        return {"n_targets": len(vals),
                "mean_precision": float(np.mean([v["precision"] for v in vals])),
                "median_precision": float(np.median([v["precision"] for v in vals])),
                "mean_enrichment": float(np.mean(enr)) if enr else None,
                "targets_enrichment_gt2": int(sum(e > 2 for e in enr)),
                "mean_pairs_mapped": float(np.mean([v["n_mapped"] for v in vals]))}

    out = {"identity_mapping": agg("identity"), "first_row_nw_mapping": agg("first_row_nw"),
           "per_target": rows}
    (REPO / "analysis" / "psicov_top_pairs_check.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "per_target"}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
